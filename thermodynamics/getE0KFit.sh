#!/bin/bash

#-----set default parameters and paths--------------------------------------------------------------------------------
inpStrAll="sc fcc bcc hcp none"; fitTypeAll="Vinet Murn Birch"; fitDef="Vinet"; inpDef="energies_Tto0K";
unitsDef1="Ang Ang^3 Bohr Bohr^3"; unitsDef2="eV meV Hartree"; unitsDef="Ang:eV"; Bmin="5"; Bdef="50";
Bmax="4500"; BderMin=".1"; BderMax="50"; BderDef="4"; minNrMeshPoints="5"; meshDef="100";
#---------------------------------------------------------------------------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -s -i -u -f -n -V -a -B -Bd -k"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage "$script             without Options: create EOS"
  printOptions "-i inpFile   input energy vs volume file-----def: \033[1m$inpDef\033[0m" \
               "-u units     units of inpFile----------------def: \033[1m$unitsDef\033[0m; supported: check -help" \
               "-s strType   structure type------------------def: \033[1mfrom path\033[0m; supported: $inpStrAll" \
               "-f fitType   type of EOS---------------------def: \033[1m$fitDef\033[0m; supported: $fitTypeAll" \
               "-n nPoints   nr of volume points in output---def: \033[1m$meshDef\033[0m"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   fit of energies vs volume curve to equation-of-state (EOS)"
  echo2 "   EOS types implemented for fitting: $fitTypeAll"
  echo2 "   energy vs volume curve is read from input file (min. $minNrMeshPoints volume points expected)"
  echo2 "   the units of the input file can be:" \
        "   1.column: $unitsDef1" \
        "   2.column: $unitsDef2"
  echo2 "   if 1.column is length then it is assumed to be the lattice constant" \
        "   and converted to a volume per atom according to structure type"
  echo2 "   if strType (-s option) is not given then it is tried to extract it" \
        "   from the path where $script is executed"
  echo2 "   following additional options can be useful in case of convergence problems" \
        "   or when debugging:" \
       "      -V startV    starting volume in Ang^3 used in fit (default: Vmin+(Vmax-Vmin)/2)" \
       "      -a starta    instead of startV a starting lattice constant (Ang) can be specified" \
       "      -B startB    starting bulk modulus in GPa used in fit (default: $Bdef GPa)" \
       "      -Bd startBd  starting B derivative used in fit (default: $BderDef)" \
       "      -k           keep mathematica output for debugging"

  exit
fi

# mathematica kernel if needed
checkAndSetMath

# get structure type and check or try to extract from path
inpStr="none"
sOp=`getOption -s`;
if [ $sOp == True ]; then
  inpStr=`getValue -s`; if [ -z "$inpStr" ]; then error "no value to -s option"; fi
  c=`echo $inpStrAll | xargs -n1 | grep "\<$inpStr\>"`
  if [ "$c" != $inpStr ]; then error "input structure wrong"; fi
else
  for i in $inpStrAll; do
    inpStr=`pwd | sed 's/.*\('"$i"'\).*/\1/'`
    if [ "$inpStr" == $i ]; then
      echo; echo;
      echo "----------------------------------------------------------------------------------------"
      echo -e "            WARN: structure type \033[1m\033[31m$inpStr\033[0m from path will be used"
      echo "----------------------------------------------------------------------------------------"
      break;
    fi
  done
fi

# if strType=hcp check if cBya file provided
cBya=1.
if [ $inpStr = hcp ]; then
  if [ ! -e cBya ]; then error "strType = hcp but no cBya file provided"; fi
  cBya=`awk 'NR==1{print $1}' cBya`
  c=`checkReal $cBya`; if [ "$c" != ok ]; then error "cBya value from cBya file wrong"; fi
fi

# get input energy file and check
file=`getOption -i`;  if [ $file  == True ]; then file=`getValue -i`; else file=$inpDef; fi
check $file

# get units of input energy file and check
units=`getOption -u`; if [ $units == True ]; then units=`getValue -u`; else units=$unitsDef; fi
volumeUnit=`echo $units | sed 's/\(.*\):\(.*\)/\1/'`
energyUnit=`echo $units | sed 's/\(.*\):\(.*\)/\2/'`
c1=`echo $unitsDef1 | xargs -n1 | awk '$0=="'$volumeUnit'"{print $0}'`
c2=`echo $unitsDef2 | xargs -n1 | awk '$0=="'$energyUnit'"{print $0}'`
if [ "$c1" != "$volumeUnit" -o "$c2" != "$energyUnit" ]; then error "units wrong"; fi

# get type of EOS fit and check
type=`getOption -f`;  if [ $type  == True ]; then type=`getValue -f`; else type=$fitDef; fi
c=`echo $fitTypeAll | xargs -n1 | grep "\<$type\>"`
if [ "$c" != $type ]; then error "EOS fitType wrong"; fi

# check if structure type is set if we have aLat in input energy file
if [ $volumeUnit == Ang -o $volumeUnit == Bohr ]; then
  if [ $inpStr == "none" ]; then error "aLat provided in $file but no structure type available"; fi
fi

# get nr of mesh points for output
mesh=`getOption -n`; if [ $mesh == True ]; then mesh=`getValue -n`; else mesh=$meshDef; fi

# get starting values for fitting
opV=`getOption -V`
opa=`getOption -a`
if [ $opV == True -a $opa == True ]; then error "options -V and -a incompatibel"; fi
if [ $opV == True -a $inpStr != none ]; then error "option -V only compatibel with structure type = none"; fi
if [ $opa == True -a $inpStr == none ]; then error "option -a incompatibel with structure type = none"; fi
Vdef=-1.;
if [ $opV == True ]; then Vdef=`getValue -V`; fi
if [ $opa == True ]; then Vdef=`getValue -a`; fi
c=`checkReal $Vdef`; if [ "$c" != ok ]; then error "value given to -V or -a option wrong"; fi

op=`getOption -B`;  if [ $op == True ]; then Bdef=`getValue -B`; fi
if [ -z "$Bdef" ]; then error "no value to -B option"; fi
op=`getOption -Bd`; if [ $op == True ]; then BderDef=`getValue -Bd`; fi
if [ -z "$BderDef" ]; then error "no value to -Bd option"; fi

# check if input energy file has enough mesh points
n=`awk 'END{print NR}' $file`; if [ $n -lt $minNrMeshPoints ]; then error "not enough mesh points in $file"; fi

# print some relevant information
if [ $inpStr == 0 ]; then inpStr="none provided"; fi
echo; echo; echo -e " fitType: \033[1m$type\033[0m  inputFile: \033[1m$file\033[0m  units: \033[1m$units\033[0m  structure: \033[1m$inpStr\033[0m"; echo

# do the fitting
echo " ... fitting ..."; echo; rm -f _tmp_math
$math >> _tmp_math << EOF
<<$path/mathematica/ALL.math;
fitToEOS["$volumeUnit","$energyUnit","$file","$inpStr",$cBya,$type,$Vdef*1.,$Bdef*1.,$Bmin*1.,$Bmax*1.,$BderDef*1.,$BderMin*1.,$BderMax*1.,$mesh]
EOF

# format ouput
E0=`awk '$1=="'E0InmeV'"{printf("%d",$2)}' _tmp_math`
aLat=`awk '$1=="'aLatInAng'"{printf("%.3f",$2)}' _tmp_math`
V0=`awk '$1=="'V0InAng3'"{printf("%.3f",$2)}' _tmp_math`
B0=`awk '$1=="'B0GPa'"{printf("%.1f",$2)}' _tmp_math`
B0d=`awk '$1=="'B0der'"{printf("%.2f",$2)}' _tmp_math`
delta=`awk '$1=="'maxDeltaInmeV'"{printf("%.2f",$2)}' _tmp_math`

# print output
echo -e " E(meV) \033[1m\033[32m$E0\033[0m  a(Ang) \033[1m\033[32m$aLat\033[0m  V(Ang^3) \033[1m\033[32m$V0\033[0m " \
        " B(GPa) \033[1m\033[32m$B0\033[0m  Bder \033[1m\033[32m$B0d\033[0m  delta(meV) \033[1m\033[31m$delta\033[0m"
echo;

# check if mathematica output should be removed or not
keep=`getOption -k`; if [ $keep != True ]; then rm _tmp_math; fi

