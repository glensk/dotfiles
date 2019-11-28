#!/bin/bash

#-----set parameters and paths------------------------
inputFileDef=EVinet; TstartDef=1; TstepDef=1;
Vs=0.99; nVDef=6; VstepDef=0.2 # angstrom^3
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -e -i -c -Tstart -Tstep -Tmax -Vstart -Vstep -nV -Vfix -f"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printNeeded  "-e elem/mass      element(s) or (effective) mass"
  printOptions "-i inpT0K         input T=0K parametrization (default: $inputFileDef)"   \
               "-c concs          if more elements then concentrations need to be given" \
               "-Tstart Tstart    start T for Fqh mesh (def: Tstart=$TstartDef)"         \
               "-Tstep  Tstep     T step for Fqh mesh (def: Tstep=$TstepDef)"            \
               "-Tmax   Tmax      max T for Fqh mesh (def: Tmax=melting T)"              \
               "-Vstart Vstart    start V for Fqh files (def: Vstart=$Vs*V0)"            \
               "-Vstep  Vstep     V step for Fqh files (def: Vstep=$VstepDef)"           \
               "-nV     nrV       nr of Fqh files along V (def: nrV=$nVDef)"             \
               "-Vfix   vol       produce Fqh file additionally at this volume"          \
               "-f                force overwriting previous output"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   calculate vibrational free energies in the Debye-Grueneisen approximation as a function" \
        "   of temperatures and at different volumes according to Moruzzi, PRB 37 (1988) 790"
  echo2 "   element(s) or effective mass (a.u.) need to be provided with -e option; if more elements are" \
        "   given then the concentration needs to be provided with the -c option; if the mass is given"   \
        "   also the maximum temperature (K) needs to be provided with -Tmax option" \
        "   e.g.:" \
        "    getDebye.sh -e Fe Co Cr Mn Ni -c .2 .2 .2 .2 .2" \
        "   or" \
        "    getDebye.sh -e 50 -Tmax 2000"
  echo2 "   EVinet or EMurn or EBirch must exist with the T=0K parametrization; can be obtained with" \
        "   one of the createFolders_murn*.sh scripts"
  echo2 "   Fqhs are created for the low T gamma=1 and for high T gamma=2/3 (see Moruzzi Ref)"
  echo2 "   if element(s) is(are) given and no -Tmax option then the max temperature is taken to be Tmelt" \
        "   which for more elements is linearly mixed from the consitutents with the concentrations"
  echo2 "   the volumes are automatically determined from Vstart=$Vs*V0, Vstep=$VstepDef, and nrV=$nVDef" \
        "   which should reasonably well span the relevant thermal expansion region; the range can be" \
        "   adjusted alternatively with the options -Vstart, -Vstep, -nV"
  echo2 "   at the end an adjusted parameters.math input file for getFqhFit.sh is created so that" \
        "   getFqhFit.sh can be started directly (see script sequence below)"
  echo2 "   typical script sequence (for emto):"                                                     \
        "    createFolders_murn_emto.sh -a           (adjust parameters.dat)"                        \
        "    createFolders_murn_emto.sh"                                                             \
        "    submit.sh run.emto.ser.cluster.workdir  (or the local version)"                         \
        "    collectEMTO.sh                          (when finished and inside NxNxNkp_XXeV folder)" \
        "    extractEnergies_emto.sh"                                                                \
        "    getE0KFit.sh -i energies_Tto0K_GGA      (cmpc01)"                                       \
        "    \033[1mgetDebye.sh\033[0m -e elem/mass                (creates Debye Fqhs)"             \
        "    getFqhFit.sh                            (cmpc01)"                                       \
        "    getThermodynamics.sh                    (collect EVinet,Fqh_fit in separate folder)"    \
        "                                            (to add Fel check createFolders_elec_emto.sh -help)"
  exit
fi

# check if previous output exists
overwrite=`getOption -f`
if [ -e log_getDebye -a $overwrite != True ]; then error "previous output exists, use -f to overwrite"; fi

# log file
echo > log_getDebye
echo " $script $options" >> log_getDebye
echo >> log_getDebye

# element(s) or mass
el=`getOption -e`; if [ $el != True ]; then error "elem or mass missing (-e elem or -e mass)"; fi
el=`getValue -e`; if [ -z "$el" ]; then error "value for -e option missing"; fi
c=`echo $el | awk '/^[0-9.]+$/{print "mass";exit};{print "elem"}'`
if [ $c == elem ]; then
  nel=`echo $el | awk '{print NF}'`
  if [ $nel -gt 1 ]; then
    # for multiple elements we run some checks concerning -e and -c option
    # we need -c option with concentrations for more than one element
    conc=`getOption -c`
    if [ $conc != True ]; then error "more elements in -e option provided but no -c option with concentrations given"; fi
    conc=`getValue -c`
    # check if number of concentrations equal to number of elements and if they add up to one
    nconc=`echo $conc | awk '{print NF}'`
    if [ $nconc -ne $nel ]; then error "number of concentrations in -c option is not equal to number of elements in -e option"; fi
    totConc=`echo $conc | awk '{c=0; for (i=1;i<=NF;i++) c=c+$i; if ((c-1)^2>0.0001) print "wrong"}'`
    if [ "$totConc" == wrong ]; then error "concentrations in -c option do not add up to one"; fi
    echo " elements:      $el" >> log_getDebye
    echo " concentration: $conc" >> log_getDebye
  else
    # for one element we set the concentration simply to one
    conc="1"
    echo " element:      $el" >> log_getDebye
  fi

  mass=0
  Tmax=0
  # now we loop over the elements (in $el) and average the masses with the concentrations (in $conc)
  for (( i=1; i<=$nel; i++ )) do
    e=`echo $el   | awk '{print $'$i'}'`
    x=`echo $conc | awk '{print $'$i'}'`
    m=`$path/getAtomicMass.sh $e`
    mass=`echo $mass $m $x | awk '{printf "%.7f",$1+$2*$3}'`
    # determine also Tmelt if no -Tmax option
    if [ `getOption -Tmax` != True ]; then
      melt=`$path/getMeltingPoint.sh $e`
      Tmax=`echo $Tmax $melt $x | awk '{printf "%.7f",$1+$2*$3}'`
    fi
  done
else
  mass=$el
  if [ `getOption -Tmax` != True ]; then error "if no elements are given -Tmax option needed"; fi
fi

echo " mass:          $mass" >> log_getDebye
if [ `getOption -Tmax` == True ]; then Tmax=`getValue -Tmax`; fi

# determine T start and step
opT=`getOption -Tstart`; if [ $opT == True ]; then Tstart=`getValue -Tstart`; else Tstart=$TstartDef; fi
if [ -z "$Tstart" ]; then error "no value to -Tstart option"; fi
opT=`getOption -Tstep`; if [ $opT == True ]; then Tstep=`getValue -Tstep`; else Tstep=$TstepDef; fi
if [ -z "$Tstep" ]; then error "no value to -Tstep option"; fi

# get T=0K energy file and check
file=`getOption -i`;  if [ $file  == True ]; then file=`getValue -i`; else file=$inputFileDef; fi
check $file

V0=`awk '{print $2}' $file`
B0=`awk '{print $3}' $file`
Bprime=`awk '{print $4}' $file`

# determine V range
opV=`getOption -Vstart`; if [ $opV == True ]; then Vstart=`getValue -Vstart`; else Vstart=`echo $Vs $V0 | awk '{printf "%.5f",$1*$2}'`; fi
if [ -z "$Vstart" ]; then error "no value to -Vstart option"; fi
opV=`getOption -Vstep`; if [ $opV == True ]; then Vstep=`getValue -Vstep`; else Vstep=$VstepDef; fi
if [ -z "$Vstep" ]; then error "no value to -Vstep option"; fi
opV=`getOption -nV`; if [ $opV == True ]; then nV=`getValue -nV`; else nV=$nVDef; fi
if [ -z "$nV" ]; then error "no value to -nV option"; fi
opV=`getOption -Vfix`; if [ $opV == True ]; then Vfix=`getValue -Vfix`; else Vfix=$V0; fi
if [ -z "$Vfix" ]; then error "no value to -Vfix option"; fi

echo "" >> log_getDebye
echo " Tstart:        $Tstart" >> log_getDebye
echo " Tstep:         $Tstep" >> log_getDebye
echo " Tmax:          $Tmax" >> log_getDebye
echo "" >> log_getDebye
echo " Vstart:        $Vstart" >> log_getDebye
echo " Vstep:         $Vstep" >> log_getDebye
echo " nV:            $nV" >> log_getDebye

echo; echo " running ..."
cp $path/utilities/Debye .
$path/fortran/debye.x $V0 $B0 $Bprime $mass $Vstart $Vstep $nV $Tstart $Tstep $Tmax $Vfix
rm Debye
echo " successful"


# generate parameters.math for getFqhFit.sh
$path/getFqhFit.sh -p FqhDebye-highT_ -f > /dev/null
sed -i -e 's|"bcc"|"aLatsAreVolumes"|' \
       -e 's|"FqhDebye-highT_"|"FqhDebye-highT_","FqhDebye-lowT_"|' parameters.math
echo " parameters.math for getFqhFit.sh (next step on cmpc) generated"

