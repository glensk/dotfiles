#!/bin/bash

#-----set parameters and paths------------------------------------
f="OUTCAR OUTCAR.gz workDir/OUTCAR"; atomFile=atoms_volume_steps
cellOut=cell; output=POSITIONs; offDef=0; corrDef=0
#-----------------------------------------------------------------

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -i -o -p -f -a -O -c -eq -r"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-i inpFile  input file to use (default: check for $f)" \
               "-o outFile  output file for POSITIONs and/or TOTAL-FORCEs (default: $output)" \
               "-p          extract only positions" \
               "-f          extract only forces" \
               "-a atom     extract only values for atom" \
               "-O off      output after offset off (default: $offDef)" \
               "-c corr     correlation length; output only every corr'th step (default: $corrDef)" \
               "-eq         write equilibrium (first) structure to equilibrium file" \
               "-r          remove mapping into original cell"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   extracts POSITIONs and TOTAL-FORCEs from an OUTCAR file for each MD step"
  echo2 "   the input OUTCAR file can be in zipped format"
  echo2 "   both POSITIONs and TOTAL-FORCEs are written by default into the same file" \
        "   (default: $output; change with -o) with the format (x,y,z: coordinates;" \
        "   Fx,Fy,Fz: forces; first index over atoms; second over MDsteps; N=nrOfAtoms;" \
        "   M=nrOfMDsteps):" \
        "      x11  y11  z11  Fx11  Fy11  Fz11" \
        "      x21  y21  z21  Fx21  Fy21  Fz21" \
        "       .    .    .     .    .     .  " \
        "       .    .    .     .    .     .  " \
        "      xN1  yN1  zN1  FxN1  FyN1  FzN1" \
        "      x12  y12  z12  Fx12  Fy12  Fz12" \
        "      x22  y22  z22  Fx22  Fy22  Fz22" \
        "       .    .    .     .    .     .  " \
        "       .    .    .     .    .     .  " \
        "      xN2  yN2  zN2  FxN2  FyN2  FzN2" \
        "       .    .    .     .    .     .  " \
        "       .    .    .     .    .     .  " \
        "       .    .    .     .    .     .  " \
        "      x1M  y1M  z1M  Fx1M  Fy1M  Fz1M" \
        "      x2M  y2M  z2M  Fx2M  Fy2M  Fz2M" \
        "       .    .    .     .    .     .  " \
        "       .    .    .     .    .     .  " \
        "      xNM  yNM  zNM  FxNM  FyNM  FzNM" 
  echo2 "   the units of the POSITIONs and TOTAL-FORCEs are as in OUTCAR (Ang and eV/Ang)"
  echo2 "   using -p or -f ony the POSITIONs or TOTAL-FORCEs can be extracted"
  echo2 "   an offset can be given ('-O off' option) to skip the first off steps" \
        "   and/or a correlation length ('-c corr' option) to skip steps in between"
  echo2 "   POSITIONs/TOTAL-FORCEs can be extracted only for a specfic atom but all" \
        "   MD steps using the '-a atom' option"
  echo2 "   in addition to the POSITIONs file a $cellOut file is created containing" \
        "   the real space cell (Ang) and a $atomFile file containing the number of" \
        "   atoms, the cell volume (Ang^3), the number of all MD steps, and the number" \
        "   of extracted MD steps (the latter differs from all MD steps if offset or" \
        "   correlation length are used)"
  echo2 "   the POSITIONs given in the OUTCAR are always mapped into the original unit" \
        "   cell; this means that an atom which has vibrated out of the cell is mapped" \
        "   across the cell to the other side which results in a discontinuous jump" \
        "   this mapping can be removed with the -r option (arbitrary cells possible)"
  exit
fi

# get and check input file
input=`getOption -i`
if [ $input == True ]; then
  file=`getValue -i`; check $file;
else
  file=0; for i in $f; do if [ -e $i ]; then file=$i; break; fi; done
  if [ $file == 0 ]; then error "no input file existing; checked for $f"; fi
fi
echo; echo " extracting POSITIONs from $file"

# check if input file is compressed
zip=`file $file | sed 's/.*gzip compressed data.*/True/'`
if [ "$zip" == True ]; then zgrep -a --text ".*" $file > _tmp_file; file=_tmp_file; fi

# get cell, nr of atoms, and volume
awk 'BEGIN{l=-10};/NIONS =/{n=$NF};/volume of cell/{v=$NF;l=NR};
     NR>l+1&&NR<l+5{print $1,$2,$3};NR==l+5{print n,v;exit}' $file > $cellOut

# separate cell and atoms and volume
l=`sed -n '4p' $cellOut`;
#echo co:$cellOut
nAt=`echo $l | awk '{print $1}'`
n=`echo $nAt | awk '{print $1+2}'`
sed -i '4d' $cellOut


# NOTE !!!!!!!
# it is not faster to implement the -p or -f option at this point by putting the procedure
# into the sed command; the reason is that the corresponding matching procedure in sed
# takes long times because it is much more involved when trying to seperate positions and forces
# NOTE !!!!!!!


# get all POSITIONs from OUTCAR; we also write just the flag POSITION into _tmp_
# in order to count the nr of steps from this
sed -n -e '/^ POSITION/,+'$n'{s/[0-9]\+/\0/p}' -e '/^ POSITION/w _tmp_' $file  | awk '{print $1,$2,$3,$4,$5,$6}' > _tmp_POSITIONs
steps=`wc -l _tmp_ | awk '{print $1}'`; rm _tmp_

# create file with nr of atoms, volume, and nr of MD steps
echo -n $l $steps > $atomFile

# check if we are using different offset and correlation length than default
off=`getOption -O`;
if [ $off == True ]; then
  off=`getValue -O`; c=`checkInteger $off`
  if [ $c != ok ]; then error "wrong or empty value for -O option"; fi
else
  off=$offDef
fi
corr=`getOption -c`;
if [ $corr == True ]; then
  corr=`getValue -c`; c=`checkInteger $corr`
  if [ $c != ok ]; then error "wrong or empty value for -c option"; fi
else
  corr=$corrDef
fi

# check if we extract only for a specific atom
atom=`getOption -a`
if [ $atom == True ]; then
  atom=`getValue -a`; c=`echo $atom | awk '$1>=1&&$1<='$nAt'{print "ok"}'`
  if [ "$c" != ok ]; then error "wrong or empty value for -a option"; fi
else
  atom=0 
fi

#check if we extract only positions or forces
pos=`getOption -p`; forc=`getOption -f`
n1=1; n2=6                                  # get them all from column 1 to column 6
if [ $pos  == True ]; then n1=1; n2=3; fi   # get only the columns 1 to 3 (==positions)
if [ $forc == True ]; then n1=4; n2=6; fi   # get only the columns 4 to 6 (==forces)

# check if we put equilibrium structure into separate file
eq=`getOption -eq`
#echo "eq:",$eq,$nAt,$n1,$n2
if [ $eq == True ]; then
  awk 'NR<='$nAt'{for (i='$n1';i<='$n2';i++) printf("%14s ",$i); printf("\n")}' _tmp_POSITIONs > equilibrium
fi

# check if we have different output file name
out=`getOption -o`;
if [ $out == True ]; then output=`getValue -o`; fi

# if offset or correlation length not default then extract corresponding coordinates
if [ $off != 0 -o $corr != 0 -o $pos == True -o $forc == True -o $atom != 0 ]; then
  awk 'BEGIN{ nAt='$nAt'; corr='$corr'; off='$off'; atom='$atom';
              n1=off*nAt; n2=n1+nAt }
       NR>n1 && NR<=n2 {if (atom==0||atom==NR-n1) {for (i='$n1';i<='$n2';i++) printf("%14s ",$i); printf("\n")}}
       NR==n2{n1=n2+corr*nAt; n2=n1+nAt}' _tmp_POSITIONs > $output
  rm _tmp_POSITIONs
  # get number of extracted steps
  extracted=`wc -l $output | awk '{print $1/'$nAt'}'`
else
  mv _tmp_POSITIONs $output
  # nr of extracted steps will be the same as all steps if no offset or correlation length used
  extracted=$steps
fi

# add nr of extracted steps
echo " $extracted" >> $atomFile

# rm temporary files
rm -f _tmp_file

# check if POSITIONs is empty
n=`awk 'END{print NR}' $output`
if [ $n == 0 ]; then rm $output; error "no POSITIONs found in $file"; fi

# remove jumps if -r option
rem=`getOption -r`
if [ $rem == True ]; then
  if [ $forc == True ]; then error "-r option not compatible with -f option"; fi
  if [ $atom != 0 ]; then nAt=1; fi
  if [ $pos == True ]; then force=0; else force=1; fi

  $path/fortran/removeMappingInMD.x $output $nAt `cat cell | xargs` $force
  mv ${output}_noJumps $output

fi

