#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -p -f -i -r"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-p         create parameters.dat and exit" \
               "-f         force template or folder creation even if existing" \
               "-i         generate only spcm input file, do not run spcm" \
               "-r         run spcm using the available input file (do not modify it)"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   generates an SQS supercell in POSCAR format"
  echo2 "   parameters.dat file must be present (generate with -p)"
  echo2 "   $script is only a wrapper for the spcm code developed" \
        "   by Andrei Ruban from KTH who should be always acknowledged"
  echo2 "   $script presently supports only the ollowing structures:" \
        "   bcc, bccPrim, fcc, fccPrim, hcp" \
        "   and for each of these only a single sublattice, only pair" \
        "   correlation functions, and fixed weights"
  echo2 "   for greater flexibility the original spcm code needs to" \
        "   be used, e.g., by using first the -i option, then modifying" \
        "   the input file in.dat and then using the -r option"
  echo2 "   $script outputs the deviations of the achieved pair" \
        "   correlation functions from the perfectly disordered ones" \
        "   the best possible value for the deviation is 0 which means" \
        "   that the pair correlation is ideally matching its target" \
        "   value; larger or smaller values show actual deviations"
  echo2 "   the larger the supercell the easier it is to achieve small" \
        "   deviations"
  echo2 "   it is difficult to give a generic critical magnitude above" \
        "   above which the corresponding correlation functions are" \
        "   critical"
  exit
fi

# check if force overwriting templates or folders
overwrite=`getOption -f`

# create parameters.dat template and exit
genPar=`getOption -p`
if [ $genPar == True ]; then
  if [ -e parameters.dat -a $overwrite != True ]; then error "parameters.dat existing; use -f to overwrite"; fi

  # write parameters.dat
  echo "
str  = bcc                # crystal structure (bcc,bccPrim,fcc,fccPrim,hcp)
sc   = 3 3 3              # supercell size
conc = 0.5 0.5            # concentrations for each species, must add up to 1
alat = xxxALATxxx         # scaling for aLat in POSCAR, can be xxxALATxxx tag
" > parameters.dat
  echo; echo "parameters.dat written"; exit
fi


if [ `getOption -r` == True ]; then
  # if -r we only run with an availble in.dat
  check in.dat
  # delete CR/LF line terminators (just in case because spcm will not work)
  sed 's|\r||g' in.dat > _tmp_in.dat
  mv _tmp_in.dat in.dat

  # get cell, supercell size, concentrations, and labels from in.dat
  cell=`grep BSX in.dat | sed 's|=| |g' | awk '{print $2,$4,$6}' | xargs`
  na=`grep "NA.*NB.*NC" in.dat | awk '{printf "%s",$2}'`
  nb=`grep "NA.*NB.*NC" in.dat | awk '{printf "%s",$4}'`
  nc=`grep "NA.*NB.*NC" in.dat | awk '{printf "%s",$6}'`
  conc=`grep -A1 "Concentrations on sublattices" in.dat | awk 'END{print}'`
  c=`echo $conc | awk '{s=0; for (i=1;i<=NF;i++) {s=s+$i}; if (s!=1) print "error"; else print "ok"}'`
  if [ "$c" == error ]; then error "concentrations do not add up to 1 (must be exactly 1)"; fi
  labels=`grep "SMB(IAT)" in.dat | awk '{for (i=2;i<=NF;i++) printf "%s ",$i}'`

  echo; echored "  WARNING:"; echo "  using cell, supercell size, concentrations from in.dat (not from parameters.dat)"; echo

else
  # -------------------- start here construction of new in.dat if not -r option ----------------
  if [ -e in.dat -a $overwrite != True ]; then
    error "previous in.dat available, use -f to overwrite"
  fi

  # read in all parameters from parameters.dat
  check parameters.dat
  str=`get str`; sc=`get sc`; conc=`get conc`; alat=`get alat`
  checkInput "$str" "$sc" "$conc" "$alat"

  # generate spcm input file
  # general header
  echo "SPC       HP......=N                                         00 May 00
JOBNAM...=out        MSGL.=  0
FOR001=spc/
FOR002=spc/
FOR004=spc/
FOR006=
FOR008=spc/
FOR009=spc/
Comment
NPRN..=  0 TEST.=  0 NCOL.=  0 STAT.=  1 TCLIM=  0 nsho.= 10" > in.dat

  # crystal structure specific information; set also number of sublattices needed below
  case "$str" in
    bccPrim) nSites=1; cell=".5 .5 -.5 -.5 .5 .5 .5 -.5 .5"
    echo "NQ3...=  1 LAT..=  3 IPRIM=  0 HIGH.=  0 NSHC.= 10 NL...=  0 NLH..=  0
A........=       1.0 B.......=       1.0 C.......=1.0
BSX......=       0.5 BSY.....=       0.5 BSZ.....=      -0.5
BSX......=      -0.5 BSY.....=       0.5 BSZ.....=       0.5
BSX......=       0.5 BSY.....=      -0.5 BSZ.....=       0.5
QX.......=       0.0 QY......=       0.0 QZ......=0.0        it..=   1" >> in.dat;;

    fccPrim) nSites=1; cell="0 .5 .5 .5 0 .5 .5 .5 0"
    echo "NQ3...=  1 LAT..=  2 IPRIM=  0 HIGH.=  0 NSHC.= 10 NL...=  0 NLH..=  0
A........=       1.0 B.......=       1.0 C.......=1.0
BSX......=       0.0 BSY.....=       0.5 BSZ.....=       0.5
BSX......=       0.5 BSY.....=       0.0 BSZ.....=       0.5
BSX......=       0.5 BSY.....=       0.5 BSZ.....=       0.0
QX.......=       0.0 QY......=       0.0 QZ......=0.0        it..=   1" >> in.dat;;

    bcc) nSites=2; cell="1 0 0 0 1 0 0 0 1"
    echo "NQ3...=  2 LAT..=  1 IPRIM=  0 HIGH.=  0 NSHC.= 10 NL...=  0 NLH..=  0
A........=       1.0 B.......=       1.0 C.......=1.0
BSX......=       1.0 BSY.....=       0.0 BSZ.....=       0.0
BSX......=       0.0 BSY.....=       1.0 BSZ.....=       0.0
BSX......=       0.0 BSY.....=       0.0 BSZ.....=       1.0
QX.......=       0.0 QY......=       0.0 QZ......=0.0        it..=   1
QX.......=       0.5 QY......=       0.5 QZ......=0.5        it..=   1" >> in.dat;;

    fcc) nSites=4; cell="1 0 0 0 1 0 0 0 1"
    echo "NQ3...=  4 LAT..=  1 IPRIM=  0 HIGH.=  0 NSHC.= 10 NL...=  0 NLH..=  0
A........=       1.0 B.......=       1.0 C.......=1.0
BSX......=       1.0 BSY.....=       0.0 BSZ.....=       0.0
BSX......=       0.0 BSY.....=       1.0 BSZ.....=       0.0
BSX......=       0.0 BSY.....=       0.0 BSZ.....=       1.0
QX.......=       0.0 QY......=       0.0 QZ......=0.0        it..=   1
QX.......=       0.0 QY......=       0.5 QZ......=0.5        it..=   1
QX.......=       0.5 QY......=       0.0 QZ......=0.5        it..=   1
QX.......=       0.5 QY......=       0.5 QZ......=0.0        it..=   1" >> in.dat;;

    hcp) nSites=2; cell="1 0 0 -0.5 0.8660254 0 0 0 1.63299316";
    echo "NQ3...=  2 LAT..=  4 IPRIM=  0 HIGH.=  0 NSHC.= 10 NL...=  0 NLH..=  0
A........=       1.0 B.......=       1.0 C.......=1.63299316
BSX......= 1.0000000 BSY.....= 0.0000000 BSZ.....= 0.0000000
BSX......=-0.5000000 BSY.....= 0.8660254 BSZ.....= 0.0000000
BSX......= 0.0000000 BSY.....= 0.0000000 BSZ.....=1.63299316
QX.......=       0.0 QY......=       0.0 QZ......=0.0        it..=   1
QX.......=       0.0 QY......=0.57735027 QZ......=0.81649658 it..=   1" >> in.dat;;

    *) error "strType not known";;
  esac


  # only Madelung matrix, not relevant for SQS?
  echo "LAMDA....=    2.5000 AMAX....=    4.5000 BMAX....=    4.5000" >> in.dat


  #supercell size related
  na=`echo $sc | awk '{printf "%3d",$1}'`
  nb=`echo $sc | awk '{printf "%3d",$2}'`
  nc=`echo $sc | awk '{printf "%3d",$3}'`
  echo "Size of the super cell
NA.......= $na NB.......= $nb NC.......= $nc  Dmax     4.5
NSDC.....=  -6 NSDS.....=  -1 NSDM.....=   1
NMAXMX...=   3 TMLIM....= 1.0" >> in.dat


  # sublattice related; nSites set above already
  sublat=`echo $nSites | awk '{for (i=1;i<=$1;i++) printf "   1"}'`
  echo "NT.......=   1
NTA(IQ)..=$sublat
NTO......=   1
NTAO(IQ).=$sublat
NQ3O.....=   1
IQO(IQ)..=$sublat" >> in.dat


  # number of atoms and concentrations
  c=`echo $conc | awk '{s=0; for (i=1;i<=NF;i++) {s=s+$i}; if (s!=1) print "error"; else print "ok"}'`
  if [ "$c" == error ]; then error "concentrations do not add up to 1 (must be exactly 1)"; fi
  natoms=`echo $conc | awk '{print NF}'`
  labels=`echo $conc | awk '{for (i=65;i<65+NF;i++) printf "   %c",i}'`
  echo "NATOM....=   $natoms
SMB(IAT).=$labels
Concentrations on sublattices:
$conc" >> in.dat


  # correlation functions and weigths
  echo "Correlation functions and weights for each pairs of elem. (A-B, A-C, ... )
Sublattice 
1
nc2  r_max
4    3.0
i    alpha           weight" >> in.dat

  # pair correlations
  for (( i=1; i<=$natoms; i++ )); do
    for (( j=$i+1; j<=$natoms; j++ )); do 

      a=`echo $i | awk '{printf "%c",64+$1}'`
      b=`echo $j | awk '{printf "%c",64+$1}'`

      echo "1     0.0          1000.0  ${a}-$b
2     0.0           800.0
3     0.0           400.0
4     0.0           200.0" >> in.dat

    done
  done

  # triplets and quadruplets all set to zero
  echo "nc3 
0
i   i1  i2  i3           <sss>          weight
nc4
0
i   i1 i2 i3 i4 i5 i6    <ssss>         weight" >> in.dat


  # temperature and steps for MC
  echo "T_i,    T_f,    delt_T
 50.0    00.0    10.0
 20" >> in.dat

fi  # -------------------- end here construction of new in.dat if not -r option ----------------


# run spcm if not -i option
if [ `getOption -i` != True ]; then
  mkdir spc
  $path/sqs/spcm.exe < in.dat > log_spcm
  rm -fr spc # out.log

  # check for Error
  c=`grep Error log_spcm`
  if [ "$c" != "" ]; then error "error in spcm, check log_spcm (supercell too small?)"; fi

  # printout pair correlations
  echo "                                       deviation from"
  echo "  pair shell      coordinate           perfect disorder (=0)"
  awk 'BEGIN{f=0;ff=0} /Optimzied correlation functions/{f=1} /Alloy components/{if (f==1) {ff=1; a=$3; b=$4}}
       f==1&&$1=="<ss>_"{if (ff==1) {ff=0; printf "  %s-%s",a,b} else printf "     ";
                         printf "   %2s %7s %7s %8s     %9s\n", $2,$3,$4,$5,$8}' out.prn
  echo
  echo "  check out.prn and -help for more information"

  echo Comment > POSCAR
  echo $alat >> POSCAR
  grep "   a" out.xyz | awk '{print $2,$3,$4}' >> POSCAR
  grep "   b" out.xyz | awk '{print $2,$3,$4}' >> POSCAR
  grep "   c" out.xyz | awk '{print $2,$3,$4}' >> POSCAR
  nall=""
  for i in $labels; do
    n=`grep ' '$i'$' out.xyz | awk 'END{print NR}'`
    nall="$nall $n"
    echo -n " $n " >> POSCAR
    grep ' '$i'$' out.xyz | awk '{print $1,$2,$3}' >> _tmp_POSCAR
  done
  echo "" >> POSCAR
  echo "Cartesian" >> POSCAR
  cat _tmp_POSCAR >> POSCAR
  rm _tmp_POSCAR

  echo; echo "  POSCAR with SQS supercell created"

  # calculate new concentrations and compare to initial ones
  concNew=`echo $nall | awk '{s=0; for (i=1;i<=NF;i++) {s=s+$i}; for (i=1;i<=NF;i++) printf "%.4f ",$i/s}'`
  concOrig=`echo $conc | awk '{for (i=1;i<=NF;i++) printf "%.4f ",$i}'`
  if [ "$concNew" != "$concOrig" ]; then
    echo; echored "  WARNING: concentrations could not be fully satisfied"
             echo "           original concentrations: $concOrig"
             echo "                new concentrations: $concNew"
  fi
fi

# acknowledgement Andrei
echo
echo "  ACKNOWLEDGMENT: if you use results for publications, acknowledge:"
echo "  Andrei Ruban, Royal Institute of Technology, Stockholm, Sweden"



