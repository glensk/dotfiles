#!/bin/bash

#-----set parameters and paths------------------------
# parameters should be ok in most cases
pathToSpecies=/data/grabowski/elk-1.3.31/species/
lmaxapw=12        # angular momentum cut-off for the APW functions
lmaxvr=12         # angular momentum cut-off for the muffin-tin density and potential
lmaxmat=10        # angular momentum cut-off for the outer-most loop in the hamiltonian and overlap matrix setup
nempty=10         # number of empty states
epsengy=0.000005  # convergence criterion for the total energy
#-----------------------------------------------------


toHar=0.036749325405189714 # eV to Rydberg
toBohr=1.88972613289636593 # Angstrom to Bohr

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -p -f"

# small help for options
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-p     create parameters.dat and exit" \
               "-f     force folder creation even if existing"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   creates input files (elk.in) for a Murnaghan calculation with elk" 
  echo2 "   elk species file are assumed to be in:" \
        "   $pathToSpecies"
  echo2 "   parameters.dat must be present"
  echo2 "   check header of this script ($script) for some defaults"
  exit
fi

genPar=`getOption -p`
if [ $genPar == True ]; then
  echo "
aLats=a1 a2 a3 a4 ...   # Angstrom

elem=XX
xc=XX         # LDA or PBE
type=XX       # fcc or bcc

rmt=XX        # Bohr (muffin tin radius)
rmtkmax=XX    # 9..12
gmax=XX       # 12..18

kp=XX
sh=XX         # true or false; if true: kp shift=(0.5 0.5 0.5)

smear=XX      # Fermi, Gauss, or Methfessel
sigma=XX      # eV
" > parameters.dat
  echo; echo "parameters.dat written"; exit
fi

# get options
overwrite=`getOption -f`;

# read in all parameters from parameters.dat
check parameters.dat
elem=`get elem`; aLats=`get aLats`; rmtkmax=`get rmtkmax`; kp=`get kp`; sh=`get sh`; xc=`get xc`;
gmax=`get gmax`; smear=`get smear`; sigma=`get sigma`; rmt=`get rmt`; type=`get type`;
checkInput "$elem" "$aLats" "$rmtkmax" "$kp" "$sh" "$xc" "$gmax" "$smear" "$sigma" "$rmt" "$type"

# check if species file available
check $pathToSpecies/$elem.in

# check if given xc ok
if [ "$xc" == LDA ]; then xc=2; else if [ "$xc" == PBE ]; then xc=20; else error "xc not known"; fi; fi

# check if given strucuture type ok
case "$type" in
fcc)
  type="0.5  0.5  0.0
  0.5  0.0  0.5
  0.0  0.5  0.5";;
bcc)
  type="0.5  0.5  -0.5
  0.5  -0.5  0.5
  -0.5  0.5  0.5";;
*)
  error "structure type $type not known";;
esac

# k-points and shift
if [ $sh == true ]; then str="-shift"; sh="0.5 0.5 0.5"; else str=""; sh="0 0 0"; fi

# Fermi smearing
case "$smear" in
Gauss)
  smear_=0;;
Methfessel)
  smear_=1;;
Fermi)
  smear_=3;;
*)
  error "smearing method $smear not known";;
esac

# change sigma to sigmaHartree (elk unit)
sigmaHartree=`echo $sigma | awk '{printf("%.8f",$1*'$toHar')}'`

# run checks first
for a in $aLats; do
  # check if folder exists and should be overwritten
  f="$rmtkmax"\RmtKmax_Rmt=$rmt\_Gmax=$gmax\_$kp\x$kp\x$kp\kp$str\_$smear${sigma}eV/$a\Ang
  if [ -e $f -a "$overwrite" != True ]; then error "folder $f exists; run \"./createFolders.sh -f\" to force overwriting"; fi
done

dir=`pwd`; rm -f jobList
# now create the folders
for a in $aLats; do
  f="$rmtkmax"\RmtKmax_Rmt=$rmt\_Gmax=$gmax\_$kp\x$kp\x$kp\kp$str\_$smear${sigma}eV/$a\Ang
  echo; echo -e "\033[1m\033[31m$f\033[0m"; rm -fr $f; mkdir -p $f; cd $f

  abohr=`echo $a | awk '{printf("%9.6f",$1*'$toBohr')}'`

  # create elk.in file (main file)
  echo "
tasks
  0

xctype
  $xc

avec
  $type

scale
  $abohr

sppath
  './'

atoms
  1
  '$elem.in'
  1
  0.0  0.0  0.0    0.0  0.0  0.0

rgkmax
  $rmtkmax

ngridk
  $kp $kp $kp

vkloff
  $sh

gmaxvr
  $gmax

stype
  $smear_

swidth
  $sigmaHartree

nempty
  $nempty

lmaxapw
  $lmaxapw

lmaxvr
  $lmaxvr

lmaxmat
  $lmaxmat

epsengy
  $epsengy
" > elk.in

  # prepare species file
  rmtstr=`echo $rmt | awk '{printf("%6.4f",$1)}'`
  sed '/sprmin, rmt/s/\(^[ ]*[^ ]*[ ]*\)[^ ]*\(.*\)/\1'"$rmtstr"'\2/' $pathToSpecies/$elem.in > $elem.in

  cd $dir
  echo `pwd`/$f >> jobList
done

