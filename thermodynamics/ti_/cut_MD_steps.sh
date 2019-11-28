#!/bin/bash

if [ $# != 2 ]; then
  echo; echo USAGE: cut.sh lambdaFOLDER MDstepToCutFrom
  exit
fi

lam=$1
n=$2

#transtr=ALL-STEPS-TRANSITION-TO-DHCP
transtr=ALL-STEPS-TRANSITION

if [ ! -e $lam ]; then
  echo ERROR: no $lam folder
  exit
fi

if [ ! -e $lam/dUdL -o ! -e $lam/OUTCAR.gz -o ! -e $lam/vasprun.xml.gz ]; then
  echo; echo ERROR: no $lam/dUdL or $lam/OUTCAR.gz or $lam/vasprun.xml.gz
  exit
fi

# backup the folder
s=`echo $lam | sed -e 's|lambda||' -e 's|/||'`
if [ -e LAMBDA${s}__$transtr ]; then
  echo; echo ERROR: backup folder exists: LAMBDA${s}__$transtr
  exit
fi
mv $lam LAMBDA${s}__$transtr
mkdir $lam
cp LAMBDA${s}__$transtr/{dUdL,OUTCAR.gz,vasprun.xml.gz,POSCAR} $lam/

# cut dUdL
awk 'NR<='$n'+1 {print}' $lam/dUdL > _tmp
mv _tmp $lam/dUdL

# cut OUTCAR
nn=`echo $n | awk '{print $1 + 2}'`
zcat $lam/OUTCAR.gz | sed '/Iteration *'$nn'/,$d' > _tmp
gzip _tmp
mv _tmp.gz $lam/OUTCAR.gz

# cut vasprun.xml
zcat $lam/vasprun.xml.gz | awk 'BEGIN{s=0} {print}; /nosekinetic/{s=s+1}; /calculation/{if (s=='$n'+1) exit}' > _tmp
gzip _tmp
mv _tmp.gz $lam/vasprun.xml.gz

echo !!!!!!!!!!!!
echo " WARNING: Remember to delete corresponding structures_uncor file to trigger an update!"
echo !!!!!!!!!!!!


