#!/bin/sh

if [ $# != 1 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE:\033[0m AH_checkScalingFactor OUTCARfile    (zipped file possible; assumes dUdL in same folder)" 1>&2
  echo -e "       AH_checkScalingFactor -a          (runs AH_checkScalingFactor for OUTCAR* in all subfolders)" 1>&2
  exit 1
fi

if [ "$1" == "-a" ]; then
  l=`du -a | grep OUTCAR | awk '{print $2}'`
  for i in $l; do
      AH_checkScalingFactor.sh $i
  done
  exit
fi

outcar=$1
dUdL=`echo $outcar | sed 's|\(.*\)[^/]*OUTCAR[^/]*|\1|'`dUdL
if [ ! -e $outcar -o ! -e $dUdL ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mERROR\033[0m: $outcar or $dUdL not available" 1>&2
  exit 1
fi
Udft=`tail -n2 $dUdL | head -n1 | awk '{print $5}'`
n=`tail -n2 $dUdL | head -n1 | awk '{print $1+1}'`
nions=`zgrep "NIONS =" $outcar | awk '{print $12}'`
fref=`zgrep "free  energy" $outcar | head -n1 | awk '{print $5}'`
which=`zgrep "free  energy" $outcar | sed -n ''$n'p' | \
       awk '{if( ((1000*($5-1*('$fref'))/'$nions'-'$Udft')^2)^(1/2) < 0.01 ) print "NIONS"; else
             if( ((1000*($5-1*('$fref'))/('$nions'-1)-'$Udft')^2)^(1/2) < 0.01 ) print "NIONS - 1"; else print "!!ERROR!! no factor matches"}'`

echo -e "$outcar        NIONS = $nions        scaling factor = \033[31m\033[1m$which\033[0m"
