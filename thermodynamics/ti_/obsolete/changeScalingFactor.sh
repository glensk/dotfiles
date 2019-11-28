#!/bin/sh

if [ $# != 1 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE:\033[0m AH_changeScalingFactor OUTCARfile    (zipped file possible; assumes dUdL in same folder)" 1>&2
  echo -e "       AH_changeScalingFactor -a          (runs AH_changeScalingFactor for OUTCAR* in all subfolders)" 1>&2
  exit 1
fi

if [ "$1" == "-a" ]; then
  l=`du -a | grep OUTCAR | awk '{print $2}'`
  for i in $l; do
      AH_changeScalingFactor.sh $i
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
             if( ((1000*($5-1*('$fref'))/('$nions'-1)-'$Udft')^2)^(1/2) < 0.01 ) print "NIONS - 1"; else print "ERROR"}'`

if [ "$which" = "ERROR" ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mERROR\033[0m: scaling factor undetermined" 1>&2
  exit 1
fi

if [ "$which" = "NIONS" ]; then
  scale=`echo $nions | awk '{printf("%.10f",$1/($1-1))}'`
  echo -e "$outcar   changing from    \033[1mNIONS\033[0m   to   \033[31m\033[1mNIONS - 1\033[0m"
else
  scale=`echo $nions | awk '{printf("%.10f",($1-1)/$1)}'`
  echo -e "$outcar   changing from    \033[1mNIONS - 1\033[0m   to   \033[31m\033[1mNIONS\033[0m"
fi
awk 'BEGIN{s='$scale'};/#/{print $0};
     !/#/{printf("%7d %9.1f %8.1f %8.1f %13.2f %9.2f %13.2f %9.2f %9.2f\n",$1,$2,$3,$4,s*$5,s*$6,s*$7,s*$8,s*$9)}' $dUdL > tmp
mv tmp $dUdL
