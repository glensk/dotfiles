#!/bin/sh

if [ $# -lt 1 -o $# -gt 3 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e " high_AvgEn.sh fitType structure       (fitType: linear or quad; structure: bcc or fcc)" 1>&2
  echo -e " high_AvgEn.sh fitType structure -g    (to add blank lines in Fah_surface for gnuplot)" 1>&2
  exit 1
fi

# check if "-g" option given
gnuplot=`echo $* | xargs -n1 | grep -e "-g" | grep "\<g\>"`

# get fit type
fit=`echo $* | xargs -n1 | grep -e lin -e quad -e const | head -n1`
if [ "$fit" == "" ]; then
  echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: fitType not known" 1>&2; exit 1
fi

# get structure
str=`echo $* | xargs -n1 | grep -e bcc -e fcc`
if [ "$str" == "" ]; then
  echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: structure not known" 1>&2; exit 1
fi
if [ $str == bcc ]; then structureFactor=2; fi
if [ $str == fcc ]; then structureFactor=4; fi

l=`ls -1d *Ang_*K | wc -l`
if [ "$l" == 0 ]; then
  echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: no *Ang_*K folders available" 1>&2; exit 1
fi

l=`ls -1d *Ang_*K`
rm -f AvgEn__*; touch AvgEn__from_fit_$fit
echo; echo -e "using fit \033[31m\033[1m$fit\033[0m"; echo
echo "aLat(Ang)   T(K)  AvgEn(meV/at) err(meV/at)"
for i in $l; do
 exts=`ls -1 $i/AvgEn__* 2> /dev/null | sed 's|.*AvgEn__\(.*\)|\1|'`
 for ext in $exts; do
   a=`echo $i | sed 's/\(.*\)Ang_\(.*\)K[\/]*/\1/'`
  v=`echo $a | awk '{printf("%.2f",$1^3/'$structureFactor')}'`
   t=`echo $i | sed 's/\(.*\)Ang_\(.*\)K[\/]*/\2/'`
   f=`grep $fit $i/AvgEn__$ext | awk '{print $2,"   ",$3}'`
   echo $t $f >> AvgEn__$ext\__$a\Ang;   sort -g AvgEn__$ext\__$a\Ang   > tmp; mv tmp AvgEn__$ext\__$a\Ang
   echo $a $f >> AvgEn__$ext\__$t\K;     sort -g AvgEn__$ext\__$t\K     > tmp; mv tmp AvgEn__$ext\__$t\K
   echo $v $f >> AvgEn__$ext\__$t\K_vol; sort -g AvgEn__$ext\__$t\K_vol > tmp; mv tmp AvgEn__$ext\__$t\K_vol
   echo $t $a $f >> AvgEn__$ext\__surface
   echo $t $v $f >> AvgEn__$ext\__surface_vol
   echo "   $a     $t     $f"
 done
done

l=`ls AvgEn__*_surface*`
for i in $l; do
  if [ "$gnuplot" == "-g" ]; then
  # add blank lines for pm3d in gnuplot
    sort -g AvgEn__$i\_surface | awk ' $1 != prev {printf "\n"; prev=$1} 
                     {print}' > tmp
    mv tmp AvgEn__$i\_surface
  fi
done

echo; echo Files written:; ls AvgEn__*;

