#!/bin/sh

if [ $# -lt 1 -o $# -gt 3 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e " get_Fah.sh fitType structure       (fitType: linear, cubic, tangens; structure: bcc or fcc)" 1>&2
  echo -e " get_Fah.sh fitType structure -g    (to add blank lines in Fah_surface for gnuplot)" 1>&2
  exit 1
fi

# check if "-g" option given
gnuplot=`echo $* | xargs -n1 | grep -e "-g" | grep "\<g\>"`

# get fit type
fit=`echo $* | xargs -n1 | grep -e lin -e cub -e tan`
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

fah=Fah
l=`ls -1d *Ang_*K | wc -l`
if [ "$l" == 0 ]; then
  echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: no *Ang_*K folders available" 1>&2; exit 1
fi

l=`ls -1d *Ang_*K`
rm -f Fah_*
echo $structureFactor > Fah_from_fit_$fit
echo; echo -e "using fit \033[31m\033[1m$fit\033[0m"; echo
echo "aLat(Ang)   T(K)  Fah(meV/at) err(meV/at)"
for i in $l; do
  if [ -e $i/$fah ]; then
    is=`grep $fit $i/$fah | wc -l`
    if [ "$is" == 1 ]; then
      a=`echo $i | sed 's/\(.*\)Ang_\(.*\)K[\/]*/\1/'`
      v=`echo $a | awk '{printf("%.2f",$1^3/'$structureFactor')}'`
      t=`echo $i | sed 's/\(.*\)Ang_\(.*\)K[\/]*/\2/'`
      f=`grep $fit $i/$fah | awk '{print $2,"   ",$3}'`
      echo $t $f >> Fah_$a\Ang
      sort -g Fah_$a\Ang > tmp; mv tmp Fah_$a\Ang
      echo $a $f >> Fah_$t\K
      echo $v $f >> Fah_$t\K_vol
      sort -g Fah_$t\K > tmp; mv tmp Fah_$t\K
      sort -g Fah_$t\K_vol > tmp; mv tmp Fah_$t\K_vol
      echo $t $a $f >> Fah_surface
      echo $t $v $f >> Fah_surface_vol
      echo "   $a     $t     $f"
    else
      echo no fit $fit in $i/$fah
    fi
  else
    echo no $fah  in $i
  fi
done

if [ "$gnuplot" == "-g" ]; then
# add blank lines for pm3d in gnuplot
  sort -g Fah_surface | awk ' $1 != prev {printf "\n"; prev=$1} 
                   {print}' > tmp
  mv tmp Fah_surface
  sort -g Fah_surface_vol | awk ' $1 != prev {printf "\n"; prev=$1} 
                   {print}' > tmp
  mv tmp Fah_surface_vol
fi

echo; echo Files written:; echo " Fah_from_fit_$fit"; echo " Fah_surface Fah_surface_vol";
echo -n " "; echo Fah_*K; echo -n " "; echo Fah_*K_vol; echo -n " "; echo Fah_*Ang
rm -f HesseRef_changed
