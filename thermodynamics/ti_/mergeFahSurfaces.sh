#!/bin/sh

if [ $# -lt 3 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e " \033[1mcombine_Fah_surfaces.sh strType combined_surface_name surface_1 surface_2 ...\033[0m     " 1>&2
  echo -e " \033[1mcombine_Fah_surfaces.sh strType combined_surface_name file_with_surfaces     \033[0m     one surface per line" 1>&2
  exit 1
fi

if [ $1 == bcc ]; then
  structureFactor=2;
else
  if [ $1 == fcc ]; then
    structureFactor=4;
  else
     if [ $1 == hcp ]; then
       if [ ! -e cBya ]; then echo -e "\033[31m\033[1mERROR\033[0m: no cBya file provided" 1>&2; exit 1; fi
       cBya=`cat cBya`
       structureFactor=`echo $cBya | awk '{printf "%.10f", 2/(sin(3.141592653589793/3)*$1)}'` # for hcp: 2/(Sin[Pi/3] cBya)
     else
       echo -e "\033[31m\033[1mERROR\033[0m: structure type $1 not known" 1>&2; exit 1;
     fi
  fi
fi

new=$2

if [ $# == 3 ]; then
  if [ `awk 'NR==1{print NF}' $3` == 1 ]; then
    surfaces=`cat $3`
  else
    surfaces=$3
  fi
else 
  surfaces=`echo $* | awk '{for (i=3;i<=NF;i++) print $i}'`
fi

c=`cat $surfaces | awk '{print $1}' | sed -n '/^0/p'`
if [ "$c" != "" ]; then echo; echo " ERROR: one of the surfaces contains leading zeros in the temperatures; remove and start again"; exit; fi

echo $surfaces | xargs -n1 > $new\_files

rm -fr tmp_Fah $new $new\_vol $new\_
for s in $surfaces; do
 if [ ! -e $s ]; then echo -e "\033[31m\033[1mERROR\033[0m: $s not existing" 1>&2; exit 1; fi
 vol=`echo $s | sed 's/.*\(_vol\)/\1/'`
 if [ "$vol" == "_vol" ]; then echo -e "\033[31m\033[1mERROR\033[0m: do not use *_vol surface" 1>&2; exit 1; fi
 cat $s >> tmp_Fah;
done
l=`awk '{print $1"_"$2}' tmp_Fah | sort -u`

mkdir $new\_
for i in $l; do
  t=`echo $i | sed 's/\(.*\)_\(.*\)/\1/'`
  a=`echo $i | sed 's/\(.*\)_\(.*\)/\2/'`
  v=`echo $a | awk '{printf("%.2f",$1^3/'$structureFactor')}'`
  f=`grep "^$t $a " tmp_Fah | awk 'BEGIN{F=0;inv=0}
                                  {Fi=$3; ei=2*$4; if (ei==0) ei=0.01; invi=1/ei^2; inv=inv+invi; F=F+invi*Fi}
                                  END{printf("%.2f %.2f",F/inv,(1/inv)^(1/2)/2)}'`
  echo "$t $a $f" >> $new
  echo "$t $v $f" >> $new\_vol
  echo "$a $f" >> $new\_/Fah_$t\K
  echo "$v $f" >> $new\_/Fah_$t\K_vol
  echo "$t $f" >> $new\_/Fah_$v\Ang
done

rm tmp_Fah
cat $new
