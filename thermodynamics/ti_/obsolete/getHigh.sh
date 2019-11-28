#!/bin/sh

if [ ! -e avg_dUdL ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mERROR\033[0m: no avg_dUdL file" 1>&2
  exit 1
fi

eps=0.0000001

if [ $# -lt 1 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e " \033[1mgetHigh.sh avg_diff__*\033[0m                  list of avg_diff__* files for upsampling; one or sequence" 1>&2
  echo -e " \033[1mgetHigh.sh fit_{const,linear,quad}__*\033[0m   list of fit__* files for upsampling; one or sequence" 1>&2
  echo -e " \033[1mgetHigh.sh file_with_fileList\033[0m           file containing list of avg_diff__* or fit_* files in 1. column" 1>&2
  echo -e "                                         avg_diff__* and fit_* files can be mixed" 1>&2
  exit 1
fi

files=`echo $* | xargs -n1 | grep -e 'avg_diff__' -e 'fit_const__' -e 'fit_linear__' -e 'fit_quad__'`
if [ "$files" == "" ]; then files=`awk '{print $1}' $1`; fi
for f in $files; do
  if [ ! -e $f ]; then echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: no $f file" 1>&2; exit 1; fi
done
lambdas=`sed '/^[ ]*#/d' avg_dUdL | awk '{print $1}'`

echo "# upsampled from:" > high_avg_dUdL
echo $files | xargs -n1 | awk '{print "#",$1}' >> high_avg_dUdL
echo "# ----" >> high_avg_dUdL
cat avg_dUdL >> high_avg_dUdL

echo "# lambda  <dUdL>(meV) err/2(meV)"
for l in $lambdas; do
  dUdL=`awk '$1=='$l' {print $2}' avg_dUdL`
  err=` awk '$1=='$l' {print $3}' avg_dUdL`
  for f in $files; do
    if [ "`awk 'NF>1&&(($1-'$l')^2)^(1/2)<'$eps'{print $1}' $f`" == "" ]; then 
      echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: no $l in $f" 1>&2; exit 1;
    fi
    dUdL=`sed '/^[ ]*#/d' $f | awk 'NF>1&&(($1-'$l')^2)^(1/2)<'$eps'{printf("%6.2f",$2+'"$dUdL"')}'`
    err=`sed '/^[ ]*#/d' $f | awk '(($1-'$l')^2)^(1/2)<'$eps'{if(NF==5)printf("%9.2f",$3+'"$err"'); if(NF==2)printf("%9.2f",'"$err"')}'`
  done
  sed 's/\(^[ ]*'$l'\)[ ]*[^ ]*[ ]*[^ ]*\([ ]*.*\)/\1    '"$dUdL $err"'\2/' high_avg_dUdL > tmp_awk; mv tmp_awk high_avg_dUdL
  echo "    $l    $dUdL $err"
done

echo  $files | xargs -n1 > up_files
echo; echo up_files written
