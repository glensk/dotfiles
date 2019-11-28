#!/bin/bash

if [ $# -gt 0 ]; then
   FileList=$1
else
   FileList=`find . -maxdepth 2 -type f -name "OUTCAR*" | sed 's|^./||g'`
   #FileList1=`ls -1 | grep OUTCAR`
   #FileList2=`ls -1 *Ang/OUTCAR*`
fi
#echo FileList:$FileList
#for i in $FileList;do
#		echo $i
#done
#echo ---

for i in $FileList; do
   #echo i:$i
   nions=`zgrep -a --text NIONS $i | awk '{print $12}'`
   disp=`zgrep -a --text -A1 "position of ions in cartesian" $i | tail -n1 | \
         awk '{print $1}'`
   s=`echo $i | sed 's/OUTCAR\.//' | sed 's/\.gz//' | sed 's|Ang.*||' | sed 's|/gz||' | sed 's|/OUTCAR||'`
   #echo i:$i s:$s
   drift=`zgrep -a --text drift $i | awk '{print $1,$2,$3,$4,$5}'`
   echo -e "$i   || $s \t$nions      $disp         $drift"
   copyto=`echo $i | sed 's|/OUTCAR.*|/|'`
   zgrep -a --text -A`expr $nions + 1` "TOTAL-FORCE" $i | tail -n+3 | \
      awk '{print $4,$5,$6}' > forces.$s
   echo $disp > disp.$s
done;
