#!/bin/sh

[ ! -e "$1" ] && echo 'please specify the file by as first argument $1' && exit -1

l1=`head -1 $1 | wc -c`
l2=`head -2 $1 | tail -1 | wc -c`
l3=`head -3 $1 | tail -1 | wc -c`
la=`echo $l1 $l2 $l3 | awk '{print ($1+$2+$3)/3}'`
#echo l1:$l1 l2:$l2 l3:$l3 la:$la
fs=`stat --format=%s $1`
#echo fs:$fs
echo $fs $la | awk '{printf "%.0f\n", $1/$2}'

#file=$1; string=$(wc -c $file); bite=${string% *}; okay=$(echo "scale=2; $bite/1024" | bc);friend=$(echo -e #"$file $okay" "kb"); echo -e "$friend"
