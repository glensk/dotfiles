#!/bin/bash
[ "$1" = "" ] && echo please specify the element first argument && exit
[ "$2" = "" ] && echo please specify the formation energy in eV as first argument && exit
tmelt=`getMeltingPoint.sh $1 -r`
kb=0.000086173423
for i in `seq 1 $tmelt`;do
echo $i | awk '{print '"$tmelt"'/$1,exp(-'"$2"'/($1*'"$kb"'))}'
done

