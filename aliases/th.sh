#!/bin/sh

hier=`pwd`
# this is for cmpc/cmmd
#echo pwd:$hier
dort=`echo $hier | sed 's|^/data/|/home/|' | sed 's|^/home/|/home/|' | sed 's|^/nas/|/home/|' | sed 's|^/Users/|/home/|'`
#echo d1: $1
#echo dort1:$dort
[ "$1" != "" ] && dort=`echo $dort | sed 's|^/home/glensk/|'"$1"'/|'`
#echo dort2:$dort
echo $dort
