#!/bin/sh
[ "$1" = "" ] && echo 'provide the cursel.def file to romve 0 for 3body interactions' && exit
[ ! -e "$1" ] && echo $1 does not exist && exit
lines=`grep -n ' 3 ' $1 | sed 's|:.*||' | xargs`
for line in $lines;do
    #echo $line
    #sed -i ''"$line"'' | 0.000E+00 ||' $1
    sed -i "$line  s| 0.000E+00 |  |" $1
done
