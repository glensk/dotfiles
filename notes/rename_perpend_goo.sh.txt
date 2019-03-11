#!/bin/sh

f=`ls -1d *`
for i in $f;do
    zw=`echo $i | sed 's|goo_||'`
    echo $i $zw
    mv $i $zw
done
