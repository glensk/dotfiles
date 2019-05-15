#!/bin/sh
f=`ls */input.nn`
hier=`pwd`
for i in $f;do
    cd $hier
    fo=`echo $i | sed 's|input.nn||'`
    echo $fo
    cd $fo
    getEnergies_byLammps.py -p . -e
    cd $hier
done

