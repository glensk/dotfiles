#!/bin/sh
f=`ls */input.nn`
hier=`pwd`
for i in $f;do
    cd $hier
    fo=`echo $i | sed 's|input.nn||'`
    echo $fo
    cd $fo

    mkdir kmc57
    cd kmc57
    getEnergies_byLammps.py -p ../ -kmc

    #getEnergies_byLammps.py -p . -e
    cd $hier
done

