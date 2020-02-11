#!/bin/sh
ff=`ls -1d n2p2_v6ag_*`
hier=`pwd`
for f in $ff;do
    cd $hier
    echo $f
    cd $f
        #getEnergies_byLammps.py -p . -sys fcc -sys_ele Al -evinet
        #getEnergies_byLammps.py -p . -sys fcc -sys_ele Al -sys_ncell 1 -ea
    cd $hier
done


