#!/bin/sh
cd $dotfiles/scripts/bin
pwd
subfolder="n2p2 qe-aiida python_thermodynamics i-pi-mc_scripts runner_scripts"
for i in $subfolder;do
    echo $i
    execs=`ls ../$i/*.sh ../$i/*.py`
    for j in $execs;do
        z=`echo $j | sed 's|../'"$i"'/||g'`
        echo "->" $j
        echo "-->" $z 
        ln -s $j $z
        echo
    done
done

