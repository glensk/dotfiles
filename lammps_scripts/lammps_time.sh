#!/bin/sh

if [ "$1" = "" ];then
    if [ -e "log.lammps" ];then 
        pfad=log.lammps
    else
        echo 'need $1 to be pfad to file'
        exit
    fi
elif [ "$1" = "all" ];then
    pfad=`find . -name "log.lammps" | sort -n`
else
pfad=$1
fi
echo $pfad

for pfadd in $pfad;do
    o=`tail -200 $pfadd | grep "^Loop time of"`
    #echo o:$o:
    if [ "$o" = "" ];then
        echo "$pfadd no timing available"
    else
    # Loop time of 10712.2 on 40 procs for 50000000 steps with 4000 atoms
    
    sec=`echo "$o" | awk '{print $4}'`
    cores=`echo "$o" | awk '{print $6}'`
    steps=`echo "$o" | awk '{print $9}'`
    atoms=`echo "$o" | awk '{print $12}'`
    hours=`echo $sec 3600 | awk '{print $1/$2}'`
    seccores_per_step=`echo $sec $cores $steps | awk '{print 1000*$1*$2/$3}'`
    seccores_per_step_atoms=`echo $sec $cores $steps $atoms | awk '{print 1000000*$1*$2/$3/$4}'`
    echo "$pfadd || $sec sec || $hours hours || steps: $steps || cores: $cores || $seccores_per_step sec*cores/(1000 steps) || $seccores_per_step_atoms sec*cores/(1000 stesp*atoms)"
    fi
done
