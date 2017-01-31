#!/bin/sh

if [ "$1" = "" ];then
    if [ -e "log.lammps" ];then 
        pfad=log.lammps
    else
        echo 'need $1 to be pfad to file'
        exit
    fi
else
pfad=$1
fi

o=`tail -200 $pfad | grep "^Loop time of"`
# Loop time of 10712.2 on 40 procs for 50000000 steps with 4000 atoms

sec=`echo "$o" | awk '{print $4}'`
cores=`echo "$o" | awk '{print $6}'`
steps=`echo "$o" | awk '{print $9}'`
atoms=`echo "$o" | awk '{print $12}'`
hours=`echo $sec 3600 | awk '{print $1/$2}'`
#echo "$sec sec || $hours hours || steps: $steps || cores: $cores"
seccores_per_step=`echo $sec $cores $steps | awk '{print 1000*$1*$2/$3}'`
echo "$seccores_per_step sec*cores/(1000 steps)"
