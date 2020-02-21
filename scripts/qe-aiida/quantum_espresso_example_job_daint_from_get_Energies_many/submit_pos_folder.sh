#!/bin/sh
ff=`ls -1d pos_*`
hier=`pwd`
for f in $ff;do
    cd $hier
    cd $f
    echo $f `pwd`
    sbatch submit.sh
    cd $hier
done


