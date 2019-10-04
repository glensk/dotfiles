#!/bin/sh
ff=`ls -1d *`
hier=`pwd`
ff=" 1   2   3   5   7  10  15  21  26  30  42  60  84 119 168 236 332 468 658 659"
[ -e "result_bill" ] && rm -f result_bill
for f in $ff;do
    #cd $hier
    echo $f

    getEnergies_byLammps.py -p ../ --units meV_pa -i train.data -wa -pe $f #-idx :6
    std=`tail -1 ene_std.npy`
    echo $f $std >> result_bill 
    mv ene_std.npy ene_std.npy_$f
    #cd $hier
done


