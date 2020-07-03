#!/bin/sh
ff=`ls -1d n2p2_v7ag_neinth20_*`
hier=`pwd`
for f in $ff;do
    cd $hier
    cd $f
    echo $f
    
    #if [ ! -e "test_auswertung/ene_diff_abs.npy" ];then
    #mkdir test_auswertung
    #cd test_auswertung
    #getEnergies_byLammps.py --units meV_pa -p .. -i ../test.data -wa
    #cd ..
    #fi
    
    #if [ ! -e "train_auswertung/ene_diff_abs.npy" ];then
    #mkdir train_auswertung
    #cd train_auswertung
    #getEnergies_byLammps.py --units meV_pa -p .. -i ../train.data -wa
    #cd ..
    #fi

    if [ ! -e "test_auswertung_dE_vs_dist/ene_diff_abs_vs_dist.npy" ];then
    mkdir test_auswertung_dE_vs_dist
    cd test_auswertung_dE_vs_dist
    getEnergies_byLammps.py --units meV_pa -p .. -i ../test.data -wa -waf
    cd ..
    fi

    cd $hier
done


