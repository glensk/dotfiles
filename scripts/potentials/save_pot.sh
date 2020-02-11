#!/bin/sh
cd /scratch/glensk/2020_retrain_almgsi/zeroth20

ff=`ls -1d random_seed_*`
hier=`pwd`
for f in $ff;do
    cd $hier
    seed=`echo $f | sed 's|random_seed_||'`
    cd $f
    fn="/home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/n2p2_v6ag_$seed"
    mkdir -p $fn
    #nr=`ls weights.012.000* | xargs -n 1 | sort -n | tail -1 | sed "s|weights.012.||" | sed "s|.out||"`

    echo $f $seed $nr
    #cp weights.012.$nr.out $fn
    #cp weights.013.$nr.out $fn
    #cp weights.014.$nr.out $fn
    #cp learning-curve.out $fn
    #cp input.nn $fn
    #cp input.data $fn
    #cp train.data $fn
    #cp test.data $fn
    cp submit_training_debug.sh $fn
    cp -r analysis_0thrun $fn
    cp -r analysis_8thrun $fn
    cp -r analysis_OUTLIERS $fn
    #cp scaling.data $fn
    cd $hier
done


