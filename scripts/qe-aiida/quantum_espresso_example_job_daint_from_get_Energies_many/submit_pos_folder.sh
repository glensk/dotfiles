#!/bin/sh
ff=`ls -1d pos_*`
hier=`pwd`
for f in $ff;do
    # add --mem=120GB
    cd $hier
    cd $f
    [ ! -e "aiida.out" ] && echo "$f does not have aiida.out" && continue
    nat=`grep "nat" aiida.in`
    jd=`grep "JOB DONE" aiida.out`
    mem=`grep "Estimated total dynamical RAM" aiida.out | awk '{print $6}'`
    nodes=`echo $mem | awk '{print $1/120}'`
    #[ "$nodes" -se "1" ] && nodes=1
    nodesround=`echo $nodes | cut -c1-2 |  sed 's|\.||' | awk '{print $1+1}'` #awk 'BEGIN { printf("%.1f\n", $1+0.5); }'`
    nodesused=`grep "MPI processes distributed on" aiida.out | awk '{print $5}'`
    nowuse=$nodesround
    [ "$nodesused" == "$nowuse" ] && nowuse=`echo $nodesused | awk '{print $1+1}'`

    #echo $f `pwd` $nat ":$jd: mem:$mem: GB"
    #if [ "$jd" == "   JOB DONE." ];then
    #    continue
    #fi
    echo $nat ":$jd: mem:$mem: GB nodes:$nodes: nodesround: $nodesround used:$nodesused: nowuse:$nowuse:" $f
    getEnergies_byLammps.py -p runner_v3ag_4998_3 -i aiida.out -fi espresso-out --units hartree -wrDFT


    #if [ "$jd" == "" ];then
    #    cp ../submit.sh .
    #    sed -i 's|#SBATCH --nodes=.*|#SBATCH --nodes='"$nowuse"'|' submit.sh
    #    sed -i "s|-nk' 4|-nk' 36|" submit.sh
    #    sed -i "s|-nk' 2|-nk' 36|" submit.sh
    #    
    #    rm aiida.out
    #    rm _*
    #    sbatch submit.sh
    #fi

    #sbatch submit.sh
    cd $hier
done


