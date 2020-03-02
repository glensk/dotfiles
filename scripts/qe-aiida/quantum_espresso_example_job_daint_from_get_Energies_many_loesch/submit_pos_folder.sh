#!/bin/sh
ff=`ls -1d pos_*`
hier=`pwd`
for f in $ff;do
    cd $hier
    cd $f
    nat=`grep "nat" aiida.in`
    jd=`grep "JOB DONE" aiida.out`

    echo $f `pwd` $nat ":$jd:"
    if [ "$jd" == "" ];then
        echo resubmit
        #sed -i 's|#SBATCH --nodes=.*|#SBATCH --nodes=2|' submit.sh
        #sed -i "s|-nk' 4|-nk' 36|" submit.sh
        rm aiida.out
        rm _*
        sbatch submit.sh



    fi

    #sbatch submit.sh
    cd $hier
done


