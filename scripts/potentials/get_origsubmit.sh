#!/bin/sh
ff=`ls -1d n2p2*`
hier=`pwd`
for f in $ff;do
    echo
    cd $hier
    cd $f
    [ -e "submit_n2p2_train.sh" ] && continue #echo "-->" $f exists && continue
    #ll=`ls README_* | wc -l`
    ll=`find . -name "README_*" | wc -l`
    [ "$ll" == "0" ] && continue
    pwd_=`cat README_* | grep "# pwd:" | awk '{print $3}' | sed 's|/potential$||'`
    [ -e "$pwd_/submit_training.sh" ] && cp $pwd_/submit_training.sh submit_n2p2_train.sh
    [ -e "$pwd_/submit_n2p2_train.sh" ] && cp $pwd_/submit_n2p2_train.sh submit_n2p2_train.sh
    
    echo "-->" $f :$ll:$ll: $pwd_
    cd $hier
done


