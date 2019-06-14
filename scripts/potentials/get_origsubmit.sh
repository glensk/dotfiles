#!/bin/sh
ff=`ls -1d n2p2*`
hier=`pwd`
for f in $ff;do
    echo
    cd $hier
    cd $f
    ll=`ls README_* | wc -l`
    pwd_=`cat README_* | grep "# pwd:"`
    
    echo $f $ll $pwd_
    cd $hier
done


