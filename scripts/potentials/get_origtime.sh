#!/bin/sh
ff=`ls -1d n2p2*`
hier=`pwd`
for f in $ff;do
    echo
    cd $hier
    cd $f
    [ -e "time.sec" ] && continue #echo "-->" $f exists && continue
    #ll=`ls README_* | wc -l`
    ll=`find . -name "README_*" | wc -l`
    [ "$ll" == "0" ] && continue
    pwd_=`cat README_* | grep "# pwd:" | awk '{print $3}' | sed 's|/potential$||'`
    [ -e "$pwd_/time.sec" ] && cp $pwd_/time.sec time.sec 
    
    echo "-->" $f :$ll:$ll: $pwd_
    cd $hier
done


