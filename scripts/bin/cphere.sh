#!/bin/sh
ff=`ls -1d *`
hier=`pwd`
for f in $ff;do
    cd $hier
    lk=`ls -l $f | awk '{print $NF}'`
    [ "$lk" = "$f" ] && echo "---> " $lk && continue
    echo $lk $f
    unlink $f
    mv $lk .
    cd $hier
done


