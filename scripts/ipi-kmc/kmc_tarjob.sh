#!/bin/sh

files=`find . -name KMC_QCACHE`
for i in $files;do
    echo $i
done

