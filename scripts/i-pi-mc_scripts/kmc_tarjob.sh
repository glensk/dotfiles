#!/bin/sh
hier=`pwd`

cd $hier
files=`find . -name KMC_QCACHE`
for i in $files;do
    folder=`echo $i | sed 's|/KMC_QCACHE$||g'`
    cd $hier
    cd $folder
    pwd
    du -sh KMC_QCACHE
    tar.gz_file_or_folder_remove_origfiles.sh KMC_QCACHE
    cd $hier
done

cd $hier
files=`find . -name log.i-pi`
for i in $files;do
    folder=`echo $i | sed 's|/log.i-pi$||g'`
    cd $hier
    cd $folder
    pwd
    tar.gz_file_or_folder_remove_origfiles.sh log*
    cd $hier
done

