#!/bin/sh
hier=`pwd`

tar="KMC_QCACHE"
files=`find . -name $tar`
for i in $files;do
    cd $hier
    folder=`echo $i | sed 's|/'"$tar"'$||g'`
    cd $hier
    cd $folder
    pwd
    #du -sh $tar 
    tar.bzip2_file_or_folder.sh $tar
    cd $hier
done
tar="KMC_ECACHE"
files=`find . -name $tar`
for i in $files;do
    cd $hier
    folder=`echo $i | sed 's|/'"$tar"'$||g'`
    cd $hier
    cd $folder
    pwd
    #du -sh $tar 
    tar.bzip2_file_or_folder.sh $tar
    cd $hier
done

# all log files
files=`find . -name log.i-pi`
for i in $files;do
    cd $hier
    folder=`echo $i | sed 's|/log.i-pi$||g'`
    cd $hier
    cd $folder
    pwd
    tar.bzip2_file_or_folder.sh log*
    cd $hier
done

