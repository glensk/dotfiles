#!/bin/sh
tmp=$HOME/tmp_tar
[ ! -e "$tmp" ] && mkdir $tmp
hier=`pwd`
files=`find $hier -name "*.tar.gz"`
#files=`find $hier -name "OUTCAR.gz"`
#files=`find $hier -name "vasprun.xml.gz"`
for i in $files;do
    echo $i
done
#exit
outof=`echo "$files" | wc -l`

count=1
for i in $files;do
    cd $hier
    [ "`ls $tmp`" != "" ] && echo "$tmp is not empty (start)" && exit
    echo
    filename=`basename $i`
    BASEDIR=$(dirname "$i")
    echo "COUONT $count out of $outof"
    count=`echo $count | awk '{print $1+1}'`
    echo "1 i       : $i"
    echo "2 filename: $filename"
    echo "3 BASEDIR : $BASEDIR"
    mv $i $tmp
    cd $tmp
    echo "4 ls      : `ls`"
    [ "`ls`" == "" ] && echo empty 1 && exit
    [ "`ls $filename`" != "$filename" ] && echo filename 1 && exit
    tar -xvf $filename && rm $filename
    echo "5         : `ls`"
    [ -e "$filename" ] && echo "$filename still exist" && exit
    [ "`ls *.tar.gz`" != "" ] && echo "*.tar.gz files exist" && exit
    tar.bzip2_file_or_folder_and_remove_originals.sh *
    echo "6 `ls`"
    mv * $BASEDIR
    cd $hier
    [ "`ls $tmp`" != "" ] && echo "$tmp is not empty (fin)" && exit
    echo "7 done"
done
