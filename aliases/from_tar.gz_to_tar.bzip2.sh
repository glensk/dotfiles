#!/bin/sh
tmp=$HOME/tmp_tar
[ ! -e "$tmp" ] && mkdir $tmp
hier=`pwd`
files=`find $hier -name "*.tar.gz"`
for i in $files;do
    cd $hier
    [ "`ls $tmp`" != "" ] && echo "$tmp is not empty (start)" && exit
    echo
    filename=`basename $i`
    BASEDIR=$(dirname "$i")
    echo "1 $i"
    echo "2 $filename"
    echo "3 $BASEDIR"
    mv $i $tmp
    cd $tmp
    echo "4 `ls`"
    tar -xvf $filename && rm $filename
    echo "5 `ls`"
    [ -e "$filename" ] && echo "$filename still exist" && exit
    [ "`ls *.tar.gz`" != "" ] && echo "*.tar.gz files exist" && exit
    tar.bzip2_file_or_folder_and_remove_originals.sh *
    echo "6 `ls`"
    mv * $BASEDIR
    cd $hier
    [ "`ls $tmp`" != "" ] && echo "$tmp is not empty (fin)" && exit
    echo "7 done"
done
