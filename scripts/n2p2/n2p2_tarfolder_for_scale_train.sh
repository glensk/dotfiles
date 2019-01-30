#!/bin/sh
weights=`find . -maxdepth 1 -name "weights.0*"`
weights_l=`find . -maxdepth 1 -name "weights.0*" | wc -l`
echo :weights:$weights_l:
if [ "$weights_l" != "0" ];then
    mkdir _weights 
    mv $weights _weights
fi

weights=`find . -maxdepth 1 -name "nnp*"`
weights_l=`find . -maxdepth 1 -name "nnp*" | wc -l`
echo :nnp:$weights_l:
if [ "$weights_l" != "0" ];then
    mkdir _nnp
    mv $weights _nnp
fi

weights=`find . -maxdepth 1 -name "sf*"`
weights_l=`find . -maxdepth 1 -name "sf*" | wc -l`
echo :sf:$weights_l:
if [ "$weights_l" != "0" ];then
    mkdir _sf
    mv $weights _sf
fi

for i in `echo "testforces testpoints trainforces trainpoints"`;do
weights=`find . -maxdepth 1 -name "$i*"`
weights_l=`find . -maxdepth 1 -name "$i*" | wc -l`
echo :$i:$weights_l:
if [ "$weights_l" != "0" ];then
    mkdir _test_train
    mv $weights _test_train
fi
done

[ -e "_weights" ] && tar.bzip2_file_or_folder.sh _weights
[ -e "_nnp" ] && tar.bzip2_file_or_folder.sh _nnp
[ -e "_sf" ] && tar.bzip2_file_or_folder.sh _sf
[ -e "_test_train" ] && tar.bzip2_file_or_folder.sh _test_train
[ -e "_nnp-train" ] && tar.bzip2_file_or_folder.sh _nnp-train
[ -e "train-log.out" ] && tar.bzip2_file_or_folder.sh train-log.out

