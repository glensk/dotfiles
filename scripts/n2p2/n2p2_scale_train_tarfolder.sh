#!/bin/sh
mkdir _weights 
mv weights.0* _weights
mkdir _nnp && mv nnp* _nnp
mkdir _sf && mv sf* _sf
mkdir _test_train
mv testforces* _test_train
mv testpoints* _test_train
mv trainforces* _test_train
mv trainpoints* _test_train
tar.bzip2_file_or_folder.sh _weights
tar.bzip2_file_or_folder.sh _test_train
tar.bzip2_file_or_folder.sh _nnp-train
tar.bzip2_file_or_folder.sh train-log.out

