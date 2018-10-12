#!/bin/sh
folder=sed-4.2
rm -rf $folder
tar -xvf sed-4.2.tar.gz
cd $folder
mkdir build && cd build
../configure --prefix=`pwd` --program-suffix=-4.2
make -j 8
sudo make install
