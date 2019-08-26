#!/bin/sh

### 6.3
export QEDIR=$HOME/sources/
cd $QEDIR
#wget https://gitlab.com/QEF/q-e/-/archive/qe-6.3/q-e-qe-6.3.tar.gz
#tar -xzvf q-e-qe-6.3.tar.gz
#mv q-e-qe-6.3 qe-6.3
#cd $QEDIR/qe-6.3
wget https://gitlab.com/QEF/q-e/-/archive/29c951ebbfc3b2248d0308102d05c4c081db891e/q-e-29c951ebbfc3b2248d0308102d05c4c081db891e.tar.gz
tar -xzvf q-e-29c951ebbfc3b2248d0308102d05c4c081db891e.tar.gz
cd q-e-29c951ebbfc3b2248d0308102d05c4c081db891e

module load intel intel-mkl intel-mpi

#CC=mpiicc FC=mpiifort  ./configure --enable-openmp --enable-shared --with-scalapack=no    # orig felix/edgar
#CC=mpiicc FC=mpiifort  ./configure --enable-openmp --enable-shared --with-scalapack=no    # works
#CC=mpiicc FC=mpiifort  ./configure --enable-openmp --enable-parallel --with-scalapack --disable-xml
CC=mpiicc FC=mpiifort  ./configure --enable-openmp --enable-parallel --with-scalapack=intel --disable-xml

make pw
