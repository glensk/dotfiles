#!/bin/sh

#####################################################################
# get n2p2 
#####################################################################
cd $HOME/Dropbox/Albert/git
[ ! -e "n2p2" ] && git clone https://github.com/CompPhysVienna/n2p2.git
cd n2p2/src
#make libnnpif-static
#make libnnptrain-static
make static  
# in my case /usr/include/eigen3/Eigen/src/Core/util/MKL_support.h has in line 56 a preprocessor command #if defined EIGEN_USE_MKL
# the EIGEN_USE_MKL is actually not acitvaed in the makefiel but is somehow, during the compilation set (which makes the compilation crash)

#####################################################################
# now the part to make lammps work with n2p2
#####################################################################

# from https://compphysvienna.github.io/n2p2/a01203.html
pwd
make libnnpif-shared
cd $HOME/sources/lammps_source_cosmo
ln -s $HOME/Dropbox/Albert/git/n2p2 lib/nnp
cd $HOME/Dropbox/Albert/git/n2p2
cp -r src/interface/LAMMPS/src/USER-NNP $HOME/sources/lammps_source_cosmo/src

cd $HOME/sources/lammps_source_cosmo/src
make yes-user-nnp
make mpi

