#!/bin/sh

# from https://compphysvienna.github.io/n2p2/a01203.html

cd $HOME/Dropbox/Albert/git
#git clone https://github.com/CompPhysVienna/n2p2.git
cd n2p2/src
pwd
make libnnpif-shared
cd $HOME/sources/lammps_source_cosmo
ln -s $HOME/Dropbox/Albert/git/n2p2 lib/nnp
cd $HOME/Dropbox/Albert/git/n2p2
cp -r src/interface/LAMMPS/src/USER-NNP $HOME/sources/lammps_source_cosmo/src

cd $HOME/sources/lammps_source_cosmo/src
make yes-user-nnp
make mpi

