#!/bin/sh

cd $HOME/Dropbox/Albert/git
#git clone https://github.com/CompPhysVienna/n2p2.git
cd n2p2/src
pwd
make libnnpif-shared
cd $HOME/sources/lammps_source_cosmo
ln -s $HOME/Dropbox/Albert/git/n2p2 lib/nnp


