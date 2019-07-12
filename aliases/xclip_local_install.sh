#!/bin/bash
 
# Script for installing xclip on systems without root access.
# xclip will be installed in $HOME/.local/bin.
# It's assumed that wget and a C/C++ compiler are installed.
 
# exit on error
set -e
 
XCLIP_VERSION=0.12
 
# create our directories
mkdir -p $HOME/.local $HOME/sources/xclip_tmp
cd $HOME/sources/xclip_tmp
 
# download source files for XClip
wget http://kent.dl.sourceforge.net/project/xclip/xclip/${XCLIP_VERSION}/xclip-${XCLIP_VERSION}.tar.gz
 
# extract files, configure, and compile
 
tar xvzf xclip-${XCLIP_VERSION}.tar.gz
cd xclip-${XCLIP_VERSION}
./configure --prefix=$HOME/.local --disable-shared
make
make install
