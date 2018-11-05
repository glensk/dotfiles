#!/bin/sh

[ ! -e "$HOME/Downloads" ] && mkdir $HOME/Downloads 
[ ! -e "$HOME/Downloads" ] && echo not folder $HOME/Downloads && exit
cd $HOME/Downloads
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh


source ~/.zshrc
conda config --add channels intel 
conda activate
conda install colorama
conda install -c conda-forge ase

#conda install -c conda-forge pyfftw=0.10.4  # is this already installed using 
#conda install argcomplete           # to get argcompletion of python scripts in bash/zsh

