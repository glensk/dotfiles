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
conda install click  # actually comes with ase
conda install -c conda-forge ase
conda install seaborn
[ "`hostname`" = "mac" ] && conda install -c conda-forge jupyter_contrib_nbextensions   # get the notebook extensions for jupyter notebooks
[ "`hostname`" = "mac" ] && jupyter contrib nbextension install --user                  # also necessary to get the notebook extensions working

#conda install -c conda-forge pyfftw=0.10.4  # is this already installed using 
#conda install argcomplete           # to get argcompletion of python scripts in bash/zsh

