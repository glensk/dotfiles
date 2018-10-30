#!/bin/sh
cd $HOME
mkdir aiida
cd aiida
git clone https://github.com/aiidateam/aiida_core
virtualenv/virtualenv.py ~/aiidapy
source ~/aiidapy/bin/activate
pip install -U setuptools pip  # Upgrade all specified packages to the newest available version
cd $HOME/aiida
pip install -e aiida_core  # -e means install in "develop mode"
verdi quicksetup  # albert.glensk@gmail.com
--> it tells me that postgresql is missing ...

############# on mac ##########
conda create -n aiida python=2.7.15
conda activate aiida
conda install postgresql
mkdir $HOME/aiida 
cd $HOME/aiida 
git clone https://github.com/aiidateam/aiida_core
pip install -e aiida_core  # was not available with conda
--> here i can not install the aiida_core also with pip ...

