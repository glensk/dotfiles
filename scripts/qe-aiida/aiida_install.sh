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
verdi quicksetup  # albert.glensk@gmail.com,Albert,Glensk,EPFL COSMO, 

Password:
Detected no known postgres setup, some information is needed to create the aiida database and grant
aiida access to it. If you feel unsure about the following parameters, first check if postgresql is
installed. If postgresql is not installed please exit and install it, then run verdi quicksetup again.
If postgresql is installed, please ask your system manager to provide you with the following parameters:
    postgres host [localhost]:
    postgres port [5432]:
    template [template1]:
    postgres super user [postgres]:
PW: Standard tow nix
%verdi profile list
%Info: configuration folder: /Users/glensk/.aiida
%Critical: no default profile configured yet, run `verdi setup`


