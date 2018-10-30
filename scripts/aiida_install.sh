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

