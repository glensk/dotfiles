#!/bin/sh

# check if BASH or ZSH
if [ "$ZSH_VERSION" != "" ];then
    #echo zsh
    SCR="$( cd "$(dirname "$0")" ; pwd -P )"
else
    #echo bash
    SCR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
fi

# define the folders which are given in this git repo 
_ipi="$SCR/i-pi-mc_scripts"
_runner="$SCR/runner_scripts"
_n2p2="$SCR/n2p2"
_aiida_o="$SCR/qe-aiida/"
_aiida_s="$SCR/qe-aiida/aiida_solutejobs_scripts/"
_aiida_a="$SCR/qe-aiida/aiida_analyze/"
_aiida_b="$SCR/qe-aiida/aiida_submitskripts/"
_lammps1="$SCR/lammps_executables"
_lammps2="$SCR/lammps_scripts"
_ase_lammps="$HOME/sources/lammps_source_cosmo"
_n2p2_lib="$HOME/Dropbox/Albert/git/n2p2/lib"

addeverywhere="$_ipi:$_n2p2:$_runner:$_aiida_o:$_aiida_s:$_aiida_a:$_aiida_b:$_lammps1:$_lammps2:$_ase_lammps/python"

export PATH="$PATH:$addeverywhere"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ase_lammps/src:$_n2p2_lib"

export LAMMPS_COMMAND="$SCR/executables/lmp_$onhost"
export IPI_COMMAND="$HOME/sources/ipi/bin/i-pi"
export N2P2_PATH="$HOME/sources/n2p2/"

if [ "$PYTHONPATH" = "" ];then
    export PYTHONPATH="$addeverywhere"              # need "export" for python
else
    export PYTHONPATH="$PYTHONPATH:$addeverywhere"  # need "export" for python
fi
export scripts=$SCR
export ESPRESSO_PSEUDO=$SCR/potentials/quantum_espresso/pseudo_SSSPefV1.1_alalloy

# for vmd
export VMDINSTALLNAME="vmd"
export VMDINSTALLBINDIR="$HOME/sources/vmd-1.9.3/bin"
export VMDINSTALLLIBRARYDIR="$HOME/sources/vmd-1.9.3/lib/lib" # lib/lib to distinguish from /lib
