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
_aiida="$HOME/sources/aiida-alloy/"
_aiida_o="$SCR/qe-aiida/"
_aiida_s="$SCR/qe-aiida/aiida_solutejobs_scripts/"
_aiida_a="$SCR/qe-aiida/aiida_analyze/"
_aiida_b="$SCR/qe-aiida/aiida_submitskripts/"
_lammps1="$SCR/lammps_executables"
_lammps2="$SCR/lammps_scripts"
_python_thermodynamics="$SCR/python_thermodynamics/"
#_ase_lammps="$HOME/sources/lammps_source_cosmo" # not necessary when added the _lammps_source/src to $LD_LIBRARY_PATH
_n2p2_lib="$HOME/sources/n2p2/lib"
_ipi_source="$HOME/sources/ipi/"
#_lammps_source="$HOME/sources/lammps_n2p2/"
_lammps_source="$HOME/sources/lammps/"
#addeverywhere="$_python_thermodynamics:$_ipi:$_n2p2:$_runner:$_aiida:$_aiida_o:$_aiida_s:$_aiida_a:$_aiida_b:$_lammps1:$_lammps2:$_ase_lammps/python:$_ipi_source"
addeverywhere="$_python_thermodynamics:$_ipi:$_n2p2:$_runner:$_aiida:$_aiida_o:$_aiida_s:$_aiida_a:$_aiida_b:$_lammps1:$_lammps2:$_ipi_source"

export PATH="$PATH:$addeverywhere"

# the $_lammps_source/src holds the liblammps.so which is necessary to use lammps in ase
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ase_lammps/src:$_n2p2_lib:$_lammps_source/src"

if [ "`hostname`" = "mac" ];then
    _n2p2_mac_predict="$HOME/miniconda2/lib"
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$_n2p2_mac_predict"
    export GSL_ROOT="/Users/glensk/miniconda2/pkgs/gsl-2.4-ha2d443c_1005" # conda install -c conda-forge gsl
    export EIGEN_ROOT="/Users/glensk/miniconda2/" # conda install -c conda-forge gsl
fi

#export LAMMPS_COMMAND="$SCR/executables/lmp_$onhost"
#if [ -e $SCR/executables/lmp_$onhost\_par ];then  # dont use "$SCR/ex..." quotes on fidis --> will not work
#    #echo 'does yes'
#    export LAMMPS_COMMAND="$SCR/executables/lmp_$onhost\_par"
#else
#    #echo 'does not'
#    export LAMMPS_COMMAND="$SCR/executables/lmp_$onhost"
#fi
# @@ on fidis this can do runner and n2p2
export LAMMPS_COMMAND="$SCR/executables/lmp_$onhost"

export IPI_COMMAND="$HOME/sources/ipi/bin/i-pi"
export N2P2_PATH="$HOME/sources/n2p2/"

if [ "$PYTHONPATH" = "" ];then
    export PYTHONPATH="$addeverywhere"              # need "export" for python
else
    export PYTHONPATH="$PYTHONPATH:$addeverywhere"  # need "export" for python
fi
export PYTHONPATH="$PYTHONPATH:$_lammps_source/python"  # lammps with python
export scripts=$SCR
export ESPRESSO_PSEUDO=$SCR/potentials/quantum_espresso/pseudo_SSSPefV1.1_alalloy

# for vmd (seems ok to define this on mac with downloaded/installed vmd)
export VMDINSTALLNAME="vmd"
export VMDINSTALLBINDIR="$HOME/sources/vmd-1.9.3/bin"
export VMDINSTALLLIBRARYDIR="$HOME/sources/vmd-1.9.3/lib/lib" # lib/lib to distinguish from /lib

export AIIDA_PATH="$HOME/aiida/.aiida"
