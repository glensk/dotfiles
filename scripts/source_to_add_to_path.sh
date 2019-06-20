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
#_ipi="$SCR/i-pi-mc_scripts"
#_runner="$SCR/runner_scripts"
#_n2p2="$SCR/n2p2"
#_aiida="$HOME/sources/aiida-alloy/"
#_aiida_o="$SCR/qe-aiida/"
#_aiida_s="$SCR/qe-aiida/aiida_solutejobs_scripts/"
#_aiida_a="$SCR/qe-aiida/aiida_analyze/"
#_aiida_b="$SCR/qe-aiida/aiida_submitskripts/"
#_lammps_scripts="$SCR/lammps_scripts"
_python_thermodynamics="$SCR/python_thermodynamics/"
#_ase_lammps="$HOME/sources/lammps_source_cosmo" # not necessary when added the _lammps_source/src to $LD_LIBRARY_PATH
_n2p2_lib="$HOME/sources/n2p2/lib"
_ipi_source="$HOME/sources/ipi/"
export LAMMPSPATH="$HOME/sources/lammps"
_bin="$SCR/bin"
_bin2="$SCR/../aliases"
#addeverywhere="$_python_thermodynamics:$_ipi:$_n2p2:$_runner:$_aiida:$_aiida_o:$_aiida_s:$_aiida_a:$_aiida_b:$_lammps_exec:$_lammps_scripts:$_ase_lammps/python:$_ipi_source"
#addeverywhere="$_python_thermodynamics:$_ipi:$_n2p2:$_runner:$_aiida:$_aiida_o:$_aiida_s:$_aiida_a:$_aiida_b:$_lammps_exec:$_lammps_scripts:$_ipi_source"
#addeverywhere="$_bin:$_aiida:$_ipi_source"  # add aiida only in aiases, since I have changed my own stuff.
addeverywhere="$_bin:$_bin2:$_ipi_source"  # add aiida only in aiases, since I have changed my own stuff.

#####################
##### PATH
#####################
export PATH="$PATH:$addeverywhere"

#####################
##### LD_LIBRARY_PATH
#####################
# the $LAMMPSPATH/src holds the liblammps.so which is necessary to use lammps in ase
MY_LD_LIBRARY_PATH="$_n2p2_lib:$LAMMPSPATH/src"
if [ "$LD_LIBRARY_PATH" = "" ];then
    export LD_LIBRARY_PATH="$MY_LD_LIBRARY_PATH"
else
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MY_LD_LIBRARY_PATH"
fi

#####################
##### PYHON_PATH
#####################
if [ "$PYTHONPATH" = "" ];then
    export PYTHONPATH="$addeverywhere"              # need "export" for python
else
    export PYTHONPATH="$PYTHONPATH:$addeverywhere"  # need "export" for python
fi
export PYTHONPATH="$PYTHONPATH:$LAMMPSPATH/python"  # lammps with python



if [ "`hostname`" = "mac" ];then
    _n2p2_mac_predict="$HOME/miniconda2/lib"
    #export PATH="$PATH:$HOME/local/bin"    # add mpic++ (MPI) which is necessary to compile n2p2 ... but does not work for compilation of n2p2
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
export scripts=$SCR
export ESPRESSO_PSEUDO=$SCR/potentials/quantum_espresso/pseudo_SSSPefV1.1_alalloy
export OMP_NUM_THREADS=1

# for vmd (seems ok to define this on mac with downloaded/installed vmd)
export VMDINSTALLNAME="vmd"
export VMDINSTALLBINDIR="$HOME/sources/vmd-1.9.3/bin"
export VMDINSTALLLIBRARYDIR="$HOME/sources/vmd-1.9.3/lib/lib" # lib/lib to distinguish from /lib

export AIIDA_PATH="$HOME/aiida/.aiida"
