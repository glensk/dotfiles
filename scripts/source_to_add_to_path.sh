# check if BASH or ZSH
if [ "$ZSH_VERSION" != "" ];then
    #echo zsh
    SCR="$( cd "$(dirname "$0")" ; pwd -P )"
else
    #echo bash
    SCR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
fi

# define the folders which are given in this git repo 
ipi="$SCR/i-pi-mc_scripts"
runner="$SCR/runner_scripts"
aiida_s="$SCR/qe-aiida/aiida_solutejobs_scripts"
aiida_a="$SCR/qe-aiida/aiida_analyze"
aiida_b="$SCR/qe-aiida/aiida_submitskripts"
lammps1="$SCR/lammps_executables"
lammps2="$SCR/lammps_scripts"
ase_lammps="$HOME/sources/lammps_source_cosmo"
n2p2="$HOME/Dropbox/Albert/git/n2p2/lib"

addeverywhere="$ipi:$runner:$aiida_s:$aiida_a:$aiida_b:$lammps1:$lammps2:$ase_lammps/python"

export PATH="$PATH:$addeverywhere"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ase_lammps/src:$n2p2"

export LAMMPS_COMMAND="$SCR/executables/lmp_$onhost"
export IPI_COMMAND="$HOME/Dropbox/Albert/git/i-pi-mc/bin/i-pi"

if [ "$PYTHONPATH" = "" ];then
    export PYTHONPATH="$addeverywhere"              # need "export" for python
else
    export PYTHONPATH="$PYTHONPATH:$addeverywhere"  # need "export" for python
fi
export scripts=$SCR
export ESPRESSO_PSEUDO=$SCR/potentials/quantum_espresso/pseudo_SSSPefV1.1_alalloy
