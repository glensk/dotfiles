if [ "$ZSH_VERSION" != "" ];then
    #echo zsh
    SCR="$( cd "$(dirname "$0")" ; pwd -P )"
else
    #echo bash
    SCR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
fi

ipi="$SCR/i-pi-mc_scripts"
runner="$SCR/runner_scripts"
aiida_s="$SCR/qe-aiida/aiida_solutejobs_scripts"
aiida_a="$SCR/qe-aiida/aiida_analyze"
aiida_b="$SCR/qe-aiida/aiida_submitskripts"
lammps1="$SCR/lammps_executables"
lammps2="$SCR/lammps_scripts"
ase_lammps="$HOME/sources/lammps_source_cosmo"

addeverywhere="$ipi:$runner:$aiida_s:$aiida_a:$aiida_b:$lammps1:$lammps2:$ase_lammps/python"

export PATH="$PATH:$addeverywhere"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$ase_lammps/src"
#export LAMMPS_COMMAND="$SCR/executables/lmp_fidis_par_2018_11_28"
export LAMMPS_COMMAND="$SCR/executables/lmp_$onhost"
#export lmp_exec="$SCR/executables/lmp_fidis_par_2018_11_28"

if [ "$PYTHONPATH" = "" ];then
    export PYTHONPATH="$addeverywhere"              # need "export" for python
else
    export PYTHONPATH="$PYTHONPATH:$addeverywhere"  # need "export" for python
fi
export scripts=$SCR
export ipi_mc=$HOME/Dropbox/Albert/scripts/i-pi-new/i-pi
export ESPRESSO_PSEUDO=$SCR/potentials/quantum_espresso/pseudo_SSSPefV1.1_alalloy
