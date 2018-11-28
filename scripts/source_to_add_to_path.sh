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
lammps1="$SCR/scripts/lammps_executables"
lammps2="$SCR/scripts/lammps_scripts"

addeverywhere="$ipi:$runner:$aiida_s:$aiida_a:$aiida_b:$lammps1:$lammps2"
export PATH="$PATH:$addeverywhere"
if [ "$PYTHONPATH" = "" ];then
    export PYTHONPATH="$addeverywhere"              # need "export" for python
else
    export PYTHONPATH="$PYTHONPATH:$addeverywhere"  # need "export" for python
fi
export scripts=$SCR
#export ipi_mc=$SCR/i-pi-mc/bin/i-pi  
#export ipi_mc=/Users/glensk/tmp_ipi_merge_ipi_kmc/i-pi/bin/i-pi
export ipi_mc=$HOME/Dropbox/Albert/scripts/i-pi-new/i-pi
export ESPRESSO_PSEUDO=$SCR/potentials/quantum_espresso/pseudo_SSSPefV1.1_alalloy
