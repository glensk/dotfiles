
if [ "$ZSH_VERSION" != "" ];then
    #echo zsh
    SCR="$( cd "$(dirname "$0")" ; pwd -P )"
else
    #echo bash
    SCR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
fi
#echo "SCR: $SCR"
#if [ -d "$SCRIPTPATH/ipi-kmc/" ] ; then
#      PATH="$PATH:$SCRIPTPATH/ipi-kmc"
#fi
#
#if [ -d "$SCRIPTPATH/runner_scripts/" ] ; then
#      PATH="$PATH:$SCRIPTPATH/runner_scripts"
#fi

#PATH="$PATH:$SCRIPTPATH/ipi-kmc:$SCRIPTSPATH/runner_scripts"
ipi="$SCR/i-pi-mc_scripts"
runner="$SCR/runner_scripts"
aiida_s="$SCR/qe-aiida/aiida_solutejobs_scripts"
aiida_a="$SCR/qe-aiida/aiida_analyze"
lammps1="$SCR/scripts/lammps_executables"
lammps2="$SCR/scripts/lammps_scripts"

PATH="$PATH:$ipi:$runner:$aiida_s:$aiida_a:$lammps1:$lammps2"
if [ "$PYTHONPATH" = "" ];then
    PYTHONPATH="$ipi:$runner:$aiida_s:$aiida_a"
else
    PYTHONPATH="$PYTHONPATH:$ipi:$runner:$aiida_s:$aiida_a"
fi
export scripts=$SCR
export ipi_mc=$SCR/i-pi-mc/bin/i-pi   
