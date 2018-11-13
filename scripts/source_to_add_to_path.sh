SCR="$( cd "$(dirname "$0")" ; pwd -P )"
#if [ -d "$SCRIPTPATH/ipi-kmc/" ] ; then
#      PATH="$PATH:$SCRIPTPATH/ipi-kmc"
#fi
#
#if [ -d "$SCRIPTPATH/runner_scripts/" ] ; then
#      PATH="$PATH:$SCRIPTPATH/runner_scripts"
#fi

#PATH="$PATH:$SCRIPTPATH/ipi-kmc:$SCRIPTSPATH/runner_scripts"
ipi="$SCR/ipi-kmc"
runner="$SCR/runner_scripts"
aiida_s="$SCR/qe-aiida/aiida_solutejobs_scripts"
aiida_a="$SCR/qe-aiida/aiida_analyze"
lammps1="$SCR/scripts/lammps_executables"
lammps2="$SCR/scripts/lammps_scripts"

PATH="$PATH:$ipi:$runner:$aiida_s:$aiida_a:$lammps1:$lammps2"

