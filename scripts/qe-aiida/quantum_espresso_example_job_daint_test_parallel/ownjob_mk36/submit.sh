#!/bin/bash
#SBATCH --no-requeue
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=mc

#module purge
module load daint-mc
module load QuantumESPRESSO
module list

export OMP_NUM_THREADS=1

echo -e "\nSLURM_JOBID        =" $SLURM_JOBID
echo -e   "SLURM_JOB_NODELIST =" $SLURM_JOB_NODELIST
echo -e   "SLURM_NNODES       =" $SLURM_NNODES
echo -e   "SLURM_SUBMIT_DIR   =" $SLURM_SUBMIT_DIR
echo -e "\nsystem info:\n" `uname -a` "\n" `head -n 1 /etc/issue` "\n"
echo -e "work dir:\n" `pwd` "\n"
echo -e "ESPRESSO_PSEUDO      =" $ESPRESSO_PSEUDO

#'srun' 'pw.x' '-nk' 2 '-in' 'aiida.in'  > 'aiida.out'
'srun' '/users/dmarchan/Install/software/QuantumESPRESSO/6.3-backports-20181003-CrayIntel-18.08/bin/pw.x' '-nk' 36 '-in' 'aiida.in'  > 'aiida.out'
rm -rf ./out
