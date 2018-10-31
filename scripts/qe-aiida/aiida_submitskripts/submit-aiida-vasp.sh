#!/bin/bash
#SBATCH --no-requeue
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --time=24:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=mc

#module purge
module load daint-mc
module load VASP
module list

echo -e "\nSLURM_JOBID        =" $SLURM_JOBID
echo -e   "SLURM_JOB_NODELIST =" $SLURM_JOB_NODELIST
echo -e   "SLURM_NNODES       =" $SLURM_NNODES
echo -e   "SLURM_SUBMIT_DIR   =" $SLURM_SUBMIT_DIR
echo -e "\nsystem info:\n" `uname -a` "\n" `head -n 1 /etc/issue` "\n"
echo -e "work dir:\n" `pwd` "\n"


'srun'  'vasp_std'   
