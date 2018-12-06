#!/bin/bash

#SBATCH --job-name=NNP-mpi
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=1
#SBATCH --ntasks 28
#SBATCH --time=00-71:50:00
#SBATCH --constraint=E5v4

set +e
source $MODULESHOME/init/bash    # necessary in the case of zsh or other init shells
module load intel intel-mpi intel-mkl gsl eigen 
export OMP_NUM_THREADS=1
touch time.out
date +"%y.%m.%d %H:%M:%S" >> time.out


#srun --hint=nomultithread --exclusive -n 14 /home/glensk/scripts/lammps/src/lmp_fidis
srun -n 21 /home/glensk/Dropbox/Albert/git/n2p2/bin/nnp-train

date +"%y.%m.%d %H:%M:%S" >> time.out
exit 0
