#!/bin/bash
#SBATCH --no-requeue
#SBATCH --job-name="aiida-103275"
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=18
#SBATCH --time=08:00:00


#SBATCH --constraint=mc

module load daint-mc

module use /users/dmarchan/Install/modules/all/
module load QuantumESPRESSO/6.3-backports-20181003-CrayIntel-18.08
export OMP_NUM_THREADS=1

'srun' '/users/dmarchan/Install/software/QuantumESPRESSO/6.3-backports-20181003-CrayIntel-18.08/bin/pw.x' '-nk' '4' '-in' 'aiida.in'  > 'aiida.out' 
