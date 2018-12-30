#!/bin/bash

#SBATCH --job-name=NNP-mpi
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=2
#SBATCH --ntasks 56
#SBATCH --time=00-71:50:00
#SBATCH --constraint=E5v4

set +e
source $MODULESHOME/init/bash    # necessary in the case of zsh or other init shells
module load intel intel-mpi intel-mkl fftw python/2.7.14
export OMP_NUM_THREADS=1
#touch time.out
#date +"%y.%m.%d %H:%M:%S" >> time.out

# sets up the internet socket for connections both for i-PI and on the lammps side
sed -i 's/<ffsocket.*/<ffsocket name="lmpserial" mode="inet">/' input-runner.xml
sed -i 's/address>.*<.addr/address>'$(hostname)'<\/addr/' input-runner.xml
sed -i 's/all ipi [^ ]*/all ipi '$(hostname)'/' in.lmp

python /home/glensk/scripts/i-pi-mc/bin/i-pi input-runner.xml &> log.i-pi &

sleep 10

for i in `seq 4`
do
      srun --hint=nomultithread --exclusive -n 14 --mem=4G /home/glensk/scripts/lammps/src/lmp_fidis < in.lmp > log.lmp$i  &
done

wait 
#date +"%y.%m.%d %H:%M:%S" >> time.out
exit 0
