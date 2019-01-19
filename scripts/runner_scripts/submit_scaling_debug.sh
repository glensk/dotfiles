#!/bin/bash
#SBATCH --job-name=NNP-mpi
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=1
#SBATCH --ntasks 28
#SBATCH --time=00-01:00:00
#SBATCH --constraint=E5v4
#SBATCH --mem=100G

set +e
# it is necessary to have all the modules which are used when compiling
export LD_LIBRARY_PATH=""
module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen
export LD_LIBRARY_PATH=$HOME/Dropbox/Albert/git/n2p2/lib:${LD_LIBRARY_PATH}
#echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH

touch time.out
date +%s >> time.out

srun -n 28 $HOME/Dropbox/Albert/git/n2p2/bin/nnp-scaling 1
mkdir _nnp && mv nnp* _nnp
mkdir _sf && mv sf* _sf
date +%s >> time.out
cat time.out | xargs | awk '{print $2-$1-10}' > time.sec

exit 0
