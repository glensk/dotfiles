#!/bin/bash
#SBATCH --no-requeue
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=00:30:00

export PATH="$PATH:/home/dmarchan/Src/ase"
source /home/dmarchan/virtualenvs/ase/bin/activate
module purge
module load intel intel-mpi intel-mkl
module load quantum-espresso
module list

export OMP_NUM_THREADS=1

QE_COMMAND="srun pw.x -in aiida.in"

echo -e "\nSLURM_JOBID        =" $SLURM_JOBID
echo -e   "SLURM_JOB_NODELIST =" $SLURM_JOB_NODELIST
echo -e   "SLURM_NNODES       =" $SLURM_NNODES
echo -e   "SLURM_SUBMIT_DIR   =" $SLURM_SUBMIT_DIR
echo -e "\nsystem info:\n" `uname -a` "\n" `head -n 1 /etc/issue` "\n"
echo -e "work dir:\n" `pwd` "\n"

function modify_setting(){
  INPUT="$1"
  SETTING_TO_CHANGE="$2"
  SETTING_NEWVALUE="$3"
  QE_TAG="$4"
  grep -q $SETTING_TO_CHANGE $INPUT \
    && sed -i s/${SETTING_TO_CHANGE}.*/"${SETTING_NEWVALUE}"/g $INPUT \
    ||  sed -i  "/${QE_TAG}/a ${SETTING_NEWVALUE}" $INPUT
}

function modify_kpoints(){
  INPUT_FILE=$1
  KMESH_L=$2

  sed -i '/K_POINTS/,+1 d' $INPUT_FILE

  echo "K_POINTS  automatic" >> $INPUT_FILE
  gen_qe_kmesh.py $INPUT_FILE $KMESH_L >> $INPUT_FILE || exit 1
}

function runqe(){
  QE_COMMAND="$1"
  RUNMODE="$2"

  FAILED_FLAGS=$(ls ./* | grep failed)
  if [ ! -z $FAILED_FLAGS ]; then
    echo $FAILED_FLAGS found exiting!
    exit 1
  fi


  echo $QE_COMMAND
  OUTPUT_FILE=aiida.out_${RUNMODE}
  $QE_COMMAND > $OUTPUT_FILE

  FINISHED=$(grep "JOB DONE" $OUTPUT_FILE)
  if [ -z "$FINISHED" ]; then
    touch ${FLAGPREFIX}.failed
    echo "Stage has failed to converge!"
    exit 1
  else
    touch ${FLAGPREFIX}.finished
  fi
}

RUNMODE=low
cp aiida.in_final aiida.in || exit 1
modify_setting aiida.in restart_mode restart_mode='from_scratch' "&CONTROL"
modify_kpoints aiida.in 40
runqe "$QE_COMMAND" "$RUNMODE"

RUNMODE=final
cp aiida.in_final aiida.in || exit 1
modify_setting aiida.in restart_mode restart_mode='restart' "&CONTROL"
modify_setting aiida.in startingwfc startingwfc='atomic+random' "&ELECTRONS"
runqe "$QE_COMMAND" "$RUNMODE"
rm -rf ./out

date
DURATION=$SECONDS
echo -e "\nelapsed time (second):" $DURATION
echo "$(($DURATION / 3600)) h : $(( ($DURATION % 3600) / 60 )) m : $(( ($DURATION % 3600) % 60 )) s"
