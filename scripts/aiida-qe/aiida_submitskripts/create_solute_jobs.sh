#!/usr/bin/env bash
rm -f aiida.pos
rm -rf /scratch/snx1600/aglensk/users/aglensk/test
#unlink calculations
cp ~/raw_data/second_test/aiida.settings .

exit
pseudo_dir="$HOME/Dropbox/Albert/scripts/dotfiles/scripts/aiida-qe/qe_pseudopotentials"
qe_pos_setting_merge="qe_pos_setting_merge.sh"              # needs to be in $PATH 
generate_structure="generate_structure.py"                  # needs to be in $PATH
virtualenvs_ase="$HOME/virtualenvs/ase/bin/activate"
SETTING_FILE=aiida.settings

DEBUG_MODE=0
PRODUCTION_MODE=0

LATTICE=4.057
SUPERCELL_SIZE="3,3,3"
MATRIX="Al"
SOLUTE="Si"

POSITION_VALUES="pure singlesolute "$(seq 1 1 1)


#SUBMIT_SCRIPT=_multistagesubmit.sh
SUBMIT_SCRIPT=submit-aiida-qe.sh 
#POSITION_FILE=../tmp/aiida.pos
POSITION_FILE=aiida.pos
INPUT_FILE=aiida.in




#########################################################################################
qe_pos_setting_merge=`which $qe_pos_setting_merge`
generate_structure=`which $generate_structure`
[ ! -f $qe_pos_setting_merge ] && echo $qe_pos_settings_merge not found && exit
[ ! -f $generate_structure ] && echo $generate_structure not found && exit

[ ! -f $virtualenvs_ase ] && echo $virtualenvs_ase not found && exit
[ ! -e $pseudo_dir ] && echo $pseudo_dir not found && exit
[ ! -f $SETTING_FILE ] && echo $SETTING_FILE not found && exit
[ ! -f $SUBMIT_SCRIPT ] && echo $SUBMIT_SCRIPT not found && exit
ntypes=`grep ntyp $SETTING_FILE | awk '{print $3}'`
elements=`grep ATOMIC_SPECIES $SETTING_FILE -A $ntypes | tail -$ntypes | awk '{print $1}'`


if [ "$1" == "debug" ]; then
  source ~/.bashrc # sbatch_debug is a custom command that must be in the bashrc
  DEBUG_MODE=1
elif [ "$1" == "production" ]; then
  PRODUCTION_MODE=1
elif [ "$1" == "dryrun" ]; then
  DRYRUN_MODE=1
else
  USAGE="USAGE: $0 <debug/dryrun/production>"
  echo $USAGE
  exit 1
fi

source $virtualenvs_ase
SCRATCH_BASE=${SCRATCH}${PWD}/calculations
mkdir -p $SCRATCH_BASE
if [ ! -e calculations ]; then
    ln -s $SCRATCH_BASE calculations
fi

for POSITION in $POSITION_VALUES;
  do
  GEN_STRUCTURE_CMD="python $generate_structure $LATTICE $SUPERCELL_SIZE $MATRIX $SOLUTE $POSITION $POSITION_FILE"
  #exit 
  #GEN_STRUCTURE_CMD=${GEN_STRUCTURE_CMD}" $LATTICE $SUPERCELL_SIZE "
  #GEN_STRUCTURE_CMD=${GEN_STRUCTURE_CMD}" $MATRIX $SOLUTE $POSITION  "
  #GEN_STRUCTURE_CMD=${GEN_STRUCTURE_CMD}" $POSITION_FILE"
  #exit
  LABEL=$($GEN_STRUCTURE_CMD)
  #exit
  WORKDIR=${PWD}/calculations/wkdir_${LABEL}
  echo Creating WORKDIR: $WORKDIR
  #exit 
  if [ -d $WORKDIR ]; then
    echo $WORKDIR already exists, skipping
    continue
  fi
  mkdir -p $WORKDIR
  #exit
  $qe_pos_setting_merge $SETTING_FILE $POSITION_FILE $WORKDIR/$INPUT_FILE || exit 1;
  rm -f $POSITION_FILE
  cp $SUBMIT_SCRIPT $WORKDIR/ || exit 1;

  cd $WORKDIR
  if [ $DEBUG_MODE == "1" ]; then
    sbatch -p debug -t 01:00:00 $SUBMIT_SCRIPT || exit 1;
  elif [ $PRODUCTION_MODE == "1" ]; then
    sbatch $SUBMIT_SCRIPT || exit 1;
  elif [ $DRYRUN_MODE == "1" ]; then
      ka=kb
    #echo "generating $WORKDIR"
  else
    echo "Script not run in a valid mode!"
    exit 1
  fi
  cd -

  done
