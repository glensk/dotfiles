#!/bin/sh
[ "$1" = "" ] && echo 'please provide $1 to be the inputxxx.data file to be calculated by DFT/aiida' && exit

### make README
readmedate=`date +"%Y-%m-%d_%H:%M:%S"`
filename="README_$readmedate.txt"
touch $filename
echo n2p2_submit_inputdata_to_aiida.sh $1 >> $filename
echo >> $filename

### scp
bn=`basename $1`
scp $1 cosmopc:/local/scratch/glensk/Dropbox/Albert/scripts/dotfiles/aiida-alloy
echo
echo "now ssh: loading structures as aiida group"
echo

ssh glensk@cosmopc15.epfl.ch "source /home/glensk/aiida/bin/activate;cd /local/scratch/glensk/Dropbox/Albert/scripts/dotfiles/aiida-alloy; ./load_runner_dataset.py -d $bn -gn $bn.group && verdi group list -A"

echo
echo "now ssh: submitting the job"
echo
ssh glensk@cosmopc15.epfl.ch "source /home/glensk/aiida/bin/activate;cd /local/scratch/glensk/Dropbox/Albert/scripts/dotfiles/aiida-alloy; ./launch_workflow_alalloy_scf.py \
    --code_node \"1\" \
    --structure_group_name \"$bn.group\" \
    --workchain_group_name \"$bn.group_calc\" \
    --base_parameter_node \"9b370584-3f56-471c-a724-dbaadf022ec5\" \
    --pseudo_familyname \"SSSP_v1.1_eff\" \
    --kptper_recipang 80 \
    --nume2bnd_ratio 0.75 \
    --max_wallclock_seconds 21600 \
    --max_active_calculations 300 \
    --sleep_interval 600"

echo retrieve: $bn.group_calc >> $filename
echo
echo "now ssh: checking job status"
echo
ssh glensk@cosmopc15.epfl.ch "source /home/glensk/aiida/bin/activate;cd /local/scratch/glensk/Dropbox/Albert/scripts/dotfiles/aiida-alloy; verdi calculation list -a -p1"

