#!/bin/sh

# for instance in /scratch/glensk/raw_data/test_3_at_1000K_for_fps

# filename="input.data.all"
# rm -f $filename  
# touch $filename 
# files=`find . -maxdepth 4 -name simulation.pos_0.xyz`
# for i in $files;do
#     xyz2runner.sh $i >> $filename
# done

rm -f symfun.output
symfun_gen.py -e Al,Si,Mg -c 16 -n 30 > symfun.output

cp $scripts/runner_scripts/input.nn .
