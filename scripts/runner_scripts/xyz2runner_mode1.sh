#!/bin/sh

# for instance from fidis: /scratch/glensk/raw_data/test_3_at_1000K_for_fps
# for instance on cosmopc: /local/scratch/glensk/runner_scratch
#
# to get it on mac:
# sudo mkdir -p /local/scratch/glensk
# sudo chown -R glensk /local
# sshfs glensk@cosmopc15.epfl.ch:/local/scratch/glensk/ /local/scratch/glensk -o reconnect -C

<<<<<<< HEAD
=======

>>>>>>> Thu Nov  8 09:48:52 CET 2018
# filename="input.data.all"
# rm -f $filename  
# touch $filename 
# files=`find . -maxdepth 4 -name simulation.pos_0.xyz`
# for i in $files;do
#     xyz2runner.sh $i >> $filename
# done

rm -f symfun.output
<<<<<<< HEAD
symfun_gen.py -e Al,Si,Mg -c 16 -n 30 > symfun.output

cp $scripts/runner_scripts/input.nn .
=======
symfun_gen.py -e Al,Si,Mg -c 12 -n 10 > symfun.output  # order is not important
symfun_gen.py -e Al,Si,Mg -c 16 -n 10 >> symfun.output  # order is not important
symfun_gen.py -e Al,Si,Mg -c 20 -n 10 >> symfun.output  # order is not important
symfun_gen.py -e Al,Si,Mg -c 16 -n 4 >> symfun.output  # order is not important
sort symfun.output | uniq > tmp
mv tmp symfun.output
rm tmp


###################################################
# get and adapt input.nn
###################################################
cp $scripts/runner_scripts/input.nn .
linebegin=`grep -n "^symfunction_short" input.nn | head -1 | sed 's|:.*||'`
lineend=`grep -n "^symfunction_short" input.nn | tail -1 | sed 's|:.*||'`
echo lb $linebegin
echo le $lineend
sed ''"$linebegin"','"$lineend"'d' input.nn > tmp
sed -i '/cutoff_type 2/r symfun.output' input.nn


frame_selector.py --prefix ni input.data.all random 113 
>>>>>>> Thu Nov  8 09:48:52 CET 2018
