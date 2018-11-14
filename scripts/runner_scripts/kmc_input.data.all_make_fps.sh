#!/bin/sh

# for instance from fidis: /scratch/glensk/raw_data/test_3_at_1000K_for_fps
# for instance on cosmopc: /local/scratch/glensk/runner_scratch
# for instance on cosmopc: /home/glensk/runner_scratch

# to get it on mac:
# sudo mkdir -p /local/scratch/glensk
# sudo chown -R glensk /local
# sshfs glensk@cosmopc15.epfl.ch:/local/scratch/glensk/ /local/scratch/glensk -o reconnect -C

# sudo mkdir -p /home/glensk/
# sudo chown -R glensk /home/

echo "here, a README.txt is missing..."
echo "scripts: $scripts"

###################################################
echo "make sure input.data.all is available"
###################################################
[ ! -e "input.data.all" ] && echo 'run fist kmc_gather_all_xyz.py on the kmc folder(s)' && exit

###################################################
echo "define executable"
###################################################
executable="$HOME/Dropbox/Albert/scripts/runner_source/RuNNer.serial.cosmopc.natascha.x"
[ ! -e "$executable" ] && echo "executable $executable does not exist; Exit" && exit

###################################################
echo "get symmetry functions"
###################################################
rm -f symfun.output
symfun_gen.py -e Al,Si,Mg -c 12 -n 10 > symfun.output  # order is not important
symfun_gen.py -e Al,Si,Mg -c 16 -n 10 >> symfun.output  # order is not important
symfun_gen.py -e Al,Si,Mg -c 20 -n 10 >> symfun.output  # order is not important
symfun_gen.py -e Al,Si,Mg -c 16 -n 4 >> symfun.output  # order is not important
sort symfun.output | uniq > tmp
mv tmp symfun.output


###################################################
echo "get and adapt input.nn"
###################################################
cp $scripts/runner_scripts/inputAlMgSi.nn input.nn
linebegin=`grep -n "^symfunction_short" input.nn | head -1 | sed 's|:.*||'`
lineend=`grep -n "^symfunction_short" input.nn | tail -1 | sed 's|:.*||'`
echo "     ... linebegin $linebegin"
echo "     ... lineend $lineend"
sed -i ''"$linebegin"','"$lineend"'d' input.nn 
sed -i '/cutoff_type 2/r symfun.output' input.nn
sed -i 's|^test_fraction .*|test_fraction 0|g' input.nn
sed -i 's|^runner_mode .*|runner_mode 1|g' input.nn
sed -i 's|^number_of_elements .*|number_of_elements 3|' input.nn
sed -i 's|^elements .*|elements Al Mg Si|' input.nn


####################################################
#echo "frame selector   --> is probably not what we need since we want to do FPS, or?"
####################################################
#frame_selector.py --prefix AlMgSi input.data.all random 113
#mv AlMgSi_selected.data input.data

###################################################
echo "start runner in mode 1 (on cosmopc)"
echo "this took 16h for 2435 structures 18:21 - 10:46 of next day (teststruct.data and debug.out)"
echo "this took  1h for  113 structures 18:10 - 18:56 of same day (teststruct.data and debug.out)"
###################################################
$executable > logfile_mode1.1&

###################################################
echo "make fps to get structures to calculate"
echo "cursel.distances.dat shows the distances"
###################################################
#CurSel.py -t 1e-3 --landmarks 40 ../runner_scratch_new_data_all/function.data ../runner_scratch_new_data_all/logfile_mode1.1
# tool 490 seconds on 40 structures from function.data of 3.2GB
#frame_selector.py input.data.all precomp cursel.landmarks  # crates input.dat_selected.data with 40 strut

CurSel.py -t 1e-3 --landmarks `grep -c begin input.data` function.data logfile_mode1.1
head -20 cursel.landmarks > cursel.landmarks.20
head -40 cursel.landmarks > cursel.landmarks.40

frame_selector.py input.data precomp cursel.landmarks.20
mv input_selected.data input_selected.data.20

frame_selector.py input.data precomp cursel.landmarks.40
mv input_selected.data input_selected.data.40

# took 22 secons on 113 structures 
# took 525 secons on 2455 structures
# this will create a file cursel.landmarks which contains indexes;
# x cursel.distances.dat (out.dat) will then give the distances
