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
echo "get environment scripts variable from click. make defulat ones"
exit 

###################################################
echo "define source folder (on fidis)"
echo "define target folder (on cosmopc, where Runner Mode 1 calculation will be done)"
echo "define executable"
###################################################
source_folder="/scratch/glensk/raw_data/test_3_at_1000K_for_fps"
target_folder="/local/scratch/glensk/runner_scratch"
executable="$HOME/Dropbox/Albert/scripts/runner_source/RuNNer.serial.cosmopc.natascha.x"

###################################################
echo "check if source folder exist $source_folder"
###################################################
[ ! -e "$source_folder" ] && echo "source_folder $source_folder does not exist (mount it); Exit" && exit

###################################################
echo "check if target folder exist $target_folder"
###################################################
[ ! -e "$target_folder" ] && echo "target_folder $target_folder does not exist (mount it); Exit" && exit

###################################################
echo "check if executable exist $executable"
###################################################
[ ! -e "$executable" ] && echo "executable $executable does not exist; Exit" && exit
cp $executable $target_folder

###################################################
echo "cd to source_folder $source_folder"
###################################################
cd $source_folder 

###################################################
echo "get input.data.all"
###################################################
filename="input.data.all"
rm -f $filename  
touch $filename 
files=`find . -maxdepth 4 -name simulation.pos_0.xyz`
for i in $files;do
    echo $i
    xyz2runner.sh $i >> $filename
done

###################################################
echo "cp $filename to target_folder $folfer_folder"
echo "cd target_folder $target_folder"
###################################################
cp $filename $target_folder
cd $target_folder

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

###################################################
echo "frame selector   --> is probably not what we need since we want to do FPS, or?"
###################################################
frame_selector.py --prefix AlMgSi input.data.all random 113
mv AlMgSi_selected.data input.data

###################################################
echo "start runner in mode 1 (on cosmopc)"
###################################################


###################################################
echo "make fps to get structures to calculate"
###################################################
CurSel.py -t 1e-3 --landmarks `grep -c begin input.data` function.data logfile_mode1.1
# took 22 secons on 113 structures 
# took 525 secons on 2455 structures
# this will create a file cursel.landmarks which contains indexes;
# x cursel.distances.dat (out.dat) will then give the distances
cp cursel.landmarks cursel.landmarks.all
# and in cursel.landmarks delete everythin after line 100 to get only first 100 structures.
frame_selector.py input.data precomp cursel.landmarks

