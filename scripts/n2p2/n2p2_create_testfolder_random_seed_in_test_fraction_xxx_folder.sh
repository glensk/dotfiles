#!/bin/sh

seeds="1 10 3049 43223 102376 1234123 1234124 2234125 5234123 8234123 9234123 32341233"


test_fraction=`pwd | sed 's|.*test_fraction_||'`
hier=`pwd`

### get actual path of subfolder
subfolder="../test_fraction_vorlage"
cd $subfolder
subfolder=`pwd`
cd $hier

for i in $seeds;do
    cd $hier
    folder="random_seed_$i"
    [ -e "$folder" ] && echo $folder exists && continue
    [ ! -e "$subfolder/input.nn" ] && echo "../input.nn does not exist;Exit" && exit
    echo $folder
    mkdir $folder
    cd $folder
    [ ! -e "$subfolder/input.nn" ] && echo "$subfolder/input.nn missing!" && exit
    [ ! -e "$subfolder/input.data" ] && echo "$subfolder/input.data missingl!" && exit
    [ ! -e "$subfolder/scaling.data" ] && echo "$subfolder/scaling.data missing!" && exit
    [ ! -e "$subfolder/submit_training_debug.sh" ] && echo "$subfolder/submit_training_debug.sh missing!" && exit
    cp $subfolder/input.nn .
    cp $subfolder/input.data .
    cp $subfolder/scaling.data .
    cp $subfolder/submit_training_debug.sh .
    sed -i 's|test_fraction.*|test_fraction '"$test_fraction"'|' input.nn
    sed -i 's|random_seed.*|random_seed '"$i"'|' input.nn
    sbatch submit_training_debug.sh
done
