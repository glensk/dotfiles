#!/bin/sh
[ "$1" = "" ] && echo 'provide $1 to be the number of random seeds' && exit

[ ! -e "input.nn" ] && echo "input.nn missing" && exit
[ ! -e "input.data" ] && echo "input.data missing" && exit
[ ! -e "scaling.data" ] && echo "scaling.data missing" && exit
[ ! -e "submit_training_debug.sh" ] && echo "submit_training_debug.sh missing" && exit

for i in `seq $1`;do
    rand=$RANDOM
    folder="random_seed_$rand"
    echo $i $rand $folder
    if [ ! -e "$folder" ];then
        mkdir $folder
        cd $folder
        cp ../input.nn .
        cp ../input.data .
        cp ../scaling.data .
        cp ../submit_training_debug.sh .
        sed -i 's|^random_seed.*|random_seed '"$rand"'|' input.nn
        sbatch submit_training_debug.sh
        cd ..
    fi
done
