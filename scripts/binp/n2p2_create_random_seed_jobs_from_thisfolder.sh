#!/bin/sh
[ "$1" = "" ] && echo 'provide $1 to be the number of random seeds' && exit
hier=`pwd`
[ ! -e "input.nn" ] && echo "input.nn missing" && exit
[ ! -e "input.data" ] && echo "input.data missing" && exit
[ ! -e "scaling.data" ] && echo "scaling.data missing" && exit
trainfile="submit_training_debug.sh"
trainfile="submit_auswertung_debug.sh"
[ ! -e "$trainfile" ] && echo "$trainfile missing" && exit

for i in `seq $1`;do
    rand=$RANDOM
    folder="random_seed_$rand"
    echo $i $rand $folder
    [ ! -e "input.nn" ] && echo input.nn missing && exit
    [ ! -e "input.data" ] && echo input.data missing && exit
    [ ! -e "scaling.data" ] && echo scaling.data missing && exit
    [ ! -e "$trainfile" ] && echo $trainfile data missing && exit
    if [ ! -e "$folder" ];then
        mkdir $folder
        cd $folder
        cp ../input.nn .
        cp ../input.data .
        cp ../scaling.data .
        cp ../$trainfile .
        sed -i 's|^random_seed.*|random_seed '"$rand"'|' input.nn
        sbatch $trainfile 
        cd ..
    fi
done

# make a README
cd $hier
ts=`date +"%Y-%m-%d_%H:%M:%S"`
fn=README_$ts.txt
touch $fn
echo `pwd` >> $fn
echo "`basename $0` $*" >> $fn
