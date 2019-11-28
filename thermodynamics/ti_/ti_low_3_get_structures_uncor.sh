#!/bin/sh


## this script creates structures_uncor_lambda$l; 
## it randomizes the structures, but
## it keeps the old known structures in same order and adds new if there are new structures
## this script assumes to be in the low run

## THIS script does not! renew the structures(_vasprun).gz files!!! it just takes the available ones! ### --> ###


## check input: 3.74 1360 0.0
## check input: 3.74 1360 0.0
[ "`echo $* | wc -w | sed 's|[ ]*||g'`" != "3" ] && echo 'input missing, e.g.: 3.74 1360 0.0 e.g: -a -a -a for all folder (or 3.74 -a 0.2)' && exit


## check and get pathes: this script assumes to be in the low run
## check and get pathes: this script assumes to be in the low run
[ "`pwd | grep high | wc -w | sed 's|[ ]*||g'`" != "0" ] && echo it seems you are in high folder && exit
low_path=`pwd`

aLats=$1
temps=$2
lambdas=$3



[ "$aLats" = "-a" ] && aLats=`ls -1d *Ang_*K | sed 's|Ang.*||' | sort | uniq`
[ "$temps" = "-a" ] && temps=`ls -1d *Ang_*K | sed 's|.*_||' | sed 's|K||g' | sed 's|/||g' | sort | uniq`
[ "$lambdas" = "-a" ] && lambdas=`ls -1d *Ang_*K/lambda*_* | sed 's|.*lambda\(.*\)_.*|\1|' | sort | uniq`


#aLats="3.62 3.66 3.7 3.74"
#temps="450 800 1100"

#echo "###########################################################"
#echo "###   aLats:   "$aLats " (-a has all possible)"
#echo "###   temps:   "$temps " (-a has all possible)"
#echo "###   lambdas: "$lambdas " (-a has all possible)"
#echo "###########################################################"



for a in $aLats; do for t in $temps; do 
  folder=$a\Ang_$t\K
  #echo "######################################################################"
  #echo "############### $folder"
  #echo "######################################################################"
  [ ! -e "$folder" ] && echo "folder $folder does not exist" && continue


for l in $lambdas; do ## schleife ueber a=aLats, t=temps, l=lambdas
  folder=$a\Ang_$t\K
  [ ! -e "$low_path/$folder" ] && echo "folder $folder does not exist" && continue
  ## falls es das lambda nicht gibt!
  [ "`find -L $low_path/$folder -maxdepth 1 -type d -name "lambda$l\_*"`" = "" ] && continue #echo $low_path/$folder does not have lambda$l && continue


  ## get all seeds
  ## get all seeds
  seeds_all=`ls -1d $low_path/$folder/lambda$l\_* | sed 's/.*lambda.*_\(.*\)[\/]*/\1/'`
  #echo seedsall: $seeds_all

  ## offset && corrlength
  ## offset && corrlength $folder = angKfolder
  [ ! -e "$low_path/$folder/avg_dUdL" ] && offset=31
  [ -e "$low_path/$folder/avg_dUdL" ] && offset=`awk '$1=='$l' {print $9+1}' $low_path/$folder/avg_dUdL` #macht aus 30 in avg_dUdL offset 31
  [ "`isnumber.sh $offset`" != "yes" ] && offset=31

  [ ! -e "$low_path/$folder/avg_dUdL" ] && corLen=15
  [ -e "$low_path/$folder/avg_dUdL" ] && corLen=`awk '/correlation length/{print $4}'  $low_path/$folder/avg_dUdL`
  [ "`isnumber.sh $corLen`" != "yes" ] && corLen=15
  #echo "   offset:$offset corLen:$corLen"

  ## create structures_uncor_lambda$l (+randomize)
  ## create structures_uncor_lambda$l (+randomize)
  for seed in $seeds_all;do
     anzahl_total=OUTCARpositions_not_available
     anzahl_uncor=OUTCARpositions_not_available
     if [ -e "`echo $low_path/$folder/lambda$l\_$seed/workDir/ion_energies`" ];then
     seed_structure=`tail -n+$offset $low_path/$folder/lambda$l\_$seed/workDir/ion_energies | awk 'NR%'"$corLen"'==1' | awk '{print $1}' | sed 's|\(.*\)|'"$seed"'_\1|'`
     anzahl_total=`cat $low_path/$folder/lambda$l\_$seed/workDir/ion_energies | wc -l | sed 's|[ ]*||g'`
     anzahl_uncor=`echo $seed_structure | wc -w | sed 's|[ ]*||g'`
     # echo "structures: $anzahl_total  || UNCOR: $anzahl_uncor"
     # echo 3: offset: $offset $low_path/$folder/lambda$l\_$seed/ion_energies
     seed_structure_all="$seed_structure_all $seed_structure"
 
 elif [ -e "`echo $low_path/$folder/lambda$l\_$seed/ion_energies`" ];then
    seed_structure=`tail -n+$offset $low_path/$folder/lambda$l\_$seed/ion_energies | awk 'NR%'"$corLen"'==1' | awk '{print $1}' | sed 's|\(.*\)|'"$seed"'_\1|'`
     anzahl_total=`cat $low_path/$folder/lambda$l\_$seed/ion_energies | wc -l | sed 's|[ ]*||g'`
     anzahl_uncor=`echo $seed_structure | wc -w | sed 's|[ ]*||g'`
     # echo "structures: $anzahl_total  || UNCOR: $anzahl_uncor"
     # echo 3: offset: $offset $low_path/$folder/lambda$l\_$seed/ion_energies
     seed_structure_all="$seed_structure_all $seed_structure"
     fi;
     echo    $low_path/$folder/lambda$lambda$l\_$seed  "||" "structures: $anzahl_total  || UNCOR: $anzahl_uncor"
  done
 
  
  
  ## get old structures_uncor_lambda$l
  ## get old structures_uncor_lambda$l
  rm -f $low_path/$folder/structures_uncor_lambda$l\_all
  #echo "$seed_structure_all"
  #exit
  #echo "$seed_structure_all" | xargs -n1 | random.pearl.sh > $low_path/$folder/structures_uncor_lambda$l\_all
  echo "$seed_structure_all" | xargs -n1 | sort -R > $low_path/$folder/structures_uncor_lambda$l\_all
  seed_structure_known=$low_path/$folder/structures_uncor_lambda$l
  seed_structure_all=$low_path/$folder/structures_uncor_lambda$l\_all

  [ ! -e "$seed_structure_known" ] && touch $seed_structure_known
  rm -f tmpx; cat $seed_structure_known $seed_structure_all > tmpx; rm -f $seed_structure_known $seed_structure_all
  perl -ne '$H{$_}++ or print' tmpx > $seed_structure_known; rm tmpx; sed -i 's|/.*||g' $seed_structure_known;
  [ "`cat $seed_structure_known | wc -w | sed 's|[ ]*||g'`" -le "1" ] && rm $seed_structure_known
done; done; done; ## schleife ueber a=aLats, t=temps, l=lambdas
