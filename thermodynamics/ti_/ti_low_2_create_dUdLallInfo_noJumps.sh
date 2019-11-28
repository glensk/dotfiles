#!/bin/sh

out=no
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'` #[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'` #[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -q -d" #[ "$out" = "yes" ] && echo options: $options


if [ "`getOption -h`" = "True" ] || [ "`getOption -help`" = "True" ]; then  # here is is ok if we dont get an argument
  usage `basename $0`
  printOptions " " \
               "-q          quiet" \
               "-d          debug mode for plotting inbetween information"
  exit
fi


###### check where you are (how deep in folder structure)
gotofolder=`pwd`
checkfolder=`pwd`
gotolam_folder=`find -L $checkfolder -maxdepth 2 -mindepth 1 -type d -name lambda'[0-9.]*_[0-9]*'`  ## keep maxdepth 2 in case you are in *Ang_*K folder
#echo lam:$gotolam_folder
gotolam_anz=`echo $gotolam_folder | wc -w | sed 's|[ ]*||g'`
[ ! -e "OUTCAR" ] && [ ! -e "OUTCAR.gz" ] && [ ! -e "vasprun.xml" ] && [ ! -e "vasprun.xml.gz" ] && [ -d "workDir" ] && gotofolder=`pwd`/workDir
[ ! -e "OUTCAR" ] && [ ! -e "OUTCAR.gz" ] && [ ! -e "vasprun.xml" ] && [ ! -e "vasprun.xml.gz" ] && [ "$gotolam_anz" != "0" ] && gotofolder=$gotolam_folder
folder=$gotofolder
angkfolders=`echo "$gotofolder" | sed 's|\(.*\)/lambda[0-9.]*_[0-9]*.*|\1|' | sort | uniq`
for angkfolder in $angkfolders;do  ## angkschleife
    cd $angkfolder
#for i in $gotofolder;do
#    echo $i
#done

qdel=`which qdel.sh`
[ "$folder" = "" ] && echo no lambda folder found && exit
hier=`pwd`


## die files jumps & jumps_avg waren bisher immer im *Ang_*K folder
rm -f jumps jumps_avg
[ "$1" != "-q" ] && echo; 
[ "$1" != "-q" ] && echo "#lambda   ||   jDist  jRadius ||   deleted      left" | awk '{printf "%-20s %-3s %-7s %-7s %-3s %-10s %-9s\n",$1,$2,$3,$4,$5,$6,$7}'
[ "$1" != "-q" ] && echo "# -----------------------------------------------------------------"
echo "# lambda  jDist jRadius" > jumps


for i in $folder;do
cd $i
[ "`getOption -d`" = "True" ] && echo ii:$i
[ ! -e dUdLallInfo ] && [ -d "workDir" ] && cd workDir

    ## check if necessary files exist
    ## check if necessary files exist
    lambda=`echo $i | sed 's|lambda\([0-9.]*\).*|\1|'`
    [ ! -e dUdLallInfo ] && [ "$1" != "-q" ] && echo "$i   || dUdLallInfo does not exist !!!!!!!!!!!!!!" | awk '{printf "%-20s %-6s %-6s %-6s %-6s %-6s %-6s\n", $1,$2,$3,$4,$5,$6,$7}' 
    [ ! -e dUdLallInfo ] && cd $hier && continue
    [ "`getOption -d`" = "True" ] && echo hh: `pwd`

    ## if dUdLallInfo_noJumps exist --> check if already done
    ## if dUdLallInfo_noJumps exist --> check if already done
    jumptonext="no"
    dUdLmax=`tail -1 dUdLallInfo | awk '{print $1}'`  # 1758
    #if [ -e "dUdLallInfo_noJumps" ];then
        #if [ -e "dUdLallInfo" ]; then


    dUdLmaxdone="0"
    [ -e dUdLallInfo_noJumps ] && dUdLmaxdone=`head -1 dUdLallInfo_noJumps | sed 's|.*###||'`
    [ "`getOption -d`" = "True" ] && echo dUdLmax: $dUdLmax  dUdLmaxdone:$dUdLmaxdone
    [ -e "dUdLallInfo_noJumps" ] && [ "$dUdLmax" = "$dUdLmaxdone" ] && [ "$1" != "-q" ] && echo "$i     || latest dUdLallInfo_max:$dUdLmax: dUdLallInfo_noJumps_maxdone:$dUdLmaxdone:" 
    [ -e "dUdLallInfo_noJumps" ] && [ "$dUdLmax" = "$dUdLmaxdone" ] && cd $hier && continue

    ## if already here and dUdLallInfo_noJumps exist it can be removed
    ## if already here and dUdLallInfo_noJumps exist it can be removed
    [ -e "dUdLallInfo_noJumps" ] && rm dUdLallInfo_noJumps

    ## part getJumps
    ## part getJumps
    #echo 1
##################################################################################################################################################################
    [ "`getOption -d`" = "True" ] && echo vor extractPOSITIONS.sh  pwd: `pwd`                                              
    echo extractPOSITIONS.sh in `pwd`
    extractPOSITIONS.sh        ## erstellt POSITIONs, cell und atoms_volume_steps aus der OUTCAR file                      
    [ "`getOption -d`" = "True" ] && echo nach extractPOSITIONS.sh                                                         
    #echo 2
    [ ! -e "cell" ] && echo "$i   || PROBLEM: getJump needs cell which does not exist! -> SURE THERE ARE NO JUMPS? i remove dUdLallInfo to be sure" && rm -f dUdLallInfo && cd $hier && continue
    [ ! -e "POSITIONs" ] && echo "$i   || PROBLEM: getJump needs POSITIONs which do not exist!  -> SURE THERE ARE NO JUMPS? i remove dUdLallInfo to be sure" && rm -f dUdLallInfo && cd $hier && continue
    [ ! -e "atoms_volume_steps" ] && echo "$i   || PROBLEM: getJump needs atoms_volume_steps which do not exist!  -> SURE THERE ARE NO JUMPS? i remove dUdLallInfo to be sure" && rm -f dUdLallInfo && cd $hier && continue               #hier!
    #echo 22222222222222222222222
    getJump.x
    # we can create the POSITIONS all the time, otherwise those use a lot of space!
    rm -f POSITIONs
    #echo 3
    [ ! -e "jumps_dist" ] || [ ! -e "jumps_radius" ] && echo "$i   || getJump_problem_!!!!!!!!!!!!" | awk '{printf "%-6s %-6s %-6s %-6s \n", $1,$2,$3,$4}'  && cd $hier && continue

    [ "`getOption -d`" = "True" ] && echo jjumps_dist jumps_radius 
    [ ! -e "jumps_dist" ] && echo prob jumps_dist missing
    [ ! -e "jumps_radius" ] && echo prob jumps_radius missing

    count=`cat jumps_dist | wc -l`                                                                 
    
    # just counts the number of jumps, not the total number of atoms which are involved in a jump 
    j=`cat jumps_dist | sed 's|0||' | sed 's|[ ]*||' | sed -e '/^\s*$/d' | wc -l | sed 's|[ ]*||'`
    #j=`awk 'BEGIN{s=0};{s=s+$1};END{print s}' jumps_dist`
    # just counts the number of jumps, not the total number of atoms which are involved in a jump 
    jr=`cat jumps_radius | sed 's|0||' | sed 's|[ ]*||' | sed -e '/^\s*$/d' | wc -l | sed 's|[ ]*||'` 
    #jr=`awk 'BEGIN{s=0};{s=s+$1};END{print s}' jumps_radius`
    if [ "$count" -gt 19 ]                                                                         
       then checkjumpexists_dist=`tail -20 jumps_dist | sort | uniq | sed 's|[ ]*||' | grep 0`          
       else checkjumpexists_dist=0                                                                      
    fi                   
                                                                              
    if [ "$count" -gt 19 ]                                                                         
       then checkjumpexists_radius=`tail -20 jumps_radius | sort | uniq | sed 's|[ ]*||' | grep 0`          
       else checkjumpexists_radius=0                                                                      
    fi
    ## in case we have a jump
    #echo 44
    #[ "$j" != "0" ] || [ "$jr" != "0" ] && ID=`OUTCAR_ID.sh` && IDcheck=`isnumber.sh $ID`
    [ "$checkjumpexists_dist" != "0" ] && [ "$checkjumpexists_radius" != "0" ]  && ID=`OUTCAR_ID.sh` && IDcheck=`isnumber.sh $ID`
    #echo ID:$ID kk:$IDcheck
    
    #which qdel

    [ "`getOption -d`" = "True" ] && echo IDcheck: $IDcheck  ID: $ID checkjumpexists_dist:$checkjumpexists_dist: checkjumpexists_radius:$checkjumpexists_radius:
    [ ! -e "$qdel" ] && echo qdel does not exist and will not WORK!!!! jobs stopped if there is a jump!
    # if there is one correct structure without a jump out of the last 20 then dont kill the job!
    [ -e "$qdel" ] && [ "$checkjumpexists_dist" != "0" ] && [ "$checkjumpexists_radius" != "0" ]  && [ "$IDcheck" = "yes" ] &&  outqdel=`$qdel $ID` && echo -e "\n\033[1;4;31m  $outqdel\033[0m"
    # not sure what jumps_dist is
    #[ -e "$qdel" ] && [ "$j" != "0" ] || [ "$jr" != "0" ] && [ "$IDcheck" = "yes" ] && $qdel $ID 
    #[[ "$j" != "0" || "$jr" != "0" ]] && [ "$IDcheck" != "yes" ] && echo "wanted to qdel, but didnt get the jobID. ID: $ID ( written to /home/$USER/JOBLIST_COULD_NOT_DELETE )"  && touch /home/$USER/JOBLIST_COULD_NOT_DELETE && echo `pwd`>> /home/$USER//JOBLIST_COULD_NOT_DELETE 
    
    ## in case we have no jump
    #jumps="$j $jr" ## $j=0 && $jr=0 --> keine jumps
    #jumps="$j $jr" ## $j=0 && $jr=0 --> keine jumps
    #echo 33
    echo $lambda $jumps >> $hier/jumps
    [ "$j" = "0" ] && [ "$jr" = "0" ] && cp dUdLallInfo dUdLallInfo_noJumps && sed -i '1 s|\(.*\)|\1  ###'"$dUdLmax"'|' dUdLallInfo_noJumps
    [ "$j" = "0" ] && [ "$jr" = "0" ] && [ "$1" != "-q" ] && echo "$i   || $j $jr" | awk '{printf "%-20s %-6s %-6s %-6s \n", $1,$2,$3,$4}' 
    [ "$j" = "0" ] && [ "$jr" = "0" ] && cd $hier && continue
    #echo 4


## part removeJumps
## part removeJumps
  rm -f jumps dUdLallInfo_noJumps; echo $j $jr > jumps
  stepmax=`cat dUdLallInfo | tail -n+2 | wc -l | sed 's|[ ]*||g'`

  [ "`getOption -d`" = "True" ] && echo stepmax:$stepmax
  rm -f tmp; echo jj 0 > tmp; awk '{print "jj",$1}' jumps_dist | head -$stepmax >> tmp
  rm -f tmp2; echo jjr 0 > tmp2; awk '{print "jjr",$1}' jumps_radius | head -$stepmax >> tmp2

  #paste dUdLallInfo tmp tmp2 | awk 'NF>1&&$NF==0||$(NF-2)==0{print $0}' | grep "jj 0\|jjr 0" | sed 's|\(.*\)jj.*|\1|' > dUdLallInfo_noJumps  && sed -i '1 s|\(.*\)|\1  ###'"$dUdLmax"'|' dUdLallInfo_noJumps  ### past dUdLallinof tmp haengt an das dUdLallInfo immer jj 0, jj 1, jj2 an jede zeile dran; (NF ist die anzahl der spalte, $NF ist die letzte spalte)
  paste dUdLallInfo tmp tmp2 | awk 'NF>1&&$NF==0||$(NF-2)==0{print $0}' | grep "jj 0\|jjr 0"  > dUdLallInfo_noJumps  && sed -i '1 s|\(.*\)|\1  ###'"$dUdLmax"'|' dUdLallInfo_noJumps  ### past dUdLallinof tmp haengt an das dUdLallInfo immer jj 0, jj 1, jj2 an jede zeile dran; (NF ist die anzahl der spalte, $NF ist die letzte spalte)
  #!!paste dUdLallInfo tmp | awk 'NF>1&&$NF==0{print $0}' | grep "jj 0"  > dUdLallInfo_noJumps  && sed -i '1 s|\(.*\)|\1  ###'"$dUdLmax"'|' dUdLallInfo_noJumps  ### past dUdLallinof tmp haengt an das dUdLallInfo immer jj 0, jj 1, jj2 an jede zeile dran; (NF ist die anzahl der spalte, $NF ist die letzte spalte)
  #rm -f tmp tmp2

  [ "`getOption -d`" = "True" ] && echo vor n
  n=`wc dUdLallInfo dUdLallInfo_noJumps | awk '$NF=="dUdLallInfo"{l1=$1};$NF=="dUdLallInfo_noJumps"{l2=$1};END{print l2-l1}'`   ## n=0 -> noJump, n=1 einJump, n=78 -> 78Jumps, ...
  
  [ "`getOption -d`" = "True" ] && echo n: $n
  left=`wc -l dUdLallInfo_noJumps | awk '{print $1-1}'`
  # wc |
  #if [ "$n" == 0 ]; then
  #    rm -f jumps_dist tmp dUdLallInfo_noJumps
  #else
  #    #awk 'NR==1{print};NR>1{printf("%7d %9.1f %8.1f %8.1f  %11.2f %9.2f  %11.2f %9.2f %9.2f\n",$1,$2,$3,0,$5,$6,$7,0,0)}' dUdLallInfo_noJumps > tmp; rm -f dUdLallInfo_noJumps; mv tmp dUdLallInfo_noJumps
  #    awk 'NR==1{print};NR>1{print $0}' dUdLallInfo_noJumps > tmp; rm -f dUdLallInfo_noJumps; mv tmp dUdLallInfo_noJumps
  #fi
  echo $lambda $n >> $hier/removed
  [ "$1" != "-q" ] && echo "$i   || $j $jr || $n $left" | awk '{printf "%-20s %-3s %7.0f %7.0f %-3s %7.0f %7.0f\n",$1,$2,$3,$4,$5,$6,$7}'


  [ "`getOption -d`" = "True" ] && echo test

#jumps=`cat jumps`
#echo "i: $i lambda:$lambda || $lambda $jumps"
#echo $lambda $jumps >> $hier/jumps



cd $hier
done

awk 'BEGIN{l=-1};l==$1{s=s+$2;s2=s2+$3};l!=$1&&l!=-1{print l,s,s2};l!=$1{l=$1;s=$2;s2=$3};END{print l,s,s2}' jumps > jumps_avg

done ### angkschleife
echo ":) `basename $0` finished, last line"
