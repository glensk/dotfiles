#!/bin/sh

out=no
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -a -c -v -s -A -ga -ol -oh";[ "$out" = "yes" ] && echo options: $options
path_ti=$path/


printhelp() {
    usage "e.g. `basename $0` 4.25Ang_1100K 4.25Ang_700K            (to evaluate named folder)"
    usage "e.g. `basename $0` low_2x2x2sc_230eV_2x2x2kp    low_2x2x2sc_230eV_3x3x3kp_NGXF120"
    usage "e.g. `basename $0` .                                     (to evaluate this current folder)"
    usage "e.g. `basename $0`                                       (to evaluate all necessary high/low/ang*K folder seen from pwd)"
    usage "e.g. `basename $0` -a -c 7200"
    printOptions " " \
        "-a         run this skript for all currently running ti_jobs (qstatn -q -r -n | awk '{print \$NF}' | grep \"low_\" | grep -v \"high\" | sed 's|lambda.*||' | sort | uniq)" \
        "-c  [sec]  run in continues mode (suggested to run in screen) every [sec] seconds (default is 7200)" \
        "-s  [meV]  set convergence criterion; jobs wiss Encorrelated error less than convergence criterion [meV] will be killed (default is 0.15 meV)" \
        "-ga [INT]  go to (low) Ang_K folder with number [INT], skip all angK folder before (this is just done once before sleep in -a mode)" \
        "-ol        do only low jobs" \
        "-oh        do only high jobs" \
        "-v         be more verbose" \
        " "
        exit

}

[ `getOption -h` = True ] && printhelp
[ "$1" = "" ] && printhelp
[ "`getOption -a`" = "True" ] && [ -e "/data/$USER" ] && cd /data/$USER
if [ "`getOption -ga`" = "True" ];then
    angKfolderstart=`getValue -ga`
    [ "`isnumber.sh $angKfolderstart`" != "yes" ] && error "please specify an INTEGER after -ga" && exit
    firsttime="yes"
    echo angKfolderstart:$angKfolderstart:
fi

########################################################################################
# check submithost
########################################################################################
hier=`pwd`
[ "`pwd`" = "$HOME" ] && echo please dont start this skript from $HOME since problems with snapshot system && exit
#echo "-----------------------------------------------------"
checkAndSetMath

# if [ "`hostname`" != "cmmc001" ];then
# if [ "`hostname`" != "cmmc002" ];then

#   if [ ! -z "$submithost" ];then
#       # submithost is defined
#       #echo submithost: $submithost
#       [ "`hostname`" = "$submithost" ] && echo please dont run this from submitshost since it does not have mathemtica && exit
#       if [ "`hostname`" != "$submithost" ];then
#           echogreen " -> this is a check if your submithost (which is $submithost) is in knownhosts and if login free ssh is possible to this host."
#           echogreen " -> if you are asked for a password: this skript will not be able to kill already converged ti jobs"
#           ssh $submithost echo 
#           fi
#       else 
#       # submithost is not defined
#       echored "submithost is not defined -> converged jobs will not be killed"
#       fi
# fi
# fi
#echo "-----------------------------------------------------"
########################################################################################
# endlosschleife
########################################################################################
angKfolder_input=""
while [ "endlos" = "endlos" ];do    ## endlosschleife
## sleep timer
## sleep timer
schlafe=0
[ `getOption -c` = True ] && schlafe=10000
[ `getOption -c` = True ] && schlafe=`getValue -c`; [ "$schlafe" = "" ] && schlafe=10000

## continuesmode
## continuesmode
continuesmode=no;[ `getOption -c` = True ] && continuesmode=yes
[ `getOption -c` = True ] && echogreen "SLEEP in between  : $schlafe [seconds] (time between consecutive runs of `basename $0`)"

## convcrit
## convcrit
#convcrit=`tcsh -c 'echo $convcrit'` # This does not work and writes a wird string in the file
convcrit=$convcrit

[ `getOption -s` = True ] && convcrit=`getValue -s`

if [ "`isnumber.sh $convcrit`" != "yes" ];then
    convcrit=0.15
    echored "variable convcrit not set; try in your .tcshrc set convcrit=0.3 or so; manually set convcrit to $convcrit"
    fi

# convcrit=0.35 # we dont want to fix it
echogreen "CONVCRIT          : $convcrit [meV] errUnc/2 (convergence criterium after which jobs are killed";echo
IDstokillall=""
## check if submithost is in key-ring, if not print warning at beginning and  exit, but make an option -i to ignore this warning (the consequence is that you will be asked for your password is a run is converged and this skript wants to kill the job)


do=`pwd`
dostart=`pwd`
for iiii in $do;do   ## das geht einfach nur in das pwd... dass koennte man sich auch sparen
cd $iiii
#/bin/echo -e "\033[31m\033[1m`pwd`\033[0m"

##################################################################################################################################################################
## get/check necessary scripts && settings
##################################################################################################################################################################
## low scripts
dUdLallInfo=ti_low_1_create_dUdLallInfo.sh
removeJumps=ti_low_2_create_dUdLallInfo_noJumps.sh  
getStrucUncor=ti_low_3_get_structures_uncor.sh
getdUdLs=ti_low_4_get_dUdLs.sh
fit=ti_low_5_fit_avg_dUdL.sh                                                ## ---->>>> this one needs mathematica
auswertung=ti_low_6_auswertung.sh
avg_dUdL_converged=`which avg_dUdL_converged.py`

## high scripts
ion_energies=ti_high_1_get_ion_energies.sh                           
fit_up=ti_high_2_highFitUp.sh                                              ## ---->>>> this one needs mathematica
auswertung=ti_low_6_auswertung.sh


qdel=0; [ -e "`which qdel.sh`" ] && qdel=`which qdel.sh`

##################################################################################################################################################################
## check if low/high scripts available 
##################################################################################################################################################################
## low
[ ! -e `which $dUdLallInfo` ] && echo $dUdLallInfo does not exist && exit -1
[ ! -e `which $removeJumps` ] && echo $removeJumps does not exist && exit -1
[ ! -e `which $getStrucUncor` ] && echo $getStrucUncor does not exist && exit -1
[ ! -e `which $getdUdLs` ] && echo $getdUdLs does not exist && exit -1
[ ! -e `which $fit` ] && echo $fit does not exist && exit -1
[ ! -e `which $auswertung` ] && echo $auswertung does not exist && exit -1        ## removes old auswertung by itself
[ ! -e `which $avg_dUdL_converged` ] && echo $avg_dUdL_converged does not exist && exit -1        ## removes old auswertung by itself

#echo jo
#cd /nas/glensk/v/pp/au/ti_bulk_fcc4/low_2x2x2sc_300eV_03x03x03kp_EDIFF1E-2_TAKE/4.17Ang_1338K
#pwd
#cat avg_dUdL_fre
#echo ----------------
#echo 
#echo $avg_dUdL_converged
#cc=`echo $convcrit | grep "[0-9][.][0-9]*"`
#cc="0.35"
#echo cc:$cc
#echo ----------------------------------------------------------
#echo -----------------------------$convcrit----------------
#echo "-----------------------------$convcrit----------------"
#echo :"$convcrit":
#echo $convcrit | awk '{printf "%.2f", $1*1.1}'
#echo ----------------------------------------------------------
#echo ">>>>>>>>>>>>>>>>>>>>>>"
##python $avg_dUdL_converged $convcrit
#python $avg_dUdL_converged $cc
#echo ">>>>>>>>>>>>>>>>>>>>>>"
##jo=`python $avg_dUdL_converged $convcrit | xargs`
#echo :$convcrit:
#echo :$cc:
#echo ppppppppppppppppppppppppppppppppppppppp
#echo "$convcrit"..
#echo "$cc"..
#echo 88888888888888888888888888888888888888888888
#[ "$convcrit" = "$cc" ] && echo GLLLLEICHEEEEEEEEEEEE
#[ "$convcrit" != "$cc" ] && echo NOOOOOOOOOOOOOOOOOOOOOOOOO
#jo=`python $avg_dUdL_converged $cc | xargs`
#echo jo:$jo
#exit
## high
[ ! -e `which $ion_energies` ] && echo $ion_energies does not exist && exit -1    ## removes old auswertung by itself  
                                                                                  ## (avg_dUdL_high_{fre,ene,eS0}, avg_dUdL_lowplushigh_{fre,ene,eS0}, AllEn, summery_highUp)

[ ! -e `which $fit_up` ] && echo $fit_up does not exist && exit -1                ## removes old auswertung by itself  
                                                                                  ## (Fah_high, Fah_lowplushigh, fit*)

[ ! -e `which $auswertung` ] && echo $auswertung does not exist && exit -1        ## removes old auswertung by itself  
                                                                                  ## (auswertung_high/ auswertung_low_plus_high/ auswertung_low_plus_high_old/}

##################################################################################################################################################################
## some settings
##################################################################################################################################################################
### complete_new: this sets, that output of this skript is delted first
complete_new=yes;[ "$1" = "no" ] && complete_new=no
all_in=`echo " $* "`
rm=no;rm_tag=`echo "$all_in" | grep -o " -rm "`;[ "$rm_tag" = " -rm " ] && rm=yes
rm=yes

h=yes;l=yes;
h_tag=`echo "$all_in" | grep -o " -h "`;[ "$h_tag" = " -h " ] && l=no ; #echo h:$h h_tag:$h_tag: 
l_tag=`echo "$all_in" | grep -o " -l "`;[ "$l_tag" = " -l " ] && h=no ; #echo l:$l l_tag:$l_tag: 


##################################################################################################################################################################
## check if input folders defined ( z.B. 3.74Ang_1100K or low_2x2x2sc_280eV_3x3x3kp_120NGXF/ low_3x3x3sc_260eV_2x2x2kp_240NGXF__high_400eV_4x4x4kp_240NGXF/ )
##################################################################################################################################################################
lowfolder_input=`echo $* | xargs -n1 | grep "^low_.*sc_.*eV.*kp" | sed 's|.*_high_.*||' | sed '/^$/d'`

highfolder_input=`echo $* | xargs -n1 | grep "^low_.*sc_.*eV.*kp_.*_high_"`
angKfolder_input=`echo $* | xargs -n1 | grep "^[0-9.]*Ang_[0-9]*K"`  ## ob dies fuer high oder fuer low ist weiss man hier noch nicht
if [ `getOption -a` = True ];then
    echogreen "Running qstat to get your currently running jobs ..... this may take a while ..."
    #angKfolder_input_1=`qstatn -q -r -n` # | awk '{print $NF}' | grep "/[0-9.]*Ang_[0-9.]*K/lambda[0-9.]*_[0-9.]*" | sed 's|lambda.*||' | sort | uniq` ## Dies schliesst highfolder aus
    angKfolder_input_1=`qstatn -r -n` # | awk '{print $NF}' | grep "/[0-9.]*Ang_[0-9.]*K/lambda[0-9.]*_[0-9.]*" | sed 's|lambda.*||' | sort | uniq` ## Dies schliesst highfolder aus
    angKfolder_input=`echo "$angKfolder_input_1" | awk '{print $NF}' | grep "/[0-9.]*Ang_[0-9.]*K/lambda[0-9.]*_[0-9.]*" | sed 's|lambda.*||' | sort | uniq` ## Dies schliesst highfolder aus
    ## make sure to include to include input_old .... but only once and not forever .... therefore better not first of all
    #angKfolder_input=`echo "$angKfolder_input" "$angKfolder_input_old" | xargs -n1 | sort | uniq`  
    ti6folder=`echo "$angKfolder_input" | sed 's|[0-9.]*Ang_[0-9.]*K/||' | sort | uniq`   # in this low folder ti_low_6....sh will be run
    lowfolder_input=`echo "$angKfolder_input" | head -1`
    lowfolder_input=$do
    allhighfolder_a=""
    ### get highfolder
    ti6folderhigh=""
    for jjj in $angKfolder_input;do
        #lowcheck=`echo $jjj | sed 's|\(.*low_.*\)/\([0-9.]*Ang_[0-9]*K\)/$|\1__high_|'`
        check_folder=`echo $jjj | sed 's|\(.*\)\(.*low_.*\)/\([0-9.]*Ang_[0-9]*K\)/$|\1|'`
        check_low=`echo $jjj | sed 's|\(.*\)\(.*low_.*\)/\([0-9.]*Ang_[0-9]*K\)/$|\2|'`
        check_low_angk=`echo $jjj | sed 's|\(.*\)\(.*low_.*\)/\([0-9.]*Ang_[0-9]*K\)/$|\3|'`
        [ "`getOption -v`" = "True" ] && echo jjj: $jjj  check_low_angk:$check_low_angk:
            check_find_high=`find -L $check_folder -maxdepth 1 -mindepth 1 -type d -name "$check_low\__high*"`
            for cl in $check_find_high;do   # loop over all possible highfolder
                [ "`getOption -v`" = "True" ] && echo cl: $cl
                check_find_high_angk=`find -L $cl -maxdepth 1 -mindepth 1 -type d -name "$check_low_angk"`
                # now we only have to look for corresponding check_angk.... other folders should not be evaluated
                for kkk in $check_find_high_angk;do # loop over all possible high angk folder
                    allhighfolder_a="$allhighfolder_a $kkk"
                    ti6folderhigh="$ti6folderhigh $cl"
                    [ "`getOption -v`" = "True" ] && echo kkk: $kkk
                done
            done

        [ "`getOption -v`" = "True" ] && echo && echo ----------------------
    done
fi
allhighfolder_a=`echo $allhighfolder_a | xargs -n1`
ti6folderhigh=`echo $ti6folderhigh | xargs -n1 | sort | uniq`
[ "`getOption -v`" = "True" ] && echo "###################################" && echo "$allhighfolder_a" && echo "############################"
## wenn angKfolder_input gegeben sind, muss man sich schon in einem lowfolder befind -Len!

##################################################################################################################################################################
## get all lowfolder  ( low_2x2x2sc_280eV_4x4x4kp_120NGXF )
## get all highfolder ( low_2x2x2sc_260eV_4x4x4kp_120NGXF__high_400eV_kp06m06m06_NGXF120)
## get all angKfolder kommt erst spaeter ( wie z.B. 3.74Ang_1100K ...)
##################################################################################################################################################################
startfolder=`pwd | sed 's|.*/||'`

## startfolder = TIfolder
[ "`echo $startfolder | grep "^ti" | wc -w | sed 's|[ ]*||g'`" = "1" ] && all_lowfolder=`ls -1d low*sc_*eV*kp* | sed 's|.*_high_.*||' | sed '/^$/d'` && all_highfolder=`ls -1d *sc_*eV*kp* | grep "^low_.*sc_.*eV.*kp_.*_high_"`
[ "`echo $startfolder | grep "^ti" | wc -w | sed 's|[ ]*||g'`" = "1" ] && [ "$lowfolder_input" != "" ] || [ "$highfolder_input" != "" ] && all_lowfolder=$lowfolder_input && all_highfolder=$highfolder_input
#[ `getOption -a` != True ] && [ "`echo $startfolder | grep "^ti" | wc -w | sed 's|[ ]*||g'`" = "1" ] && [ "$angKfolder_input" != "" ] && echo you are in ti folder but defined angK folder && exit -1 ## DIES IST ZU SPEZIELL

## startfolder = lowfolder       or        startfolder = highfolder
[ "`echo $startfolder | grep "^low" | wc -w | sed 's|[ ]*||g'`" = "1" ] && all_lowfolder=`pwd` && all_highfolder=""
[ "`echo $startfolder | grep "^low.*_high_" | wc -w | sed 's|[ ]*||g'`" = "1" ] && all_lowfolder="" && all_highfolder=`pwd`

## startfolder = angKfolder -------------------------------- >> V
## angKfolder  = lowfolder or highfolder?  -------------------------------------------------------------- >> V --------------------------------------------- >> V
[ "`echo $startfolder | grep "^[0-9.]*Ang_[0-9]*K" | wc -w | sed 's|[ ]*||g'`" = "1" ] && [ "`pwd | grep "low_.*" | wc -w | sed 's|[ ]*||g'`" = "1" ] && [ "`pwd | grep "low_.*_high_" | wc -w | sed 's|[ ]*||g'`" = "0" ] && all_lowfolder=`pwd`   && all_highfolder=""
[ "`echo $startfolder | grep "^[0-9.]*Ang_[0-9]*K" | wc -w | sed 's|[ ]*||g'`" = "1" ] && [ "`pwd | grep "low_.*" | wc -w | sed 's|[ ]*||g'`" = "1" ] && [ "`pwd | grep "low_.*_high_" | wc -w | sed 's|[ ]*||g'`" = "1" ] && all_lowfolder=""  && all_highfolder=`pwd`
[ "`echo $startfolder | grep "^[0-9.]*Ang_[0-9]*K" | wc -w | sed 's|[ ]*||g'`" = "1" ] && [ "$all_highfolder" != "" ] && echo "start for highfolders in lowfolder not in *Ang*K folder" && exit -1

## if selected just high or just low
[ "$l" = "no" ] && all_lowfolder=""
[ "$h" = "no" ] && all_highfolder=""


### if naming is unknown
[ "$all_lowfolder" = "" ] && [ "$all_highfolder" = "" ] && all_lowfolder=`find -L $dostart -maxdepth 3 -mindepth 1 -type d -name lambda'[0-9.]*_[0-9]*' | grep "[0-9.]*Ang_[0-9]*K/lambda[0-9.]*_[0-9]*" | sed 's|\(.*\)lambda[0-9.]*_[0-9]*.*|\1|' | sort | uniq`
[ "$lowfolder_input" != "" ] && all_lowfolder=$lowfolder_input
[ "`getOption -a`" = "True" ] && [ "$allhighfolder_a" != "" ] && all_highfolder=$allhighfolder_a

## format output
## echo to screen
/bin/echo -e "\033[31m\033[1m   lowfolder_input::\033[0m" $lowfolder_input  ## you can manually select what to evaluate
/bin/echo -e "\033[31m\033[1m high folder_input::\033[0m" $highfolder_input ## you can manually select what to evaluate
/bin/echo -e "\033[31m\033[1m angK folder_input: (in this folders `basename $0` is run) :\033[0m" "\n$angKfolder_input"
/bin/echo -e "\033[31m\033[1m     startfolder::\033[0m" $startfolder
echo ""
/bin/echo -e "\033[31m\033[1m      alllowfolder: (in this folders `basename $0` is run) :\033[0m \n$all_lowfolder"
echo ""
/bin/echo -e "\033[31m\033[1m     allhighfolder: (in this folders `basename $0` is run) :\033[0m \n$all_highfolder"
echo ""
/bin/echo -e "\033[31m\033[1m     ti6folder    : (in this folders $auswertung is run) :\033[0m \n$ti6folder"
echo 
/bin/echo -e "\033[31m\033[1m     ti6folderhigh    : (in this folders $auswertung is run) :\033[0m \n$ti6folderhigh"
echo 
echo xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx:
echo "rm:   $rm"
echo "high: $h"
echo "low:  $l"
echo xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx





##############################################################################################################################################################t#
##############################################################################################################################################################t#
##############################################################################################################################################################t#
################################################################################################################################################################
################################################################################################################################################################
## schleife ueber all lowfolder  ( low_2x2x2sc_280eV_4x4x4kp_120NGXF )
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
############################################################################
## Print all lowfolder to screen
############################################################################
if [ "`getOption -oh`" = "True" ];then
    echo "NOT RUNNING LOWFOLDER DUE TO OPTION -oh"
else
  if [ "$all_lowfolder" = "." ] || [ "`echo $all_lowfolder | wc -w | sed 's|[ ]*||g'`" -lt "2" ];then noprint=noprint; else 
  /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
  /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
  for lowfolder in $all_lowfolder;do /bin/echo -e "\033[31m\033[1m### lowfolder: $lowfolder\033[0m $host"; done
  /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
  /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
  fi
hier=`pwd`
anz_lowfolder=`echo $all_lowfolder | wc -w | sed 's|[ ]*||g'`
nr_lowfolder=0
[ "$all_lowfolder" = "" ] && all_lowfolder=$lowfolder_input
for lowfolder in $all_lowfolder;do ## schleife lowfolder
  nr_lowfolder=` expr $nr_lowfolder + 1 `
  echo cd lowfolder::: cd $lowfolder
  cd $lowfolder                ## lowfolder ist der entsprechende lowfolder
  lowfolder=`pwd`
  if [ "$all_lowfolder" = "." ] || [ "`echo $all_lowfolder | wc -w | sed 's|[ ]*||g'`" -lt "2" ];then noprint=noprint; else 
  /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
  /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
  /bin/echo -e "\033[31m\033[1m#### lowfolder [$nr_lowfolder/$anz_lowfolder]: $lowfolder\033[0m $host"
  /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
  /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
  fi


  ############################################################################
  ## get all_angKfolder  ( 3.74Ang_1360K )
  ############################################################################
  [ "`pwd | sed 's|.*/||' | grep "[0-9.]*Ang_[0-9]*K" | wc -w | sed 's|[ ]*||g'`"  = "1" ] && all_angKfolder=. ## dis checkd ob man gerade in einem ang*K folder ist
  [ "`pwd | sed 's|.*/||' | grep "[0-9.]*Ang_[0-9]*K" | wc -w | sed 's|[ ]*||g'`" != "1" ] && all_angKfolder=`find . -maxdepth 1 -mindepth 1 -type d -name "*Ang_*K" | sed 's|./||' | sed 's|\([0-9.]*\)Ang_\([0-9]*\)K.*|\2 \1Ang_\2K|' | sort -n | awk '{print $2}'`
  [ "$angKfolder_input" != "" ] && [ "$all_angKfolder" != "." ] && all_angKfolder=$angKfolder_input
  [ "`getOption -a`" = "True" ] && all_angKfolder=$angKfolder_input

  ############################################################################
  ## Print all angKfolder to screen
  ############################################################################
  if [ "$all_angKfolder" != "." ];then 
  /bin/echo -e "\033[34m\033[1m##############################################################################################################################\033[0m"
  for angKfolder in $all_angKfolder;do \/bin/echo -e "\033[34m\033[1m### angKfolder: $angKfolder in [$nr_lowfolder/$anz_lowfolder]\033[0m $host";done
  /bin/echo -e "\033[34m\033[1m##############################################################################################################################\033[0m"
  fi
  anz_angKfolder=`echo $all_angKfolder | wc -w | sed 's|[ ]*||g'`
  nr_angKfolder=0



  ##########################################################################
  #### schleife ueber alle AngK folder
  ##########################################################################
  IDswhichreallytokill=""
  zuvorangK=`pwd`
  for angKfolder in $all_angKfolder;do ## schleife angKfolder
      cd $zuvorangK
      cd $angKfolder     ## angK ist der entsprechende angK folder
      angKfolder=`pwd`
      ang=`echo $angKfolder | sed 's|.*/||' | sed 's|Ang.*||'`
      K=`echo $angKfolder | sed 's|.*/||' | sed 's|.*_||' | sed 's|K||'`
      nr_angKfolder=` expr $nr_angKfolder + 1 `

      [ "$firsttime" = "yes" ] && [ "`getOption -ga`" = "True" ] && [ "$nr_angKfolder" -lt "$angKfolderstart" ] && echo "nr_angKfolder:$nr_angKfolder: angKfolderstart:$angKfolderstart:" && continue
      firsttime="no"
          
  ##########################################################################
  #### low no mathematica (possible at cmmd010 001)
  ##########################################################################
      /bin/echo -e "\033[34m\033[1m##############################################################################################################################\033[0m"
      /bin/echo -e "\033[34m\033[1m### angKfolder($nr_angKfolder/$anz_angKfolder) in [$nr_lowfolder/$anz_lowfolder]: $angKfolder \033[0m $host"; 
      /bin/echo -e "\033[34m\033[1m##############################################################################################################################\033[0m"

      ## last test before skript
      ## last test before skript
      [ "`pwd | sed 's|.*/||' | grep "[0-9.]*Ang_[0-9]*K" | wc -w | sed 's|[ ]*||g'`" != "1" ] && echo "`pwd` doesnt seem to be *ang*K folder" && continue

      ## RUN SKRIPTS
      ## RUN SKRIPTS
      
      echo "##############################################################################################################################"
      echo "### angKfolder($nr_angKfolder/$anz_angKfolder) in [$nr_lowfolder/$anz_lowfolder]: `which $dUdLallInfo` $host"
      /bin/echo -e "\033[31m\033[1m### $angKfolder\033[0m"
      echo "##############################################################################################################################"
      for i in `ls lambda[0-9.]*_[0-9]* -1d`; do cd $i;
          #echo `pwd`;
          $dUdLallInfo; 
          cd $angKfolder; done
      
      echo "##############################################################################################################################"
      echo "### angKfolder00($nr_angKfolder/$anz_angKfolder) in [$nr_lowfolder/$anz_lowfolder]: `which $removeJumps` THIS qdels jobs if jump"            ## /dUdL_create_dUdLallInfo_noJumps.sh -q
      /bin/echo -e "\033[31m\033[1m### $angKfolder\033[0m $host"
      echo "##############################################################################################################################"
      [ "`echo "$*" | grep -o "\-n" | wc -w | sed 's|[ ]*||g'`" = "0" ] && $removeJumps #-d ##-q ## -q suppresses standard output of $removeJumps but not the errors
##    [ "`echo "$*" | grep -o "\-n" | wc -w | sed 's|[ ]*||g'`" = "0" ] && $removeJumps ##-q ## -q suppresses standard output of $removeJumps but not the errors

      echo "##############################################################################################################################"
      echo "### angKfolder11($nr_angKfolder/$anz_angKfolder) in [$nr_lowfolder/$anz_lowfolder]: `which $getStrucUncor` ang:$ang K:$K"         ## /fah_high_1get_structures_uncor.sh
      /bin/echo -e "\033[31m\033[1m### $angKfolder\033[0m $host"
      echo "##############################################################################################################################"
      cd $angKfolder
      #[ "$lowfolder" != "$angKfolder" ] && cd $lowfolder
      #[ "$lowfolder"  = "$angKfolder" ] && cd ..
      cd ..
      $getStrucUncor $ang $K -a
##    $getStrucUncor $ang $K -a
      cd $angKfolder

      echo "##############################################################################################################################"
      echo "### angKfolder22($nr_angKfolder/$anz_angKfolder) in [$nr_lowfolder/$anz_lowfolder]: `which $getdUdLs`" ## /home/glensk/scripts/ThermodynamicIntegration/AH_getdUdLs.sh
        /bin/echo -e "\033[31m\033[1m### $angKfolder\033[0m $host"
      echo "##############################################################################################################################"
      [ "$rm" = "yes" ] && rm -f */getAvgAndDev*
      [ "$rm" = "yes" ] && rm -f dUdL_all fit* avg_dUdL_*   # dudl all contains info about seed and step
      $getdUdLs -ns  ## -ns removes: dUdLallInfo_noJumps not found        ;
##    [ "$rm" = "yes" ] && rm -f */getAvgAndDev*
##    [ "$rm" = "yes" ] && rm -f dUdL_all fit* avg_dUdL_*
##    $getdUdLs -ns  ## -ns removes: dUdLallInfo_noJumps not found        ;
      ##########################################################################
      #### low !!mathematica!! (NOT possible at cmmd010 001)
      ##########################################################################

      echo "##############################################################################################################################"
      echo "### angKfolder33($nr_angKfolder/$anz_angKfolder): in [$nr_lowfolder/$anz_lowfolder] `which $fit` "  ## fah_low_fit_avg_dUdL.sh
        /bin/echo -e "\033[31m\033[1m### $angKfolder\033[0m $host"
      echo "##############################################################################################################################"
##    [ "$rm" = "yes" ] && rm -f Fah Fah.old
##    $fit
      [ "$rm" = "yes" ] && rm -f Fah Fah.old
      $fit

    if [ -e "avg_dUdL_fre" ];then
        ##########################################################################
        #### check if some lambdas are already converged and could be killed, convcrit is set at the beginning of this skript
        ##########################################################################
        # convcrit local

        if [ -e "convcrit" ];then
            convcrit_take=`cat convcrit`
            [ "`isnumber.sh $convcrit_take`" != "yes" ] && convcrit_take=$convcrit
            echo "CONFCRIT FILE FOUND: convcrit: $convcrit_take"
        else
            if [ -e "../convcrit" ];then
                convcrit_take=`cat ../convcrit`
                [ "`isnumber.sh $convcrit_take`" != "yes" ] && convcrit_take=$convcrit
                echo "CONFCRIT FILE FOUND: ../convcrit: $convcrit_take"
            else 
                convcrit_take=$convcrit
            fi
        fi
        #convcrit_take=0.35
        echo "----> convcrit_take:$convcrit_take<----"
        
        lambdastokill=`python $avg_dUdL_converged $convcrit_take | xargs`

        #lambdastokill1=`grep -v "^#" avg_dUdL_fre | sed 's|^[ ]*||g'`
        #echo lambdastokill1:
        #echo "$lambdastokill1"

        #lambdastokill2=`grep -v "^#" avg_dUdL_fre | sed 's|^[ ]*||g' | awk '{print $1,$3,'"`echo $convcrit_take`"'}'`
        #echo lambdastokill2:
        #echo "$lambdastokill2"


        ##lambdastokill3=`grep -v "^#" avg_dUdL_fre | sed 's|^[ ]*||g' | awk '{print $1,$3}' | sed 's|\(.*\)|\1 '"$convcrit_take"'|' | awk '{print $1,$2,$3}'`
        ##echo lambdastokill3:
        ##echo "$lambdastokill3"

        #lambdastokill4=`echo "$lambdastokill2" | awk '{print $1,$2,$3,$2+$1+$3}'`  #'($2<=$3){print $1};($3>$2){print "no"}'`
        #echo lambdastokill4:
        #echo "$lambdastokill4"

        # lambdastokill=`grep -v "^#" avg_dUdL_fre | sed 's|^[ ]*||g' | awk '{print $1,$3}' | sed 's|\(.*\)|\1 '"$convcrit_take"'|' | awk '{print $1,$2,$3}' | awk '$2<=$3{print $1}'`
        #lambdastokill=`grep -v "^#" avg_dUdL_fre | sed 's|^[ ]*||g' | awk '$3<='"$convcrit_take"'{print $1}' | xargs`
        echo lambdastokill:$lambdastokill

        for lamtok in $lambdastokill;do
            IDstokill=`OUTCAR_ID.sh lambda$lamtok\_* | xargs`
            echo IDstokill:$IDstokill:
		## check if we have enough steps in those ID's
		for id in `echo $IDstokill`;do
			idtake=""
			echo id:$id:
			steps=`cat dUdL_all | awk '{if ($NF=='"$id"') print $10}'`
			echo steps:$steps 
			[ "$steps" != "" ] && [ "$steps" -ge "301" ] && idtake=$id
			IDswhichreallytokill="$IDswhichreallytokill $idtake"
		done

            IDstokillall=`echo $IDstokillall $IDstokill | xargs`                

            #echo IDstokillall:$IDstokillall
            #seeds=`ls -1d lambda$lamtok\_* | sed 's|lambda||' | sed 's|.*_||' | xargs`
            #echo seedstokill: $seeds
            done
            #[ "`echo $IDstokillall | wc -w | sed 's|[ ]*||g'`" != "0" ] && [ "$qdel" != "0" ] && $qdel $IDstokillall
        fi
        echo "IDs which could be killed since those are converged (due to convcrit)":$IDstokillall
  cd $lowfolder
  done ## schleife angKfolder

  if [ "`hostname`" == "cmmc001" -o "`hostname`" == "cmmc002" ]; then
    IDstokillall=`echo "$IDstokillall" | grep -o "[0-9]*" | xargs`
    echo 
    echo
    echo ----------------------------------------------
    echo "IDstokillallreal (due to convcrit and steps > 300): " $IDswhichreallytokill            
    [ "`echo $IDswhichreallytokill | wc -w | sed 's|[ ]*||g'`" != "0" ] && [ "$qdel" != "0" ] && $qdel $IDswhichreallytokill  
    echo ----------------------------------------------
    echo ok
  fi

  ausw_low=0
  [ "`pwd | sed 's|.*/||' | grep "^low" | wc -w | sed 's|[ ]*||g'`" = "1" ] && ausw_low=1
  checkfolder=`pwd`
  checkangk=`find -L $checkfolder -maxdepth 1 -mindepth 1 -type d -name '[0-9.]*Ang_*[0-9]*K' | wc -l | sed 's|[ ]*||g'`
  [ "$checkangk" != "0" ] && ausw_low=1
  checklam=`find -L $checkfolder -maxdepth 1 -mindepth 1 -type d -name lambda'[0-9.]*_[0-9]*' | wc -l | sed 's|[ ]*||g'`
  [ "$checklam" != "0" ] && cd .. && ausw_low=1
  if [ "$ausw_low" = "1" ];then
  echo "##############################################################################################################################"
  echo "### angKfolder($nr_angKfolder/$anz_angKfolder): in [$nr_lowfolder/$anz_lowfolder] `which $auswertung` $host"
  echo "##############################################################################################################################"
  ausw_low=1
  $auswertung
  fi


  cd $hier ## in den ausgangsfolder
done ## schleife all_lowfolder
fi



















#################################################################################################################################################################
#################################################################################################################################################################
#################################################################################################################################################################
#################################################################################################################################################################
#################################################################################################################################################################
#### schleife ueber alle high folder (low_2x2x2sc_280eV_3x3x3kp_120NGXF__high_400eV_6x6x6kpm0_120NGXF_ED1E-3/)
#################################################################################################################################################################
#################################################################################################################################################################
#################################################################################################################################################################
#################################################################################################################################################################
#################################################################################################################################################################
############################################################################
## Print all high to screen
############################################################################
if [ "`getOption -ol`" = "True" ];then
    echo "NOT RUNNING HIGHFOLDER DUE TO OPTION -ol"
else
    if [ "$all_highfolder" = "." ] || [ "`echo $all_$highfolder | wc -w | sed 's|[ ]*||g'`" = "1" ];then noprint=noprint; else ## else print
    /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
    /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
    for highfolder in $all_highfolder;do /bin/echo -e "\033[31m\033[1m### highfolder: $highfolder\033[0m"; done
    /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
    /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
    fi

vorhigh=`pwd`
anz_highfolder=`echo $all_highfolder | wc -w | sed 's|[ ]*||g'`; [ "$anz_highfolder" = "0" ] && anz_highfolder=1
nr_highfolder=0
for highfolder in $all_highfolder;do
    [ "$highfolder" = "." ] && highfolder=`pwd`
    nr_highfolder=` expr $nr_highfolder + 1 `
    /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
    /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
    /bin/echo -e "\033[31m\033[1m#### highfolder [$nr_highfolder/$anz_highfolder]: $highfolder\033[0m"
    /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"
    /bin/echo -e "\033[31m\033[1m##############################################################################################################################\033[0m"

    cd $highfolder
    rm -f highUp_summary

    ############################################################################
    ## get all_angKfolderhigh  ( 3.74Ang_1360K )
    ############################################################################
    [ "`pwd | sed 's|.*/||' | grep "[0-9.]*Ang_[0-9]*K" | wc -w | sed 's|[ ]*||g'`"  = "1" ] && all_angKfolderhigh=. ## dis checkd ob man gerade in einem ang*K folder ist
    [ "`pwd | sed 's|.*/||' | grep "[0-9.]*Ang_[0-9]*K" | wc -w | sed 's|[ ]*||g'`" != "1" ] && all_angKfolderhigh=`ls -1d *Ang_*K | sed 's|\([0-9.]*\)Ang_\([0-9]*\)K.*|\2 \1Ang_\2K|' | sort -n | awk '{print $2}'`
    all_angKfolderhigh="-a"
    [ "$angKfolder_input" != "" ] && [ "$all_angKfolderhigh" != "." ] && all_angKfolderhigh=$angKfolder_input
    [ "`getOption -a`" = "True" ] && all_angKfolderhigh="-c"
    [ "`getOption -v`" = "True" ] && echo jo: `pwd`
    [ "`getOption -v`" = "True" ] && echo ja:$all_angKfolderhigh:


    /bin/echo -e "\033[34m\033[1m##############################################################################################################\033[0m"
    /bin/echo -e "\033[34m\033[1m## [$nr_highfolder/$anz_highfolder]: $ion_energies ( `pwd` )\033[0m"
    /bin/echo -e "\033[31m\033[1m#### highfolder [$nr_highfolder/$anz_highfolder]: $highfolder\033[0m"
    /bin/echo -e "\033[34m\033[1m##############################################################################################################\033[0m"
    # removes old stuff on its own
    # ti_high_1_get_ion_energies.sh
    $ion_energies $all_angKfolderhigh
    
    /bin/echo -e "\033[34m\033[1m##############################################################################################################\033[0m"
    /bin/echo -e "\033[34m\033[1m## [$nr_highfolder/$anz_highfolder]: $fit_up ( `pwd` )\033[0m"
    /bin/echo -e "\033[31m\033[1m#### highfolder [$nr_highfolder/$anz_highfolder]: $highfolder\033[0m"
    /bin/echo -e "\033[34m\033[1m##############################################################################################################\033[0m"
    ## removes old stuff on its own
    # ti_high_2_highFitUp.sh
    $fit_up $all_angKfolderhigh
    
    # mache generell die auswertung in highfolder
    $auswertung
    if [ "`echo $all_angKfolderhigh | wc -w`" -ge "2" ];then
      /bin/echo -e "\033[34m\033[1m##############################################################################################################\033[0m"
      /bin/echo -e "\033[34m\033[1m## [$nr_highfolder/$anz_highfolder]: $auswertung ( `pwd` )\033[0m"
      /bin/echo -e "\033[31m\033[1m#### highfolder [$nr_highfolder/$anz_highfolder]: $highfolder\033[0m"
      /bin/echo -e "\033[34m\033[1m##############################################################################################################\033[0m"
      ## removes old stuff on its own
      ## auswertung kann immer gemacht werden wenn high angKfolder erneuert werden
      [ "`getOption -a`" != "True" ] && $auswertung
      fi


    cd $vorhigh
    done ## schleife ueber all_highfolder

    if [ "`getOption -a`" = "True" ];then
        for kkk in $ti6folderhigh;do
            cd $vorhigh
            cd $kkk
            $auswertung
            cd $vorhigh
            done
        fi


fi

cd $dostart
done ## schleife ueber seperat definierte folder


echo ""
echo yes

# ================== Papers to cite ============
echo
echored "==================================== Papers to cite =========================================================================="
echored "When using UP-TILD, pls. cite:"
echo "B. Grabowski, L. Ismer, T. Hickel, and J. Neugebauer, Phys. Rev. B 79, 134106 (2009)."
echored "When using TU-TILD, pls. cite:"
echo "A. I. Duff, T. Davey, D. Korbmacher, A. Glensk, B. Grabowski, J. Neugebauer, and M. W. Finnis, Phys. Rev. B 91, 214311 (2015)."
echo "A. I. Duff, M. W. Finnis, P. Maugis, B. J. Thijsse, and M. H. F. Sluiter, Comput. Phys. Commun. 196, 439 (2015)."
echored "==================================== Papers to cite =========================================================================="
echo
# ================== Papers to cite ============

[ "$continuesmode" != "yes" ] && echo "no continues mode" && exit
if [ "$continuesmode" = "yes" ];then
    zuvor=`pwd`
    for jjj in $ti6folder;do
        echo "FOLDER: ---> make low auswertung in $jjj"
        cd $jjj
        $auswertung
        cd $zuvor
    done
    cd $zuvor
    /bin/echo -e "\033[32m\033[1mvor  sleep ($schlafe)sec `date`\033[0m"
    rmEqw.sh
    sleep $schlafe
    /bin/echo -e "\033[32m\033[1mnach sleep ($schlafe)sec `date`\033[0m"
fi
done   # endlosschleife
