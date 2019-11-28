#!/bin/bash
    
    [ ! -e /data/$USER ] && echo you need to have /data/$USER mounted && exit
    
    submithost_fallback=cmmd002.mpie.de
    [ ! -z "$submithost" ] && submithost=$submithost   ## this set $submithost if it is defined as evn variable
    [ -z "$submithost" ] && submithost=$submithost_fallback   ## this set $submithost if not defied as env_variable
    #echo ss:$submithost

    ###### ANMERKUNG ##################
    ### Wenn der intelcompiler geleaden wird wenn man auf den $submithost einloggt sollte dass so gehen
    ###### ANMERKUNG ##################
    echo "... make_Fx_files_equal_length.sh ..."
    make_Fx_files_equal_length.sh
    echo "... make_Fx_files_equal_length.sh .................................................................... DONE "
    echo 
    ##############################################
    ## GetGibbsEnergyOfFormation.sh
    ##############################################
    time=`date +%s`; 
    tmpfolder=/data/$USER/tmp_results_$time
    tmpfolder_cmmd=/home/$USER/tmp_results_$time
    mkdir $tmpfolder
    cp * $tmpfolder 2> /dev/null

   # ssh $submithost "cd $tmpfolder_cmmd; /home/grabowski/Thermodynamics/getGibbsEnergyOfFormation.sh"
   # ssh $submithost "cd $tmpfolder_cmmd;module load intelcompiler/11.1; /home/glensk/scripts/Thermodynamics/getGibbsEnergyOfFormation.sh"
    echo "... go to $submitfolder and do `basename $0` ..."
    ssh $submithost "cd $tmpfolder_cmmd;module load intelcompiler/11.1; /home/$USER/Thermodynamics/getGibbsEnergyOfFormation.sh"
    echo "... go to $submitfolder and do `basename $0` ............................................................. DONE"
    cp -R $tmpfolder/* . 2> /dev/null
    rm -rf $tmpfolder
    ##############################################
    ## GetGibbsEnergyOfFormation.sh
    ##############################################
