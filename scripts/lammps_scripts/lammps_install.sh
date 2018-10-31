#!/bin/sh

folder_lammps_sources=$HOME/sources
folder_lammps_save=$scripts/lammps_executables
folder_lammps_makefiles=$scripts/lammps_makefiles
makeversion="serial"
makefile=""
if [ "$hostname" = "fidis" ] && makeversion="fidis" && makefile=$folder_lammps_makefiles/fidis_deneb_2018-10-31/


[ "$makeversion" != "serial" ] && [ ! -e "$makefile" ] && echo makefile $makefile not found && exit
echo "-----------------------------------------------------------------------------------"
echo "folder_lammps_sources: $folder_lammps_sources"
echo "folder_lammps_save: $folder_lammps_save"
echo "makeversion: $makeversion (either serial or fidis)"
echo "makefile: $makefile (either empty(==serial) or path to makefile)"
echo "-----------------------------------------------------------------------------------"

#@ if [ ! -e "$folder_lammps_sources" ];then
#@     read -p "Should I crate the folder $folder_lammps_sources ? [yes y no n] " yn
#@     case $yn in
#@         [Yy]* ) mkdir $folder_lammps_sources; break;;
#@         [Nn]* ) exit;;
#@         * ) echo "Please answer yes or no.";;
#@     esac
#@ fi
#@ 
#@ cd $folder_lammps_sources
#@ pwd
#@ [ ! -e "lammps_source_cosmo" ] && git clone https://github.com/cosmo-epfl/lammps.git lammps_source_cosmo && echo `date +"%Y_%m_%d"` > ANMERKUNG.txt
#@ cd lammps_source_cosmo
#@ git checkout runner-lammps
#@ cd src
#@ pwd
#@ make yes-CLASS2 yes-KSPACE yes-MANYBODY yes-MISC yes-MOLECULE yes-REPLICA yes-RIGID yes-USER-MISC
#@ make yes-USER-RUNNER
#@ sed -i 's|^#define MAXNEIGH.*|#define MAXNEIGH 500|' pair_runner.h
#@ rm -f lmp_serial
[ "$hostname" = "fidis" ] && 
pwd
exit
make $makeversion | tee -a make_serial_out_`date +"%Y_%m_%d"`
pwd
[ ! -e "lmp_serial" ] && echo "lmp_serial was NOT CREATED! THE COMPILATION FAILED!"
if [ -e "lmp_serial" ];then
    echo "lmp_serial SUCCESSFULLY compiled!"
    if [ -e "$folder_lammps_save" ];then
        echo mv lmp_serial $folder_lammps_save/lmp_serial_`hostname`_runner_`date +"%Y_%m_%d"`
        mv lmp_serial $folder_lammps_save/lmp_serial_`hostname`_runner_`date +"%Y_%m_%d"`
    fi
fi

# mac:      done    (serial) WORKS, energy checked
# cosmopc:  done    (serial) WORKS, energy checked 
# fidis:            (parallel since makefile exists)
# still to anser the question weather lammps sources sould be kept or not ... no ... just delete it!
