#!/bin/sh
#############################################################
# Start editing here ########################################
#############################################################
folder_lammps_sources=$HOME/sources
folder_lammps_save=$scripts/lammps_executables
#folder_lammps_change=$scripts/lammps_scripts/change_src/

makeversion="serial"  # for n2p2 with serial the -DNOMPI needs to be enabled
makefile=""
[ "`hostname`" = "fidis" ] && makeversion="fidis" && makefile=$scripts/lammps_makefiles/fidis_deneb_2018-10-31/MINE
#############################################################
# Stop editing here #########################################
#############################################################


#[ ! -e "$folder_lammps_change" ] && echo "$folder_lammps_change does not exist" && exit
#[ ! -e "$folder_lammps_change/dump_xyz.cpp" ] && echo "$folder_lammps_change/dump_xyz.cpp does not exist" && exit
#[ ! -e "$folder_lammps_change/dump_xyz.h" ] && echo "$folder_lammps_change/dump_xyz.h does not exist" && exit

[ "$makeversion" != "serial" ] && [ ! -e "$makefile" ] && echo makefile $makefile not found && exit
echo "-----------------------------------------------------------------------------------"
echo "folder_lammps_sources: $folder_lammps_sources"
echo "folder_lammps_save: $folder_lammps_save"
echo "makeversion: $makeversion (either serial or fidis)"
echo "makefile: $makefile (either empty(==serial) or path to makefile)"
echo "-----------------------------------------------------------------------------------"

if [ ! -e "$folder_lammps_sources" ];then
    read -p "Should I crate the folder $folder_lammps_sources ? [yes y no n] " yn
    case $yn in
        [Yy]* ) mkdir $folder_lammps_sources; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
fi

cd $folder_lammps_sources
pwd
[ ! -e "lammps_source_cosmo" ] && git clone https://github.com/cosmo-epfl/lammps.git lammps_source_cosmo && echo `date +"%Y_%m_%d"` > ANMERKUNG.txt
cd lammps_source_cosmo
git checkout runner-lammps
cd src
#cp $folder_lammps_change/dump_xyz.cpp .
#cp $folder_lammps_change/dump_xyz.h .
pwd
make yes-CLASS2 yes-KSPACE yes-MANYBODY yes-MISC yes-MOLECULE yes-REPLICA yes-RIGID yes-USER-MISC
make yes-USER-RUNNER
make yes-user-nnp   # this will not work if n2pc has not been installed 
sed -i 's|^#define MAXNEIGH.*|#define MAXNEIGH 500|' pair_runner.h
rm -f lmp_$makeversion
[ "`hostname`" = "fidis" ] && cp -r $makefile MAKE && source $MODULESHOME/init/bash && module load intel intel-mpi intel-mkl fftw python/2.7.14
pwd
echo "now make"
make $makeversion | tee -a make_$makeversion\_out_`date +"%Y_%m_%d"`
make mode=shlib $makeversion
pwd
[ ! -e "lmp_$makeversion" ] && echo "lmp_$makeversion was NOT CREATED! THE COMPILATION FAILED!"
if [ -e "lmp_$makeversion" ];then
    echo "lmp_$makeversion SUCCESSFULLY compiled!"
    
    # move the makefile if possible  
    if [ -e "$folder_lammps_save" ];then
        echo mv lmp_$makeversion $folder_lammps_save/lmp_$makeversion\_`hostname`_runner_`date +"%Y_%m_%d"`
        mv lmp_$makeversion $folder_lammps_save/lmp_$makeversion\_`hostname`_runner_`date +"%Y_%m_%d"`
    fi
fi

# mac:      done    (serial) WORKS, energy checked
# cosmopc:  done    (serial) WORKS, energy checked 
# fidis:    done    (parallel since makefile exists), check weahther it works by redoing one oldjob
#
# in general, dont remove the ~/sources/lammps_source_cosmo folder 
