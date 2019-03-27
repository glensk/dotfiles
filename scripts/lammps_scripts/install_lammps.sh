#!/bin/sh
#############################################################
# Start editing here ########################################
#############################################################
lammps="original" # "original" or "cosmo"
#lammps="cosmo" # "original" or "cosmo"

n2p2_folder=$HOME/Dropbox/Albert/git/n2p2  # is no =""

installfolder=$HOME/sources         # where to clone the lammps folder 
move_exec_to=$scripts/lammps_executables  # if exec should be moved



makefile=""         # "" means no special/own makefile
makeversion="fidis"   # "serial" or "mpi" or "fidis"; 
                      # for serial and n2p2 enable -DNOMPI
[ "$1" != "" ] && makeversion=$1

lammpsfolder=lammps_$lammps\_$makeversion
#############################################################
# Stop editing here #########################################
#############################################################
[ "$makeversion" == "fidis" ] && [ "$makeversion" == "mpi" ] && echo "makeversion has to be fidis or mpi" && exit

[ "`hostname`" = "fidis" ] && [ "$makeversion" == "fidis" ] && makefile=$dotfiles/scripts/lammps_makefiles/fidis_deneb_2018-10-31/MINE

[ "`hostname`" = "fidis" ] && [ "$makeversion" == "fidis" ] && [ ! -e "$makefile" ] && echo "makefile $makefile does not exist" && exit 

src=$installfolder/$lammpsfolder/src
src_=$installfolder/$lammpsfolder
[ "$lammps" == "cosmo" ] && gitfrom="https://github.com/cosmo-epfl/lammps.git"
[ "$lammps" == "original" ] && gitfrom="https://github.com/lammps/lammps.git"
[ "$n2p2_folder" != "" ] && [ ! -e "$n2p2_folder" ] && echo n2p2 folder $n2p2_folder not found && exit

echo "----------------------------------------------------------------------"
echo "lammps       : $lammps (original or cosmo)"
echo "installfolder: $installfolder"
echo "lammpsfolder : $lammpsfolder"
echo "move_exec_to : $move_exec_to"
echo "makeversion  : $makeversion (mpi or fidis or serial)"
echo "makefile     : $makefile (in case of custome makefile)"
echo "n2p2_folder  : $n2p2_folder"
echo "src          : $src"
echo "----------------------------------------------------------------------"

#[ -e "$src" ] && echo "$src folder exists;Exit!" && exit
#### check if installfolder exists
if [ ! -e "$installfolder" ];then
    read -p "Should I crate the folder $installfolder ? [yes y no n] " yn
    case $yn in
        [Yy]* ) mkdir $installfolder; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
fi

#### git clone
cd $installfolder
echo pwd: `pwd`
[ ! -e "$lammpsfolder" ] && git clone $gitfrom $lammpsfolder 
cd $lammpsfolder 
echo `date +"%Y_%m_%d"` > ANMERKUNG.txt 
[ "$lammps" = "cosmo" ] && git checkout runner-lammps

cd $src
make clean-all
make yes-CLASS2 yes-KSPACE yes-MANYBODY yes-MISC yes-MOLECULE yes-REPLICA yes-RIGID yes-USER-MISC yes-MEAM

if [ "$n2p2_folder" != "" ];then
    cd $src_
    echo ""
    echo "get n2p2 lib/nnp and USER-NNP to src"
    ln -s $n2p2_folder lib/nnp
    cp -r $n2p2_folder/src/interface/LAMMPS/src/USER-NNP src
    cd $src
    make yes-user-nnp
fi

if [ "$lammps" = "cosmo" ];then 
    cd $src
    echo""
    echo "make yes-USER-RUNNER"
    make yes-USER-RUNNER  
    sed -i 's|^#define MAXNEIGH.*|#define MAXNEIGH 500|' pair_runner.h
fi


#### load modules on fidis
if [ "`hostname`" = "fidis" ];then
    cd $src
    #conda deactivate
    [ "$makeversion" == "fidis" ] && cp -r $makefile MAKE # !!! copy the makefile
    echo "now loading the modules ..."
    source $MODULESHOME/init/bash
    module purge
    module load intel
    module load intel-mpi
    module load intel-mkl
    module load fftw
    module load python/2.7.14
    module load gsl   # necessary to have it similarly to n2p2
    module load eigen # necessary to have it similarly to n2p2
fi


#### mow make
cd $src
echo "now make"
pwd
rm -f lmp_$makeversion
make $makeversion | tee -a make_$makeversion\_out_`date +"%Y_%m_%d"`
#make mode=shlib $makeversion
#make mode=shlib fidis
pwd
exit



[ ! -e "lmp_$makeversion" ] && echo "lmp_$makeversion was NOT CREATED! THE COMPILATION FAILED!"
if [ -e "lmp_$makeversion" ];then
    echo "lmp_$makeversion SUCCESSFULLY compiled!"
    
    # move the makefile if possible  
    if [ -e "$move_exec_to" ];then
        echo mv lmp_$makeversion $move_exec_to/lmp_$makeversion\_`hostname`_runner_`date +"%Y_%m_%d"`
        mv lmp_$makeversion $move_exec_to/lmp_$makeversion\_`hostname`_runner_`date +"%Y_%m_%d"`
    fi
fi

# mac:      done    (serial) WORKS, energy checked
# cosmopc:  done    (serial) WORKS, energy checked 
# fidis:    done    (parallel since makefile exists), check weahther it works by redoing one oldjob
#
# in general, dont remove the ~/sources/lammps_source_cosmo folder 
