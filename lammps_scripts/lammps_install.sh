#!/bin/sh

folder_lammps_sources=$HOME/sources
#[ "$1" != "" ] && folder_lammps_sources=$1
#echo folder_lammps_sources $folder_lammps_sources
#exit

[ ! -e "$folder_lammps_sources" ] && read yn && mkdir $folder_lammps_sources
cd $folder_lammps_sources
pwd
[ ! -e "lammps_source_cosmo" ] && git clone https://github.com/cosmo-epfl/lammps.git lammps_source_cosmo && echo `date +"%Y_%m_%d"` > ANMERKUNG.txt
cd lammps_source_cosmo
git checkout runner-lammps
cd src
pwd
make yes-CLASS2 yes-KSPACE yes-MANYBODY yes-MISC yes-MOLECULE yes-REPLICA yes-RIGID yes-USER-MISC
make yes-USER-RUNNER
sed -i 's|^#define MAXNEIGH.*|#define MAXNEIGH 500|' pair_runner.h
make serial
