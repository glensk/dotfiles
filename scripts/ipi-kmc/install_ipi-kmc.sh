#!/bin/sh

# folder where to save i-pi-kmc
folder_sources=$HOME/sources

# define particular other older for some users
[ ! -z "$scripts" ] && folder_sources=$scripts


echo "-----------------------------------------------------------------------------------"
echo "folder_sources: $folder_sources"
echo "-----------------------------------------------------------------------------------"

if [ ! -e "$folder_sources" ];then
    read -p "Should I crate the folder $folder_sources ? [yes y no n] " yn
    case $yn in
        [Yy]* ) mkdir $folder_sources; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
fi

[ -e "$folder_sources/i-pi-mc" ] && echo "$folder_sources/i-pi-mc does already exist; Exit" && exit

cd $folder_sources
[ ! -e "ipi-mc" ] && echo "git clone https://github.com/ceriottm/i-pi-mc" && git clone https://github.com/ceriottm/i-pi-mc i-pi-mc
cd i-pi-mc
git checkout kmc-al6xxx
