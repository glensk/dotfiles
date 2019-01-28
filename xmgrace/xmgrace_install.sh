#!/bin/sh
[ -z "$scripts" ] && echo source source_to_add_to_path.sh first && exit
[ ! -z "$scripts" ] && target_folder=$scripts
[ "$USER" = "glensk" ] && target_folder="$HOME/Dropbox/Albert/git"
ipi_folder_name="i-pi-mc"
echo "-----------------------------------------------------------------------------------"
echo "target_folder: $target_folder"
echo
echo "-----------------------------------------------------------------------------------"

if [ ! -e "$target_folder" ];then
    read -p "Should I crate the folder $target_folder ? [yes y no n] " yn
    case $yn in
        [Yy]* ) mkdir $target_folder; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
fi

[ -e "$target_folder/$ipi_folder_name" ] && echo "$target_folder/$ipi_folder_name does already exist; Exit" && exit

cd $target_folder
[ ! -e "$ipi_folder_name" ] && echo "git clone https://github.com/ceriottm/i-pi-mc" && git clone --depth 1 -b kmc-al6xxx https://github.com/ceriottm/i-pi-mc $ipi_folder_name  # 33 mb
cd $ipi_folder_name
git checkout kmc-al6xxx






git clone --depth 1 https://github.com/fxcoudert/xmgrace 
