#!/bin/sh

install_ctags="yes"

goto=$HOME/scripts/dotfiles/vim
[ ! -e "$goto" ] && echo goto: $goto does not exist && exit
cd $goto



if [ "$install_ctags" = "yes" ] ;then
cd $goto
folder=ctags-5.8
ctags_installfolder=installfolder
[ ! -e "$folder" ] && tar zxf ctags-5.8.tar.gz
cd $folder
mkdir -p $ctags_installfolder
cd $ctags_installfolder
./../configure --prefix=$goto/$folder/$ctags_installfolder 
make && make install
fi
