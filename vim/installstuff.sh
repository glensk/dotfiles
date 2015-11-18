#!/bin/sh

install_ctags="yes"

goto=$HOME/scripts/dotfiles/vim
[ ! -e "$goto" ] && echo goto: $goto does not exist && exit
cd $goto


#rm -rf vim-pathogen/ autoload bundle
#mkdir -p autoload bundle
#git clone git://github.com/tpope/vim-pathogen/
#cp vim-pathogen/autoload/* autoload
#rm -rf vim-pathogen/

cd bundle
#git clone git://github.com/altercation/vim-colors-solarized.git
#git clone https://github.com/vim-scripts/restore_view.vim.git
#git clone git://github.com/vim-scripts/Efficient-python-folding.git
#git clone https://github.com/scrooloose/syntastic.git 
#git clone git://github.com/scrooloose/nerdtree.git
#git clone https://github.com/jcf/vim-latex.git
#git clone https://github.com/kshenoy/vim-signature.git


#git clone git://github.com/klen/python-mode.git
#git clone git://github.com/vim-scripts/taglist.vim.git
#git clone git://github.com/scrooloose/nerdcommenter.git
#git submodule add https://github.com/jcf/vim-latex.git dotfiles/vim/bundle/vim-latex

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
