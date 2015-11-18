#!/bin/sh

dotfiles=$HOME/Dropbox/scripts/dotfiles
[ ! -e "$dotfiles" ] && echo Folder $dotfiles does not exist && exit

cd $dotfiles

echo "##########"
git init 
echo "##########"


echo "#############"
echo "# git add ..."
echo "#############"
git add bash 
git add bin
git add commands 
git add cron
git add git
git add ipython
#git add ipython-qt-console-solarized-syntax-highlighting/  # has git!
git add screen
git add ssh
git add subversion
git add tcsh

#######################################
#git add vim                 # has git!
#######################################
git add vim/startup
git add vim/spell
git add vim/syntax
git add vim/thesaurus

git add vim/ctags
git add vim/ctags-5.8.tar.gz
git add vim/install_ctags.sh
git add vim/installstuff.sh
git add vim/py.template
git add vim/sh.template
git add vim/test.py
git add vim/vim-shortcuts.png
git add vim/vim_cheatsheet
git add vim/vimrc
git add vim/vimrc_cool
git add vim/vimrc_gerard
git add vim/vimrc_my_old
git add vim/vimrc_my_oldest
git add vim/vimrc_physik_kassel
git add vim/vimrcbackup
git add vim/vimrcbackup2
git add vim/ANMERKUNG_vim_on_mac
git add vim/ANMERKUNG_reinstall_vundle.sh
git add vim/bundle_manual
git add vim/bundle_update_manually.sh
git add vim/colors

echo "#############"
echo "# git submodule add ..."
echo "#############"
git submodule add   https://github.com/gmarik/Vundle.vim.git	vim/bundle/Vundle.vim/
git submodule add   https://github.com/kien/ctrlp.vim.git	vim/bundle/ctrlp.vim/
git submodule add   https://github.com/sjl/gundo.vim	vim/bundle/gundo.vim/
git submodule add   https://github.com/shiblon/latex-makefile.git vim/bundle/latex-makefile/
git submodule add	https://github.com/yegappan/mru.git vim/bundle/mru/
git submodule add	https://github.com/Shougo/neocomplete.vim vim/bundle/neocomplete.vim/
git submodule add	https://github.com/kien/rainbow_parentheses.vim.git vim/bundle/rainbow_parentheses.vim/
git submodule add	https://github.com/ervandew/supertab.git vim/bundle/supertab/
git submodule add	https://github.com/AndrewRadev/switch.vim.git vim/bundle/switch.vim/
git submodule add	https://github.com/scrooloose/syntastic.git vim/bundle/syntastic/
git submodule add	https://github.com/majutsushi/tagbar vim/bundle/tagbar/
git submodule add	https://github.com/tomtom/tcomment_vim.git vim/bundle/tcomment_vim/
git submodule add	https://github.com/bling/vim-airline.git vim/bundle/vim-airline/
git submodule add	https://github.com/altercation/vim-colors-solarized.git vim/bundle/vim-colors-solarized/
git submodule add	https://github.com/Lokaltog/vim-easymotion.git vim/bundle/vim-easymotion/
git submodule add	https://github.com/henrik/vim-indexed-search.git vim/bundle/vim-indexed-search/
git submodule add	https://github.com/gerw/vim-latex-suite.git vim/bundle/vim-latex-suite/
git submodule add	https://github.com/jcf/vim-latex.git vim/bundle/vim-latex/
git submodule add	https://github.com/kshenoy/vim-signature.git vim/bundle/vim-signature/
git submodule add	https://github.com/tpope/vim-surround.git vim/bundle/vim-surround/



dotfiles=$HOME/Dropbox/scripts/dotfiles
[ ! -e "$dotfiles" ] && echo Folder $dotfiles does not exist && exit
cd $dotfiles

echo "#############"
echo "# cd folder && git submodule init/update  ..."
echo "#############"
cd $dotfiles && cd vim/bundle/Vundle.vim/              && git submodule init && git submodule update && cd $dotfiles 
cd $dotfiles && cd vim/bundle/ctrlp.vim/               && git submodule init && git submodule update && cd $dotfiles  
cd $dotfiles && cd vim/bundle/gundo.vim/               && git submodule init && git submodule update && cd $dotfiles 
cd $dotfiles && cd vim/bundle/latex-makefile/          && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/mru/                     && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/neocomplete.vim/         && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/rainbow_parentheses.vim/ && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/supertab/                && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/switch.vim/              && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/syntastic/               && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/tagbar/                  && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/tcomment_vim/            && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-airline/             && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-colors-solarized/    && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-easymotion/          && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-indexed-search/      && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-latex-suite/         && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-latex/               && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-signature/           && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-surround/            && git submodule init && git submodule update && cd $dotfiles

# the following give (untracked content) in git status  (after manual install seem to be ok)
git submodule add	https://github.com/vim-scripts/taglist.vim.git vim/bundle/taglist.vim/
git submodule add	https://github.com/kana/vim-fakeclip vim/bundle/vim-fakeclip/
git submodule add  https://github.com/gmarik/vundle.git vim/bundle/vundle
cd $dotfiles && cd vim/bundle/taglist.vim/             && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vim-fakeclip/            && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd vim/bundle/vundle/                  && git submodule init && git submodule update && cd $dotfiles 


git add xmgrace
git add xmodmap

echo "#######################################"
echo "#git add zsh                 # has git!"
echo "#######################################"
git add zsh/ANMERKUNG_git
git add zsh/zshrc
git add zsh/zshrc.save1
git add zsh/zshrc.save2
#git add zsh/zsh-history-substring-search
#git add zsh/zsh-syntax-highlighting

echo "#######################################"
echo "#git add submodules 1"
echo "#######################################"
#git submodule add git://github.com/robbyrussell/oh-my-zsh.git zsh/oh-my-zsh/
git submodule add git://github.com/zsh-users/zsh-syntax-highlighting.git zsh/zsh-syntax-highlighting/
git submodule add git://github.com/zsh-users/zsh-history-substring-search.git zsh/zsh-history-substring-search
#cd zsh/oh-my-zsh/ && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd zsh/zsh-syntax-highlighting && git submodule init && git submodule update && cd $dotfiles
cd $dotfiles && cd zsh/zsh-history-substring-search && git submodule init && git submodule update && cd $dotfiles

echo "#######################################"
echo "#git add submodules 2"
echo "#######################################"
git submodule add https://github.com/rupa/z z
cd $dotfiles && cd z && git submodule init && git submodule update && cd $dotfiles

echo "#######################################"
echo "#git add bin" 
echo "#######################################"
git add bin
git add .gitignore
git add git_setup_dotfiles.sh
git add terminal_colors

echo "#######################################"
echo "#git commit" 
echo "#######################################"
git commit -m "first commit"
git remote add origin https://github.com/glensk/dotfiles.git
git push -u origin master
