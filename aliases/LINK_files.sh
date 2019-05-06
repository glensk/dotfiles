#!/bin/sh

echo 'Dropbox (DB) folder (required) is $HOME/Dropbox (either real or linked)'
echo 'If Dropbox does not exist: mkdir -p $HOME/Dropbox/Albert/scripts and checkout dotfiles from git'
echo 'Google Drive (GD) (not required) is $HOME/google_drive (either real or linked)'
echo 
host=`hostname`
onhost=`echo $host | sed 's/\..*$//'`
mylaptop="mac"

onmac='false'
oncmmc='false'
oncmmd='false'
ondaint='false'
oncmpc='false'
oncosmopc='false'
DB='false'  # Dropbox is installed?
GD="false"  # google_drive is installed?

[ "`echo $onhost | grep -o mac`"     = "mac"     ] && onmac="true"     && DB="true"  && GD="false" 
[ "`echo $onhost | grep -o cmmc`"    = "cmmc"    ] && oncmmc="true"    && DB="false" && GD="false"
[ "`echo $onhost | grep -o cmmd`"    = "cmmd"    ] && oncmmd="true"    && DB="false" && GD="false"
[ "`echo $onhost | grep -o cmpc`"    = "cmpc"    ] && oncmpc="true"    && DB="false" && GD="false"
[ "`echo $onhost | grep -o daint`"   = "daint"   ] && ondaint="true"   && DB="false" && GD="false"
[ "`echo $onhost | grep -o cosmopc`" = "cosmopc" ] && oncosmopc="true" && DB="true"  && GD="true"
[ "`echo $onhost | grep -o fidis`"   = "fidis"   ] && oncfidis="true"  && DB="false" && GD="false"


echo "host      :$host:"
echo "onhost    :$onhost:"
echo
echo "onmac     :$onmac:" 
echo "oncmmc    :$oncmmc:" 
echo "oncmmd    :$oncmmd:"
echo "oncmpc    :$oncmpc:"
echo "oncosmopc :$oncosmopc:"
echo
echo "DB        :$DB"
echo "GD        :$GD"
echo

[ "$oncosmopc" = "true" ] && [ ! -e "$HOME/Dropbox" ]      && ln -s /local/scratch/glensk/Dropbox $HOME/Dropbox 
[ "$oncosmopc" = "true" ] && [ ! -e "$HOME/google_drive" ] && ln -s /local/scratch/glensk/google_drive $HOME/google_drive


[ "$DB" = "true" ] && echo checking if Dropbox folder exists .... && [ ! -e "$HOME/Dropbox" ] && echo $HOME/Dropbox does not exist && exit
[ "$GD" = "true" ] && echo checking if google_drive exists folder .... && [ ! -e "$HOME/google_drive" ] && echo $HOME/google_drive does not exist && exit

createfolder () {
    echo "please create once the folder $1"
    echo
    echo "I ASSUME THAT THE DROPBOX CLIENT IS NOT INSTALLED/RUNNING ON THIS MACHINE"
    echo "can i create folder $1 ? (use y for yes or anything else for no)"
    read yn
    if [ "$yn" = "y" ];then
        echo "making directory $1";mkdir -p $1
    else
        echo "you have chosen no, therefore exit";exit
    fi
}

DROPBOX=$HOME/Dropbox
DROPBOXME=$DROPBOX/Albert
DROPBOXMESCRIPTS=$DROPBOX/Albert/scripts
#DROPBOXMEGOOGLEDRIVE=$DROPBOXME/google_drive   # it is nice to have one in the other so taht all my cosmo data are saved on dropbox automatically which they otherwise would not 
#################################################################################################
# where are the dotfiles (use from Dropbox when possible (is faster) otherwise from google_drive
#################################################################################################
dotfiles=$DROPBOXMESCRIPTS/dotfiles
[ ! -e "$dotfiles" ] && echo "dotfiles folder $dotfiles does not exist" && exit
echo "dotfiles are in $dotfiles"

####################################################################
# checking Dropbox folder
####################################################################
echo "checking if Dropbox folder exist ... (will be created here if it not exists)"
if [ "$DB" = "false" ];then
    [ ! -e "$DROPBOX" ]              && createfolder $DROPBOX              
    [ ! -e "$DROPBOXME" ]            && createfolder $DROPBOXME            
    #[ ! -e "$DROPBOXMEGOOGLEDRIVE" ] && createfolder $DROPBOXMEGOOGLEDRIVE 
    [ ! -e "$DROPBOXMESCRIPTS" ]     && createfolder $DROPBOXMESCRIPTS
fi

#echo "checking if Dropbox folder exist ... (will exit here if it not exists)"
[ ! -e "$DROPBOX" ]              && echo "no DROPBOX              folder $DROPBOX"              && exit
[ ! -e "$DROPBOXME" ]            && echo "no DROPBOXME            folder $DROPBOXME"            && exit
#[ ! -e "$DROPBOXMEGOOGLEDRIVE" ] && echo "no DROPBOXMEGOOGLEDRIVE folder $DROPBOXMEGOOGLEDRIVE" && exit
[ ! -e "$DROPBOXMESCRIPTS" ]     && echo "no DROPBOXMESCRIPTS     folder $DROPBOXMEGOOGLEDRIVE" && exit
[ ! -e "$dotfiles" ]             && echo "no dotfiles             folder $dotfiles"             && exit
echo "                                 ... successfull"
cd $dotfiles
chmod u+x $dotfiles/generalrc/*
chmod u+x $dotfiles/aliases/*
chmod u+x $dotfiles/scripts/lammps_scripts/*


####################################################################
# checking dotfiles
####################################################################
echo
echo "# checking dotfiles ... ###########################################################"
[ ! -e "$dotfiles" ] && echo "$dotfiles does not exist!" && exit
echo "                    ... successfull, dotfiles from: $dotfiles"


#####################################################################
# this was for cmmc
#####################################################################
if [ -e "/data/glensk" ];then
    # this exists at mac and cmpc
    hier=`pwd`
    cd /data/glensk
    [ -h Thermodynamics ] && unlink Thermodynamics
    [ ! -e Thermodynamics ] && [ -e Dropbox/scripts/Thermodynamics ] && ln -s Dropbox/scripts/Thermodynamics Thermodynamics
    [ ! -e Thermodynamics ] && [ -e Dropbox/Thermodynamics ] && ln -s Dropbox/Thermodynamics Thermodynamics
    [ -e Dropbox ] && echo "linking Thermodynamics (data)    ...";echo
    [ -e Dropbox ] && echo "                                 ... successfull";echo
    cd $hier
fi

echo
echo "# linking ... #####################################################################"
linkdropbox_home () {
    if [ "$DB" = "true" ];then     # dropbox path is there
        echo "trying to link $1 ..."
        if [ -e "$DROPBOXME/$1" ];then           # the true dropbox path is there
            if [ ! -e "$HOME/$1" ];then
                echo "making new link     linking $HOME/$1 ... from ... $DROPBOXME/$1" && ln -s $DROPBOXME/$1 $HOME/$1
            else
                echo "(did already exist) linking $HOME/$1 ... from ... `readlink -f $HOME/$1`"
            fi;fi;fi
}
linkdropbox_home proj 
linkdropbox_home scripts
linkdropbox_home v
linkdropbox_home Thermodynamics


echo
echo "# unlink tcsh,bash,zsh ############################################################"
unlinktcsh() {
    [ -h "$HOME/.tcshrc" ] && unlink $HOME/.tcshrc
    [ -h "$HOME/.tcshrc.complete" ] && unlink $HOME/.tcshrc.complete
    [ -h "$HOME/.tcshrc.alias" ] && unlink $HOME/.tcshrc.alias
    [ -h "$HOME/.tcshrc.bindkey" ] && unlink $HOME/.tcshrc.bindkey
    [ -h "$HOME/.tcshrc.set" ] && unlink $HOME/.tcshrc.set
    [ -h "$HOME/.tcshrc.local" ] && unlink $HOME/.tcshrc.local
}
unlinkbash() {
    [ -h "$HOME/.bashrc" ] && unlink $HOME/.bashrc
    [ -h "$HOME/.bash_profile" ] && unlink $HOME/.bash_profile
    [ -h "$HOME/.bash_alias" ] && unlink $HOME/.bash_alias
    [ -h "$HOME/.inputrc" ] && unlink $HOME/.inputrc
}
unlinkzsh() {
    [ -h "$HOME/.zshrc" ] && unlink $HOME/.zshrc
    [ -h "$HOME/.oh-my-zsh" ] && unlink $HOME/.oh-my-zsh
    [ -h "$HOME/.bash_alias" ] && unlink $HOME/.bash_alias
}
unlinktcsh
unlinkbash
unlinkzsh

echo
echo "# link tcsh,bash,zsh,generalrc ####################################################"
file=$HOME/.tcshrc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tcsh/tcshrc $file
file=$HOME/.bashrc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/bash/bashrc $file
file=$HOME/.bash_profile
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo mv $file $file.save && mv $file $file.save
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/bash/bash_profile $file
file=$HOME/.zshrc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/zsh/zshrc $file
file=$HOME/.generalrc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/generalrc/generalrc_.sh $file

if [ "$host" = "$mylaptop" ];then 
echo
echo "# on mylaptop (local) #############################################################"
    file=$HOME/.gvimrc
    [ -h "$file" ] && unlink $file 
    [ -f "$file" ] && echo $file rm && rm -rf $file 
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/vim/startup/gvimrc $file
    
    file=$HOME/.ideavimrc
    [ -h "$file" ] && unlink $file 
    [ -f "$file" ] && echo $file rm && rm -rf $file 
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/pycharm/ideavimrc $file
fi

echo
echo "# vim/bundle #################################################################"
$dotfiles/vim/vim_install_plugins.py

echo
echo "# LINK EVERYTHING #################################################################"
loadeverywhere () {
file=$1
from=$2
echo "LINK $file from $from check:$3:" | sed 's|'"$HOME"'/|~/|' | sed 's| from.*Albert/scripts/| from $|' | awk '{printf "%4s %32s %5s %53s %10s\n",$1,$2,$3,$4,$5}'
[ -h "$file" ] && unlink $file 
[ -f "$file" ] && echo $file rm && rm -rf $file 
[ ! -e "$file" ] && ln -s $from $file
[ ! -e "$file" ] && echo "$file was not linked!" && exit
if [ "$3" = "checklinkdir" ];then
	if [ -L "$file" ];then
	if [ -d "$file" ];then
	#echo "$file is a symlink to a directory"
	k=k # pass since this is ok
	else
	    echo "$file NOT a symlink to a directory (check what it is and evlt. delete it)"
	    exit
	fi
    fi
fi
}


loadeverywhere $HOME/.zlogin                        $dotfiles/zsh/zlogin
loadeverywhere $HOME/.iterm2_shell_integration.zsh  $dotfiles/zsh/iterm2_shell_integration.zsh
loadeverywhere $HOME/.vim                           $dotfiles/vim/             "checklinkdir"
loadeverywhere $HOME/.ctags                         $dotfiles/vim/ctags             
loadeverywhere $HOME/.Xmodmap                       $dotfiles/xmodmap/Xmodmap 
loadeverywhere $HOME/.dir_colors                    $dotfiles/terminal_colors

file=$HOME/.ipython/profile_default/ipython_config.py
[ -h $file ] && unlink $file
[ -e $file ] && echo $file rm && rm -rf $file
[ ! -f $HOME/.ipython/profile_default ] && mkdir -p $HOME/.ipython/profile_default
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/ipython/profile_default/ipython_config.py $file

loadeverywhere $HOME/.gitignore                     $dotfiles/git/gitignore_global
loadeverywhere $HOME/.gitconfig                     $dotfiles/git/gitconfig
[ ! -e "$HOME/.subversion" ] && mkdir $HOME/.subversion 
loadeverywhere $HOME/.subversion/config             $dotfiles/subversion/config
loadeverywhere $HOME/.pyiron                        $dotfiles/other_dotfiles/pyiron
loadeverywhere $HOME/.vimrc                         $dotfiles/vim/vimrc
loadeverywhere $HOME/.screenrc                      $dotfiles/screen/screenrc
loadeverywhere $HOME/.Xresources                    $dotfiles/screen/Xresources
loadeverywhere $HOME/.xmgracerc                     $dotfiles/xmgrace/xmgracerc
loadeverywhere $HOME/.tmux.conf                     $dotfiles/tmux/tmux.conf

[ ! -e "$HOME/.ssh" ] && mkdir $HOME/.ssh 
loadeverywhere $HOME/.ssh/config                    $dotfiles/ssh/config 
[ ! -e "$dotfiles/ssh/known_hosts" ] && touch $dotfiles/ssh/known_hosts

if [ -e "$HOME/.kde/share/apps/konsole" ];then
    file=$HOME/.kde/share/apps/konsole
    [ -h "$file" ] && unlink $file 
    [ -d "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/konsole/ $file
    fi
if [ -e "$HOME/.kde/share/config" ];then
    file=$HOME/.kde/share/config/yakuakerc
    [ -h "$file" ] && unlink $file
    [ -f "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/yakuake/yakuakerc $file
    fi
if [ -e "$HOME/.config/terminator" ]; then
    file=$HOME/.config/terminator/config
    [ -h "$file" ] && unlink $file
    [ -f "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/terminator-solarized/config $file
    fi
chmod 600 ~/.ssh/config

echo
echo "# autojump ########################################################################"
echo "# autojump .... (if autojump is not working enable this in the LINK_files.sh skript"
if [ ! -e "autojump" ];then
    echo "installing autojump"
    cd $dotfiles
    git clone https://github.com/wting/autojump.git
    cd $dotfiles/autojump 
    ./install.py 
    cd $dotfiles
fi

echo
echo "# zsh-history-substring-search ####################################################"
if [ ! -e "$dotfiles/zsh/zsh-history-substring-search/zsh-history-substring-search.zsh" ];then
    echo installing zsh-history-substring-search since it is not available
    cd $dotfiles/zsh
    rm -rf zsh-history-substring-search 
    git clone https://github.com/zsh-users/zsh-history-substring-search.git
fi

echo
echo "# zsh-syntax-highlighting #########################################################"
if [ ! -e "$dotfiles/zsh/zsh-syntax-highlighting/zsh-syntax-highlighting.zsh" ];then
    echo installing zsh-syntax-highlighting since it is not available
    cd $dotfiles/zsh
    rm -rf zsh-syntax-highlighting
    git clone https://github.com/zsh-users/zsh-syntax-highlighting.git
fi

echo
echo "# make git not upload *.pyc files #########################################################"
cd $dotfiles
git rm --cached *.pyc
find . -name '*.pyc' | xargs -n 1 git rm --cached

# set rights for id_rsa
chmod 400 $HOME/.ssh/id_rsa

echo
echo "set up passwords to be able to sx and rx without password"
echo "### pw is normal (normal zwei ohne sterne 1234)"
ssh-copy-id -i $HOME/.ssh/id_rsa.pub aglensk@ela.cscs.ch

