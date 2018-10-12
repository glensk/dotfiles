#!/bin/sh

echo 'This script assumes that Dropbox folder is        @ $HOME/Dropbox'
echo '                          if it is somewhere else, link it to $HOME/Dropbox'
echo '                          if Dropbox does not exist, copy it there and possibly use git for that'
echo 'This script assumes that Google Drive folder is   @ $HOME/google_drive'
echo 'This script assumes that Google Drive folder is   @ $HOME/Dropbox/Albert/google_drive (on mac, to save Google drive files to Dropbox)'
echo 
host=`hostname`
onhost=`echo $host | sed 's/\..*$//'`
mylaptop="mac"

oncmmc='false'
oncmmd='false'
oncmpc='false'
oncosmopc='false'
DB='false'  # Dropbox is installed?
GD="false"  # google_drive is installed?

[ "`echo $onhost | grep -o mac`" = "mac" ]          && oncmmc="true"    && DB="true"  && GD="true" 
[ "`echo $onhost | grep -o cmmc`" = "cmmc" ]        && oncmmc="true"    && DB="false" && GD="false"
[ "`echo $onhost | grep -o cmmd`" = "cmmd" ]        && oncmmd="true"    && DB="false" && GD="false"
[ "`echo $onhost | grep -o cmpc`" = "cmpc" ]        && oncmpc="true"    && DB="false" && GD="false"
[ "`echo $onhost | grep -o cosmopc`" = "cosmopc" ]  && oncosmopc="true" && DB="true"  && GD="true"

echo "host      :$host:"
echo "onhost    :$onhost:"
echo
echo "oncmmc    :$oncmmc:" 
echo "oncmmd    :$oncmmd:"
echo "oncmpc    :$oncmpc:"
echo "oncosmopc :$oncosmopc:"
echo
echo
echo "dropboxinstalled:$dropboxinstalled:"
echo "googledriveinstalled:$googledriveinstalled:"
echo



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

echo
echo
HOMEPATH=$HOME              # can be /home/glensk or /Users/glensk
[ "$oncosmopc" = "true" ] && [ ! -e "$HOME/Dropbox" ] && ln -s /local/scratch/glensk/Dropbox $HOME/Dropbox 

# This should be kept fix since Dropbox should always be in $HOME/Dropbox, even if only link
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
echo "checking if Dropbox folder exist ... (can be created here if it not exists)"
if [ "$DB" = "false" ];then
    [ ! -e "$DROPBOX" ]              && createfolder $DROPBOX              
    [ ! -e "$DROPBOXME" ]            && createfolder $DROPBOXME            
    #[ ! -e "$DROPBOXMEGOOGLEDRIVE" ] && createfolder $DROPBOXMEGOOGLEDRIVE 
    [ ! -e "$DROPBOXMESCRIPTS" ]     && createfolder $DROPBOXMESCRIPTS
fi

echo "checking if Dropbox folder exist ... (will exit here if it not exists)"
[ ! -e "$DROPBOX" ]              && echo "no DROPBOX              folder $DROPBOX"              && exit
[ ! -e "$DROPBOXME" ]            && echo "no DROPBOXME            folder $DROPBOXME"            && exit
#[ ! -e "$DROPBOXMEGOOGLEDRIVE" ] && echo "no DROPBOXMEGOOGLEDRIVE folder $DROPBOXMEGOOGLEDRIVE" && exit
[ ! -e "$DROPBOXMESCRIPTS" ]     && echo "no DROPBOXMESCRIPTS     folder $DROPBOXMEGOOGLEDRIVE" && exit
[ ! -e "$dotfiles" ]             && echo "no dotfiles             folder $dotfiles"             && exit
echo "                                 ... successfull"
echo
cd $dotfiles
chmod u+x $dotfiles/generalrc/*
chmod u+x $dotfiles/bin/*


####################################################################
# checking dotfiles
####################################################################
echo "checking dotfiles                ..."
[ ! -e "$dotfiles" ] && echo "$dotfiles does not exist!" && exit
echo "                                 ... successfull, dotfiles from: $dotfiles"
echo


####################################################################
## Google Drive
####################################################################
#echo "checking Google_Drive            ..."
#[ "$host" = "mac" ] && [ ! -e "$GOOGLEDRIVEtarget" ] && [ -e "$DROPBOXMEGOOGLEDRIVE" ] && ln -s $DROPBOXMEGOOGLEDRIVE $GOOGLEDRIVEtarget
#[ ! -e "$GOOGLEDRIVEtarget" ] && echo $GOOGLEDRIVEtarget folder does not exist && exit
#[ -e "$GOOGLEDRIVEtarget" ] && echo "                                 ... successfully linked to (or existing) ~/google_drive"
#echo

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

####################################################################
# 
####################################################################
linkdropbox_home () {
    if [ "$DB" = "true" ];then     # dropbox path is there
        echo "trying to link $1 ..."
        if [ -e "$DROPBOXME/$1" ];then           # the true dropbox path is there
            if [ ! -e "$HOMEPATH/$1" ];then
                echo "making new link     linking $HOMEPATH/$1 ... from ... $DROPBOXME/$1" && ln -s $DROPBOXME/$1 $HOMEPATH/$1
            else
                echo "(did already exist) linking $HOMEPATH/$1 ... from ... `readlink -f $HOMEPATH/$1`"
            fi;fi;fi
}

echo "linking ..."
linkdropbox_home proj 
linkdropbox_home scripts
linkdropbox_home v
linkdropbox_home Thermodynamics


unlinktcsh() {
echo "# unlink tcsh                                                #"
    [ -h "$HOMEPATH/.tcshrc" ] && unlink $HOMEPATH/.tcshrc
    [ -h "$HOMEPATH/.tcshrc.complete" ] && unlink $HOMEPATH/.tcshrc.complete
    [ -h "$HOMEPATH/.tcshrc.alias" ] && unlink $HOMEPATH/.tcshrc.alias
    [ -h "$HOMEPATH/.tcshrc.bindkey" ] && unlink $HOMEPATH/.tcshrc.bindkey
    [ -h "$HOMEPATH/.tcshrc.set" ] && unlink $HOMEPATH/.tcshrc.set
    [ -h "$HOMEPATH/.tcshrc.local" ] && unlink $HOMEPATH/.tcshrc.local
}
unlinkbash() {
echo "# unlink bash                                               #"
    [ -h "$HOMEPATH/.bashrc" ] && unlink $HOMEPATH/.bashrc
    [ -h "$HOMEPATH/.bash_profile" ] && unlink $HOMEPATH/.bash_profile
    [ -h "$HOMEPATH/.bash_alias" ] && unlink $HOMEPATH/.bash_alias
    [ -h "$HOMEPATH/.inputrc" ] && unlink $HOMEPATH/.inputrc
}
unlinkzsh() {
echo "# unlink zsh                                                #"
    [ -h "$HOMEPATH/.zshrc" ] && unlink $HOMEPATH/.zshrc
    [ -h "$HOMEPATH/.oh-my-zsh" ] && unlink $HOMEPATH/.oh-my-zsh
    [ -h "$HOMEPATH/.bash_alias" ] && unlink $HOMEPATH/.bash_alias
}
echo
echo
echo
echo "#############################################################"
echo "# unlink tcsh,bash,zsh                                      #"
echo "#############################################################"
unlinktcsh
unlinkbash
unlinkzsh

echo "#############################################################"
echo "# link tcsh,bash,zsh                                        #"
echo "#############################################################"
file=$HOMEPATH/.tcshrc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tcsh/tcshrc $file
file=$HOMEPATH/.bashrc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/bash/bashrc $file
file=$HOMEPATH/.zshrc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/zsh/zshrc $file
echo
echo
echo

if [ "$host" = "$mylaptop" ];then 
echo "###########################################################################"
echo "#                 on mylaptop (local)                                     #"
echo "###########################################################################"
    file=$HOMEPATH/.gvimrc
    [ -h "$file" ] && unlink $file 
    [ -f "$file" ] && echo $file rm && rm -rf $file 
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/vim/startup/gvimrc $file
    
    file=$HOMEPATH/.ideavimrc
    [ -h "$file" ] && unlink $file 
    [ -f "$file" ] && echo $file rm && rm -rf $file 
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/pycharm/ideavimrc $file

fi


echo "###########################################################################"
echo "#                 LOAD EVERYWHERE                                         #"
echo "###########################################################################"
loadeverywhere () {
file=$1
from=$2
echo "LOAD EVERYWHERE $file"
echo "             from $from"
[ -h "$file" ] && unlink $file 
[ -f "$file" ] && echo $file rm && rm -rf $file 
[ ! -e "$file" ] && ln -s $from $file
[ ! -e "$file" ] && echo "$file was not linked!" && exit
if [ "$3" = "checklinkdir" ];then
	if [[ -L "$file" && -d "$file" ]];then
	echo "$file is a symlink to a directory"
	else
	echo "$file NOT a symlink to a directory (check what it is and evlt. delete it)"
	exit
	fi
fi
}

loadeverywhere $HOMEPATH/.zlogin                        $dotfiles/zsh/zlogin
loadeverywhere $HOMEPATH/.iterm2_shell_integration.zsh  $dotfiles/zsh/iterm2_shell_integration.zsh
loadeverywhere $HOMEPATH/.vim                           $dotfiles/vim/             "checklinkdir"
loadeverywhere $HOMEPATH/.ctags                         $dotfiles/vim/ctags             
loadeverywhere $HOMEPATH/.Xmodmap                       $dotfiles/xmodmap/Xmodmap 
loadeverywhere $HOMEPATH/.dir_colors                    $dotfiles/terminal_colors

file=$HOMEPATH/.ipython/profile_default/ipython_config.py
[ -h $file ] && unlink $file
[ -e $file ] && echo $file rm && rm -rf $file
[ ! -f $HOMEPATH/.ipython/profile_default ] && mkdir -p $HOMEPATH/.ipython/profile_default
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/ipython/profile_default/ipython_config.py $file

loadeverywhere $HOMEPATH/.gitignore                     $dotfiles/git/gitignore_global
loadeverywhere $HOMEPATH/.gitconfig                     $dotfiles/git/gitconfig
[ ! -e "$HOMEPATH/.subversion" ] &&  loadeverywhere $HOMEPATH/.subversion             	$dotfiles/subversion/      "checklinkdir"
loadeverywhere $HOMEPATH/.subversion/config             $dotfiles/subversion/config
loadeverywhere $HOMEPATH/.pyiron                        $dotfiles/other_dotfiles/pyiron
loadeverywhere $HOMEPATH/.vimrc                         $dotfiles/vim/vimrc
loadeverywhere $HOMEPATH/.screenrc                      $dotfiles/screen/screenrc
loadeverywhere $HOMEPATH/.Xresources                    $dotfiles/screen/Xresources
loadeverywhere $HOMEPATH/.xmgracerc                     $dotfiles/xmgrace/xmgracerc
loadeverywhere $HOMEPATH/.tmux.conf                     $dotfiles/tmux/tmux.conf

[ ! -e "$HOMEPATH/.ssh" ] && mkdir $HOMEPATH/.ssh 
loadeverywhere $HOMEPATH/.ssh/config                    $dotfiles/ssh/config 
[ ! -e "$dotfiles/ssh/known_hosts" ] && touch $dotfiles/ssh/known_hosts

if [ -e "$HOMEPATH/.kde/share/apps/konsole" ];then
    file=$HOMEPATH/.kde/share/apps/konsole
    [ -h "$file" ] && unlink $file 
    [ -d "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/konsole/ $file
    fi
if [ -e "$HOMEPATH/.kde/share/config" ];then
    file=$HOMEPATH/.kde/share/config/yakuakerc
    [ -h "$file" ] && unlink $file
    [ -f "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/yakuake/yakuakerc $file
    fi
if [ -e "$HOMEPATH/.config/terminator" ]; then
    file=$HOMEPATH/.config/terminator/config
    [ -h "$file" ] && unlink $file
    [ -f "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/terminator-solarized/config $file
    fi
chmod 600 ~/.ssh/config

echo "###########################################################################"
echo "#                         autojump                                        #"
echo "###########################################################################"
echo "autojump .... (if autojump is not working enable this in the LINK_files.sh skript"
#git clone https://github.com/wting/autojump.git
#cd $dotfiles/autojump && ./install.py && cd $dotfiles
#echo "installing autojump .... successfull"

echo "###########################################################################"
echo "#                 zsh-history-substring-search                            #"
echo "###########################################################################"
if [ ! -e "$dotfiles/zsh/zsh-history-substring-search/zsh-history-substring-search.zsh" ];then
    echo installing zsh-history-substring-search since it is not available
    cd $dotfiles/zsh
    rm -rf zsh-history-substring-search 
    git clone https://github.com/zsh-users/zsh-history-substring-search.git
fi

