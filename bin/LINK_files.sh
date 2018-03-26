#!/bin/sh

host=`hostname`
mylaptop="mac"
cd $HOME/Dropbox/Albert/scripts/dotfiles

linktcsh="false"
linkbash="false"
linkbash="true"
linkbash="false"
linkzsh="false"
oncmmc='false'
oncmmd='false'
oncmpc='false'

onhost=`echo $host | sed 's/\..*$//'`
echo host  :$host:
echo onhost:$onhost:

[ "`echo $onhost | grep -o cmmc`" = "cmmc" ] && oncmmc="true"
[ "`echo $onhost | grep -o cmmd`" = "cmmd" ] && oncmmd="true"
[ "`echo $onhost | grep -o cmpc`" = "cmpc" ] && oncmpc="true"

[ "$host" = "$mylaptop" ] && linkzsh="true"
[ "$host" = "$mylaptop" ] && linktcsh="true"
[ "$host" = "$mylaptop" ] && linkbash="true"

[ "$oncmmc" = "true" ]    && linkzsh="true"
#[ "$oncmmc" = "true" ]    && linktcsh="true"

#[ "$oncmmd" = "true" ]    && linktcsh="true"
[ "$oncmmd" = "true" ]    && linkzsh="true"

[ "$oncmpc" = "true" ]    && linkzsh="true"

echo
echo oncmmc:$oncmmc:
echo oncmmd:$oncmmd:
echo oncmpc:$oncmpc:
echo
echo linktcsh :$linktcsh:
echo linkbash :$linkbash:
echo linkzsh  :$linkzsh:
echo
[ "$linkcsh" = "false" ] && [ "$linkcsh" = "false" ] && [ "$linkcsh" = "false" ] && echo "all, {tcsh,bash,zsh}link are false, WHICH SYSTEM?" && exit

####################################################################
## checking if Dropbox folder exist
####################################################################
echo "checking if Dropbox folder exist ..."
HOMEPATH=$HOME   # can be /home/glensk or /Users/glensk
DB=$HOME/Dropbox
ME=$DB/Albert
DATAPATH=/data/$USER
dotfiles=$ME/scripts/dotfiles
[ ! -e "$ME/scripts" ] && echo "$ME/scripts does not exist!" && exit
[ ! -e "$dotfiles" ] && echo "$dotfiles does not exist!" && exit
echo "                                 ... successfull"
# @cmmd: folder exist, created by hand
# @cmpc: folder is linked to /data/glensk/Dropbox by hand
# @ mac: folder exists, created by dropbox.app

####################################################################
## Dropbox
####################################################################
echo linking Dropbox ...;
[ ! -e "$DB" ] && [ "$oncmpc" = "true" ] && [ -e "/data/glensk/Dropbox" ] && ln -s /data/glensk/Dropbox ~/Dropbox 
[ ! -e "$DB" ] && echo $DB folder does not exist && exit
[ -e "$DB" ] && echo "                                 ... successfull"

####################################################################
# scripts (this i want to link to Dropbox/Albert/scripts not to $HOMEPATH/Dropbox/Albert/scripts)
####################################################################
echo "linking scripts ..."
hier=`pwd`
cd $HOMEPATH
[ -h scripts ] && unlink scripts
[ -e scripts ] && echo "$HOMEPATH/scripts seems to exist and not be a link!" && exit
ln -s $ME/scripts scripts
[ -h scripts ] && echo "                                 ... successfull"
cd $hier
   
####################################################################
# Thermodynamics
####################################################################
echo "linking Thermodynamics ..."
hier=`pwd`
cd $HOMEPATH
[ -h Thermodynamics ] && unlink Thermodynamics
[ -e Thermodynamics ] && echo "$HOMEPATH/Thermodynamics seems to exist and not be a link!" && exit
ln -s $ME/Thermodynamics Thermodynamics
[ -h Thermodynamics ] && echo "                                 ... successfull"
cd $hier



#@ echo HOMEPATH:$HOMEPATH
#@     if [ -e "$HOMEPATH/Thermodynamics" ];then
#@         echo "$HOMEPATH/Thermodynamics DOES already exist!!! can not be linked!!!";echo
#@     else
#@         if [ ! -e "$HOMEPATH/Albert/Dropbox/Thermodynamics" ];then
#@             echo $HOMEPATH/Albert/Dropbox/Thermodynamics not found
#@         else
#@             # here we know the we dont have the Thermodynamics folder linked but it is available
#@             cd $HOMEPATH   # /home/glensk @ cmpc
#@             ln -s $HOMEPATH/Albert/Dropbox/Thermodynamics $HOMEPATH/Thermodynamics
#@             [ ! -e "$HOMEPATH/Thermodynamics" ] && echo "!!!!!!!!!!!!!!! Thermodynamics could not be linked1 !!!!!!!!!!!!!!"
#@             [ -e "$HOMEPATH/Thermodynamics" ] && echo "Thermodynamics was correctly linked"
#@             #[ -e Dropbox/scripts/Thermodynamics ] && [ ! -e Thermodynamics ] && unlink Thermodynamics && ln -s Dropbox/scripts/Thermodynamics Thermodynamics
#@             #[ -e Dropbox/Thermodynamics ] && [ ! -e Thermodynamics ] && unlink Thermodynamics && ln -s Dropbox/Thermodynamics Thermodynamics
#@             #[ ! -e Thermodynamics ] && echo "!!!!!!!!!!!!!!! Thermodynamics could not be linked1 !!!!!!!!!!!!!!"
#@             cd $hier
#@         
#@         fi
#@     fi
#@ #fi
if [ -e "/data/glensk" ];then
    echo "linking Thermodynamics (data)...";echo
    # this exists at mac and cmpc
    hier=`pwd`
    cd /data/glensk
    [ -h Thermodynamics ] && unlink Thermodynamics
    [ ! -e Thermodynamics ] && [ -e Dropbox/scripts/Thermodynamics ] && ln -s Dropbox/scripts/Thermodynamics Thermodynamics
    [ ! -e Thermodynamics ] && [ -e Dropbox/Thermodynamics ] && ln -s Dropbox/Thermodynamics Thermodynamics
    cd $hier
fi

####################################################################
# 
####################################################################
echo "linking diverses ..."
[ ! -e "$ME/proj" ] && [ -e "$ME/proj" ] && echo "linking ~/proj" && ln -s $ME/proj $HOMEPATH/proj

[ ! -e "$HOMEPATH/Documents" ] && [ -e "$ME/Documents" ] && echo "~/linking ~/Documents" && ln -s $ME/Documents $HOMEPATH/Documents

[ -e "$HOMEPATH/scripts/dotfiles/commands" ] && [ ! -e "$HOMEPATH/scripts/commands" ] && echo "linking ~/scripts/commands" && ln -s $dotfiles/commands $HOMEPATH/scripts/commands








unlinktcsh() {
echo "#############################################################"
echo "# unlink tcsh                                                #"
echo "#############################################################"
    [ -h "$HOMEPATH/.tcshrc" ] && unlink $HOMEPATH/.tcshrc
    [ -h "$HOMEPATH/.tcshrc.complete" ] && unlink $HOMEPATH/.tcshrc.complete
    [ -h "$HOMEPATH/.tcshrc.alias" ] && unlink $HOMEPATH/.tcshrc.alias
    [ -h "$HOMEPATH/.tcshrc.bindkey" ] && unlink $HOMEPATH/.tcshrc.bindkey
    [ -h "$HOMEPATH/.tcshrc.set" ] && unlink $HOMEPATH/.tcshrc.set
    [ -h "$HOMEPATH/.tcshrc.local" ] && unlink $HOMEPATH/.tcshrc.local
}
unlinkbash() {
echo "#############################################################"
echo "# unlink bash                                               #"
echo "#############################################################"
    [ -h "$HOMEPATH/.bashrc" ] && unlink $HOMEPATH/.bashrc
    [ -h "$HOMEPATH/.bash_profile" ] && unlink $HOMEPATH/.bash_profile
    [ -h "$HOMEPATH/.bash_alias" ] && unlink $HOMEPATH/.bash_alias
    [ -h "$HOMEPATH/.inputrc" ] && unlink $HOMEPATH/.inputrc
}
unlinkzsh() {
echo "#############################################################"
echo "# unlink zsh                                                #"
echo "#############################################################"
    [ -h "$HOMEPATH/.zshrc" ] && unlink $HOMEPATH/.zshrc
    [ -h "$HOMEPATH/.oh-my-zsh" ] && unlink $HOMEPATH/.oh-my-zsh
    [ -h "$HOMEPATH/.bash_alias" ] && unlink $HOMEPATH/.bash_alias
}
echo
echo
echo
unlinktcsh
unlinkbash
unlinkzsh
echo
echo
echo

if [ "$linktcsh" = "true" ];then
echo "#############################################################"
echo "# link tcsh                                                 #"
echo "#############################################################"
    file=$HOMEPATH/.tcshrc
    [ -h "$file" ] && unlink $file
    [ -f "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tcsh/tcshrc $file

    #file=$HOMEPATH/.tcshrc.complete
    #[ -h $file ] && unlink $file
    #[ -e $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tcsh/tcshrc.complete $file
    
    #file=$HOMEPATH/.tcshrc.alias
    #[ -h $file ] && unlink $file
    #[ -e $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tcsh/tcshrc.alias $file
    
    #file=$HOMEPATH/.tcshrc.bindkey
    #[ -h $file ] && unlink $file
    #[ -e $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tcsh/tcshrc.bindkey $file
    
    #file=$HOMEPATH/.tcshrc.set
    #[ -h $file ] && unlink $file
    #[ -e $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tcsh/tcshrc.set $file
    
    #file=$HOMEPATH/.tcshrc.local
    #[ -h $file ] && unlink $file
    #[ -e $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tcsh/tcshrc.local $file
    
    #file=$HOMEPATH/.git-completion.tcsh    # seems unused in tcsh
    #[ -h $file ] && unlink $file
    #[ -e $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/git/git-completion.tcsh $file
fi


if [ "$linkbash" = "true" ];then
echo "#############################################################"
echo "# link bash                                                 #"
echo "#############################################################"
    file=$HOMEPATH/.bashrc
    [ -h "$file" ] && unlink $file
    [ -f "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/bash/bashrc $file
    
    #file=$HOMEPATH/.bash_profile
    #[ -h "$file" ] && unlink $file
    #[ -f "$file" ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/bash/bash_profile $file
    #
    #file=$HOMEPATH/.bash_alias
    #[ -h "$file" ] && unlink $file
    #[ -f "$file" ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/bash/bash_alias $file

    #file=$HOMEPATH/.git-completion.bash   # seems unused in bash
    #[ -h $file ] && unlink $file
    #[ -e $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/git/git-completion.bash $file
    
    #file=$HOMEPATH/.git-completion.sh      # seems unused in bash
    #[ -h $file ] && unlink $file
    #[ -e $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/git/git-completion.bash $file

    #file=$HOMEPATH/.inputrc      #only used in bashrc
    #[ -h "$file" ] && unlink $file 
    #[ -f "$file" ] && echo $file rm && rm -rf $file 
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/bash/inputrc $file
fi

if [ "$linkzsh" = "true" ];then
echo "#############################################################"
echo "# link zsh                                                  #"
echo "#############################################################"
    file=$HOMEPATH/.zshrc
    [ -h "$file" ] && unlink $file
    [ -f "$file" ] && echo $file rm && rm -rf $file
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/zsh/zshrc $file

    #file=$HOMEPATH/.oh-my-zsh
    #[ -h $file ] && unlink $file
    #[ -d $file ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/zsh/oh-my-zsh $file

    #file=$HOMEPATH/.bash_alias    # here are also the zsh aliases stored
    #[ -h "$file" ] && unlink $file
    #[ -f "$file" ] && echo $file rm && rm -rf $file
    #[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/bash/bash_alias $file
fi

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
echo "#                 load everywhere                                         #"
echo "###########################################################################"
file=$HOMEPATH/.zlogin  # can do this since I use zsh everywhere
[ -h "$file" ] && unlink $file 
[ -f "$file" ] && echo $file rm && rm -rf $file 
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/zsh/zlogin $file

file=$HOMEPATH/.iterm2_shell_integration.zsh  # can do this since I use zsh everywhere
[ -h "$file" ] && unlink $file 
[ -f "$file" ] && echo $file rm && rm -rf $file 
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/zsh/iterm2_shell_integration.zsh $file

file=$HOMEPATH/.vim
[ -h $file ] && unlink $file
[ -d $file ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/vim/ $file


file=$HOMEPATH/.ctags
[ -h $file ] && unlink $file
[ -d $file ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/vim/ctags $file

file=$HOMEPATH/.Xmodmap
[ -h $file ] && unlink $file
[ -e $file ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/xmodmap/Xmodmap $file

file=$HOMEPATH/.dir_colors
[ -h $file ] && unlink $file
[ -e $file ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/terminal_colors $file
echo 'nnnnnnnn'
file=$HOMEPATH/.ipython/profile_default/ipython_config.py
[ -h $file ] && unlink $file
[ -e $file ] && echo $file rm && rm -rf $file
[ ! -f $HOMEPATH/.ipython/profile_default ] && mkdir -p $HOMEPATH/.ipython/profile_default
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/ipython/profile_default/ipython_config.py $file

#file=$HOMEPATH/.autojump
#[ -h $file ] && unlink $file
#[ -e $file ] && echo $file rm && rm -rf $file
#[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/autojump $file

file=$HOMEPATH/.gitignore
[ -h $file ] && unlink $file
[ -e $file ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/git/gitignore_global $file

file=$HOMEPATH/.gitconfig
[ -h $file ] && unlink $file
[ -e $file ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/git/gitconfig $file

file=$HOMEPATH/.subversion/config
[ -h $file ] && unlink $file
[ -e $file ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/subversion/config $file

file=$HOMEPATH/.vimrc
[ -h "$file" ] && unlink $file 
[ -f "$file" ] && echo $file rm && rm -rf $file 
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/vim/vimrc $file

file=$HOMEPATH/.screenrc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/screen/screenrc $file

file=$HOMEPATH/.Xresources
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/screen/Xresources $file

file=$HOMEPATH/.xmgracerc
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/xmgrace/xmgracerc $file

file=$HOMEPATH/.tmux.conf
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/tmux/tmux.conf $file

[ ! -e "$HOMEPATH/.ssh" ] && mkdir $HOMEPATH/.ssh 
file=$HOMEPATH/.ssh/config
[ -h "$file" ] && unlink $file
[ -f "$file" ] && echo $file rm && rm -rf $file
[ ! -e "$file" ] && echo $file link && ln -s $dotfiles/ssh/config $file


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
    [ ! -e "$file" ] && echo $file link && ln -s $dotfiles/terminator/config $file
    fi
chmod 600 ~/.ssh/config
cd $dotfiles/autojump && ./install.py && cd $dotfiles
