#!/bin/sh

host=`hostname`
mylaptop="mac"
cd $HOME/Dropbox/scripts/dotfiles

linktcsh="false"
linkbash="false"
linkzsh="false"
oncmmc='false'
oncmmd='false'
oncmpc='false'

onhost=`echo $host | sed 's/\..*$//'`
echo host  :$host:
echo onhost:$onhost:
echo

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

#file=/home/$USER/v
#[ -e $file ] && echo $file unlink && unlink $file
#[ ! -e "$file" ] && echo $file link && ln -s /nas/glensk/v /home/$USER
#echo ""
####################################################################
##  Dropbox
####################################################################
#file=/home/$USER/doc 
#[ -e $file ] && echo $file unlink && unlink $file
#[ ! -e "$file" ] && echo $file link && ln -s /home/$USER/Dropbox/doc /home/$USER/doc
#echo ""


########################## dotfiles ################################
## .vim
## .vimrc
## .tcshrc
## .screenrc
## .xmgracerc
## .ctags
## konsole
########################## dotfiles ################################



####################################################################
## Dropbox
####################################################################
echo;echo linking Dropbox ...;
HOMEPATH=$HOME   # can be /home/glensk or /Users/glensk
DATAPATH=/data/$USER
# @cmmd: folder exist, created by hand
# @cmpc: folder is linked to /data/glensk/Dropbox by hand
# @ mac: folder exists, created by dropbox.app

[ ! -e "$HOMEPATH/Dropbox" ] && [ "$oncmpc" = "true" ] && [ -e "/data/glensk/Dropbox" ] && ln -s /data/glensk/Dropbox ~/Dropbox 
[ ! -e "$HOMEPATH/Dropbox" ] && echo $HOMEPATH/Drobpox folder does not exist && exit




####################################################################
# scripts (this i want to link to Dropbox/scripts not to $HOMEPATH/Dropbox/scripts)
####################################################################
if [ -e "$HOMEPATH/scripts" ];then
    # $HOMEPATH @ mac:  /Users/glensk
    # $HOMEPATH @ cmpc: /home/glensk
    # $HOMEPATH @ cmmd: /home/glensk
    # scripts exist (in every of the above mentioned homefolders the scripts folder should exist)
     echo;echo linking scripts ...;echo
     hier=`pwd`
     cd $HOMEPATH
     unlink scripts 
     ln -s Dropbox/scripts scripts
     cd $hier
else
    # scripts does not exist
    if [ -e "$HOMEPATH/Dropbox/scripts" ];then
        echo;echo linking scripts ...;echo 
        hier=`pwd`
        cd $HOMEPATH
        #unlink scripts 
        ln -s Dropbox/scripts scripts
        cd $hier
    else
        echo;echo "!!!!!!!!!!!!!!! scripts could not be linked !!!!!!!!!!!!!!";echo
    fi
fi
   
####################################################################
# Thermodynamics
####################################################################
echo HOMEPATH:$HOMEPATH
#if [ -e "$HOMEPATH/Dropbox/scripts/Thermodynamics" ];then
#    echo "linking Thermodynamics (home)...";echo
#    [ -h $HOMEPATH/Thermodynamics ] && echo unlink $HOMEPATH/Thermodynamics
#    hier=`pwd`
#    cd $HOMEPATH   # /home/glensk @ cmpc
#    [ -e Dropbox/scripts/Thermodynamics ] && [ ! -e Thermodynamics ] && ln -s Dropbox/scripts/Thermodynamics Thermodynamics
#    [ -e Dropbox/Thermodynamics ] && [ ! -e Thermodynamics ] && ln -s Dropbox/Thermodynamics Thermodynamics
#    [ ! -e Thermodynamics ] && echo "!!!!!!!!!!!!!!! Thermodynamics could not be linked1 !!!!!!!!!!!!!!"
#    cd $hier
#else
    if [ -e "$HOMEPATH/Thermodynamics" ];then
        echo "$HOMEPATH/Thermodynamics DOES already exist!!! can not be linked!!!";echo
    else
        if [ ! -e "$HOMEPATH/Dropbox/Thermodynamics" ];then
            echo $HOMEPATH/Dropbox/Thermodynamics not found
        else
            # here we know the we dont have the Thermodynamics folder linked but it is available
            cd $HOMEPATH   # /home/glensk @ cmpc
            ln -s $HOMEPATH/Dropbox/Thermodynamics $HOMEPATH/Thermodynamics
            [ ! -e "$HOMEPATH/Thermodynamics" ] && echo "!!!!!!!!!!!!!!! Thermodynamics could not be linked1 !!!!!!!!!!!!!!"
            [ -e "$HOMEPATH/Thermodynamics" ] && echo "Thermodynamics was correctly linked"
            #[ -e Dropbox/scripts/Thermodynamics ] && [ ! -e Thermodynamics ] && unlink Thermodynamics && ln -s Dropbox/scripts/Thermodynamics Thermodynamics
            #[ -e Dropbox/Thermodynamics ] && [ ! -e Thermodynamics ] && unlink Thermodynamics && ln -s Dropbox/Thermodynamics Thermodynamics
            #[ ! -e Thermodynamics ] && echo "!!!!!!!!!!!!!!! Thermodynamics could not be linked1 !!!!!!!!!!!!!!"
            cd $hier
        
        fi
    fi
#fi
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

[ ! -e "$HOMEPATH/proj" ] && [ -e "$HOMEPATH/Dropbox/proj" ] && echo "linking ~/proj" && ln -s $HOMEPATH/Dropbox/proj $HOMEPATH/proj

[ ! -e "$HOMEPATH/Documents" ] && [ -e "$HOMEPATH/Dropbox/Documents" ] && echo "~/linking ~/Documents" && ln -s $HOMEPATH/Dropbox/Documents $HOMEPATH/Documents

[ -e "$HOMEPATH/scripts/dotfiles/commands" ] && [ ! -e "$HOMEPATH/scripts/commands" ] && echo "linking ~/scripts/commands" && ln -s $HOMEPATH/scripts/dotfiles/commands $HOMEPATH/scripts/commands





dotfiles=$HOMEPATH/Dropbox/scripts/dotfiles



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
fi


echo "###########################################################################"
echo "#                 load everywhere                                         #"
echo "###########################################################################"
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
