#!/bin/bash

aliastcshtobash () {
		echo alias ${1}=\'`echo "${2}" | sed "s:':'\\\\\\\\'':"`\'
}

mkalias () {
    # for bash $2 need to be in ""
    eval `aliastcshtobash $1 "$2"`
}

#echo $currenthost
#echo $whichalias
#echo $currentshell

if [ "$currenthost" = "" ];then
[ "`echo $host | grep -o cmpc`" = "cmpc" ] && currenthost=cmpc
[ "`echo $host | grep -o cmmd`" = "cmmd" ] && currenthost=cmmd
[ "`echo $host | grep -o cmmc`" = "cmmc" ] && currenthost=cmmc
fi
aliasfile=$HOME/Dropbox/scripts/dotfiles/generalrc/generalrc_alias_$currenthost
#source $HOME/Dropbox/scripts/dotfiles/generalrc_alias_.sh $whichalias $currentshell 
source $HOME/Dropbox/scripts/dotfiles/generalrc/generalrc_alias_.sh mkalias $currentshell 
rm -f $aliasfile
#alias | grep -v "^-" | grep -v "^../=" | sed 's|\(.*\)|alias \1|' > $aliasfile
alias | grep -v "^-" | grep -v "^../=" > $aliasfile

#[ "$1" = "" ] && echo '$1 needs to be the aliasfilepath' && exit
#echo $1
#echo ----------
#echo `eval alias`
