#!/bin/bash

aliastcshtobash () {
		echo alias ${1}=\'`echo "${2}" | sed "s:':'\\\\\\\\'':"`\'
}

mkalias () {
    # for bash $2 need to be in ""
    eval `aliastcshtobash $1 "$2"`
}


#echo $whichalias
#echo $currentshell

#[ "$onhost" = "UNKNOWN" ] && echo "onhost is UNKNOWN --> EXIT!" && exit
#[ "$onhost" = "" ] && echo "onhost empty (not found) --> EXIT!" && exit
onhost=$myhost
#echo xxx $myhost
echo "alias   : $onhost"

aliasfile=$dotfiles/generalrc/generalrc_alias_$onhost
echo "aliasfile: $aliasfile"
#source $dotfiles/generalrc_alias_.sh $whichalias $currentshell 
source $dotfiles/generalrc/generalrc_alias_.sh mkalias $currentshell 
rm -f $aliasfile
#alias | grep -v "^-" | grep -v "^../=" | sed 's|\(.*\)|alias \1|' > $aliasfile
alias | grep -v "^-" | grep -v "^../=" > $aliasfile

#[ "$1" = "" ] && echo '$1 needs to be the aliasfilepath' && exit
#echo $1
#echo ----------
#echo `eval alias`
