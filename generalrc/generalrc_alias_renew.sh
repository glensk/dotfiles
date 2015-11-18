#!/bin/bash

aliastcshtobash () {
		echo alias ${1}=\'`echo "${2}" | sed "s:':'\\\\\\\\'':"`\'
}

mkalias () {
    # for bash $2 need to be in ""
    eval `aliastcshtobash $1 "$2"`
}

echo $currenthost
echo $whichalias
echo $currentshell
aliasfile=$generalrc/generalrc_alias_$currenthost
source $generalrc/generalrc_alias_.sh $whichalias $currentshell 
rm -f $aliasfile
#alias | grep -v "^-" | grep -v "^../=" | sed 's|\(.*\)|alias \1|' > $aliasfile
alias | grep -v "^-" | grep -v "^../=" > $aliasfile

#[ "$1" = "" ] && echo '$1 needs to be the aliasfilepath' && exit
#echo $1
#echo ----------
#echo `eval alias`
