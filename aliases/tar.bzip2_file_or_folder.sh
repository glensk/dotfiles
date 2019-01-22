#!/bin/sh
[ "$1" = "" ] && echo 'need $1 to the the foldername' && exit
[ -e "$1.tar.bzip2" ] && echo $1.tar.bzip2 does already exist && exit
tar --use-compress-program=lbzip2 -cvf "$1.tar.bzip2" $1

###########################################################
# if lbzip2 not available
###########################################################
#cd $dotfiles/sources

