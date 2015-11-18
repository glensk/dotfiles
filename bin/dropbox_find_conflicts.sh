#!/bin/sh
#find ~/Dropbox | grep onflic
find $HOME/Dropbox/ -name "*onflic*" ! -name "dropbox_find_conflicts.sh"

dd=`date +%s`
cp $0 $HOME/$dd\_dr
#[ "$1" = "-rf" ] && i
find $HOME/Dropbox/ -name "*onflic*" ! -name "dropbox_find_conflicts.sh" -exec rm {} \;
