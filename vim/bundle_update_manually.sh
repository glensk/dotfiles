#!/bin/sh

hier=$HOME/Dropbox/scripts/dotfiles/vim/bundle
cd $hier
folder=`ls -1d */`
for i in $folder;do
    cd $hier
    cd $i
    echo 
    echo
    echo $i
    git status -u
done

