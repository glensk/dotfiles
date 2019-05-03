#!/bin/sh
if [ ! -e ".git" ];then
echo normal dotfile require old pw with one
echo dotfiles $dotfiles
cd $dotfiles
else
echo current folder `pwd`
fi

git config credential.helper store
git push `git config --get remote.origin.url`

echo aiida-alloy $HOME/sources/aiida-alloy
cd $HOME/sources/aiida-alloy

git config credential.helper store
git push `git config --get remote.origin.url`
