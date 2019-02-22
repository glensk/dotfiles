#!/bin/sh
if [ ! -e ".git" ];then
echo dotfiles $dotfiles
cd $dotfiles
else
echo current folder `pwd`
fi

git config credential.helper store
git push `git config --get remote.origin.url`

echo aiida-alloy $_aiida
cd $_aiida

git config credential.helper store
git push `git config --get remote.origin.url`
