#!/bin/sh
echo dotfiles $dotfiles
cd $dotfiles
git config credential.helper store
git push `git config --get remote.origin.url`

