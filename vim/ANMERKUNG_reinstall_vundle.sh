#!/bin/sh
echo in case deinstalled vundel
cd bundle
rm -rf vundle Vundle.vim
git clone https://github.com/gmarik/Vundle.vim.git ~/.vim/bundle/Vundle.vim
