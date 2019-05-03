#!/bin/sh

mkdir -p ~/local/
cd ~/local
git clone https://github.com/vim/vim.git
cd vim/src
./configure --with-features=huge
make
