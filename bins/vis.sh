#!/bin/sh
echo '$1       :' $1
echo 'which $1 :' `which $1`
echo "which vim:" `which vim`
echo "EDITOR   :" $EDITOR
echo "MYVIM    :" $MYVIM
echo "will run :" $MYVIM `which $1`
echo 
[ "$MYVIM" == "" ] && echo MYVIM is not defined && exit
#$HOME/sources/nvim/bin/nvim `which $1`
eval $MYVIM `which $1`
