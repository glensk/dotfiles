#!/bin/sh

[ "$1" = "" ] && echo 'provide $1 to be the inputfile' && exit
[ ! -f "$1" ] && echo file $1 not found && exit


$HOME/Dropbox/scripts/dotfiles/aliases/cpdf-binaries/OSX-Intel/cpdf $1 -split -chunk 2 -o $time\_out%%%.pdf

