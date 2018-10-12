#!/bin/sh

# from http://atomsk.univ-lille1.fr/dl.php
# to compile : go to src and type make atomsk (this was successful on mac with gfortran)

# I expect to be in sources folder
sources=$HOME/Dropbox/scripts/dotfiles/sources
[ "`pwd`" != "$sources" ] && echo "not in sources folder" && exit
echo "OK, I am in sources folder :)"
file=atomsk_b0.8.3
tar -xvf $file.tar.gz
cd $file/src
make atomsk
cp atomsk $sources\_bin/atoms_$host
cd $sources\_bin
ln -s atomsk_$host atomsk
cd $sources
#rm -rf $file

#### ALSO ADD $sources/atomsk_*/ to .gitignore
