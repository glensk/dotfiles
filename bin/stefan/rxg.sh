#!/bin/sh

scp aglen@gate.rzg.mpg.de:~/.exchange/archive.tar.bzip2 .
tar --use-compress-program=lbzip2 -xvf archive.tar.bzip2
rm archive.tar.bzip2

