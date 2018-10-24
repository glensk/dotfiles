#!/bin/sh
[ "$1" = "" ] && echo 'need $1 to the the foldername' && exit
[ -e "$1.tar.gz" ] && echo $1.tar.gz does already exist && exit
tar -zcvf "$1".tar.gz "$1"

