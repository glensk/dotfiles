#!/bin/sh
[ "$1" = "" ] && echo 'need $1 to the the foldername' && exit
[ -e "$1.zpaq" ] && echo $1.zpaq does already exist && exit
zpaq a "$1".zpaq "$1" -m5

