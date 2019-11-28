#!/bin/bash

set | grep BASH_SOURCE
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
$path/fortran/dos.x $*

