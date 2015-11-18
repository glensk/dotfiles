#!/bin/sh

hier=`pwd`
dort=`echo $hier | sed 's|^/data/|/data/|' | sed 's|^/home/|/data/|' | sed 's|^/nas/|/data/|' | sed 's|^/Users/|/data/|'`
echo $dort
