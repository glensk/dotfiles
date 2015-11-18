#!/bin/sh

hier=`pwd`
dort=`echo $hier | sed 's|^/data/|/nas/|' | sed 's|^/home/|/nas/|' | sed 's|^/nas/|/nas/|' | sed 's|^/Users/|/nas/|'`
echo $dort
