#!/bin/sh

hier=`pwd`
dort=`echo $hier | sed 's|^/data/|/home/|' | sed 's|^/home/|/home/|' | sed 's|^/nas/|/home/|' | sed 's|^/Users/|/home/|'`

echo $dort
