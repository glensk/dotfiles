#!/bin/sh
dropboxscript=$HOME/Dropbox/scripts/bin/dropbox_.py

[ ! -e $dropboxscript ] && echo $dropboxscript does not exist && exit

excludefolder="
$HOME/Dropbox/scripts/mac_tools
$HOME/Dropbox/scripts/mathematica
$HOME/Dropbox/
"
excludefolder="
$HOME/Dropbox/Handy
"

for i in $excludefolder;do
    [ -e $i ] && echo $dropboxscript exclude add $i && $dropboxscript exclude add $i
    [ ! -e $i ] && echo $i does not exist
    echo ""
done
