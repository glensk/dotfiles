#!/bin/bash

pid=`ps -A -o pid,comm | grep "/Applications.*Dropbox" | awk '{print $1}'`
[ "`echo $pid | wc -w`" != "1" ] && echo not 1
echo pid:$pid
#cpulimit -p "$pid" -z -l 7

[ "$1" = "" ] && cpulimit -p "$pid" -z -l 30 > /dev/null &   # this works jo
[ "$1" != "" ] && cpulimit -p "$pid" -z -l $1 > /dev/null &   # this works jo
#### (skript > /dev/null ) > & /dev/null &
