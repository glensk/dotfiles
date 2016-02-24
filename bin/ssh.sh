#!/bin/bash

tab-color() {
echo -ne "\033]6;1;bg;red;brightness;$1\a"
echo -ne "\033]6;1;bg;green;brightness;$2\a"
echo -ne "\033]6;1;bg;blue;brightness;$3\a"
}

tab-reset() {
echo -ne "\033]6;1;bg;*;default\a"
}

# tab-reset
# tab-color 255 0 0
# tab-color 0 255 0
# tab-color 0 0 255
# tab-reset

black="0 0 0"
red="255 0 0"
green="0 255 0"
blue="0 0 255"
turquoise="64 224 208"   # seems too blue, needs to get greener
mediumturquoise="72 209 204"
zomp="57 167 142"

goto=""
#echo myhost:$myhost:
#[ "`echo $* | grep -o mac`" = "mac" ] && goto=$green
[ "`echo $* | grep -o cmpc`" = "cmpc" ] && goto=$red
[ "`echo $* | grep -o cmmd`" = "cmmd" ] && goto=$blue
[ "`echo $* | grep -o cmmc`" = "cmmc" ] && goto=$turquoise


currenthost=""
[ "`echo $host | grep -o cmpc`" = "cmpc" ] && currenthost=$red
[ "`echo $host | grep -o cmmd`" = "cmmd" ] && currenthost=$blue
[ "`echo $host | grep -o cmmc`" = "cmmc" ] && currenthost=$turquoise
#[ "`echo $host | grep -o cmdft`" = "cmdft" ] && currenthost=$blue
[ "`echo $host | grep -o mac`" = "mac" ] && currenthost=$black


#echo goto:$goto:
[ "$goto" = "" ] && tab-reset
[ "$goto" != "" ] && tab-color $goto
ssh $*

########################################################################
## revert color depends on hostname current
#echo currenthost:$currenthost:
[ "$currenthost" = "" ] && tab-reset
[ "$currenthost" != "" ] && tab-color $currenthost
