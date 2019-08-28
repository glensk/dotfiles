#!/bin/bash

tab-color() {
echo -ne "\033]6;1;bg;red;brightness;$1\a"
echo -ne "\033]6;1;bg;green;brightness;$2\a"
echo -ne "\033]6;1;bg;blue;brightness;$3\a"
}

tab-reset() {
echo -ne "\033]6;1;bg;*;default\a"
}

## tab-reset
## tab-color 255 0 0
## tab-color 0 255 0
## tab-color 0 0 255
## tab-reset
#
#black="0 0 0"
#red="255 0 0"
#green="0 255 0"
#blue="0 0 255"
#turquoise="64 224 208"   # seems too blue, needs to get greener
#mediumturquoise="72 209 204"
#zomp="57 167 142"
#magenta="255 0 255"
#orange="255 165 0"
##@
##@goto=""
##@
##@
##@#echo myhost:$myhost:
##@#[ "`echo $* | grep -o mac`" = "mac" ] && goto=$green
##@[ "`echo $* | grep -o cmpc`" = "cmpc" ] && goto=$red
##@[ "`echo $* | grep -o cmmd`" = "cmmd" ] && goto=$blue
##@[ "`echo $* | grep -o cmmc`" = "cmmc" ] && goto=$turquoise
##@
##@[ "`echo $* | grep -o daint`" = "daint" ] && goto=$blue
##@[ "`echo $* | grep -o cosmopc`" = "cosmopc" ] && goto=$red
##@[ "`echo $* | grep -o fidis`" = "fidis" ] && goto=$turquoise
##@[ "`echo $* | grep -o helvetios`" = "helvetios" ] && goto=$orange
##@
##@
#currenthost=""
#[ "`echo $host | grep -o cmpc`" = "cmpc" ] && currenthost=$red
#[ "$myhost" = "mac" ] && currenthost=$magenta
#[ "`echo $host | grep -o cmmd`" = "cmmd" ] && currenthost=$blue
#[ "`echo $host | grep -o cmmc`" = "cmmc" ] && currenthost=$turquoise
##[ "`echo $host | grep -o cmdft`" = "cmdft" ] && currenthost=$blue
#
#[ "`echo $host | grep -o mac`" = "mac" ] && currenthost=$magenta
#
#
##@#echo goto:$goto:
##@[ "$goto" = "" ] && tab-reset
##@[ "$goto" != "" ] && tab-color $goto
#echo "currenthost:$currenthost:  myhost:$myhost mycolor:$mycolor rgb:$mycolor_rgb"
#echo "when comming back will change to back:$mycolor_rgb"
ssh $* 

########################################################################
## revert color depends on hostname current
#echo currenthost:$currenthost:
#[ "$currenthost" = "" ] && tab-reset
#[ "$currenthost" != "" ] && tab-color $currenthost

#echo "myhost back to:$myhost color:`tab-color-host $myhost`"
#[ "$currenthost" = "" ] && tab-reset
[ "$mycolor_rgb" = "" ] && tab-reset
#[ "$currenthost" != "" ] && tab-color $currenthost
#echo "back:$mycolor_rgb:"
#echo "curr:$currenthost:"
#[ "$currenthost" != "" ] && tab-color $mycolor_rgb
[ "$mycolor_rgb" != "" ] && tab-color $mycolor_rgb
