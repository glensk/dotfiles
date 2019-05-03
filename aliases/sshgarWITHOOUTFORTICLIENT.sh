#!/bin/sh
to="cmmc002"
[ "$1" != "" ] && to="cmmc001"
echo to:$to

a=64 b=244 c=208
echo "\033]6;1;bg;red;brightness;$a\a"
echo "\033]6;1;bg;green;brightness;$b\a"
echo "\033]6;1;bg;blue;brightness;$c\a"

ssh -Y -X -t -i $HOME/.ssh/gate aglen@gate.rzg.mpg.de ssh -X -Y  aglen@$to.bc.rzg.mpg.de
