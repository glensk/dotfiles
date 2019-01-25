#!/bin/bash

x=`echo $1 | sed 's/\(.*\):\(.*\)/\1/' | awk '/^[0-9]+$/{print $1;exit};{print "error"}'`
y=`echo $1 | sed 's/\(.*\):\(.*\)/\2/' | awk '/^[0-9]+$/{print $1;exit};{print "error"}'`
if [ "$x" == error -o "$y" == error ]; then
  echo; echo USAGE: xmgraceblock.sh x:y files
  exit
fi
rm -fr ${HOME}/_tmp_xmgrace
inp=`echo $* | xargs -n1 | tail -n+2 | xargs`
mkdir ${HOME}/_tmp_xmgrace/
for i in $inp; do
   new=${HOME}/_tmp_xmgrace/_tmp_xmgrace_`echo $i | sed 's|/||g'`
   sed '/#/d' $i | awk 'NR==1{col=NF}; NF==col{if ('$x'=='0') print NR,$'$y'; else print $'$x',$'$y'} NF==0{print}' > $new
done
#xmgrace -geom 1000x800 -viewport .11 .27 1.08 0.98 ${HOME}/_tmp_xmgrace/_tmp_xmgrace_*
#DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -param ~/.xmgracerc ${HOME}/_tmp_xmgrace/_tmp_xmgrace_*
DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -param $dotfiles/xmgrace/tpl_log_log.par ${HOME}/_tmp_xmgrace/_tmp_xmgrace_*
rm -fr ${HOME}/_tmp_xmgrace
