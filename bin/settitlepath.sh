#!/bin/bash
prin="yesno"
lmax=20
get() {
    echo $1
#sed -n '$s/.*\(...............\)$/\1/p'
}

out=`pwd | sed 's|.*/||'`;outl=`echo $out | wc -c | sed 's|[ ]*||'`; #outm=`echo $out | get`
[ "$prin" = "yes" ] && echo out:$out: #outm:$outm:
[ "$outl" -gt "$lmax" ] && echo -ne "\033]0;$out\007" && exit

###
# get only last 16 characters!!!!
###
###############################################################################################
## check if longer name to tab is ok
###############################################################################################
## here we have less $lmax characters 
out1=`pwd | sed 's|.*/\(.*\)/\(.*\)$|\1/\2|'`;out1l=`echo $out1 | wc -c | sed 's|[ ]*||'`; #out1m=`echo $out1 | get`
[ "$prin" = "yes" ] && echo out1:$out1: #outm:$out1:
[ "$out1l" -gt "$lmax" ] && echo -ne "\033]0;$out\007" && exit


## here we have less $lmax characters
out2=`pwd | sed 's|.*/\(.*\)/\(.*\)/\(.*\)$|\1/\2/\3|'`;out2l=`echo $out2 | wc -c | sed 's|[ ]*||'`; #out2m=`echo $out2 | get`
[ "$prin" = "yes" ] && echo out2:$out3: #outm:$out1:
[ "$out2l" -gt "$lmax" ] && echo -ne "\033]0;$out1\007" && exit


## here we have less $lmax characters
out3=`pwd | sed 's|.*/\(.*\)/\(.*\)/\(.*\)/\(.*\)$|\1/\2/\3/\4|'`;out3l=`echo $out3 | wc -c | sed 's|[ ]*||'`; #out3m=`echo $out3 | cut -c$cutt-`
[ "$out3l" -gt "$lmax" ] && echo -ne "\033]0;$out2\007" && exit

## if we got to this point we can just plot last out
echo -ne "\033]0;$out3\007"
