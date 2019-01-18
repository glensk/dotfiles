#!/bin/sh

#echo $@
rm -f $HOME/.tmp_xmg_command
rm -f $HOME/.tmp_xmg_to_delete*
touch tta
nr=0
for folder in $@;do
    #echo f1: $folder
    folder=`echo $folder | sed 's|KMC_AL6XXX$||'`
    folder=`echo $folder | sed 's|simlation.out$||'`
    folder=`echo $folder | sed 's|simulation.out$||'`
    #echo f2: $folder
    file=$folder/KMC_AL6XXX
    [ ! -e "$file" ] && echo $file does not exist && exit
        # make tmp file
        file_tmp=$HOME/.tmp_xmg_to_delete$nr
        [ -e "$file_tmp" ] && rm $file_tmp
        awk '{print $1*2.4e-17,$3}' $file > $file_tmp
        ((nr++))
    #echo " -block $file_tmp -bxy 1:3"    # when block 
    #echo " -block $file_tmp -bxy 1:3" >> $HOME/.tmp_xmg_command  # when block
done
#commanda=`cat $HOME/.tmp__xmg_command | xargs`

#DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -block $@ -bxy 1:3
DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch  $HOME/.tmp_xmg_to_delete*
