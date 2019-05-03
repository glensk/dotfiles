#!/bin/sh
volume="/data/$USER"
check=`mount |grep $volume`

umount $volume 
umount -f $volume 
info=`pgrep -lf sshfs`
echo sshfsmounts:
echo "$info"

#echo check:$check
if mount|grep $volume; then
    diskutil unmount force $volume
    if mount|grep $volume;then
        echo "\033[31m\033[1mSTILL MOUNTED!!!!!!!!!!!!!!!!!!!!!!\033[0m"
    else
        echo "\033[38m\033[1m   -> DONE -> unmounted\033[0m"
    fi 
else
        echo "\033[38m\033[1m   -> DONE -> unmounted\033[0m"
fi
