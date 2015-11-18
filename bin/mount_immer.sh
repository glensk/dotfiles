#!/bin/sh
####################### daten von aussen ######################################
#base=`basename $0 | sed 's|mount.sh||'`
[ "$1" = "" ] && echo '$1 needs to be the volume (i.e. nas)' && exit
base=$1    # e.g. nas
[ "$2" = "" ] && echo '$2 needs to be the host to connect to (i.e. cmpc34)' && exit
myhost=$2  # e.g. cmpc34
volume="/$base/$USER"
#echo base in mount_immer:$base:
#echo volume:$volume:
####################### daten von aussen ######################################


echo "############ unmounting /$base/$USER ##############################################"
check=`mount |grep $volume`
umount $volume 
umount -f $volume 
info=`pgrep -lf sshfs`
#echo sshfsmounts:$info:
#echo check:$check:  # heck:glensk@cmpc02:/nas/glensk on /nas/glensk (osxfusefs, nodev, nosuid, synchronous, mounted by glensk):

if mount|grep $volume; then
    diskutil unmount force $volume
    if mount|grep $volume;then
        echo "   \033[1;4;31m->STILL MOUNTED!!!!!!!!!\033[0m"
    else
        echo "   \033[1;4;32m->DONE->unmounted\033[0m"
    fi 
else
        echo "   \033[1;4;32m->DONE->unmounted\033[0m"
fi
echo "############ unmounting /$base/$USER ##############################################"
echo 
echo 
echo "############ mounting /$base/$USER ################################################"
#sshfs $USER@$myhost:/$base/$USER /$base/$USER -o reconnect -C -o workaround=all,transform_symlinks,BatchMode=yes
echo USER:$USER:
echo myhost:$myhost:
echo base:$base:
sshfs $USER@$myhost:/$base/$USER /$base/$USER #-o reconnect -C -o workaround=all,transform_symlinks,BatchMode=yes
echo "############ mounting /nas/glensk done ###########################################"
echo 
echo
echo "############ ls /$base/$USER ######################################################"
sleep 1
ls --color /$base/$USER 1>&2
echo "############ ls /$base/$USER done #################################################"
