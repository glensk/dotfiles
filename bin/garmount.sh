#!/bin/sh
#########################################################
## in case executed on cmpc
#########################################################
[ "`hostname | \grep -o cmpc`" = "cmpc" ] && echo dont umount nas ... this is not cecessary && echo "unmounting would be fusermount -u /nas/glensk/cmmc" && echo "now mount nas" && sshfs -o idmap=user,uid=10010,gid=10010 aglen@cmmc001.bc.rzg.mpg.de:/u/aglen/ /nas/glensk/cmmc/ -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc",transform_symlinks,BatchMode=yes

# to umount
#fusermount -u /nas/glensk/cmmc

#########################################################
## in case executed on mac
#########################################################
if [ "`hostname | \grep -o mac`" = "mac" ];then
    # umount and then mount /nas/glensk
    /Users/glensk/Dropbox/scripts/dotfiles/bin/nasmount.sh
    #echo myhost:$myhost: # works: cmpc34
    # login to cmpc and mount cmmc
    echo "########### ---->> in"
    echo "This seems currently to work!"
    #ssh glensk@$myhost "cd /nas/glensk;./garmount.sh"  # this seems currently to work and loads /nas/glensk/cmmc also from mac
    ssh glensk@$myhost "garmount.sh"  # this seems currently to work and loads /nas/glensk/cmmc also from mac
    # the next command is only necessary to map /nas/glensk/cmmc onto /u/aglen 
    sshfs glensk@`echo $myhost.mpie.de`:/nas/glensk/cmmc /u/aglen/ -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc",transform_symlinks,BatchMode=yes
    sshfs glensk@`echo $myhost.mpie.de`:/nas/glensk/cmmc /cmmc/u/aglen/ -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc",transform_symlinks,BatchMode=yes
    #ssh glensk@$myhost "fusermount -u /nas/glensk/cmmc && sshfs -o idmap=user,uid=10010,gid=10010 aglen@cmmc001.bc.rzg.mpg.de:/u/aglen/ /nas/glensk/cmmc/ -o reconnect -C -o workaround=all,Ciphers=\"blowfish-cbc\",transform_symlinks,BatchMode=yes"
    echo "########### ----->> out"
    # aternative
    # ssh -f glensk@cmpc34 -L 2222:cmmc001.bc.rzg.mpg.de:22 -N
    # sshfs -p 2222 aglen@cmmc001.bc.rzg.mpg.de:/u/aglen/ /u/aglen
    # sshfs -p 2222 aglen@mac:/u/aglen/ /u/aglen
    # sshfs -p 2222 -o idmap=user,uid=10010,gid=10010 aglen@cmmc001.bc.rzg.mpg.de:/u/aglen/ /u/aglen -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc", transform_symlinks,BatchMode=yes 
fi


#ips=`ifconfig -a | grep -o "inet 172."`
#ips_check=`echo "$ips" | wc -w | sed 's|[ ]*||g'`
#check_172=`echo "$ips" | grep "^172."`
#check_ok=`echo $check_172 | wc -w | sed 's|[ ]*||g'`
##echo was once: 172.16.13.208 @rosegger
##echo was once: 172.21.2.166  @mpie
##echo was once: 172.21.2.147  @mpie
##echo check_172:$check_172
##echo check_ok:$check_ok
##echo ips:  $ips
##echo ipsC: $ips_check
##exit
#
#    #[ -e $HOME/v/pp ] && echo "OK LINKED"
#    #[ ! -e $HOME/v/pp ] && echo "NOT LINKED"
#[ "$ips_check" = "0" ] && echo "not connected to MPIE: wrong inet adress -> cannot mount /nas/glensk ips:$ips:  ips_check:$ips_check:" && exit
#sshfs glensk@`echo $myhost`:/nas/glensk/home_cmmc /u/aglen/ -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc",transform_symlinks,BatchMode=yes
#
#volume="/u/aglen/"
#check=`mount |grep $volume`
##echo check:$check
#if mount|grep $volume; then
#    echo "\033[38m\033[1m   -> DONE -> mounted\033[0m"
#    if [ -e $HOME/v/pp ];then
#        echo "\033[38m\033[1m   -> DONE -> linked\033[0m"
#        #ln -s /nas/glensk/v ~/v   ## enough is this is done once
#    else
#        echo try again ...
#        sshfs glensk@`echo $myhost.mpie.de`:/nas/glensk/home_cmmc /u/aglen/ -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc",transform_symlinks,BatchMode=yes
#        if [ -e "/u/aglen/v/" ];then
#            echo "\033[38m\033[1m   -> DONE -> linked\033[0m"
#        else
#            garumount.sh
#            sshfs glensk@`echo $myhost.mpie.de`:/nas/glensk/home_cmmc /u/aglen/ -o reconnect -C -o workaround=all,Ciphers="blowfish-cbc",transform_symlinks,BatchMode=yes
#
#            if [ -e "/u/alglen/v" ];then
#                echo "\033[38m\033[1m   -> DONE -> linked\033[0m"
#            else
#                echo "\033[31m\033[1m DID NOT LINK PROPERLY!!!!!!!!!!!!!!!!!!!!!!!!!\033[0m"
#            fi
#        fi
#    fi
