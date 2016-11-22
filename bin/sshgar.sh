#!/bin/sh

#echo myhost:$myhost
#echo hostname:$hostname
#echo hostname:`hostname`
    
##############################################################
# mac -> cmmc001 
##############################################################
if [ "`hostname`" = "mac" ];then
    #########################################################
    # from within mpie or loged in via sslvpn/forticlient
    #########################################################
    #echo host:`hostname`
    ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 -i $HOME/.ssh/garching2 aglen@cmmc001.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh /u/aglen` ] && cd `th.sh /u/aglen`;ls -F --color=auto --group-directories-first --show-control-chars --hide=\"Icon?\"; zsh"
    #########################################################
    # from anywhere in the world (no forticlient or sslvpn)
    #########################################################
    #ssh -t aglen@gate.rzg.mpg.de ssh aglen@cmmc001.bc.rzg.mpg.de

fi

##############################################################
# cmpc -> cmmc001
# oncmmc first bash is loaded
##############################################################
if [[ "`hostname`" = "$myhost" || "`hostname`" = "$myhost.mpie.de" ]];then
    ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 aglen@cmmc001.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh /u/aglen` ] && cd `th.sh /u/aglen`;ls -F --color=auto --group-directories-first --show-control-chars --hide=\"Icon?\"; zsh"
fi


if [[ "`hostname`" != "$myhost" && "`hostname`" != "$myhost.mpie.de" && "`hostname`" != "mac" ]];then
    echo host $host is neither mac nor $myhost
fi

##############################################################
# old login: from mac -> cmpc -> cmmc 
##############################################################
#ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -Y -X -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 \"cd `/Users/glensk/Dropbox/scripts/dotfiles/bin/th.sh /u/aglen`;zsh\""
