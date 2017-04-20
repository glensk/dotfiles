#!/bin/sh

#echo myhost:$myhost
#echo hostname:$hostname
#echo hostname:`hostname`
to="cmmc002"
[ "$1" != "" ] && to="cmmc001"
echo to:$to
##############################################################
# mac -> cmmc001 
##############################################################
if [ "`hostname`" = "mac" ];then
    #########################################################
    # from within mpie or loged in via sslvpn/forticlient
    #########################################################
    #echo host:`hostname`
    # with ls
    #ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 -i $HOME/.ssh/garching2 aglen@$to.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh /u/aglen` ] && cd `th.sh /u/aglen`;ls -F --color=auto --group-directories-first --show-control-chars --hide=\"Icon?\"; zsh"
    # whithout ls
    ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 -i $HOME/.ssh/garching2 aglen@$to.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh /u/aglen` ] && cd `th.sh /u/aglen`; zsh"
    #########################################################
    # from anywhere in the world (no forticlient or sslvpn)
    #########################################################
    #ssh -t aglen@gate.rzg.mpg.de ssh aglen@$to.bc.rzg.mpg.de

fi

##############################################################
# cmpc -> cmmc001
# oncmmc first bash is loaded
##############################################################
if [[ "`hostname`" = "$myhost" || "`hostname`" = "$myhost.mpie.de" ]];then
    # with ls
    #ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 aglen@$to.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh /u/aglen` ] && cd `th.sh /u/aglen`;ls -F --color=auto --group-directories-first --show-control-chars --hide=\"Icon?\"; zsh"
    # without ls
    ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 aglen@$to.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh /u/aglen` ] && cd `th.sh /u/aglen`; zsh"
fi


if [[ "`hostname`" != "$myhost" && "`hostname`" != "$myhost.mpie.de" && "`hostname`" != "mac" ]];then
    echo host $host is neither mac nor $myhost
fi

##############################################################
# old login: from mac -> cmpc -> cmmc 
##############################################################
#ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -Y -X -t aglen@$to.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 \"cd `/Users/glensk/Dropbox/scripts/dotfiles/bin/th.sh /u/aglen`;zsh\""
