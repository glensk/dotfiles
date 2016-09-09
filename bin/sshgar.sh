#!/bin/sh

echo myhost:$myhost
echo hostname:$hostname
echo hostname:`hostname`
    
##############################################################
# cmpc
##############################################################
if [ "`hostname`" = "mac" ];then
    ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 aglen@cmmc001.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh` ] && cd `th.sh`; zsh"
fi

##############################################################
# cmpc
##############################################################
if [[ "`hostname`" = "$myhost" || "`hostname`" = "$myhost.mpie.de" ]];then
    ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 aglen@cmmc001.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh` ] && cd `th.sh`; zsh"
fi


if [[ "`hostname`" != "$myhost" && "`hostname`" != "mac" ]];then
    echo host $host is neither mac nor $myhost
fi

#alias gar='sshimmer -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh` ] && cd `th.sh`; tcsh"'
#ssh -t glensk@$myhost "cd ~/Dropbox;pwd;gar;cd ~/Dropbox"
#ssh -t glensk@$myhost "gar"
#ssh -t glensk@$myhost "[ -e `th.sh` ] && cd `th.sh`;pwd; gar;tcsh"
#ssh -t glensk@$myhost "[ -e `th.sh` ] && cd `th.sh`;pwd; $gar;zsh"
#if [ "`hostname`" = "mac" ];then
#ssh.sh -t glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80"
#fi


# da jetzt das login eingerichtet wurde (das direkte) kann die if schleife weggelassen werden
#ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -Y -X -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 \"cd `/Users/glensk/Dropbox/scripts/dotfiles/bin/th.sh /u/aglen`;zsh\""
#else
    #ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh` ] && cd `th.sh`; zsh"




# ssh.sh -t glensk@$myhost.mpie.de "[ -e `th.sh` ] && cd `th.sh`; zsh" works to logon to cmpc
#
# ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 glensk@$myhost.mpie.de "[ -e `th.sh` ] && cd `th.sh`;zsh"
#
# %ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 \"zsh\""
#
#%ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 \"cd /u/aglen/Dropbox/scripts/dotfiles/bin/;zsh\""
#ssh.sh -t -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200 glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 \"cd `/Users/glensk/Dropbox/scripts/dotfiles/bin/th.sh /u/aglen`;zsh\""

