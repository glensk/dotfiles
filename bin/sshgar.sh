#!/bin/sh

#alias gar='sshimmer -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh` ] && cd `th.sh`; tcsh"'
#ssh -t glensk@$myhost "cd ~/Dropbox;pwd;gar;cd ~/Dropbox"
#ssh -t glensk@$myhost "gar"
#ssh -t glensk@$myhost "[ -e `th.sh` ] && cd `th.sh`;pwd; gar;tcsh"
#ssh -t glensk@$myhost "[ -e `th.sh` ] && cd `th.sh`;pwd; $gar;zsh"
if [ "`hostname`" = "mac" ];then
ssh.sh -t glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80"
fi
