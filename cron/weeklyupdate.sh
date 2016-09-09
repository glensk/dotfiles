#!/bin/zsh
host=`hostname`
echo "start:$(date)" >> $HOME/Dropbox/scripts/dotfiles/cron/runtime_$host
doit() {
echo ""
echo ""
echo ""
echo ""
echo ""
echo "##############################################################################"
echo "##############################################################################"
echo "##############    $*" 
echo "##############################################################################"
echo "##############################################################################"
$* #conda update conda -y
echo "##############################################################################"
echo "##############################################################################"
echo "##############    $* done!!!!!"
echo "##############################################################################"
echo "##############################################################################"
echo ""
echo ""
echo ""
echo ""
echo ""
}

doit conda update conda -y
doit conda update anaconda -y
if [ "$host" = "mac" ];then
    doit brew update
    doit brew upgrade --all
    doit brew cleanup
    doit find -L ~/Dropbox -type l -exec unlink {} \;  # this removes dead symlinks wich make high cpuusage through opendirectoryd
    doit diskutil repairPermissions /
    #doit diskutil verifyvolume /
    #doit diskutil repairvolume /
    fi

echo "stop :$(date)" >> $HOME/Dropbox/scripts/dotfiles/cron/runtime_$host
