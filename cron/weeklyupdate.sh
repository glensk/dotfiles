#!/bin/zsh
echo "start:$(date)" >> /Users/glensk/Dropbox/scripts/dotfiles/cron/runtime
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
doit brew update
doit brew upgrade --all
doit brew cleanup
doit diskutil repairPermissions /
#doit diskutil verifyvolume /
#doit diskutil repairvolume /


echo "stop :$(date)" >> /Users/glensk/Dropbox/scripts/dotfiles/cron/runtime
