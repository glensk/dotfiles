#!/bin/bash

if [ "$1" = "-h" ] || [ "$1" = "-help" ] || [ "$1" = "--help" ];then #|| "-help" || "--help" ]];then
    echo "OPTTIONS:"
    echo "-h, -help, --help         print this help message and exit"
    echo "-gus                      update submodules"
    echo "-oo                       only update others"
    exit
fi


    
#echo "#### cd $home/scripts #################";cd $home/scripts | echo "";\\
#echo "#### now in `pwd` #####################";echo;\\
#echo "#### cd $HOME/scripts #################";cd $HOME/scripts | echo "";\\
#echo "#### now in `pwd` #####################";echo;\\
#echor : \033[1;31m\!*\033[m
#dot=$HOME/Dropbox/scripts/dotfiles/
dot=$dotfiles

cd $dot
[ "$1" = "-aa" ] && cd $HOME/sources/aiida-alloy


echo `pwd`
one=`pwd | sed 's|^/cmmc/||'`
one=`echo $one/`

if [ "$1" != "-aa" ];then 
if [ "$one" != "$dot" ];then
    echo "pwd              : `pwd`"
    echo "but you are in   : $one" 
    echo "you have to be in: $dot"
    exit
fi
fi

if [ "$1" != "-oo" ];then
echo "#######################################"
echo "#### git add -A -v ####################"
echo "#######################################"
echo;git add -A -v;echo;

echo "#######################################"
echo "#### git commit -a -m "`date`" ########"
echo "#######################################"
echo;git commit -a -m "`date`";echo;


echo "#######################################"
echo "#### git pull (here it often wants a message)"
echo "#######################################"
# git pull is basically two actions at once: git fetch followed by a git merge
echo;git pull;echo;

echo "#######################################"
echo "#### git push #########################"
echo "#######################################"
echo;git push;



if [ "$1" = "-gus" ];then
echo "#######################################"
echo "#### to update git submodules #########"
echo "#### to run git update submodules use -gus option with this script"
echo "#######################################"
echo "git submodule update --recursive --remote"
git submodule update --recursive --remote
echo "you will now see (dong git status) ...(new commits)"
echo "git submodule update"
git submodule update
git commit -am "Updated submodules to latest"
git push
fi
fi



#echo myhost:$myhost:
#echo mylaptop:$mylaptop:
#if [ "`hostname`" = "$mylaptop" ];then
#    echo
#    echo
#    echo
#    echo
#    echo
#    echo "#######################################"
#    echo "#### update other hosts ###############"
#    echo "#######################################"
#
#
#    #@ echo "#######################################"
#    #@ echo "#### update cmpc (and cmmd) ###########"
#    #@ echo "#######################################"
#    #@ #ssh.sh -t glensk@$myhost.mpie.de "cd /home/glensk/Dropbox/scripts/dotfiles;./bin/gitdown.sh;./generalrc/generalrc_alias_renew.sh;./ssh/change_rights_of_config_to_600.sh"
#    #@ ssh.sh -t glensk@$myhost.mpie.de "cd /home/glensk/Dropbox/Albert/scripts/dotfiles;./bin/gitdown.sh;./generalrc/generalrc_alias_renew.sh;./ssh/change_rights_of_config_to_600.sh"
#    
#    #@ echo "#######################################"
#    #@ echo "#### update cmpc (and cmmd) DONE ######"
#    #@ echo "#######################################"
#    #@ echo
#    #@ echo
#    #@ echo
#    #@ echo
#    #@ echo
#    #echo "#######################################"
#    #echo "#### update cmmc ######################"
#    #echo "#######################################"
#    #ssh.sh -t glensk@$myhost.mpie.de "/home/glensk/Dropbox/Albert/scripts/dotfiles/bin/ssh.sh -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 \"cd /cmmc/u/aglen/Dropbox/Albert/scripts/dotfiles;module load git;./bin/gitdown.sh\""
#    #echo
#    #echo
#    #echo
#fi
