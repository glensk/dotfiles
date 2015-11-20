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
dot=$HOME/Dropbox/scripts/dotfiles
[ "`pwd`" != "$dot" ] && echo you have to be in $dot directory and you are not && exit

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
echo "#### git pull #########################"
echo "#######################################"
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
if [ "`hostname`" = "$mylaptop" ];then
    echo
    echo
    echo
    echo
    echo
    echo "#######################################"
    echo "#### update other hosts ###############"
    echo "#######################################"
    echo "#######################################"
    echo "#### update cmpc (and cmmd) ###########"
    echo "#######################################"
    ssh.sh -t glensk@$myhost.mpie.de "cd /home/glensk/Dropbox/scripts/dotfiles;./bin/gitdown.sh;.generalrc/generalrc_alias_renew.sh"
    echo "#######################################"
    echo "#### update cmpc (and cmmd) DONE ######"
    echo "#######################################"
    echo
    echo
    echo
    echo
    echo
    #echo "#######################################"
    #echo "#### update cmmc ######################"
    #echo "#######################################"
    #ssh.sh -t glensk@$myhost.mpie.de "/home/glensk/Dropbox/scripts/dotfiles/bin/ssh.sh -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 \"cd /cmmc/u/aglen/Dropbox/scripts/dotfiles;module load git;./bin/gitdown.sh\""
    #echo
    #echo
    #echo
fi
