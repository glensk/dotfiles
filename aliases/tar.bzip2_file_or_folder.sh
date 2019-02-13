#!/bin/sh
check_remove=`echo $0 | sed 's|.*aliases/||' | grep "and_remove_originals"`
check_sendx=`echo $0 | sed 's|.*aliases/||' | grep "for_sendx"`
verbose="true"
[ "$check_sendx" != "" ] && verbose="false"
#echo "check_remove:$check_remove:"
#echo "check_sendx:$check_sendx:"
#exit
tar=`command -v gtar`
[ "$tar" == "" ] && tar=`command -v tar`
lbzip2=`command -v lbzip2`
if ([ "$lbzip2" = "" ] && [ ! -e "$HOME/.local/bin/lbzip2" ]);then
    echo "##################### installing lbzip2 #################"
    echo "lbzip2:$lbzip2: (if empty will be installed)"
    echo "##################### installing lbzip2 #################"
    hier=`pwd`
    [ ! -e "$HOME/sources" ] && mkdir $HOME/sources
    [ ! -e "$HOME/.local" ] && mkdir $HOME/.local
    cd $HOME/sources
    wget http://archive.lbzip2.org/lbzip2-2.5.tar.gz
    tar -xvf lbzip2-2.5.tar.gz
    cd lbzip2-2.5/
    ./configure --prefix="$HOME/.local"
    make
    make install
    cd $hier
fi
lbzip2=`command -v lbzip2`
if ([ "$lbzip2" = "" ] && [ -e "$HOME/.local/bin/lbzip2" ]);then
    echo lbzip2 $lbzip2
    echo lbzip2 exists but $HOME/.local/bin is not in your PATH
    echo make sure that $HOME/.local/bin is in your PaTH
    exit
fi
#echo ll $lbzip2
#exit
#############################################################
# tar.bzip2_file_or_folder.sh weights/ --> DOES NOT WORK
# tar.bzip2_file_or_folder.sh weights  --> DOES WORK
# --> to make both work
#############################################################


#############################################################
# DEFINE WHAT TO SEND (allpaths, can be files and or folder)
#############################################################
allpaths=$@
allpaths=`echo $@ | sed 's|/ | |g'`
[ "$allpaths" = "" ] && echo 'need $1 to the the foldername' && exit
[ "$verbose" == "true" ] && echo "allpath to tar: $allpaths"

#############################################################
# DEFINE how to name the tarred archive  
#############################################################
saveas=""
[ "`echo $* | wc -w`" = "1" ] && saveas=`echo "$1" | sed 's|/$||'`
if [ "`echo $* | wc -w`" != "1" ];then
    saveas=`echo $* | xargs -n1 | sed -e '1{h;d;}' -e 'G;s,\(.*\).*\n\1.*,\1,;h;$!d' | sed 's/\.$//g'`
    [ "$saveas" = "" ] && saveas="$1_and_other_files"
    [ "`echo $saveas | wc -c`" = "0" ] && saveas="$1_and_other_files"
fi
[ "$check_sendx" != "" ] && saveas="sendx"
[ "$verbose" == "true" ] && echo "saveas        : $saveas.tar.bzip2"
[ -e "$saveas.tar.bzip2" ] && echo $saveas does already exist && exit
################# make the tar
addcommand="";[ "$check_remove" != "" ] && addcommand=" --remove-files "


[ "$verbose" == "true" ] && echo "addcommand    : $addcommand (empty means original files are not deleted)"
$tar $addcommand --use-compress-program=lbzip2 -cvf "$saveas.tar.bzip2" $allpaths

#############################################################
# in case lbzip2 is not available or can not be compiled
#############################################################
#tar --remove-files -zcvf $saveas.tar.gz $*

if [ -e "$saveas.tar.bzip2" ];then
    [ "$verbose" == "true" ] && echo "success       : $saveas.tar.bzip2 has been created."
else
    [ "$verbose" == "true" ] && echo "ERROR  !!!!!  : $saveas.tar.bzip2 has NOT been created."
fi

###########################################################
# if lbzip2 not available
###########################################################
#cd $dotfiles/sources

