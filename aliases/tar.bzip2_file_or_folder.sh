#!/bin/sh
lbzip2=`command -v lbzip2`
#exit
if ([ "$lbzip2" = "" ] && [ ! -e "$HOME/.local/bin/lbzip2" ]);then
    echo lbzip2 $lbzip2
    echo "insalling lbzip2"
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

################# all files/folders
allpaths=$@
allpaths=`echo $@ | sed 's|/ | |g'`
[ "$allpaths" = "" ] && echo 'need $1 to the the foldername' && exit
echo "allpath to tar: $allpaths"
################# name of the archive
saveas=""
[ "`echo $* | wc -w`" = "1" ] && saveas=$1
if [ "`echo $* | wc -w`" != "1" ];then
    saveas=`echo $* | xargs -n1 | sed -e '1{h;d;}' -e 'G;s,\(.*\).*\n\1.*,\1,;h;$!d' | sed 's/\.$//g'`
    [ "$saveas" = "" ] && saveas="$1_and_other_files"
    [ "`echo $saveas | wc -c`" = "0" ] && saveas="$1_and_other_files"
fi
echo "saveas        : $saveas.tar.bzip2"
[ -e "$saveas.tar.bzip2" ] && echo $saveas does already exist && exit
#exit
################# make the tar
tar --remove-files --use-compress-program=lbzip2 -cvf "$saveas.tar.bzip2" $allpaths
#############################################################
# in case lbzip2 is not available or can not be compiled
#############################################################
#tar --remove-files -zcvf $saveas.tar.gz $*

if [ -e "$saveas.tar.bzip2" ];then
    echo "success       : $saveas.tar.bzip2 has been created."
else
    echo "ERROR  !!!!!  : $saveas.tar.bzip2 has NOT been created."
fi

###########################################################
# if lbzip2 not available
###########################################################
#cd $dotfiles/sources

