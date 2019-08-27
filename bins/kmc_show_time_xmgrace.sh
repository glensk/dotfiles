#!/bin/sh

echo :$@:
[ "$1" = "" ] && echo "please provide a folder/file! (or path is not mounted)" && exit
rm -f $HOME/.tmp_xmg_command
rm -f $HOME/.tmp_xmg_to_delete*
touch tta
nr=0
echo "---------------------"
for folderin in $@;do
    echo folder: $folderin
    [ ! -e "$folderin" ] && echo $folderin does not exist && exit
done
echo "---------------------"
echo
for folderin in $@;do
    #echo f1: $folder
    folder=`echo $folderin | sed 's|KMC_AL6XXX$||'`
    folder=`echo $folderin | sed 's|simlation.out$||'`
    folder=`echo $folderin | sed 's|simulation.out$||'`
    echo folder1 $folder  # has to be the seedxxx folder
    # now remove seedxxx (just in case there were some)
    folder=`echo $folder | sed 's|/seed.*|/|g'`
    echo folder2 $folder  # has to be the seedxxx folder
    # now re introduce folder
    seedfolder=`ls -1d $folder/seed*`
    for folder in $seedfolder;do
        #if [ ! -e "$folder" 
        #echo f2: $folder
        file=$folder/KMC_AL6XXX
        if [ ! -e "$file" ];then
           a=1 
        fi
        [ ! -e "$file" ] && echo $file does not exist && exit

            # make tmp file
            file_tmp=$HOME/.tmp_xmg_to_delete$nr
            [ -e "$file_tmp" ] && rm $file_tmp
            awk '{print $1*2.4e-17,$3}' $file > $file_tmp
            ((nr++))
        #echo " -block $file_tmp -bxy 1:3"    # when block 
        #echo " -block $file_tmp -bxy 1:3" >> $HOME/.tmp_xmg_command  # when block
    done
done
#commanda=`cat $HOME/.tmp__xmg_command | xargs`

#DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -block $@ -bxy 1:3
DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch  $HOME/.tmp_xmg_to_delete*
