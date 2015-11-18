#!/bin/bash
commitmessage=$1
[ "$1" = "" ] && [ "`pwd | sed 's|.*/||'`" = "scripts" ] && commitmessage=`date`

commit () {
    [ "$1" = "" ] && echo no commit messge && exit
    svn commit -m "$1"
}

excludefolder="dotfiles/vim/view ka/kafolder"

################################################################################################################################
# svn status
################################################################################################################################
status1=`svn status -u | sed 's|^Status against revision:.*||'`
status=`echo "$status1" | grep -v "^.[ \t]*\..*" | sed 's|^[ ]*||g'`
echo -e "\033[1;4;32m------ svn status --------------------------------------------------------\033[0m"
echo "$status"
echo -e "\033[1;4;32m------ svn status done ---------------------------------------------------\033[0m"
[ "$status" = "" ] && echo && echo -e "\033[1;4;32m------- NOTHING CHANGED :) DONE ------------\033[0m"
echo



################################################################################################################################
# look at changed files
################################################################################################################################
echo -e "\033[1;4;32m------ look at changed files ---------------------------------------------\033[0m"
# conflicts
conflicts=`echo "$status" | grep -v "^[ ]" | grep "^.[ ]*\*" | awk '{print $NF}'`
conflicts_nr=`echo "$conflicts" | wc -w | sed 's|^[ ]*||'`
    ################################################################################################################################
    # look at conflicted files
    ################################################################################################################################
    theirs_full=""
    for i in $conflicts;do
        for j in $excludefolder;do
        echo "conflicting file: $i"
        if [ "`echo $i | grep "^$j"`" = "$j" ];then
            addex="yes"
        else
            addex="no"
        fi

        echo "conflicting file: $i      |||| addex: $addex"
        [ "$addex" = "yes" ] && theirs_full="$theirs_full $i" && conflicts_nr=`echo $conflicts_nr | awk '{print $1-1}'`
        # svn update dotfiles/vim/view/~=+.tcshrc.alias= --accept theirs-full
    done
    done
    if [ "$theirs_full" != "" ];then
        echo -e "\033[1;4;32m----UPDATING conflicting files -------\033[0m"
        echo --------theirs_full:-------
        echo $theirs_full
        echo --------theirs_full:-------
        svn update $theirs_full --accept theirs-full
        echo -e "\033[1;4;32m----UPDATING conflicting files DONE---\033[0m"
    fi

[ "$conflicts_nr" != "0" ] && echo "conflicts:$conflicts_nr:  files: $conflicts" && exit

allelines=5
ex="no"
# all entries
all=`echo "$status" | awk '{print $NF}'`
allnr=`echo "$all" | wc -w | sed 's|^[ ]*||'`

[ "$allnr" = "0" ] && echo && echo -e "\033[1;4;32m------- NOTHING CHANGED :) DONE ------------\033[0m"
echo -e "\033[1;0;31m                sum:$allnr files \033[0m" # (=allnr (number of all lines during svn status -u) to get in the end)"

# entries with ^? to add
add=`echo "$status" | grep "^?" | awk '{print $NF}'`
addnr=`echo "$add" | wc -w | sed 's|[ ]*||'`
[ "$addnr" -ge "1" ] && echo "$add" | xargs svn add   # add all new (?) files
sum=$addnr
echo -e "1/$allelines  ?:$addnr:\tsum:$sum"
[ "$allnr" = "$sum" ] && ex="yes" #&& echo -e "\033[1;4;32mEverything added now commit $commitmessage \033[0m" && commit "$commitmessage" && echo -e "\033[1;4;32mdone\033[0m"
[ "$sum" -gt "$allnr" ] && echo "ERROR: sum:$sum: greater allnr:$allnr: something is wrong" && exit

# entries with ^! to remove
remove=`echo "$status" | grep "^!" | awk '{print $NF}'`
removenr=`echo "$remove" | wc -w | sed 's|[ ]*||g'`
#echo -------------
#echo $remove | xargs
#echo ---------------
[ "$removenr" -ge "1" ] && echo "$remove" | xargs svn rm --force     # removes all files which have been deleted (!)
sum=`expr $sum + $removenr`
echo -e "2/$allelines  !:$removenr:\tsum:$sum"
[ "$allnr" = "$sum" ] && ex="yes" #echo -e "\033[1;4;32mEverything added now commit $commitmessage \033[0m" && commit "$commitmessage" && echo -e "\033[1;4;32mdone\033[0m" && exit

# entries with ^D which were deleted
D=`echo "$status" | grep "^D" | awk '{print $NF}'`   # deleted files
Dnr=`echo "$D" | wc -w | sed 's|[ ]*||g'`
sum=`expr $sum + $Dnr`
echo -e "3/$allelines  D:$Dnr:\tsum:$sum"
[ "$allnr" = "$sum" ] && ex="yes" #echo -e "\033[1;4;32mEverything added/modified now commit $commitmessage \033[0m" && commit "$commitmessage" && echo -e "\033[1;4;32mdone\033[0m" && exit


# entries with ^M which were modified
M=`echo "$status" | grep "^M" | awk '{print $NF}'`   # modified files
Mnr=`echo "$M" | wc -w | sed 's|[ ]*||g'`
sum=`expr $sum + $Mnr`
echo -e "4/$allelines  M:$Mnr:\tsum:$sum"
[ "$allnr" = "$sum" ] && ex="yes" #echo -e "\033[1;4;32mEverything added/modified now commit $commitmessage \033[0m" && commit "$commitmessage" && echo -e "\033[1;4;32mdone\033[0m" && exit

# entries with ^A for addition
A=`echo "$status" | grep "^A" | awk '{print $NF}'`   # added files
Anr=`echo "$A" | wc -w | sed 's|[ ]*||g'`
sum=`expr $sum + $Anr`;sumnoupd=$sum
echo -e "5/$allelines  A:$Anr:\tsum:$sum"
echo -e "\033[1;4;32m------ look at changed files done ----------------------------------------\033[0m"
echo




################################################################################################################################
# UPLOAD (svn commit)
################################################################################################################################
echo -e "\033[1;4;32m------ UPLOADING changed files (svn commit) ------------------------------\033[0m"
#[ "$allnr" = "$sum" -o "$ex" = "yes" ] && echo -e "\033[1;4;32mEverything added/modified now commit $commitmessage \033[0m" && commit "$commitmessage" && echo -e "\033[1;4;32mdone\033[0m" && exit
[ "$allnr" = "$sum" -o "$ex" = "yes" ] && commit "$commitmessage"
[ "$allnr" = "$sum" -o "$ex" = "yes" ] && echo -e "\033[1;4;32m------ UPLOADING changed files done --------------------------------------\033[0m" && exit
echo 

################################################################################################################################
# N: new files on server which need to be updated (svn update)
################################################################################################################################
echo -e "\033[1;4;32m------ DOWNLOADING (svn update) -----------------------------------------\033[0m"
IFSbefore=$IFS
IFS='
'
N=""
Conf=""
for line in `echo "$status" | grep "^*"`;do
    #echo line:$line:
    one=`echo "$line" | awk '{print $1}'`
    two=`echo "$line" | awk '{print $2}'`
    three=`echo "$line" | awk '{print $3}'`
    four=`echo "$line" | awk '{print $4}'`
    last=`echo "$line" | awk '{print $NF}'`
    #echo "one:$one: two:$two: three:$three: four:$four:"
    add="no"
    [ "$add" = "no" ] && [ "$one" = '*' ] && [ "$three" = "" ] && [ "$four" = "" ] && add=$two    # if the file is only on server, not locally
    [ "$add" = "no" ] && [ "$one" = '*' ] && [ "`isnumber.sh $two`" = "yes" ] && [ "$four" = "" ] && add=$three       # if file in on server and on specific revision nr
    [ "$add" != "no" ] && N="$N $add"
    [ "$add" =  "no" ] && Conf="$Conf $last"
    #echo "one:$one: two:$two: three:$three: four:$four: |||| add:$add: |||  N:$N:"
    #exit
done
IFS=$IFSbefore
#echo N:"${N}":
Nnr=`echo "$N" | wc -w`
sum=`expr $sum + $Nnr`
echo -e "*:$Nnr\tsum:$sum"
#echo sum:$sum:
#echo allnr:$allnr:
[ "$sum" -gt "$allnr" ] && echo "ERROR: sum:$sum: greater allnr:$allnr: something is wrong" && exit
echo -e "\033[1;4;32m-------------------------- svn update ------------------------------------\033[0m"
[ "$N" != "" ] && svn update $N
[ "$sumnoupd" != "0" ] && echo -e "\033[1;4;32m-------------------------- svn update done -------------------------------\033[0m"
[ "$sumnoupd" = "0" ] && echo -e "\033[1;4;32m-------------------------- svn update done, no newer local files, DONE -----\033[0m" && exit
[ "$allnr" = "$sum" ] && echo -e "\033[1;4;32m  Everything added/modified now commit $commitmessage \033[0m" && commit "$commitmessage" && echo -e "\033[1;4;32mdone\033[0m" && exit
