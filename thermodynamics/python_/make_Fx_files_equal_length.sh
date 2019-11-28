#!/bin/sh

out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -p -i -k -o -a -d -c";[ "$out" = "yes" ] && echo options: $options


##################################################
# join or gjoin?
#################################################
lom=`getlinuxormac`
gjoin="----------";
[ "$lom" = "mac" ] && gjoin=`which gjoin`
join=`which join`
[ -e "$gjoin" ] && join=$gjoin


filesdef=`ls -1d F{el,qh,ah}\_[db]_[0-9]* 2> /dev/null`
filesbulk=`ls -1d F{el,qh,ah} 2> /dev/null`
files=`echo "$filesdef" "$filesbulk"`
echo ... files: $files
for file1 in $files;do
    for file2 in $files;do
        [ "$file1" = "$file2" ] && continue
        wcf1=`wc -l $file1 | awk '{print $1}'`
        wcf2=`wc -l $file2 | awk '{print $1}'`
        cf1=`cat $file1 | awk '{print NF}' | sort | uniq`
        cf2=`cat $file2 | awk '{print NF}' | sort | uniq`
        #echo $cf1 | wc -w | sed 's|[ ]*||g'
        #echo $cf2 | wc -w | sed 's|[ ]*||g'
        if [ "`echo $cf1 | wc -w | sed 's|[ ]*||g'`" != "1" ] || [ "`echo $cf2 | wc -w | sed 's|[ ]*||g'`" != "1" ];then
            echo "it seems file1:$file1: contains different rows:cf1:$cf1:"
            echo "it seems file2:$file2: contains different rows:cf2:$cf2:"
            echo "file1: $file1    wcf1:$wcf1   ch1:$cf1:"
            echo "file2: $file2    wcf2:$wcf2   ch2:$cf2:"
            exit
            fi
        [ "`echo $cf2 | wc -w | sed 's|[ ]*||g'`" != "1" ] && echo "it seems $file2 contains different rows:$cf2:" && exit 
        rm -f $file1\_tmp  $file2\_tmp
        $join --nocheck-order -1 1 -2 1 $file1 $file2 | cut -d" " -f1-$cf1 > $file1\_tmp;
        $join --nocheck-order -1 1 -2 1 $file2 $file1 | cut -d" " -f1-$cf2 > $file2\_tmp; 
         
        mv $file1\_tmp $file1
        mv $file2\_tmp $file2      
done
done
