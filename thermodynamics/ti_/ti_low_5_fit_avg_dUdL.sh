#!/bin/sh


########################################################################
## check for mathematica
########################################################################
#path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`

path=`echo \`readlink -f "$0"\` | sed 's|'\`basename $0\`'||'`
scriptpath=$path; path=$path/../

script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -s -i -u -f -n -V -a -B -Bd -k"

# mathematica kernel if needed
checkAndSetMath

host=`hostname`

########################################################################
## check for fitting script
########################################################################
fitting=`which low_fit_avg_dUdL.txt`
#echo fitting: 
[ ! -e "$fitting" ] && echo fitting script does not exist && exit


########################################################################
## check for input files
########################################################################


[ "`pwd | sed 's|.*/||' | grep "Ang_.*K" | wc -w | sed 's|[ ]*||g'`" != "1" ] && echo "`pwd` you have to be in AngK folder" && exit

list="fre ene eS0"
output=Fah
[ -e "$output" ] && rm $output 


for var in $list;do
[ ! -e "avg_dUdL_$var" ] && echo "`pwd` : avg_dUdL_$var does not exist" && continue
[ "`cat avg_dUdL_$var | wc -l | sed 's|[ ]*||g'`" -le "3" ] && echo "`pwd` : avg_dUdL_$var has no input" && continue
## check lambda values

lambdas=`cat avg_dUdL_$var | tail -n+4 | awk '{print $1}' | xargs`
[ "`echo $lambdas | wc -l | sed 's|[ ]*||g'`" = "1" ] && [ "$lambdas" = "1.0" ] && echo "just lambda 1.0, exit and nofit created!" && exit
lambda_greater_half=no
for i in $lambdas;do
[ "`echo 0.4999 $i | awk '($1<=$2){print "yes"};($1>$2){print "no"}'`" = "yes" ] && lambda_greater_half=yes && break
done
[ "$lambda_greater_half" = "no" ] && echo -e "`pwd` \033[31m\033[1m PROBLEM \033[0m:no lambda greater 0.5" && exit

######################################################################
## this part if at least on lambda >= 0.5
######################################################################
cp avg_dUdL_$var avg_in
$math < $fitting
sed 's|\(.*\)|'"$var"'_\1|' Fahout.dat >> $output
echo "" >> $output; 
rm -f Fahout.dat
[ -e "fit_tangens" ] && mv fit_tangens fit_tangens_$var
[ -e "fit_cubic" ] && mv fit_cubic fit_cubic_$var

done

[ ! -e "$output" ] && echo "`pwd` : $output does not exist" && exit
awk '{printf "%-20s %7.2f %7.2f %7.2f\n",$1,$2,$3,$4}' $output >> $output.tmp; mv $output.tmp $output
rm avg_in
