#!/bin/bash
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -w -c -r -t -a -d -wh"
if [ "`getOption -d`" = "True" ];then 
    echo path:$path
    echo script:$script
    echo options:$options
fi
file=$path/utilities/db_melting_points_CRC_Handbook
file_we=$path/utilities/db_melting_points_webelements
file_calphad=$path/utilities/db_CALPHAD_PURE5_SGTE.TDB
printdebuginfo=false
[ ! -e "$file" ] && echo $file not found && exit


#echo path:$path
#echo file:$file

elements_CRC=" Ac Ag Al Am Ar As At Au B Ba Be Bi Bk Br C Ca Cd Ce Cf Cl Cm Co Cr Cs Cu Dy Er Es Eu F Fe Fm Fr Ga Gd Ge H He Hf Hg Ho I In Ir K Kr La Li Lr Lu Md Mg Mn Mo N Na Nb Nd Ne Ni No Np O Os P Pa Pb Pd Pm Po Pr Pt Pu Ra Rb Re Rh Rn Ru S Sb Sc Se Si Sm Sn Sr Ta Tb Tc Te Th Ti Tl Tm U V W Xe Y Yb Zn Zr "

elements_calphad=" \
Ac Ag Al Am Ar As At Au Bb Ba Be Bi Bb Cc Ca Cd Ce Cf Cl \
Cm Co Cr Cs Cu D1 Dt D2 Dy Er Es Eu Ff Fe Fm Fr Ga Gd Ge \
Hh He Hf Hg Ho Ii In Ir Kk Kr La Li Lu Mg Mn Mo Nn Na Nb \
Nd Ne Ni Np Oo Os Pp Pa Pb Pd Pm Po Pr Pt Pu Ra Rb Re Rh \
Rn Ru Ss Sb Sc Se Si Sm Sn Sr T1 T2 Ta Tb Tc Te Th Ti Tl \
Tm Uu Vv Ww Xe Yy Yb Zn Zr "

errorexit() {
echo ""
echo "        This script will return the Melting point of an Element."
echo ""
echo "        The data is taken from $file" 
echo "        Ref: Haynes, William M. [Hrsg.] : CRC handbook of chemistry and physics."
echo "             Boca Raton, Fla.[u.a.] : CRC Press, Taylor & Francis Group, 2007."
echo " "
echo "Usage: `basename $0` [ELEMENT] [OPTIONS] (ELEMENT is not case sensitive)"
echo "  e.g. `basename $0` ag"
echo ""
echo "[ELEMENT]:"
echo "$elements_CRC"
echo ""
echo "[OPTIONS]:"
echo "       -w:    swich to melting point extracted from http://www.webelements.com"
echo "       -c:    swich to melting point extracted from the CALPHAD PURE5 SGTE database"
echo "       -r:    round the temperature up if necessary (to next higher integer)"
echo "       -t:    print real name of the element and exit"
echo "       -h:    print this help"
echo "       -a:    list all values and corresponding database without rounding"
echo "       -d:    debug modus, print additional infos"
echo "       -wh:   additional information for Tmelts from http://www.webelements.com"
echo ""
exit
}
#echo hallo :`echo $*`:
#echo hallo :`echo "$*"|wc -w`:
#### checks
#[ "`echo $* | wc -w | sed 's| ||g'`" -ge "4" ] && echo "you can use maximal two options like -c -r but not more:$*:" && exit   # getMeltingpoint.sh Pb -c -r will have 3 arguments
#options=""; [ "`echo $* | wc -w`" -ge "2" ] && options=`echo $* | grep -o "\-[a-z]"` # | sed 's|-||g'`
#options_known="-c -r -h -a -t -d"
#optionsout=""
#echo op:$options:
#for i in $options;do
#    grepfor=`echo "$i" | sed 's|-||'`  # sollte nur w r oder h sein
#    echo i:"$i": grepfor:$grepfor: wc:`echo "$i" | wc -c`: getminus:`echo "$i" | cut -c 1`:
#    [ "`echo "$i" | cut -c 1`" != "-" ] && echo OPTIONS have to start with a minus "-" && echo "OPTION \"$i\" is not known; try `basename $0` -h for help" && exit
#    [ "`echo "$i" | wc -c | sed 's| ||g'`" != "3" ] && echo "OPTION \"$i\" is not known; try `basename $0` -h for help" && exit  ## 3 is correct bc wc adds 1
#    [ "`echo $options_known | grep -o "$grepfor" | wc -w`" != "1" ] && echo "option \"$i\" is not known; try `basename $0` -h for help" && exit
#    optionsout="$optionsout $grepfor"
#    #echo i:$i: grepfor:$grepfor: opout:$optionsout:
#done
#[ "`echo $1 | grep "\-h" | wc -w | sed 's| ||g'`" != "0" ] && errorexit
#[ "`echo $optionsout | grep -o "h" | wc -w | sed 's| ||g'`" != "0" ] && errorexit

[ "`getOption -h`" = "True" ] && errorexit
[ "`getOption -help`" = "True" ] && errorexit

if [ "`getOption -wh`" = "True" ]; then
  echo
  echo " there are differences between Tmelts from the CRC Handbook and webelements"
  echo " differences for the d-elements (others are the same):
   CRC       web       delta
   Ir 2719   Ir 2739   -20
   Lu 1936   Lu 1925   11
   Mo 2895   Mo 2896   -1
   Re 3458   Re 3459   -1
   Rh 2236   Rh 2237   -1
   Ru 2606   Ru 2607   -1
   Ti 1943   Ti 1941    2
   W 3687    W 3695    -8
   Y 1795    Y 1799    -4
   Zr 2127   Zr 2128   -1"
  echo
  echored " the more recent CRC data should be preferred"
  exit
fi

### skript begin
element=`echo $1 | tr "[:lower:]" "[:upper:]"`
[ "`getOption -d`" = "True" ] && echo element:$element
check=`echo "$elements_CRC" | grep -i -o " $element " | sed 's| ||g'`

[ "`getOption -d`" = "True" ] && echo check:$check 
[ "`echo $check | wc -w | sed 's| ||g'`" != "1" ] && echo element \""$1"\" is not is not an available ELEMENT:"$elements_CRC" && echo "check:$check: is not 1 word; counted `echo $check | wc -w` words!" && exit



####################################
# get the melting pint
####################################
#echo element:$element:

## a) from CRC Handbook file
meltpoint=`grep -i "^$element " $file | awk '{print $2}'`
[ "`getOption -d`" = "True" ] && echo CRC_Handbook:$meltpoint

## b) from webelements
name_of_element=`grep -i "^[0-9]* $element " $file_we | awk '{print $3}'`
[ "`getOption -d`" = "True" ] && echo name_of_element:$name_of_element
[ "`getOption -t`" = "True" ] && echo "$name_of_element" && exit

webelements=`grep -i -A2 "^[0-9]* $element " $file_we | grep "^WEL" | awk '$3=="K"{print $2}'`  # this is necessary to have a first estimate
[ "`getOption -d`" = "True" ] && echo webelements:$webelements
#echo 3 web:$webelements
calphad=`sed -n '/FUNCT GLIQ'"$element"'/,/\!/p' $file_calphad | grep " Y" | grep -o "[0-9.]* Y" | awk '{print sqrt(($1-'"$webelements"')^2), $1}' |sort -n | head -1 | awk '{print $2}'`

[ "`getOption -d`" = "True" ] && echo calphadcheck:$calphad
if [ "$calphad" = "" ];then
    calphad=`sed -n '/FUNCT GHSER'"$element"'/,/\!/p' $file_calphad | grep " Y" | grep -o "[0-9.]* Y" | awk '{print sqrt(($1-'"$webelements"')^2), $1}' |sort -n | head -1 | awk '{print $2}'`
fi
[ "`getOption -d`" = "True" ] && echo calphad:$calphad
[ "$calphad" = "" ] && echo not found in CALPHAD db && exit
#echo calphad:$calphad

########################################
# decide which output to print
########################################
out=$meltpoint
[ "`getOption -w`" = "True" ] && out=$webelements
[ "`getOption -c`" = "True" ] && out=$calphad
[ "`getOption -d`" = "True" ] && echo vor runden:$out
[ "`getOption -r`" = "True" ] && out=`echo $out | awk '{printf "%.0f\n", $1+0.4999999999}'`
[ "`getOption -r`" = "True" ] && meltpoint=`echo $meltpoint | awk '{printf "%.0f\n", $1+0.4999999999}'`
[ "`getOption -r`" = "True" ] && webelements=`echo $webelements | awk '{printf "%.0f\n", $1+0.4999999999}'`
[ "`getOption -r`" = "True" ] && calphad=`echo $calphad | awk '{printf "%.0f\n", $1+0.4999999999}'`
[ "`getOption -d`" = "True" ] && echo nach runden:$out
[ "`getOption -w`" = "True" ] && out="$out  <-- check -wh option for possible issues with webelements"

if [ "`getOption -a`" = "True" ];then
  echo -e "$meltpoint   \t:  CRC Handbook"
  echo -e "$webelements   \t:  webelements  <-- check -wh option for possible issues with webelements"
  echo -e "$calphad   \t:  calphad"
  exit
fi

echo $out

#echo differences foudn in -- ir np pa th pu
#echo $webelements $calphad | awk '{print $1-$2}'
