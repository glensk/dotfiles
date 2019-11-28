#!/bin/sh
out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -v -qm";[ "$out" = "yes" ] && echo options: $options

optionstmp=$options
options=`echo $options | sed 's| -qm||' | sed 's| -v||'`
nr_opt=`echo $options | wc -w | sed 's|[ ]*||'`
if [ `getOption -h` = True ] || [ `getOption -help` = True ] || [ "$nr_opt" = "0" ] || [ "$nr_opt" != "1" -a "$nr_opt" != "2" ]; then
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e " \033[1m`basename $0` folderToOldRef folderToNewRef\033[0m     folders must contain old and new Fqh_fromExactFreqs_classical_*" 1>&2
  echo -e " \033[1m`basename $0` fileWithFolders\033[0m                   folderToOldRef in 1. line; folderToNewRef in 2. line" 1>&2
  printOptions " " \
               "-qm         swich from classical to quantenmechanical (Fqh_fromExactFreqs_classical_* -> Fqh_fromExactFreqs_*)" \
               "            NOTE: -qm is just for testing purposes; since we run classical MD's we also need classical changeHesseRef" \
               "-v          verbose mode"
  exit
fi
options=$optionstmp

################
#  -qm
################
classical=classical_
[ "`getOption -qm`" = "True" ] && classical=""

l=`ls Fah_[0-9.]*K | wc -l | sed 's|[ ]*||g'`
if [ "$l" == 0 ]; then
  echo "" 1>&2; echo -e "  \033[31m\033[1mERROR\033[0m: no Fah_* files available" 1>&2; exit 1
fi

################
# check if HssseRef was changed before
###############
if [ -e HesseRef_changed ]; then
  echo 1>&2; echo -e "  \033[31m\033[1mERROR\033[0m: reference has been changed already; to force changing again delete HesseRef_changed file" 1>&2; exit 1
fi

if [ $# == 1 ]; then old=`sed -n '1p' $1`; new=`sed -n '2p' $1`; else old=$1; new=$2; fi
if [ ! -d $old -o ! -d $new ]; then
  echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: old or new ref folder not existing" 1>&2; exit 1
fi

check() {
  if [ ! -e $1 ]; then echo "" 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: file $1 not existing" 1>&2; exit 1; fi
}

Fit() {
        checkAndSetMath
        t=$1;file=$2;checkf=$3
        [ $# != 3 ] && echo 'need 2 input' && exit
        [ "`getOption -v`" = "True" ] && echo math:$math:
        rm  -f tmpfileon tmpfileon.dat tmpfileonout.dat tttk
        cp  $2 tmpfileon
        sed -i 's|\t| |g' tmpfileon
        sed -i '/^ *$/d' tmpfileon   # necessary otherwise {} problems when inport
        #echo "------------ start fit -------------"
        #cat tmpfileon
        #echo "-------------done ------------------"
#        $math << EOF
$math >>/dev/null << EOF
    data=Import["tmpfileon","Table"];
    (*Print["--1--",data//Length," ",data//Dimensinos];*)
    datain=data;
    get = $t;

    Print["--2--"];
    kk = data[[1 ;; -1, 1]];

    Print["--3--"];
    (*kk = If[Last[kk] == {}, Delete[kk, -1]];*)

    Print["--4--"];
    (*Print["kk: ",kk];*)
    
    xdata = Select[Select[kk, # <= (get + 3) &], # >= (get - 3) &];
    (*Print["xdata: ",kk];*)
    datafit=Table[Select[data, #[[1]] == xdata[[i]] &][[1]], {i, 1, xdata // Length}];
    (* in a test the linear version gave 0.0005 meV/atom error in the extrapolation*)
    If[Length[xdata] < 2, {PRINT["ERROR:",Length[xdata]],Quit[]}];
    parabola = Fit[datafit, {1, x, x^2}, x];
    kb=parabola /. x -> get;
    dataout=Partition[{datain, {get, kb}} // Flatten, 2] // Sort;
    Print["len1: ",Length[datain]];
    Print["len2: ",Length[dataout]];
    Export["tmpfileonout.dat",dataout];
EOF
        [ ! -e tmpfileonout.dat ] && echo "file tmpfileonout.dat does not exst" && exit
        #echo "------------ stop fit -------------"
        outcheck=`awk 'BEGIN{t="'$t'"; s=0};(($1-t)^2)^(1/2)<.0001{s++};END{print s}' tmpfileonout.dat`
        outF=`awk 'BEGIN{t="'$t'"; tt=-1000};(($1-t)^2)^(1/2)<.0001&&(($1-t)^2)^(1/2)<((t-tt)^2)^(1/2){tt=$1;f=$2};END{print f}' tmpfileonout.dat`
        rm -f tttk
        [ "$checkf" = "check" ] && echo $outcheck > tttk
        [ "$checkf" = "F" ] && echo $outF > tttk
        rm -f tmpfileon tmpfileon.dat tmpfileonout.dat
}

check Fah_from_fit_*
structureFactor=`awk 'NR==1{print $1}' Fah_from_fit_*`

l=`ls Fah_[0-9.]*K`
echo; echo -n " running checks ... ";echo
#for i in $l; do
#  t=`echo $i | sed 's|Fah_\(.*\)K|\1|'`
#  aLats=`awk '{print $1}' $i`
#  for a in $aLats; do
#    echogreen " ---> t:$t a:$a"
#    check $old/Fqh_fromExactFreqs_$classical$a; check $new/Fqh_fromExactFreqs_$classical$a
#    oldF=`awk 'BEGIN{t="'$t'"; s=0};(($1-t)^2)^(1/2)<.0001{s++};END{print s}' $old/Fqh_fromExactFreqs_$classical$a`
#    newF=`awk 'BEGIN{t="'$t'"; s=0};(($1-t)^2)^(1/2)<.0001{s++};END{print s}' $new/Fqh_fromExactFreqs_$classical$a`
#    
#    ## if new value does not exist 
#    makefit="no"
#    [ $oldF == 0 ] && Fit $t $old/Fqh_fromExactFreqs_$classical$a check && oldF=`cat tttk` && makefit="yes"
#    [ $newF == 0 ] && Fit $t $new/Fqh_fromExactFreqs_$classical$a check && newF=`cat tttk` && makefit="yes"
#    rm -f tttk 
#    [ "$oldF" = "0" ] && echored "ERROR: T=$t K not existing in $old/Fqh_fromExactFreqs_$classical$a" && exit
#    [ "$newF" = "0" ] && echored "ERROR: T=$t K not existing in $new/Fqh_fromExactFreqs_$classical$a" && exit
#    
#    #echo oldF:$oldF
#    #echo newF:$newF
#  done
#done
echo ok; echo

rm -f Fah_*K_vol Fah_*Ang Fah_surface*
for i in $l; do
  t=`echo $i | sed 's|Fah_\(.*\)K|\1|'`
  aLats=`awk '{print $1}' $i`
  echo "$i"
  echo "old":$old
  echo "new":$new
  echo "classical:$classical"
  echo "aLat   oldFah  newFah"
  for a in $aLats; do
    echogreen " ---> t:$t a:$a"
    [ "`getOption -qm`" = "False" ] && [ ! -e "$old/Fqh_fromExactFreqs_$classical$a" ] && echo old:$old:in:`pwd` && getSingleSpeciesPhonons.sh -folder $old -Fe -c -a $a
    [ "`getOption -qm`" = "False" ] && [ ! -e "$new/Fqh_fromExactFreqs_$classical$a" ] && echo new:$new:in:`pwd` && getSingleSpeciesPhonons.sh -folder $new -Fe -c -a $a
    check $old/Fqh_fromExactFreqs_$classical$a; check $new/Fqh_fromExactFreqs_$classical$a
    
    #checks
    oldFc=`awk 'BEGIN{t="'$t'"; s=0};(($1-t)^2)^(1/2)<.0001{s++};END{print s}' $old/Fqh_fromExactFreqs_$classical$a`
    oldoldFc=$oldFc
    newFc=`awk 'BEGIN{t="'$t'"; s=0};(($1-t)^2)^(1/2)<.0001{s++};END{print s}' $new/Fqh_fromExactFreqs_$classical$a`
    newnewFc=$newFc
    
    echo oldFc:$oldFc: newFc:$newFc:
    ## get values
    if [ $oldFc == 0 ];then  ## value has to be interpolated
        Fit $t $old/Fqh_fromExactFreqs_$classical$a check && oldFc=`cat tttk` && makefit="yes"
        [ "$oldFc" = "0" ] && echored "ERROR: T=$t K not existing in $old/Fqh_fromExactFreqs_$classical$a" && exit
        rm -f tttk 
        Fit $t $old/Fqh_fromExactFreqs_$classical$a F && oldF=`cat tttk` && makefit="yes"
        #echo oldF:$oldF:
        rm -f tttk 
    else
        oldF=`awk 'BEGIN{t="'$t'"; tt=-1000};(($1-t)^2)^(1/2)<.0001&&(($1-t)^2)^(1/2)<((t-tt)^2)^(1/2){tt=$1;f=$2};END{print f}' $old/Fqh_fromExactFreqs_$classical$a`
    fi
    [ "`checkReal $oldF`" != "ok" ] && echored "oldF: $oldF : is not a number oldFc:$oldFc: oldoldFc:$oldoldFc: newFc:$newFc: t:$t:" && exit


    if [ $newFc == 0 ];then
        Fit $t $new/Fqh_fromExactFreqs_$classical$a check && newFc=`cat tttk` && makefit="yes"
        [ "$newFc" = "0" ] && echored "ERROR: T=$t K not existing in $new/Fqh_fromExactFreqs_$classical$a" && exit
        rm -f tttk 
        Fit $t $new/Fqh_fromExactFreqs_$classical$a F && newF=`cat tttk` && makefit="yes"
        rm -f tttk 
    else
        newF=`awk 'BEGIN{t="'$t'"; tt=-1000};(($1-t)^2)^(1/2)<.0001&&(($1-t)^2)^(1/2)<((t-tt)^2)^(1/2){tt=$1;f=$2};END{print f}' $new/Fqh_fromExactFreqs_$classical$a`
    fi
    [ "`checkReal $newF`" != "ok" ] && echored "newF: $newF : is not a number " && exit


	echo oldF:$oldF newF:$newF
    delta=`echo $oldF $newF | awk '{printf("%.8f",$1-$2)}'`
    oldFah=`awk '$1=='$a'{printf("%6.2f",$2)}' $i`
    awk '$1=='$a'{printf("%s %.2f %s\n",$1,$2+'$delta',$3)};$1!='$a'{print $0}' $i > tmp; mv tmp $i
    newFah=`awk '$1=='$a'{printf("%6.2f",$2)}' $i`
    aa=`echo $a | awk '{printf("%.2f",$1)}'`
    echo -e "$aa  \033[1m$oldFah  \033[31m$newFah\033[0m"


    err=`awk '$1=='$a'{print  $3}' $i`
    v=`echo $a | awk '{printf("%.2f",$1^3/'$structureFactor')}'`
    echo $t $newFah $err >> Fah_$a\Ang
    echo $v $newFah $err >> Fah_$t\K_vol
    echo $t $a $newFah $err >> Fah_surface
    echo $t $v $newFah $err >> Fah_surface_vol
  done
  echo
done

echo $old > changeHesseRef_folders
echo $new >> changeHesseRef_folders
touch HesseRef_changed
echogreen "HesseRef_changed and changeHesseRef_folders written (Fah_xxxAng Fah_xxxK Fah_xxxK_vol Fah_surface Fah_surface_vol)"

