#!/bin/bash

out=no #yes #(print additional info for debugging when yes)
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo path: $path
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`;[ "$out" = "yes" ] && echo script: $script
options=$*; . $path/../utilities/functions.include; checkOptions "-h -help -v -vol -s -f";[ "$out" = "yes" ] && echo options: $options

if [ `getOption -h` = True ] || [ `getOption -help` = True ]; then
  usage `basename $0`
  printOptions " " \
               "-f              fitted surface; default:Fah_surface_fit or Fqh_fromExactFreqs_fit_order2"\
               "-vol [VOLUME]   voume to fit"\
               "-s [real]       scaling factor"\
               "-e              only show volumes and factor and exit"\
               "-v              verbose mode"
  exit
fi

#### define modules
#isnumber() {
#for isnumber in $*;do
#        number=`echo "$isnumber" | grep -o "[-]*[+]*[0-9.]*[0-9]*"`
#        [ "$number" != "$isnumber" ] && echo no1 && exit
#        [ "`echo $number | grep "[0-9]" | wc -w`" -lt "1" ] && echo no2 && exit
#        [ "`echo "$isnumber" | grep -o "[0-9.]*" | wc -w`" -gt "1" ] && echo no3 && exit
#done
#echo yes
#}


[ "`getOption -f`" = "True" ] && surface=`getValue -f`


### Fah_surface_fit
surface=Fah_surface_fit
if [ -f "$surface" ];then
alats=`find . -maxdepth 1 -mindepth 0 -type f -name "Fah*Ang" | sed 's|./Fah_||g' | sed 's|Ang||g' | xargs`
sc=`grep "^sc" fit.input | sed 's|.*=||' | sed 's|[ ]*||' | sed 's|;.*||'`
atoms=`grep "^nAtoms" fit.input | sed 's|.*=||' | sed 's|[ ]*||' | sed 's|;.*||'`
sf=`grep "^structureFactor" fit.input | sed 's|.*=||' | sed 's|[ ]*||' | sed 's|;.*||'`
volumes=""
makefile="Fah_"
fi


#### Fqh_fromExactFreqs_fit_order2 
[ ! -e "$surface" ] && surface=Fqh_fromExactFreqs_fit_order2
if [ -f "$surface" ];then
if [ "$alats" = "" ];then
alats=`find . -maxdepth 1 -mindepth 0 -type f -name "Fqh_fromExactFreqs_[0-9.]*" | sed 's|./Fqh_fromExactFreqs_||g' | xargs`
atoms=`cat fitFqh.input | grep "^s =" | awk '{print $3}' | sed 's|[ ]*||' | sed 's|;.*||'`
sc=`grep "^sc" fitFqh.input | sed 's|.*=||' | sed 's|[ ]*||' | sed 's|;.*||'`
sf=`grep "^structureFactor" fitFqh.input | sed 's|.*=||' | sed 's|[ ]*||' | sed 's|;.*||'`
makefile="Fqh_fromExactFreqs_"
fi
fi


### fitFqh.input
if [ ! -e "$surface" ];then
    echo 1
if [ -e "fitFqh.input" ];then
    echo 2
if [ "$alats" = "" ];then
    echo 3
base=`grep baseName fitFqh.input | sed 's|.*{"\(.*\)"}.*|\1|'`
alats=`find . -maxdepth 1 -mindepth 0 -type f -name "$base[0-9.]*" | sed 's|./'"$base"'||g' | xargs`
atoms=`cat fitFqh.input | grep "^s =" | awk '{print $3}' | sed 's|[ ]*||' | sed 's|;.*||'`
sc=`grep "^sc" fitFqh.input | sed 's|.*=||' | sed 's|[ ]*||' | sed 's|;.*||'`
sf=`grep "^structureFactor" fitFqh.input | sed 's|.*=||' | sed 's|[ ]*||' | sed 's|;.*||'`
makefile=$base
surface=$base\fit_order2
fi;fi;fi


[ "`getOption -v`" = "True" ] && echo alats:$alats:
[ "`getOption -v`" = "True" ] && echo atoms:$atoms:
[ "`getOption -v`" = "True" ] && echo sc:$sc:
[ "`getOption -v`" = "True" ] && echo sf:$sf:
[ "`getOption -v`" = "True" ] && echo makefile:$makefile:
[ "`getOption -v`" = "True" ] && echo surface:$surface:


for a in $alats;do
    vol=`echo "($a*$sc)^3" | bc -l`
    factor=`echo "1/$atoms" | bc -l`
    [ "`getOption -v`" = "True" ] && echo $a $vol $factor
    [ "$sf" != "1" ] && vol=`echo "$vol/$sf" | bc -l`
    [ "`getOption -v`" = "True" ] && echo $a $vol $factor
    volumes="$volumes $vol"
    volumes_alats="$volumes_alats $vol\_$a"
done

[ "`getOption -v`" = "Treu" ] && exit



# get from commandline
[ "`getOption -vol`" = "True" ] && volumes=`getValue -vol`
[ "`getOption -s`" = "True" ] && factor=`getValue -s`


# checks
echo ------------------------------------------------------
echo -----------  `basename $0` -------------
echo ------------------------------------------------------
for a in $alats;do
    vol=`echo "($a*$sc)^3" | bc -l`
    [ "$sf" != "1" ] && vol=`echo "$vol/$sf" | bc -l`
    factor=`echo "1/$atoms" | bc -l`
    #echo $a $volume $factor

    ## checks
    [ ! -f $surface ] && echo "surface :$surface: does not exist" && exit
    [ "`isnumber.sh $vol`" != "yes" ] && echored "volume $vol is not a number" && exit
    [ "`isnumber.sh $factor`" != "yes" ] && echored "factor $factor is not a number" && exit
    echo "surface:$surface: factor:$factor: `echogreen "alat:$a: vol:$vol:"`"
done
echo ------------------------------------------------------

#### make fromfit folder
fitfolder=fromfit
rm -rf $fitfolder
zuvor=`pwd`
for a in $alats;do
    vol=`echo "($a*$sc)^3" | bc -l`

    [ "$sf" != "1" ] && vol=`echo "$vol/$sf" | bc -l`

    factor=`echo "1/$atoms" | bc -l`
    
    ## checks
    [ "`isnumber.sh $vol`" != "yes" ] && echored "volume $vol is not a number" && continue
    [ "`isnumber.sh $factor`" != "yes" ] && echored "factor $factor is not a number" && continue

    ### get sur
    lines=`wc -l $surface | awk '{print $1+1}'`
    echo "# vol:$vol lines:$lines"
    mkdir -p $fitfolder
    fit=$fitfolder/$makefile$a
    rm -f $fit
    for i in `seq 1 $lines`;do  # this runs over temperatures
        
        line=`awk 'NR=='"$i"'' $surface`
        order=`echo $line | wc -w | sed 's|[ ]*||'`
    
        [ "$order" = "4" ] && y=`echo $vol $line | awk '{print $3,$4,$5}'    | awk '{print $1+$2*'"$vol"'+$3*'"$vol"'*'"$vol"'}'`
        [ "$order" = "5" ] && y=`echo $vol $line | awk '{print $3,$4,$5,$6}' | awk '{print $1+$2*'"$vol"'+$3*'"$vol"'*'"$vol"'+$4*'"$vol"'*'"$vol"'*'"$vol"'}'`
        
        [ "`getOption -v`" = "True" ] && echo temperature_i:$i order:$order: vol:$vol: line:$line  wc:`echo $line | wc -w`: factor:$factor:  y:$y:
        [ ! -e "$fit" ] && touch $fit

        #echo $y $line | awk '{print $2,$1*'"$factor"'}' | tee -a Fvib_$vol\_fitted
        echo $y $line | awk '{printf "%.0f %.10f\n", $2,$1*'"$factor"'}' >> $fit
    done
done


cd $zuvor
if [ "$makefile" != "Fah_" ];then
for vola in $volumes_alats;do
    vol=`echo "$vola" | sed 's|_.*||' | grep -o "[0-9.]*"`
    a=`echo "$vola" | sed 's|.*_||'`
    from=$makefile$a
    [ ! -e "$from" ] && from=$from\Ang
    fit=$fitfolder/$makefile$a
    
    #echo "$vola --> :$vol:  :$a:"
    if [ -e "$from" ] && [ -e $fit ];then    #Fvib keine ""
        wc1=`wc -l $from | awk '{print $1}'`
        wc2=`wc -l $fit | awk '{print $1}'`
        if [ "$wc1" = "$wc2" ];then
            echogreen "$from and $fit found"
            file=$fitfolder/correction_add_$a
            rm -f $file
            paste $from $fit | awk '{print $1,$2-$4}' > $file
        else
            echored "$from :$wc1: and $fit :$wc2: differnet length"
        fi
    else
        echored "$from or $fit NOT found"
    fi
done
fi
