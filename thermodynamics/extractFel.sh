#!/bin/bash

#-----set default parameters and paths--------------------------------------------------------------------------------
order=3
lowestSigma=0.01
highestSigma=0.06
highestSigmaallowed=0.5     ## maximal sigmavalue before throwing an exit (check of consistency)
minNSigma=6
kB=11604.5059554883246
#---------------------------------------------------------------------------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -d -o"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage "$script"
  printOptions "-d      disable the sigma constraints" \
               "-o ORD  change order for ETto0K extrapolation (default: $order)"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   the script extracts electronic free energies from a set of OUTCARs"
  echo2 "   the OUTCARs are expected to be named OUTCAR.aLat_sigma{,.gz} where aLat"          \
        "   is the lattice constant and sigma is the smearing parameter"
  echo2 "   the number of aLats is arbitrary but the sigmas must lie on a reasonable"         \
        "   mesh with certain criteria: enough mesh points (>=$minNSigma), lowest provided"   \
        "   sigma smaller than or equal to $lowestSigma, highest sigma larger than or"        \
        "   equal to $highestSigma (sigma values in eV)"
  echo2 "   if these conditions are not satisfied an error is issued; this can be"            \
        "   avoided using the -d option to disable the constraints; BUT be careful"           \
        "   the constraint are meant to ensure good accuracy which might fail otherwise"
  echo2 "   $script first extracts the sigma->0 energies from the OUTCARs and"                \
        "   fits a polynomial of order $order (or determined by -o option) to get a highly"   \
        "   accurate extrapolation to T=0K (=ETto0K)"
  echo2 "   in a second step the electronic free energies are extracted and ETto0K is"        \
        "   subtracted to get the pure T dependent contribution"
  echo2 "   typical script sequence:"                                                         \
        "     createFolders_Fel.sh"                                                           \
        "     collectOUTCARs.sh          (one folder higher than createFolders_Fel.sh)"       \
        "     \033[1mextractFel.sh\033[0m"                                                    \
        "     getFelFit.sh"                                                                   \
        "     getThermodynamics.sh       (in separate folder)"
  exit
fi

# if -d option we disable all constraints by setting them to specific values
dOp=`getOption -d`
if [ $dOp == True ]; then
  lowestSigma=10000
  highestSigma=-10000
  highestSigmaallowed=10000
  minNSigma=1
fi

# do we change extrapolation order with -o ORD?
oOp=`getOption -o`
if [ $oOp == True ]; then
  order=`getValue -o`
  if [ "$order" == "" ]; then error "no value provided to -o option"; fi
  c=`checkInteger $order`
  if [ "$c" != ok ]; then error "value provided to -o option not an integer"; fi
fi

# mathematica kernel if needed
checkAndSetMath

l=`ls OUTCAR.*_* 2> /dev/null | wc -l`
if [ "$l" == 0 ]; then error "no OUTCAR.*_* files available"; fi

rm -f ETto0K_*Ang
l=`ls OUTCAR.*_*`

echo; echo extracting data, this may take a while ...
for i in $l; do
     a=`echo $i | sed 's/OUTCAR\.\([.0-9]*\).*/\1/'`
     s=`echo $i | sed 's/OUTCAR\..*_\([0-9]*\.[0-9]*\)[.gz]*/\1/'`
     e=`zgrep "energy  w" $i | awk 'END{print $NF}'`
     volume=`zgrep "volume of cell :" $i | awk '{print $5}' | tail -1`
     element=`zgrep POTCAR $i | sort | uniq | awk '{print $3}' | xargs | sed 's| |-|g'`

     if [ "`echo $s $highestSigmaallowed | awk '$1>=$2{print "g"}'`" = "g" ];then
      error "your sigma of $s ($i) corresponds to a temperature of `echo $s | awk '{print $1*'"$kB"'}'` K is that correct?
       consider renaming your files using: mmv -r \"OUTCAR.*_*.gz\" \"OUTCAR.#2_#1.gz\""
     fi
     lattype="-";
     ltg=`zgrep "LATTYP:" $i`
     [ "`echo $ltg | sed 's|.*LATTYP: Found a face centered cubic cell.*|fcc|'`" = "fcc" ] && lattype=fcc  ## important not to use echo "$ltg" here but echo $ltg
     [ "`echo $ltg | sed 's|.*LATTYP: Found a body centered cubic cell.*|fcc|'`" = "fcc" ] && lattype=bcc  ## important not to use echo "$ltg" here but echo $ltg
     echo $s $e $lattype $element $volume >> ETto0K_$a\Ang
     #echo $a $s $e $lattype $volume $element
     #echo ETto0K_$a\Ang; exit
done
lattype=`cat ETto0K_*Ang | awk '{print $3}' | sort | uniq`
element=`cat ETto0K_*Ang | awk '{print $4}' | sort | uniq`


rm -f Fel_*Ang Fel*eV extractFel.log
echo "# lattype: $lattype element: $element" > extractFel.log
echo; echo "# using $order. order DOS based extrapolation of ETto0K" >> extractFel.log
echo -e "\033[1m`tail -n1 extractFel.log`\033[0m"
echo "# aLat(A) ETto0K(eV) MeanDev(meV) MaxDev(meV)" >> extractFel.log
tail -n1 extractFel.log
l=`ls ETto0K_*Ang`

for i in $l; do ## loop over all OUTCARS
    a=`echo $i | sed 's/ETto0K_\(.*\)Ang/\1/'`
    mes=`awk 'BEGIN{lowest=10;highest=0}
              $1<lowest{lowest=$1}
              $1>highest{highest=$1}
              END{     if (lowest>'$lowestSigma') printf "lowest sigma (%s) larger than allowed ('$lowestSigma'); check -help",lowest;
                  else if (highest<'$highestSigma') printf "highest sigma (%s) smaller than allowed ('$highestSigma'); check -help",highest;
                  else if (NR<'$minNSigma') printf "number of sigma values (%s) smaller than allowed ('$minNSigma'); check -help",NR;
                  else print "ok"}' $i`
    if [ "$mes" != "ok" ]; then error "for aLat $a $mes"; fi
    cat $i | sort -g > _tmp; cp _tmp $i
    $math >> _tmp_log << EOF
        e0=Import["_tmp","Table"];
        
        basis[T_]:=Table[T^i,{i,0,$order}];
        ff[e_,T_]:=1/(Exp[e/T]+1)
        sel[e_,T_]:=(1 - ff[e, T]) Log[1 - ff[e, T]] + ff[e, T] Log[ff[e, T]]
        int[a_,e_,T_]=Integrate[a sel[e,T],e];
        Sel[a_,T_]:=int[a,10,T]
        Fel[a_,T_]:=T Sel[a,T]
        fitFunc[T_, coef_] := coef[[1]] + Fel[ Sum[coef[[1+k]] basis[T][[k]], {k, basis[T] // Length}], T]
        squares[coef_] := Sum[(e0[[j,2]] - fitFunc[e0[[j,1]], coef])^2, {j, e0// Length}]
        
        startcoef = Table[0.0001, {1+(basis[T] // Length)}];
        coef = Array[c, {1+(basis[T] // Length)}];
        
        fit = FindMinimum[ squares[coef], {coef // Flatten, startcoef // Flatten}\[Transpose]][[2]];
        last=Transpose[e0][[1,-1]];
        fitted=Table[{T,fitFunc[T,coef]/.fit},{T,0.0001,last,last/100}];
        delta=Abs[Transpose[e0][[2]]-(fitFunc[T,coef]/.fit/.T->Transpose[e0][[1]])];
        mean=Mean[delta];
        max=Max[delta];
        zero=fitFunc[0.0001,coef]/.fit;
        output={zero,1000mean,1000max};
        Export["_tmp_output",output,"Table"];
        Export["_tmp_fitted",fitted,"Table"];
EOF
    # check whether mathematica has ran properly
    c=`grep "Mathematica cannot find a valid password" _tmp_log`
    if [ -n "$c" ]; then error "mathematica licence unavailable"; fi
    rm -f _tmp_log

    echo >> $i; cat _tmp_fitted >> $i
    zero=`head -n1 _tmp_output`
    ll=`ls OUTCAR.$a\_*`
    for ii in $ll; do
        f=`zgrep "free  e" $ii | awk 'END{printf("%.6f",1000*($(NF-1)-1*('"$zero"')))}'`
        s=`echo $ii | sed 's/OUTCAR\..*_\([0-9]*\.[0-9]*\)[.gz]*/\1/'`
        T=`echo $s | awk '{printf("%.8f",$1*'$kB')}'`
        echo $T $f >> Fel_$a\Ang
        echo $a $f >> Fel_$s\eV
    done
  
    cat Fel_$a\Ang | sort -g > _tmp; mv _tmp Fel_$a\Ang
    echo "$a   `cat _tmp_output | xargs | awk '{printf("%.6f   %.3f   %.3f",$1,$2,$3)}'`" >> extractFel.log
    tail -n1 extractFel.log
done

rm _tmp_output _tmp_fitted

