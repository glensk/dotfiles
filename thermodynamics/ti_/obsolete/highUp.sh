#!/bin/sh

#-----set parameters and paths------------------------
refPath=reference_0K
#-----------------------------------------------------


# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
path=$path/../../
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -f -d"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-f folder(s)   to run on given folder(s)" \
               "-d             diffs only"
  echo "Note:    default is to run for all *Ang_*K folders" 1>&2 
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   some details ..." \
        "   some more details ..."
  echo2 "   yet more ..."
  echo2 "   and so on..." \
        "   ..."
  exit
fi

checkAndSetMath
check diffs_to_take $refPath

diffsOnly=`getOption -d`; folder=`getOption -f`
if [ $folder == True ]; then l=`getValue -f`; else l=`ls -1d *Ang_*K`; fi

diffsToMake=`awk '{print $1"___"$2}' diffs_to_take`
dir=`pwd`
for i in $l; do
 if [ -d "$i" ]; then
  echo; echo -n $i; cd $i
  a=`echo $i | sed 's|\(.*\)Ang_.*|\1|'`
  if [ $diffsOnly == False ]; then
    rm -f lambda*/ion_energies_el_{freeEnergies,energies,T0K}__*eV_*kp*
    ll=`ls lambda*/*_*/{OUTCAR_*.gz,*/OUTCAR.gz} 2> /dev/null`
    run=`ls lambda*/*_*/*/OUTCAR 2> /dev/null`
    last=-1
    echo -n "  l="
    rm -f tmp_lambda*_new tmp_lambda*_all tmp_lambda*_ext
    for ii in $ll; do
      lambda=`echo $ii | sed 's|lambda\([^/]*\)/\([^/]*\)_\([^/]*\)/.*|\1|'`
      if [ $last != $lambda ]; then
        echo -n " $lambda"
        last=$lambda
      fi
      seed=`echo $ii | sed 's|lambda\([^/]*\)/\([^/]*\)_\([^/]*\)/.*|\2|'`
      structure=`echo $ii | sed 's|lambda\([^/]*\)/\([^/]*\)_\([^/]*\)/.*|\3|'`
      outcar=`echo $ii | sed 's|.*/\([^/]*\)|\1|'`
      if [ "$outcar" == "OUTCAR.gz" ]; then
        ext=`echo $ii | sed 's|.*/\([^/]*\)/OUTCAR\.gz|\1|'`
        iiNew=lambda$lambda/$seed\_$structure/OUTCAR_$ext.gz
        mv $ii $iiNew; ii=$iiNew
        rm -fr lambda$lambda/$seed\_$structure/$ext
        refOUTCAR=$dir/$refPath/$ext/$a\Ang/OUTCAR.gz
        if [ ! -e $refOUTCAR ]; then continue; fi
        echo $ii >> tmp_lambda$lambda\_new
      else
        ext=`echo $ii | sed 's/.*OUTCAR_\(.*\)\.gz/\1/'`
        refOUTCAR=$dir/$refPath/$ext/$a\Ang/OUTCAR.gz
        if [ ! -e $refOUTCAR ]; then continue; fi
      fi
      echo $ext >> tmp_lambda$lambda\_ext
      echo $ii >> tmp_lambda$lambda\_all
      nions=`zgrep "NIONS" $refOUTCAR | awk '{print $12}'`

      ref=`zgrep "free  en" $refOUTCAR | awk '{print $5}'`
      free=`zgrep "free  en" $ii | awk '{printf("%.2f",1000*($5-('$ref'))/'$nions')}'`

      ref=`zgrep "energy  w" $refOUTCAR | awk '{print $7}'`
      T0=`zgrep "energy  w" $ii | awk '{printf("%.2f",1000*($7-('$ref'))/'$nions')}'`

      ref=`zgrep "energy  w" $refOUTCAR | awk '{print $4}'`
      ener=`zgrep "energy  w" $ii | awk '{printf("%.2f",1000*($4-('$ref'))/'$nions')}'`

      n=`awk '$3=='$seed'&&$4=='$structure'{print $1}' lambda$lambda/ion_energies_el_freeEnergies__reference`
      if [ "$n" == "" ]; then
        echo 1>&2;
        echo -e "\033[31m\033[1mERROR\033[0m: no matching seed $seed and str $structure in lambda$lambda/ion_energies_el_freeEnergies__reference" 1>&2;
        exit
      fi

      file=lambda$lambda/ion_energies_el_freeEnergies__$ext; echo $n $free >> $file; sort -g $file > tmp; mv tmp $file
      file=lambda$lambda/ion_energies_el_T0K__$ext;          echo $n $T0   >> $file; sort -g $file > tmp; mv tmp $file
      file=lambda$lambda/ion_energies_el_energies__$ext;     echo $n $ener >> $file; sort -g $file > tmp; mv tmp $file
    done
    ll=`ls tmp_lambda*_all`
    echo; echo "   lambda   all   new  running  all_per_extension"
    for ii in $ll; do
      lambda=`echo $ii | sed 's|tmp_lambda\(.*\)_all|\1|'`
      all=`cat tmp_lambda$lambda\_all | wc -l`
      ext=`cat tmp_lambda$lambda\_ext | xargs -n1 | sort -u | xargs`
      extStr=""
      for iii in $ext; do
        num=`grep $iii tmp_lambda$lambda\_all | wc -l`
        extStr="$extStr  \033[31m\033[1m$num\033[0m__$iii"
      done
      if [ -e tmp_lambda$lambda\_new ]; then
        new=`cat tmp_lambda$lambda\_new | wc -l`
      else
        new=0
      fi
      running=`echo $run | xargs -n1 | grep $lambda | awk 'END{print NR}' `
      echo -e "     $lambda     $all     $new     $running   $extStr"
    done

    rm -f tmp_lambda*_new tmp_lambda*_all tmp_lambda*_ext
  fi
  echo -e "   calculating diffs ..."

  ll=`ls -1d lambda*`
  rm -f avg_diff__*
  for ii in $ll; do
    rm -f $ii/diff__* $ii/avg_diff__*
    lambda=`echo $ii | sed 's|.*lambda\([^/]*\).*|\1|'`
    lll=`ls $ii/ion_energies_el_freeEnergies__*`
    for iii in $diffsToMake; do
      ext1=`echo $iii | sed 's|\(.*\)___\(.*\)|\1|'`
      ext2=`echo $iii | sed 's|\(.*\)___\(.*\)|\2|'`
      if [ -e $ii/ion_energies_el_freeEnergies__$ext1 -a -e $ii/ion_energies_el_freeEnergies__$ext2 ]; then
        nn=`awk '{print $1}' $ii/ion_energies_el_freeEnergies__$ext1`
        sum=0;count=0
        for a in $nn; do
          f1=`grep "^$a " $ii/ion_energies_el_freeEnergies__$ext1 | awk '{print $2}'`
          f2=`grep "^$a " $ii/ion_energies_el_freeEnergies__$ext2 | awk '{print $2}'`
          if [[ "$f2" != "" && "$f1" != "" ]]; then
            diff=`echo $f1 $f2 | awk '{print $1-1*($2)}'`
            sum=`echo $diff $sum | awk '{print $1+$2}'`
            count=`expr $count + 1`
            avg=`echo $sum $count | awk '{printf("%.2f",$1/$2)}'`
            echo $a $diff >> $ii/diff__$ext1\__$ext2
            echo $a $avg >> $ii/avg_diff__$ext1\__$ext2
          fi
        done
        cp $ii/diff__$ext1\__$ext2 tmp_math
        $math >> /dev/null << EOF
dat=Import["tmp_math","Table"]
dat=Transpose[dat][[2]];
std=StandardDeviation[dat];
err=std/Sqrt[Length[dat]];
Export["tmp_stdDev",{{Mean[dat],std,err,err/2},{Null,Null,Null,Null}},"Table","FieldSeparators"->" "];
EOF
        mean=`head -n1 tmp_stdDev | awk '{printf("%.2f",$1)}'`
        std=`head -n1 tmp_stdDev | awk '{printf("%.2f",$2)}'`
        err=`head -n1 tmp_stdDev | awk '{printf("%.2f",$3)}'`
        errHalf=`head -n1 tmp_stdDev | awk '{printf("%.2f",$4)}'`
        if [ ! -e avg_diff__$ext1\__$ext2 ]; then
          echo "# lambda avg(meV/at) err/2(meV/at)    stdDev(meV/at) err(meV/at)" > avg_diff__$ext1\__$ext2
        fi
        echo "     $lambda     $mean     $errHalf                $std         $err" >> avg_diff__$ext1\__$ext2
      fi
    done
  done
  rm tmp_stdDev tmp_math
  cd $dir
 else
  echo; echo no directory $i; echo
 fi
done

