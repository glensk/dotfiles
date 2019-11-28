#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions
  echo "Note:    EATOM is not included!" 1>&2 
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   extract energies per atom from vasp Murn calculation"
  echo2 "   results are assumed in OUTCAR.aLat.gz files"
  echo2 "   T->0K energies and free energies are extracted"
  echo2 "   EATOM is not included!"
  exit
fi

eVTomeV=1000; eVToRyd=0.0734986508052261761

echoFile () { echo "$1" >> log_extractEnergies
}


# get OUTCARS
l1=`find . -maxdepth 1 -mindepth 1 -name "OUTCAR.*.gz" | sed 's|^./||'`;geta="l1";l=$l1
[ "`echo $l1`" = "" ] && l2=`find . -maxdepth 2 -mindepth 2 -name "OUTCAR.gz" | sed 's|^./||'` && geta="l2" && l=$l2
[ "`echo $l1`" = "" ] && [ "`echo $l2`" = "" ] && l3=`find . -maxdepth 3 -mindepth 3 -name "OUTCAR.gz" | sed 's|^./||'` && geta="l3" && l=$l3
[ "`echo $l`" = "" ] && error "no OUTCAR.*.gz or */OUTCAR.gz files existing." && exit

# check lattice constant of every OUTCAR
for i in $l;do
    [ "$geta" = "l1" ] && a=`echo $i | sed 's|.*OUTCAR\.\(.*\)\.gz|\1|'`
    [ "$geta" = "l2" ] && a=`echo $i | sed 's|Ang/OUTCAR.*||'`
    [ "$geta" = "l3" ] && a=`echo $i | sed 's|.*/\(.*\)Ang/OUTCAR.*|\1|'`
    [ "$a" = "" ] && error "lattice constant not found for $i" && exit
    echo i:$i geta:$geta: a:$a 
done
smearFirst=none; min=100000;

rm -f tmp energies_Tto0K free_energies log_extractEnergies; touch log_extractEnergies
for i in $l; do
  zgrep ".*" $i | grep -e ISMEAR -e NIONS -e "energy  w" -e "energy " -e "free  e" > _tmp_OUTCAR
  smear=`grep ISMEAR _tmp_OUTCAR | awk '{print $3}' | sed 's/\(.*\);/\1/'`
  if [ "$smear" != -1 ]; then str="  (no free_energies file)"; str2=" and"; else str=""; str2=", free_energies, and"; fi
  if [ $smearFirst == none ]; then smearFirst=$smear; echoFile; echoFile "smearing method $smear$str"; fi
  if [ "$smear" != "$smearFirst" ]; then error "smearing method inconsistent in $i"; fi

    [ "$geta" = "l1" ] && a=`echo $i | sed 's|.*OUTCAR\.\(.*\)\.gz|\1|'`
    [ "$geta" = "l2" ] && a=`echo $i | sed 's|Ang/OUTCAR.*||'`
    [ "$geta" = "l3" ] && a=`echo $i | sed 's|.*/\(.*\)Ang/OUTCAR.*|\1|'`

  NN=`awk '/NIONS/{print $NF}' _tmp_OUTCAR`
  e=`grep "energy  w" _tmp_OUTCAR | awk 'END{print $NF}'`
  d=`grep "energy w" _tmp_OUTCAR | tail -n2 | awk 'NR==1{e=$NF};{printf("%d",(($NF-e)^2)^(1/2)*'$eVTomeV')}'`
  if [ $d != 0 -a $d != 00 ]; then strd="WARNING:_last_delta=${d}_meV"; else strd=""; fi
  if [ $smear == '-1' ]; then f=`grep "free  e" _tmp_OUTCAR | awk '{print $(NF-1)}'`; fi
  min=`echo $e | awk '{if ($1<'$min') m=$1; else m='$min'; printf("%.6f",m)}'`
  E=`echo $e | awk '{printf("%.8f",$1/'$NN')}'`
  f=`echo $f | awk '{printf("%.8f",$1/'$NN')}'`
  echo $a $e $strd >> tmp
  echo $a $E >> energies_Tto0K; echo $a $f >> free_energies
done

echoFile; echoFile "aLat(Ang)  energyT->0K(Ry)    energyT->0K(eV)   shifted(meV)"
cat tmp | sort -g | awk '{printf("  %6.3f     %10.6f       %10.6f     %6.1f   %s   %s\n",
    $1,$2*'$eVToRyd'/'$NN', $2/'$NN',($2-(1*'$min'))*'$eVTomeV'/'$NN',$3,$4)}' >> log_extractEnergies
echoFile; echoFile "files written (eV): energies_Tto0K$str2 log_extractEnergies"
cat log_extractEnergies; rm _tmp_OUTCAR
cat energies_Tto0K | sort -g > tmp; mv tmp energies_Tto0K
cat free_energies | sort -g > tmp; mv tmp free_energies
if [ "$smear" != "-1" ]; then rm free_energies; fi

