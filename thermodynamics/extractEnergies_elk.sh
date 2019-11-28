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
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   extract energies from elk Murn calculation"
  echo2 "   results are assumed in *Ang/TOTENERGY.OUT files"
  echo2 "   T->0K energies and free energies are extracted"
  exit
fi

hartreeToeV=27.2113838546484104
hartreeTomeV=27211.3838546484104

echoFile () { echo "$1" >> $dir/log_extractEnergies
}

l=`ls -1d *Ang/ 2> /dev/null | wc -l`
if [ $l == 0 ]; then error "no *Ang folders existing"; fi

smearFirst=none; min=100000
l=`ls -1d *Ang/`; dir=`pwd`
rm -f tmp energies_Tto0K free_energies log_extractEnergies; touch log_extractEnergies
for i in $l; do
  cd $i
  if [ ! -e TOTENERGY.OUT ]; then error "no TOTENERGY.OUT file existing in $i"; fi
  if [ ! -e INFO.OUT ]; then error "no INFO.OUT file existing in $i"; fi
  smear=`grep -A1 'Smearing type' INFO.OUT | awk 'END{print $1}'`
  if [ "$smear" == Gaussian -o "$smear" == 'Methfessel-Paxton' ]; then
    str="  (no free_energies file)"; str2=" and"; else str=""; str2=", free_energies, and";
  fi
  if [ $smearFirst == none ]; then smearFirst=$smear; echoFile; echoFile "smearing method $smear$str"; fi
  if [ "$smear" != "$smearFirst" ]; then error "smearing method inconsistent in $i"; fi
  a=`echo $i | sed 's|\(.*\)Ang/|\1|'`
  e=`cat TOTENERGY.OUT | awk 'END{print $NF}'`
  d=`cat TOTENERGY.OUT | tail -n2 | awk 'NR==1{e=$NF};{printf("%d",(($NF-e)^2)^(1/2)*'$hartreeTomeV')}'`
  if [ $d != 0 -a $d != 00 ]; then strd="WARNING:_last_delta=${d}_meV"; else strd=""; fi
  ts=`grep -e "electron entropic" INFO.OUT | awk 'END{print $NF}'`
  f=""
  if [ $smear == 'Fermi-Dirac' ]; then f=$e; e=`echo $f $ts | awk '{printf("%.6f",$1-(1*$2/2))}'`; fi
  min=`echo $e | awk '{if ($1<'$min') m=$1; else m='$min'; printf("%.6f",m)}'`
  E=`echo $e | awk '{printf("%.8f",$1*'$hartreeToeV')}'`
  f=`echo $f | awk '{printf("%.8f",$1*'$hartreeToeV')}'`
  echo $a $e $E $strd >> $dir/tmp
  echo $a $E >> $dir/energies_Tto0K; echo $a $f >> $dir/free_energies
  cd $dir
done

echoFile; echoFile "aLat(Ang)  energyT->0K(Ry)    energyT->0K(eV)   shifted(meV)"
awk '{printf("  %s     %s     %s     %6.1f   %s   %s\n",$1,$2,$3,($2-(1*'$min'))*'$hartreeTomeV',$4,$5)}' tmp >> log_extractEnergies
echoFile; echoFile "files written (eV): energies_Tto0K$str2 log_extractEnergies"
cat log_extractEnergies; rm tmp

