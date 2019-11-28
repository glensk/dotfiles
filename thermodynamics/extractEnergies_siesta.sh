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
  echo2 "   extract energies per atom from siesta murn calculation"
  echo2 "   results are assumed in log.aLat.gz files"
  echo2 "   T->0K energies and free energies are extracted"
  echo2 "   typical script sequence:"                     \
        "     createFolders_murn.sh -a ..."               \
        "     getALatRangeMurn.sh -e ..."                 \
        "     createFolders_murn.sh"                      \
        "     \033[1mextractEnergies_siesta.sh\033[0m"    \
        "     getE0KFit.sh           (cmpc01)"            \
        "     getThermodynamics.sh   (needs Fqh additionally)"
  exit
fi

eVTomeV=1000; eVToRyd=0.0734986508052261761

echoFile () { echo "$1" >> log_extractEnergies
}

l=`ls log.*.gz 2> /dev/null | wc -l`
if [ $l == 0 ]; then error "no log.*.gz files existing"; fi

smearFirst=none; min=100000;
l=`ls log.*.gz`;
rm -f tmp energies_Tto0K free_energies log_extractEnergies; touch log_extractEnergies
for i in $l; do
  zgrep ".*" $i | sed -n -e '/siesta: iscf/,/^$/p' -e '/OccupationFunction/p' -e '/Number of atoms/p' -e '/Etot    =/p' -e '/FreeEng =/p' > _tmp_log
  smear=`grep OccupationFunction _tmp_log | awk '{print $NF}'`
  if [ "$smear" != FD ]; then str="  (no free_energies file)"; str2=" and"; else str=""; str2=", free_energies, and"; fi
  if [ $smearFirst == none ]; then smearFirst=$smear; echoFile; echoFile "smearing method $smear$str"; fi
  if [ "$smear" != "$smearFirst" ]; then error "smearing method inconsistent in $i"; fi
  a=`echo $i | sed 's|log\.\(.*\)\.gz|\1|'`
  NN=`awk '/initatomlists: Number of atoms/{at=$(NF-2)}; END{print at}' _tmp_log`
  e=`grep "Etot" _tmp_log | awk 'END{print $NF}'`
  d=`sed -n -e '/siesta: iscf/,/^$/p' _tmp_log | awk 'NF==7{eold=e; e=$4};END{printf("%d",((e-eold)^2)^(1/2)*'$eVTomeV')}'`
  if [ $d != 0 -a $d != "00" ]; then strd="WARNING:_last_delta=${d}_meV"; else strd=""; fi
  if [ $smear == 'FD' ]; then f=`grep "FreeEng" _tmp_log | awk 'END{print $NF}'`; fi
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
cat log_extractEnergies; rm _tmp_log
cat energies_Tto0K | sort -g > tmp; mv tmp energies_Tto0K
cat free_energies | sort -g > tmp; mv tmp free_energies
if [ "$smear" != FD ]; then rm free_energies; fi

