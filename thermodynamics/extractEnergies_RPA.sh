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
  echo2 "   extract energies per atom from vasp RPA Murn calculation"
  echo2 "   results are assumed in *Ang/{01_DFT,02_EXX,03_diag,04_ACFDT,05_HFc}/OUTCAR{,.gz} files"
  echo2 "   T->0K energies are only extracted"
  echo2 "   EATOM is not included!"
  exit
fi


l=`ls -1d *Ang/`

eVTomeV=1000; eVToRyd=0.0734986508052261761; min=100000;

echoFile () { echo "$1" >> log_extractEnergies
}

rm -f tmp energies_Tto0K* log_extractEnergies exchange_HF; touch log_extractEnergies
for i in $l; do
  a=`echo $i | sed 's|Ang.*||'`
  NN=`zcat -f $i/01_DFT/OUTCAR   2> /dev/null | awk '/NIONS/{print $NF}'`
  e0=`zcat -f $i/01_DFT/OUTCAR   2> /dev/null | awk '/energy  w/{print $NF}'`
  e1=`zcat -f $i/02_EXX/OUTCAR   2> /dev/null | awk '/energy  w/{print $NF}'`
  eHF=`zcat -f $i/02_EXX/OUTCAR  2> /dev/null | awk '/-exchange      EXHF   =/{print $NF}'`
  eno=`zcat -f $i/00_noXC/OUTCAR 2> /dev/null | awk '/energy  w/{print $NF}'`
  ano=`zcat -f $i/00_noXC/OUTCAR 2> /dev/null | awk '/atomic energy  EATOM  =/{print $NF}' | tail -n1`
  aHF=`zcat -f $i/02_EXX/OUTCAR  2> /dev/null | awk '/atomic energy  EATOM  =/{print $NF}' | tail -n1`
  e2=`zcat -f $i/04_ACFDT/OUTCAR 2> /dev/null | awk '/converged v/{print $(NF-1)}'`
  e3=`zcat -f $i/05_HFc/OUTCAR   2> /dev/null | awk '/HF-corr/{print $NF}'`

  str=""
  if [ `echo $eno | awk '{print NF}'` != 1 ]; then str=$str"__problemIn00_noXC"; fi
  if [ `echo $NN  | awk '{print NF}'` != 1 ]; then str=$str"__problemIn01_DFT"; fi
  if [ `echo $e1  | awk '{print NF}'` != 1 ]; then str=$str"__problemIn02_EXX"; fi
  if [ `echo $e2  | awk '{print NF}'` != 1 ]; then str=$str"__problemIn04_ACFDT"; fi
  if [ `echo $e3  | awk '{print NF}'` != 1 ]; then str=$str"__problemIn05_HFc"; fi

  e=`echo $e1 $e2 $e3 | awk '{printf "%12.6f",$1+1.*($2)+1.*($3)}'`
  min=`echo $e | awk '{if ($1<'$min') m=$1; else m='$min'; printf("%.6f",m)}'`
  E=`echo $e | awk '{printf("%.8f",$1/'$NN')}'`
  E0=`echo $e0 | awk '{printf("%.8f",$1/'$NN')}'`
  EEXX=`echo $e1 | awk '{printf("%.8f",$1/'$NN')}'`
  EHF=`echo $eHF | awk '{printf("%.8f",$1/'$NN')}'`
  ECORR=`echo $e2 | awk '{printf("%.8f",$1/'$NN')}'`
  hf=`echo $e1 $aHF $eno $ano | awk '{printf("%.8f",(($1-$2)-($3-$4))/'$NN')}'`
  echo $a $e $str >> tmp
  echo $a $E >> energies_Tto0K
  echo $a $E0 >> energies_Tto0K_DFT
  echo $a $EHF >> energies_Tto0K_HF
  echo $a $EEXX >> energies_Tto0K_EXX
  echo $a $ECORR >> energies_Tto0K_CORR
  echo $a $hf >> exchange_HF
done

echoFile; echoFile "aLat(Ang)  energyT->0K(Ry)    energyT->0K(eV)   shifted(meV)"
cat tmp | sort -g | awk '{printf("  %6.3f     %10.6f       %10.6f     %6.1f   %s   %s\n",
    $1,$2*'$eVToRyd'/'$NN', $2/'$NN',($2-(1*'$min'))*'$eVTomeV'/'$NN',$3,$4)}' >> log_extractEnergies
echoFile; echoFile "files written (eV): energies_Tto0K$str2 log_extractEnergies"
cat log_extractEnergies;
cat energies_Tto0K | sort -g > tmp; mv tmp energies_Tto0K
cat energies_Tto0K_DFT | sort -g > tmp; mv tmp energies_Tto0K_DFT
cat energies_Tto0K_HF | sort -g > tmp; mv tmp energies_Tto0K_HF
cat energies_Tto0K_EXX | sort -g > tmp; mv tmp energies_Tto0K_EXX
cat energies_Tto0K_CORR | sort -g > tmp; mv tmp energies_Tto0K_CORR
cat exchange_HF | sort -g > tmp; mv tmp exchange_HF

paste energies_Tto0K_EXX energies_Tto0K_DFT | awk '{printf "%s  %.8f\n",$1,$2-$4}' > energies_Tto0K_EXX-DFT


