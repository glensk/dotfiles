#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -n"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-n   do not clean"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   extract energies from wien2k Murn or Rmt-convergence calculation"
  echo2 "   results are assumed in *Ang/base/base.scf files or in"             \
        "   *RmtKmax_Rmt=*/base/base.scf in case of Rmt-convergence"
  echo2 "   T->0K energies and free energies are extracted"
  echo2 "   by default the wien2k output is cleaned using clean_wien.sh (only" \
        "   important files are kept) to prevent cleaning use the -n option"
  exit
fi

rydbergToeV=13.6056919282781491
rydbergTomeV=13605.6919282781491

echoFile () { echo "$1" >> $dir/log_extractEnergies
}

l=`ls -1d [0-9.]*Ang/ 2> /dev/null | wc -l`
l2=`ls -1d Rmt=[0-9.]*_[0-9.]*Ang/ 2> /dev/null | wc -l`
if [ $l == 0 -a $l2 == 0 ]; then error "no [0-9.]*Ang/ nor Rmt=[0-9.]*_[0-9.]*Ang/ folders existing"; fi

# check if -n option is given (no cleaning)
noClean=`getOption -n`

smearFirst=none; min=100000
if [ $l == 0 ]; then
  l=`ls -1d Rmt=[0-9.]*_[0-9.]*Ang/`
  mode=Rmt
else
  l=`ls -1d [0-9.]*Ang/`;
  mode=Ang
fi
dir=`pwd`; ww="";
rm -f tmp energies_Tto0K* free_energies* delta_energies_Tto0K delta_free_energies log_extractEnergies; touch log_extractEnergies
for i in $l; do
  cd $i
  ll=`ls -1d */ 2> /dev/null | wc -l`
  if [ $ll != 1 ]; then error "no subfolder or too many in folder $i"; fi
  base=`ls -1d */ | sed 's|/||'`
  cd $base/
  smear=`awk 'NR==3{print $1}' $base.in2`
  if [ "$smear" != TETRA -a "$smear" != GAUSS -a "$smear" != TEMP -a "$smear" != TEMPS ]; then
    error "smearing method $smear not supported";
  fi
  if [ \( "$smear" == TETRA -o "$smear" == GAUSS \) -a $mode == True ]; then
      str="  (no free_energies file)";  str2=" and"; else str=""; str2=", free_energies, and";
  fi
  if [ \( "$smear" == TETRA -o "$smear" == GAUSS \) -a $mode != True ]; then
      str="  (no free_energies_*RmtKmax files)";  str2=" and"; else str=""; str2=", free_energies_*RmtKmax, and";
  fi
  if [ $smearFirst == none ]; then smearFirst=$smear; echoFile; echoFile "smearing method $smear$str"; fi
  if [ "$smear" != "$smearFirst" ]; then error "smearing method inconsistent in $i"; fi
  if [ $mode == Ang ]; then
    a=`echo $i | sed 's|\(.*\)Ang/|\1|' | awk '{printf("%.2f",$1)}'`
  else
    rmt=`echo $i | sed 's|Rmt=\(.*\)_\(.*\)Ang/|\1|' | awk '{printf("%.3f",$1)}'`
    a=`  echo $i | sed 's|Rmt=\(.*\)_\(.*\)Ang/|\2|' | awk '{printf("%.2f",$1)}'`
  fi
  e=`grep :ENE $base.scf | awk 'END{print $NF}'`
  d=`grep :ENE $base.scf | tail -n2 | awk 'NR==1{e=$NF};{printf("%d",(($NF-e)^2)^(1/2)*'$rydbergTomeV')}'`
  if [ $d != 0 -a $d != 00 ]; then strd="WARNING:_last_delta=${d}_meV"; else strd=""; fi
  warn=`grep :ENE $base.scf | tail -n1 | awk '$3=="*WARNING**" {print $3}'`
  
  # clean only if -n option not given and if no warnings
  if [ $noClean != True -a "$strd" == "" -a "$warn" == "" ]; then $path/clean_wien.sh > /dev/null; fi
  if [ "$strd" != "" -o "$warn" != "" ]; then ww=WARNINGS; fi

  # if TEMP we get -(T*S)/2; uf TEMPS we get -(T*S)
  ts=`grep -e "-(T\*S)" $base.scf | awk 'END{print $NF}'`
  f=""
  if [ $smear == TEMP  ]; then f=`echo $e $ts | awk '{printf("%.6f",$1+(1*$2))}'`; fi
  if [ $smear == TEMPS ]; then f=$e; e=`echo $f $ts | awk '{printf("%.6f",$1-(1*$2/2))}'`; fi
  min=`echo $e | awk '{if ($1<'$min') m=$1; else m='$min'; printf("%.6f",m)}'`
  E=`echo $e | awk '{printf("%.8f",$1*'$rydbergToeV')}'`
  f=`echo $f | awk '{printf("%.8f",$1*'$rydbergToeV')}'`
  charge=`grep :NEC $base.scf | tail -n3 | awk 'BEGIN{delta=0};
    {deltaNew=$(NF-2)-$(NF-1); if (deltaNew>delta) delta=deltaNew;} END {printf("%10.5f",delta)}'`
  if [ $mode == Ang ]; then
    echo $a $e $E $charge $warn $strd >> $dir/tmp
    echo $a $E >> $dir/energies_Tto0K; echo $a $f >> $dir/free_energies
  else
    echo $rmt $E $charge >> $dir/energies_Tto0K_$a; echo $rmt $f $charge >> $dir/free_energies_$a
  fi
  cd $dir
done

echoFile;
if [ $mode == Ang ]; then
  echoFile "aLat(Ang)  energyT->0K(Ry)    energyT->0K(eV)   shifted(meV)  MissingCharge"
  cat tmp | sort -g | awk '{printf("  %s     %s     %s     %6.1f    %11s   %s   %s\n",$1,$2,$3,($2-(1*'$min'))*'$rydbergTomeV',$4,$5,$6)}' >> log_extractEnergies
  echoFile; echoFile "files written (eV): energies_Tto0K$str2 log_extractEnergies"
  if [ "$ww" == WARNINGS ]; then echoFile; echoFile "NOTE: output left for rerun due to warnings"; fi
else
  paste energies_Tto0K_* | 
    awk 'NR==1{for (i=1;i<=NF/3-1;i++) f[i]=($(3*i-1)-1*($((i+1)*3-1)))*1000}
         {printf("%s  ",$1); for (i=1;i<=NF/3-1;i++) printf("%6.2f ",($(3*i-1)-1*($((i+1)*3-1)))*1000-f[i])
          delta=0; for (i=1;i<=NF/3-1;i++) if ($(3*i)>delta) delta=$(3*i); printf ("  %17.5f\n",delta)}' > delta_energies_Tto0K
  paste free_energies_* | 
    awk 'NR==1{for (i=1;i<=NF/3-1;i++) f[i]=($(3*i-1)-1*($((i+1)*3-1)))*1000}
         {printf("%s  ",$1); for (i=1;i<=NF/3-1;i++) printf("%6.2f ",($(3*i-1)-1*($((i+1)*3-1)))*1000-f[i])
          delta=0; for (i=1;i<=NF/3-1;i++) if ($(3*i)>delta) delta=$(3*i); printf ("  %17.5f\n",delta)}' > delta_free_energies
  echoFile; echoFile "delta_energies_Tto0K:"; echoFile; 
  awk 'NR==1{print "  Rmt(Bohr) Delta(s)(meV)    MissingCharge"};{print "    "$0}' delta_energies_Tto0K >> log_extractEnergies
fi
cat log_extractEnergies
l=`ls energies_Tto0K* free_energies*`
for i in $l; do cat $i | sort -g > tmp; mv tmp $i; done


