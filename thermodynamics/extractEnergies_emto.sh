#!/bin/bash

# following 3 lines must always be present
path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
script=`set | grep BASH_SOURCE | sed 's|.*\".*/\([^/]*\)\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
options=$*; . $path/utilities/functions.include; checkOptions "-h -help -v -l"

# small help for options
# space between option and value must be 1
# space between option and Comment or value and Comment must be >2
h=`getOption -h`
if [ $h == True ]; then
  usage $script
  printOptions "-v   write energies_Tto0K_volume and free_energies_volume additionally" \
               "-l   extract also LDA energies"
  exit
fi

# detailed help
help=`getOption -help`
if [ $help == True ]; then
  details $script
  echo2 "   extract energies per atom from EMTO Murn calculation"
  echo2 "   results are assumed in OUT-EMTO.[0-9.]*{,.gz} files"
  echo2 "   T->0K energies and free energies are extracted"
  echo2 "   the -v option produces in addition to the energy files with aLats"  \
        "   energy files which use volumes; this is helpful if later using"     \
        "   getE0KFit.sh for fitting hcp runs"
  exit
fi

pi=3.14159265358979312; eVTomeV=1000; bohr3Toang3=0.148184709267640019; rydToeV=13.6056917292532837

echoFile () { echo "$1" >> log_extractEnergies
}

l=`ls OUT-EMTO.[0-9.]*[0-9.]{,.gz} 2> /dev/null | wc -l`
if [ $l == 0 ]; then error "no OUT-EMTO.[0-9.]*{,.gz} files existing"; fi

# check if -v option for additional energy files vs volume given
vOp=`getOption -v`

rm -f log_extractEnergies; touch log_extractEnergies

#check if -l option given
lOp=`getOption -l`
if [ $lOp == True ]; then xcs="LDA GGA"; else xcs="GGA"; fi

for xc in $xcs; do
  echoFile; echoFile "==== $xc ===="
  smearFirst=none; min=100000;
  l=`ls OUT-EMTO.[0-9.]*[0-9.]{,.gz} 2> /dev/null`;
  rm -f tmp energies_Tto0K_$xc free_energies_$xc
  if [ $vOp == True ]; then rm -f energies_Tto0K_volume_$xc free_energies_volume_$xc; fi
  for i in $l; do
    zgrep ".*" $i | grep -e ZMSH -e 'T*S' -e "TOT-$xc" -e Iteration -e TS_el > _tmp
    smear=`grep ZMSH _tmp | awk '{print $3}'`
    if [ "$smear" = "E" -o "$smear" = "e" ]; then str="  (no free_energies file)"; str2=" and"; else str=""; str2=", free_energies, and"; fi
    if [ $smearFirst == none ]; then smearFirst=$smear; echoFile; echoFile "smearing method $smear$str"; fi
    if [ "$smear" != "$smearFirst" ]; then error "smearing method inconsistent in $i"; fi
    e=`grep TOT-$xc _tmp | awk '{printf "%.6f",$4*'$rydToeV'}'`    # energy per site
    ts=`grep 'TS_el' _tmp | awk '{printf "%.6f",$4*'$rydToeV'}'` # T*S per site
    a=`echo $i | sed 's|OUT-EMTO\.\(.*\)|\1|' | sed 's|\.gz||'`
    d=`grep Iteration _tmp | awk 'END{printf("%d",$NF*'$rydToeV'*'$eVTomeV')}'`
    if [ $d != 0 -a $d != 00 ]; then strd="WARNING:_last_delta=${d}_meV"; else strd=""; fi
    if [ $smear = 'F' -o $smear = 'f' ]; then f=$e; e=`echo $f $ts | awk '{printf "%.6f",$1+0.5*($2)}'`; fi
    min=`echo $e | awk '{if ($1<'$min') m=$1; else m='$min'; printf("%.6f",m)}'`
    echo $a $e $strd >> _tmp_energies
    echo $a $e >> energies_Tto0K_$xc
    if [ $smear = 'F' -o $smear = 'f' ]; then echo $a $f >> free_energies_$xc; fi
    if [ $vOp == True ]; then
      v=`grep TOT-$xc _tmp | awk '{printf "%.5f", '$bohr3Toang3'*4*'$pi'*($(NF-1))^3/3}'` # RMT(bohr) to volume(ang^3)
      echo $v $e >> energies_Tto0K_volume_$xc
      if [ $smear = 'F' -o $smear = 'f' ]; then echo $v $f >> free_energies_volume_$xc; fi
    fi
  done

  echoFile; echoFile "vol(Ang^3)  energyT->0K(Ry)    energyT->0K(eV)   shifted(meV)"
  cat _tmp_energies | sort -g | awk '{printf("  %6.3f     %10.6f       %10.6f     %6.1f   %s   %s\n",
      $1,$2/'$rydToeV', $2,($2-(1*'$min'))*'$eVTomeV',$3,$4)}' >> log_extractEnergies
  echoFile; echoFile "files written (eV): energies_Tto0K$str2 log_extractEnergies"
  cat log_extractEnergies; rm _tmp_energies
  cat energies_Tto0K_$xc | sort -g > _tmp; mv _tmp energies_Tto0K_$xc
  if [ $smear = 'F' -o $smear = 'f' ]; then cat free_energies_$xc | sort -g > _tmp; mv _tmp free_energies_$xc; fi
done

