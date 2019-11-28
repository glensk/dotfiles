#!/bin/sh

if [ $# -lt 1 ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e " highClean.sh -a                (to run on all *Ang_*K folders)" 1>&2
  echo -e " highClean.sh *Ang*K_folders    (to run on given folders)" 1>&2
  exit 1
fi
if [ "$1" == "-a" ]; then l=`ls -1d *Ang_*K`; else l="$*"; fi

echo; dir=`pwd`; c=0; s=0
for i in $l; do
 if [ -d "$i" ]; then
  echo -n $i; cd $i
  ll=`ls -1d lambda*/*_*/`
  for ii in $ll; do
    n=`ls -1 $ii/ | wc -l | sed 's|[ ]*||g'`
    p=`ls -1 $ii/`
    if [ "$n" == 1 -a "$p" == "POSCAR" ]; then 
      c=`expr $c + 1`; rm -fr $ii
    fi
  done
  echo " $c removed"; s=`expr $s + $c`; c=0; cd $dir
 else
  echo; echo no directory $i; echo
 fi
done
echo "sum $s removed"
