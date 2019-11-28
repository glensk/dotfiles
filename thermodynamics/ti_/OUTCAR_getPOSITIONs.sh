#!/bin/sh
## DIESES SKRIPT WIRD DIREKT IM JOBFOLDER AUSGEFUEHRT

#################################################################################
###### DIESER TEIL braucht noch die OUTCAR, alternativ! aus der vasprun.xml erstellen
#################################################################################

## braucht numer of atoms
## braucht
## braucht anzahl der steps (sollte man eher aus der structures{_vasprun}.gz rausholen

pfad=`XXXCAR_link-to-basename.sh $0 $1`

if [ -e "$pfad" ];then
### from OUTCAR
### from OUTCAR
nions=`OUTCAR_number_of_atoms.sh $pfad`
f=$pfad
zgrep -a --text ".*" $f | awk 'BEGIN{s=-1;l=-10};/NIONS =/{n=$12};/volume of cell/{v=$5;s++;l=NR};
                  s==0&&NR>l+1&&NR<l+5{print $1,$2,$3};END{print n,v,s}' > cell
l=`tail -n1 cell`  ## sowas wie 31 418.51 4  ## numatoms volumecell anzahlion steps
head -n3 cell > tmp; mv tmp cell
a=`echo $l | awk '{print $1+1}'`  ## numatoms+1
rm -f POSITIONs
echo $l > POSITIONs


else
### from vasprun
### from vasprun
[ "$1" = "" ] && pfadzurfile=""
[ "$1" != "" ] && [ -d = "$1" ] && pfadzurfile=$1/
[ "$1" != "" ] && [ -d != "$1" ] && pfadzurfile=`echo $1 | sed 's|OUTCAR.*||'`
[ "$pfadzurfile" = "" ] && pfadzurfile=.
      #echo pfadzurfile: $pfadzurfile
pfad=`find $pfadzurfile -maxdepth 1 -mindepth 1 -type f -name "vasprun.xml*"`
[ "`echo $pfad | wc -w | sed 's|[ ]*||g'`" != "1" ] && echo "`pwd` PROBLEM: fund 0 or more than 1 vasprun.xml file" && exit
nions=`zgrep -a --text "<atoms>" $pfad | sed 's|<atoms>||' | sed 's|</atoms>||' | sed 's|^[ \t]*||;s|[ \t]*$||'`
volume=`zgrep -a --text '"volume"' $pfad | head -1 | sed 's|.*">||' | sed 's|<.*||' | sed 's|^[ \t]*||;s|[ \t]*$||'`
cell=`zgrep -a --text -A 3 'name="basis"' $pfad | head -4 | tail -3 | sed 's|<v>||' | sed 's|</v>||' | sed 's|^[ \t]*||;s|[ \t]*$||'`
rm -f cell; echo "$cell" > cell
rm -f POSITIONs
echo $nions $volume 0 > POSITIONs
fi
#echo 1;cat POSITIONs

#################################################################################
###### DIESER TEIL WIRD AUS DER structures{_vasprun}.gz ERSTELLT
#################################################################################
structures=structures.gz
[ ! -e "$structures" ] && structures=structures_vasprun.gz; [ ! -e "$structures" ] && echo "`pwd` structures not found! " && exit
if [ -e "$structures" ];then
zgrep -a --text ".*" $structures | awk '{print $2,$3,$4}' >> POSITIONs 
lines=`zgrep -a --text ".*" $structures | wc -l | sed 's|[ ]*||g' | awk '{print $1}'`
steps=` echo "$lines/$nions" | bc`
#echo steps: $steps
#echo nions: $nions
#echo lines: $lines
sed -i '1 s|\([0-9]*\)\([ ]*\)\([0-9.]*\).*|\1 \3 '"$steps"'|' POSITIONs
exit
fi


## if both are not available
#zgrep -a --text -A$a POSITION $f | grep -v -e -- -e POSITION | awk '{print $1,$2,$3}' >> POSITIONs
echo PROBLEM: no files with structures found!
rm -f POSITIONs

