#!/bin/bash

path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
file=$path/utilities/atomic_weights_and_ionic_compositions_NIST

error () { echo 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: $1" 1>&2; exit 1;
}

if [ ! -e $file ]; then error "reference file not available:\n $file"; fi

if [ $# != 1 -o "$1" == "-h" ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e "   \033[1mgetAtomicMass.sh ElementName\033[0m" 1>&2
  echo "" 1>&2
  echo "" 1>&2
  echo -e "   Notes: - atomic weights are extracted from file:" 1>&2
  echo -e "            `ls $file`" 1>&2
  echo "" 1>&2
  echo -e "          - `grep Reference $file`" 1>&2
  exit 1
fi

elem=$1
m=`grep -i -A 4 "^Atomic Symbol = $elem" ~/Thermodynamics/utilities/atomic_weights_and_ionic_compositions_NIST | grep "Standard Atomic Weight =" | awk '{print $5}' | sort | uniq`

#m=`awk 'BEGIN{l=-10;m=0};
#        /Atomic Symbol =/{if ($NF=="'$elem'") l=NR};
#        l+4==NR{if (m!=0&&m!=$NF) m=-1; else m=$NF; l=-10};
#        END{print m}' $file`

if [ $m == 0 ];  then error "element $elem not existing in database"; fi
if [ $m == -1 ]; then error "isotopes have different mass"; fi

# removing last digit in round brackets (inaccuracy interval)
# and square brackets which are present for some elements
echo $m | sed -e 's/(.*)//' -e 's/\[//' -e 's/\]//'

