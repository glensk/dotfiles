#!/bin/bash

tol=0.0001 # Angstrom

if [ $# -gt 1 -o "$1" == "-h" ]; then
  echo "" 1>&2
  echo -e "\033[31m\033[1mUSAGE\033[0m:" 1>&2
  echo -e "   \033[1mgetDisplacements.sh\033[0m      optional input: tolerance (Ang)" 1>&2
  echo "" 1>&2
  echo "" 1>&2
  echo -e "   Note: - default for tolerance is $tol (Ang)" 1>&2
  exit 1
fi

path=`set | grep BASH_SOURCE | sed 's|.*\"\(.*\)/[^/]*\".*|\1|' | sed '/BASH_SOURCE=/s/.*/\./'`
error () { echo 1>&2; echo -e "\033[31m\033[1mERROR\033[0m: $1" 1>&2; exit 1;
}

# check if input files available
l="cell reduced_coords"
for i in $l; do if [ ! -e $i ]; then error "input file $i not existing"; fi; done

# check and import cell
check=`awk 'BEGIN{n=3};NF!=3{n=-3};END {if (n==3&&NR==3) print "ok"; else print "error"}' cell`;
if [ $check == "error" ]; then error "wrong # of rows or columns in cell file"; fi
A=`awk 'NR==1{print}' cell`; B=`awk 'NR==2{print}' cell`; C=`awk 'NR==3{print}' cell`;

# vector lengths and angles
pf=$path/fortran/
a=`$pf/norm.x $A`; b=`$pf/norm.x $B`; c=`$pf/norm.x $C`
aa=`$pf/dotProduct.x $B $C`; bb=`$pf/dotProduct.x $A $C`; cc=`$pf/dotProduct.x $A $B`;

check=`awk 'BEGIN{n=4};NF!=4{n=-4};END{if (n==-4) print "ok"}' reduced_coords`
if [ "$check" == "error" ]; then error "reduced_coords inconsistent"; fi

# construct sgroup input
echo P > sgroup.dat
echo "$a $b $c $aa $bb $cc" >> sgroup.dat
echo >> sgroup.dat
echo `awk 'END{print NR}' reduced_coords` >> sgroup.dat
awk '{printf("%s %s %s\n%s",$2,$3,$4,$1)}' reduced_coords >> sgroup.dat

# run sgroup
$path/sgroup/sgroup sgroup.dat > log_sgroup

