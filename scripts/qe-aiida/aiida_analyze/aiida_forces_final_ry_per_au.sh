#!/bin/sh

#You have: rydberg/bohrradius
#You want: eV/angstrom
#  * 25.711032
#  / 0.038893811
#You have:

atoms=`head -100 $1 | grep "number of atoms/cell" | awk '{print $5}'`
linebegin=`grep -n "^     Forces acting on atoms (cartesian axes, Ry/au):" $1 | tail -1 | sed 's|:.*||' | awk '{print $1+2}'`
lineend=`echo $linebegin $atoms | awk '{print $1+$2-1}'`
#echo atoms: $atoms, $line $lineend
sed -n ''"$linebegin"','"$lineend"'p' $1 | awk '{print $7,$8,$9}'
#sed -n ''"$linebegin"','"$lineend"'p' aiida.out | awk '{print $7*25.711032,$8*25.711032,$9*25.711032}'  # in eV
#sed -n ''"$linebegin"','"$lineend"'p' aiida.out | awk '{print sqrt(($7*25.711032)^2),sqrt(($8*25.711032)^2),sqrt(($9*25.711032)^2)}' 


#sed -n ''"$linebegin"','"$lineend"'p' $1 | awk '{print sqrt(($7)^2),sqrt(($8)^2),sqrt(($9)^2)}' 
