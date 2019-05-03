#!/bin/sh
atoms=`head -100 $1 | grep "number of atoms/cell" | awk '{print $5}'`

#line=`grep -n "^     Forces acting on atoms (cartesian axes, Ry/au):" aiida.out | tail -1 | sed 's|:.*||' | awk '{print $1+2}'`
#lineend=`echo $line $atoms | awk '{print $1+$2-1}'`
##echo atoms: $atoms, $line $lineend
##sed -n ''"$line"','"$lineend"'p' aiida.out | awk '{print $7,$8,$9}'
##sed -n ''"$line"','"$lineend"'p' aiida.out | awk '{print $7*25.711032,$8*25.711032,$9*25.711032}'  # in eV
##sed -n ''"$line"','"$lineend"'p' aiida.out | awk '{print sqrt(($7*25.711032)^2),sqrt(($8*25.711032)^2),sqrt(($9*25.711032)^2)}' 
#sed -n ''"$line"','"$lineend"'p' aiida.out | awk '{print sqrt(($7)^2),sqrt(($8)^2),sqrt(($9)^2)}' 

tac $1 | grep -m 1 "ATOMIC_POSITIONS (angstrom)" -B $atoms | tac | tail -$atoms

