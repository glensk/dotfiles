#!/bin/sh
atoms=`head -100 aiida.out | grep "number of atoms/cell" | awk '{print $5}'`
echo $atoms

