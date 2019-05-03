#!/bin/sh
atoms=`head -100 $1 | grep "number of atoms/cell" | awk '{print $5}'`
echo $atoms

