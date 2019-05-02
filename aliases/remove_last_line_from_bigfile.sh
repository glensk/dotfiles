#!/bin/sh
filename="lammps_pos.txt"
# tail -n1 log.lammps | wc -c  === 36
# stat --format=%s log.lammps === 14237564928
# seek === 14237564892
# dd --help  === Copy a file, converting and formatting according to the operands.
dd if=/dev/null of=$filename bs=1 seek=$(echo $(stat --format=%s $filename ) - $( tail -n1 $filename | wc -c) | bc )
