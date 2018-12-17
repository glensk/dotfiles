#!/usr/bin/env python
import myutils as my
from ase.io import read as ase_read

atomsa = ase_read('out.runner',index=':',format='runner')
print('ta :',type(atomsa),len(atomsa))
atomsb = ase_read('out.runner.1',index=':',format='runner')
print('ta 1',type(atomsb),len(atomsb))
print()
print('----------------')
for i in atomsa:
    my.lammps_ext_calc(i)
#print('----------------')
#for i in atomsb:
#    my.lammps_ext_calc(i)
#
