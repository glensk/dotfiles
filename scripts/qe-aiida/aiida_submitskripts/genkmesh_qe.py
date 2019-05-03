#!/usr/bin/env python
import ase
import ase.io
import numpy as np
import sys

def get_kmesh_size(structure, kmesh_l):
    reci_cell = structure.get_reciprocal_cell()
    kmesh = [np.ceil(kmesh_l * np.linalg.norm(reci_cell[i]))
             for i in range(len(reci_cell))]
    return kmesh

input_file = str(sys.argv[1])
kmesh_l = int(sys.argv[2])

structure = ase.io.read(input_file, format="espresso-in")
kmesh = get_kmesh_size(structure, kmesh_l)

print "{} {} {} 0 0 0".format(
        int(kmesh[0]),
        int(kmesh[1]),
        int(kmesh[2])
        )
