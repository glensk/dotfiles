from __future__ import print_function

import re
import numpy as np

from ase.atoms import Atoms
from ase.calculators.calculator import all_properties, Calculator
from ase.calculators.singlepoint import SinglePointCalculator
from ase.parallel import paropen
from ase.utils import basestring
from ase import units

def read_runner(fileobj, index=-1):
    """
    Read from a file in RuNNer format

    index is the frame to read, default is last frame (index=-1).
    """
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    if not isinstance(index, int) and not isinstance(index, slice):
        raise TypeError('Index argument is neither slice nor integer!')

    # If possible, build a partial index up to the last frame required
    last_frame = None
    if isinstance(index, int) and index >= 0:
        last_frame = index
    elif isinstance(index, slice):
        if index.stop is not None and index.stop >= 0:
            last_frame = index.stop

    # scan through file to find where the frames start
    fileobj.seek(0)
    frames = []
    natoms = 0
    while fileobj:
        line = fileobj.readline()
        if line.strip() == '':
            break
        elif line.strip() == 'begin':
            frame_pos = fileobj.tell()
        elif line.split()[0] == 'atom':
            natoms += 1
        elif line.strip() == 'end':
            frames.append((frame_pos, natoms))
            natoms = 0
            if last_frame is not None and len(frames) > last_frame:
                break

    if isinstance(index, int):
        if index < 0:
            tmpsnp = len(frames) + index
            trbl = range(tmpsnp, tmpsnp + 1, 1)
        else:
            trbl = range(index, index + 1, 1)
    elif isinstance(index, slice):
        start = index.start
        stop = index.stop
        step = index.step

        if start is None:
            start = 0
        elif start < 0:
            start = len(frames) + start

        if step is None:
            step = 1

        if stop is None:
            stop = len(frames)
        elif stop < 0:
            stop = len(frames) + stop

        trbl = range(start, stop, step)
        if step < 0:
            trbl.reverse()

    for index in trbl:
        frame_pos, natoms = frames[index]
        fileobj.seek(frame_pos)

        # comment line
        line = fileobj.readline() 
        # Is there any valuable info to strip from the comment line?
        # Is the comment line always present?
        # info = key_val_str_to_dict(line)

        info = {}
        arrays = {}

        cell = [] # Now we have to read 3 lines to extract information about the lattice

        for ln in range(3):
            line = fileobj.readline()
            vals = line.split()
            cell.append(vals[1:4])
        cell = np.array(cell,dtype='float')*units.Bohr # This should convert to Angstrom

        positions = []
        symbols = []
        forces = []
        numbers = []
        for ln in range(natoms):

            line = fileobj.readline()
            vals = line.split()
            coords = tuple([float(val) for val in vals[1:4]])
            tmpforces = tuple([float(val) for val in vals[7:]])

            symbols.append(vals[4])
            positions.append(coords)
            forces.append(tmpforces)


        try:
            positions = np.array(positions)*units.Bohr
            forces = np.array(forces)*units.Hartree/units.Bohr # Transform from Ha/bohr to eV/angstrom

        except TypeError:
            raise IOError('Badly formatted data, ' +
                          'or end of file reached before end of frame')


        line = fileobj.readline()
        energy = float(line.split()[1])
        line = fileobj.readline()
        charge = float(line.split()[1])

        arrays['forces'] = forces
        info['nrg'] = energy

        structure = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)

        calc = SinglePointCalculator(structure, energy=energy, forces=forces)
        structure.set_calculator(calc)


        #atoms = Atoms(symbols=symbols,
        #              positions=positions,
        #              cell=cell,
        #              pbc=pbc,
        #              info=info)

        yield structure
