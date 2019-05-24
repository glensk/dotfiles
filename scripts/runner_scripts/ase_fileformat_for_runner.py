from __future__ import print_function

import re
import numpy as np
import sys

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
        # info = key_val_str_to_dict(line)
        #print('9999999 88888',line)
        #sys.exit('876999')

        #info = {}
        info = {'comment': line.strip()}
        arrays = {}

        cell = [] # Now we have to read 3 lines to extract information about the lattice

        for ln in range(3):
            line = fileobj.readline()
            vals = line.split()
            cell.append(vals[1:4])
        cell = np.array(cell,dtype='float')*units.Bohr # This converts to Angstrom, which is the default unit in ASE

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
        energy = float(line.split()[1])*units.Hartree
        line = fileobj.readline()
        charge = float(line.split()[1])

        arrays['forces'] = forces
        info['nrg'] = energy

        structure = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True,info=info)

        calc = SinglePointCalculator(structure, energy=energy, forces=forces)
        structure.set_calculator(calc)


        #atoms = Atoms(symbols=symbols,
        #              positions=positions,
        #              cell=cell,
        #              pbc=pbc,
        #              info=info)

        yield structure

def write_runner(fileobj,images,comment=None,append=False,setenergy_eV=False,setforces_ase_units=False,wrap=True):
    """
    Write output in runner format. Written quickly with no regard to the form
    """

    if isinstance(fileobj, basestring):
        mode = 'w'
        if append:
            mode = 'a'
        fileobj = open(fileobj,mode)

    if hasattr(images, 'get_positions'):
        images = [images]

    for idximage,atoms in enumerate(images):
        nat = atoms.get_number_of_atoms()
        #print('idximage',idximage,nat,setforces_ase_units,setenergy_eV)

        ####################################################
        # forces
        ####################################################
        if type(setforces_ase_units) == bool:
            if setforces_ase_units == False:
                try: forces = atoms.get_forces()/units.Hartree * units.Bohr
                except:
                    sys.stderr.write("No forces found, setting them to 0\n")
                    forces = np.zeros((nat,3))
        else:
            if type(setforces_ase_units) == str:
                forces = np.zeros((nat,3))
            else:
                forces = setforces_ase_units/units.Hartree * units.Bohr

        ####################################################
        # energy
        ####################################################
        if setenergy_eV == False:
            try: energy = atoms.get_potential_energy()/units.Hartree
            except:
                sys.stderr.write("No energy found, setting it to 0\n")
                energy = 0
        else:
            #print('uh',units.Hartree)
            energy = setenergy_eV/units.Hartree

        #atoms.wrap()
        fileobj.write('begin\n')
        if comment is None:
            #print('!-===1',atoms.info)
            #print('!-===2',atoms.info.get('comment'))
            #print('!-===2',atoms.info.get('comment').split())
            #sys.exit()
            if 'comment' in atoms.info:
                listcheck = atoms.info.get('comment').split()
                if len(listcheck) > 1 and listcheck[0] == 'comment':
                    fileobj.write(atoms.info.get('comment') + "\n")
                else:
                    fileobj.write('comment '+atoms.info.get('comment') + "\n")
            else:
                fileobj.write('comment ' + str(atoms.get_number_of_atoms())  + ' atoms, species ' + str(atoms.get_chemical_formula()) + "\n")
        else:
            fileobj.write('comment ' + comment)

        cell = atoms.get_cell()/units.Bohr
        for idx in range(3):
            #fileobj.write("lattice " + str(cell[idx][0]) + " " + str(cell[idx][1]) + " " + str(cell[idx][2]) + "\n")
            fileobj.write('lattice  %16.10f %16.10f %16.10f\n' % (cell[idx][0],cell[idx][1],cell[idx][2]))

        nat = atoms.get_number_of_atoms()
        positions = atoms.get_positions(wrap=False)/units.Bohr  # Wrapping takes huge amout of time
        for idx in range(nat):
            fileobj.write('atom  %16.10f %16.10f %16.10f  %3s  %1.1f %1.1f   %15.8e %15.8e %15.8e\n' % (positions[idx,0],positions[idx,1],positions[idx,2],atoms[idx].symbol,0.0,0.0,forces[idx][0],forces[idx][1],forces[idx][2]))



        fileobj.write('energy %.15e\n' % energy)
        fileobj.write('charge 0.0000\n')
        fileobj.write('end\n')
