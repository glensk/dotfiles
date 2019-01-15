"""The functions below are for reference only.
We use the implementation from extxyz module, which is backwards
compatible with standard XYZ format."""
import sys
from ase.atoms import Atoms
from ase.io.extxyz import read_extxyz as read_xyz, write_extxyz as write_xyz

__all__ = ['read_ipi', 'write_ipi']


#def simple_read_xyz(fileobj, index):
def read_ipi(fileobj, index):
    lines = fileobj.readlines()
    natoms = int(lines[0])
    cell_str = lines[1]
    #print('nat',natoms)
    #print('cel',cell_str)
    #print('cel',type(cell_str))
    lst = cell_str.split()[2:8]
    cell_out = [float(i) for i in lst]
    #print('out',out)
    #print('cell_out:')
    #print(cell_out)
    pbc = (True, True, True)
    #print('xxx',cell_str[2])

    #sys.exit()
    nimages = len(lines) // (natoms + 2)
    for i in range(*index.indices(nimages)):
        symbols = []
        positions = []
        n = i * (natoms + 2) + 2
        for line in lines[n:n + natoms]:
            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
        #atoms[0].set_cell(atoms[0].get_cell_lengths_and_angles()))
        yield Atoms(symbols=symbols, positions=positions, cell=cell_out, pbc=pbc)


#def simple_write_xyz(fileobj, images, comment=''):
def write_ipi(fileobj, images, comment=''):
    symbols = images[0].get_chemical_symbols()
    laa = images[0].get_cell_lengths_and_angles()
    comment = '# CELL(abcABC):   '+str(laa[0])+"  "+str(laa[1])+"  "+str(laa[2])+"  "+str(round(laa[3]))+"  "+str(round(laa[4]))+"  "+str(round(laa[5]))+" Step: 4  Bead: 0 positions{angstrom} cell{angstrom}"

    natoms = len(symbols)
    for atoms in images:
        fileobj.write('%d\n%s\n' % (natoms, comment))
        for s, (x, y, z) in zip(symbols, atoms.positions):
            fileobj.write('%-2s %16.8f %16.8f %16.8f\n' % (s, x, y, z))
