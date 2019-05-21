"""The functions below are for reference only.
We use the implementation from extxyz module, which is backwards
compatible with standard XYZ format."""
import sys
from ase.atoms import Atoms
from ase.io.extxyz import read_extxyz as read_xyz, write_extxyz as write_xyz
from numpy.linalg import norm
import numpy as np

__all__ = ['read_ipi', 'write_ipi']

def is_upper_triangular(arr, atol=1e-8):
    """test for upper triangular matrix based on numpy"""
    # must be (n x n) matrix
    assert len(arr.shape)==2
    assert arr.shape[0] == arr.shape[1]
    return np.allclose(np.tril(arr, k=-1), 0., atol=atol)

def convert_cell(cell,pos):
    """
    Convert a parallelepipedal (forming right hand basis)
    to lower triangular matrix LAMMPS can accept. This
    function transposes cell matrix so the bases are column vectors
    """
    cell = np.matrix.transpose(cell)

    if not is_upper_triangular(cell):
        # rotate bases into triangular matrix
        tri_mat = np.zeros((3, 3))
        A = cell[:, 0]
        B = cell[:, 1]
        C = cell[:, 2]
        tri_mat[0, 0] = norm(A)
        Ahat = A / norm(A)
        AxBhat = np.cross(A, B) / norm(np.cross(A, B))
        tri_mat[0, 1] = np.dot(B, Ahat)
        tri_mat[1, 1] = norm(np.cross(Ahat, B))
        tri_mat[0, 2] = np.dot(C, Ahat)
        tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
        tri_mat[2, 2] = norm(np.dot(C, AxBhat))

        # create and save the transformation for coordinates
        volume = np.linalg.det(cell)
        trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
        trans /= volume
        coord_transform = np.dot(tri_mat, trans)
        pos = np.dot(coord_transform, pos.transpose())
        pos = pos.transpose()

        return tri_mat, pos
    else:
        return cell, pos

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
    if len(images) != 1:
         sys.exit('you can only give write_ipi one frame at a time!')
    frame = images[0].copy()

    symbols = frame.get_chemical_symbols()
    laa = frame.get_cell_lengths_and_angles()
    comment = '# CELL(abcABC):   '+str(laa[0])+"  "+str(laa[1])+"  "+str(laa[2])+"  "+str(round(laa[3]))+"  "+str(round(laa[4]))+"  "+str(round(laa[5]))+" Step: 4  Bead: 0 positions{angstrom} cell{angstrom}"
    newcell, newpos = convert_cell(frame.cell, frame.positions)
    frame.set_cell(newcell)
    frame.set_positions(newpos)

    natoms = len(symbols)
    #for atoms in images:
    fileobj.write('%d\n%s\n' % (natoms, comment))
    for s, (x, y, z) in zip(symbols, frame.positions):
        fileobj.write('%-2s %16.8f %16.8f %16.8f\n' % (s, x, y, z))
