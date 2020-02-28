"""The functions below are for reference only.
We use the implementation from extxyz module, which is backwards
compatible with standard XYZ format."""
import sys
from ase.atoms import Atoms
from ase.io.extxyz import read_extxyz as read_xyz, write_extxyz as write_xyz
import numpy as np
from numpy.linalg import norm
#from copy import copy,deepcopy


__all__ = ['read_quippy', 'write_quippy']

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
    cell = np.matrix.transpose(cell[:])

    if not is_upper_triangular(cell):
        #print('quippy cell is not upper_triangular')
        #print('quippy cell',cell)
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
def read_quippy(fileobj, index):
    lines = fileobj.readlines()
    natoms = int(lines[0])
    cell_str = lines[1]
    #print('nat',natoms)
    #print('quippy cel_str',cell_str)
    #print('cel',type(cell_str))
    #lst = cell_str.split()[1:9]
    lst = cell_str.split('"')[1].split()
    #print('quippy lst xx:',lst)
    #print('quippy lst:',lst)
    cell_out = [float(i) for i in lst]
    #print('out',out)
    #print('quippy cell_out:',cell_out)
    cell_out_f = np.array(cell_out).reshape((3, 3))
    #print('quippy cell_out_f:',cell_out_f)
    #from ase.geometry import cell_to_cellpar
    #cx = cell_to_cellpar(cell_out_f)
    #print('cx',cx)
    #sys.exit('not yet')
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
        yield Atoms(symbols=symbols, positions=positions, cell=cell_out_f, pbc=pbc)


#def simple_write_xyz(fileobj, images, comment=''):
def write_quippy(fileobj, images, comment=''):
    if len(images) != 1:
        sys.exit('you can only give write_quippy one frame at a time!')
    frame = images[0].copy()
    #laa = images[0].get_cell_lengths_and_angles()
    #print('oldpos')
    #print((images[0].positions)[:5])
    #print('--------')
    #print((frame.positions)[:5])
    #print('oldcell quippy:',frame.cell)
    newcell, newpos = convert_cell(frame.cell, frame.positions)
    #print('newcell quippy:',newcell)
    laa = np.matrix.transpose(newcell)
    #print('newpos')
    #print(newpos[:5])
    frame.set_cell(newcell)
    frame.set_positions(newpos)

    # Lattice="14.34366105636912 0 0 7.171830528184559 12.421974858089195 0 7.171830528184559 4.140658286029731 11.71155021051156"
    comment = 'Lattice="'+str(laa[0,0])+"  "+str(laa[0,1])+"  "+str(laa[0,2])+"  "+str(laa[1,0])+"  "+str(laa[1,1])+"  "+str(laa[1,2])+"  "+str(laa[2,0])+"  "+str(laa[2,1])+"  "+str(laa[2,2])+'"'
    #print('aa',(frame.get_positions())[:5])
    symbols = frame.get_chemical_symbols()
    natoms = len(symbols)
    #for atoms in frame:
    fileobj.write('%d\n%s\n' % (natoms, comment))
    for s, (x, y, z) in zip(symbols, frame.positions):
        fileobj.write('%-2s %16.8f %16.8f %16.8f\n' % (s, x, y, z))
