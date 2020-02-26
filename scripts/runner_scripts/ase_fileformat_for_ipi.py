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
    #print('cell',cell,type(cell))
    #print('cell2',cell[0])
    #print('cell3',cell[:])
    #print('cell4',type(cell[:]))
    cell = np.matrix.transpose(cell[:])

    if not is_upper_triangular(cell):
        #print('ipi cell is not upper_triangular')
        #print('ipi cell:',cell)
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
        #print('ipi cell is !! upper_triangular')
        #print('ipi cell:',cell)
        return cell, pos

#def simple_read_xyz(fileobj, index):
def read_ipi(fileobj, index):
    lines = fileobj.readlines()
    natoms = int(lines[0])
    cell_str = lines[1]
    #print('nat',natoms)
    #print('ipi cell_str:',cell_str)
    lst = cell_str.split()[2:8]
    #print('ipi lst:',lst)
    lst2 = cell_str.split()[8:]
    #print('ipi lst2:',lst2)
    conv = False
    for i in lst2:
        #print(i)
        if "positions" in i:
            #print('--p->',i)
            if "atomic_unit" in i:
                conv = 0.52917721
            elif "angstrom" in i:
                conv = 1.
        if "cell" in i:
            #print('--c->',i)
            if "atomic_unit" in i:
                conv = 0.52917721
            elif "angstrom" in i:
                conv = 1.
    if conv == False:
        print('ipi inputfile:',fileobj)
        sys.exit("ipi.py: Error: did not recognise units of ipi inputfile")
    cell_out = [float(i) for i in lst] # [14.34366105636912, 14.34366105636912, 14.34366105636912, 60.0, 60.0, 60.0]
    #print('out',out)
    #print('cell_out:',cell_out)
    cell_out[0] = cell_out[0]*conv
    cell_out[1] = cell_out[1]*conv
    cell_out[2] = cell_out[2]*conv
    #print(cell_out[0],type(cell_out[0]))
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
            positions.append([float(x)*conv, float(y)*conv, float(z)*conv])
        #atoms[0].set_cell(atoms[0].get_cell_lengths_and_angles()))
        yield Atoms(symbols=symbols, positions=positions, cell=cell_out, pbc=pbc)


#def simple_write_xyz(fileobj, images, comment=''):
def write_ipi(fileobj, images, comment='',write_x_reference=False): #,convert_ang_to_bohrradius=False):
    #print('convert_ang_to_bohrradius',convert_ang_to_bohrradius)
    if len(images) != 1:
         sys.exit('you can only give write_ipi one frame at a time!')
    frame = images[0].copy()
    print('frame.pos',frame.positions)

    symbols = frame.get_chemical_symbols()
    laa = frame.get_cell_lengths_and_angles()
    #if convert_ang_to_bohrradius:
    #    atb = 1.8897261
    #    laa[0:3] = laa[0:3]*atb
    comment = '# CELL(abcABC):   '+str(laa[0])+"  "+str(laa[1])+"  "+str(laa[2])+"  "+str(round(laa[3]))+"  "+str(round(laa[4]))+"  "+str(round(laa[5]))+" Step: 4  Bead: 0 positions{angstrom} cell{angstrom}"
    if True:  # for cu_theta_prime (checked with adp) this seems not to do what it is supposed to do.  But without it, this seems to be better
        # nope, acutally this is correct for theta_prime, just that x_refecrece for harmonic has to be taken from this positions.
        newcell, newpos = convert_cell(frame.cell, frame.positions)
        frame.set_cell(newcell)
        frame.set_positions(newpos)
        if type(write_x_reference) == str:
            ang_to_bohr = 1.8897261
            np.savetxt(write_x_reference,frame.positions.flatten()*ang_to_bohr,newline=" ",fmt="%3.15f")


    natoms = len(symbols)
    #for atoms in images:
    fileobj.write('%d\n%s\n' % (natoms, comment))
    for s, (x, y, z) in zip(symbols, frame.positions):
        fileobj.write('%-2s %16.8f %16.8f %16.8f\n' % (s, x, y, z))
    return
