#!/usr/bin/env python
import sys,os
import argparse
import numpy as np
from numpy.linalg import norm
import ase
from ase.io import read,write


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

def write_lammps_data(filename, atoms, atom_types, comment=None, cutoff=None,
                      molecule_ids=None, charges=None, units='metal'):

    if isinstance(filename, str):
        fh = open(filename, 'w')
    else:
        fh = filename

    if comment is None:
        comment = 'lammpslib autogenerated data file'
    fh.write(comment.strip() + '\n\n')

    fh.write('{0} atoms\n'.format(len(atoms)))
    fh.write('{0} atom types\n'.format(len(atom_types)))


    fh.write('\n')
    cell, coord_transform = convert_cell(atoms.get_cell())
    fh.write('{0:16.8e} {1:16.8e} xlo xhi\n'.format(0.0, cell[0, 0]))
    fh.write('{0:16.8e} {1:16.8e} ylo yhi\n'.format(0.0, cell[1, 1]))
    fh.write('{0:16.8e} {1:16.8e} zlo zhi\n'.format(0.0, cell[2, 2]))
    fh.write('{0:16.8e} {1:16.8e} {2:16.8e} xy xz yz\n'
             ''.format(cell[0, 1], cell[0, 2], cell[1, 2]))

    fh.write('\nMasses\n\n')
    sym_mass = {}
    masses = atoms.get_masses()
    symbols = atoms.get_chemical_symbols()
    for sym in atom_types:
        for i in range(len(atoms)):
            if symbols[i] == sym:
                sym_mass[sym] = masses[i]
                break
            else:
                sym_mass[sym] = (atomic_masses[chemical_symbols.index(sym)] /
                                 unit_convert("mass", units))

    for (sym, typ) in sorted(list(atom_types.items()), key=operator.itemgetter(1)):
        fh.write('{0} {1}\n'.format(typ, sym_mass[sym]))

    fh.write('\nAtoms # full\n\n')
    if molecule_ids is None:
        molecule_ids = np.zeros(len(atoms), dtype=int)
    if charges is None:
        charges = atoms.get_initial_charges()
    for i, (sym, mol, q, pos) in enumerate(
            zip(symbols, molecule_ids, charges, atoms.get_positions())):
        typ = atom_types[sym]
        fh.write('{0} {1} {2} {3:16.8e} {4:16.8e} {5:16.8e} {6:16.8e}\n'
                 .format(i + 1, mol, typ, q, pos[0], pos[1], pos[2]))

    if isinstance(filename, str):
        fh.close()




def main(filename, fileformat_in=False,fileformat_out=False,outputfile=False):
    #if outputfile == False:
    #    outputfile = filename+'.lmp'
    frame = read(filename,format=fileformat_in)
    otherlist = ['lmp', 'lmp.runner','ipi']
    if fileformat_out not in otherlist:
        write(outputfile,frame,format=fileformat_out)
        print('written '+outputfile)
        sys.exit()
    elif fileformat_out is 'lmp':
        save_ase_object_as_lmp(frame,outputfile,comment=filename,runner=False)
    elif fileformat_out is 'lmp.runner':
        save_ase_object_as_lmp_runner(frame,outputfile,comment=filename)


def save_ase_object_as_lmp_runner(frame,outputfile,comment=""):
    save_ase_object_as_lmp(frame,outputfile,comment=comment,runner=True)


def save_ase_object_as_lmp(frame,outputfile,comment="",runner=False):
    newcell, newpos = convert_cell(frame.cell, frame.positions)
    frame.set_cell(newcell)
    frame.set_positions(newpos)

    species = frame.get_atomic_numbers()
    unique_species, i_unique_species = np.unique(species, return_index=True)
    masses = frame.get_masses()
    positions = frame.get_positions()

    lammps_indices = species.copy()
    for ii in range(len(unique_species)):
        lammps_indices[species == unique_species[ii]] = ii + 1

    ((hxx, hxy, hxz), (hyx, hyy, hyz) , (hzx, hzy, hzz)) = frame.get_cell()


    fout = open(outputfile,'w')
    #print('jo',outputfile)

    fout.write("LAMMPS positionsfile from: "+str(comment)+"\n")
    fout.write("\n")
    fout.write(str(len(species))+" atoms\n")
    if runner is False:
        fout.write(str(len(unique_species))+" atom types\n")
    elif runner is True:
        # in principle this can be automatically obtained from the nn potential,
        # for now however I am lazy
        fout.write("3 atom types\n")
    else:
        sys.exit('variable runner is neigher True nor False; Exit.')
    fout.write("\n")
    fout.write("0.0 "+str(hxx)+" xlo xhi\n")
    fout.write("0.0 "+str(hyy)+" ylo yhi\n")
    fout.write("0.0 "+str(hzz)+" zlo zhi\n")
    fout.write(str(hxy)+" "+str(hxz)+" "+str(hyz)+" xy xz yz\n")
    fout.write("\n")
    fout.write("Masses\n")
    fout.write("\n")
    if runner is False:
        for ii in range(len(unique_species)):
            #print(lammps_indices[i_unique_species[ii]], masses[i_unique_species[ii]])
            fout.write(str(lammps_indices[i_unique_species[ii]])+" "+str(masses[i_unique_species[ii]])+"\n")
    elif runner is True:
        # in principle this can be automatically obtained from the nn potential,
        # for now however I am lazy
        fout.write("1 24.305\n")
        fout.write("2 26.9815385\n")
        fout.write("3 28.085\n")
    fout.write("\n")
    fout.write("Atoms\n")
    fout.write("\n")
    for ii in range(len(species)):
        fout.write(("%5d %5d %s" % (ii + 1, lammps_indices[ii], "  ".join(map(str, positions[ii]))) )+"\n")
        #fout.write(("%5d %5d %.10f %.10f %.10f" % (ii + 1, lammps_indices[ii], "  ".join(map(str, positions[ii]))) )+"\n")
    fout.write("\n")

    fout.close()



if __name__ == "__main__":
    #p = argparse.ArgumentParser(description=pp.pprint(x),
    p = argparse.ArgumentParser(description='',
            formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    string='''e.g. convert_fileformats.py PathToInputfile.out --formatin 'espresso-out' --formatout lmp'''
    parser = argparse.ArgumentParser(description=string)
    parser.add_argument("infile", type = str, help = "The name of the file")
    parser.add_argument("--showformats",'-sf', action='store_true', default=False, help = "show the possible in/outputformats (not showing lammps) and exit")
    parser.add_argument("--formatin",'-fi', type=str,default="lmp", help="Format of the input file, this can be all of the listed using --showformats, e.g. espresso-out,xyz,vasp, ...; default:lmp")
    parser.add_argument("--formatout",'-fo', type=str,default="lmp", help="Format of the output file, this can be most of the listed using --showformats and additionally lmp for LAMMPS, e.g. lmp,xyz,vasp, ...; default:lmp")
    parser.add_argument("--outfilename", type=str,default=False, help="if nothing is provided --> default=infile+.lmp")


    args = parser.parse_args()
    if args.showformats:
        x = ase.io.formats.all_formats
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(x)
        sys.exit()

    if args.outfilename == False:
        #args.outfilename = args.infile+'.'+args.formatout
        args.outfilename = os.getcwd()+'/'+os.path.basename(args.infile)+'.'+args.formatout
    else:
        args.outfilename = os.getcwd()+'/'+args.outfilename+"."+args.formatout  #args.outfilename #+'.'+args.formatout

    print()
    print('infile       :',args.infile)
    print('infileformat :',args.formatin)
    print('outfileformat:',args.formatout)
    print('outfilename  :',args.outfilename)
    print()
    sys.exit(main(args.infile, args.formatin, args.formatout,args.outfilename))
