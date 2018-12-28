#!/usr/bin/env python
from __future__ import print_function
import inspect
import sys,os,copy
import argparse
import numpy as np
from numpy.linalg import norm
import ase
from ase.io import read,write
import myutils as my

#my.exit_if_not_python3()  # for what??

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

def write_lammps_data(infile, atoms, atom_types, comment=None, cutoff=None,
                      molecule_ids=None, charges=None, units='metal'):

    if isinstance(infile, str):
        fh = open(infile, 'w')
    else:
        fh = infile

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


def show_ase_frame_or_frames_content(frame_or_frames):
    #print('frame5',frame_or_frames.cell)
    #print('frame5',frame_or_frames[:].cell)
    #print('tt',type(frame_or_frames))
    #print('t?',type(frame_or_frames) == list)
    #for idx,i  in enumerate(frame_or_frames):
    #    print('idx,i',idx,i)
    if type(frame_or_frames) != list:
        frame_or_frames = [frame_or_frames]
    show_first = 3

    for idx,i  in enumerate(frame_or_frames):
        print("frame:",idx)
        print("------------")
        print("cell:")
        print(frame_or_frames[0].cell)
        #print(frame_or_frames[0].get_cell())
        print()
        print("positions:")
        print(frame_or_frames[idx].positions[:show_first])
        print()
        print("...")
        print(frame_or_frames[0].positions[-show_first:])
        print()
        print("elements:",frame_or_frames[0].get_chemical_symbols()[:show_first],"...",frame_or_frames[0].get_chemical_symbols()[-show_first:])
        print()
        print('kk',frame_or_frames[idx].get_cell())
    return

def read_in_file_or_files_and_make_ase_object(infile,formatin=False,verbose=False):
    ''' reads in a file (or several files) and creates out of it one
    ase object containing one (or several) frames
    '''

    # get all known formats from ase; all known formats can be read.
    known_formats = ase_get_known_formats()

    if formatin in known_formats:
        print('formatin         :',formatin,"(known by ase by default)")
        frame_or_frames = read(infile,':',format=formatin) # all structures
        print('infile (read in) :',infile,"(successfully)")
        print('frames           :',len(frame_or_frames))
        if len(frame_or_frames) == 0:
            sys.exit("NO FRAMES FOUND!")



    #elif formatin == 'ipi':
    #    print('fi ipi')
    #    #frame_or_frames = read(infile,':',format='extxyz')
    #    frame_or_frames = read(infile,':',format='xyz')
    #    print('f',frame_or_frames)
    #    print('f',type(frame_or_frames))
    #    #sys.exit('jo')
    #    #print('fc',frame_or_frames.cell)
    #    #print('fp',frame_or_frames.positions)
    #    #with open(outfilename, 'r') as file:
    #    #    # read a list of lines into data
    #    #    data = file.readlines()
    #    #print('d0',data[0])
    #    #print('d1',data[1])
    #    #print('d2',data[2])
    else:
        sys.exit('formatin '+formatin+' not in the list of known formats! Exit.')


    if args.verbose == True:
        show_ase_frame_or_frames_content(frame_or_frames)

    #sys.exit('jojoa')
    if len(frame_or_frames) == 1:
        return [frame_or_frames]
    else:
        return frame_or_frames

def convert_file(infile, formatin=False,formatout=False,outfilename=False,args=False):

    ###########################################################################
    # read the inputfile
    ###########################################################################
    frame_or_frames = read_in_file_or_files_and_make_ase_object(infile=infile,formatin=formatin,verbose=args.verbose)
    #sys.exit('nach readin')
    ###########################################################################
    # write the outputfile
    ###########################################################################
    print('formatout        :',formatout)

    #print('nowframes (1): ',len(frame_or_frames))
    if args.write_particular_frame != -1: # -1 is the default
        #frame_or_frames = frame_or_frames[args.write_particular_frame]
        frame_or_frames = [frame_or_frames[args.write_particular_frame]]

    #print('nowframes (2): ',len(frame_or_frames))
    print('frames writing   :',len([frame_or_frames]),"(-wf argument)")
    #print('nowframes: ',len(frame_or_frames))

    otherlist = ['lmp', 'lmp.runner','ipi']
    #otherlist = ['lmp', 'lmp.runner']

    known_formats = ase_get_known_formats()

    if len(frame_or_frames) > 1:
        my.get_from_prompt_Yy_orexit("Do you want to write "+str(len(frame_or_frames))+" structures to drive? [Yy]")

    for idx,frameone in enumerate(frame_or_frames):

        if len(frame_or_frames) == 1:
            idx = ""
        #print('idx:',idx,type(idx))
        outfilename = get_outfilename(args,idx)
        #print('aaa',outfilename)
        if formatout in otherlist:
            print('convert_fileformats.py: formatout        :',formatout,"(not default fileformat, but known - lmp or lmp.runner -)")
            if formatout == 'lmp':
                save_ase_object_as_lmp(frameone,outfilename,comment=infile,runner=False)
            if formatout == 'lmp.runner':
                save_ase_object_as_lmp_runner(frameone,outfilename,comment=infile)
            if formatout == 'ipi':
                save_ase_object_as_ipi_format(frameone,outfilename)

        elif formatout in known_formats:
            print('formatout        :',formatout,"known by default")
            save_ase_object_in_ase_format(frameone,outfilename,formatout)
        else:
            sys.exit('unknown output format '+formatout)

    my.create_READMEtxt(os.getcwd())
    return


def save_ase_object_in_ase_format(ase_object,outfilename,formatout):
    ase.io.write(outfilename,ase_object,format=formatout)
    if formatout == 'espresso-in':
        f = open(outfilename,"r")
        lines = f.readlines()
        f.close()
        print('change some stuff')
        settings = os.environ['scripts']+"/qe-aiida/aiida_submitskripts/aiida.in.top"
    print('written (2)      : '+outfilename)


def save_ase_object_as_ipi_format(frame,outfilename):
    ase.io.write(outfilename,frame,format='xyz')
    laa = frame.get_cell_lengths_and_angles()
    with open(outfilename, 'r') as file:
        # read a list of lines into data
        data = file.readlines()

    data[1] = '# CELL(abcABC):   '+str(laa[0])+"  "+str(laa[1])+"  "+str(laa[2])+"  "+str(round(laa[3]))+"  "+str(round(laa[4]))+"  "+str(round(laa[5]))+"  Step: 4  Bead: 0 positions{angstrom}  cell{angstrom}"+'\n'

    # and write everything back
    with open(outfilename, 'w') as file:
        file.writelines( data )
    print('written (3)      : '+outfilename)
    return


def save_ase_object_as_lmp_runner(frame,outfilename,comment=""):
    save_ase_object_as_lmp(frame,outfilename,comment=comment,runner=True)
    return


def save_ase_object_as_lmp(frame,outfilename,comment="",runner=False):
    if type(frame) != list:
        frame = [frame]
    # important not to change frame (aseobject), i dont know why but when this is called
    # externally it changes external frame (aseobject)
    #print('frame2',frame.cell)
    framework = copy.deepcopy(frame)
    #print('frame3',frame.cell)
    #print('frame4',framework.cell)
    #framework = frame
    #print('----------------------------------------xxxxxxxxxxxxx')
    #show_ase_frame_or_frames_content(frame)
    #print('----------------------------------------xxxxxxxxxxxxx')
    #show_ase_frame_or_frames_content(framework)
    #sys.exit()
    #print('fw',framework.positions)
    framework = framework[0]
    #sys.exit()

    # this converts from ase to lammps; newcell,newpos give same energies?
    # I would actually need something which reverts this
    newcell, newpos = convert_cell(framework.cell, framework.positions)
    framework.set_cell(newcell)
    framework.set_positions(newpos)

    species = framework.get_atomic_numbers()
    #print('species',species)
    unique_species, i_unique_species = np.unique(species, return_index=True)
    #print('unique_species',unique_species)
    masses = framework.get_masses()
    positions = framework.get_positions()

    lammps_indices = species.copy()
    for ii in range(len(unique_species)):
        lammps_indices[species == unique_species[ii]] = ii + 1

    ((hxx, hxy, hxz), (hyx, hyy, hyz) , (hzx, hzy, hzz)) = framework.get_cell()


    fout = open(outfilename,'w')

    fout.write("LAMMPS positionsfile from: "+str(comment)+"\n")
    fout.write("\n")
    fout.write(str(len(species))+" atoms\n")

    if runner is False:  # "normal"  lammps format
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
        fout.write("1 24.305\n")        # Mg
        fout.write("2 26.9815385\n")    # Al
        fout.write("3 28.085\n")        # Si
        species_runner = species.copy()
        species_runner[species_runner == 12] = 1
        species_runner[species_runner == 13] = 2
        species_runner[species_runner == 14] = 3
        lammps_indices = species_runner
        #print('sp',species_runner)

    fout.write("\n")
    fout.write("Atoms\n")
    fout.write("\n")
    # take care here for runner / n2p2 if only 2 species (say Mg,Si) for Mg,Al,Si potential
    for ii in range(len(species)):  # if only Mg and Si get is, there will be only type 1 and 2 and not 1 and 3!
        fout.write(("%5d %5d %s" % (ii + 1, lammps_indices[ii], "  ".join(map(str, positions[ii]))) )+"\n")
        #fout.write(("%5d %5d %.10f %.10f %.10f" % (ii + 1, lammps_indices[ii], "  ".join(map(str, positions[ii]))) )+"\n")
    fout.write("\n")

    fout.close()
    print('written (4)      : '+outfilename)
    return

def check_if_ase_knows_runner():
    pass

def ase_get_known_formats(show=False,add_missing_formats=False,verbose=True):
    ''' adds formats runner and lammps-runner to ase '''
    known_formats = []
    x = ase.io.formats.all_formats
    for i in x:
        known_formats.append(i)

    if show:
        import pprint
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(x)

    if add_missing_formats:
        if verbose:
            print('cc',ase.io.__file__)
        formatspy = os.path.dirname(ase.io.__file__)+"/formats.py"
        if verbose:
            print('formatspy',formatspy)

        scripts = my.scripts()
        missing1 = scripts+"/runner_scripts/ase_fileformat_for_runner.py"
        missing2 = scripts+"/runner_scripts/ase_fileformat_for_lammpsrunner.py"
        missing3 = scripts+"/runner_scripts/ase_fileformat_for_lammpsdata.py"
        from shutil import copyfile
        if verbose:
            print('copying files to',os.path.dirname(ase.io.__file__))
        copyfile(missing1,os.path.dirname(ase.io.__file__)+"/runner.py")
        copyfile(missing2,os.path.dirname(ase.io.__file__)+"/lammpsrunner.py")
        copyfile(missing3,os.path.dirname(ase.io.__file__)+"/lammpsdata.py")

        if 'runner' in x:
            print('missing format are already added in formats.py (of ase).')
        else:
            print('adapting ase formats.py .... ')
            if not os.path.isfile(formatspy):
                print('formatspy',formatspy)
                sys.exit('did not find '+str(formatspy))

            print('now changing formatspy')

            f = open(formatspy, "r")
            contents = f.readlines()
            f.close()
            insert=0
            insert2=0
            for idx,i in enumerate(contents):
                #print('i',idx,i)
                #print("|"+i[:20]+"|")
                if i[:20] == "    'abinit': ('ABIN":
                    insert = idx
                if i[:30] == "    'lammps-data': 'lammpsdata":
                    insert2 = idx

            contents.insert(insert, "    'runner': ('Runner input file', '+F'),\n")
            contents.insert(insert, "    'lammps-runner': ('LAMMPS data input file for n2p2 or runner', '1F'),\n")
            contents.insert(insert2, "    'lammps-runner': 'lammpsrunner',\n")
            print('insert',insert)

            f = open(formatspy, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()


    return known_formats

if __name__ == "__main__":
    #p = argparse.ArgumentParser(description=pp.pprint(x),
    p = argparse.ArgumentParser(description='',
            formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    string='''
    e.g.
    convert_fileformats.py PathToInputfile.out --formatin 'espresso-out' --formatout lmp
    convert_fileformats.py -fi runner -fo espresso-in input_selected.data.20
    convert_fileformats.py -fi runner -fo espresso-in input_selected.data.1


    '''
    parser = argparse.ArgumentParser(description=string)
    parser.add_argument("infile", type = str, help = "The name of the file")
    parser.add_argument("--showformats",'-sf', action='store_true', default=False, help = "show the possible in/outputformats (not showing lammps) and exit")
    parser.add_argument("--make_ase_runner",'-ar', action='store_true', default=False, help = "adapt ase so that it can read/write runner format.")
    parser.add_argument("--formatin",'-fi', type=str,default="lmp", help="Format of the input file, this can be all of the listed using --showformats, e.g. espresso-out,xyz,vasp, ...; default:lmp")
    parser.add_argument("--formatout",'-fo', type=str,default="lmp", help="Format of the output file, this can be most of the listed using --showformats and additionally lmp for LAMMPS, e.g. lmp,xyz,vasp, ...; default:lmp")
    parser.add_argument("--outfilename", type=str,default=False, help="if nothing is provided --> default=infile+.lmp or when split option: infile+NUMBER+.lmp")
    parser.add_argument("--write_particular_frame",'-wf', type=int,default=-1, help="when writing of files, specify here if a particular frame should be written; (mind: 0 is frame one, 1 is frame 2, ...); default: all frames are written")
    parser.add_argument('-v','--verbose', help='verbose', action='count', default=False)


    args = parser.parse_args()
    if args.showformats or args.make_ase_runner:
        known_formats = ase_get_known_formats(show=True,add_missing_formats=args.make_ase_runner)
        sys.exit()

    def get_outfilename(args,frame=""):
        if frame == "":
            addidx =""
        elif type(frame) is int:
            addidx = "_"+str(frame)

        if args.outfilename == False:
            outfilename = os.getcwd()+'/'+os.path.basename(args.infile)+addidx+'.'+args.formatout
        else:
            outfilename = os.getcwd()+'/'+args.outfilename+addidx+"."+args.formatout  #args.outfilename #+'.'+args.formatout
        return outfilename

    outfilename = get_outfilename(args,frame="")

    sys.exit(convert_file(args.infile, args.formatin, args.formatout,outfilename,args))
