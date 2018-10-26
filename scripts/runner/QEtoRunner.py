#!/usr/bin/env python

import numpy as np
import ase
import argparse
from ase.io import read, write
import os,sys

#print(sys.argv)
#print(len(sys.argv))

def convertQEtoRunner(fileOut,pathRead,indx=-1):
#ASE reads energies from QE in eV
    aToBr, ryToHa=1.88973,1./(13.605691932782346*2.)
    atoms=ase.io.read(filename=pathRead,index=indx, format='espresso-out')
#if-condition is TRUE when $atoms is trajectory and contains more than 1 frame
    if type(indx) is str:
        for frame in range(len(atoms)):
            fileOut.write("begin\ncomment\n")
            for i in range(3):
                acell=atoms[frame].get_cell()[i]*aToBr
                fileOut.write("lattice %.5f %.5f %.5f\n" % (acell[0], acell[1], acell[2]))
            for i in range(atoms[frame].get_number_of_atoms()):
                atCor=atoms[frame].get_positions()[i]*aToBr
                atFor=atoms[frame].get_forces()[i]*ryToHa/aToBr
                potEn=atoms[frame].get_potential_energy()*ryToHa
                element=atoms[frame].get_chemical_symbols()[i]
                fileOut.write("atom   %.6f    %.6f   %.6f %s  0.0   0.0  %.10f  %.10f  %.10f\n"  % (atCor[0], atCor[1],atCor[2],element, atFor[0],atFor[1],atFor[2]) )
            fileOut.write("energy %.5f\n" % (potEn))
            fileOut.write("charge 0\nend\n")
#if-condition is FALSE when $atoms is single frame
    else:
        fileOut.write("begin\ncomment\n")
        for i in range(3):
            acell=atoms.get_cell()[i]*aToBr
            fileOut.write("lattice %.5f %.5f %.5f\n" % (acell[0], acell[1], acell[2]))
        for i in range(atoms.get_number_of_atoms()):
            atCor=atoms.get_positions()[i]*aToBr
            atFor=atoms.get_forces()[i]*ryToHa/aToBr
            potEn=atoms.get_potential_energy()*ryToHa
            element=atoms.get_chemical_symbols()[i]
            fileOut.write("atom   %.6f    %.6f   %.6f %s  0.0   0.0  %.10f  %.10f  %.10f\n"  % (atCor[0], atCor[1],atCor[2],element, atFor[0],atFor[1],atFor[2]) )
        fileOut.write("energy %.5f\n" % (potEn))
        fileOut.write("charge 0\nend\n")
    #print("QE: atPos = Angstrom, E = Ry, forces = Ry/Bohr\nRunner: atPos = Bohr, E = Ha, forces = Ha/Bohr")A

#fileO=open(pathW+filenameW,"w")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('-i', '--inputfolder', type = str, required="True", help = 'Path to the folder which contain Quantum Espresso input')
    parser.add_argument('-o', '--outputfilename', type = str, default="input.data.all", help = 'Path to the folder which contain Quantum Espresso input')
    args = parser.parse_args()
    #print('yes')
    #print(args.inputfolder)
    pathR=args.inputfolder
    pathW=os.getcwd()
    index=-1   # takes last structures of every file
    index="1:5" # takes from each file the 1:5 structure
    index=":"  # takes every structure of every file

    #print(os.getcwd())

    filenameW="out.runner"
    filenameW="input.data"
    filenameW=args.outputfilename
    filename=pathW+"/"+filenameW
    if os.path.isfile(filename):
	sys.exit(filename+" does exist!")
    #print(filename)
    #sys.exit()
    fileO=open(filename,"w")
    for subdir, dirs, files in os.walk(pathR):
    #for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if file.endswith(".out"):
                pathToFileR=os.path.join(subdir, file)
                print(pathToFileR)
                convertQEtoRunner(fileO, pathToFileR,index)
                #convertQEtoRunner(fileO, pathToFileR)
    fileO.close()
    print("written to",filename)
