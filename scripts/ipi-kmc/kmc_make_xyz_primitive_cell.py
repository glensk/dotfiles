#!/usr/bin/env python

from ase.lattice.cubic import FaceCenteredCubic
import ase
import numpy as np
import sys
import random
import argparse

def help(p = None):
    description = '''
    description
    '''
    p = argparse.ArgumentParser(description=description,formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    p.add_argument('-sc',   '--sc',default=8,type=int,help="supercell size")
    p.add_argument('-nmg',   '--nmg',default=0,type=int,help="number of Mg atoms")
    p.add_argument('-nsi',   '--nsi',default=0,type=int,help="number of Si atoms")
    p.add_argument('-nvac', '--nvac',default=0,type=int,help="number of vacancies")
    p.add_argument('-alat', '--alat',default=4.057,type=float,help="lattice constant")
    return p

p = help()
args = p.parse_args()
print('args',args)

alat=args.alat
sc=args.sc
nsi=args.nsi
nmg=args.nmg
nva=args.nvac


#if False: # 6x6x6 sc daniele
#    alat = 4.06
#    sc = 6
#    nsi = 6
#    nmg = 6
#    nva = 2
#
#if False:  # 8x8x8 structure daniele
#    alat = 4.057
#    sc = 8
#    nsi = 12
#    nmg = 12
#    nva = 2
#
#if False: # 8x8x8 structure own
#    alat = 4.057
#    sc = 8
#    nsi = 6        # 1.17 %
#    nmg = 5        # 1 %
#    nva = 1        # 0.001953125
#
#if True: # 8x8x8 structure own
#    alat = 4.057
#    sc = 8
#    nsi = 6        # 1.17 %
#    nmg = 6        # 1 %
#    nva = 1        # 0.001953125
#
#
#if True: # 10x10x10 structure own
#    alat = 4.057
#    sc = 10
#    nsi = 12        # 1.17 %
#    nmg = 12        # 1 %
#    nva = 2        # 0.001953125
#
#if False: # 10x10x10 structure own
#    alat = 4.057
#    sc = 10
#    nsi = 6        # 1.17 %
#    nmg = 5        # 1 %
#    nva = 1        # 0.001953125
#
#if False: # 16x16x16 structure own
#    alat = 4.057
#    sc = 16
#    nsi = 6        # 1.17 %
#    nmg = 5        # 1 %
#    nva = 1        # 0.001953125



atom = ase.build.bulk('Al',crystalstructure='fcc',a=alat)
atomsc = atom.repeat(sc)

number_of_atoms = atomsc.get_number_of_atoms()
print('nat:',number_of_atoms)
nal = number_of_atoms - nsi - nmg
print('nal:',nal)
print('nmg:',nmg)
print('nsi:',nsi)
print('nva:',nva)
print('number_of_atoms (before removing vac)',number_of_atoms)

delrandom = False
if delrandom:
    randomlist= np.arange(number_of_atoms - nva)
    random.shuffle(randomlist)

print()
if True: # 8x8x8 supercell of daniele
    for i in np.arange(nsi):
        atomsc[i].symbol = 'Si'
    for i in np.arange(nsi,nmg+nsi):
        atomsc[i].symbol = 'Mg'
    #for i in np.arange(nmg+nsi,nmg+nsi+nva):
    #    del atomsc[i]
    for i in np.arange(nva):
        del atomsc[-1]

if False: # 6x6x6 supercell of daniele
    for i in np.arange(nsi):
        atomsc[i].symbol = 'Si'
    for i in np.arange(nsi,nmg+nsi):
        atomsc[i].symbol = 'Mg'

if False: # for my random supercells
    #print(randomlist[:20])
    for i in randomlist[0:nmg]:
        print('i Mg',i)
        atomsc[i].symbol = 'Mg'
    for i in randomlist[nmg:nmg+nsi]:
        print('i Si',i)
        atomsc[i].symbol = 'Si'
    #print()
    #print('rrr',randomlist[nmg+nsi:nmg+nsi+nva])
    #print('rrr',np.sort(randomlist[nmg+nsi:nmg+nsi+nva]))
    #print('rrr',np.sort(randomlist[nmg+nsi:nmg+nsi+nva])[::-1])
    #for i in randomlist[nmg+nsi:nmg+nsi+nva]:
    for i in np.sort(randomlist[nmg+nsi:nmg+nsi+nva])[::-1]:  # example: [759 703 455 234 158]
        print('i Va',i)
        #print(atomsc[i-1])
        #print(atomsc[i])
        #print(atomsc[i+1])
        del atomsc[i]
        #print(atomsc[i-1])
        #print(atomsc[i])
        #print(atomsc[i+1])
        #print()

print()
number_of_atoms = atomsc.get_number_of_atoms()
print('nat:',number_of_atoms)
nal = number_of_atoms - nsi - nmg
print('nal:',nal)
print('nmg:',nmg)
print('nsi:',nsi)
print('nva:',nva)
print('number_of_atoms (after removing vac)',number_of_atoms)


#for i in np.arange(nmg+nsi,nmg+nsi+nva):
#    del atomsc[i]

#from ase.visualize import view
#atom
#view(atom)

#atom2 = atom.repeat(2)
#atom2
#atom2.positions
#atom2.cell
#view(atom2)

#atomxyz = atom.repeat((2,2,4))
# WRITE XYZ FILE
#sys.exit()

print(str(sc)+"x"+str(sc)+"x"+str(sc)+"_alat"+str(alat))
scstr = str(sc)+"x"+str(sc)+"x"+str(sc)+"_alat"+str(alat)+"_"+str(nal)+"al_"+str(nsi)+"si_"+str(nmg)+"mg_"+str(nva)+"va_"+str(number_of_atoms)+"atoms"
if delrandom is True:
    scstr = scstr + "_random"
#print('kk',scstr)
filenameoutxyz = 'al'+scstr+'.xyz'
ase.io.write(filenameoutxyz,atomsc)

#filenameoutlammps = 'al'+scstr+'.lammpsdata'
#ase.io.write(filenameoutlammps,atomsc,format='lammps-dump')
#ase.io.write(filenameoutlammps,atomsc)



# WRITE XYZ FILE for i-pi
filenameoutipi = 'al'+scstr+'.ipi'
ase.io.write(filenameoutipi,atomsc,format='xyz')
laa = atomsc.get_cell_lengths_and_angles()
print('atomsc.get_cell_lengths_and_angles',laa)

with open(filenameoutipi, 'r') as file:
    # read a list of lines into data
    data = file.readlines()
print('kkk',str(round(laa[3])))
# now change the 2nd line, note that you have to add a newline
data[1] = '# CELL(abcABC):   '+str(laa[0])+"  "+str(laa[1])+"  "+str(laa[2])+"  "+str(round(laa[3]))+"  "+str(round(laa[4]))+"  "+str(round(laa[5]))+"  Step:           4  Bead:       0 positions{angstrom}  cell{angstrom}"+'\n'

# and write everything back
with open(filenameoutipi, 'w') as file:
    file.writelines( data )


from subprocess import call
file_to_exec = '/Users/glensk/Dropbox/Albert/scripts/dotfiles/scripts/ipi-kmc/kmc_xyz_or_ase_to_lammpsinput.py'
call(["python", file_to_exec,filenameoutxyz ])
