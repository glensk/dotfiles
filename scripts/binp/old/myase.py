#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,shutil
from socket import gethostname

## import ase stuff
import ase
from ase.build import bulk as ase_build_bulk
from ase.constraints import StrainFilter
from ase.neighborlist import NeighborList, neighbor_list, NewPrimitiveNeighborList
from ase.constraints import ExpCellFilter
from ase.spacegroup import crystal
from ase.constraints import StrainFilter
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.optimize import BFGS
from ase.optimize import LBFGS
from ase.optimize import FIRE
from ase.optimize import GPMin
from ase.optimize.basin import BasinHopping
from ase.optimize.minimahopping import MinimaHopping
from ase import units as aseunits

from myutils import mypot,pot_all

try:
    from ase.calculators.lammpslib import LAMMPSlib
except ImportError:
    print("ERROR when importing LAMMPSlib ... possibly you have to change your (conda/aiida) environment")

##################################################################################
## ase funcions
##################################################################################
def get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0,cubic=False,create_fake_vacancy=False,whichcell="fcc",normal_ordering=True):
    """Creating bulk systems.

        Crystal structure and lattice constant(s) will be guessed if not
        provided.

        name: str
            Chemical symbol or symbols as in 'MgO' or 'NaCl'.
        crystalstructure: str
            Must be one of sc, fcc, bcc, hcp, diamond, zincblende,
            rocksalt, cesiumchloride, fluorite or wurtzite.
        a: float
            Lattice constant.
        c: float
            Lattice constant.
        covera: float
            c/a ratio used for hcp.  Default is ideal ratio: sqrt(8/3).
        u: float
            Internal coordinate for Wurtzite structure.
        orthorhombic: bool
            Construct orthorhombic unit cell instead of primitive cell
            which is the default.
        cubic: bool
            Construct cubic unit cell if possible.
    """
    if whichcell == "fcc":
        atom = ase_build_bulk('Al',crystalstructure='fcc',a=a0,cubic=cubic)
    elif whichcell == "hcp":
        a = 3.21
        c = 5.21
        atom = crystal('Mg', [(1./3., 2./3., 3./4.)], spacegroup=194, cellpar=[a, a, c, 90, 90, 120])
    elif whichcell == "dc":
        a = 5.47
        atom = crystal('Si', [(0,0,0)], spacegroup=227, cellpar=[a, a, a, 90, 90, 90])
    else:
        sys.exti("whichcell has to be in fcc or hcp")

    atomsc = atom.repeat(ncell)
    number_of_atoms = atomsc.get_number_of_atoms()
    nal = number_of_atoms - nsi - nmg

    #for i in np.arange(nmg):
    #    atomsc[i].symbol = 'Mg'
    #for i in np.arange(nmg,nmg+nsi):
    #    atomsc[i].symbol = 'Si'

    # This is the order which ipi kmc expects
    if type(normal_ordering) == bool:
        for i in np.arange(nsi):
            atomsc[i].symbol = 'Si'
        for i in np.arange(nsi,nmg+nsi):
            atomsc[i].symbol = 'Mg'

        if create_fake_vacancy == False:
            for i in np.arange(nvac):
                del atomsc[-1]
        elif create_fake_vacancy == True:
            startsubst = -1
            for i in np.arange(nvac):
                #print('startsubst',startsubst)
                atomsc[startsubst].symbol = 'V'
                startsubst -= 1
        else:
            sys.exit("create_fake_vacancy has to be True or False")
    elif type(normal_ordering) == str:
        normal_ordering_element = normal_ordering.split("_")[0]
        normal_ordering_pos     = normal_ordering.split("_")[1]

        for i in np.arange(nvac):
            atomsc[i].symbol = 'V'
        #print('normal_ordering_element',normal_ordering_element)
        #symb = 'Si'
        #symb = 'Mg'
        symb = normal_ordering_element
        if normal_ordering_pos == '0':
            return atomsc


        if normal_ordering_pos == '1':
            atomsc[1].symbol = symb
        if normal_ordering_pos == '2':
            atomsc[5].symbol = symb
        if normal_ordering_pos == '3':
            atomsc[25].symbol = symb
        if normal_ordering_pos == '4':
            atomsc[1].symbol = symb
            atomsc[5].symbol = symb
        if normal_ordering_pos == '5':
            atomsc[1].symbol = symb
            atomsc[5].symbol = symb
            atomsc[25].symbol = symb
    else:
        sys.exit('random_ordering has to be of type bool or str!')
    #number_of_atoms = atomsc.get_number_of_atoms()
    #nal = number_of_atoms - nsi - nmg
    #ase.io.write('kalmp',atomsc,format='lammps-dump')
    return atomsc

def ase_get_neighborlist(frame,atomnr=0,cutoff=3.,skin=0.1):
    NN = NeighborList([cutoff/2.]*frame.get_number_of_atoms(),skin=skin,self_interaction=False,bothways=True,primitive=NewPrimitiveNeighborList)
    #print('NN')
    NN.update(frame)
    #print(NN.get_connectivity_matrix())
    #for idx,i in enumerate(NN.get_connectivity_matrix()):
    #    print('idx',idx,i)
    NN_indices, offsets = NN.get_neighbors(atomnr)
    #print('NN idx:',np.sort(NN_indices))
    #sys.exit()
    return np.sort(NN_indices)

def count_amount_1NN_around_vacancies(filename,cutoffa=3.,cutoffb=4.5,skin=0.1,format='ipi',vac_symbol="V",save_every = 10,filename_save="KMC_analyze_akmc_ext"):
    print()
    print("########### count_amount_1NN_around_vacancies(...) #######################")
    print('reading',os.path.abspath(filename),'...')
    frames = ase_read(filename,index=":",format=format)
    print('reading',os.path.abspath(filename),'done.')

    structures = len(frames)
    structures = structures - 1 # just to make the lenght equal between
    print('structures',structures)

    all_vac_idx = ([atom.index for atom in frames[0] if atom.symbol == vac_symbol])
    print('all_vac_idx',all_vac_idx)

    if cutoffa == False:
        print(frames[0].cell)
        print(frames[0].get_number_of_atoms())
        #a0 = frames[0].cell[0,0]/5.*sqrt(2.)
        nndist = cutoffa = 2.95
    if cutoffb == False:
        a0 = cutoffb = 4.29

    if cutoffa == False:sys.exit('cutoffa')
    if cutoffb == False:sys.exit('cutoffb')
    print('cutoffa',cutoffa)
    print('cutoffb',cutoffb)

    filename_analyze_all = []
    al_mg_si_all = []
    for vac_nr,vac_idx in enumerate(all_vac_idx):
        #filename_analyze = filename +  ".1NN.al_mg_si_vac_"+str(vac_nr)+".dat"
        filename_analyze = filename_save+"_"+str(vac_nr)
        filename_analyze_all.append(filename_analyze)
        print('filename_analyze',filename_analyze)

        if os.path.isfile(filename_analyze):
            al_mg_si = np.loadtxt(filename_analyze)
            al_mg_si_all.append(al_mg_si)
        else:
            al_mg_si = np.zeros((structures,7))
            al_mg_si_all.append(al_mg_si)
    print()
    print('filename_analyze_all',filename_analyze_all)
    print()

    def test_anz_nn(al_mg_si_all,vac_nr,step,exit=False):
        #if al_mg_si_all[vac_nr][step][0] == 0:
        #    return True
        anz_1NN = np.sum(al_mg_si_all[vac_nr][step][1:4])
        anz_2NN = np.sum(al_mg_si_all[vac_nr][step][4:7])
        #print(al_mg_si[step][1:4])
        #print('sum',np.sum(al_mg_si[step][1:4]))
        testanz = True
        do_continue = True
        printed = False
        if testanz:
            if anz_1NN != 12:
                #print('step:',str(step).ljust(6),"not 12 1NN but "+str(anz_1NN),'exit',exit)
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_a',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'exit',exit,"not 12 1NN but "+str(anz_1NN),'-> REDO')
                printed = True
                do_continue = False
                if exit == True: sys.exit("ERROR")
            if anz_2NN != 6 and printed == False:
                #print('step:',str(step).ljust(6),"not  6 2NN but "+str(anz_2NN),'exit',exit)
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_b',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'exit',exit,"not  6 2NN but "+str(anz_2NN),'-> REDO')
                printed = True
                do_continue = False
                if exit == True: sys.exit("ERROR")
        if printed == False:
            print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_c',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'continue',do_continue)
        return do_continue

    for step in np.arange(structures):
        all_vac_idx = ([atom.index for atom in frames[step] if atom.symbol == vac_symbol])
        #print('step',step,'all_vac_idx',all_vac_idx)
        for vac_nr,vac_idx in enumerate(all_vac_idx):
            if al_mg_si_all[vac_nr][step,0] != 0:
                do_continue = test_anz_nn(al_mg_si_all,vac_nr,step,exit=False)
                if do_continue == True:
                    continue

            NN_1_indices, NN_2_indices = ase_get_neighborlist_1NN_2NN(frames[step],atomnr=vac_idx,cutoffa=cutoffa,cutoffb=cutoffb,skin=skin)
            NN_1_sym = [atom.symbol for atom in frames[step] if atom.index in NN_1_indices]
            NN_2_sym = [atom.symbol for atom in frames[step] if atom.index in NN_2_indices]
            NN_1_al = NN_1_sym.count("Al")
            NN_1_mg = NN_1_sym.count("Mg")
            NN_1_si = NN_1_sym.count("Si")
            NN_2_al = NN_2_sym.count("Al")
            NN_2_mg = NN_2_sym.count("Mg")
            NN_2_si = NN_2_sym.count("Si")
            al_mg_si_all[vac_nr][step,1] = NN_1_al
            al_mg_si_all[vac_nr][step,2] = NN_1_mg
            al_mg_si_all[vac_nr][step,3] = NN_1_si
            al_mg_si_all[vac_nr][step,4] = NN_2_al
            al_mg_si_all[vac_nr][step,5] = NN_2_mg
            al_mg_si_all[vac_nr][step,6] = NN_2_si
            al_mg_si_all[vac_nr][step,0] = step
            anz_1NN = np.sum(al_mg_si_all[vac_nr][step][1:4])
            anz_2NN = np.sum(al_mg_si_all[vac_nr][step][4:7])
            #str_1NN = str(NN_1_al).ljust(3)+str(NN_1_mg).ljust(3)+str(NN_1_si).ljust(3)
            #str_2NN = str(NN_2_al).ljust(3)+str(NN_2_mg).ljust(3)+str(NN_2_si).ljust(3)
            #print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, "||",str(NN_1_al).ljust(3),NN_1_mg,NN_1_si,"||",NN_2_al,NN_2_mg,NN_2_si)
            #print('-------',filename_analyze_all[vac_nr])
            #print(al_mg_si_all[vac_nr])
            if anz_1NN != 12 or anz_2NN != 6:
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'ERROR!')
                sys.exit("ERROR see above")
            print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'OK')
            #do_continue = test_anz_nn(al_mg_si_all,vac_nr,step,exit=True)

            if step > 0 and step in np.arange(structures)[::save_every]:
                np.savetxt(filename_analyze_all[vac_nr],al_mg_si_all[vac_nr],fmt='%i')
                print('saving',os.path.abspath(filename_analyze_all[vac_nr]),'at step',step)

    # save everything in the very end
    np.savetxt(filename_analyze_all[vac_nr],al_mg_si_all[vac_nr],fmt='%i')
    print('saving (very end)',os.path.abspath(filename_analyze_all[vac_nr]),'at step',step)
    return

    def test_anz_nn(al_mg_si_all,vac_nr,step,exit=False):
        #if al_mg_si_all[vac_nr][step][0] == 0:
        #    return True
        anz_1NN = np.sum(al_mg_si_all[vac_nr][step][1:4])
        anz_2NN = np.sum(al_mg_si_all[vac_nr][step][4:7])
        #print(al_mg_si[step][1:4])
        #print('sum',np.sum(al_mg_si[step][1:4]))
        testanz = True
        do_continue = True
        printed = False
        if testanz:
            if anz_1NN != 12:
                #print('step:',str(step).ljust(6),"not 12 1NN but "+str(anz_1NN),'exit',exit)
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_a',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'exit',exit,"not 12 1NN but "+str(anz_1NN),'-> REDO')
                printed = True
                do_continue = False
                if exit == True: sys.exit("ERROR")
            if anz_2NN != 6 and printed == False:
                #print('step:',str(step).ljust(6),"not  6 2NN but "+str(anz_2NN),'exit',exit)
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_b',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'exit',exit,"not  6 2NN but "+str(anz_2NN),'-> REDO')
                printed = True
                do_continue = False
                if exit == True: sys.exit("ERROR")
        if printed == False:
            print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_c',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'continue',do_continue)
        return do_continue

    for step in np.arange(structures):
        all_vac_idx = ([atom.index for atom in frames[step] if atom.symbol == vac_symbol])
        #print('step',step,'all_vac_idx',all_vac_idx)
        for vac_nr,vac_idx in enumerate(all_vac_idx):
            if al_mg_si_all[vac_nr][step,0] != 0:
                do_continue = test_anz_nn(al_mg_si_all,vac_nr,step,exit=False)
                if do_continue == True:
                    continue

            NN_1_indices, NN_2_indices = ase_get_neighborlist_1NN_2NN(frames[step],atomnr=vac_idx,cutoffa=cutoffa,cutoffb=cutoffb,skin=skin)
            NN_1_sym = [atom.symbol for atom in frames[step] if atom.index in NN_1_indices]
            NN_2_sym = [atom.symbol for atom in frames[step] if atom.index in NN_2_indices]
            NN_1_al = NN_1_sym.count("Al")
            NN_1_mg = NN_1_sym.count("Mg")
            NN_1_si = NN_1_sym.count("Si")
            NN_2_al = NN_2_sym.count("Al")
            NN_2_mg = NN_2_sym.count("Mg")
            NN_2_si = NN_2_sym.count("Si")
            al_mg_si_all[vac_nr][step,1] = NN_1_al
            al_mg_si_all[vac_nr][step,2] = NN_1_mg
            al_mg_si_all[vac_nr][step,3] = NN_1_si
            al_mg_si_all[vac_nr][step,4] = NN_2_al
            al_mg_si_all[vac_nr][step,5] = NN_2_mg
            al_mg_si_all[vac_nr][step,6] = NN_2_si
            al_mg_si_all[vac_nr][step,0] = step
            anz_1NN = np.sum(al_mg_si_all[vac_nr][step][1:4])
            anz_2NN = np.sum(al_mg_si_all[vac_nr][step][4:7])
            #str_1NN = str(NN_1_al).ljust(3)+str(NN_1_mg).ljust(3)+str(NN_1_si).ljust(3)
            #str_2NN = str(NN_2_al).ljust(3)+str(NN_2_mg).ljust(3)+str(NN_2_si).ljust(3)
            #print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, "||",str(NN_1_al).ljust(3),NN_1_mg,NN_1_si,"||",NN_2_al,NN_2_mg,NN_2_si)
            #print('-------',filename_analyze_all[vac_nr])
            #print(al_mg_si_all[vac_nr])
            if anz_1NN != 12 or anz_2NN != 6:
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'ERROR!')
                sys.exit("ERROR see above")
            print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'OK')
            #do_continue = test_anz_nn(al_mg_si_all,vac_nr,step,exit=True)

            if step > 0 and step in np.arange(structures)[::save_every]:
                np.savetxt(filename_analyze_all[vac_nr],al_mg_si_all[vac_nr],fmt='%i')
                print('saving',os.path.abspath(filename_analyze_all[vac_nr]),'at step',step)

    # save everything in the very end
    np.savetxt(filename_analyze_all[vac_nr],al_mg_si_all[vac_nr],fmt='%i')
    print('saving (very end)',os.path.abspath(filename_analyze_all[vac_nr]),'at step',step)
    return

def ase_get_neighborlist_1NN_2NN(frame,atomnr=0,cutoffa=3.,cutoffb=4.5,skin=0.1):
    NN_1_indices       = ase_get_neighborlist(frame,atomnr=atomnr,cutoff=cutoffa,skin=skin)
    NN_1_2_indices_tmp = ase_get_neighborlist(frame,atomnr=atomnr,cutoff=cutoffb,skin=skin)
    NN_2_indices       = np.sort(np.array(diff(NN_1_2_indices_tmp,NN_1_indices)))
    return NN_1_indices, NN_2_indices

def ase_get_neighborlist_1NN_2NN(frame,atomnr=0,cutoffa=3.,cutoffb=4.5,skin=0.1):
    NN_1_indices       = ase_get_neighborlist(frame,atomnr=atomnr,cutoff=cutoffa,skin=skin)
    NN_1_2_indices_tmp = ase_get_neighborlist(frame,atomnr=atomnr,cutoff=cutoffb,skin=skin)
    NN_2_indices       = np.sort(np.array(diff(NN_1_2_indices_tmp,NN_1_indices)))
    return NN_1_indices, NN_2_indices

def ase_get_unique_frames(frames):
    '''
    this function only takes care of exactly same frames;
    structures which are close by another function will be necessary;
    '''
    framesout = deepcopy(frames)
    length = len(frames)
    #for idx,midx in enumerate(tdqm(range(len(frames))[::-1])):
    for idx,midx in enumerate(range(len(frames))[::-1]):
        progress(idx, length , status='')
        isin=False
        if frames[midx] in frames[:midx]:
            isin=True
            del framesout[midx]
        #print(idx,midx,frames[midx].positions[0,0],isin)
    print('returning framesout ... (unique)')
    return framesout

def ase_enepot(atoms,units='eV',verbose=False):
    ''' units: eV, eV_pa, hartree, hartree_pa
        check before if calculator is attached
    '''

    #print('now in ene')
    #print('ac',atoms.cell)
    try:
        # in the case of "DFT"                 , it just retrieves the energy  -> get_stress() can NOT be obtained.
        # in the case of of an ace calculations, it calculates the energy      -> get_stress() CAN     be obtained.
        #print('before')
        ene = atoms.get_potential_energy()
        #print('ene:',ene,atoms.info)
        #uuid = atoms.info["comment"]
        #print('uuid',uuid)
        #stress = atoms.get_stress()
        #print('stress:',stress)
    except: # RuntimeError:
        #print("had runtime error, e.g. ther cant be an energy in the POSCAR")
        ene = 0.
        #stress = False
    if verbose > 1:
        print('ene eV',ene,"(not per atom)")
    units_split = units.split("_")
    #print('us',units_split,units_split[1])
    if units_split[0].lower() == 'ev':
        pass
    elif units_split[0].lower() == 'mev':
        ene = ene*1000.
    elif units_split[0] == "hartree" or units_split[0] == "Hartree":
        ene = ene/aseunits.Hartree

    if len(units_split) == 2:
        if units_split[1] == 'pa':
            ene = ene/atoms.get_number_of_atoms()
        else:
            sys.exit("energy can not have this units (ending must be pa, eV_pa or hartree_pa)")

    return ene

def ase_get_chemical_symbols_to_number_of_species(atoms):
    symbols = atoms.get_chemical_symbols()
    #numat = atoms.get_number_of_atoms()
    #print('symbols',symbols)
    #print('numat',numat)

    uniquesym = set(atoms.get_chemical_symbols())
    d = {}
    for i in uniquesym:
        #print(i,symbols.count(i),numat)
        d[i] = symbols.count(i)

    def dcheck(element):
        if element in d.keys():
            #print(element+" exists")
            pass
        else:
            d[element] = 0

    dcheck("Mg")
    dcheck("Si")
    dcheck("Al")
    return d

def ase_get_chemical_symbols_to_conz(atoms):
    symbols = atoms.get_chemical_symbols()
    numat = atoms.get_number_of_atoms()
    #print('symbols',symbols)
    #print('numat',numat)

    uniquesym = set(atoms.get_chemical_symbols())
    # print("uniquesym",uniquesym) # --> set(['Mg', 'Al'])
    d = {}
    for i in uniquesym:
        #print(i,symbols.count(i),numat)
        #print('ii',i,'---',float(symbols.count(i))/float(numat))
        d[i] = float(symbols.count(i))/float(numat)
    #print('kk',d)
    def dcheck(element):
        if element in d.keys():
            #print(element+" exists")
            pass
        else:
            d[element] = 0.0
    #print('dd',d)
    dcheck("Mg")
    dcheck("Si")
    dcheck("Al")
    return d

class ase_calculate_ene( object ):
    '''
    - should be evaluated just once to initialize and NOT FORE EVERY CALCULATION
      (due to the fact that mypot would be called unnecessarily)
    ase_calculate_ene (ace) class which holds lammps commands to be executed
    if only pot is defined, static calculation.
    '''
    def __init__(self,
            pot,
            potpath,
            use_different_epoch=False,
            units=False,
            geopt=False,
            kmc=False,
            verbose=False,
            temp=False,
            elastic=False,
            ):

        #self.pot = pot
        self.LAMMPS_COMMAND = False
        self.potpath        = potpath
        self.mypot          = False
        self.units          = units.lower()
        self.geopt          = geopt          # so far only for ene object.
        self.elastic        = elastic
        self.elastic_relax  = True
        self.nsteps         = 0
        self.verbose        = verbose
        self.atoms          = False          # ase atoms object (frame)
        if self.verbose:
            print('>> ase_calculate_ene: initializing mypot .... to self.pot')
        self.pot            = mypot(pot,self.potpath,use_different_epoch = use_different_epoch,verbose = self.verbose)

        #####################
        # for the calculator
        #####################
        self.calculator = "lammps"
        self.lmpcmd     = False         # in case we run through ase (also needs lmpcmd) or external lammps
        self.atom_types = False     # for nn pot
        self.keep_alive = True

        #print('init')
        # case of MD or KMC
        self.kmc  = kmc
        self.temp = temp

        #self.eos = [ False, False, False, False] # e0_meV/pa, v0_ang^3/pa, B0, B0der]
        return


    def print_variables_ase(self,text=""):
        if self.verbose > 1:
            tt = 'ase_calculate_ene.'
            print()
            print(text,tt+'pot.pot      (1) :',self.pot.pot)    # : n2p2_v2ag
            print(text,tt+'units        (1) :',self.units)      # : ev
            print(text,tt+'geopt        (1) :',self.geopt)      # : False
            print(text,tt+'elastic      (1) :',self.elastic)    # : False
            print(text,tt+'nsteps       (1) :',self.nsteps)     # : 0
            print(text,tt+'lmpcmd       (1) :',self.lmpcmd)     # : False
            print(text,tt+'kmc          (1) :',self.kmc)        # : False
            print(text,tt+'temp         (1) :',self.temp)       # : False
            print(text,tt+'verbose      (1) :',self.verbose)       # : False
            print(text,tt+'LAMMPS_COMMAND1) :',self.LAMMPS_COMMAND)    # : n2p2_v2ag
            print()
        return


    def lammps_command_potential_n2p2(self):
        #units_giulio_ene = "0.0367493254"
        ase_units_ene    = "0.03674932247495664" # 1./ase.units.Hartree
        #units_giulio_bohr = "1.8897261328"
        ase_units_bohr    = "1.8897261258369282" # 1./ase.units.Bohr

        command = [
        # showewsum 1 showew yes resetew no maxew 1000000
        'variable nnpDir string \"'+self.pot.potpath_work+'\"',
        "pair_style nnp dir ${nnpDir} showew no resetew yes maxew 100000000 cflength "+ase_units_bohr+" cfenergy "+ase_units_ene,
        "pair_coeff * * "+str(self.pot.potcutoff),
        #"#write_data ./pos.data # would this be the final struct?"
        ]
        return command

    def lammps_command_potential_runner(self):
        command = [
        # comment
        #"# thermo 1 # for geopt",
        'variable nnpDir string \"'+self.pot.potpath_work+'\"',
        "pair_style runner dir ${nnpDir} showewsum 1 showew yes resetew no maxew 1000000",
        #"# pair_coeff * * 7.937658735",
        #"pair_coeff * *  14.937658735"
        "pair_coeff * * "+str(self.pot.potcutoff)
        ]
        return command

    def lammps_command_masses(self):
        command = [
                "mass 1 24.305",
                "mass 2 26.9815385",
                "mass 3 28.0855",
                ]
        return command

    def pot_get_and_ase_lmp_cmd(self,kmc=False,temp=False,nsteps=0,ffsocket='inet',address=False):
        ''' geoopt (geometry optimization) is added / or not in
            lammps_write_inputfile(); here only the potential is set.
            ffsocket: ipi ffsocket [ "unix" or "inet" ]
        '''
        if self.verbose > 1:
            print('PPh potDONE:',self.pot.potDONE)
        if self.pot.potDONE == False:
            if self.verbose:
                print("PPi self.pot.get()")
            self.pot.get()

        self.kmc = kmc
        self.temp = temp
        self.nsteps = nsteps
        self.ffsocket = ffsocket
        if self.ffsocket not in [ "unix", "inet" ]:
            print('ffsocket:',ffsocket)
            sys.exit('ffsocket has to be "unix" or "inet"; Exit!')



        if self.verbose > 2:
            tt = 'PPj pot_get_and_ase_lmp_cmd_A '
            print(tt+'pot.pot           :',self.pot.pot)
            print(tt+'pot.potpath       :',self.pot.potpath)
            print(tt+'pot.potpath_work  :',self.pot.potpath_work)
            print(tt+'pot.pottype       :',self.pot.pottype)
            print()

        #sys.exit()
        # this depends only on the potential which is already defined
        # so should be easy to make this general.
        self.lmpcmd = [ "########## lmpcmd.begin #############" ]
        self.lmpcmd = self.lmpcmd + self.lammps_command_masses()

        if self.pot.pottype == "n2p2":
            # showewsum 1 showew yes resetew no maxew 1000000
            self.lmpcmd = self.lmpcmd + self.lammps_command_potential_n2p2()
            self.atom_types = {'Mg':1,'Al':2,'Si':3}

        elif self.pot.pottype == "runner":
            self.lmpcmd = self.lmpcmd + self.lammps_command_potential_runner()
            self.atom_types = {'Mg':1,'Al':2,'Si':3}
        else:
            sys.exit('pot '+str(self.pot.pot)+' not found! (X)')

        if self.kmc:
            if self.ffsocket == "unix": add = "unix"
            if self.ffsocket == "inet": add = ""
            if address == False:
                address = gethostname()
            self.lmpcmd = self.lmpcmd + [
                "",
                "timestep 0.001   # timestep (ps)",
                "velocity all create "+str(self.temp)+" 4928459",  # create initial velocities 4928459 is random seed for velocity initialization"
                "thermo 1   # screen output interval (timesteps)",
                "fix 1 all ipi "+str(address)+" 12345 "+str(add),
                ]
                # with n2p2 in the parallel version, inet is not working
                # "fix 1 all ipi fidis 12345",     # for fidis job
                # "fix 1 all ipi mac 77776 unix",  # for mac job

        self.lmpcmd = self.lmpcmd + [ "########## lmpcmd.end  #############" ]
        if self.verbose > 1:
            print('HERE THE lmpcmd I got',self.lmpcmd)
        self.print_variables_ase("pot_get_and_ase_lmp_cmd_FIN")
        self.LAMMPS_COMMAND = get_LAMMPS_executable(exit=True) #,verbose=self.verbose)
        return

    def define_wrapped_self_atoms(self,atoms=False):
        if atoms == False:
            atoms = self.atoms
        else:
            atoms = atoms
        atoms.wrap()
        #self.atomsin = deepcopy(atoms)
        self.atoms = atoms
        return self.atoms

    def ene_allfix(self,atoms=False):
        atoms = self.define_wrapped_self_atoms(atoms)
        keep_alive = False
        asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        atoms.set_calculator(asecalcLAMMPS)
        ene = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        return ene

    def ene_new(self,atoms=False,
            atomrelax=False,
            volumerelax=False,
            cellshaperelax=False,
            print_minimization_to_screen=False,minimizer="LGBFGS"):
        ''' atoms is an ase object
            if don_change_atomsobject is chosen,
        '''
        ## now the atoms object is not changed
        #atoms = atomsin.copy()
        atoms = self.define_wrapped_self_atoms(atoms)

        if atomrelax == False and self.geopt == False:
            pass
        if atomrelax == True and self.geopt == True:
            pass
        if atomrelax == True and self.geopt == False:
            pass  # realx wins, since locally set explicitely
        if atomrelax == False and self.geopt == True:
            atomrelax = True
        #print('svb',self.verbose)
        if self.verbose > 2:
            print()
            print("#####################################################")
            print('##--lmpcmd:',self.lmpcmd)
            print('##--atom_types',self.atom_types)
            print('##--geopt',self.geopt)
            print('##--atomrelax',atomrelax)
            print('##--minimizer',minimizer)
            print("#####################################################")
            show_ase_atoms_content(atoms,showfirst=10,comment = "START ASE INTERNAL CALUCLATION !!!")
            print()

        ################################
        ## set up the calculator
        ## (by attaching the calculator to atoms)
        ## this is valid for relaxations and static calcs
        ################################
        keep_alive = False
        if atomrelax == False: keep_alive = False
        if atomrelax == True:  keep_alive = True
        asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        atoms.set_calculator(asecalcLAMMPS)

        ### attach to atoms to relax the cell
        constraint = False
        if volumerelax == False and atomrelax == False and cellshaperelax == False:
            constraint = False
        elif volumerelax == True and atomrelax == False and cellshaperelax == False:
            print('before murn')
            vinet = self.get_murn(atoms,verbose=False,return_minimum_volume_frame=True)
            print('after murn')
            #constraint = ExpCellFilter(atoms, hydrostatic_strain=True) does not work
        elif atomrelax == False and cellshaperelax == True:
            constraint = StrainFilter(atoms)  # this relaxes the cell shape & the volume while keeping atomic positions fixed
        else:
            sys.exit('still to define')

        #if cellrelax == True and atomrelax == False:
        #    constraint = StrainFilter(atoms)  # this relaxes the cell shape & the volume while keeping atomic positions fixed
        #    ## in this case it does not work out
        #    sys.exit('This gives a segmentation fault (coredump) when cellrelax == True and atomrelax == True ... in this case do only one, then the other one')
        #elif cellrelax == True and atomrelax == True:
        #    ## when doint both this is recommended
        #    constraint = ExpCellFilter(atoms)  # Modify the supercell and the atom positions.

        ## atomrelax = False and cellrelax = False works
        ## atomrelax = True  and cellrelax = False works
        ## atomrelax = True  and cellrelax = True  works
        ## atomrelax = False and cellrelax = True  NOPE
        if atomrelax == True or cellshaperelax == True or volumerelax == True:
            #if atomrelax == False and cellrelax == True:
            #    print('NOOOOOOOOOOOOOOOWWWWWWWWWWWWWWWWW'*3)
            ################################################
            ### with geometry optimization
            ################################################
            # in case of a verbose run:
            #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, log_file='./xlolg.lammps.log',tmp_dir="./",keep_alive=True,atom_types=self.atom_types)

            # in case  of a non verbose run
            #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
            #atoms.set_calculator(asecalcLAMMPS)
            #from ase.io.trajectory import Trajectory
            #traj = Trajectory('ka', mode='w',atoms=atoms)
            #opt = BFGS(atoms,trajectory="ni.traj")
            #opt.run(steps=20)
            minimizer_choices = [ 'BFGS', 'LGBFGS', 'FIRE', 'GPMin', 'bh', 'mh' ]
            if minimizer not in minimizer_choices:
                print("your minimizer",minimizer)
                print("available:",minimizer_choices)
                sys.exit("choose one of the proper minimizer choices")
            if minimizer == 'mh' and atomrelax == True and cellrelax == True:
                sys.exit("use minimahopping (mh) with either with atomrelax or with cellrelax; or use both but with other minimizer;  but not with both at a time, rather run those sequentially o")

            if print_minimization_to_screen:
                print("AAA nat:",atoms.get_number_of_atoms())
                print("AAA pos:",atoms.get_positions()[:4])
                print("AAA for:",atoms.get_forces()[:4])
                print("AAA fmx:",abs(atoms.get_forces()).max())
                print("AAA vol:",atoms.get_volume())
                print("AAA vpa:",atoms.get_volume()/atoms.get_number_of_atoms())

            ## the syntax apparently vareis a bit depending
            ## if constraint or not
            use = atoms
            if type(constraint) != bool:
                use = constraint

            if print_minimization_to_screen:
                print("minimizer:",minimizer)
                logfile="-" # output to screen
            else:
                logfile="tmp" # output to file and not to screen


            #if print_minimization_to_screen:
            #print('BBB logfile',logfile)

            if minimizer == 'BFGS':
                opt1 = BFGS(use,logfile=logfile) #,trajectory="ni.traj")
            elif minimizer == 'LGBFGS':
                opt1 = LBFGS(use,logfile=logfile) #,trajectory="test.traj")
            elif minimizer == 'GPMin':
                opt1 = GPMin(use,logfile=logfile) #,trajectory="test.traj")
            elif minimizer == 'FIRE':
                opt1 = FIRE(use,logfile=logfile) #,trajectory="ni.traj")
            elif minimizer == 'bh':
                kB = 1.38064852e-23
                kB = 1.6021765e-19
                opt1 = BasinHopping(atoms=use, # the system to optimize
                  temperature=1*kB, # 'temperature' to overcome barriers
                  dr=0.5,      # maximal stepwidth
                  optimizer=LBFGS, # optimizer to find local minima
                  fmax=0.1,      # maximal force for the optimizer
                  logfile=logfile)
            elif minimizer == 'mh':
                if os.path.isfile("tmp"):
                    os.remove("tmp")
                opt1 = MinimaHopping(atoms=use,logfile=logfile)
                opt1(totalsteps=10)
            #print('startrun....')
            maxsteps = 200
            if minimizer == 'bh':
                maxsteps = 3

            ######################################################
            ## MINIMIZE
            ######################################################
            if minimizer not in ['mh','bh']: # in all cases but
                #opt1.run(steps=maxsteps,fmax=0.005)
                #opt1.run(steps=maxsteps,fmax=0.0001)
                #print('start')
                opt1.run(fmax=0.0001)
                #print('maxsteps                ',maxsteps,type(maxsteps))
                #print('opt1.get_number_of_steps',opt1.get_number_of_steps(),type(opt1.get_number_of_steps()))

                if maxsteps == opt1.get_number_of_steps():
                    print('DID NOT CONVErGE IN '+str(maxsteps)+' number of minimizer steps!')
                    #print('---- cell -----')
                    #print(atoms.get_cell())
                    #print('---- positions -----')
                    #print(atoms.get_positions())
                    return np.nan

        if print_minimization_to_screen:
            print('UUU atomrelax:',atomrelax)
            print('UUU cellrelax:',cellrelax)
            print("UUU nat:",atoms.get_number_of_atoms())
            print("UUU pos[:4]:",atoms.get_positions()[:4])
            print("UUU for[:4]:",atoms.get_forces()[:4])
            print("UUU fmx:",abs(atoms.get_forces()).max())
            print("UUU vol:",atoms.get_volume())
            print("UUU vpa:",atoms.get_volume()/atoms.get_number_of_atoms())

        if self.verbose > 1:
            print('ZZ done2')
            print('ZZ self.units',self.units)
        ######################################################
        # calculate the energy
        ######################################################
        #print('atxxx',atoms)
        ene = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        if print_minimization_to_screen:
            print('atoms')
            print(atoms.get_positions()[:3])
        if self.verbose > 1:
            print('ZZ ene:',ene,self.units)
        if not print_minimization_to_screen and os.path.isfile("tmp"):
            os.remove("tmp")
        #sys.exit()
        #ene = atoms.get_total_energy()
        #if self.verbose:
        #    print('ene',ene)
        #return ene,ene/atoms.get_number_of_atoms()*1000.
        if self.verbose > 2:
            show_ase_atoms_content(atoms,showfirst=10,comment="FINISHED ASE INTERNAL CALUCLATION")
            print()
            print()

        #print('forces out',atomrelax,cellrelax)
        #print(atoms.get_forces()[:3])
        return ene

    def ene(self,atoms=False,
            atomrelax=False,
            cellrelax=False,
            cellshaperelax=False,
            cellvolumerelax=False,
            print_minimization_to_screen=False,
            minimizer="LGBFGS",
            debug=False):
        ''' atoms is an ase object
            if don_change_atomsobject is chosen,
        '''
        if debug:
            print_minimization_to_screen=True
            print('777 degub is on for calculation of ene')
        #unique_elements = list(set(atoms.get_chemical_symbols()))
        ## now the atoms object is not changed
        #atoms = atomsin.copy()
        atoms = self.define_wrapped_self_atoms(atoms)

        if atomrelax == False and self.geopt == False:
            pass
        if atomrelax == True and self.geopt == True:
            pass
        if atomrelax == True and self.geopt == False:
            pass  # realx wins, since locally set explicitely
        if atomrelax == False and self.geopt == True:
            atomrelax = True
        #print('svb',self.verbose)
        if self.verbose > 2 or debug:
            print()
            print("#####################################################")
            print('##--lmpcmd:',self.lmpcmd)
            print('##--atom_types',self.atom_types)
            print('##--geopt',self.geopt)
            print('##--atomrelax',atomrelax)
            print('##--minimizer',minimizer)
            print("#####################################################")
            show_ase_atoms_content(atoms,showfirst=10,comment = "START ASE INTERNAL CALUCLATION !!!")
            print()

        ################################
        ## set up the calculator
        ## (by attaching the calculator to atoms)
        ## this is valid for relaxations and static calcs
        ################################
        if atomrelax == False: keep_alive = False
        if atomrelax == True:  keep_alive = True
        self.keep_alive = keep_alive
        if debug:
            print('ATTACHING CALCULATOR!!')
        self.get_calculator(atoms)

        ### attach to atoms to relax the cell
        if debug:
            print('attaching constraints')
        constraint = False
        if cellrelax == True: # and atomrelax == False:
            constraint = StrainFilter(atoms)  # this relaxes the cell shape & the volume while keeping atomic positions fixed
            ## in this case it does not work out
            #sys.exit('This gives a segmentation fault (coredump) when cellrelax == True and atomrelax == True ... in this case do only one, then the other one')
        elif cellrelax == True and atomrelax == True:
            ## when doint both this is recommended
            constraint = ExpCellFilter(atoms)  # Modify the supercell and the atom positions.

        ## atomrelax = False and cellrelax = False works
        ## atomrelax = True  and cellrelax = False works
        ## atomrelax = True  and cellrelax = True  works
        ## atomrelax = False and cellrelax = True  NOPE
        if atomrelax == True or cellrelax == True:
            if atomrelax == False and cellrelax == True:
                print('NOOOOOOOOOOOOOOOWWWWWWWWWWWWWWWWW'*3)
            ################################################
            ### with geometry optimization
            ################################################
            # in case of a verbose run:
            #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, log_file='./xlolg.lammps.log',tmp_dir="./",keep_alive=True,atom_types=self.atom_types)

            # in case  of a non verbose run
            #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
            #atoms.set_calculator(asecalcLAMMPS)
            #from ase.io.trajectory import Trajectory
            #traj = Trajectory('ka', mode='w',atoms=atoms)
            #opt = BFGS(atoms,trajectory="ni.traj")
            #opt.run(steps=20)
            minimizer_choices = [ 'BFGS', 'LGBFGS', 'FIRE', 'GPMin', 'bh', 'mh' ]
            if minimizer not in minimizer_choices:
                print("your minimizer",minimizer)
                print("available:",minimizer_choices)
                sys.exit("choose one of the proper minimizer choices")
            if minimizer == 'mh' and atomrelax == True and cellrelax == True:
                sys.exit("use minimahopping (mh) with either with atomrelax or with cellrelax; or use both but with other minimizer;  but not with both at a time, rather run those sequentially o")

            if print_minimization_to_screen:
                print("AAA nat:",atoms.get_number_of_atoms())
                print("AAA pos:",atoms.get_positions()[:4])
                print("AAA for:",atoms.get_forces()[:4])
                print("AAA fmx:",abs(atoms.get_forces()).max())
                print("AAA vol:",atoms.get_volume())
                print("AAA vpa:",atoms.get_volume()/atoms.get_number_of_atoms())

            ## the syntax apparently vareis a bit depending
            ## if constraint or not
            atoms_or_constraint = atoms
            if type(constraint) != bool:
                atoms_or_constraint = constraint

            if print_minimization_to_screen:
                print("minimizer:",minimizer)
                logfile="-" # output to screen
            else:
                logfile="tmp" # output to file and not to screen


            #if print_minimization_to_screen:
            #print('BBB logfile',logfile)

            if minimizer == 'BFGS':
                opt1 = BFGS(atoms_or_constraint,logfile=logfile) #,trajectory="ni.traj")
            elif minimizer == 'LGBFGS':
                opt1 = LBFGS(atoms_or_constraint,logfile=logfile) #,trajectory="test.traj")
            elif minimizer == 'GPMin':
                opt1 = GPMin(atoms_or_constraint,logfile=logfile) #,trajectory="test.traj")
            elif minimizer == 'FIRE':
                opt1 = FIRE(atoms_or_constraint,logfile=logfile) #,trajectory="ni.traj")
            elif minimizer == 'bh':
                kB = 1.38064852e-23
                kB = 1.6021765e-19
                opt1 = BasinHopping(atoms=atoms_or_constraint, # the system to optimize
                  temperature=1*kB, # 'temperature' to overcome barriers
                  dr=0.5,      # maximal stepwidth
                  optimizer=LBFGS, # optimizer to find local minima
                  fmax=0.1,      # maximal force for the optimizer
                  logfile=logfile)
            elif minimizer == 'mh':
                if os.path.isfile("tmp"):
                    os.remove("tmp")
                opt1 = MinimaHopping(atoms=atoms_or_constraint,logfile=logfile)
                opt1(totalsteps=10)
            #print('startrun....')
            maxsteps = 200
            if minimizer == 'bh':
                maxsteps = 3

            ######################################################
            ## MINIMIZE
            ######################################################
            if minimizer not in ['mh','bh']: # in all cases but
                #opt1.run(steps=maxsteps,fmax=0.005)
                #opt1.run(steps=maxsteps,fmax=0.0001)
                #print('start')
                opt1.run(fmax=0.0001)
                #print('maxsteps                ',maxsteps,type(maxsteps))
                #print('opt1.get_number_of_steps',opt1.get_number_of_steps(),type(opt1.get_number_of_steps()))

                if maxsteps == opt1.get_number_of_steps():
                    print('DID NOT CONVErGE IN '+str(maxsteps)+' number of minimizer steps!')
                    #print('---- cell -----')
                    #print(atoms.get_cell())
                    #print('---- positions -----')
                    #print(atoms.get_positions())
                    return np.nan

        if print_minimization_to_screen:
            print('UUA atomrelax:',atomrelax)
            print('UUA cellrelax:',cellrelax)
            print("UUA nat:",atoms.get_number_of_atoms())
            print("UUA pos[:4]:",atoms.get_positions()[:4])
            #print('for',atoms.get_forces.__module__)
            #print('for',atoms.get_forces.__globals__)
            print("UUA for:",atoms.get_forces()[:4])
            print("UUA fmx:",abs(atoms.get_forces()).max())
            print("UUA vol:",atoms.get_volume())
            print("UUA vpa:",atoms.get_volume()/atoms.get_number_of_atoms())

        if self.verbose > 1:
            print('ZZ done243')

        if self.verbose > 1:
            print('ZZ done236')
            print('ZZ self.units',self.units)
        ######################################################
        # calculate the energy
        ######################################################
        #print('atxxx',atoms)
        ene = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        #print('jo')
        if print_minimization_to_screen:
            print('atoms')
            print(atoms.get_positions()[:3])
        if self.verbose > 1:
            print('ZZ ene:',ene,self.units)
        if not print_minimization_to_screen and os.path.isfile("tmp"):
            os.remove("tmp")
        #sys.exit()
        #ene = atoms.get_total_energy()
        #if self.verbose:
        #    print('ene',ene)
        #return ene,ene/atoms.get_number_of_atoms()*1000.
        if self.verbose > 2:
            show_ase_atoms_content(atoms,showfirst=10,comment="FINISHED ASE INTERNAL CALUCLATION")
            print()
            print()

        #print('forces out',atomrelax,cellrelax)
        #print(atoms.get_forces()[:3])
        return ene

    def stress(self,atoms=False):
        atoms = self.define_wrapped_self_atoms(atoms)
        try:
            stress = atoms.get_stress()
        except:
            stress = False
        #print('stress',stress)
        return stress

    def get_v0(self,atomsin=False):
        ''' the function will never change the atomsobject '''
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        atoms_murn = atomsin.copy()
        atoms_murn.wrap()

        keep_alive = False
        atomrelax = False
        if atomrelax == False: keep_alive = False
        if atomrelax == True:  keep_alive = True
        asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        atoms_murn.set_calculator(asecalcLAMMPS)

        ### relax the atoms_murn first to the equilibrium
        self.ene(atoms_murn,cellrelax=True,atomrelax=True)
        return atoms_murn.get_volume()

    def get_v0_pa(self,atomsin=False):
        ''' the function will never change the atomsobject '''
        return self.get_v0(atomsin=atomsin)/atomsin.get_number_of_atoms()

    def get_elastic_external(self,atomsin=False,verbose=False,text=False,get_all_constants=False):
        ''' the function will never change the atomsobject '''
        #print('######## get_elastic_external #############')
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        frame = atomsin.copy()
        frame.wrap()
        #print('3')

        self.get_calculator(frame)  # to be able to calculate stress
        #print('4x')
        #print('stress original frame:',frame.get_stress())
        #print('stress original frame:',frame.get_stress())
        #print('5x')
        if verbose:
            print('frame cell::',frame.get_cell())
            print('verbose   ::',verbose)
        if self.elastic_relax == True:
            self.ase_relax_cellshape_and_volume_only(frame,verbose=verbose)
        #print('stress !relaxed! frame :',frame.get_stress())
        #print('volume !relaxed! frame :',frame.get_volume())
        #print('ene    !relaxed! frame :',frame.get_potential_energy())
        print('### get_elastic_external: !RELAXED! stress(max),vol,ene',abs(frame.get_stress()).max(),frame.get_volume(),frame.get_potential_energy())
        #print('frame cell',frame.get_cell())
        #ase_write('pos.runner',frame,format='runner')
        #print('-------------- lammps_ext_calc -----------')
        ene_pot_lmp = lammps_ext_calc(frame,self,get_elastic_constants=get_all_constants)
        #print('ene_pot_lmp...kk',ene_pot_lmp)
        #sys.exit('88')

        if verbose:
            if type(text) != bool:
                text = printred(text)
                #if get_all_constants == True:
                ene_pot_lmp = ene_pot_lmp.replace('Elastic Constant ', 'Elastic Constant '+text+" ")
                #else:
                #    ene_pot_lmp = ene_pot_lmp.replace('Elastic Constant ', 'Elastic Constant '+text+" ")
                #ene_pot_lmp = ene_pot_lmp.replace('C44all =', printred('C44all ='))
            print("relaxed? "+str(self.elastic_relax)+";",ene_pot_lmp)
        return

    def get_elastic(self,atomsin=False,verbose=False):
        ''' the function will never change the atomsobject '''

        if atomsin == False:
            sys.exit('need to define atoms in this case XX')

        atoms_h = atomsin.copy()
        atoms_h.wrap()

        keep_alive = False
        atomrelax = False
        if atomrelax == False: keep_alive = False
        if atomrelax == True:  keep_alive = True
        #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        #atoms_h.set_calculator(asecalcLAMMPS)  # wird in parcalc.py gesetzt
        #/home/glensk/miniconda2/lib/python2.7/site-packages/parcalc/parcalc.py


        #### load the elastic stuff
        #from parcalc import ClusterVasp, ParCalculate
        #from elastic import get_pressure, BMEOS, get_strain
        #from elastic import get_elementary_deformations, scan_volumes
        #from elastic import get_BM_EOS, get_elastic_tensor

        print('############## from elastic ################')
        if False:
            print('ene   ',self.ene(atoms_h))
            print('stress1',atoms_h.get_stress())
        self.ase_relax_cellshape_and_volume_only(atoms_h,verbose=False)
        if False:
            print('stress2',atoms_h.get_stress())
            print('cell',atoms_h.get_cell())

        cell_ref = (atoms_h.copy()).get_cell()
        atoms_work = atoms_h.copy()
        self.ase_relax_cellshape_and_volume_only(atoms_work,verbose=False)

        for i in [0.98,0.99,1.00,1.01, 1.02, 1.03]:
            cell_work = cell_ref.copy()
            #print('cell_ref',atoms_h.get_cell())
            cell_work[0,0] = cell_ref[0,0]*i
            #print('cell_ref',cell_ref)
            atoms_work.set_cell(cell_work,scale_atoms=True)
            #print('cell_ref?',atoms_h.get_cell())
            if False:
                print('strain',i,'stress?',atoms_work.get_stress())
            #print('cell_ref',cell_ref)



        ################################################################
        # from elastic
        # http://wolf.ifj.edu.pl/elastic/lib-usage.html
        ################################################################
        try:
            from elastic.elastic import get_cart_deformed_cell, get_lattice_type, get_elementary_deformations
        except ImportError:
            return
        from elastic import get_pressure, BMEOS, get_strain
        from elastic import get_BM_EOS, get_elastic_tensor
        from parcalc import ParCalculate

        sym = get_lattice_type(atoms_h)
        print('sym',sym)
        # Create 10 deformation points on the a axis
        systems = []
        ss=[]
        for d in np.linspace(-0.2,0.2,11):
            # get_cart_deformed_cell:
            # The axis is specified as follows: 0,1,2 = x,y,z ;
            # sheers: 3,4,5 = yz, xz, xy.
            # d: The size of the deformation is in percent and degrees, respectively.
            struct = get_cart_deformed_cell(atoms_h, axis=0, size=d)
            if verbose:
                print()
            strc = struct.get_cell()
            if verbose:
                print('d      ',d)
                print('struct :',strc[0],strc[1],strc[2])
            strca = atoms_h.get_cell()
            if verbose:
                print('atoms_h:',strca[0],strca[1],strca[2])
            stress = struct.get_stress()
            strain = get_strain(struct, atoms_h)
            pressure = get_pressure(stress)
            if verbose:
                print("stress :",stress)
                print("strain :",strain)
                print('pressure:',pressure)
            ss.append([strain, stress])
            systems.append(struct)

        def myparcalc():
            systems_all = get_elementary_deformations(atoms_h, n=5, d=0.33)
            if type(systems_all) != type([]) :
                sysl=[systems_all]
                #print('11')
            else:
                sysl=systems_all
                #print('22')

            res = []
            for n,s in enumerate(sysl):
                s.get_potential_energy()
                res.append([n,s])
            return [r for ns,s in enumerate(sysl) for nr,r in res if nr==ns]
        res = myparcalc()
        Cij, Bij = get_elastic_tensor(atoms_h, systems=res)
        print("Cij (GPa):", Cij/aseunits.GPa)

        ss=np.array(ss)
        lo=min(ss[:,0,0])
        hi=max(ss[:,0,0])
        mi=(lo+hi)/2
        wi=(hi-lo)/2
        xa=np.linspace(mi-1.1*wi,mi+1.1*wi, 50)

        # Now fit the polynomials to the data to get elastic constants
        # C11 component
        f=np.polyfit(ss[:,0,0],ss[:,1,0],3)
        c11=f[-2]/aseunits.GPa
        #print('ffff')
        #print(f)
        #print()
        #print(f[-2])

        # C12 component
        f=np.polyfit(ss[:,0,0],ss[:,1,1],3)
        c12=f[-2]/aseunits.GPa

        #np.savetxt('c11.dat',np.transpose([ss[:,0,0],ss[:,1,0]]))
        print('C11 = %.3f GPa, C12 = %.3f GPa => K= %.3f GPa' % (
                    c11, c12, (c11+2*c12)/3))

        ################################################################
        # daniels manual way
        ################################################################
        from daniel_strainstuff import _gen_strainmatrix, _apply_strain
        from daniel_strainstuff import _2lammpslattice
        from daniel_lmprun import _find_compliance_viaenergy
        from daniel_lmprun import run_fcc, _print_compliance_components
        from lammps import lammps
        lmp = lammps()

        print('############## daniels way  ################')
        for d in np.linspace(-0.2,0.2,5):
            struct = get_cart_deformed_cell(atoms_h, axis=3, size=d)
            if verbose:
                print()
            strc = struct.get_cell()
            if verbose:
                print('d      ',d)
                print('struct :',strc[0],strc[1],strc[2])
            strca = atoms_h.get_cell()
            if verbose:
                print('atoms_h:',strca[0],strca[1],strca[2])
            stress = struct.get_stress()
            strain = get_strain(struct, atoms_h)
            pressure = get_pressure(stress)
            if verbose:
                print("stress :",stress)
                print("strain :",strain)

        def my_get_cart_deformed_cell(base_cryst, size=1,verbose=False,vol=False):
            from ase.atoms import Atoms
            cryst = Atoms(base_cryst)
            uc = base_cryst.get_cell()
            if vol != False:
                uc = base_cryst.get_cell()*vol
            s = size/100.0
            L = np.diag(np.ones(3))
            #L = L * 0.9997
            if verbose:
                print(L)
                print()
            if False: # not volume conserving
                L[1, 2] += s/2.
                L[2, 1] += s/2.
            if True: # volume conserving
                L[0, 1] += s/2.
                L[1, 0] += s/2.
                L[2, 2] += (s**2.)/(4.-s**2.)
            if verbose:
                print(L)
                print()
            uc = np.dot(uc, L)
            cryst.set_cell(uc, scale_atoms=True)
            return cryst


        print()
        print("########### now only one deformed cell ###########")
        print('atoms_h.get_cell()     :')
        print(atoms_h.get_cell())
        print('atoms_h.get_stress()   :')
        print(atoms_h.get_stress())
        volfact = 1.0111  # C44 36.86712322917455  vol (68.43265399591313) = 17.108163499  d = 0.557 ang^3
        volfact = 1.0051  # C44 39.770796753113444 vol (67.22160401031846) = 16.805401003  d = 0.254 ang^3
        volfact = 1.0011  # C44 41.41759291379788  vol (66.42222758795769) = 16.605556897  d = 0.055 ang^3
        volfact = 1.0001  # C44 41.793101000472575 vol (66.2233786205119)  = 16.555844655  d = 0.005 ang^3
        volfact = 1.0000  # C44 41.82985611457602  vol (66.20351557966632) = 16.550878895  d = 0.000 ang^3
        volfact = 0.995477787394 # C44 43.3418676  vol (65.30941200550778) = 16.550878895  d = 0.000 ang^3
        volfact = 0.             # C44 42.3672012  vol (65.90412920697601) = 16.476032301  d = 0.075 ang^3

        # DFT                                      vol (65.905655194)      = 16.476413798491  (first converged)
        # DFT                                      vol (65.904129207)      = 16.476032301744  (beset converge)
        V0DFT = 16.476413798491 #  first converged
        V0DFT = 16.476032301744 #  beset converge
        a0DFT = (4.*V0DFT)**(1./3.)
        volfact = ((V0DFT*4.)/66.20351557966632)
        print('volfact',volfact)

        #cryst.set_cell(np.diag(np.ones(3))*a0, scale_atoms=True)
        #cell = cryst.get_cell()
        #cryst = my_get_cart_deformed_cell(cryst, size=0.2)


        ### This is to be at volume of DFT! (JUST IF YOU WANT TO CHECK HWO LARGE THE
        ### ERROR WOULD BE IF WE USE THIS (DFT) VOLUME
        if False:
            if volfact >= 1.0:
                atoms_h.set_cell(atoms_h.get_cell()*volfact, scale_atoms=True)
            else:
                atoms_h.set_cell(np.diag(np.ones(3))*a0DFT, scale_atoms=True)

        e0 = atoms_h.get_potential_energy()
        V0 = atoms_h.get_volume()
        print('atoms_h.get_cell()     :')
        print(atoms_h.get_cell())
        print('atoms_h.get_potential():',e0)
        print('atoms_h.V0',V0)
        print()
        print('volum0',V0,'   s 0               ene0: e0',e0)
        print('----------------------------------------------------------------------------------------')
        #for d in np.linspace(-0.1,0.1,4):
        points = 10
        ene_vs_strain = np.zeros((points,2))
        ene_vs_strain_wo = np.zeros((points,2))
        for idx,d in enumerate(np.linspace(-10.,10.,points)):
            sd = 0.2
            sd = d
            s = sd/100.

            cryst = my_get_cart_deformed_cell(atoms_h, size=sd,vol=False)
            cell = cryst.get_cell()
            scheck =strain= cell[0,1]/cell[0,0]*2.
            print('straincheck',s,scheck)
            ene_vs_strain[idx,0] = s
            ene_vs_strain[idx,0] = scheck
            ene_vs_strain_wo[idx,0] = scheck
            stress = cryst.get_stress()
            enecryst_eV_cell = cryst.get_potential_energy()  # eV for whole cell
            ene_vs_strain[idx,1] = (enecryst_eV_cell-e0)*1000./cryst.get_number_of_atoms()
            ene_vs_strain_wo[idx,1] = (enecryst_eV_cell)*1000./cryst.get_number_of_atoms()
            vol = cryst.get_volume()

            if True:
                if False:
                    print()
                    print('cryst.get_cell()     :')
                    print(cryst.get_cell())
                    print('cryst.get_stress()   :')
                    print(cryst.get_stress())
                if False:
                    print('cryst.energy:',enecryst_eV_cell)
                    print('cryst.get_volume()')
            C44 = (enecryst_eV_cell - e0)/vol*(2./(strain**2.))
            ase_write("out_c_check_vol_cons_widerange_DFTV0.runner",cryst,format='runner',append=True)
            #print(d,'de',de2)
            atb = ANGSTROM_TO_BOHRRADIUS = 1./aseunits.Bohr
            print("volume",str(vol).ljust(20),'s',str(round(s,5)).ljust(10),'ene',enecryst_eV_cell,'c44:',str(round(C44/aseunits.GPa,2)).ljust(5),cell[0]*atb,cell[2]*atb) #stress)
            print("volume",str(vol).ljust(20),'s',str(round(s,5)).ljust(10),'ene',enecryst_eV_cell,'c44:',str(round(C44/aseunits.GPa,2)).ljust(5),cell[0],cell[2]) #stress)
        np.savetxt("elastic_ene.dat",np.array([C44]))
        np.savetxt("ene_vs_strain_NN.dat",ene_vs_strain)
        np.savetxt("ene_vs_strain_NN_wo.dat",ene_vs_strain_wo)
        sys.exit()
        print()
        cryst = my_get_cart_deformed_cell(atoms_h, size=0.2,vol=False)
        print(cryst.get_cell())
        print('stress                ',stress)
        print('stress/2              ',stress/2.)
        print('stress/aeunits.GPa2/2.',stress/aseunits.GPa/2.)
        print('stress/aeunits.GPa2   ',stress/aseunits.GPa)
        print('stress',1000*stress/aseunits.GPa/2.)
        print('st C44',1000*stress[3]/aseunits.GPa/2.)
        print('st C44',1000*stress[5]/aseunits.GPa/2.)
        print()
        print()
        print()
        print("########### get murn structures ###########")
        sys.exit()
        print('linsp',np.linspace(-0.03,0.03,9))
        for d in np.linspace(-0.03,0.03,9):
            cryst.set_cell(atoms_h.get_cell()*(1.+d), scale_atoms=True)
            #print('ah--',atoms_h.get_cell()/4.)
            #print('volh',atoms_h.get_volume()/4.)
            print('volc',cryst.get_volume()/4.,cryst.get_potential_energy())
            #ase_write("out_murn.runner",cryst,format='runner',append=True)

        V0 = 16.476413798491 #  first converged
        V0 = 16.476032301744 #  beset converge
        a0 = (4.*V0)**(1./3.)
        print("V0",V0)
        print("a0",a0)
        cryst.set_cell(np.diag(np.ones(3))*a0, scale_atoms=True)
        print(np.diag(np.ones(3)))
        print(np.diag(np.ones(3))*a0)
        cell = cryst.get_cell()
        print('-->cryst.cell:',cryst.get_cell())
        print()
        #for d in [0,0.01,-0.01]: #np.linspace(-0.01,0.01,3):
        for d in np.linspace(-0.01,0.01,3):
            print()
            cryst.set_cell(cell*(1.+d), scale_atoms=True)
            print('cell0',cryst.get_cell())
            print('cell1',cell)
            print('-->d',d,'    volc',cryst.get_volume()/4.,cryst.get_potential_energy())
            #ase_write("out_c_check.runner",cryst,format='runner',append=True)
        print()
        print("##########")
        cryst.set_cell(cell, scale_atoms=True)
        cryst = my_get_cart_deformed_cell(cryst, size=0.2)
        print('cell0 ',cryst.get_cell())
        print('cell0v',cryst.get_volume()/4.)
        ase_write("out_c_check_vol_cons.runner",cryst,format='runner',append=True)
        return


    def get_calculator(self,atoms):
        if self.verbose > 1:
            for i in self.lmpcmd:
                print("lmpcmds    :",i)
            print("atom_types :",self.atom_types)
            print("keep_alive :",self.keep_alive)
        if self.calculator == "lammps":
            asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=self.keep_alive)
            atoms.set_calculator(asecalcLAMMPS)
        return

    def ase_relax_atomic_positions_only(self,atoms,fmax=0.0001,verbose=False):
        ''' The strain filter is for optimizing the unit cell while keeping scaled positions fixed. '''
        self.keep_alive = True
        self.get_calculator(atoms)

        if verbose:
            print('relax atomic positions; stress:',atoms.get_stress(),"volume per atom:",ase_vpa(atoms))

        logfile="-" # output to screen
        logfile="tmp" # output to file and not to screen
        opt = LBFGS(atoms,logfile=logfile)
        opt.run(fmax=fmax)
        if os.path.isfile("tmp"):
            os.remove("tmp")
        if verbose:
            print('relax atomic positions; stress:',atoms.get_stress(),"volume per atom:",ase_vpa(atoms))
        return

    def ase_relax_cellshape_and_volume_only(self,atoms,verbose=False):
        ''' The strain filter is for optimizing the unit cell while keeping scaled positions fixed. '''
        self.keep_alive = True
        #print('1')
        self.get_calculator(atoms)
        #print('2')

        if verbose: self.check_frame('ase_relax_cellshape_and_volume_only in',frame=atoms,verbose=verbose)
        #print('3')

        sf = StrainFilter(atoms)
        logfile="-" # output to screen
        logfile="tmp" # output to file and not to screen
        opt = BFGS(sf,logfile=logfile)
        opt.run(0.005)
        if os.path.isfile("tmp"):
            os.remove("tmp")
        if verbose: self.check_frame('ase_relax_cellshape_and_volume_only out',frame=atoms)
        return

    def get_murn(self,atomsin=False,verbose=False,
            return_minimum_volume_frame=False,
            return_frame_with_volume_per_atom=False,
            atomrelax=False,
            write_energies=False,
            get_to_minvol_first=True):
        ''' the murn will never change the atomsobject
        return_frame_with_volume_per_atom : volume can be specified and he frame scaled
        '''
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        if return_minimum_volume_frame == True or type(return_frame_with_volume_per_atom) != bool:
            atoms_murn = atomsin
        else:
            atoms_murn = atomsin.copy()

        atoms_murn.wrap()

        # probably not necessary since this is set up when necessary
        self.keep_alive = True
        self.get_calculator(atoms_murn)

        ### relax the atoms_murn first to the equilibrium
        if verbose: self.check_frame('get_murn 1 in',frame=atoms_murn)
        self.ase_relax_cellshape_and_volume_only(atoms_murn,verbose=verbose)
        if verbose: self.check_frame('get_murn 2 atfer cellshape relax',frame=atoms_murn)

        if atomrelax:
            self.ase_relax_atomic_positions_only(atoms_murn,fmax=0.0001,verbose=False)
            if verbose: self.check_frame('get_murn 2 atfer atomrelax only',frame=atoms_murn)

        dvol_rel=[0.97,0.975,0.98,0.985,0.99,0.995,0.998,1.0,1.002,1.005,1.01,1.015,1.02,1.025,1.03]
        #dvol_rel = np.arange(0.97,1.03,0.001)
        #dvol_rel = np.arange(0.995,1.005,0.0003)
        vol_pa = np.zeros(len(dvol_rel))
        ene_pa = np.zeros(len(dvol_rel))


        cell_ref = atoms_murn.get_cell()
        nat = atoms_murn.get_number_of_atoms()

        atoms_murn_loop = atoms_murn.copy()

        for idx,i in enumerate(dvol_rel):
            if verbose > 2:
                print()
                print('000 idx:',idx,'i:',i)
            atoms_murn_loop.set_cell(cell_ref*i,scale_atoms=True)
            if verbose > 2:
                print('111 cell',atoms_murn_loop.get_cell())
            #print('111 cell',atoms_murn_loop.get_cell())
            #print('111 cell',atoms_murn_loop.get_cell()[0,0]/atoms_murn_loop.get_cell()[1,1])
            vol=atoms_murn_loop.get_volume()
            if verbose > 2:
                print('222 vol',atoms_murn_loop.get_volume())
            #ene = self.ene(atoms_murn_loop)                         # works
            #ene = self.ene_new(atoms_murn_loop)                         # works
            ene = self.ene_allfix(atoms_murn_loop)                       # works
            #print('ams3',atoms_murn_loop.get_stress(),ase_vpa(atoms_murn_loop))
            if verbose > 2:
                stress = atoms_murn_loop.get_stress()[:3]
                cell = atoms_murn_loop.get_cell()
                pos = atoms_murn_loop.get_positions()[1]/cell[0,0]

                if cell[0,1] == cell[0,2] == cell[1,0] == cell[1,2] == cell[2,0] == cell[2,1] == 0:
                    if cell[0,0] == cell[1,1] == cell[2,2]:
                        if round(stress[0],6) == round(stress[1],6) == round(stress[2],6):
                            print('murn cell++  ',round(cell[0,0],6),pos)
                        else:
                            print('murn cell--  ',cell[0,0],stress)
                    else:
                        print('murn cell---  ',cell[0,0],cell[1,1],cell[2,2])
            #ene = ase_enepot(atoms_murn_loop) #,units=ace.units)    # core dump
            #ene = atoms_murn_loop.get_potential_energy()            # core dump
            #ene = ase_enepot(atoms_murn_loop,units=self.units,verbose=self.verbose)  # core dump
            if verbose > 2:
                print('333 ene',ene) #,ene2)
            vol_pa[idx] = vol/nat
            ene_pa[idx] = ene/nat
            if verbose > 2:
                print('idx:',str(idx).ljust(3),'i:',str(i).ljust(10),'vol:',str(vol).ljust(10),'ene:',ene)
            if verbose:
                stress = atoms_murn_loop.get_stress()[:3]
                print('i',str(i).ljust(5),'vol/nat',str(round(vol/nat,7)).ljust(10),'ene/nat',str(ene/nat).ljust(19),stress)

        if verbose: self.check_frame('get_murn 3 atfer loop           ',frame=atoms_murn)
        if verbose > 1:
            stress = atoms_murn.get_stress()[:3]
            cell = atoms_murn.get_cell()
            pos = atoms_murn.get_positions()[1]/cell[0,0]
            print("---1-->>",pos)
        if write_energies:
            print('we1')
            if type(write_energies) == bool:
                print('we2')
                write_energies = "energies.dat"
            np.savetxt(write_energies,np.transpose([vol_pa,ene_pa]))
        if verbose > 1:
            print('loop done')
        vinet = eos()
        data=np.transpose([vol_pa,ene_pa])
        vinet.fit_to_energy_vs_volume_data(datax=vol_pa,datay=ene_pa)
        #self.eos = vinet.parameters
        if verbose > 1:
            print('pars',vinet.parameters)
        if verbose: self.check_frame('get_murn 4 before min vol ret   ',frame=atoms_murn)
        #if return_minimum_volume_frame == True or type(return_frame_with_volume_per_atom) != bool:
        #    print('adapt atoms in')
        #    atomsin = atoms_murn.copy()
        #        # old
        #        #volume_in  = ase_vpa(atomsin)
        #        #volume_out = vinet.parameters[1]
        #        #if type(return_frame_with_volume_per_atom) != bool:
        #        #    volume_out = return_frame_with_volume_per_atom
        #        #volume_scale = (volume_out/volume_in)**(1./3.)
        #        ##print('volume_in',volume_in)
        #        ##print('volume_out',volume_out)
        #        ##print('scale',volume_scale)
        #        #atomsin.set_cell(atomsin.get_cell()*volume_scale,scale_atoms=True)

        if verbose: self.check_frame('get_murn 5 after  min vol ret   ',frame=atoms_murn)



        if verbose > 2:
            stress = atoms_murn.get_stress()[:3]
            cell = atoms_murn.get_cell()
            pos = atoms_murn.get_positions()[1]/cell[0,0]
            print("---2-->>",pos)
        return vinet.parameters

    def get_fh(self,atomsin=False,disp=0.03,debug=False,try_readfile=False,atomrelax=True):
        ''' the function will never change the atomsobject
        atomrelax: in most cases it is desirable to relax the atomic positions
        in few cases we might be tempted to assess the free energy for particular positions (e.g. to check weather the DFT equilibrium position is the stable position of a NN)
        '''
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        atoms_h = atomsin.copy()
        atoms_h.wrap()

        T0shift_ev_atom = self.ene(atoms_h)/atoms_h.get_number_of_atoms()
        #print('--> T0shift_ev_atom',T0shift_ev_atom)
        #sys.exit('T0shift_ev_atom')

        #if return_units == "mev_pa":
        #    return_mult = 1.
        #elif return_units == "ev_cell":
        #    return_mult = atomsin.get_number_of_atoms()/1000.
        #else:
        #    raise NameError('return_units conversion not yet set for '+return_units)

        if try_readfile:
            if os.path.isfile(try_readfile+"_hessematrix"):
                if debug: print('tryread 1x .....',try_readfile+"_hessematrix")
                hessematrix = hesse.read_Hessematrix(try_readfile+"_hessematrix")
                if debug: print('hesse 2x ...',hessematrix)
                if debug: print('hesse 3x ...',atoms_h.get_chemical_symbols())
                hes = hesse.hesseclass(listin=atoms_h.get_chemical_symbols(),H=hessematrix,show_negative_eigenvalues = False, Tmax=1000, T0shift_ev_atom = T0shift_ev_atom)
                if hes.has_negative_eigenvalues == True:
                    print("Negative Eigenvalues",hes.freqs)

                if debug: print('done 1x ...')

                #free_ene      = hes.ene_atom
                #print('fe pa',free_ene[:3])
                #print('fe pa',free_ene[-3:])
                #free_ene_cell = hes.ene_cell
                #print()
                #print('fe cell',free_ene_cell[:3])
                #print('fe cell',free_ene_cell[-3:])
                #print()
                #print('fe cell ev',hes.ene_cell_only_ev)
                #print()
                #print('ka',hes.ene_cell_only_ev_T0shifted[:3])
                #print('ka',hes.ene_cell_only_ev_T0shifted[-3:])
                #sys.exit()
                #get = np.loadtxt(try_readfile)
                #return np.transpose([get[:,0],get[:,1]*return_mult])
                return hes


        self.keep_alive = False
        self.get_calculator(atoms_h)

        if debug:
            print("###########################################")
            #`print("forces harmonic 3:",atoms_h.get_forces()[:3])
            #print("###########################################")
            print('in atoms_h str',atoms_h.get_stress())
            maxforce = np.abs(atoms_h.get_forces()).max()
            print('in maxforce in',maxforce)

        #ene = self.ene(atoms_h,cellrelax=True,atomrelax=True,print_minimization_to_screen=debug)
        if atomrelax == True:
            self.ase_relax_atomic_positions_only(atoms_h,fmax=0.0001,verbose=False)
            print('harmonic cell stress',atoms_h.get_stress())
        maxforce = np.abs(atoms_h.get_forces()).max()
        if debug:
            print("###########################################")
            #print("forces harmonic 4:",atoms_h.get_forces()[:3])
            #print("###########################################")
            print('out atoms_h str',atoms_h.get_stress())
            print('out maxforce out',maxforce)

        if atomrelax == True and maxforce > 0.0001:
            print("forces harmonic 4:",atoms_h.get_forces()[:3])
            print('maxforce',maxforce)
            sys.exit('maxforce is too large')

        if maxforce > 0.0001:
            # from this it needs to be deduced that NEGATIVE EIGENVALUES
            # calculating the hesse makes a segmentation fault with ace lammps
            return False
        if debug:
            print("###########################################")
            print('atoms_h.get_cell() before',atoms_h.get_cell())

        nat = atoms_h.get_number_of_atoms()
        if nat < 20:
            atoms_h *= (2,2,2)
        if debug:
            nat = atoms_h.get_number_of_atoms()
            print("###########################################")
            print('!!!!!!! nat:',nat)
            print("forces harmonic 2:",atoms_h.get_forces()[:3])
            print("stress:",atoms_h.get_stress())
            print("forces max harmonic 2:",abs(atoms_h.get_forces()).max())
            print('atoms_h.get_cell() after mult',atoms_h.get_cell())
            print("###########################################")
        pos0 = atoms_h.get_positions()
        hessematrix=np.zeros((pos0.shape[0]*3,pos0.shape[0]*3))

        ### schleife ueber alle atome, 1..32

        for iidx,i in enumerate(pos0): # loop over all atoms
            progress(iidx,len(pos0),status=try_readfile)
            #print(iidx,"/",pos0.shape[0]) # loop over xyz 1..3
            for jidx,j in enumerate(i):
                pos1 = np.copy(pos0)
                pos1[iidx,jidx] = pos1[iidx,jidx]+disp

                atoms_h.set_positions(pos1)
                fah = atoms_h.get_forces()
                hessematrix[iidx*3+jidx] = fah.reshape((1,pos0.shape[0]*        3))/(-disp)

        #np.savetxt("HesseMatrix.dat",hessematrix/97.173617,fmt="%.13f")
        if debug:
            print("get free energy ...")
            print('get_chemical_symbols()',atoms_h.get_chemical_symbols())
            print()

        hes = hesse.hesseclass(listin=atoms_h.get_chemical_symbols(),H=hessematrix,show_negative_eigenvalues = True, Tmax=1000, T0shift_ev_atom = T0shift_ev_atom)
        print('eigenvalues',hes.freqs)
        #try:
        #    #free_ene = (hes.ene_atom[300]-hes.ene_atom[0])[1]
        #    free_ene      = hes.ene_atom
        #    free_ene_cell = hes.ene_cell
        #except IndexError:
        #    free_ene = "UNSTABLE"
        #if debug and type(free_ene) != str:
        #    hes.write_ene_atom()
        if try_readfile:
            if True: #hes.has_negative_eigenvalues == False:
                # write in any case, if it has negative eigenvalues so be it
                hes.write_hessematrix(try_readfile+"_hessematrix")
                if hes.has_negative_eigenvalues == False:
                    # only write for ground state
                    hes.write_ene_atom(try_readfile+"_per_atom")
                    hes.write_ene_cell(try_readfile+"_per_cell")
            #np.savetxt(try_readfile,free_ene)

        #print('k T)shift   ',T0shift_ev_atom)
        #print('hes.ene_cell',hes.ene_cell)
        #print('hes.ene_atom',hes.ene_atom)
        #get = np.loadtxt(try_readfile)
        #return np.transpose([get[:,0],get[:,1]*return_mult])
        #return free_ene*return_mult
        return hes


    def check_frame(self,text,frame=False,verbose=True,setupcalc=True):
        #print('a')
        if type(text) != str:
            raise TypeError("need a text!")
        #print('b')

        if setupcalc == True:
            #print('c',verbose)
            self.keep_alive = True
            self.get_calculator(frame)
        #print('d',verbose)

        if verbose:
            print('check_frame::')
            print('check_frame::',frame.get_stress())
        check_stress_max = round(abs(frame.get_stress()).max(),5)
        check_vpa = round(ase_vpa(frame),2)  # only 2 digits can be nicely fitted
        check_force_max = round(abs(frame.get_forces()).max(),5)
        if verbose:
            print(text.ljust(42),'fm,sm,curr vol ',[check_force_max,check_stress_max,check_vpa]) # forces max , stress max
        return ['fm,sm,curr vol ',check_force_max,check_stress_max,check_vpa]



    def submit_aiida(self,atomsin=False):
        ## o) move this from the ase calass to something separate
        ## a) create inputxxx.data
        ## b) follow kmc_submit_inputdata_to_aiida.sh
        return

def ase_vpa(atoms):
    return atoms.get_volume()/atoms.get_number_of_atoms()
def ase_epa(atoms):
    return atoms.get_potential_energy()/atoms.get_number_of_atoms()
def ase_mepa(atoms):
    return atoms.get_potential_energy()/atoms.get_number_of_atoms()*1000.
def ase_fmax(atoms):
    return abs(atoms.get_forces()).max()

def ase_repeat_structure(atoms,repeat):
    atomsc = atoms.repeat(repeat)
    cell_o = atoms.get_cell()
    cell_n = cell_o * repeat


    pos_o = atoms.get_positions()
    for_o = atoms.get_forces()
    nat_o = atoms.get_number_of_atoms()
    #print()
    #print(pos_o)
    #print()
    #print('sh old',pos_o.shape)
    #print('nat_o',nat_o)
    nat_n = nat_o*repeat**3
    #print('nat_n',nat_n)
    pos_n = np.zeros((nat_n,3))
    for_n = np.zeros((nat_n,3))
    #print('sh new',pos_n.shape)
    ## create translation matrix
    N1 = N2 = N3 = repeat
    var1 = np.mgrid[ 0:N1, 0:N2, 0:N3 ]
    ntrans1 = var1[0,:,:].flatten(1)
    ntrans2 = var1[1,:,:].flatten(1)
    ntrans3 = var1[2,:,:].flatten(1)
    transmat = np.zeros( [ len(ntrans1),3] )
    transmat[:,0] = ntrans1
    transmat[:,1] = ntrans2
    transmat[:,2] = ntrans3
    #print('tr',transmat)
    sclattpos = np.dot( transmat, cell_o )

    for ij, sclatt in enumerate(sclattpos):
        #print('-> ij',ij,sclatt)
        for idx,i in enumerate(pos_o):
            #print('ij',ij,'i',i,type(ij),type(i))
            #print('-> pos_o',i,'--->',sclatt+i)
            pos_n[ij*nat_o+idx] = sclatt + pos_o[idx]
            for_n[ij*nat_o+idx] = for_o[idx]
            #print('i',i,pos_n[ij+idx])

    #print('pos c ??')
    atomsc.set_positions(pos_n)
    #print(atomsc.get_positions()[:5])
    #sys.exit()
    return atomsc,for_n

class ase_get_known_formats_class():
    """
    helptext
    """
    def __init__(self,verbose = False):
        self.formatspy              = os.path.dirname(ase.io.__file__)+"/formats.py"
        self.all_known_formats      = []
        self.all_known_formats_ase  = False
        self.my_formats_shall       = [ 'runner',   'lammps-runner', 'lammps-data'   ,'ipi'   , 'quippy'    ]
        self.my_formats_filenames   = [ "runner.py","lammpsrunner.py","lammpsdata.py","ipi.py", "quippy.py" ]
        self.my_formats_is          = []
        self.verbose                = verbose
        self.needcopy_since_missing = False
        return

    def get_all_known_formats(self):
        if len(self.all_known_formats) == 0:
            self.all_known_formats_ase = ase.io.formats.all_formats
            for i in self.all_known_formats_ase:
                self.all_known_formats.append(i)
        return

    def check_if_format_in_know_formats(self,typ):
        if typ in self.all_known_formats:
            if self.verbose: print(">> formats.py knows", typ)
            return True
        else:
            self.verbose = True
            if self.verbose: print(">> ERROR, formats.py does not know",typ)
            return False

    def copy(self):
        scripts = os.environ['scripts']
        from_ = scripts+"/runner_scripts/ase_fileformat_for_"
        to = os.path.dirname(ase.io.__file__)+"/"
        for ff in self.my_formats_filenames:
            #print('copying ',from_+ff,'to',to+ff)
            if self.verbose:
                print('copying ',from_+ff)
                print('                    to',to+ff)
            shutil.copyfile(from_+ff,to+ff)

    def adapt_formatspy(self,writeformatspy = False):
        # check if formatspy exist
        if not os.path.isfile(self.formatspy):
            print('formatspy',self.formatspy)
            sys.exit('did not find '+str(self.formatspy))

        f = open(self.formatspy, "r")
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

        for fo in [ 'runner', 'ipi', 'quippy' ]:
            if fo in self.all_known_formats_ase:
                if self.verbose:
                    print(fo.ljust(14)+'format are already added in formats.py (of ase).')
            else:
                print(fo.ljust(14)+'format NOT KNOWN in formats.py (of ase, WILL BE ADDED).')
                contents.insert(insert, "    '"+fo+"': ('"+fo+" input file', '+F'),\n")
                writeformatspy = True


        if 'lammps-runner' in self.all_known_formats_ase:
            if self.verbose:
                print('lammps-runner format are already added in formats.py (of ase).')
        else:
            contents.insert(insert, "    'lammps-runner': ('LAMMPS data input file for n2p2 or runner', '1F'),\n")
            contents.insert(insert2,"    'lammps-runner': 'lammpsrunner',\n")
            writeformatspy = True

        if writeformatspy == True:
            print('! adapting formatspy ... ! (seems to be save in any case)')
            #print('insert',insert)

            f = open(self.formatspy, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()
            sys.exit('! since adapting formatspy need to exit here!')
        else:
            if self.verbose:
                print('everything was already in',self.formatspy)
        return

    def check_if_default_formats_known(self,copy_and_adapt_formatspy_anyhow=False):
        self.get_all_known_formats()
        if self.verbose: print(">> formats.py", self.formatspy)
        for i in self.my_formats_shall:
            if self.check_if_format_in_know_formats(i) == False: self.needcopy_since_missing = True

        if self.needcopy_since_missing == True or copy_and_adapt_formatspy_anyhow == True: # only in this case copy stuff and add to formats.py
            self.verbose = True
            self.copy()
            self.adapt_formatspy(writeformatspy = copy_and_adapt_formatspy_anyhow)
        return

def convert_energy(ene,units_in_,units_out_,frame,verbose=False):
    known = [ "hartree_pa", "ev_pa", "mev_pa", 'ev', "hartree" ]
    units_in  = units_in_.lower()
    units_out = units_out_.lower()
    if verbose:
        print('units_in :',units_in)
        print('units_out:',units_out)
    if units_in not in known or units_out not in known:
        print('units_in :',units_in)
        print('units_out:',units_out)
        print('known    :',known)
        sys.exit("units_in or units_out not in known")

    units_in_split  = units_in.split("_")
    units_out_split = units_out.split("_")
    #print('us in ',units_in_split,units_in_split[0])
    #print('us out',units_out_split,units_out_split[0])
    #print('ene in',ene)

    if units_in_split[0] == units_out_split[0]:
        pass
    elif units_in_split[0] == 'ev' and units_out_split[0] == 'hartree':
        if verbose:
            print('111',aseunits.Hartree)
        ene = ene / aseunits.Hartree
    elif units_in_split[0] == 'hartree' and units_out_split[0] == 'ev':
        ene = ene * aseunits.Hartree
        if verbose:
            print('2222',aseunits.Hartree)
    else: sys.exit('conversion from to now yet implemented')
    if verbose:
        print('ene (1)',ene)

    if len(units_in_split) == 2 and len(units_out_split) == 1:
        ene = ene * frame.get_number_of_atoms()
    elif len(units_in_split) == 1 and len(units_out_split) == 2:
        ene = ene / frame.get_number_of_atoms()
    if verbose:
        print('ene (2)',ene)
    return ene

def show_ase_atoms_content(atoms,showfirst=10,comment = ""):
    print()
    print("#####################################################")
    print("# "+comment+" BEGIN #######")
    print("#####################################################")
    print('## atoms')
    print(atoms)
    print('## atoms.get_number_of_atoms()',atoms.get_number_of_atoms())
    #print(atoms.positions)
    print(atoms.get_positions()[:showfirst])
    print('## elements get_chemical_symbols()')
    print(atoms.get_chemical_symbols()) #[:showfirst])
    print(list(set(atoms.get_chemical_symbols())))
    print('## atoms.cell')
    print(atoms.cell)
    print('## aa.get_cell_lengths_and_angles()')
    print(atoms.get_cell_lengths_and_angles())
    print('##atom.numbers',atoms.numbers)
    print("#####################################################")
    print("# "+comment+" END #######")
    print("#####################################################")
    print()
    return


if __name__ == '__main__':
    pass
