#!/usr/bin/env python

# general libraries
import os
import sys
import math
import numpy as np
import copy
import glob
#import argcomplete, argparse
from argparse import ArgumentDefaultsHelpFormatter
import argparse
import shutil
import imp
import numpy.core.arrayprint as arrayprint
import contextlib
import datetime
#from numba import jit
from scipy.ndimage        import map_coordinates


# non general libraries
import crystal_generator
reload(crystal_generator)
import utils
reload(utils)
import pot_parametrize
reload(pot_parametrize)
import hesse
#print list(sys.modules.keys())
#print dir(utils)

import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 25)
pd.set_option('display.max_colwidth', 9)

# TODO: (Ueberlegungen und Ideen):
# - ein atom0 hat einen ersten nachbarn1 und einen zweiten nachbarn2; das ganze sytem schwingt;
#   beide atome haben unterschiedliche funktionen der wechselwirkung mit atom 0;
#   in abhaengigkeit von dem Abstand von nachbar1 und nachbar2 von atom0 ist die frage
#   ob die zweite nachbar wechselwirkung gaendert werden sollte zur ersten nachbar wechselwirkung wenn
#   einmal der nachbar2 naeher an atom0 ist als nachbar1; auch ist hier die frage inwieweit
#   der tausch der funkionen einen platzwechsel beguenstigt (oder nicht).
# - aus der MD sollte man mal eine lineare kette von atomen rausnehmen und dann deren
#   schwingungen untersuchen. Alle atome welche nicht in der kette sind einfach mal rausnehmen.
#   Es waere interessant zu sehen ob man nur aus den positionen der atome direkt die langwelligen
#   schwingungen sehen kann.
# - what is also interesting is to ask weather there might be some kind of coupling between the
#   to / ti modes with the longitudinal mode. Could this be estimated from T=0K or would MD be necessary?
# - is it possible to construct (parts of) phonon dispersion curves directly form the measured
#   forces?


@contextlib.contextmanager
def printoptions(strip_zeros=True, **kwargs):
    origcall = arrayprint.FloatFormat.__call__
    def __call__(self, x, strip_zeros=strip_zeros):
        return origcall.__call__(self, x, strip_zeros)
    arrayprint.FloatFormat.__call__ = __call__
    original = np.get_printoptions()
    np.set_printoptions(**kwargs)
    yield
    np.set_printoptions(**original)
    arrayprint.FloatFormat.__call__ = origcall

#@jit
def inversepot(r,alpha,B,e):
    '''
    returnes energy in [eV]
    alpha (no units)
    B (Angstrom)
    e (eV)
    in alfes paper:
        alpha = 6.7
        B = 1.85
        e = 1.
    '''
    return 4.*e*((B/r)**alpha)

#@jit
def inversepot_derivative(r,alpha,B,e):
    ''' derivative of inersepot
    returnes [eV/angstrom]
    in alfes paper:
        alpha = 6.7
        B = 1.85
        e = 1.
    '''
    return (-4.*alpha*B*e*((B/r)**(alpha-1)))/r**2

def LJ(r,eps,rm):
    ''' Lennard Jones
    returnes energy in [eV]
    r: distance of atoms [angstrom]
    rm: nearest neighbor distance [angstrom]
    eps: depth of potential well in [eV]

    '''
    return eps*( (rm/r)**12. - 2.*(rm/r)**6. )

def LJ_derivative(r,eps,rm):
    ''' derivative of Lennard Jones
    returnes forces in [eV/angstrom]
    r: distance of atoms [angstrom]
    rm: nearest neighbor distance [angstrom]
    eps: depth of potential well in [eV]

    '''
    return eps*( -12.*rm**12./r**13. + 12.*rm**6./r**7. )
#@jit
def Morse(r,De,aa,re):
    ''' Morse potential
    returnes energy in [eV]
    r: distance of atoms [angstrom]
    re: nearest neighbor distance [angstrom]
    '''
    return De*(1.-np.exp(-aa*(r-re)))**2
#@jit
def Morse_derivative(r,De,aa,re):
    ''' derivative of Morse potential
    returnes [eV/angstrom]
    '''
    return 2.*aa*De*np.exp(-aa*(r-re))*(1.-np.exp(-aa*(r-re)))

#@jit
def mc1(r,De,aa,re,B,A):
    ''' Energy: x(mathematica) -> r-re(python) '''
    #return Morse(r,De,aa,re)+A*np.exp((-3.+B)*(r-re))+\
    #        (-(6.*(-7.+B)/(-3.+B)**5.)+(6.*(-7.+B)*(r-re)/(-3.+B)**4.)+\
    #        (3.*(-7.+B)*(r-re)**2./(-3.+B)**3.)+((-7.+B)*(r-re)**3./(-3.+B)**2.)+(r-re)**4./(-3.+B))
    return Morse(r,De,aa,re)+A*((r-re)**3.)*np.exp(-(3.-B)*(r-re))
    #+\(-(6.*(-7.+B)/(-3.+B)**5.)+(6.*(-7.+B)*(r-re)/(-3.+B)**4.)+\
    #(3.*(-7.+B)*(r-re)**2./(-3.+B)**3.)+((-7.+B)*(r-re)**3./(-3.+B)**2.)+(r-re)**4./(-3.+B))

#@jit
def mc1_derivative(r,De,aa,re,B,A):
    ''' force '''
    #return Morse_derivative(r,De,aa,re) + A*(r-re)**3.*(1.+(r-re))*np.exp(-(3.-B)*(r-re))
    return Morse_derivative(r,De,aa,re) + \
            A*np.exp((-3.+B)*(r-re))* (r-re)**2 * (3.+(r-re)*(1.+B+(-3.+B)*(r-re)))

def disp_to_longvec(vec0,x,y=False,z=False):
    # in the first case x is the disp vector
    # in the second case x is just the x coordinate of the disp
    if y == False and z == False:
        return vec0-x
    else:
        return vec0-np.array([x,y,z])

def longvec_to_disp(vec0,longvec):
    return vec0-longvec

class parameters_init_class( object ):
    '''Class which savec parameters'''
    def __init__( self ):
        self.u1nn_pottype       = False
        self.u1nn_potparam      = False
        self.u1nn_potadd        = False
        self.u1nn_potaddparam   = False

        self.verbose            = False
        self.save_vecs_to_file_for_DOS = False


class pot_energy_forces_class( object ):
    '''Class to calculate energies anf forces from simple pairwise interaction'''
    def __init__(self):
        ''' necessary:  - cell or cellfile
                        - coordfile_cart or coord_cart or coofile_rrel or coord_rrel
                        - coordfile0_cart or coord0_cart or coofile0_rrel or coord0_rrel
        '''
        self.verbose            = False

        self.coordfile_cart = False
        self.coord_cart = False
        self.coordfile_rrel = False
        self.coord_rrel = False

        self.coordfile0_cart = False
        self.coord0_cart = False
        self.coordfile0_rrel = False
        self.coord0_rrel = False

        self.cellfile = False
        self.cell = False

        self.crystal0neverchange = False  #crystal instance of the undisplaces structure
        self.crystalneverchange = False   #crystal instance of the current positions
        self.numberofatoms = False
        self.loop_over_all_atoms_with_1nn_pot = False
        self.mix_uref_uharmonic = False
        #self.params = False                 # has to be loaded in from other script

        self.dudlpkl = 'dUdL.pkl'    # filename for pkl; default dUdL.pkl
        self.measure_time = False
        self.steps_up_to = False
        self.write_vecs_to_disc = False
        self.write_everything_to_pkl = False

        self.params = parameters_init_class()

        self.load_parameters_from_file = True

        self.energy_DFT = False
        self.ss = False

        ################### vecs for plotting ##############################
        self.longvec_all_mapped_to_first_quadrant = False


    ######################################################################################
    # DETERMINE POSITIONS AT THE BEGINNING
    ######################################################################################
    def get_or_load_parameters(self):
        ''' loads parameters form calculate_energy_and_forces file '''
        if self.load_parameters_from_file == False:
            # in this case you have to define before:
            #self.params = parameters_init_class()
            # self.params.u1nn_pottype = ....
            pass
        else:
            # here we load all parameters from calculate_energy_and_forces file
            calculate_energy_and_forces = imp.load_source('calculate_energy_and_forces', os.getcwd()+'/calculate_energy_and_forces')
            self.params = calculate_energy_and_forces.parameters()

            try:
                self.mix_uref_uharmonic = self.params.mix_uref_uharmonic
            except AttributeError:
                self.mix_uref_uharmonic = False
        self.datapkl = "2x2x2sc_data_3x3x3kp.pkl"
        if os.path.isfile(self.datapkl):
            self.pot_parametrize = pot_parametrize.forcesneighbors()
            self.pot_parametrize.loadDatapd(self.datapkl)

        self.shells = self.determine_shells()
        return

    def get_or_load_positions(self):
        ''' loads in positions in either cartesian or reduced coordinates
        cartesian coordinaes are read in from the file: cartesian_coords (if no filename specified or array given directly)
        relative coordinates are read in from the file: EqCoords_direct (if no filename specified or array given directly)
        you can also define self.coord_{cart,rrel} and self.cell and self.coord0_{cart,rrel}
        '''
        if type(self.coordfile_cart) == bool and type(self.coord_cart) == bool \
            and type(self.coordfile_rrel) == bool and type(self.coord_rrel) == bool:
                if os.path.isfile("cartesian_coords") == True:
                    self.coordfile_cart = "cartesian_coords"
        if type(self.coordfile_cart) == bool and type(self.coord_cart) == bool \
            and type(self.coordfile_rrel) == bool and type(self.coord_rrel) == bool:
                sys.exit("you need the displaced structure \
                        e.gl. cartesian_coords")

        self.get_or_load_cell()
        self.get_or_load_crystal0()

        # We have to ensure we dont read in the EqCoords file every time;
        if self.verbose > 1:
            print ""
            print "0: cartesian coords: ---------------------"
            print "1:","self.coordfile_cart:",  self.coordfile_cart
            print "2:","self.coord_cart    :",  self.coord_cart
            print "3:","self.cellfile      :",  self.cellfile
            print "4:","self.cell          :",  self.cell
            print "5:","self.coordfile_rrel:",  self.coordfile_rrel
            print "6:","self.coord_rrel    :",  self.coord_rrel
        self.crystal = crystal_generator.crystal()
        self.crystal.load_positions_cell(coordfile_cart = self.coordfile_cart, coord_cart = self.coord_cart,
            cellfile = self.cellfile, cell = self.cell,
            coordfile_rrel = self.coordfile_rrel, coord_rrel = self.coord_rrel)


        ##############################################################################
        # here we need: self.crystal.remove_mapping_into_originalcell()
        ##############################################################################
        if self.verbose > 3:
            print "self.crystal.rcar"
            print self.crystal.rcar
            print ""
            print "self.crystal0.rcar"
            print self.crystal0.rcar
        u = self.crystal.rcar-self.crystal0.rcar
        #print "uorig:",u
        for ind,i in enumerate(u):
            if  u[ind][0] > self.crystal.cellvec[0,0]/2.:
                u[ind][0] = u[ind][0] - self.crystal.cellvec[0,0]
            if  u[ind][1] > self.crystal.cellvec[0,0]/2.:
                u[ind][1] = u[ind][1] - self.crystal.cellvec[0,0]
            if  u[ind][2] > self.crystal.cellvec[0,0]/2.:
                u[ind][2] = u[ind][2] - self.crystal.cellvec[0,0]
        if self.verbose > 3:
            print "u"
            print u
        self.crystal.rcar = u + self.crystal0.rcar
        self.crystal.update_rrel_from_rcar()
        self.crystal.update_xyzcar_from_rcar()
        self.crystal.update_xyzrel_from_xyzcar()
        if self.verbose > 3:
            print "self.crystal.rcar"
            print self.crystal.rcar


        if self.verbose > 1:
            print ""
            print "0: EqCoordsDirect: ---------------------"
            print "1:","self.coordfile0_cart:", self.coordfile0_cart
            print "2:","self.coord0_cart    :", self.coord0_cart
            print "3:","self.cellfile       :", self.cellfile
            print "4:","self.cell           :", self.cell
            print "5:","self.coordfile0_rrel:", self.coordfile0_rrel
            print "6:","self.coord0_rrel    :", self.coord0_rrel

        # copy class instance to not have to load in againg every time neighbors were found
        # (when neighbors are found relative positions (and car coords) are changed)
        self.crystal0neverchange = copy.deepcopy(self.crystal0)
        self.crystalneverchange = copy.deepcopy(self.crystal)
        self.numberofatoms = self.crystal.rcar.shape[0]


        # This ist just for information to show!
        self.print_neighbors(copy.deepcopy(self.crystal0))
        return

    def get_or_load_crystal0(self, EqCoords_direct_path = False, cellpath = False):
        ''' initializes crystal0, (EqCoords_direct) crystal structure

            considered are:     - self.coord0_cart
                                - self.coord0_rrel
                                - self.cell

            if coordinates are saved in files:
                                - self.coordfile0_cart
                                - self.coordfile0_rrel
                                - self.cellfile

            the path to this files can also be directly given to this function
                    EqCoords_direct_path
                    cellpath
        '''
        if type(self.coordfile0_cart) == bool and type(self.coord0_cart) == bool \
            and type(self.coordfile0_rrel) == bool and type(self.coord0_rrel) == bool:
                if os.path.isfile("EqCoords_direct") == True:
                    self.coordfile0_rrel = "EqCoords_direct"
        if type(self.coordfile0_cart) == bool and type(self.coord0_cart) == bool \
            and type(self.coordfile0_rrel) == bool and type(self.coord0_rrel) == bool:
                sys.exit("you need the undisplaced structure \
                        (e.gl. EqCoords_direct) as reference to get 1NN atoms")
        self.crystal0 = crystal_generator.crystal()

        self.crystal0.load_positions_cell(coordfile_cart = self.coordfile0_cart, coord_cart = self.coord0_cart,
            cellfile = self.cellfile, cell = self.cell,
            coordfile_rrel = self.coordfile0_rrel, coord_rrel = self.coord0_rrel)
        return

    def get_or_load_cell(self):
        if type(self.cell) == bool and type(self.cellfile) == bool:
            if os.path.isfile("cell") == True:
                self.cellfile = "cell"
        if type(self.cell) == bool and type(self.cellfile) == bool:
            if os.path.isfile("POSCAR") == True:
                utils.run2("rm -f cell; POSCAR_cell_cartesian.sh > cell") # to create cell
            if os.path.isfile("cell") == True:
                self.cellfile = "cell"
        if type(self.cell) == bool and type(self.cellfile) == bool:
            if os.path.isfile("OUTCAR") == True or os.path.isfile("OUTCAR.gz"):
                utils.run2("rm -f cell; OUTCAR_cell-last-cartesian-ARRAY.sh > cell") # to create cell
            if os.path.isfile("cell") == True:
                self.cellfile = "cell"
        if type(self.cell) == bool and type(self.cellfile) == bool:
            sys.exit("you need the cell file")
        return

    def repeat_supercell(self):
        ''' repeats rcar to sc, necessary for 2NN and furter out NN
            before the whole cell can be repeated, the mapping to the original cell
            has to be remove (in other words: the atoms have to be close to their origial
            undisplaced positions, otherwise crystal_generator.center_atoms_around_atom
            does not work properly '''
        self.sc0 = crystal_generator.supercell()
        self.sc = crystal_generator.supercell()
        nsc = 2  # for now on we will just double the supercell, in general maybe factor
                 # 3 or greater necessary
                 # mit 2 koennen wir bei fcc (cubic cell) bis zum 7ten nachbarn gehen
                 #
        self.sc0.create_supercell( self.crystal0, nsc, nsc, nsc, newsorting = True )
        self.sc.create_supercell(  self.crystal, nsc, nsc, nsc , newsorting = True )

        self.print_neighbors(copy.deepcopy(self.crystal0),copy.deepcopy(self.sc0))

        check0 = utils.remove_duplicates_of_2d_array(self.sc0.rcar)
        if check0.shape != self.sc0.rcar.shape:
            print "shape:",check0.shape,self.sc.rcar.shape
            sys.exit("there were duplicates in self.sc0.rcar")
        check = utils.remove_duplicates_of_2d_array(self.sc.rcar)
        if check0.shape != self.sc.rcar.shape:
            print "shape:",check0.shape,self.sc.rcar.shape
            sys.exit("there were duplicates in self.sc.rcar")
        if np.array_equal(self.crystal.rcar,self.sc.rcar[:self.numberofatoms]) != True:
            print self.crystal.rcar
            print "-------------"
            print self.sc.rcar[:self.numberofatoms]
            sys.exit("atoms have swiched 1")
        self.sc0neverchange = copy.deepcopy(self.sc0)
        self.scneverchange = copy.deepcopy(self.sc)
        return

    ######################################################################################
    # PRINT INFORMATION / SAVE INFORMATION TO HARDDRIVE
    ######################################################################################
    def print_neighbors(self, copy_of_crystal0_instance,copy_of_sc0_instance = False):
        ''' shows information of how many atoms are in certain shell '''
        if self.verbose:
            NNabst = copy_of_crystal0_instance.get_NNlist(0, 1,
            cell = copy_of_crystal0_instance.cellvec,
            coord_rrel = copy_of_crystal0_instance.rrel,
            return_NNdist = True)
            NNabstsmall = copy.copy(NNabst)
            show = copy_of_crystal0_instance
            if type(copy_of_sc0_instance) != bool:  # if we do have a supercell
                NNabst = copy_of_sc0_instance.get_NNlist(0, 1,
                cell = copy_of_sc0_instance.cellvec,
                coord_rrel = copy_of_sc0_instance.rrel,
                return_NNdist = True)
                show = copy_of_sc0_instance
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!
            # HERE CHANGE not to get sc numbers!!!!!!!!!!!!!!!!!!!!!!!!!!

                self._info_nnlist_sc = np.zeros(len(NNabstsmall))
                self._info_nnlist_sc[:] = np.nan
                self._info_nnatoms_sc = np.zeros(len(NNabstsmall))
                self._info_nnatoms_sc[:] = np.nan
            if type(copy_of_sc0_instance) == bool:  # if we do not have a supercell
                self._info_nnlist_pc = np.zeros(len(NNabstsmall))
                self._info_nnlist_pc[:] = np.nan
                self._info_nnatoms_pc = np.zeros(len(NNabstsmall))
                self._info_nnatoms_pc[:] = np.nan
            for nnatom,nnabst in enumerate(NNabstsmall):
                if nnatom == 0.0:
                    continue
                #print "nnatom>>:",nnatom
                nnlist = show.get_NNlist(0, nnatom,
                        cell = show.cellvec,
                        coord_rrel = show.rrel)
                        #return_NNdist = True)
                #print "nnlist;",nnlist
                abstand = np.linalg.norm(show.rcar[nnlist[0]])
                abstandstr = utils.number_to_string_of_certain_length(abstand, 5, 10)

                printelements = 15
                if len(nnlist) > printelements:
                    if self.verbose > 1:
                        print nnatom,"\t",len(nnlist),"\t",abstandstr,nnlist[:printelements]," ..."
                else:
                    if self.verbose > 1:
                        print nnatom,"\t",len(nnlist),"\t",abstandstr,nnlist
                if type(copy_of_sc0_instance) != bool:  # if we do have a supercell
                    self._info_nnlist_sc[nnatom-1] = int(nnatom)
                    self._info_nnatoms_sc[nnatom-1] = int(len(nnlist))
                if type(copy_of_sc0_instance) == bool:  # if we do not have a supercell
                    self._info_nnlist_pc[nnatom-1] = int(nnatom)
                    self._info_nnatoms_pc[nnatom-1] = int(len(nnlist))
            if self.verbose > 1:
                print "-------------------"*4

        if self.verbose > 1:
            print ""
            print "-------------------"*4
            print "NNabst:",NNabst,NNabstsmall
            print "len(NNabst):",len(NNabst)
            print "NN     atoms       distance atomindex"
            print "-------------------"*4

            # zeige alle moeglichen NN in dieser zelle && welche schalen sind voll drin in sc?
            # 6 in der 2x2x2fcc superzelle
            # 10 in der 3x3x3fcc superzelle
            #
            # in 3x3x3sc:
            # NNabst: [ 0.      2.9204  4.13    5.0582  5.8407  6.5301  7.1534  7.7265  8.7611
            #   9.6857]
            #   len: 10
            # 1     12  [ 0.     2.065  2.065]  full (gesampled in quer)
            # 2     6   [ 0.    0.    4.13   ]  full (gesampled in xdir)
            # 3     24  [ 4.13   2.065  2.065]  full (gesampled in 3NNdir)
            # 4     12  [ 0.    4.13  4.13   ]  full (gesampled in quer)
            # --------------------------------
            # 5     12  [ 0.     2.065  6.195]  DIFFRENT IN 2x2x2 and 3x3x3?
            # 6     8   [ 4.13  4.13  4.13   ]
            # 7     24  [ 4.13   2.065  6.195]
            # 8     3   [ 0.     6.195  6.195]
            # 9     6   [ 4.13   6.195  6.195]
            #
            # in 2x2x2sc:
            # NNabst: [ 0.      2.9204  4.13    5.0582  5.8407  7.1534]  # 7.1534 is the atom exactly in the middle of the supercell
            # len: 6
            # 1     12  (12) [ 0.     2.065  2.065]  full
            # 2     3   ( 6) [ 0.    0.    4.13   ]  NOTFULL (50%) (gesampled in xdir)
            # 3     12  (24) [ 4.13   2.065  2.065]  NOTFULL (50%) (gesampled in 3NNdir)
            # 4     3   (12) [ 0.    4.13  4.13   ]  NOTFULL (25%) (gesampled in quer)
            # -------------------------------------
            # 5     1   [ 4.13  4.13  4.13        ]  NOTFULL  (gesampled in 111dir)  DIFFRENT IN 2x2x2 and 3x3x3?

            # generell
            # SHOW1
            # ----------------------------------------------------------------------------
            # 1     12     2.81428 [   8   16   24  148  156  650  666  798 3209 3217 3349 3851]
            # 2     6         3.98 [   1    2    4  132  642 3201]
            # 3     24     4.87448 [  12   18   25  140  150  157  654  658  667  782  790  799 3213 3219 3225]  ...
            # 4     12     5.62857 [   3    5    6  133  134  643  646  774 3203 3205 3333 3843]
            # 5     24     6.29293 [  9  10  17  20  26  28 144 149 152 158 648 651 664 670 794]  ...
            # 6     8      6.89356 [   7  135  647  775 3207 3335 3847 3975]
            # 7     48      7.4459 [ 13  14  19  22  27  29 141 142 146 151 153 159 652 655 659]  ...
            # ------ up to here with nsc = 2
            # 8     6         7.96 [  32  128  160  640  800 3200]


            # for 2x2x2sc:
            # 1     12     2.82136 [ 8  9 10 11 16 17 20 21 24 26 28 30]
            # 2     3         3.99 [1 2 4]
            # 3     12     4.88673 [12 13 14 15 18 19 22 23 25 27 29 31]
            # 4     3      5.64271 [3 5 6]
            # 5     1      6.91088 [7]
            #
            # NN     atoms       distance atomindex
            # ----------------------------------------------------------------------------
            # 1     12     2.82136 [  8  16  24  52  60  74  90 126 137 145 181 203]
            # 2     6         3.99 [  1   2   4  36  66 129]
            # 3     24     4.88673 [ 12  18  25  44  54  61  78  82  91 110 118 127 141 147 153]  ...
            # 4     12     5.64271 [  3   5   6  37  38  67  70 102 131 133 165 195]
            # 5     24     6.30874 [  9  10  17  20  26  28  48  53  56  62  72  75  88  94 122]  ...
            # ----------------------------------------------------------------------------
            # array([[ 0.   ,  1.995,  1.995],  8
            #        [ 1.995,  0.   ,  1.995],  16
            #        [ 1.995,  1.995,  0.   ],  24
            #        [-1.995,  0.   ,  1.995],  52
            #        [-1.995,  1.995,  0.   ],  60
            #        [ 0.   , -1.995,  1.995],  74
            #        [ 1.995, -1.995,  0.   ],  90
            #        [-1.995, -1.995,  0.   ],  126
            #        [ 0.   ,  1.995, -1.995],  137
            #        [ 1.995,  0.   , -1.995],  145
            #        [-1.995,  0.   , -1.995],  181
            #        [ 0.   , -1.995, -1.995]]) 203


            # this will only hold now for fcc,bcc where we think of the 0th atom in the 'origin'

        return

    def print_parametrization(self):
        ''' print which parameters were used '''
        if not self.verbose: # self.verbose == False or 0
            return
        @contextlib.contextmanager
        def printoptions(strip_zeros=True, **kwargs):
            origcall = arrayprint.FloatFormat.__call__
            def __call__(self, x, strip_zeros=strip_zeros):
                return origcall.__call__(self, x, strip_zeros)
            arrayprint.FloatFormat.__call__ = __call__
            original = np.get_printoptions()
            np.set_printoptions(**kwargs)
            yield
            np.set_printoptions(**original)
            arrayprint.FloatFormat.__call__ = origcall

        with printoptions(precision=3, suppress=True):
            for shell in self.shells:
                print "lon"+str(shell),utils.string_add_spaces(eval("self.params.u"+str(shell)+"nn_pottype"),4,False),utils.string_to_list(eval("self.params.u"+str(shell)+"nn_potparam"))
                for xyz in [ 'x', 'y', 'z' ]:
                    try:
                        if type(eval("self.params.u1nn_top"+xyz)) != bool:
                            print "top"+xyz+":",utils.string_add_spaces(eval("self.params.u1nn_top"+xyz),4,False)
                    except AttributeError:
                        pass

                #if type(self.params.u1nn_topx) != bool:
                #    print "topx:",utils.string_add_spaces(self.params.u1nn_topx,4,False)
                #if type(self.params.u1nn_topy) != bool:
                #    print "topy:",utils.string_add_spaces(self.params.u1nn_topy,4,False)
                #if type(self.params.u1nn_topz) != bool:
                #    print "topz:",utils.string_add_spaces(self.params.u1nn_topz,4,False)
        #print "---------------------------------------------------------------------------"
        return

    def print_energies(self):
        rou=4
        #print "energy / energy (mev):",round(energy,rou),"\t",round(energy*1000/(numberofatoms-1),rou)
        print 2*"---------------------------------------------------------------------------"
        #print "energylong           :",round(np.sum(energylong)/2.,rou)    ,"\t",round(np.sum(energylong)/2.*1000/(numberofatoms-1)           ,rou)
        #print "energytrantopx       :",round(np.sum(energytrantopx)/2.,rou),"\t",round(np.sum(energytrantopx)/2.*1000/(numberofatoms-1)       ,rou)
        #print "energytrantopy       :",round(np.sum(energytrantopy)/2.,rou),"\t",round(np.sum(energytrantopy)/2.*1000/(numberofatoms-1)       ,rou)
        #print "energytrantopz       :",round(np.sum(energytrantopz)/2.,rou),"\t",round(np.sum(energytrantopz)/2.*1000/(numberofatoms-1)       ,rou)
        #print "hrestev/hrestmev     :",round(hrestev,rou),"\t",round(hrestmev,rou)
        print utils.printred("ENERGY (mev/atom):"+\
                str(round(self.energymev,rou))),\
                \
                utils.printgreen(
                "energylong :"+\
                str(round(self.energymevlong       ,rou))\
                ),\
                utils.printblue(
                "top{x,y,z}:"+\
                str(round(self.energymevtopx     ,rou))+" "+\
                str(round(self.energymevtopy     ,rou))+" "+\
                str(round(self.energymevtopz     ,rou))+" ",\
                ),\
                \
                utils.printyellow(
                "tip{x,y,z}:"+\
                str(round(self.energymevtipx    ,rou))+" "+\
                str(round(self.energymevtipy    ,rou))+" "+\
                str(round(self.energymevtipz    ,rou))+" ",\
                )#,\
                #\
                #utils.printyellow(
                #"tip{x,y,z}:",\
                #round(np.sum(energytrantipx)/2.*1000/(numberofatoms-1)       ,rou),\
                #round(np.sum(energytrantipy)/2.*1000/(numberofatoms-1)       ,rou),\
                #round(np.sum(energytrantipz)/2.*1000/(numberofatoms-1)       ,rou)\
                #)
        if os.path.isfile("ene_free_last") == True:
            print "noa:",self.numberofatoms
            self.energy_DFT = np.loadtxt("ene_free_last")*1000/(self.numberofatoms-1)
            print "DFT:",self.energy_DFT
        print utils.printred("ENERGY (DFT/atom):")
        print 2*"---------------------------------------------------------------------------"
        return

    def print_forces(self):
        ''' shows forces and energyies '''
        if self.verbose < 1:
            return
        neighborlist = self.crystal0.get_NNlist(0, 1,
                coord_cart = False, #crystal0.rcar, # hopefully thos do not \
                cell = self.crystal0.cellvec,
                coord_rrel = self.crystal0neverchange.rrel,
                return_NNdist = False, return_d_NNidx = True)
        #print "neighborlist:",neighborlist # neighborlist: [ 0.  2.  2.  4.  2.  4.  4.  5.  1.  1.  1.  1.  3.  3.  3.  3.  1.  1.  3.  3.  1.  1.  3.  3.  1.  3.  1.  3.  1.  3.  1.  3.]
        np.set_printoptions(threshold=np.nan)  # print the whole array
        np.set_printoptions(linewidth=240)    # print only 6 digist after .
        np.set_printoptions(precision=3)    # print only 6 digist after .
        #np.set_printoptions(precision=1)    # print only 6 digist after .
        if os.path.isfile("forces_OUTCAR"):
            fab=np.loadtxt("forces_OUTCAR")
        else:
            fab = np.zeros((self.numberofatoms,3))

        printforces = np.zeros((self.numberofatoms,23))
        faktorshow = 100.
        #faktorshow = 1

        printforces[:,0]  = faktorshow*(self.forceslong[:,0]                                   )
        printforces[:,1]  = faktorshow*(self.forceslong[:,1]                                   )
        printforces[:,2]  = faktorshow*(self.forceslong[:,2]                                   )

        printforces[:,3]  = faktorshow*(self.forcestrantopx[:,0]+self.forcestrantopy[:,0]+self.forcestrantopz[:,0]                                   )
        printforces[:,4]  = faktorshow*(self.forcestrantopx[:,1]+self.forcestrantopy[:,1]+self.forcestrantopz[:,1]                                   )
        printforces[:,5]  = faktorshow*(self.forcestrantopx[:,2]+self.forcestrantopy[:,2]+self.forcestrantopz[:,2]                                   )

        #printforces[:,6]  = faktorshow*(forcestrantopy[:,0]                                   )
        #printforces[:,7]  = faktorshow*(forcestrantopy[:,1]                                   )
        #printforces[:,8]  = faktorshow*(forcestrantopy[:,2]                                   )

        printforces[:,6]  = faktorshow*(self.forcestrantipx[:,0]+self.forcestrantipy[:,0]+self.forcestrantipz[:,0]                                   )
        printforces[:,7]  = faktorshow*(self.forcestrantipx[:,1]+self.forcestrantipy[:,1]+self.forcestrantipz[:,1]                                   )
        printforces[:,8]  = faktorshow*(self.forcestrantipx[:,2]+self.forcestrantipy[:,2]+self.forcestrantipz[:,2]                                   )

        #printforces[:,6]  = faktorshow*(forcestrantipx[:,0]                                   )
        #printforces[:,7]  = faktorshow*(forcestrantipx[:,1]                                   )
        #printforces[:,8]  = faktorshow*(forcestrantipx[:,2]                                   )

        #printforces[:,9]  = faktorshow*(forcestrantopz[:,0]                                   )
        #printforces[:,10]  = faktorshow*(forcestrantopz[:,1]                                   )
        #printforces[:,11]  = faktorshow*(forcestrantopz[:,2]                                   )

        printforces[:,9]   = faktorshow*(self.forces[:,0]                                   )
        printforces[:,10]  = faktorshow*(self.forces[:,1]                                   )
        printforces[:,11]  = faktorshow*(self.forces[:,2]                                   )

        #printforces[:,9]  = faktorshow*(forces[:,0]                                       )
        #printforces[:,10] = faktorshow*(forces[:,1]                                       )
        #printforces[:,11] = faktorshow*(forces[:,2]                                       )

        #printforces[:,9]  = faktorshow*(forcestrantip[:,0]                                   )
        #printforces[:,10] = faktorshow*(forcestrantip[:,1]                                   )
        #printforces[:,11] = faktorshow*(forcestrantip[:,2]                                   )

        printforces[:,12] = np.arange(self.numberofatoms)
        #VASP forces
        printforces[:,13] = faktorshow*(fab[:,0]                                          )
        printforces[:,14] = faktorshow*(fab[:,1]                                          )
        printforces[:,15] = faktorshow*(fab[:,2]                                          )


        #printforces[:,13] = faktorshow*(hrestf[:,0]                                          )
        #printforces[:,14] = faktorshow*(hrestf[:,1]                                          )
        #printforces[:,15] = faktorshow*(hrestf[:,2]                                          )

        printforces[:,16] = neighborlist
        printforces[:,17] = faktorshow*(self.forces[:,0] - fab[:,0]                            )
        printforces[:,18] = faktorshow*(self.forces[:,1] - fab[:,1]                            )
        printforces[:,19] = faktorshow*(self.forces[:,2] - fab[:,2]                            )

        # positions
        printforces[:,20] = self.crystal0neverchange.rcar[:,0]
        printforces[:,21] = self.crystal0neverchange.rcar[:,1]
        printforces[:,22] = self.crystal0neverchange.rcar[:,2]

        if fab.shape[0] != self.numberofatoms:
            printforces[:] = np.nan
            print "wrong OUTCAR (nuber of atoms not correct"

        ###################################################################################
        # print parametrization
        ###################################################################################
        with printoptions(precision=2, suppress=True):
            self.print_parametrization()

        ###################################################################################
        # print fores
        ###################################################################################
        emp=" "*19+"|"
        emps=" "*20+"|"
        empss=" "*18+"|"
        #a  = " lon:"+emp+"tox:"+emps+"toy:"+empss+"toz:"+empss+"NUMAT: "+"HES:"+emp+"NEI:   "+"DIF:"+empss+"pos:"
        a  = " lon:"+emp+"tox:"+emp+"tix:"+emp+"SUM:"+emp+"NUMAT: "+"VASP:"+emp+"SHELL: "+"|DIFF:"+empss+"pos:"
        #a  = "forcestotal:                     forceslong:                   forcestran_out:               forcestran_in:                   VASP:                          DIFF:"
        #a = "forcestotal:           forceslong:           forcestan:          VASP:                 DIFF:                HesseRest:"
        b =  "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"

        # print heading:
        #  lon:                   |tox:                   |tix:                   |SUM:                   |NUMAT: VASP:                   |SHELL: |DIFF:                  |pos:
        print b
        print a
        print b

        def printfield(printforces,start,stop,color):
            #print "start:",start,"stop:",stop,"shape:",printforces[printforces[:,16].argsort()].shape
            #print "--------------"
            #A = printforces[printforces[:,16].argsort()][start:stop]
            #A = printforces[printforces[:,16].argsort()][start:stop]
            A = printforces[printforces[:,16].argsort()][int(start):int(stop)]
            B = A[np.lexsort((A[:, 22], A[:, 21], A[:,20], A[:,16]))]
            #utils.printarray(B,color=color,decimals=2)
            utils.printarray(B,color=color,decimals=[2,2,2,2,2,2,2,2,2,1,1,1,0,1,1,1,0,1,1,1,2,2,2,2],prepend=['','','','','','','','','','|','','','|','','','','|','|','','','|'])

        # print forces on first atom
        printfield(printforces,0,1,utils.printgreen)

        # print forces on other atoms
        print_shells = self.shells
        if self.ss != False:
            print_shells = np.arange(0,self.ss)
        print "print_shells:", self.ss,print_shells
        for shell in print_shells:
            colorshell = [ 'yellow', 'red', 'blue', 'green', 'yellow' ]
            #print utils.printred("hallo")
            #print eval("utils.print"+str(colorshell[shell])+"(\"hallo\")")
            #printfield(printforces,np.sum(self._info_nnatoms_pc[:shell-1])+1,np.sum(self._info_nnatoms_pc[:shell])+1,utils.printred)
            printfield(printforces,np.sum(self._info_nnatoms_pc[:shell-1])+1,np.sum(self._info_nnatoms_pc[:shell])+1,eval("utils.print"+str(colorshell[shell])))

        self.print_energies()
        return

    def print_info_loop1(self):
        if self.verbose > 2:
        #if (abs(self.elong) > 1e-6 and self.verbose == 3) or self.verbose == 4: # elong is not yet defined
            #print "self.verbose",self.verbose
            #sys.exit()
            potaddparams_write = self.u_potaddparam
            if self.u_potadd == 'spline':
                potaddparams_write = self.u_potaddparam.get_coeffs()[[0,-1]]
            print "---------------------------------------------------------------"
            print utils.printred("++INDI"),utils.printred(utils.number_to_string_of_certain_length(self.indi,0,3)),\
            utils.printred("  shell:"),utils.printred(self.shell),"pottype:",utils.string_add_spaces(self.u_pottype,4,False),\
            utils.printgreen("nndist:"),utils.printgreen(utils.number_to_string_of_certain_length(self.nndisteq,5,7))#,\
            #" params:",self.u_potparam," potadd:",self.u_potadd," potaddparams:",potaddparams_write
        return

    def print_info_loop2_lon_vec(self):
        if self.verbose > 2:
            ene_print = utils.number_to_string_of_certain_length(self.elong,4,6)
            eneadd_print = utils.number_to_string_of_certain_length(self.elongadd,4,6)
            #if (abs(self.elong) > 1e-6 and self.verbose == 3) or self.verbose == 4:
            #    print '>>indi',utils.number_to_string_of_certain_length(self.indi,0,3),\
            #    "  shell:",self.shell,'indj',utils.string_add_spaces(self.indj,3),\
            #    "|lon("+utils.string_add_spaces(self._u_pottype_curr,4)+")|:",\
            #    utils.printgreen(utils.number_to_string_of_certain_length(self.longvecnorm,5,7)),\
            #    utils.printblue("ene:"),\
            #    utils.printblue(ene_print),\
            #    "posi:",self.sc.rcar[self.indi],"posj:",self.posj,"|i-j|-|lon|:",\
            #    np.linalg.norm(self.sc.rcar[self.indi]-self.posj)-np.linalg.norm(self.longvec)\

            #if (abs(self.elong) > 1e-6 and self.verbose == 3) or self.verbose == 4:
            #    print '>>indi',utils.number_to_string_of_certain_length(self.indi,0,3),\
            #    "  shell:",self.shell,'indj',utils.string_add_spaces(self.indj,3),\
            #    "|lon("+utils.string_add_spaces(self._u_pottype_curr,4)+")|:",\
            #    utils.printgreen(utils.number_to_string_of_certain_length(self.longvecnorm,5,7)),\
            #    utils.printblue("ene:"),\
            #    utils.printblue(ene_print),\
            #    "posi:",self.sc.rcar[self.indi],"posj:",self.posj,"flong:",\
            #    self.flong,"flongadd:",self.flongadd

            if (abs(self.elong) > 1e-6 and self.verbose == 3) or self.verbose == 4:
                print '>>indi',utils.printyellow(utils.number_to_string_of_certain_length(self.indi,0,3)),\
                "  shell:",self.shell,'indj',utils.printpink(utils.string_add_spaces(self.indj,3)),\
                "|lon("+utils.string_add_spaces(self._u_pottype_curr,4)+")|:",\
                utils.printgreen(utils.number_to_string_of_certain_length(self.longvecnorm,5,7)),\
                utils.printblue("ene:"),\
                utils.printblue(ene_print),utils.printblue(eneadd_print),\
                "posj:",self.posj,"flong:",\
                self.flong,"flongadd:",self.flongadd
        return

    def save_information_of_run_initialize(self):
        ''' saves every longvec, toxvec, .... '''
        try:
            if self.params.save_vecs_to_file_for_DOS != True:
                return
        except AttributeError:
            return


        ###########################################
        # tox
        ###########################################
        if os.path.isfile("tvecoutall.dat") == True:
            self._info_tvecoutall = np.loadtxt("tvecoutall.dat")
        else:
            self._info_tvecoutall = np.array([[0.0],[0.0]])
        if self._info_tvecoutall.shape[0] == 2:
            self._info_tvecoutall = np.array([[0.0],[0.0]])
        self._info_tvecoutall = self._info_tvecoutall.reshape((len(self._info_tvecoutall),1))


        ###########################################
        # tix
        ###########################################
        if os.path.isfile("tvecinsall.dat") == True:
            self._info_tvecinsall = np.loadtxt("tvecinsall.dat")
        else:
            self._info_tvecinsall = np.array([[0.0],[0.0]])
        if self._info_tvecinsall.shape[0] == 2:
            self._info_tvecinsall = np.array([[0.0],[0.0]])
        self._info_tvecinsall = self._info_tvecinsall.reshape((len(self._info_tvecinsall),1))

        ###########################################
        # lon / lonx / lony / lonz
        ###########################################
        if os.path.isfile("tveclonall.dat") == True:
            self._info_tveclonall = np.loadtxt("tveclonall.dat")
            self._info_tveclonallx = np.loadtxt("tveclonallx.dat")
        else:
            self._info_tveclonall = np.array([[0.0],[0.0]])
            self._info_tveclonallx = np.array([[0.0],[0.0]])
        if self._info_tveclonall.shape[0] == 2:
            self._info_tveclonall = np.array([[0.0],[0.0]])
        if self._info_tveclonallx.shape[0] == 2:
            self._info_tveclonallx = np.array([[0.0],[0.0]])
        self._info_tveclonall = self._info_tveclonall.reshape((len(self._info_tveclonall),1))
        self._info_tveclonallx = self._info_tveclonallx.reshape((len(self._info_tveclonallx),1))

        ###########################################
        # lonangle
        ###########################################
        if os.path.isfile("tveclonallangle.dat") == True:
            self._info_tveclonallangle = np.loadtxt("tveclonallangle.dat")
        else:
            self._info_tveclonallangle= np.array([[0.0],[0.0]])
        if self._info_tveclonallangle.shape[0] == 2:
            self._info_tveclonallangle = np.array([[0.0],[0.0]])
        self._info_tveclonallangle = self._info_tveclonallangle.reshape((len(self._info_tveclonallangle),1))

        ###########################################
        # lonangle vs norm
        ###########################################
        self._info_tveclonallanglevsnorm = False
        self._info_lonvec = False
        return

    def save_information_of_run_add(self,contr = False):
        ''' adds data to tvecoutallARRAY, tvecinsallARRAY, tveclonallARRAY '''
        try:
            if self.params.save_vecs_to_file_for_DOS != True:
                return
        except AttributeError:
            return

        if contr == 'lon':
            if self.longvecnorm > 0.00001:
                self._info_tveclonall = np.concatenate((self._info_tveclonall,np.array([[self.longvecnorm]])),axis=0)
                self._info_tveclonallx = np.concatenate((self._info_tveclonallx,np.array([[self.longvec[0]-self.vec0[0]]])),axis=0)
                self._info_tveclonallx = np.concatenate((self._info_tveclonallx,np.array([[self.longvec[1]-self.vec0[1]]])),axis=0)
                self._info_tveclonallx = np.concatenate((self._info_tveclonallx,np.array([[self.longvec[2]-self.vec0[2]]])),axis=0)
                self._info_tveclonallangle = np.concatenate((self._info_tveclonallangle,np.array([[self.longvecangle]])),axis=0)

                self._info_tveclonallanglevsnorm = utils.append_row_to_2d_array( \
                        inarray = self._info_tveclonallanglevsnorm, \
                        addrow=[self.longvecnorm,self.longvecangle])
                self._info_lonvec = utils.append_row_to_2d_array(\
                        inarray = self._info_lonvec, \
                        addrow = list([self.longvec-self.vec0]))

        #if contr == 'to':
            #if np.linalg.norm(self.tvecout) > 0.00001:

            # to self.tvecoutnorm was not defined (AttributeError)
            #if self.tvecoutnorm > 0.00001:
            #    #print "shape;",tvecoutall.shape
            #    tvecoutall = np.concatenate((self._info_tvecoutall,np.array([[self.tvecoutnorm]])),axis=0)
            #    tvecoutall = np.concatenate((self._info_tvecoutall,np.array([[-self.tvecoutnorm]])),axis=0)
        return

    def save_information_of_run_save(self):
        ''' saves tvecoutall.dat, tvecinsall.dat, tveclonall.dat '''
        try:
            if self.params.save_vecs_to_file_for_DOS != True:
                return
        except AttributeError:
            return

        np.savetxt("tveclonall.dat",self.tveclonall,fmt="%.6f")
        np.savetxt("tvecoutall.dat",self.tvecoutall,fmt="%.6f")
        np.savetxt("tvecinsall.dat",self.tvecinsall,fmt="%.6f")

        if False:
            if os.path.isfile("energy_long") != True:
                open("energy_long", 'a').close()
            utils.write_inputdata("energy_long","step "+str(round(energymevlong,4)))
            if os.path.isfile("energy_tox") != True:
                open("energy_tox", 'a').close()
            utils.write_inputdata("energy_tox","step "+str(round(energymevtopx,4)))

    ######################################################################################
    # CHECKS
    ######################################################################################
    def check_forces_sum(self):
        ''' checks if sum of forces is 0 '''
        def fehlermelungforces(tol=None):
            print "tolerance:",str(tol)
            print "sum:",np.sum(self.forces)
            print "sum forceslong:",np.sum(self.forceslong)
            print "sum forcestrantopx:",np.sum(self.forcestrantopx)
            print "sum forcestrantopy:",np.sum(self.forcestrantopy)
            print "sum forcestrantopz:",np.sum(self.forcestrantopz)
            print "sum forcestrantipx:",np.sum(self.forcestrantipx)
            print "sum forcestrantipy:",np.sum(self.forcestrantipy)
            print "topx:",np.sum(self.forcestrantopxcheck)
            print "topy:",np.sum(self.forcestrantopycheck)
            print "topz:",np.sum(self.forcestrantopzcheck)
            print "tipx:",np.sum(self.forcestrantipxcheck)
            print "tipy:",np.sum(self.forcestrantipycheck)
            print "tipz:",np.sum(self.forcestrantipzcheck)
            for idxftix,ftix in enumerate(self.forcestrantipxcheck):
                if np.sum(ftix) != 0.0:
                    print "idxftix:",idxftix,np.sum(ftix)
                    print ftix
            sys.exit("sum of forces is not 0")

        def fehlermelungenergy():
            print "self.energylong:",self.energylong,    " sum:",np.sum(self.energylong),     " meV:",self.energymevlong," meV/2:",self.energymevlong/2.
            print "self.energytopx:",self.energytrantopx," sum:",np.sum(self.energytrantopx), " meV:",self.energymevtopx," meV/2:",self.energymevtopx/2.
            print "self.energytopy:",self.energytrantopy," sum:",np.sum(self.energytrantopy), " meV:",self.energymevtopy," meV/2:",self.energymevtopy/2.
            print "self.energytopz:",self.energytrantopz," sum:",np.sum(self.energytrantopz), " meV:",self.energymevtopz," meV/2:",self.energymevtopz/2.
            print "self.energytipx:",self.energytrantipx," sum:",np.sum(self.energytrantipx), " meV:",self.energymevtipx," meV/2:",self.energymevtipx/2.
            print "self.energytipy:",self.energytrantipy," sum:",np.sum(self.energytrantipy), " meV:",self.energymevtipy," meV/2:",self.energymevtipy/2.
            print "self.energytipz:",self.energytrantipz," sum:",np.sum(self.energytrantipz), " meV:",self.energymevtipz," meV/2:",self.energymevtipz/2.
            print "-----------------------------------------------------------------------"
            print "SUM (mev) (all/2):",np.sum(self.energymevlong)/2,"+",np.sum(self.energymevtopx)/2.,"+",np.sum(self.energymevtopy)/2.,"+ ... =",self.energymev
            #self.energymevlong = np.sum(self.energylong)    /2.*1000/(self.numberofatoms-1)
            #self.energymevtopx = np.sum(self.energytrantopx)/2.*1000/(self.numberofatoms-1)
            #self.energymevtopy = np.sum(self.energytrantopy)/2.*1000/(self.numberofatoms-1)
            #self.energymevtopz = np.sum(self.energytrantopz)/2.*1000/(self.numberofatoms-1)
            #self.energymevtipx = np.sum(self.energytrantipx)/2.*1000/(self.numberofatoms-1)
            #self.energymevtipy = np.sum(self.energytrantipy)/2.*1000/(self.numberofatoms-1)
            #self.energymevtipz = np.sum(self.energytrantipz)/2.*1000/(self.numberofatoms-1)

            #self.energymev     = np.sum(self.energy    )    *1000/(self.numberofatoms-1)  # is already devided by 2

            #self.energy =  \
            #        np.sum(self.energylong)/2. +\
            #        np.sum(self.energytrantopx)/2. +\
            #        np.sum(self.energytrantopy)/2. +\
            #        np.sum(self.energytrantopz)/2. +\
            #        np.sum(self.energytrantipx)/2. +\
            #        np.sum(self.energytrantipy)/2. +\
            #        np.sum(self.energytrantipz)/2.

            #self.energymev     = np.sum(self.energy    )    *1000/(self.numberofatoms-1)  # is already devided by 2
            #self.energymevlong = np.sum(self.energylong)    /2.*1000/(self.numberofatoms-1)
            #self.energymevtopx = np.sum(self.energytrantopx)/2.*1000/(self.numberofatoms-1)
            #self.energymevtopy = np.sum(self.energytrantopy)/2.*1000/(self.numberofatoms-1)
            #self.energymevtopz = np.sum(self.energytrantopz)/2.*1000/(self.numberofatoms-1)
            #self.energymevtipx = np.sum(self.energytrantipx)/2.*1000/(self.numberofatoms-1)
            #self.energymevtipy = np.sum(self.energytrantipy)/2.*1000/(self.numberofatoms-1)
            #self.energymevtipz = np.sum(self.energytrantipz)/2.*1000/(self.numberofatoms-1)


        tolerance = 0.0001
        if abs(np.sum(self.forceslong)) >= tolerance:
            print "fehler,forceslong",abs(np.sum(self.forceslong))
            fehlermelungforces(tolerance)
        if abs(np.sum(self.forces)) >= tolerance:
            print "fehler,forces"
            fehlermelungforces(tolerance)


        if abs(np.sum(self.forcestrantopx)) >= tolerance:
            print "fehler,forcestrantopx"
            fehlermelungforces(tolerance)
        if abs(np.sum(self.forcestrantopy)+np.sum(self.forcestrantopz)) >= 0.000001:
            print "fehler,forcestrantop{y,z}",abs(np.sum(self.forcestrantopx)-np.sum(self.forcestrantopy))
            fehlermelungforces(tolerance)

        if abs(np.sum(self.forcestrantipx)) >= tolerance:
            print "fehler,forcestrantipx"
            fehlermelungforces(tolerance)
        if abs(np.sum(self.forcestrantipy)) >= tolerance:
            print "fehler,forcestrantipy"
            fehlermelungforces(tolerance)



        if self.energy < -1.0:
            print "energy:",self.energy,"in meV:",self.energymev
            self.print_forces()
            fehlermelungenergy()
            sys.exit("Potential energy is negative! (probably there is a jump / change of poitions from equilibrium structure)")
        if self.verbose > 2:
            fehlermelungenergy()
        return

    def checks_for_fcc(self):
        ''' checks if found neighbors ... are correct '''
        listcheck = [
            [ 1,self.numberofatoms,12 ],
            [ 2,self.numberofatoms, 6 ],
            [ 3,self.numberofatoms,24 ]]
        fccatoms = [ 32, 108, 256 ]
        if self.numberofatoms not in fccatoms:
            sys.exit("this does not seem to be fcc, inclde check ....or dele this")

        for listchk in listcheck:
            if self.shell == listchk[0]:
                if self.positions.shape[0] == listchk[1] and \
                        self.NNlist.shape[0] != listchk[2]:
                    sys.exit("did not find "+str(listchk[2])+" NN atoms in \
                            shell "+str(listchk[0])+" (assumed fcc structure)")
        return

    ######################################################################################
    # ENERGY FORCES CALCULATION
    ######################################################################################
    def get_or_load_parameters_lon_to_ti_depending_on_shell(self):
        ''' loads the parameters specified in ./calculate_energy_and_forces depending on shell'''
        ################################
        # long
        ################################
        self.u_pottype      = eval("self.params.u"+str(self.shell)+"nn_pottype")
        self.u_potparam     = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_potparam"))



        ################################
        # longadd
        ################################
        try:
            self.u_potadd       = eval("self.params.u"+str(self.shell)+"nn_potadd")
        except AttributeError:
            self.u_potadd       = False

        if self.u_potadd == 'spline':
            try:
                self.pot_parametrize
            except AttributeError:
                self.pot_parametrize = pot_parametrize.forcesneighbors()
                #self.pot_parametrize.loadDatapd("/Users/glensk/Dropbox/Understand_distributions/displacements_dense/Ir/2x2x2sc_data_3x3x3kp.pkl")
                self.pot_parametrize.loadDatapd(eval("self.params.u"+str(self.shell)+"nn_potaddpath"))
            #print "||||aaaa;",self.pot_parametrize.data[1].lon.quer.dos.morse
            self.u_potaddparam = eval('self.pot_parametrize.fit[self.shell].'+self.params.u1nn_potaddparam)
            #print "kkkkkk:",self.u_potaddparam
        else:
            try:
                self.u_potaddparam  = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_potaddparam"))
            except AttributeError:
                self.u_potaddparam  = False

        if type(self.u_pottype) != bool and type(self.u_potparam) == bool:
            sys.exit("define u_potparam")
        if type(self.u_pottype) == bool and type(self.u_potparam) != bool:
            sys.exit("define u_pottype")
        if type(self.u_potadd) != bool and type(self.u_potaddparam) == bool:
            sys.exit("define u_potaddparam")
        if type(self.u_potadd) == bool and type(self.u_potaddparam) != bool:
            sys.exit("define u_potadd")

        ############################################################################
        ############################################################################
        # qubus corrections
        ############################################################################
        ############################################################################
        try:
            self.u_qubus        = eval("self.params.u"+str(self.shell)+"nn_qubus")
        except AttributeError:
            self.u_qubus       = False

        if self.u_qubus       != False:
            try:
                #self.fmlx
                self.Fx
            except AttributeError:
                #self.pot_parametrize = pot_parametrize.forcesneighbors()

            #def load_qubus_npz(filename):
                folder = self.u_qubus
                print utils.printred("loading "+folder+"qubus.F.npz")
                var = np.load(folder+'qubus.F.npz')
                self.x = var['x']
                self.y = var['y']
                self.z = var['z']

                #self.fmlx = var['fmlx']
                #self.fmly = var['fmly']
                #self.fmlz = var['fmlz']
                #self.emlx = var['emlx']
                #self.emly = var['emly']
                #self.emlz = var['emlz']

                self.Fx = var['Fx']
                self.Fy = var['Fy']
                self.Fz = var['Fz']
                self.Ex = var['Ex']
                self.Ey = var['Ey']
                self.Ez = var['Ez']

                print utils.printred("loading "+folder+"qubus.fml.npz")
                var = np.load(folder+'qubus.fml.npz')
                self.fmlx = var['fmlx']
                self.fmly = var['fmly']
                self.fmlz = var['fmlz']
                self.emlx = var['emlx']
                self.emly = var['emly']
                self.emlz = var['emlz']

        ################################################################
        # pot2 is not necessary
        ################################################################
        # ################################
        # # long2
        # ################################
        # try:
        #     self.u_pot2type      = eval("self.params.u"+str(self.shell)+"nn_pot2type")
        # except AttributeError:
        #     self.u_pot2type      = False

        # try:
        #     self.u_pot2param     = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_pot2param"))
        # except AttributeError:
        #     self.u_pot2param     = False

        # ################################
        # # long2add
        # ################################
        # try:
        #     self.u_pot2add       = eval("self.params.u"+str(self.shell)+"nn_pot2add")
        # except AttributeError:
        #     self.u_pot2add       = False

        # try:
        #     self.u_pot2addparam  = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_pot2addparam"))
        # except AttributeError:
        #     self.u_pot2addparam  = False

        # if type(self.u_pot2type) != bool and type(self.u_pot2param) == bool:
        #     sys.exit("define u_pot2param")
        # if type(self.u_pot2type) == bool and type(self.u_pot2param) != bool:
        #     sys.exit("define u_pot2type")
        # if type(self.u_pot2add) != bool and type(self.u_pot2addparam) == bool:
        #     sys.exit("define u_pot2addparam")
        # if type(self.u_pot2add) == bool and type(self.u_pot2addparam) != bool:
        #     sys.exit("define u_pot2add")


        ################################
        # tox
        ################################
        try:
            self.u_topx = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_topx"))
        except AttributeError:
            self.u_topx = False
        try:
            self.u_topy = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_topy"))
        except AttributeError:
            self.u_topy = False
        try:
            self.u_topz = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_topz"))
        except AttributeError:
            self.u_topz = False

        if  type(self.u_topx) == bool and \
            type(self.u_topy) == bool and \
            type(self.u_topz) == bool:
            self.u_top = False
        else:
            self.u_top = True

        ################################
        # tix
        ################################
        try:
            self.u_tipx = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_tipx"))
        except AttributeError:
            self.u_tipx = False
        try:
            self.u_tipy = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_tipy"))
        except AttributeError:
            self.u_tipy = False
        try:
            self.u_tipz = utils.string_to_list(eval("self.params.u"+str(self.shell)+"nn_tipz"))
        except AttributeError:
            self.u_tipz = False

        if  type(self.u_tipx) == bool and \
            type(self.u_tipy) == bool and \
            type(self.u_tipz) == bool:
            self.u_tip = False
        else:
            self.u_tip = True

        if self.u_top == False and self.u_tip == False:
            self.u_toti = False
        else:
            self.u_toti = True
        return

    def declare_empty_arrays(self):
        ''' empty arrays ...
            eigentlich muesste man dies fuer jede schale seperat definieren '''
        #self.positions = np.copy(self.crystal.rcar)
        #self.positions = np.copy(self.sc.rcar[:self.numberofatoms])
        self.positions = np.copy(self.scneverchange.rcar[:self.numberofatoms])

        #self.energylong =          np.zeros((self.positions.shape[0],1))
        self.energylong =          np.zeros(self.positions.shape[0])
        self.forceslong =          np.zeros((self.positions.shape[0],3))

        #self.energytrantopx =      np.zeros((self.positions.shape[0],1))
        self.energytrantopx =      np.zeros(self.positions.shape[0])
        self.forcestrantopx =      np.zeros((self.positions.shape[0],3))
        #self.energytrantopy =      np.zeros((self.positions.shape[0],1))
        self.energytrantopy =      np.zeros(self.positions.shape[0])
        self.forcestrantopy =      np.zeros((self.positions.shape[0],3))
        #jself.energytrantopz =      np.zeros((self.positions.shape[0],1))
        self.energytrantopz =      np.zeros(self.positions.shape[0])
        self.forcestrantopz =      np.zeros((self.positions.shape[0],3))

        #self.energytrantipx =      np.zeros((self.positions.shape[0],1))
        self.energytrantipx =      np.zeros(self.positions.shape[0])
        self.forcestrantipx =      np.zeros((self.positions.shape[0],3))
        #self.energytrantipy =      np.zeros((self.positions.shape[0],1))
        self.energytrantipy =      np.zeros(self.positions.shape[0])
        self.forcestrantipy =      np.zeros((self.positions.shape[0],3))
        #self.energytrantipz =      np.zeros((self.positions.shape[0],1))
        self.energytrantipz =      np.zeros(self.positions.shape[0])
        self.forcestrantipz =      np.zeros((self.positions.shape[0],3))

        self.forcestrantopxcheck = np.zeros((self.positions.shape[0],self.positions.shape[0],3))
        self.forcestrantopycheck = np.zeros((self.positions.shape[0],self.positions.shape[0],3))
        self.forcestrantopzcheck = np.zeros((self.positions.shape[0],self.positions.shape[0],3))
        self.forcestrantipxcheck = np.zeros((self.positions.shape[0],self.positions.shape[0],3))
        self.forcestrantipycheck = np.zeros((self.positions.shape[0],self.positions.shape[0],3))
        self.forcestrantipzcheck = np.zeros((self.positions.shape[0],self.positions.shape[0],3))

        # sublists for every shell (first element just stays empty self.longvec_all[1] ist data
        self.longvec_all = [ False for i in np.arange(len(self.shells)+1) ]
        self.longvecnorm_all = [ np.array([]) for i in np.arange(len(self.shells)+1) ]
        return

    def determine_shells(self):
        ''' determine over which shells of nearest neighbours to loop '''
        self.shells = [1]
        checkshells = range(2,6)  # check up to 5th shell
        for i in checkshells:
            try:
                if type(eval("self.params.u"+str(i)+"nn_potparam")) != bool:
                    self.shells = range(1,i+1) # range(1,3) = [1,2]
            except AttributeError:
                return self.shells
        return self.shells


    def get_lon_vec(self):
        ''' gets longvec, longvecnorm, self.vec0, self.vec0norm '''
        self.longvec = np.copy(self.posj)
        #self.longvecnorm = np.copy(np.linalg.norm(self.posj))
        self.vec0 = self.sc0.rcar[self.indj]
        self.disp = self.vec0 - self.longvec
        self.longvecangle = utils.anglevec(self.longvec,self.vec0)
        #self.vec0norm = np.copy(np.linalg.norm(self.vec0))

        # append longvec and longvecnorm to lists
        #self.longvec_all[self.shell] = utils.append_row_to_2d_array(inarray = self.longvec_all[self.shell],
        #        addrow=[self.longvec[0],self.longvec[1],self.longvec[2]])
        self.longvecnorm_all[self.shell] = np.append(self.longvecnorm_all[self.shell],self.longvecnorm)
        return

    def get_lon_energy_forces(self):
        ''' returnes the energy and the fores between two atoms '''
        if  type(self.u_pottype) == bool:
            return
        #if  type(self.u_pottype) == bool and \
        #    type(self.u_pot2type) == bool:
        #        return

        # adjust potential if pot2param are defined
        self._u_pottype_curr      = self.u_pottype
        self._u_potparam_curr     = self.u_potparam
        self._u_potadd_curr       = self.u_potadd
        self._u_potaddparam_curr  = self.u_potaddparam
        ##################################################################################
        # pot2 is not necessary
        ##################################################################################
        #if type(self.u_pot2param) != bool:
        #    if self.u_pot2type != 'm' and self.u_pot2type != 'mc1' \
        #            and self.u_pot2type != 'poly':
        #        print "u_pot2:",self.u_pot2type
        #        sys.exit("u_pot2type has to be m or mc1 or poly")
        #    if self.longvecnorm > self.nndisteq:
        #        self._u_pottype_curr      = self.u_pot2type
        #        self._u_potparam_curr     = self.u_pot2param
        #        self._u_potadd_curr       = self.u_pot2add
        #        self._u_potaddparam_curr  = self.u_pot2addparam




        #print "self.u_pottype:----------",self.u_pottype,self.u_potparam
        #self.longvecproj = utils.project_vector(self.longvec,self.vec0)
        # projecting the longvec onthe the original longvec is not good!
        #print "||:",self.longvec,self.longvecproj,np.linalg.norm(self.longvec),np.linalg.norm(self.longvecproj)
        self.elong, self.flong = self.getefvec(self.longvec,self._u_potparam_curr,pot = self._u_pottype_curr,vecnorm = self.longvecnorm )
        #self.elong, self.flong = self.getefvec(self.longvecproj,self._u_potparam_curr,pot = self._u_pottype_curr,vecnorm = np.linalg.norm(self.longvecproj))
        #print "yo",self.u_potadd,self.u_potaddparam
        self.elongadd, self.flongadd = self.getefvec(self.longvec,self._u_potaddparam_curr,pot = self._u_potadd_curr, vecnorm = self.longvecnorm)

        ###############################################################################
        #$xmgrace -block  /Users/glensk/Dropbox/Understand_distributions/ti/Ir/31__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD_YES_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda1.0/tests/test1_qubus_correction/dUdLref  -bxy 1:7 -block /Users/glensk/Dropbox/Understand_distributions/ti/Ir/31__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD_YES_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda1.0/dUdL -bxy 1:7

        if self.u_qubus  != False:
            # wenn positiv self.vec0 und self.longvec bekommt man
            # die kraft auf das atom welches man anschaut
            #
            # wenn negativ self.vec0 und self.longvec bekommt man
            # die kraft auf das ausgelenkte atom
            self.Fq,self.Eq = self.map_vec_back_to_first_quadrant(
                    vec0_curr = -self.vec0,
                    longvec_curr = -self.longvec)
            #print "vec;",vec,self.vec0,self.vec0-vec
            #vecd = self.vec0-self.longvec
            #fq,eq = self.get_ef_rest(vecd[0],vecd[1],vecd[2])
            #print vecd,"eq,fq:",eq,fq,"fvec:",fvec
            if self.verbose > 2:
                if abs(np.linalg.norm(self.Fq)) > 0.0000000001:
                    print ""
                    print "        ",utils.printpink(self.indj),"|||||",self.vec0,self.longvec,"============> Fq:",self.Fq,"Eq:",utils.printblue(self.Eq)
        #if self.shell == 2 or self.shell == 3:
        #    self.elong = self.elong/2.
        #    self.elongadd = self.elongadd/2.
        #    self.flong = self.flong/2.
        #    self.flongadd = self.flongadd/2.
        #print "kk1:",self.energylong[self.indi]
        #print "kk2:",self.elong
        #print "kk3:",self.elongadd
        #sys.exit()
        #self.oldway = False
        #if self.oldway == True:
        if self.u_qubus != False:
            #sys.exit("mit? qubus_correction?")
            #print "-->>>>",self.indi,"self.forceslong[self.indi]",self.forceslong[self.indi]
            self.forceslong[self.indi] = self.forceslong[self.indi] + self.Fq
            self.energylong[self.indi] = self.energylong[self.indi] + np.linalg.norm(self.Eq)
        else:
            #sys.exit("ohne qubus_correction?")
            self.forceslong[self.indi] = self.forceslong[self.indi] + self.flong + self.flongadd
            self.energylong[self.indi] = self.energylong[self.indi] + self.elong + self.elongadd

        return

    def getefvec(self, vec, params, pot = False, vecnorm = False):
        ''' returnes energies forces '''
        if type(params) == bool:
            return 0.0,np.array([0.0,0.0,0.0])
        if type(pot) == bool:
            print "pot:",pot
            print "self.indi:",self.indi
            print "self.indj:",self.indj
            print "self.shell:",self.shell,"self.nndisteq:",self.nndisteq
            sys.exit(" you have to define the pot as m, mc1, poly, ...")
        if vecnorm == False:
            vecnorm=np.linalg.norm(vec)
        if abs(vecnorm) <= 0.000000000000000000000001:
            return 0.0,np.array([0.0,0.0,0.0])
        #print "vn:",vecnorm,params,pot
        enorm,fnorm = self.getef(vecnorm, params, pot)
        fvec = vec/vecnorm*fnorm


        return enorm, fvec

    def getef(self, vecnorm, params, pot = False):
        ''' returnes energies forces '''
        if type(pot) == bool:
            sys.exit(" you have to define the pot as m, mc1, poly, ...")
        #print "params:",params,type(params)
        if type(params) == bool:
            return 0.0,0.0
        #print "    vecnorm:",vecnorm
        #print "    params:",params
        #print "    pot:",pot
        if pot == 'poly': # np.poly is assumed
            #f = np.polynomial.Polynomial(params)(vecnorm)
            #fNOWWRONG = np.poly1d(params[::-1])(vecnorm)  # assume poly1d to be quicker
            #  np.poly1d([1, 2, 3]) == 1 x^2 + 2 x + 3
            #  np.polynomial.Polynomials([1,2,3]) == 1 + 2x + 3x^2
            #eNOWWRONG = np.polyint(np.poly1d(params[::-1]))(vecnorm)
            e = np.poly1d(params[::-1])(vecnorm)
            f = np.polyder(np.poly1d(params[::-1]))(vecnorm)
            return e,f
        #if pot == False: # or pot == '135' or pot == '1357' or pot == '13579':
        #    if pot == '135'     and len(params) != 3: sys.exit('135 but not 3')
        #    if pot == '1357'    and len(params) != 4: sys.exit('1357 but not 4')
        #    if pot == '13579'   and len(params) != 5: sys.exit('1359 but not 5')
        #    if len(params) == 3:
        #        e =     Trans135(vecnorm,*params)
        #        f = Trans135_der(vecnorm,*params)
        #    elif len(params) == 4:
        #        e =     Trans1357(vecnorm,*params)
        #        f = Trans1357_der(vecnorm,*params)
        #    elif len(params) == 5:
        #        #print "in 13579:",vecnorm,params
        #        e =     Trans13579(vecnorm,*params)
        #        f = Trans13579_der(vecnorm,*params)
        #    return e,f
        function = None
        functionder = None
        if pot == 'i' or pot == 'inverse':
            function = inversepot
            functionder = inversepot_derivative
        if pot == 'l':
            function = LJ
            functionder = LJ_derivative
        if pot == 'm' or pot == 'morse':
            function = Morse
            functionder = Morse_derivative
        if pot == 'mc1':
            function = mc1
            functionder = mc1_derivative
        if pot == 'spline':
            #import pot_parametrize
            #ka = pot_parametrize.forcesneighbors()
            #ka.loadDatapd("2x2x2sc_data_3x3x3kp.pkl")
            #print "self,shell:",self.shell
            # for fcc we need/want quer for long parametrization
            #spl = ka.fit[ self.shell,'lon','quer','dos','morse'].ix['splinedata']
            spl = params
            e_antiderivative_achsenabschnitt = spl.antiderivative()(self.nndisteq)
            e = spl.antiderivative()(vecnorm) - e_antiderivative_achsenabschnitt   # das hier muss ein - sein
            #f = spl.derivative()(vecnorm)
            f = spl(vecnorm)                                # hier wes es ueberhaupt nicht gut wenn eine restkraft bei self.nndisteq auftritt
            #print "--->2:",params,"splinef(",vecnorm,"):",f
            #print "4,e,f:",e,f
            return e,f
        if function == None:
            print "pot:",pot,type(pot)
            sys.exit("pot not recognized")
        #e =    function(longvecnorm,*params)
        #f = functionder(longvecnorm,*params)
        #print "vecnorm:",vecnorm
        #print "params:",params
        e =    function(vecnorm,*params)
        f = functionder(vecnorm,*params)
        #print "--->1:",params,"(",vecnorm,"):",f
        return e,f

    def get_ef_rest(self,x,y=False,z=False):
        ''' '''
        if type(y) == bool and type(z) == bool:
            x,y,z = x[0],x[1],x[2]
        scale = (self.x.size-1)/2./self.x.max()
        #print "self.x:",self.x
        #print "self.x.max()",self.x.max(),self.x.size,self.x.size,(self.x.size-1)/2,scale
        def xyz_to_map_coords(x):
            ''' x = -1.3 ---> out = 0
                x = 0    ---> out = 13
                x = 1.3  ---> out = 26 '''
            return (x+self.x.max())*scale

        coords = np.array([[xyz_to_map_coords(x), xyz_to_map_coords(y),xyz_to_map_coords(z)]])
        coords = coords.T
        #zi = scipy.ndimage.map_coordinates(q.fall, coords, order=2, mode='nearest')

        #fx = map_coordinates(self.fmlx, coords, order=2, mode='nearest')[0]
        #fy = map_coordinates(self.fmly, coords, order=2, mode='nearest')[0]
        #fz = map_coordinates(self.fmlz, coords, order=2, mode='nearest')[0]
        #ex = map_coordinates(self.emlx, coords, order=2, mode='nearest')[0]
        #ey = map_coordinates(self.emly, coords, order=2, mode='nearest')[0]
        #ez = map_coordinates(self.emlz, coords, order=2, mode='nearest')[0]

        fx = map_coordinates(self.Fx, coords, order=2, mode='nearest')[0]
        fy = map_coordinates(self.Fy, coords, order=2, mode='nearest')[0]
        fz = map_coordinates(self.Fz, coords, order=2, mode='nearest')[0]
        ex = map_coordinates(self.Ex, coords, order=2, mode='nearest')[0]
        ey = map_coordinates(self.Ey, coords, order=2, mode='nearest')[0]
        ez = map_coordinates(self.Ez, coords, order=2, mode='nearest')[0]
        return [fx,fy,fz],[ex,ey,ez]

    def map_vec_back_to_first_quadrant(self,vec0_curr,longvec_curr=False,disp_curr=False):
        '''
        ##########################################################################
        # in the end we need:
        # vec0_curr
        # vec_curr  == longvec (correspongs to undisplaced vec0_curr)
        # disp_curr
        #
        # vec0_orig == sollte immer [1.955, 1.955, 0 ] sein
        # disp_orig
        # vec_orig
        ###########################################################################

        vec0_curr: is the vector in the undisplaced reference e.g. [1.955, 1.955, 0.0]
        vecs0 = np.array([ [ 1., 1.,0.], [ -1., 1.,0.], [-1.,-1.,0.], [ 1.,-1.,0.],
                       [ 1., 0.,1.], [ -1., 0.,1.], [-1.,0.,-1.], [ 1.,0.,-1.],
                       [ 0., 1.,1.], [ 0., -1.,1.], [0., -1.,-1.], [ 0., 1.,-1.],
        ])

        Rorig is the original xyz refernze frame with basis [1,0,0]
                                                            [0,1,0]
                                                            [0,0,1]
        R is the reference frme of the current longitudinal vector. Whith this function we
        convert this reference frame to Rorig, get there the energies and the forces, and map
        those back to the current referenze frame
        '''

        Rorig = np.eye(3)
        Rorig[0] = [1.0,0.0,0.0]
        Rorig[1] = [0.0,1.0,0.0]
        Rorig[2] = [0.0,0.0,1.0]
        R = np.zeros((3,3))

        # erste reihe
        R[0,0] = vec0_curr[0]
        if R[0,0] == 0.0:
            R[0,1] = vec0_curr[1]
            R[1,2] = vec0_curr[2]
            R[2] = np.cross(R[0],R[1])

        # zweite reihe
        else:  # R[0,0] ist bereits gefuellt mit einer 1 oder -1
            R[1,1] = vec0_curr[1]
            if R[1,1] != 0.0:
                R[2] = np.cross(R[0],R[1])
            if R[1,1] == 0.0:
                R[1,2] = vec0_curr[2]
                R[2] = np.cross(R[0],R[1])

        ####################### normalize matrix ########
        R = R / np.linalg.norm(R, axis=-1)[:, np.newaxis]

        ######## map vec0_curr to orig #######################
        vec0_orig = np.dot(R,vec0_curr) # sollte immer [1.955,1.955,0.0] sein
        if vec0_orig[2] != 0.0:
            sys.exit("vec0_orig[2] != 0.0")
        if vec0_orig[0] < 0.0:
            sys.exit("vec0_orig[0] < 0.0")
        if vec0_orig[1] < 0.0:
            sys.exit("vec0_orig[1] < 0.0")

        ######## get longvec_curr disp_curr  ###############################
        if type(longvec_curr) == bool and type(disp_curr) == bool:
            sys.exit("either longvec_vec or disp_curr, not both")

        if type(disp_curr) != bool:
            longvec_curr =  disp_to_longvec(vec0_curr,disp_curr)
        else:
            disp_curr = longvec_to_disp(vec0_curr,longvec_curr)

        ######## get longvec_orig disp_orig  ###############################
        longvec_orig = np.dot(R,longvec_curr)
        disp_orig = np.dot(R,disp_curr)

        ####### get forces in original position ############################
        ef_orig = self.get_ef_rest(disp_orig)
        ef_curr_f = np.dot(R.T,ef_orig[0])  # forces
        ef_curr_e = np.dot(R.T,ef_orig[1])  # energies


        # bei der energie macht nur die ef_orig[1] sinn, da positiv,
        # NICHT jecho die zurueckgemappte ef_curr_e
        #print "{disp,longvec,vec0,disp}_curr",disp_curr,"|",longvec_curr,"|",vec0_curr
        #print "{disp,longvec,vec0,disp}_orig",disp_orig,"|",longvec_orig,"|",vec0_orig
        #print "ef_{orig,curr}[0]",ef_orig[0],utils.printred(ef_curr_f[0],ef_curr_f[1],ef_curr_f[2])
        #print "",type(ef_orig[1]),ef_orig[1],"--->",np.sum(ef_orig[1])
        #print "",type(ef_curr_e),ef_curr_e,"-->",np.sum(ef_curr_e)
        #print "",type(ef_curr_f),ef_curr_f
        #print ""
        #print ""
        return ef_curr_f,ef_curr_e

    def sum_up_forces(self):
        ''' takes the sum of all forces '''
        self.energy =  \
                np.sum(self.energylong)/2. +\
                np.sum(self.energytrantopx)/2. +\
                np.sum(self.energytrantopy)/2. +\
                np.sum(self.energytrantopz)/2. +\
                np.sum(self.energytrantipx)/2. +\
                np.sum(self.energytrantipy)/2. +\
                np.sum(self.energytrantipz)/2.

        self.energymev     = np.sum(self.energy    )    *1000/(self.numberofatoms-1)  # is already devided by 2
        self.energymevlong = np.sum(self.energylong)    /2.*1000/(self.numberofatoms-1)
        self.energymevtopx = np.sum(self.energytrantopx)/2.*1000/(self.numberofatoms-1)
        self.energymevtopy = np.sum(self.energytrantopy)/2.*1000/(self.numberofatoms-1)
        self.energymevtopz = np.sum(self.energytrantopz)/2.*1000/(self.numberofatoms-1)
        self.energymevtipx = np.sum(self.energytrantipx)/2.*1000/(self.numberofatoms-1)
        self.energymevtipy = np.sum(self.energytrantipy)/2.*1000/(self.numberofatoms-1)
        self.energymevtipz = np.sum(self.energytrantipz)/2.*1000/(self.numberofatoms-1)


        # make a check if sum of forces is zero!
        self.forces = self.forceslong + \
                self.forcestrantopx + self.forcestrantopy + self.forcestrantopz + \
                self.forcestrantipx + self.forcestrantipy + self.forcestrantipz
        return

    def get_to_ti_direction_vecs_fcc(self):
        self.topdirectionx = np.array([0.0,0.0,0.0])  # those are the directions in which the forces act on to / ti part of vector
        self.topdirectiony = np.array([0.0,0.0,0.0])
        self.topdirectionz = np.array([0.0,0.0,0.0])
        self.topdirectionxnorm = 0.0  # those are the directions in which the forces act on to / ti part of vector
        self.topdirectionynorm = 0.0
        self.topdirectionznorm = 0.0

        self.tipdirection  = np.array([0.0,0.0,0.0])
        self.tipdirectionx = np.array([0.0,0.0,0.0])
        self.tipdirectiony = np.array([0.0,0.0,0.0])
        self.tipdirectionz = np.array([0.0,0.0,0.0])
        self.tipdirectionnorm  = 0.0  # those are the directions in which the forces act on to / ti part of vector
        self.tipdirectionxnorm = 0.0  # those are the directions in which the forces act on to / ti part of vector
        self.tipdirectionynorm = 0.0
        self.tipdirectionznorm = 0.0

        if self.vec0[0] == 0.0 and self.shell == 1:  # dies gilt nur fuer die ersten Nachbarn, fuer die 2NN
                            # muss dies anders definiert werden
                            # hier sollte man sich mehr gedanken ueber symmetrien machen
                            # fuer die 2NN werden die ti und to gleich sein
            # z.b. self.vec0 = [0.0,  2.065,  2.065] 1NN
            # z.b. self.vec0 = [0.0,  2.065, -2.065] 1NN
            # z.b. self.vec0 = [0.0, -2.065,  2.065] 1NN
            # z.b. self.vec0 = [0.0, -2.065, -2.065] 1NN
            #topdirectionx = np.array([vec0[2]/abs(vec0[1]),0.0,0.0])  # only the direction seems to be important not the sign
            # wenn g = (a,b) ist dann ist der senkrechte vektor (-b,a) order (b,-a)
            self.topdirebene    = np.array([0.0, -self.vec0[2],self.vec0[1]])
            self.topdirebeneout = np.cross( self.vec0, self.topdirebene )

            self.topdirectionx = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
            self.topdirectiony = np.array([0.0,-self.vec0[1]/abs(self.vec0[1]),0.0])
            self.topdirectionz = np.array([0.0,0.0,-self.vec0[2]/abs(self.vec0[2])])

            # for tip directions it is probably best to project on self.vec0? or on 1/0/0 .... this should only be a matter of koordinate transformation.
            # therefore lets start with x/y directions as done for the tox toy directions
            #
            # tipdirection{x,y,z} and tipdirection{long,senk} are both possible bases, this is just a matter of coordinate transformation
            self.tipdirectionx = np.array([0.0,-self.vec0[1]/abs(self.vec0[1]),0.0])  # evtl auch hier -vec0[xyz]/abs(vec0)
            self.tipdirectiony = np.array([0.0,0.0,-self.vec0[2]/abs(self.vec0[2])])
            self.tipdirectionz = np.array([0.0,0.0,0.0])  # dieser vector wird nicht gebraucht da er senkrecht auf tipdirection steht

            self.topdirectionxnorm = 1.0
            self.topdirectionynorm = 1.0
            self.topdirectionznorm = 1.0
            self.tipdirectionxnorm = 1.0
            self.tipdirectionynorm = 1.0
            self.tipdirectionznorm = 1.0

            self.tipdirectionlong = self.vec0
            self.tipdirectionsenk = np.array([0.0,self.vec0[1]/abs(self.vec0[1]),-self.vec0[2]/abs(self.vec0[2])])

            self.tipdirection = np.cross(self.vec0,self.topdirectionx)


        if self.vec0[1] == 0.0 and self.shell == 1:  # dies gilt nur fu
            # z.b. self.vec0 = [ 2.065, 0.0,  2.065] 1NN
            # z.b. self.vec0 = [ 2.065, 0.0, -2.065] 1NN
            # z.b. self.vec0 = [-2.065, 0.0,  2.065] 1NN
            # z.b. self.vec0 = [-2.065, 0.0, -2.065] 1NN
            # wenn g = (a,b) ist dann ist der senkrechte vektor (-b,a) order (b,-a)
            self.topdirebene    = np.array([-self.vec0[2],0.0,self.vec0[0]])
            self.topdirebeneout = np.cross( self.vec0, self.topdirebene )

            self.topdirectionx = np.array([0.0,1.0,0.0])
            self.topdirectiony = np.array([-self.vec0[0]/abs(self.vec0[0]),0.0,0.0])
            self.topdirectionz = np.array([0.0,0.0,-self.vec0[2]/abs(self.vec0[2])])
            #tipdirection_out = np.array([vec0[0],0.0,-vec0[2]])
            #tipdirection_in = self.vec0
            self.tipdirectionx = np.array([-self.vec0[0]/abs(self.vec0[0]),0.0,0.0])  # evtl auch hier -vec0[xyz]/abs(vec0)
            self.tipdirectiony = np.array([0.0,0.0,1.0])
            self.tipdirectionz = np.array([0.0,1.0,0.0])  # dieser vector wird nicht gebraucht da er senkrecht auf tipdirection steht
            self.topdirectionxnorm = 1.0
            self.topdirectionynorm = 1.0
            self.topdirectionznorm = 1.0
            self.tipdirectionxnorm = 1.0
            self.tipdirectionynorm = 1.0
            self.tipdirectionznorm = 1.0

            self.tipdirectionlong = self.vec0
            # rechtsdrehend
            if self.vec0[0] > 0.0 and self.vec0[2] < 0.0: self.tipdirectionsenk = np.array([-abs(self.vec0[0]),0.0,-abs(self.vec0[2])])
            if self.vec0[0] < 0.0 and self.vec0[2] > 0.0: self.tipdirectionsenk = np.array([ abs(self.vec0[0]),0.0, abs(self.vec0[2])])
            if self.vec0[0] < 0.0 and self.vec0[2] < 0.0: self.tipdirectionsenk = np.array([-abs(self.vec0[0]),0.0, abs(self.vec0[2])])
            if self.vec0[0] > 0.0 and self.vec0[2] > 0.0: self.tipdirectionsenk = np.array([ abs(self.vec0[0]),0.0,-abs(self.vec0[2])])
            #tipdirectionsenk = np.array([vec0[0]/abs(vec0[0]),0.0,-vec0[2]/abs(vec0[2])])
            self.tipdirection = np.cross(self.vec0,self.topdirectionx)


        if self.vec0[2] == 0.0 and self.shell == 1:
            # z.b. self.vec0 = [ 2.065,  2.065, 0.0] 1NN
            # z.b. self.vec0 = [ 2.065, -2.065, 0.0] 1NN
            # z.b. self.vec0 = [-2.065,  2.065, 0.0] 1NN
            # z.b. self.vec0 = [-2.065, -2.065, 0.0] 1NN
            # wenn g = (a,b) ist dann ist der senkrechte vektor (-b,a) order (b,-a)
            self.topdirebene    = np.array([-self.vec0[1],self.vec0[0],0.0])
            self.topdirebeneout = np.cross( self.vec0, self.topdirebene )


            self.topdirectionx = np.array([0.0,0.0,1.0])
            #topdirectiony = np.array([0.0,1.0,0.0])
            #topdirectionz = np.array([1.0,0.0,0.0])
            self.topdirectiony = np.array([0.0,-self.vec0[1]/abs(self.vec0[1]),0.0])
            self.topdirectionz = np.array([-self.vec0[0]/abs(self.vec0[0]),0.0,0.0])

            self.tipdirectionx = np.array([-self.vec0[0]/abs(self.vec0[0]),0.0,0.0])  # evtl auch hier -vec0[xyz]/abs(vec0)
            self.tipdirectiony = np.array([0.0,-self.vec0[1]/abs(self.vec0[1]),0.0])  # evtl auch hier -vec0[xyz]/abs(vec0)
            self.tipdirectionz = np.array([0.0,0.0,1.0])  # dieser vector wird nicht gebraucht da er senkrecht auf tipdirection steht
            self.topdirectionxnorm = 1.0
            self.topdirectionynorm = 1.0
            self.topdirectionznorm = 1.0
            self.tipdirectionxnorm = 1.0
            self.tipdirectionynorm = 1.0
            self.tipdirectionznorm = 1.0


            # bin von atom 87 ausgegangen: [-2.065,  2.065, 0.0] ist ref
            self.tipdirectionlong = self.vec0
            self.tipdirectionsenk = np.array([-self.vec0[0]/abs(self.vec0[0]),self.vec0[1]/abs(self.vec0[1]),0.0])

            self.tipdirection = np.cross(self.vec0,self.topdirectionx)

        #print "vec0:",self.vec0,"Senkrechtebene:",self.topdirebene,"Senkebeneout:",self.topdirebeneout


        # original reference frame (orf) : x,y,z  === {1,0,0},{0,1,0},{0,0,1}
        # vec0(orf): [0,1,1]      R2 = [[ 0.,  1.,  0.],[0.,  0.,  1.],[ 1.,  0.,  0.]]
        #                         R2 = [[ 0.,  0.,  1.],[0.,  1.,  0.],[ 1.,  0.,  0.]]
        #             x y z       np.dot(R2,[0.,1.,1.]) == [ 1.,  1.,  0.]

        # vec0(orf): [0,-1,1]     R2 = [[ 0.,  -1.,  0.],[0.,  0.,  1.],[ 1.,  0.,  0.]]
        # vec0(orf): [0,1,-1]     R2 = [[ 0.,   1.,  0.],[0.,  0., -1.],[ 1.,  0.,  0.]]
        # vec0(orf): [0,-1,-1]    R2 = [[ 0.,  -1.,  0.],[0.,  0., -1.],[ 1.,  0.,  0.]]





        # vec0(orf): [1,0,1]      R2 = [[ 1.,  0.,  0.],[0.,  0.,  1.],[ 0.,  1.,  0.]]
        #                         R2 = [[ 0.,  0.,  1.],[1.,  0.,  0.],[ 0.,  1.,  0.]]
        #                         np.dot(R2,[1.,0.,1.])  == [ 1.,  1.,  0.]

        # vec0(orf): [-1,0,1]     R2 = [[ -1.,  0.,  0.],[0.,  0.,  1.],[ 0.,  -1.,  0.]]
        #                         R2 = [[ -1.,  0.,  0.],[0.,  0.,  1.],[ 0.,   1.,  0.]]
        #                         R2 = [[ 0.,  0.,  1.],[-1.,  0.,  0.],[ 0.,  1.,  0.]]
        #                         R2 = [[ 0.,  0.,  1.],[-1.,  0.,  0.],[ 0.,  -1.,  0.]]
        #                         np.dot(R2,[-1.,0.,1.])  == [ 1.,  1.,  0.]

        #                         R2 = [[ 1.,  0.,  0.],[0.,  0.,  -1.],[ 0.,   1.,  0.]]
        #                         R2 = [[ 1.,  0.,  0.],[0.,  0.,  -1.],[ 0.,  -1.,  0.]]
        # vec0(orf): [1,0,-1]     R2 = [[ 0.,  0.,  -1.],[1.,  0.,  0.],[ 0.,  -1.,  0.]]
        #                         R2 = [[ 0.,  0.,  -1.],[1.,  0.,  0.],[ 0.,   1.,  0.]]

        # vec0(orf): [-1,0,-1]    R2 = [[ -1.,  0.,  0.],[0.,  0.,  -1.],[ 0.,  -1.,  0.]]



        # current reference frame  : xachese =
        # vec0(orf): [0,1,1]
        #             z x y

        # da wo der vec0 0 ist, kommt auch der entsprechende vektor hin.


        # wenn z == 0 dann ist es einfach, dann x und y switchen und vorzeichen lassen vom vec0
        # vec0(orf): [ 1, 1,0]
        #              x y  z

        #                         R2 = [[ -1.,  0.,  0.],[0.,  1.,  0.],[ 0., 0.,   1.]]
        # vec0(orf): [-1, 1,0]    R2 = [[ 0.,  1.,  0.],[-1.,  0.,  0.],[ 0.,  0.,  1.]]
        #             -x  x z     np.dot(R2,[-1.,1.,0.])  == [ 1.,  1.,  0.]


        # vec0(orf): [-1,-1,0]    R2 = [[ -1.,  0.,  0.],[0.,  -1.,  0.],[ 0., 0.,  1.]]
        #                         R2 = [[ 0.,  -1.,  0.],[-1.,  0.,  0.],[ 0., 0.,  1.]]
        #             -x -y z     np.dot(R2,[-1.,-1.,0.])  == [ 1.,  1.,  0.]

        #                         R2 = [[ 1.,  0.,  0.],[0.,  -1.,  0.],[ 0., 0.,  1.]]
        # vec0(orf): [ 1,-1,0]    R2 = [[ 0.,  -1.,  0.],[1.,  0.,  0.],[ 0.,  0., 1.]]
        #              x -y z     np.dot(R2,[1.,-1.,0.])  == [ 1.,  1.,  0.]



        # @ # # original reference frame (original longvec for qubus)
        # @ # a = orig_tvec = np.array([1.,-1.,0])/np.linalg.norm(np.array([1.,-1.,0]))
        # @ # b = orig_vec0 = np.array([1.,1.,0])/np.linalg.norm(np.array([1.,1.,0]))
        # @ # c = np.cross( orig_tvec, orig_vec0 )


        # @ # if self.vec0[2] == 0.0 and self.shell == 1:
        # @ #     pass


        # @ # if abs(self.longvecnorm - self.nndisteq) > 0.0000001:
        # @ #     print "a:",a
        # @ #     print "b:",b
        # @ #     print "c:",c

        # @ # #print "a:",a,b,c
        # @ # R = np.eye(3)
        # @ # R[0,:] = a
        # @ # R[1,:] = b
        # @ # R[2,:] = c

        # @ # # current reference frame ( simpy xy  reference frame)
        # @ # a2 = self.topdirebene/np.linalg.norm(self.topdirebene)
        # @ # b2 = self.vec0/np.linalg.norm(self.vec0)
        # @ # c2 = np.cross( a2,b2 )

        # @ # a2 = np.array([1.,0.,0.])
        # @ # b2 = np.array([0.,1.,0.])
        # @ # c2 = np.cross( a2,b2 )

        # @ # if abs(self.longvecnorm - self.nndisteq) > 0.0000001:
        # @ #     print "a2:",a2
        # @ #     print "b2:",b2
        # @ #     print "c2:",c2
        # @ # #R2 = np.array(a2,b2,c2)
        # @ # R2 = np.eye(3)
        # @ # R2[0,:] = a2
        # @ # R2[1,:] = b2
        # @ # R2[2,:] = c2
        # @ # print utils.printblue(np.dot(R2, self.longvec))

        # @ # self.longvec_orig = np.dot( R.T, np.dot(R2, self.longvec) )
        # @ # self.vec0_orig = np.dot( R.T, np.dot(R2, self.vec0) )
        # @ # #f_o = np.dot( R.T, np.dot(R2, tvec_curr) )
        # @ # if np.linalg.norm(self.longvec) != np.linalg.norm(self.vec0):
        # @ #     if abs(self.longvecnorm - self.nndisteq) > 0.0000001:
        # @ #         print utils.printred("longvec_curr:",self.longvec,"in orgin:",self.longvec_orig)
        # @ # else:
        # @ #     if abs(self.longvecnorm - self.nndisteq) > 0.0000001:
        # @ #         print "longvec_curr:",self.longvec,"in orgin:",self.longvec_orig
        # @ # if abs(self.longvecnorm - self.nndisteq) > 0.0000001:
        # @ #     print "vec0_curr:",self.vec0,"in origin:",self.vec0_orig
        # @ #     print ""
        # @ # if self.vec0_orig[2] != 0.0:
        # @ #     sys.exit("self.vec0_orig[2] != 0.0")
        # @ # if abs(self.vec0_orig[1] - 1.995) > 0.05:
        # @ #     print self.vec0_orig[1],self.vec0_orig[1]-1.955
        # @ #     sys.exit("self.vec0_orig[1] != 1.995")
        # @ # if abs(self.vec0_orig[0] - 1.995) > 0.05:
        # @ #     print self.vec0_orig[0]
        # @ #     sys.exit("self.vec0_orig[0] != 1.995")

        # @ # if self.indi == 0 and self.indj == 24:
        # @ #     sys.exit()
        if self.vec0[0] == 0.0 and self.shell == 2: # vec in y-z ebene
            self.topdirectionx = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
            if self.vec0[1] == 0.0 and self.shell == 2: # vec in y-z ebene
                self.topdirectiony = np.array([0.0,1.0,0.0])  # minus xdir  will also be correct
            if self.vec0[2] == 0.0 and self.shell == 2: # vec in y-z ebene
                self.topdirectiony = np.array([0.0,0.0,1.0])  # minus xdir  will also be correct
            self.topdirectionxnorm = 1.0
            self.topdirectionynorm = 1.0
            self.topdirectionznorm = 1.0
            self.tipdirectionxnorm = 1.0
            self.tipdirectionynorm = 1.0
            self.tipdirectionznorm = 1.0

        if self.vec0[1] == 0.0 and self.shell == 2: # vec in x-z ebene
            self.topdirectionx = np.array([0.0,1.0,0.0])  # minus xdir  will also be correct
            if self.vec0[0] == 0.0 and self.shell == 2: # vec in y-z ebene
                self.topdirectiony = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
            if self.vec0[2] == 0.0 and self.shell == 2: # vec in y-z ebene
                self.topdirectiony = np.array([0.0,0.0,1.0])  # minus xdir  will also be correct
            self.topdirectionxnorm = 1.0
            self.topdirectionynorm = 1.0
            self.topdirectionznorm = 1.0
            self.tipdirectionxnorm = 1.0
            self.tipdirectionynorm = 1.0
            self.tipdirectionznorm = 1.0

        if self.vec0[2] == 0.0 and self.shell == 2: # vec in x-y ebene
            self.topdirectionx = np.array([0.0,0.0,1.0])  # minus xdir  will also be correct
            if self.vec0[0] == 0.0 and self.shell == 2: # vec in y-z ebene
                self.topdirectiony = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
            if self.vec0[1] == 0.0 and self.shell == 2: # vec in y-z ebene
                self.topdirectiony = np.array([0.0,1.0,0.0])  # minus xdir  will also be correct
            self.topdirectionxnorm = 1.0
            self.topdirectionynorm = 1.0
            self.topdirectionznorm = 1.0
            self.tipdirectionxnorm = 1.0
            self.tipdirectionynorm = 1.0
            self.tipdirectionznorm = 1.0

        if self.shell == 1:
            if self.topdirectionxnorm == 0.0 or self.topdirectionynorm == 0.0 or self.topdirectionznorm == 0.0 or \
                self.tipdirectionxnorm == 0.0 or self.tipdirectionynorm == 0.0: # or np.linalg.norm(tipdirectionz) == 0.0:
                #print "nn:",nn,"indi:",indi,"indj:",indj,"posj:",posj,"|| first ||"
                print "topdirectionx:",self.topdirectionx
                print "topdirectiony:",self.topdirectiony
                print "topdirectionz:",self.topdirectionz
                print "tipdirectionx:",self.tipdirectionx
                print "tipdirectiony:",self.tipdirectiony
                print "tipdirectionz:",self.tipdirectionz
                sys.exit("t{o,i}pdirection{x,y,z} not found")

        if self.shell == 2:
            if self.topdirectionxnorm == 0.0 or self.topdirectionynorm == 0.0: # or \
                #np.linalg.norm(tipdirectionx) == 0.0 or np.linalg.norm(tipdirectiony) == 0.0: # or np.linalg.norm(tipdirectionz) == 0.0:
                #print "nn:",nn,"indi:",indi,"indj:",indj,"posj:",posj,"|| second ||"
                print "topdirectionx:",self.topdirectionx
                print "topdirectiony:",self.topdirectiony
                print "topdirectionz:",self.topdirectionz
                print "tipdirectionx:",self.tipdirectionx
                print "tipdirectiony:",self.tipdirectiony
                print "tipdirectionz:",self.tipdirectionz
                sys.exit("t{o,i}pdirection{x,y,z} not found")


        ##################################################################################
        # checks
        ##################################################################################
        if self.vec0[0] == 0.0 and self.shell > 2:
            sys.exit("not yet defined 33")

        if  type(self.u_topx) != bool or \
            type(self.u_topy) != bool or \
            type(self.u_topz) != bool:
                if self.topdirectionxnorm == 0.0 or self.topdirectionynorm == 0.0 or self.topdirectionynorm == 0.0:
                    sys.exit("not yet defined to{xyz}")

        if  type(self.u_tipx) != bool or \
            type(self.u_tipy) != bool or \
            type(self.u_tipz) != bool:
                if self.tipdirectionxnorm == 0.0 or self.tipdirectionynorm == 0.0 or self.tipdirectionynorm == 0.0:
                    sys.exit("not yet defined ti{xyz}")
        return

    def get_to_vec(self):
        ''' this we should only need if we do have transversal forces included '''
        if  type(self.u_topx) == bool and \
            type(self.u_topy) == bool and \
            type(self.u_topz) == bool and \
            type(self.u_tipx) == bool and \
            type(self.u_tipy) == bool and \
            type(self.u_tipz) == bool:
                return
        self.tvec = self.longvec - self.vec0
        self.tvecoutdir = utils.project_vector(self.tvec,self.topdirectionx)   # [ 0.3, 0.0, 0.0 ]
        self.rejection = utils.reject_vector(self.longvec, self.vec0)
        self.projection = utils.project_vector(self.longvec, self.vec0)
        self.tvecout = utils.project_vector(self.rejection,self.topdirectionx)   # hier ist es egal ob man tvecoutdir oder tvecoutdir*2 schreibet, nur die richtung ist von interesse
        self.tvecoutnorm = np.linalg.norm(self.tvecout)
        self.tipdirection = np.cross(self.vec0,self.topdirectionx)
        return

    def rotate_vec0_to_vec0_original_inplane_and_getforce_current(self, vec0_curr,longvec,vec0_orig_formapping,mappedvec_orig_dir,forcefuncx,forcefuncy,forcefuncz,
            tvecwirkdir = False, atomlist = False, atomi = False, atomj = False,bbbtext = False):
        '''
        forcefunz will come later when this function is extended to 3D (probably the rotation can alredy handle this)
        tvec will have to be rejected on mappedvec_orig_dir, then forcesfunc{x,y,z} can be applied

        @param {vec0,longvec}_{curr,orig}   : The current and original ve0 (equilibrium vector) and longvec (current distrubed vector) (correct length of vec0 is also important)
        @type  {vec0,longvec}_{curr,orig}   : numpy array, len 3
        @vec0_orig_formapping               : for to{x,y,z}: np.array([0.0, -2.065, -2.065])     for ti{x,y,z}: np.array([2.065,2.065,0]) (always looking from the considered
                                              atom to the 0-Atom and keeping the coordinate system as it is
        @mappedvecdir                       : is the direction in which the original forces were mapped e.g. [1.0, 0.0, 0.0] or [ 1.0, 1.0, 0.0 ]
                                              make sure mappedvec_orig_dir shows toward the posivive direction and not! [ -1.0, 0.0, 0.0 ] the sign is important to show the rotation
                                              which was / is the positive direction
        @mappedvecdir                       : numpy array, len 3

        @tvecwirkdir                        : the direction in which the force should act (for to{x,y,z} this will always be the direction out of the plane, for ti{x,y,z} this will
                                              be a direction (senkrecht/prallel) in the plane
                                              or a direction in x or y direction (in the plane)

                                               - a projection does not care about the direction it is projected on: So we only get the "correct" sign if the vector as such has the "correct" direction.

                                                 [209]utils.project_vector(np.array([ 1. , 0.,  0.]),np.array([ 0.402534  ,0.,        0.402534]))
                                                 Out[209]: array([ 0.5,  0. ,  0.5])

                                                 19:18:35 glensk@mac /Users/glensk/Dropbox/Understand_distributions/Pt/2x2x2sc_quer_3x3x3kp
                                                 [210]utils.project_vector(np.array([ 1. , 0.,  0.]),np.array([ -0.402534  ,0.,        -0.402534]))
                                                 Out[210]: array([ 0.5,  0. ,  0.5])

                                                 [211]utils.project_vector(np.array([ -1. , 0.,  0.]),np.array([ -0.402534  ,0.,        -0.402534]))
                                                 Out[211]: array([-0.5,  0. , -0.5])

                                                 19:22:31 glensk@mac /Users/glensk/Dropbox/Understand_distributions/Pt/2x2x2sc_quer_3x3x3kp
                                                 [212]utils.project_vector(np.array([ -1. , 0.,  0.]),np.array([ 0.402534  ,0.,        0.402534]))
                                                 Out[212]: array([-0.5,  0. , -0.5])

                                                 solange tvecwirkdir nur eine richtung hat [1.0, 0.0, 0.0] und wie bei tox die Kraftfunktion symmetrisch ist, wird durch dass zurueckdrehen immer
                                                 das richtige vorzeichen erlangt. (zumindest glaube ich das momentan)
                                                # tvec_o will always be positive: for tox the rotation decides then if the real sign stays positive or negative
                                                                                    for tix we have either have (pos,pos,0) or (neg,neg,0) nothing else for tvec_o

                                                bei tix haut dass noch nicht hin! tvec_o im folgenden fall ist falsch:
                                                    atomi:16 atomj:0
                                                    bbbtext:TI{x,y,z}
                                                    ('TI{x,y,z} atomi:16 atomj:0 tvec_curr:', array([ 0.659744,  0.      , -0.659744]), 'tvec_o:', array([ 0.659744,  0.659744,  0.      ]))
                                                    ('tvec_orig1;', array([ 0.219915,  0.219915, -0.879659]), 'tvec_orig:', array([ 0.29322 ,  0.806354, -0.366525]))
                                                      vec0_curr: [-2.065  0.    -2.065]
                                                      vec0_orig_formapping: [-2.065  2.065  0.   ]
                                                      longvec: [-1.065  0.    -2.065] 2.32345647689
                                                      tvecfull: [ 1.  0.  0.]
                                                      tvecwirkdir: [ 0.659744  0.       -0.659744]
                                                      tvec_curr_dir: [ 0.5  0.  -0.5] 16 || [20, 16, 26]
                                                      tvec_curr: [ 0.659744  0.       -0.659744]     (tvec_curr = utils.project_vector(tvecfull,tvecwirkdir))
                                                      longvec_curr: [-1.405256  0.       -2.724744]
                                                      tvec_o = np.dot( R.T, np.dot(R2, tvec_curr) )
                                                      R2: [[ 0.707107  0.       -0.707107]
                                                     [-0.707107  0.       -0.707107]
                                                     [ 0.        1.        0.      ]]
                                                      np.dot(R2, tvec_curr): [ 0.933019  0.        0.      ]
                                                      np.dot( R.T, np.dot(R2, tvec_curr) ): [ 0.659744  0.659744  0.      ]
                                                      tvec_o: [ 0.659744  0.659744  0.      ]
                                                      mappedvec_orig_dir: [ 1.  1.  0.]
                                                      perpendicular_fromvec_tovec: [ 0.  0.  0.]
                                                      forcefuncx: [0.0, 0.0, 0.0006989016595, -0.006947637395, -0.008869118256, 0.0021385945798, 0.0013974540059, -3.2447587e-05]
                                                      forcefuncy: [0.0, 0.0, 0.0001951898033, 0.0155931071382, 0.0117455886566, 0.0016183055654, -0.001648220302, -0.002814611506, -0.00090722266, 0.0007605163137, 0.000362115408]
                                                      forcefuncz: False
                                                    fx_,fy_,fz_: -0.031772789092 0.0667925403173 0.0
                                                    fxorig,fyorig,fzorig: [-0.031773  0.        0.      ] [ 0.        0.066793  0.      ] [ 0.  0.  0.]
                                                    fx    ,fy    ,fz    : [-0.031773  0.       -0.      ] [ 0.        0.       -0.066793] [ 0.  0.  0.]




        '''
        ex = ey = ez = 0.0
        fxorig = np.array([0.0,0.0,0.0])
        fyorig = np.array([0.0,0.0,0.0])
        fzorig = np.array([0.0,0.0,0.0])
        fx     = np.array([0.0,0.0,0.0])
        fy     = np.array([0.0,0.0,0.0])
        fz     = np.array([0.0,0.0,0.0])
        fx_    = np.array([0.0,0.0,0.0])
        fy_    = np.array([0.0,0.0,0.0])
        fz_    = np.array([0.0,0.0,0.0])
        tvec_mapped = np.array([0.0,0.0,0.0])
        longvec_orig = np.array([0.0,0.0,0.0])
        tvec_orig = np.array([0.0,0.0,0.0])

        if type(forcefuncx) == bool and type(forcefuncy) == bool and type(forcefuncz) == bool:
            return ex,ey,ez,fx,fy,fz


        if type(tvecwirkdir)  != np.ndarray:
            print "tvecwirkdir:",tvecwirkdir
            print "type(tvecwirkdir):",type(tvecwirkdir)
            sys.exit("tvecwirkdir is not a np.ndarray")



        ######################################################################
        # tvecfull, longvec -> tvec_curr and longvec_curr
        ######################################################################
        tvecfull  = longvec - vec0_curr

        # von tvec brauchen wir nur den teil der in wirkdir wirkt :) (==tvec_curr, wobei curr fuer current considered part steht)  (je nachdem ob wir gerade to oder tis tip anschauen)
        # in der ebene ist dies ganz einfach: siehe thoughts.pptx; wenn wir aus der ebene rausgehen, ...., dann haben wir wieder eine Winkelabhaengigkeit!
        # bei ti{x,y,z} vermeiden wir dies (vorlaeufig) indem wir einfach auf die senkrechte richtung rejizieren
        # alles was bei ti{x,y,z} aus der ebene geht , geht in tox teile
        # bei to{x,y,z} .... muessen wir das genau so machen, alles, also tvec, sollte aufgespalten werden, am besten von vorneherin in to und ti anteile!
        tvec_curr_dir = utils.project_vector(tvecfull,tvecwirkdir)
        rejection = utils.reject_vector(longvec, vec0_curr)  # dies ist so nur ok wenn man die rejektion auf quer wuenscht (fuer tox funktioniert es auch, siehe pptx), jedoch wollen wir etwas anderes
        tvec_curr = utils.project_vector(rejection,tvecwirkdir)
        tvec_curr = utils.project_vector(tvecfull,tvecwirkdir)

        longvec_curr =  tvec_curr + vec0_curr



        if np.linalg.norm(tvec_curr) == 0.0:
            return ex,ey,ez,fx,fy,fz
            #return fxorig,fyorig,fzorig,fx,fy,fz,longvec_orig,vec0_orig_formapping,tvec_mapped,tvec_orig,fx_,fy_,fz_

        # TODO: spaeter auch mal alles anstatt mit rejection / alles einfach mal mit projections versuchen!
        #print ""
        #print "tvec_curr_dir:",tvec_curr_dir,np.linalg.norm(tvec_curr_dir)
        #print "tvec_curr", tvec_curr,np.linalg.norm(tvec_curr)
        #print "longvec_curr",longvec_curr
        import warnings
        warnings.filterwarnings('error')
        try:
            tvec_curr_dir = utils.project_vector(tvecfull,tvecwirkdir)
            tvec_curr = utils.project_vector(rejection,tvecwirkdir)
            tvec_curr = utils.project_vector(tvecfull,tvecwirkdir)
        except RuntimeWarning:
            print "  XXtvec_curr_dir:",tvec_curr_dir
            print "  XXtvec_curr", tvec_curr
            print "  XXlongvec_curr",longvec_curr
            sys.exit()


        fromvec     = vec0_curr
        tovec       = vec0_orig_formapping

        # TODO (much later):
        # define out of plane and inplane direction, this we can not do internally in this module but have to define externally ant then give to this module only the tvec in the
        # correct inplane / out of plane part. Practically what "outofplane" means will depend on the CURRENT position of the respective atoms and not on the undisplaced vec0 positions.
        # We could consider this later on!




        #Rotmat      = utils.R_2vect(vec0_curr,vec0_orig_formapping)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS
        #Rotmatback  = utils.R_2vect(vec0_orig_formapping,vec0_curr)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS
        #                                             # if mappedvec_orig_dir is perpendicular to the plane spaned by fromvec(==vec0_curr),tovec -> no problem
        #                                             # otherwise: further checking necessary
        #                                             # in other words: we have to make sure that tvec points in same direction as originally!
        #Rotmatt      = utils.R_2vect(tvec_curr,mappedvec_orig_dir)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS
        #Rotmattback  = utils.R_2vect(mappedvec_orig_dir,tvec_curr)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS


        ## check if mappedvec_orig_dir is parallel to the perpendicular vector spanned by {from,to}vec


        ## after applying rotation on ... there should be a vector like (x,0,0) or (0,y,0) or (0,z,0)
        #tvec_orig1 = np.dot(Rotmat,tvec_curr)
        #longvec_orig1 = np.dot(Rotmat,longvec_curr)

        ## with Rotmat2 we need to make sure that it only rotates around vec0_orig_formapping; what happend previously is taat in order to ensure that tvec points in the correct direction
        ## vor the -1/-1 direction Rotvec was simply movec back (negative einheitsmatrix) this can only happen if direction vec0 is in the same direction as vec0_orig_formapping and both
        ## go in exactly opposite directions. To avoid this one needs to make sure that Rotmat2 only rotates around vec0_orig_formapping. Typically Rotmat2 will only be a rotation around
        ## vec0_orig_formapping by 180 degrees or none.
        ## Dies hier funktionier noch nicht 100% ig

        ##Rotmat2      = utils.R_2vect(tvec_orig1,mappedvec_orig_dir,fixed_rotation_axis=vec0_orig_formapping)  # the length of the vectors does not play a role!
        ##Rotmatback2  = utils.R_2vect(mappedvec_orig_dir,tvec_orig1,fixed_rotation_axis=vec0_orig_formapping)

        #Rotmat2      = utils.R_2vect(tvec_curr,mappedvec_orig_dir)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS
        #Rotmatback2  = utils.R_2vect(mappedvec_orig_dir,tvec_curr)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS

        #tvec_orig = np.dot(Rotmat2,tvec_orig1)
        #longvec_orig = np.dot(Rotmat2,longvec_orig1)

        ####################################################################################################################
        # bei der rotation back muessen wir auch sicherstellen dass der tvec_curr wieder hergestellt wird
        # also brauchen wir auch da wieder 2 rotationen zurueck!!!!!
        ####################################################################################################################





        # TODO: What we want in principle is to map the Force as a function of x-y (and z) position of an atom. Since the next nearest neighbors are so important we will get with
        # this a) the angular dependence and b) most probably the forces in arbitrary position correct.

        #######################################################################################################################################################
        # wenn hier tvec_curr senkrecht auf vec0_curr steht (sollte es) sollte es auch mit matrizenmultiplikation gehen (vorgesschlagen von gerard)
        # was machen wir wenn tvec_curr nicht senkrecht auf vec0_curr steht?
        # R == refernez!
        # R is the matrix which transposes a vector from the x,y,z refernce frame onto the reference
        a = orig_xvec = mappedvec_orig_dir/np.linalg.norm(mappedvec_orig_dir)
        b = orig_vec0 = vec0_orig_formapping/np.linalg.norm(vec0_orig_formapping)
        c = np.cross( orig_xvec,orig_vec0 )
        #print "a:",a,b,c
        R = np.eye(3)
        R[0,:] = a
        R[1,:] = b
        R[2,:] = c


        # R2 is the vecor we currently look at
        # R2 is the matrix which transposes the current vector (or current coordinate system) a2,b2,c2 into x,y,z reference frame
        a2 = tvec_curr/np.linalg.norm(tvec_curr)
        b2 = vec0_curr/np.linalg.norm(vec0_curr)
        c2 = np.cross( a2,b2 )
        #R2 = np.array(a2,b2,c2)
        R2 = np.eye(3)
        R2[0,:] = a2
        R2[1,:] = b2
        R2[2,:] = c2
        #print R2
        #print np.dot(R2, tvec_curr)
        tvec_o = np.dot( R.T, np.dot(R2, tvec_curr) )
        f_o    = np.dot( R.T, np.dot(R2, tvec_curr) )
        #######################################################################################################################################################
        # check if tvec_o is parallel to mappedvec_orig_dir
        perpendicular_fromvec_tovec = np.cross(tvec_o,mappedvec_orig_dir)


        ###############################################################################
        # old
        ###############################################################################
        # anwenden der forcefuncx, forcefuncy
        # now we want to reject the tvec_orig onto the originally mapped vector
        #rejection_orig = utils.reject_vector(longvec_orig, vec0_orig_formapping)
        #tvec_mapped = utils.project_vector(rejection_orig,mappedvec_orig_dir)   # hier ist es egal ob man tvecoutdir oder tvecoutdir*2 schreibet, nur die richtung ist von interesse
        #ex, fx_ = getef(math.copysign(np.linalg.norm(tvec_mapped),tvec_mapped[0]),forcefuncx,pot=False)
        #fxorig = np.array([fx_,0.0,0.0])
        #ey, fy_ = getef(math.copysign(np.linalg.norm(tvec_mapped),tvec_mapped[1]),forcefuncy,pot=False)
        #fyorig = np.array([0.0,fy_,0.0])
        #ez, fz_ = getef(math.copysign(np.linalg.norm(tvec_mapped),tvec_mapped[2]),forcefuncz,pot=False)
        #fzorig = np.array([0.0,0.0,fz_])

        # tvec_o will always be positive: in case of tix
        ex, fx_ = self.getef(math.copysign(np.linalg.norm(tvec_o),tvec_o[0]),forcefuncx,pot='poly')
        fxorig = np.array([fx_,0.0,0.0])
        ey, fy_ = self.getef(math.copysign(np.linalg.norm(tvec_o),tvec_o[1]),forcefuncy,pot='poly')
        fyorig = np.array([0.0,fy_,0.0])
        ez, fz_ = self.getef(math.copysign(np.linalg.norm(tvec_o),tvec_o[2]),forcefuncz,pot='poly')
        fzorig = np.array([0.0,0.0,fz_])

        fx = np.dot( R2.T, np.dot(R, fxorig) )
        fy = np.dot( R2.T, np.dot(R, fyorig) )
        fz = np.dot( R2.T, np.dot(R, fzorig) )

        def print_vecs_to_screen():
            print utils.printred("atomi:"+str(atomi)+" atomj:"+str(atomj))
            print utils.printred("bbbtext:"+str(bbbtext))
            print utils.printyellow(str(bbbtext)+" atomi:"+str(atomi)+" atomj:"+str(atomj)+" tvec_curr:",tvec_curr,"tvec_o:",tvec_o)
            print utils.printyellow("tvec_orig1;",tvec_orig1,"tvec_orig:",tvec_orig) #,fx_,fy_,fz_)
            print "  vec0_curr:",vec0_curr
            print "  vec0_orig_formapping:",vec0_orig_formapping
            print "  longvec:",longvec,np.linalg.norm(longvec)
            print "  tvecwirkdir:",tvecwirkdir
            print ""
            print "  tvecfull  = longvec - vec0_curr:",tvecfull
            print "  tvec_curr_dir = utils.project_vector(tvecfull,tvecwirkdir):",tvec_curr_dir
            print "  rejection = utils.reject_vector(longvec, vec0_curr):",rejection
            print "  tvec_curr = utils.project_vector(rejection,tvecwirkdir):",tvec_curr
            print "  longvec_curr:",longvec_curr
            print "  tvec_o = np.dot( R.T, np.dot(R2, tvec_curr) ):",tvec_o
            print " "
            print "  R2:",R2
            print "  np.dot(R2, tvec_curr):",np.dot(R2, tvec_curr)
            print "  np.dot( R.T, np.dot(R2, tvec_curr) ):",np.dot( R.T, np.dot(R2, tvec_curr) )
            print "  tvec_o:",tvec_o
            print "  mappedvec_orig_dir:",mappedvec_orig_dir
            print "  perpendicular_fromvec_tovec:",perpendicular_fromvec_tovec
            print "  forcefuncx:",forcefuncx
            print "  forcefuncy:",forcefuncy
            print "  forcefuncz:",forcefuncz
            print "fx_,fy_,fz_:",fx_,fy_,fz_
            print "fxorig,fyorig,fzorig:",fxorig,fyorig,fzorig
            print "fx    ,fy    ,fz    :",fx    ,fy    ,fz
            print ""
            if abs(np.linalg.norm(perpendicular_fromvec_tovec)) >= 1e-9:
               print utils.printred("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
               print "perpendicular_fromvec_tovec = np.cross(tvec_o,mappedvec_orig_dir)",perpendicular_fromvec_tovec
               print "np.linalg.norm(perpendicular_fromvec_tovec)",np.linalg.norm(perpendicular_fromvec_tovec)
               print "tvec_o:",tvec_o
               print "mappedvec_orig_dir:",mappedvec_orig_dir
               print "tvec_orig1:",tvec_orig1," == np.dot(Rotmat,tvec_curr)","  tvec_curr:",tvec_curr
               print utils.printred("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
               sys.exit("not 0:")

        # second check in case tvec_o had to be gained by rotation
        # is currently not used
        perpendicular_fromvec_tovec = np.cross(tvec_o,mappedvec_orig_dir)

        if type(atomlist) != bool and (atomi in atomlist or atomj in atomlist):
            print_vecs_to_screen()
        return ex,ey,ez,fx,fy,fz

    def get_to_ti_energy_forces(self, contr = False):
        if type(contr) == False:
            sys.exit("contr not in [ \'to\', \'ti' ]")
        if contr not in [ 'to', 'ti' ]:
            sys.exit("contr not in [ \'to\', \'ti' ]")
        self.ftopx = np.array([0.0,0.0,0.0])
        self.ftopy = np.array([0.0,0.0,0.0])
        self.ftopz = np.array([0.0,0.0,0.0])
        self.etopx = np.array([0.0])
        self.etopy = np.array([0.0])
        self.etopz = np.array([0.0])

        self.ftipx = np.array([0.0,0.0,0.0])
        self.ftipy = np.array([0.0,0.0,0.0])
        self.ftipz = np.array([0.0,0.0,0.0])
        self.etipx = np.array([0.0])
        self.etipy = np.array([0.0])
        self.etipz = np.array([0.0])

        if  type(eval("self.u_"+contr+"px")) == bool and \
            type(eval("self.u_"+contr+"py")) == bool and \
            type(eval("self.u_"+contr+"pz")) == bool:
                return

        #if  type(self.u_topx) == bool and \
        #    type(self.u_topy) == bool and \
        #    type(self.u_topz) == bool and \
        #    type(self.u_tipx) == bool and \
        #    type(self.u_tipy) == bool and \
        #    type(self.u_tipz) == bool:
        #        return
        #sys.exit('das ist der teil fuer TOX; hier sollte ich nicht hin!')

        divide = False
        if self.numberofatoms == 32:
            divide = 2.
        if self.numberofatoms == 108:
            divide = 3.
        if type(divide) == bool:
            sys.exit("sc size unknown")

        cellvecdirnorm = -(self.crystal0.cellvec[0][0]/divide)/np.sqrt(2)
        if self.shell == 1:
            #vec0_orig_formapping = np.array([0.0, -2.065, -2.065])
            cellvecdirnorm = -(self.crystal0.cellvec[0][0]/divide)/np.sqrt(2)
            if contr == 'to':
                self.vec0_orig_formapping = np.array([0.0, cellvecdirnorm, cellvecdirnorm])   # verktor zu atom 8
                self.mappedvec_orig_dir = np.array([1.0, 0.0, 0.0]) # hier entlang haben wir ausgelenkt um die kraft auf atom 8 (in x richtung) zu mappen
            if contr == 'ti':
                self.vec0_orig_formapping = np.array([-cellvecdirnorm, cellvecdirnorm, 0.0]) # vektor zu atom
                self.mappedvec_orig_dir = np.array([1.0, 1.0, 0.0])   # hier entlang haben wir ausgelenkt um die kraft auf atom 60 bzw 29 bzw 27 zu mappen
        if self.shell == 2:
            #vec0_orig_formapping = np.array([0.0, -4.13, 0.0])
            self.vec0_orig_formapping = np.array([0.0, -self.crystal0.cellvec[0][0]/divide, 0.0])
            self.mappedvec_orig_dir = np.array([1.0, 0.0, 0.0])
            if contr == 'ti':
                sys.exit('not yet defined')
        if self.shell > 2:
            sys.exit("not yet 34")

        #atomlist = [8]
        #atomlist = False # only for output to screen
        #if contr == 'to':
        #    self.etopx,self.etopy,self.etopz,self.ftopx,self.ftopy,self.ftopz = \
        #        self.rotate_vec0_to_vec0_original_inplane_and_getforce_current(\
        #        self.vec0,\
        #        self.longvec,\
        #        self.vec0_orig_formapping,\
        #        self.mappedvec_orig_dir,\
        #        self.u_topx,\
        #        self.u_topy,\
        #        self.u_topz,\
        #        tvecwirkdir = self.topdirectionx,\
        #        atomlist=atomlist, \
        #        atomi=self.indi,\
        #        atomj=self.indj, \
        #        bbbtext="TO{X,Y,Z}")


        def project_vec_into_original_cell(\
                # current (long)vec
                vec,\
                # original reference frame
                vec0_orig,vec0_orig_outdirection,\
                # current reference frame
                vec0_curr,vec0_curr_outdirection,\
                forcefunc,\
                        ):
            '''
            vec                     : is the current longvec you want to map into the original cell
            original reference frame: vec0_orig and vec0_orig_outdirection have to be orthogonal to each other and live inside the original cell
            current reference frame :


            '''
            #####################################
            # span the original reference frame
            #####################################
            a = orig_xvec = vec0_orig_outdirection/np.linalg.norm(vec0_orig_outdirection)
            b = orig_vec0 = vec0_orig/np.linalg.norm(vec0_orig)
            c = np.cross( orig_xvec,orig_vec0 )
            #print "a:",a,b,c
            R = np.eye(3)
            R[0,:] = a
            R[1,:] = b
            R[2,:] = c

            #####################################
            # span the original reference frame
            #####################################
            # R2 is the vecor we currently look at
            # R2 is the matrix which transposes the current vector (or current coordinate system) a2,b2,c2 into x,y,z reference frame
            a2 = vec0_curr_outdirection/np.linalg.norm(vec0_curr_outdirection)
            b2 = vec0_curr/np.linalg.norm(vec0_curr)
            c2 = np.cross( a2,b2 )
            R2 = np.eye(3)
            R2[0,:] = a2
            R2[1,:] = b2
            R2[2,:] = c2
            vec_orig = np.dot( R.T, np.dot(R2, vec) )   # vec == longvec



            #####################################
            # get forces in original cell
            #####################################
            ecurr, fnorm = self.getef(math.copysign(np.linalg.norm(vec_orig),vec_orig[0]),forcefunc,pot='poly')
            # normalize vec0_orig_outdirection === wirdirection der Kraft
            vec0_orig_outdirection_normalized = vec0_orig_outdirection/np.linalg.norm(vec0_orig_outdirection)
            forig = vec0_orig_outdirection_normalized * fnorm
            #forig = np.array([0.0,0.0,fnorm])  # zeigt hier entlang tvecwirkdir

            #####################################
            # map back into current cell
            #####################################
            forig = np.dot( R2.T, np.dot(R, forig) )
            return vec_orig, forig, ecurr




        cellvecdirnorm = -(self.crystal0.cellvec[0][0]/divide)/np.sqrt(2)
        # tox
        self.vec0_orig_formapping = np.array([0.0, cellvecdirnorm, cellvecdirnorm])   # verktor zu atom 8
        self.mappedvec_orig_dir = np.array([1.0, 0.0, 0.0]) # hier entlang haben wir ausgelenkt um die kraft auf atom 8 (in x richtung) zu mappen

        # tox Funktioniert, richtung und betrag
        vec0_orig = np.array([0.0, self.crystal0.cellvec[0][0]/divide/2, self.crystal0.cellvec[0][0]/divide/2])
        vec0_orig_outdirection = np.array([1.0, 0.0, 0.0])


        ######################### tox in andere richtungen zu definieren funktioniert nicht ! Die Richtung wird falsch!
        # @ # tox does not work: betrag und richgung ok, vorzeichen falsch
        # @ vec0_orig = np.array([self.crystal0.cellvec[0][0]/divide/2, self.crystal0.cellvec[0][0]/divide/2,0.0])
        # @ vec0_orig_outdirection = np.array([0.0, 0.0, 1.0])

        # @ # tox does not work: betrag und richgung ok, vorzeichen falsch
        # @ vec0_orig = np.array([-self.crystal0.cellvec[0][0]/divide/2, -self.crystal0.cellvec[0][0]/divide/2,0.0])
        # @ vec0_orig_outdirection = np.array([0.0, 0.0,1.0])



        vec0_curr = self.vec0
        self.vec0_curr = self.vec0



        #if contr == 'to':
        #    self.vec_orig_cell, self.forces_to, ene_to = project_vec_into_original_cell(\
        #        #vec = utils.project_vector(self.longvec - vec0_curr,self.topdirectionx),
        #        vec = utils.project_vector(self.longvec,self.topdirectionx),
        #        vec0_orig = vec0_orig,
        #        vec0_orig_outdirection = vec0_orig_outdirection,
        #        vec0_curr = vec0_curr,
        #        vec0_curr_outdirection = self.topdirectionx,
        #        forcefunc = self.u_topx)
        #    self.ftopx = self.forces_to
        #    self.etopx = ene_to

        def get_new_local_basis(longvec,longvec0,topdirection):
            ''' longvec0 und topdirection spannnen die basis auf
                longvec ist der atuelle longvec
            '''
            l  = longvec
            l0 = longvec0
            t1 = topdirection
            t2 = np.cross(l0,t1)

            t11_unnorm = np.cross(l,t2)
            t22_unnorm = np.cross(l,t11_unnorm)

            # die * -1. braucht man fuer die to sachen damit die vorzeichen der vektoren stimmen
            if np.linalg.norm(t11_unnorm) == 0.0:
                t11 = np.array([0.0,0.0,0.0])
            else:
                t11 = t11_unnorm/np.linalg.norm(t11_unnorm) * -1.

            if np.linalg.norm(t22_unnorm) == 0.0:
                t22 = np.array([0.0,0.0,0.0])
            else:
                t22 = t22_unnorm/np.linalg.norm(t22_unnorm) * -1.

            return t11,t22



        #l = self.longvec

        #l0 = self.vec0_curr
        #t1 = self.topdirectionx
        #t2 = np.cross(l0,t1)

        #t11_unnorm = np.cross(l,t2)
        #t22_unnorm = np.cross(l,t11_unnorm)

        #t11 = t11_unnorm/np.linalg.norm(t11_unnorm)
        #t22 = t22_unnorm/np.linalg.norm(t22_unnorm)

        #t11, t22 = get_new_local_basis(l = self.longvec, longvec0 = self.vec0_curr, topdirectionx = self.topdirectionx)

        #print "topdirx:",t1,"topnew:",t11

        if contr == 'to':
            vec_to = utils.project_vector(self.longvec,self.topdirectionx) # hat aber die richtige richtung
            #print "vec_ti_:",vec_ti_
            vec_to_length = np.linalg.norm(vec_to)
            #print vec_to,vec_to_length
            #print "vec_ti_norm:",vec_ti_norm

            if vec_to_length == 0.0:
                vec_to_normalized_corr_direction = np.array([0.0,0.0,0.0])
            else:
                vec_to_normalized_corr_direction = vec_to/vec_to_length

            #print "vec_ti_normalized:",vec_ti_normalized
            ene_to, fnorm_to = self.getef(vec_to_length,self.u_topx,pot='poly')
            ene_to_a = ene_to
            forces_to = vec_to_normalized_corr_direction * fnorm_to

            #vec_to_ = utils.project_vector(vec_to,self.topdirectionx)

            t11_norm, t22_norm = get_new_local_basis(longvec = self.longvec, longvec0 = self.vec0_curr, topdirection = vec_to_normalized_corr_direction)


            # wenn ich es richtig sehe gibt es hier 2 Moeglichkeiten:
            # a) nehme die "alte" Kraft und "alte" energie und wende die an mit dem neuen tvec_vector -> energies bleiben gleich -> Kraefte (Phasenraum) aendert sich
            # b) pojeziere den alten transversalen vektor auf den neuen --> energien werden was kleiner, Kraefte auch
            # X) Zusaetzliche Frage: warum zeigen die to kraefte in die falsche richtung?

            # a)
            forces_to_orthbasis_a = t11_norm * fnorm_to

            # b)
            #print vec_to, t11_norm
            vec_orthbasis_b = utils.project_vector(vec_to,t11_norm)
            ene_to_b, fnorm_to_b = self.getef(np.linalg.norm(vec_orthbasis_b),self.u_topx,pot='poly')
            forces_to_orthbasis_b = t11_norm * fnorm_to_b
            #print vec_to, vec_orthbasis_b

            #print "topdir+richtung:",vec_to,vec_to_normalized_corr_direction,"topnew:",t11,"---> oldF:",forces_to,"newF:",forces_to_orthbasis_a


            #############################################################
            # check that the greatest value is on the same index
            #############################################################
            values = abs(forces_to)
            lvalues = list(values)
            old_i = lvalues.index(max(lvalues))

            values_n = abs(forces_to_orthbasis_a)
            lvalues_n = list(values_n)
            new_i = lvalues_n.index(max(lvalues_n))
            #print values,lvalues, old_i,new_i
            if old_i != new_i:
                sys.exit("greatest indizes (value positions) in old and new basis differ")

            #############################################################
            # check that the greatest value are of same sign
            #############################################################
            if abs(forces_to[old_i]) > 0.0001:
                if np.sign(forces_to[old_i]) != np.sign(forces_to_orthbasis_a[new_i]):
                    print "old:",forces_to
                    print "new:",forces_to_orthbasis_a
                    sys.exit('signes of old biggest value and new biggest value differ old:'+str(forces_to)+"  new:"+str(c_to_orthbasis_a))

            #############################################################
            # check how much energy deviates between a and b
            #############################################################
            aaa = round((ene_to_a-ene_to_b)/(ene_to_a)*100.,0)
            bbb = round((ene_to_a-ene_to_b)/(ene_to_b)*100.,0)

            if False:
                print "",self.longvec,vec_to_normalized_corr_direction,"topnew_norm:",t11_norm,"---> oldF:",forces_to,"newF_a:",forces_to_orthbasis_a,"_b:",forces_to_orthbasis_b,aaa,bbb


            #print "forces_ti:",forces_ti
            self.etopx = ene_to
            self.ftopx = forces_to

            self.etopx = ene_to_a
            self.ftopx = forces_to_orthbasis_a

            #self.etopx = ene_to_b
            #self.ftopx = forces_to_orthbasis_b


            #print "self.ftipx:",self.ftipx

        if contr == 'ti':
            vec_ti_ = utils.project_vector(self.longvec,self.tipdirection)
            #print "vec_ti_:",vec_ti_
            vec_ti_norm = np.linalg.norm(vec_ti_)
            #print "vec_ti_norm:",vec_ti_norm

            if vec_ti_norm == 0.0:
                vec_ti_normalized = np.array([0.0,0.0,0.0])
            else:
                vec_ti_normalized = vec_ti_/vec_ti_norm

            #print "vec_ti_normalized:",vec_ti_normalized
            ene_ti, fnorm_ti = self.getef(vec_ti_norm,self.u_tipx,pot='poly')
            forces_ti = vec_ti_normalized * fnorm_ti * -1.0
            #print "forces_ti:",forces_ti
            self.etipx = ene_ti
            self.ftipx = forces_ti
            #print "self.ftipx:",self.ftipx


            #ene_ti, fnorm_ti = self.getef(math.copysign(np.linalg.norm(vec_orig),vec_orig[0]),forcefunc,pot='poly')
            #self.tix_vec_orig_cell, self.forces_ti, ene_ti = project_vec_into_original_cell(\
            #        #vec = utils.project_vector(self.longvec - vec0_curr,self.topdirectionx),
            #        vec = utils.project_vector(self.longvec,self.tipdirection),
            #        vec0_orig = vec0_orig,
            #        vec0_orig_outdirection = vec0_orig_outdirection,
            #        vec0_curr = vec0_curr,
            #        vec0_curr_outdirection = self.topdirectionx,
            #        forcefunc = self.u_topx)
            # try instead of    vec = self.longvec:
            #                   vec = tvec_curr = utils.project_vector(tvecfull,tvecwirkdir)
            #                   with tvecfull  = longvec - vec0_curr
            #                   and  tvecwirkdir = self.topdirectionx
            #print "vec0_curr:",vec0_curr,"vec0_orig:",vec0_orig,"|| vec_curr:",self.longvec,"vec_orig:",self.vec_orig_cell,
            #print "vec0_curr:",vec0_curr,"tipdir:",self.tipdirection,"|| vec_curr:",self.longvec,"vec_ti_:",vec_ti_,"|||vec_orig:",self.vec_orig_cell,"|to:",self.forces_to,"|ti:",self.forces_ti
            #print "vec0_curr:",vec0_curr,"tipdir:",self.tipdirection,"|| vec_curr:",self.longvec,"|||vec_orig:",self.vec_orig_cell,"vec_ti_:",vec_ti_,"norm:",vec_ti_norm,"ene_ti:",ene_ti,"fnorm_ti:",fnorm_ti,"forces_ti:",forces_ti
            #print "|||vec_orig:",self.vec_orig_cell,"Fto:",self.forces_to,"ti_norm:",vec_ti_norm,"fnorm_ti:",fnorm_ti,"forces_ti:",forces_ti

            #self.longvec_all_mapped_to_first_quadrant = utils.append_row_to_2d_array(inarray = self.longvec_all_mapped_to_first_quadrant, addrow = self.vec_orig_cell)
            #if self.ftopx[0] != self.forces_to[0]:
            #    sys.exit("prob")
            #if self.ftopx[1] != self.forces_to[1]:
            #    sys.exit("prob")
            #if self.ftopx[2] != self.forces_to[2]:
            #    sys.exit("prob")


            #if vec0_curr[0] == 0.0:
            #    self.longvec[0]

        #print self.ftipx
        #if np.isnan(self.ftipx[0]) == True:
        #    print "isnan! "
        #    sys.exit()
        self.forcestrantopx[self.indi] = self.forcestrantopx[self.indi] + self.ftopx + self.ftipx
        self.energytrantopx[self.indi] = self.energytrantopx[self.indi] + self.etopx + self.etipx
        self.forcestrantopy[self.indi] = self.forcestrantopy[self.indi] + self.ftopy + self.ftipy
        self.energytrantopy[self.indi] = self.energytrantopy[self.indi] + self.etopy + self.etipy
        self.forcestrantopz[self.indi] = self.forcestrantopz[self.indi] + self.ftopz + self.ftipz
        self.energytrantopz[self.indi] = self.energytrantopz[self.indi] + self.etopz + self.etipz
        return

    def loop_over_atoms(self):
        ''' outer loop over atoms '''
        if self.verbose > 1:
            print "self.shells:",self.shells

        ##################################################################################
        # loop over all atoms in the original cell == self.positions (not the created supercell)
        # NNlist_all are the nnindizes of all coresponding cells
        ##################################################################################
        for self.indi,self.posi in enumerate(self.positions):
            NNlist_all, nndist_all = self.sc0.get_NNlist(
                        self.indi,
                        self.shells,
                        cell = self.sc0.cellvec,
                        coord_rrel = self.sc0neverchange.rrel,
                        return_result_and_NNdist = True)

            self.sc.center_atoms_around_atom(self.indi,coord_cart=self.sc.rcar,cell=self.sc.cellvec)
            self.longvecnorm_all_outerloop = np.linalg.norm(self.sc.rcar,axis=1)

            ##############################################################################
            # loop over all shells of interest
            ##############################################################################
            for self.shell in self.shells:

                ##########################################################################
                # get nearest neighbors of this atom
                # To make this faster ensure that we know the coord_cart and coord_rrel
                # we dont want to load those in every time we are in the inner loop
                ##########################################################################
                self.NNlist = NNlist_all[self.shell-1]
                self.nndisteq = nndist_all[self.shell]

                self.checks_for_fcc()
                self.get_or_load_parameters_lon_to_ti_depending_on_shell() # it depends on indj and nndist if we load poly or mc1 if lon2_pottype is defined!
                self.print_info_loop1()

                ##############################################################################
                # inner loop
                ##############################################################################
                #print "self.NNlist (inner loop):",self.NNlist
                for self.indj in self.NNlist:
                    self.posj = self.sc.rcar[self.indj]
                    self.longvecnorm = self.longvecnorm_all_outerloop[self.indj]

                    self.get_lon_vec()
                    #if abs(self.longvecnormcheck -self.longvecnorm ) >= 0.0000001:
                    #    sys.exit("diff in longvec")
                    self.get_lon_energy_forces()
                    self.print_info_loop2_lon_vec()
                    self.save_information_of_run_add(contr = 'lon')

                    self.get_to_ti_direction_vecs_fcc()
                    if self.u_toti == True:
                        self.get_to_ti_direction_vecs_fcc()
                        self.get_to_vec()
                        self.get_to_ti_energy_forces(contr = 'to')
                        self.get_to_ti_energy_forces(contr = 'ti')
                    self.save_information_of_run_add(contr = 'to')
        return

    def get_harmonic_energy_forces(self):
        h = hesse.read_Hessematrix("HesseMatrix_sphinx")
        fHAR,eHAR = hesse.get_energy_forces(
            pot='h',
            pot2=False,
            #potparam=args.potparam,
            #hessefile = "HesseMatrix_sphinx",
            h = h,
            #hessefile1nn = args.inputfile1nn,
            #coordfile_cart = "cartesian_coords",
            coord_cart = posaktuell,
            only_return_forces_harmonic = True,
            verbose=self.verbose
            )
        return

    def pot_run(self):
        ''' executes the necessary functions to calculate energies/forces/.... '''
        self.get_or_load_parameters()
        self.get_or_load_positions()
        self.repeat_supercell()
        self.declare_empty_arrays()
        self.save_information_of_run_initialize()
        self.loop_over_atoms()
        if self.mix_uref_uharmonic == True:
            self.get_harmonic_energy_forces()  # is executed only if necessary
        self.sum_up_forces()
        self.check_forces_sum()
        self.print_forces()
        return

    #######################################################################################
    # ADDITIONAL ANALYSIS SCRIPTS
    #######################################################################################
    def from_OUTCAR_create_everything_necessary(self):
        ''' cell, ... '''
        if os.path.isfile("cell") != True:
            if os.path.isfile("OUTCAR") == True or os.path.isfile("OUTCAR.gz") == True:
                print "cell"
                utils.run2("OUTCAR_cell-last-cartesian-ARRAY.sh > cell")
        if os.path.isfile("dUdL") != True:
            if os.path.isfile("forces_OUTCAR") != True:
                if os.path.isfile("OUTCAR") == True or os.path.isfile("OUTCAR.gz") == True:
                    print "forces_OUTCAR"
                    utils.run2("rm -f forces_OUTCAR; OUTCAR_forces-last-ARRAY.sh > forces_OUTCAR")
            if os.path.isfile("cartesian_coords") != True:
                if os.path.isfile("OUTCAR") == True or os.path.isfile("OUTCAR.gz") == True:
                    print "cartesian_coords"
                    utils.run2("rm -f cartesian_coords; OUTCAR_positions-last-ARRAY.sh > cartesian_coords")
        if os.path.isfile("dUdL") == True:
            if os.path.isfile("POSITIONs") != True:
                if os.path.isfile("OUTCAR") == True or os.path.isfile("OUTCAR.gz") == True:
                    print "POSITIONs"
                    utils.run2("extractPOSITIONS.sh")
        if os.path.isfile("EqCoords_direct") != True:
            print "EqCoords_direct"
            if os.path.isfile("../EqCoords_direct") == True:
                eq = np.loadtxt("../EqCoords_direct")
                if eq.shape[0] == self.numberofatoms:
                    shutil.copyfile("../EqCoords_direct","EqCoords_direct")
        if os.path.isfile("EqCoords_direct_") == True:
            atoms = np.loadtxt("EqCoords_direct_").shape[0]
        if os.path.isfile("EqCoords_direct") != True:
            print utils.printred("EqCoords_direct still missing ....")
            sys.exit()

            checkfor = False
            if atoms == 32:
                checkfor = "/Users/glensk/Thermodynamics/utilities/fcc/EqCoords_direct_fcc_2x2x2sc"
            if atoms == 108:
                checkfor = "/Users/glensk/Thermodynamics/utilities/fcc/EqCoords_direct_fcc_3x3x3sc"
            if type(checkfor) != bool:
                if os.path.isfile(checkfor) == True:
                    shutil.copyfile(checkfor,"EqCoords_direct")


        if os.path.isfile("EqCoords_direct") != True:
            sys.exit("you need to get EqCoords_direct")
        check0 = np.loadtxt("EqCoords_direct")


        # if we rund dUdLharmonic we still dont need cartesian_coords

        #if os.path.isfile("cartesian_coords") != True:
        #    sys.exit("you need to get cartesian_coords")
        if os.path.isfile("cartesian_coords") == True:
            try:
                check1 = np.loadtxt("cartesian_coords")
            except ValueError:
                sys.exit("cartesian coords file corrupt; check for XXCAR")
            if check0.shape != check1.shape:
                print "check0:",check0
                print "check1:",check1
                sys.exit("check0.shape != check1.shape")
        return

    def getpos_from_POSITIONs(self,schritte,pos,atoms):
        ''' returnes from POITIONs file (pos) the only coordinates of a certain given step '''
        if type(schritte) == np.ndarray:
            for i in schritte:
                print pos[atoms*schritte:atoms*schritte+atoms]
        return pos[atoms*schritte:atoms*schritte+atoms]


    def dudlrefharmonic(self, dudlnew = False, dudlnewharmonic = False, schritte_up_to = False):
        ''' creates dUdLref and or dUdLharmonic'''
        dudlposcm = False
        sb = False
        se = False
        verbose = False
        if os.path.isfile("POSITIONs") != True:
            self.from_OUTCAR_create_everything_necessary()
        if os.path.isfile("POSITIONs") != True:
            sys.exit("need POSITIONs")
        if os.path.isfile("dUdL") != True:
            self.from_OUTCAR_create_everything_necessary()
        if os.path.isfile("dUdL") != True:
            sys.exit("need dUdL")
        if os.path.isfile("cell") != True:
            self.from_OUTCAR_create_everything_necessary()
        if os.path.isfile("cell") != True:
            sys.exit("need cell")
        cell = np.loadtxt("cell")
        if os.path.isfile("EqCoords_direct") != True:
            self.from_OUTCAR_create_everything_necessary()
        if os.path.isfile("EqCoords_direct") != True:
            sys.exit("need EqCoords_direct")
        pos0rel = np.loadtxt("EqCoords_direct")
        pos = np.loadtxt("POSITIONs")
        listdudl = np.loadtxt("dUdL")
        eqpos_rrel = np.loadtxt("EqCoords_direct")

        h = False
        if dudlnewharmonic != False:
            h = hesse.read_Hessematrix("HesseMatrix_sphinx")


        ########################################################################
        # from here we dont have to load anything
        ########################################################################
        listdudlnew = np.copy(listdudl)
        listdudlnew[:,5] = np.nan  # Uref
        listdudlnew[:,6] = np.nan  # dUdL
        listdudlnewharmonic = copy.deepcopy(listdudlnew)
        schritte_dudl = listdudl.shape[0]+1

        atoms = eqpos_rrel.shape[0]
        print "atoms:",atoms
        schritte = pos.shape[0]/atoms
        print "schritte:",schritte
        out = np.zeros((schritte,atoms))

        fhar_vs_steps = False
        fdft_vs_steps = False
        fref_vs_steps = False
        eref_vs_steps = False
        elon_vs_steps = False
        etox_vs_steps = False

        self.get_or_load_parameters()
        lonvecs = [ 'lonvecnorm'+str(i) for i in self.shells]
        lonvecsmean = [ 'lonvecnormmean'+str(i) for i in self.shells]

        self.df = pd.DataFrame(columns = [ 'step', 'u', 'uref', 'dudlref', 'uhar', 'dudlhar', 'pos', 'f', 'fref', 'fhar' ]+lonvecs+lonvecsmean).set_index(['step'])
        self.df = self.df.astype(object)
        if os.path.isfile(self.dudlpkl) == True:
            df = pd.read_pickle(self.dudlpkl)
            if df.shape[0] > 3:
                self.df = df

        self.coord0_rrel = pos0rel
        self.cell = cell

        self.longvec_all_allsteps = [ False for i in np.arange(len(self.shells)+1) ]
        self.longvecnorm_all_allsteps = [ np.array([]) for i in np.arange(len(self.shells)+1) ]
        self.longvecnorm_all_max_vs_steps = [ False for i in np.arange(len(self.shells)+1) ]
        self.longvecnorm_all_min_vs_steps = [ False for i in np.arange(len(self.shells)+1) ]

        #schritte = 200
        print "steps_up_to:",self.steps_up_to,type(self.steps_up_to)
        if type(self.steps_up_to) != bool:
            schritte = self.steps_up_to

        self._info_lonvec_all = False
        self._info_lonvec = False
        for i in np.arange(schritte):
            if schritte_up_to != False:
                if schritte_up_to + 1 > i:
                    return
            #print utils.printred("schritte:"+str(i))
            if self.df.shape[0] > i+1:
                continue
            if i > 0 and type(sb) != bool and i < sb:
                continue
            if i > 0 and type(se) != bool and i > se:
                continue
            posforceaktuell = self.getpos_from_POSITIONs(i,pos,atoms)
            posaktuell = posforceaktuell[:,:3]  # those are the positions generated by the run
            fDFT = posforceaktuell[:,3:]
            print "CURRENT:",os.getcwd()
            #old
            #self.df.loc[i] = [ np.nan in range(self.df.shape[1]-1) ]
            #self.df.loc[i] = [False]*(self.df.shape[1]-1)
            self.df.loc[i] = False
            self.df = self.df.astype(object)
            self.df.loc[i]['step'] = int(i)
            if self.write_everything_to_pkl == True:
                self.df.loc[i]['pos'] = posaktuell
                self.df.loc[i]['f'] = fDFT

            self.coord_cart = posaktuell

            # here we can actually access anything from our class
            #import cProfile
            #cProfile.run('self.pot_run()')
            #############################################################################
            # calculate energy
            #############################################################################
            self.pot_run()
            #self.pot_run()
            #cProfile.run('self.pot_run()')
            forces = self.forces
            energymevlong = self.energymevlong
            energymevtopx = self.energymevtopx

            #print self._info_tveclonallanglevsnorm
            #print "---------------------->>>>",self._info_tveclonallanglevsnorm.shape
            #ka = np.loadtxt("a_vs_b.dat")
            #kb = np.concatenate((ka,self._info_tveclonallanglevsnorm),axis=0)
            #np.savetxt("a_vs_b.dat",kb)
            if type(self._info_lonvec_all) != bool:
                print "self._info_lonvec_all.shape",self._info_lonvec_all.shape
            else:
                print "self._info_lonvec_all.shape",type(self._info_lonvec_all)
            #print "self._info_lonvec.shape:",self._info_lonvec.shape

            self._info_lonvec_all = utils.append_row_to_2d_array(\
                    self._info_lonvec_all,self._info_lonvec)
            #self._info_lonvec_all = np.concatenate((self._info_lonvec_all,self._info_lonvec),axis=0)
            #ka = np.loadtxt("lonvecallxyz.dat")
            #print "kkk:",self._info_tveclonallx.shape
            #print "kk2:",self._info_tveclonallx[0].shape
            #print "kk2:",self._info_tveclonallx[1].shape
            #print "b  :",ka.shape
            #kb = np.concatenate((ka,self._info_tveclonallx.flatten()),axis=0)
            #np.savetxt("lonvecallxyz.dat",kb)

            if dudlnewharmonic != False:
                # eHAR we could in principle just take from the dUdL file
                fHAR,eHAR = hesse.get_energy_forces(
                    pot='h',
                    pot2=False,
                    #potparam=args.potparam,
                    #hessefile = "HesseMatrix_sphinx",
                    h = h,
                    #hessefile1nn = args.inputfile1nn,
                    #coordfile_cart = "cartesian_coords",
                    coord_cart = posaktuell,
                    only_return_forces_harmonic = True,
                    verbose=verbose
                    )
                eHAR = eHAR*1000/(atoms-1)
                self.df.loc[i]['uhar'] = eHAR
                if self.write_everything_to_pkl == True:
                    self.df.loc[i]['fhar'] = fHAR
                fhar_vs_steps = utils.append_row_to_2d_array(fhar_vs_steps, np.insert(fHAR.flatten(),0,i).flatten())
            energymev = self.energy*1000/(atoms-1)
            eNEW = self.energymev
            print "listdudl.shape:",listdudl.shape,"i:",i,"energymev:",self.energymev,energymev
            self.df.loc[i]['uref'] = self.energymev
            if self.write_everything_to_pkl == True:
                self.df.loc[i]['fref'] = self.forces

            eDFT  = listdudl[i-1][4]
            #eREF = listdudl[i-1][5]
            #print "eDFT:",eDFT
            #print "eHAR:",eHAR
            #print "eNEW:",eNEW
            ##print "eREF:",eREF
            #print "energymev:",energymev
            if i == 0:
                self.df.loc[i]['u'] = 0.0
                self.df.loc[i]['uref'] = 0.0
                self.df.loc[i]['dudlref'] = 0.0
            else:
                self.df.loc[i]['u'] = eDFT


            self.df.loc[i]['dudlref'] =  self.df.loc[i]['u'] - self.df.loc[i]['uref']

            if dudlnewharmonic != False:
                self.df.loc[i]['dudlhar'] =  self.df.loc[i]['u'] - self.df.loc[i]['uhar']
            else:
                self.df.loc[i]['dudlhar'] = 0

            print utils.printred(">>>>>>"*4+" DOS_POSITIONS_auswerten.py "+"<<<<<<"*4)
            print utils.printred(">>>>>>>",i,"energymev:",energymev,"eDFT:",eDFT)  #,"eREF:",eREF)
            print utils.printred(">>>>>>"*4+" DOS_POSITIONS_auswerten.py "+"<<<<<<"*4)
            listdudlnew[i-1][5] = energymev
            listdudlnew[i-1][6] = eDFT - energymev
            listdudlnew[i-1][7] = 0 #listdudlnew[:,6].mean()
            listdudlnew[i-1][8] = 0
            listdudlsave = listdudlnew[:i]
            ##############################################################################
            # here we save dUdLref
            ##############################################################################
            if dudlnew != False:
                np.savetxt("dUdLref",listdudlsave,fmt="%7.2f%10.1f%9.1f%9.1f%14.2f%16.2f%14.2f%10.2f%10.2f",header=" step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset")

            ##############################################################################
            # here we save dUdLharmonic
            ##############################################################################
            if dudlnewharmonic != False:
                listdudlnewharmonic[i-1][5] = eHAR
                listdudlnewharmonic[i-1][6] = eDFT - eHAR
                listdudlnewharmonic[i-1][7] = 0 #listdudlnewharmonic[:,6].mean()
                listdudlnewharmonic[i-1][8] = 0
                listdudlsaveharmonic = listdudlnewharmonic[:i]
                np.savetxt("dUdLharmonic",listdudlsaveharmonic,fmt="%7.2f%10.1f%9.1f%9.1f%14.2f%16.2f%14.2f%10.2f%10.2f",header=" step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset")

            eref_vs_steps = utils.append_row_to_2d_array(eref_vs_steps, np.array([i,energymevlong+energymevtopx]))
            elon_vs_steps = utils.append_row_to_2d_array(elon_vs_steps, np.array([i,energymevlong]))
            etox_vs_steps = utils.append_row_to_2d_array(etox_vs_steps, np.array([i,energymevtopx]))
            fdft_vs_steps = utils.append_row_to_2d_array(fdft_vs_steps, np.insert(fDFT.flatten(),0,i).flatten())
            fref_vs_steps = utils.append_row_to_2d_array(fref_vs_steps, np.insert(self.forces.flatten(),0,i).flatten())


            # for DOS
            #self.longvec_all_allsteps = [ False for i in np.arange(len(shells)+1) ]
            #self.longvecnorm_all_allsteps = [ np.array([]) for i in np.arange(len(shells)+1) ]
            #for shell in np.arange(len(self.shells)):
            for shell in self.shells:
                # shell has to be increased by 1 (there is no 0 shell
                #self.longvec_all_allsteps[shell] = np.vstack((self.longvec_all_allsteps[shell],self.longvec_all[shell]))
                self.longvecnorm_all_allsteps[shell] = np.append(self.longvecnorm_all_allsteps[shell],self.longvecnorm_all[shell])


                self.longvecnorm_all_max_vs_steps[shell] = utils.append_row_to_2d_array(self.longvecnorm_all_max_vs_steps[shell],np.array([i,self.longvecnorm_all[shell].max()]))
                self.longvecnorm_all_min_vs_steps[shell] = utils.append_row_to_2d_array(self.longvecnorm_all_min_vs_steps[shell],np.array([i,self.longvecnorm_all[shell].min()]))
                self.df.loc[i]['lonvecnorm'+str(shell)] = self.longvecnorm_all[shell]
                self.df.loc[i]['lonvecnormmean'+str(shell)] = self.longvecnorm_all[shell].mean()

            ###############################################################################
            # save all data every 10 steps
            # - longvecnorm_1
            # - longvecnorm_1_DOS
            # - longvecnorm_1_max_vs_steps
            # - longvecnorm_1_min_vs_steps
            # - elon_vs_steps
            # - etox_vs_steps
            #
            ###############################################################################
            if i in np.arange(schritte)[::50]:
                if i == 0:
                    continue
                self.df.to_pickle(self.dudlpkl)
                #np.savetxt("kkk50",self.longvec_all_mapped_to_first_quadrant)
                #print self.longvec_all_mapped_to_first_quadrant.shape,"jo"

            if i == schritte-1:
                self.df.to_pickle(self.dudlpkl)

            if self.write_vecs_to_disc == True or self.params.save_vecs_to_file_for_DOS == True:
                if i in np.arange(schritte)[::50]:
                    if i == 0:
                        continue
                    #for shell in np.arange(len(shells)):
                    for shell in self.shells:
                        np.savetxt("longvecnorm_"+str(shell),self.longvecnorm_all_allsteps[shell])
                        dos = utils.getDOS_of_1d_array(self.longvecnorm_all_allsteps[shell])
                        np.savetxt("longvecnorm_"+str(shell)+"_DOS",dos)
                        np.savetxt("longvecnorm_"+str(shell)+"_max_vs_steps",self.longvecnorm_all_max_vs_steps[shell],fmt="%.0f %.6f")
                        np.savetxt("longvecnorm_"+str(shell)+"_min_vs_steps",self.longvecnorm_all_min_vs_steps[shell],fmt="%.0f %.6f")

                    np.savetxt("lonvec_all",self._info_lonvec_all,fmt="%.6f")
                    np.savetxt("elon_vs_steps",elon_vs_steps,fmt="%.0f %.6f")
                    np.savetxt("etox_vs_steps",etox_vs_steps,fmt="%.0f %.6f")
                    np.savetxt("fdft_vs_steps",fdft_vs_steps,fmt="%.6f")
                    np.savetxt("fref_vs_steps",fref_vs_steps,fmt="%.6f")
                    #np.savetxt("fref_vs_steps",fref_vs_steps,fmt="%.6f")
                    ka = fdft_vs_steps[:,1:].flatten()
                    print "shape.ka:",ka.shape
                    if ka.shape[0] == 0:
                        print "YOU HAVE TO DO AN extractPOSITIONs.sh to get forces in POSITIONS file"
                    kb = fref_vs_steps[:,1:].flatten()
                    #print "shape.kb:",kb.shape
                    fref_vs_fdft = np.zeros((ka.shape[0],2))
                    print "shape.kb:",kb.shape
                    fref_vs_fdft[:,0] = ka
                    fref_vs_fdft[:,1] = kb
                    np.savetxt("fref_vs_fdft.dat",fref_vs_fdft,fmt="%.6f %.6f")

                    minimum = ka.min();
                    if kb.min() < minimum:
                        minimum = kb.min()
                    maximum = ka.max();
                    if kb.max() > maximum:
                        maximum = kb.max()
                    np.savetxt("fdft_vs_fdft.dat",np.array([[minimum*1.2,minimum*1.2],[maximum*1.2,maximum*1.2]]),fmt="%.6f %.6f")
                    if dudlnewharmonic != False:
                        kc = fhar_vs_steps[:,1:].flatten()
                        fhar_vs_fdft = np.zeros((kc.shape[0],2))
                        fhar_vs_fdft[:,0] = ka
                        fhar_vs_fdft[:,1] = kc
                        np.savetxt("fhar_vs_fdft.dat",fhar_vs_fdft,fmt="%.6f %.6f")

        if self.write_vecs_to_disc == True or self.params.save_vecs_to_file_for_DOS == True:
            if type(self._info_lonvec_all) != bool:
                np.savetxt("lonvec_all",self._info_lonvec_all,fmt="%.6f")
            if type(elon_vs_steps) != bool:
                np.savetxt("elon_vs_steps",elon_vs_steps,fmt="%.6f")
            if type(etox_vs_steps) != bool:
                np.savetxt("etox_vs_steps",etox_vs_steps,fmt="%.6f")
            if type(fdft_vs_steps) != bool:
                np.savetxt("fdft_vs_steps",fdft_vs_steps,fmt="%.6f")
            if type(fref_vs_steps) != bool:
                np.savetxt("fref_vs_steps",fref_vs_steps,fmt="%.6f")
        return

    def eval_dudl_pkl(self, dudlpkl = False, getdos = False, verbose = False):
        ''' take filename (usually dUdL.pkl), load it and evaluate '''
        if dudlpkl == False:
            dudlpkl = self.dudlpkl
        else:
            dudlpkl = dudlpkl
        if verbose == True:
            print "dudlpkl:",dudlpkl
        if os.path.isfile(dudlpkl) == True:
            df = pd.read_pickle(dudlpkl)
            if df.shape[0] > 3:
                self.df = df
        else:
            sys.exit("dudlpkl "+dudlpkl+" not found")

        if getdos:
            # to get the pot we need the element
            element = utils.run2("OUTCAR_elements.sh").rstrip()
            print "element:",element
            import my_atom
            import getDOS
            reload(getDOS)
            import pot_parametrize
            self.element = my_atom.atom([element])
            print "self.element:",self.element
            #self.df.loc[:,['lonvecnormmean1','dudlhar']].plot(x='lonvecnormmean1',y='dudlhar',style='.')
            #np.savetxt("uval",self.df.u.values[1:])
            #np.savetxt("urefharval",self.df.urefhar.values[1:])
            self.udos = getDOS.get_dos(self.df.u.values[1:])
            self.upot = pot_parametrize.dos_to_pot(self.udos, self.element.melting_rounded[0])
            #np.savetxt("udos",udos)
            if self.df.uhar.values[1] != False:
                self.udoshar = getDOS.get_dos(self.df.uhar.values[1:])
                self.upothar = pot_parametrize.dos_to_pot(self.udoshar, self.element.melting_rounded)
        #np.savetxt("udoshar",udoshar)
        self.dudlhar_mean = False
        self.dudlref_mean = False
        self.dudlhar_std = False
        self.dudlref_std = False
        for i in np.arange(2,self.df.dudlhar.values[1:].shape[0]):
            addrow = [ i-1, self.df.dudlhar.values[1:i].mean()]
            self.dudlhar_mean = utils.append_row_to_2d_array(inarray = self.dudlhar_mean, addrow=addrow)
            addrow = [ i-1, self.df.dudlref.values[1:i].mean()]
            self.dudlref_mean = utils.append_row_to_2d_array(inarray = self.dudlref_mean, addrow=addrow)
            addrow = [ i-1, self.df.dudlhar.values[1:i].std()]
            self.dudlhar_std = utils.append_row_to_2d_array(inarray = self.dudlhar_std, addrow=addrow)
            addrow = [ i-1, self.df.dudlref.values[1:i].std()]
            self.dudlref_std = utils.append_row_to_2d_array(inarray = self.dudlref_std, addrow=addrow)
        #print dudlpkl+" evaluated"
        #print ""
        return

    def plot_dudl_pkl(self, whichplot, dudlpkls = False):
        import matplotlib.pyplot as plt
        plt.clf()
        ##################################################################################
        # get pkl filename is not specified
        ##################################################################################
        if dudlpkls == False:
            if type(self.dudlpkl) == str:
                dudlpkls_all = [ self.dudlpkl ]
            if type(self.dudlpkl) == list:
                dudlpkls_all = self.dudlpkl
        else:
            dudlpkls_all = dudlpkls
        if 'dUdL_harmonic.pkl' in dudlpkls_all: dudlpkls_all.remove('dUdL_harmonic.pkl')
        dudlpkls_all = utils.lsn(dudlpkls_all)
        #print "dudlpkls",dudlpkls

        ##################################################################################
        # get info / plots
        ##################################################################################
        if whichplot == 'std_conv':
            #dudlpkls_all = [ "dUdL_17.7_2.12_1.0.pkl" ]
            #dudlpkls_all = glob.glob("*.pkl")
            if len(dudlpkls_all) == 1:
                as_ = np.arange(0.5,1.01,0.01)
                bs_ = np.arange(0.0,0.6,0.01)
            else:
                as_ = [ 1.00 ]  # ip
                bs_ = [ 0.00 ]  # harmonic
            self.results_compare = pd.DataFrame(columns = [ 'a', 'b', 'alpha', 'B', 'std' ])

            if os.path.isfile("dUdL_harmonic.pkl") == True:
                pot.eval_dudl_pkl("dUdL_harmonic.pkl")
                harmdudl = "dUdL_harmonic.pkl"
            else:
                sys.exit("harmdudl ../dUdL_harmonic.pkl not found")
            self.dfharm = copy.copy(self.df)
            print utils.printgreen("found harmonic "+str(harmdudl)+" file")



        self.ip = False
        #minimum = [a,b,c,d,std]
        minimum = [0,0,0,0, 1000000000000000000000000000000000000]
        for i in dudlpkls_all:
            #print ""
            print "I:",i
            self.eval_dudl_pkl(dudlpkl = i)
            if whichplot == 'ene_dos':
                plt.plot(self.udos[:,0],self.udos[:,1],'r-')
                plt.plot(self.udoshar[:,0],self.udoshar[:,1],'k-')
                #plt.plot(self.upot[:,0],self.upot[:,1],'r-')
                #plt.plot(self.upothar[:,0],self.upothar[:,1],'k-')
                plt.xlabel("vibrational energy (meV)")
                plt.ylabel("DOS (arb. units)")
            if whichplot == 'std_conv_plot':
                np.savetxt("std_"+i,self.dudlref_std)

            if whichplot == 'std_conv':
                for a in as_:
                    print "a:",a,"minimum:",minimum
                    for b in bs_:
                        aip = self.df.uref*a
                        bhar = self.dfharm.uhar[:self.df.uref.shape[0]]*b    # due to definition in dario alfes paper
                        uref =  aip+bhar  # as defined in alfes paper
                        self.dudlip = self.dfharm.u[:self.df.uref.shape[0]] - uref
                        plt.plot(self.dudlhar_std[:,0],self.dudlhar_std[:,1],'k-',label='harmonic')
                        c = float(i.split("_")[1])
                        d = float(i.split("_")[2])
                        std = round(self.dudlip[:1000].std(),2)
                        #print "a:",a,"b:",b,float(i.split("_")[1]),float(i.split("_")[2]),str(round(self.dudlip[:1000].std(),2))
                        if std < minimum[4]:
                            minimum = [ round(a,2), round(b,2), round(c,2), round(d,2), round(std,2) ]


                        addrow = [a,b,float(i.split("_")[1]),float(i.split("_")[2]),float(round(self.dudlip[:1000].std(),2))]
                        self.ip = utils.append_row_to_2d_array(inarray = self.ip, addrow=addrow)

                        self.results_compare.loc[self.results_compare.shape[0]] = [ a, b, float(i.split("_")[1]),float(i.split("_")[2]),round(self.dudlip[:1000].std(),2) ]
                        plt.plot(self.dudlref_std[:,0],self.dudlref_std[:,1],label='ref '+i)
                        plt.xlabel("steps")
                        plt.ylabel("< dUdL.std > (meV/atom)")
                print "minimum:",minimum
            if whichplot == 'dudl_conv':
                plt.plot(self.dudlhar_mean[:,0],self.dudlhar_mean[:,1],'k-',label='harmonic')
                plt.plot(self.dudlref_mean[:,0],self.dudlref_mean[:,1],'r-',label='harmonic')
                plt.xlabel("steps")
                plt.ylabel("< dUdL > (meV/atom)")
            if whichplot == 'ene_vs_disp':
                ka = pd.read_pickle("dUdL.pkl")
                #ka.sort(columns='lonvecnormmean1').loc[:,['lonvecnormmean1','dudlhar']].plot(x='lonvecnormmean1',y='dudlhar',style='.')
                ka.loc[:,['lonvecnormmean1','dudlhar']].plot(x='lonvecnormmean1',y='dudlhar',style='.')
                a = ka.loc[:,['lonvecnormmean1','dudlhar']].values
                np.savetxt("ene_vs_disp_harmonic",a)
                b = np.array([a[:,0] - a[0,0],a[:,1]]).transpose()
                np.savetxt("ene_vs_disp_harmonic_shifted",b)

                ka.loc[:,['lonvecnormmean1','dudlref']].plot(x='lonvecnormmean1',y='dudlref',style='.')
                a = ka.loc[:,['lonvecnormmean1','dudlref']].values
                np.savetxt("ene_vs_disp_ref",a)
                b = np.array([a[:,0] - a[0,0],a[:,1]]).transpose()
                np.savetxt("ene_vs_disp_ref_shifted",b)

                #import getDOS
                #np.savetxt("uval",ka.u.values[1:])
                #np.savetxt("urefharval",ka.urefhar.values[1:])
                #udos = getDOS.get_dos(ka.u.values[1:])
                #np.savetxt("udos",udos)
                #udoshar = getDOS.get_dos(ka.urefhar.values[1:])
                #np.savetxt("udoshar",udoshar)

            if whichplot != 'std_conv_plot':
                #legend = plt.legend(loc='upper center', shadow=True)
                legend = plt.legend(loc='best',fancybox=True)
                legend.get_frame().set_alpha(0.5)
                plt.grid(True)
                plt.ion()

        if whichplot != 'std_conv_plot':
            print ""
            v1 = np.unique(pot.ip[:,2])
            for a in v1:
                print "a:",a
                data = self.ip[np.nonzero(self.ip[:,2] == a)][:,3:]
                np.savetxt("auswertung_"+str(a),data)
        return

    def auswertung_lonvec_all(self):
        self.auswertung_DOS2d()
        sys.exit()
        self.al = np.loadtxt("lonvec_all")
        self.al[:,2] = np.absolute(self.al[:,2])

        def lonvec_all_z(von,bis):
            self.alv01 = self.al[np.nonzero((self.al[:,2] > von) & (self.al[:,2] < bis))]
            np.savetxt("lonvec_all_"+str(von)+"_"+str(bis),self.alv01)
            np.savetxt("lonvec_all_"+str(von)+"_"+str(bis)+"_forDOS",self.alv01[:,[0,1]].flatten())
            utils.run2("getDOS.sh lonvec_all_"+str(von)+"_"+str(bis)+"_forDOS; mv DOS DOS_"+str(von)+"_"+str(bis)) # to create cell

        lonvec_all_z(0.0,0.1)
        lonvec_all_z(0.1,0.2)
        lonvec_all_z(0.2,0.3)
        lonvec_all_z(0.3,0.4)
        lonvec_all_z(0.4,0.5)
        lonvec_all_z(0.5,0.6)
        lonvec_all_z(0.6,0.7)
        lonvec_all_z(0.7,0.8)
        lonvec_all_z(0.8,0.9)
        lonvec_all_z(0.9,2.0)
        return

    def auswertung_DOS2d(self):
        filename = "DOS2d_0.0_0.1"
        filename = "./ti/Ir/31__PTS_dosall_from_2x2x2sc__LON_quer_3x3x3kp_all_morse_LONADD_YES_LON2_NONE_LON2ADD__NO_TOX_NONE_TOY_NONE_TOZ_NONE/lambda1.0/tests/longecxyz/DOS2d_0.0_0.1"
        self.al = np.loadtxt(filename)
        von = 0.01
        bis = 0.013
        self.alv01 = self.al[np.nonzero((self.al[:,2] > von) & (self.al[:,2] < bis))]
        from math import atan2
        print self.alv01[:,[0,1]]
        self.alv01[:,[0,1]].sort(key=lambda c:atan2(c[0], c[1]))
        np.savetxt(filename+"__"+str(von)+"_"+str(bis),self.alv01[:,[0,1]])

    def help(self):
        '''
        to masure time:

            python -m cProfile -o prof /Users/glensk/Thermodynamics/python_thermodynamics/pot_energy_forces.py -dudlref -s 20

            pstats.Stats('prof').strip_dirs().sort_stats("cumulative").print_stats(30)'''

        p = argparse.ArgumentParser(description='''help string''')

        p.add_argument('-harmonic', action='store_true',
                help='get harmonic forces and energies', default=False)
        p.add_argument('-dudlpklall', action='store_true',
                help='write_everything_to_pkl add forces{DFT,ref,harmonic} poitions to pkl file', default=False)
        p.add_argument('-dudlref', action='count',
                help='go thround POSITIONs and create dUdLref', default=False)
        p.add_argument('-dudlharmonic', action='count',
                help='go thround POSITIONs and create dUdLharmonic (+fhexplode)', default=False)
        p.add_argument('-auswertung_lonvec_all', action='count',
                help='evaluates lonvec_all file', default=False)
        p.add_argument('-ci', action='count',
                help='create all necessary inputfiles from OUTCAR', default=False)
        p.add_argument('-wtd', action='store_true',
                help='write_vecs_to_disc like {elon,etox,fdft,fref}_vs_step longvecnorm_{1,2,3,4}_{max_vs_steps,min_vs_steps,DOS}', default=False)

        p.add_argument('-dudlpkl',
                help='change pkl filename of dUdL file; default: dUdL.pkl',nargs='+', default=self.dudlpkl)
        p.add_argument('-evalpkl', action='store_true',
                help='get dos, pot, from (dUdL).pkl file', default=False)
        p.add_argument('-p', choices=['std_conv', 'std_conv_plot', 'dudl_vs_longvecnorm', 'ene_vs_disp', 'ene_dos', 'dudl_conv'],
                help='plot from pkl file one of the above choices', default=False)
        p.add_argument('-s', '--steps',
                help='dont run through all steps of dUdL file but rahter up to (int)', default=False)

        p.add_argument('-plot2', action='store_true',
                help='plot energies as a function of displacement', default=False)

        p.add_argument('-t', '--measure_time', action='store_true',
                help='measure and show timings of single steps of execution', default=False)
        p.add_argument('-ss', type=int,choices=[1,2,3,4,5],help='define number of shells for which forces will be shown, default: for the shells for which force functions are defined',default=False)
        p.add_argument('-v', '--verbose',action='count',
                help='verbosity level: v:verbose, vv:even more verbose, vvv, vvvv ...', default=False)
        #argcomplete.autocomplete(p)
        args = p.parse_args()

        self.verbose = args.verbose
        self.dudlpkl = args.dudlpkl
        self.measure_time = args.measure_time
        self.steps_up_to = args.steps
        self.write_vecs_to_disc = args.wtd
        self.write_everything_to_pkl = args.dudlpklall
        self.ss = args.ss

        if self.steps_up_to != False:
            self.steps_up_to = int(self.steps_up_to)
        print "self.pkl:",self.dudlpkl

        if args.auswertung_lonvec_all:
            self.auswertung_lonvec_all()
            sys.exit()

        if args.ci:
            self.from_OUTCAR_create_everything_necessary()
            sys.exit()

        if args.dudlref or args.dudlharmonic:
            self.dudlrefharmonic(args.dudlref,args.dudlharmonic)
            sys.exit()

        if args.evalpkl:
            self.eval_dudl_pkl()
            sys.exit()

        if args.p:
            print "args.p:",args.p
            self.plot_dudl_pkl(args.p)
            sys.exit()

        return args


if __name__ == '__main__':

    pot = pot_energy_forces_class()
    args = pot.help()  # this eventually starts the calculatino if option is given
    pot.pot_run()

    if args.dudlref != False or args.dudlharmonic != False or args.p != False or args.ci != False:
        pass
    else:
        np.savetxt("forces",pot.forces)
        np.savetxt("energy",np.array([pot.energy]))
