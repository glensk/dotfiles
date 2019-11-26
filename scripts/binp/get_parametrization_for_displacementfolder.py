#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from ase import Atoms
from scipy.fftpack import fft
from scipy.fftpack import ifft
from lmfit import Model
from scipy.optimize import curve_fit
import glob,copy
import sys,os,argparse
import utils_rename as utils
import my_atom
import myutils as my
from subprocess import call
from numpy import linalg as LA

# Matplotlib for additional customization
from matplotlib import pyplot as plt

import pot_parametrize
import pot_energy_forces
import hesse


def getindex(array,findvec):
    findvec[0] = round(findvec[0],13)
    findvec[1] = round(findvec[1],13)
    findvec[2] = round(findvec[2],13)
    try:
        #return np.where((array[:,0] == findvec[0]) & (array[:,1] == findvec[1]) & (array[:,2] == findvec[2]))[0][0]
        return np.where((np.round(array[:,0],13) == findvec[0]) & (np.round(array[:,1],13) == findvec[1]) & (np.round(array[:,2],13) == findvec[2]))[0][0]
    except IndexError:
        print("ERROR IN getindex() (DEFINITION)")
        print("array: (==the positions)")
        print(array)
        sys.exit()

def restrict_xrange_of_array(array,xmin=-np.inf,xmax=np.inf):
    array = array[array[:,0] > xmin]
    array = array[array[:,0] < xmax]
    return array

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
       return v
    return v / norm

def getforce(distvec,force):
    '''the force should be a scalar'''
    return - normalize(np.array(distvec))*force



class get_one_disp():
    def __init__(self,element,disp):
        self.sc = 5
        self.element = element
        self.pos = False
        self.forces = False
        self.disp = disp
        self.elements = ["Pd","Pb","Al","Ir","Cu","Rh","Pt","Ag","Au"]
        self.path="/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/"


        #################################################################
        # for one particular displacement
        #################################################################
        if True: # for one particular displacement
            self.folder=glob.glob(
                self.path+str(self.element)+"/"+str(self.sc)+"x"+str(self.sc)+"x"+str(self.sc)+
                "sc_quer*/*Ang_"+str(self.disp)+'/')

            #print("--aaa",self.folder)

            if len(self.folder) == 0:
                sys.exit("No folder found with element \""+self.element+"\".")
            if len(self.folder) != 1:
                for i in self.folder:
                    print('folder:',i)
                sys.exit("more then one folder")
            else:
                self.folder = self.folder[0]

            if self.element not in self.elements:
                sys.exit("The element "+self.element+" is not known!")

            self.pos = np.loadtxt(self.folder+'pos')
            self.forces = np.loadtxt(self.folder+'forces')
            self.alat = float(self.folder.split("vasp4/")[1].split("Ang")[0])
            #print(self.folder)
            #print(self.element,self.alat)

            self.idx_1_1_0 = getindex(self.pos/self.alat,np.array([1,1,0]))
            self.idx_05_05_0 = getindex(self.pos/self.alat,np.array([0.5,0.5,0]))

            self.p_1_1_0   = self.pos[self.idx_1_1_0]
            self.p_05_05_0 = self.pos[self.idx_05_05_0]

            self.f_1_1_0   = self.forces[self.idx_1_1_0]
            self.f_05_05_0 = self.forces[self.idx_05_05_0]
        return


class get_all_disps():
    def __init__(self,
            element=False,
            sc=False,
            alat=False,
            dofor=False,
            verbose=False,
            verbose2=False,
            only_return_parametrization_file=False,

            folder_alat_lattice_T0K = False,
            shift_parametrization_to_alat = False,
            correct_for_110_forces = False
            ):
        '''
        required:
            element = ["Al"]
            sc = 5

        NOT necessary:
            alat
            dofor = "/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400"
            verbose
        '''
        if folder_alat_lattice_T0K != False:
            if element == "Al": dofor = "Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400"
            else:
                sys.exit("No T0K folder jet defined for "+element)

        if dofor == ".": dofor = os.getcwd()
        if type(dofor) != bool and not os.path.isdir(dofor):
            dofor = "/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/"+dofor
            if not os.path.isdir(dofor):
                sys.exit("dofor "+dofor+" does not exist (1)!")

        self.sc = sc
        self.element = element
        self.alat = alat
        self.dofor = dofor
        self.verbose = verbose
        self.verbose2 = verbose2
        self.only_return_parametrization_file = only_return_parametrization_file
        self.folder_alat_lattice_T0K = folder_alat_lattice_T0K
        self.shift_parametrization_to_alat = shift_parametrization_to_alat
        self.correct_for_110_forces = correct_for_110_forces

        if self.shift_parametrization_to_alat != False:
            if type(self.shift_parametrization_to_alat) not in [np.float64,float]:
                print("self.shift_parametrization_to_alat",self.shift_parametrization_to_alat)
                print("type(self.shift_parametrization_to_alat)",type(self.shift_parametrization_to_alat))
                sys.exit("shift_parametrization_to_alat is set, but not a float!")


        if self.dofor:
            self.element = False
            self.sc = False
            self.alat = False

        self.pos = False
        self.forces = False
        self.elements = ["Pd","Pb","Al","Ir","Cu","Rh","Pt","Ag","Au"]
        self.path="/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/"

        ### make sure that element is defined: either dofor or given
        if self.element == False and self.dofor == False:
            sys.exit("need either element or dofor")

        if self.element == False:
            element = []
            check_for_element = self.dofor.split("/")
            for e in self.elements:
                if e in check_for_element:
                    element.append(e)
            if len(element) == 1:
                self.element = element[0]
            else:
                print('element(s):',element)
                sys.exit("found more or none elements")


        ### element is defined!
        self.rmin_distmax = my_atom.rmin[self.element]         # which range to consider
        self.rmax_distmax = my_atom.rmax[self.element]         # which range to consider
        self.alatT0K      = my_atom.alatT0K[self.element]
        self.alatTmelt    = my_atom.alatTmelt[self.element]
	self.nndistT0K    = self.alatT0K/np.sqrt(2.)
	self.nndistTmelt  = self.alatTmelt/np.sqrt(2.)
        if self.verbose:
            print("self.sc           (1):",self.sc)
            print("self.dofor        (1):",self.dofor)
            print("self.element      (1):",self.element)
            print("self.rmin_distmax (1):",self.rmin_distmax)
            print("self.rmax_distmax (1):",self.rmax_distmax)
            print("self.alatT0K      (1):",self.alatT0K)
            print("self.alatTmelt    (1):",self.alatTmelt)
            print("self.alat         (1):",self.alat     )
            print("self.shift_parametrization_to_alat (1):",self.shift_parametrization_to_alat)

        if self.dofor != False:
            # in this case this will be defined further down
            self.alat = False
            self.nndist = False
            self.sc = False
        else:
            if self.alat == False: self.alat = self.alatTmelt
            self.nndist = self.alat/np.sqrt(2.)
            if self.sc == False:
                sys.exit('please define supercell')


        if self.element == "Ir" and self.dofor == False:
            # here for the correct Tmelt lattice constant
            self.sc = 3
            take_from = "3x3x3sc_3.99Ang_quer_4x4x4kp_vasp4"


        if self.verbose:
            print("self.sc           (2):",self.sc)
            print("self.alat         (2):",self.alat)
        #print('self.nndist',self.nndist)
        #print()


        #################################################################
        # for one particular displacement
        # if dofor    : check for displacemetns
        # if NOT dofor: try to find the folder with the displacements
        # Creates:
        # --> self.dofor (folder with the displacements)
        # --> self.folder_alldisp (list with folder to every disp)
        #################################################################
        ####### get or check self.folder_alldisp
        def print_dofor_self_dofor(dofor,selfdofor,verbose,add="X"):
            if verbose:
                print('dofor             ('+add+'):',dofor)
                print('self.dofor        ('+add+'):',selfdofor)

        if self.verbose:
            print_dofor_self_dofor(dofor,self.dofor,self.verbose,add="X")

        if self.dofor == False:
            self.folder_alldisp = glob.glob(
                self.path+str(self.element)+"/"+str(self.sc)+"x"+str(self.sc)+"x"+str(self.sc)+
                "sc_quer*/"+str(self.alat)+"Ang_*/")
            if self.element == "Ir":
                # /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Ir/3x3x3sc_3.99Ang_quer_4x4x4kp_vasp4
                # 3x3x3sc_3.99Ang_quer_4x4x4kp_vasp4
                self.folder_alldisp = glob.glob(
                    self.path+str(self.element)+"/"+take_from+"/"+str(self.alat)+"Ang_0.*/")
                #print("Ir,!!")
            if len(self.folder_alldisp) == 0:
                self.folder_alldisp = glob.glob(
                self.path+str(self.element)+"/"+str(self.sc)+"x"+str(self.sc)+"x"+str(self.sc)+
                "sc_"+str(self.alat)+"Ang_quer*/")
        else:
            #print('dofor:',self.dofor)
            self.folder_alldisp = glob.glob(self.dofor+"/*Ang_*/")
        if self.verbose:
            print_dofor_self_dofor(dofor,self.dofor,self.verbose,add="Y")


        self.folder_alldisp = utils.list_sorted(self.folder_alldisp)
        if self.verbose:
            print('---------------------------------- all found displacements -------------')
            for fia in self.folder_alldisp:
                print(fia)
            print('---------------------------------- all found displacements -------------')
        if self.verbose:
            print_dofor_self_dofor(dofor,self.dofor,self.verbose,add="Z")

        ####### get or check self.dofor
        recheck_dofor = True
	if type(self.dofor) == str == type(dofor):
            if self.dofor == dofor:
                recheck_dofor = False

	verbose1 = True
	if verbose1 == True:
		print('xx sfd',self.folder_alldisp)

        if recheck_dofor == True:
            str_list = utils.common_prefix(self.folder_alldisp).split("/")
	    if verbose1 == True:
            	print('xx sfd',str_list)
            str_list = filter(None, str_list)
	    if verbose1 == True:
            	print('xx sdf',str_list)
            str_list= "/"+"/".join(str_list[:-1])
	    if verbose1 == True:
            	print('xx sdf',str_list)
            dofor = str_list
	    if verbose1 == True:
            	print('xx dof',dofor)

        if self.verbose:
            print_dofor_self_dofor(dofor,self.dofor,self.verbose,add="A")

        if self.dofor == False:
            self.dofor = os.path.abspath(dofor)
        else:
            if os.path.abspath(self.dofor) == os.path.abspath(dofor):
                self.dofor = os.path.abspath(dofor)
            else:
                print('dofor     ',dofor)
                print('self.dofor',self.dofor)
                print_dofor_self_dofor(dofor,self.dofor,add="B")
                sys.exit('self.dofor is not dofor; Exit')

        print_dofor_self_dofor(dofor,self.dofor,self.verbose,add="C")
        return

    def parametrize_it(self):
        if self.shift_parametrization_to_alat == False:  # the normal case
            parametrizationfile_morse = self.dofor+"/disp_fit.parameters_morse.dat"
        if self.shift_parametrization_to_alat != False:
            parametrizationfile_morse = self.dofor+"/disp_fit.parameters_morse_alat_"+str(self.shift_parametrization_to_alat)+".dat"
        if self.correct_for_110_forces != False and self.shift_parametrization_to_alat != False:
            sys.exit('write this few lines of code')
        if self.correct_for_110_forces != False:
            parametrizationfile_morse = self.dofor+'/disp_fit_corr_110.parameters_morse.dat'


        ## return parametrization file in case it exists already
        if self.only_return_parametrization_file == True and os.path.isfile(parametrizationfile_morse):
            if self.verbose:
                print('parametrizationfile_morse:',parametrizationfile_morse)
            return parametrizationfile_morse

        if len(self.folder_alldisp) == 0:
            sys.exit('no folders found!')

        ###########################################################################
        # Read in the forces
        # create 4.14Ang_0.6Ang/pos
        # create 4.14Ang_0.6Ang/forces
        # create 4.14Ang_0.6Ang/POSITIONs
        # create POSITIONs  ( which combines all positions )
        ###########################################################################
        self.disps = []
        self.disp_vs_force      = np.zeros((len(self.folder_alldisp)*2,2))
        self.disp_vs_force_110  = np.zeros((len(self.folder_alldisp)*2,2))

        self.dist_0_05_05_at0   = np.zeros((len(self.folder_alldisp),3))
        self.dist_0_05_05_at1   = np.zeros((len(self.folder_alldisp),3))

        self.force_0_05_05      = np.zeros((len(self.folder_alldisp),3))
        self.force_05_45_0      = np.zeros((len(self.folder_alldisp),3))
        self.force_1_1_1        = np.zeros((len(self.folder_alldisp),3))
        self.force_1_1_0        = np.zeros((len(self.folder_alldisp),3))
        self.force_45_45_0      = np.zeros((len(self.folder_alldisp),3))
        self.force_0_05_05_rest = np.zeros((len(self.folder_alldisp),3))
        self.force_05_05_0_rest = np.zeros((len(self.folder_alldisp),3))

        da = len(self.folder_alldisp)
        if self.verbose:
            print()
            for idx,i in enumerate(self.folder_alldisp):
                if idx == 0 and self.sc == False:
                    if os.path.isfile(self.folder_alldisp[0]+"/POSITIONs"):
                        nat, cryst, self.sc, self.alat, cell, pos0 = my.try_to_get_from_positions_only__nat_cryst_sc_cell(1,1,or_folder_to_positions=self.folder_alldisp[0]+"/POSITIONs")
                        sc = self.sc
                        alat = self.alat
                        if self.verbose:
                            print('xX nat',nat)
                            print('xX cryst',cryst)
                            print('xX self.sc',self.sc)
                            print('xX self.alat',self.alat)
                            print('xX cell',cell)

                        #positions_forces = np.loadtxt(self.folder_alldisp[0]+"/POSITIONs")
                        #positions = positions_forces[:,[0,1,2]]
                        #forces = positions_forces[:,[3,4,5]]
                        #nat, cryst, self.sc, self.alat, cell, pos0 = my.try_to_get_from_positions_only__nat_cryst_sc_cell(positions,forces,verbose=False)

            print()

        for idx,i in enumerate(self.folder_alldisp):
            if idx == 0 and self.verbose:
                print("folder_alldisp[i] (W):",i)
            disp = float(i.split("Ang_")[-1].split("/")[0])
            self.disps.append(disp)
            pos,forces,cell = my.folder_get_pos_forces_cell(i)
            sc, alat = my.try_to_get_alat_and_sc_from_cell_and_positions(cell,pos)

            ### get self.sc in case it is necessary
            if self.sc == False:
                self.sc = sc
            else:
                if self.sc != sc:
                    print('sc',sc)
                    print('self.sc',self.sc)
                    sys.exit('sc not self.sc (77); Exit')

            if self.alat == False:
                self.alat = alat
            else:
                if self.alat != alat:
                    print('alat',alat)
                    print('self.alat',self.alat)
                    sys.exit('alat not self.alat (77); Exit')

            self.nndist = self.alat/np.sqrt(2.)

            if self.verbose > 1:
                print('self.sc  ',self.sc)
                print('self.alat',self.alat)

            if idx == 0 and self.verbose == True:
                print("self.dofor        (3):",self.dofor)
                print("self.element      (3):",self.element)
                print("self.rmin_distmax (3):",self.rmin_distmax)
                print("self.rmax_distmax (3):",self.rmax_distmax)
                print("self.alatT0K      (3):",self.alatT0K)
                print("self.alatTmelt    (3):",self.alatTmelt)
                print("self.sc           (3):",self.sc)
                print("self.alat         (3):",self.alat)

            #print('100',forces)
            #print('pos')
            #print(pos)


            if idx == 0:
                self.force_all          = np.zeros((len(self.folder_alldisp),4*self.sc**3,3))
                self.pos_all            = np.zeros((len(self.folder_alldisp),4*self.sc**3,3))

            self.force_all[idx] = forces
            pos_ = my.center_positions_around_0(pos,self.sc,self.alat)
            self.pos_all[idx] = pos_
            #for idx,i in enumerate(pos_):
            #    print('idx',idx,pos[idx],pos_[idx])

            if disp == 0.0:
                atoms = Atoms(self.element+str(len(pos_)),pos_,cell=cell,pbc=[1,1,1])
                #for idx,i in enumerate(atoms.positions):
                #    print('idx',idx,atoms.get_chemical_symbols()[idx],atoms.positions[idx])
                #print(atoms.cell)
                #print(atoms.pbc)
                NN1,NN2,NN3 = my.get_NN1_NN2_NN3_from_positions(atoms,alat,atomnr=0,verbose=args.verbose)


            ##########################################################################################
            # This is so for for a random displacement, whichever was loaded last
            ##########################################################################################
            def dosearch(search,str,pos,selfalat):
                try:
                    #idx_45_45_0 = getindex(pos/self.alat,np.array([round(self.sc-0.5,12),round(self.sc-0.5,12),0]))
                    idx_45_45_0 = getindex(pos/selfalat,np.array(search))
                    #print('idx_-05_-05_0 (the index ot the atom at -05_-05_0)',idx_45_45_0)
                except IndexError:
                    print("ERROR, alat for "+str,selfalat)
                    print("search:",np.array(search)) #round(a,14)
                    print("pos",pos/selfalat)
                    sys.exit()
                return idx_45_45_0

            idx_05_05_0 = dosearch([0.5,0.5,0],"idx_05_05_0",pos,self.alat)
            p_05_05_0   = pos[idx_05_05_0]
            f_05_05_0   = forces[idx_05_05_0]

            idx_1_1_0   = dosearch([1.0,1.0,0],"idx_1_1_0",pos,self.alat)
            p_1_1_0     = pos[idx_1_1_0]
            f_1_1_0     = forces[idx_1_1_0]
            self.force_1_1_0[idx] = f_1_1_0

            idx_45_45_0 = dosearch([round(self.sc-0.5,12),round(self.sc-0.5,12),0],"idx_45_45_0",pos,self.alat)
            p_45_45_0   =    pos[idx_45_45_0]
            f_45_45_0   = forces[idx_45_45_0]
            self.force_45_45_0[idx] = f_45_45_0

            idx_40_40_0 = dosearch([round(self.sc-1.0,12),round(self.sc-1.0,12),0],"idx_40_40_0",pos,self.alat)
            p_40_40_0   =    pos[idx_40_40_0]
            f_40_40_0   = forces[idx_40_40_0]

            idx_45_0_45 = dosearch([round(self.sc-0.5,12),0,round(self.sc-0.5,12)],"idx_45_0_45",pos,self.alat)
            p_45_0_45   =    pos[idx_45_0_45]
            f_45_0_45   = forces[idx_45_0_45]

            idx_0_45_45 = dosearch([0, round(self.sc-0.5,12),round(self.sc-0.5,12)],"idx_0_45_45",pos,self.alat)
            p_0_45_45   =    pos[idx_0_45_45]
            f_0_45_45   = forces[idx_0_45_45]

            idx_05_45_0 = dosearch([0.5,round(self.sc-0.5,12),0],"idx_05_45_0",pos,self.alat)
            p_05_45_0   =    pos[idx_05_45_0]
            f_05_45_0   = forces[idx_05_45_0]
            self.force_05_45_0[idx] = f_05_45_0

            idx_45_05_0 = dosearch([round(self.sc-0.5,12),0.5,0],"idx_45_05_0",pos,self.alat)
            p_45_05_0   =    pos[idx_45_05_0]
            f_45_05_0   = forces[idx_45_05_0]

            if self.verbose > 1:
                print("idx_45_45_0",idx_45_45_0)
                print('idx',disp,'repulsive force:',f_05_05_0,'attractive force:',f_45_45_0)

            idx_0_05_05 = getindex(pos/self.alat,np.array([0,0.5,0.5]))
            p_0_05_05 = pos[idx_0_05_05]
            f_0_05_05 = forces[idx_0_05_05]
            self.force_0_05_05[idx] = f_0_05_05


            d_0_05_05_at0 = p_0_05_05-pos[0]
            d_0_05_05_at1 = p_0_05_05 - p_05_05_0

            self.dist_0_05_05_at0[idx] = d_0_05_05_at0
            self.dist_0_05_05_at1[idx] = d_0_05_05_at1
            if self.verbose2:
                print('idx',disp,'f_0_05_05:',f_0_05_05,'attractive force:',d_0_05_05_at0, LA.norm(d_0_05_05_at0))




            if disp == 0.7:
                if self.verbose:
		    print('doing this for disp = ',disp)
                self.pos = pos
                self.forces = forces

                self.idx_05_05_0 = idx_05_05_0
                self.idx_45_45_0 = idx_45_45_0
                self.idx_0_05_05 = idx_0_05_05

                self.f_0_05_05__01 = f_0_05_05

                self.vec_to_0_05_05__01 = p_0_05_05 - pos[0]
                self.dist_to_0_05_05__01 = LA.norm(p_0_05_05 - pos[0])
                self.forces__01 = forces
                self.pos__01 = pos
		self.f_110 = LA.norm(f_1_1_0)*np.sign(f_1_1_0[0])

                if self.verbose:
                    print('force on 05_05_0:',disp,f_05_05_0,LA.norm(f_05_05_0))
                    print('force on 1_1_0  :',disp,f_1_1_0,LA.norm(f_1_1_0)*np.sign(f_1_1_0[0]))
                    print('force on 0_05_05:',disp,LA.norm(f_0_05_05),f_0_05_05)
                self.ratio0510  = (LA.norm(f_1_1_0)*np.sign(f_1_1_0[0]))/LA.norm(f_05_05_0)
                self.ratio0510 = round(self.ratio0510,3)
		add = "ATTRACTIVE_to_"
                if self.verbose:
                    print('ratio f[0.5,0.5,0]/f[1,1,0]',self.ratio0510)
                #np.sign(ele.forces[0][0])

            if disp < 0.5:
                #
                #print('i',disp,LA.norm(p_05_05_0-pos[0]))
                #print('i',dist_att_rest-p_45_45_0,p_05_05_0,LA.norm(p_45_45_0-pos[0]))
                #print('i',disp, pos[0],p_0_05_05, p_0_05_05 - pos[0],self.f_0_05_05,LA.norm(p_0_05_05 - pos[0]))
                pass


            ##############################################################################
            # put together disp vs forces at atom [0.5,0.5,0
            ##############################################################################
            self.disp_vs_force[idx,0]       = LA.norm(p_05_05_0-pos[0]) # repulsive dist
            self.disp_vs_force[idx,1]       = -np.sign(f_05_05_0[0])*LA.norm(f_05_05_0) # repulsive forces
            self.disp_vs_force[idx+da,0]    = self.nndist + disp # attractive dist
            self.disp_vs_force[idx+da,1]    = np.sign(f_45_45_0[0])*LA.norm(f_45_45_0) # attractive forces

            ##############################################################################
            # put together disp vs forces at atom 1_1_0 to substract from force @ 05_05_0 as a correction
            ##############################################################################
            self.disp_vs_force_110[idx,0]    = LA.norm(p_05_05_0-pos[0])  # distance
            self.disp_vs_force_110[idx,1]    = -np.sign(f_1_1_0[0])*LA.norm(f_1_1_0)           # forces repulsive ... in this case the sign is different
            self.disp_vs_force_110[idx+da,0] = self.nndist + disp
            self.disp_vs_force_110[idx+da,1] = np.sign(f_40_40_0[0])*LA.norm(f_40_40_0)

            ##############################################################################
            # put together disp vs forces at atom 05_45_0
            ##############################################################################




    #######################################################################################
    #######################################################################################
    # here all folders are defined
    #######################################################################################
    #######################################################################################
        self.disp_vs_force = restrict_xrange_of_array(array = self.disp_vs_force,
                         xmin = self.nndist - self.rmin_distmax,
                         xmax = self.nndist + self.rmax_distmax)
        self.disp_vs_force_110 = restrict_xrange_of_array(array = self.disp_vs_force_110,
                         xmin = self.nndist - self.rmin_distmax,
                         xmax = self.nndist + self.rmax_distmax)

        ##############################################################################
        # remove duplicates
        ##############################################################################
        if self.verbose > 2:
            print(self.disp_vs_force)
        self.disp_vs_force                      = my.remove_duplicates_in_numpy_xy_array_and_sort(self.disp_vs_force,roundto=10)
        self.disp_vs_force_110                  = my.remove_duplicates_in_numpy_xy_array_and_sort(self.disp_vs_force_110,roundto=10)
        self.disp_vs_force_corrected_for_110    = np.copy(self.disp_vs_force)
        self.disp_vs_force_corrected_for_110[:,1] += self.disp_vs_force_110[:,1]
        if self.verbose > 2:
            print('---------------------------')
            print('self.disp_vs_force')
            print('---------------------------')
            print(self.disp_vs_force)
            print('---------------------------')
        #print('??')
        #print('corr')
        #print(self.disp_vs_force_corrected_for_110)
        #sys.exit()


        #print("#########################################################################")
        #print("# fits on [0.5,0.5,0.0]")
        #print("#########################################################################")
        fit     = pot_parametrize.fit_to_func(self.disp_vs_force,function='morse',fixzeroat=self.nndist)
        self.disp_vs_force_rel = np.copy(self.disp_vs_force)
        self.disp_vs_force_rel[:,0] = self.disp_vs_force_rel[:,0]/self.alat
        #print('disp_vs_force_rel')
        #print(disp_vs_force_rel)
        #print()
        #print(self.nndist/self.alat)
        fit_rel           = pot_parametrize.fit_to_func(self.disp_vs_force_rel,function='morse',fixzeroat=self.nndist/self.alat)
        fit_rel_repulsive = pot_parametrize.fit_to_func(self.disp_vs_force_rel,function='ma',fixzeroat=self.nndist/self.alat)
        try:
            fitmc1 = pot_parametrize.fit_to_func(self.disp_vs_force,function='mc1',fixzeroat=self.nndist)
        except TypeError:
            fitmc1 = False
        fit_corrected_for_110 = pot_parametrize.fit_to_func(self.disp_vs_force_corrected_for_110,function='morse',fixzeroat=self.nndist)
        print('# fits on [0.5,0.5,0.0]: fit.parameters                   morse',fit.parameters)
        print('# fits on [0.5,0.5,0.0]: fit_rel.parameters               morse',fit_rel.parameters)
        print('# fits on [0.5,0.5,0.0]: fit_rel_repulsive.parameters     morse',fit_rel_repulsive.parameters)
        print('# fits on [0.5,0.5,0.0]: fit_corrected_for_110.parameters morse',fit_corrected_for_110.parameters)
	#for i in fit.parameters:
	#    print(i)
	#print()
	#for i in fit_rel.parameters:
	#    print(i)
        #sys.exit('111111111112233')
        fit_on_05_05_0 = np.zeros((len(fit.fit),7))
        fit_on_05_05_0[:,0] = fit.fit[:,0]
        dist__hydrogen   = (1./7.2)*fit.fit[:,0]*10**-10
        dist__           = (1./1.0)*fit.fit[:,0]*10**-10
        fit_on_05_05_0[:,1] = self.disp_vs_force[:,1]   # actual (VASP) forces
        fit_on_05_05_0[:,2] = fit.fit[:,1]              # fitted morse forces
        #print('dist__',dist__)
        k__ = (0.9*10**10)
        e__ = (1.6*10**(-19))
        dist__meter = dist__*10**(-10)   # angstrom to meter
        newton_to_mev_per_angstrom = 6.2415091**11
        q1__ = ((1.6*78.)*10**(-19)) #78
        q1__ = ((1.6*2.8)*10**(-19)) #78
        fit_on_05_05_0[:,6] =  (-k__*((q1__**2) /(dist__**2))*newton_to_mev_per_angstrom)+12.              # columb forces

        # get screened potential (https://en.wikipedia.org/wiki/Electric-field_screening)
        #for idx,i in enumerate(dist__):
        #    print('d',str(idx).ljust(3),str(dist__[idx]).ljust(25),'force VASP',str(fit_on_05_05_0[idx,1]).ljust(30),"columb",fit_on_05_05_0[idx,6])
        if type(fitmc1) != bool:
            fit_on_05_05_0[:,3] = fitmc1.fit[:,1]
        fit_on_05_05_0[:,4] = fit_on_05_05_0[:,1] - fit.fit[:,1]
        if type(fitmc1) != bool:
            fit_on_05_05_0[:,5] = fit_on_05_05_0[:,1] - fitmc1.fit[:,1]

        if type(fitmc1) != bool:
            np.savetxt(self.dofor+'/disp_fit.parameters_mc1.dat',fitmc1.parameters)
        np.savetxt(self.dofor+'/disp_fit.parameters_morse.dat',fit.parameters)
        np.savetxt(self.dofor+'/disp_fit_corr_110.parameters_morse.dat',fit_corrected_for_110.parameters)
        np.savetxt(self.dofor+'/disp_vs_forces.dat',self.disp_vs_force)
        np.savetxt(self.dofor+'/disp_vs_forces_rel.dat',self.disp_vs_force_rel)
        print('saved',self.dofor+'/disp_vs_forces.dat')
        print('saved',self.dofor+'/disp_vs_forces_rel.dat')
        print('saved',self.dofor+'/disp_fit.parameters_morse.dat')
        if type(fitmc1) != bool:
            np.savetxt(self.dofor+'/disp_vs_fitted_mc1.dat',fitmc1.fit)
        np.savetxt(self.dofor+'/disp_vs_fitted_morse.dat',fit.fit)
        print('saved',self.dofor+'/disp_vs_fitted_morse.dat')
        np.savetxt(self.dofor+'/disp_vs_fitted_morse_rel.dat',fit_rel.fit)
        print('saved',self.dofor+'/disp_vs_fitted_morse_rel.dat')
        np.savetxt(self.dofor+'/disp_vs_fitted_morse_rel_repulsive.dat',fit_rel_repulsive.fit)
        print('saved',self.dofor+'/disp_vs_fitted_morse_rel_repulsive.dat')
        f1 = hesse.Morse_repulsive_derivative_to_normaldata(fit_rel_repulsive.datax, *fit_rel_repulsive.parameters)
        f2 = hesse.Morse_repulsive_derivative(fit_rel_repulsive.datax, *fit_rel_repulsive.parameters)
        np.savetxt(self.dofor+'/disp_vs_f1.dat',np.array([fit_rel_repulsive.datax,f1]).T)
        np.savetxt(self.dofor+'/disp_vs_f2.dat',np.array([fit_rel_repulsive.datax,f2]).T)

        ##############################################################################
        # shift to other alat
        ##############################################################################
        if self.shift_parametrization_to_alat:
            nndist_shifted = self.shift_parametrization_to_alat/np.sqrt(2.)
            if self.verbose:
                print('---> old alat (self.nndist     ):',self.nndist)
                print('--->',self.disp_vs_force)
                print('---> new alat (nndist_shifted):',nndist_shifted)

            force_at_nndist_otheralat = hesse.Morse_derivative(nndist_shifted, *fit.parameters)
            if self.verbose:
                print('---> force_at_nndist_otheralat:',force_at_nndist_otheralat)
            disp_vs_force_yshift = self.disp_vs_force.copy()  # shifted to simulate a distance @ lattice constants @Tmelt
            disp_vs_force_yshift[:,1] = disp_vs_force_yshift[:,1] - force_at_nndist_otheralat # shifted to simulate a distance @ lattice constants @Tmelt
            if self.verbose:
                print('--->',disp_vs_force_yshift)
            fit_shifted = pot_parametrize.fit_to_func(disp_vs_force_yshift,function='morse',fixzeroat=nndist_shifted)

            filenameout = self.dofor+'/disp_vs_forces_shifted_to_alat_'+str(self.shift_parametrization_to_alat)+'.dat'
            np.savetxt(filenameout,disp_vs_force_yshift)

        addto = self.dofor+'/disp_vs_forces.dat'
        with open(addto, "a") as myfile:
                myfile.write("\n")
                myfile.write(str(self.disp_vs_force[:,0].min())+" 0\n")  # horizontal line
                myfile.write(str(self.disp_vs_force[:,0].max())+" 0\n")  # horizontal line
                myfile.write("\n")
                myfile.write(str(self.nndist)+" "+str(self.disp_vs_force[:,1].max())+"\n")  # vertical line, should meet morse at zero force
                myfile.write(str(self.nndist)+" "+str(-self.disp_vs_force[:,1].max())+"\n") # vertical line, should meet morse at zero force
                myfile.write("\n")
                myfile.write(str(self.nndistTmelt)+" "+str(self.disp_vs_force[:,1].max())+"\n") # vertical line @nndist @Tmelt
                myfile.write(str(self.nndistTmelt)+" "+str(-self.disp_vs_force[:,1].max())+"\n") # vertical line @nndist @Tmelt
                if self.shift_parametrization_to_alat:
                    myfile.write("\n")
                    myfile.write(str(nndist_shifted - 0.1)+" "+str(force_at_nndist_otheralat)+"\n") # vertical line @nndist @Tmelt
                    myfile.write(str(nndist_shifted + 0.1)+" "+str(force_at_nndist_otheralat)+"\n") # vertical line @nndist @Tmelt

        if self.shift_parametrization_to_alat != False:
            #print('fs',fit_shifted.parameters)
            #sys.exit()
            np.savetxt(parametrizationfile_morse,fit_shifted.parameters)
            if self.only_return_parametrization_file == True:
                return parametrizationfile_morse

        if True:
            #print("############################################################")
            #print("# parametrize Michaels Model ueber alle nachbarn (a,b,c,d)")
            #print("############################################################")
            #############################################################
            # Al (as a rule of thumb: ene std/ 2 = error in free energy)
            #############################################################
            # --- 7 displacements -----------------------------------------------
            # displacements_/Al/2x2x2sc_4.13Ang_quer_3x3x3kp/shmall_subset/POSITIONs (== parametrization)  ( the convergence criteria seem fine, not sure why convergence so bad)
            # LA no tox       7 displacements (                    ) (amc -e Al -v -pm 7_poly           )       ps: 0.99618  ene std: 2.2041 meV/atom  for std: 0.22204
            # LA no tox       7 displacements (                    ) (amc -e Al -v -pm 7_morse          )       ps: 0.99772  ene std: 1.3189 meV/atom  for std: 0.33527

            # displacements_/Al/3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400  (==parametrization)
            # LA no tox       7 displacements (orig from disp foler) (amc -e Al -v -pm 7_morse          )       ps: 0.99835  ene std: 0.2973 meV/atom  for std: 0.01297
            # LA no tox       7 displacements                        (amc -e Al -v -pm 7_morseprl -t1o 0)       pc: 0.99723  ene std: 0.4165 meV/atom  for std: 0.01604
            # LA with tox     7 displacements (orig from disp foler) (amc -e Al -v -pm 7_morse -t1 -0.065)      ps: 0.99927  ene std: 0.0164 meV/atom  for std: 0.00619 a_mor:1.5256625969 D_mor:0.213311082
            # LA with tox     7 displacements                        (amc -e Al -v -pm 7_morseprl       )       ps: 0.99875  ene std: 0.0883 meV/atom  for std: 0.00810
	    # POLY            7 displacements                        (amc -e Al -v -pm 7_poly   ) (cd 0.88      ps: 0.99935  ene std: 0.0190 meV/atom  for std: 0.00588  ** best (=POLY)

            # --- SUM_run1/first_20000 ------------------------------------------
            # LA no   tox     SUM_run1/first_20000 displacements     (amc -e Al -v -pm 20000_morse)             ps: 0.99061  ene std: 3.7131 meV/atom  for std: 0.11227
            # LA no   tox     SUM_run1/first_20000 displacements     (amc -e Al -v -pm 20000_morseprl -t1o 0)   ps: 0.98824  ene std: 4.5969 meV/atom  for std: 0.14161
            # LA with tox     SUM_run1/first_20000 displacements     (amc -e Al -v -pm 20000_morse -t1 -0.065)  ps: 0.99616  ene std: 0.8286 meV/atom  for std: 0.05077
            # LA with tox     SUM_run1/first_20000 displacements     (amc -e Al -v -pm 20000_morseprl)          ps: 0.99539  ene std: 0.8035 meV/atom  for std: 0.06043
            # POLY            SUM_run1/first_20000 displacements     (amc -e Al -v -pm 20000_poly ) (  cd)      ps: 0.99638  ene std: 0.7584 meV/atom  for std: 0.04931  ** best in std_ene (=POLY)
            # POLY            SUM_run1/first_20000 displacements     (amc -e Al -v -pm 20000_poly ) (abcd)      ps: 0.99632  ene std: 0.7897 meV/atom  for std: 0.04978

            # --- 2000_PTS_dosall -#####--SAME      LATTICE CONSTANT ------------
            # POLY             7 displacements             (amc -e Al -pm 7_poly for 4.13 MD                    ps: 0.99855  ene std: 0.0365
            # POLY            2000 0__PTS                  (amc -e Al -pm 7_poly --f_md_2x2x2_30                ps: 0.99637  ene std: 2.0677

            # --- 2000_PTS_dosall -#####--DIFFERENT LATTICE CONSTANT ------------
            # LA no   tox     2000 0__PTS_dosall_ ...      (amc -e Al --f_md_2x2x2_30)                          ps: 0.98858  ene std: 8.8278 meV/atom  for std: 0.15698
            # LA with tox     2000 0__PTS_dosall_ ...      (amc -e Al --f_md_2x2x2_30 -t1 -0.072)               ps: 0.99679  ene std: 2.3736 meV/atom  for std: 0.06278
            # POLY            2000 0__PTS_dosall_ ...      (amc -e Al --f_md_2x2x2_30 -pm prl2015_polycorralat) ps: 0.99634  ene std: 2.0745 meV/atom  for std: 0.05248
            # POLY            2000 0__PTS_dosall_ ...      (amc -e Al --f_md_2x2x2_30 -pm 20000_poly)           ps: 0.99637  ene std: 2.1900 meV/atom  for std: 0.06278  ** best in std_ene (=POLY)
            #

            #############################################################
            # Pt (as a rule of thumb: ene std/ 2 = error in free energy)
            #############################################################
            # in /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Pt/3x3x3sc_4.1Ang_quer_2x2x2kp_vasp4:
            # LA no tox       7 displacements              (amc -e Pt -f .)                                     ps: 0.99755  ene std: 0.6370  meV/atom  for std: 0.03135 (a_mor=1.849;D_mor=0.25)
            # LA with tox     7 displacements              (amc -e Pt -f . -t1 -0.135)                          ps: 0.99884  ene std: 0.0186  meV/atom  for std: 0.01746
            #                                              (== amc -e Pt -f . -t1 -0.135 -alat_mor 4.1 -a_mor 1.8496052614 -D_mor 0.2504607705)
            # LA with tox     7 displacements   (amc -e Pt -f . -t1 -0.135 -alat_mor 4.1 -a_mor 1.8496 -D_mor 0.25)
            #                                                                                                   ps: 0.99883  ene std: 0.0186  meV/atom  for std: 0.01754
            # POLY            7 displacements         (amc -e Pt -pm 7_poly)    (cd 0.88)                       ps: 0.99825  ene std: 0.4035  meV/atom  for std: 0.02261
            # POLY            7 displacements         (amc -e Pt -pm 7_poly)    (abcd 0.88)                     ps: 0.99870  ene std: 0.3360  meV/atom  for std: 0.02021
            # POLY            7 displacements         (amc -e Pt -pm 7_poly -t1o 0.05)    (abcd 0.88)           ps: 0.99879  ene std: 0.1132  meV/atom  for std: 0.01788
            # POLY            7 displacements         (amc -e Pt -pm 7_poly -t1o 0.041)   (abcd 0.88)           ps: 0.99881  ene std: 0.1512  meV/atom  for std: 0.01788
            # POLY            7 displacements         (amc -e Pt -pm 7_poly -t1o -0.06)   (abcd) (rcut=0.84)    ps: 0.99892  ene std: 0.0410  meV/atom  for std: 0.01704  ** best (=POLY+TOX)  -> cant reproduce t1o seems wrong



            # LA with tox     2000 0__PTS_dosall_ ..(amc -e Pt --f_md_2x2x2_30 -t1 -0.135 -alat_mor 4.1 -a_mor 1.8496 -D_mor 0.25)
            #                                                                                                   ps: 0.98739  ene std: 9.8001  meV/atom  for std: 0.20689
            # LA no   tox     2000 0__PTS_dosall_ ..(amc -e Pt --f_md_2x2x2_30)                                 ps: 0.97421  ene std: 31.7570 meV/atom  for std: 0.42125
            # POLY            2000 0__PTS_dosall ...(amc -e Pt -pm 7_poly --f_md_2x2x2_30) (params_cd rcut0.88) ps: 0.98477  ene std: 15.3497 meV/atom  for std: 0.29601
            # POLY            2000 0__PTS_dosall ...(amc -e Pt -pm 7_poly --f_md_2x2x2_30) (params_cd rcut0.84) ps: 0.98582  ene std: 14.6261 meV/atom  for std: 0.27421
            # POLY            2000 0__PTS_dosall ...(amc -e Pt -pm 7_poly --f_md_2x2x2_30) (params_abcd)        ps: 0.98584  ene std: 13.5977 meV/atom  for std: 0.21300
            # POLY            2000 0__PTS_dosall ...(amc -e Pt -pm 7_poly --f_md_2x2x2_30 -t1o 0.041) (par_abcd)ps: 0.98671  ene std: 9.6987  meV/atom  for std: 0.22360  ** best (=POLY+TOX)
            # POLY            2000 0__PTS_dosall ...(amc -e Pt -pm 7_poly --f_md_2x2x2_30 -t1o 0.05 ) (par_abcd)ps: 0.98644  ene std: 9.8821  meV/atom  for std: 0.23063
            # POLY            2000 0__PTS_dosall ...(amc -e Pt -pm 7_poly --f_md_2x2x2_30 -t1o -0.06) (rcut0.84)ps: 0.98868  ene std: 8.2412  meV/atom  for std: 0.21188  ** best (=POLY+TOX)

            #############################################################
            # Ir (as a rule of thumb: ene std/ 2 = error in free energy)
            #############################################################
            # POLY         2000 displacements  (amc  Ir -pm 7_poly --f_md_2x2x2_30 -t1o 0.08)(  cd; rcut=084) ps: 0.98372  ene std:  10.0525 meV/atom  for std:
            # POLY            7 displacements  (amc  Ir -pm 7_poly                          )(  cd; rcut=084) ps: 0.98817  ene std:  10.0525 meV/atom  for std:

            # POLY            7 displacements  (amc  Ir -pm 7_poly)                  (abcd; rcut=088)         ps: 0.98512  ene std:  2.4026  meV/atom  for std: 0.11722
            # POLY            7 displacements  (amc  Ir -pm 7_poly)                  (  cd; rcut=088)         ps: 0.98887  ene std:  2.1500  meV/atom  for std: 0.09742
            # POLY            7 displacements  (amc  Ir -pm 7_poly)                  (  cd; rcut=084)         ps: 0.98817  ene std:  2.2560  meV/atom  for std: 0.09952
            # POLY            7 displacements  (amc  Ir -pm 7_poly -t1o 0.26)        (  cd; rcut=088)         ps: 0.99255  ene std:  0.9503  meV/atom  for std: 0.06706
            # LA              7 displacements  (amc  Ir -pm 7_morse         )        (              )         ps: 0.98480  ene std:  2.2345  meV/atom  for std: 0.09378
            # LA + tox        7 displacements  (amc  Ir -pm 7_morse -t1o 0.3)        (              )         ps: 0.98914  ene std:  0.8450  meV/atom  for std: 0.07454
            #
            # LA + tox     2000 displacements  (amc  Ir -pm 7_morse -t1o 0.3 --f_md_2x2x2_30                  ps: 0.96823  ene std: 22.9657  meV/atom  for std: 0.46523
            # POLY + tox   2000 displacements  (amc  Ir -pm 7_poly  -t1o 0.26 --f_md_2x2x2_30                 ps: 0.97609  ene std: 17.8781  meV/atom  for std: 0.40539
            # POLY         2000 displacements  (amc  Ir -pm 7_poly            --f_md_2x2x2_30                 ps: 0.98246  ene std: 13.4180  meV/atom  for std: 0.34881 ** best (=POLY WITHOUT TOX)
            # POLY + tox   2000 displacements  (amc  Ir -pm 7_poly  -t1o 0.07 --f_md_2x2x2_30                 ps: 0.98389  ene std:  9.7953  meV/atom  for std: 0.31164

            print('self.dofor',self.dofor)
            params_abcd,params_cd = my.get_michaels_paramerization(pos_all=self.pos_all,force_all=self.force_all,NN1=NN1,alat=alat,atoms=atoms,rcut=0.88,save_parametrization=self.dofor)
            print("# parametrize Michaels Model ueber alle nachbarn (a,b,c,d)",params_abcd)
            print("# parametrize Michaels Model ueber alle nachbarn (    c,d)",params_cd)
            params_abcd,params_cd = my.get_michaels_paramerization(pos_all=self.pos_all,force_all=self.force_all,NN1=NN1,alat=alat,atoms=atoms,rcut=0.84,save_parametrization=self.dofor)
            print("# parametrize Michaels Model ueber alle nachbarn (a,b,c,d)",params_abcd)
            print("# parametrize Michaels Model ueber alle nachbarn (    c,d)",params_cd)
            params_abcd,params_cd = my.get_michaels_paramerization(pos_all=self.pos_all,force_all=self.force_all,NN1=NN1,alat=alat,atoms=atoms,rcut=0.85355,save_parametrization=self.dofor)
            print("# parametrize Michaels Model ueber alle nachbarn (a,b,c,d)",params_abcd)
            print("# parametrize Michaels Model ueber alle nachbarn (    c,d)",params_cd)
            sys.exit('778866')

        ##############################################################################
        # force on 0_05_05 (previously tox)
        # Weight 1Morse 05_05_0
        # Weight 4Morse 0_05_05
        ##############################################################################
        print("#########################################################################")
	print("# xx@Forces on [0.5,0.5,0.0]")
        print("#########################################################################")
	print("  dist  VASP   morse  mc1    dmorse  dmc1  Columb")
        for i in fit_on_05_05_0:
            print(i)
        print("#########################################################################")
        print("# fits on [0.0,0.5,0.5] previously tox (params_morse_normal)")
        print("#########################################################################")
        if True:
            print('dist pos[0] to pos_0_05_05 ')
            params_morse_normal     = fit.parameters
            print('params_morse_normal',params_morse_normal)
            for idx,i in enumerate(self.dist_0_05_05_at0):
                dist_norm                           = np.around(LA.norm(i),5)
                fm                                  = np.round(hesse.Morse_derivative(LA.norm(i), *params_morse_normal),4)
                forces_morse                        = hesse.getefvec(i,params_morse_normal,pot = 'm')
                self.force_0_05_05_rest[idx]        = self.force_0_05_05[idx] + forces_morse[1]
                vasp_force                          = self.force_0_05_05[idx] #.ljust(22)
                vasp_force_norm                     = str(np.around(-LA.norm(abs(vasp_force)),decimals=3)).ljust(6)
                vasp_minus_morse                    = vasp_force+forces_morse[1]
                print(str(i).ljust(22),str(dist_norm).ljust(7),'F_0_05_05_vasp',str(vasp_force).ljust(22),"norm:",vasp_force_norm,'F_morse',str(-forces_morse[1]).ljust(22),'norm',str(fm).ljust(7),'vasp-morse',vasp_minus_morse,'norm',str(np.around(LA.norm(vasp_minus_morse),3)).ljust(7))
            #print(self.dist_0_05_05_at0)
            print()
            print('force on 0_05_05       remaining forces on 0_05_05          this should be a sum of')
            print('                       (after substr. morse)                a*vec at_000 + b*vec at_05050')
            for idx,i in enumerate(self.dist_0_05_05_at0):
                print(self.force_0_05_05[idx],self.force_0_05_05_rest[idx],self.dist_0_05_05_at0[idx],self.dist_0_05_05_at1[idx])
            print()
            for idx,i in enumerate(self.dist_0_05_05_at0):
                print(self.force_0_05_05[idx],abs(self.force_0_05_05_rest[idx]),abs(self.force_0_05_05_rest[idx]).sum()*4.)

            from lmfit import Model,minimize
            print('alat',self.alat,self.alat/np.sqrt(2.))
            morsemodel = Model(hesse.Morse_derivative)
            morsemodel.set_param_hint('re' ,value=self.alat/np.sqrt(2.), vary=False)
            morsemodel.set_param_hint('De' ,value=0.25, min=0.01, max=0.5)
            morsemodel.set_param_hint('aa' ,value=1.5,  min=1.0,  max=3.0)


            r1 = self.disp_vs_force[:da,0][::-1]
            y1 = self.disp_vs_force[:da,1][::-1]
            #for idx,i in enumerate(r1):print(r1[idx],y1[idx])
            result = morsemodel.fit(y1, r=r1)
            print('vgl (obtained by fit)',params_morse_normal)
            print('vgl (only repulsive )',np.round([result.best_values.get('De'),result.best_values.get('aa'),result.best_values.get('re')],3))

            r1a = self.disp_vs_force[(da-1):,0]
            y1a = self.disp_vs_force[(da-1):,1]

            rall = np.concatenate((r1, r1a), axis=None)
            yall = np.concatenate((y1, y1a), axis=None)
            #for idx,i in enumerate(rall):print(rall[idx],yall[idx])
            result = morsemodel.fit(yall, r=rall)
            print('vgl (all forces     )',np.round([result.best_values.get('De'),result.best_values.get('aa'),result.best_values.get('re')],3))

            raddx =  LA.norm(self.dist_0_05_05_at0,axis=1)
            raddy = -LA.norm(abs(self.force_0_05_05),axis=1)
            np.savetxt(self.dofor+'/disp_vs_forces_weighted.dat',np.transpose([raddx,raddy]))
            for i in [1,2,3,4]:
                rall = np.concatenate((rall, raddx), axis=None)
                yall = np.concatenate((yall, raddy), axis=None)
                #for idx,i in enumerate(rall):print(rall[idx],yall[idx])
                result = morsemodel.fit(yall, r=rall)
                print('vgl (all forces '+str(i)+'   )',np.round([result.best_values.get('De'),result.best_values.get('aa'),result.best_values.get('re')],3))
            params_weighted = [result.best_values.get('De'),result.best_values.get('aa'),result.best_values.get('re')]
            np.savetxt(self.dofor+'/disp_fit.parameters_morse_weighted.dat',params_weighted)
            print('saved in',self.dofor)
            print("STILL DID NOT INCLUDE THE ATTRACTIVE PART TIMES 4")
            for idx,i in enumerate(self.dist_0_05_05_at0):
                fm = hesse.Morse_derivative(LA.norm(i), *params_weighted)
                fv = hesse.getefvec(i,params_morse_normal,pot = 'm')
                self.force_0_05_05_rest[idx] = self.force_0_05_05[idx] + fv[1]
                print(i,LA.norm(i),'force_0_05_05 vasp_full',self.force_0_05_05[idx],"norm",-LA.norm(abs(self.force_0_05_05[idx])),'force morse',fm,'fv',-fv[1]) #,LA.norm(fv[1]))
            print('this parametrization, even for T=0K displacements, is not optimal!, try to add attractive forces on 0_45_45')


        if True:
            print("###################################################")
            print("# parametrize Morse from morse_michaels_plot_3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400.txt")
            print("###################################################")
            path = "/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/morse_michaels_plot_3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400.txt"
            xymorse = np.loadtxt(path)
            morsemodel = Model(hesse.Morse_derivative)
            morsemodel_var = [ 'De','aa','re' ]
            morsemodel.set_param_hint('re' ,value=1./np.sqrt(2.), vary=False)
            morsemodel.set_param_hint('De' ,value=0.25)
            morsemodel.set_param_hint('aa' ,value=1.5)
            result = morsemodel.fit(xymorse[:,1], r=xymorse[:,0])
            parameters_morse = []
            for i in morsemodel_var:
                coef = result.best_values.get(i)
                print(i,"Morse_derivative:",coef)
                parameters_morse.append(coef)
            np.savetxt("xy"+add+"_morse.dat",np.array([xymorse[:,0],result.best_fit]).T)
            print('saved',"xy"+add+"_morse.dat")

            # now for unscaled positions
            morsemodel.set_param_hint('re' ,value=self.alat/np.sqrt(2.), vary=False)
            morsemodel.set_param_hint('De' ,value=0.25)
            morsemodel.set_param_hint('aa' ,value=1.5)
            result = morsemodel.fit(xymorse[:,1], r=xymorse[:,0]*self.alat)
            parameters_morse = []
            for i in morsemodel_var:
                coef = result.best_values.get(i)
                print("for positions in angstrom:",i,"Morse_derivative:",coef)
                parameters_morse.append(coef)

            # now for unscaled positions to obtain michaels parametrization
            morsemodel = Model(hesse.Michael_polynomial_for_morsevalues_positive)
            morsemodel_function_var     = [ 'a','b','c','d','e' ]

        if True:
            print("###################################################")
            print("# parametrize michaels plot: Michaels_model_michaels_plot_3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400")
            print("###################################################")
            from lmfit import Model,minimize
            path ="/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/Michaels_model_michaels_plot_3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400_2.txt"
            xy = np.loadtxt(path)
            xmin = xy[:,0].min()
            xmax = xy[:,0].max()
            print('xmin',xmin)
            print('xmax',xmax)
            xred_dense = np.arange(xmin,xmax,0.001)
            #print('xy')
            #print(xy)
            print(hesse.__file__)

            ################## define the fitting function
            for allvar in [True,False]:
                print()
                print()
                mmodel_function     = hesse.Michael_poly_der
                mmodel_function_var     = [ 'a','b','c','d' ]
                mmodel_function_der     = hesse.Michael_poly_der_der
                mmodel_function_var_der = [ 'a','b','c','d' ]
                add = 'allvar____'

                if allvar == False:
                    mmodel_function     = hesse.Michael_poly_der_noA_noB
                    mmodel_function_var     = [ 'c','d' ]
                    mmodel_function_der     = hesse.Michael_poly_der_der_noA_noB
                    mmodel_function_var_der = [ 'c','d' ]
                    add = "allvar_noB"

                #######################
                # hier die kraefte
                #######################
                mmodel = Model(mmodel_function)
                for i in mmodel_function_var: mmodel.set_param_hint(i ,value=1) #, min=0.01, max=0.5)
                result = mmodel.fit(xy[:,1], r=xy[:,0])
                parameters = []
                for i in mmodel_function_var:
                    coef = result.best_values.get(i)
                    print(i,"Forces xx:",coef)
                    parameters.append(coef)
                print('parameters    ',parameters)
                name = "xy"+add+"_michael_from_plot.dat"
                np.savetxt(name,np.array([xy[:,0],result.best_fit]).T)
                print('saved AAA FORCES',name)

                parameters_der = []
                for i in mmodel_function_var_der:
                    coef = result.best_values.get(i)
                    print(i,"derivative:",coef)
                    parameters_der.append(coef)
                print('parameters_der',parameters_der)
                if add == "allvar_noB":
                    parameters_mi = parameters
                    mmodel_key = 'Michael_poly_der_noA_noB'

                np.savetxt(name+"_der.dat",np.array([xred_dense,mmodel_function_der(xred_dense,*parameters_der)]).T)
                print('saved AAA FORCES derivative',name+"_der.dat")

        if False:
            print("###################################################")
            print("# parametrize ueber 7 auslenkungen entlang [110]")
            print("###################################################")
            print(self.disp_vs_force)
            xy = self.disp_vs_force
            xy[:,0] = xy[:,0]/self.alat
            print(self.disp_vs_force)
            sys.exit()


        if False:
            print("#########################################################################")
            print("# fits on [0.5,-0.5,0.0]")
            print("#########################################################################")
            for idx,i in enumerate(self.force_05_45_0):
                dd = self.pos_all[idx,idx_05_45_0]-self.pos_all[idx,0]
                dd = self.pos_all[idx,0]-self.pos_all[idx,idx_05_45_0]
                dist_norm                           = np.around(LA.norm(dd),5)
                fm                                  = np.round(hesse.Morse_derivative(LA.norm(dd), *params_morse_normal),4)
                forces_morse                        = hesse.getefvec(dd,params_morse_normal,pot = 'm')
                print('idx',idx,i,"##""##",i[0],i[1],i[2],"##f",self.force_all[idx,idx_05_45_0],"##p",self.pos_all[idx,idx_05_45_0],self.pos_all[idx,0],dd,dist_norm,'for',fm,'ro',forces_morse)




        if False:
            print("############################################################")
            print("# parametrize Michaels Model only along 110 (axial) (a,b,c,d)")
            print("############################################################")
            params_axial_abcd,params_axial_cd = my.get_michaels_paramerization(pos_all=self.pos_all,force_all=self.force_all,NN1=NN1,alat=alat,atoms=atoms,parametrize_only_idx=[idx_05_05_0,idx_45_45_0],rcut=0.88)
            mmodel_axial_abcd = hesse.Michael_poly_der
            print('params_axial_abcd',params_axial_abcd)
            print('params_axial_cd  ',params_axial_cd)
            np.savetxt("xy_CHECK_axial_abcd.dat",np.array([xred_dense,-1*mmodel_axial_abcd(xred_dense,*params_axial_abcd)]).T)
            np.savetxt("xy_CHECK_axial_cd.dat",np.array([xred_dense,-1*mmodel_axial_abcd(xred_dense,*params_axial_cd)]).T)


        if False:
            print("############################################################")
            print("# parametrize first Michael_poly_der numerical")
            print("# parametrize Morse over all neighbors")
            print("############################################################")
            params_m,params_m_,params_m__= my.get_michaels_paramerization(pos_all=self.pos_all,force_all=self.force_all,NN1=NN1,alat=alat,atoms=atoms,parametrize_only_idx=[idx_05_05_0,idx_45_45_0],function=hesse.Michael_poly_der)
            print('first only axial forces')
            mmodel_axial_abcd = hesse.Michael_poly_der
            print('params_axial_abcd  ',params_axial_abcd)
            print('params_axial_abcd__',params_axial_cd)
            np.savetxt("xy_CHECK_axial.dat",np.array([xred_dense,-1*mmodel_axial_abcd(xred_dense,*params_axial_abcd)]).T)
            np.savetxt("xy_CHECK_axial__.dat",np.array([xred_dense,-1*mmodel_axial_abcd(xred_dense,*params_axial_cd)]).T)
        sys.exit('1234556')

        print()
        print("#########################################################################")
        print("# fits ueber alle ersten nachbarn vom [0.0,0.0,0.0] atom")
        print("# first look at atoms in (001) plane")
        print("#########################################################################")
        import inspect, re
        def varname(p):
            for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
                m = re.search(r'\bvarname\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)', line)
                if m:
                    return m.group(1)

        indexes_NN              = [ idx_05_05_0, idx_05_45_0, idx_45_05_0, idx_45_45_0]
        indexes_NN_forces_sign  = [ -1.,         1          , 1          , 1 ]

        # repulsive and attractive
        indexes_NN              = [ idx_05_05_0, idx_45_45_0]
        indexes_NN_forces_sign  = [ -1.,         1. ]

        # repulsive and attractive and tox 0,0.5,0.5
        indexes_NN              = [ idx_05_05_0, idx_45_45_0, idx_0_05_05 ]
        indexes_NN_forces_sign  = [ -1.,         1.,          1 ]

        # only repulsive
        #indexes_NN              = [ idx_05_05_0]
        #indexes_NN_forces_sign  = [ -1.       ]

        # all atoms in the 001 plane
        indexes_NN              = [ idx_05_05_0, idx_45_45_0, idx_05_45_0, idx_45_05_0, idx_0_05_05, idx_0_45_45 ]
        indexes_NN_forces_sign  = [ -1.,         1          , 1          , 1          , 1          , 1]

        # all atoms in the 001 plane
        #indexes_NN              = [ idx_05_05_0,  idx_0_05_05]
        #indexes_NN_forces_sign  = [ -1.        , 1          ]
        x = []
        y = []
        ymi = []
        print('params_morse_normal    :',params_morse_normal)
        print('params_morse_normal_rel:',fit_rel.parameters)
        print('parameters             :',parameters)
        print('mmodel_function        :',mmodel_function)
        print('mmodel_function_der    :',mmodel_function_der)
        print('mmodel_key             :',mmodel_key)


        for idx_nn, indexname_nn in enumerate(indexes_NN):
            print()
            for disp_idx in np.arange(len(self.pos_all)):
                D  = self.pos_all[disp_idx,0]-self.pos_all[disp_idx,indexname_nn]
                dist_norm     = np.around(LA.norm(D),5)
                dist_norm_rel = dist_norm/self.alat
                x += [LA.norm(D)]
                Fv = self.force_all[disp_idx,indexname_nn]
                fv = indexes_NN_forces_sign[idx_nn]*LA.norm(Fv)   # affected by sign
                y += [fv]

                #### use a particular model
                if True: # do the test
                    fmo_test = hesse.Morse_derivative(LA.norm(D/self.alat), *fit_rel.parameters)
                    Fmo_test = hesse.getefvec(D/self.alat,fit_rel.parameters,pot = 'm')[1]

                ### use michaels model by summing over all neighbors
                #print('f_all')
                Fmo = -my.get_forces_on_atom_by_considering_all_its_1NNs(fit_rel.parameters,fit_rel.model,atoms,alat,atomnr=indexname_nn,pos=self.pos_all[disp_idx]/self.alat)
                Fma = -my.get_forces_on_atom_by_considering_all_its_1NNs(fit_rel_repulsive.parameters,hesse.Morse_repulsive_derivative,atoms,alat,atomnr=indexname_nn,pos=self.pos_all[disp_idx]/self.alat)
                f_all    = my.get_forces_on_atom_by_considering_all_its_1NNs(params_cd,mmodel_abcd,atoms,alat,atomnr=indexname_nn,pos=self.pos_all[disp_idx]/self.alat)
                f_all_ax = my.get_forces_on_atom_by_considering_all_its_1NNs(params_axial_cd,mmodel_abcd,atoms,alat,atomnr=indexname_nn,pos=self.pos_all[disp_idx]/self.alat)


                fmieq  = mmodel_function(1/np.sqrt(2.),*parameters)
                fmipd  = mmodel_function(dist_norm_rel,*parameters)
                fmi = fmipd - fmieq
                #fmaeq  = mmodel_function(1/np.sqrt(2.),*parameters)
                #fmapd  = mmodel_function(dist_norm_rel,*parameters)
                if False:
                    fal = mmodel_abcd(dist_norm_rel,*par_full) - mmodel_function(1/np.sqrt(2.),*par_full)

                fmi_test = hesse.Michael_poly_der_noA_noB_axial_subtract(dist_norm_rel,*parameters)
                ymi += [fmi]
                D0 = self.pos_all[0,0]-self.pos_all[0,indexname_nn]
                Fmieq = hesse.getefvec(D0/self.alat,parameters,pot = mmodel_key,paramsder=parameters_der)[1]
                Fmid  = hesse.getefvec(D/self.alat,parameters,pot = mmodel_key,paramsder=parameters_der)[1]
                #Fmieq = hesse.getefvec(D0/self.alat,parameters,pot = mmodel_key)[1] #,paramsder=parameters_der)[1]
                #Fmid  = hesse.getefvec(D/self.alat,parameters,pot = mmodel_key)[1] #,paramsder=parameters_der)[1]
                Fmi   = Fmid - Fmieq
                #if disp_idx == 0 or disp_idx == 7:
                #    print(
                #        '|fmipd|:',str(np.round(fmipd,5)).ljust(8),
                #        '-|fmieq|:',str(np.round(fmieq,5)).ljust(8),
                #        'Fv:',str(Fv).ljust(22),
                #        'Fmieq:',str(Fmieq).ljust(22),
                #        'Fmid:',str(Fmid).ljust(22),
                #        'Fmi:',str(Fmi).ljust(22)
                #            )
                print(
                        str(indexname_nn).ljust(3),
                        str(self.pos_all[disp_idx,indexname_nn]).ljust(19),
                        '|d|',str(np.round(dist_norm,4)).ljust(6),str(np.round(dist_norm_rel,4)).ljust(6),
                        'D:',str(D).ljust(23),
                        'D:',str(D/self.alat).ljust(23),
                        '||',
                        #'|fv!|',str(np.round(fv,5)).ljust(8),
                        #'|fmo_t|:',str(np.round(fmo_test,5)).ljust(8),
                        #'|fmi|:',str(np.round(fmi,5)).ljust(8),
                        #'|fmi_test|:',str(np.round(fmi_test,5)).ljust(8),
                        #'|fmipd|:',str(np.round(fmipd,5)).ljust(8),
                        #'-|fmieq|:',str(np.round(fmieq,5)).ljust(8),
                        #'Fv:',str(Fv).ljust(22),
                        #'Fmo:',str(Fmo).ljust(22),
                        #'Fma:',str(Fma).ljust(22),
                        #'Fmo_t:',str(Fmo_test).ljust(22),
                        #'Fmi:',str(Fmi).ljust(22),
                        #'Fmid:',str(Fmid).ljust(22),
                        #'Fmieq:',str(Fmieq).ljust(22),
                        #'Fmo:' ,str(Fmo).ljust(22),
                        'Fd_Fmo:' ,str(Fv - Fmo).ljust(22),
                        'Fd_Fma:' ,str(Fv - Fma).ljust(22),
                        'Fd_mi_ax:' ,str(Fv - Fmi).ljust(22),
                        'dF_all:'   ,str(Fv -f_all).ljust(22),
                        'dF_all_ax:',str(Fv -f_all_ax).ljust(22),
                        #'Fmieq:',str(Fmieq).ljust(22),
                        #'-> r=0.71',np.round(fmp,5),
                        #'r=|d|:',np.round(fmpd,5),
                        #'fmi:',np.round(fmfmii,5)
                        )
                if disp_idx == 7 and False:
                    # this need basically to loop over all interacting particles (here: first NN)
                    # Fv = sum over all forces
                    vec = D/self.alat               # given for every particle
                    vecnorm = r = np.linalg.norm(vec)   # given for every particle
                    vec_over_vecnorm = vec/vecnorm
                    fnorm = hesse.Michael_poly_der_noA_noB(r,*parameters)  # this gives only the force, I would like to get this
                    #fnorm =
                    force = vec_over_vecnorm * fnorm
                    print("@@!!",
                        'v_o_v:',str(vec_over_vecnorm).ljust(22),
                        'Fv:',str(Fv).ljust(22),
                        'D:',str(vec).ljust(23),
                        'Fmid:',str(Fmid).ljust(22),
                        'fnorm:',str(fnorm).ljust(22),
                        'force:',str(force).ljust(22),
                        " 1. get a function that for a given ATOM, gives all its 1NN (NN=a,b,c,d,...) (or its 2NN, 3NN, ...), basically, D for particular neighbors.",
                        " 2. write down the equation to solve: FORCE_VASP(ATOM)__x = SUM_over_all_neighbors[v_o_v__NN__x*Forcefunction(r)]",
                        "                                      FORCE_VASP(ATOM)__y = SUM_over_all_neighbors[v_o_v__NN__y*Forcefunction(r)]",
                        "                                      FORCE_VASP(ATOM)__z = SUM_over_all_neighbors[v_o_v__NN__z*Forcefunction(r)]",
                        " 3. This gives, for every displacement (=7), for every NN (=12), three (=x,y,z) equations.",
                        " was ich jetzt loesen muss ist ein gleichungssystem von vielen vektoren D and die VASP force, mache , basically, D for particular neighbors."
                        )
        print(      print())
        x = np.array(x)
        y = np.array(y)
        xrel = x/(self.alat/np.sqrt(2.))/np.sqrt(2)
        print('x',x)
        print('xrel',xrel)
        print('y',y)

        take_xrel = True
        add = ""
        if take_xrel == True:
            x = xrel
            add = "_xrel"

        np.savetxt("xy_from_VASP"+add+".dat",np.array(self.disp_vs_force_rel))
        np.savetxt("xy_chosen_displacements_"+add+".dat",np.array([x,y]).T)
        print('written',"xy_from_VASP"+add+".dat")
        print('written',"xy_chosen_displacements_"+add+".dat")



        if True:
            print("###################################################")
            print("# parametrize michaels function from original vasp forces (only axial)")
            print("###################################################")
            print('x',x)
            print('y',y)
            mmodel_function = hesse.Michael_poly_der_noA_noB_axial_subtract
            mmodel_function_resulting = hesse.Michael_poly_der_noA_noB
            mmodel = Model(mmodel_function)
            print('mmodel:',mmodel)
            mmodel_function_var = [ 'c','d' ]
            for i in mmodel_function_var: mmodel.set_param_hint(i ,value=1) #, min=0.01, max=0.5)
            result = mmodel.fit(y, r=x)
            parameters_mi = []
            for i in mmodel_function_var:
                coef = result.best_values.get(i)
                print(i,"Forces xx:",coef)
                parameters_mi.append(coef)

            np.savetxt("xy_chosen_displacements_"+add+"fit.dat",np.array([xred_dense,mmodel_function(xred_dense,*parameters_mi)]).T)
            np.savetxt("xy_chosen_displacements_"+add+"fit_resulting.dat",np.array([xred_dense,mmodel_function_resulting(xred_dense,*parameters_mi)]).T)
        #np.savetxt("xy_from_mich"+add+".dat",np.array([x,ymi]).T)
        #print(self.disp_vs_force)





        sys.exit()
        if True:
            print("###################################################")
            print("# parametrize Michael_polynomial_for_morsevalues_positive (x oder xrel)")
            print("###################################################")
            mmodel = Model(hesse.MP_for_morsevalues_positive)
            mmodel.set_param_hint('req' ,value=self.alat/np.sqrt(2.), vary=False)
            if take_xrel == True:
                mmodel.set_param_hint('req' ,value=1./np.sqrt(2.), vary=False)
            mmodel.set_param_hint('a' ,value=1) #, min=0.01, max=0.5)
            mmodel.set_param_hint('b' ,value=1) #,  min=1.0,  max=3.0)
            mmodel.set_param_hint('c' ,value=1) #,  min=1.0,  max=3.0)
            mmodel.set_param_hint('d' ,value=1) #,  min=1.0,  max=3.0)
            mmodel.set_param_hint('e' ,value=1) #,  min=1.0,  max=3.0)
            #mmodel.set_param_hint('req' ,value=1/np.sqrt(2.)) #,  min=1.0,  max=3.0)
            result = mmodel.fit(y, r=x)
            print('Michaels Model for polynomial positive',np.round([
                result.best_values.get('a'),
                result.best_values.get('b'),
                result.best_values.get('c'),
                result.best_values.get('d'),
                result.best_values.get('e'),
                result.best_values.get('req'),
                ],3))
            #print("forces with this model @eq:")
            #print('-->',result.best_values)
            np.savetxt("xy"+add+"_michael_for_polynomial_pos.dat",np.array([x,result.best_fit]).T)
            print('result.best_fit',result.best_fit)


        if True:
            print("###################################################")
            print("# parametrize MP (x oder xrel)")
            print("###################################################")
            mmodel = Model(hesse.MP)
            mmodel.set_param_hint('req' ,value=self.alat/np.sqrt(2.), vary=False)
            if take_xrel == True:
                mmodel.set_param_hint('req' ,value=1./np.sqrt(2.), vary=False)
            mmodel.set_param_hint('a' ,value=1) #, min=0.01, max=0.5)
            mmodel.set_param_hint('b' ,value=1) #,  min=1.0,  max=3.0)
            mmodel.set_param_hint('c' ,value=1) #,  min=1.0,  max=3.0)
            mmodel.set_param_hint('d' ,value=1) #,  min=1.0,  max=3.0)
            mmodel.set_param_hint('e' ,value=1) #,  min=1.0,  max=3.0)
            result = mmodel.fit(y, r=x)
            print('Michaels Model for rel x',np.round([
                result.best_values.get('a'),
                result.best_values.get('b'),
                result.best_values.get('c'),
                result.best_values.get('d'),
                result.best_values.get('e')
                ],3))
            np.savetxt("xy"+add+"_michael_polynomial.dat",np.array([x,result.best_fit]).T)


        sys.exit('78786')
        print()
        print('idx_05_05_0',idx_05_05_0,'diff pos:',self.pos_all[0,idx_05_05_0])
        print('idx_05_45_0',idx_05_45_0,'diff pos:',self.pos_all[0,idx_05_45_0])
        print('idx_45_05_0',idx_45_05_0,'diff pos:',self.pos_all[0,idx_45_05_0])
        print('idx_45_45_0',idx_45_45_0,'diff pos:',self.pos_all[0,idx_45_45_0])
        ##sys.exit()
        #print('force on 0_05_05 after subtracting mores')
        #print('fp',*fit.parameters)
        #fp = hesse.Morse_derivative(2.9, *fit.parameters)
        #print('fp 2.9 ',fp)
        #fp = hesse.Morse_derivative(2.95, *fit.parameters)
        #print('fp 2.95',fp)
        #print('at @ 05_05_0 -> pushes 0_05_05 in x and z')
        #print('at @ 00_00_0 -> pushes 0_05_05 in y and z and in sec approx in x y ')
        #sys.exit()


        ### here some stoff for disp 0.1
        print("dist_to_0_05_05",self.dist_to_0_05_05__01,"nndist:",self.nndist,"T=0K dist:",self.nndistT0K)
        #print("vec_to_0_05_05",self.vec_to_0_05_05__01)
        if type(fitmc1) != bool:
            force_repulsive_mc1 = pot_energy_forces.mc1_derivative(self.dist_to_0_05_05__01,*fitmc1.parameters)
        force_repulsive_morse = pot_energy_forces.Morse_derivative(self.dist_to_0_05_05__01,*fit.parameters)
       	print()
	print('----- now calcs --------')
        if type(fitmc1) != bool:
	    print("force_repulsive_mc1  ",force_repulsive_mc1)
	print("force_repulsive_morse",force_repulsive_morse)

	morseforce_on_0_05_05 = getforce(self.vec_to_0_05_05__01,force_repulsive_morse)
        if type(fitmc1) != bool:
            mc1force_on_0_05_05   = getforce(self.vec_to_0_05_05__01,force_repulsive_mc1)

        remaining_on_0_05_05_morse = self.f_0_05_05__01 - morseforce_on_0_05_05
        if type(fitmc1) != bool:
            remaining_on_0_05_05_mc1   = self.f_0_05_05__01 - mc1force_on_0_05_05

        print('morse forces on 0_05_05',morseforce_on_0_05_05)
        if type(fitmc1) != bool:
            print('mc1 forces on 0_05_05  ',mc1force_on_0_05_05)
	print()
        print('remaining on 0_05_05 (VASP - morse)',remaining_on_0_05_05_morse)
        if type(fitmc1) != bool:
            print('remaining on 0_05_05 (VASP - mc1)  ',remaining_on_0_05_05_mc1)
	print("")






        #print('normalize(morseforce_on_0_05_05)',normalize(morseforce_on_0_05_05))
        #print('normalize(remaining_on_0_05_05)',normalize(remaining_on_0_05_05))
        #parallelism = np.abs(np.dot(normalize(morseforce_on_0_05_05),normalize(remaining_on_0_05_05)))

        #print("parallelism",parallelism)
        print("LA.norm(remaining_on_0_05_05_morse",LA.norm(remaining_on_0_05_05_morse))
        if type(fitmc1) != bool:
            print("LA.norm(remaining_on_0_05_05_mc1  ",LA.norm(remaining_on_0_05_05_mc1))



        #self.ratio0510  = (LA.norm(f_1_1_0)*np.sign(f_1_1_0[0]))/LA.norm(f_05_05_0)
        #self.ratio0510 = round(self.ratio0510,3)

	self.ratio00505_morse  = round((LA.norm(remaining_on_0_05_05_morse)*np.sign(remaining_on_0_05_05_morse[0]))/LA.norm(f_05_05_0),3)
        if type(fitmc1) != bool:
	    self.ratio00505_mc1  = round((LA.norm(remaining_on_0_05_05_mc1)*np.sign(remaining_on_0_05_05_mc1[0]))/LA.norm(f_05_05_0),3)

	print('ratio remaining morse f[0,0.5,0.5]/f[1,1,0]',self.ratio00505_morse)
        if type(fitmc1) != bool:
	    print('ratio remaining mc1   f[0,0.5,0.5]/f[1,1,0]',self.ratio00505_mc1)
	print()
	print("@Forces on [0.5,0.5,0.0]")
	print("    dist  VASP   morse  mc1    dmorse  dmc1  Columb")
        print(fit_on_05_05_0)
	print()
        if type(fitmc1) != bool:
	    print(round(self.f_110,3),round(LA.norm(remaining_on_0_05_05_morse),3),round(LA.norm(remaining_on_0_05_05_mc1),3))
	print()
        if args.shift_parametrization_to_alat == False:
            if self.correct_for_110_forces == False:
                print("DONE alat {alat:4.3f} ||| alat_mor {alat_mor:4.3f} a_mor {a_mor:7.5f}  D_mor {D_mor:7.5f}".format(alat=self.alat,alat_mor=fit.parameters[2]*np.sqrt(2.),a_mor=fit.parameters[1],D_mor=fit.parameters[0]))
            elif self.correct_for_110_forces != False:
                print("DONE alat {alat:4.3f} ||| alat_mor {alat_mor:4.3f} a_mor {a_mor:7.5f}  D_mor {D_mor:7.5f}".format(alat=self.alat,alat_mor=fit_corrected_for_110.parameters[2]*np.sqrt(2.),a_mor=fit_corrected_for_110.parameters[1],D_mor=fit_corrected_for_110.parameters[0]))
        else:
            print("DONE alat {alat:4.3f} ||| alat_mor {alat_mor:4.3f} a_mor {a_mor:7.5f}  D_mor {D_mor:7.5f}".format(alat=self.alat,alat_mor=fit_shifted.parameters[2]*np.sqrt(2.),a_mor=fit_shifted.parameters[1],D_mor=fit_shifted.parameters[0]))

        if self.only_return_parametrization_file == True and os.path.isfile(parametrizationfile_morse):
            return parametrizationfile_morse
        else:
            return
	#print()
	#print("force on 0_05_05 given T=0K nndist: repulsive (r) if dist_to_0_05_05 < T=0K dist")
	#print("                       V")
	#print("f   1_1_0 055mo 055mc1   remain on 0_05_05morse remain on 0_05_05 mc1 ")
	#print("Pd ")
	#print("Pb ")
	#print("Al  [-0.016 -0.034 -0.03 ]  [-0.018 -0.029 -0.023]*")
	#print("Ir ")
	#print("Cu ")
	#print("Rh ")
	#print("Pt ")
	#print("Ag ")
	#print("Au ")
	#print()
        #print("- Todo: use the obtained morse and mc1 parametrization using for Al:")
        #print("        a) ")


        #####################################
        # rename -gnp (get_new_parametrization) to get_5x5x5_parametrization_at_tmelt -gp5x5x5Tmelt
        # CURRENT AIM -> TRY TO GO FROM PARAMETRIZATION from 4.04 to 4.14 !!  (TRY FIRST FOR AL, then for CU)
        # a) argparse_my_code: use -gnp to check weather the parametrization exists in the folder, if yes load it, if no, make it first and then load it
        # b) argparse_my_code: make sure same results as before --> delete stuff in /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Analyze_forces
        # c) look with Al /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400 parametrization what force exists when a distance of 4.13.
        # d) use force from c) to substract from the repulsive AND attractive forces of parametrization with 4.13
        # e) reparametrize the morse given this forces. (and check the std when using in the MD with 4.13)
        # f) if that helps, to this also for Cu
        # -) now one could also take into account that atom at 05_05_0 pulls on atom at 0_05_05, anyhow (only relevatn for tox)
        #
        # DISP: argparse_md_code.py -e Al -gnp -waf -rp /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400/POSITIONs --alat 4.04 --alat_mor 4.13 -v
        # DISP: argparse_md_code.py -e Al -gnp -waf -rp /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400/POSITIONs --alat 4.04 --alat_mor 4.05 -v
        # DISP:                                                             ==> Al std:    0.0 dudl/2:    0.0 alat_mor: 4.13 alat: 4.04 delta: 1.02   (std force 0.014)
        # DISP:                                                             ==> Al std:    0.0 dudl/2:    0.0 alat_mor: 4.04 alat: 4.04 delta: 1.00   (std force 0.035)
        # MD: argparse_md_code.py -md2 -e Al -gnp -dbl -waf                 ==> Al std:    6.5 dudl/2:    5.3 alat_mor: 4.13 alat: 4.13 delta: 1.00   (std force 0.110)
        # MD: argparse_md_code.py -md2 -e Al -gnp_alat_mor_T0K -dbl -waf    ==> Al std:    4.0 dudl/2:    1.9 alat_mor: 4.05 alat: 4.13 delta: 0.98   (std force 0.104)
        #  .... for alat of 4.13 (MD)  -> a_mor 4.05 is better
        #  .... for alat of 4.04 (DISP)-> a_mor 4.13 is better  -> check for the DISP weather another parametrization with 4.05 could be better (should)
        #                              -> ok, this is clear, since the parametrization used (-gnp) was done for a 4.13 alat and does not fit for 4.04 at all!
        #                              -> NEXT: parametrize DISP from 4.04 and check std for 4.04 (DISP) and 4.13 (DISP & MD)
        #                              -> UNDERSTAND HOW TO TRANSLATE PARAMETRIZATION FROM 4.04 to 4.13 and vice versa! (HOPEFULLY THE MORSE PARAMETERS WILL BE THE SAME OR CLOSE AFTER TRANSFORMATION
        #                              -> 1. morse parametrisation @4.04
        #                              -> 2. morse parametrisation @4.13 BUT the force of 4.04parametrization betwenn atom at 1_1_0 and 05_05_0 has to be substracted from 05_05_0
        # -> DISP is better with alat_mor = 4.13
        # -> MD is better for alat_mor = 4.05
        # -----> however checked for different alats DISP vs MD
        #print("          - for argparse_md_code.py -e Al -gnp -waf -rp /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400/POSITIONs --alat 4.04 --alat_mor 4.{04,13} -v   --> make sure that no std Energy std (DFT-LA)  :  4.026 meV/atom (works in general lon and tox) IS NOT WRITTEN! (since it is not available); then -> make it available and compare it!)")



        #print("          - fir morse without alat_mor.")
        #print("          - do argparse_md_code.py -md2 -ea -gnp_alat_mor_T0K -dbl -vm2 but get rid of the angle (e.e. for Cu) by adjusting the other morse parameter from T=0K")

        #print("          - try for all elements the T=0K parametrization for the Tmelt lattice constant.")
        #print("            for al: first or second try: self.displacements = /Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400/")
        #print("            for al: second try: use the forces obtained with the 4.14 (Tmelt) displacements using the 5x5x5sc but fixing the a_mor to 4.045 angstom and ")
        #print("                    writing force_on_05_05_0 == Morse(attractive_to__1_1_0) + Morse(XXX_to_disp_from_0_0_0) where XXX is attr and rep depending on dist between 0_0_0displaced and 05_05_0; the Morse(attractive_to__1_1_0) will alwasy want to move the atom towards atom at 1_1_0 ")
        #print("            for al: best would be having the murn T=0K instead of the phonon one.")
        #print("          - for Ir the 5sc parametrisation uses T=0K lattice constant. check parametrisation with alat @ Tmelt.")
        #print("          - check weather (for all elements) with smaller lattice constant (maybe from T=0K) better fit with smaller std in MD can be made, start this for Al!")
        #print("          - check to what end the difference in fitted a_mor vs MD alat make a difference in std (change in parametrization)")
        #print("          - look at these NEW (-gnp) parametrizations as functions of displacemetn vasp vs fitted")
        #print("          - note the corresponding \"f   1_1_0 055mo 055mc1   remain on 0_05_05morse remain on 0_05_05 mc1\" lines")
        #print("          - forces_vs_dft_ANGLE_FROM_TOX.agr: show that negative forces are from tox (probably), those show best how angle will be")
        #print("          - think how this can be cured")
        #print("          - check weather with smaller lattice constant (maybe from T=0K) better fit with smaller std in MD can be made")
        #print("          - make forces_vs_dft_ANGLE_FROM_TOX.agr for other elements (from 5x5x5 supercell) save in Analyze_forces and make small readme or better: make there the corresponding folder for the elements.")
        #print("          - make forces_vs_dft_ANGLE_FROM_TOX.agr in a 2x2x2 supercell and the 3x3x3 supercell for al")
        #print("          - ")
        #print("        b) the used displacements to get std mit mc1")
        #print("        c) a 2x2x2 md to get std mit morse")
        #print("        d) a 2x2x2 md to get std mit mc1")
        #print("        e) a 4x4x4 md to get std mit morse")
        #print("        f) a 4x4x4 md to get std mit mc1")
        #print("        g) do everything up to 2x2x2 md for every other element")
        #print("        h) now try to play with the parametrization of the displacements.")
        #print()
	#print("- Die besten fits sind die wo die remaining y komponente klein ist!")
	#print("- Die schlechtesten y komponenten (morse) sind fuer Al, Ir, Rh, Au")
	#print("- Waere da auch noch Ir dabei dann waere es sicher dass die y komponente ...")
	#print("  ... mit dem fit klein gemacht werden muss um einen guten fit zu bekommen")
	#print("- Bei Ir haben wir jedoch die groessten kraefte auf 110.")
	#print("-> a) BENUTZE(mc1)/MACHE parametrisierungen wo die verbleibende ...")
	#print("      ... y komponente klein ist!")
	#print("-> aa) dies bedeutet dass der run (MD) mit mc1 schon viel besser sein muesste fuer mc1 als fuer morse fuer Al, Ir, Pt, Au fuer Cu hingegen koennte Morse besser sein.")
	#print("-> aa) die restkraefte auf das atom bei [0,0.5,0.5] (d2) sind 4 mal so ...")
	#print("       wichgit wie die restkraefte auf das atom bei [0.5,0.5,0] (d1)")
	#print("-> ab) zu loesen ist: M = Morse : M(d1) & 4*M(d2) wobei die kraefte auf d2 so zu waehlen sind dass die y komponente auf [0,0.5,0.5] verschwindet.")

	#print("-> aa) dies hiesse dass die repulsive flanke gaendert wird (schlechter gemacht wird) -> schaue ob dies verbessert werden kann wenn man")
	#print("-> b) checke diese parametrisierung (MD) vs originale morse parameter.")
	#print("- Wenn die verbleigende kraft auf 0_05_05 in y 0 ist, sollte man eine Wechselwirking vom atom bei 0.5,0.5,0 auf das atom 0,0.5,0.5 machen. (3body), d.h. wenn man bei dem atom 0,0,0 ist ist ein paerchen die atome [0.5,0.5,0] und [ 0,0.5,0.5] bei welchen dann (da beide repulsiv zum atom 0,0,0 istn) eine relativ starke restliche kraft bleibt.")




#for element in elements:
#    ele = get_one_disp(element,0.1)
#    ratio = round(ele.f_1_1_0[0]/ele.f_05_05_0[0]*514.)
#    rationorm = round(np.linalg.norm(ele.f_1_1_0)/np.linalg.norm(ele.f_05_05_0),2)*514.
#    #remaining_on_0_05_05 = np.loadtxt(path+'remaining_on_0_05_05.'+element+str(ele.alat)+'.dat')
#    rationorm2 = round(np.linalg.norm(remaining_on_0_05_05)/np.linalg.norm(ele.f_1_1_0),2)*514.
#    print(ele.element,ele.alat,ratio,rationorm,rationorm2,"scale tox? f110[0]:",ele.forces[ele.idx_1_1_0][0],'remain 0_05_05')#,#remaining_on_0_05_05,'f05050',ele.forces[ele.idx_05_05_0][0])
#print()

if __name__ == '__main__':
    def help(p = None):
        string = '''
        export base="/Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/"
        e.g. get_parametrization_for_displacementfolder.py -f . -v   # being in a particular forlder ....
        e.g. get_parametrization_for_displacementfolder.py -e Al -sc 5 -v     # to make the parametrization only for Al (in the 5x5x5 supercell)
        e.g. get_parametrization_for_displacementfolder.py -e Al -sc 5 -rp    # just to get the path to the morse parametrization
        e.g. get_parametrization_for_displacementfolder.py -ea   -sc 5        # to make the parametrization in every 5x5x5 folder
        e.g. get_parametrization_for_displacementfolder.py -f $base/Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400
        e.g. get_parametrization_for_displacementfolder.py -f       Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400   # works as well
        e.g. get_parametrization_for_displacementfolder.py -f . -spta 4.09 -rp
        e.g. get_parametrization_for_displacementfolder.py -e Al --folder_alat_lattice_T0K    # searches a particular folder for the parametrization
        e.g. get_parametrization_for_displacementfolder.py -e Al --folder_alat_lattice_T0K -v -spta 4.13   # to get corresponding shifted parametrization at other alat
        e.g. get_parametrization_for_displacementfolder.py -e Al -f /Users/glensk/Dropbox/Albert/Understanding_distributions/displacements_/Al/3x3x3sc_4.14Ang_quer_10x10x10kp_vasp4_ENCUT400/ -v





        '''
        p = argparse.ArgumentParser(description=string,
                formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
        p.add_argument('-e',   '--element', required=False,
                help='Element used to set the mass', type=str, default=False)
        p.add_argument('-ea',   '--element_all' , default=False,required=False,action='store_true',
                help='Go over all available fcc elements')
        p.add_argument('-sc',   '--supercell', '-N' , required=False,
                help='supercell size (default = 2) -> 2*2*2*4 = 32 atoms for fcc', type=int, default=False)
        p.add_argument('-f',   '--folder', required=False,help='Define folder where parametrisation should be made', type=str, default=False)
        p.add_argument('-f_alat_lat_T0K',   '--folder_alat_lattice_T0K', required=False, action='store_true',
                help='get a parametrisation which has been done @T=0K alat', default=False)

        p.add_argument('-rp','--only_return_parametrization_file',required=False,action='store_true',default=False,help='returns the path to the parametrization file')


        p.add_argument('-spta', '--shift_parametrization_to_alat', required=False, type=float, default=False, help="at another atom distance, the morse would have a static force which cancels out in equilibrium. Find the static force for the chosed alat_lattice and shift the forces to parametrize correspondingly.")
        p.add_argument('-cf110', '--correct_for_110_forces', required=False, type=float, default=False, help="correct forces on 05_05_0 atom by the forces on 1_1_0 atom.")


        p.add_argument('-v','--verbose',help='verbose', action='count', default=False)
        p.add_argument('-v2','--verbose2',help='verbose2', action='count', default=False)
        return p
    p = help()
    args = p.parse_args()

    ### get the elements
    if args.element:
        elements = [args.element]
    else:
        elements = [False]

    if args.element_all:
        elements = ["Pd","Pb","Al","Ir","Cu","Rh","Pt","Ag","Au"]



    for element in elements:
        ele = get_all_disps(
                element=element,
                sc=args.supercell,
                dofor=args.folder,
                verbose=args.verbose,
                verbose2=args.verbose2,
                only_return_parametrization_file=args.only_return_parametrization_file,

                folder_alat_lattice_T0K = args.folder_alat_lattice_T0K,
                shift_parametrization_to_alat = args.shift_parametrization_to_alat,
                correct_for_110_forces = args.correct_for_110_forces,
                )
        ele.parametrize_it()

