#!/usr/bin/env python

from __future__ import print_function
import numpy as np
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

            self.p_1_1_0 = self.pos[self.idx_1_1_0]
            self.p_05_05_0 = self.pos[self.idx_05_05_0]

            self.f_1_1_0 = self.forces[self.idx_1_1_0]
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
        str_list = utils.common_prefix(self.folder_alldisp).split("/")
        #print(str_list)
        str_list = filter(None, str_list)
        #print(str_list)
        str_list= "/"+"/".join(str_list[:-1])
        #print(str_list)
        dofor = str_list


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
        disp_vs_force1          = np.zeros((len(self.folder_alldisp),2))
        disp_vs_force2          = np.zeros((len(self.folder_alldisp),2))
        disp_vs_force1_110      = np.zeros((len(self.folder_alldisp),2))
        disp_vs_force2_110      = np.zeros((len(self.folder_alldisp),2))
        self.disp_vs_force      = np.zeros((len(self.folder_alldisp)*2,2))
        self.disp_vs_force_110  = np.zeros((len(self.folder_alldisp)*2,2))
        self.dist_0_05_05_at0   = np.zeros((len(self.folder_alldisp),3))
        self.dist_0_05_05_at1   = np.zeros((len(self.folder_alldisp),3))

        self.force_0_05_05      = np.zeros((len(self.folder_alldisp),3))
        self.force_45_45_0      = np.zeros((len(self.folder_alldisp),3))
        self.force_0_05_05_rest = np.zeros((len(self.folder_alldisp),3))
        self.force_05_05_0_rest = np.zeros((len(self.folder_alldisp),3))

        da = len(self.folder_alldisp)
        if self.verbose:
            print()
            for idx,i in enumerate(self.folder_alldisp):
                print(idx,i)
            print()

        create_dofor_POSITIONs = False
        create_dofor_u_OUTCAR = False

        if not os.path.isfile(self.dofor+"/POSITIONs"):
            create_dofor_POSITIONs = True
            open(self.dofor+"/POSITIONs", 'a').close()
        if not os.path.isfile(self.dofor+"/u_OUTCAR"):
            create_dofor_u_OUTCAR = True
            open(self.dofor+"/u_OUTCAR", 'a').close()

        for idx,i in enumerate(self.folder_alldisp):
            if idx == 0 and self.verbose:
                print("folder_alldisp[i] (W):",i)
            disp = float(i.split("Ang_")[-1].split("/")[0])
            #try:
            #    print(i.split("Ang_"))
            #    sys.exit()
            #except ValueError:
            #    print(i.split("Ang_"))
            #    sys.exit()
            self.disps.append(disp)
            if not os.path.isfile(i+'/pos'):
                with my.cd(i):
                    call(["OUTCAR_positions-last-ARRAY.sh > pos"],shell=True)
            if not os.path.isfile(i+'/forces'):
                with my.cd(i):
                    call(["OUTCAR_forces-last-ARRAY.sh > forces"],shell=True)
            if not os.path.isfile(i+'/POSITIONs'):
                with my.cd(i):
                    call(["extractPOSITIONS.sh"],shell=True)
            if not os.path.isfile(i+'/u_OUTCAR'):
                with my.cd(i):
                    call(["OUTCAR_ene-potential_energy_without_first_substracted.sh > u_OUTCAR"],shell=True)

            if create_dofor_POSITIONs == True:
                with my.cd(i):
                    call(["cat POSITIONs >> ../POSITIONs"],shell=True)
            if create_dofor_u_OUTCAR == True:
                with my.cd(i):
                    call(["cat u_OUTCAR >> ../u_OUTCAR"],shell=True)
            #alat = float(i.split("vasp4/")[1].split("Ang")[0])

            if self.verbose > 1:
                print('loadtxt (pos):',i+'/pos')
                print('loadtxt (forces):',i+'/forces')
            pos = np.loadtxt(i+'/pos')

            ### get self.sc in case it is necessary
            if self.sc == False:
                #print('ps',pos.shape[0])
                #print('ps',pos.shape[0]/4)
                #print('ps',(pos.shape[0]/4)**(1./3.))
                #print('ps',round((pos.shape[0]/4)**(1./3.),10))
                #print('ps',int(round((pos.shape[0]/4)**(1./3.),10)))
                self.sc = int(round((pos.shape[0]/4)**(1./3.),10))
            if self.verbose > 1:
                print('self.sc',self.sc)
                #print(pos[:10])
                #print(pos.flatten()[:10])
                #print(np.sort(pos.flatten())[:10])
                #print(np.unique(np.sort(pos.flatten()))[:10])

            if self.alat == False:
                if pos[0][2] == pos[1][0] == pos[1][1] == 0:
                    alat = pos[1][2]
                #alat = np.unique(np.sort(pos.flatten()))[2]
                if alat > 6 or alat < 3:
                    print('alat',alat)
                    sys.exit("this does not seem to be a correct alat")
                self.alat = alat
                self.nndist = self.alat/np.sqrt(2.)
            if self.verbose > 1:
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


            forces = np.loadtxt(i+'/forces')
            #if self.verbose:
            #    print('pos')
            #    print(pos[-1])
            #    print("forces")
            #    print(forces[-1])


            ##########################################################################################
            # This is so for for a random displacement, whichever was loaded last
            ##########################################################################################
            idx_05_05_0 = getindex(pos/self.alat,np.array([0.5,0.5,0]))
            #print('idx_05_05_0 (the index ot the atom at 05_05_0)',idx_05_05_0)
            p_05_05_0 = pos[idx_05_05_0]
            f_05_05_0 = forces[idx_05_05_0]

            idx_1_1_0 = getindex(pos/self.alat,np.array([1.,1.,0]))
            p_1_1_0 = pos[idx_1_1_0]
            f_1_1_0 = forces[idx_1_1_0]

            if self.verbose > 1:
                print("idx_05_05_0",idx_05_05_0)
                print("p_05_05_0",p_05_05_0)
                print("f_05_05_0",f_05_05_0)

            #def dosearch():


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

            idx_45_45_0 = dosearch([round(self.sc-0.5,12),round(self.sc-0.5,12),0],"idx_45_45_0",pos,self.alat)
            p_45_45_0   =    pos[idx_45_45_0]
            f_45_45_0   = forces[idx_45_45_0]

            idx_40_40_0 = dosearch([round(self.sc-1.0,12),round(self.sc-1.0,12),0],"idx_40_40_0",pos,self.alat)
            p_40_40_0   =    pos[idx_40_40_0]
            f_40_40_0   = forces[idx_40_40_0]

            idx_45_0_45 = dosearch([round(self.sc-0.5,12),0,round(self.sc-0.5,12)],"idx_45_0_45",pos,self.alat)
            p_45_0_45   =    pos[idx_45_0_45]
            f_45_0_45   = forces[idx_45_0_45]

            if self.verbose > 1:
                print("idx_45_45_0",idx_45_45_0)
                print('idx',disp,'repulsive force:',f_05_05_0,'attractive force:',f_45_45_0)

            idx_0_05_05 = getindex(pos/self.alat,np.array([0,0.5,0.5]))
            p_0_05_05 = pos[idx_0_05_05]
            f_0_05_05 = forces[idx_0_05_05]

            d_0_05_05_at0 = p_0_05_05-pos[0]
            d_0_05_05_at1 = p_0_05_05 - p_05_05_0

            self.dist_0_05_05_at0[idx] = d_0_05_05_at0
            self.dist_0_05_05_at1[idx] = d_0_05_05_at1
            self.force_0_05_05[idx] = f_0_05_05
            self.force_45_45_0[idx] = f_45_45_0
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
		add = "ATTRACTIVE! to "
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
            # put together disp vs forces
            ##############################################################################
            disp_vs_force1[idx,0] = LA.norm(p_05_05_0-pos[0])  # distance
            disp_vs_force1[idx,1] = - np.sign(f_05_05_0[0])*LA.norm(f_05_05_0)       # forces repulsive
            #print('disp_vs_force1 == repulsive')
            #print(disp_vs_force1)

            disp_vs_force2[idx,0] = self.nndist + disp         #
            disp_vs_force2[idx,1] = np.sign(f_45_45_0[0])*LA.norm(f_45_45_0)
            #print('disp_vs_force2 == attractive')
            #print(disp_vs_force2)

            self.disp_vs_force[idx,0] = disp_vs_force1[idx,0]
            self.disp_vs_force[idx,1] = disp_vs_force1[idx,1]
            self.disp_vs_force[idx+da,0] = disp_vs_force2[idx,0]
            self.disp_vs_force[idx+da,1] = disp_vs_force2[idx,1]
            #print('disp_vs_force == all')
            #print(disp_vs_force)
            #sys.exit

            ##############################################################################
            # put together disp vs forces at atom 1_1_0 to substract from force @ 05_05_0 as a correction
            ##############################################################################
            disp_vs_force1_110[idx,0] = LA.norm(p_05_05_0-pos[0])  # distance
            disp_vs_force1_110[idx,1] = - np.sign(f_1_1_0[0])*LA.norm(f_1_1_0)           # forces repulsive ... in this case the sign is different
            disp_vs_force2_110[idx,0] = self.nndist + disp         #
            disp_vs_force2_110[idx,1] = np.sign(f_40_40_0[0])*LA.norm(f_40_40_0)
            self.disp_vs_force_110[idx,0] = disp_vs_force1_110[idx,0]
            self.disp_vs_force_110[idx,1] = disp_vs_force1_110[idx,1]
            self.disp_vs_force_110[idx+da,0] = disp_vs_force2_110[idx,0]
            self.disp_vs_force_110[idx+da,1] = disp_vs_force2_110[idx,1]


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
        self.disp_vs_force = my.remove_duplicates_in_numpy_xy_array_and_sort(self.disp_vs_force,roundto=10)
        self.disp_vs_force_110 = my.remove_duplicates_in_numpy_xy_array_and_sort(self.disp_vs_force_110,roundto=10)
        self.disp_vs_force_corrected_for_110 = np.copy(self.disp_vs_force)
        self.disp_vs_force_corrected_for_110[:,1] += self.disp_vs_force_110[:,1]
        if self.verbose > 2:
            print('---------------------------')
            print('self.disp_vs_force')
            print('---------------------------')
            print(self.disp_vs_force)
            print('---------------------------')
        #print('??')
        #print(self.disp_vs_force_110)
        #print('corr')
        #print(self.disp_vs_force_corrected_for_110)
        #sys.exit()


        print("#########################################################################")
        print("# fits on [0.5,0.5,0.0]")
        print("#########################################################################")
        fit = pot_parametrize.fit_to_func(self.disp_vs_force,function='morse',fixzeroat=self.nndist)
        try:
            fitmc1 = pot_parametrize.fit_to_func(self.disp_vs_force,function='mc1',fixzeroat=self.nndist)
        except TypeError:
            fitmc1 = False
        fit_corrected_for_110 = pot_parametrize.fit_to_func(self.disp_vs_force_corrected_for_110,function='morse',fixzeroat=self.nndist)
        print('fit.parameters',fit.parameters)

        fit_on_05_05_0 = np.zeros((len(fit.fit),7))
        fit_on_05_05_0[:,0] = fit.fit[:,0]
        dist__hydrogen   = (1./7.2)*fit.fit[:,0]*10**-10
        dist__           = (1./1.0)*fit.fit[:,0]*10**-10
        fit_on_05_05_0[:,1] = self.disp_vs_force[:,1]   # actual (VASP) forces
        fit_on_05_05_0[:,2] = fit.fit[:,1]              # fitted morse forces
        print('dist__',dist__)
        for idx,i in enumerate(dist__):
            print('d',idx,dist__[idx],'force VASP',fit_on_05_05_0[idx,1])
        k__ = (0.9*10**10)
        e__ = (1.6*10**(-19))
        dist__meter = dist__*10**(-10)   # angstrom to meter
        newton_to_mev_per_angstrom = 6.2415091**11
        q1__ = ((1.6*78.)*10**(-19)) #78
        q1__ = ((1.6*2.8)*10**(-19)) #78
        fit_on_05_05_0[:,6] =  (-k__*((q1__**2) /(dist__**2))*newton_to_mev_per_angstrom)+12.              # columb forces

        # get screened potential (https://en.wikipedia.org/wiki/Electric-field_screening)
        for idx,i in enumerate(dist__):
            print('d',idx,dist__[idx],fit_on_05_05_0[idx,6])
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
        if type(fitmc1) != bool:
            np.savetxt(self.dofor+'/disp_vs_fitted_mc1.dat',fitmc1.fit)
        np.savetxt(self.dofor+'/disp_vs_fitted_morse.dat',fit.fit)


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

        ##############################################################################
        # force on 0_05_05 (previously tox)
        # Weight 1Morse 05_05_0
        # Weight 4Morse 0_05_05
        ##############################################################################
	print("xx@Forces on [0.5,0.5,0.0]")
	print("xx    dist  VASP   morse  mc1    dmorse  dmc1  Columb")
        print("xx",fit_on_05_05_0)
        print("#########################################################################")
        print("# fits on [0.0,0.5,0.5] previously tox")
        print("#########################################################################")
        if True:
            print('dist between pos[0] and pos_0_05_05')
            params = fit.parameters
            print('params',params)
            for idx,i in enumerate(self.dist_0_05_05_at0):
                dist_norm                           = np.around(LA.norm(i),5)


                fm                                  = np.round(hesse.Morse_derivative(LA.norm(i), *params),4)
                forces_morse                        = hesse.getefvec(i,params,pot = 'm')
                self.force_0_05_05_rest[idx]        = self.force_0_05_05[idx] + forces_morse[1]
                vasp_force                          = self.force_0_05_05[idx] #.ljust(22)
                vasp_force_norm                     = str(np.around(-LA.norm(abs(vasp_force)),decimals=3)).ljust(6)
                vasp_minus_morse                    = vasp_force+forces_morse[1]
                print(str(i).ljust(22),str(dist_norm).ljust(7),'F_0_05_05_vasp',str(vasp_force).ljust(22),"norm:",vasp_force_norm,'F_morse',str(-forces_morse[1]).ljust(22),'norm',str(fm).ljust(7),'vasp-morse',vasp_minus_morse,'norm',str(np.around(LA.norm(vasp_minus_morse),3)).ljust(7))
            #print(self.dist_0_05_05_at0)
            sys.exit('23')
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
            print(disp_vs_force1)
            morsemodel = Model(hesse.Morse_derivative)
            morsemodel.set_param_hint('re' ,value=self.alat/np.sqrt(2.), vary=False)
            morsemodel.set_param_hint('De' ,value=0.25, min=0.01, max=0.5)
            morsemodel.set_param_hint('aa' ,value=1.5,  min=1.0,  max=3.0)


            r1 = disp_vs_force1[:,0]
            y1 = disp_vs_force1[:,1]
            #for idx,i in enumerate(r1):print(r1[idx],y1[idx])
            result = morsemodel.fit(y1, r=r1)
            print('vgl (obtained by fit)',params)
            print('vgl (only repulsive )',np.round([result.best_values.get('De'),result.best_values.get('aa'),result.best_values.get('re')],3))

            r1a = disp_vs_force2[:,0]
            y1a = disp_vs_force2[:,1]
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
                fv = hesse.getefvec(i,params,pot = 'm')
                self.force_0_05_05_rest[idx] = self.force_0_05_05[idx] + fv[1]
                print(i,LA.norm(i),'force_0_05_05 vasp_full',self.force_0_05_05[idx],"norm",-LA.norm(abs(self.force_0_05_05[idx])),'force morse',fm,'fv',-fv[1]) #,LA.norm(fv[1]))
            print('this parametrization, even for T=0K displacements, is not optimal!, try to add attractive forces on 0_45_45')


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
        e.g. get_parametrization_for_displacementfolder.py -e Al -sc 5 -v     # to make the parametrization only for Al (in the 5x5x5 supercell)
        e.g. get_parametrization_for_displacementfolder.py -e Al -sc 5 -rp    # just to get the path to the morse parametrization
        e.g. get_parametrization_for_displacementfolder.py -ea   -sc 5        # to make the parametrization in every 5x5x5 folder
        e.g. get_parametrization_for_displacementfolder.py -f $base/Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400
        e.g. get_parametrization_for_displacementfolder.py -f       Al/3x3x3sc_4.04Ang_quer_10x10x10kp_vasp4_ENCUT400   # works as well
        e.g. get_parametrization_for_displacementfolder.py -f . -spta 4.09 -rp
        e.g. get_parametrization_for_displacementfolder.py -e Al --folder_alat_lattice_T0K    # searches a particular folder for the parametrization
        e.g. get_parametrization_for_displacementfolder.py -e Al --folder_alat_lattice_T0K -v -spta 4.13   # to get corresponding shifted parametrization at other alat




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

