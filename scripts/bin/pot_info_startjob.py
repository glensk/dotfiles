#!/usr/bin/env python

####################################
# examples:
####################################
# run ~/Thermodynamics/python_thermodynamics/pot_info_startjob.py -v -e al -dispdirection quer -a 4.07 -kp 10x10x10kp -sp 3 -stradd "_vasp4_ENCUT400" -pnnforcesnpz

# run /Users/glensk/Thermodynamics/python_thermodynamics/pot_info_startjob.py  -e al -dispdirection quer -a 4.13 -v -pnnforcesnpz
#
# run /Users/glensk/Thermodynamics/python_thermodynamics/pot_info_startjob.py  -e al -dispdirection quer -a 4.16 -ppkldata
#
#
# run  /Users/glensk/Thermodynamics/python_thermodynamics/pot_info_startjob.py  -e al -a 4.14 -sp 5 -kp 2x2x2kp -dispdirection quer -pnnforcesnpz -v -stradd "_vasp4_ENCUT250_GGA_2_atoms_displaced"
#
#
# pot.disp.alle.atoms = 500
# pot.disp.alle.struct = 'fcc'
# pot.forparametrization_sc = 5
# pot.disp.alle.a = 4.14
# pot.alat --> 4.13
#
# pot.disp.alle.p.shape -> (1800, 500, 3)
# pot.disp.alle.pr.shape -> (1800, 500, 3)
# AM ANFANG SIND ERSTMAL alle.p == alle.pr
#
#
#
#
#
# In [39]: pot.disp.alle.sc
# Out[39]:
# array([[ 8.26,  0.  ,  0.  ],
#        [ 0.  ,  8.26,  0.  ],
#        [ 0.  ,  0.  ,  8.26]])
#
# In [41]: pot.disp.quer.{p,f}.shape   # from cartesian_coords/ forces_OUTCAR vasp  normal cell
# Out[41]: (200, 32, 3)
#
# pot.disp.quer.e{1-200} # gives the ene_free_last
#
# In [58]: pot.disp.quer.{P,F}.shape # positions in big repeated cell
# Out[58]: (200, 256, 3)
#
#
# In [41]: pot.disp.quer.pr.shape  # removes mapping from positions
# Out[41]: (200, 32, 3)

# In [41]: pot.disp.quer.Pr.shape  # removes mapping from positions in big cell
# Out[41]: (200, 256, 3)
#
# In [74]: pot.disp.quer.P0.shape  # undisplaced structure big cell
# Out[74]: (256, 3)
#
#
# In [69]: pot.disp.quer.dvec[1]   # displacement of first atom
# Out[69]: array([ 0.007,  0.007,  0.   ])
#
# In [70]: pot.disp.quer.dnorm[1]  # displacementnorm of first atom
# Out[70]: 0.0099984898859777818
#
#
#
# In [45]: pot.disp.quer.p[0]  -> undisplaced cartesian positions
# In [45]: pot.disp.quer.p[1]  -> first atome displaced ...
#
# In [45]: pot.disp.quer.f[0]  -> vasp forces
# In [45]: pot.disp.quer.f[1]  -> vasp forces for for first atome displaced ...


# pot.disp.alle.data[1,'lon','quer','alle'].file  # path to file with 1NN long forces
#
# instead of 'alle' use 'dos'A
#
# pot.disp.quer.f[:,0]  gives all forces on first atom
# pot.disp.quer.p[:,0]  gives all position of fist atom




########################################################
# Morse parametrization
########################################################
# --> for 4.13 supercell Al (2x2x2 supercell)
# In [107]: pot.disp.alle.data[1,'lon','quer','dos'].morse
# Out[107]: array([ 0.282,  1.427,  2.92 ])
#
# In [108]: pot.disp.alle.data[1,'lon','quer','dos'].morse[0]
# Out[108]: 0.28215356970188138
#
# In [109]: pot.disp.alle.data[1,'lon','quer','dos'].morse[1]
# Out[109]: 1.427037696092982
#
# In [110]: pot.disp.alle.data[1,'lon','quer','dos'].morse[2]
# Out[110]: 2.9203510063004408

# --> for 4.16 supercell Al (2x2x2 supercell)
# In [119]: pot.disp.alle.data[1,'lon','quer','dos'].morse[0]
# Out[119]: 0.24630777411338961
#
# In [120]: pot.disp.alle.data[1,'lon','quer','dos'].morse[1]
# Out[120]: 1.4565490150540261
#
# In [121]: pot.disp.alle.data[1,'lon','quer','dos'].morse[2]
# Out[121]: 2.9415642097360379
#

import os
import sys
import copy
import numpy as np
import glob
import shutil
import argparse  # PYTHON_ARGCOMPLETE_OK
from argparse import RawTextHelpFormatter
import utils
reload(utils)
import pot_parametrize
reload(pot_parametrize)
import pot_energy_forces
reload(pot_energy_forces)
import hesse
reload(hesse)
import DOS_POSITIONS_auswerten
reload(DOS_POSITIONS_auswerten)
import qubus_interpolate
reload(qubus_interpolate)

import pot_energy_forces
reload(pot_energy_forces)



#def effectlontox(foldernrs=args.effectlontox, verbose = False):
def effectlontox(foldernrs, verbose = False):
    ''' comment '''
    if verbose:
        print utils.printblue("get effect of TOX ...")
    ##################################################
    # get effect of {TOX,TOY,TOZ,(TIX,TIY,TIZ)}
    ##################################################
    # 0. OK get both jobnumbers e.g. 3 && 0
    # 0. define the folder of interest f1 f2
    # 0. try to find out what is diffrent of both filenames
    # 1. get dUdLvs_disp (xdir/quer) for the two folder
    # 2. make sure that x values of dUdLvs_disp are the same
    # 3. make difference of of both files
    # 4. mae folder parametrization_3_minus_0_effect_tox
    print "     foldernrs:",foldernrs
    f1 = glob.glob(str(foldernrs[0])+"_*")
    f2 = glob.glob(str(foldernrs[1])+"_*")
    if len(f1) != 1:
        print "     f1:",f1
        sys.exit("len(f1) is not 1")
    if len(f2) != 1:
        print "     f2:",f2
        sys.exit("len(f2) is not 1")
    f1 = f1[0]
    f2 = f2[0]
    f1q_path = f1+"/disp_dudl/quer/dUdLvs_disp"
    f1x_path = f1+"/disp_dudl/xdir/dUdLvs_disp"
    f2q_path = f2+"/disp_dudl/quer/dUdLvs_disp"
    f2x_path = f2+"/disp_dudl/xdir/dUdLvs_disp"
    if os.path.isfile(f1q_path) != True or \
       os.path.isfile(f1x_path) != True or \
       os.path.isfile(f2q_path) != True or \
       os.path.isfile(f2x_path) != True:
           sys.exit("f1q or f1x or f2q or f2x not available")

    f1q = np.loadtxt(f1q_path)
    f1x = np.loadtxt(f1x_path)
    f2q = np.loadtxt(f2q_path)
    f2x = np.loadtxt(f2x_path)
    if np.array_equal(f1q[:,0],f2q[:,0]) != True:
        sys.exit("11quit")
    if np.array_equal(f1x[:,0],f2x[:,0]) != True:
        sys.exit("12quit")
    fdq = np.array([f1q[:,0],f2q[:,1]-f1q[:,1]]).transpose()
    fdx = np.array([f1x[:,0],f2x[:,1]-f1x[:,1]]).transpose()

    print "     f1:",f1
    print "     f2:",f2
    f1lon = f1.split("_LON_")[-1].split("_LON2_")[0]
    f2lon = f2.split("_LON_")[-1].split("_LON2_")[0]
    f1lon2 = f1.split("_LON2_")[-1].split("_TOX_")[0]
    f2lon2 = f2.split("_LON2_")[-1].split("_TOX_")[0]
    f1tox = f1.split("_TOX_")[-1].split("_TOY_")[0]
    f2tox = f2.split("_TOX_")[-1].split("_TOY_")[0]
    f1toy = f1.split("_TOY_")[-1].split("_TOZ_")[0]
    f2toy = f2.split("_TOY_")[-1].split("_TOZ_")[0]
    f1toz = f1.split("_TOY_")[-1].split("_TOZ_")[-1]
    f2toz = f2.split("_TOY_")[-1].split("_TOZ_")[-1]
    print "     f1lon:",f1lon
    print "     f2lon:",f2lon
    print "     f1lon2:",f1lon2
    print "     f2lon2:",f2lon2
    print "     f1tox:",f1tox
    print "     f2tox:",f2tox
    print "     f1toy:",f1toy
    print "     f2toy:",f2toy
    print "     f1toz:",f1toz
    print "     f2toz:",f2toz

    ##### get diffs:
    diff_lon = None
    diff_lon2 = None
    diff_tox = None
    diff_toy = None
    diff_toz = None
    diff_lon_str  = ""
    diff_lon2_str = ""
    diff_tox_str  = ""
    diff_toy_str  = ""
    diff_toz_str  = ""
    if f1lon != f2lon:
        diff_lon_str = "LON_"+f1lon+"_vs_"+f2lon #+"_"
    if f1lon2 != f2lon2:
        diff_lon2_str = "LON2_"+f1lon2+"_vs_"+f2lon2 #+"_"
    if f1tox != f2tox:
        diff_tox_str = "TOX_"+f1tox+"_vs_"+f2tox #+"_"
    if f1toy != f2toy:
        diff_toy_str = "TOY_"+f1toy+"_vs_"+f2toy #+"_"
    if f1toz != f2toz:
        diff_toz_str = "TOZ_"+f1toz+"_vs_"+f2toz #+"_"

    diffstr = "DIFF_"+str(foldernrs[0])+"_"+str(foldernrs[1])+"_"+diff_lon_str+diff_lon2_str+diff_tox_str+diff_toy_str+diff_toz_str
    print "     diffstr:",diffstr
    if os.path.isdir("diffs") != True:
        os.makedirs("diffs")
    np.savetxt("diffs/"+diffstr+"_quer",fdq)
    np.savetxt("diffs/"+diffstr+"_xdir",fdx)
    if pot.verbose:
        print utils.printblue("     get effect of TOX complete ...")
    return

def start_dispdudl(nndist = False):
    ''' creates dudlnew and dudl_vs_disp '''
    if type(nndist) == False:
        sys.exit("please provide nndist")
    folder_dispdudl = os.getcwd()
    #utils.run2("DOS_POSITIONS_auswerten.py -dudlposc")
    DOS_POSITIONS_auswerten.dudlposc()
    dudlnewfile = folder_dispdudl+"/dUdLnew"
    if os.path.isfile(dudlnewfile) != True:
        sys.exit("dudlnew was not created! in "+str(folder_dispdudl))
    print "dudlnewfile:",dudlnewfile
    dudlnew = np.loadtxt(dudlnewfile)

    dudl_vs_disp = np.array([dudlnew[:,0]+nndist,dudlnew[:,6]]).transpose()
    np.savetxt(folder_dispdudl+"/dUdLvs_disp",dudl_vs_disp,fmt="%.3f %.4f")
    return dudlnew, dudl_vs_disp

def create_job_parametrization(jobpath_parametrization = False, lon = False, lon2 = False, tox = False, toy = False, toz = False, \
        lonallp = False, toxallp = False, toyallp = False, tozallp = False, \
        verbose = False):
    if os.path.isdir(jobpath_parametrization) != True:
        os.makedirs(jobpath_parametrization)
    if type(lon) != bool:
        shutil.copy2(lon,jobpath_parametrization)
    if type(lon2) != bool:
        shutil.copy2(lon2,jobpath_parametrization)
    if type(tox) != bool:
        shutil.copy2(tox,jobpath_parametrization)
    if type(toy) != bool:
        shutil.copy2(toy,jobpath_parametrization)
    if type(toz) != bool:
        shutil.copy2(toz,jobpath_parametrization)
    if type(lonallp) != bool:
        shutil.copy2(lonallp,jobpath_parametrization)
    if type(toxallp) != bool:
        shutil.copy2(toxallp,jobpath_parametrization)
    if type(toyallp) != bool:
        shutil.copy2(toyallp,jobpath_parametrization)
    if type(tozallp) != bool:
        shutil.copy2(tozallp,jobpath_parametrization)


    if verbose:
        print " lon  copied:",lon
        print " lon2 copied:",lon2
        print " tox  copied:",tox
        print " toy  copied:",toy
        print " toz  copied:",toz
        print " lonallp    :",lonallp
        print " toxallp    :",toxallp
        print " toyallp    :",toyallp
        print " tozallp    :",tozallp
    return


def get_parametrizationfiles_oldpath(jobpath, folder_displacement_direction_all):
    ''' returnes lonfile lon2file tox toy toz file which were used in tdi run for parametrization '''
    i = jobpath
    #print "i:",i
    lon = i.split("__LON_")
    #print "lon:",lon,len(lon)
    #print "xxx;",lon[-1].split("_")
    #print "xxx;",lon[-1].split("_")[0]
    lonfrom = lon[-1].split("_")[0]
    #print "lonfrom:",lonfrom
    lon2 =  lon[-1].split("_LON2_")
    lon2from = lon2[-1].split("_")[0]
    if lonfrom != lon2from:
        print "lon     :",lon
        print "lon2    :",lon2
        print "lonfrom :",lonfrom
        print "lon2from:",lon2from
        sys.exit("lonfrom != lon2from")
    #print "lon2:",lon2,len(lon2)
    tox =  lon2[-1].split("_TOX_")
    toxfrom = tox[-1].split("_")[0]
    #print "tox:",tox
    #print "tox:",tox[0]
    #print "tox:",tox[0].split("kp_")
    toy =  tox[-1].split("_TOY_")
    #print "toy:",toy
    toz =  toy[-1].split("_TOZ_")
    #print "toz:",toz
    lon = lon2[0].split("kp_")[-1]
    lon2 = tox[0].split("kp_")[-1]
    tox = toy[0].split("kp_")[-1]
    toy = toz[0].split("kp_")[-1]
    toz = toz[1].split("kp_")[-1].split("/")[0]
    #print "-->lon: ",lon
    #print "-->lon2:",lon2
    #print "-->tox: ",tox
    #print "-->toy: ",toy
    #print "-->toz: ",toz
    #print "lonfrom:",lonfrom
    #print "toxrom:",toxfrom
    lonfolder = False
    tofolder = False
    for j in folder_displacement_direction_all:
        check = j.split("/")[-1]
        if lonfrom in check:
            lonfolder = j
        if toxfrom in check:
            tofolder = j
    #print folder_displacement_direction_all[0].split("/")[-1]
    #print folder_displacement_direction_all[1].split("/")[-1]
    #lonfolder = folder_displacement_direction_all[0]
    #tofolder = folder_displacement_direction_all[1]
    if type(lonfolder) == bool:
        sys.exit("lonfolderprob")
    #if type(toxfolder) == bool:
    #        sys.exit("toxfolderprob")

    #print "lonfolder;",lonfolder
    #print "tofolder;",tofolder
    lonf = glob.glob(lonfolder+"/nnforces/lon*"+lon+"*")
    lon2f =  glob.glob(lonfolder+"/nnforces/lon*"+lon2+"*")
    lonallf = glob.glob(lonfolder+"/nnforces/lon_1nn_all")
    if len(lonallf) != 1:
        sys.exit("len lonallf not 1")
    lonallf = lonallf[0]
    if os.path.isfile(lonallf) != True:
        sys.exit("lonallf does not exist")
    if len(lonf) != 1 or len(lon2f) != 1:
        sys.exit("len(lonf) or len(lon2f)  not 1")
    lonf = lonf[0]
    lon2f = lon2f[0]

    def get_toxyzf(string='tox',toxyz=tox, infolder = False):
        if toxyz == "NONE":
            return False, False
        toxf =  glob.glob(infolder+"/nnforces/"+string+"*"+toxyz+"*")
        toxallf =  glob.glob(infolder+"/nnforces/"+string+"_1nn_all")
        if len(toxallf) != 1:
            sys.exit("len(toxallf) not 1")
        toxallf = toxallf[0]
        if os.path.isfile(toxallf) != True:
            sys.exit("toxyzallf does not exist")

        #print "TOXYZF:",toxf
        #print "toxyz:",toxyz
        if len(toxf) > 1:
            sys.exit("len tox 1")
        if toxyz != "NONE" and len(toxf) != 1:
            syslexti("len tox 2")
        if len(toxf) == 1:
            toxf = toxf[0]
        return toxf,toxallf

    toxf, toxallf = get_toxyzf(string='tox',toxyz=tox,infolder = tofolder)
    toyf, toyallf = get_toxyzf(string='toy',toxyz=toy,infolder = tofolder)
    tozf, tozallf = get_toxyzf(string='toz',toxyz=toz,infolder = tofolder)
    #print ""
    #print "------"
    #print " lonf:",lonf, len(lonf)
    #print " lon2f:",lon2f, len(lon2f)
    #print " toxf:",toxf
    #print " toyf:",toyf
    #print " tozf:",tozf
    return lonf, lon2f, toxf, toyf, tozf, lonallf, toxallf, toyallf, tozallf


class pot( object ):
    '''
        General class to get/create pot
    '''
    def __init__( self ):
        '''
        '''
        #############################################################################################
        # Start editing here
        #############################################################################################
        self.element                    = False     # string like Al or Ir
        self.sc_create                  = False
        self.forparametrization_sc      = False
        self.forparametrization_stringadd = ""
        self.alat                       = False  # for pot_parametrize
        self.folder_base                = "/Users/glensk/Dropbox/Understanding_distributions/"
        self.folder_base_jobvorlage     = self.folder_base+"/jobvorlage_all/"
        self.utilsfolder                = "/Users/glensk/Thermodynamics/utilities/"
        self.folder_job_tdi_base        = "/u/aglen/Understanding_distributions/ti/"  # Where are jobs created
        self.folder_job_tdi_base        = "/Users/glensk/Dropbox/Understand_distributions/ti_2/"  # Where are jobs created
        self.calculate_energy_forces_input = []
        self.jobname_param              = ""
        self.folder_displacement_from   = False   # "displacements" oder "displacements_dense"
        self.folder_displacement_from = "displacements_dense"
        self.folder_displacement_from = "displacements_"
        self.structure                  = "fcc"                     # necessary for EqCoords
        self.steps                      = 1000
        #self.piarametrizatioin_folder = False
        self.verbose                    = False
        self.forparametrization_displacement = False  # means get all (xdir, quer, 3nnd, 4nnd)
        self.parametrize_nnforcesnpz    = False
        self.parametrize_data_pkl       = False
        self.parametrize_fml_npz        = False
        self.analyze                    = False
        self.analyze_qubus              = False
        self.create_qubus_npz           = False
        self.dont_load_pkl_npz          = False
        #############################################################################################
        # Stop editing here
        #############################################################################################
        if os.path.isdir(self.folder_base) != True:
            sys.exit("ERROR: folder_base (in pot_info_startjob.py) does not exist:"+str(self.folder_base))
        if os.path.isdir(self.folder_base+"/"+self.folder_displacement_from) != True:
            sys.exit("ERROR: folder_displacements_from does not exist:"+str(self.folder_displacement_from))

    def init_variables( self, kpstringoptional = False ):
        print ""
        print utils.printblue(self.element+" init_variables for createjob ...")

        ##################################################################################################
        # get alat for the sc we are going to create
        ##################################################################################################
        if type(self.alat) == bool:
            self.alat = self.get_element_2x2x2sc_info(self.element, self.sc_create)[0]
        if type(self.alat) == bool:
            print "self.sc_create:",self.sc_create
            print "self.alat:",self.alat
            print "self.element:",self.element
            sys.exit("did not find self.alat v max; self.alat:False")
        if self.verbose == True:
            print "     self.alat:",self.alat
        self.scalat = self.sc_create*self.alat
        self.cell = np.array([[self.scalat,0.0,0.0],[0.0,self.scalat,0.0],[0.0,0.0,self.scalat]])
        if self.verbose == True:
            print "     self.scalat:",self.scalat
            #print "     self.cell:",self.cell

        self.nndist = self.alat/np.sqrt(2.)
        if self.verbose:
            print "     self.nndist:",self.nndist

        ##################################################################################################
        # check if self.jobvorlage_createfolder is available
        ##################################################################################################
        self.sc_create_string = str(self.sc_create)+"x"+str(self.sc_create)+"x"+str(self.sc_create)
        self.jobvorlage_createfolder = self.folder_base_jobvorlage+"/"+self.sc_create_string+"sc_"+self.element+"_"+str(self.alat)
        if self.verbose:
            print "     self.jobvorlage_createfolder   :",self.jobvorlage_createfolder
        if os.path.isdir(self.jobvorlage_createfolder) != True:
            sys.exit(self.jobvorlage_createfolder + " self.jobvorlage_createfolder does not exist 3 !")


        ###########################################################################################
        # general checks
        ###########################################################################################

        utils.isfile(self.jobvorlage_createfolder+"/HesseMatrix_sphinx")
        utils.isfile(self.jobvorlage_createfolder+"/INCAR")
        utils.isfile(self.jobvorlage_createfolder+"/KPOINTS")

        self.POTCAR = "/Users/glensk/Thermodynamics/vasp_potentials/PAW-GGA-PBE_vasp4.6__DO-NOT-USE--TAKE-NEW-POTCARs-FROM-vasp5.2-INSTEAD/"+self.element+"/POTCAR"
        utils.isfile(self.POTCAR)

        self.eqcoordsfile_job = self.utilsfolder+"/"+self.structure+"/EqCoords_direct_"+self.structure+"_"+self.sc_create_string+"sc"
        utils.isfile(self.eqcoordsfile_job)

        # init self.param
        self.param = pot_parametrize.forcesneighbors()


        #self.eqcoordsfile_getpot = self.utilsfolder+"/"+self.structure+"/EqCoords_direct_"+self.structure+"_"+self.forparametrization_sc_string+"sc"
        #if os.path.isfile(self.eqcoordsfile_getpot) != True:
        #    sys.exit(self.eqcoordsfile_getpot+" self.eqcoordsfile_getpot does not exist 2 !")

        #if type(self.folder_displacement_from) == bool:
        #    sys.exit("please define self.folder_displacement_from, e.g. \"displacements_dense\"")

        #self.folder_displacement = self.folder_base+"/"+self.folder_displacement_from+"/"+self.element
        ##self.folder_displacement = self.folder_base+"/displacements/"+self.element
        #if os.path.isdir(self.folder_displacement) != True:
        #    sys.exit("self.folder_displacement: "+self.folder_displacement +" does not exist 1 !")

        #if type(kpstringoptional) == bool:
        #    searchfor1 = self.folder_displacement+"/"+self.forparametrization_sc_string+"sc_quer_*"
        #    searchfor2 = self.folder_displacement+"/"+self.forparametrization_sc_string+"sc_xdir_*"
        #    self.folder_displacement_direction = glob.glob(searchfor1) + glob.glob(searchfor2)
        #else:
        #    searchfor1 =self.folder_displacement+"/"+self.forparametrization_sc_string+"sc_quer_"+kpstringoptional+"*"
        #    searchfor2 =self.folder_displacement+"/"+self.forparametrization_sc_string+"sc_xdir_"+kpstringoptional+"*"
        #    self.folder_displacement_direction = \
        #        glob.glob(searchfor1) + \
        #        glob.glob(searchfor2)

        #    if len(self.folder_displacement_direction) == 0:
        #        searchfor1 = self.folder_displacement+"/"+self.forparametrization_sc_string+"sc_quer_"+"*"
        #        searchfor2 = self.folder_displacement+"/"+self.forparametrization_sc_string+"sc_xdir_"+"*"
        #        self.folder_displacement_direction = \
        #        glob.glob(searchfor1) + \
        #        glob.glob(searchfor2)

        #if len(self.folder_displacement_direction) != 2:
        #    print "searched for:",searchfor1
        #    print "            :",searchfor2
        #    print "found self.folder_displacement_direction:",self.folder_displacement_direction
        #    sys.exit("no dispfolder found in self.folder_displacement: "+self.folder_displacement)
        #self.folder_displacement_direction_all = copy.copy(self.folder_displacement_direction)
        #for jkl in self.folder_displacement_direction_all:
        #    if os.path.isfile(jkl+'/EqCoords_direct') != True:
        #        shutil.copyfile(self.eqcoordsfile_getpot,jkl+'/EqCoords_direct')
        #if os.path.isfile(self.folder_displacement+'/EqCoords_direct') != True:
        #    shutil.copyfile(self.eqcoordsfile_getpot,self.folder_displacement+'/EqCoords_direct')
        #if self.verbose:
        #    print "     self.forparametrization_sc_string:",self.forparametrization_sc_string
        #    print "     self.eqcoordsfile_getpot:",self.eqcoordsfile_getpot
        #    print "     self.folder_displacement:",self.folder_displacement
        #    print "     self.folder_displacement_direction_all:"
        #    for jkl in self.folder_displacement_direction_all:
        #        print "             -> ",jkl

        # get self.alat, scalat


        # define DOS, DOScut, is of course nice, but what happens if nnforces.npz not available?
        # this should only be loaded when nnforzes.npz is available!
        #self.param = pot_parametrize.forcesneighbors()
        #self.param.loadforces(a = self.alat)
        #self.param.infodos(jobvorlage = self.jobvorlage_createfolder, verbose = self.verbose)


        ###########################################################################################
        # POSITIONs  --> in self.param_quer.P or self.param_xdir.P
        ###########################################################################################
        if False:
            for pos in self.folder_displacement_direction_all:
                #print "pos:",pos
                positions = pos+"/POSITIONs"
                #if os.path.isfile(positions) != True:
                    #sys.exit(positions+" positionsfile does not exist")


        ###########################################################################################
        # get DOSlonpy
        ###########################################################################################
        if False:
            print "kk:",self.jobvorlage_createfolder_vecnormlon_file
            data = np.loadtxt(self.jobvorlage_createfolder_vecnormlon_file)[2:]
            #data = np.loadtxt("/Users/glensk/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/tveclonall.dat")
            sys.exit()
            dist_space = np.linspace( min(data), max(data), 500 )
            fact = 1.0
            def my_kde_bandwidth(obj, fac=fact):
                """We use Scott's Rule, multiplied by a constant factor."""
                return np.power(obj.n, -1./(obj.d+4)) * fac

            from scipy.stats.kde import gaussian_kde
            kde = gaussian_kde( data ,bw_method=my_kde_bandwidth)
            #kde = gaussian_kde( data )
            np.savetxt(self.jobvorlage_createfolder+"/DOS_dfn_1.0_py"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())
            #np.savetxt(self.jobvorlage_createfolder+"/DOSlonpy"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())
            #np.savetxt("/Users/glensk/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/t"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())


        ###########################################################################################
        # get lon xidr/ quer und diff xdir - quer ( dann auch alles im Dcut bereich)
        ###########################################################################################
        if False:
            self.disp_lon_all_quer_path = self.folder_displacement_direction_all[0]+"/nnforces/lon_1nn_all"
            self.disp_lon_all_xdir_path = self.folder_displacement_direction_all[1]+"/nnforces/lon_1nn_all"
            if os.path.isfile(self.disp_lon_all_quer_path) != True:
                sys.exit(self.disp_lon_all_quer_path+" lonfile missing")
            if os.path.isfile(self.disp_lon_all_xdir_path) != True:
                sys.exit(self.disp_lon_all_xdir_path+" lonfile missing")
            self.disp_lon_all_quer = np.loadtxt(self.disp_lon_all_quer_path)
            self.disp_lon_all_xdir = np.loadtxt(self.disp_lon_all_xdir_path)
            if self.disp_lon_all_quer.shape != self.disp_lon_all_xdir.shape:
                print "self.disp_lon_all_quer.shape:",self.disp_lon_all_quer.shape
                print "self.disp_lon_all_xdir.shape:",self.disp_lon_all_xdir.shape
                sys.exit("self.disp_lon_all_quer.shape != self.disp_lon_all_xdir.shape")
            f1,f2 = utils.return_fuctions_on_same_grid(self.disp_lon_all_quer, self.disp_lon_all_xdir)
            self.disp_lon_all_diff = np.array([f1[:,0],f1[:,1]-f2[:,1]]).transpose()
            self.disp_lon_all_diff_file = self.folder_displacement_direction_all[0]+"/nnforces/lon_1nn_all_quer_min_xdir"
            if self.verbose:
                print "     self.disp_lon_all_diff_file:",self.disp_lon_all_diff_file
            #print self.disp_lon_all_diff
            if False:
                np.savetxt(self.disp_lon_all_diff_file,self.disp_lon_all_diff)

            # create diffDOScut
            self.disp_lon_all_diffcut = utils.cut_function_at_DOS(self.disp_lon_all_diff,self.jobvorlage_createfolder_DOSloncut)
            self.disp_lon_all_diffcut_file = self.disp_lon_all_diff_file+"_cut"
            self.disp_lon_all_diffcutshifted = np.array([self.disp_lon_all_diffcut[:,0]-self.nndist,self.disp_lon_all_diffcut[:,1]]).transpose()
            self.disp_lon_all_diffcutshifted_file = self.disp_lon_all_diffcut_file+ "shifted0"


        if False:
            np.savetxt(self.disp_lon_all_diffcut_file,self.disp_lon_all_diffcut)
            np.savetxt(self.disp_lon_all_diffcutshifted_file,self.disp_lon_all_diffcutshifted)

        print utils.printblue("     init_variables complete ...")
        print ""
        return

        # init of parametrization, do not load forces yet!!
        lastpartdata = self.folder_displacement_direction[0].split("/")[-1].replace("_quer_","_data_").replace("_xdir_","_data_")
        lastpartfit = self.folder_displacement_direction[0].split("/")[-1].replace("_quer_","_fit__").replace("_xdir_","_fit__")

        self.param.pkl = "/".join(self.folder_displacement_direction[0].split("/")[:-1])+'/'+lastpartdata+'.pkl'
        #self.param.fit_file =  "/".join(self.folder_displacement_direction[0].split("/")[:-1])+'/'+lastpartfit+'.pkl'
        #print self.param.fit_file
        if os.path.isfile(self.param.pkl) == True: # and os.path.isfile(self.param.fit_file):
            print utils.printred("loading "+self.param.pkl)
            self.param.loadDatapd(self.param.pkl)
            #self.param.loadFitpd(self.param.fit_file)
            print utils.printgreen("     loaded data & fit")
        else:
            print utils.printred("     need to parametrize displacements!")
        print utils.printblue("     init_variables complete")
        print ""
        return

    def get_eqcoords_path(self, structure, sc, copy_to_folder = False):
        eqcoordsfile_path = self.utilsfolder+"/"+structure+"/EqCoords_direct_"+structure+"_"+str(sc)+"x"+str(sc)+"x"+str(sc)+"sc"
        utils.isfile(eqcoordsfile_path)
        if type(copy_to_folder) == str:
            if os.path.isdir(copy_to_folder) == True:
                shutil.copy2(eqcoordsfile_path,copy_to_folder+"/EqCoords_direct")
            else:
                sys.exit(copy_to_folder+" does not exist!")
        return eqcoordsfile_path

    def get_alats( self, folder ):
        alats = glob.glob(folder+"/[0-9.]*Ang_*")
        alatsstringlist = []
        for a in alats:
            #print a.split("Ang")
            disp = a.split(folder)[-1]
            alat = disp.split("Ang")[0]
            #disps = disp.split("Ang")[1]
            alatsstring = alat.split("/")[-1]
            alatsstringlist.extend([alatsstring])
        self.alatsstringlist = list(set(alatsstringlist))
        return self.alatsstringlist

    #def load_dataframe(self, sc):
    #    ''' currently unused '''
    #    print ""
    #    print utils.printblue("load_dataframe ...")
    #    forparametrization_sc_string = str(sc)+"x"+str(sc)+"x"+str(sc)+"sc"
    #    if self.verbose:
    #        print "     forparametrization_sc_string:",forparametrization_sc_string
    #    forparametrization_kpoint_string = self.get_element_2x2x2sc_info(self.element, sc)[1]
    #    if self.verbose:
    #        print "     forparametrization_kpoint_string:",forparametrization_kpoint_string

    #    search = self.folder_displacement_element+'/'+forparametrization_sc_string+"_*_"+forparametrization_kpoint_string
    #    if self.verbose:
    #        print "     search:",search
    #    folder_displacement_direction_all =  glob.glob(search)
    #    if self.verbose:
    #        for folder in folder_displacement_direction_all:
    #            print utils.printgreen("     folder    : -->>>> "+folder) #,"type(DOScut):",type(DOScut)
    #            if os.path.isfile(folder+"/EqCoords_direct") != True:
    #                self.get_eqcoords_path(self.structure, sc, folder)
    #    # if we let the parametrization run we want to update those files!
    #    #if os.path.isfile(self.param.pkl) == True and os.path.isfile(self.param.fit_file):
    #    #    self.param.loadDatapd(self.param.pkl)
    #    #    self.param.loadFitpd(self.param.fit_file)
    #    #    print utils.printblue("     pot_parametrize_function complete...")
    #    #    return
    #    print utils.printblue("     load_dataframe complete...")
    #    return

    def pot_load_parametrize_function_get_folder( self ):
        ''' returns self.folder_displacement_direction_all which is the folder for the parametrization
            e.g. /Users/glensk/Understand_distributions//displacements_dense/Ir/2x2x2sc_xdir_3x3x3kp
            the direction is taken from self.forparametrization_displacement '''
        print ""
        print utils.printblue("     pot_load_parametrize_function_get_folder ...")

        self.forparametrization_sc_string = str(self.forparametrization_sc)+"x"+str(self.forparametrization_sc)+"x"+str(self.forparametrization_sc)+"sc"
        if self.verbose:
            print "     self.forparametrization_sc_string     :",self.forparametrization_sc_string
        if type(self.forparametrization_kpoint_string) == bool:
            self.forparametrization_kpoint_string = self.get_element_2x2x2sc_info(self.element, self.forparametrization_sc)[1]
        if self.verbose:
            print "     self.forparametrization_kpoint_string :",self.forparametrization_kpoint_string
            print "     self.forparametrization_stringadd     :",self.forparametrization_stringadd

        #self.forparametrization_displacement = False  # means get all (xdir, quer, 3nnd, 4nnd)

        if self.forparametrization_displacement == False:
            forparametrization_displacement = [ "*" ]
        else:
            forparametrization_displacement = self.forparametrization_displacement
        if self.verbose:
            print "     forparametrization_displacement       :",forparametrization_displacement
        #forparametrization_displacement = '{xdir,quer}'

        #                       Ir                                  2x2x2sc                            */{xdir,quer,3nnd,4nnd}             3x3x3kp
        self.folder_displacement_direction_all = []
        for i in forparametrization_displacement: # {"*",quer,xdir,midd,3nnd,...}
            # this is the search pattern !!!!!!
            search = self.folder_displacement_element+'/'+self.forparametrization_sc_string+"_"+str(self.alat)+"Ang_"+i+"_"+self.forparametrization_kpoint_string+self.forparametrization_stringadd
            if self.verbose:
                print "     search                                :",search
            self.folder_displacement_direction_all = self.folder_displacement_direction_all + glob.glob(search)
            #if self.verbose:
            #    print "     self.folder_displacement_direction_all:",self.folder_displacement_direction_all
        #print ""
        self.folder_displacement_direction_all = utils.sort_list_to_sortlist(self.folder_displacement_direction_all,['quer','xdir'])
        if len(self.folder_displacement_direction_all) == 0:
            sys.exit("pot.folder_displacement_direction_all is an empty list! there semms to be no folder "+str(search)+" !")

        if self.verbose:
            for i in self.folder_displacement_direction_all:
                print "     self.folder_displacement_direction_all:",i
        # if we let the parametrization run we want to update those files!
        #if os.path.isfile(self.param.pkl) == True and os.path.isfile(self.param.fit_file):
        #    self.param.loadDatapd(self.param.pkl)
        #    self.param.loadFitpd(self.param.fit_file)
        #    print utils.printblue("     pot_parametrize_function complete...")
        #    return
        return

    def pot_load_parametrize_function( self, parametrize = False ):
        ''' this function loads and/or parametrizes the T=0K displacements;
        it goes through displacements folders (xdir,quer) and fits forces to lon,tox functions;
        if parametrize == Fals forces are just loaded from pkl file;
        if parametrize == True the parametrization is performed to T=0K is performed'''
        print ""
        print utils.printblue("pot_load_parametrize_function ... 0")

        ###################################################################################
        # get some info about the job
        ###################################################################################
        self.pot_load_parametrize_function_get_folder()

        for folder in self.folder_displacement_direction_all:
            print utils.printgreen("     folder    : -->>>> "+folder) #,"type(DOScut):",type(DOScut)
        print ""

        ###################################################################################
        # got through xdir,quer,3nnd,4nnd and load forces stuff (and parametrize if wished)
        ###################################################################################
        # folder    : -->>>> /Users/glensk/Understand_distributions//displacements_dense/Ir/2x2x2sc_quer_3x3x3kp
        # folder    : -->>>> /Users/glensk/Understand_distributions//displacements_dense/Ir/2x2x2sc_xdir_3x3x3kp

        print utils.printblue("pot_load_parametrize_function ... 1")
        if type(self.alat) == bool:
            self.alatsstringlist = self.get_alats(folder)
        else:
            self.alatsstringlist = [ self.alat ]

        self.disp = pot_parametrize.forcesneighbors_all_disp()
        self.disp.folder_base = self.folder_base
        self.disp.folder_displacement_direction_all = self.folder_displacement_direction_all
        self.disp.alatsstringlist = self.alatsstringlist
        self.disp.pkl = self.folder_displacement_element+'/'+self.forparametrization_sc_string+"_"+str(self.alat)+"Ang_data_"+self.forparametrization_kpoint_string+self.forparametrization_stringadd+'.pkl'
        self.disp.npz = self.folder_displacement_element+'/'+self.forparametrization_sc_string+"_"+str(self.alat)+"Ang_data_"+self.forparametrization_kpoint_string+self.forparametrization_stringadd+'.npz'


        print utils.printblue("pot_load_parametrize_function ... 2")
        self.disp.parametrize_nnforcesnpz = self.parametrize_nnforcesnpz
        self.disp.parametrize_data_pkl = self.parametrize_data_pkl
        self.disp.parametrize_fml_npz = self.parametrize_fml_npz
        self.disp.create_qubus_npz = self.create_qubus_npz

        print utils.printblue("pot_load_parametrize_function ... 3")
        if self.verbose == True:
            print "     self.dont_load_pkl_npz:",self.dont_load_pkl_npz
        if self.dont_load_pkl_npz == True:
            pass
        else: # if self.dont_load_pkl_npz == False
            self.disp.alle.struct = pot.structure
            self.disp.alle.verbose = self.verbose
            self.disp.load_all_pkl_npz()  # form pot_parametrize.py
            # self.disp.load_all_pkl_npz() loads in line 2331: self.alle.creatennforcesnpz(folder = folder, a = a)

        print utils.printblue("pot_load_parametrize_function ... 4")
        if self.analyze == True:
            self.disp.analyze_restforces()
        if self.analyze_qubus == True or self.create_qubus_npz:
            self.disp.analyze_qubus()

        print utils.printblue("pot_load_parametrize_function ... 5")
        #self.disp.parametrize_tox()

        print utils.printblue("     pot_load_parametrize_function complete...")
        print ""
        return

    def plot_lon(self, shell = 1, out = 'd'):
        ''' plot the lon datapoints for a certain shell
            out:{\'a\' == all,\'d\' == dos, \'dc\' == doscut '''
        import matplotlib.pyplot as plt
        #from matplotlib import rc
        #plt.rc('text', usetex=True)
        #plt.rc('font', family='serif')
        plt.clf()
        xliml = self.param.nndist[shell]-1.0
        xlimr = self.param.nndist[shell]+1.0
        if type(self.param.DOSlon[shell]) != bool:
            xliml = self.param.DOSlon[shell][:,0].min()
            xlimr = self.param.DOSlon[shell][:,0].max()
        if out == 'a':
            param_xdir = self.param_xdir.lonall
            param_quer = self.param_quer.lonall
        if out == 'd':
            param_xdir = self.param_xdir.lonallDOS
            param_quer = self.param_quer.lonallDOS
        if out == 'dc':
            param_xdir = self.param_xdir.lonallDOScut
            param_quer = self.param_quer.lonallDOScut

        plotfunc = param_xdir
        plt.plot(plotfunc[shell][:,0],plotfunc[shell][:,1],'b.-',label='xdir')
        x1,x2,y1,y2 = utils.plot_find_ylim_from_xlim(plotfunc[shell],xliml,xlimr)
        plotfunc = param_quer
        plt.plot(plotfunc[shell][:,0],plotfunc[shell][:,1],'r.-',label='quer')
        x3,x4,y3,y4 = utils.plot_find_ylim_from_xlim(plotfunc[shell],xliml,xlimr)

        plt.xlim(np.array([x1,x3]).min(),np.array(x2,x4).max())
        plt.ylim(np.array([y1,y3]).min(),np.array(y2,y4).max())
        dy = np.array(y2,y4).max() - np.array([y1,y3]).min()
        if type(self.param.DOSlon[shell]) != bool:
            dosy = self.param.DOSlon[shell][:,1].max()
            #print "dy:",dy
            #print "dosy:",dosy
            plt.plot(self.param.DOSlon[shell][:,0],
                (self.param.DOSlon[shell][:,1])/dosy*dy+np.array([y1,y3]).min(),
                'g--',label="DOS")
        plt.legend(loc='best')
        plt.plot(np.array([self.param_quer.nndist[shell]]),np.array([0.0]),'ro',ms=10)
        plt.xlabel('Internuclear distance (\AA)', fontsize=14)
        plt.ylabel('Internuclear force (eV/\AA)')
        plt.title(str(shell)+' neighbor (Element: '+str(self.param.atom.symbol[0])+')')
        plt.grid(True)
        plt.ion()
        return

    #def get_parametrization_folder(self, disp = False, stringadd = False):
    #    ''' returnes something like:
    #            2x2x2sc_quer_3x3x3kp or
    #            2x2x2sc_xdir_3x3x3kp
    #            '''
    #    if type(disp) == bool or type (stringadd) == bool:
    #        return False
    #    #if self.verbose:
    #    #    print ""
    #    #    print "get_parametrization_folder ..."
    #    searchfor = self.folder_displacement+"/"+self.forparametrization_sc_string+"sc_"+disp+"_"+stringadd
    #    #if self.verbose:
    #    #    print "     searchfor:",searchfor
    #    dispfolder = glob.glob(searchfor)
    #    if len(dispfolder) != 1:
    #        print ""
    #        print utils.printred("get_parametrization_folder ...")
    #        print "-----------ERROR since len(dispfolder) != 1 ----------------"
    #        print 'searchfor                = self.folder_displacement+"/"+self.forparametrization_sc_string+"sc_"+disp+"_"+stringadd'
    #        print "self.folder_displacement =",self.folder_displacement
    #        print "self.forparametrization_sc_string    =",self.forparametrization_sc_string
    #        print "disp (xdir/quer)         =",disp
    #        print "searchfor                =",searchfor
    #        print ""

    #        print "self.sc_string_arams:",self.forparametrization_sc_string
    #        print "stringadd:",stringadd
    #        print "dispfolder:",dispfolder
    #        if len(dispfolder) == 0:
    #            sys.exit("expected one dispfolder, got none")
    #        if len(dispfolder) != 0:
    #            sys.exit("expected one dispfolder, got several")
    #    self.folder_displacement_direction = dispfolder[0]
    #    #print "dispfolder:",dispfolder
    #    if os.path.isdir(self.folder_displacement_direction) != True:
    #        sys.exit(self.folder_displacement_direction+" not found self.folder_displacement_direction")
    #    #print "ddd:",self.folder_displacement_direction
    #    return self.folder_displacement_direction


    def get_parametrization_parameters_from_file(self,filename): #,contr, parametrization_string):
        ''' contr an be 'lon', 'lon2', 'tox', 'toy', 'toz'
            parametrization_string: all_poly_9thorder , all_rlv_0.25_mc1 ,  all_poly_1storder

        '''
        #print 'filename;',filename,type(filename)
        if type(filename) == bool:
            return filename,filename
        #print "still in"
        par_lon_file = filename
        ## get parametrization from filename of file
        par_lon_params = par_lon_file.split("/")
        pars1 = par_lon_params[-1]
        #print "pars1:",pars1
        pars2 = pars1.split("___")
        if len(pars2) != 2:
            print "pars2:",pars2
            sys.exit("expected two parts of pars2!")
        par_lon_params = pars2[-1]
        #print "dd:",par_lon_params
        #print "par_lon_params:",par_lon_params
        par_lon_pot = False
        if pars1.split("_")[0] == "lon":
        #if contr == 'lon':
            #if self.verbose:
            #    print "     parstr:",parametrization_string
            if "_poly" in pars1:
                par_lon_pot = "poly"
            if "_morse" in pars1:
                par_lon_pot = "m"
            if "_mc1" in pars1:
                par_lon_pot = "mc1"
            if par_lon_pot == False:
                sys.exit("par_lon_pot not found")
        if pars1.split("_")[0] == "tox" or pars1.split("_")[0] == "toy" or pars1.split("_")[0] == "toz" or pars1.split("_")[0] == "tix" or pars1.split("_")[0] == "tiy" or pars1.split("_")[0] == "tiz":
            if "_poly" in pars1:
                par_lon_pot = "poly"
        return par_lon_params,par_lon_pot

    def pars_from_pd(self, *kw):
        ''' gets a parametrization from pandas dataframe '''
        print utils.printblue("     pars_from_pd ...")
        #pot.pars_from_pd(1, 'lonadd', 'from_2x2x2sc', 'get_kpstring', True, 'quer', 'dos', 'morse', 'polybeste' )
        self.shell      = kw[0]     # {1,2,3,4}
        self.controrig  = kw[1]     # 'lon', 'lonadd', 'tox', 'toy', 'tix', ...
        self.scstring   = kw[2]     # 'from_2x2x2sc', 'from_3x3x3sc'
        self.kpstring   = kw[3]     # 'get_kpstring' or '10x10x10kp'
        self.runornot   = kw[4]     # True or False
        self.disp       = kw[5]     # 'quer', 'xdir'
        self.rang       = kw[6]     # 'dos', 'doscut', 'all'
        self.func       = kw[7]     # 'morse', 'mc1', 'poly'
        self.add        = False     # NoEntry or 'polybeste'


        ########################################################################################
        # checks regarding functions
        ########################################################################################
        if self.controrig == 'lon' or self.controrig == 'lonadd':
            if self.shell == 1 and self.disp != 'quer':
                sys.exit("use quer disp for shell 1 lon contribution, xdir might have jumps")
            if self.shell == 2 and self.disp != 'xdir':
                sys.exit("use xdir disp for shell 2 lon contribution, quer might have jumps:")

        ########################################################################################
        # checks of input parameters
        ########################################################################################
        if 'add' in self.controrig:
            if len(kw) != 9: sys.exit("pars_from_pd() needs 9 keywords for {lon,tox,...}_add, given: "+str(kw))
            self.add = kw[8]
        else:
            if len(kw) != 8: sys.exit("pars_from_pd() needs 8 keywords, given: "+str(kw))
            self.add = ''

        if False in kw: return  # kw[2] meistens

        if self.runornot != True: sys.exit('kw[2] has to be False or True')

        if self.scstring != 'from_2x2x2sc' and self.scstring != 'from_3x3x3sc':sys.exit('kw[6] is neither \"from_2x2x2sc\" nor \"from_3x3x3sc\"')
        if self.scstring == 'from_2x2x2sc': self.scstring = "2x2x2sc"
        if self.scstring == 'from_3x3x3sc': self.scstring = "3x3x3sc"
        if self.kpstring == 'get_kpstring':
            if self.scstring == "2x2x2sc": self.kpstring = self.get_element_2x2x2sc_info(self.element, sc=2)[1]   # gets 3x3x3kp string
            if self.scstring == "3x3x3sc": self.kpstring = self.get_element_2x2x2sc_info(self.element, sc=3)[1]   # gets 3x3x3kp string
        self.param.pkl = self.folder_displacement_element+'/'+self.scstring+"_data_"+self.kpstring+'.pkl'
        #self.param.fit_file  = self.folder_displacement_element+'/'+self.scstring+"_fit__"+self.kpstring+'.pkl'
        #print "picklefile_data:",self.param.pkl
        #print "picklefile_fit :",self.param.fit_file
        self.param.loadDatapd(self.param.pkl)
        #self.param.loadFitpd(self.param.fit_file)

        self.dataframe = self.param.data
        self.contr = self.controrig
        if 'add' in kw[1]:
            self.dataframe = self.param.fit
            self.contr = self.controrig.split('add')[0]
        #print "self.contr2",self.contr
        #print "self.u"+str(self.shell)+"nn_"+self.contr, kw[1]

        #####################################################################################################
        # calculate_energy_and_forces
        #####################################################################################################
        if self.controrig == 'lon':
            append = '        self.u'+str(self.shell)+'nn_pottype          = \''+self.func+'\''
            self.calculate_energy_forces_input.append(append)
            ka = self.dataframe[self.shell,self.contr,self.disp,self.rang].ix[self.func]
            def boolkakb(kx):
                print "kab:",kx,type(kx),"== pot.dataframe["+str(self.shell)+",\""+str(self.contr)+"\",\""+str(self.disp)+"\",\""+str(self.rang)+"\"] is not defined! (bool)"
                print "rerun AUSWERTUNG with -pnnforcesnpz -ppkldata"
                sys.exit("e.g: run ~/Thermodynamics/python_thermodynamics/pot_info_startjob.py -e ni -ps 2 -vvv -pnnforcesnpz -ppkldata")
            if type(ka) == bool:
                boolkakb(ka)
            kb = utils.list_to_string(self.dataframe[self.shell,self.contr,self.disp,self.rang].ix[self.func])
            if type(kb) == bool:
                boolkakb(kb)
            append = '        self.u'+str(self.shell)+'nn_potparam         = \''+utils.list_to_string(self.dataframe[self.shell,self.contr,self.disp,self.rang].ix[self.func])+'\''
            self.calculate_energy_forces_input.append(append)
        if self.controrig == 'lonadd':
            append = '        self.u'+str(self.shell)+'nn_potadd           = \''+'poly'+'\''   # we have to specify which function we do use, the add stuff will always be poly (or spline)
            self.calculate_energy_forces_input.append(append)
            append = '        self.u'+str(self.shell)+'nn_potaddparam      = \''+utils.list_to_string(self.dataframe[self.shell,self.contr,self.disp,self.rang,self.func].ix['polybeste'])+'\''
            self.calculate_energy_forces_input.append(append)

        #####################################################################################################
        # jobname
        #####################################################################################################
        if 'add' in kw[1]:
            #self.jobname_param = self.jobname_param + str(self.shell)+"_"+self.controrig+"_"+self.scstring+"_"+self.disp+"_"+self.rang+"_"+self.func+"_"+self.add+"___"
            self.jobname_param = self.jobname_param + str(self.shell)+"_"+self.controrig+"_"+self.add+"___"
        else:
            self.jobname_param = self.jobname_param + str(self.shell)+"_"+self.controrig+"_"+self.scstring+"_"+self.disp+"_"+self.rang+"_"+self.func+"___"
        print utils.printblue("     pars_from_pd ... done")
        return

    def get_parametrization(self, shell = False, contr = False, disp = False, stringadd = False, parametrization_string = False):
        ''' disp: 'quer', 'xdir'                           # here the parametrization cames from
            stringadd: 3x3x3kp_vasp4, 10x10x10kp, ...      # here the parametrization cames from
            parametrization_string: neg_rlv_0.15_mc1       # here the parametrization cames from
            contr = [ 'lon', 'lonadd', 'lon2', 'lon2add', 'tox', 'toy', 'toz' ]
            '''
        def print_indata():
            print "shell                    :",shell
            print "contr                    :",contr
            print "disp                     :",disp
            print "stringadd                :",stringadd
            print "parametrization_string   :",parametrization_string
            return

        if self.verbose >=2:
            print_indata()

        if type(contr) == bool or type(shell) == bool:
            return
        shellstr = str(shell)

        shell_choices = [ 1,2,3,4 ]
        contr_choices = [ 'lon', 'lonadd', 'lon2', 'lon2add', 'tox', 'toy', 'toz' ]
        disp_choices = [ 'xdir', 'quer' ]
        parametrization_string_choices = [ 'doscut_morse', 'dos_mores', 'doscut_mc1',
                'dos_mc1' ]
        if contr not in contr_choices:
            sys.exit("contr: \""+str(contr)+"\" but contr has to be one of "+str(contr_choices))
        if shell not in shell_choices:
            sys.exit("shell: \""+str(shell)+"\" but shell has to be one of "+str(shell_choices))
        if disp not in disp_choices:
            sys.exit("disp: \""+str(disp)+"\" but disph has to be one of "+str(disp_choices))
        if parametrization_string not in parametrization_string_choices:
            sys.exit("parametrization_string: \""+str(parametrization_string)+"\" but disph has to be one of "+str(parametrization_string_choices))

        self._par_param_folder =  self.get_parametrization_folder(disp = disp, stringadd = stringadd)  # stringadd
        self.fitfolder = self._par_param_folder+"/nnforces/"
        if os.path.isdir(self.fitfolder) != True:
            sys.exit("self.fitfolder "+self.fitfolder+" not found")
        if self.verbose >=2:
            print ""
            print "self._par_param_folder        :",self._par_param_folder
            print "self.fitfolder           :",self.fitfolder

        #############################################################################################
        if contr in [ 'lonadd' , 'lon2add' ]:
            if contr == 'lonadd':
                contr = 'lon'
                searchfor =     self.fitfolder+"/deltasfit/"+contr+"_"+shellstr+"nn_"+parametrization_string+"*"
                searchforallp = self.fitfolder+"/deltas/"+contr+"_"+shellstr+"nn_"+parametrization_string+"*"
                contr = 'lonadd'
            if contr == 'lon2add':
                contr = 'lon'
                searchfor =     self.fitfolder+"/deltasfit/"+contr+"_"+shellstr+"nn_"+parametrization_string+"*"
                searchforallp = self.fitfolder+"/deltas/"+contr+"_"+shellstr+"nn_"+parametrization_string+"*"
                contr = 'lon2add'

        #############################################################################################
        # for par lon file
        #############################################################################################
        if contr in [ 'lon' , 'lon2', 'tox', 'toy', 'toz' ]:
            searchfor = self.fitfolder+"/"+contr+"_"+shellstr+"nn_"+parametrization_string+"*"


        self._par_file = glob.glob(searchfor)
        def printtmpstuff():
            print "@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@"*3
            print "parametrization_folder:",self._par_param_folder
            print_indata()
            print "-----------------------------------------------"
            print "searchfor:",searchfor
            print "disp:",disp
            print "stringadd:",stringadd
            print "self._par_file",self._par_file
            print "@@@@@@@@@ ERROR @@@@@@@@@@@@@@@@"*3
            sys.exit("expected one file!")
        if len(self._par_file) != 1:
            printtmpstuff()

        self._par_file = self._par_file[0]
        print "self._par_file:",self._par_file

        par_lon_params, par_lon_pot = self.get_parametrization_parameters_from_file(self._par_file)
        if contr in [ 'lonadd' , 'lon2add' ]:
            par_lon_pot = 'poly'

        if contr == 'lon':
            self.par_lon_file[shell] = self._par_file
            self.par_lon_param_folder[shell] = self._par_param_folder

        return par_lon_params,self._par_file,self._par_param_folder,par_lon_pot

    def check_calculate_energy_and_forces_parameters_with_current(self,calculate_energy_and_forces_path):
        import imp
        calculate_energy_and_forces = imp.load_source('calculate_energy_and_forces',calculate_energy_and_forces_path)
        from calculate_energy_and_forces import \
                u1nn_pot, \
                u1nn_potparam, \
                u1nn_pot2, \
                u1nn_pot2param, \
                u1nn_topx, \
                u1nn_topy, \
                u1nn_topz, \
                u1nn_tipx, \
                u1nn_tipy, \
                u1nn_tipz, \
                u2nn_pot, \
                u2nn_potparam, \
                u2nn_pot2, \
                u2nn_potparam2, \
                u2nn_topx, \
                u2nn_topy, \
                u2nn_topz, \
                u2nn_tipx, \
                u2nn_tipy, \
                u2nn_tipz, \
                u3nn_pot, \
                u3nn_potparam, \
                u3nn_pot2, \
                u3nn_potparam2, \
                save_vecs_to_file_for_DOS

        if u1nn_potparam != self.par_lon_params:
            sys.exit("55 lon")
        if u1nn_pot2param != self.par_lon2_params:
            sys.exit("55 lon2")
        if u1nn_topx != self.par_tox_params:
            sys.exit("55 tox")
        if u1nn_topy != self.par_toy_params:
            sys.exit("55 toy")
        if u1nn_topz != self.par_toz_params:
            sys.exit("55 toz")
        return


    def UNUSED_STUFF_DELETE_LATER(self):
        def params_init(self):
            import pandas as pd
            data = {'year': [2010, 2011, 2012, 2011, 2012, 2010, 2011, 2012],
                    'team': ['Bears', 'Bears', 'Bears', 'Packers', 'Packers', 'Lions', 'Lions', 'Lions'],
                    'wins': [11, 8, 10, 15, 11, 6, 10, 4],
                    'losses': [5, 8, 6, 1, 5, 10, 6, 12]}
            football = pd.DataFrame(data, columns=['year', 'team', 'wins', 'losses'])
            self.left = pd.DataFrame({'key': ['foo', 'foo'], 'lval': [1, 2]})
            print "left:"
            print self.left
            index = pd.date_range('1/1/2000', periods=8)
            print "index:"
            print index,type(index)
            s = pd.Series(np.random.randn(5), index=['a', 'b', 'c', 'd', 'e'])
            df = pd.DataFrame(np.random.randn(8, 3), index=index,columns=['A', 'B', 'C'])
            print "s:"
            print s
            print "df;"
            print df
            print "-----------------------"
            self.df = pd.DataFrame({
                'one'   : pd.Series(np.random.randn(3), index=['a', 'b', 'c']),
                'two'   : pd.Series(np.random.randn(4), index=['a', 'b', 'c', 'd']),
                'three' : pd.Series(np.random.randn(3), index=['b', 'c', 'd'])})
            indexx = [ 'shell', 'pot', 'params', 'file', 'file_orig', 'disp',  ]  # parametrization_string = doscut_morse, nn_dos_mc1, best would be [ all, dos, doscut ] and include xdir/quer
            #self.df = pd.DataFrame({
            #    'lon'       : pd.Series(['m','0.0_1.0',1, 'filepath'], index=indexx),
            #    'lon_add'   : pd.Series(np.random.randn(4), index=indexx),
            #    'lon2'      : pd.Series(np.random.randn(4), index=indexx),
            #    'lon2_add'  : pd.Series(np.random.randn(4), index=indexx)
            #    })
            print self.df


            shellsmax = 4+1  # +1 da 0 sell always False
            self.par_lon_shell       = [ shell for shell in range(shellsmax) ]  # mostly for checking purpuses
            self.par_lon_pot         = [ False for shell in range(shellsmax) ]
            self.par_lon_params      = [ False for shell in range(shellsmax) ]
            self.par_lon_file        = [ False for shell in range(shellsmax) ]
            self.par_lon_file_orig   = [ False for shell in range(shellsmax) ]
            self.par_lon_disp        = [ False for shell in range(shellsmax) ]

            self.par_lon_stringadd   = [ False for shell in range(shellsmax) ]
            self.par_lon_param_string = [ False for shell in range(shellsmax) ]
            self.par_lon_param_folder = [ False for shell in range(shellsmax) ]
            #
            self.par_lon_add_pot         = [ False for shell in range(shellsmax) ]
            self.par_lon_add_params      = [ False for shell in range(shellsmax) ]
            self.par_lon_add_file        = [ False for shell in range(shellsmax) ]
            self.par_lon_add_file_orig   = [ False for shell in range(shellsmax) ]
            self.par_lon_add_shell       = [ shell for shell in range(shellsmax) ]  # mostly for checking purpuses
            self.par_lon_add_disp        = [ False for shell in range(shellsmax) ]
            self.par_lon_add_stringadd   = [ False for shell in range(shellsmax) ]
            self.par_lon_add_param_string = [ False for shell in range(shellsmax) ]
            self.par_lon_add_param_folder = [ False for shell in range(shellsmax) ]

            self.par_lon2_pot         = [ False for shell in range(shellsmax) ]
            self.par_lon2_params      = [ False for shell in range(shellsmax) ]
            self.par_lon2_file        = [ False for shell in range(shellsmax) ]
            self.par_lon2_file_orig   = [ False for shell in range(shellsmax) ]
            self.par_lon2_shell       = [ shell for shell in range(shellsmax) ]  # mostly for checking purpuses
            self.par_lon2_disp        = [ False for shell in range(shellsmax) ]
            self.par_lon2_stringadd   = [ False for shell in range(shellsmax) ]
            self.par_lon2_param_string = [ False for shell in range(shellsmax) ]
            self.par_lon2_param_folder = [ False for shell in range(shellsmax) ]

            self.par_lon2_add_pot         = [ False for shell in range(shellsmax) ]
            self.par_lon2_add_params      = [ False for shell in range(shellsmax) ]
            self.par_lon2_add_file        = [ False for shell in range(shellsmax) ]
            self.par_lon2_add_file_orig   = [ False for shell in range(shellsmax) ]
            self.par_lon2_add_shell       = [ shell for shell in range(shellsmax) ]  # mostly for checking purpuses
            self.par_lon2_add_disp        = [ False for shell in range(shellsmax) ]
            self.par_lon2_add_param_string = [ False for shell in range(shellsmax) ]
            self.par_lon2_add_param_folder = [ False for shell in range(shellsmax) ]

            self.par_tox_add_pot         = [ False for shell in range(shellsmax) ]
            self.par_tox_add_params      = [ False for shell in range(shellsmax) ]
            self.par_tox_add_file        = [ False for shell in range(shellsmax) ]
            self.par_tox_add_file_orig   = [ False for shell in range(shellsmax) ]
            self.par_tox_add_shell       = [ shell for shell in range(shellsmax) ]  # mostly for checking purpuses
            self.par_tox_add_disp        = [ False for shell in range(shellsmax) ]
            self.par_tox_add_param_string = [ False for shell in range(shellsmax) ]
            self.par_tox_add_param_folder = [ False for shell in range(shellsmax) ]

            self.par_toy_add_pot         = [ False for shell in range(shellsmax) ]
            self.par_toy_add_params      = [ False for shell in range(shellsmax) ]
            self.par_toy_add_file        = [ False for shell in range(shellsmax) ]
            self.par_toy_add_file_orig   = [ False for shell in range(shellsmax) ]
            self.par_toy_add_shell       = [ shell for shell in range(shellsmax) ]  # mostly for checking purpuses
            self.par_toy_add_disp        = [ False for shell in range(shellsmax) ]
            self.par_toy_add_param_string = [ False for shell in range(shellsmax) ]

            self.par_toz_add_pot         = [ False for shell in range(shellsmax) ]
            self.par_toz_add_params      = [ False for shell in range(shellsmax) ]
            self.par_toz_add_file        = [ False for shell in range(shellsmax) ]
            self.par_toz_add_file_orig   = [ False for shell in range(shellsmax) ]
            self.par_toz_add_shell       = [ shell for shell in range(shellsmax) ]  # mostly for checking purpuses
            self.par_toz_add_disp        = [ False for shell in range(shellsmax) ]
            self.par_toz_add_param_string = [ False for shell in range(shellsmax) ]
            return
        print "yo"
        def get_par_lon(self, shell = False, disp = False, stringadd = False, parametrization_string = False):
            if self.verbose:
                print utils.printblue("get_par_lon ...")
            self.lon_shell = shell
            self.lon_disp = disp
            self.lon_stringadd = stringadd
            self.lon_parametrization_string = parametrization_string

            # in case there is a LON2 this will be changed in get_par_lon2r_lon2() function
            self.lon2_shell = False
            self.lon2_disp = False
            self.lon2_stringadd = False
            self.lon2_parametrization_string = False
            self.par_lon_params,self.par_lon_file,self._par_param_folder, self.par_lon_pot, self.par_lon_file_orig = \
                    self.get_parametrization(shell = shell, disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = 'lon')

        def get_par_lon_add(self, shell = False, disp = False, stringadd = False, parametrization_string = False):
            ''' disp: 'quer', 'xdir'                           # here the parametrization cames from
                stringadd: 3x3x3kp_vasp4, 10x10x10kp, ...      # here the parametrization cames from
                parametrization_string: neg_rlv_0.15_mc1       # here the parametrization cames from
                '''
            if self.verbose:
                print utils.printblue("get_par_lon_add ...")
            self.lon_add_shell = shell
            self.lon_add_disp = disp
            self.lon_add_stringadd = stringadd
            self.lon_add_parametrization_string = parametrization_string

            # in case there is a LON2 this will be changed in get_par_lon2r_lon2() function
            #`self.lon2_disp = False
            #`self.lon2_stringadd = False
            #`self.lon2_parametrization_string = False
            #print "disp:",disp
            #print "stringadd:",stringadd
            #print "parametrization_string:",parametrization_string
            self.par_lon_add_params,self.par_lon_add_file,self._par_param_folder, self.par_lon_add_pot, self.par_lon_add_file_orig = False, False, False, False, False
            if type(shell) == bool: return
            if type(disp) == bool and type(stringadd) == bool and type(parametrization_string) == bool:
                return
            self.par_lon_add_params,self.par_lon_add_file,self._par_param_folder, self.par_lon_add_pot, self.par_lon_add_file_orig = \
                self.get_parametrization(shell = shell, disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = 'lonadd')

        def get_par_lon2_add(self, shell = False, disp = False, stringadd = False, parametrization_string = False):
            ''' disp: 'quer', 'xdir'                           # here the parametrization cames from
                stringadd: 3x3x3kp_vasp4, 10x10x10kp, ...      # here the parametrization cames from
                parametrization_string: neg_rlv_0.15_mc1       # here the parametrization cames from
                '''
            if self.verbose:
                print utils.printblue("get_par_lon2_add ...")
            self.lon2_add_shell = shell
            self.lon2_add_disp = disp
            self.lon2_add_stringadd = stringadd
            self.lon2_add_parametrization_string = parametrization_string

            # in case there is a LON2 this will be changed in get_par_lon2r_lon2() function
            #`self.lon2_disp = False
            #`self.lon2_stringadd = False
            #`self.lon2_parametrization_string = False
            #print "disp:",disp
            #print "stringadd:",stringadd
            #print "parametrization_string:",parametrization_string
            self.par_lon2_add_params,self.par_lon2_add_file,self._par_param_folder, self.par_lon2_add_pot, self.par_lon2_add_file_orig = False, False, False, False, False
            if type(shell) == bool: return
            if type(disp) == bool and type(stringadd) == bool and type(parametrization_string) == bool:
                self.par_lon2_add_params,self.par_lon2_add_file,self._par_param_folder, self.par_lon2_add_pot, self.par_lon2_add_file_orig = False, False, False, False, False
                return
            else:
                if self.par_lon2_pot == 'poly':
                    sys.exit("you cant have self.ar_lon2_pot == 'poly' and add another poly (not yet)")
                self.par_lon2_add_params,self.par_lon2_add_file,self._par_param_folder, self.par_lon2_add_pot, self.par_lon2_add_file_orig = \
                    self.get_parametrization(shell = shell, disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = 'lon2add')

        def get_par_lon2(self, shell = False, disp = False, stringadd = False, parametrization_string = False):
            if self.verbose:
                print utils.printblue("get_par_lon2 ...")
            self.lon2_disp = disp
            self.lon2_stringadd = stringadd
            self.lon2_parametrization_string = parametrization_string
            if type(disp) == bool:
                self.par_lon2_params = False
                self.par_lon2_file = False
                self.par_lon2_file_orig = False
                self.par_lon2_pot = False
                self._par_param_folder = False
            else:
                self.par_lon2_params, self.par_lon2_file, self._par_param_folder, self.par_lon2_pot, self.par_lon2_file_orig = \
                        self.get_parametrization(shell = shell, disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = 'lon')

        def get_par_tox(self, shell = False, disp = False, stringadd = False, parametrization_string = False, contr = 'tox'):
            if self.verbose:
                print utils.printblue("get_par_tox ...")
            self.tox_disp = disp
            self.tox_stringadd = stringadd
            self.tox_parametrization_string = parametrization_string
            if type(disp) == bool:
                self.par_tox_params = False
                self.par_tox_file = False
                self.par_tox_file_orig = False
                self._par_param_folder = False
            else:
                self.par_tox_params,self.par_tox_file,self._par_param_folder, self.par_tox_pot, self.par_tox_file_orig = \
                        self.get_parametrization(shell = shell, disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = contr )

        def get_par_toy(self, shell = False, disp = False, stringadd = False, parametrization_string = False, contr = 'toy'):
            if self.verbose:
                print utils.printblue("get_par_toy ...")
            self.toy_disp = disp
            self.toy_stringadd = stringadd
            self.toy_parametrization_string = parametrization_string
            if type(disp) == bool:
                self.par_toy_params = False
                self.par_toy_file = False
                self.par_toy_file_orig = False
                self._par_param_folder = False
            else:
                self.par_toy_params,self.par_toy_file,self._par_param_folder, self.par_toy_pot, self.par_toy_file_orig = \
                        self.get_parametrization(shell = shell, disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = contr )

        def get_par_toz(self, shell = False, disp = False, stringadd = False, parametrization_string = False, contr = 'toz'):
            if self.verbose:
                print utils.printblue("get_par_toz ...")
            self.toz_disp = disp
            self.toz_stringadd = stringadd
            self.toz_parametrization_string = parametrization_string
            if type(disp) == bool:
                self.par_toz_params = False
                self.par_toz_file = False
                self.par_toz_file_orig = False
                self._par_param_folder = False
            else:
                self.par_toz_params,self.par_toz_file,self._par_param_folder, self.par_toz_pot, self.par_toz_file_orig = \
                        self.get_parametrization(shell = shell, disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = contr )

    def print_parametrizations(self):
        print ""
        print utils.printblue("print_parametrizations ...")
        print "lon :",self.par_lon_pot,self.par_lon_params
        print "lonadd:",self.par_lon_add_pot,self.par_lon_add_params
        print "lon2:",self.par_lon2_pot,self.par_lon2_params
        #print "lon2add:",self.par_lon2_pot,self.par_lon2_params
        if type(self.par_tox_params) != bool:
            print "tox :",self.par_tox_pot,self.par_tox_params
        if type(self.par_toy_params) != bool:
            print "toy :",self.par_toy_pot,self.par_toy_params
        if type(self.par_toz_params) != bool:
            print "toz :",self.par_toz_pot,self.par_toz_params
        print ""
        print " lono   :",self.par_lon_file_orig
        print " lon    :",self.par_lon_pot,self.par_lon_file
        print ""
        print " lonaddo:",self.par_lon_add_file_orig
        print " lonadd :",self.par_lon_add_pot,self.par_lon_add_file
        print ""
        print " lon2   :",self.par_lon2_pot,self.par_lon2_file
        if type(self.par_tox_file) != bool:
            print " tox :",self.par_tox_file
        if type(self.par_toy_file) != bool:
            print " toy :",self.par_toy_file
        if type(self.par_toz_file) != bool:
            print " toz :",self.par_toz_file
        if type(self.par_tox_file) != bool:
            print " toxo:",self.par_tox_file_orig
        if type(self.par_toy_file) != bool:
            print " toyo:",self.par_toy_file_orig
        if type(self.par_toz_file) != bool:
            print " tozo:",self.par_toz_file_orig
        print utils.printblue("     print_parametrizations complete ...")
        print ""
        return

    def create_calculate_energy_and_forces(self):
        filevorlage = self.folder_base_jobvorlage+"/calculate_energy_and_forces"
        filevorlagetmp = filevorlage+"_tmp"

        self.dataka = utils.file2dat(filevorlage)
        #self.dataka.insert(30,[ i+"\n" for i in self.calculate_energy_forces_input] )
        self.dataka.insert(30,''+'\n' )
        for idx,i in enumerate(self.calculate_energy_forces_input):
            #print "idx:",idx
            #self.dataka.insert(idx+31,self.calculate_energy_forces_input[0]+'\n' )
            self.dataka.insert(idx+30,i+'\n' )
        utils.dat2file(filevorlagetmp,self.dataka)
        return filevorlagetmp

    def get_job_number(self):
        ''' creates self.jobnumber and self.jobnumberstr '''
        jobnumber = ""
        self.jobnumber = False
        alljobs = glob.glob(self.folder_job_tdi+"[0-9]*")
        alljobsnrs = [ i.split(self.folder_job_tdi)[1].split("_")[0] for i in alljobs ]
        alljobsnrs = list(set(alljobsnrs))
        if self.verbose:
            print "     alljobsnrs: ...",alljobsnrs[-5:]
        if len(alljobsnrs) == 0:
            self.jobnumber = 0
            self.jobnumberstr = '0___'
            return self.jobnumber
        for i in np.arange(99):
            if str(i) in alljobsnrs:
                continue
            else:
                self.jobnumber = i
                if self.jobnumber < 10.:
                    self.jobnumberstr = str(i)+"___"
                if 99. > self.jobnumber >= 10.:
                    self.jobnumberstr = str(i)+"__"
                if self.jobnumber >= 100.:
                    self.jobnumberstr = str(i)+"_"
                break
        return self.jobnumber

    #def get_job_name(self):
    #    lonpart = "from_"+self.forparametrization_sc_string+"sc__LON_"+self.lon_disp+"_"+self.lon_stringadd+"_"+self.lon_parametrization_string
    #    if self.par_lon2_params == False:
    #        lon2part = "_LON2_NONE"
    #    else:
    #        lon2part = "_LON2_"+self.lon2_disp+"_"+self.lon2_stringadd+"_"+self.lon2_parametrization_string

    #    if type(self.tox_disp) == bool or type(self.tox_stringadd) == bool or type(self.tox_parametrization_string) == bool:
    #        toxpart = "_TOX_NONE"
    #    else:
    #        toxpart = "_TOX_"+self.tox_disp+"_"+self.tox_stringadd+"_"+self.tox_parametrization_string

    #    if type(self.toy_disp) == bool or type(self.toy_stringadd) == bool or type(self.toy_parametrization_string) == bool:
    #        toypart = "_TOY_NONE"
    #    else:
    #        toypart = "_TOY_"+self.toy_parametrization_string

    #    if type(self.toz_disp) == bool or type(self.toz_stringadd) == bool or type(self.toz_parametrization_string) == bool:
    #        tozpart = "_TOZ_NONE"
    #    else:
    #        tozpart = "_TOZ_"+self.toz_parametrization_string

    #    if type(self.par_lon_add_params) == bool:
    #        lonaddpart = "_LONADD__NO"
    #    else:
    #        lonaddpart = "_LONADD_YES"

    #    if type(self.par_lon2_add_params) == bool:
    #        lon2addpart = "_LON2ADD__NO"
    #    else:
    #        lon2addpart = "_LON2ADD_YES"

    #    doscut = "PTS_allpts_"
    #    if self.DOSparametrization == True:
    #        doscut = "PTS_dosall_"
    #    if self.DOScutparametrization == True:
    #        doscut = "PTS_doscut_"

    #    self.jobnamestr = doscut + lonpart + lonaddpart + lon2part + lon2addpart + toxpart + toypart + tozpart
    #    self.jobnamestr = self.jobnamestr.replace("order","")
    #    self.jobnamestr = self.jobnamestr.replace("_fit","")
    #    return self.jobnamestr

    def get_job_info(self,ls = False, execute = False, jobnumber = False, verbose = False):
        ''' if at some point we will calculate ohter alats, jobnumber will have to also take care of this
            currently no mater which jobnumber, the alat ist the same '''
        print utils.printblue("get_job_info ...")

        # check that u is mounted
        if os.path.isdir(self.folder_job_tdi_base) != True:
            sys.exit(self.folder_job_tdi_base+" not found self.folder_job_tdi_base --> garmount")


        if self.verbose == True:
            print "     self.folder_job_tdi_base:",self.folder_job_tdi_base

        # get path for job
        if self.verbose == True:
            print "     self.folder_job_tdi:",self.folder_job_tdi


        ###########################################################################################
        ###########################################################################################
        # in case we got a jubnumber, load old job info
        ###########################################################################################
        ###########################################################################################
        self.jobnew = True
        if type(jobnumber) != bool:
            self.jobnew = False
            self.jobnumber = int(jobnumber)
            if self.verbose == True:
                print "     self.jobnumber:",self.jobnumber

            if type(self.jobnumber) != int:
                print 'self.jobnumber:',self.jobnumber,"type:",type(self.jobnumber)
                sys.exit("self.jobnumber has to be an integer")

        if self.verbose == True:
            print utils.printyellow("     self.jobnew:"+str(self.jobnew))

        ##########################################################################################################
        # EXISTING JOB (OLD JOB)
        ##########################################################################################################
        if self.jobnew != True:
            # get old path
            search = self.folder_job_tdi+"/"+str(self.jobnumber)+"_*"
            self.jobpath = glob.glob(search)
            if len(self.jobpath) != 1:
                print "search:",search
                print "self.jobpath:",self.jobpath
                sys.exit("not found self.jobpath:")
            self.jobpath = self.jobpath[0]
            if self.verbose == True:
                print utils.printgreen("     self.jobpath:"+self.jobpath)

            self.jobnamestr = ""

            self.par_lon_params,self.par_lon_file,self._par_param_folder, self.par_lon_pot, self.par_lon_file_orig =     "", "", "", "", ""
            self.par_lon2_params,self.par_lon2_file,self._par_param_folder, self.par_lon2_pot, self.par_lon2_file_orig = "", "", "", "", ""
            self.par_tox_params,self.par_tox_file,self._par_param_folder, self.par_tox_pot, self.par_tox_file_orig =     "", "", "", "", ""
            self.par_toy_params,self.par_toy_file,self._par_param_folder, self.par_toy_pot, self.par_toy_file_orig =     "", "", "", "", ""
            self.par_toz_params,self.par_toz_file,self._par_param_folder, self.par_toz_pot, self.par_toz_file_orig =     "", "", "", "", ""

            # get old parametrization
            self.par_lon_file, self.par_lon2_file, self.par_tox_file, self.par_toy_file, self.par_toz_file, \
            self.par_lon_file_orig, self.par_tox_file_orig, self.par_toy_file_orig, self.par_toz_file_orig = \
             get_parametrizationfiles_oldpath(self.jobpath,self.folder_displacement_direction_all)
            #    lonfdot, toxfdot, toyfdot, tozfdot = get_parametrizationfiles_oldpath(self.jobpath,self.folder_displacement_direction_all)


            self.par_lon_params, self.par_lon_pot = self.get_parametrization_parameters_from_file(self.par_lon_file)
            self.par_lon2_params, self.par_lon2_pot = self.get_parametrization_parameters_from_file(self.par_lon2_file)
            self.par_tox_params, self.par_tox_pot = self.get_parametrization_parameters_from_file(self.par_tox_file)
            self.par_toy_params, self.par_toy_pot = self.get_parametrization_parameters_from_file(self.par_toy_file)
            self.par_toz_params, self.par_toz_pot = self.get_parametrization_parameters_from_file(self.par_toz_file)

            ##############################################################
            # do a check that old parametrization is current one
            ##############################################################
            calculate_energy_and_forces_path =  self.jobpath+"/calculate_energy_and_forces"
            if os.path.isfile(calculate_energy_and_forces_path) != True:
                sys.exit(calculate_energy_and_forces_path+" calculate_energy_and_forces_path not found!")
            self.check_calculate_energy_and_forces_parameters_with_current(calculate_energy_and_forces_path)


        ##########################################################################################################
        # NEW JOB
        ##########################################################################################################
        if self.jobnew == True:

            # get job number
            self.jobnumber = self.get_job_number()
            if self.verbose == True:
                print "     self.jobnumber:",self.jobnumber


            # get job name
            self.jobname = (self.jobnumberstr+self.jobname_param)[:-3]  # last 3 chars are ___
            if self.verbose == True:
                print "     self.jobname:",self.jobname

            self.jobpath = self.folder_job_tdi+self.jobname
            if self.verbose == True:
                print "     self.jobpath:",self.jobpath

            if os.path.isdir(self.jobpath)  == True:
                sys.exit("self.jobpath does already exist "+self.jobpath)


            # create file calculate_energy_and_forces_tmp
            self.calculate_energy_and_forces_pathtmp = self.create_calculate_energy_and_forces()
            if self.verbose == True:
                print "     self.calculate_energy_and_forces_pathtmp:",self.calculate_energy_and_forces_pathtmp

        print utils.printblue("     get_job_info complete ...")
        print ""
        return

    def create_job_dispdudl(self):
        ''' is a job from the POSITIONs of the displacement folder '''
        if self.verbose:
            print ""
            print utils.printblue("create_job_dispdudl ...")


        self.create_self_jobpath()

        for disp in [ "quer", "xdir" ]:
            #parametrization_folder =  self.get_parametrization_folder(disp = disp, stringadd = self.kpstring)
            print self.get_element_2x2x2sc_info(self.element, self.sc_create)[1]
            parametrization_folder = self.folder_displacement_element+'/'+self.sc_create_string+"sc_"+disp+"_"+self.get_element_2x2x2sc_info(self.element, self.sc_create)[1]
            if os.path.isdir(parametrization_folder) != True:
                sys.exit("parametrization_folder not found: "+parametrization_folder)
            print "     parametrization_folder:",parametrization_folder
            POSITIONs = parametrization_folder+"/POSITIONs_short"
            dudl = parametrization_folder+"/dudl_short"
            print "     POSITIONs:",POSITIONs
            if os.path.isfile(POSITIONs) != True:
                sys.exit("evaluate disp, no POSITIONs file found")
            print "     dudl:",dudl
            if os.path.isfile(dudl) != True:
                sys.exit("evaluate disp, no dudl file found")

            # make a tmp folder to calculate this
            if os.path.isdir(self.jobpath+"/disp_dudl") != True:
                os.makedirs(self.jobpath+"/disp_dudl")

            folder_dispdudl = self.jobpath+"/disp_dudl/"+disp
            if os.path.isdir(folder_dispdudl) != True:
                os.makedirs(folder_dispdudl)
            print "     folder_dispdudl:",folder_dispdudl


            sys.exit()
            # necessary:
            # ----------
            # calculate_energy_and_forces OK
            # dUdL OK
            # EqCoords_direct  OK
            # POSITIONs OK
            # cell
            # --> DOS_POSITIONS_auswerten.py -dudlposc
            shutil.copy2(self.calculate_energy_and_forces_pathtmp,folder_dispdudl+"/calculate_energy_and_forces")
            shutil.copy2(self.eqcoordsfile_job,folder_dispdudl+"/EqCoords_direct")
            shutil.copy2(POSITIONs,folder_dispdudl+"/POSITIONs")
            shutil.copy2(dudl,folder_dispdudl+"/dUdL")
            shutil.copy2(self.jobvorlage_createfolder_info.jobvorlage_DOSlon_file,folder_dispdudl+"/DOSlon")
            shutil.copy2(self.jobvorlage_createfolder_info.jobvorlage_DOSlon_filecut,folder_dispdudl+"/DOSloncut")
            np.savetxt(folder_dispdudl+"/cell",self.cell)
            hier=os.getcwd()
            print "------------1"
            os.chdir(folder_dispdudl)
            print "------------2"
            print "------------3"


            self.dudlnew, self.dudl_vs_disp = start_dispdudl(nndist = self.nndist)
            os.chdir(hier)

        if self.verbose:
            print utils.printblue("     create_job_dispdudl complete ...")
            print ""
        return

    def create_self_jobpath(self):
        # create jobpath  (.../ti/Cu/1___from_2x2x2sc__LON_xdir...)
        if os.path.isdir(self.jobpath) != True:
            os.makedirs(self.jobpath)
            print "     created  self.jobpath:",self.jobpath

        #print ""
        #print "self.jobpath :",self.jobpath
        if os.path.isdir(self.jobpath) != True:
            sys.exit("creation of self.jobpath faied! "+str(self.jobpath))
        return

    def create_job_tdi(self,ls = False, verbose = False):
        print ""
        print utils.printblue("create_job_tdi ...")
        # check if all files are available for job creation

        if os.path.isdir(self.jobpath) == True:
            print utils.printred('      job exists: '+self.jobpath)
            return 'job exists: '+self.jobpath


        #########################################################
        # create jobpath  (creates self.jobpath)
        # e.g. /u/aglen/Understand_distributions/ti//Pt/3___from_2x2x2sc__LON_xdir_3x3x3kp_neg_\
        #           rlv_0.25_mc1_LON2_xdir_3x3x3kp_pos_poly9_TOX_xdir_3x3x3kp_all_poly_1st_TOY_all_poly_1st_TOZ_all_poly_1st
        #########################################################
        self.create_self_jobpath()


        POSCARVORLAGE = self.folder_base_jobvorlage+"/POSCAR_"+self.sc_create_string+"sc"
        POSCARVORLAGETMP = self.folder_base_jobvorlage+"/POSCAR_"+self.sc_create_string+"sc_tmp"
        print "+self.sc_create_string:",self.sc_create_string
        print "POSCARVORLAGETMP:",POSCARVORLAGETMP
        sys.exit()
        utils.isfile(POSCARVORLAGE)
        shutil.copy2(POSCARVORLAGE,POSCARVORLAGETMP)

        #print "self.scalat:",self.scalat,str(self.scalat)
        utils.sed(filename = POSCARVORLAGETMP,stringsearch='xxxSCxxx',stringreplace=str(self.scalat))

        #print "self.eq:",self.eqcoordsfile_job

        jobList = self.jobpath+"/jobList"

        shutil.copy2(self.calculate_energy_and_forces_pathtmp,self.jobpath+"/calculate_energy_and_forces")

        #########################################################
        # create_job_dispdudl()
        #########################################################
        self.create_job_dispdudl()
        sys.exit()

        #########################################################
        # create_job_dispparametrization()
        #########################################################
        self.jobpath_parametrization = self.jobpath+"/parametrization"
        create_job_parametrization(jobpath_parametrization = self.jobpath_parametrization, \
                lon = self.par_lon_file, \
                lon2 = self.par_lon2_file, \
                tox = self.par_tox_file, \
                toy = self.par_toy_file, \
                toz = self.par_toz_file, \
                lonallp = self.par_lon_file_orig,\
                toxallp = self.par_tox_file_orig,\
                toyallp = self.par_toy_file_orig,\
                tozallp = self.par_toz_file_orig\
                )

        for l in ls:
            jobpathl = self.jobpath+"/lambda"+str(l)
            print utils.printgreen("    jobpathl: "+jobpathl)
            os.makedirs(jobpathl)
            shutil.copy2(self.calculate_energy_and_forces_pathtmp,jobpathl+"/calculate_energy_and_forces")
            shutil.copy2("/Users/glensk/Thermodynamics/python_thermodynamics/hesse.py",jobpathl+"/hesse.py")
            #shutil.copy2(self.folder_base_jobvorlage+"/hesse.py",jobpathl+"/hesse.py")
            shutil.copy2(self.eqcoordsfile_job,jobpathl+"/EqCoords_direct")
            shutil.copy2(self.jobvorlage_createfolder+"/HesseMatrix_sphinx",jobpathl+"/HesseMatrix_sphinx")
            shutil.copy2(self.jobvorlage_createfolder+"/INCAR",jobpathl+"/INCAR")
            #print "a:"
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR LAMBDA "+str(l)+"\'", dont_raise_exceptino = False)
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR REF_TYPE extern"+"\'", dont_raise_exceptino = False)
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR PRE_EQ_N -400"+"\'", dont_raise_exceptino = False)
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR NSW "+str(self.steps)+"\'", dont_raise_exceptino = False)
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR NPAR 4"+"\'", dont_raise_exceptino = False)
            #print "b:"
            #shutil.copy2(self.jobvorlage_createfolder+"/POSCAR",jobpathl+"/POSCAR")
            shutil.copy2(POSCARVORLAGETMP,jobpathl+"/POSCAR")
            shutil.copy2(self.jobvorlage_createfolder+"/KPOINTS",jobpathl+"/KPOINTS")
            shutil.copy2(self.POTCAR,jobpathl+"/POTCAR")

            # append to jobList
            if os.path.isfile(jobList) != True:
                open(jobList, 'a').close()
            with open(jobList, "a") as myfile:
                myfile.write(jobpathl+"\n")

        return

    def get_element_2x2x2sc_info(self,element, sc):
        # max volume; kpoint in ti;
        if sc == 2:            # parameters for parametrization are taken from here (2x2x2sc disp quer xdir)
            if element == 'Al': return [ 4.13 , "3x3x3kp", "250" , '934' ]
            if element == 'Pb': return [ 5.13 , "4x4x4kp", "300" , '601' ]
            if element == 'Cu': return [ 3.75 , "3x3x3kp", "260" , '1360' ]
            if element == 'Rh': return [ 3.98 , "2x2x2kp", "270" , '2237' ]
            if element == 'Rh': return [ 3.98 , "6x6x6kp", "270" , '2237' ]
            if element == 'Pd': return [ 4.1  , "3x3x3kp", "300" , '1830' ]
            if element == 'Ag': return [ 4.31 , "3x3x3kp", "225" , '1235' ]
            if element == 'Ir': return [ 3.99 , "3x3x3kp", "300" , '2739' ]
            #if element == 'Ir': return [ 3.99 , "6x6x6kp", "300" , '2739' ]
            if element == 'Pt': return [ 4.1  , "3x3x3kp", "200" , '2042' ]
            if element == 'Au': return [ 4.25 , "3x3x3kp", "300" , '1338' ]
            if element == 'Ni': return [ 3.62 , "3x3x3kp", "300" , '1728' ]
        if sc == 3:         # parameters for parametrization are taken from here (2x2x2sc disp quer xdir)
            if element == 'Al': return [ 4.13 , "2x2x2kp", "250" , '934' ]
            if element == 'Pb': return [ 5.13 , "4x4x4kp", "300" , '601' ]
            if element == 'Cu': return [ 3.75 , "2x2x2kp", "260" , '1360' ]
            if element == 'Rh': return [ 3.98 , "2x2x2kp", "270" , '2237' ]
            if element == 'Pd': return [ 4.1  , "2x2x2kp", "300" , '1830' ]
            if element == 'Ag': return [ 4.31 , "2x2x2kp", "225" , '1235' ]
            if element == 'Ir': return [ 3.99 , "2x2x2kp", "300" , '2739' ]
            if element == 'Pt': return [ 4.1  , "2x2x2kp", "200" , '2042' ]
            if element == 'Au': return [ 4.25 , "2x2x2kp", "300" , '1338' ]
            if element == 'Ni': return [ 3.62 , "2x2x2kp", "300" , '1728' ]
        print "sc:",sc
        print "element:",element
        sys.exit("not found sc or element; sc (where ENCUT and TMELT are predefined):"+str(sc)+":   element:"+str(element)+":")

    def printelement(self):
        print ""
        print utils.printred("#########################################################")
        print utils.printred("# "+ pot.element)
        print utils.printred("#########################################################")
        return

    # unimportant or just once used scripts
    def repair_disp_dudl_folder():
        ''' repair disp_dudl folder '''
        #self.folder_job_tdi: /u/aglen/Understand_distributions/ti//Au/
        for disp in [ "quer", "xdir" ]:
            #to_rep_str = pot.folder_job_tdi+"/[2-5]*/disp_dudl/"+disp
            to_rep_str = pot.folder_job_tdi+"/[0]*/disp_dudl/"+disp
            to_rep_folder = glob.glob(to_rep_str)
            parametrization_folder =  pot.get_parametrization_folder(disp = disp, stringadd = pot.kpstring)
            positions = parametrization_folder+"/POSITIONs"
            dudl = parametrization_folder+"/dUdL"
            if os.path.isfile(positions) != True:
                sys.exit("pos do not exist")
            if os.path.isfile(dudl) != True:
                sys.exit("dudl do not exist")
            print ""
            print disp,positions
            for folder_dispdudl in to_rep_folder:
                print "folder_dispdudl:",folder_dispdudl,pot.nndist
                shutil.copy2(positions,folder_dispdudl+"/POSITIONs")
                shutil.copy2(dudl,folder_dispdudl+"/dUdL")
                np.savetxt(folder_dispdudl+"/cell",pot.cell)
                shutil.copy2(pot.eqcoordsfile_job,folder_dispdudl+"/EqCoords_direct")

                hier=os.getcwd()
                print "------------1"
                os.chdir(folder_dispdudl)
                print "------------2"
                print "------------3"
                DOS_POSITIONS_auswerten.dudlposc()
                print "------------4"
                #utils.run2("DOS_POSITIONS_auswerten.py -dudlposc")
                os.chdir(hier)

                dudlnewfile = folder_dispdudl+"/dUdLnew"
                if os.path.isfile(dudlnewfile) != True:
                    sys.exit("dudlnew was not created! in "+str(folder_dispdudl))
                pot.dudlnew = np.loadtxt(dudlnewfile)
                print pot.dudlnew
                dudl_vs_disp = np.array([pot.dudlnew[:,0]+pot.nndist,pot.dudlnew[:,6]]).transpose()
                np.savetxt(folder_dispdudl+"/dUdLvs_disp",dudl_vs_disp,fmt="%.3f %.4f")
        return

    def get_elements(self, elements = False):
        self.elements_all = [ "Al", "Pb", "Cu", "Rh", "Pd", "Ag", "Ir", "Pt", "Au", "Ni" ]

        if elements != False:  # z.B. al
            self.elements = [i.title() for i in elements]
        else:
            self.elements = self.elements_all

        print "elements:",self.elements
        return


if __name__ == '__main__':
    # help
    help_string = '''
    This is the main skript which controls: pot_parametrize.py, pot_...
    Start with, e.g., run ~/Thermodynamics/python_thermodynamics/pot_info_startjob.py -v -e ir

    run ~/Thermodynamics/python_thermodynamics/pot_info_startjob.py -v -e ir -dispdirection quer -a 3.99
    run ~/Thermodynamics/python_thermodynamics/pot_info_startjob.py -v -e al -sp 3 -dispdirection quer -a 4.07 -kp 10x10x10kp -stradd "_vasp4_ENCUT400"
    run  /Users/glensk/Thermodynamics/python_thermodynamics/pot_info_startjob.py  -e al -a 4.14 -sp 5 -kp 2x2x2kp -dispdirection quer -pnnforcesnpz -v -stradd "_vasp4_ENCUT250_GGA_2_atoms_displaced"


    und danach:
    pot.disp.alle.data[1,'lon','quer','dos'].morse[0]
    pot.disp.alle.data[1,'lon','quer','dos'].morse[1]
    pot.disp.alle.data[1,'lon','quer','dos'].morse[2]


    The prametrization is performed by creating a sequence of files:
        1) nnforces.npz         --> Forces, positions from displacements from VASP
        2) xxx___data___.pkl    --> parametrization long
        3) FML.npz              --> xdirforce - longforce
                                --> querforce - longforce

                                '''
    p = argparse.ArgumentParser(description=help_string,formatter_class=RawTextHelpFormatter)
    p.add_argument('-e',   nargs="+", choices=['al', 'pb', 'cu', 'rh', 'pd', 'ag', 'ir', 'pt', 'au', 'ni' ], help='specify element', default=False)
    p.add_argument('-a', default=False, type=float, help='define the alat for parametrization')
    p.add_argument('-sp', default=2, type=int,
            help='supercell for parametrization; default = 2; needs sc (int) as argument')
    p.add_argument('-kp', default=False, type=str, help='String of kpoint parametrization, e.g. 3x3x3kp; default=False')
    p.add_argument('-stradd', default="", type=str, help='String to add to make parametrizatoinfolder unique; default=""')
    p.add_argument('-dispdirection',   nargs="+", choices=['quer', 'xdir', 'midd', '3nnd', '4nnd' ],
        help='direction of displacements when parametrizing; default == False == all == {quer+xdir+midd+3nnd+4nnd}',
        default=False)




    p.add_argument('-dnl', default=False, action = 'store_true',
        help='do not load all {nnforces,FML}.npz data.pkl files (default = True))')

    p.add_argument('-pnnforcesnpz', default=False, action = 'store_true',
        help='1. parametrization of new nnforces.npz files (containing the forces at T=0K); this also creates new global xxx_data_xxx.pkl file; default = False')
    p.add_argument('-ppkldata', default=False, action = 'store_true',
        help='(2.) parametrizations of xxx_data_xxx.pkl files (containing the parametrized long functions at T=0K); ALSO RUN THIS IF NO FML.npz files')
    p.add_argument('-pnpzfml', default=False, action = 'store_true',
        help='parametrizations of new FML.npz file and ? analyze.pkl files (containing the parametrization fml,tox,tix functions at T=0K); default = False\n\n\n')



    p.add_argument('-analyze', default=False, action = 'store_true',
        help='print analyze/{Element} files')
    p.add_argument('-pnpzqubus', default=False, action = 'store_true',
        help='create qubus.{fml,F}.npz')

    p.add_argument('-analyze_qubus', default=False, action = 'store_true',
        help='print {Element}/analyze_qubus files')



    p.add_argument('-j', default=False, action='store_true',
        help='create job')
    p.add_argument('-je', default=False, type=int, choices=[1,2,3],
        help='executable Nr.\n'
        '1: executable 1\n'
        '2: executable 2\n'
        )
    p.add_argument('-jd', default=False, action='store_true',
        help='create job from die POSITIONs of the displacements_dense file (== dudl from T=0K (quer/xdir) displacements)\n\n\n')






    #p.add_argument('-l',
    #            help='load all data', action='store_true', default=False)
    p.add_argument('-o', '--o',
                help='get parametrization of old job; just give the jobnumber (int)', action='store_true', default=False)
    p.add_argument('-r', '--r',
                help='r', action='store_true', default=False)
    p.add_argument('-effectlontox', '--effectlontox',
                type=int, nargs=2, help='difference between two parametrizations (need tow jobnumbers as input)', default=False)
    p.add_argument('-dudl_vs_disp', action='store_true', default=False,
                help='start dUdLvs_disp here (will be created in current xdir/quer folder)')
    p.add_argument('-v', '--verbose',action='count',
                help='verbosity level: v:verbose, vv:even more verbose', default=False)
    args = p.parse_args()


    hier=os.getcwd()
    pot = pot()

    pot.forparametrization_kpoint_string = args.kp
    pot.forparametrization_stringadd = args.stradd
    pot.alat = args.a
    pot.verbose = args.verbose
    pot.get_elements(args.e)

    # initial values, can be changed
    pot.sc_create = 2                               # create a 2x2x2sc or a 3x3x3sc ? (is alredy used in pot.init_variables)
    pot.folder_displacement_from = "displacements"
    pot.folder_displacement_from = "displacements_dense"
    pot.folder_displacement_from = "displacements_"
    #pot.folder_displacement_from = "displacements_dense_6x6x6kp"

    pot.forparametrization_sc = args.sp              # parameters for parametrization are taken from here (2x2x2sc disp quer xdir)
    pot.forparametrization_displacement = args.dispdirection  # (xdir, quer, 3nnd, 4nnd, False = * = all)
    pot.parametrize_nnforcesnpz = args.pnnforcesnpz
    pot.parametrize_data_pkl = args.ppkldata
    pot.parametrize_fml_npz = args.pnpzfml
    pot.analyze = args.analyze
    pot.analyze_qubus = args.analyze_qubus
    pot.create_qubus_npz = args.pnpzqubus
    pot.dont_load_pkl_npz = args.dnl

    if pot.verbose == True:
        print "     pot.folder_base                         :",pot.folder_base
        print "     pot.folder_displacement_from            :",pot.folder_displacement_from
        print "     pot.forparametrization_sc               :",pot.forparametrization_sc
        print "     pot.forparametrization_displacement     :",pot.forparametrization_displacement
        print "     pot.parametrize_nnforcesnpz             :",pot.parametrize_nnforcesnpz
        print "     pot.parametrize_data_pkl                :",pot.parametrize_data_pkl
        print "     pot.parametrize_fml_npz                 :",pot.parametrize_fml_npz
        print "     pot.analyze                             :",pot.analyze
        print "     pot.analyze_qubus                       :",pot.analyze_qubus
        print "     pot.create_qubus_npz                    :",pot.create_qubus_npz
        print "     pot.dont_load_pkl_npz                   :",pot.dont_load_pkl_npz
        print "     pot.elements                            :",pot.elements
        print "     pot.element                             :",pot.element, "will be defined in first loop: for el in pot.elements:"
        print "     pot.forparametrization_kpoint_string    :",pot.forparametrization_kpoint_string

    ######################################################################################
    # load all data / fits in 2x2x2sc and 3x3x3sc if exists
    # perform this always
    # currently this does not make sense since we overwrite all variables
    ######################################################################################
    print utils.printred("-------------------------------- 1 -----------------------------")
    #for el in pot.elements:
    #    pot.element    = el             # Rh, Al, ...
    #    pot.init_variables()
    #    pot.load_dataframe(2)
    #    pot.load_dataframe(3)
    #    # a) load all dataframes (data/fit) in a list  (2x2x2sc and 3x3x3sc)
    #    #       - make function to loadpklfile(sc = 2)
    #    # b) concatenate all dataframes

    print utils.printred("-------------------------------- 2 -----------------------------")
    for el in pot.elements:
        pot.element    = el             # Rh, Al, ...
        if pot.verbose == True:
            print "     pot.folder_base                 :",pot.folder_base
            print "     pot.folder_displacement_from    :",pot.folder_displacement_from
        pot.folder_displacement_element = pot.folder_base+'/'+pot.folder_displacement_from+'/'+pot.element+"/"
        pot.folder_job_tdi = pot.folder_job_tdi_base +"/"+str(pot.element)+"/"
        if pot.verbose == True:
            print "     pot.folder_displacement_element :",pot.folder_displacement_element
            print "     pot.folder_job_tdi              :",pot.folder_job_tdi


        #pot.kpstring = pot.get_element_2x2x2sc_info(pot.element)[1]   # gets 3x3x3kp string

        #pot.printelement()
        ##pot.kpstring = "3x3x3kp_ISYM0"
        #pot.init_variables( kpstringoptional = pot.kpstring)
        #print utils.printred("-------------------------------- 3 -----------------------------")
        # This should only be loaded if jobs are being created
        #pot.init_variables()  # if this is only related for creating jobs, is it necessary for auswertung?

        ##################################################################################
        ###
        ##################################################################################
        #pot.DOSparametrization = True   # always leave true otherwise all points are taken
        #pot.DOScutparametrization = False   # if this is turned on the parametrization is even smaller

        print utils.printred("-------------------------------- 4 -----------------------------")
        if args.dudl_vs_disp:
            pot.dudlnew, pot.dudl_vs_disp = start_dispdudl(nndist = pot.nndist)
            sys.exit()

        print utils.printred("-------------------------------- 5 -----------------------------")
        # this simply loads the parametrization if args.pc == False
        pot.pot_load_parametrize_function()
        sys.exit()

        #if args.l:
        #    pot.load_all_data()
        #pot.param.loadDatapd(pot.param.pkl)
        #pot.param.loadFitpd(pot.param.fit_file)
        print utils.printred("-------------------------------- 6 -----------------------------")
        ############################################################################
        print utils.printblue("pars_from_pd ...")
        pot.pars_from_pd(1, 'lon'   , 'from_2x2x2sc', 'get_kpstring', True, 'quer', 'dos', 'morse' )
        #pot.pars_from_pd(1, 'lonadd', 'from_2x2x2sc', 'get_kpstring', True, 'quer', 'dos', 'morse', 'polybeste' )
        #pot.pars_from_pd(1, 'tox'   , 'from_2x2x2sc', 'get_kpstring', False, 'xdir', 'all', 'poly' )

        #pot.pars_from_pd(2, 'lon',    'from_3x3x3sc', 'get_kpstring', True, 'xdir', 'dos', 'morse' )
        #pot.pars_from_pd(2, 'lonadd', 'from_3x3x3sc', 'get_kpstring', True, 'xdir', 'dos', 'morse', 'polybeste' )
        #pot.pars_from_pd(2, 'tox'   , 'from_3x3x3sc', 'get_kpstring', False, 'xdir', 'all', 'poly' )


        for i in pot.calculate_energy_forces_input:
            print i
        print utils.printblue("     pars_from_pd complete...")
        print ""
        ############################################################################

        if args.effectlontox:
            effectlontox(foldernrs, verbose = pot.verbose)



        if args.r:
            repair_disp_dudl_folder()


        if args.j:
            pot.get_job_info(jobnumber = args.o)  # has to be before pot.pot_parametrize_function() but after pot.get_par_{lon,lon2,tox,toy,toz}
            pot.steps = 1000
            pot.create_job_tdi(ls=[0.0,1.0], \
                    execute=False,\
                    verbose = args.verbose \
                    )
            pot.execute_job_tdi(args.je)

        if args.jd:
            ''' is generally greated with args.j, this is only to redo the dispdudl folders '''
            pot.create_job_dispdudl( \
                    verbose = args.verbose \
                    )

        if args.je:
            #execute="/u/aglen/start.4.6.28.tid.quick.sh",\
            execute = args.je
            if os.path.isfile(execute) != True:
                sys.exit("execute not found: "+str(execute))
            utils.run2(command = "tcsh -c \'"+"ssh aglen@cmmc002.bc.rzg.mpg.de \"cd "+self.jobpath+";"+execute+";echo jo\""+"\'", dont_raise_exceptino = False)



        # make moduel or simple list to get the one parametrization
    os.chdir(hier)


        # Cu auslenkung 0.0: -117.23976675
        # Cu: tdi energy unausgelenkt:
        #         free energy    TOTEN  =      -117.23979717 eV
        #         energy without entropy =     -116.81457379  energy(sigma->0) =     -117.02718548
        # Cu: tdi zweiter schirtt (erster schritt im dUdL file)
        #     free  energy   TOTEN  =      -110.300832 eV
        #     energy  without entropy=     -109.900028  energy(sigma->0) =     -110.100430
        # Cu: tdi dUdL file:
        # $head dUdL
        # #  step   time(fs)  temp(K) average       U(meV/at)    Uref          dUdL   average    offset
        #       1      10.0   1764.9      0.0       223.84    181.37        42.46      0.00      0.00
        #       2      20.0   1658.3      0.0       208.65    180.63        28.02      0.00      0.00
        #
        # [129](-110.300832--117.23979717)*1000/31.
        # Out[129]: 223.83758612903236
        #
        # [130](-110.300832--117.23979717)*1000/32.
        # Out[130]: 216.8426615625001



        # TODO:
        #   2) change startscript so that it
        #       - creates calculate_energy_forces_new
        #       - gets and incudes the parametrization for the 2NN / 3NN / 4NN
