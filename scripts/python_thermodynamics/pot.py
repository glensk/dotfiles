#!/usr/bin/env python

import os
import sys
import copy
import numpy as np
import glob
import shutil
import argparse
import utils
reload(utils)
import pot_parametrize
reload(pot_parametrize)
import hesse
reload(hesse)
import DOS_POSITIONS_auswerten
reload(DOS_POSITIONS_auswerten)


p = argparse.ArgumentParser(description='''help string''')
p.add_argument('-e',   '--element', nargs="+", choices=['al', 'pb', 'cu', 'rh', 'pd', 'ag', 'ir', 'pt', 'au' ], help='specify element', default=False)
p.add_argument('-j', default=False, action='store_true',
    help='create job')
p.add_argument('-jd', default=False, action='store_true',
    help='create job for dudl from T=0K (quer/xdir) displacements')
p.add_argument('-e', default=False, action='store_true',
    help='evaluate parametrizations')
p.add_argument('-v', '--verbose',
            help='verbose', action='store_true', default=False)
p.add_argument('-k', '--k',
            help='k', action='store_true', default=False)
p.add_argument('-r', '--r',
            help='r', action='store_true', default=False)
p.add_argument('-o', '--o',
            type=str, help='old job nr', default=False)
args = p.parse_args()


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
        self.element        = False
        self.sc_forparams   = False
        self.alat           = False  # for evaluate_pot
        self.folder_base     = "/Users/glensk/Dropbox/Understand_distributions/"
        self.folder_base_jobvorlage = self.folder_base+"/jobvorlage_all/"
        self.utilsfolder    = "/Users/glensk/Thermodynamics/utilities/"
        self.folder_job_tdi_base      = "/u/aglen/Dropbox/Understand_distributions/ti/"

        self.dispstrigpos   = [ "lon", "tox" ]
        self.structure      = "fcc"                     # necessary for EqCoords
        self.steps          = 2000
        #self.piarametrizatioin_folder = False
        #############################################################################################
        # Stop editing here
        #############################################################################################


    def init_variables( self, kpstringoptional = False ):
        if self.verbose:
            #print ""
            print utils.printblue("init_variables ...")
        if type(self.sc_forparams) == bool:
            sys.exit("init_variables: please define self.sc_forparams")
        self.sc_string_params = str(self.sc_forparams)+"x"+str(self.sc_forparams)+"x"+str(self.sc_forparams)
        self.eqcoordsfile_getpot = self.utilsfolder+"/"+self.structure+"/EqCoords_direct_"+self.structure+"_"+self.sc_string_params+"sc"
        if os.path.isfile(self.eqcoordsfile_getpot) != True:
            sys.exit(self.eqcoordsfile_getpot+" self.eqcoordsfile_getpot does not exist 2 !")

        self.folder_displacement = self.folder_base+"/displacements/"+self.element
        if os.path.isdir(self.folder_displacement) != True:
            sys.exit("self.folder_displacement: "+self.folder_displacement +" does not exist 1 !")

        if type(kpstringoptional) == bool:
            self.folder_displacement_direction = \
                glob.glob(self.folder_displacement+"/"+self.sc_string_params+"sc_quer_*") + \
                glob.glob(self.folder_displacement+"/"+self.sc_string_params+"sc_xdir_*")
        else:
            self.folder_displacement_direction = \
                glob.glob(self.folder_displacement+"/"+self.sc_string_params+"sc_quer_"+kpstringoptional+"*") + \
                glob.glob(self.folder_displacement+"/"+self.sc_string_params+"sc_xdir_"+kpstringoptional+"*")
        #for i in self.folder_displacement_direction:
        #    print i
        if len(self.folder_displacement_direction) == 0:
            print "self.folder_displacement:",self.folder_displacement
            print "self.folder_displacement_direction:",self.folder_displacement_direction
            sys.exit("no dispfolder found in self.folder_displacement: "+self.folder_displacement)
        self.folder_job_tdi = self.folder_job_tdi_base +"/"+str(self.element)+"/"
        self.folder_displacement_direction_all = copy.copy(self.folder_displacement_direction)
        if self.verbose:
            print "     self.sc_string_params:",self.sc_string_params
            print "     self.eqcoordsfile_getpot:",self.eqcoordsfile_getpot
            print "     self.folder_displacement:",self.folder_displacement
            print "     self.folder_displacement_direction_all:"
            for jkl in self.folder_displacement_direction_all:
                print "             -> ",jkl
            print utils.printblue("     init_variables complete")
            print ""
        return

        # get alat, scalat
        alat = False  # only necessary for self.jobvorlage
        if self.sc_create == 2:
            alat = pot.get_element_2x2x2sc_info(self.element)[0]
        if alat == False:
            print "self.sc_create:",self.sc_create
            print "alat:",alat
            print "self.element:",self.element
            sys.exit("did not find alat v max; alat:False")
        if self.verbose == True:
            print "     alat:",alat
        self.scalat = self.sc_create*alat
        self.cell = np.array([[self.scalat,0.0,0.0],[0.0,self.scalat,0.0],[0.0,0.0,self.scalat]])
        if self.verbose == True:
            print "     self.scalat:",self.scalat
            #print "     self.cell:",self.cell

        self.nndist = alat/np.sqrt(2.)
        if self.verbose:
            print "     self.nndist:",self.nndist

        # check if self.jobvorlage is available
        self.sc_string_create = str(self.sc_create)+"x"+str(self.sc_create)+"x"+str(self.sc_create)
        self.jobvorlage = self.folder_base_jobvorlage+"/"+self.sc_string_create+"sc_"+self.element+"_"+str(alat)
        if self.verbose:
            print "     self.jobvorlage   :",self.jobvorlage
        if os.path.isdir(self.jobvorlage) != True:
            sys.exit(self.jobvorlage + " self.jobvorlage does not exist 3 !")


        ######################################################################
        # tveclonall.dat
        self.jobvorlage_vecnormlon_file = self.jobvorlage+"/tveclonall.dat"   # somewhat shorter
        self.jobvorlage_vecnormlon_file = self.jobvorlage+"/atoms_1nn_all_f/dfn_1.0"
        if os.path.isfile(self.jobvorlage_vecnormlon_file) != True:
            sys.exit(self.jobvorlage_vecnormlon_file + " self.jobvorlage_vecnormlon_file does not exist 3 !")
        if self.verbose:
            print "     self.jobvorlage_vecnormlon_file   :",self.jobvorlage_vecnormlon_file

        ######################################################################
        # DOSlon  DOSDOSDOS
        self.jobvorlage_DOSlon_file = self.jobvorlage+"/DOSlon"
        self.jobvorlage_DOSlon_file = self.jobvorlage+"/DOS_dfn_1.0_py3.0"
        self.jobvorlage_DOSlon_filecut = self.jobvorlage+"/DOS_dfn_1.0_py3.0cut"
        self.jobvorlage_DOSlon_filecutshiftednndist = self.jobvorlage+"/DOS_dfn_1.0_py3.0cutshiftednndist"
        if os.path.isfile(self.jobvorlage_DOSlon_file) != True:
            sys.exit(self.jobvorlage_DOSlon_file + " self.jobvorlage_DOSlon_file does not exist 3 !")
        self.jobvorlage_DOSlon = np.loadtxt(self.jobvorlage_DOSlon_file)
        self.jobvorlage_DOSloncut = np.loadtxt(self.jobvorlage_DOSlon_filecut)
        self.jobvorlage_DOSlon_max = self.jobvorlage_DOSlon[:,0].max()
        self.jobvorlage_DOSlon_min = self.jobvorlage_DOSlon[:,0].min()
        self.jobvorlage_DOSloncut_max = self.jobvorlage_DOSloncut[:,0].max()
        self.jobvorlage_DOSloncut_min = self.jobvorlage_DOSloncut[:,0].min()
        self.jobvorlage_DOSloncutshiftednndist = np.array([self.jobvorlage_DOSloncut[:,0]-self.nndist,self.jobvorlage_DOSloncut[:,1]]).transpose()
        if False:
            np.savetxt(self.jobvorlage_DOSlon_filecutshiftednndist,self.jobvorlage_DOSloncutshiftednndist)
        if self.verbose:
            print "     self.jobvorlage_DOSlon: OK"
            print "     self.jobvorlage_DOSloncut: OK"
            print "     self.jobvorlage_DOSlon_file   :",self.jobvorlage_DOSlon_file,self.jobvorlage_DOSlon_min,self.jobvorlage_DOSlon_max
            print "     self.jobvorlage_DOSlon_filecut:",self.jobvorlage_DOSlon_filecut,self.jobvorlage_DOSloncut_min,self.jobvorlage_DOSloncut_max

        ###########################################################################################
        # POSITIONs
        ###########################################################################################
        #print "|||||||||||||||",self.folder_displacement_direction_all
        for pos in self.folder_displacement_direction_all:
            #print "pos:",pos
            positions = pos+"/POSITIONs"
            if os.path.isfile(positions) != True:
                sys.exit(positions+" positionsfile does not exist")


        ###########################################################################################
        # get DOSlonpy
        ###########################################################################################
        if False:
            data = np.loadtxt(self.jobvorlage_vecnormlon_file)[2:]
            #data = np.loadtxt("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/tveclonall.dat")
            dist_space = np.linspace( min(data), max(data), 500 )
            fact = 1.0
            def my_kde_bandwidth(obj, fac=fact):
                """We use Scott's Rule, multiplied by a constant factor."""
                return np.power(obj.n, -1./(obj.d+4)) * fac

            from scipy.stats.kde import gaussian_kde
            kde = gaussian_kde( data ,bw_method=my_kde_bandwidth)
            #kde = gaussian_kde( data )
            np.savetxt(self.jobvorlage+"/DOS_dfn_1.0_py"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())
            #np.savetxt(self.jobvorlage+"/DOSlonpy"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())
            #np.savetxt("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/t"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())

        ###########################################################################################
        # get DOS  (make DOS_dfn_1.0_py3.0cut  &&  DOS_dfn_1.0_py3.0cutconservative )
        ###########################################################################################
        if False:
            # takes time to load in, just load in if necessary
            #self.jobvorlage_vecnormlon = np.loadtxt(self.jobvorlage_vecnormlon_file)[2:]
            self.jobvorlage_vecnormlon = np.loadtxt(self.jobvorlage_vecnormlon_file)
            lonmin = round(self.jobvorlage_vecnormlon.min()-0.01,2)
            lonmax = round(self.jobvorlage_vecnormlon.max()+0.01,2)
            #lonmin = self.jobvorlage_vecnormlon.min()
            #lonmax = self.jobvorlage_vecnormlon.max()
            print "lonmin/max:",lonmin,lonmax
            #print self.jobvorlage_DOSlon
            #print "smaller----"
            smaller = np.nonzero(self.jobvorlage_DOSlon[:,0]<=lonmax)[0]
            greater = np.nonzero(self.jobvorlage_DOSlon[:,0]>=lonmin)[0]
            intersect = utils.schnittmenge(smaller,greater)

            greaterconservative = np.nonzero(self.jobvorlage_DOSlon[:,1]>=0.003)[0]
            #print "greaterconservative:"
            #print greaterconservative
            union = utils.union(greaterconservative,intersect)

            #print "greater:"
            #print greater
            #print "intersect:"
            #print "last--"
            #print np.array([smaller,self.jobvorlage_DOSlon[:,0][smaller],self.jobvorlage_DOSlon[:,1][smaller]]).transpose()
            #print np.array([intersect,self.jobvorlage_DOSlon[:,0][intersect],self.jobvorlage_DOSlon[:,1][intersect]]).transpose()
            #DOSloncut = np.array([self.jobvorlage_DOSlon[:,0][intersect],self.jobvorlage_DOSlon[:,1][intersect]]).transpose()
            DOSloncut = np.array([self.jobvorlage_DOSlon[:,0][greaterconservative],self.jobvorlage_DOSlon[:,1][greaterconservative]]).transpose()
            np.savetxt(self.jobvorlage_DOSlon_file+"cut",DOSloncut)
            DOSloncutshifted = np.array([self.DOSloncut[:,0]-self.nndist,DOSloncut[:,1]]).transpose()
            np.savetxt(self.jobvorlage_DOSlon_file+"cutshifted",DOSloncutshifted)
            #DOSlonconservative = np.array([self.jobvorlage_DOSlon[:,0][union],self.jobvorlage_DOSlon[:,1][union]]).transpose()
            #np.savetxt(self.jobvorlage_DOSlon_file+"cutconservative",DOSlonconservative)


        ###########################################################################################
        # get DOS shifted to max point
        ###########################################################################################
        if False:
            maximum = self.jobvorlage_DOSlon[:,1].max()
            print "maximum:",maximum
            maximumidx = np.nonzero(self.jobvorlage_DOSlon[:,1]==maximum)[0][0]
            print "maximumidx:",maximumidx
            maximum_xvalue = self.jobvorlage_DOSlon[maximumidx][0]
            print "maximum_xvalue:",maximum_xvalue
            DOSlonshifted0 = np.array([self.jobvorlage_DOSlon[:,0]-maximum_xvalue,self.jobvorlage_DOSlon[:,1]]).transpose()
            np.savetxt(self.jobvorlage_DOSlon_file+"shifted0",DOSlonshifted0)

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
            self.disp_lon_all_diffcut = utils.cut_function_at_DOS(self.disp_lon_all_diff,self.jobvorlage_DOSloncut)
            self.disp_lon_all_diffcut_file = self.disp_lon_all_diff_file+"_cut"
            self.disp_lon_all_diffcutshifted = np.array([self.disp_lon_all_diffcut[:,0]-self.nndist,self.disp_lon_all_diffcut[:,1]]).transpose()
            self.disp_lon_all_diffcutshifted_file = self.disp_lon_all_diffcut_file+ "shifted0"


        if False:
            np.savetxt(self.disp_lon_all_diffcut_file,self.disp_lon_all_diffcut)
            np.savetxt(self.disp_lon_all_diffcutshifted_file,self.disp_lon_all_diffcutshifted)


        ###########################################################################################
        # general checks
        ###########################################################################################

        utils.isfile(self.jobvorlage+"/HesseMatrix_sphinx")
        utils.isfile(self.jobvorlage+"/INCAR")
        utils.isfile(self.jobvorlage+"/KPOINTS")

        self.POTCAR = "/Users/glensk/Thermodynamics/vasp_potentials/PAW-GGA-PBE_vasp4.6__DO-NOT-USE--TAKE-NEW-POTCARs-FROM-vasp5.2-INSTEAD/"+self.element+"/POTCAR"
        utils.isfile(self.POTCAR)

        self.eqcoordsfile_job = self.utilsfolder+"/"+self.structure+"/EqCoords_direct_"+self.structure+"_"+self.sc_string_create+"sc"
        utils.isfile(self.eqcoordsfile_job)



        if self.verbose:
            print utils.printblue("     init_variables complete")
            print ""
        return




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


    def pot_parametrize_function( self, onefolder = False ):
        ''' goes through displacements (xdir,quer) and fits forces to lon,tox function '''
        ''' goes through displacements (xdir,quer) and fits forces to lon,tox function '''
        if type(onefolder) != bool:
            self.folder_displacement_direction = [ onefolder ]
        if self.verbose:
            print ""
            print utils.printblue("evaluate_pot ...")
            print "     self.folder_displacement_direction:",self.folder_displacement_direction

        hier=os.getcwd()
        for folder in self.folder_displacement_direction:
            print utils.printgreen("     folder:"+folder)
            if os.path.isfile(folder+"/EqCoords_direct") != True:
                shutil.copyfile(self.eqcoordsfile_getpot,folder+"/EqCoords_direct")

            if type(self.alat) == bool:
                self.alatsstringlist = self.get_alats(folder)
            else:
                self.alatsstringlist = [ self.alat ]

            os.chdir(folder)
            for a in self.alatsstringlist:
                self.n = pot_parametrize.forcesneighbors()
                self.n.a = float(a)
                self.n.verbose = False
                print "## 11 n.loadforces:"
                self.n.DOScut = False
                if self.DOSparametrization:
                    self.n.nnforcesaddname = "DOS"    # this is the dos which take all the points from the DOS of the MD but not more!
                    self.n.DOScut = self.jobvorlage_DOSlon
                if self.DOScutparametrization:
                    self.n.nnforcesaddname= "DOScut"  # this is the dos with slightly less points
                    self.n.DOScut = self.jobvorlage_DOSloncut
                self.n.loadforces() #DOScut = True, DOScutarray = self.jobvorlage_DOSloncut)
                print "## 22  n.getforces:, self.n.nnforcesaddname:",self.n.nnforcesaddname,self.n.DOScut.shape
                self.n.getforces()
                print "## 33 ---"
        os.chdir(hier)
        print utils.printblue("evaluate_pot complete...")
        return

    def get_parametrization_folder(self, disp = False, stringadd = False):
        if type(disp) == bool or type (stringadd) == bool:
            return False
        #if self.verbose:
        #    print ""
        #    print "get_parametrization_folder ..."
        searchfor = self.folder_displacement+"/"+self.sc_string_params+"sc_"+disp+"_"+stringadd
        #if self.verbose:
        #    print "     searchfor:",searchfor
        dispfolder = glob.glob(searchfor)
        if len(dispfolder) != 1:
            print ""
            print utils.printred("get_parametrization_folder ...")
            print "-----------ERROR----------------"
            print "sfd:",self.folder_displacement
            print "self.sc_string_arams:",self.sc_string_params
            print "disp:",disp
            print "stringadd:",stringadd
            print "searchfor:",searchfor
            print "dispfolder:",dispfolder
            if len(dispfolder) == 0:
                sys.exit("expected one dispfolder, got none")
            if len(dispfolder) != 0:
                sys.exit("expected one dispfolder, got several")
        self.folder_displacement_direction = dispfolder[0]
        #print "dispfolder:",dispfolder
        if os.path.isdir(self.folder_displacement_direction) != True:
            sys.exit(self.folder_displacement_direction+" not found self.folder_displacement_direction")
        #print "ddd:",self.folder_displacement_direction
        return self.folder_displacement_direction


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




    def get_parametrization(self, disp = False, stringadd = False, parametrization_string = False, contr = False):
        if type(contr) == bool:
            print "contr:",contr
            sys.exit("which contribution?")
        parametrization_folder =  self.get_parametrization_folder(disp = disp, stringadd = stringadd)
        searchfor = parametrization_folder+"/nnforces/"+contr+"_1nn_"+parametrization_string+"*"
        searchforallp = parametrization_folder+"/nnforces/"+contr+"_1nn_all"
        par_lon_file = glob.glob(searchfor)
        par_lon_fileall = glob.glob(searchforallp)
        if len(par_lon_fileall) != 1:
            sys.exit(searchforallp+" not found")
        par_lon_fileall = par_lon_fileall[0]
        if os.path.isfile(par_lon_fileall) != True:
            print "searchforall:",searchforall
            sys.exit("original file not found")
        if len(par_lon_file) != 1:
            print "parametrization_folder:",parametrization_folder
            print "searchfor:",searchfor
            print "disp:",disp
            print "stringadd:",stringadd
            print "par_lon_file:",par_lon_file
            sys.exit("expected one file!")
        par_lon_file = par_lon_file[0]
        #print par_lon_file
        #par_lon_params, par_lon_pot_check = self.get_parametrization_parameters_from_file(par_lon_file)
        par_lon_params, par_lon_pot = self.get_parametrization_parameters_from_file(par_lon_file)
        #print "par_lon_params:",par_lon_params
        #par_lon_pot = False
        #if contr == 'lon':
        #    #if self.verbose:
        #    #    print "     parstr:",parametrization_string
        #    if "_poly" in parametrization_string:
        #        par_lon_pot = "poly"
        #    if "_morse" in parametrization_string:
        #        par_lon_pot = "m"
        #    if "_mc1" in parametrization_string:
        #        par_lon_pot = "mc1"
        #    if par_lon_pot == False:
        #        sys.exit("par_lon_pot not found")
        #if par_lon_pot != par_lon_pot_check:
        #    sys.exit("pot check 55")
        return par_lon_params,par_lon_file,parametrization_folder,par_lon_pot,par_lon_fileall


    def check_calculate_energy_and_forces_parameters_with_current(self,calculate_energy_and_forces_path):
        import imp
        calculate_energy_and_forces = imp.load_source('calculate_energy_and_forces',calculate_energy_and_forces_path)
        from calculate_energy_and_forces import \
                u1nn_pot, \
                u1nn_potparam, \
                u1nn_pot2, \
                u1nn_potparam2, \
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
        if u1nn_potparam2 != self.par_lon2_params:
            sys.exit("55 lon2")
        if u1nn_topx != self.par_tox_params:
            sys.exit("55 tox")
        if u1nn_topy != self.par_toy_params:
            sys.exit("55 toy")
        if u1nn_topz != self.par_toz_params:
            sys.exit("55 toz")
        return


    def get_par_lon(self, disp = False, stringadd = False, parametrization_string = False):
        if self.verbose:
            print utils.printblue("get_par_lon ...")
        self.lon_disp = disp
        self.lon_stringadd = stringadd
        self.lon_parametrization_string = parametrization_string

        # in case there is a LON2 this will be changed in get_par_lon2r_lon2() function
        self.lon2_disp = False
        self.lon2_stringadd = False
        self.lon2_parametrization_string = False
        self.par_lon_params,self.par_lon_file,self.parametrization_folder, self.par_lon_pot, self.par_lon_file_orig = \
                self.get_parametrization(disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = 'lon')


    def get_par_lon_add(self, disp = False, stringadd = False, parametrization_string = False):
        ''' disp: 'quer', 'xdir'                           # here the parametrization cames from
            stringadd: 3x3x3kp_vasp4, 10x10x10kp, ...      # here the parametrization cames from
            parametrization_string: neg_rlv_0.15_mc1       # here the parametrization cames from
            '''
        if self.verbose:
            print utils.printblue("get_par_lon_add ...")
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
        if type(disp) == bool and type(stringadd) == bool and type(parametrization_string) == bool:
            self.par_lon_add_params,self.par_lon_add_file,self.parametrization_folder, self.par_lon_add_pot, self.par_lon_add_file_orig = False, False, False, False, False
            return
        else:
            self.par_lon_add_params,self.par_lon_add_file,self.parametrization_folder, self.par_lon_add_pot, self.par_lon_add_file_orig = \
                self.get_parametrization(disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = 'lonadd')


    def get_par_lon2_add(self, disp = False, stringadd = False, parametrization_string = False):
        ''' disp: 'quer', 'xdir'                           # here the parametrization cames from
            stringadd: 3x3x3kp_vasp4, 10x10x10kp, ...      # here the parametrization cames from
            parametrization_string: neg_rlv_0.15_mc1       # here the parametrization cames from
            '''
        if self.verbose:
            print utils.printblue("get_par_lon2_add ...")
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
        if type(disp) == bool and type(stringadd) == bool and type(parametrization_string) == bool:
            self.par_lon2_add_params,self.par_lon2_add_file,self.parametrization_folder, self.par_lon2_add_pot, self.par_lon2_add_file_orig = False, False, False, False, False
            return
        else:
            if self.par_lon2_pot == 'poly':
                sys.exit("you cant have self.ar_lon2_pot == 'poly' and add another poly (not yet)")
            self.par_lon2_add_params,self.par_lon2_add_file,self.parametrization_folder, self.par_lon2_add_pot, self.par_lon2_add_file_orig = \
                self.get_parametrization(disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = 'lon2add')



    def get_par_lon2(self, disp = False, stringadd = False, parametrization_string = False):
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
            self.parametrization_folder = False
        else:
            self.par_lon2_params, self.par_lon2_file, self.parametrization_folder, self.par_lon2_pot, self.par_lon2_file_orig = \
                    self.get_parametrization(disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = 'lon')

    def get_par_tox(self, disp = False, stringadd = False, parametrization_string = False, contr = 'tox'):
        if self.verbose:
            print utils.printblue("get_par_tox ...")
        self.tox_disp = disp
        self.tox_stringadd = stringadd
        self.tox_parametrization_string = parametrization_string
        if type(disp) == bool:
            self.par_tox_params = False
            self.par_tox_file = False
            self.par_tox_file_orig = False
            self.parametrization_folder = False
        else:
            self.par_tox_params,self.par_tox_file,self.parametrization_folder, self.par_tox_pot, self.par_tox_file_orig = \
                    self.get_parametrization(disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = contr )

    def get_par_toy(self, disp = False, stringadd = False, parametrization_string = False, contr = 'toy'):
        if self.verbose:
            print utils.printblue("get_par_toy ...")
        self.toy_disp = disp
        self.toy_stringadd = stringadd
        self.toy_parametrization_string = parametrization_string
        if type(disp) == bool:
            self.par_toy_params = False
            self.par_toy_file = False
            self.par_toy_file_orig = False
            self.parametrization_folder = False
        else:
            self.par_toy_params,self.par_toy_file,self.parametrization_folder, self.par_toy_pot, self.par_toy_file_orig = \
                    self.get_parametrization(disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = contr )

    def get_par_toz(self, disp = False, stringadd = False, parametrization_string = False, contr = 'toz'):
        if self.verbose:
            print utils.printblue("get_par_toz ...")
        self.toz_disp = disp
        self.toz_stringadd = stringadd
        self.toz_parametrization_string = parametrization_string
        if type(disp) == bool:
            self.par_toz_params = False
            self.par_toz_file = False
            self.par_toz_file_orig = False
            self.parametrization_folder = False
        else:
            self.par_toz_params,self.par_toz_file,self.parametrization_folder, self.par_toz_pot, self.par_toz_file_orig = \
                    self.get_parametrization(disp = disp, stringadd = stringadd, parametrization_string = parametrization_string, contr = contr )

    def print_parametrizations(self):
        print ""
        print utils.printblue("print_parametrizations ...")
        print "lon :",self.par_lon_pot,self.par_lon_params
        print "lonadd:",self.par_lon_add_pot,self.par_lon_add_params
        print "lon2:",self.par_lon2_pot,self.par_lon2_params
        if type(self.par_tox_params) != bool:
            print "tox :",self.par_tox_pot,self.par_tox_params
        if type(self.par_toy_params) != bool:
            print "toy :",self.par_toy_pot,self.par_toy_params
        if type(self.par_toz_params) != bool:
            print "toz :",self.par_toz_pot,self.par_toz_params
        print ""
        print " lono:",self.par_lon_file_orig
        print " lon :",self.par_lon_pot,self.par_lon_file
        print " lon2:",self.par_lon2_pot,self.par_lon2_file
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
        if os.path.isfile(filevorlage) != True:
            sys.exit(filevorlage+" not found filevorlage")
        shutil.copy2(filevorlage,filevorlagetmp)
        utils.sed(filename = filevorlagetmp,stringsearch='^LONPOT=None.*',stringreplace='LONPOT=\''+self.par_lon_pot+'\'')
        utils.sed(filename = filevorlagetmp,stringsearch='^LONPARAMETRIZATION=None.*',stringreplace='LONPARAMETRIZATION=\''+self.par_lon_params+'\'')


        if type(self.par_lon_add_pot) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^LONPOTADD=None.*',stringreplace='LONPOTADD='+str(self.par_lon_add_pot)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^LONPOTADD=None.*',stringreplace='LONPOTADD=\''+str(self.par_lon_add_pot)+'\'')


        if type(self.par_lon_add_params) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^LONADDPARAMETRIZATION=None.*',stringreplace='LONADDPARAMETRIZATION='+str(self.par_lon_add_params)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^LONADDPARAMETRIZATION=None.*',stringreplace='LONADDPARAMETRIZATION=\''+str(self.par_lon_add_params)+'\'')


        if type(self.par_lon2_add_pot) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^LON2POTADD=None.*',stringreplace='LON2POTADD='+str(self.par_lon2_add_pot)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^LON2POTADD=None.*',stringreplace='LON2POTADD=\''+str(self.par_lon2_add_pot)+'\'')

        if type(self.par_lon2_add_params) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^LON2ADDPARAMETRIZATION=None.*',stringreplace='LON2ADDPARAMETRIZATION='+str(self.par_lon2_add_params)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^LON2ADDPARAMETRIZATION=None.*',stringreplace='LON2ADDPARAMETRIZATION=\''+str(self.par_lon2_add_params)+'\'')


        if type(self.par_lon2_pot) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^LON2POT=None.*',stringreplace='LON2POT='+str(self.par_lon2_pot)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^LON2POT=None.*',stringreplace='LON2POT=\''+str(self.par_lon2_pot)+'\'')

        if type(self.par_lon2_params) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^LON2PARAMETRIZATION=None.*',stringreplace='LON2PARAMETRIZATION='+str(self.par_lon2_params)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^LON2PARAMETRIZATION=None.*',stringreplace='LON2PARAMETRIZATION=\''+str(self.par_lon2_params)+'\'')


        if type(self.par_tox_params) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^TOXPARAMETRIZATION=None.*',stringreplace='TOXPARAMETRIZATION='+str(self.par_tox_params)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^TOXPARAMETRIZATION=None.*',stringreplace='TOXPARAMETRIZATION=\''+self.par_tox_params+'\'')

        if type(self.par_toy_params) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^TOYPARAMETRIZATION=None.*',stringreplace='TOYPARAMETRIZATION='+str(self.par_toy_params)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^TOYPARAMETRIZATION=None.*',stringreplace='TOYPARAMETRIZATION=\''+self.par_toy_params+'\'')

        if type(self.par_toz_params) == bool:
            utils.sed(filename = filevorlagetmp,stringsearch='^TOZPARAMETRIZATION=None.*',stringreplace='TOZPARAMETRIZATION='+str(self.par_toz_params)+'')
        else:
            utils.sed(filename = filevorlagetmp,stringsearch='^TOZPARAMETRIZATION=None.*',stringreplace='TOZPARAMETRIZATION=\''+self.par_toz_params+'\'')

        return filevorlagetmp


    def get_job_number(self):
        jobnumber = ""
        self.jobnumber = False
        alljobs = glob.glob(self.folder_job_tdi+"[0-9]*")
        alljobsnrs = [ i.split(self.folder_job_tdi)[1].split("_")[0] for i in alljobs ]
        alljobsnrs = list(set(alljobsnrs))
        if self.verbose:
            print "     alljobsnrs:",alljobsnrs
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


    def get_job_name(self):
        lonpart = "from_"+self.sc_string_params+"sc__LON_"+self.lon_disp+"_"+self.lon_stringadd+"_"+self.lon_parametrization_string
        if self.par_lon2_params == False:
            lon2part = "_LON2_NONE"
        else:
            lon2part = "_LON2_"+self.lon2_disp+"_"+self.lon2_stringadd+"_"+self.lon2_parametrization_string

        if type(self.tox_disp) == bool or type(self.tox_stringadd) == bool or type(self.tox_parametrization_string) == bool:
            toxpart = "_TOX_NONE"
        else:
            toxpart = "_TOX_"+self.tox_disp+"_"+self.tox_stringadd+"_"+self.tox_parametrization_string

        if type(self.toy_disp) == bool or type(self.toy_stringadd) == bool or type(self.toy_parametrization_string) == bool:
            toypart = "_TOY_NONE"
        else:
            toypart = "_TOY_"+self.toy_parametrization_string

        if type(self.toz_disp) == bool or type(self.toz_stringadd) == bool or type(self.toz_parametrization_string) == bool:
            tozpart = "_TOZ_NONE"
        else:
            tozpart = "_TOZ_"+self.toz_parametrization_string

        if type(self.par_lon_add_params) == bool:
            lonaddpart = "_LONADD__NO"
        else:
            lonaddpart = "_LONADD_YES"

        if type(self.par_lon2_add_params) == bool:
            lon2addpart = "_LON2ADD__NO"
        else:
            lon2addpart = "_LON2ADD_YES"

        doscut = "PTS_allpts_"
        if self.DOSparametrization == True:
            doscut = "PTS_dosall_"
        if self.DOScutparametrization == True:
            doscut = "PTS_doscut_"

        self.jobnamestr = doscut + lonpart + lonaddpart + lon2part + lon2addpart + toxpart + toypart + tozpart
        self.jobnamestr = self.jobnamestr.replace("order","")
        self.jobnamestr = self.jobnamestr.replace("_fit","")
        return self.jobnamestr

    def get_job_info(self,ls = False, execute = False, jobnumber = False, verbose = False):
        ''' if at some point we will calculate ohter alats, jobnumber will have to also take care of this
            currently no mater which jobnumber, the alat ist the same '''
        if self.verbose:
            print ""
            print utils.printblue("get_job_info ...")


        # get alat, scalat
        alat = False  # only necessary for self.jobvorlage
        if self.sc_create == 2:
            alat = pot.get_element_2x2x2sc_info(self.element)[0]
        if alat == False:
            print "self.sc_create:",self.sc_create
            print "alat:",alat
            print "self.element:",self.element
            sys.exit("did not find alat v max; alat:False")
        if verbose == True:
            print ""
            print "v--------------- get_job_info -------------------v"
            print "alat:",alat
        self.scalat = self.sc_create*alat
        self.cell = np.array([[self.scalat,0.0,0.0],[0.0,self.scalat,0.0],[0.0,0.0,self.scalat]])
        if verbose == True:
            print "self.scalat:",self.scalat
            print "self.cell:",self.cell

        self.nndist = alat/np.sqrt(2.)
        if self.verbose:
            print "     alat:",alat
            print "     self.nndist:",self.nndist

        # check if self.jobvorlage is available
        self.sc_string_create = str(self.sc_create)+"x"+str(self.sc_create)+"x"+str(self.sc_create)
        self.jobvorlage = self.folder_base_jobvorlage+"/"+self.sc_string_create+"sc_"+self.element+"_"+str(alat)
        if self.verbose:
            print "     self.jobvorlage   :",self.jobvorlage
        if os.path.isdir(self.jobvorlage) != True:
            sys.exit(self.jobvorlage + " self.jobvorlage does not exist 3 !")

        ######################################################################
        # tveclonall.dat
        self.jobvorlage_vecnormlon_file = self.jobvorlage+"/tveclonall.dat"   # somewhat shorter
        self.jobvorlage_vecnormlon_file = self.jobvorlage+"/atoms_1nn_all_f/dfn_1.0"
        if os.path.isfile(self.jobvorlage_vecnormlon_file) != True:
            sys.exit(self.jobvorlage_vecnormlon_file + " self.jobvorlage_vecnormlon_file does not exist 3 !")
        if self.verbose:
            print "     self.jobvorlage_vecnormlon_file   :",self.jobvorlage_vecnormlon_file

        ######################################################################
        # DOSlon
        self.jobvorlage_DOSlon_file = self.jobvorlage+"/DOSlon"
        self.jobvorlage_DOSlon_file = self.jobvorlage+"/DOS_dfn_1.0_py3.0"
        self.jobvorlage_DOSlon_filecut = self.jobvorlage+"/DOS_dfn_1.0_py3.0cut"
        if os.path.isfile(self.jobvorlage_DOSlon_file) != True:
            sys.exit(self.jobvorlage_DOSlon_file + " self.jobvorlage_DOSlon_file does not exist 3 !")
        self.jobvorlage_DOSlon = np.loadtxt(self.jobvorlage_DOSlon_file)
        self.jobvorlage_DOSloncut = np.loadtxt(self.jobvorlage_DOSlon_filecut)
        self.jobvorlage_DOSlon_max = self.jobvorlage_DOSlon[:,0].max()
        self.jobvorlage_DOSlon_min = self.jobvorlage_DOSlon[:,0].min()
        self.jobvorlage_DOSloncut_max = self.jobvorlage_DOSloncut[:,0].max()
        self.jobvorlage_DOSloncut_min = self.jobvorlage_DOSloncut[:,0].min()
        if self.verbose:
            print "     self.jobvorlage_DOSlon_file   :",self.jobvorlage_DOSlon_file,self.jobvorlage_DOSlon_min,self.jobvorlage_DOSlon_max
            print "     self.jobvorlage_DOSlon_filecut:",self.jobvorlage_DOSlon_filecut,self.jobvorlage_DOSloncut_min,self.jobvorlage_DOSloncut_max




        ###########################################################################################
        # POSITIONs
        ###########################################################################################
        #print "|||||||||||||||",self.folder_displacement_direction_all
        for pos in self.folder_displacement_direction_all:
            #print "pos:",pos
            positions = pos+"/POSITIONs"
            if os.path.isfile(positions) != True:
                sys.exit(positions+" positionsfile does not exist")


        ###########################################################################################
        # get DOSlonpy
        ###########################################################################################
        if False:
            data = np.loadtxt(self.jobvorlage_vecnormlon_file)[2:]
            #data = np.loadtxt("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/tveclonall.dat")
            dist_space = np.linspace( min(data), max(data), 500 )
            fact = 1.0
            def my_kde_bandwidth(obj, fac=fact):
                """We use Scott's Rule, multiplied by a constant factor."""
                return np.power(obj.n, -1./(obj.d+4)) * fac

            from scipy.stats.kde import gaussian_kde
            kde = gaussian_kde( data ,bw_method=my_kde_bandwidth)
            #kde = gaussian_kde( data )
            np.savetxt(self.jobvorlage+"/DOS_dfn_1.0_py"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())
            #np.savetxt(self.jobvorlage+"/DOSlonpy"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())
            #np.savetxt("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/t"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())

        ###########################################################################################
        # get DOS  (make DOS_dfn_1.0_py3.0cut  &&  DOS_dfn_1.0_py3.0cutconservative )
        ###########################################################################################
        if False:
            # takes time to load in, just load in if necessary
            #self.jobvorlage_vecnormlon = np.loadtxt(self.jobvorlage_vecnormlon_file)[2:]
            self.jobvorlage_vecnormlon = np.loadtxt(self.jobvorlage_vecnormlon_file)
            lonmin = round(self.jobvorlage_vecnormlon.min()-0.01,2)
            lonmax = round(self.jobvorlage_vecnormlon.max()+0.01,2)
            #lonmin = self.jobvorlage_vecnormlon.min()
            #lonmax = self.jobvorlage_vecnormlon.max()
            print "lonmin/max:",lonmin,lonmax
            #print self.jobvorlage_DOSlon
            #print "smaller----"
            smaller = np.nonzero(self.jobvorlage_DOSlon[:,0]<=lonmax)[0]
            greater = np.nonzero(self.jobvorlage_DOSlon[:,0]>=lonmin)[0]
            intersect = utils.schnittmenge(smaller,greater)

            greaterconservative = np.nonzero(self.jobvorlage_DOSlon[:,1]>=0.003)[0]
            #print "greaterconservative:"
            #print greaterconservative
            union = utils.union(greaterconservative,intersect)

            #print "greater:"
            #print greater
            #print "intersect:"
            #print "last--"
            #print np.array([smaller,self.jobvorlage_DOSlon[:,0][smaller],self.jobvorlage_DOSlon[:,1][smaller]]).transpose()
            #print np.array([intersect,self.jobvorlage_DOSlon[:,0][intersect],self.jobvorlage_DOSlon[:,1][intersect]]).transpose()
            #DOSloncut = np.array([self.jobvorlage_DOSlon[:,0][intersect],self.jobvorlage_DOSlon[:,1][intersect]]).transpose()
            DOSloncut = np.array([self.jobvorlage_DOSlon[:,0][greaterconservative],self.jobvorlage_DOSlon[:,1][greaterconservative]]).transpose()
            np.savetxt(self.jobvorlage_DOSlon_file+"cut",DOSloncut)
            #DOSlonconservative = np.array([self.jobvorlage_DOSlon[:,0][union],self.jobvorlage_DOSlon[:,1][union]]).transpose()
            #np.savetxt(self.jobvorlage_DOSlon_file+"cutconservative",DOSlonconservative)


        ###########################################################################################
        # get DOS shifted to max point
        ###########################################################################################
        if False:
            maximum = self.jobvorlage_DOSlon[:,1].max()
            print "maximum:",maximum
            maximumidx = np.nonzero(self.jobvorlage_DOSlon[:,1]==maximum)[0][0]
            print "maximumidx:",maximumidx
            maximum_xvalue = self.jobvorlage_DOSlon[maximumidx][0]
            print "maximum_xvalue:",maximum_xvalue
            DOSlonshifted0 = np.array([self.jobvorlage_DOSlon[:,0]-maximum_xvalue,self.jobvorlage_DOSlon[:,1]]).transpose()
            np.savetxt(self.jobvorlage_DOSlon_file+"shifted0",DOSlonshifted0)



        ###########################################################################################
        # general checks
        ###########################################################################################

        utils.isfile(self.jobvorlage+"/HesseMatrix_sphinx")
        utils.isfile(self.jobvorlage+"/INCAR")
        utils.isfile(self.jobvorlage+"/KPOINTS")

        self.POTCAR = "/Users/glensk/Thermodynamics/vasp_potentials/PAW-GGA-PBE_vasp4.6__DO-NOT-USE--TAKE-NEW-POTCARs-FROM-vasp5.2-INSTEAD/"+self.element+"/POTCAR"
        utils.isfile(self.POTCAR)

        # check that u is mounted
        if os.path.isdir(self.folder_job_tdi_base) != True:
            sys.exit(self.folder_job_tdi_base+" not found self.folder_job_tdi_base")


        self.eqcoordsfile_job = self.utilsfolder+"/"+self.structure+"/EqCoords_direct_"+self.structure+"_"+self.sc_string_create+"sc"
        utils.isfile(self.eqcoordsfile_job)

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
        # EXISTING (OLD) JOB
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

            self.par_lon_params,self.par_lon_file,self.parametrization_folder, self.par_lon_pot, self.par_lon_file_orig =     "", "", "", "", ""
            self.par_lon2_params,self.par_lon2_file,self.parametrization_folder, self.par_lon2_pot, self.par_lon2_file_orig = "", "", "", "", ""
            self.par_tox_params,self.par_tox_file,self.parametrization_folder, self.par_tox_pot, self.par_tox_file_orig =     "", "", "", "", ""
            self.par_toy_params,self.par_toy_file,self.parametrization_folder, self.par_toy_pot, self.par_toy_file_orig =     "", "", "", "", ""
            self.par_toz_params,self.par_toz_file,self.parametrization_folder, self.par_toz_pot, self.par_toz_file_orig =     "", "", "", "", ""

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
        # NEWJOB
        ##########################################################################################################
        if self.jobnew == True:
            # create file calculate_energy_and_forces_tmp
            self.calculate_energy_and_forces_pathtmp = self.create_calculate_energy_and_forces()
            if self.verbose == True:
                print "     self.calculate_energy_and_forces_pathtmp:",self.calculate_energy_and_forces_pathtmp

            # get job number
            self.jobnumber = self.get_job_number()
            if self.verbose == True:
                print "     self.jobnumber:",self.jobnumber


            # get job name
            self.jobnamestr = self.get_job_name()
            if self.verbose == True:
                print "     self.jobnamestr:",self.jobnamestr

            # check if self.jobnamestr does already exist
            # if self.jobpath does already exist take the jobname which already exists
            print 'xxx', self.folder_job_tdi+self.jobnumberstr+self.jobnamestr
            check1 = self.folder_job_tdi+"*"+self.jobnamestr
            check2 = glob.glob(check1)
            if len(check2) != 0:
                if len(check2) == 1:
                    self.jobpath = check2[0]
                    print utils.printred("     self.jobpath: "+self.jobpath+" does already exist (== self.jobpath)")
                if len(check2) > 1:
                    print "check1:",check1
                    print "check2:",check2
                    sys.exit("found more then one folder with same name (see check1 and check2)")


            # create new jobname since it does not exist yet
            if len(check2) == 0:
                self.jobpath = self.folder_job_tdi+self.jobnumberstr+self.jobnamestr
            if self.verbose == True:
                print utils.printgreen("     self.jobpath:"+self.jobpath)


        if self.verbose:
            print utils.printblue("     get_job_info complete ...")
            print ""
        return




    def create_job_dispdudl(self):
        if self.verbose:
            print ""
            print utils.printblue("create_job_dispdudl ...")
        alat = self.scalat/self.sc_create
        self.nndist = alat/np.sqrt(2.)
        if self.verbose:
            print "     alat:",alat
            print "     self.nndist:",self.nndist


        self.create_self_jobpath()

        for disp in [ "quer", "xdir" ]:
            parametrization_folder =  self.get_parametrization_folder(disp = disp, stringadd = self.kpstring)
            print "     parametrization_folder:",parametrization_folder
            POSITIONs = parametrization_folder+"/POSITIONs"
            dudl = parametrization_folder+"/dudl"
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
            shutil.copy2(self.jobvorlage_DOSlon_file,folder_dispdudl+"/DOSlon")
            shutil.copy2(self.jobvorlage_DOSlon_filecut,folder_dispdudl+"/DOSloncut")
            np.savetxt(folder_dispdudl+"/cell",self.cell)
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
            print "dudlnewfile:",dudlnewfile
            self.dudlnew = np.loadtxt(dudlnewfile)
            print self.dudlnew
            dudl_vs_disp = np.array([self.dudlnew[:,0]+self.nndist,self.dudlnew[:,6]]).transpose()
            np.savetxt(folder_dispdudl+"/dUdLvs_disp",dudl_vs_disp,fmt="%.3f %.4f")

        if self.verbose:
            print utils.printblue("     create_job_dispdudl complete ...")
            print ""
        return


    def create_self_jobpath(self):
        # create jobpath  (.../ti/Cu/1___from_2x2x2sc__LON_xdir...)
        if os.path.isdir(self.jobpath) != True:
            print "     creating self.jobpath:",self.jobpath
            os.makedirs(self.jobpath)
            print "     created  self.jobpath:",self.jobpath

        #print ""
        #print "self.jobpath :",self.jobpath
        if os.path.isdir(self.jobpath) != True:
            sys.exit("creation of self.jobpath faied! "+str(self.jobpath))
        return


    def create_job_tdi(self,ls = False, execute = False, verbose = False):
        print ""
        print utils.printblue("create_job_tdi ...")
        # check if all files are available for job creation

        if os.path.isdir(self.jobpath) == True:
            print utils.printred('      job exists: '+self.jobpath)
            return 'job exists: '+self.jobpath

        if type(execute) == bool:
            sys.exit("please define execute")
        if os.path.isfile(execute) != True:
            sys.exit("execute not found: "+str(execute))

        #########################################################
        # create jobpath  (creates self.jobpath)
        # e.g. /u/aglen/Dropbox/Understand_distributions/ti//Pt/3___from_2x2x2sc__LON_xdir_3x3x3kp_neg_\
        #           rlv_0.25_mc1_LON2_xdir_3x3x3kp_pos_poly9_TOX_xdir_3x3x3kp_all_poly_1st_TOY_all_poly_1st_TOZ_all_poly_1st
        #########################################################
        self.create_self_jobpath()


        POSCARVORLAGE = self.folder_base_jobvorlage+"/POSCAR_"+self.sc_string_create+"sc"
        POSCARVORLAGETMP = self.folder_base_jobvorlage+"/POSCAR_"+self.sc_string_create+"sc_tmp"
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
            shutil.copy2(self.folder_base_jobvorlage+"/hesse.py",jobpathl+"/hesse.py")
            shutil.copy2(self.eqcoordsfile_job,jobpathl+"/EqCoords_direct")
            shutil.copy2(self.jobvorlage+"/HesseMatrix_sphinx",jobpathl+"/HesseMatrix_sphinx")
            shutil.copy2(self.jobvorlage+"/INCAR",jobpathl+"/INCAR")
            #print "a:"
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR LAMBDA "+str(l)+"\'", dont_raise_exceptino = False)
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR REF_TYPE extern"+"\'", dont_raise_exceptino = False)
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR PRE_EQ_N -400"+"\'", dont_raise_exceptino = False)
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR NSW 2000"+"\'", dont_raise_exceptino = False)
            utils.run2(command = "tcsh -c \'"+"INCAR_change.sh "+jobpathl+"/INCAR NPAR 4"+"\'", dont_raise_exceptino = False)
            #print "b:"
            #shutil.copy2(self.jobvorlage+"/POSCAR",jobpathl+"/POSCAR")
            shutil.copy2(POSCARVORLAGETMP,jobpathl+"/POSCAR")
            shutil.copy2(self.jobvorlage+"/KPOINTS",jobpathl+"/KPOINTS")
            shutil.copy2(self.POTCAR,jobpathl+"/POTCAR")

            # append to jobList
            if os.path.isfile(jobList) != True:
                open(jobList, 'a').close()
            with open(jobList, "a") as myfile:
                myfile.write(jobpathl+"\n")

        #print "execute:",execute
        print ""
        utils.run2(command = "tcsh -c \'"+"ssh aglen@cmmc002.bc.rzg.mpg.de \"cd "+self.jobpath+";"+execute+";echo jo\""+"\'", dont_raise_exceptino = False)
        return

    def get_element_2x2x2sc_info(self,element):
        # max volume; kpoint in ti;
        maxv_Al = [ 4.13 , "3x3x3kp", "250" , '934' ]
        maxv_Pb = [ 5.13 , "4x4x4kp", "300" , '601' ]
        maxv_Cu = [ 3.75 , "3x3x3kp", "260" , '1360' ]
        maxv_Rh = [ 3.98 , "2x2x2kp", "270" , '2237' ]
        maxv_Pd = [ 4.1  , "3x3x3kp", "300" , '1830' ]
        maxv_Ag = [ 4.31 , "3x3x3kp", "225" , '1235' ]
        maxv_Ir = [ 3.99 , "3x3x3kp", "300" , '2739' ]
        maxv_Pt = [ 4.1  , "3x3x3kp", "200" , '2042' ]
        maxv_Au = [ 4.25 , "3x3x3kp", "300" , '1338' ]
        returninfo = "maxv_"+element
        return eval(returninfo)

    def printelement(self):
        print ""
        print utils.printred("#########################################################")
        print utils.printred("# "+ pot.element)
        print utils.printred("#########################################################")
        return

    # unimportant or just once used scripts
    def repair_disp_dudl_folder():
        ''' repair disp_dudl folder '''
        #self.folder_job_tdi: /u/aglen/Dropbox/Understand_distributions/ti//Au/
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



if __name__ == '__main__':
    hier=os.getcwd()
    pot = pot()
    pot.verbose = args.verbose

    if False:
        data = np.loadtxt("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/tveclonall.dat")
        dist_space = np.linspace( min(data), max(data), 500 )
        fact = 1.5
        def my_kde_bandwidth(obj, fac=fact):
            """We use Scott's Rule, multiplied by a constant factor."""
            return np.power(obj.n, -1./(obj.d+4)) * fac

        from scipy.stats.kde import gaussian_kde
        kde = gaussian_kde( data ,bw_method=my_kde_bandwidth)
        np.savetxt("/Users/glensk/Dropbox/Understand_distributions/jobvorlage_all/2x2x2sc_Pt_4.1/tests/t"+str(fact),np.array([dist_space, kde(dist_space)]).transpose())


    # oben
    elements = [ "Al", "Pb", "Cu", "Rh", "Pd", "Ag", "Ir", "Pt", "Au" ]
    elements = [ "Al",       "Cu", "Rh", "Pd", "Ag", "Ir", "Pt", "Au" ]
    elements = [                                           "Pt", "Au" ]
    elements = [ "Al", "Pb", "Cu", "Rh", "Pd", "Ag", "Ir", "Pt", "Au" ]
    #elements = [             "Cu", "Rh", "Pd", "Ag", "Ir", "Pt", "Au" ]
    #elements = [   "Au" ]
    #elements = [ "Al", "Pb", "Cu", "Rh", "Pd", "Ag", "Ir", "Pt", "Au" ]
    #elements = [                               "Ag", "Ir", "Pt", "Au" ]

    if args.element != False:  # z.B. al
        elements = [i.title() for i in args.element]


    for el in elements:
        pot.element    = el             # Rh, Al, ...
        pot.sc_forparams = 2            # parameters for parametrization are taken from here (2x2x2sc disp quer xdir)
        pot.sc_create = 2               # create a 2x2x2sc or a 3x3x3sc ?
        pot.printelement()
        pot.kpstring = pot.get_element_2x2x2sc_info(pot.element)[1]
        pot.init_variables( kpstringoptional = pot.kpstring)
        #pot.kpstring = "10x10x10"
        #pot.kpstring = "3x3x3kp_vasp4"
        #pot.kpstring = "3x3x3kp"
        #pot.kpstring = "10x10x10kp_vasp4"
        #pot.kpstring = "3x3x3kp_ISYM0"

        ##################################################################################
        ###
        ##################################################################################
        pot.DOSparametrization = True   # always leave true otherwise all points are taken
        pot.DOScutparametrization = False   # if this is turned on the parametrization is even smaller

        if args.sdh:
            pot.dudlnew, pot.dudl_vs_disp = start_dispdudl(nndist = pot.nndist)
            sys.exit()

        if args.p:
            pot.alat = False                                    # for pot_parametrize_function
            pot.pot_parametrize_function()

        ##################################################
        #pot.get_par_lon(    "xdir", pot.kpstring, "all_rlv_0.2_mc1" )              # parameters for parametrization
        #pot.get_par_lon(    "xdir", pot.kpstring, "all_rlv_0.25_mc1" )              # parameters for parametrization
        #pot.get_par_lon(    "quer", pot.kpstring, "all_rlv_0.25_mc1" )              # parameters for parametrization
        #pot.get_par_lon_add("quer", pot.kpstring, "neg_rlv_0.35_mc1" )              # parameters for parametrization
        #pot.get_par_lon("quer", pot.kpstring, "all_mc1" )              # parameters for parametrization
        pot.get_par_lon(    "quer", pot.kpstring, "all_mc1" )              # parameters for parametrization
        pot.get_par_lon_add("quer", pot.kpstring, "all_mc1" )              # parameters for parametrization
        #pot.get_par_lon(    "quer", pot.kpstring, "all_morse" )              # parameters for parametrization
        #pot.get_par_lon_add("quer", pot.kpstring, "all_morse" )              # parameters for parametrization
        #pot.get_par_lon_add(False )              # parameters for parametrization
        #pot.get_par_lon(   "quer", pot.kpstring, "neg_rlv_0.3_mc1" )               # parameters for parametrization
        #pot.get_par_lon(   "quer", pot.kpstring, "neg_rlv_0.25_mc1" )              # parameters for parametrization
        #pot.get_par_lon(   "xdir", pot.kpstring, "neg_rlv_0.25_mc1" )              # parameters for parametrization
        #pot.get_par_lon(   "quer", pot.kpstring, "neg_rlv_0.35_mc1" )              # parameters for parametrization
        pot.get_par_lon2(  False  )                                                # parameters for parametrization
        #pot.get_par_lon2(   "xdir", pot.kpstring, "pos_poly9" )                     # parameters for parametrization
        #pot.get_par_lon2(   "quer", pot.kpstring, "pos_poly6" )                     # parameters for parametrization
        pot.get_par_lon2_add(False )                     # parameters for parametrization
        #pot.get_par_lon2_add("quer", pot.kpstring, "pos_poly9" )                     # parameters for parametrization

        #pot.get_par_tox(False )                     # parameters for parametrization
        pot.get_par_tox(    "xdir", pot.kpstring, "all_poly1" )                     # parameters for parametrization
        #pot.get_par_toy(    "xdir", pot.kpstring, "all_poly_1storder" )            # parameters for parametrization
        #pot.get_par_toz(    "xdir", pot.kpstring, "all_poly_1storder" )            # parameters for parametrization
        #pot.get_par_toz(    False  )                                               # parameters for parametrization
        #pot.get_par_toz(    "xdir", pot.kpstring, "all_poly_1storder" )            # parameters for parametrization
        #pot.get_par_tox(    "xdir", pot.kpstring, "all_poly_9thorder" )            # parameters for parametrization
        #pot.get_par_toy(    "xdir", pot.kpstring, "all_poly_9thorder" )            # parameters for parametrization
        #pot.get_par_toz(    "xdir", pot.kpstring, "all_poly_9thorder" )            # parameters for parametrization
        #pot.get_par_tox(False)                # parameters for parametrization
        pot.get_par_toy(False)                # parameters for parametrization
        #pot.get_par_tox(    "xdir", pot.kpstring, "all_poly_1storder" )            # parameters for parametrization
        #pot.get_par_toy(    "xdir", pot.kpstring, "all_poly_1storder" )            # parameters for parametrization
        pot.get_par_toz(False)                # parameters for parametrization

        #pot.get_par_tox(    "xdir", pot.kpstring, "all_poly_9thorder" )            # parameters for parametrization
        #pot.get_par_toy(    "xdir", pot.kpstring, "all_poly_9thorder" )            # parameters for parametrization
        #pot.get_par_toz(    "xdir", pot.kpstring, "all_poly_9thorder" )            # parameters for parametrization
        if args.o:
            print "<<<<<<<<<<<<<<<<<<<<<:",args.o
        else:
            print "<<<<<<<<<<<<<<<<<<<<<:",args.o

        pot.get_job_info(jobnumber = args.o)  # has to be before pot.pot_parametrize_function() but after pot.get_par_{lon,lon2,tox,toy,toz}
        pot.print_parametrizations()

        ##################################################
        # get parametrizations
        if args.e:
            pot.alat = False                                    # for evaluate_pot
            pot.evaluate_pot()
        ##################################################

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

        ##################################################
        # create jobs
        #if args.j or args.jd:
        #pot.kpstring = "3x3x3kp_vasp4"
        #pot.kpstring = "3x3x3kp"
        #pot.kpstring = "10x10x10kp_vasp4"
        #pot.kpstring = "3x3x3kp_ISYM0"
        pot.get_par_lon(   "quer", pot.kpstring, "all_rlv_0.15_mc1" )              # parameters for parametrization
        #pot.get_par_lon(   "quer", pot.kpstring, "all_rlv_0.25_allpos_morse_" )    # parameters for parametrization
        #pot.get_par_lon(   "quer", pot.kpstring, "all_rlv_0.25_allpos_mc1_" )      # parameters for parametrization
        #pot.get_par_lon(   "quer", pot.kpstring, "neg_rlv_0.3_mc1" )               # parameters for parametrization
        #pot.get_par_lon(   "quer", pot.kpstring, "neg_rlv_0.25_mc1" )               # parameters for parametrization
        #pot.get_par_lon(   "xdir", pot.kpstring, "neg_rlv_0.25_mc1" )               # parameters for parametrization
        #pot.get_par_lon(   "quer", pot.kpstring, "neg_rlv_0.35_mc1" )              # parameters for parametrization
        #pot.get_par_lon2(  False  )                                                # parameters for parametrization
        #pot.get_par_lon2(  "xdir", pot.kpstring, "pos_poly9_fit" )                                                  # parameters for parametrization
        pot.get_par_lon2(  "quer", pot.kpstring, "pos_poly9_fit" )                                                  # parameters for parametrization

        pot.get_par_tox(    "xdir", pot.kpstring, "all_poly_1storder" )                # parameters for parametrization
        #pot.get_par_toy(    "xdir", pot.kpstring, "all_poly_1storder" )                # parameters for parametrization
        #pot.get_par_toz(    "xdir", pot.kpstring, "all_poly_1storder" )                # parameters for parametrization
        #pot.get_par_toz(    False  )                # parameters for parametrization
        #pot.get_par_toz(    "xdir", pot.kpstring, "all_poly_1storder" )                # parameters for parametrization
        #pot.get_par_tox(    "xdir", pot.kpstring, "all_poly_9thorder" )                # parameters for parametrization
        #pot.get_par_toy(    "xdir", pot.kpstring, "all_poly_9thorder" )                # parameters for parametrization
        #pot.get_par_toz(    "xdir", pot.kpstring, "all_poly_9thorder" )                # parameters for parametrization
        #pot.get_par_tox(False)                # parameters for parametrization
        pot.get_par_toy(False)                # parameters for parametrization
        #pot.get_par_tox(    "xdir", pot.kpstring, "all_poly_1storder" )                # parameters for parametrization
        #pot.get_par_toy(    "xdir", pot.kpstring, "all_poly_1storder" )                # parameters for parametrization
        pot.get_par_toz(False)                # parameters for parametrization

        #pot.get_par_tox(    "xdir", pot.kpstring, "all_poly_9thorder" )                # parameters for parametrization
        #pot.get_par_toy(    "xdir", pot.kpstring, "all_poly_9thorder" )                # parameters for parametrization
        #pot.get_par_toz(    "xdir", pot.kpstring, "all_poly_9thorder" )                # parameters for parametrization
        if args.o:
            print "<<<<<<<<<<<<<<<<<<<<<:",args.o
        else:
            print "<<<<<<<<<<<<<<<<<<<<<:",args.o

        pot.get_job_info(jobnumber = args.o)
        pot.print_parametrizations()
        ##################################################


        ##################################################
        # get alat to do # currently we just have one
        if args.k:
            if False:
                pass
            print "redo folder parametrizations:"
            to_rep_str = pot.folder_job_tdi+"/[0-8]*/"
            to_rep_folder = glob.glob(to_rep_str)
            for i in to_rep_folder:  # ( == self.jobpath )
                lonf, lon2f, toxf, toyf, tozf, lonfdot, toxfdot, toyfdot, tozfdot = get_parametrizationfiles_oldpath(i,pot.folder_displacement_direction_all)
                print "    lonf:",lonf
                print "    lon2f:",lon2f
                print "    toxf:",toxf


                pot.jobpath_parametrization = i+"/parametrization"
                print "pot.jobpath_parametrization:",pot.jobpath_parametrization
                create_job_parametrization(jobpath_parametrization = pot.jobpath_parametrization, \
                        lon = lonf,\
                        lon2 =lon2f, \
                        tox = toxf,\
                        toy = toyf,\
                        toz = tozf, \
                        lonallp = lonfdot,\
                        toxallp = toxfdot,\
                        toyallp = toyfdot,\
                        tozallp = tozfdot\
                        )
                print "------"
                print ""



        if args.r:
            repair_disp_dudl_folder()

        if args.j:
            pot.create_job_tdi(ls=[0.0,1.0], \
                    execute="/u/aglen/start.4.6.28.tid.quick.sh",\
                    verbose = args.verbose \
                    )

        if args.jd:
            pot.create_job_dispdudl( \
                    verbose = args.verbose \
                    )


        # make moduel or simple list to get the one parametrization
    os.chdir(hier)
