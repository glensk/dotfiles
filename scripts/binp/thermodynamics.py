#!/usr/bin/env python
import numpy as np
from scipy.optimize import minimize_scalar
import sympy
import sys
import os
import glob
import warnings
import argparse
import shutil
import errno


import feos
import fqh
import fah
import fel
import utils
import hesse as h
reload(feos)    # necessary when working in ipython and changing code on the run
reload(fqh)     # necessary when working in ipython and changing code on the run
reload(fah)     # necessary when working in ipython and changing code on the run
reload(fel)     # necessary when working in ipython and changing code on the run
reload(utils)   # necessary when working in ipython and changing code on the run
reload(h)       # necessary when working in ipython and changing code on the run


def g_to_conc(filename = None,revert = False):
    kB = 0.0861734229648141309
    if filename == None:
        sys.exit('please provide a filename')
    g = np.loadtxt(filename)
    if revert == False:
        pass
    else:
        gtmp = np.copy(g)
        g[:,0] = gtmp[:,1]
        g[:,1] = gtmp[:,0]
    out = np.copy(g)
    out[:,1] = np.exp(-g[:,1]/(kB*g[:,0]))
    np.savetxt("concentration_from_"+str(filename),out,fmt="%.0f %.12f")
    return out

def conc_to_g(filename = None, revert = False):
    kB = 0.0861734229648141309
    if filename == None:
        sys.exit('please provide a filename')
    c = np.loadtxt(filename)
    if revert == False:
        pass
    else:
        ctmp = np.copy(c)
        c[:,0] = ctmp[:,1]
        c[:,1] = ctmp[:,0]

    if c[0,0] == 0:
        c[0,0] = 1
    out = np.copy(c)
    out[:,1] = -kB*c[:,0]*np.log(c[:,1])
    np.savetxt("gibbs_formation_from_"+str(filename),out,fmt="%.0f %.12f")
    return out

def _printred(var):
    ENDC = '\033[0m'
    red = '\033[31m'
    print red + str(var) + ENDC

def _printgreen(var):
    ENDC = '\033[0m'
    red = '\033[32m'
    print red + str(var) + ENDC

def write_inputdata(filename = None, text = None):
    if filename == None:
        print "No output written since no filename"
        return
    if text == None:
        text = ""
    if os.path.isfile != "True":
        with open(filename, "a") as myfile:
                myfile.write(text+"\n")
    return

def _unused_read_inputdata_not_finished(filename = None):
    if filename == None:
        print "No output written since no filename"
        return
    print "filename:",filename
    if os.path.isfile(filename) == True:
        # write 1 si ... part of sysargv in
        with open(filename) as f:
            content = f.readlines()
        content = [line.strip() for line in content]
        for i in content:
            print i.split()
        print "content",content
        print ""
        print ""

        # write EVinet and Fqh
        # write empty space
        # this is repeated for every part: al-si, bulk-refernce, si
        # copy corresponding files in there if they do not exist else... rename(1/2)?
        pass # write 1 si ... part of sysargv in


def copyanything(src, dst):
    try:
        shutil.copytree(src, dst)
    except OSError as exc: # python >2.5
        if exc.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else: raise

class surface:
    def __init_(self, surf = None):
        if surf == None:
            sys.exit("This claas needs an Free Energy surface as input")
        return



class f:
    def __init__(self, args = None):
        self.pressure = 0.000101325   # GPa
        self.pressure = 0.0001   # GPa
        self.outputfolder = None

        self.fixv = None
        self.fixv_file = args.fixv  # default = None

        if self.fixv_file != None:
            self.fixv = np.loadtxt(args.fixv)
            #print self.fixv


        self._verbose = False
        self.__verbose = False

        self.temperatures = np.arange(5.)  # more then 3 necessary for later evaluation
        self.contributions = []
        self._contributions = []
        self._v = False
        self._d = False
        self._eos = False
        self._fqh = False
        self._fel = False
        self._fah = False
        self._prefix = ""
        self._suffix = ""

        if args:
            self._verbose = args.v
            self._d = args.d
            self.pressure = args.P
            self._fqh = args.fqh
            self._fel = args.fel
            self._fah = args.fah
            self._eos = args.eos
            self._prefix = args.prefix
            self._suffix = args.suffix

        if self._verbose:
            print "args:",args

        ############################################
        # set pressure
        ############################################
        GPaTomeVAng3=6.2415097
        self._pressure = self.pressure * GPaTomeVAng3

        ###################################################################################
        # import feos (and read Evinet) (needs at least to exist)
        ###################################################################################
        self.feos = feos.eos()
        self.feos.import_parameters_data(filename = self._eos)
        if self._verbose:
            print "self.feos.inputfile_parameters:",self.feos.inputfile_parameters
            print "self.feos.parameters:",self.feos.parameters
        if self.feos.parameters == None:
            sys.exit("EVient or similar has to exist")
        v0 = self.feos.parameters[1]
        e0 = self.feos.parameters[0]
        self.contributions.append(self.feos._outfile)
        self._contributions.append(os.path.basename(self.feos._outfile))



        ###################################################################################
        # import self.fqh.surface (+ ensure temperatures from 0 to Tmax in 1K stpes)
        ###################################################################################
        self.fqh = fqh.qh()
        self.fqh._verbose = self._verbose
        if self._verbose:
            print "self._fqh:",self._fqh
        self.fqh.surface = self.fqh.import_surface(surface_filename = self._fqh, interpolate = True, sysexit = False)
        if type(self.fqh.surface) == np.ndarray:   # here if we just have one surface
            if self.fqh.surface_filename:
                self.contributions.append(self.fqh.surface_filename)
                self._contributions.append(os.path.basename(self.fqh.surface_filename))
                self.temperatures = self.fqh.surface[:,0]
                temperatures_chek = np.arange(self.temperatures.max()+1)
                if np.array_equal(self.temperatures,temperatures_chek) != True:
                    print "self.temperatures:",self.temperatures
                    print "temperatures_chek:",temperatures_chek
                    sys.exit("temperatures are not in 1K steps or internal error")
            if self._verbose:
                if type(self.fqh.surface) != np.ndarray:
                    _printred("Fqh surface NOT found!")

        ###################################################################################
        # import ah (only if fqh surface found) (+ ensures same temperatures in 1 K steps)
        ###################################################################################
        if type(self.fqh.surface) == np.ndarray:
            self.fah = fah.ah()
            self.fah.surface = self.fah.import_surface(surface_filename = self._fah, interpolate = True, sysexit = False)
            if self.fah.surface_filename:
                self.contributions.append(self.fah.surface_filename)
                self._contributions.append(os.path.basename(self.fah.surface_filename))

            # make the surfaces have equal length
            if type(self.fah.surface) == np.ndarray:
                length = np.array([len(self.fqh.surface[:,0]), len(self.fah.surface[:,0])]).min()
                self.fqh.surface = self.fqh.surface[:length]
                self.fah.surface = self.fah.surface[:length]
                self.temperatures = self.temperatures[:length]
                temperatures_chek = np.arange(self.fah.surface[:,0].max()+1)
                if (self.temperatures == temperatures_chek).all() != True:
                    print "self.temperatures2:",self.temperatures
                    print "temperatures_chek2:",temperatures_chek
                    sys.exit("temperatures are not in 1K steps or internal error")

            if self._verbose:
                if type(self.fah.surface) != np.ndarray:
                    _printred("Fah surface NOT found!")




        ###################################################################################
        # import self.fel.surface (+ ensure temperatures from 0 to Tmax in 1K stpes)
        ###################################################################################
        if type(self.fqh.surface) == np.ndarray:
            self.fel = fel.el()
            print "fel in:",self._fel
            self.fel.surface = self.fel.import_surface(surface_filename = self._fel, interpolate = True, sysexit = False)
            print "fel:",self.fel.surface
            print "fel:",self.fel.surface_filename

            if self.fel.surface_filename:
                print ":::",fel
                print "has surface fel>>>>"
                print self.fel.surface
                print "has surface fel <<<<<<"

                if self.fel.surface_filename:
                    print "in2"
                    self.contributions.append(self.fel.surface_filename)
                    self._contributions.append(os.path.basename(self.fel.surface_filename))

                # make the surfaces have equal length
                if type(self.fel.surface) == np.ndarray:
                    print "in3"
                    length = np.array([len(self.fqh.surface[:,0]), len(self.fel.surface[:,0])]).min()
                    self.fqh.surface = self.fqh.surface[:length]
                    self.fel.surface = self.fel.surface[:length]
                    self.temperatures = self.temperatures[:length]
                    temperatures_chek = np.arange(self.fel.surface[:,0].max()+1)
                    if (self.temperatures == temperatures_chek).all() != True:
                        print "self.temperatures3:",self.temperatures
                        print "temperatures_chek3:",temperatures_chek
                        sys.exit("3:temperatures are not in 1K steps or internal error")

                if self._verbose:
                    if type(self.fel.surface) != np.ndarray:
                        _printred("Fel surface NOT found!")



        ############################################
        # printout verbose help
        ############################################
        if self._verbose:
            print "saveqh:",self.fqh.surface_filename
            print "saveah:",self.fah.surface_filename
            print "contributions:    ",self.contributions

        #print  utils.get_commom_part_of_filenames(liste = self.contributions)
        _printgreen("--------------------------------------------------------------------")
        print "contributions:    ",self._contributions
        print "pressure:         ",self.pressure,"GPa"
        print "temperature range:",int(self.temperatures.min()),"...",int(self.temperatures.max()),"Kelvin"
        _printgreen("--------------------------------------------------------------------")

        ############################################
        # outputfolder .qh.ah
        ############################################
        add_qh = ""
        if self.fqh.surface_filename != None:
            add_qh = ".qh"
        add_el = ""
        if self.fel.surface_filename != None:
            add_el = ".el"
        add_ah = ""
        if self.fah.surface_filename != None:
            add_ah = ".ah"

        ############################################
        # outputfolder prefix suffix
        ############################################
        #add_common = utils.common_prefix(self.contributions, common_suffix = True)
        prefix = self._prefix
        if prefix != "":
            prefix = self._prefix+"_"

        suffix = self._suffix
        if suffix != "":
            suffix = "_"+self._suffix
        ############################################
        # outputfolder fqh alats
        ############################################
        fqhname = os.path.basename(self.fqh.surface_filename)
        #print fqhname
        #print fqhname.split("d_order__")
        #print fqhname.split("d_order__")[1]
        #alats = fqhname.split("d_order__")[1].split("_e")[0].split("_v")[0]
        #print "fqhname:",fqhname
        alats = fqhname.split("d_order__")
        #print "len:",len(alats),type(alats)
        if len(alats) > 1:
            alats = fqhname.split("d_order__")[1]
        else:
            alats = ""

        #print "a:",alats
        if alats == fqhname or alats == "":
            alats = ""
        else:
            alats = "_"+alats

        ############################################
        # outputfolder
        ############################################
        if self._verbose:
            print "prefix:",prefix
            print "add_qh:",add_qh
            print "add_el:",add_el
            print "add_ah:",add_ah
            print "fqname:",fqhname
            print "alats :",alats
        self.outputfolder = prefix+"output_"+str(self.pressure)+"GPa_eos"+add_qh+add_el+add_ah+alats+suffix+"/"


        if os.path.isdir(self.outputfolder) == True:
            if self._d == True:
                pass
            else:
                sys.exit("outputfolder \""+str(self.outputfolder)+"\" exists; use -d option to delete")

        if self._verbose:
            print "write add_qh       :",add_qh
            print "write add_ah       :",add_ah
            print "write self.pressure:", self.pressure
            print "self.outputfolder       :",self.outputfolder


        ##################################################################################
        ##################################################################################
        ##################################################################################
        # sum up the surface
        ##################################################################################
        ##################################################################################
        ##################################################################################
        self.volume_expansion = np.empty((len(self.temperatures),2))
        self.gibbs_energy = np.empty((len(self.temperatures),2))

        #V = sympy.Symbol('V')
        warnings.filterwarnings('error')

        if self._verbose:
            import time
            t1 = time.time()
            _printred("start 1 ... expansion ...")

        ##################################################################################
        # get volume for 0 pressure
        ##################################################################################

        ###########################################################################
        # define surface (eos, eos+qh, eos+qh+ah)
        ###########################################################################
        #print self.surf_qh(0)


        #def surf(V):
        #    return self.surf_eos(V)
        def surf_eos(V):  # not this! we need to add pressure
            return self.surf_eos(V)

        #if type(self.fqh.surface) == np.ndarray:
        #    def surf(V):
        #        return self.surf_eos(V) + self.surf_qh(ind,V)
        #    def surf_eos_qh(V):
        #        return self.surf_eos(V) + self.surf_qh(ind,V)
        #if type(self.fqh.surface) == np.ndarray:
        #    def surf(V):
        #        return surf(V) + surf_qh(ind,V)


        #if type(self.fel.surface) == np.ndarray:
        #    def surf(V):
        #        return self.surf(ind,V) + self.surf_el(ind,V)

        #if type(self.fah.surface) == np.ndarray:
        #    def surf(V):
        #        return surf(ind,V) + surf_ah(ind,V)
        #
        fqh_is = False
        fel_is = False
        fah_is = False
        if type(self.fqh.surface) == np.ndarray: fqh_is = True
        if type(self.fel.surface) == np.ndarray: fel_is = True
        if type(self.fah.surface) == np.ndarray: fah_is = True
        print "fqh_is:",fqh_is
        print "fel_is:",fel_is
        print "fah_is:",fah_is

        if fqh_is == True and fel_is == False and fah_is == False:
            def surf(V):
                return self.surf_eos(V) + self.surf_qh(ind,V)
            def surf_eos_qh(V):
                return self.surf_eos(V) + self.surf_qh(ind,V)
            print "------------------- feos -----------------------"
            print "------------------- fqh -----------------------"


        if fqh_is == True and fel_is == True and fah_is == False:
            def surf(V):
                return self.surf_eos(V) + self.surf_qh(ind,V) + self.surf_el(ind,V)
            print "------------------- feos -----------------------"
            print "------------------- fqh -----------------------"
            print "------------------- fel -----------------------"

        if fqh_is == True and fel_is == False and fah_is == True:
            def surf(V):
                return self.surf_eos(V) + self.surf_qh(ind,V) + self.surf_ah(ind,V)
            print "------------------- feos -----------------------"
            print "------------------- fqh -----------------------"
            print "------------------- fah -----------------------"

        if fqh_is == True and fel_is == True and fah_is == True:
            def surf(V):
                return self.surf_eos(V) + self.surf_qh(ind,V) + self.surf_el(ind,V) + self.surf_ah(ind,V)
            print "------------------- feos -----------------------"
            print "------------------- fqh -----------------------"
            print "------------------- fel -----------------------"
            print "------------------- fah -----------------------"

        if self._verbose:
            print self.temperatures[-1]
            if self.fqh.temperatures != None:
                print self.fqh.temperatures[-1]


        ###########################################################################
        # loop over temperature if no external volume given
        ###########################################################################
        if self.fixv == None:
            print "getting gibbs energy at corresponding pressure of "+str(self.pressure)+"GPa"
            v0 = False
            for ind, temp in enumerate(self.temperatures):  # not self.fqh.temperatures
                #print "surf:",surf,type(surf)
                #vol_all = minimize_scalar(surf)  # [ind] is changed in surf_eos...
                try:
                    vol_all = minimize_scalar(surf)  # [ind] is changed in surf_eos...
                    if self._verbose == 2:
                        print "ind,",ind,"temp:",temp,"vol_all:",vol_all
                    vol = vol_all.x  # [ind] is changed in surf_eos...
                    if self._verbose == 2:
                        print "ind,",ind,"temp:",temp,"vol_all:",vol_all
                except RuntimeWarning:
                    _printred("exc caught at "+str(temp)+"K")
                    print "ind,",ind,"temp:",temp,"vol_all:",vol_all
                    self.volume_expansion = self.volume_expansion[:ind]
                    self.gibbs_energy = self.gibbs_energy[:ind]
                    break

                if v0 == False: # in vorbereitung auf einen 3D plot
                    v0 = vol
                    volrange = np.arange(vol-v0/30.,vol+v0/30.,.1)
                    enerange = [surf(volr) for volr in volrange]

                    _printred("writing: e_vs_v_"+str(int(temp))+"K.dat")
                    #print np.transpose([np.array(volrange),np.array(enerange)])
                    enerange_all = np.empty([len(self.temperatures),len(volrange)])
                    enerange_all[:] = np.NAN
                    np.savetxt("e_vs_v_"+str(int(temp))+"K.dat",np.transpose([np.array(volrange),np.array(enerange)]))

                    #y_temp_empty = temps_save
                    #z_ener_empty = np.copy(x_vols_empty)
                    #empty_array = np.empty([temps_save.size,np.arange(vol-v0/10.,vol+v0/10.,0.1).size,np.arange(vol-v0/10.,vol+v0/10.,0.1).size])
                    #empty_array[:] = np.NAN
                enerange = [surf(volr) for volr in volrange]
                enerange_all[ind,:] = enerange

                self.volume_expansion[ind] = [temp, vol]
                self.gibbs_energy[ind] = [temp, surf(vol)]


        if self.fixv != None:
            print "getting gibbs energy at corresponding volume "
            print self.fixv
            print "--"
            for ind, temp in enumerate(self.temperatures):  # not self.fqh.temperatures
                vol = None
                #print "::",temp
                if temp in self.fixv[:,0]:
                    tempidx = np.where(self.fixv[:,0]==temp)[0][0]
                    vol = self.fixv[tempidx,1]
                    #print "::",temp,tempidx,vol

                if vol != None:
                    self.volume_expansion[ind] = [temp, vol]
                    self.gibbs_energy[ind] = [temp, surf(vol)]
                else:
                    self.volume_expansion[ind] = [temp, np.NAN]
                    self.gibbs_energy[ind] = [temp, np.NAN]

        print "DONE"

        # save e_vs_v_all.dat
        #np.savetxt("e_vs_v_enes_"+str(self.pressure)+"GPa.dat",enerange_all)
        #np.savetxt("e_vs_v_vols_"+str(self.pressure)+"GPa.dat",volrange)
        #np.savetxt("e_vs_v_temp_"+str(self.pressure)+"GPa.dat",self.temperatures)

        ###########################################################
        # get gibbs energy at T=0K (only 0-point vibrations)
        ###########################################################
        # we have to use the correct volume including 0-point vibrations
        self.feos.parameters[0] = 0
        ind=0
        ene0 = surf(self.volume_expansion[0,1])
        ene0shifted = self.gibbs_energy[0,1]
        #print "ene0:",ene0,"ene0shifted:",ene0shifted
        self.gibbs_energy_vib = np.copy(self.gibbs_energy)
        self.gibbs_energy_vib[:,1] = self.gibbs_energy[:,1]-self.gibbs_energy[0,1]+ene0

        # reset parameters correctly
        self.feos.parameters[0] = e0




        if self._verbose:
            t2 = time.time()
            _printred("stop  1 ... expansion ... "+str(t2-t1))


        ##################################################################################
        # in case we want to get it quicker we can simply do everty 2nd or third
        # temperature and afterwards spline fit
        ##################################################################################
        #import hesse
        #self.gibbs_energy_fit = hesse.fqh_interpolate(self.gibbs_energy[::2])

        ##################################################################################
        # fit/interpolate stuff ONLY IF NECESSARY!
        ##################################################################################
        # gibbs_energy_fit
        print "interpolating"
        check_temprange = np.unique(np.diff(self.gibbs_energy[:,0]))
        print check_temprange
        self.gibbs_energy_fit = self.gibbs_energy
        if len(check_temprange) != 1:
            if check_temprange[0] != 1.:
                self.gibbs_energy_fit = h.fqh_interpolate(self.gibbs_energy)
        self.gibbs_energy_fit = self.gibbs_energy_fit[~np.isnan(self.gibbs_energy_fit).any(axis=1)]
        self.volume_expansion = self.volume_expansion[~np.isnan(self.volume_expansion).any(axis=1)]
        print "fitted"
        self.entropy = np.copy(self.gibbs_energy_fit)
        self.entropy[:,1] = np.gradient(self.gibbs_energy_fit[:,1])*-11.604506  # mev/K -> kB   (clip(0) makes every < 0 values to 0; this happens usually only close to 0
        self.entropy[:,1] = self.entropy[:,1].clip(0)
        self.entropy[-1,1] = self.entropy[-2][1]+self.entropy[-2][1]-self.entropy[-3][1]


        self.heat_capacity_isobaric = np.copy(self.entropy)
        self.heat_capacity_isobaric[:,1] = np.gradient(self.entropy[:,1])*self.gibbs_energy_fit[:,0]
        self.heat_capacity_isobaric[:,1] = self.heat_capacity_isobaric[:,1].clip(0)
        self.heat_capacity_isobaric[-2,1] = self.heat_capacity_isobaric[-3][1]+self.heat_capacity_isobaric[-3][1]-self.heat_capacity_isobaric[-4][1]
        self.heat_capacity_isobaric[-1,1] = self.heat_capacity_isobaric[-2][1]+self.heat_capacity_isobaric[-2][1]-self.heat_capacity_isobaric[-3][1]


        ##################################################################################
        # write stuff
        ##################################################################################
        self.write()

    def surf_eos(self, V):
        return feos.vinet(V, *self.feos.parameters)+self._pressure*V # here we want numpy

    def surf_qh(self, ind, V): # this is general for 2nd and 3rd order surfaces
        order1 = self.fqh.surface[ind][1] \
                +self.fqh.surface[ind][2]*V
        _qh_order = self.fqh.surface.shape[1] - 2
        if _qh_order == 1:
            return order1
        order2 = order1 \
                +self.fqh.surface[ind][3]*V*V
        if _qh_order == 2:
            return order2
        order3 = order2 \
                +self.fqh.surface[ind][4]*V*V*V
        if _qh_order == 3:
            return order3

    def surf_ah(self, ind, V): #this is general for 2nd and 3rd order surfaces
        order1 = self.fah.surface[ind][1] \
                +self.fah.surface[ind][2]*V
        _ah_order = self.fah.surface.shape[1] - 2
        if _ah_order == 1:
            return order1
        order2 = order1 \
                +self.fah.surface[ind][3]*V*V
        if _ah_order == 2:
            return order2
        order3 = order2 \
                +self.fah.surface[ind][4]*V*V*V
        if _ah_order == 3:
            return order3

    def surf_el(self, ind, V): #this is general for 2nd and 3rd order surfaces
        order1 = self.fel.surface[ind][1] \
                +self.fel.surface[ind][2]*V
        _el_order = self.fel.surface.shape[1] - 2
        if _el_order == 1:
            return order1
        order2 = order1 \
                +self.fel.surface[ind][3]*V*V
        if _el_order == 2:
            return order2
        order3 = order2 \
                +self.fel.surface[ind][4]*V*V*V
        if _el_order == 3:
            return order3

    ##################################################################################
    # write thermodynamic properties
    ##################################################################################
    def write(self):
        if os.path.isdir(self.outputfolder) == True:
            if self._d == True:
                shutil.rmtree(self.outputfolder)
            else:
                sys.exit("outputfolder \""+str(self.outputfolder)+"\" exists; use -d option to delete")
        if os.path.isdir(self.outputfolder) != True:
            os.makedirs(self.outputfolder)


        np.savetxt(self.outputfolder+"/volume_expansion",self.volume_expansion, fmt="%.0f %.12f")
        np.savetxt(self.outputfolder+"/gibbs_energy_orig",self.gibbs_energy, fmt="%.0f %.12f")
        np.savetxt(self.outputfolder+"/gibbs_energy_vib",self.gibbs_energy_vib, fmt="%.0f %.12f")
        np.savetxt(self.outputfolder+"/gibbs_energy",self.gibbs_energy_fit, fmt="%.0f %.12f")
        np.savetxt(self.outputfolder+"/entropy",self.entropy, fmt="%.0f %.12f")
        np.savetxt(self.outputfolder+"/heat_capacity_isobaric",self.heat_capacity_isobaric, fmt="%.0f %.12f")


        ##################################################################################
        # write inputdata
        ##################################################################################
        inputfolder=self.outputfolder+"/inputdata/"
        if os.path.isdir(inputfolder) != True:
            os.makedirs(inputfolder)

        filename = inputfolder+"folder"
        _printgreen("created:  "+self.outputfolder)

        f = self.feos.inputfile_parameters
        #print "F:",f
        if f:
            shutil.copyfile(f, inputfolder+os.path.basename(f))
            write_inputdata(filename, "EVinet: "+os.path.abspath(f))

        f  = self.fqh.surface_filename
        #print "F:",f
        if f:
            shutil.copyfile(f, inputfolder+os.path.basename(f))
            write_inputdata(filename, "Fqh: "+os.path.abspath(f))

        f  = self.fah.surface_filename
        #print "F:",f
        if f:
            shutil.copyfile(f, inputfolder+os.path.basename(f))
            write_inputdata(filename, "Fah: "+f)

            cf = glob.glob(os.path.dirname(f)+"/fit.input*")
            print "############# cf",cf
            if len(cf) == 1:
                print "############# cf",cf
                cf = cf[0]
                print "cf:",cf
                print "to:",inputfolder+os.path.basename(cf)+"ah"
                shutil.copyfile(cf, inputfolder+os.path.basename(cf)+"ah")

            cf2 = glob.glob(os.path.dirname(f)+"/Fah_surface")
            print "############# cf2",cf2
            if len(cf2) == 1:
                print "############# cf2"
                cf2 = cf2[0]
                print "cf:",cf
                print "cf2:",cf2
                print "to:",inputfolder+os.path.basename(cf2)+"ah"
                shutil.copyfile(cf2, inputfolder+os.path.basename(cf2)+"ah")

        write_inputdata(filename, "")
        write_inputdata(filename, "started: "+" ".join(sys.argv))

class fformmurn:
    def __init__(self, args = None):
        '''
        args == [ "eos_qh_d_31_output_0.0001GPa",  "-31/32" "eos_qh_b_32_output_0.0001GPa" ]
        '''
        print args


class fform:
    def __init__(self, args = None):
        '''
        args == [ "eos_qh_d_31_output_0.0001GPa",  "-31/32" "eos_qh_b_32_output_0.0001GPa" ]
        '''
        self._prefix = ""

        if not args:
            sys.exit("you need to specify the inputfolder")
        self._d = args.d
        self._prefix = args.prefix
        self._suffix = args.suffix

        if self._prefix == None:
            self._prefix = ""
        if self._prefix == False:
            self._prefix = ""
        if self._prefix != "":
            self._prefix = self._prefix + "_"
        if self._suffix == None:
            self._suffix = ""
        if self._suffix == False:
            self._suffix = ""
        if self._suffix != "":
            self._suffix = "_" + self._suffix



        self.args = args.fform
        self._verbose = args.v

        if len(self.args) == 1:
            self.args = self.args[0].split()
            print "self.args:",self.args



        ##############################################
        # path and n && checks
        ##############################################
        self.path_string = self.args[1::2]
        self.n_string = self.args[0::2]
        self.path = []
        self.f_in = []
        self.f_int = []
        self.n = []
        if len(self.path_string) != len(self.n_string):
            sys.exit("self.n and self.path and f have different length, do: NUMBER path NUMBER path ...")
        nsp=utils.NumericStringParser()
        for i in self.n_string:
            self.n.append(utils.string_to_mathematical_expression(i))

        ##########################################
        # self.outputfolder
        ##########################################
        self.path_string_basename = self.path_string[:]
        for ind,string in enumerate(self.path_string):
            if self.path_string[ind][-1] == "/":
                self.path_string[ind] = self.path_string[ind][:-1]
            self.path_string_basename[ind] = os.path.basename(self.path_string[ind])
        #self.path_string[ind] = self.path_string[ind].replace("/$","")
        if self._verbose:
            _printgreen("----------------------------------")
            print "path_string  :",self.path_string
            print "path_stringbn:",self.path_string_basename
            print "len:",len(self.path_string),type(self.path_string)
            print "len:",len(self.path_string_basename),type(self.path_string_basename)
            _printgreen("----------------------------------")
        self.path_prefix = utils.common_prefix(self.path_string_basename)
        self.path_suffix = utils.common_prefix(self.path_string_basename, common_suffix = True)
        if self.path_prefix != "":
            self.path_prefix = "_"+self.path_prefix
        if self.path_suffix != "":
            self.path_suffix = "_"+self.path_suffix


        if self._verbose:
            print "path_string:",self.path_string
            print "prefix     :",self.path_prefix
            print "suffix     :",self.path_suffix
            print "self.prefix:",self._prefix
        self.outputfolder = "formation"+self._prefix+self.path_prefix+self.path_suffix+self._suffix

        if os.path.isdir(self.outputfolder) == True:
            if self._d == True:
                pass
            else:
                sys.exit("outputfolder \""+str(self.outputfolder)+"\" exists; use -d option to delete")



        ##########################################
        # loop through inputfolder
        ##########################################
        print ""
        inputdata_folder = ""
        for ind,i in enumerate(self.path_string):
            #print "i:",i
            fin = i+"/gibbs_energy"
            if os.path.isfile(fin) != True:
                fin = i+"/Gibbs_energy"
                if os.path.isfile(fin) != True:
                    sys.exit(fin+" does not exist ...")
            self.path.append(fin)
            _printgreen("importing "+fin)
            farr = np.loadtxt(fin)
            farr_int = h.fqh_interpolate(farr)
            self.f_in.append(farr)
            self.f_int.append(farr_int)
            ##############################################################################
            # copy inputfolder to formation folder
            ##############################################################################
            fif = i+"/inputdata/folder"
            inputdata_folder = inputdata_folder + str(self.n_string[ind])+" "+ str(os.path.abspath(self.path_string[ind]))+"\n"
            if os.path.isfile(fif) == True:
                #_printred("--------------------------------")
                # write 1 si ... part of sysargv in
                with open(fif) as f:
                    content = f.readlines()
                for s in content:
                    if 'started: ' not in s:
                        inputdata_folder = inputdata_folder + s
            else:
                inputdata_folder = inputdata_folder + os.path.abspath(fin)+ '\n'
            #print ">> inp:",inputdata_folder
            #print ""

        if self._verbose:
            print "self.path",self.path
            print "self.n:",self.n


        print "getting temperature intersection ..."
        tempslist = []
        for ind, i in enumerate(self.path):
            tempslist.append(self.f_int[ind][:,0])
        self.temps = utils.schnittmenge_several_lists(tempslist, verbose = False)


        print "getting equally lenghy fform.data and fform.f (only intersecting temperatures)"
        self.data = np.zeros((len(self.temps),len(self.path_string)+1))
        self.f = np.zeros((len(self.temps),len(self.path_string)))
        self.fn = np.zeros((len(self.temps),len(self.path_string)))
        self.fform = np.zeros((len(self.temps),2))
        self.data[:,0] = self.temps
        self.fform[:,0] = self.temps

        for ind,i in enumerate(self.f_int):
            self.f[:,ind] = self.f_int[ind][:,1][:len(self.temps)]
            self.fn[:,ind] = self.f_int[ind][:,1][:len(self.temps)]*self.n[ind]
            self.data[:,ind+1] = self.f_int[ind][:,1][:len(self.temps)]
            if ind == 0:
                self.fform[:,1] = self.fn[:,ind]
            else:
                self.fform[:,1] = self.fn[:,ind] + self.fform[:,1]
            print ind+1,"/",len(self.f_in),"   | n = ",self.n[ind]
        self.plot = self.fform[:,0],self.fform[:,1]


        ###############################################################
        # missing stuff == conc
        ###############################################################
        warnings.filterwarnings("ignore")
        #fform = np.transpose([temps, fform])
        kB = 0.0861734229648141309
        self.concentration = np.copy(self.fform)
        self.concentration[:,1] = np.exp(-self.fform[:,1]/(kB*self.fform[:,0]))
        self.t_vs_concentration = np.copy(self.fform)
        self.t_vs_concentration[:,0] = self.concentration[:,1]
        self.t_vs_concentration[:,1] = self.concentration[:,0]
        self.g_vs_c = np.copy(self.fform)
        self.g_vs_c[:,0] = self.concentration[:,1]


        ###############################################################
        # missing stuff == entropy
        ###############################################################
        self.entropy = np.copy(self.fform)
        self.entropy[:,1] = np.gradient(self.fform[:,1])*-11.604506  # mev/K -> kB   (clip(0) makes every < 0 values to 0; this happens usually only close to 0
        self.entropy[:,1] = self.entropy[:,1].clip(0)
        self.entropy[-1,1] = self.entropy[-2][1]+self.entropy[-2][1]-self.entropy[-3][1]

        ###################################################################################
        # write stuff
        ###################################################################################
        if os.path.isdir(self.outputfolder) == True:
            if self._d == True:
                shutil.rmtree(self.outputfolder)
            else:
                sys.exit("outputfolder \""+str(self.outputfolder)+"\" exists; use -d option to delete")
        if os.path.isdir(self.outputfolder) != True:
            os.makedirs(self.outputfolder)

        np.savetxt(self.outputfolder+"/gibbs_formation",self.fform, fmt="%.0f %.12f")
        np.savetxt(self.outputfolder+"/concentration",self.concentration, fmt="%.0f %.12f")
        np.savetxt(self.outputfolder+"/t_vs_concentration",self.t_vs_concentration, fmt="%.12f %.0f")
        np.savetxt(self.outputfolder+"/g_vs_c",self.g_vs_c, fmt="%.12f %.0f")
        np.savetxt(self.outputfolder+"/entropy",self.entropy, fmt="%.0f %.12f")

        ################################################################
        # write inputdata
        ################################################################
        inputfolder=self.outputfolder+"/inputdata/"
        if os.path.isdir(inputfolder) != True:
            os.makedirs(inputfolder)
        _printgreen("created "+self.outputfolder+"/")


        ################################################################
        # write inputdata folder
        ################################################################
        filename = inputfolder+"folder"

        sysargv = sys.argv[0]+" -fform \""+" ".join(self.args)+"\""
        if self.path_suffix != "":
            sysargv = sysargv + " -s "+args.suffix
        if self.path_prefix != "":
            sysargv = sysargv + " -p "+args.prefix
        write_inputdata(filename, "")
        #write_inputdata(filename, " ".join(sys.argv).replace("fform ","fform \""))
        write_inputdata(filename, inputdata_folder)
        write_inputdata(filename, "")
        write_inputdata(filename, "started: "+sysargv)

        ################################################################
        # write inputdata {0_,1_, ...} save sources of fformation
        ################################################################
        for ind,i in enumerate(self.path_string):
            #print "path_str:",self.path_string[ind]
            if os.path.isdir(self.path_string[ind]+'/inputdata'):
                copyanything(self.path_string[ind]+'/inputdata',inputfolder+str(ind)+"_"+str(os.path.basename(self.path_string[ind])))
            elif os.path.isfile(self.path_string[ind]+'/gibbs_energy'):
                shutil.copyfile(self.path_string[ind]+'/gibbs_energy',inputfolder+str(ind)+"_gibbs_energy")
            elif os.path.isfile(self.path_string[ind]+'/Gibbs_energy'):
                shutil.copyfile(self.path_string[ind]+'/Gibbs_energy',inputfolder+str(ind)+"_Gibbs_energy")

class changehesse:
    def __init__(self, args = None):
        if not args:
            sys.exit("you need to specify the inputfolder")

        self.old = args.old
        self.new = args.new
        self.fah_data = args.fah_data
        self._verbose = args.verbose



if __name__ == '__main__':
    def help(p = None):
        string = """
        possible inputfiles (which are also recognized automatically):
        EVinet(_...)
        Fqh(_...)
        Fah(_...)
        Fel(_...) yet to implement
        getThermodynamics help text"""

        if p == None:
            from argparse import ArgumentDefaultsHelpFormatter
            p = argparse.ArgumentParser(description=string,
                    formatter_class=ArgumentDefaultsHelpFormatter)
        p.add_argument('-P',
            help='define pressure in GPa',
            type=float, default=0.0001)
        p.add_argument('-d',
            help='delete outputfolder if it exists',action='store_true', default=False)
        p.add_argument('-eos',type=str,
                help='set eos parameters filename (e.g. -eos EVinet_b_31)', default=None)
        p.add_argument('-fqh',type=str,
                help='set fqh surface filename (e.g. -fqh Fqh_b_32)', default=None)
        p.add_argument('-fel',type=str,
                help='set fel surface filename (e.g. -fel Fel_b_32)', default=None)
        p.add_argument('-fah',type=str,
                help='set fah surface filename (e.g. -fah Fah_b_32)', default=None)
        p.add_argument('-p', '--prefix',type=str,
                help='prefix of the created outputfolder', default="")
        p.add_argument('-s', '--suffix',type=str,
                help='suffix of the created outputfolder (currently only with fform)', default="")
        p.add_argument('-fform', nargs='+',
                help='make formation energy of outptfolder; e.g. -fform \"1 eos_qh_d_31_output_0.0001GPa -31/32 eos_qh_b_32_output_0.0001GPa\"', default=None, type=str)
        p.add_argument('-fformmurn', nargs='+',
                help='make formation energy of outptfolder; e.g. -fform \"1 eos_qh_d_31_output_0.0001GPa -31/32 eos_qh_b_32_output_0.0001GPa\"', default=None, type=str)
        p.add_argument('-old',type=str,
                help='for changehesse: set path to old Fqh_Helmholz files', default=None)
        p.add_argument('-new',type=str,
                help='for changehesse: set path to new Fqh_Helmholz files', default=None)
        p.add_argument('-fixv',type=str,
                help='path to external volume expansion where surface is fittet to', default=None)
        p.add_argument('-cg',type=str,
                help='concentration(2nd row) vs temperature(1st row) to gibbs_formation(2nd row) vs temperature(1st row)', default=None)
        p.add_argument('-cg_',type=str,
                help='concentration(1st row) vs temperature(2st row) to gibbs_formation(2nd row) vs temperature(1st row)', default=None)

        p.add_argument('-gc',type=str,
                help='gibbs_formation(2nd row) vs temperature(1st row) to concentration(2nd row) vs temperature(1st row)', default=None)
        p.add_argument('-gc_',type=str,
                help='gibbs_formation(1nd row) vs temperature(2st row) to concentration(2nd row) vs temperature(1st row)', default=None)
        #p.add_argument('-v',
        #    help='write volume_expansion_[pressure]GPa',
        #    action='store_true', default=False)
        p.add_argument('-v', action='count',
                help='verbose')
        return p
    #sys.exit()
    p = help()  # this gives the possibility to change some __init__ settings
    args = p.parse_args()

    if args.cg != None:
        conc_to_g(args.cg)
        sys.exit()

    if args.cg_ != None:
        conc_to_g(args.cg_, revert = True)
        sys.exit()

    if args.gc != None:
        g_to_conc(args.gc)
        sys.exit()

    if args.gc_ != None:
        g_to_conc(args.gc_, revert = True)
        sys.exit()

    if args.fform != None:
        fform = fform(args)
    else:
        f = f(args)


