#!/usr/bin/env python

import argparse
import os
import sys
import glob
import utils
import shutil
import numpy as np
import pylab
from scipy.integrate import quad
import hesse as h

np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=6)    # print only 6 digist after .

def _printred(var):
    ENDC = '\033[0m'
    red = '\033[31m'
    print red + str(var) + ENDC

def _printgreen(var):
    ENDC = '\033[0m'
    red = '\033[32m'
    print red + str(var) + ENDC

def find_surface_filename(sysexit = True):
    ''' '''
    possible_filename = [ "Fah", "Fah_b_*", "Fah_d_*", "Fah_*" ]
    for i in possible_filename:
        filecheck = glob.glob(i)
        if filecheck:
            if len(filecheck) == 1:
                surface_filename = filecheck[0]
                return surface_filename

    # if we came to this point we have to exit
    if sysexit:
        sys.exit("Did not found any Fah surface file in current folder")
    else:
        #print("Did not found any Fah surface file in current folder")
        return None

def get_dudlmeanfit(lambdas_array, energies_array, return_fit = False, return_parameters = False):
    ''' Fits tangence function to dudl vs lambda data.
    input:
        v: the volume you want to fit
        t: the temperature you want to fit
    How can we in general give n parameters to any function without having to type parameters[0], parameters[1], parameters[2], .....?
    '''

    def tanFunction(x, a0, a1, a2, a3):
        return -a0*np.tan(np.pi*((1-a1)*x+a2+0.5))+a3
    def tanFunctionguess(energies_array):
        return [ energies_array[0],0.5,0.3,-0.6]

    #lambdas_array = self.l
    #energies_array = self.dudlmean[v, t, :]

    if np.isnan(lambdas_array).any() == True or np.isnan(lambdas_array).any() == True:
        a = np.array([lambdas_array, energies_array])  # this still has nans
        print "a||:---",a
        if np.isnan(a).any() == True:
            print "WHY:"
            lambdas_array, energies_array = np.ma.compress_rows(np.ma.fix_invalid(a.T)).T  #no nas
    #print "--> compressed:",lambdas_array, energies_array
    if len(energies_array) is 0: return np.nan

    from scipy.optimize import curve_fit
    parameters, evinet_covariance = \
        curve_fit(tanFunction, lambdas_array, energies_array, tanFunctionguess(energies_array), maxfev=30000)
    #jprint ";;;energies_r:",energies_array
    #print ";;;parameters:",parameters
    #print ";;;"

    if parameters[1] < 0:
        print energies_array
        _printred("parameter 2 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 2 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 2 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 2 is sammer 0 we might have a probleme in fitting")
    if parameters[2] < 0:
        print energies_array
        _printred("parameter 3 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 3 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 3 is sammer 0 we might have a probleme in fitting")
        _printred("parameter 3 is sammer 0 we might have a probleme in fitting")

    # fit
    points = 101
    fit_x_min = 0
    fit_x_max = 1

    fit_x_values = np.linspace(fit_x_min, fit_x_max, points)
    fit_y_values = tanFunction(fit_x_values, parameters[0],parameters[1],parameters[2],parameters[3])

    # delta to fit
    fitdeltas = tanFunction(lambdas_array, parameters[0],parameters[1],parameters[2],parameters[3]) - energies_array

    # delta to fit max
    fitdeltasmax = abs(fitdeltas).max()

    fah = quad(tanFunction, 0, 1, args=(parameters[0],parameters[1],parameters[2],parameters[3]))[0]
    if return_fit == True:
        return fah, [fit_x_values,fit_y_values]
    if return_parameters == True:
        return fah, [fit_x_values,fit_y_values], parameters
    return fah #,fitdeltasmax

def write_Fah_surface(a, temps, dudlmeanfit, filename = 'Fah_surface'):
    print "a    :",a
    print "temps:",temps
    anz = np.sum(np.isnan(dudlmeanfit))
    out  = np.empty([a.size * temps.size - anz, 3])

    out[:]      = np.NAN
    i = 0
    for ia,a in enumerate(a):
        for it,t in enumerate(temps):
            #print ">",t,a,dudlmeanfit[ia,it]
            if np.isnan(dudlmeanfit[ia,it]) == True:
                continue
            out[i][0] = t
            out[i][1] = a
            out[i][2] = np.around(dudlmeanfit[ia,it], decimals=2)
            i = i+1
    print "filename:",filename
    print out
    np.savetxt(filename,out,fmt="%.1f %.7f %.3f")
    return


class ah():
    ''' anharmonic contribution to Gibbs free energy '''
    def __init__(self):
        self._verbose       = False
        self._filestring                        = None
        self._pwd                               = os.getcwd()
        self._dudl_cutoff_first_structures      = 50
        self._dudl_correlation_length           = 10

        self.files          = None
        self.surface        = None
        self.surface_expr   = None
        self.surface_filename = None
        self.a              = None      # array([ 4.05, 4.09, 4.13])
        self.t              = None      # array([ 250.,  500.,  700.,  934.])
        self.l              = None      # array([ 0.  ,  0.15,  0.5 ,  0.85,  1.  ])
        self.s              = None      #
        self.e              = None      #
        self.dudl           = None      # --> same as dUdL(_noJumps) file
        self.dudlmean       = None      # --> see below (mean dUdL for every lambda and temp)
        self.dudlmeanfit    = None      # --> see below same as avg_dUdL_fre
        self.dudlste        = None      # mean dUdl ste for every lambda
        self.dudlstd        = None      # mean dUdl std for every lambda
        self.data           = None      #
        self.ah             = None      # anharmonic corrections at certain temp and vol
        self.ahste          = None      # anharmonic ste         at certain temp and vol
        self.ahstd          = None      # anharmonic std         at certain temp and vol

        # [158]f.dudlmean (f.dudlstd, f.dudlste)
        # Out[158]:
        # array([[[  0.51,   0.27,   0.01,  -0.53,  -0.64],     <- 250  <- 4.05
        #         [  2.88,   1.66,  -0.16,  -1.57,  -2.45],     <- 500  <- 4.05
        #         [  6.09,   3.81,   0.12,  -2.74,  -4.43],     <- 700  <- 4.05
        #         [ 10.79,   6.62,   0.58,  -5.46,  -8.57]],    <- 934  <- 4.05
        #
        #        [[  0.73,   0.6 ,  -0.09,  -0.5 ,  -0.92],     <- 250  <- 4.05
        #         [  3.76,   2.34,   0.02,  -1.92,  -2.71],     <- 500  <- 4.05
        #         [  7.06,   4.63,   0.62,  -3.43,  -5.23],     <- 700  <- 4.05
        #         [ 13.32,   8.45,   1.25,  -5.57, -10.17]]     <- 934  <- 4.05
        #
        #        [[  1.28,   0.61,   0.04,  -0.57,  -0.83],     <- 250  <- 4.05
        #         [  4.88,   3.73,   0.36,  -2.06,  -3.07],     <- 500  <- 4.05
        #         [  9.27,   5.84,   1.09,  -3.36,  -5.84],     <- 700  <- 4.05
        #         [ 17.44,  11.05,   1.96,  -5.83, -11.71]]])   <- 934  <- 4.05

        # [163]f.dudlmeanfit
        # Out[163]:
        # array([[-0.064615, -0.02395 ,  0.390583,  0.621287],          <-4.05
        #        [-0.019864,  0.158969,  0.634918,  1.375083],          <-4.05
        #        [ 0.04227 ,  0.647443,  1.216228,  2.375211]])         <-4.05


        return

    def initialize_files_data(self):
        """ defines:    self.files
                        self.a
                        self.l
                        self.s
                        self.data filled with NaNs """
        print "initialize_files_data ..."
        self.files = utils.lsn(self._filestring)
        self.filesinfo = []

        atls = np.array([ utils.string_to_num_list(string) for string in self.files ])
        if self._verbose:
            print "atls:",atls
            print "atls.shape:",atls.shape
            print "atls.shape[1]:",atls.shape[1]

        # wir setzen a(alts) und t(temperatures) voraus, alles andere optional
        # oder kann spaeter hinzugefuegt werden
        self.a      = np.unique(atls[:,0])
        self.t      = np.unique(atls[:,1])
        if atls.shape[1] <= 2:
            return
        if atls.shape[1] > 2:
            self.l      = np.unique(atls[:,2])

        # determine maximun amount of seeds (seeds are not distinguishable)
        nr_s = 0
        for a in self.a:
            array_a = atls[np.nonzero(atls[:,0] == a)[0]]
            #print "a:",array_a
            for t in self.t:
                array_t = array_a[np.nonzero(array_a[:,1] == t)[0]]
                if len(array_t) is 0: continue # in case a temp/vol missing
                #print "t:",array_t,"len_t:",len(array_t)
                for l in self.l:
                    array_l = array_t[np.nonzero(array_t[:,2] == l)[0]]
                    if len(array_l) > nr_s: nr_s= len(array_l)

        self.s      = np.arange(nr_s)
        self.data   = np.empty((
            len(self.a),
            len(self.t),
            len(self.l),
            len(self.s)), dtype=object)  # necessary for arrays of unequal length
        self.data[:] = np.NAN
        return

    def import_dudl_files_data(self):
        """ Imports all files to self.data; Working with this data (e.g.
        extracting dudls) should be separated; Doing so this can be used for
        sphinx and vasp """
        def get_inds(a,t,l):
            for s in self.s:
                #print "a,t,l,s:",a,t,l,s
                if np.isnan(np.any(self.data[a,t,l,s])):
                    return s

        for file in self.files:
            dudl = pylab.loadtxt(file)
            atls = utils.string_to_num_list(file)
            a = atls[0]
            t = atls[1]
            l = atls[2]
            s = atls[3]

            inda = np.nonzero(self.a == a)[0][0]  # index
            indt = np.nonzero(self.t == t)[0][0]  # index
            indl = np.nonzero(self.l == l)[0][0]  # index
            inds = get_inds(inda,indt,indl)
            print "file:",file,a,t,l,s,"(ind:",inda,indt,indl,inds,")"
            self.filesinfo.append([file,"ind:",[inda,indt,indl,inds]])
            self.data[inda, indt, indl, inds] = dudl
        return

    def import_avg_dudl_data(self, avg_dudl_filenames = "*Ang_*K/avg_dUdL_lowplushigh_fre", verbose = False):
        '''
        can import:
         - "*Ang_*K/avg_dUdL_lowplushigh_fre"   (in high folder)
         - "*Ang_*K/avg_dUdL_fre"               (in low  folder) ????

        '''
        #################################################################################
        # get alat, temperatures, lambdas
        #################################################################################
        avg_dudl_files = utils.lsn(avg_dudl_filenames)
        self.atls = np.array([ utils.string_to_num_list(string) for string in avg_dudl_files ])
        self.a      = np.unique(self.atls[:,0])
        self.t      = np.unique(self.atls[:,1])

        self.l = np.array([])  # here we want to load all data! not exclude all lambda 0.5 if one lambda 0.5 is missing
        for i in avg_dudl_files:
            #print "i:",i
            avg_dudl = np.loadtxt(i)
            #print "avg_dudl:",avg_dudl
            #print "avg_dudl.shape:",avg_dudl.shape
            if len(avg_dudl.shape) == 1:
                continue
            #print "avg_dudl.shape:",len(avg_dudl.shape)
            l_tmp = avg_dudl[:,0]
            #print "l_tmp:",l_tmp
            for ll in l_tmp:
                self.l = np.append(self.l,ll)
            self.l      = np.unique(self.l)
            #print avg_dudl
            #print avg_dudl.shape
            #print "l:",l
        if verbose:
            print "a:",self.a
            print "t:",self.t
            print "l:",self.l
        #if len(self.l) == 9:
        #    print type(self.l)
        #    self.l = np.array([0.0, 0.15, 0.5, 0.85, 1.0])

        self.dudlmean     = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudllen      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlste      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlstd      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlmeanlen  = np.empty([self.a.size, self.t.size])
        self.dudlmeanfit  = np.empty([self.a.size, self.t.size])
        self.dudlmeanfit_tanfit  = np.empty([self.a.size, self.t.size,101,2])

        self.dudlmean[:]     = np.NAN
        self.dudllen[:]      = np.NAN
        self.dudlste[:]      = np.NAN
        self.dudlstd[:]      = np.NAN
        self.dudlmeanlen[:]  = np.NAN
        self.dudlmeanfit[:]  = np.NAN

        #self.dudlavg = np.copy(self.dudl)

        np.set_printoptions(suppress=True)   # display arrays withou 000000
        np.set_printoptions(precision=6)    # print only 6 digist after .

        #################################################################################
        # got trhount avg_dUdL files and fit
        #################################################################################
        for i in avg_dudl_files:
            if verbose:
                print "------------------>i:",i
            atls = utils.string_to_num_list(i)
            a = atls[0]
            t = atls[1]
            inda = np.nonzero(a == self.a)[0]
            indt = np.nonzero(t == self.t)[0]
            avg_dudl = np.loadtxt(i)
            #print alll
            #print avg_dudl,"|",len(avg_dudl),"|",avg_dudl.shape,len(avg_dudl.shape)
            if len(avg_dudl) == 0:
                continue # empty
            if len(avg_dudl.shape) == 1:
                avg_dudl = np.array([avg_dudl])

            for lline in avg_dudl:
                if verbose:
                    print "lline:",lline
                l = lline[0]
                indl = np.nonzero(l == self.l)[0]
                if verbose:
                    print "a,inda:",a,inda
                    print "t,indt:",t,indt
                    print "l,indl:",l,indl
                self.dudlmean[inda,indt,indl] = lline[1]
                self.dudlste[inda,indt,indl] = lline[2] #*2   # so we now have ste and not ste/2
                self.dudlstd[inda,indt,indl] = lline[3]

            ## get energy vs temperatrues data
            lambdas_array = avg_dudl[:,0]
            #energies_array = self.dudlmean[inda, indt, :]
            energies_array = avg_dudl[:,1]

            ## only make a fit if we have > 3 lambda values
            if len(lambdas_array) <= 3:
                if verbose:
                    print "only following lambdas:",lambdas_array," therefore no fit!"
                if len(lambdas_array) ==1:
                    if lambdas_array[0] == 0.5:
                        #print "yes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                        self.dudlmeanfit[inda,indt] = energies_array[0]
                ## but if only lambda 0.5 one could take this value instead...
                continue

            if verbose:
                print "La;",lambdas_array
                print "LA;",energies_array
            if self._verbose:
                print "i:",i,lambdas_array,type(lambdas_array),energies_array[0],type(energies_array)
            fit = get_dudlmeanfit(lambdas_array, energies_array, return_fit = True)
            self.dudlmeanfit[inda,indt] = fit[0]
            self.dudlmeanfit_tanfit[inda,indt] = np.transpose(fit[1])
        return


    def save_data(self):
        print "saving data ..."
        np.savez_compressed( "dudl",                \
                files       = self.files,           \
                data        = self.data,            \
                dudl        = self.dudl,            \
                a           = self.a,               \
                t           = self.t,               \
                l           = self.l,               \
                s           = self.s,               \
                dudlmean    = self.dudlmean,        \
                dudlmeanfit = self.dudlmeanfit,     \
                dudlste     = self.dudlste,         \
                dudlstd     = self.dudlstd,         \
                ah          = self.ah,              \
                ahste       = self.ahste,           \
                ahstd       = self.ahstd,           \
                _filestring = self._filestring,     \
                _pwd        = self._pwd             \
                )
        return

    def load_data(self):
        ''' '''
        print "loading data ..."
        var = np.load( "dudl" + ".npz" )

        self.files          = var['files']
        self.data           = var['data']
        self.a              = var['a']
        self.t              = var['t']
        self.l              = var['l']
        self.s              = var['s']
        self.dudlmean       = var['dudlmean']
        self.dudlmeanfit    = var['dudlmeanfit']
        self.dudlste        = var['dudlste']
        self.dudlstd        = var['dudlstd']
        self.ah             = var['ah']
        self.ahste          = var['ahste']
        self.ahstd          = var['ahstd']
        return

    def get_dudls_and_averages_from_self_data(self):
        ''' 1st index: volume / lattice constant    (self.a)
            2nd index: temperature                  (self.t)
            3rd index: lambda                       (self.s)'''
        print "building dudls ..."
        self.dudl = np.empty((  self.a.size, \
                                self.t.size, \
                                self.l.size, \
                                self.s.size), dtype=object)
        self.dudl[:] = np.NAN

        self.dudlmean     = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudllen      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlste      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlstd      = np.empty([self.a.size, self.t.size, self.l.size])
        self.dudlmeanlen  = np.empty([self.a.size, self.t.size])
        self.dudlmeanfit  = np.empty([self.a.size, self.t.size])

        self.dudlmean[:]     = np.NAN
        self.dudllen[:]      = np.NAN
        self.dudlste[:]      = np.NAN
        self.dudlstd[:]      = np.NAN
        self.dudlmeanlen[:]  = np.NAN
        self.dudlmeanfit[:]  = np.NAN

        #self.dudlavg = np.copy(self.dudl)

        np.set_printoptions(suppress=True)   # display arrays withou 000000
        np.set_printoptions(precision=6)    # print only 6 digist after .

        for a in range(self.a.size):
            for t in range(self.t.size):
                for l in range(self.l.size):
                    data = np.array([])
                    for s in range(self.s.size):
                        if np.isnan(np.any(f.data[a,t,l,s])):
                            # isnan -> no information
                            #print a,t,l,s,"NONE"
                            pass
                        else:
                            # This is the only input file specifit line
                            # sphinx
                            self.dudl[a,t,l,s] = self.data[a,t,l,s][:,7]

                            # remove first steps of equilibration
                            cut = self._dudl_cutoff_first_structures
                            self.dudl[a,t,l,s] = self.dudl[a,t,l,s][cut:]

                            # concat all dudls
                            data = np.concatenate((data,self.dudl[a,t,l,s]))
                            #print a,t,l,s,"--> DATA"

                            # get dudlavg
                            #self.dudlavg[a,t,l,s] = \
                            #[self.dudl[a,t,l,s][0:i].mean() for i in \
                            #np.arange(len(self.dudl[a,t,l,s]))+1] # sphinx

                    #print a,t,l,s,"||", self.a[a], self.t[t], self.l[l]
                    #print "data:",data
                    corl = self._dudl_correlation_length
                    if data is not np.array([]):
                        self.dudlmean[a,t,l]      = data.mean()
                        self.dudlste[a,t,l]   = \
                            data[0::corl].std()/np.sqrt(len(data[0::corl]))
                        self.dudlstd[a,t,l]   = data.std()
                        self.dudllen[a,t,l]   = len(data)
                #print "a:",a,"t:",t
                # a and t are just an index
                # building dudls ...
                # a: 0 t: 0
                # a: 0 t: 1
                # a: 0 t: 2
                # a: 0 t: 3
                # a: 0 t: 4
        self.dudlmeanfit = duelmean_to_dudlmeanfit(self.dudlmean, self.l)

                #lambdas_array = self.l
                #energies_array = self.dudlmean[a, t, :]
                ##self.dudlmeanfit[a,t] = self.get_dudlmeanfit(a, t)
                #print "ii:",lambdas_array,type(lambdas_array),energies_array,type(energies_array)
                #self.dudlmeanfit[a,t] = get_dudlmeanfit(lambdas_array, energies_array)
                ##print "sdmf:",self.dudlmeanfit[a,t]
        return


    def dudlmean_to_dudlmeanfit(dudlmean = None, lambdas = None):
        if dudlmean == None:
            sys.exit("expecting dudlmean to be an array[a,t,l] of energies")
        if type(dudlmean) != np.ndarray:
            sys.exit("expecting dudlmean to be an array[a,t,l] of energies")
        if lambdas == None:
            sys.exit("expecting lambdas to be an array; e.g. [0.0, 0.15, 0.5, 0.85, 1.0]")
        if len(lambdas) != dudlmean.shape[-1]:
            sys.exit("len(lambdas) != dudlmean.shpae[-1] which it has to be!")

        dudlmeanfit  = np.empty([dudlmean.shape[0], dudlmean.shape[1]])
        dudlmeanfit[:]  = np.NAN
        for inda in np.arange(dudlmean.shape[0]):
            for indt in np.arange(dudlmean.shape[1]):
                energies = dudlmean[inda, indt, :]
                dudlmeanfit[a,t] = get_dudlmeanfit(lambdas, energies)
        return dudlmeanfit



    def get_dudl_vs_temp_folder(self, verbose = False):
        ''' needs:
            - self.dudlmean
            - self.a
            - self.t
            - self.dudldste
            '''
        def quadFunction(x, a0):
            #return a0*x**2.718
            return a0*x**2.
            #return -1.+np.exp(a0*x)
            #return -1.+x*np.exp(a0)
            #return -1.+a0*np.exp(a0*x)

        folder = "auswertung_dudl_vs_temp"

        # remove folder if existing and create new one
        if os.path.isdir(folder):
           shutil.rmtree(folder)
        os.makedirs(folder)

        self.dudlmean_from_quadfit = np.copy(self.dudlmean)

        #for inda,a in np.arange(len(self.a)):
        for inda,a in enumerate(self.a):
            #for indl,l in np.arange(len(self.l)):
            for indl,l in enumerate(self.l):
                if verbose:
                    print ""
                    print "a:",a," l:",l
                #print "a:",self.a[a],"l:",self.l[l]
                x = self.t
                y = self.dudlmean[inda,:,indl]
                dy = self.dudlste[inda,:,indl]
                dellist = np.array([])
                for iiind,ii in enumerate(y):
                    #print "Y:",iiind,y[iiind]
                    if np.isnan(y[iiind]) == True:
                        dellist = np.append(dellist,iiind)
                        #print "K:",iiind
                #print "dellst:",dellist,len(dellist)
                #print ">x:",x
                #print ">y:",y
                #print "?1:",y[~np.isnan(y)]   #~ reverts np.isnan
                #print "?2:",~np.isnan(y)   #~ reverts np.isnan
                #print "?x:",x
                #print "?y:",y
                y = y[~np.isnan(y)]   #~ reverts np.isnan
                # NO THIS IS HERE WRONG: !x = x[~np.isnan(y)]   #~ reverts np.isnan
                if len(dellist) > 0:
                    x = np.delete(x,dellist)
                    dy = np.delete(dy,dellist)
                    ## NO THIS IS HERE WRONG: y = np.delete(x,dellist)
                #print "!x:",x
                #print "!y:",y

                if len(x) != len(y) != len(dy):
                    sys.exit("len(x):"+str(len(x))+" is not len y:"+str(len(y)))
                if verbose:
                    print "x:",x
                    print "y:",y
                if len(x) != len(y):
                    print "len x != len y !!!!!!!!! ERROR"
                    continue
                if len(x) == 0:
                    continue
                #print "dy:",dy
                #print ""
                filebase = folder+"/"+str(a)+"_"+str(l)
                np.savetxt(filebase+"_data",np.transpose([x,y,dy]), fmt="%.0f %.3f %.3f")

                from scipy.optimize import curve_fit
                parameters, covariance = \
                    curve_fit(quadFunction, x, y, maxfev=33000)

                #print "pars:",parameters
                #print "pars:",parameters[0]
                # fit
                points = 101
                fit_x_min = 0
                fit_x_max = self.t.max()

                fit_x_values = np.linspace(fit_x_min, fit_x_max, points)
                fit_y_values = quadFunction(fit_x_values, *parameters)
                np.savetxt(filebase+"_fit",np.transpose([fit_x_values,fit_y_values]), fmt="%.0f %.6f")

                for indt,t in enumerate(self.t):
                    self.dudlmean_from_quadfit[inda,indt,indl] = quadFunction(t,*parameters)

        # make Fqh surface from quad fit
        self.dudlmean_from_quadfit_samenans = np.copy(self.dudlmean_from_quadfit)
        self.dudlmean_from_quadfit_samenansfit = np.copy(self.dudlmeanfit)
        self.dudlmean_from_quadfit_samenansfit[:] = np.NAN
        for a in np.arange(self.dudlmean.shape[0]):
            for t in np.arange(self.dudlmean.shape[1]):
                for e in np.arange(self.dudlmean.shape[2]):
                    print "a:",a," t:",t,self.dudlmean[a,t,e]
                    if np.isnan(self.dudlmean[a,t,e]) == True:
                        self.dudlmean_from_quadfit_samenans[a,t,e] = np.nan
                #print "???",self.l, self.dudlmean_from_quadfit_samenans[a,t]
                #print "???", get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans[a,t])
                self.dudlmean_from_quadfit_samenansfit[a,t] = get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans[a,t])
        for a in np.arange(self.dudlmeanfit.shape[0]):
            for t in np.arange(self.dudlmeanfit.shape[1]):
                if np.isnan(self.dudlmeanfit[a,t]) == True:
                    self.dudlmean_from_quadfit_samenansfit[a,t] = np.nan

        # make Fqh surface from quad fit only highest temperature
        self.dudlmean_from_quadfit_samenans_highestT = np.copy(self.dudlmean)
        self.dudlmean_from_quadfit_samenans_highestTfit = np.copy(self.dudlmeanfit)
        self.dudlmean_from_quadfit_samenans_highestTfit[:] = np.NAN
        for a in np.arange(self.dudlmeanfit.shape[0]):
            for t in np.arange(self.dudlmean.shape[1]):
                for e in np.arange(self.dudlmean.shape[2]):
                    if t != self.dudlmean.shape[1]-1:
                        self.dudlmean_from_quadfit_samenans_highestT[a,t,e] = np.nan
                print "??? a:",a,"t:",t,self.l, self.dudlmean_from_quadfit_samenans_highestT[a,t]
                print "??? a:",a,"t:",t, get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans_highestT[a,t])
                print "!!!", np.isnan(self.dudlmean_from_quadfit_samenans_highestT[a,t])
                if np.isnan(get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans_highestT[a,t])):
                    print "isnan!!!!!!!!!!!!!"
                self.dudlmean_from_quadfit_samenans_highestTfit[a,t] = get_dudlmeanfit(self.l, self.dudlmean_from_quadfit_samenans_highestT[a,t])

        return






    def plot_dudl_vs_lambda(self, a, t, l):
        ''' TODO: NOT READY YET
            1st index: a
            2nd index: t
            3rd index: l '''
        energies_array = None
        if utils.is_int(a) is True:
            a = int(a)
        elif a is ":":
            a = eval(":")
            #energies_array = self.dudlmean[:]
        else: sys.exit("index a:"+str(a)+" not correct")

        energies_array = self.dudlmean[a]

        #if utils.is_int(t) is True:
        #    energies_array = energies_array[int(t)]
        #elif t is ":":
        #    energies_array = energies_array[:]

        print "energies_array:",energies_array
        return


    def import_fah_data(self, filename = None):
        ''' '''
        self.surface = None
        possible_filename = [ "Fah_surface" ]
        if filename == None:
            for i in possible_filename:
                filecheck = glob.glob(i)
                if filecheck:
                    if len(filecheck) == 1:
                        _printgreen("importing fah data : "+filecheck[0]+" "+"."*20)
                        self.surface = np.loadtxt(filecheck[0])
                        break

        if self.surface == None:
            sys.exit("fah surface not found!")
        _printgreen("importing fah data ............................... DONE")
        print ""
        return self.surface


    def import_surface(self, surface_filename = None, interpolate = False, sysexit = True):
        ''' '''
        surface = None
        if type(surface_filename) == bool or surface_filename == None:
            surface_filename = find_surface_filename(sysexit = sysexit)

        # check surface_filename
        #print "000:",surface_filename
        if type(surface_filename) == bool or surface_filename == None:
            if sysexit:
                sys.exit("Did not found any Fah surface file in current folder")
            else:
                return False

        if os.path.isfile(surface_filename) != True:
            if sysexit:
                sys.exit("Did not found any Fah surface file in current folder")
            else:
                return False

        self.surface_filename = surface_filename
        surface = np.loadtxt(surface_filename)
        if interpolate:
            surface = h.fqh_surface_interpolate(surface)
        return surface

    def import_surface_expr(self, surface_filename = None, surface = None, sysexit = True):
        if surface == None:
            surface = self.import_surface(surface_filename = surface_filename, sysexit = sysexit)
        try:
            if surface == False:
                return False
        except ValueError:
            pass # in this case we have an array and everything is correct
        self.surface_expr = np.array([])
        import sympy
        V = sympy.Symbol('V')
        for ind, temp in enumerate(np.arange(surface.shape[0])):
            add =    surface[ind][1] \
                    +surface[ind][2]*V \
                    +surface[ind][3]*V*V
            if surface.shape[1] == 4:
                self.surface_expr = np.append(self.surface_expr, add)
            elif surface.shape[1] == 5:
                self.surface_expr = np.append(self.surface_expr, add \
                    +surface[ind][4]*V*V*V)
        return self.surface_expr

if __name__ == '__main__':
    f = ah()


    p = argparse.ArgumentParser(description='''help string''')
    p.add_argument('-ws', default=False, action='store_true',
        help='Write Fqh_surface (Fah_surface and Fah_surface_x2) start this in high folder')
    p.add_argument('-lla', default=False, action='store_true',
        help='load low  avg_dUdL_fre files'),
    p.add_argument('-v', '--verbose',
                help='verbose', action='store_true', default=False)
    args = p.parse_args()

    if args.ws:
        if args.verbose:
            print "###################### 1"
        f.import_avg_dudl_data(avg_dudl_filenames = "*Ang_*K/avg_dUdL_lowplushigh_fre", verbose = args.verbose)

        if args.verbose:
            print "###################### 2"
        f.get_dudl_vs_temp_folder(verbose = args.verbose)
        if args.verbose:
            print "###################### 3"
        write_Fah_surface(a = f.a, temps = f.t, dudlmeanfit = f.dudlmeanfit)
        write_Fah_surface(a = f.a, temps = f.t, dudlmeanfit = f.dudlmean_from_quadfit_samenansfit, filename = 'Fah_surface_x2')


    if args.lla:
        f.import_avg_dudl_data(avg_dudl_filenames = "*Ang_*K/avg_dUdL_fre", verbose = False)

    #read_from_files = True
    #if read_from_files:
    #    f = ah()
    #    f._filestring             = "4.16Ang_*K/lambda*/moldyn.dat*"
    #    f.initialize_files_data()
    #    f.import_dudl_files_data()
    #    f.save_data()
    #    f.get_dudls_and_averages_from_self_data()
    #    f.save_data()
    #else:
    #    pass
    #    #f.load_data()


    #f.get_dudls_and_averages_from_self_data()


    #a=1
    #plot(f.l,f.dudlmean[a,0],f.l,f.dudlmean[a,1],f.l,f.dudlmean[a,2],f.l,f.dudlmean[a,3],f.l,f.dudlmean[a,4])


    ## Fah vs volume
    #a=0
    #plot(f.a,f.dudlmeanfit[a])
    #

    ## Fah vs temperature
    #t=0
    #plot(f.t,f.dudlmeanfit[:,t]
