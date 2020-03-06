#!/usr/bin/env python

#import pylab
import numpy as np
import glob
import sys
import re
import os
import scipy as sp
from scipy.interpolate import griddata
import myutils as my

# for plots
from mpl_toolkits.mplot3d import *
import matplotlib
#import matplotlib.pyplot as plt
from random import random, seed
from matplotlib import cm
import imp

# own stuff
#import hesse as h
#reload(h)    # necessary when working in ipython and changing code on the run
# to show position of module
# print u.__file__

np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=6)    # print only 6 digist after .

# functions with _name are not seen when imported
def _printred(var):
    ENDC = '\033[0m'
    red = '\033[31m'
    print(red + str(var) + ENDC)

def _printgreen(var):
    ENDC = '\033[0m'
    red = '\033[32m'
    print(red + str(var) + ENDC)


def _string_to_num_list(string, tostring = False):
    """ Turn a string into a list of int/float chunks.
        "4.Ang_300K_.2m-9m=.77k" -> [4.0, 300.0, 0.2, 9.0, 0.77] """

    def tryfloat(stringtest):
        try:
            return float(stringtest)
        except:
            return stringtest
    import re
    out = [ tryfloat(c) for c in re.findall(r"[-+]?\d*\.\d+|\d+", string)]
    if tostring:
        return [ str(i) for i in out]
    return out

def _lsn(searchstring):
    """ Returns files sorted numerically
        searchstring: e.g. "*Ang_*K"            """
    import glob
    files = glob.glob(searchstring)
    return sorted((files), key=_string_to_num_list)

def _common_prefix(strings, common_suffix = False):
    """ Find the longest string that is a prefix (or suffix) of all the strings of the input list.
    """
    if not strings:
        return ''

    if common_suffix == True:
        strings = [i[::-1] for i in strings]

    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    if common_suffix == True:
        return prefix[::-1]
    else:
        return prefix


class qh(object):
    def __init__(self, args = None):
        self._verbose = False
        self._ignore_volumerange_error = False
        self._oqh = False
        self._filessearchstring = ""
        self._filessearchstring_common = ""
        self._filessearchstring_suffix = ""
        self._alat_fcc = False
        self._alat_bcc = False
        self._alat_dfcc = False
        self._atoms = False
        self._multiply_energy = 1.
        self._multiply_volume = 1.
        self._divide_volume = 1.
        self._divide_energy = 1.
        self._plots = False
        self._verbose = False
        self.__verbose = False
        self._wqh1 = None
        self._wqh2 = None
        self._wqh3 = None
        self.directory = None
        self.files = None
        self.surface_filename = None
        self.data = None

        self.temperatures = None
        self.volumes = np.array([])

        self.surface_2nd = None
        self.surface_3rd = None
        self.surface_2nd_expr = None
        self.surface_3rd_expr = None
        self.surface = None
        self.surface_expr = None

        if args: # dont use: if args != None
            self.files = args.inputfiles
            self._filessearchstring = args.filessearchstring
            self._divide_energy = args.divide_energy
            self._multiply_energy = args.multiply_energy
            self._multiply_volume = args.multiply_volume
            self._divide_volume = args.divide_volume
            self._alat_fcc = args.fcc
            self._alat_bcc = args.bcc
            self._alat_dfcc = args.dfcc
            self._verbose = args.verbose
            self._wqh1 = args.wqh1
            self._wqh2 = args.wqh2
            self._wqh3 = args.wqh3
            self._plots = args.plots
            self._oqh = args.oqh
            self._ignore_volumerange_error = args.ie
            if self._verbose == True:
                print("args:",args)

        if self._multiply_energy != 1.:
            self._multiply_energy = float(self._multiply_energy)
        if self._multiply_volume != 1.:
            self._multiply_volume = float(self._multiply_volume)
        if self._divide_energy != 1.:
            self._divide_energy = float(self._divide_energy)
        if self._divide_volume != 1.:
            self._divide_volume = float(self._divide_volume)



    def help(self, p = None):
        string = '''
        Import Fqh files. Fqh filename need to contain the volume per atom (e.g. Fqh_fromMeshFreqs_16.*).
        Fqh file 1st column --> Temperature [K]; 2nd column --> energy [meV];
        '''
        import argparse
        if p == None:
            p = argparse.ArgumentParser(description=string)

        p.add_argument('-s',   '--filessearchstring',
            help='input files search string (wildcards as *.? are allowed) e.g. Fqh_fromMeshFreqs_*)',
            type=str, default=None)
        p.add_argument('-i','--inputfiles',nargs='+',
            help='Use given files as Fqh inputfiles; You can use wildcards (*.{}); e.g. Fqh_4.*;', default=None)
        p.add_argument('-me', '--multiply_energy',
            help='multiplies the 2nd column of inputfile by [int]',type=int, default=1.)
        p.add_argument('-de', '--divide_energy',
            help='divedes the 2nd column of inputfile by [int]',type=int, default=1.)
        p.add_argument('-mv', '--multiply_volume',
            help='multiplies volume by [integer]',type=float, default=1.)
        p.add_argument('-dv', '--divide_volume',
            help='devides volume by [integer]',type=float, default=1.)
        #p.add_argument('-m', '--multiply',
        #    help='multiply (scale) energy and volume by [integer]',type=float, default=1.)
        p.add_argument('-fcc',
            help='inputfilename has fcc alat instead volume (alat**3/4 to get volume)',
            action='store_true', default=False)
        p.add_argument('-bcc',
            help='inputfilename has bcc alat instead volume (alat**3/2 to get volume)',
            action='store_true', default=False)
        p.add_argument('-dfcc',
            help='inputfilename has doublefcc alat instead volume (alat**3/8 to get volume)', action='store_true', default=False)
        p.add_argument('-wqh1', default=False, action='store_true',
            help='Write Fqh_surface_1st_order')
        p.add_argument('-wqh2', default=False, action='store_true',
            help='Write Fqh_surface_2nd_order')
        p.add_argument('-wqh3', default=False, action='store_true',
            help='Write Fqh_surface_3rd_order')


        p.add_argument('-ih', default=False, nargs='?',
            help='define input helmholtz file for fitting')

        p.add_argument('-iqh', default=False, nargs='?',
            help='define input surface')
        p.add_argument('-oqh', default=False, nargs='?',
            help='if surface was scaled use this falg to save outpuf surface')

        p.add_argument('-ov', default=False, type=float,
            help='output Helmholz free energy from surface for given volume')
        p.add_argument('-ofcca', default=False, type=float,
            help='output Helmholz free energy from surface for given fcc alat')
        p.add_argument('-ob', default="Fqh_fitted_", type=str,
            help='basenaem for fitted Helmholz free energy file')

        p.add_argument('-p','--plots', default=False, type=int,
            help='Write analysis plots')
        p.add_argument('-v', '--verbose',
            help='be more verbose', action='store_true', default=False)
        p.add_argument('-q', '--quickanalysis',
            help='quick analysis fqh  & evinet', action='store_true', default=False)
        p.add_argument('-ie',
            help='ignore error due to volume range',
            action='store_true', default=False)
        return p

    def import_fqh_files(self, scale=None, vol_min_error=10, vol_max_error=35):
        '''
        e.g.
            files = "Fqh_fromMeshFreqs_*"

        description:
            imports a list of files;
            Expects 1st column to be the temperature [K]
            Expects 2nd column to be the energies [meV]
            Expects the volume per atom in [Angstrom^3] to be given in the filename
        '''

        possible_filenames = [ \
        "Fvib_fromExactFreqs_1_[0-9]*", \
        "Fqh_fromExactFreqs_[0-9]*", \
        "Fqh_fromMeshFreqs_[0-9]*", \
        "Fqh_fromMeshfreqs_[0-9]*", \
        "Fqh_fromMeshFreq_[0-9]*", \
        "Fqh_fromMeshfreq_[0-9]*", \
        "Fqh_fromMsehFreq_[0-9]*",\
        "Fqh_fromMsehfreq_[0-9]*",\
        "Fqh_[0-9]*", \
        "FvibSupercell_perAtom_[0-9]*", \
        #"FvibSupercell_perAtom_[0-9]*", \
        "FvibSupercell_[0-9]*" \
        "Fqh_*" \
        ]

        ############################################################################
        # get self.files (inputfiles)
        ############################################################################
        print('self.files (0)',self.files)
        if self.files == None:
            if self._filessearchstring != None: # specified "Fqh_from_xyz*"
                self.files = _lsn(self._filessearchstring)

            print('self.files (1)',self.files)
            if self._filessearchstring == None: # no searchstring specified ("Fqh_from_xyz*")
                #_printred("You didn't sepcify files to import. Trying default names of Fqh files:")
                for files in possible_filenames:
                    filenames = _lsn(files)
                    print("filenames:",filenames,"len:",len(filenames))
                    if filenames == None:
                        sys.exit("No filenames specified")
                    if len(filenames) is not 0:
                        self._filessearchstring = files
                        #_printgreen("   found : "+files)
                        _printgreen("import_fqh_files: .............................................")
                        self.files = filenames
                        break

                if self.files == None:
                    print("No Fqh files found; trying Yaml files;")
                    filenames = glob.glob("thermal_properties*.yaml")
                    import yaml
                    #with open(filenames[0]) as f:
                    #    stream = yaml.load(f)
                    if len(filenames) == 0:
                        sys.exit("No Fqh_xxx files or thermal_properties*.yaml files found")
                    stream = open(filenames[0],'r')
                    ka = yaml.load(stream)
                    print("done")
                    #stream = open('thermal_properties-1.yaml','r')
                    #print yaml.load(stream)

        #########################################################################
        # check if files exist
        #########################################################################

        for i in self.files:
            print('A',i)
        self.files = [f for f in self.files if '_fit_' not in f]
        for i in self.files:
            print('B',i)
        self.files = [f for f in self.files if '_Surface_' not in f]
        for i in self.files:
            print('C',i)
        for i in self.files:
            print(i)
        for i in self.files:
            if os.path.isfile(i) != True:
                sys.exit(i+" does not exist!")

        #########################################################################
        # get common name of files for ...surface_2nd_order
        # get common name of files for ...surface_3rd_order
        #########################################################################
        ka = _common_prefix(self.files)
        import re
        ka = re.sub('\.$','',ka)
        ka = re.sub('[0-9]*$','',ka)
        ka = re.sub('\.$','',ka)
        self._filessearchstring_common = ka
        #print "self._filessearchstring_common:", self._filessearchstring_common, len(self._filessearchstring_common), i[len(self._filessearchstring_common):]

        #########################################################################
        # get suffix of files ( volumes )
        #########################################################################
        self._filessearchstring_suffix = ""
        print('self.files',self.files)
        for i in self.files:
            #print('i',i)
            suffix = i[len(self._filessearchstring_common):]
            #print('suffix',suffix,len(suffix),'files',len(self.files),suffix[:10])
            pos_dot = len(suffix)
            for pos_dot_,j in enumerate(suffix):
                if j == ".":
                    pos_dot = pos_dot_
                    break
            #print('pos_dot',pos_dot)
            digits = 4
            if pos_dot+digits < len(suffix):
                suffix = suffix[:pos_dot+digits+1]
            #print('suffix',suffix,len(suffix),'files',len(self.files),suffix[:10])
            #sys.exit('3')
            self._filessearchstring_suffix = self._filessearchstring_suffix + "_" + suffix

        if self.files == None:
            sys.exit("No Fqh files found (1)")

        if len(self.files) == 0:
            sys.exit("No Fqh files found (2)")

        # get self.directory
        if len(self.files) >= 1:
            self.directory = "/".join(os.path.realpath(self.files[0]).split("/")[:-1])

        ############################################
        # loop throuth inputfiles
        ############################################
        self._datain = None
        self._temperaturesin = None
        _printred( "make check of 0-point vibrations vs volume")
        print("file:                            volume:                 T=0K Energy:")
        print("________________________________________________________________________")

        for i in range(len(self.files)):
            file = self.files[i]

            # get volume string and volume
            string_in = re.findall(r"[-+]?\d*\.\d+|\d+", file)
            #print "string_in:",string_in,"type:",type(string_in),len(string_in)
            #print "string_in:", string_in, len(string_in)
            if type(string_in) is str:
                string_in = string_in
            elif type(string_in) is list:
                string_in = string_in[-1]
            else:
                sys.exit(
                    "Found several (or none) real numbers in string but wanted one.\n\
                Please check the filename.\n\
                filename          : " + str(file) + "\n\
                    strings_real      : " + str(string_in) + "\n\
                    type(strings_real): " + str(type(string_in)) + "\n\
                    len(string_in)  : " + str(len(string_in))) #+\

            # get volume, checking if vlaue makes sence: see further down
            string_vol = float(string_in)
            if self._alat_bcc == True:
                vol = string_vol**3/2
            elif self._alat_fcc == True:
                vol = string_vol**3/4
            elif self._alat_dfcc == True:
                vol = string_vol**3/8
            else:
                vol = string_vol

            vol = vol*self._multiply_volume/self._divide_volume
            # self.data (import fqh files)
            if self._verbose:
                print("_multiply_volume",self._multiply_volume)
                print("_divide_volume",self._divide_volume)
                print("file:",file)
            ############################################
            # import files
            ############################################
            #data = pylab.loadtxt(file)
            data = np.loadtxt(file)
            ############################################
            # check imported files
            ############################################
            if len(data.shape) != 2:
                print("files:",self.files)
                print("------------------")
                print(data)
                print("------------------")
                sys.exit("(1) check filenames, did you want to import all this files?")
            if data.shape[1] != 2:
                with open(file, 'r') as f:
                      first_line = f.readline()
                if first_line.split() == ['#Input', 'files', 'are', '[phonon.sxb];']:
                    data = data[:,[0,1]]
            if data.shape[1] != 2:
                print("data.shape:",data.shape)
                print("files:",self.files)
                print("------------------")
                print(data)
                print("------------------")
                sys.exit("(2) check filenames, did you want to import all this files?")



            #print "danach:",file
            tempin = data[:,0]
            enein = data[:,1]*self._multiply_energy/self._divide_energy
            data[:,1] = data[:,1]*self._multiply_energy/self._divide_energy

            # self.volumes (append)
            self.volumes = np.append(self.volumes,[vol])
            print("[{one}, \t\t{two}\t\t{three} <-- check VOLUME with OUTCAR!]".format(one=file, two=vol, three=enein[0]))

            # self._datain (append)
            if self._datain == None:
                self._datain = np.empty(len(self.files),dtype=object)
                self._datain[:] = np.NAN
            self._datain[i] = data


            # self._datain (append)
            #if self._datain == None:
            #    self._datain = np.array(data, dtype=object)
            #    #self._datain = np.empty(len(self.files),dtype=object)
            #    #self._datain[:] = np.NAN
            ##self._datain[i] = data
            #else:
            #    pass
            #    #self._datain.append(data)

            #print('self._temperaturesin:', self._temperaturesin)
            #print('tempin              :', tempin)
            if self._temperaturesin is None:
                self._temperaturesin = np.unique(tempin)
            self._temperaturesin = np.intersect1d(self._temperaturesin,tempin)


        #print "danach:",file
        # self.temperatures
        temp_min = 0.
        temp_max = self._temperaturesin.max()
        self.temperatures = np.arange(temp_min, temp_max + 1.)

        # self.data
        #print("len self.files:",len(self.files))
        #print("temp_max",temp_max)
        self.data = np.empty((len(self.files),int(temp_max) + 1))
        #print "dim:",self.data.shape
        self.data[:] = np.NAN


        # create mask, get indizes indes of list
        # to get indexed elemetnes:
        # a = np.array([1,2,3,4,5,6,7,8,9])
        # a[[0,1,6]] --> [1,2,7]
        # to get index???  ... geht google A==B
        # _printgreen("create surface ...")
        if self._verbose:
            print("lenL:",len(self._datain))
            print("self._datain",self._datain[0])
        for volindex in np.arange(len(self._datain)):
            #print "v:",volindex
            tempin = self._datain[volindex][:,0]
            enein = self._datain[volindex][:,1]
            #print "volindex:",volindex
            #print "[tempin]:",tempin.shape

            #print "tempin:",tempin
            #print "tempin[0]:",tempin[0]
            #print "tempin[1]:",tempin[1]
            #print "tempin[-1]:",tempin[-1]
            #print "tempin[-2]:",tempin[-2]
            #print "data.shape:",self._datain[volindex].shape
            #print "datain[:,0,-1]:",self._datain[volindex][:,0][-1]
            #print "datain[:,0,-2]:",self._datain[volindex][:,0][-2]
            #print "datain[:,1,-1]:",self._datain[volindex][:,1][-1]
            #print "datain[:,1,-2]:",self._datain[volindex][:,1][-2]
            #print('enein',enein)
            #print('tempin',tempin)
            #print('tempin',tempin.astype(int))
            self.data[volindex,[tempin.astype(int)]] = enein

            # fill T = 0K (make equal to T = 1K or T = 2 K)
            if np.isnan(self.data[volindex,0]):
                if np.isnan(self.data[volindex,1]) == False:
                    self.data[volindex,0] = self.data[volindex,1]

            if np.isnan(self.data[volindex,0]):
                if np.isnan(self.data[volindex,2]) == False:
                    self.data[volindex,0] = self.data[volindex,2]

            if np.isnan(self.data[volindex,0]):
                sys.exit("No energies at T=0K found for volume:"\
                        +str(self.volumes[volindex]))

        # check shape of data and tempteratures
        if self.data.shape[1] != self.temperatures.shape[0]:
            sys.exit("shape of data and temperatures different")


        #def uniquearr(a):
        #    order = np.lexsort(a.T)
        #    a = a[order]
        #    diff = np.diff(a, axis=0)
        #    ui = np.ones(len(a), 'bool')
        #    ui[1:] = (diff != 0).any(axis=1)
        #    return a[ui]

        # check volume
        volmin = self.volumes.min()
        if vol_min_error*self._multiply_volume/self._divide_volume < volmin < vol_max_error*self._multiply_volume/self._divide_volume:
            pass
        else:
            e1 = "The minimum volume imported, "+str(volmin)+" [Angstrom^3], "
            e2 = "seems souspiciously big/small to be scaled to volume per atom;"
            e3 = "Typical volumes per atom are between "+str(vol_min_error)+" and "
            e4 = str(vol_max_error)+" [Angstrom^3]. Are you sure this is scaled per "
            e5 = "atom? In case this is a volume per atom you can indrease/decrease "
            e6 = "the vol_max_error vol_min_error option. || "
            e8 = str(vol_min_error*self._multiply_volume/self._divide_volume)+" < "+str(volmin*self._multiply_volume/self._divide_volume)+" < "+str(vol_max_error)
            e9 = "To ignore this use the -ie option"
            if self._ignore_volumerange_error == True:
                pass
            else:
                del self.volumes
                del self.data
                del self.temperatures
                sys.exit(e1 + e2 + e3 + e4 + e5 + e6+ e8 + e9)
        # _printgreen("removing nans ...")
        a = np.transpose(self.data)
        self.data = np.transpose(a[~np.isnan(a).any(1)])
        self.temperatures = self.temperatures[~np.isnan(a).any(1)]
        _printgreen("import_fqh_files: ............................................. DONE")
        print("")
        return


    def fit_surface(self, volume=None):
        ''' energy vs volume: polynomial fits of order 2 and 3 for every temperature
        defines:
        - surf_fit_{2n,3r}d:       [[t1, a1, b1, c1, ...],
                                    [t2, a2, b2, c2, ...],
                                    ...
                                    [tn, an, bn, cn, ...]]


        - surf_coefs_{2n,3r}d:     [[c1, b1, a1, ...],
                                    [c2, b2, a2, ...],
                                    ...
                                    [cn, bn, an, ...]]

               where n runs over all temperatures
               and a, b, c, d, .... are the fitting coeffitients
                 order n=2  : a + bx + cx^2
                 order n=3  : a + bx + cx^2 + dx^3

               ...
        output alternatives when set output to:
           "vin"        : array of volumes on which fit is performed
           "enein"      : array of energies on which fit is performed
           "enefit"     : array of fitted energies at vin
           "enediff"    : enein - enefit
           "enediffmax" : (enein - enefit).max

           returnes:
           ([     v1 ,      v2 ,      v3 , ...,      vn ] == "vin",
            [   y(v1),    y(v2),    y(v3), ...,    y(vn)] == "enein",
            [yfit(v1), yfit(v2), yfit(v3), ..., yfit(vn)] == "enefit",
                   3        2
           ([-252 x + 2577 x - 9185 x + 1.106e+04 ],
           [ d, c, b, a],
           [y(v1) - yfit(v1), y(v2) - yfit(v2), ... ],
           [ max y(v) - yfit(v) ])
        '''

        _printgreen("fit_surface ...................................................")

        # define surface
        np.set_printoptions(precision=12)   # to make the output exact
        self.surf_coefs_1st = np.empty((len(self.temperatures),2))  # -wf1
        self.surf_coefs_2nd = np.empty((len(self.temperatures),3))  # -wf2
        self.surf_coefs_3rd = np.empty((len(self.temperatures),4))  # -wf3

        self.surface_1st  = np.empty((len(self.temperatures),3))
        self.surface      = np.empty((len(self.temperatures),4))
        self.surf_out_3rd = np.empty((len(self.temperatures),5))

        self.surf_function = []

        # define fits
        self.data_fit = np.empty((len(self.volumes),len(self.temperatures)))
        fitpoints = 100
        fitgrid = np.linspace(self.volumes.min(),self.volumes.max(),fitpoints + 1)
        self._fitgrid = fitgrid
        self._data_fit = np.empty((fitpoints + 1,len(self.temperatures)))


        if self.data.shape[0] > 3:
            self.surf_function_3rd = []
            self.data_fit_3rd = np.empty((len(self.volumes),len(self.temperatures)))
            self._data_fit_3rd = np.empty((fitpoints + 1,len(self.temperatures)))

        for i,t in enumerate(self.temperatures):
            energies = self.data[:,i]
            if self.__verbose:
                print("temperature:", t)
                print("vol        :", self.volumes)
                print("energies   :", energies)

            # polynomial fit
            self.surf_coefs_1st[i] = np.polyfit(self.volumes, energies, 1)
            self.surf_coefs_2nd[i] = np.polyfit(self.volumes, energies, 2)
            self.surface_1st[i]    = np.append(self.surf_coefs_1st[i],[t])[::-1]
            self.surface[i]        = np.append(self.surf_coefs_2nd[i],[t])[::-1]
            self.surf_function.append(np.poly1d(self.surf_coefs_2nd[i]))


            #bounds = [10, 12]
            #crit_points = bounds + [x for x in p.deriv().r if x.imag == 0 and bounds[0] < x.real < bounds[1]]

            #crit = c.deriv().r
            #r_crit = crit[crit.imag==0].real
            #test = c.deriv(2)(r_crit)
            #x_min = r_crit[test>0]
            #print('i',i,crit_points)

            if self.data.shape[0] > 3:
                self.surf_coefs_3rd[i] = np.polyfit(self.volumes, energies, 3)
                self.surf_out_3rd[i]   = np.append(self.surf_coefs_3rd[i],[t])[::-1]
                self.surf_function_3rd.append(np.poly1d(self.surf_coefs_3rd[i]))

                #fit_fn = np.poly1d(self.surf_coefs_3rd)
                #ymin = sp.optimize.minimize(fit_fn,11)
                #c = np.poly1d(self.surf_coefs_3rd[i])
                #crit = c.deriv().r
                #v_ = crit[0]
                #print('i',i,ymin)

            #polycoefs: [  -252.00857010244    2576.739955929857  -9185.475148705054
            #         11059.807076691339]
            #function:    3        2
            #       -252 x + 2577 x - 9185 x + 1.106e+04

            # fits
            self.data_fit[:,i] = self.surf_function[i](self.volumes)
            self._data_fit[:,i] = self.surf_function[i](fitgrid)

            if self.data.shape[0] > 3:
                self.data_fit_3rd[:,i] = self.surf_function_3rd[i](self.volumes)
                self._data_fit_3rd[:,i] = self.surf_function_3rd[i](fitgrid)

            # TODO: c) write surface
            # TODO: d) read surface
        _printred("fit_surface ... here we might consider just to load the before fitted surf")
        _printgreen("fit_surface ................................................... DONE")
        print("")
        if self._wqh1: self.write_fqh_surface_1st()
        if self._wqh2: self.write_fqh_surface_2nd()
        if type(self._plots) != bool:  # Then it is a number
            self.plot_energy_vs_volume_qh(self._plots)
        if self.data.shape[0] > 3:
            if self._wqh3: self.write_fqh_surface_3rd()
        return

    def helmholtz_free_energy(self, v = False, fccalat = False, save = False):
        ''' get from the fqh surface a certain fqh for a certain volume'''
        if v == False and fccalat == False:
            sys.exit("define v or fccalat")
        if v != False and fccalat != False:
            sys.exit("define either v or fccalat, not both")
        out = np.empty([len(self.temperatures),2])
        out[:] = np.nan
        if v != False:
            stringout = str(v)
        if fccalat != False:
            v = (fccalat**3.)/4.
            stringout = str(fccalat)

        for i,t in enumerate(self.temperatures):
            out[i,0] = t
            out[i,1] = self.surf_function[i](v)
        if save != False:
            outputfile = "Fqh_fitted_"+stringout
            if save != True:
                outputfile = save+stringout
            np.savetxt(outputfile,out,fmt="%.0f %.12f")
            print((outputfile+" written ..."))
            return
        return out


    def find_surface_filenames(self, sysexit = True):
        possible_filename = [ "Fqh", "Fqh_surface", "Fqh*urface_2nd_order", "Fqh_b_*", "Fqh_d_*", "Fvib*", "Fqh_*fit*order*", "Fqh_*"]
        all_found = []
        for i in possible_filename:
            if self._verbose:
                print("iii:",i)
            filecheck = glob.glob(i)
            if filecheck:
                for j in filecheck:
                    #print "ldsasdf:",os.path.abspath(j)
                    all_found.append(os.path.abspath(j))

        # output
        if all_found != []:
            if len(all_found) == 1:
                self.surface_filename = all_found[0]
                return self.surface_filename
            else:
                return all_found

        # if we came to this point we have to exit
        if sysexit:
            sys.exit("Did not found any Fqh surface file in current folder")
        else:
            #print("Did not found any Fqh surface file in current folder")
            return None

    def import_surface(self, surface_filename = None, interpolate = True, sysexit = True):
        self.surface = None
        if type(surface_filename) == bool or surface_filename == None:
            surface_filename = self.find_surface_filenames(sysexit = sysexit)
            if self._verbose:
                print("Fqh surface_filename:",surface_filename)

        if type(surface_filename) == bool or surface_filename == None:
            sysexit = "Did not find any Fqh surface file in current folder"

        if type(surface_filename) == str:
            loopover = [surface_filename]
            self.surface_filename = os.path.abspath(surface_filename)

        if type(surface_filename) == list:
            loopover = surface_filename

        if type(surface_filename) != str:
            if type(surface_filename) != list:
                if sysexit == True:
                    sys.exit("No Fqh file found")
                else:
                    return None

        all_surfaces = []
        self.surface = None
        if self._verbose:
            print("loopover:",loopover)
        for ind, surface_filename_ in enumerate(loopover):
            # print "ind:",ind,"surface_filename_:",surface_filename_
            # check if file exists
            if surface_filename_ == None:
                sys.exit("Fqh surface file not given!")
            if os.path.isfile(surface_filename_) != True:
                sys.exit("Fqh surface "+str(surface_filename_)+" does not exist!")

            if self.surface == None:  #case of just one singele surface
                self.surface = np.loadtxt(surface_filename_,dtype='float64')
                #self.surface = np.loadtxt(surface_filename_)
                if self.surface.shape[1] == 4:
                    self.surface_2nd = self.surface
                elif self.surface.shape[1] == 5:
                    self.surface_3rd = self.surface
                self.temperatures = self.surface[:,0]
            else:                     # case of several surfaces
                print("ds")
                self.surface = [self.surface, np.loadtxt(surface_filename_)]



            ## case of several surfaces
            #if type(surface_filename) == list:
            #    print "aappend"
            #    all_surfaces = np.append(all_surfaces, np.loadtxt(surface_filename_))

        #print "as:",all_surfaces
        #self.surface = all_surfaces
        self.surface = np.float64(self.surface)

        print(self.surface)
        print("---")
        # e32 = e1 * 32 ->
        # v32 = v1 * 32 -> v1 = v32/32
        # e1  = a1 + b1v1 + c1v1^2 + d1v1^3
        # if we now multiply energy and volume by 32:
        # e32 = 32 * (a1 + b1(v32/32) + c1(v32/32)^2 + d1(v32/32)^3)
        #        ^  multiply energy
        #                        ^            ^               ^ multiply volume
        #                 self.surface[:,2]
        #                               self.surface[:,3]
        #                                           self.surface[:,4]
        if self._multiply_energy:
            self.surface[:,1:] = self.surface[:,1:]*self._multiply_energy

        #print "ss:",self.surface.shape[1]
        if self._multiply_volume:
            if self.surface.shape[1] > 5:
                sys.exit("add if higher order of surface then 3")
            if self.surface.shape[1] >= 3:
                self.surface[:,2] = self.surface[:,2]*(1./self._multiply_volume)
            if self.surface.shape[1] >= 4:
                self.surface[:,3] = self.surface[:,3]*(1./self._multiply_volume)**2.
            if self.surface.shape[1] == 5:
                self.surface[:,4] = self.surface[:,4]*(1./self._multiply_volume)**3.

            print("mv:",self._multiply_volume)
            print("me:",self._multiply_energy)


        if interpolate == True:
            from . import hesse as h
            imp.reload(h)    # necessary when working in ipython and changing code on the run
            self.surface = h.fqh_surface_interpolate(self.surface)
            if self._verbose:
                print("interpolate:")
                print("ss:",self.surface)
        #print "ska:",self.surface_interpolated
        self.temperatures = self.surface[:,0]
        self.surf_coefs_2nd = np.fliplr(self.surface[:,1:])
        self.surf_coefs_1st = np.fliplr(self.surface_1st[:,1:])

        self.surf_function = []
        for i,t in enumerate(self.temperatures):
            self.surf_function.append(np.poly1d(self.surf_coefs_2nd[i]))
        #print self.surface

        if self._oqh != False:
            if self.surface.shape[1] == 5:
                self.write_fqh_surface_3rd(surface = self.surface)
            if self.surface.shape[1] < 5:
                self.write_fqh_surface_2nd()


        return self.surface



    def import_surface_expr_sympy(self, surface_filename = None, surface = None, sysexit = True):
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

    def load_data(self, filename):
        ''' .. '''
        tmpdata = np.load(filename + ".npz")
        self.data = tmpdata['data']
        tmpdata.close()
        return


    def save_data(self, filename):
        ''' .. '''
        np.savez_compressed(filename, data=self.data)
        return

    def write_fqh_surface_1st(self, filename = None):
        ''' docstring here '''
        if filename == None:
            print('self._filessearchstring_common',self._filessearchstring_common)
            print('self._filessearchstring_suffix',self._filessearchstring_suffix)
            filename = self._filessearchstring_common+"Surface_1st_order_"+self._filessearchstring_suffix+"_e"+str(int(self._multiply_energy))+"_v"+str(int(self._multiply_volume))
        np.savetxt(filename, self.surface_1st,fmt="%.0f %.16f %.16f")
        _printgreen(filename+" written")
        return


    def write_fqh_surface_2nd(self, filename = None):
        ''' docstring here '''
        if filename == None:
            print('self._filessearchstring_common',self._filessearchstring_common)
            print('self._filessearchstring_suffix',self._filessearchstring_suffix)
            filename = self._filessearchstring_common+"Surface_2nd_order_"+self._filessearchstring_suffix+"_e"+str(int(self._multiply_energy))+"_v"+str(int(self._multiply_volume))
        np.savetxt(filename, self.surface,fmt="%.0f %.16f %.16f %.16f")
        _printgreen(filename+" written")
        return

    def write_fqh_surface_3rd(self, filename = None, surface = None):
        ''' docstring here '''
        if filename == None:
            filename = self._filessearchstring_common+"Surface_3rd_order_"+self._filessearchstring_suffix+"_e"+str(int(self._multiply_energy))+"_v"+str(int(self._multiply_volume))
        if surface == None:
            np.savetxt(filename, self.surf_out_3rd,fmt="%.0f %.12f %.12f %.12f %.12f")
        else:
            np.savetxt(filename, surface,fmt="%.0f %.16f %.16f %.16f %.16f")
        #np.savetxt(filename, self.surf_out_3rd)
        _printgreen(filename+" written")
        return


        t0 = self.temperatures[0]
        tm = self.temperatures[-1]
        x = self.volumes
        y = self.data[:,0]-self.data_fit[:,0]
        save = np.transpose([x,y])
        np.savetxt("delta_"+str(temp)+"K.dat",save)

    def plot_energy_vs_volume_qh(self, *tempin_all):
        fignum = 2
        import math
        round = 5
        roundnum = 10**round
        plt.ion();


        if type(self._plots) == bool: #then it is a number
            plt.clf();
        print("all:")

        for tempin in tempin_all:
            idx = min(list(range(len(self.temperatures))), \
                    key=lambda i: abs(self.temperatures[i]-tempin))

            temp = int(self.temperatures[idx])

            # labletext
            a0 = str(math.ceil(self.surf_coefs_2nd[idx][0]*roundnum)/roundnum)
            a1 = str(math.ceil(self.surf_coefs_2nd[idx][1]*roundnum)/roundnum)
            a2 = str(math.ceil(self.surf_coefs_2nd[idx][2]*roundnum)/roundnum)

            # labletext second label
            labeltext2 = ""
            for vind,i in enumerate(self.data[:,idx]):
                ene = str(math.ceil(i*roundnum)/roundnum)
                vol = str(math.ceil(self.volumes[vind]*roundnum)/roundnum)

                labeltext2 = labeltext2+"("+vol+"/"+ene+")"+"\n"

            if fignum == 1:
                plt.plot(   self.volumes,   self.data[:,idx],'o', \
                        label=labeltext2
                        )
                plt.plot(   self._fitgrid,  self._data_fit[:,idx],'-',\
                        label=str(temp)+" K: "+a0+"x^2 + "+a1+"x + "+a2
                            )
            if fignum == 2:
                x = self.volumes
                y = self.data[:,idx]-self.data_fit[:,idx]
                xmesh = self._fitgrid
                ymesh = self._data_fit[:,idx]

                ydata = self.data[:,idx]
                if type(self._plots) != bool: #then it is a number
                    fit_delta = np.transpose([x,y])
                    fit_delta = fit_delta[fit_delta[:,0].argsort()]
                    print("writing plot_fit_delta_"+str(temp)+"K.dat")
                    np.savetxt("plot_fit_delta_"+str(temp)+"K.dat",fit_delta)
                    print("writing plot_data_"+str(temp)+"K.dat")
                    np.savetxt("plot_data_"+str(temp)+"K.dat",np.transpose([x,ydata]))
                    print("writing plot_data_fit_"+str(temp)+"K.dat")
                    np.savetxt("plot_data_fit_"+str(temp)+"K.dat",np.transpose([xmesh,ymesh]))
                    self._plots = False # to be able to plot in ipython
                    return

                plt.plot(   self.volumes,   self.data[:,idx]-self.data_fit[:,idx],'o-', \
                        label=labeltext2+"  "+str(temp)+" K"
                        )
            # fit 3rd order
            if self.data.shape[0] > 3:
                plt.plot(   self._fitgrid,  self._data_fit_3rd[:,idx],'-',\
                        label=str(temp)+" K (fit 3rd order)")
            if self.data.shape[0] > 3:
                plt.plot(   self.volumes,   self.data_fit_3rd[:,idx],'.',)

        plt.xlim(self.volumes.min()-.1,self.volumes.max()+.1)
        plt.grid(True)  # Trunes on grid
        plt.title = 'Energy vs. volume quasiharmonic'
        leg  = plt.legend(loc='best', fancybox=True)
        #leg  = plt.legend(loc='best', fancybox=True, bbox_to_anchor=(1.05, 1))
        #leg  = plt.legend(loc=2, fancybox=True,mode="expand", bbox_to_anchor=(1.05, 1), borderaxespad=0.)

        # second label
        #labeltext2 = "v: "
        #for i in self.data[:,idx]:
        #    labeltext2 = labeltext2+" | "+str(math.ceil(i*roundnum)/roundnum)
        #leg2 = plt.legend([labeltext2],loc='best', fancybox=True)  # This removes the label before
        #from matplotlib.pyplot import gca
        #gca().add_artist(leg) # add l1 as a separate artist to the axes

        # make label fancy transparent
        print("self._plots:",self._plots)
        leg.get_frame().set_alpha(0.5)
        print("self._plots:",self._plots)
        return


    def plot_surface(self, fignum=1, points_in_plot=50, plotpoints=False):
        ''' hier muss rein welche data / datavol '''
        plt.ion()                                    # to free ipython shell
        fig = plt.figure(fignum)
        fig.clf()
        ax = fig.gca(projection='3d')               # to work in 3d
        plt.hold(True)

        temp, vol = np.meshgrid(self.temperatures, self.volumes)

        lines_n = len(self.temperatures) / points_in_plot
        if lines_n <= 1:
            lines_n = 1
        #print "lines_n:",lines_n

        # get every volume line seperate
        for ind, volume in enumerate(self.volumes):
            print("volume:", volume)

            x = np.repeat(self.volumes[ind],len(self.temperatures))
            y = self.temperatures
            z = self.data[ind]
            #ax.plot_surface(temp, vol, fqh, cmap=cm.hot);    # plot a 3d surface plot
            ax.plot(x, y, z, label='parametric curve')
            if plotpoints:
                ax.scatter(x, y, z);                        # plot a 3d scatter plot

        ax.set_xlabel('Volume [Ang^3]')
        ax.set_ylabel('Temperature [K]')
        ax.set_zlabel('Fqh')
        return


    def plot_delta_to_fit(self, draw_lines = 20):
        from mpl_toolkits.mplot3d import axes3d
        from matplotlib import cm

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        lines = int(len(self.temperatures)/draw_lines/2)
        print("l:",lines)
        data = self.data[:,::lines]
        temperatures = self.temperatures[::lines]
        data_fit = self.data_fit[:,::lines]
        volumes = self.volumes
        X = np.copy(data)
        Y = np.copy(data)
        Z = abs(data - data_fit)

        for i in np.arange(len(data)):
            X[i] = temperatures

        for i in np.arange(len(temperatures)):
            Y[:,i] = volumes

        #ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
        ax.plot_wireframe(X, Y, Z, rstride=3, cstride=3)
        #cset = ax.contour(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
        #cset = ax.contour(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
        #cset = ax.contour(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)

        ax.set_xlabel('Temperature (K)')
        ax.set_xlim(0, temperatures.max())
        ax.set_ylabel('Volume (Angstrom^3)')
        ax.set_ylim(volumes.min(), volumes.max())
        ax.set_zlabel('delta to 2nd order fit (meV)')
        ax.set_zlim(Z.min(), Z.max())
        plt.show()
        return

    def __str__(self):
        return self.data




if __name__ == '__main__':
    # %head Fqh_32at_cell_per_atom_Surface_1st_order__16.6315_16.8822_17.1354_17.3912_17.6494_17.9103_18.1737_18.4396_e1_v1
    # 0 114.425456223140 -4.550048118149
    # 1 114.425456223140 -4.550048118149
    # 2 114.425456223140 -4.550048118149
    # 3 114.425456223140 -4.550048118149
    # --> (ouput in mev/atom), example for al
    # 114.425456223140-4.550048118149*16.631514793360534 = 38.751263636
    # 114.425456223140-4.550048118149*18.43965773972576  = 30.524126225

    p = qh().help()
    args = p.parse_args()
    my.create_READMEtxt(os.getcwd())
    qh = qh(args)

    if args.quickanalysis == 1:
        Fqh_surf_1 = "Fqh_32at_cell_per_atom_Surface_1st_order__16.6315_16.8822_17.1354_17.3912_17.6494_17.9103_18.1737_18.4396_e1_v1"
        Fqh_surf_1 = "Fqh_81at_cell_per_atom_Surface_1st_order__16.1854_16.4293_16.6758_16.9247_17.1760_17.4299_17.6862_17.9450_e1_v1"
        volumes = np.arange(16.0,17.3,0.01)
        fqh_ene = np.zeros(len(volumes))
        evinet = np.loadtxt('../evinet/EVinet_1')
        import feos
        vi = feos.eos()
        T = 0
        vi.e0    = evinet[0]
        vi.v0    = evinet[1]
        vi.b0    = evinet[2]
        vi.b0der = evinet[3]
        vi.parameters = [vi.e0, vi.v0, vi.b0, vi.b0der]
        evinet_ene = feos.vinet(volumes, *vi.parameters)
        print('evient_ene',evinet_ene)
        np.savetxt("evinet_ene_"+str(T)+"K",np.transpose([volumes,evinet_ene]))

        for T in [0,100,300,500]:
            fqh = np.loadtxt(Fqh_surf_1)
            print(fqh)
            #Temperature = idx = fqh[:,0]
            for vidx,v in enumerate(volumes):
                fqh_ene[vidx] = fqh[T,1] + fqh[T,2]*v
            np.savetxt("fqh_ene_"+str(T)+"K",np.transpose([volumes,fqh_ene]))
            np.savetxt("sum_ene_"+str(T)+"K",np.transpose([volumes,fqh_ene+evinet_ene]))
        sys.exit()



    if args.ih != False:
        print(args.ih)
        from . import hesse
        out=hesse.fqh_interpolate(fqh=np.loadtxt(args.ih))
        np.savetxt("out.fitted",out)
        print("saved out.iftted")
        sys.exit()

    if args.iqh != False:
        print("importing",args.iqh)
        qh.import_surface(surface_filename = args.iqh, interpolate = False, sysexit = True)

    # import files
    if args.iqh == False:
        print("importing Fqh files")
        qh.import_fqh_files()   #(files="Fqh_fromExact[Ff]reqs_[0-9.]*")
        qh.fit_surface()

    if args.ov != False or args.ofcca != False:
        qh.helmholtz_free_energy(v = args.ov, fccalat=args.ofcca, save = args.ob)

