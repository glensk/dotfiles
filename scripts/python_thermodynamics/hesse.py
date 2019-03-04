#!/usr/bin/env python

# /Users/glensk/Dropbox/proj/current_parabola_to_morse/Al_displacements_2x2x2sc_quer/4.04Ang_0.3_quer_wirklich
# hesse.py al -ene -p m --potparam 0.300 1.478 2.856
# hesse.py al -ene -p mc1 --potparam 0.283 1.48 2.856 0.749 0.372
# ~/Thermodynamics/python_thermodynamics/hesse.py al -eneext -v
# /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2015.08_phonon_lifetimes_2_for_munich/calculations/allq_Al_2_4.13_1000000steps_timestep0.001_totaltime10000ps_900K_dump_10_test_manual_short/check_one_disp/coords_pyton_forces_1

#from . import my_atom as atom
import my_atom as atom
import os
import sys
import math
import numpy as np
import argparse
from argparse import ArgumentDefaultsHelpFormatter
import textwrap
import shutil
# temporarily outcommented .... lets see weather this is a good idea
#try:
#    from . import crystal_generator
#try:
#    from . import utils
sys.dont_write_bytecode = True   # Important: We dont want bytecode, this could possibly save "old" parameters wich are not the current onse (but is probably faster)
np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=6)    # print only 6 digist after .
import numpy.core.arrayprint as arrayprint
import contextlib

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

# functions with _name are not seen when imported
def _printred(var):
    ''' help text '''
    ENDC = '\033[0m'
    red = '\033[31m'
    print(red + str(var) + ENDC)


def _printgreen(var):
    '''
    print green to screen
    '''
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

def fvib_classical(T, freqs):
    '''
    help text
    '''
    kB = 8.6173423*10**(-2)
    return kB * T * np.log(1 - np.exp(-(freqs/(kB*T)) ))

def fit_fcc_alat_mean_freqs():
    ''' return second order polynomial fit for omegas(freqs) vs volume '''
    data = np.loadtxt("mean_freqs")
    data[:,0] = data[:,0]**3./4.  # here we assume fcc alat
    np.savetxt("mean_freqsx",data)

    coefs_freqs_vs_volume = np.polyfit(data[:,0], data[:,1], 2)

    def quadfunc(x,a,b,c):
        #return a+b*x*c*x**2  NEVER LIKE THIS, scipy.optimize needs : c*x**2+b*x+a
        return c*x**2+b*x+a
    from scipy import optimize
    coefs_freqs_vs_volume_, covariance = \
            optimize.curve_fit(quadfunc,data[:,0],data[:,1])

    def meanFreqsFit(V):
        ''' return a Frequence for a certain volume '''
        return np.poly1d(coefs_freqs_vs_volume)(V)
    #funct_freqs_vs_volume = np.poly1d(coefs_freqs_vs_volume) # as a function of volme
    t=934
    volmesh = np.linspace(data[:,0].min(),data[:,0].max(),num=100)
    #print "volmesh:",volmesh
    freqs_fitted = np.zeros((len(volmesh),2))
    freqs_fitted_ = np.zeros((len(volmesh),2))
    #print "freqs_fitted:",freqs_fitted
    freqs_fitted[:,0] = volmesh
    freqs_fitted_[:,0] = volmesh
    for ind,vol in enumerate(volmesh):
        #freqs_fitted[ind,1] = funct_freqs_vs_volume(vol)
        #freqs_fitted_[ind,1] = quadfunc(vol, *coefs_freqs_vs_volume_)

        # bei einem bestimmtem volume eine bestimmte frequenz VVVVVV
        #freqs_fitted[ind,1] = fvib_classical(t,funct_freqs_vs_volume(vol))

        # was ich hier brauche ist ein 3d fit, zurzeit mache ich alles bei einer fixen temperatur
        freqs_fitted[ind,1] = fvib_classical(t, meanFreqsFit(vol))
        freqs_fitted_[ind,1] = fvib_classical(t,quadfunc(vol,*coefs_freqs_vs_volume_))
        #print "t: ",t," vol: ",vol," freqs_fit: ",freqs_fitted[ind,1]
    #print ">>>",quadfunc(volmesh, *coefs_freqs_vs_volume_)
    #print ">>>",freqs_fitted
    np.savetxt("FvibSave.dat",freqs_fitted)
    np.savetxt("FvibSave_.dat",freqs_fitted_)
    print(volmesh)
    print("t:",t)
    freq1 = data[0][1]
    freq2 = data[1][1]
    freq3 = data[2][1]
    print("freq1:",freq1)
    print("freq2:",freq2)
    print("freq3:",freq3)
    print("out:",fvib_classical(t,freq1))
    print("out:",fvib_classical(t,freq2))
    print("out:",fvib_classical(t,freq3))
    return coefs_freqs_vs_volume

def fqh_cell(T, ExactFreqs = None):
    '''
    harmonic free energy for whole supercell
    '''
    kB=0.086173423
    return np.sum(ExactFreqs/2+(kB*T*np.log(1-np.exp(-ExactFreqs/(kB * T)))) )

def fqh_atom(T, ExactFreqs = None):
    '''
    harmonic free energy per atom, how many ExactFreqs should we have here?:
        It should be 3N-3! (The three 0 freqs shoud be excluded) ?'''
    return fqh_cell(T = T, ExactFreqs = ExactFreqs)/(len(ExactFreqs)/3.)

def fqh_fit(fqh = None, maxfev = 16000, maxdiffallowed = 0.001, verbose = False, print_diff = True):
    """ this script fits a fqh (as a function of T); it can work well (Al) but also not
    so well (Si), see correspondig print_diff"""
    #fqh = np.loadtxt("Fqh_fromMeshFreqs_4.12")
    if verbose:
        print("fqh:",fqh)
    if fqh == None:
        error = "This script needs fqh as input: \
               [[1, 28.8], \
                [4, 26.6], \
                [7, 24.4], \
                 ..., \
                [934, -654.4]]"
        sys.exit(error)
    guess = [ 1., 2.5,
            29., 28., 27., 26., 25., 24., 23., 22., 21., 20., 15.,
            31.]
    def fqh_fit_fun(T, a, scaling,
            omega1, omega2, omega3, omega4, omega5, omega6, omega7, omega8, omega9, omega10, omega11,
            omegal):
        kB=0.086173423
        return (a +
               omega1/2+(kB*T*np.log(1-np.exp(-omega1/(kB * T)))) + \
               omega2/2+(kB*T*np.log(1-np.exp(-omega2/(kB * T)))) + \
               omega3/2+(kB*T*np.log(1-np.exp(-omega3/(kB * T)))) + \
               omega4/2+(kB*T*np.log(1-np.exp(-omega4/(kB * T)))) + \
               omega5/2+(kB*T*np.log(1-np.exp(-omega5/(kB * T)))) + \
               omega6/2+(kB*T*np.log(1-np.exp(-omega6/(kB * T)))) + \
               omega7/2+(kB*T*np.log(1-np.exp(-omega7/(kB * T)))) + \
               omega8/2+(kB*T*np.log(1-np.exp(-omega8/(kB * T)))) + \
               omega9/2+(kB*T*np.log(1-np.exp(-omega9/(kB * T)))) + \
               omega10/2+(kB*T*np.log(1-np.exp(-omega10/(kB * T)))) + \
               omega11/2+(kB*T*np.log(1-np.exp(-omega11/(kB * T)))) + \
               omegal/2 +(kB*T*np.log(1-np.exp(-omegal/( kB * T)))))/scaling
    from scipy import optimize
    print(fqh[:,0].min())
    print(fqh[:,0].max())
    #freqs, covar = optimize.curve_fit(fqh_fit_fun, fqh[:,0], fqh[:,1], guess, maxfev = maxfev)
    freqs, covar = optimize.curve_fit(fqh_fit_fun, fqh[:,0], fqh[:,1], guess, maxfev = maxfev)
    if verbose:
        print("freqs:",freqs)
    fqh_fit = fqh_fit_fun(fqh[:,0], *freqs)
    fqh_diff = fqh_fit - fqh[:,1]
    if verbose:
        print("fqh_fit:",fqh_fit)
        print("fqh_diff:",fqh_diff)
    maxdiff = abs(fqh_diff).max()
    if maxdiff > maxdiffallowed:
        sys.exit("maximal difference is too large: "+str(maxdiff))
    tmp_max = fqh[:,0][-1]
    print("tmp_max:",tmp_max)
    all_tmp = np.arange(0.,int(tmp_max)+1)
    all_fqh = fqh_fit_fun(all_tmp, *freqs)
    all_fqh_out = np.empty((len(all_tmp),2))
    all_fqh_out[:,0] = all_tmp
    all_fqh_out[:,1] = all_fqh
    if print_diff:
        return fqh_diff
    return all_fqh_out

def fqh_interpolate(fqh = None):
    import copy
    write_out = copy.deepcopy(fqh)
    print("fqh1:",fqh,type(fqh))
    print("write_out1:",write_out,type(write_out))
    if type(fqh) == str:
        if os.path.isfile(fqh) != True:
            sys.exit(fqh+" does not exist!")
        fqh=np.loadtxt(fqh)
        print("write_out2:",write_out,type(write_out))

    if fqh[0][0] != 0:
        if fqh[0][0] <= 5:
            for i in np.arange(0.,fqh[0][0])[::-1]:
                insert = [i,fqh[0][1]]
                fqh = np.insert(fqh,0,insert,axis=0)
        else:
            sys.exit("no energies up to 5 K")

    # spline fit
    from scipy.interpolate import interp1d
    x = fqh[:,0]
    y = fqh[:,1]
    tempall = np.arange(0.,x.max()+1)
    f = interp1d(x, y)
    f2 = interp1d(x, y, kind='cubic')
    print("tempall:",tempall)

    # otput
    out_f = np.empty((len(tempall),2))
    out_f[:,0] = tempall
    out_f[:,1] = f(tempall)

    out_f2 = np.empty((len(tempall),2))
    out_f2[:,0] = tempall
    out_f2[:,1] = f2(tempall)

    out_f2_diff = f2(fqh[:,0])-fqh[:,1]
    print("write_out2:",write_out,type(write_out))
    print("fqh:",fqh)
    if type(write_out) != False:
        np.savetxt(write_out+"_interpolated",out_f2,fmt='%.1f %.17f')
    return out_f2

def fqh_surface_interpolate(fqh = None):
    if fqh == None:
        sys.exit("need fqh surface as input")
    temps = np.arange(fqh[-1][0]+1)
    out = np.zeros((len(temps),fqh.shape[1]))
    out[:,0] = temps
    if np.array_equal(temps,fqh[:,0]):
        return fqh
    print("interpolating surface ...")
    print("temps, is    :",fqh[:,0][:7],"...",fqh[:,0][-7:])
    print("temps, should:",temps[:7],"...",temps[-7:])
    for col in np.arange(fqh.shape[1]-1):   # -1 since first column is the temperature
        #out[:,col+1] = fqh_interpolate(fqh[:,[0,col+1]])
        #print "------------------"
        #print fqh_interpolate(np.array(fqh[:,[0,col+1]]))
        out[:,col+1] = fqh_interpolate(np.array(fqh[:,[0,col+1]]))[:,-1]
        #print "------------------"
    return out

def read_Hessematrix(hessematrixfile = "HesseMatrix_sphinx" ):
    hessematrix = np.loadtxt(hessematrixfile)*97.173617
    #You have:hartree/bohrradius^2
    #You want:eV/angstrom^2
    #    *97.173617
    #    /0.010290859
    return hessematrix

def qh_forces(dpos = None, h = None ):
    '''
    dpos = deltas of positions in cartesian coords
    h = hessematrix
    '''
    #print dpos
    #print '---'
    #print h
    #print '---'
    f = np.dot(-h, np.ndarray.flatten(dpos))
    return f.reshape(dpos.shape)

def qh_energy_cell(dpos = None, h = None ):
    f = qh_forces( dpos = dpos, h = h )
    e = -np.dot(np.ndarray.flatten(dpos),np.ndarray.flatten(f))/2.
    if e < 1e-8 and e > -1e-8:
        e = 0
    return e

def qh_energy_atom(dpos = None, h = None ):
    ecell = qh_energy_cell( dpos = dpos, h = h )
    #print "ecell:",ecell
    atoms = len(np.ndarray.flatten(dpos))/3.-1.
    #print "atoms",atoms
    e = ecell/atoms*1000.
    if e < 1e-8 and e > -1e-8:
        e = 0
    return e

def plot_morse(parameters):
    '''
    help text
    '''
    #fitx = np.linspace(datax.min(),datax.max(), num=100)
    #fity = Morse_derivative(fitx, parameters[0], parameters[1], NN)
    fitx = np.linspace(2.0,5.0, num=100)
    fity = Morse_derivative(fitx, parameters[0], parameters[1],parameters[2])
    np.savetxt(
            "fit_morse_"+str(parameters[0])[:5]
            +"_"+str(parameters[1])[:5]
            +"_"+str(parameters[2])[:5],
                np.transpose([fitx, fity]),
                fmt='%.18f',
                delimiter='  ')   # X is an array

def get_hessematrix_lj_morse(
        pot = False,   # choices h (hesse), l (lennard jones), m (Morse)
        potparam = False,  # only used in case of pot == "l" or pot == "m"

        coordfile0_cart = False,
        coord0_cart = False,
        coordfile0_rrel = False,
        coord0_rrel = False,

        cellfile = False,
        cell = False,
        disp=0.01):
    if pot == "l" or pot == "m":
        pass
    else:
        sys.exit("you can only create a Hessian Matrix for Lennard Jones or Morse")

    if type(coordfile0_cart) == bool and type(coord0_cart) == bool \
        and type(coordfile0_rrel) == bool and type(coord0_rrel) == bool:
            if os.path.isfile("EqCoords_direct") == True:
                coordfile0_rrel = "EqCoords_direct"
    if type(coordfile0_cart) == bool and type(coord0_cart) == bool \
        and type(coordfile0_rrel) == bool and type(coord0_rrel) == bool:
            sys.exit("you need the undisplaced structure \
                    (e.gl. EqCoords_direct) as reference to get 1NN atoms (1)")
    if type(cell) == bool and type(cellfile) == bool:
        if os.path.isfile("cell") == True:
            cellfile = "cell"
    if type(cell) == bool and type(cellfile) == bool:
        if os.path.isfile("POSCAR") == True:
            utils.run2("rm -f cell; POSCAR_cell_cartesian.sh > cell") # to create cell
        if os.path.isfile("cell") == True:
            cellfile = "cell"
    if type(cell) == bool and type(cellfile) == bool:
        sys.exit("you need the cell file")

    crystal0 = crystal_generator.crystal()
    crystal0.load_positions_cell(coordfile_cart = coordfile0_cart, coord_cart = coord0_cart,
            cellfile = cellfile, cell = cell, coordfile_rrel = coordfile0_rrel,
            coord_rrel = coord0_rrel)

    pos0 = np.copy(crystal0.rcar)
    hessematrix=np.zeros((pos0.shape[0]*3,pos0.shape[0]*3))
    stradd=""
    for pppidx,ppp in enumerate(potparam):
        if pppidx == 0:
            stradd = stradd+str(potparam[pppidx])
        else:
            stradd = stradd+"_"+str(potparam[pppidx])
    hessefilename =  "HesseMatrix_"+str(crystal0.cellvec[0,0])+"Ang_"+pot+"_"+stradd
    if os.path.isfile(hessefilename) == True:
        sys.exit(hessefilename+" does already exist")
    print("hessefilename:",hessefilename)
    a = crystal0.cellvec[0]
    b = crystal0.cellvec[1]
    c = crystal0.cellvec[2]
    volume = np.dot(a,np.cross(b,c))
    if volume <= 1.0:
        sys.exit("volume of cell is "+str(volume)+"!")
    print("volume:",volume)


    for iidx,i in enumerate(pos0):   # schleife ueber atome, 1..32
        print(iidx,"/",pos0.shape[0]) # schleife ueber xyz 1..3
        for jidx,j in enumerate(i):
            pos1 = np.copy(pos0)
            pos1[iidx,jidx] = pos1[iidx,jidx]+disp
            eahmev,eah,fah = get_energy_forces(
                pot=pot,
                potparam=potparam,
                coord_cart = pos1,
                coord0_cart = pos0,
                cell=crystal0.cellvec,
                printresult=False)
            hessematrix[iidx*3+jidx] = fah.reshape((1,pos0.shape[0]*3))/(-disp)
    np.savetxt(hessefilename,hessematrix/97.173617,fmt="%.13f")
    return hessematrix


##################################################################
# help
##################################################################
def help(p = None):
    string = '''
    examples:
    hesse.py al -i HesseMatrix_4.13 -wf (defines eigenfrequencies and Free Energy from HesseMatrix)
    hesse.py al -i HesseMatrix_4.13 -p h -ene (calculates energy and forces from Hessematrix and cartesian_coords)
    hesse.py al -eneext -v  (calculates forces and energies from external calculate_enery_forces file)
    hesse.py al -i FORCE_CONSTANTS -p h -v -wh  (to create HesseMatrix_sphinx_al)
    hesse.py al -i HesseMatrix_sphinx_al -p h -v -ene (to get the energy)



    curently not working:
    hesse.py al -ene -p m --potparam 0.279300 1.432673 2.920351
    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    p.add_argument('elements', metavar='elements', nargs='+',
    #p.add_argument('elements', metavar='elements', nargs='?',
        help='list of elements, e.g. Mg; AL cU; si 2 u 7; u o 8 6')
    p.add_argument('-i',   '--inputfile',
       help='specify name of inputfile', type=str, default=False) #"Hessematrix_sphinx")
    p.add_argument('-i1nn',   '--inputfile1nn',
       help='specify name of inputfile of HesseMatrix with 1stNN',
       type=str, default=False) #"Hessematrix_sphinx_1NN")
    p.add_argument('-l', action='store_true',
       help='print output only to temperatrue of lowest melting \
       elemetn if several exist', default=False)
    p.add_argument('--fqh_interpolate', action='store_true', default=False,
        help='interpolate file which must be given with -i argument (e.g. --fqh_interpolate -i Gibbs_energy')

    # print to screen
    p.add_argument('-f', action='store_true',
       help='print ExactFreqs', default=False)
    p.add_argument('-fm', action='store_true',
       help='print mean_ExactFreqs', default=False)
    p.add_argument('-fitmf', action='store_true',
       help='fit fcc alat mean freqs', default=False)

    # write
    p.add_argument('-wf', nargs='?', type=str, const="ExactFreqs",
        help='write ExactFreqs file (default: ExactFreqs)') #,
    p.add_argument('-e', action='store_true',
       help='print Fqh_fromExactFreqs', default=False)
    p.add_argument('-we', nargs='?', type=str, const="Fqh_fromExactFreqs",
        help='write Fqh_fromExactFreqs file (default: ExactFreqs)') #,
    p.add_argument('-wh', nargs='?', type=str, const="HesseMatrix_sphinx",
        help='write HesseMatrix_sphinx file (default: HeseMatrix_sphinx)') #,
    p.add_argument('-w', action='store_true',default=False,
        help='write all output files above')



    #############################################
    # morse lj
    #############################################
    p.add_argument('-p', '--pot',choices=['h', 'm', 'l', 'mc1', 'i', '135', '1357', '13579',
    '1357911', 'extern', 'externorig' ], default=False, help=textwrap.dedent('''\
            define potential
            h=hesse (needs Hessematrix with -i option)
            m=Morse
            l=Lennard-Jones
            mc1=Morsec2
            i=inversepotential (Dario Alfe))'''))
    p.add_argument('-p2', '--pot2', default=False, nargs='+',
            help='define potential 2 for positive (attractive) longitudinal \
                    part e.g.: mc1 0.447778 1.051045 2.920351 0.684528 -0.24431 \
                    (parameters defined in potparam go to repulsive part')
    p.add_argument('-pp', '--potparam',nargs='+', default=False, help=textwrap.dedent('''\
            define parameter for potential (is to be used with -p option)
            e.g. hesse.py al -ene -p m --potparam 0.279300 1.432673 2.920351   (4.13/np.sqrt(2) = 2.9203510063004412)
            e.g. hesse.py al -ene -p m --potparam 0.300 1.478 2.856
            e.g. hesse.py al -ene -p mc1 --potparam 0.283 1.48 2.856 0.749 0.372'''))
    p.add_argument('-ene', action='store_true', default=False,
        help='calculate energy and forces for given pot')
    p.add_argument('-eneext', action='store_true', default=False,
        help='calculate energy and forces by using current calculate energy and forces (sets pot to -extern')
    p.add_argument('-eneexto', action='store_true', default=False,
        help='calculate energy and forces by using calculate energy and forces from ~/Thermodynamics/python_thermodynamics foler (sets pot to -externorig')
    p.add_argument('-tm', action='store_true', default=False,
        help='get transitional forces of 1NN from hessematrix which must be available')
    p.add_argument('-hesse1nn', nargs='+', default=False,
        help='create Hessematrix with only first nearest neighbor atoms')
    p.add_argument('-hrest1', nargs='+', default=False,
        help='calculate energy and forces from HesseMatrix from all but 1NN atoms (specify the Hessematrix and Hessematrix1NN as parameters')
    p.add_argument('-getforcepot', nargs='+', default=False,
        help='calculate forces vector for a given pot and potparam, the forcesvector is given as 3 numbers')




    p.add_argument('-plotm',nargs='+', default=False, help='plot morse derivatvie for given parameters')

    p.add_argument('--fitpotf',
            help='fit forces to morse, 1st arg: filename; 2nd arg: NN dist; define pot with -p option', nargs='+')
    p.add_argument('--polyfit',
            help='fit forces to polynomial first argument defines order', nargs='+')
    p.add_argument('--ffd','-ffd',nargs='+',default=False,
            help='get forces from displacements')
            #help='get forces from displacements', action='store_true', default=False)
    p.add_argument('-morseh', nargs='+',
       help='create HesseMatrix with morse parameters')
    p.add_argument('-ljf',
       help='fit forces to lennard-jones', type=str)
    p.add_argument('-ljh', nargs='+',
       help='create HesseMatrix with lj parameters')
    #p.add_argument('-fitfqh', nargs='?', type=str, const="Fqh_fromExactFreqs",
    #    help='fit Fqh file') #,

    # verbose help
    p.add_argument('-v','--verbose',
            help='verbose', action='store_true', default=False)

    return p


class hesseclass( object ):
    '''
    defines everything related to the harmonic approximation.
    defines eigenfrequencies and Free Energy from HesseMatrix
    '''
    def __init__( self, args = False , listin = False, H = False, show_negative_eigenvalues = True):
        '''
        units of h (hessematrix): [eV/Angstrom^2]
        if HesseMatrix is imported units are expected in [hartree/bohrradius^2]
        '''
        #def __init__( self, listin = None, H = None, verbose = None ):
        # im __init__ sollte nie ein systemabbrch erfolgen (keine checks)
        # im __init__ darf ein systemabbruch erfolgen wenn infiputvars in gesetz werde koennen
        # help really should be external

        self.H = H
        self.M = None
        self.freqs = False
        self.freqsunsorted = False
        self.ene_atom = False
        self.ene_cell = None

        self._verbose = None
        self.__verbose = None
        self.show_negative_eigenvalues = show_negative_eigenvalues
        self.listin = listin      # listin   --> 'Al', 'Si1Al31'
        self.inputfile = "Hessematrix_sphinx"
        self._inputfile_default = "Hessematrix_sphinx"
        self.inputfile_all = None
        self._filename_addstring = ""
        self._l = False
        self._fm = False
        self.inputfile = False
        self._H_filename_units_possible = ["hartree/bohrradius^2","eV/angstrom^2"]
        self._H_filename_units = "hartree/bohrradius^2"
        self.writeoutput = False

        if args:
            self._verbose = args.verbose
            self.inputfile = args.inputfile
            self.listin = args.elements
            self._l = args.l
            self._fm = args.fm
            self.writeoutput = True  # when called from shell

        if self._verbose:
            print("args:",args)

        #############################################################
        ## listin defines the atoms which are used with the hessematrix
        #############################################################
        # exit if necessary
        m1="Call hesse class with"
        m2="    hesse=h.hesseclass(['al'],              h=\"filename\")"
        m3="or  hesse=h.hesseclass(['al', 'si', 31, 1], h=np.array([...]))"
        if self.listin == False:
            sys.exit(m1+"\n"+m2+"\n"+m3)
        #print('listin',self.listin)
        self.atomdata = atom.atom(self.listin)
        #print('ad',self.atomdata)
        self.mass = self.atomdata.mass
        #print('mass',self.mass)
        self._temp = int((self.atomdata.melting+1.0).max())  # +1 since int rounds down
        if self._l:
            self._temp = int((self.atomdata.melting+1.0).min())  # +1 since int rounds down
        if self._verbose:
            print("self.mass        :",self.mass)
            print("self._temp       :",self._temp)
            print("self.atomdata    :",self.atomdata)


        ####################################################################
        ### check if inputfile is ok, if only one Hesse... exists, take it
        ####################################################################
        if type(self.H) == bool:
            self.import_H()

        self.atoms = self.H.shape[0]/3

        ############################################################
        # self.m -- Mass Matrix
        ############################################################
        self.get_mass_matrix()
        #print('self.M',self.M)


        ############################################################
        # self.freqs -- array
        ############################################################
        try:
            self.freqs = self.get_freqs()
        except: # catch *all* exceptions
            self.freqs = sys.exc_info()
        if self.__verbose:
            print("self.freqs:",self.freqs)
        if self._fm:
            print(self.freqsmean)


        ############################################################
        # self.freene -- array
        ############################################################
        try:
            self.ene_atom = self.get_ene_atom()
        except: # catch *all* exceptions
            self.ene_atom = sys.exc_info()

        try:
            self.ene_cell = self.get_ene_cell()
        except: # catch *all* exceptions
            self.ene_cell = sys.exc_info()

        if self._verbose:
            print('self.write() now ....')

        if self.writeoutput:
            self.write()
        return


    def import_H(self):
        print("self.inputfile:",self.inputfile)
        searchhessefile = False
        if type(self.inputfile) == bool:
            searchhessefile = True
        else:
            if os.path.isfile(self.inputfile) != True:
                searchhessefile = True

        #if os.path.isfile(self.inputfile) != True:
        print("searchhessefile",searchhessefile,"(False is good)")
        if searchhessefile == True:
            #_printred("Inputfile "+str(self.inputfile)+" does not exist, looking for Hesse* files")
            import glob
            all_h_files = glob.glob("Hesse*")
            if len(all_h_files) == 1:
                self.inputfile = all_h_files[0]
            else:
                print("no Hesse* file found ...")
                if os.path.isfile("FORCE_CONSTANTS"):
                    self.inputfile = "FORCE_CONSTANTS"
                elif len(all_h_files) >= 0:
                    m1='Found '+str(len(all_h_files))+' Hesse* files: '
                    m2=', '.join(all_h_files)
                    m3='. Specify using the -i option which one to use.'
                    self.inputfile_all = all_h_files
                    #sys.exit(m1+m2+m3)
                else:
                    m1='Inputfile \"'+args.inputfile+"\" does not exist. "
                    m2='Please define an inputfile using -i FILENAME'
                    sys.exit(m1+m2)
        if self._verbose:
            print("self.inputfile",self.inputfile)
            print("self.inputfile_all",self.inputfile_all)
        #for counter in self.inputfile_all:
        self.H = None
        #self.inputfile = counter
        print("self.inputfile:",self.inputfile)
        if self.inputfile == False:
            sys.exit("no inputfile found")


        #############################################################
        ## self._filename_addstring
        #############################################################
        if self.inputfile != self._inputfile_default:
            self._filename_addstring = _string_to_num_list(self.inputfile, tostring = True)

        if self._verbose:
            print("################################################")
            print("self._filename_addstring1 :",self._filename_addstring)
            print(".join(self.listin)      1 :","_".join(self.listin))
        self._filename_addstring = "_".join(self._filename_addstring)
        if len(self._filename_addstring) > 1:
            self._filename_addstring = ""
            if self._verbose:
                print("self._filename_addstring1 :",self._filename_addstring)
        if self._filename_addstring == "":
            self._filename_addstring = ""
            if self._verbose:
                print("self._filename_addstring2 :",self._filename_addstring)
                print(".join(self.listin)      2 :","_".join(self.listin))

        #if "_".join(self.listin) == self._filename_addstring:
        #    self._filename_addstring = ""
        else:
            self._filename_addstring = "_"+self._filename_addstring
            if self._verbose:
                print("self._filename_addstring3 :",self._filename_addstring)
                print(".join(self.listin)      3 :","_".join(self.listin))
        # ################################################
        # if -i = HesseMatrix_sphinx_al_3.96
        # self._filename_addstring1 : ['3.96']
        # .join(self.listin)      1 : al
        # self._filename_addstring3 : _3.96
        # .join(self.listin)      3 : al
        # ################################################
        self._filename_addstring = "_"+ "_".join(self.listin) + self._filename_addstring


        #self._filename_addstring = "_"+ "_".join(self.listin)
        if self._verbose:
            print("################################################")
            print(".join(self.listin)       :","_".join(self.listin))
            print("self._filename_addstring :",self._filename_addstring)
            print("################################################")


        if self.inputfile == None: # --> need to import hessematrix
            sys.exit(m1+"\n"+m2+"\n"+m3)
        if type(self.inputfile) == np.ndarray:
            self.H = self.inputfile

        if type(self.inputfile) == str:
           #print "six:",self.inputfile, self.H
           # check weather file exists
           if os.path.isfile(self.inputfile) != True:
               sys.exit('Inputfile \"'+self.inputfile+"\" does not exist.")

           if os.path.isfile(self.inputfile) == True:
               # check for FORCE_CONSTANTS (phonopy)
               from itertools import islice
               with open(self.inputfile) as myfile:
                   head=list(islice(myfile,3))
                   if len(_string_to_num_list(head[0])) == 1 and \
                   len(_string_to_num_list(head[1])) == 2 and \
                   len(_string_to_num_list(head[2])) == 3:
                       self._H_filename_units = "eV/angstrom^2"
                       self.H = self.phonopy_FORCE_CONSTANTS_to_h(filename = self.inputfile)
                       if self._verbose:
                            print("SH:",self.H,type(self.H))
               #print("kkk")
               #print("kkk",self.H)
               # otherwise Hessematrix
               #if type(self.H) == NoneType or type(self.H) == bool:
               if self.H == None or type(self.H) == None or type(self.H) == bool:
                   if self._H_filename_units == "hartree/bohrradius^2":
                       self.H = np.loadtxt(self.inputfile)*97.173617
                   if self._H_filename_units == "eV/angstrom^2":
                       self.H = np.loadtxt(self.inputfile)


        ############################################################
        # check length of self.mass vs dimensions of h
        ############################################################
        #print("HHH",self.H)
        #if self.H == None:
        #    sys.exit('11')
        if len(self.H.shape) != 2:
            sys.exit("Your Hessematirxi is not a 2D array")
        if self.H.shape[0] != self.H.shape[1]:
            sys.exit("Your Hessematirxi is not a square matrix (same number of rows and columns)")
        return


    def get_mass_matrix(self):
        if len(self.mass) == 1:
            #print('sm',self.mass)
            #print('sm',self.mass*self.mass)
            #print('sa',self.atoms)
            #print('sa',self.atoms*3)
            self.M = np.tile(self.mass*self.mass,(int(self.atoms*3),int(self.atoms*3)))
            if self._verbose:
                print("self.M.shape:",self.M.shape)

        elif len(self.mass) == self.atoms:
            self.M = np.zeros((self.atoms*3,self.atoms*3))
            for zeile,mass1 in enumerate(self.mass):
                for spalte,mass2 in enumerate(self.mass):
                    self.M[zeile*3,spalte*3:spalte*3+3] = mass1*mass2
                    self.M[zeile*3+1,spalte*3:spalte*3+3] = mass1*mass2
                    self.M[zeile*3+2,spalte*3:spalte*3+3] = mass1*mass2
        else: # len(self.mass) != 1:
            if len(self.mass) != self.atoms:
                m1="Your Hessematrix and amount of atoms have different dimensions. \n"
                m2="Hessematrix: "+str(self.H.shape)+" (="+str(self.atoms)+" atoms)\n"
                m3="atoms in: "+str(self.atomdata.symbol)+" (having mass:"+str(self.mass)+")"
                sys.exit(m1+m2+m3)

        if self._verbose:
            print("self.atoms:",self.atoms)
            print("self.mass:",self.mass)
            print("len(self.mass):",len(self.mass))
        if self.__verbose:
            print("M:",self.M)
        return

    def write(self):
        ############################################################
        # write stuff out to file
        ############################################################
        if args.f:
            for i in h.freqs:
                print(i)
            #sys.exit(0)

        # -wf
        if args.wf or args.w:
            self.write_freqs(args.wf)

        # -e
        #if args.e:
        #    for i in h.ene_atom:
        #        print i[0],i[1]

        # -we
        if args.we or args.w:
            self.write_ene_atom(args.we)

        # -wh
        #print('args.wh, argsw',args.wh,args.w)

        if args.wh or args.w:
            print('writing ....')
            self.write_hessematrix(args.wh)
            #np.savetxt(args.wh, self.H/97.173617)

        # -fm
        if args.fm:
            #np.savetxt(filename + self._filename_addstring, self.freqs, fmt='%.12f')
            filename = 'ExactFreqs'
            if self._verbose:
                print("filename                 :",filename)
                print("self._filename_addstring :",self._filename_addstring)
                print("self.freqs:",self.freqs)
            np.savetxt(filename + self._filename_addstring+"_mean",
                    np.array([self.freqs.mean()]), fmt='%.12f')
            print(filename + self._filename_addstring+"_mean written")

            #self.write_hessematrix(args.wh)
            #np.savetxt(args.wh, self.H/97.173617)

        ############################################################
        #
        ############################################################
        def read_POSITIONS(self):  ## -> get std
            pass

        self.structure_undisplaced = None
        self.structure_displaced = None


    def read_Hessematrix( self , hessematrixfile = "HesseMatrix_sphinx" ):
        # FORCES in VASP are in eV/Angstrom
        # Hessematrix is in eV/Angstrom^2 == FORCES/displacemeners
        # /glensk/Dropbox/proj/current_cu_paper/version4_Cu_PRL/ForResubmission_PSS_4/'
        #You have:hartree/bohrradius^2
        #You want:eV/angstrom^2
        #    *97.173617
        #    /0.010290859
        self.H = np.loadtxt(hessematrixfile)*97.173617
        # push calc of eigenfreqs
        return self.H


    def H_plot(self, temp = 0.05):
        ''' MatrixPlot of the HesseMatirxi '''
        import matplotlib
        import matplotlib.pyplot as plt
        plt.clf()
        plt.imshow(self.H)
        plt.clim([-temp,temp])

    def M_plot(self, temp = 0.05):
        ''' MatrixPlot of the HesseMatirxi '''
        import matplotlib
        import matplotlib.pyplot as plt
        plt.clf()
        plt.imshow(self.M, cmap='RdBu',vmin=self.m.min(), vmax=self.m.max())
        plt.colorbar()   # makes the bar showing the vales
        print("a colormap would be great ....")
        #plt.clim([-temp,temp])

    def get_freqs(self, tol = 1e-5): #1e-8):
        ''' returnes Exact Freqs in [meV] of supercell without the three zero frequencies '''
        from numpy import linalg

        # You have: hbar(hartree/(bohrradius^2 u))^(1/2)
        # You want: meV
        #    * 637.33912
        dynMatToFreq = 637.339113585464

        # You have: hbar(eV/(angstrom^2 u))^(1/2)
        # You want: meV
        #    * 64.654148
        #    / 0.015466912
        dynMatToFreq = 64.654148

        d = self.H/np.sqrt(self.M)
        eigenvalues,   eigenvectors = linalg.eig(d)
        eigenvalues_realpart = np.sort(eigenvalues.real)
        eigenvalues_realpart_unsorted = eigenvalues.real
        self.eigenvectors = eigenvectors
        if self.__verbose:
            print("-------------- eigenvectors ----------------------")
            print(eigenvectors)
            print("-------------- eigenvalues -----------------------")
            print(eigenvalues)
            print("-------------- eigenvalues_realpart --------------")
            print(eigenvalues_realpart)
            print("--------------------------------------------------")

        #### Chop vales smaller 1e-8 (this is little)
        from numpy import ma
        eigenvalues_out = ma.masked_inside(eigenvalues_realpart,-tol, tol)
        eigenvalues_out_unsorted = ma.masked_inside(eigenvalues_realpart_unsorted,-tol, tol)
        ev =  eigenvalues_out.filled(0)
        evunsorted =  eigenvalues_out_unsorted.filled(0)
        if self.__verbose:
            print("-------------- ev -----------------------")
            print(ev)
            print("-----------------------------------------")

        ### check if we have negative parts
        if ev[0] < 0:
            if self.show_negative_eigenvalues:
                _printred("NEGATIVE EIGENVALUES!")
                print("ev[0]:",ev[0])
                print("ev:",ev)
            self.freqsNEGATIVE = ev*dynMatToFreq
            if self.show_negative_eigenvalues:
                np.savetxt("ExactFreqs_NEGATIVE",self.freqsNEGATIVE)
                print("np.sqrt(ev)*dynMatToFreq:",np.sqrt(ev)*dynMatToFreq)
            sys.exit("ERROR: Negative Eigenvalues : "+str(ev[0]))

        a = np.sqrt(ev)*dynMatToFreq
        aunsorted = np.sqrt(evunsorted)*dynMatToFreq


        ### check if first 3 eigenvalues are 0
        if a[0] == -0. : a[0] = 0
        if a[1] == -0. : a[1] = 0
        if a[2] == -0. : a[2] = 0
        if a[0] != 0 : my.exit("a[0] "+str(a[0])+" is not 0")
        if a[1] != 0 : my.exit("a[1] "+str(a[1])+" is not 0")
        if a[2] != 0 : my.exit("a[2] "+str(a[2])+" is not 0")

        out = np.nan_to_num(a[3:])
        outunsorted = np.nan_to_num(aunsorted)
        #self.freqs = out[::-1]
        self.freqs = out
        self.freqsunsorted = outunsorted

        if self.__verbose:
            print("-------------- self.freqs -----------------------")
            print(self.freqs)
            print("-------------------------------------------------")
        self.freqsmean = np.mean(self.freqs)


        return self.freqs

    def get_ene_atom(self):
        ''' returns free energy in [meV] as a function of T in [K] per atom '''
        Tmax = int(self._temp)  # 7K
        Ts = list(range(Tmax+1))      # 0..7K
        out = np.zeros((len(Ts),2))  # 8
        out[:] = np.nan
        out[:,0] = Ts

        for T in Ts:
            TT = T
            if T == 0: TT = 1  ## just to avoid durch 0 teilen
            out[T,1] = fqh_atom(TT, self.freqs )

        self.ene_atom = out
        return self.ene_atom

    def get_ene_cell(self):
        ''' returns free energy in [meV] as a function of T in [K] per supercell '''
        out = np.copy(self.get_ene_atom())
        out[:,1] = out[:,1]*float(self.atoms)
        return out

    def write_freqs(self, filename = 'ExactFreqs'):
        if filename == None:
            filename = 'ExactFreqs'
        if type(self.freqs) != bool:
            if self._verbose:
                print("filename                 :",filename)
                print("self._filename_addstring :",self._filename_addstring)
                print("self.freqs:",self.freqs)
            np.savetxt(filename + self._filename_addstring, self.freqs, fmt='%.12f')
            np.savetxt(filename + self._filename_addstring+"_unsorted", self.freqsunsorted, fmt='%.12f')

    def write_ene_atom(self, filename = "Fqh_fromExactFreqs_peratom" ):
        if filename == None:
            filename = "Fqh_fromExactFreqs"
        if self._verbose:
            print("--------------- ka ------------------------")
            print("self.filename",filename)
            print("self._filename_addstring:",self._filename_addstring)
            print("self.listin:","_".join(self.listin))
            print("---------------- ka done ------------------")

        if type(self.ene_atom) != bool:
            np.savetxt(filename + self._filename_addstring, self.ene_atom,fmt='%.0f %.12f')

    def write_ene_cell(self, filename = "Fqh_fromExactFreqs_percell" ):
        if filename == None:
            filename = "Fqh_fromExactFreqs_percell"
        if filename == "Fqh":
            filename = filename + self._filename_addstring
        if self._verbose:
            print("------------------ kb ---------------------")
            print("self.filename",filename)
            print("self._filename_addstring:",self._filename_addstring)
            print("---------------------------------------")
        if self.ene_cell != None:
            np.savetxt(filename, self.ene_cell,fmt='%.12f')


    def write_hessematrix(self, filename = "HesseMatrix_sphinx" ):
        if filename == None:
            filename = "HesseMatrix_sphinx"
        if filename == "HesseMatrix_sphinx":
            filename = filename + self._filename_addstring
        if self._verbose:
            print("----------------- kc ----------------------")
            print("self.filename",filename)
            print("self._filename_addstring:",self._filename_addstring)
            print('self.ene_cell',self.ene_cell)
            print("---------------------------------------")
        if type(self.ene_cell) == None or type(self.ene_cell) != bool:
            np.savetxt(filename, self.H/97.173617)
    #####################################################
    # cooperate with phonopy
    #####################################################

    def phonopy_FORCE_CONSTANTS_to_h(self, filename = "FORCE_CONSTANTS"):
        if os.path.isfile(filename) != True:
            sys.exit("file "+filename+" does not exist; Exit")
        print("Importing phonopy "+filename)
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fc = np.genfromtxt(filename, skip_header=2,invalid_raise=False)
        fc_lines = fc.shape[0]
        if fc.shape[1] != 3:
            sys.exit(filename+" has wrong dimensions")
        atoms = int(np.sqrt(fc_lines/3))
        h_lines = atoms*3

        if atoms*atoms*3 != fc_lines:
            sys.exit("number of fc_lines wrong in "+filename)

        h = np.zeros((atoms*3,atoms*3))
        h[:] = np.nan

        for i in np.arange(atoms): # i = 0..31 = line of hessematrix
            fc_3hesselines = np.array(np.vsplit(fc,atoms))[i]
            for j in np.arange(atoms): # j = 0..31 = atom
                #print "i_line:",i_line,"j_atom:",j_atom
                fc_3hesselines = fc[i*h_lines:i*h_lines+h_lines]
                h[i*3:i*3+3][:,j*3:j*3+3] = fc_3hesselines[j*3:j*3+3]
        #print "atoms:",atoms,"fc_lines:",fc_lines,"h_lines:",h_lines
        #return fc, ph, fc_3hesselines
        return h



def inversepot(r,alpha,B,e):
    '''
    returnes energy in [eV]
    alpha (no units)
    B (Angstrom)
    e (meV)
    in alfes paper:
        alpha = 6.7
        B = 1.85
        e = 1.
    '''
    return 4.*e*((B/r)**alpha)

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

def Morse(r,De,aa,re):
    ''' Morse potential
    returnes energy in [eV]
    r: distance of atoms [angstrom]
    re: nearest neighbor distance [angstrom]
    (a - b)^2 = (a - b)(a - b) = a^2 - 2ab + b^2
    De*(1.-np.exp(-aa*(r-re)))**2  ==
    De^2 - 2*De*np.exp(-aa*(r-re)) + np.exp(-aa*(r-re))^2

    subst: r=(r-re)
    subst: a=-aa
    derivative of np.exp(a*r)^2 == 2a np.exp(2*a*r)
    derivative of np.exp(a*r)   ==  a np.exp(  a*r)

    '''
    return De*(1.-np.exp(-aa*(r-re)))**2

def Morse_derivative(r,De,aa,re):
    ''' derivative of Morse potential
    returnes [eV/angstrom]
    import hesse
    hesse.Morse_derivative(2.82035196530504,0.279300,1.432673,4.13/np.sqrt(2))

    2.*aa*De*np.exp(-aa*(r-re))*(1.-np.exp(-aa*(r-re)))

    subst: np.exp(-aa*(r-re)) == dum
    2.*aa*De*dum*(1.-dum)

    '''
    return 2.*aa*De*np.exp(-aa*(r-re))*(1.-np.exp(-aa*(r-re)))

def mc1(r,De,aa,re,B,A):
    ''' Energy: x(mathematica) -> r-re(python) '''
    #return Morse(r,De,aa,re)+A*np.exp((-3.+B)*(r-re))+\
    #        (-(6.*(-7.+B)/(-3.+B)**5.)+(6.*(-7.+B)*(r-re)/(-3.+B)**4.)+\
    #        (3.*(-7.+B)*(r-re)**2./(-3.+B)**3.)+((-7.+B)*(r-re)**3./(-3.+B)**2.)+(r-re)**4./(-3.+B))

    #+\(-(6.*(-7.+B)/(-3.+B)**5.)+(6.*(-7.+B)*(r-re)/(-3.+B)**4.)+\
    #(3.*(-7.+B)*(r-re)**2./(-3.+B)**3.)+((-7.+B)*(r-re)**3./(-3.+B)**2.)+(r-re)**4./(-3.+B))
    return Morse(r,De,aa,re)+A*((r-re)**3.)*np.exp(-(3.-B)*(r-re))

def mc1_derivative(r,De,aa,re,B,A):
    ''' force '''
    #return Morse_derivative(r,De,aa,re) + A*(r-re)**3.*(1.+(r-re))*np.exp(-(3.-B)*(r-re))
    return Morse_derivative(r,De,aa,re) + \
            A*np.exp((-3.+B)*(r-re))* (r-re)**2 * (3.+(r-re)*(1.+B+(-3.+B)*(r-re)))

def Trans135(r,a,b,c):
    return (a*r**2.)/2.+(b*r**4.)/4.+(c*r**6.)/6.

def Trans135_der(r,a,b,c):
    return  a*r + b*r**3. + c*r**5.

def Trans1357(r,a,b,c,d):
    return  (a*r**2.)/2.+ \
            (b*r**4.)/4.+ \
            (c*r**6.)/6.+ \
            (d*r**8.)/8.

def Trans13579(r,a,b,c,d,e):
    return  (a*r**2.)/2.+ \
            (b*r**4.)/4.+ \
            (c*r**6.)/6.+ \
            (d*r**8.)/8.+ \
            (e*r**10.)/10.

def Trans13579_der(r,a,b,c,d,e):
    return  a*r + \
            b*r**3. + \
            c*r**5. + \
            d*r**7. + \
            e*r**9.

def Trans1357_der(r,a,b,c,d):
    return  a*r + \
            b*r**3. + \
            c*r**5. + \
            d*r**7.

def Trans1357911(r,a,b,c,d,e,f):
    return  (a*r**2.)/2.+ \
            (b*r**4.)/4.+ \
            (c*r**6.)/6.+ \
            (d*r**8.)/8.+ \
            (e*r**10.)/10. + \
            (f*r**12.)/12.

def Trans1357911_der(r,a,b,c,d,e,f):
    return  a*r + \
            b*r**3. + \
            c*r**5. + \
            d*r**7. + \
            e*r**9. + \
            f*r**11.




def hesse_to_hesse_1NN(
        hfilename = False,
        h = False,

        coordfile0_cart = False,
        coord0_cart = False,
        coordfile0_rrel = False,
        coord0_rrel = False,

        cellfile = False,
        cell = False):
    verbose = False
    if verbose:
        print('---- hesse_to_hesse_1NN ----')
        print("hfilename:",hfilename)
        print("cellfile:",cellfile)
        print("cell:",cell)


    if type(hfilename) == bool and type(h) == bool:
        sys.exit("please provied h or hfilename")
    if type(hfilename) != bool and type(h) != bool:
        sys.exit("please provied h or hfilename, not both")
    if type(hfilename) != bool:
        h = np.loadtxt(hfilename)*97.173617
        hout1NN = np.zeros((h.shape[0],h.shape[0]))
        hout1NNno = np.copy(h)

    crystal0 = crystal_generator.crystal()
    crystal0.load_positions_cell(coordfile_cart = coordfile0_cart, coord_cart = coord0_cart,
            cellfile = cellfile, cell = cell,
            coordfile_rrel = coordfile0_rrel, coord_rrel = coord0_rrel)
    positions = np.copy(crystal0.rcar)
    if verbose:
        print("---crystal0.cell---")
        print(crystal0.cellvec)
        print('---- positions ----')
        print(positions)
    if hout1NN.shape[0]/3 != positions.shape[0]:
        sys.exit("h and pos have diffrent size")


    ######################################################
    for atomidx,pos in enumerate(positions):
        NNlist1 = crystal0.get_NNlist(atomidx, 1,
                coordfile_cart = coordfile0_cart, coord_cart = coord0_cart,
                cell = crystal0.cellvec,
                coordfile_rrel = coordfile0_rrel, coord_rrel = coord0_rrel, return_NNdist = False)
        print("atomidx:",atomidx,NNlist1)
        ###############################
        # go through every 1st NN atom
        ###############################
        for atomjdx in NNlist1:
            hout1NN[atomidx*3+0][atomjdx*3+0] = h[atomidx*3+0][atomjdx*3+0]
            hout1NN[atomidx*3+0][atomjdx*3+1] = h[atomidx*3+0][atomjdx*3+1]
            hout1NN[atomidx*3+0][atomjdx*3+2] = h[atomidx*3+0][atomjdx*3+2]

            hout1NN[atomidx*3+1][atomjdx*3+0] = h[atomidx*3+1][atomjdx*3+0]
            hout1NN[atomidx*3+1][atomjdx*3+1] = h[atomidx*3+1][atomjdx*3+1]
            hout1NN[atomidx*3+1][atomjdx*3+2] = h[atomidx*3+1][atomjdx*3+2]

            hout1NN[atomidx*3+2][atomjdx*3+0] = h[atomidx*3+2][atomjdx*3+0]
            hout1NN[atomidx*3+2][atomjdx*3+1] = h[atomidx*3+2][atomjdx*3+1]
            hout1NN[atomidx*3+2][atomjdx*3+2] = h[atomidx*3+2][atomjdx*3+2]


            hout1NNno[atomidx*3+0][atomjdx*3+0] = 0.0
            hout1NNno[atomidx*3+0][atomjdx*3+1] = 0.0
            hout1NNno[atomidx*3+0][atomjdx*3+2] = 0.0

            hout1NNno[atomidx*3+1][atomjdx*3+0] = 0.0
            hout1NNno[atomidx*3+1][atomjdx*3+1] = 0.0
            hout1NNno[atomidx*3+1][atomjdx*3+2] = 0.0

            hout1NNno[atomidx*3+2][atomjdx*3+0] = 0.0
            hout1NNno[atomidx*3+2][atomjdx*3+1] = 0.0
            hout1NNno[atomidx*3+2][atomjdx*3+2] = 0.0

    for atomxyz in np.arange(h.shape[0])[0::3]:   # 1 zeile 4 zeile ... der HesseMatrix
        xf = hout1NN[atomxyz][0::3].sum()
        yf = hout1NN[atomxyz][1::3].sum()
        zf = hout1NN[atomxyz][2::3].sum()
        hout1NN[atomxyz][atomxyz+0] = -xf
        hout1NN[atomxyz][atomxyz+1] = -yf
        hout1NN[atomxyz][atomxyz+2] = -zf

        hout1NNno[atomxyz][atomxyz+0] = 0.0
        hout1NNno[atomxyz][atomxyz+1] = 0.0
        hout1NNno[atomxyz][atomxyz+2] = 0.0
        xfno = hout1NNno[atomxyz][0::3].sum()
        yfno = hout1NNno[atomxyz][1::3].sum()
        zfno = hout1NNno[atomxyz][2::3].sum()
        hout1NNno[atomxyz][atomxyz+0] = -xfno
        hout1NNno[atomxyz][atomxyz+1] = -yfno
        hout1NNno[atomxyz][atomxyz+2] = -zfno

    for atomxyz in np.arange(h.shape[0])[1::3]:  # 2 zeile 5 zeile ...
        xf = hout1NN[atomxyz][0::3].sum()
        yf = hout1NN[atomxyz][1::3].sum()
        zf = hout1NN[atomxyz][2::3].sum()
        hout1NN[atomxyz][atomxyz+0-1] = -xf
        hout1NN[atomxyz][atomxyz+1-1] = -yf
        hout1NN[atomxyz][atomxyz+2-1] = -zf

        hout1NNno[atomxyz][atomxyz+0-1] = 0.0
        hout1NNno[atomxyz][atomxyz+1-1] = 0.0
        hout1NNno[atomxyz][atomxyz+2-1] = 0.0
        xfno = hout1NNno[atomxyz][0::3].sum()
        yfno = hout1NNno[atomxyz][1::3].sum()
        zfno = hout1NNno[atomxyz][2::3].sum()
        hout1NNno[atomxyz][atomxyz+0-1] = -xfno
        hout1NNno[atomxyz][atomxyz+1-1] = -yfno
        hout1NNno[atomxyz][atomxyz+2-1] = -zfno

    for atomxyz in np.arange(h.shape[0])[2::3]:  # 3 zeile 6 zeile ...
        xf = hout1NN[atomxyz][0::3].sum()
        yf = hout1NN[atomxyz][1::3].sum()
        zf = hout1NN[atomxyz][2::3].sum()
        hout1NN[atomxyz][atomxyz+0-2] = -xf
        hout1NN[atomxyz][atomxyz+1-2] = -yf
        hout1NN[atomxyz][atomxyz+2-2] = -zf

        hout1NNno[atomxyz][atomxyz+0-2] = 0.0
        hout1NNno[atomxyz][atomxyz+1-2] = 0.0
        hout1NNno[atomxyz][atomxyz+2-2] = 0.0
        xfno = hout1NNno[atomxyz][0::3].sum()
        yfno = hout1NNno[atomxyz][1::3].sum()
        zfno = hout1NNno[atomxyz][2::3].sum()
        hout1NNno[atomxyz][atomxyz+0-2] = -xfno
        hout1NNno[atomxyz][atomxyz+1-2] = -yfno
        hout1NNno[atomxyz][atomxyz+2-2] = -zfno

    hout1NN = hout1NN/97.173617
    hout1NNno = hout1NNno/97.173617
    np.savetxt(hfilename+"_1NN",hout1NN,fmt="%.14f")
    np.savetxt(hfilename+"_-1NN",hout1NNno,fmt="%.14f")
    return hout1NN,hout1NNno

def get_energy_forces(
        pot = False,   # choices h (hesse), l (lennard jones), m (Morse)
        potparam = False,  # only used in case of pot == "l" or pot == "m"
        pot2 = False,
        hessefile = False,
        h = False,
        hessefile1nn = False,
        h1nn = False,

        coordfile_cart = False,
        coord_cart = False,
        coordfile_rrel = False,
        coord_rrel = False,

        coordfile0_cart = False,
        coord0_cart = False,
        coordfile0_rrel = False,
        coord0_rrel = False,

        #ddposfile = False,
        #pos = False,
        #posfile0 = False,
        #pos0 = False,

        cellfile = False,
        cell = False,

        add_transverse = False,
        hrest1 = False,
        returndu = False,
        printresult=False,
        only_return_forces_harmonic = False,

        verbose=False

        ):
    '''
    calculate energies and forces
    pot choices: h (hesse), l (lennard jones), m (Morse)
    please provide a hessefile or h as a numpy matrix

    '''
    if verbose == True:
        print()
        print()
        print("STRAT get_energy_forces() !!")
        print("verbose      :",verbose)
        #print ">hessefile:",hessefile
    #############################################
    # get positions
    #############################################
    if type(coordfile0_cart) == bool and type(coord0_cart) == bool \
        and type(coordfile0_rrel) == bool and type(coord0_rrel) == bool:
            if os.path.isfile("EqCoords_direct") == True:
                coordfile0_rrel = "EqCoords_direct"
    if type(coordfile0_cart) == bool and type(coord0_cart) == bool \
        and type(coordfile0_rrel) == bool and type(coord0_rrel) == bool:
            sys.exit("you need the undisplaced structure (e.gl. EqCoords_direct) \
                    as reference to get 1NN atoms (2)")
    if type(cell) == bool and type(cellfile) == bool:
        if os.path.isfile("cell") == True:
            cellfile = "cell"
    if type(cell) == bool and type(cellfile) == bool:
        if os.path.isfile("POSCAR") == True:
            utils.run2("rm -f cell; POSCAR_cell_cartesian.sh > cell") # to create cell
        if os.path.isfile("cell") == True:
            cellfile = "cell"
    if type(cell) == bool and type(cellfile) == bool:
        sys.exit("you need the cell file")

    crystal = crystal_generator.crystal()
    crystal0 = crystal_generator.crystal()
    #print 'cell',cell,type(cell)
    #print 'cellfile',cellfile
    crystal.load_positions_cell(coordfile_cart = coordfile_cart, coord_cart = coord_cart,
            cellfile = cellfile, cell = cell,
            coordfile_rrel = coordfile_rrel, coord_rrel = coord_rrel)
    crystal0.load_positions_cell(coordfile_cart = coordfile0_cart, coord_cart = coord0_cart,
            cellfile = cellfile, cell = cell,
            coordfile_rrel = coordfile0_rrel, coord_rrel = coord0_rrel)

    #print "<< pot:     ",pot
    #print "<< potparam:",potparam
    #############################################
    # get pot
    #############################################
    fehlermeldung = "choose one of pot choices: \"h\" (hesse), \"l\" (lennard jones), \
            \"m\" (Morse) \"mc1\" (mc1), \"i\" (inverse potential) \
            or \"extern\" (calculate_energy_and_forces file)"
    if type(pot) == bool:
        sys.exit(fehlermeldung)
    if pot == "l" or pot == "h" or pot == "m" or pot == 'mc1' or pot == 'i' or pot == 'extern' or pot == 'externorig':
        pass
    else:
        sys.exit(fehlermeldung)
    if pot == "l" or pot == "m" or pot == 'mc1' or pot == 'i':
        if type(potparam) == bool:
            sys.exit("plase provide parameters for lennard jones or morse potential in \
                    potparam variable; e.g. potparam = [ 0.1, 2,86 ] (lennard jones) or \
                    potparam = [0.234631, 1.48197, 2.86] (morse)")

    ##############################################################################
    # potential
    ##############################################################################
    #################
    # harmonic
    #################
    if pot == "h" or hrest1 == True:
        if type(hessefile) == bool and type(h) == bool:
            if os.path.isfile("HesseMatrix_sphinx"):
                hessefile = "HesseMatrix_sphinx"
        if type(hessefile) == bool and type(h) == bool:
            sys.exit("please provide a hessefile or h as a numpy matrix")
        if type(hessefile) != bool and type(h) != bool:
            sys.exit("please provide a hessefile or h as a numpy matrix")
            sys.exit("please provide either a hessefile or h as a numpy matrix, not both")

        # h
        if hessefile != False:
            h = read_Hessematrix(hessematrixfile = hessefile)
        else:
            h = h

        if hrest1 == True:
            hessefile1NN = hessefile+"_1NN"
            if os.path.isfile(hessefile1NN) != True:
                sys.exit(hessefile1NN+" does not exist")

            # else file exists
            h1NN = read_Hessematrix(hessematrixfile = hessefile1NN)

        # map atoms back onto original positions
        if crystal.rcar.shape[0] == 32 or crystal.rcar.shape[0] == 108:
            pass
        else:
            sys.exit("The mapping back into the cell currently only works for fcc")
        #print "=============crystal.rcar"
        #print crystal.rcar
        #print "=============crystal0.rcar"
        #print crystal0.rcar
        #print "============="
        u = crystal.rcar-crystal0.rcar
        #print "uorig:",u
        for ind,i in enumerate(u):
            if  u[ind][0] > crystal.cellvec[0,0]/2:
                u[ind][0] = u[ind][0] - crystal.cellvec[0,0]
            if  u[ind][1] > crystal.cellvec[0,0]/2:
                u[ind][1] = u[ind][1] - crystal.cellvec[0,0]
            if  u[ind][2] > crystal.cellvec[0,0]/2:
                u[ind][2] = u[ind][2] - crystal.cellvec[0,0]
        #print "u out:",u
        if returndu:
            return u



        if hrest1 == True:
            fh1NN = qh_forces(u,h1NN)
            eh1NN = qh_energy_cell(u,h1NN)
            ehmev1NN = eh1NN*1000/(crystal.rcar.shape[0]-1)
        #print "vor fh:",
        print("u:",u)
        print('h:',h)
        fh = qh_forces(u,h)
        eh = qh_energy_cell(u,h)
        ehmev = eh*1000/(crystal.rcar.shape[0]-1)
        print("fh:",fh)
        print("eh:",eh,"ehmev:",ehmev)
        if only_return_forces_harmonic:
            return fh,eh
        np.savetxt("forces_harmonic",fh)




        if hrest1 == True:
            fhrest = fh - fh1NN
            ehrest = eh - eh1NN
            ehmevrest = ehmev - ehmev1NN

        if printresult and hrest1 == True:
            print("ehrest  [eV/cell]:",ehrest,"  ehrest  [meV/atom]:",ehmevrest)
            print(fhrest)
            print("")
        if printresult:
            print("eh  [eV/cell]:",eh,"  eh  [meV/atom]:",ehmev)
            print(fh)
        return ehmev,eh,fh

    #################
    # Morse, Lennard Jones
    #################
    if pot == "i":
        usepot=[ 'inversepot', potparam[0], potparam[1], potparam[2] ]
    if pot == "l":
        usepot=[ 'LJ', potparam[0], potparam[1] ]
    if pot == "m":
        usepot=[ 'Morse', potparam[0], potparam[1], potparam[2] ]
    if pot == "mc1":
        usepot=[ 'mc1', potparam[0], potparam[1], potparam[2], potparam[3], potparam[4] ]
    if pot == "extern":
        usepot = [ 'extern' ]
    if pot == "externorig":
        usepot = [ 'externorig' ]
    if pot == "l" or pot == "m" or pot == 'mc1' or pot == 'i' or pot == 'extern' or pot == 'externorig':
        if verbose:
            print("usepot       :",usepot,type(usepot))
            #print "cryst0:",crystal0.rcar
            #print "---"
            #print crystal0.cellvec
            print("hh:",hessefile,hessefile1nn)


        eah0, fah0 = get_energy_forces_pairpot(
                coord_cart = crystal0.rcar,
                cell = crystal0.cellvec,
                NNshell = 1,
                usepot=usepot,
                #pot2=pot2,
                add_transverse = add_transverse,
                hessefile = hessefile,
                hessefile1nn = hessefile1nn,
                return_NNdist = False,
                verbose = verbose)


        #print eah0, fah0
        #print "---------------------------"
        #print "------ NOW DISPLACED STRUCTURE ---------------------"
        eahausgelenkt, fah = get_energy_forces_pairpot(
                coord_cart = crystal.rcar,
                cell = crystal.cellvec,
                NNshell = 1,
                usepot=usepot,
                #pot2=pot2,
                add_transverse = add_transverse,
                hessefile = hessefile,
                hessefile1nn = hessefile1nn,
                verbose = verbose)
        #print "eahausgelenkt",eahausgelenkt,eahausgelenkt*1000/(crystal.rcar.shape[0]-1)
        #print "eah0:",eah0
        eah = eahausgelenkt - eah0
        eahmev = eah*1000/(crystal.rcar.shape[0]-1)
        if printresult:
            print("eah [eV/cell]:",eah,"  eah [meV/atom]:",eahmev)
            print(fah)
        return eahmev,eah,fah

def getefvec(vec,params,pot = False):
    vecnorm=np.linalg.norm(vec)
    if abs(vecnorm) <= 0.000000000000000000000001:
        return 0.0,np.array([0.0,0.0,0.0])
    #print "vn:",vecnorm,params,pot
    enorm,fnorm = getef(vecnorm,params,pot)
    fvec = vec/vecnorm*fnorm
    return enorm, fvec

def getef(vecnorm,params,pot = False):
    '''
    help text
    '''
    #print "params:",params,type(params)
    if type(params) == bool:
        return 0.0,0.0
    #print "    vecnorm:",vecnorm
    #print "    params:",params
    #print "    pot:",pot
    if pot == False or pot == 'poly': # np.poly is assumed
        #f = np.polynomial.Polynomial(params)(vecnorm)
        #fNOWWRONG = np.poly1d(params[::-1])(vecnorm)  # i assume poly1d to be quicker also never tested
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
    if pot == 'i':
        function = inversepot
        functionder = inversepot_derivative
    if pot == 'l':
        function = LJ
        functionder = LJ_derivative
    if pot == 'm':
        function = Morse
        functionder = Morse_derivative
    if pot == 'mc1':
        function = mc1
        functionder = mc1_derivative
    if function == None:
        sys.exit("pot not recognized")
    #e =    function(longvecnorm,*params)
    #f = functionder(longvecnorm,*params)
    #print "vecnorm:",vecnorm
    #print "params:",params
    e =    function(vecnorm,*params)
    f = functionder(vecnorm,*params)
    return e,f

def rotate_vec0_to_vec0_original_inplane_and_getforce_current(vec0_curr,longvec,vec0_orig_formapping,mappedvec_orig_dir,forcefuncx,forcefuncy,forcefuncz,
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
        print("tvecwirkdir:",tvecwirkdir)
        print("type(tvecwirkdir):",type(tvecwirkdir))
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
        print("  XXtvec_curr_dir:",tvec_curr_dir)
        print("  XXtvec_curr", tvec_curr)
        print("  XXlongvec_curr",longvec_curr)
        sys.exit()


    fromvec     = vec0_curr
    tovec       = vec0_orig_formapping

    # TODO (much later):
    # define out of plane and inplane direction, this we can not do internally in this module but have to define externally ant then give to this module only the tvec in the
    # correct inplane / out of plane part. Practically what "outofplane" means will depend on the CURRENT position of the respective atoms and not on the undisplaced vec0 positions.
    # We could consider this later on!




    Rotmat      = utils.R_2vect(vec0_curr,vec0_orig_formapping)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS
    Rotmatback  = utils.R_2vect(vec0_orig_formapping,vec0_curr)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS
                                                 # if mappedvec_orig_dir is perpendicular to the plane spaned by fromvec(==vec0_curr),tovec -> no problem
                                                 # otherwise: further checking necessary
                                                 # in other words: we have to make sure that tvec points in same direction as originally!
    Rotmatt      = utils.R_2vect(tvec_curr,mappedvec_orig_dir)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS
    Rotmattback  = utils.R_2vect(mappedvec_orig_dir,tvec_curr)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS


    # check if mappedvec_orig_dir is parallel to the perpendicular vector spanned by {from,to}vec
    #perpendicular_fromvec_tovec = np.dot(np.array(fromvec),np.array(tovec))


    # after applying rotation on ... there should be a vector like (x,0,0) or (0,y,0) or (0,z,0)
    tvec_orig1 = np.dot(Rotmat,tvec_curr)
    longvec_orig1 = np.dot(Rotmat,longvec_curr)

    # with Rotmat2 we need to make sure that it only rotates around vec0_orig_formapping; what happend previously is taat in order to ensure that tvec points in the correct direction
    # vor the -1/-1 direction Rotvec was simply movec back (negative einheitsmatrix) this can only happen if direction vec0 is in the same direction as vec0_orig_formapping and both
    # go in exactly opposite directions. To avoid this one needs to make sure that Rotmat2 only rotates around vec0_orig_formapping. Typically Rotmat2 will only be a rotation around
    # vec0_orig_formapping by 180 degrees or none.
    # Dies hier funktionier noch nicht 100% ig

    #Rotmat2      = utils.R_2vect(tvec_orig1,mappedvec_orig_dir,fixed_rotation_axis=vec0_orig_formapping)  # the length of the vectors does not play a role!
    #Rotmatback2  = utils.R_2vect(mappedvec_orig_dir,tvec_orig1,fixed_rotation_axis=vec0_orig_formapping)

    Rotmat2      = utils.R_2vect(tvec_curr,mappedvec_orig_dir)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS
    Rotmatback2  = utils.R_2vect(mappedvec_orig_dir,tvec_curr)   # the directon which is perpendicular to tovec and dromvec WILL NOT BE DESCRIBED BY THIS

    tvec_orig = np.dot(Rotmat2,tvec_orig1)
    longvec_orig = np.dot(Rotmat2,longvec_orig1)

    ####################################################################################################################
    # bei der rotation back muessen wir auch sicherstellen dass der tvec_curr wieder hergestellt wird
    # also brauchen wir auch da wieder 2 rotationen zurueck!!!!!
    ####################################################################################################################



    # now we want to reject the tvec_orig onto the originally mapped vector
    rejection_orig = utils.reject_vector(longvec_orig, vec0_orig_formapping)
    tvec_mapped = utils.project_vector(rejection_orig,mappedvec_orig_dir)   # hier ist es egal ob man tvecoutdir oder tvecoutdir*2 schreibet, nur die richtung ist von interesse


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
    f_o = np.dot( R.T, np.dot(R2, tvec_curr) )
    #######################################################################################################################################################
    # check if tvec_o is parallel to mappedvec_orig_dir
    perpendicular_fromvec_tovec = np.cross(tvec_o,mappedvec_orig_dir)



    # anwenden der forcefuncx, forcefuncy
    #ex, fx_ = getef(math.copysign(np.linalg.norm(tvec_mapped),tvec_mapped[0]),forcefuncx,pot=False)
    #fxorig = np.array([fx_,0.0,0.0])
    #ey, fy_ = getef(math.copysign(np.linalg.norm(tvec_mapped),tvec_mapped[1]),forcefuncy,pot=False)
    #fyorig = np.array([0.0,fy_,0.0])
    #ez, fz_ = getef(math.copysign(np.linalg.norm(tvec_mapped),tvec_mapped[2]),forcefuncz,pot=False)
    #fzorig = np.array([0.0,0.0,fz_])

    # tvec_o will always be positive: in case of tix
    print("tvec_o:",tvec_o)
    print("forcefuncx:",forcefuncx)
    ex, fx_ = getef(math.copysign(np.linalg.norm(tvec_o),tvec_o[0]),forcefuncx,pot=False)
    print("fx_:",fx_)
    fxorig = np.array([fx_,0.0,0.0])
    ey, fy_ = getef(math.copysign(np.linalg.norm(tvec_o),tvec_o[1]),forcefuncy,pot=False)
    fyorig = np.array([0.0,fy_,0.0])
    ez, fz_ = getef(math.copysign(np.linalg.norm(tvec_o),tvec_o[2]),forcefuncz,pot=False)
    fzorig = np.array([0.0,0.0,fz_])

    fx = np.dot( R2.T, np.dot(R, fxorig) )
    fy = np.dot( R2.T, np.dot(R, fyorig) )
    fz = np.dot( R2.T, np.dot(R, fzorig) )

    print("fx:",fx)
    def print_vecs_to_screen():
        print(utils.printred("atomi:"+str(atomi)+" atomj:"+str(atomj)))
        print(utils.printred("bbbtext:"+str(bbbtext)))
        print(utils.printyellow(str(bbbtext)+" atomi:"+str(atomi)+" atomj:"+str(atomj)+" tvec_curr:",tvec_curr,"tvec_o:",tvec_o))
        print(utils.printyellow("tvec_orig1;",tvec_orig1,"tvec_orig:",tvec_orig)) #,fx_,fy_,fz_)
        print("  vec0_curr:",vec0_curr)
        print("  vec0_orig_formapping:",vec0_orig_formapping)
        print("  longvec:",longvec,np.linalg.norm(longvec))
        print("  tvecwirkdir:",tvecwirkdir)
        print("")
        print("  tvecfull  = longvec - vec0_curr:",tvecfull)
        print("  tvec_curr_dir = utils.project_vector(tvecfull,tvecwirkdir):",tvec_curr_dir)
        print("  rejection = utils.reject_vector(longvec, vec0_curr):",rejection)
        print("  tvec_curr = utils.project_vector(rejection,tvecwirkdir):",tvec_curr)
        print("  longvec_curr:",longvec_curr)
        print("  tvec_o = np.dot( R.T, np.dot(R2, tvec_curr) ):",tvec_o)
        print(" ")
        print("  R2:",R2)
        print("  np.dot(R2, tvec_curr):",np.dot(R2, tvec_curr))
        print("  np.dot( R.T, np.dot(R2, tvec_curr) ):",np.dot( R.T, np.dot(R2, tvec_curr) ))
        print("  tvec_o:",tvec_o)
        print("  mappedvec_orig_dir:",mappedvec_orig_dir)
        print("  perpendicular_fromvec_tovec:",perpendicular_fromvec_tovec)
        print("  forcefuncx:",forcefuncx)
        print("  forcefuncy:",forcefuncy)
        print("  forcefuncz:",forcefuncz)
        #print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        #print utils.printred("---------------------")
        #print "  tvec_orig1 (determined through rotation):",tvec_orig1
        #print utils.printred("---------------------")
        #print "  tvec_orig (determined through rotation):",tvec_orig
        #print utils.printred("---------------------")
        #print "  longvec_orig (determind through rotation):",longvec_orig
        #print "  mappedvec_orig_dir:",mappedvec_orig_dir
        #print "  np.cross(tvec_orig,mappedvec_orig_dir):",np.cross(tvec_orig,mappedvec_orig_dir)
        #print "  np.linalg.norm(np.cross(tvec_orig,mappedvec_orig_dir)):", np.linalg.norm(np.cross(tvec_orig,mappedvec_orig_dir))
        #print "  vec0_orig_formapping:",vec0_orig_formapping
        #print "  Rotmat:",Rotmat
        #print "  Rotmat2:",Rotmat2
        #print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
        #print "tvecwirkdir:",tvecwirkdir
        #print "tvec_curr:",tvec_curr
        #print "Rotmat:",Rotmat
        print("fx_,fy_,fz_:",fx_,fy_,fz_)
        print("fxorig,fyorig,fzorig:",fxorig,fyorig,fzorig)
        print("fx    ,fy    ,fz    :",fx    ,fy    ,fz)
        print("")
        #if True:
        #   sys.exit("not 0:")
        if abs(np.linalg.norm(perpendicular_fromvec_tovec)) >= 1e-9:
           print(utils.printred("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
           print("perpendicular_fromvec_tovec = np.cross(tvec_o,mappedvec_orig_dir)",perpendicular_fromvec_tovec)
           print("np.linalg.norm(perpendicular_fromvec_tovec)",np.linalg.norm(perpendicular_fromvec_tovec))
           print("tvec_o:",tvec_o)
           print("mappedvec_orig_dir:",mappedvec_orig_dir)
           print("tvec_orig1:",tvec_orig1," == np.dot(Rotmat,tvec_curr)","  tvec_curr:",tvec_curr)
           print(utils.printred("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"))
           sys.exit("not 0:")

    #if np.linalg.norm(perpendicular_fromvec_tovec) != 0.0:
    if False: #np.linalg.norm(perpendicular_fromvec_tovec) != 0.0:
       ex, fx_ = getef(math.copysign(np.linalg.norm(tvec_orig),tvec_orig[0]),forcefuncx,pot=False)
       fxorig = np.array([fx_,0.0,0.0])
       ey, fy_ = getef(math.copysign(np.linalg.norm(tvec_orig),tvec_orig[1]),forcefuncy,pot=False)
       fyorig = np.array([0.0,fy_,0.0])
       ez, fz_ = getef(math.copysign(np.linalg.norm(tvec_orig),tvec_orig[2]),forcefuncz,pot=False)
       fzorig = np.array([0.0,0.0,fz_])

       tvec_o = tvec_orig

       # Rotate back ( the forces, do the same also for positions to check if we ended up in tvec_curr )
       fx2 = np.dot(Rotmatback2,fxorig)
       fy2 = np.dot(Rotmatback2,fyorig)
       fz2 = np.dot(Rotmatback2,fzorig)
       fx  = np.dot(Rotmatback,fx2)
       fy  = np.dot(Rotmatback,fy2)
       fz  = np.dot(Rotmatback,fz2)
       check1 = np.dot(Rotmatback2,tvec_orig)
       check2 = np.dot(Rotmatback,check1)

       #print utils.printyellow("!!!!!!!!","bbbtext:",bbbtext,"atomi:",atomi,"tvec_curr:",tvec_curr)

    # second check in case tvec_o had to be gained by rotation
    perpendicular_fromvec_tovec = np.cross(tvec_o,mappedvec_orig_dir)

    #if np.linalg.norm(perpendicular_fromvec_tovec) != 0.0:
    if False: #np.linalg.norm(perpendicular_fromvec_tovec) != 0.0:
        print("perpendicular_fromvec_tovec = np.cross(tvec_o,mappedvec_orig_dir)")
        print(perpendicular_fromvec_tovec)
        print("tvec_o:",tvec_o)
        print("mappedvec_orig_dir:",mappedvec_orig_dir)
        print("tvec_orig1:",tvec_orig1," == np.dot(Rotmat,tvec_curr)","  tvec_curr:",tvec_curr)
        print('------')
        print("vec0_curr:",vec0_curr)
        print("vec0_orig_formapping:",vec0_orig_formapping)
        print("Rotmat:",Rotmat)
        print('------')
        print("tvec_curr:",tvec_curr)
        print("mappedvec_orig_dir:",mappedvec_orig_dir)
        print("Rotmat2:",Rotmat2)
        print("np.dot(Rotmat,Rotmat2):",np.dot(Rotmat,Rotmat2))
        sys.exit("not 0:")



    ##return fxorig,fyorig,fzorig,fx,fy,fz,longvec_orig,vec0_orig_formapping,tvec_mapped,tvec_orig,fx_,fy_,fz_
    #if False:
    #    if type(atomlist) != bool:
    #        if type(atomlist) == list:
    #            for i in atomlist:
    #                print "i:",i,"longvec_orig:",longvec_orig
    #                print "i:",i,"tvec_orig:",tvec_orig
    #                print "i:",i,"tvec_curr:",tvec_curr
    #                print "i:",i,"tvec_back:",check2
    #                print "i:",i,"fxorig,fyorig,fzorig  :",fxorig,fyorig,fzorig
    #                print "i:",i,"fx,fy,fz              :",fx,fy,fz
    #                print "i:",i,"Rotmat:",Rotmat
    #                print "i:",i,"Rotmat2:",Rotmat2

    # with the first roation we can never make sure that tvec points in the right direction, also it should show in the proper direction (make this sure by checking
    # versus the  mappedvec_orig_dir --> mappedvec_orig_dir has to be parallel! to tvec_orig


    #if abs(np.linalg.norm(np.cross(tvec_o,mappedvec_orig_dir))) >= 1e-8 or type(atomlist) != bool and atom in atomlist:

    if type(atomlist) != bool and (atomi in atomlist or atomj in atomlist):
        print_vecs_to_screen()

    return ex,ey,ez,fx,fy,fz

def get_energy_forces_pairpot(
        coordfile_cart = False,
        coord_cart = False,
        coordfile_rrel = False,
        coord_rrel = False,

        coordfile0_cart = False,
        coord0_cart = False,
        coordfile0_rrel = False,
        coord0_rrel = False,

        cellfile = False,
        cell = False,

        NNshell = 1,
        usepot = False,
        #pot2 = False,
        add_transverse = True, #False,

        h = False,      # to add 2NN to pairpot
        h1nn = False,   # to add 2NN to pairpot
        hessefile = False,      # to add 2NN to pairpot
        hessefile1nn = False,   # to add 2NN to pairpot
        verbose = False,
        return_NNdist = False,
        save_vecs_to_file_for_DOS = False):
    '''
    run e.g.: hesse.get_energy_forces_pairpot(coordfile_cart="cartesian_coords",
    usepot=[ 'LJ', 0.1, 2.92035])  no. seems not to work :(
    choose as potential: "LJ", "Morse" or "inversepot"
    '''
    ## BEGINPROG
    if verbose:
        print()
        print("START get_energy_forces_pairpot() !!!")
        print("verbose      :",verbose)
        print("usepot       :",usepot)
    print("hessefile    :",hessefile)
    print("hessefile1nn :",hessefile1nn)
    #############################################
    # get positions
    #############################################
    if type(hessefile) != bool and type(h) != bool:
        sys.exit("either h or hessefile, not both")
    if type(hessefile1nn) != bool and type(h1nn) != bool:
        sys.exit("either h1nn or hessefile1nn, not both")
    if type(hessefile) != bool:
        if os.path.isfile(hessefile) != True:
            sys.exit(hessefile+" not existing")
        if os.path.isfile(hessefile1nn) != True:
            sys.exit(hessefile1nn+" not existing")
    if type(coordfile_cart) == bool and type(coord_cart) == bool \
        and type(coordfile_rrel) == bool and type(coord_rrel) == bool:
            if os.path.isfile("cartesian_coords") == True:
                coordfile_cart = "cartesian_coords"
    if type(coordfile_cart) == bool and type(coord_cart) == bool \
        and type(coordfile_rrel) == bool and type(coord_rrel) == bool:
            sys.exit("you need the displaced structure \
                    e.gl. cartesian_coords")

    if type(coordfile0_cart) == bool and type(coord0_cart) == bool \
        and type(coordfile0_rrel) == bool and type(coord0_rrel) == bool:
            if os.path.isfile("EqCoords_direct") == True:
                coordfile0_rrel = "EqCoords_direct"
    if type(coordfile0_cart) == bool and type(coord0_cart) == bool \
        and type(coordfile0_rrel) == bool and type(coord0_rrel) == bool:
            sys.exit("you need the undisplaced structure \
                    (e.gl. EqCoords_direct) as reference to get 1NN atoms (3)")

    if type(cell) == bool and type(cellfile) == bool:
        if os.path.isfile("cell") == True:
            cellfile = "cell"
    if type(cell) == bool and type(cellfile) == bool:
        if os.path.isfile("POSCAR") == True:
            utils.run2("rm -f cell; POSCAR_cell_cartesian.sh > cell") # to create cell
        if os.path.isfile("cell") == True:
            cellfile = "cell"
    if type(cell) == bool and type(cellfile) == bool:
        sys.exit("you need the cell file")


    # We have to ensure we dont read in the EqCoords file every time;
    #print "0: cartesian coords: ---------------------"
    #print "1:",coordfile_cart
    #print "2:",coord_cart
    #print "3:",cellfile
    #print "4:",cell
    #print "5:",coordfile_rrel
    #print "6:",coord_rrel
    crystal = crystal_generator.crystal()
    crystal0 = crystal_generator.crystal()
    crystal.load_positions_cell(coordfile_cart = coordfile_cart, coord_cart = coord_cart,
            cellfile = cellfile, cell = cell,
            coordfile_rrel = coordfile_rrel, coord_rrel = coord_rrel)
    #print "0: EqCoordsDirect: ---------------------"
    #print "1:",coordfile0_cart
    #print "2:",coord0_cart
    #print "3:",cellfile
    #print "4:",cell
    #print "5:",coordfile0_rrel
    #print "6:",coord0_rrel
    crystal0.load_positions_cell(coordfile_cart = coordfile0_cart, coord_cart = coord0_cart,
            cellfile = cellfile, cell = cell,
            coordfile_rrel = coordfile0_rrel, coord_rrel = coord0_rrel)

    #print "0: EqCoordsDirect: ---------------------"
    # copy class instance to not have to load in againg every time neighbors were found
    # (when neighbors are found relative positions (and car coords) are changed)
    import copy
    crystal0neverchange = copy.deepcopy(crystal0)
    crystalneverchange = copy.deepcopy(crystal)
    numberofatoms = crystal.rcar.shape[0]



    # zeige alle moeglichen NN in dieser zelle && welche schalen sind voll drin in sc?
    # 6 in der 2x2x2fcc superzelle
    # 10 in der 3x3x3fcc superzelle
    #
    # in 3x3x3sz:
    # NNabst: [ 0.      2.9204  4.13    5.0582  5.8407  6.5301  7.1534  7.7265  8.7611
    #   9.6857]
    #   len: 10
    # 1     12  [ 0.     2.065  2.065]  full (gesampled in quer)
    # 2     6   [ 0.    0.    4.13]     full (gesampled in xdir)
    # 3     24  [ 4.13   2.065  2.065]  full (gesampled in 3NNdir)
    # 4     12  [ 0.    4.13  4.13]     full (gesampled in quer)
    # 5     12  [ 0.     2.065  6.195]
    # 6     8   [ 4.13  4.13  4.13]
    # 7     24  [ 4.13   2.065  6.195]
    # 8     3   [ 0.     6.195  6.195]
    # 9     6   [ 4.13   6.195  6.195]
    #
    # in 2x2x2sz:
    # NNabst: [ 0.      2.9204  4.13    5.0582  5.8407  7.1534]
    # len: 6
    # 1     12  (12) [ 0.     2.065  2.065]  full
    # 2     3   ( 6) [ 0.    0.    4.13]     NOTFULL  (gesampled in xdir)
    # 3     12  (24) [ 4.13   2.065  2.065]  NOTFULL  (gesampled in 3NNdir)
    # 4     3   (12) [ 0.    4.13  4.13]     NOTFULL  (gesampled in quer)
    # 5     1   [ 4.13  4.13  4.13]     NOTFULL  (gesampled in 111dir)

    # this will only hold now for fcc,bcc where we think of the 0th atom in the 'origin'
    #
    # This ist just for information to show!
    if verbose:
        NNabst = crystal0.get_NNlist(0, 1,
               cell = crystal0.cellvec,
               coord_rrel = crystal0neverchange.rrel,
               return_NNdist = True)
        #NNlist = crystal0.get_NNlist(0, 1,
        #       cell = cell,
        #       coord_rrel = crystal0neverchange.rrel,
        #       return_NNdist = False)
        #
        print("NNabst:",NNabst)
        print("len:",len(NNabst))
        for nnabstidx,nnatom in enumerate(NNabst):
            if nnabstidx == 0.0:
                continue
            nnlist = crystal0.get_NNlist(0, nnabstidx,
                    cell = crystal0.cellvec,
                    coord_rrel = crystal0neverchange.rrel,
                    return_NNdist = return_NNdist)
            print(nnabstidx,"\t",len(nnlist),"\t",crystal0neverchange.rcar[nnlist[0]])


    ################################################################
    ## Beginning of 2NN stuff
    ################################################################
    ## # create 32 atom supercell
    ## unit = crystal_generator.unit_cell()
    ## unit.load_fcc_cell( 1.0 )
    ## sc1 =  crystal_generator.supercell()
    ## sc1.create_supercell( unit, 2, 2, 2 )

    ## # copy things
    ## sc1.cellvec = crystal.cellvec
    ## sc1.rcar = crystal.rcar



    ## sc = crystal_generator.supercell()
    ## sc.create_supercell( sc1, 2, 2, 2 )
    ## #sc.create_supercell( crystal, 2, 2, 2 )

    ## print crystal.cellvec
    ## print crystal.rcar
    ## print sc.rrel.shape
    ## print sc.rcar


    #############################################
    # get pot and check pot and pot parameters
    #############################################

    fehlermeldung="you have to define the potential e.g. [ 'LJ', 2.8567, 1.2 ] or \
            [ 'inversepot', 6.7, 1.85, 1.0 ]"
    if usepot == False:
        print("usepot == False")
        sys.exit(fehlermeldung)
    if type(usepot) != list:
        print("usepot:",usepot,"is not a list")
        sys.exit(fehlermeldung)

    # add if new pot
    if usepot[0] != "LJ" and usepot[0] != "inversepot" and usepot[0] != "Morse" \
            and usepot[0] != "mc1" and usepot[0] != "extern" and usepot[0] != 'externorig':
        sys.exit(fehlermeldung)
    if usepot[0] == "LJ":
        if len(usepot) != 3:
            sys.exit(fehlermeldung)
    if usepot[0] == "inversepot" or usepot[0] == "Morse":
        if len(usepot) != 4:
            sys.exit(fehlermeldung)
    if usepot[0] == "mc1":
        if len(usepot) != 6:
            sys.exit(fehlermeldung)

    print("usepot[0]    :",usepot[0])
    ##################################################
    # get pot from calculate_energy_and_forces
    ##################################################
    if usepot[0] == "extern" or usepot[0] == 'externorig':
        if usepot[0] == 'externorig':  # in this case we want the calculate_energy_and_forces from the Thermodynamics folder
            import imp
            calculate_energy_and_forces = imp.load_source('calculate_energy_and_forces', os.environ['HOME']+'/Thermodynamics/python_thermodynamics/calculate_energy_and_forces')
        if usepot[0] == 'extern':  # in this case we want the calculate_energy_and_forces from the Thermodynamics folder
            import imp
            calculate_energy_and_forces = imp.load_source('calculate_energy_and_forces', os.getcwd()+'/calculate_energy_and_forces')
        from calculate_energy_and_forces import \
                u1nn_pot, \
                u1nn_potparam, \
                u1nn_potadd, \
                u1nn_potaddparam, \
                u1nn_pot2, \
                u1nn_pot2param, \
                u1nn_pot2add, \
                u1nn_pot2addparam, \
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



        u1nn_pot  # stays as it is defined (string)
        u1nn_potparam = utils.string_to_list(u1nn_potparam)
        u1nn_potadd  # stays as it is defined (string)
        u1nn_potaddparam = utils.string_to_list(u1nn_potaddparam)
        u1nn_pot2 = u1nn_pot2  #tays as it is defined (string)
        u1nn_pot2param = utils.string_to_list(u1nn_pot2param)
        u1nn_pot2add = u1nn_pot2  #tays as it is defined (string)
        u1nn_pot2addparam = utils.string_to_list(u1nn_pot2addparam)

        u1nn_topx = utils.string_to_list(u1nn_topx)
        u1nn_topy = utils.string_to_list(u1nn_topy)
        u1nn_topz = utils.string_to_list(u1nn_topz)
        u1nn_tipx = utils.string_to_list(u1nn_tipx)
        u1nn_tipy = utils.string_to_list(u1nn_tipy)
        u1nn_tipz = utils.string_to_list(u1nn_tipz)
        #pot2 = u1nn_pot2
        #potparam = u1nn_potparam
        u2nn_pot  # stays as it is defined (string)
        u2nn_potparam = utils.string_to_list(u2nn_potparam)
        u2nn_pot2 = u2nn_pot2 #string m or mc1
        u2nn_potparam2 = utils.string_to_list(u2nn_potparam2) # string m or mc1

        u2nn_topx = utils.string_to_list(u2nn_topx)
        u2nn_topy = utils.string_to_list(u2nn_topy)
        u2nn_topz = utils.string_to_list(u2nn_topz)
        u2nn_tipx = utils.string_to_list(u2nn_tipx)
        u2nn_tipy = utils.string_to_list(u2nn_tipy)
        u2nn_tipz = utils.string_to_list(u2nn_tipz)

        u3nn_pot  # stays as it is defined (string)
        u3nn_potparam = utils.string_to_list(u3nn_potparam)
        u3nn_pot2 = u2nn_pot2 #string m or mc1
        u3nn_potparam2 = utils.string_to_list(u3nn_potparam2) # string m or mc1


        def potprint(u1nn_pot2param,types = False):
            '''
            print pot to screen
            '''
            contrib = [ \
             'u1nn_pot' \
            ,"u1nn_potparam"       \
            ,"u1nn_potadd"       \
            ,"u1nn_potaddparam"       \
            ,"u1nn_pot2"           \
            ,"u1nn_topx"           \
            ,"u1nn_topy"           \
            ,"u1nn_topz"           \
            ,"u1nn_tipx"           \
            ,"u1nn_tipy"           \
            ,"u1nn_tipz"           \
            ,"u2nn_pot"            \
            ,"u2nn_potparam"       \
            ,"u2nn_pot2"           \
            ,"u2nn_topx"           \
            ,"u2nn_topy"           \
            ,"u2nn_topz"           \
            ,"u2nn_tipx"           \
            ,"u2nn_tipy"           \
            ,"u2nn_tipz"           \
            ,"u3nn_pot"            \
            ,"u3nn_potparam"       \
            ,"u3nn_pot2"           \
            #,"u3nn_topx"           \
            ]
            contribmust = [
            #'u1nn_pot' \
            #,"u1nn_potparam"       \
            #,"u1nn_pot2"           \
             "u1nn_topx"           \
            ,"u1nn_topy"           \
            ,"u1nn_topz"           \
            ,"u1nn_tipx"           \
            ,"u1nn_tipy"           \
            #,"u1nn_tipz"           \
            #,"u2nn_pot"            \
            #,"u2nn_potparam"       \
            #,"u2nn_pot2"           \
            #,"u2nn_topx"           \
            #,"u2nn_topy"           \
            #,"u2nn_topz"           \
            #,"u3nn_pot"            \
            #,"u3nn_potparam"       \
            ]
            with printoptions(precision=3, suppress=True):
                #print "u1nn_pot   :",eval("u1nn_pot"),"\tu1nn_potparam:",eval("u1nn_potparam"),"  u1nn_pot2 :",eval("u1nn_pot2"),"\tu1nn_pot2param:",eval("u1nn_pot2param")
                print("u1nn_pot   :",u1nn_pot,"\tu1nn_potparam:",np.array(u1nn_potparam),"  u1nn_pot2 :",u1nn_pot2,"\tu1nn_pot2param:",np.array(u1nn_pot2param))
                print("u1nn_potadd:",u1nn_potadd,"\tu1nn_potaddparam:",np.array(u1nn_potaddparam)) #,"  u1nn_pot2add :",u1nn_pot2add,"\tu1nn_pot2addparam:",u1nn_pot2addparam
                print("u2nn_pot   :",eval("u2nn_pot"),"\tu2nn_potparam:",eval("u2nn_potparam"))
            for i in contribmust:
                if type(eval(i)) == bool and i not in contribmust:
                    continue
                else:
                    print(i,":",eval(i))
            if False:
                print("kkk---------------")
                print("u1nn_pot             :",u1nn_pot,"u1nn_potparam        :",u1nn_potparam)
                print("u1nn_pot2            :",u1nn_pot2,"u2nn_potparam        :",u2nn_potparam)
                print("u1nn_topx            :",u1nn_topx)
                print("u1nn_topy            :",u1nn_topy)
                print("u1nn_topz            :",u1nn_topz)
                print("u1nn_tipx            :",u1nn_tipx)
                print("u1nn_tipy            :",u1nn_tipy)
                print("u1nn_tipz            :",u1nn_tipz)
                print("")
                print("u2nn_pot             :",u2nn_pot)
                print("u2nn_potparam        :",u2nn_potparam)
                print("u2nn_pot2            :",u2nn_pot2)
                print("u2nn_topx            :",u2nn_topx)
                print("u2nn_topy            :",u2nn_topy)
                print("u2nn_topz            :",u2nn_topz)
                print("u2nn_tipx            :",u2nn_tipx)
                print("u2nn_tipy            :",u2nn_tipy)
                print("u2nn_tipz            :",u2nn_tipz)
                print("u3nn_pot             :",u2nn_pot)
                #if types:
            if False:
                print("")
                print("u1nn_pot             :",type(u1nn_pot              ))
                print("u1nn_potparam        :",type(u1nn_potparam         ))
                print("u1nn_pot2            :",type(u1nn_pot2             ))
                print("u1nn_topx            :",type(u1nn_topx             ))
                print("u1nn_topy            :",type(u1nn_topy             ))
                print("u1nn_tipx            :",type(u1nn_tipx             ))
                print("u1nn_tipy            :",type(u1nn_tipy             ))
                print("")
                print("u2nn_pot             :",type(u2nn_pot              ))
                print("u2nn_potparam        :",type(u2nn_potparam         ))
                print("u2nn_pot2            :",type(u2nn_pot2             ))
                print("u2nn_topx            :",type(u2nn_topx             ))
                print("u2nn_topy            :",type(u2nn_topy             ))
                print("u2nn_tipx            :",type(u2nn_tipx             ))
                print("u2nn_tipy            :",type(u2nn_tipy             ))
                print("u3nn_pot             :",type(u3nn_pot              ))
                print("u3nn_potparam        :",type(u3nn_potparam         ))
                print("u3nn_pot2            :",type(u3nn_pot2             ))
            print("---------------------------------------------------------------------------")
            return


    # add pot to pos2
    #if type(pot2) != bool:
    #    pot2[:0] = [u1nn_pot]



    positions = np.copy(crystal.rcar)
    #print "pos:"
    #for idx,i in enumerate(positions):
    #    print idx,i


    # Dies wollen wir eigentlich fuer jede nachbarschale definieren
    energylong = np.zeros((positions.shape[0],1))
    forceslong = np.zeros((positions.shape[0],3))

    energytrantopx = np.zeros((positions.shape[0],1))
    forcestrantopx = np.zeros((positions.shape[0],3))
    energytrantopy = np.zeros((positions.shape[0],1))
    forcestrantopy = np.zeros((positions.shape[0],3))
    energytrantopz = np.zeros((positions.shape[0],1))
    forcestrantopz = np.zeros((positions.shape[0],3))

    energytrantipx = np.zeros((positions.shape[0],1))
    forcestrantipx = np.zeros((positions.shape[0],3))
    energytrantipy = np.zeros((positions.shape[0],1))
    forcestrantipy = np.zeros((positions.shape[0],3))
    energytrantipz = np.zeros((positions.shape[0],1))
    forcestrantipz = np.zeros((positions.shape[0],3))

    forcestrantopxcheck = np.zeros((positions.shape[0],positions.shape[0],3))
    forcestrantopycheck = np.zeros((positions.shape[0],positions.shape[0],3))
    forcestrantopzcheck = np.zeros((positions.shape[0],positions.shape[0],3))
    forcestrantipxcheck = np.zeros((positions.shape[0],positions.shape[0],3))
    forcestrantipycheck = np.zeros((positions.shape[0],positions.shape[0],3))
    forcestrantipzcheck = np.zeros((positions.shape[0],positions.shape[0],3))

    ###############################################################################
    # Diese schleife laeuft einmal fuer die undisplaced structure durch und dann fuer die
    # displaced structure
    # dies ist definiert in get_energy_forces
    ###############################################################################
    #print ""
    #print ""

    ################### save stuff if necessary
    if save_vecs_to_file_for_DOS == True:
        if os.path.isfile("tvecoutall.dat") == True:
            tvecoutall = np.loadtxt("tvecoutall.dat")
        else:
            tvecoutall = np.array([[0.0],[0.0]])
        if tvecoutall.shape[0] == 2:
            tvecoutall = np.array([[0.0],[0.0]])
        tvecoutall = tvecoutall.reshape((len(tvecoutall),1))


        if os.path.isfile("tvecinsall.dat") == True:
            tvecinsall = np.loadtxt("tvecinsall.dat")
        else:
            tvecinsall = np.array([[0.0],[0.0]])
        if tvecinsall.shape[0] == 2:
            tvecinsall = np.array([[0.0],[0.0]])
        tvecinsall = tvecinsall.reshape((len(tvecinsall),1))

        if os.path.isfile("tveclonall.dat") == True:
            tveclonall = np.loadtxt("tveclonall.dat")
        else:
            tveclonall = np.array([[0.0],[0.0]])
        if tveclonall.shape[0] == 2:
            tveclonall = np.array([[0.0],[0.0]])
        tveclonall = tveclonall.reshape((len(tveclonall),1))

        #if os.path.isfile("tveclonallf.dat") == True:
        #    tveclonallf = np.loadtxt("tveclonallf.dat")
        #else:
        #    tveclonallf = np.array([[0.0],[0.0]])
        #if tveclonallf.shape[0] == 2:
        #    tveclonallf = np.array([[0.0],[0.0]])
        #tveclonallf = tveclonallf.reshape((len(tveclonallf),1))

        #if os.path.isfile("tveclonalle.dat") == True:
        #    tveclonalle = np.loadtxt("tveclonalle.dat")
        #else:
        #    tveclonalle = np.array([[0.0],[0.0]])
        #if tveclonalle.shape[0] == 2:
        #    tveclonalle = np.array([[0.0],[0.0]])
        #tveclonalle = tveclonalle.reshape((len(tveclonalle),1))

        if os.path.isfile("tvecdiffall.dat") == True:
            tvecdiffall = np.loadtxt("tvecdiffall.dat")
        else:
            tvecdiffall= np.array([[0.0],[0.0]])
        if tvecdiffall.shape[0] == 2:
            tvecdiffall= np.array([[0.0],[0.0]])
        tvecdiffall= tvecdiffall.reshape((len(tvecdiffall),1))
        print("tvecdiffall1:",tvecdiffall,"save_vecs_to_file_for_DOS 1:",save_vecs_to_file_for_DOS)



    ####################################################################################
    # aeusere schleife ueber alle atome
    ####################################################################################

    for indi,posi in enumerate(positions):
        # get nearest neighbors of this atom
        # To make this faster we have to ensure that here we know the coord_cart and coord_rrel
        # we dont want to load those in every time
        NNlist1 = crystal0.get_NNlist(indi, NNshell,
                cell = crystal0.cellvec,
                coord_rrel = crystal0neverchange.rrel,
                return_NNdist = return_NNdist)
        NNlist2 = crystal0.get_NNlist(indi, 2,
                cell = crystal0.cellvec,
                coord_rrel = crystal0neverchange.rrel,
                return_NNdist = return_NNdist)
        NNlist3 = crystal0.get_NNlist(indi, 3,
                cell = crystal0.cellvec,
                coord_rrel = crystal0neverchange.rrel,
                return_NNdist = return_NNdist)
        NNlist4 = crystal0.get_NNlist(indi, 4,
                cell = crystal0.cellvec,
                coord_rrel = crystal0neverchange.rrel,
                return_NNdist = return_NNdist)


        if positions.shape[0] == 32 and NNlist1.shape[0] != 12:
            sys.exit("did not find 12 NN atoms in presumably fcc structure")
        if positions.shape[0] == 108 and NNlist1.shape[0] != 12:
            sys.exit("did not find 12 NN atoms in presumably fcc structure")
        if positions.shape[0] == 108 and NNlist2.shape[0] != 6:
            sys.exit("did not find 6 2NN atoms in presumably fcc structure")
        if positions.shape[0] == 108 and NNlist3.shape[0] != 24:
            sys.exit("did not find 24 3NN atoms in presumably fcc structure")
        crystal.center_atoms_around_atom(indi,coord_cart=crystal.rcar,cell=crystal.cellvec)
        #print crystal.rcar
        ########################################################################
        # schleife ueber jedes atom
        # wobei hier nur gehandelt wird wenn es sich um erste nachbarn handelt
        ########################################################################
        NNlist = np.append(NNlist1,NNlist2)  # erste und zweite Nachbarn
        NNlist = np.append(NNlist,NNlist3)   # erste und zweite Nachbarn + dritte Nachbarn
        NNlist = NNlist1                     # nur die ersten Nachbarn
        if verbose:
            print("NNlist       :",NNlist)
        ####################################################################################
        # innere schleife ueber naechste nachbar atome (je nachdem welche 1NN oder 1NN+2NN ...)
        ####################################################################################
        for indj in NNlist:

            # bestimme potential erst hier da pot pot2 wechseln koennen
            posj = crystal.rcar[indj]
            print("indi",indi,'>indj:',indj,"posi",posi,"posj",posj,"norm",np.linalg.norm(posj))
            u_pot          = None
            u_potparam     = None
            u_potadd       = None
            u_potaddparam  = None
            u_pot2         = None
            u_pot2param    = None
            u_pot2add      = None
            u_pot2addparam = None
            u_topx         = None
            u_topy         = None
            u_topz         = None
            u_tipx         = None
            u_tipy         = None
            u_tipz         = None
            nn             = None

            if indj in NNlist1:
                u_pot          = u1nn_pot
                u_potparam     = u1nn_potparam
                u_potadd       = u1nn_potadd
                u_potaddparam  = u1nn_potaddparam
                u_pot2         = u1nn_pot2
                u_pot2param    = u1nn_pot2param
                u_pot2add      = u1nn_pot2add
                u_pot2addparam = u1nn_pot2addparam
                u_topx         = u1nn_topx
                u_topy         = u1nn_topy
                u_topz         = u1nn_topz
                u_tipx         = u1nn_tipx
                u_tipy         = u1nn_tipy
                u_tipz         = u1nn_tipz
                nn             = 1

            if indj in NNlist2:
                u_pot          = u2nn_pot
                u_potparam     = u2nn_potparam
                u_potadd       = ""
                u_potaddparam  = ""
                u_pot2         = u2nn_pot2
                u_pot2param    = u2nn_potparam2
                u_pot2add      = ""
                u_pot2addparam = ""
                u_topx         = u2nn_topx
                u_topy         = u2nn_topy
                u_topz         = u2nn_topz
                u_tipx         = u2nn_tipx
                u_tipy         = u2nn_tipy
                u_tipz         = u2nn_tipz
                nn             = 2

            if indj in NNlist3:
                u_pot          = u3nn_pot
                u_potparam     = u3nn_potparam
                u_potadd       = ""
                u_potaddparam  = ""
                u_pot2         = u3nn_pot2
                u_pot2param    = u3nn_potparam2
                u_pot2add      = ""
                u_pot2addparam = ""
                u_topx         = False #u2nn_topx
                u_topy         = False #u2nn_topy
                u_topz         = False #u2nn_topz
                u_tipx         = False #u2nn_tipx
                u_tipy         = False #u2nn_tipy
                u_tipz         = False #u2nn_tipz
                nn             = 3

            # der folgende check sollte doch nicht jedes mal gemacht werden ....
            if u_pot          == None: sys.exit('pot not defined 1')
            if u_potparam     == None: sys.exti('pot not defined 2')
            #print "u_pot2:",u_pot2
            if u_pot2         == None: sys.exti('pot not defined 3')
            if u_pot2param    == None: sys.exti('pot not defined 2')
            if u_topx         == None: sys.exti('pot not defined 4')
            if u_topy         == None: sys.exti('pot not defined 5')
            if u_topz         == None: sys.exti('pot not defined 6')
            if u_tipx         == None: sys.exti('pot not defined 7')
            if u_tipy         == None: sys.exti('pot not defined 8')
            if u_tipz         == None: sys.exti('pot not defined 9')
            if nn             == None: sys.exti('pot not defined nn')
            #print "u_potadd:",u_potadd
            if u_potadd       == None: sys.exit('pot not defined 11')
            if u_potaddparam  == None: sys.exti('pot not defined 12')
            if u_pot2add         == None: sys.exti('pot not defined 13')
            if u_pot2addparam    == None: sys.exti('pot not defined 14')

            if type(u_pot) == bool and \
                type(u_potparam) == bool and \
                type(u_pot2) == bool and \
                type(u_pot2param) == bool and \
                type(u_topx) == bool and \
                type(u_topy) == bool and \
                type(u_topz) == bool and \
                type(u_tipx) == bool and \
                type(u_tipy) == bool and \
                type(u_tipz) == bool:
                continue

            #print "nn:",nn,"indi:",indi,"indj:",indj,"u_pot:",u_pot,"u_potparam:",u_potparam


            ##################################################################
            ##################################################################
            ################# LONGVEC ########################################
            # braucht umbedingt         : posj, indj, u_pot, u_potparam
            # braucht nicht umbedingt   : u_pot2, upot2param,
            ##################################################################
            ##################################################################
            longvec = np.copy(posj)
            longvecnorm = np.copy(np.linalg.norm(posj))
            vec0 = crystal0.rcar[indj]
            vec0norm = np.copy(np.linalg.norm(vec0))
            #print longvecnorm-vec0norm,posj,vec0,posj[0],posj[1],posj[2],vec0[0],vec0[1],vec0[2] #longvecnorm,vec0norm

            if type(u_pot) != bool or \
                type(u_potparam) != bool or \
                type(u_pot2) != bool or \
                type(u_pot2param) != bool:


                #if type(u_pot2) != bool or type(u_pot2param) != bool:  # u_pot2 == m or mc1
                if type(u_pot2param) != bool:  # u_pot2 == m or mc1
                    if u_pot2 != 'm' and u_pot2 != 'mc1' and u_pot2 != 'poly':
                        print("u_pot2:",u_pot2)
                        sys.exit("u_pot2 has to be m or mc1 or poly")
                    #print "u_pot2;",u_pot2
                    nndist = float(u_potparam[2])  # here [2] is the third object since it is a numpy array (not a list)
                    #print "nn:",nndist
                    if longvecnorm > nndist:
                        u_pot = u_pot2
                        u_potparam = u_pot2param
                        u_potadd = u_pot2add
                        u_potaddparam = u_pot2addparam
                        #getfrompot2 = True
                # np.linalg.norm(vec0 - longvec) != np.linalg.norm(longvec - vec0)
                # The reason is that out EqCoords are not exact: the positions are
                # 0.0000000000
                # 0.1666666667
                # 0.3333333333
                # 0.5000000000
                # 0.6666666667
                # 0.8333333333
                # This in the end "only" gives an accuracy of "only" when doing XXX) -8.26000601251e-10 and 4.13000300625e-10; The difference is ~4e-10
                # THEREFORE: DONT USE THE FOLLOWING PART!
                # if longvecnorm == vec0norm:
                #    elong = 0.0
                #    flongnorm = 0.0
                #    flong = np.array([0.0, 0.0, 0.0])
                #else:
                # DONT USE THE PART BEFORE!

                elong, flong = getefvec(longvec,u_potparam,pot = u_pot)
                elongadd, flongadd = getefvec(longvec,u_potaddparam,pot = u_potadd)

                forceslong[indi] = forceslong[indi] + flong + flongadd
                energylong[indi] = energylong[indi] + elong + elongadd
                #if indi == 1:
                #    sys.exit()
                #if elong > 0.1:
                #    print "longvec:",longvec,"norm:",np.linalg.norm(longvec)
                #    print "u_potparam:",u_potparam
                #    print "u_pot:",u_pot
                #    print "elong, flong:",elong,flong
                #    print "elongadd,flongadd:",elongadd,flongadd
                #    sys.exit()

                if save_vecs_to_file_for_DOS == True:
                    if np.linalg.norm(longvec) > 0.00001:
                        #print "shape;",tveclonall.shape
                        tveclonall = np.concatenate((tveclonall,np.array([[np.linalg.norm(longvec)]])),axis=0)
                        #tveclonallf = np.concatenate((tveclonallf,np.array([[np.linalg.norm(forceslong[indi])]])),axis=0)
                        #tveclonalle = np.concatenate((tveclonalle,np.array([[np.linalg.norm(elong)]])),axis=0)


                #if abs(elong) > 0.00001:
                #print "indi:",indi,"nn:",nn,"indj:",indj,"longvecnorm:",longvecnorm,"longvec:",longvec,"flongnorm:",flongnorm,"elong:",elong
                #if indi == 2:
                #    sys.exit()


            #############################################################
            #############################################################
            ############## TVEC #########################################
            #############################################################
            #############################################################

            ###########################################################
            # for fcc:
            ###########################################################
            # am besten waere hier wenn man tvecout:
            #   1) so hindrehte wie die kraefte gemappt wruden
            #   2) den entsprechenden kraftvektor berechnet
            #   3) die drehung von 1) wieder rueckgaengig macht
            #   e.g.:
            #   for first NN tox we always want to rotate the current tout vector on [1,0,0]
            #   it will be necessary to rotate from {-+1,0,0}   - > 1,0,0
            #                                       {  0,-+1,0} - > 1,0,0  (z stays, xy swap)
            #                                       {  0,0,-+1} - > 1,0,0  (y stays, xz swap)
            #   x was already mapped to correct vector with topdirectino wich is now topdirectionx
            # in the end this are all just rotation matrizes and can be implemented
            # generally for every neighbor for every direction
            if  type(u_topx) != bool or \
                type(u_topy) != bool or \
                type(u_topz) != bool or \
                type(u_tipx) != bool or \
                type(u_tipy) != bool or \
                type(u_tipz) != bool:

                ##################################################################
                ############### TVEC ###############
                tvec = longvec - vec0
                ############### TVEC ###############
                ##################################################################

                topdirectionx = np.array([0.0,0.0,0.0])  # those are the directions in which the forces act on to / ti part of vector
                topdirectiony = np.array([0.0,0.0,0.0])
                topdirectionz = np.array([0.0,0.0,0.0])

                tipdirectionx = np.array([0.0,0.0,0.0])
                tipdirectiony = np.array([0.0,0.0,0.0])
                tipdirectionz = np.array([0.0,0.0,0.0])

                if vec0[0] == 0.0 and nn == 1:  # dies gilt nur fuer die ersten Nachbarn, fuer die 2NN
                                    # muss dies anders definiert werden
                                    # hier sollte man sich mehr gedanken ueber symmetrien machen
                                    # fuer die 2NN werden die ti und to gleich sein
                    # z.b. vec0 = [0.0,  2.065,  2.065] 1NN
                    # z.b. vec0 = [0.0,  2.065, -2.065] 1NN
                    # z.b. vec0 = [0.0, -2.065,  2.065] 1NN
                    # z.b. vec0 = [0.0, -2.065, -2.065] 1NN
                    #topdirectionx = np.array([vec0[2]/abs(vec0[1]),0.0,0.0])  # only the direction seems to be important not the sign
                    topdirectionx = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
                    topdirectiony = np.array([0.0,-vec0[1]/abs(vec0[1]),0.0])
                    topdirectionz = np.array([0.0,0.0,-vec0[2]/abs(vec0[2])])

                    # for tip directions it is probably best to project on vec0? or on 1/0/0 .... this should only be a matter of koordinate transformation.
                    # therefore lets start with x/y directions as done for the tox toy directions
                    #
                    # tipdirection{x,y,z} and tipdirection{long,senk} are both possible bases, this is just a matter of coordinate transformation
                    tipdirectionx = np.array([0.0,-vec0[1]/abs(vec0[1]),0.0])  # evtl auch hier -vec0[xyz]/abs(vec0)
                    tipdirectiony = np.array([0.0,0.0,-vec0[2]/abs(vec0[2])])
                    tipdirectionz = np.array([0.0,0.0,0.0])  # dieser vector wird nicht gebraucht da er senkrecht auf tipdirection steht

                    tipdirectionlong = vec0
                    tipdirectionsenk = np.array([0.0,vec0[1]/abs(vec0[1]),-vec0[2]/abs(vec0[2])])



                if vec0[1] == 0.0 and nn == 1:  # dies gilt nur fu
                    # z.b. vec0 = [ 2.065, 0.0,  2.065] 1NN
                    # z.b. vec0 = [ 2.065, 0.0, -2.065] 1NN
                    # z.b. vec0 = [-2.065, 0.0,  2.065] 1NN
                    # z.b. vec0 = [-2.065, 0.0, -2.065] 1NN
                    topdirectionx = np.array([0.0,1.0,0.0])
                    topdirectiony = np.array([-vec0[0]/abs(vec0[0]),0.0,0.0])
                    topdirectionz = np.array([0.0,0.0,-vec0[2]/abs(vec0[2])])
                    #tipdirection_out = np.array([vec0[0],0.0,-vec0[2]])
                    #tipdirection_in = vec0
                    tipdirectionx = np.array([-vec0[0]/abs(vec0[0]),0.0,0.0])  # evtl auch hier -vec0[xyz]/abs(vec0)
                    tipdirectiony = np.array([0.0,0.0,1.0])
                    tipdirectionz = np.array([0.0,1.0,0.0])  # dieser vector wird nicht gebraucht da er senkrecht auf tipdirection steht

                    tipdirectionlong = vec0
                    # rechtsdrehend
                    if vec0[0] > 0.0 and vec0[2] < 0.0: tipdirectionsenk = np.array([-abs(vec0[0]),0.0,-abs(vec0[2])])
                    if vec0[0] < 0.0 and vec0[2] > 0.0: tipdirectionsenk = np.array([ abs(vec0[0]),0.0, abs(vec0[2])])
                    if vec0[0] < 0.0 and vec0[2] < 0.0: tipdirectionsenk = np.array([-abs(vec0[0]),0.0, abs(vec0[2])])
                    if vec0[0] > 0.0 and vec0[2] > 0.0: tipdirectionsenk = np.array([ abs(vec0[0]),0.0,-abs(vec0[2])])
                    #tipdirectionsenk = np.array([vec0[0]/abs(vec0[0]),0.0,-vec0[2]/abs(vec0[2])])

                if vec0[2] == 0.0 and nn == 1:
                    # z.b. vec0 = [ 2.065,  2.065, 0.0] 1NN
                    # z.b. vec0 = [ 2.065, -2.065, 0.0] 1NN
                    # z.b. vec0 = [-2.065,  2.065, 0.0] 1NN
                    # z.b. vec0 = [-2.065, -2.065, 0.0] 1NN
                    topdirectionx = np.array([0.0,0.0,1.0])
                    #topdirectiony = np.array([0.0,1.0,0.0])
                    #topdirectionz = np.array([1.0,0.0,0.0])
                    topdirectiony = np.array([0.0,-vec0[1]/abs(vec0[1]),0.0])
                    topdirectionz = np.array([-vec0[0]/abs(vec0[0]),0.0,0.0])

                    tipdirectionx = np.array([-vec0[0]/abs(vec0[0]),0.0,0.0])  # evtl auch hier -vec0[xyz]/abs(vec0)
                    tipdirectiony = np.array([0.0,-vec0[1]/abs(vec0[1]),0.0])  # evtl auch hier -vec0[xyz]/abs(vec0)
                    tipdirectionz = np.array([0.0,0.0,1.0])  # dieser vector wird nicht gebraucht da er senkrecht auf tipdirection steht


                    # bin von atom 87 ausgegangen: [-2.065,  2.065, 0.0] ist ref
                    tipdirectionlong = vec0
                    tipdirectionsenk = np.array([-vec0[0]/abs(vec0[0]),vec0[1]/abs(vec0[1]),0.0])

                    #tipdirection_out = np.array([vec0[0],-vec0[1],0.0])
                    #tipdirection_in = vec0

                # lets think about the 2NN:
                # first of all the longitudinal part should work without problems without tox OK
                # second: tox
                # VEC1
                #if nn == 2:
                #    topdirectionx = np.array([1.0,0.0,0.0])
                #    topdirectiony = np.array([0.0,1.0,0.0])
                #    topdirectionz = np.array([0.0,0.0,1.0])
                #if nn == 3:  # brauchen wir erst mal nicht
                #    topdirectionx = np.array([1.0,0.0,0.0])
                #    topdirectiony = np.array([0.0,1.0,0.0])
                #    topdirectionz = np.array([0.0,0.0,1.0])

                if vec0[0] == 0.0 and nn == 2: # vec in y-z ebene
                    topdirectionx = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
                    if vec0[1] == 0.0 and nn == 2: # vec in y-z ebene
                        topdirectiony = np.array([0.0,1.0,0.0])  # minus xdir  will also be correct
                    if vec0[2] == 0.0 and nn == 2: # vec in y-z ebene
                        topdirectiony = np.array([0.0,0.0,1.0])  # minus xdir  will also be correct

                if vec0[1] == 0.0 and nn == 2: # vec in x-z ebene
                    topdirectionx = np.array([0.0,1.0,0.0])  # minus xdir  will also be correct
                    if vec0[0] == 0.0 and nn == 2: # vec in y-z ebene
                        topdirectiony = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
                    if vec0[2] == 0.0 and nn == 2: # vec in y-z ebene
                        topdirectiony = np.array([0.0,0.0,1.0])  # minus xdir  will also be correct

                if vec0[2] == 0.0 and nn == 2: # vec in x-y ebene
                    topdirectionx = np.array([0.0,0.0,1.0])  # minus xdir  will also be correct
                    if vec0[0] == 0.0 and nn == 2: # vec in y-z ebene
                        topdirectiony = np.array([1.0,0.0,0.0])  # minus xdir  will also be correct
                    if vec0[1] == 0.0 and nn == 2: # vec in y-z ebene
                        topdirectiony = np.array([0.0,1.0,0.0])  # minus xdir  will also be correct

                if nn == 1:
                    if np.linalg.norm(topdirectionx) == 0.0 or np.linalg.norm(topdirectiony) == 0.0 or np.linalg.norm(topdirectionz) == 0.0 or \
                        np.linalg.norm(tipdirectionx) == 0.0 or np.linalg.norm(tipdirectiony) == 0.0: # or np.linalg.norm(tipdirectionz) == 0.0:
                        print("nn:",nn,"indi:",indi,"indj:",indj,"posj:",posj,"|| first ||")
                        print("topdirectionx:",topdirectionx)
                        print("topdirectiony:",topdirectiony)
                        print("topdirectionz:",topdirectionz)
                        print("tipdirectionx:",tipdirectionx)
                        print("tipdirectiony:",tipdirectiony)
                        print("tipdirectionz:",tipdirectionz)
                        sys.exit("t{o,i}pdirection{x,y,z} not found")

                if nn == 2:
                    if np.linalg.norm(topdirectionx) == 0.0 or np.linalg.norm(topdirectiony) == 0.0: # or \
                        #np.linalg.norm(tipdirectionx) == 0.0 or np.linalg.norm(tipdirectiony) == 0.0: # or np.linalg.norm(tipdirectionz) == 0.0:
                        print("nn:",nn,"indi:",indi,"indj:",indj,"posj:",posj,"|| second ||")
                        print("topdirectionx:",topdirectionx)
                        print("topdirectiony:",topdirectiony)
                        print("topdirectionz:",topdirectionz)
                        print("tipdirectionx:",tipdirectionx)
                        print("tipdirectiony:",tipdirectiony)
                        print("tipdirectionz:",tipdirectionz)
                        sys.exit("t{o,i}pdirection{x,y,z} not found")


                ###############################################################
                ###############################################################
                ################# tvecout  ####################################
                ###############################################################
                ###############################################################
                # oldway do first of all not delete
                # tvecout ist immer senkrecht auf longvec in fcc
                # direction of tvecoutricx will always be correct since tvec is only projected
                tvecoutdir = utils.project_vector(tvec,topdirectionx)   # [ 0.3, 0.0, 0.0 ]
                #tvecout = utils.project_vector(utils.reject_vector(longvec, vec0),tvecoutdir)   # hier ist es egal ob man tvecoutdir oder tvecoutdir*2 schreibet, nur die richtung ist von interesse
                # sollte das gleiche liefern
                # the rejection (rejected vector) is the projection of tvec (yes tvec, also it only gets longvec, but longvec+vec0=tvec) ... the projection of tvec on vec0;
                # BUT: if tvec and vec0 have are perpendicular, --> rejection = tvec
                # with the only difference that we take the rejection which somewhat longer/shorter in length and not the projection
                rejection = utils.reject_vector(longvec, vec0)
                projection = utils.project_vector(longvec, vec0)
                tvecout = utils.project_vector(rejection,topdirectionx)   # hier ist es egal ob man tvecoutdir oder tvecoutdir*2 schreibet, nur die richtung ist von interesse

                if save_vecs_to_file_for_DOS == True:
                    if np.linalg.norm(tvecout) > 0.00001:
                        #print "shape;",tvecoutall.shape
                        tvecoutall = np.concatenate((tvecoutall,np.array([[np.linalg.norm(tvecout)]])),axis=0)
                        tvecoutall = np.concatenate((tvecoutall,np.array([[-np.linalg.norm(tvecout)]])),axis=0)

                # # ######################################################################################
                # # # rejection = utils.reject_vector(longvec, vec0)
                # # # utils.project_vector(rejection,AXIS) is always necessary if we pant to project the part of longvec on AXIS (where AXIS is perpendicular to vec0)
                # # ######################################################################################
                # # tvecoutnorm = np.linalg.norm(tvecout)

                # # tvecoutx = tvecout
                # # tvecouty = tvecoutnorm*topdirectiony
                # # tvecoutz = tvecoutnorm*topdirectionz

                # # ####################################################################
                # # # TOX, TOY, TOZ
                # # # this does not have to be done if atoms are at NN positions, we will nevertheless do it to get crrect forces to 0
                # # ####################################################################
                # # # to spilt this in x,y,z is absolutely redundant; tvecout{x,y,z} are the same vector showing in different direction; this is just
                # # # a poor mans implementation to get the forces back in the correct direction
                # # # oldway
                # # etopx,ftopx = getefvec(tvecoutx,u_topx,pot=False)
                # # etopy,ftopy = getefvec(tvecouty,u_topy,pot=False)
                # # etopz,ftopz = getefvec(tvecoutz,u_topz,pot=False)
                # # ftopxx = np.copy(ftopx)
                # # ftopyy = np.copy(ftopy)
                # # ftopzz = np.copy(ftopz)
                # # #print ">>> indi:",indi,indj,"posi:",posi,"posj:",posj,"vec0:",vec0,"longvec:",longvec,np.linalg.norm(longvec),"flong:",flong
                # # #print "  OLD: tvec:",tvec," -> tvecoutx,tvecouty,tvecoutz:",tvecoutx,tvecouty,tvecoutz
                # # #print "  OLD: ftopx,ftopy,ftopz:",ftopx,ftopy,ftopz


                #print "!!!!!!!!fold:",ftopxx,ftopyy,ftopzz
                atomlist=[8]
                atomlist = False
                #if indi == 27: atomlist = [27]
                #if indi == 35: atomlist = [35]
                #atomlist = [0, 27, 29, 33, 35,     54, 56, 72, 74,    81, 87, 99, 105]
                #bbbtext = ">>> indi:",indi,indj,"posi:",posi,"posj:",posj,"vec0:",vec0,"longvec:",longvec,np.linalg.norm(longvec),"flong:",flong,"tvec:",tvec
                if nn == 1:
                    vec0_orig_formapping = np.array([0.0, -2.065, -2.065])
                    mappedvec_orig_dir = np.array([1.0, 0.0, 0.0])
                    etopx,etopy,etopz,ftopx,ftopy,ftopz = \
                        rotate_vec0_to_vec0_original_inplane_and_getforce_current(vec0,longvec,vec0_orig_formapping,mappedvec_orig_dir,u_topx,u_topy,u_topz,tvecwirkdir = topdirectionx,atomlist=atomlist, atomi=indi,atomj=indj, bbbtext="TO{X,Y,Z}")
                if nn == 2:
                    vec0_orig_formapping = np.array([0.0, -4.13, 0.0])
                    mappedvec_orig_dir = np.array([1.0, 0.0, 0.0])
                    etopx,etopy,etopz,ftopx,ftopy,ftopz = \
                        rotate_vec0_to_vec0_original_inplane_and_getforce_current(vec0,longvec,vec0_orig_formapping,mappedvec_orig_dir,u_topx,u_topy,u_topz,tvecwirkdir = topdirectionx,atomlist=atomlist, atomi=indi,atomj=indj, bbbtext="TO{X,Y,Z}")


                #print "fx,fy,fz         :",fx,fy,fz
                #print "ftopx,ftopy,ftopz:",ftopx,ftopy,ftopz
                #print ""
                #print "||",fx+fy+fz,ftopx+ftopy+ftopz

                # # check if old way and new way are equivalent
                # if abs((fx+fy+fz)[2] - (ftopx+ftopy+ftopz)[2]) >= 0.000001:
                #     print fx+fy+fz
                #     print ftopx+ftopy+ftopz
                #     sys.exit("x not 0")
                # if abs((fx+fy+fz)[0] - (ftopx+ftopy+ftopz)[0]) >= 0.000001:
                #     print fx+fy+fz
                #     print ftopx+ftopy+ftopz
                #     sys.exit("x not 0")
                # if abs((fx+fy+fz)[1] - (ftopx+ftopy+ftopz)[1]) >= 0.000001:
                #     print fx+fy+fz
                #     print ftopx+ftopy+ftopz
                #     sys.exit("x not 0")

                #print "tvecoutx:",tvecoutx
                #print "tvecouty:",tvecouty
                #print "tvecoutz:",tvecoutz



                ###############################################################
                ###############################################################
                ################# tvecinsenk  #################################
                ###############################################################
                ###############################################################
                ftipx = np.array([0.0,0.0,0.0])
                ftipy = np.array([0.0,0.0,0.0])
                ftipz = np.array([0.0,0.0,0.0])
                etipx = np.array([0.0])
                etipy = np.array([0.0])
                etipz = np.array([0.0])

                mappedvec_orig_dir = np.array([1.0, 1.0, 0.0])   # if mapped on quer/parallel
                tvecinsenkdir = utils.project_vector(tvec,tipdirectionsenk)

                #mappedvec_orig_dir = np.array([1.0, 0.0, 0.0])   # if mapped on quer/parallel
                #tvecinsenkdir = utils.project_vector(tvec,tipdirectionx)

                tvecinsenk = utils.project_vector(utils.reject_vector(longvec, vec0),tvecinsenkdir) # this should remove the forces on 87 and 99 in quer direction
                #tvecinsenk = np.array([math.copysign(tvecinsenk[0],tvecinsenkdir[0]),math.copysign(tvecinsenk[1],tvecinsenkdir[1]),math.copysign(tvecinsenk[2],tvecinsenkdir[2])])

                # lets try to map everythin on xy direction:
                #mappedvec_orig_dir = np.array([1.0, 0.0, 0.0])   # if mapped on x,y
                #tvecinsenkdir = utils.project_vector(tvec,tipdirectionx)
                #tvecinsenk = utils.project_vector(utils.reject_vector(longvec, vec0),tvecinsenkdir) # this should remove the forces on 87 and 99 in quer direction
                ################################################
                # Anharmonic approximation a)
                # if we want to stay in "close" to harmonic approximation we can just have one ti mode (ohne long, one to and one ti - corresponding to sample x,y,z direction in harmonic approx) / ==
                # the nice thing is that we could show explicitly what the nonlinearity of forces can do
                #
                #
                # Anharmonic approximation b)
                # This clearly needs angular terms. Can we argue to use those in an "anharmonic" approximation?
                #
                if np.linalg.norm(tvecinsenk) > 0.00001:
                    #print "tvecinsenk:",tvecinsenk,np.linalg.norm(tvecinsenk)
                    pass
                if save_vecs_to_file_for_DOS == True:
                    if np.linalg.norm(tvecinsenk) > 0.00001:
                        tvecinsall = np.concatenate((tvecinsall,np.array([[np.linalg.norm(tvecinsenk)]])),axis=0)
                        tvecinsall = np.concatenate((tvecinsall,np.array([[-np.linalg.norm(tvecinsenk)]])),axis=0)


                # from here we need only to proceed if np.linalg.norm(tvecinsenkdir) is != 0.0  (in this case the displacement simply has no inplane quer direction)
                atomlist=[8,16]
                atomlist=False
                atomlist=[26]
                atomlist=False
                atomlist=[20,16,26]
                atomlist=False

                if np.linalg.norm(tvecinsenk) != 0.0:  # keep this to avoid problems with rotation matrix
                    vec0_orig_formapping = np.array([-2.065, 2.065, 0.0])   # welches atom haben wir uns angeschaut? 87
                    etipx,etipy,etipz,ftipx,ftipy,ftipz = \
                        rotate_vec0_to_vec0_original_inplane_and_getforce_current(vec0,longvec,vec0_orig_formapping,mappedvec_orig_dir,u_tipx,u_tipy,u_tipz,tvecwirkdir = tvecinsenk,atomlist=atomlist, atomi=indi,atomj=indj,bbbtext="TI{x,y,z}")


                ###############################################################
                ###############################################################
                ################# tvecinpar ###################################
                ###############################################################
                ###############################################################
                tvecinpar = projection - vec0

                # the following two lines only work if tvecinsenk would only show in one direction with correct norm like: [ 1.2, 0.0, 0.0 ] or [ 0.0, 1.2, 0.0 ] or [ 0.0, 0.0, 1.2 ]
                # as is the case for tvecout{x,y,z}
                #etipsenkx, ftipsenkx = getefvec(tvecinsenk,u_tipx,pot=False)   # x is here the taken on basis of 1NN dir
                #etipsenky, ftipsenky = getefvec(tvecinsenk,u_tipy,pot=False)   # x is here the taken on basis of 1NN dir

                #tvecinsenkx = tvecinsenk
                #tvecinsenky = tvecinsenknorm*tipdirectiony
                #tvecinsenkz = tvecinsenknorm*tipdirectionz

                #tvecin = tvec - tvecout
                #tvecinnorm = np.linalg.norm(tvecin)
                #tvecout = tvecout*-1  # currently tvecout shows towards plane but we need the opposite direction for the forces to be intuitive



                ####################################################################
                # TIX, TIY, TIZ
                ####################################################################
                # makes RuntimeWarning: invalid value encountered in divide
                #etip_outx =     Trans135(tvecin_outnorm,-0.0539196,0.0,0.0) # 6x6x6kp
                #ftip_outx = Trans135_der(tvecin_outnorm,-0.0539196,0.0,0.0) # 6x6x6kp
                #if abs(np.linalg.norm(ftopx)) <= 0.0000000001:
                #    ftip_outx = np.array([0.0, 0.0, 0.0])
                #else:
                #    ftip_outx = tvecin_out/tvecin_outnorm*ftip_outx


                #if abs(np.linalg.norm(ftopx)) <= 0.0000000001:
                #    ftip_outy = np.array([0.0, 0.0, 0.0])
                #else:
                #    ftip_outy = tvecin_out/tvecin_outnorm*ftip_outy


                ####################################################################################
                # MMM
                # MMM
                # Bei Anna: mache paar auslenkungen in -y ,y ,-z, z und anschauen kraefte auf die 1NN 2NN, sind symmetrien korrekt?
                # Es ist der hammer, durch die endiche anzahl der in den EqCoords machen wir in den 1NN tox atomen einen Fehler von 0.00000001 in den Kraeften
                #
                # ftop_y macht es zurzeit schlechter welche moeglichkeiten:
                #       a)  die implementation mit vec0transout ist nicht korrekt fuer andere
                #           ebenen (ueberpruefbar durch x auslenkung in (0,1,0),(0,-1,0)...
                #           und anschauen des atoms auf welches die out Kraft wirkt)
                #       b)  fuer andere als das 9te atom ist diese implementation zwar richtig
                #           aber macht ohne das einfuehren der inplane moden keinen sinn
                #       c)  go to a bigger supercell to avoid effects of repetitions in supercell
                ####################################################################################
                # /Users/glensk/Dropbox/Understand_distributions/Al/Al_displacements/Al_displacements_3x3x3sc_xdir_4x4x4kp/4.13Ang_0.3
                #
                # die restkraft ist noch   [-0.23    0.11    0.11] auf at   [ 0.      2.07    2.07]
                # sollte aber ganz auf null sein. Es stimmt noch etwas nicht mit toy und toz
                # VORGEHEN:
                #   check if np.roll always works:
                #       do displacements in -x +y -y +z -z direction
                #
                #
                #

                #print "ftopx,y,z:",ftopx,ftopy,ftopz
                #print "ftipx,y,z:",ftipx,ftipy,ftipz
                #print "tiex,tiey,z:",tiex,tiey,tiez
                if np.abs(ftopx).max() > 0.00000000001:
                    print("ftopx,y,z:",ftopx,ftopy,ftopz)
                    #sys.exit()

                forcestrantopx[indi] = forcestrantopx[indi] + ftopx
                energytrantopx[indi] = energytrantopx[indi] + etopx
                forcestrantopy[indi] = forcestrantopy[indi] + ftopy
                energytrantopy[indi] = energytrantopy[indi] + etopy
                forcestrantopz[indi] = forcestrantopz[indi] + ftopz
                energytrantopz[indi] = energytrantopz[indi] + etopz


                forcestrantipx[indi] = forcestrantipx[indi] + ftipx
                energytrantipx[indi] = energytrantipx[indi] + etipx
                forcestrantipy[indi] = forcestrantipy[indi] + ftipy
                energytrantipy[indi] = energytrantipy[indi] + etipy
                forcestrantipz[indi] = forcestrantipz[indi] + ftipz
                energytrantipz[indi] = energytrantipz[indi] + etipz

                # checks in case sum of forces is not 0
                if np.linalg.norm(ftopx) != 0.0 or np.linalg.norm(ftopy) != 0.0 or np.linalg.norm(ftopz) != 0.0:
                    pass
                    #print utils.printgreen(indi,indj,"FTOPX:",ftopx,ftopy,ftopz)
                forcestrantopxcheck[indi,indj] = ftopx
                forcestrantopycheck[indi,indj] = ftopy
                forcestrantopzcheck[indi,indj] = ftopz
                forcestrantipxcheck[indi,indj] = ftipx
                forcestrantipycheck[indi,indj] = ftipy
                forcestrantipzcheck[indi,indj] = ftipz

                if True: #False: #verbose == True:
                    # 3x3x3sc 1NN tox atoms:                  2NN tox atoms: hmmm. sind das to or ti atome? koennen beides sein!
                    # tox atoms in xdir: ( no tix in xdir )
                    # 27.   [ 0.      2.07    2.07]           1. [  0.      0.      4.13]
                    # 29.   [ 0.      2.07   10.32]
                    # 33.   [ 0.     10.32    2.07]
                    # 35.   [ 0.     10.32   10.32]

                    # tix atoms in x and quer: (tox komponent in quer auslenkung)
                    # 54.   [ 2.07    0.      2.07]
                    # 56.   [ 2.07    0.     10.32]
                    # 72.   [ 10.32    0.      2.07]
                    # 74.   [ 10.32    0.     10.32]
                    #
                    # tix atoms in x and quer: ( no tox component in x and quer)
                    # 81.   [  2.07    2.07    0.]
                    # 87.   [  2.07   10.32    0.]
                    # 99.   [ 10.32    2.07    0.]
                    # 105.  [ 10.32   10.32    0.]
                    #   how to map tix forces: (two good possibilities in fcc)
                    #       a)  x displacement, get force of 81 (or better 87 to keep the atom
                    #           of the quer displacement), substract long, ( forces easy,
                    #           but the dispvector (x axis of force) in (1,-1,0) has to be
                    #           projected or rejected on (1,-1,0) direction
                    #       b)  DONE; quer displacement, get force on 87, substract long,
                    #           ( forces easy,
                    #           dispvector (x axis of force) in (1,1,0) direction directly there.
                    #           -> b) seems to be easier, check if consistent with a)
                    #
                    sc3x3x3_1nn = np.array([0, 27, 29, 33, 35,     54, 56, 72, 74,    81, 87, 99, 105])
                    sc3x3x3_1nn_xz = np.array([54, 56, 72, 74])
                    sc3x3x3_1nn_to = np.array([0, 27, 29, 33, 35])
                    #sc3x3x3_1nn_to = np.array([0, 27, 29])
                    #sc3x3x3_1nn = np.array([0, 87, 99])
                    sc3x3x3_1nn_ti_xy = np.array([0, 87])
                    sc3x3x3_2nn = np.array([0, 3, 1, 6, 2])  # 18, 9 are the once only
                    sc3x3x3_3nn = np.array([0, 36])  # 18, 9 are the once only
                    #affected by the longvec
                    #sc3x3x3_1nn = np.array([0, 27, 54])
                    # VEC2
                    sc2x2x2_1nn = np.array([0, 16, 17, 24, 26])
                    sc2x2x2_1nn = np.array([0, 24])
                    sc2x2x2_1nn = np.array([0, 8])
                    sc2x2x2_1nn = np.array([8])
                    showvec = sc3x3x3_1nn
                    showvec = sc2x2x2_1nn
                    if False: #True:
                        if indi in showvec:
                            if indj in showvec:
                                #if indi == 8 or indi == 24 or indi == 26:

                                #print "   indi:",indi,indj,posi,"  indj(long):",indj,posj,
                                #"vec0:",vec0,"LONG:",flong,"|tvec:",tvec,tvecout,tvecin_out
                                chop = 0.00000001
                                if abs(np.linalg.norm(tvec)) <= chop: # and abs(np.linalg.norm(tvecout)) <= chop: # and abs(np.linalg.norm(tvecin)) <= chop:
                                    pass
                                else:
                                    print(">>> indi:",indi,indj,"posi:",posi,"posj:",posj,"vec0:",vec0,"longvec:",longvec,np.linalg.norm(longvec),"flong:",flong,"tvec:",tvec,"tvecinsenk:",tvecinsenk) # ,"outx:",tvecoutx
                                    #print "fold:",ftopxx,ftopyy,ftopzz
                                    #print "fneo:",fxorig,fyorig,fzorig
                                    print("ftop:",ftopx,ftopy,ftopz)
                                    print("ftip:",ftipx,ftipy,ftipz)
                                    print("")
                                    print("")
                    #showinfo()

    ##############################################################################
    # ADD HESSE 2nd neighbor !!!!!!!!!! DONT USE IT!
    # INTRODUCES FORCES ON UNDISPLACED STRUCTURE  !!!!!!!!!!!!!!!!!!!!!
    # !!!!!!!!!! DONT USE IT! INTRODUCES FORCES ON UNDISPLACED STRUCTURE  !!!!!!!!!!!!!!!!!!!!!
    # !!!!!!!!!! DONT USE IT! INTRODUCES FORCES ON UNDISPLACED STRUCTURE  !!!!!!!!!!!!!!!!!!!!!
    # for this we need the Hessematrix and Hessematrix_1NN
    ##############################################################################
    hrestmev, hrestev, hrestf = 0.0, 0.0, np.zeros((positions.shape[0],3))
    if type(hessefile) != bool and type(hessefile1nn) != bool:
        hrestmev, hrestev, hrestf = get_energy_forces_hesse_rest(hessefile = hessefile,
        hessefile1nn = hessefile1nn)
        print("############# HRESTF:")
        print(hrestf)


    ##############################################################################
    # SUM UP FORCES / energies
    ##############################################################################
    energy =  \
            np.sum(energylong)/2. +\
            np.sum(energytrantopx)/2. +\
            np.sum(energytrantopy)/2. +\
            np.sum(energytrantopz)/2. +\
            np.sum(energytrantipx)/2. +\
            np.sum(energytrantipy)/2. +\
            np.sum(energytrantipz)/2.
    print("||||||||||||||||",energy)
    print("||||||long||||||",energylong)
    print("||||||tox ene|||",energytrantopx)
    print("||||||tox force|",forcestrantopx)

    energymevlong = np.sum(energylong)/2.*1000/(numberofatoms-1)
    energymevtopx = np.sum(energytrantopx)/2.*1000/(numberofatoms-1)
    energymevtopy = np.sum(energytrantopx)/2.*1000/(numberofatoms-1)
    energymevtopz = np.sum(energytrantopx)/2.*1000/(numberofatoms-1)
    energymevtipx = np.sum(energytrantipx)/2.*1000/(numberofatoms-1)
    energymevtipy = np.sum(energytrantipx)/2.*1000/(numberofatoms-1)
    energymevtipz = np.sum(energytrantipx)/2.*1000/(numberofatoms-1)


    # make a check if sum of forces is zero!
    forces = forceslong + forcestrantopx + forcestrantopy + forcestrantopz \
                        + forcestrantipx + forcestrantipy + forcestrantipz

    if type(hessefile) != bool and type(hessefile1nn) != bool:
        forces = forces + hrestf
        energy = energy + hrestev




    #######################################################################################
    # SHOW RESULTS SHOW
    # welchen fehler machen wir: alleine in den ersten Nachbarn:
    # der Longitudinale teil in Al hat einen feher von ca:
    #   (1) - 0.016 wenn wir den d = 0.3 bereich waehlen (fitted with morse)
    #   (2) - 0.004 wenn wir den d = 0.2 bereich waehlen (fitted with morse)
    #   (3) - 0.0026 im repulsiven teil wenn nur dieser   fitted with morse mc1
    #   (4) - 0.00008 im attraktiven teil wenn nur dieser fitted with morse mc1
    #       - um den tox teil zu berechnen nehmen wir den (4) beitrag, von daher sollte der
    #         fehler auf den tox teil klein sein.
    #######################################################################################
    if verbose == True:  # SHOW RESULTS
        # this just gives us the neighborlist, we need this for printing
        neighborlist = crystal0.get_NNlist(0, 1,
                coord_cart = False, #crystal0.rcar, # hopefully thos do not \
                cell = crystal0.cellvec,
                coord_rrel = crystal0neverchange.rrel,
                return_NNdist = return_NNdist, return_d_NNidx = True)

        #print "numberofatoms:",numberofatoms
        np.set_printoptions(threshold=np.nan)  # print the whole array
        np.set_printoptions(linewidth=240)    # print only 6 digist after .
        np.set_printoptions(precision=3)    # print only 6 digist after .
        #np.set_printoptions(precision=1)    # print only 6 digist after .
        if os.path.isfile("forces_OUTCAR"):
            fab=np.loadtxt("forces_OUTCAR")
        else:
            fab = np.zeros((numberofatoms,3))

        printforces = np.zeros((numberofatoms,23))
        faktorshow = 100.

        # lon:
        printforces[:,0]  = faktorshow*(forceslong[:,0]                                   )
        printforces[:,1]  = faktorshow*(forceslong[:,1]                                   )
        printforces[:,2]  = faktorshow*(forceslong[:,2]                                   )

        # |tox:
        printforces[:,3]  = faktorshow*(forcestrantopx[:,0]+forcestrantopy[:,0]+forcestrantopz[:,0]                                   )
        printforces[:,4]  = faktorshow*(forcestrantopx[:,1]+forcestrantopy[:,1]+forcestrantopz[:,1]                                   )
        printforces[:,5]  = faktorshow*(forcestrantopx[:,2]+forcestrantopy[:,2]+forcestrantopz[:,2]                                   )

        #printforces[:,6]  = faktorshow*(forcestrantopy[:,0]                                   )
        #printforces[:,7]  = faktorshow*(forcestrantopy[:,1]                                   )
        #printforces[:,8]  = faktorshow*(forcestrantopy[:,2]                                   )

        # |tix:
        printforces[:,6]  = faktorshow*(forcestrantipx[:,0]+forcestrantipy[:,0]+forcestrantipz[:,0]                                   )
        printforces[:,7]  = faktorshow*(forcestrantipx[:,1]+forcestrantipy[:,1]+forcestrantipz[:,1]                                   )
        printforces[:,8]  = faktorshow*(forcestrantipx[:,2]+forcestrantipy[:,2]+forcestrantipz[:,2]                                   )

        #printforces[:,6]  = faktorshow*(forcestrantipx[:,0]                                   )
        #printforces[:,7]  = faktorshow*(forcestrantipx[:,1]                                   )
        #printforces[:,8]  = faktorshow*(forcestrantipx[:,2]                                   )

        #printforces[:,9]  = faktorshow*(forcestrantopz[:,0]                                   )
        #printforces[:,10]  = faktorshow*(forcestrantopz[:,1]                                   )
        #printforces[:,11]  = faktorshow*(forcestrantopz[:,2]                                   )


        # |SUM:
        printforces[:,9]   = faktorshow*(forces[:,0]                                   )
        printforces[:,10]  = faktorshow*(forces[:,1]                                   )
        printforces[:,11]  = faktorshow*(forces[:,2]                                   )

        #printforces[:,9]  = faktorshow*(forces[:,0]                                       )
        #printforces[:,10] = faktorshow*(forces[:,1]                                       )
        #printforces[:,11] = faktorshow*(forces[:,2]                                       )

        #printforces[:,9]  = faktorshow*(forcestrantip[:,0]                                   )
        #printforces[:,10] = faktorshow*(forcestrantip[:,1]                                   )
        #printforces[:,11] = faktorshow*(forcestrantip[:,2]                                   )

        # |NUMAT: VASP:
        printforces[:,12] = np.arange(numberofatoms)
        #VASP forces
        printforces[:,13] = faktorshow*(fab[:,0]                                          )
        printforces[:,14] = faktorshow*(fab[:,1]                                          )
        printforces[:,15] = faktorshow*(fab[:,2]                                          )


        #printforces[:,13] = faktorshow*(hrestf[:,0]                                          )
        #printforces[:,14] = faktorshow*(hrestf[:,1]                                          )
        #printforces[:,15] = faktorshow*(hrestf[:,2]                                          )

        # |SHELL:
        printforces[:,16] = neighborlist
        # |DIFF:
        printforces[:,17] = faktorshow*(forces[:,0] - fab[:,0]                            )
        printforces[:,18] = faktorshow*(forces[:,1] - fab[:,1]                            )
        printforces[:,19] = faktorshow*(forces[:,2] - fab[:,2]                            )

        # |pos:
        printforces[:,20] = crystal0neverchange.rcar[:,0]
        printforces[:,21] = crystal0neverchange.rcar[:,1]
        printforces[:,22] = crystal0neverchange.rcar[:,2]


        if fab.shape[0] != numberofatoms:
            printforces[:] = np.nan
            print("wrong OUTCAR (nuber of atoms not correct")
        print("###############",fab.shape)
        #VASP forces
        #printforces[:,20] = faktorshow*(fab[:,0]                                          )
        #printforces[:,21] = faktorshow*(fab[:,1]                                          )
        #printforces[:,22] = faktorshow*(fab[:,2]                                          )

        #hesse.py al -p mc1 -pp 0.290265 1.292275 2.920351 -0.209085 -0.152073 -ene -i HesseMatrix_4.13 -i1nn HesseMatrix_4.13_1NN
        #hesse.py al -p mc1 -pp 0.290265 1.292275 2.920351 -0.209085 -0.152073 -ene

        # 1 # hesse.py al -p mc1 -pp 0.290264 1.292274 2.920351 -0.20908 -0.15207 -ene
        # 2 # hesse.py al -p mc1 -pp 0.207988 1.551485 2.920351 2.804549 0.033907 -ene
        # 3 # hesse.py al -p mc1 -pp 0.221583 1.520757 2.920351 8.541780 0.196864 -ene -p2 mc1 0.447778 1.051045 2.920351 0.684528 -0.24431   # fitted to 6x6x6kp
        # 4 # hesse.py al -p mc1 -pp 0.221583 1.520757 2.920351 8.541780 0.196864 -ene -i HesseMatrix_4.13 -i1nn HesseMatrix_4.13_1NN -p2 mc1 0.447778 1.051045 2.920351 0.684528 -0.24431

        # dies ist nuer tar test fuer dq=0.3
        # 1 # DIFFxdir: 0.020247
        # 2 # DIFFxdir: -0.007237
        # 3 # DIFFxdir: 0.004227
        # 1 # DIFFquer: 0.017539
        # 2 # DIFFquer: 0.000673
        # 3 # DIFFquer: 0.008501
        emp=" "*19+"|"
        emps=" "*20+"|"
        empss=" "*18+"|"
        #a  = " lon:"+emp+"tox:"+emps+"toy:"+empss+"toz:"+empss+"NUMAT: "+"HES:"+emp+"NEI:   "+"DIF:"+empss+"pos:"
        a  = " lon:"+emp+"tox:"+emp+"tix:"+emp+"SUM:"+emp+"NUMAT: "+"VASP:"+emp+"SHELL: "+"|DIFF:"+empss+"pos:"
        #a  = "forcestotal:                     forceslong:                   forcestran_out:               forcestran_in:                   VASP:                          DIFF:"
        #a = "forcestotal:           forceslong:           forcestan:          VASP:                 DIFF:                HesseRest:"
        b =  "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        print(b)   # ------------------------------------
        print(a)   # lon:                   |tox: ...
        print(b)   # ------------------------------------
        def printfield(printforces,start,stop,color):
            A = printforces[printforces[:,16].argsort()][start:stop]
            B = A[np.lexsort((A[:, 22], A[:, 21], A[:,20], A[:,16]))]
            #utils.printarray(B,color=color,decimals=2)
            # decimals defines
            #utils.printarray(B,color=color,decimals=[2,2,2,2,2,2,2,2,2,1,1,1,0,1,1,1,0,1,1,1,2,2,2,2],prepend=['','','','','','','','','','|','','','|','','','','|','|','','','|'])
            posaccuracy = 3
            utils.printarray(B,color=color,decimals=[2,2,2,2,2,2,2,2,2,1,1,1,0,   2,2,2,0,   1,1,1,posaccuracy,posaccuracy,posaccuracy,2],prepend=['','','','','','','','','','|','','','|','','','','|','|','','','|'])

        printfield(printforces,1,32,color=utils.printred)   # 32 defines how many lines to write out !
        #printfield(printforces,13,19,utils.printblue)
        #printfield(printforces,19,43,utils.printgreen)
        printfield(printforces,0,1,utils.printred)

        checkexit=abs((forces-fab)).max()
        print("|||||||||",checkexit," == forces max diff, save_vecs_to_file_for_DOS 2:",save_vecs_to_file_for_DOS)
        if save_vecs_to_file_for_DOS == True:
            tvecdiffall = np.concatenate((tvecdiffall,np.array([[checkexit]])),axis=0)
            print("tvecdiffall2:",tvecdiffall)

        #if checkexit > 5e-4:
        #    print "checkexit:",checkexit
        #    sys.exit("ERROR checkexit too large")


        def eneprint():
            rou=4
            #print "energy / energy (mev):",round(energy,rou),"\t",round(energy*1000/(numberofatoms-1),rou)
            print(2*"---------------------------------------------------------------------------")
            #print "energylong           :",round(np.sum(energylong)/2.,rou)    ,"\t",round(np.sum(energylong)/2.*1000/(numberofatoms-1)           ,rou)
            #print "energytrantopx       :",round(np.sum(energytrantopx)/2.,rou),"\t",round(np.sum(energytrantopx)/2.*1000/(numberofatoms-1)       ,rou)
            #print "energytrantopy       :",round(np.sum(energytrantopy)/2.,rou),"\t",round(np.sum(energytrantopy)/2.*1000/(numberofatoms-1)       ,rou)
            #print "energytrantopz       :",round(np.sum(energytrantopz)/2.,rou),"\t",round(np.sum(energytrantopz)/2.*1000/(numberofatoms-1)       ,rou)
            #print "hrestev/hrestmev     :",round(hrestev,rou),"\t",round(hrestmev,rou)
            print(utils.printred("ENERGY (mev/atom):"+\
                    str(round(energy*1000/(numberofatoms-1),rou))),\
                    \
                    utils.printgreen(
                    "energylong :"+\
                    str(round(np.sum(energylong)/2.*1000/(numberofatoms-1)       ,rou))\
                    ),\
                    utils.printblue(
                    "top{x,y,z}:"+\
                    str(round(np.sum(energytrantopx)/2.*1000/(numberofatoms-1)       ,rou))+" "+\
                    str(round(np.sum(energytrantopy)/2.*1000/(numberofatoms-1)       ,rou))+" "+\
                    str(round(np.sum(energytrantopz)/2.*1000/(numberofatoms-1)       ,rou))+" ",\
                    ),\
                    \
                    utils.printyellow(
                    "tip{x,y,z}:"+\
                    str(round(np.sum(energytrantipx)/2.*1000/(numberofatoms-1)       ,rou))+" "+\
                    str(round(np.sum(energytrantipy)/2.*1000/(numberofatoms-1)       ,rou))+" "+\
                    str(round(np.sum(energytrantipz)/2.*1000/(numberofatoms-1)       ,rou))+" ",\
                    ))#,\
                    #\
                    #utils.printyellow(
                    #"tip{x,y,z}:",\
                    #round(np.sum(energytrantipx)/2.*1000/(numberofatoms-1)       ,rou),\
                    #round(np.sum(energytrantipy)/2.*1000/(numberofatoms-1)       ,rou),\
                    #round(np.sum(energytrantipz)/2.*1000/(numberofatoms-1)       ,rou)\
                    #)
            print(2*"---------------------------------------------------------------------------")
        eneprint()
        print("--kk")
        with printoptions(precision=2, suppress=True):
            potprint(u1nn_pot2param)

    #######################################################################################
    # check forces sum
    #######################################################################################
    def fehlermelungforces():
        print("sum:",np.sum(forces))
        print("sum forceslong:",np.sum(forceslong))
        print("sum forcestrantopx:",np.sum(forcestrantopx))
        print("sum forcestrantopy:",np.sum(forcestrantopy))
        print("sum forcestrantopz:",np.sum(forcestrantopz))
        print("sum forcestrantipx:",np.sum(forcestrantipx))
        print("sum forcestrantipy:",np.sum(forcestrantipy))
        print("topx:",np.sum(forcestrantopxcheck))
        print("topy:",np.sum(forcestrantopycheck))
        print("topz:",np.sum(forcestrantopzcheck))
        print("tipx:",np.sum(forcestrantipxcheck))
        print("tipy:",np.sum(forcestrantipycheck))
        print("tipz:",np.sum(forcestrantipzcheck))
        for idxftix,ftix in enumerate(forcestrantipxcheck):
            if np.sum(ftix) != 0.0:
                print("idxftix:",idxftix,np.sum(ftix))
                print(ftix)
        sys.exit("sum of forces is not 0")
    if abs(np.sum(forceslong)) >= 0.000001:
        print("fehler,forceslong")
        fehlermelungforces()
    if abs(np.sum(forces)) >= 0.000001:
        print("fehler,forces")
        fehlermelungforces()


    if abs(np.sum(forcestrantopx)) >= 0.000001:
        print("fehler,forcestrantopx")
        fehlermelungforces()
    if abs(np.sum(forcestrantopy)+np.sum(forcestrantopz)) >= 0.000001:
        print("fehler,forcestrantop{y,z}",abs(np.sum(forcestrantopx)-np.sum(forcestrantopy)))
        fehlermelungforces()

    if abs(np.sum(forcestrantipx)) >= 0.000001:
        print("fehler,forcestrantipx")
        fehlermelungforces()
    if abs(np.sum(forcestrantipy)) >= 0.000001:
        print("fehler,forcestrantipy")
        fehlermelungforces()


    if energy < -1.0:
        print("energy:",energy)
        sys.exit("Potential energy is negative! (probably there is a jump / change of poitions from equilibrium structure)")
    np.savetxt("forces",forces)
    np.savetxt("energy",np.array([energy]))
    #print "hf:",hessefile,hessefile1nn,hrestev
    #ENDPROG
    if save_vecs_to_file_for_DOS == True:
        np.savetxt("tvecoutall.dat",tvecoutall,fmt="%.6f")
        np.savetxt("tvecinsall.dat",tvecinsall,fmt="%.6f")
        np.savetxt("tveclonall.dat",tveclonall,fmt="%.6f")  # == longvec
        np.savetxt("tvecdiffall.dat",tvecdiffall,fmt="%.6f")  # == diff between forces_OUTCAR and forces
        print("tvecdiffall3:",tvecdiffall,"save_vecs_to_file_for_DOS 3:",save_vecs_to_file_for_DOS)

        #np.savetxt("tveclonallf.dat",tveclonallf,fmt="%.6f")
        #np.savetxt("tveclonalle.dat",tveclonalle,fmt="%.6f")
        if os.path.isfile("energy_long") != True:
            open("energy_long", 'a').close()
        utils.write_inputdata("energy_long","step "+str(round(energy*1000/(numberofatoms-1),4)))
        if os.path.isfile("energy_tox") != True:
            open("energy_tox", 'a').close()
        utils.write_inputdata("energy_tox","step "+str(round(np.sum(energytrantopx)/2.*1000/(numberofatoms-1),4)))

    #return energy,forces   # energy is the energy per cell in eV
    # END END END
    return energy,forces,energymevlong,energymevtopx   # energy is the energy per cell in eV

def get_energy_forces_hesse_rest(hessefile = False, hessefile1nn = False):
    if type(hessefile) == bool:
        sys.exit("please specify hessefile")
    if type(hessefile1nn) == bool:
        sys.exit("please specify hessefile1nn")
    if os.path.isfile('EqCoords_direct') != True:
        sys.exit('EqCoords_direct'+" does not exist")
    if os.path.isfile(hessefile) != True:
        sys.exit(hessefile+" does not exist")


    if os.path.isfile(hessefile1nn) != True:
        sys.exit(hessefile1nn+" does not exist")

    hfullmev, hfulle, hfullf = get_energy_forces(
            pot='h',
            potparam=False,
            hessefile = hessefile,
            coordfile_cart = "cartesian_coords",
            printresult = False
            )
    h1NNmev, h1NNe, h1NNf = get_energy_forces(
            pot='h',
            potparam=False,
            hessefile = hessefile1nn,
            coordfile_cart = "cartesian_coords",
            printresult = False
            )
    hrestmev = hfullmev - h1NNmev
    hreste = hfulle - h1NNe
    hrestf = hfullf - h1NNf
    return hrestmev, hreste, hrestf

def getforcepot(longvec,pot,potparam,pot2 = False):
    usepot=['nan','nan','nan','nan','nan','nan','nan']
    if pot == 'm': usepot[0] = "Morse"
    if pot == 'l': usepot[0] = "LJ"
    if pot == 'i': usepot[0] = "inversepot"
    if pot == 'mc1': usepot[0] = "mc1"
    for idx,p in enumerate(potparam):
        usepot[idx+1] = float(p)

    longvecnorm=np.linalg.norm(longvec)
    getfrompot2 = False
    if type(pot2) != bool:
        if usepot[0] == "Morse" or usepot[0] == "mc1":
            NNdist = float(pot2[3])
            #print "NNdist:",NNdist,longvecnorm
            if longvecnorm > NNdist:
                getfrompot2 = True
        else:
            sys.exit("pot2 only m or mc1")

    if usepot[0] == "inversepot":
        elong = inversepot(longvecnorm,usepot[1],usepot[2],usepot[3])
        flongnorm = inversepot_derivative(longvecnorm,usepot[1],usepot[2],usepot[3])

    if usepot[0] == "Morse":
        if getfrompot2 == False:
            elong = Morse(longvecnorm,usepot[1],usepot[2],usepot[3])
            flongnorm = Morse_derivative(longvecnorm,usepot[1],usepot[2],usepot[3])
        else:
            if pot2[0] != 'm':
                sys.exit("not m")
            elong =                Morse(longvecnorm,float(pot2[1]),float(pot2[2]),float(pot2[3]))
            flongnorm = Morse_derivative(longvecnorm,float(pot2[1]),float(pot2[2]),float(pot2[3]))


    if usepot[0] == "mc1" or usepot[0] == 'mc1':
        #print "usepot:",usepot
        if getfrompot2 == False:
            elong = mc1(longvecnorm,usepot[1],usepot[2],usepot[3],usepot[4],usepot[5])
            flongnorm = mc1_derivative(longvecnorm,usepot[1],usepot[2],usepot[3],usepot[4],usepot[5])
        else:
            if pot2[0] != 'mc1':
                print("pot2:",pot2)
                sys.exit("not mc1")
            elong =                mc1(longvecnorm,float(pot2[1]),float(pot2[2]),float(pot2[3]),float(pot2[4]),float(pot2[5]))
            flongnorm = mc1_derivative(longvecnorm,float(pot2[1]),float(pot2[2]),float(pot2[3]),float(pot2[4]),float(pot2[5]))


    if usepot[0] == "LJ":
        elong = LJ(longvecnorm,usepot[1],usepot[2])
        flongnorm = LJ_derivative(longvecnorm,usepot[1],usepot[2])

    # flong ist die longitudinale Kraft
    flong = longvec/longvecnorm*flongnorm   # der volle vektor == voller longitudinale vektor
    if flong[0] == -0.0:
        flong[0] = 0.0
    if flong[1] == -0.0:
        flong[1] = 0.0
    if flong[2] == -0.0:
        flong[2] = 0.0
    #print np.around(flong[0],decimals=7),np.around(flong[1],decimals=7),np.around(flong[2],decimals=7)
    return flong

if __name__ == '__main__':
    p = help()  # this gives the possibility to change some __init__ settings
    args = p.parse_args()

    if args.potparam:
        if not args.pot:
            sys.exit("potparam but no pot defined!")
        #print "ap:",args.potparam,len(args.potparam)
        #if args.pot != "m" or args.pot != "l" or args.pot != "h":   # this is taken care of in the help by "choices"
        #    sys.exit("pot has to be m or l or h")
        #if len(args.potparam) == 3:  # this is done in skript
        #    potparam = [ 0, 0, 0]
        #    potparam[0] = float(args.potparam[0])
        #    potparam[1] = float(args.potparam[1])
        #    potparam[2] = float(args.potparam[2])
        #    args.potparam = potparam

    if args.getforcepot:
        #print "potparam:",args.potparam,"||",args.getforcepot
        longvec = np.array([np.nan, np.nan, np.nan])
        if len(args.getforcepot) != 3:
            sys.exit("please give the vector as an argument ")
        for idx,number in enumerate(args.getforcepot):
            #print idx,number
            longvec[idx] = float(number)
        #print "--"
        getforcepot(longvec=longvec,pot=args.pot,potparam=args.potparam,pot2=args.pot2)
        sys.exit()

    if args.ene:
        hessefile = False
        if args.pot == "h" or args.hrest1 == True:
            hessefile = args.inputfile
        if os.path.isfile("cartesian_coords") != True:
            sys.exit("please provide cartesian_coords file")
        if args.verbose:
            print("verbose      :",args.verbose)
            print("pot          :",args.pot)
            print("potparam     :",args.potparam)
            print("hessefile    :",hessefile)
            print("trans from f :",args.tm)
            print("args.hrest1  :",args.hrest1)
        get_energy_forces(
                pot=args.pot,
                pot2=args.pot2,
                potparam=args.potparam,
                hessefile = args.inputfile,
                hessefile1nn = args.inputfile1nn,
                coordfile_cart = "cartesian_coords",
                add_transverse = args.tm,
                hrest1 = args.hrest1,
                verbose=args.verbose
                )
        sys.exit()

    if args.eneext:
        if os.path.isfile("cartesian_coords") != True:
            sys.exit("please provide cartesian_coords file")
        get_energy_forces_pairpot(usepot=['extern'],verbose=args.verbose)
        sys.exit()

    if args.eneexto:
        if os.path.isfile("cartesian_coords") != True:
            sys.exit("please provide cartesian_coords file")
        get_energy_forces_pairpot(usepot=['externorig'],verbose=args.verbose)
        sys.exit()

    if args.hesse1nn:
        print("args.hesse1nn[0]:",args.hesse1nn[0])
        print("args.hesse1nn[1]:",args.hesse1nn[1])
        ka = hesse_to_hesse_1NN(
            hfilename = args.hesse1nn[0],
            cellfile = args.hesse1nn[1],
            coordfile0_rrel = "EqCoords_direct")
        print(ka)
        sys.exit()

    if args.hrest1:
        ka = get_energy_forces_hesse_rest(hessefile = args.hrest1[0], hessefile1nn = args.hrest1[1])
        print(ka)
        sys.exit()

    if args.plotm:
        if len(args.plotm) != 3:
            sys.exit("len morseh is not correct")
        potparam = [ 0, 0, 0]
        potparam[0] = float(args.plotm[0])
        potparam[1] = float(args.plotm[1])
        potparam[2] = float(args.plotm[2])
        plot_morse(potparam)
        sys.exit()

    if args.fitpotf:
        print(args.fitpotf)
        if args.pot == False:
            sys.exit("please specify a pot with -p or --pot")
        NN = False
        if len(args.fitpotf) == 2:
            NN = args.fitpotf[1]
        get_fit_forces_to_pot(filename=args.fitpotf[0], foldername = os.getcwd(), NN = NN, pot=args.pot)
        sys.exit()

    if args.polyfit:
        #if len(args.polyfit) != 2:
        #    sys.exit("first argument has to be the filename, second the order")
        #if os.path.isfile(args.polyfit[0]) != True:
        #    print "filename:",args.polyfit[0]," does not exist"
        #    sys.exit("filename does not exist")
        fitdelta = polyfit(*args.polyfit)
        sys.exit()

    if args.morseh:
        if len(args.morseh) != 3:
            sys.exti("len morseh is not correct")
        potparam = [ 0, 0, 0]
        potparam[0] = float(args.morseh[0])
        potparam[1] = float(args.morseh[1])
        potparam[2] = float(args.morseh[2])
        get_hessematrix_lj_morse(pot='m',potparam=potparam)
        sys.exit()

    if args.ljh:
        if len(args.morseh) != 2:
            sys.exti("len ljh is not correct")
        potparam = [ 0, 0]
        potparam[0] = float(args.morseh[0])
        potparam[1] = float(args.morseh[1])
        get_hessematrix_lj_morse(pot='l',potparam=potparam)
        sys.exit()

    #if args.fitfqh:
    #    print args.fitfqh
    #    sys.exit()
    if args.fitmf:
        fit_fcc_alat_mean_freqs()
        sys.exit()

    print('args;',args)
    if args.ffd:
        import shutil
        import glob
        sc = "2x2x2"
        elements = [ "Al", "Pb", "Cu", "Rh", "Pd", "Au", "Ir", "Pt", "Au" ]
        elements = [ "Cu" ]
        xq = "xdir"
        xq = "quer"
        xq = "xdir"
        hier = os.getcwd()
        for el in elements:
            searchfor = os.getcwd()+"/"+el+"/2x2x2sc_"+xq+"_[1-4]*kp/"
            searchfor = os.getcwd()+"/"+el+"/3x3x3sc_"+xq+"_8x8x8kp/"
            searchfor = os.getcwd()+"/"+el+"/2x2x2sc_"+xq+"_3x3x3kp_vasp4/"
            searchfor = "/Users/glensk/Dropbox/Understand_distributions/Cu/2x2x2sc_quer_3x3x3kp_vasp4"
            searchfor = "/Users/glensk/Dropbox/Understand_distributions/Cu/2x2x2sc_xdir_3x3x3kp_vasp4"
            dirquercheck = glob.glob(searchfor)
            dirquer = False
            if len(dirquercheck) == 0:
                print(el," no folder found; searchfor:",searchfor)
                continue
            if len(dirquercheck) == 1:
                dirquer = dirquercheck[0]
            else:
                print(dirquercheck)
                sys.exit("more then one folder found")
            print(el,dirquer)
            #if False:
            if True:
                os.chdir(dirquer)




                shutil.copyfile("/Users/glensk/Dropbox/Understand_distributions/EqCoords_direct",os.getcwd()+"/EqCoords_direct")
                n = forcesneighbors()
                allfolder=glob.glob("*Ang_*")
                def is_substr(find, data):
                    if len(data) < 1 and len(find) < 1:
                        return False
                    for i in range(len(data)):
                        if find not in data[i]:
                            return False
                    return True

                def long_substr(data):
                    substr = ''
                    if len(data) > 1 and len(data[0]) > 0:
                        for i in range(len(data[0])):
                            for j in range(len(data[0])-i+1):
                                if j > len(substr) and is_substr(data[0][i:i+j], data):
                                    substr = data[0][i:i+j]
                    return substr
                n.u = 4.1
                print("--1--")
                print(allfolder)
                print("--2--")
                print(long_substr(allfolder))
                print("--3--")
                print(long_substr(allfolder).split("Ang"))
                print("--4--")
                n.a = float(long_substr(allfolder).split("Ang")[0])
                n.verbose = False
                print("## n.loadforces:")
                n.loadforces()
                print("## n.getforces:")
                n.getforces()


                os.chdir(hier)
        sys.exit()

    if args.fqh_interpolate:
        print(args.inputfile)
        fqh_interpolate(fqh = args.inputfile)
        sys.exit()


    if True:
        h = hesseclass(args)
        print(h.H)
        print()
        print(h.M)
        if args.verbose:
            print("verbose      :",args.verbose)
            print("pot          :",args.pot)
            print("potparam     :",args.potparam)
            print("hessefile    :",args.inputfile)
            print("trans from f :",args.tm)
            print("args.hrest1  :",args.hrest1)
        get_energy_forces(
                pot=args.pot,
                pot2=args.pot2,
                potparam=args.potparam,
                hessefile = args.inputfile,
                hessefile1nn = args.inputfile1nn,
                coordfile_cart = "cartesian_coords",
                add_transverse = args.tm,
                hrest1 = args.hrest1,
                verbose=args.verbose
                )
        sys.exit()

    if args.eneext:
        if os.path.isfile("cartesian_coords") != True:
            sys.exit("please provide cartesian_coords file")
        get_energy_forces_pairpot(usepot=['extern'],verbose=args.verbose)
        sys.exit()

    if args.eneexto:
        if os.path.isfile("cartesian_coords") != True:
            sys.exit("please provide cartesian_coords file")
        get_energy_forces_pairpot(usepot=['externorig'],verbose=args.verbose)
        sys.exit()



##hesse = hesseclass(args.elements,h_filename=args.inputfile,verbose=args.verbose)

    ##hesse = hesseclass(args.elements)
    ##crystal.read_Hessematrix()
    ##pos1 = crystal.read_rcar_positions()
    ##pos1r = crystal.rrel
    ##pos2 = crystal.read_rrel_positions()
    ##pos2r = crystal.rrel
    #### map back
    ##diff = pos1 - pos2
    ##diffr = pos1r - pos2r
    ##for i in np.arange(3):
    ##    a = diff[:,i]
    ##    b = diffr[:,i]
    ##    maska = a > crystal.cellvec[0,0]/2.
    ##    maskb = b > 0.5
    ##    a[maska] += -crystal.cellvec[0,0]
    ##    b[maskb] += -1.

    ##diffrr   = np.dot( diffr, crystal.cellvec )
    ##h = crystal.read_Hessematrix()
    ##f = crystal.qh_forces(diff, h)
    ##e = crystal.qh_energy_atom(diff,h)
    ##np.savetxt("energy",np.array([e]),fmt='%.12f')
    ##np.savetxt("forces",np.array(f),fmt='%.12f')
