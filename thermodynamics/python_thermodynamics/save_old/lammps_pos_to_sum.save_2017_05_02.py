#!/usr/bin/env python

# lammps_pos_to_sum.py bcc 6 3.253 4 4 4 50   # dominiques job

# /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2016.09_phonon_lifetimes_3_nach_elternzeit/tutild_2016_12_02_10sc_300K_4.04_2_hoch_5_steps

# 2*10^5 schritte --> ein qpunkt hat 19MB
# 1*10^6 schritte --> ein qpunkt hat 200MB  --> ein ast geht locker (20 qpunkte sind 4GB)


###########################################
# Timing of this skript a single qpoint for maximal length
###########################################
# started using:
# python -m cProfile -s cumtime ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py -N 10 -a 4.14 -dt 40 -q 8 8 20 -ps -calclast
###########################################

##########################################################################################
# run argparse first to make it quick
##########################################################################################
#def multi_run_wrapper(args):
#    return add(*args)
#
#def currentfunc_helper(args):
#    return get_lifetime_for_particular_mdsteps(*args)
#
#def add(x,y):
#    return x+y
import argparse
import textwrap
#from argparse import ArgumentDefaultsHelpFormatter

def lorenz(w,w0,a,b):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe
        b breite'''
    return a/((w-w0)**2.+b**2.)

def lorenzbesser(w,w0,a,b):
    ''' w: frequenzabhaengigkeit
        w0: die frequenz
        a hoehe
        b breite'''
    return a/((w**2-w0**2)**2.+w**2*b**2.)

def number_of_atoms_to_struct(atoms):
    fcc = []
    bcc = []
    for i in np.arange(100)+1:
        #print i,"fcc:",4*(i**3),"\tbcc:",2*(i**3)
        fcc.append(4*(i**3))
        bcc.append(2*(i**3))
    #print fcc
    #print bcc
    #print "4 in bcc:",4 in bcc
    #print "4 in fcc:",4 in fcc
    for i,k in enumerate(np.arange(100)+1):
        #print "fcc["+str(i)+"] in bcc:",fcc[i] in bcc
        pass
    if atoms in fcc:
        sc = N = int(round(float((atoms/float(4))**(1/3.)),0))
        return "fcc",sc
    if atoms in bcc:
        sc = N = int(round(float((atoms/float(2))**(1/3.)),0))
        return "bcc",sc

def get_inputvariables_from_calculation(positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in", inN=False,inalat=False,indt=False,inqvec=False,verbose=False,args=False):
    ''' get the lammps input variables
        gets:
            - atoms (number of atoms e.g. 32)
            - sc = N (supercellsize e.g. 2)
            - alat (in Angstrom e.g. 4.14)
            - supercelllength (in Angstrom e.g. 8.28)
            - dt (timestep of dumping of positions)

    '''
    positionsfile = glob.glob(positionsfilename)
    if len(positionsfile) != 1:
        #print "ERROR: positionsfile:",positionsfile
        return False

    infile = glob.glob(infilefilename)  # in_file_dynamics.in
    #print "infile:",infile
    if len(infile) != 1:
        print "ERROR: infile:",infile
        return False

    print "JOB           : LAMMPS"
    atoms1 = os.popen('wc -l '+positionsfile[0] + "| awk '{print $1}'").read()
    atoms = int(atoms1)-8
    atoms2 = int(os.popen('head -n 2 '+positionsfile[0] + "| tail -1 | awk '{print $1}'").read().rstrip())
    #print "atoms2:",atoms2
    if atoms != atoms2:
        print "positionsfile:",positionsfile[0]
        print "ERROR: atoms: (wc -l positionsfile -8)",atoms,"atoms2 (head -n 2 positionsfile):",atoms2,"from",positionsfile[0]
        return False


    supercelllength = float(os.popen('head -n 4 '+positionsfile[0] + "| tail -1 | awk '{print $2}'").read().rstrip())
    if verbose:
        print "supercelllength:",supercelllength
    #sc = N = int(round(float((atoms/float(args.usestruct))**(1/3.)),0))
    args.structure, args.supercell = number_of_atoms_to_struct(atoms)
    sc = N = args.supercell
    alat = a = supercelllength/N
    if verbose:
        print "alat:",alat

    # check if structure and atoms make sense
    if type(args.usestruct) != bool:
        atoms3 = args.usestruct * sc**3.
        if atoms != atoms3:
            print "args:",args
            print "usestruct:",args.usestruct, "4=fcc, 2=bcc use option: -s {bcc,fcc}"
            print "ERROR: atoms:",atoms,"atoms3:",atoms3, "maybe you have bcc instad of fcc or viec versa? define this!"
            return False

    stepslammps = int(os.popen('grep "^run " '+infile[0] + "| awk '{print $2}'").read().rstrip())
    timestepfs = int(os.popen('grep "^timestep " '+infile[0] + "| awk '{print $2*1000}'").read().rstrip())
    dt = dump = os.popen('grep "^dump dump1" '+infile[0] + "| awk '{print $5}'").read().rstrip()
    try:
        dt = dump = int(os.popen('grep "^dump dump1" '+infile[0] + "| awk '{print $5}'").read().rstrip())*timestepfs
    except ValueError:
        dt = False

    if verbose:
        print "timestepfs:",timestepfs
        print "atoms:",atoms,type(atoms)
        print "sc:",sc
        print "stepslammps:",stepslammps
        print "dump:",dump
    stepswritten = stepslammps/dt+1
    if verbose:
        print "stepswritten:",stepswritten
    lineswritten = stepswritten*atoms
    if verbose:
        print "lineswritten:",lineswritten

    #########################################
    # check against invariables
    #########################################
    if type(args.folder) == bool:  # if False --> check
        if type(inN) != bool:
            if inN != N:
                print "args.folder:",args.folder
                print "inN:",inN
                print "N  :",N
                sys.exit("ERROR: inN not equal N!")
    if type(inalat) != bool:
        if inalat != alat:
            print "inalat:",inalat
            print "alat  :",alat
            print "os.getcwd():",os.getcwd()
            if args.folder == False:
                sys.exit("ERROR: inalat not equal alat!")
            else:
                print "now changing from inalat:",inalat,"to alat  :",alat
                args.alat = alat
                inalat = alat
    if type(indt) != bool:
        if indt != dt:
            print "indt:",indt
            print "dt  :",dt
            sys.exit("ERROR: indt not equal dt!")

    #########################################
    # check against invariables
    #########################################
    qvec = inqvec
    if type(inqvec) == bool:
        qvec = [ 'all', 'all', 'all' ]


    if args.supercell == False: args.supercell = N
    elif args.supercell != N: sys.exit("args.supercell, N: "+str(args.supercell)+" "+str(N))


    if args.dt == False: args.dt = dt
    elif args.dt != dt: sys.exit("args.dt, dt: "+str(args.dt)+" "+str(dt))

    if args.alat == False: args.alat = alat
    elif args.alat != alat: sys.exit("args.alat, alat: "+str(args.alat)+" "+str(alat))

    if args.qvec == False: args.qvec = qvec
    return

def get_inputvariables_from_vasp_job():
    vaspinput = False
    if os.path.isfile("INCAR") and os.path.isfile("POSCAR"): # and os.path.isfile("POSITIONs"):
        vaspinput = True
    if vaspinput:
        print "JOB           : VASP"
    else:
        return
    try:
        atoms = int(float(os.popen("head -n 7 POSCAR | tail -1 | awk '{print $1}'").read().rstrip()))
    except ValueError:
        atoms = int(float(os.popen("head -n 6 POSCAR | tail -1 | awk '{print $1}'").read().rstrip()))
    dt = int(float(os.popen("grep POTIM INCAR | awk '{print $3}'").read().rstrip()))
    scale = float(os.popen("head -n 2 POSCAR | tail -1 | awk '{print $1}'").read().rstrip())
    sc11 = float(os.popen("head -n 3 POSCAR | tail -1 | awk '{print $1}'").read().rstrip())
    #print "atoms:",atoms
    #print "scale:",scale
    #print "sc11:",sc11
    args.dt = dt
    args.alat = scale*sc11


    args.structure, args.supercell = number_of_atoms_to_struct(atoms)
    args.alat = scale*sc11/args.supercell
    return

def help(p = None):
    string = '''
    - use this skript without options if in a lammps job (in_file_dynamics.in and positions*.dat available!)
      otherwise define required options manually!

    - required input is: (if no means of parsing from in_file_dynamics.in and positions*)
                    - N = sc = n
                    - dt
                    - a = alat
                    - q = qvec

    - examples:
        lammps_pos_to_sum.py (use without options if in a lammps job)
        lammps_pos_to_sum.py -N 10 -a 4.13 -dt 40 -q l 0 0
        lammps_pos_to_sum.py -N 3 -a 4.14 -dt 1 -q lt 0 0 -lt
        lammps_pos_to_sum.py -N 3 -a 4.14 -dt 1 -q l l l -ps -fftpy
        lammps_pos_to_sum.py -sc 10 -dt 40 -a 4.14 -q all 0 0

    - how to convert trj_lammps.out -> trj_lammpsnew.out:
             grep "^1 " trj_lammps.out | awk \'{print $2,$3,$4}\' > trj_lammpsnew.out

    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)

    p.add_argument('-s',   '--structure',choices=[ 'fcc', 'bcc' ],
       help='which structure was calculated? Currently only fcc and bcc are supported. Simple cubic should work too.', type=str, default='fcc')
    p.add_argument('-f',   '--folder',
       help='which folder to evaluate; expansion symbols are allowed but enclose in quotes e.g. "run_*" or "run_[1-5]"', type=str, default=False)
    p.add_argument('-u',   '--usestruct',help=argparse.SUPPRESS, type=int, default=4) # usestruct is 4 for fcc (4 atoms in cubic cell) and 2 for bcc (2 atoms in cubic cell)
    p.add_argument('-n', '-sc','--supercell' , required=False,
       help='supercell: 2x2x2 -> 2', type=int, default=False)
    p.add_argument('--N', required=False,
       help=argparse.SUPPRESS, type=int, default=False)
    p.add_argument('-dt',   '--dt' , required=False,
       help='timestep of your calculation in femto seconds.', type=float, default=False)
    p.add_argument('-a',   '--alat' , required=False,
       help='alat of your qubic cell in angstromr; e.g. 4.13', type=float, default=False)
    #p.add_argument('-ka', '--ka',action='append',nargs=3,metavar=('q1','q2','q3'),required=False,type=str,help='help:')
    p.add_argument('-q',   '--qvec' , required=False, action='append',nargs=3, type=str,
            help=textwrap.dedent('''\
       q1 q2 q3;
       Use 3 0 0      (three integers) to specivy a particular q-vector;
       use l 0 0      speify all in Longitudinal [ i 0 0 ] direction;
       use t 0 0      speify all in Transversal  [ 2N i 0 ] direction;
       use lt 0 0     speify all in Long and Trans   [ 2N i 0 ] direction;

       use l l l      speify all in Longitudinal [ i i i ] direction;
       use t t t      speify all in Transversal  [ ? ? ? ] direction;

       use l l 0      speify all in Longitudinal [ i i 0 ] direction;
       use t1 t1 0    speify all in Transversal1  [ i i 2N ] direction;
       use t2 t2 0    speify all in Transversal1  [ 2N+i 2N-i 0 ] direction;

       use all 0 0    if one of qvec is "all" then all branches are evaluated in Longitudinal and Transversal direction;
       use all l l    if one of qvec is "all" then all branches are evaluated in Longitudinal and Transversal direction;
       '''))

    p.add_argument('-cl','--create_lammps_inputfile', action='store_true', default=False,
        help='create in_file_dynamics.in including defined qpoints')
    p.add_argument('-space_fft','--space_fft', action='store_true', default=False,
        help='get space_fft_x_x_x.npy from xaa_, xab_, ... folder')

    p.add_argument('-lt','--lifetimes', action='store_true', default=False,
            help='print lifetimes in the specified q-vecotor direction (and write to file); (units: THz );')
    p.add_argument('-freqs','--freqs', action='store_true', default=False,
            help='print frequencies in the specified q-vecotor direction (and write to file); (units: THz);')

    #        help='average several lifetimes/frequencies files); (units: THz );')
    p.add_argument('--lifetimesaverage','-ltsum' , required=False, type=str, nargs='*',default=False,
       help='average several lifetimes/frequencies files')

    p.add_argument('-showeq','--showeq', action='store_true', default=False,
        help='print equivalent qpoint to screen;')
    p.add_argument('-showlammps','--showlammpsinput', action='store_true', default=False,
        help='print input for lammps inputfile for this paritcular qpoint(s) to screen;')

    p.add_argument('-fftpy','--make_space_fft_from_py', action='store_true', default=False,
        help='create space_fft_x_x_x.npy file using pandas read csv;')
    p.add_argument('-fftlammpslog','--make_space_fft_from_lammpslog', action='store_true', default=False,
        help='create space_fft_x_x_x.npy lammps.log file using numpy genfromtxt;')
    p.add_argument('-fftc','--make_space_fft_from_c', action='store_true', default=False,
        help='create space_fft_x_x_x.npy file using c++ skript;')
    p.add_argument('-ps','--make_power_spectrum', action='store_true', default=False,
        help='create ps_x_x_x_MDSTEPS.dat file and evaluate error due to MD/smoothing/background;')
    p.add_argument('-calclast','--calclast', action='store_true', default=False,
        help='calculate powerspectrum only for all (maximum number of) steps and not intermediate step values to estimate the MD error;')
    p.add_argument('-write_full_ps', '-write_ps','--write_full_ps', action='store_true', default=False,
        help='Write the full powerspectrum to hd. Files can be relatively larg, 15-150 MB')
    p.add_argument('-ser', '-seriell','--seriell','-serial', action='store_true', default=False,
        help='calculate everything seriell instead of parallell(==delfault); default = False')
    p.add_argument('-write_smooth_in', '-write_si','--write_smooth_in','-si' , action='store_true', default=False,
        help='Write the smoothed powerspectrum which is used to find the linewidths.')
    p.add_argument('-i',   '--lammps_infile',
       help='name of the lammps inputfile. can also be read from file infile.infilefilename', type=str, default='in_file_dynamcs.in')


    p.add_argument('-smoothing',   '--smoothing' , required=False,
       help='wite to hd a particular smoothing ; e.g. -mdsteps 30000 -smoothing 0.001', type=float, default=False)
    p.add_argument('--mdsteps','-mdsteps' , required=False, type=int, nargs='+',default=False,
       help='calculate powerspectrum for particular mdsteplength; e.g. -mdsteps 30000 60000')
    p.add_argument('-e',   '--exit' ,  action='store_true', default=False,
       help='exit after collecting job information')
    p.add_argument('-ee',   '--exclude_existing' ,  action='store_true', default=False,
       help='exclude existing space_fft_xxx files from qpoints to be used')
    p.add_argument('-v','--verbose',
            help='verbose', action='count', default=False)
    p.add_argument('-atoms','--atoms', required=False,
       help=argparse.SUPPRESS, type=int, default=False)

    return p
p = help()  # this gives the possibility to change some __init__ settings
args = p.parse_args()

##########################################################################################
# argparse first done
##########################################################################################
def importt(args):
    ''' import statements
    # timing: 3.73059201241 for pandas
    # timing: 1.33269906044 for scipy.fftpack
    '''
    global time
    import time
    anfang= time.time()
    start = time.time()
    global fft
    global ifft
    from scipy.fftpack import fft,ifft
    end = time.time()
    if args.verbose > 1:
        print "scipy.fftpack     imported in",str((end-start)),"sec."
    start = time.time()
    global sys
    import sys
    end = time.time()
    if args.verbose > 1:
        print "sys               imported in",str((end-start)),"sec."
    start = time.time()
    global os
    import os
    end = time.time()
    if args.verbose > 1:
        print "os                imported in",str((end-start)),"sec."
    start = time.time()
    global math
    import math
    end = time.time()
    if args.verbose > 1:
        print "math              imported in",str((end-start)),"sec."
    start = time.time()
    global glob
    import glob
    end = time.time()
    if args.verbose > 1:
        print "glob              imported in",str((end-start)),"sec."
    start = time.time()
    global np
    import numpy as np
    end = time.time()
    if args.verbose > 1:
        print "numpy             imported in",str((end-start)),"sec."
    start = time.time()
    global permutations
    from itertools import permutations
    global multiprocessing
    import multiprocessing
    end = time.time()
    if args.verbose > 1:
        print "itertools         imported in",str((end-start)),"sec."
    start = time.time()
    #from pandas import set_option
    end = time.time()
    if args.verbose > 1:
        print "pandas set_option imported in",str((end-start)),"sec."
    if args.make_space_fft_from_lammpslog or args.make_space_fft_from_py:
        start = time.time()
        global read_csv
        global set_option
        from pandas import read_csv,set_option
        set_option('max_columns', 0)
        end = time.time()
        if args.verbose > 1:
            print "pandas read csv   imported in",str((end-start)),"sec."
    if args.verbose > 1:
        print "-----------------------------------------------------"
    #import pyfftw
    #import h5py
    #import timeit

    end = time.time()
    # OPTIONS
    if args.verbose > 1:
        print "all python module imported in",str((end-anfang)),"sec."
        print "-----------------------------------------------------"
    return

importt(args)
np.set_printoptions(suppress=True)   # display arrays withou 000000
np.set_printoptions(precision=6,linewidth=150)

##########################################################################################
# general functions and helperfunctions ##################################################
##########################################################################################
def printblue(text):
    return ('\x1b[4;34;40m' + text + '\x1b[0m')
def printorange(text):
    return ('\x1b[4;33;40m' + text + '\x1b[0m')
def printgreen(text):
    return ('\x1b[4;32;40m' + text + '\x1b[0m')
def printred(text):
    return ('\x1b[4;31;40m' + text + '\x1b[0m')

def parnpsavetxt(filename,data,fmt=False):
    #np.savetxt(filenameout,ps[:len(ps)/factor_psout])  # saves only half the powerspektrum (left peak)
    if type(fmt) == bool:
        p = multiprocessing.Process(target=np.savetxt, args=[filename,data])
    else:
        p = multiprocessing.Process(target=np.savetxt, args=[filename,data,fmt])
    p.start()
    return

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def get_lifetimessummary_from_several_folders(base="run_*",lifetimestring="2.0"):
    if args.lifetimesaverage == []:
        pass
        #print "[]",args.lifetimesaverage
    else:
        base = args.lifetimesaverage[0]
        #print "[[[",args.lifetimesaverage
    alatstring = str(args.alat)
    lifetimestring = str(args.dt)
    np.set_printoptions(precision=4)
    np.set_printoptions(suppress=True)
    searchforall=["l_0_0", "t_0_0", "l_l_0", "t1_t1_0", "t2_t2_0" , "l_l_l", "t_t_t" ]
    forfiles= [ "lifetimesgood_", "freqsgood_" ]
    for searchfor in searchforall:
        for forefile in forfiles:
            #print
            #print
            #print
            filename = forefile+searchfor+"_alat"+alatstring+"_dt"+lifetimestring+"_THz.dat"
            searchstring = base+filename
            #print searchfor,searchstring
            found=glob.glob(searchstring)
            if len(found) == 0:
                print "base:",base
                print "searchstring:",searchstring
                print ""
                print "--- this are the details ---"
                print "forefile:",forefile
                print "searchfor:",searchfor
                print "alatstring:",alatstring
                print "--> lifetimestring:",lifetimestring
                sys.exit("no files found")
            if len(found) >=1:
                #print "found!!!!!!!!!!!!!",len(found),found
                for idx,i in enumerate(found):
                    print "$$",idx,i
                for idx,i in enumerate(found):
                    #data=np.loadtxt(i)
                    data=np.genfromtxt(i)
                    #print "idx:",idx,i,data
                    if idx == 0:
                        #print "idx:",idx,"data.shape:",data.shape,len(found),len(data.shape)
                        if len(data.shape) == 2:
                            datacumulated = np.zeros((len(found),data.shape[0],data.shape[1]))
                        if len(data.shape) == 1:
                            datacumulated = np.zeros((len(found),data.shape[0]))
                        #print "tt;",datacumulated.shape
                    datacumulated[idx]=data
                #print "mean ++++++++++++++++++++++++++++",searchfor
                #print "shape:",datacumulated.shape,len(datacumulated.shape)
                mean = np.mean(datacumulated,axis=0)
                max = np.max(datacumulated,axis=0)
                min = np.min(datacumulated,axis=0)
                print "kk:",len(datacumulated.shape)
                if len(datacumulated.shape) > 2: # more than one line
                  mean[:,5] = 0
                  max[:,5] = 0
                  min[:,5] = 0
                else:  # only one line
                  mean[5] = 0
                  max[5] = 0
                  min[5] = 0
                #print mean
                #print "std  ++++++++++++++++++++++++++++",searchfor
                std = np.std(datacumulated,axis=0)
                #print std
                #print "ste  ++++++++++++++++++++++++++++",searchfor,"################",len(found)
                ste = np.std(datacumulated,axis=0)/np.sqrt(len(found))
                #print ste


                #print "ste  part +++++++++++++++++++++++",searchfor
                #toget = np.std(datacumulated,axis=0)/np.sqrt(len(found))
                #print "len(ste.shape):",len(ste.shape)
                if len(ste.shape) == 2:
                    #print "------> mean[:,2],ste[:,1]:",mean[:,2],ste[:,1]
                    #mean[:,2] = mean[:,6]+ (ste[:,1]+std[:,1])/2.  # error up
                    #mean[:,3] = mean[:,7]+ (ste[:,1]+std[:,1])/2.  # error down
                    mean[:,2] = max[:,1]-mean[:,1] # error up
                    mean[:,3] = -(min[:,1]-mean[:,1])# error down
                if len(ste.shape) == 1:
                    #print "------>   mean[2],ste[1]  :",mean[2],ste[1]
                    #mean[2] = mean[2]+ (ste[1]+std[1])/2.       # error up
                    #mean[3] = mean[3]+ (ste[1]+std[1])/2.       # error down
                    mean[2] = max[1]-mean[1] # error up
                    mean[3] = -(min[1]-mean[1]) # error down
                #print "mean  again  +++++++++++++++++++++",searchfor
                #print mean
                np.savetxt(filename,mean,fmt="%.5f")
                if len(ste.shape) == 1: # more than one line
                    #np.savetxt(filename,np.transpose(mean),fmt="%.5f")
                    #print "mean:",mean
                    np.savetxt(filename,np.array([mean]),fmt="%.5f")
                print
        pass
    return

def qpointstring(q1,q2=False,q3=False):
    ''' creates string for the corresponding qpoint for files to write out'''
    #print "q1:",q1,type(q1)
    #print "q2:",q2,type(q2)
    #print "q3:",q3,type(q3)
    qpointstr=False
    if type(q2)==bool and type(q3)==bool:
        #print "1Xx"
        if type(q1)==np.ndarray or type(q1) == list:
            #print "2Xx"
            if len(q1) == 3:
                qpointstr = str(q1[0])+"_"+str(q1[1])+"_"+str(q1[2])
    else:
        qpointstr = str(q1)+"_"+str(q2)+"_"+str(q3)
    #print qpointstr
    #sys.exit()
    return qpointstr

def get_all_qpoints(listeqpoints,args=False,verbose=False):
    '''
    [ 1 0 0 ] # particular qpoint
    [ 3 3 3 ] # particular qpoint

    [ l 0 0 ] # 100-L, get all long in this direction
    [ l l 0 ] # 110-L, get all long in this direction
    [ l l l ] # 111-L, get all long in this direction

    [ t 0 0 ] # 100-T, get all trans in this direction
    [ t t t ] # 111-T, get all trans in this direction
    [ t1 t1 0 ] # 110-T1, get all trans1 in this direction
    [ t2 t2 0 ] # 110-T2, get all trans2 in this direction
    '''
    N = args.supercell
    qpoints_all =[]
    appropratestring = [ 'l', 't', 't1', 't2', 'all', 'lt', 'tnew', 'gar' ]
    #print listeqpoints,type(listeqpoints)
    # DEFAULT QPOINTS
    if listeqpoints is None:
        listeqpoints = [['all','all','all']]
    for N1N2N3 in listeqpoints:
        #print "iiiiiiiii:",N1N2N3
        N1 = N1N2N3[0]
        N2 = N1N2N3[1]
        N3 = N1N2N3[2]

        #print "type(N1):",N1,type(N1)
        #print "str(N1) :",N1,int(N1)
        #print "type(N2):",N2,type(N2)
        #print "str(N2) :",N2,int(N2)
        #print "type(N3):",N3,type(N3)
        liste = [N1,N2,N3]
        for idx,i in enumerate(liste):
            #print "i:",i,type(i)
            if type(i) != str and type(i) != int and type(i) != float:
                sys.exit("ERROR: qpoint needs to be a string, integer or float!")


            try:
                liste[idx]=int(i)
            except ValueError:
                if i not in appropratestring:
                    sys.exit("ERROR: qpoint as string should one of \"l,t,t1,t2\" but is \""+i+"\"")

        N1 = liste[0]
        N2 = liste[1]
        N3 = liste[2]

        if type(N1) != str and type(N2) != str and type(N3) != str:
            #qpoints_all = [[N1,N2,N3]] # we have an particular qpoint
            qpoints_all.append([N1,N2,N3])

        # in case however we had a string in there
        N = int(N)  # supercellsize
        all = False
        l100 = False
        t100 = False
        t100new = False
        l110 = False
        t1_110 = False
        t2_110 = False
        l111 = False
        t111 = False
        gar = False
        #gar = True
        if N1 == 'all' or N2 == 'all' or N3 == 'all': all = True            # [ All ]
        if N1 == 'gar' or N2 == 'gar' or N3 == 'gar': gar = True            # [ gar ]


        # [ 1 0 0 ]
        if (N1 == 'l' and N2 == 0. and N3 == 0.): l100 = True               # [ N 0 0 ] (L)
        if (N1 == 't' and N2 == 0. and N3 == 0.): t100 = True               # [ N 0 0 ] (T)
        if (N1 == 'tnew' and N2 == 0. and N3 == 0.): t100new = True               # [ N 0 0 ] (T)
        if (N1 == 'lt' and N2 == 0. and N3 == 0.): l100 = True;t100 = True  # [ N 0 0 ] (L+T)


        # [ 1 1 0 ]
        if (N1 == 'l' and N2 == 'l' and N3 == 0.): l110 = True               # [ N N 0 ] (L)
        if (N1 == 't1' and N2 == 't1' and N3 == 0.): t1_110 = True               # [ N N 0 ] (L)
        if (N1 == 't2' and N2 == 't2' and N3 == 0.): t2_110 = True               # [ N N 0 ] (L)

        # [ 1 1 1 ]
        if (N1 == 'l' and N2 == 'l' and N3 == 'l'): l111 = True               # [ N N N ] (L)
        if (N1 == 't' and N2 == 't' and N3 == 't'): t111 = True               # [ N N N ] (L)
        if (N1 == 'lt' and N2 == 'lt' and N3 == 'lt'): l111 = True;t111 = True # [ N N N ] (L)

        def check_if_to_add(args,qpoint,qpointsall,verbose=False):
            ''' add or not '''
            #print "args.exclude_existing:",args.exclude_existing
            if args.exclude_existing == False:
                #print "in False"
                if qpoint not in qpoints_all:
                    if verbose: print qpoint,args.exclude_existing,"->",True,len(qpointsall)
                    return True
                else:
                    if verbose: print qpoint,args.exclude_existing,"->",False,len(qpointsall)
                    return False
            elif args.exclude_existing == True:
                check1 = "space_fft_"+qpointstring(qpoint)+".npy"
                check = os.path.isfile(check1)
                # if file exists, dont add
                if check:
                    if verbose: print qpoint,args.exclude_existing,"->",check,"->",False,len(qpointsall)
                    return False
                else:
                    if verbose: print qpoint,args.exclude_existing,"->",check,"->",True,len(qpointsall)
                    return True
            else:
                sys.exit("Error 789")
            return

        def get_q_length(q1,q2,q3,N,struct=False):
            ka = np.array([int(q1),int(q2),int(q3)])/float(N)
            #if struct=='fcc' or st
            if True:
                ka2=ka%2
                ka22 = ka2[np.nonzero(ka2)[0]]
                ka3=ka%1
                ka4 = ka3[np.nonzero(ka3)[0]]
                #print "kkk:",ka,ka2,ka22,ka3,ka4
                qlen = ka22.min()  # this is however not working for ttt path
                # for t_t_t
                if q1==q2 and q3 != 0 and q3 != q1 and q3 != 2*N and struct=='fcc':
                    qlen = ka4.min()
                #if in1==in2==in3=='t':
                #    qlen = ka4.min()

                return qlen
            #elif struct=='bcc':
            #    print "kkk:",ka,ka%2

        symeq= []




        def qpoints_all_append(qp,add="",verbose=verbose):
            cr = check_qpoints_for_crossing(qp, args.supercell, args.structure)
            co=""
            if cr:
                co=u'\u2713'
                co="CROSSING"
            #print "cr:",cr
            qpoints_all.append(qp)
            i = qp  # i[0], i[1], i[2]
            ############ get length of qpoint
            #print "qp:",qp
            ja = 3
            symeq.append(equivalent_qs(i).shape[0])
            if verbose:
                print "    "+add.ljust(13),"[",str(i[0]).ljust(ja),str(i[1]).ljust(ja),str(i[2]).ljust(ja),"]",qpointstring(i).ljust(10),"("+str(equivalent_qs(i).shape[0])+")".ljust(5),'\t',str(len(qpoints_all)).ljust(4),"("+str(np.sum(symeq))+")".ljust(3),'\t',"",str(get_q_length(q1=i[0],q2=i[1],q3=i[2],N=N,struct=args.structure)),'\t',co
            return
        #print "agg:",args.showeq
            #if args.showlammpsinput or args.create_lammps_inputfile: print equivalent_qs(i,show_lammps_inputfile=True,write_lammps_inputfile=args.create_lammps_inputfile,N=args.supercell,alat=args.alat,dt=args.dt)

        #print "all:",all
        if verbose:
            print "qpoints_all  :     [ qpoint      ] qpstr     symeq     Nr. symeqtot     Reduced wave vector"
            print "-------------------------------------------------------------------------------------------"
        standard = True
        if gar == True:
            standard = False
        if standard:
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 0 0 ] >>>>>>>>>>>>>>>>>>>>>>
            if l100 or all == True:                                         # [ N 0 0 ] (L)
                for i in np.arange(N)+1:
                    qp = [i,0,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (L)")
                if verbose: print
            if t100 or all == True:                                         # [ N 0 0 ] (T)
                for i in np.arange(N)+1:
                    qp = [2*N,i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (T)")
                if verbose: print
            #if t100new or all == True:                                         # [ N 0 0 ] (T)
            #    for i in np.arange(N)+1:
            #        qpoints_all.append([2*N,i,2*N])
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 0 0 ] >>>>>>>>>>>>>>>>>>>>>>

            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 0 ] >>>>>>>>>>>>>>>>>>>>>>
            if args.structure == 'fcc': NN = N
            if args.structure == 'bcc': NN = N/2
            if l110 == True or all == True:      # [ N N 0 ] (L)
                for i in np.arange(NN)+1:
                    qp = [i,i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (L)")
                if verbose: print
            if t2_110 == True or all == True:    # [ N N 0 ] (T2)
                for i in np.arange(NN)+1:
                    qp = [i,i,2*N]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (T2)")
                if verbose: print
            if t1_110 == True or all == True:    # [ N N 0 ] (T1)
                for i in np.arange(NN)+1:
                    qp = [2*N+i,2*N-i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (T1)")
                if verbose: print
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 0 ] >>>>>>>>>>>>>>>>>>>>>>


            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 1 ] >>>>>>>>>>>>>>>>>>>>>>
            if l111 == True or all == True:     # [ N N N ] (L)
                if args.structure == 'fcc':
                    for i in np.arange(N/2)+1:
                        qp = [i,i,i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (L)")
                    if verbose: print
                if args.structure == 'bcc':
                    for i in np.arange(N)+1:
                        qp = [i,i,i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (L)")
                    if verbose: print
            if t111 == True or all == True:     # [ N N N ] (T)
                if args.structure == 'fcc':
                    for i in (np.arange(N/2)+N/2)[::-1]:
                        qp = [i,i,2*N-i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (T)")
                    if verbose: print
                if args.structure == 'bcc':
                    for i in (np.arange(N)+1): #[::-1]:
                        qp = [i,i,2*N-i]
                        if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (T)")
                    if verbose: print
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 1 ] >>>>>>>>>>>>>>>>>>>>>>

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ gar ] >>>>>>>>>>>>>>>>>>>>>>


        #print "gar1:",gar,len(qpoints_all)
        #if gar == True:
        if gar:
            if l100 or all == True:                                         # [ N 0 0 ] (L)
                for i in np.arange(N)+1:                    # [ N 0 0 ] (L)
                    qp = [2*N+i,0,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (L)")
                if verbose: print
            if t100 or all == True:                                         # [ N 0 0 ] (L)
                for i in np.arange(N)+1:                    # [ N 0 0 ] (T)
                    qp = [2*N,2*N,i]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N 0 0 (T)")
                if verbose: print
            if l110 == True or all == True:      # [ N N 0 ] (L)
                for i in np.arange(N)+1:                    # [ N N 0 ] (L)
                    qp = [2*N+i,2*N+i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (L)")
                if verbose: print
            if t2_110 == True or all == True:    # [ N N 0 ] (T2)
                for i in np.arange(N)+1:                    # [ N N 0 ] (T2)
                    qp = [i,i,-2*N]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (T2)")
                if verbose: print
            if t1_110 == True or all == True:    # [ N N 0 ] (T2)
                for i in np.arange(N)+1:                    # [ N N 0 ] (T1)
                    qp = [2*N+i,2*N-i,0]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N 0 (T1)")
                if verbose: print
            if l111 == True or all == True:     # [ N N N ] (L)
                for i in np.arange(N/2)+1:                  # [ N N N ] (L)
                    qp = [N+i,N+i,N+i]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (L)")
                if verbose: print
            if t111 == True or all == True:     # [ N N N ] (L)
                for i in np.arange(N/2)+1:                  # [ N N N ] (T)
                    qp = [N-i,N-i,N+i]
                    if check_if_to_add(args,qp,qpoints_all): qpoints_all_append(qp,"N N N (T)")
                if verbose: print
        #print "gar2:",gar,len(qpoints_all)
        if len(qpoints_all) == 0:
            sys.exit("ERROR: this path is not known!")
        #print "iiiiiiiii:",N1N2N3,"done"
    #print "iftrue:",qpoints_all[0] == qpoints_all[1],qpoints_all[0],qpoints_all[1]

    #print "qpoints_all  : [ qpoint ] \tqpstr\teq\tsum                "
    if verbose:
        print "-------------------------------------------------------------------------------------------"
    qpstring_all = [];tmp=0;
    for idx,i in enumerate(qpoints_all):
        tmp+=equivalent_qs(i).shape[0]
        qpstring_all.append(qpointstring(i))
        #if args.verbose > 1:
        #    print "              ",i,'\t',qpointstring(i),'\t',equivalent_qs(i).shape[0],'\t',tmp
        #else:
        #    if idx < 3 or idx > len(qpoints_all)-4:
        #        print "              ",i,'\t',qpointstring(i),'\t',equivalent_qs(i).shape[0],'\t',tmp
        #    if idx >= 3 and len(qpoints_all)-4 >= 6 and idx == 4:
        #        print "                  ...            "
        #    if idx >= 3 and len(qpoints_all)-4 >= 6 and idx == 5:
        #        print "                  ...            "
        #    if idx >= 3 and len(qpoints_all)-4 >= 6 and idx == 6:
        #        print "                  ...            "
        #qpstring_all.append(qpointstring(i))
    if verbose:
        print qpstring_all
        print "-------------------------------------------------------------------------------------------"
    if args.showeq:
        if verbose:
            print "-------------------------------------------------------------------------------------------"
        for idx,i in enumerate(qpoints_all):
            if args.showeq: print equivalent_qs(i,show_without_ipi=True)
        if verbose:
            print "-------------------------------------------------------------------------------------------"
    return qpoints_all

def print_and_save_lifetimesgood(args = False,verbose=True,printtoscreen=True): #alat, in1, in2, in3, dt):
    ''' keyword is
            - lifetimesgood l 0 0
            - lifetimesgood t2 t2 t2 mev
            - lifetimesgood l l l meV
            - freqsgood t1 t1 0 meV
    '''
    #if args.qvec[0]!='all' and args.qvec[1]!='all' and args.qvec[2]!='all':
    #    return # at least one needs to be all
    print "-------------------------------------------------------"
    print "---------------- def print_and_save_lifetimesgood() ---"
    print "-------------------------------------------------------"
    dispdir_all = [['l','0','0'],['t','0','0'],['l','l','0'],['t1','t1','0'],['t2','t2','0'],['l','l','l'],['t','t','t']]
    for dispdir in dispdir_all:
        in1, in2, in3 = dispdir[0],dispdir[1],dispdir[2]
        #qpoints_all = get_all_qpoints(['all','all','all'],args)
        #qpoints_all = get_all_qpoints(args.qvec,args)
        #print "in1;",in1
        qpoints_all = get_all_qpoints([[in1, in2, in3]],args)  # dont change this, it is correct
        for keyword in [ 'freqsgood', 'lifetimesgood' ]:
            #for keyword2 in [ 'THz', 'meV' ]:
            for keyword2 in [ 'THz' ]:
                if keyword2 == 'THz':units = 'THz'
                if keyword2 == 'meV':units = 'meV'
                #keyword2 = ""
                #print "----len(sys.argv):",len(sys.argv)
                #if len(sys.argv) > 9:
                #    keyword2 = sys.argv[9]  # {mev,meV}
                #print "----keyword2     :",keyword2



                mdsteps = []
                THz_to_meV = 1.0
                if keyword2 == "mev" or keyword2 == "meV":
                    #print "ioooooooooooo"
                    THz_to_meV = 4.1356673


                #files =  glob.glob("ps*.dat."+keyword+".dat")
                #print files
                #print
                files = []

                #for i in qpoints_all:
                #    print i,qpointstring(i)
                #print
                #print "keyword:",keyword
                #print

                for i in qpoints_all:
                    appd = glob.glob("ps"+qpointstring(i)+"*.dat."+keyword+".dat")
                    if len(appd) > 0:
                        files.append(glob.glob("ps"+qpointstring(i)+"*.dat."+keyword+".dat")[0])
                    else:
                        print("ERROR: no files foud with name: ps"+qpointstring(i)+"*.dat."+keyword+".dat")
                        return



                #print files
                #sys.exit()
                for file in files:
                    #print file
                    f1 = file.split("_")
                    f2 = f1[3].split(".dat")[0]
                    #print "f2:",f2,float(f2),int(float(f2))
                    mdsteps.append(int(float(f2)))  # keep int(float())

                mdstepsmax = np.sort(np.array(list(set(mdsteps))))[-1]
                if printtoscreen:
                    print "-------------------------------------------------------"
                    print "[",in1,in2,in3,"]",keyword,"("+str(units)+str(")")," || mdstepsmax:",mdstepsmax
                import re
                numbers = re.compile(r'(\d+)')
                def numericalSort(value):
                    parts = numbers.split(value)
                    parts[1::2] = map(int, parts[1::2])
                    return parts

                #print mdstepsmax
                files = []
                for i in qpoints_all:
                    files.append(glob.glob("ps"+qpointstring(i)+"_"+str(mdstepsmax)+".dat."+keyword+".dat")[0])
                #files =  sorted(glob.glob("*"+str(mdstepsmax)+".dat."+keyword+".dat"), key=numericalSort)

                #for file in files:
                #    print file
                def make_phonon_disperion_pic(in1,in2,in3):
                    ''' convert qpointstrin to info for picture '''

                    # 1 0 0
                    if in2 == "0" and in3 == "0":   # (T)
                        if in1 == "t":
                            out = q2;faktor=1.;rd=0;sc=-1;  # rd=0 and sc=-1 cancells out

                        if in1 == "l":              # (L)
                            out = q1;faktor=1.;rd=0;sc=-1.; # rd=0 and sc=-1 cancells out

                    # 1 1 0
                    if in1 == "l" and in2 == "l" and in3 == "0":
                        out = q1;faktor=1.;rd=0.;sc=-1.;
                    if in1 == "t1" and in2 == "t1" and in3 == "0":
                        out = q1;faktor=1.;rd=(2*args.supercell);sc=-(1./args.supercell);      # faktor =3 and rd=3 and sc=1 scheint zu stimmen (faktor 3 da bis zu [ 30 10 0 ]
                    if in1 == "t2" and in2 == "t2" and in3 == "0":
                        #out = q1;faktor=1.;rd=-1.;sc=-1;      # faktor =1 and rd=1 and sc=1 scheint zu stimmen
                        out = q2;faktor=1.;rd=0;sc=-1;      # faktor =1 and rd=1 and sc=1 scheint zu stimmen


                    # 1 1 1
                    # print 'ag;',args.supercell % 2 == 0 # odd -> False
                    # print 'agg',args.supercell % 2 == 1 # odd -> True
                    # 9_9_11, 8_8_12 .. 5_5_15 for N = 10
                    if in1 == "l" and in2 == "l" and in3 == "l":  # rd=0 and sc=-1 cancells out
                        # 1_1_1 bis 5_5_5 for N = 10
                        out = q3;faktor=1.;rd=0.;sc=-1.;
                        # out = 1..5 (L)
                        #print "out l:",out
                        if args.supercell % 2 == 1: # odd -> True
                            sc = sc*0.5

                    if in1 == "t" and in2 == "t" and in3 == "t":  # rd=0 and sc=-1 cancells out
                        #out = q3;faktor=1.;rd=args.supercell;sc=-2./args.supercell;
                        out = int(q3)-args.supercell;faktor=1.;rd=0.;sc=-1.;
                        #print "out t:",out
                        # out:11..15 (T)
                        if args.supercell % 2 == 1: # odd -> True
                            sc = sc*0.5

                    # qlen = sc*(rd-(float(out)*faktor/float(qmax)))
                    return out,faktor,rd,sc

                qmax = 0;
                for file in files:
                    #print  file.split("ps")[1].split("_") #[0]
                    q1 = file.split("ps")[1].split("_")[0]
                    q2 = file.split("ps")[1].split("_")[1]
                    q3 = file.split("ps")[1].split("_")[2]
                    #print q1,q2,q3

                    # get qmax
                    out,faktor,rd,sc = make_phonon_disperion_pic(in1,in2,in3)
                    if int(out) > qmax:
                        qmax = int(out)
                    if in1 == "t1" and in2 == "t1" and in3 == "0":
                        qmax=1
                    #if in1 == "t" and in2 == "t" and in3 == "t":
                    #    qmax=1

                    #print "out:",out,"qmax:",qmax
                    #sc*(rd-(float(out)*faktor/float(qmax)))

                #print "keyword:-----------",keyword
                #print "keyword2-----------:",keyword2
                if printtoscreen:
                    print "-------------------------------------------------------"
                f = open(str(keyword)+"_"+str(in1)+"_"+str(in2)+"_"+str(in3)+"_alat"+str(args.alat)+"_dt"+str(args.dt)+"_"+str(units)+'.dat', 'wb')
                f.write("#q\twidth\terr+\terr-\tsigma\t999999999___qpoint\terr_unt\terr_ob\terr_md\n")
                #outlifetimes = np.array([[mdsteps, ltgood,error_pos,error_neg,sigmagood,999999,error_gaussian_lt_oben_ist_max,error_gaussian_lt_unten_ist_min,lifetimestocheckerrormd[-1]]])
                # idx                         0        1       2         3         4       5          6                                  7                              8

                for file in files:
                    q1 = file.split("ps")[1].split("_")[0]
                    #print file,q1
                    q2 = file.split("ps")[1].split("_")[1]
                    q3 = file.split("ps")[1].split("_")[2]
                    lastline = np.loadtxt(file)
                    #print "file    :",file
                    #print "lastline:",lastline
                    #print "q1:",q1,"q2:",q2,"q3:",q3


                    out,faktor,rd,sc = make_phonon_disperion_pic(in1,in2,in3)


                    #print "out:",out,"qmax:",qmax,"faktor:",faktor,"rd:",rd,"sc:",sc,"qmax:",qmax

                    #print lastline

                    roundto=4
                    #print  str(round(sc*(rd-(float(out)*faktor/float(qmax))),4)),"\t",round(lastline[1]*THz_to_meV,roundto),"\t",round(lastline[2]*THz_to_meV,roundto),"\t",round(lastline[3]*THz_to_meV,roundto),"\t",round(lastline[4],7)
                    #qlen = sc*(rd-(float(out)*faktor/float(qmax)))


                    def get_q_length(q1,q2,q3,N):
                        ka = np.array([int(q1),int(q2),int(q3)])/float(N)
                        ka2=ka%2
                        ka22 = ka2[np.nonzero(ka2)[0]]
                        ka3=ka%1
                        ka4 = ka3[np.nonzero(ka3)[0]]
                        qlen = ka22.min()  # this is however not working for ttt path
                        # for t_t_t
                        if q1==q2 and q3 != 0 and q3 != q1 and q3 != 2*N and args.structure=='fcc':
                            qlen = ka4.min()
                        #if in1==in2==in3=='t':
                        #    qlen = ka4.min()
                        return qlen

                    #ka = np.array([int(q1),int(q2),int(q3)])/float(N)
                    #ka2=ka%2
                    #ka22 = ka2[np.nonzero(ka2)[0]]
                    #ka3=ka%1
                    #ka4 = ka3[np.nonzero(ka3)[0]]
                    #qlen = ka22.min()  # this is however not working for ttt path
                    #if in1==in2==in3=='t':
                    #    qlen = ka4.min()
                    #print "N1:",N1
                    factplot = 1
                    if args.structure == 'bcc':
                        #print "99999999999999999999999999999",q1,q2,q3,type(q3),in1,in2,in3
                        # This scaling is for dominiques plots!
                        factplot = ""
                        factshift = ""
                        if in1 == 'l' or in1 == 't':
                            if in2 == in3 == '0':
                                factplot = 146
                                factshift = 0
                        if in1 == in2 == in3:
                            if in1 == 'l' or in1 == 't':
                                factplot = -1.*(396 - 146)
                                factshift = 396
                        if in1 == in2 and in3 == '0':
                            factshift = 396
                            factplot = (498 - 396)*2.
                        #print "-->",factplot,factshift


                    qlen = get_q_length(int(q1),int(q2),int(q3),N)
                    if args.structure == 'bcc':
                        qlen = factshift + factplot*qlen
                    #ka = str(qlen)+"\t"+\
                    ka = str(round(qlen,4)).ljust(8)+"\t"+\
                            str(round(lastline[1]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[2]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[3]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[4],7))+"\t"+\
                            str(999999999)+"___"+str(q1)+"_"+str(q2)+"_"+str(q3)+"\t"+\
                            str(round(lastline[6]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[7]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[8]*THz_to_meV,roundto))+"\t"+\
                            "\n"
                    ktoscreen = "\t["+str(q1)+" "+str(q2)+" "+str(q3)+"] qmax:"+str(qmax)+"\n"
                    if printtoscreen:
                        print ka.rstrip()+ktoscreen.rstrip()
                    f.write(ka)
                f.close()
    return

def equivalent_qs(q,show_without_ipi=False,show_lammps_inputfile=False,write_lammps_inputfile=False,N=False,alat=False,dt=False):
    ''' - finds all the equivalent q vectors for the given q vector
        - q can be a list or a numpy array e.b. [2,0,0] or np.array([2,0,0])
    '''
    #print "type(q):",type(q)
    if type(q) == list:
        q = np.array(q)
    if type(q) != np.ndarray:
        sys.exit("q needs to be a list or a numpy array!")
    if len(q) != 3:
        sys.exit("q needs to be of length 3")

    ergebnis=[];
    dum1=np.array(list(permutations(q)))
    for i in range(0,dum1.shape[0]):
        for m in range(8):
            dum4=m;
            dum31=(-1)**(dum4%2);
            dum4=(dum4-(dum4%2))/2;
            dum32=(-1)**(dum4%2);
            dum4=(dum4-(dum4%2))/2;
            dum33=(-1)**(dum4%2);
            dum3=(dum31, dum32, dum33)
            ergebnis.append((dum3*dum1[i]))
    a=np.array(ergebnis)
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)

    equivalents = a[idx]
    #print equivalents
    #print equivalents*1j*2*np.pi

    # fixedq: [ 1 0 0 ] and not general  (== old)
    #if q[1] != 0:
    #    sys.exit("currently for [N,0,0] implemented")
    #if q[2] != 0:
    #    sys.exit("currently for [N,0,0] implemented")
    #N = q[0]
    ## aus qvector [N,0,0] wird:
    #equivalents_fixq = np.array([[N,0,0],[0,N,0],[0,0,N],[0,0,-N],[0,-N,0],[-N,0,0]])
    #return equivalents_fixq*1j*2*np.pi
    #return equivalents*1j*2*np.pi   # has + - entries

    ########################################################
    # remove all i = -i
    ########################################################
    out = []
    for i in list(equivalents):
        i=list(i)
        j=[i[0]*-1,i[1]*-1,i[2]*-1]

        #print "OUT:",out
        #print "i:",i,type(i),"j:",j,"(i in out):",i in out,"(j in out):",j in out,"OUT:",out,type(i)
        if (i in out) or (j in out):
            #print "exists"
            continue
        else:
            if len(out) == 0:
                out = [i]
            else:
                out.append(i)
            #print i,"-------------",out #i in out,"_______________",out
        #print
    out= np.array(out)
    if show_without_ipi == True:
        return out


    if show_lammps_inputfile == True or write_lammps_inputfile:
        if type(N) == bool:
            sys.exit('ERROR: you need to specify N if show_lammps_inputfile')
        print out.shape  # z.b. 12,3
        outstr = []
        maxvalue = 0
        qs_computes = []
        if write_lammps_inputfile:
            runline = os.popen('grep "^run" in_file_dynamics.in').read().rstrip()
            os.popen("sed -i 's|^run.*||' in_file_dynamics.in").read().rstrip()
            print "runline:",runline
            target = open("in_file_dynamics.in", 'a')
        for idx,sym_qp in enumerate(out):
            # x
            if sym_qp[0] == 0: o1 = ""
            if sym_qp[0] != 0:
                o1 = str(sym_qp[0])+"*x+";
                if maxvalue < sym_qp[0]: maxvalue = sym_qp[0]
            # y
            if sym_qp[1] == 0: o2 = ""
            if sym_qp[1] != 0:
                o2 = str(sym_qp[1])+"*y+"
                if maxvalue < sym_qp[1]: maxvalue = sym_qp[1]

            # z
            if sym_qp[2] == 0: o3 = ""
            if sym_qp[2] != 0:
                o3 = str(sym_qp[2])+"*z+"
                if maxvalue < sym_qp[2]: maxvalue = sym_qp[2]

            #print idx,sym_qp
            #for xyzidx,xyz in enumerate(sym_qp):

            #if out[idx][idx] == 0:
            #    outstr.append(" ")
            #else:
            #    outstr.append(str(out[idx][idx]))
            out=(o1+o2+o3).replace("+-", "-")[:-1]
            #out2='(('+out+')/${alat}*2.0*PI/'+str(maxvalue)+')'
            #out2='(('+out+')/${alat}*2.0*PI/'+str(N)+')'
            out2='(('+out+')/'+str(alat)+'*2.0*PI/'+str(N)+')'
            #print idx,sym_qp

            ########################################################################
            # schreibe variablen
            ########################################################################
            realcom = [ 'r','i' ]
            for idx2,i in enumerate([ 'cos', 'sin' ]):
                variablename = qpointstring(q)+'_'+str(idx+1)+str(realcom[idx2])
                write1 = "variable "+variablename+" atom "+i+out2
                print write1
                if write_lammps_inputfile:
                    target.write(write1+"\n")
                qs_computes.append(variablename)
        ########################################################################
        # schreibe fuer jede variable ein compute
        ########################################################################
        for variablename in qs_computes:
            write2 = "compute c"+variablename+" all reduce sum v_"+variablename
            if write_lammps_inputfile:
                target.write(write2+"\n")
        for variablename in qs_computes:
            write3 = "variable "+variablename+"_ equal c_c"+variablename
            if write_lammps_inputfile:
                target.write(write3+"\n")
        print
        write4 = 'fix '+qpointstring(q)+' all print '+str(dt)+' "'+ ' '.join(['${'+s+'_}' for s in qs_computes]) + '" append space_fft_'+qpointstring(q)+".dat screen no"
        if write_lammps_inputfile:
            target.write(write4+"\n")
            target.write(runline+"\n")
            target.close()
        print
        print 'thermo_style custom '+ ' '.join(['c_c'+s for s in qs_computes])
            #print idx,sym_qp,"outstr:","cos"+out2
            #print "           outstr:","sin"+out2
        #print "cos("+outstr[0]+"*x+",outstr[1],"*y+",outstr[2],"*z)"
        return
    return out*1j*2*np.pi

def print_warning(text):
    print bcolors.OKGREEN + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.OKGREEN + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.WARNING + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.FAIL + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.FAIL + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.FAIL + text +bcolors.ENDC
    print bcolors.FAIL + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.FAIL + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.WARNING + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.OKGREEN + "##########################################################################################################################################################################################" +bcolors.ENDC
    print bcolors.OKGREEN + "##########################################################################################################################################################################################" +bcolors.ENDC

def mysmooth(l,N):
    if (N % 2 == 0): #even
        sys.exit("N has to be odd")
    if N <=2:
        sys.exit("make N > 2")

    d = (N-1)/2
    #print "d:",d
    #print "l:",l
    #print
    out = np.copy(l)
    for idx,i in enumerate(l):
        neg=idx-d
        if neg<0:
            neg=0
        #print idx,i,"window:",l[neg:idx+d+1]
        out[idx] = np.mean(l[neg:idx+d+1])
        #if idx >= d and idx <= len(l)-d:
        #    print idx,i,"window:",l[idx-d:idx+d+1]
        #else:
        #    print idx,i,"window:"
    return out

##########################################################################################
# functions to create space_fft ##########################################################
##########################################################################################
# needs read csv
def get_space_fft_from_positions(filename, qpoints_all = False, lines_of_one_mdstep = False, mdsteps_per_chunk = 1, scale_data = 1,chunkitornot = False,args=args):
    ''' reads a file in chunks to save memory and creates the sum (for a certain list of
        qpoints) which can be directly fourier transformed.

    - per chunk, amount of lines given by 'lines_of_one_mdstep' (=atoms_in_supercell) times
      mdsteps_per_chunk

    - chunkitornot = {False,'chunk','full'}  where full reads full file and chunk reads in in chunks.

    - example:
        filename = pos_lammps
        qpoint = np.array([1,0,0])
        lines_of_one_mdstep = 32    # == number of atoms; for a 2x2x2 fcc supercell
        mdsteps_per_chunk = 100     # this would read 32*100 lines per chunk
        scale_data = 4.13*2         # 4.13 Angstrom * 2; this is necessary to get DIRECT coords
    '''
    #qpoints_all=np.array([[1,0,0],[2,0,0],[3,0,0],[4,0,0]])
    #qpoints_all=np.array([[2,0,0]])
    #qpoints_all=np.array([[2,0,0],[3,0,0]])
    if os.path.isfile(filename) != True: sys.exit("ERROR:"+str(filenmae)+" does not exist")
    print "filename:",filename
    print "scale_data:",scale_data,type(scale_data)
    if type(lines_of_one_mdstep) == bool: sys.exit("ERROR: lines_of_one_mdstep needs to be a number!")


    ############################################################
    # start doint this for al list of qpoints
    # this would only work only when 2*10**8/amountqpoint (say 20) = 1*10**7 is the max number of timesteps
    ############################################################
    qs_all = dict()
    sum_all = dict()
    sum_all_NEW = dict()

    ############################################################
    # start reading in the positions file
    ############################################################
    import time
    start1 = time.time()
    start = time.time()
    if filename == "POSITIONs" or filename == 'pos' or filename == 'trj_lammpsnew.out' or filename == 'trj_lammpsnew.out_noJumps_small':
        print "start reader using c engine!"
        chunksize=lines_of_one_mdstep*mdsteps_per_chunk
        #chunksize=4000*1000
        if chunksize < 1: sys.exit("ERROR: chunksize needs to be greater 0!")
        print "lines_of_one_mdstep:",lines_of_one_mdstep
        print "mdsteps_per_chunk  :",mdsteps_per_chunk
        print "resuling chunksize :",chunksize
        filesize = os.path.getsize(filename)
        print "filesize:",filesize
        if filesize < 1*10**10: # 10 GB
            print "--> smaller 1*10**10 (10 GB) -> no splitting necessary"
            def file_len(fname):
                with open(fname) as f:
                    for i, l in enumerate(f):
                        pass
                return i + 1
            len_file = file_len(filename)
            print "len_file (lines):",len_file
            len_steps = len_file/args.atoms
            print "len_steps:",len_steps
            # cat POSITIONs| awk '{print $1,$2,$3}' > tmp_xyz
            # mv tmp_xyz > POSITIONs
            if os.path.isfile("tmp_xyzz"):
                os.remove("tmp_xyzz")
            os.popen("cat POSITIONs| awk '{print $1,$2,$3}' > tmp_xyzz; mv tmp_xyzz POSITIONs").read()
        else:
            print "WARNING: this POSITIONS/pos/trj_lammpsnew.out/trj_lammpsnew.out_noJumps_small file ist greater than 10GB, consider the c skript using -fftc since this might take a while."
        reader = read_csv(filename, sep=' ', header=None,engine='c',chunksize=chunksize) # 3.7 sec
        print "start reader finished"

    else:
        sys.exit('ERROR: dont know this filename:'+str(filename))


    ############################################################
    # look for existing qpoint files
    ############################################################
    maxavail = np.zeros(len(qpoints_all))
    print "qpoints_all:",qpoints_all
    for idx,qpoint in enumerate(qpoints_all):
        qpointstr = qpointstring(qpoint)
        print "qpoint   :",qpoint,"qpointstr:",qpointstr,"mdsteps_per_chunk:",str(mdsteps_per_chunk)
        #qpointfiles = glob.glob("sum_all_new_"+qpointstr+"__*__"+str(mdsteps_per_chunk)+".npy")
        qpointfiles = glob.glob("space_fft_"+qpointstr+"__*__"+str(mdsteps_per_chunk)+".npy")

        maxavailqpoint = 0
        for j in qpointfiles: # for particular qpoint
            #check = int(j.split("sum_all_new_"+qpointstr+"__")[1].split("__"+str(mdsteps_per_chunk)+".npy")[0])
            check = int(j.split("space_fft_"+qpointstr+"__")[1].split("__"+str(mdsteps_per_chunk)+".npy")[0])
            if check > maxavailqpoint:
                maxavailqpoint=check
        maxavail[idx] = maxavailqpoint
        #print qpointstr, maxavailqpoint,maxavail,int(maxavail.min())
    maxavail = int(maxavail.min())
    print "maxavail:",maxavail
    if maxavail == 0:
        maxavail = -1;  # just so that nothing goes wrong on step 0
    print "maxavail:",maxavail

    ############################################################
    ############################################################
    # start iterating over infile
    ############################################################
    ############################################################
    print "starting iterating over infile ..."
    nextcheck = mdsteps_per_chunk
    xxx_ideal_equilibrium = np.zeros((N**3*args.usestruct,3))
    xxx_ideal_equilibrium = False
    for chunk,pos in enumerate(reader):
        #print "chunk:",chunk,"pos:",pos
        #sys.exit()
        if chunk < maxavail:
            if chunk == mdsteps_per_chunk or chunk == nextcheck:
                nextcheck = chunk + mdsteps_per_chunk
            continue
        #print pos.as_matrix()
        #print "scale_data:",scale_data,type(scale_data)
        ### MAIN !!!!!!!!!!!!
        xxx = pos.as_matrix()/float(scale_data)  # keep the float even if you have already a float
        xxx = np.reshape(xxx,((-1,N**3*args.usestruct,3)))  # shape: (4000, 256, 3) for last step it could be not 3 but 2 oder 1 (32,3)


        ############################### get transversal to first BZ
        ############################### get transversal to first BZ
        if type(xxx_ideal_equilibrium) == bool:
            xxx_ideal_equilibrium = np.round([(xxx[0]*2*N+0.5)%(2*N)-0.5])/2.
        ### durch diese definition sind die auslenkungen in der 4er zelle 4mal zu gros und in der 2er zelle 2mal zu gross, hat dies einen effekt?
        xxx_dr = ((xxx*2*N+0.5)%(2*N)-0.5)/2. - xxx_ideal_equilibrium

        # frage an michael:

        # fuer qpoint [1 0 0]  -> equivalent.qs = [1,0,0] [0,1,0] [0,0,1]
        sigma = np.array([0,1,0])
        sigma = np.array([0,-1,0])
        sigma = np.array([0,0,1])
        sigma = np.array([0,0,-1])
        xxx_phase = np.sum(sigma*xxx_dr,axis=2)
        print "xxx_phase.shape:",xxx_phase.shape
        np.save("space_fft_phase"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),xxx_phase)

        # wie speichert man die? man muss ueber alle phasen druebergehen , aber dass kann man schon hier amchen!
        #@ np.save("space_fft_phase"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),sum_all[idx])

        #np.save("xxx",xxx)
        #np.save("xxx_ideal_equilibrium",xxx_ideal_equilibrium)
        #np.save("xxx_dr",xxx_dr)
        #print "jo"
        #sys.exit()

        mdstepbegin = chunk*mdsteps_per_chunk
        mdstepend = chunk*mdsteps_per_chunk+xxx.shape[0]
        print "chunk (=index):",chunk,"step:",mdstepbegin,"-",mdstepend,"nextcheck:",nextcheck

        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste [ 1,0,0] [2,0,0] ] [3,0,0]
            print "qpoint:",qpoint,len(qpoints_all)
            qs_all[idx] = equivalent_qs(qpoint) # [1 0 0] [0 1 0] [0 0 1]
            sum_all[idx] = np.zeros((len(qs_all[idx]),mdstepend-mdstepbegin),dtype=np.complex64)  # np.complex_32 would also be an option!

            sum_all_NEW[idx] = np.zeros((len(qs_all[idx])*2,mdstepend-mdstepbegin),dtype=np.complex64)  # np.complex_32 would also be an option!
            for ind_qs_pol,qs_pol in enumerate(np.array([[0,1,0],[0,0,1]])): # here qs is the polarisation
                print "qpoint:",qpoint,",qs_new (polarization):",ind_qs_pol,qs_pol
                sum_all_NEW[idx][ind_qs_pol] = np.sum(xxx_phase*np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs_pol]),axis=2)),axis=1)
            print "done"
            sys.exit()
            for ind_qs,qs in enumerate(qs_all[idx]):  # fuer alle symmetrieequivalenten qpoints also fuer [1 0 0] haben wir [1 0 0] [0 1 0] [0 0 1]
                print "qs:",qs
                # calculate via exp when dumping every chunk (amount number of steps)
                ### MAIN !!!!!!!!!!!!!!!!!!!  (inner sum to make one column per md step, outer sum to make one complex number per mdstep)
                eulerforexp = False  # --> gibt beides exakt die gleiche streuamplitude

                # check q und -q: also in /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_32_q_minus_q_auswirkung_timeinversion
                #np.save("xxx",xxx)
                #sys.exit()
                # z.B. qs = np.array([0,0,1])*2*np.pi*1j
                # array([ 0.+0.j        ,  0.+0.j        ,  0.+6.28318531j])
                # z.B. qs = np.array([0,0,-1])*2*np.pi*1j
                # array([ 0.+0.j        ,  0.+0.j        , -0.-6.28318531j])
                #
                #
                # In [19]: qs = np.array([0,0,1])*2*np.pi*1j
                #
                # In [20]: ampl = np.sum(np.exp(np.sum((xxx*qs),axis=2)),axis=1);print ampl;np.save("100ampl",ampl)
                # [ 0.89848075-0.16619753j  0.90832258-0.17126346j  0.91808158-0.176336j   ...,
                #  -0.15524963-0.10834658j -0.15839359-0.1109407j  -0.16157893-0.11330545j]
                #
                # In [21]: qs = np.array([0,0,-1])*2*np.pi*1j
                #
                # In [22]: ampl = np.sum(np.exp(np.sum((xxx*qs),axis=2)),axis=1);print ampl;np.save("100amplmin",ampl)
                # [ 0.89848075+0.16619753j  0.90832258+0.17126346j  0.91808158+0.176336j   ...,
                #  -0.15524963+0.10834658j -0.15839359+0.1109407j  -0.16157893+0.11330545j]
                #
                #print "qpoint:",qpoint
                #print "qs_all[idx][ind_qs]:",qs_all[idx][ind_qs]
                #sys.exit()

                #sum_all_NEW[idx][ind_qs] = np.sum(xxx_ph*np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2)),axis=1)
                #sum_all_NEW[idx][ind_qs] = xxx_ph*np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2))
                #sum_all_NEW[idx][ind_qs] = np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2))

                #kkk = np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2))
                #print "kkk.shap:",kkk.shape,kkk
                #print
                #print "xxx.phase.shape:",xxx_phase.shape,xxx_phase
                #print
                #ab = xxx_phase*kkk
                #print "---->",ab.shape,ab
                #print
                #ac = np.sum(ab,axis=1)
                #print ".....> 2:",ac.shape,ac
                #sys.exit()
                #sum_all_NEW[idx] = np.sum(xxx_phase*np.exp(np.sum((xxx*float(N)*qs_all[idx][ind_qs]),axis=2)),axis=1)

                if eulerforexp == False:
                    sum_all[idx][ind_qs] = np.sum(np.exp(np.sum((xxx*qs_all[idx][ind_qs]),axis=2)),axis=1)

                # calculate via exp when dumping only at the very end
                # sum_all[idx][ind_qs,mdstepbegin : mdstepend] = np.sum(np.exp(np.sum((xxx*qs_all[idx][ind_qs]),axis=2)),axis=1
                # OUTPUT TO SCREEN WHEN NECESSARY
                #print "sum_all:",chunk,idx,ind_qs,sum_all[idx][ind_qs].shape,"<-- anzahl mdsteps"
                #pass

                # calculate via eulers formula (sin, cos) when dumping every chunk (amount number of steps)
                ### MAIN !!!!!!!!!!!!!!!!!!!  (first sum to make one column per md step, second sum to make one complex number per mdstep)
                if eulerforexp == True:
                    a = np.sum((xxx*qs_all[idx][ind_qs]),axis=2)  # np.sum is over all atoms of one particular mdstep to give one column (per atom one number (complex))
                    sum_all[idx][ind_qs].real = np.sum(np.cos(a.imag),axis=1)  # correct ralpart due to euler
                    sum_all[idx][ind_qs].imag = np.sum(np.sin(a.imag),axis=1)  # correct imagpart to euler
                #print "realp. :",chunk,idx,ind_qs,sum_all[idx][ind_qs].shape,"<-- anzahl mdsteps"
                #print "imagp. :",chunk,idx,ind_qs,sum_all[idx][ind_qs].shape,"<-- anzahl mdsteps"
                #print
                #sys.exit()
        #print "----------------",chunk,"start",sum_all[idx].shape
        #print sum_all
        #print "----------------",chunk,"fin",sum_all[idx].shape
        #if chunk == 2:
        #    sys.exit()


        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
            qs_all[idx] = equivalent_qs(qpoint)
            qpointstr = qpointstring(qpoint)
            #np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),sum_all[idx])
            np.save("space_fft_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),sum_all[idx])

        # for every 10-100 chunks, put all qpoints together and remove the last ones.
        # e.g. chunk  0-9  put all together to 9
        # e.g. chunk 10-19 put all together to 19
        if chunk == mdsteps_per_chunk or chunk == nextcheck:

            loadall = np.arange(nextcheck-mdsteps_per_chunk,nextcheck+1)
            #print loadall
            for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
                qpointstr = qpointstring(qpoint)
                #a = np.load("sum_all_new_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
                a = np.load("space_fft_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
                for b in loadall[1:]:
                    #b = np.load("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                    b = np.load("space_fft_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                    a = np.concatenate((a, b), axis=1)
                # this will overwrite the last chunkfile
                #np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
                np.save("space_fft_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
                for b in loadall[:-1]:  # delete all but last cunkfile
                    #os.remove("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                    os.remove("space_fft_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")

            nextcheck = chunk + mdsteps_per_chunk
        #if chunk == 10:
        #    sys.exit()
    # get last steps
    if chunk < nextcheck:
        loadall = np.arange(nextcheck-mdsteps_per_chunk,chunk+1)
        print loadall
        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
            qpointstr = qpointstring(qpoint)
            #a = np.load("sum_all_new_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
            a = np.load("space_fft_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
            for b in loadall[1:]:
                #b = np.load("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                b = np.load("space_fft_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                a = np.concatenate((a, b), axis=1)
            # this will overwrite the last chunkfile
            #np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
            np.save("space_fft_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
            for b in loadall[:-1]:  # delete all but last cunkfile
                #os.remove("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                os.remove("space_fft_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")


    print "DONE! mdstepend:",mdstepend,"------------------------------------------------"
    end = time.time()
    print "#1) read in data :",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    print "ranaming files:"
    for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
        qpointstr = qpointstring(qpoint)
        #file = glob.glob("sum_all_new_"+qpointstr+"__*")
        file = glob.glob("space_fft_"+qpointstr+"__*")
        print len(file),file
        if len(file) == 1:
            #os.rename(file[0], 'sum_all_new_'+qpointstr+'.npy')
            os.rename(file[0], 'space_fft_'+qpointstr+'.npy')

    #for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
    #    qpointstr = qpointstring(qpoint)
    #    sum_all[idx] = sum_all[idx][:,:mdstepend]
    #    #np.save("sum_all_new_"+str(qpoint[0])+str(qpoint[1])+str(qpoint[2]),sum_all)
    #    np.save("sum_all_new_"+qpointstr,sum_all[idx][:,:mdstepend])
    #    #np.savetxt("sum_all.txt",sum_all)
    #    print "sum_all[idx].shape:",sum_all[idx].shape
    #    #for ind_qs,qs in enumerate(qs_all):
    #    #    np.savetxt('sum'+str(qpoint[0])+str(qpoint[1])+str(qpoint[2])+str(ind_qs)+'_'+str(mdstepend)+'.out',sum)
    return

def get_space_fft_prepare_lammps_log(filename):
    ''' this reads directly from lammps.log; comments above and below the output should be no problem '''
    print "starting get_space_fft_from_lammps_log ..."
    #if type(qpoint) == bool:
    #    sys.exit("ERROR qpoint is bool")
    print "filename     :",filename
    skipfirst=False
    from itertools import islice
    with open(filename) as myfile:
        head = list(islice(myfile, 2000))
        for idx,i in enumerate(head):
            if "Memory usage per processor =" in i:
                skipfirst = idx+2
    skipfirst = int(skipfirst)

    print "skipfirst    :",skipfirst
    if type(skipfirst)== bool:
        sys.exit("ERROR: did not find skipfirst")


    skiplast = os.popen('tail -200 '+filename+' | grep -n "Loop time of" | sed \'s|:.*||\'').read()
    skiplast = skiplast.rstrip()
    print "skiplast1    :",skiplast #,":",type(skiplast)
    dothisall = True

    if skipfirst == 0 and skiplast == "" and os.path.isfile(filename+'.head') and os.path.isfile(filename+'.tail'):
        dothisall = False


    if dothisall:
        if skiplast == "":
            sys.exit("ERROR: skiplast not found")
        skiplast = 200 - int(skiplast) + 1
        print "skiplast2    :"+str(skiplast)+":",type(skiplast)
        if type(skiplast) == bool:
            sys.exit("ERROR: did not find skiplast")

        oben_chars = os.popen('head -'+str(skipfirst)+' '+filename+' | wc -c').read()  # writes log.lammps.head
        unten_chars = os.popen('tail -'+str(skiplast)+' '+filename+' | wc -c').read()  # writes log.lammps.head
        #print "oben_chars1  :",oben_chars
        #print "unten_chars1 :",unten_chars
        print "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
        count_oben=oben_chars.rstrip('\n')
        count_unten=unten_chars.rstrip('\n')
        print "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
        print "count_oben   :",count_oben
        print "count_unten  :",count_unten


        # write head (wird ueberschrieben)
        if skipfirst > 0:
            print "write ",filename,".head"
            os.system('head -'+str(skipfirst)+' '+filename+' > '+filename+'.head')  # writes log.lammps.head
        # write tail (wird ueberschrieben)
            print "write ",filename,".tail"
        os.system('tail -'+str(skiplast)+' '+filename+' > '+filename+'.tail')  # writes log.lammps.head



        # remove from the top (works great)
        # this was taken from: http://stackoverflow.com/questions/17330188/remove-first-n-lines-of-a-file-in-place-in-unix-command-line
        print "remove from the top ..."
        os.system('dd if="'+filename+'" bs="'+str(count_oben)+'" skip=1 of="'+filename+'" conv=notrunc')
        print "remove from the top truncate ..."
        os.system('truncate -s "-'+str(count_oben)+'" "'+filename+'"')

        # remove from the bottom
        # dd if=/dev/null of=log.lammps bs=1 seek=$(echo $(stat --format=%s log.lammps ) - $( tail -n1 log.lammps | wc -c) | bc ) #which removes from the bottom
        #os.system('dd if=/dev/null of='+filename+' bs=1 seek=$(echo $(stat --format=%s log.lammps ) - $( '+str(count_unten)') | bc )') #which removes from the bottom
        #for i in np.arange(skiplast):
        #    os.system('dd if=/dev/null of=log.lammps bs=1 seek=$(echo $(stat --format=%s log.lammps ) - $( tail -n1 log.lammps | wc -c) | bc )') #which removes from the bottom
        #
        # faster without the loop
        print "remove from the bottom ..."
        os.system('dd if=/dev/null of=log.lammps bs=1 seek=$(echo $(stat --format=%s log.lammps ) - $( tail -n'+str(skiplast)+' log.lammps | wc -c) | bc )') #which removes from the bottom

    print "get qpoints_read_in ..."
    qpoints_read_in = os.popen('tail -1 '+filename+'.head').read()
    qpoints_read_in = qpoints_read_in.rstrip()
    #print qpoints_read_in
    header = qpoints_read_in
    header = header.split(" ")

    print "nun muss gruperit werden um die space_fft_x_x_x zu speichern.",len(header)
    print header

    list_of_qpoints = []
    for i in header:
        if len(i.split("c_c")) == 2: # nur wenn es das gibt wird es interesant
            qpoint_string_in="_".join((i.split("c_c")[1]).split("_")[0:3])  # gives the ['10', '10', '20'] for every qpoint
            list_of_qpoints.append(qpoint_string_in)

    #print "list_of_qpoints:",list_of_qpoints
    fromm_too = []
    fromm_too_name = []
    for idx,i in enumerate(list_of_qpoints):
        #print idx,i
        if idx==0:
            current = i
            fromm = idx
            fromm_too_name.append(i)
            #print "from",idx
        if idx > 0:
            if i == current:
                continue
            else:
                current = i
                too = idx #-1
                #print "to",idx #-1
                fromm_too.append([fromm,too])
                fromm = idx
                fromm_too_name.append(i)
                #print "from",idx

    #print "to",len(list_of_qpoints) #-1
    fromm_too.append([fromm,len(list_of_qpoints)])

    #print
    #for idx,i in enumerate(list_of_qpoints):
    #    print idx,i
    #for idx,i in enumerate(fromm_too):
    #    print idx,i,fromm_too_name[idx]

    #for i in fromm_too:
    #    print "-->",i[0],i[1],"qp:",header[i[0]],"lll:",list_of_qpoints[i[0]]
    #    filename="space_fft_"+list_of_qpoints[i[0]]+".npy"
    #    filecontent=sum_all_new[:,i[0]:i[1]]
    #    np.savetxt(filename,filecontent)
    #print "lll:|",len(list_of_qpoints)
    #print "kkk:|",len(header)
    #print "???:|",header[98]
    for i in np.arange(len(list_of_qpoints),len(header)):
        #print i,header[i]
        fromm_too.append([i,i+1])
        fromm_too_name.append(header[i])
    #print
    #for idx,i in enumerate(fromm_too):
    #    print idx,i,fromm_too_name[idx]

    return fromm_too, fromm_too_name

# needs read csv
def get_space_fft_from_lammps_log(filename):
    ''' this reads directly from lammps.log; comments above and below the output should be no problem '''

    from_to, from_to_name = get_space_fft_prepare_lammps_log(filename)
    print
    for idx,i in enumerate(from_to):
        print idx,i,from_to_name[idx]
    print
    print "pandas read csv ... c engine (15 GB log.lammps takes 29.15% of the cmmc's memory)"
    reader = read_csv(filename, sep='\s+',header=None,engine='c')
    print "save as matrix"
    reader = reader.as_matrix()
    print "save space_fft ... lines",reader.shape
    lines = reader.shape[0]
    print "ll:",lines
    print
    for idx,i in enumerate(from_to):
        print idx,i,from_to_name[idx]
        filename="space_fft_"+from_to_name[idx]+".npy"
        filecontent=np.transpose(reader[:,i[0]:i[1]])
        columns = filecontent.shape[0]
        print "cc:",columns,filecontent.shape
        if columns > 1:
            print "columns > 1"
            sum_all_new = np.zeros((columns/2,lines),dtype=np.complex64)  # 4*10^8 seems to be max for cluster
            for j in np.arange(columns/2):
                sum_all_new[j] = filecontent[2*j] +1j * filecontent[2*j+1]  #out[:,0] - 1j * out[:,1]
                np.save(filename,sum_all_new)
        else:
            print "columns = 1"
            np.save(filename,filecontent)
            #sum_all_new[5,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,0] - 1j * out[:,1]
        #np.save(filename,filecontent.astype(np.complex64))


    if False:
        sys.exit()
        #print sum_all_new,sum_all_new.shape
        #print
        #tempin = sum_all_new[:,-1]

        ### analysis of hov many datapots should be excluded
        ### -> REUSLT: NOT excluding the first 4 points results in an error of 0.4Kelvin
        ###
        ##@cutoffmax = 40
        ##@out = np.zeros((cutoffmax,3))

        ##@for cutoff in np.arange(cutoffmax):
        ##@    temp = np.copy(tempin[cutoff:])


        ##@    #np.savetxt("temp.dat",temp)
        ##@    temp_average = np.copy(temp)
        ##@    temp_std = np.copy(temp)
        ##@    temp_ste = np.copy(temp)
        ##@    print "Temperature average is",temp.mean(),"Kelvin"
        ##@    for idx,i in enumerate(temp):
        ##@        #print idx,i,temp[:(idx+1)].mean()
        ##@        temp_average[idx]=temp[:(idx+1)].mean()
        ##@        temp_std[idx]=temp[:(idx+1)].std()
        ##@        temp_ste[idx]=temp[:(idx+1)].std()/np.sqrt(idx+1)
        ##@        #if idx == 3:
        ##@        #    sys.exit()
        ##@    #np.savetxt("temp_average.dat",temp_average)
        ##@    #np.savetxt("temp_std.dat",temp_std)
        ##@    #np.savetxt("temp_ste.dat",temp_ste)
        ##@    print cutoff,temp_ste[-1]
        ##@    out[cutoff][0] = cutoff
        ##@    out[cutoff][1] = temp.mean()
        ##@    out[cutoff][2] = temp_ste[-1]
        ##@np.savetxt("out.dat",out)
        ##@print
        ##@print

        #cutoff=4
        #temp = np.copy(tempin[cutoff:])


        ##np.savetxt("temp.dat",temp)
        #temp_average = np.copy(temp)
        #temp_std = np.copy(temp)
        #temp_ste = np.copy(temp)
        #for idx,i in enumerate(temp):
        #    #print idx,i,temp[:(idx+1)].mean()
        #    temp_average[idx]=temp[:(idx+1)].mean()
        #    temp_std[idx]=temp[:(idx+1)].std()
        #    temp_ste[idx]=temp[:(idx+1)].std()/np.sqrt(idx+1)
        #    #if idx == 3:
        #    #    sys.exit()
        #np.savetxt("temp_average.dat",temp_average)
        #np.savetxt("temp_std.dat",temp_std)
        #np.savetxt("temp_ste.dat",temp_ste)
        #print "Temperature average is",temp.mean(),"Kelvin +-",temp_ste[-1]
        #print



        #sys.exit()
        ##sum_all_new = np.zeros((columns,2*10**8),dtype=np.complex64)  # 4*10^8 seems to be max for cluster
    return

def get_space_fft_from_lammps_log_old_but_in_chunks(filename, columns = 8, qpoint = False, has_temperature_column = True):
    ''' - this is for the [N,0,0] qpoint
        - the result (sum_all_new) can however be evaluated without knowledge of N
          if lammps could finc c_1 c_2 this could be automated (at least the skipfirst)
          can python also find the end (of a possibly huge file?)  make it fixed first
    '''
    print "starting lammppace_fft_from_lammps_log_old_but_in_chunks ..."
    if type(qpoint) == bool:
        sys.exit("ERROR qpoint is bool")
    skipfirst = None
    skiplast = None

    skipfirst = None
    skiplast = 27
    print "fi:",filename
    from itertools import islice
    with open(filename) as myfile:
        head = list(islice(myfile, 2000))
        for idx,i in enumerate(head):
            #if "c_1 c_2" in i:
            if "Memory usage per processor =" in i:
                #print idx,i
                skipfirst = idx+1

    print "skipfirst",skipfirst

    checkline = False
    # chck line i+1 to know how many elements
    if skipfirst != None:
        with open(filename) as myfile:
            head = list(islice(myfile, 200))
            for idx,i in enumerate(head):
                if idx == skipfirst+1:
                    checkline = i

    if type(checkline) != bool:
        #print "checkline:",checkline
        checkline = checkline.split(" ") #.remove("' '")
        checkline[-1] = checkline[-1].strip()  # remove '\n'
        checkline = filter(None, checkline)
        #print "checkline:",checkline,
        elements = len(checkline)
        if has_temperature_column == True:
            elements = elements - 1
        # chck line i+1 to know how many elements



    print "read data from ",filename," at some point, if this file gets too large, a reading in chunks may become necessary"
    #print "with 108082514 lines this gives a memory Erroro --> chunk!"

    ###################################
    # new approach read in chunks
    ###################################
    print "read data from lammpsfiel with frequencies done!"
    #sum_all_new = np.zeros((columns/2,2*10**8),dtype=np.complex64)  # 4*10^8 seems to be max for cluster
    sum_all_new = np.zeros((columns,2*10**8),dtype=np.complex64)  # 4*10^8 seems to be max for cluster
    chunksize=100000
    chunksize=3;


    #reader = read csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize,error_bad_lines=False)
    #reader = read csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize)
    #reader = read csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,skipfooter=30)

    # reader for log_for_sum_all
    print "columns:",columns
    print "columns/2:",columns/2,np.arange(columns/2)
    print "sum_all_new.shape",sum_all_new.shape
    #sys.exit()
    for chunk,pos in enumerate(reader):
        #print "chunk    :",chunk
        #print "pos      :",pos
        #print "pos.shape:",pos.shape
        #print "pos.as_matrix():"#,pos.as_matrix()
        #print "pos"


        #!#sys.exit()
        out = pos.as_matrix()
        #print "read line nr:",chunk*chunksize

        #print "kkk:",out
        #print "out.shape:",out.shape,"out[:,0]:",out[:,0]
        #print "chunk:",chunk,"chunksize:",chunksize,"out.shape:",out.shape,"out.shape[0]:",out.shape[0],"|||",chunk*chunksize,chunk*chunksize+out.shape[0]
        #
        #
        # sum_all[ind_qs,mdstepbegin : mdstepend] = np.sum(np.exp(np.sum((xxx*qs_all[ind_qs]),axis=2)),axis=1)
        # np.exp(x+iy) wird in lammps eingesetzt wobei wir aber ueber alle positions gehen.
        # cos(x/${alat}/${N}*2.0*${N}*PI)        == realteil
        # sin(x/${alat}/${N}*2.0*${N}*PI)        == imaginaerteil
        #       wobei /${alat}/${N} das scaling (in direkte coordinaten) ist
        #       *2.0*${N}*PI ist das equivalent_qs (bekommt nur N und macht daraus die qpunkte)
        #

        # ths should be the general case probably
        for idx in np.arange(columns/2):  #  array([0, 1, 2, 3]) for dom and  array([0, 1, 2]) for fcc [N 0 0]
            #print "idx:",idx
            sum_all_new[idx,chunk*chunksize:chunk*chunksize+out.shape[0]]         = out[:,idx*2] + 1j * out[:,idx*2+1]
            #print "idx+col:",idx+columns/2
            sum_all_new[idx+columns/2,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,idx*2] - 1j * out[:,idx*2+1]

        # the next 3 lines word for fcc
        #sum_all_new[0,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,0] + 1j * out[:,1]
        #sum_all_new[1,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,2] + 1j * out[:,3]
        #sum_all_new[2,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,4] + 1j * out[:,5]
        ## (for dominique) sum_all_new[3,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,6] + 1j * out[:,7]

        #####################################################################################################
        ### die naechsten 3 sind eigentlich redundant und sollten eigentlich
        ### erst spaeter bei der berechnung des powerspektrums
        ### beruecksichtigt werden.
        ### DIS GILT FUER DIE LONGITUDINALEN, JEDOCH AUCH FUER DIE TRANSVERSALEN???
        #####################################################################################################
        #sum_all_new[5,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,0] - 1j * out[:,1]
        #sum_all_new[4,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,2] - 1j * out[:,3]
        #sum_all_new[3,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,4] - 1j * out[:,5]

    #print "chunklast!:",chunk,"out.shape[0]",out.shape[0],"--------->",chunk*chunksize+out.shape[0]
    sum_all_new = sum_all_new[:,:(chunk*chunksize+out.shape[0])]
    #print "sum_all_new:"
    #print sum_all_new
    #sys.exit()

    if abs(sum_all_new[:,0].max()) > 10**-13:
        print "first line:",sum_all_new[:,0]
        #sys.exit("first line seems to be too high")
    #qs1 = data[...,0] + 1j * data[...,1]
    #qs2 = data[...,2] + 1j * data[...,3]
    #qs3 = data[...,4] + 1j * data[...,5]
    #qs4 = data[...,4] + 1j * -data[...,5]
    #qs5 = data[...,2] + 1j * -data[...,3]
    #qs6 = data[...,0] + 1j * -data[...,1]

    #print "data.shape:",data.shape
    #print data
    #print "-----------"
    ##sum_all_new = np.zeros((6,data.shape[0]),dtype=np.complex64)
    #i=0
    #for vorzeichen in np.array([+1.0,-1.0]):
    #    for qspair in np.array([[0,1],[2,3],[4,5]]):
    #        #print "vz:",vorzeichen,"qspair:",qspair
    #        #qs = data[...,qspair[0]] + 1j * vorzeichen * data[...,qspair[1]]
    #        #qs = 1j*Data[...,1]; result += Data[...,0]
    #        qs = 1j * vorzeichen * data[...,qspair[1]];qs += data[...,qspair[0]]  # saves memory, no intermediate result
    #        print "qs:",qs,qs.shape,i
    #        sum_all_new[i] = qs
    #        i=i+1
    #print "sum_all_new[0]:",sum_all_new[0]
    #np.save("sum_all_new_zuvor",sum_all_new)
    #sum_all_new = sum_all_new[~np.isnan(sum_all_new).any(axis=1)]
    #print "sum_all_new[0]:",sum_all_new[0]

    print "np.save(sum_all_new)"
    qpointstr = qpointstring(qpoint)
    #np.save("sum_all_new_"+str(qpoint[0])+str(qpoint[1])+str(qpoint[2]),sum_all)
    #np.save("sum_all_new_"+qpointstr,sum_all_new)
    np.save("space_fft_"+qpointstr,sum_all_new)

    #np.save("sum_all_new",sum_all_new)
    return sum_all_new

def get_space_fft_from_xaa_xab_folder(qpstring):
    ''' qpstring is a list of all qpoints which are to be evaluated
    e.b. [ '1_0_0', '2_0_0', ...] '''
    #print "qpstring:",qpstring
    #print "len(qpstring):",len(qpstring),type(qpstring)
    if type(qpstring) == str:
        qpstring = [ qpstring ]
    #def sum_them_all(j):  # where j is one the qpstirng
    foldername="x[a-z][a-z]_"
    filenames_vor_qp="sum_all_new_"
    filenames_nach_qp=".dat"

    folder=sorted(glob.glob(foldername))

    files=sorted(glob.glob(folder[0]+"/"+filenames_vor_qp+"*"+filenames_nach_qp))
    #qpstring= []
    ##print "xx:",qpstring
    #print "----------------folder-------------------"
    #print "len(folder):",len(folder),qpstring

    #for num,i in enumerate(folder):
    #    print num,i
    #print
    #print "----------------files for folder[0]------"
    #for i in files:
    #    print i
    #for i in files:
    #    qpstring.append(i.split(filenames_vor_qp)[1].split(filenames_nach_qp)[0])
    #print qpstring

    #######################################################
    ## start filtering for specific branch
    #######################################################
    ## currently only evaluate [1 1 1] L branch
    #qpstringsubset = []
    #for i in qpstring:
    #    #print "possible:",i
    #    alle = i.split("_")
    #    #print alle,alle[0]==alle[1]==alle[2]

    #    # 1 1 1 (L)
    #    #if alle[0]==alle[1]==alle[2]:
    #    #    qpstringsubset.append(i)

    #    # 1 1 0 (L)
    #    #if alle[0]==alle[1] and alle[2] == '0':
    #    #    qpstringsubset.append(i)

    #    ## 1 1 0 (T1)
    #    #if alle[0]==alle[1] and alle[2] == '20':
    #    #    qpstringsubset.append(i)


    #    ## 1 0 0 (T)
    #    #if alle[0]=='20' and (int(alle[1]) <= 10 and int(alle[1]) >=1) and alle[2]=='0':
    #    #    qpstringsubset.append(i)

    #    #if alle[0]==alle[2]=='20':  # 1 0 0 Tnew
    #    #    qpstringsubset.append(i)
    #    print "alle:",alle[0],alle[1],alle[2],i,"==============",int(alle[0])==20+int(i)
    #    if int(alle[0])==20+int(i) and int(alle[1])==20 -int(i) and alle[2]=='0':
    #        qpstringsubset.append(i)


    #print
    #print "subset:",qpstringsubset
    #qpstring=qpstringsubset
    #print
    #print "out:",qpstring
    #print
    #for idx,i in enumerate(qpstring):
    #    print "in--->>>:",idx,i

    ######################################################
    # start filtering for specific branch ... DONE
    ######################################################

    #theList= range(6)
    #N=2
    #print theList
    #subList = [theList[n:n+N] for n in range(0, len(theList), N)]
    #print subList
    #
    ##print grouper(2,np.arange(6))


    #from multiprocessing import Process
    #def doit(j):
    #    print "all stuff to do with j",j
    #return
    for idxj,j in enumerate(qpstring):  # this needs to be parllelized
        #print "j:",j
        #sys.exit()
    #def sum_them_all(j):
        if os.path.isfile('space_fft_'+j+'.npy'):
            print 'space_fft_'+j,"does already exist...continue"
            continue
        else:
            print 'space_fft_'+j+'.npy',"does not exist --> creating it!"


        #print "idxj:",idxj,"j:",j
        #if j != "1_0_0":
        #    continue
        for idxi,i in enumerate(folder):
            #strout="idxj:",idxj,"j:",j,"         idxi:",idxi,"i:",i
            file=i+"/"+filenames_vor_qp+j+filenames_nach_qp
            #print "importing "+file
            data=np.loadtxt(file)
            #print file,"data.shape:",data.shape,len(data.shape)
            if len(data.shape) == 1:
                print data,type(data)
                data = np.array([data])
                #print "--.>",data,type(data)
            #print file,"data.shape:",data.shape,data.shape[0],data.shape[1]
            #print "importing "+file+" done!"

            shape = data.shape
            lines=shape[0]
            columns=shape[1]
            theList=range(columns)
            N=2
            groupby = [theList[n:n+N] for n in range(0, len(theList), N)]
            sum_all_curr = np.zeros((columns/2,lines),dtype=np.complex64)
            if idxi==0: #  and idxj==0:
                print "IDXI +++ 0000000000000000000000"
                sum_all_new = np.zeros((columns/2,lines),dtype=np.complex64)
                #print "shape:",shape
                #print "groupby:",groupby

            for idk,k in enumerate(groupby):
                #print idk,k
                #print data[:,k[0]],data[:,k[0]].shape
                #print data[:,k[1]],data[:,k[1]].shape
                #print sum_all_curr[:,idk].shape
                #print
                #print
                sum_all_curr[idk].real = data[:,k[0]]
                sum_all_curr[idk].imag = data[:,k[1]]
            if idxi==0: # and idxj==0:
                print "idxj:",idxj,"j:",j,"idxi:",idxi,"i:",i,file+" b aaa fft.shape:",sum_all_new.shape,"curr.shape:",sum_all_curr.shape
                sum_all_new = sum_all_curr
                print "idxj:",idxj,"j:",j,"idxi:",idxi,"i:",i,file+" a aaa fft.shape:",sum_all_new.shape,"curr.shape:",sum_all_curr.shape
            else:
                print "idxj:",idxj,"j:",j,"idxi:",idxi,"i:",i,file+" b bbb fft.shape:",sum_all_new.shape,"curr.shape:",sum_all_curr.shape
                sum_all_new = np.hstack((sum_all_new,sum_all_curr))
                print "idxj:",idxj,"j:",j,"idxi:",idxi,"i:",i,file+" a bbb fft.shape:",sum_all_new.shape,"curr.shape:",sum_all_curr.shape
                #print "saa:",sum_all_new
            print
        np.save("space_fft_"+j,sum_all_new)
    return


##########################################################################################
# functions to create space_fft from lammps outputfiles ##################################
##########################################################################################
def lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(f1,f2,rm=False):
    if f1 == f2:
        print "both files equal therefore return"
        return True
    c1 = os.popen('tail -1 '+f1+" | awk '{print $NF}'").read().rstrip().rstrip()
    c2 = os.popen('tail -1 '+f2+" | awk '{print $NF}'").read().rstrip().rstrip()
    c3 = os.popen('tail -1 '+f1+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c4 = os.popen('tail -1 '+f2+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    print "c1,c2:",c1,c2
    print "c3,c4:",c3,c4
    if c1 == c2 and c3==c4:
        print "c1 == c2 and c3 == c4"
        if type(rm) != bool:
            if rm == f1 or rm == f2:
                if rm == f1:
                    print "removing:",rm,"since same content in",f2
                if rm == f2:
                    print "removing:",rm,"since same content in",f1
                os.remove(rm)
        return True

    else:
        print f1,"and",f2,"are not equal!"
        return False

def lammps_check_if_trj_lammps_and_trj_lammpsnew_equal_also_check_head(f1,f2,rm=False):
    if f1 == f2:
        print "both files equal therefore return"
        return True
    c1 = os.popen('tail -1 '+f1+" | awk '{print $NF}'").read().rstrip().rstrip()
    c2 = os.popen('tail -1 '+f2+" | awk '{print $NF}'").read().rstrip().rstrip()
    c3 = os.popen('tail -1 '+f1+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c4 = os.popen('tail -1 '+f2+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()

    c5 = os.popen('head -1 '+f1+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c6 = os.popen('head -1 '+f2+" | awk '{print $(NF-1)}'").read().rstrip().rstrip()
    c7 = os.popen('head -1 '+f1+" | awk '{print $(NF)}'").read().rstrip().rstrip()
    c8 = os.popen('head -1 '+f2+" | awk '{print $(NF)}'").read().rstrip().rstrip()
    print "c1,c2:",c1,c2
    print "c3,c4:",c3,c4
    if c1 == c2 and c3==c4 and c5 == c6 and c7 == c8:
        print "c1 == c2 and c3 == c4 and c5 == c6 and c7 == c8"
        if type(rm) != bool:
            if rm == f1 or rm == f2:
                if rm == f1:
                    print "removing:",rm,"since same content in",f2
                if rm == f2:
                    print "removing:",rm,"since same content in",f1
                os.remove(rm)
        return True

    else:
        print f1,"and",f2,"are not equal!"
        return False

def lammps_split_lammps_trajectory_in_xaa_files(filenamelammps="trj_lammps.out",filenamepos="trj_lammpsnew.out",positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in",linesokformemory = 800000000,qpoints_all=False,args=False,verbose=False):
    ''' splits lammps trajectory file in xaa, xab, xac ...'''
    print bcolors.OKGREEN + "lammps_split_lammps_trajectory_in_xaa_files" +bcolors.ENDC
    print
    print
    print "das splitten und so weiter sollte alles auf einen knoten submittet werden!"
    print "daher: exit after grep und split!"
    print
    print
    positionsfile = glob.glob(positionsfilename)
    if len(positionsfile) != 1:
        print "positionsfile:",positionsfile
        sys.exit("Not one positionsfile but "+str(len(positionsfile)))

    print "positionsfile:",positionsfile

    infile = glob.glob(infilefilename)
    print "infile:",infile
    if len(infile) != 1:
        print "infile:",infile
        sys.exit("Not one infile but "+str(len(infile)))

    atoms1 = os.popen('wc -l '+positionsfile[0] + "| awk '{print $1}'").read()
    atoms = int(atoms1)-8
    atoms2 = int(os.popen('head -n 2 '+positionsfile[0] + "| tail -1 | awk '{print $1}'").read().rstrip())
    print "atoms2:",atoms2
    if atoms != atoms2:
        print "atoms:",atoms,"atoms2:",atoms2
        sys.exit("atoms not atoms2")

    supercelllength = float(os.popen('head -n 4 '+positionsfile[0] + "| tail -1 | awk '{print $2}'").read().rstrip())
    print "supercelllength:",supercelllength
    print "structure:",args.structure,"usestruct:",args.usestruct
    sc = N = round(float((atoms/args.usestruct)**(1/3.)),0)
    alat = supercelllength/N
    print "alat:",alat
    stepslammps = int(os.popen('grep "^run " '+infile[0] + "| awk '{print $2}'").read().rstrip())
    dump = int(os.popen('grep "^dump dump1" '+infile[0] + "| awk '{print $5}'").read().rstrip())

    print "atoms:",atoms,type(atoms)
    print "sc:",sc
    print "stepslammps:",stepslammps
    print "dump:",dump
    stepswritten = stepslammps/dump+1
    print "stepswritten:",stepswritten
    lineswritten = stepswritten*atoms
    print "lineswritten:",lineswritten,"in trj_lammpsnew.out"
    print "linesokformemory:",linesokformemory
    nsteps = linesokformemory / atoms
    nsteps = int(round(linesokformemory / atoms/10,0)*10)
    print "nsteps:",nsteps,"per xaa file"
    from string import ascii_lowercase
    stepsremain = stepswritten

    print "######################################################################################"
    print "# grep (now sed)"
    print "# it would be best to submit this to the que on one core (grep and sed)"
    print "# it is assumed thtat this was already done by the que, if not, comment in"
    print "######################################################################################"
    if os.path.isfile(filenamelammps):
        print filenamelammps,"exists!"
        if not os.path.isfile(filenamepos):
            print filenamepos,"does not exist! --> comment in the next lines!!"
            #print "here we have to ensure that we are in a screen session"
            #import socket
            #hostname=(socket.gethostname())
            #if hostname == 'cmmc001' or hostname == 'cmmc002':
            #    checkscreen = os.popen('echo $STY').read().rstrip()
            #    if checkscreen == "":
            #        onscreen = False
            #        print "onscreen:",onscreen
            #        sys.exit("ERROR: YOU WANT TO GREP A POSSIBLY LARGE FILE and should do this in a screen session!")
            #    else:
            #        onscreen = True
            #        print "onscreen:",onscreen,os.getcwd()
            #        print "10 GB take approx. 2min 20 sec ..."
            #        print "100 GB take approx. 1400 sec = 23 min ..."
            #        #os.popen('grep "^1 " '+filenamelammps + "| awk '{print $2,$3,$4}' > trj_lammpsnew.out").read()
            #        print "hier nun der submit fuer sed und split welcher dann wieder lammps_pos_to_sum startet und am besten gleich noch die c++ jobs startet"
            #        print "sed und split brauchen nur die nsteps, das ist einfach"
            #        print "hier sollte nur ein ping gegeben werden wenn sed/split/ gebraucht wird da sonst alles nicht da"
            #        os.popen("sed -n 's/^1 //p' "+filenamelammps + " > trj_lammpsnew.out").read()
            #        if os.path.isfile(filenamelammps) and os.path.isfile(filenamepos):
            #            lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(filenamelammps,filenamepos,rm=filenamelammps)
    else:
        print filenamelammps,"does not exist anymore, presumably grep already done!"

    if os.path.isfile(filenamelammps) and os.path.isfile(filenamepos):
        lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(filenamelammps,filenamepos,rm=filenamelammps)

    print
    print
    splitanzahl=-1
    print "      stepstot","\t","steps  ","\t","stepsremain"
    print "---------------------------------------------------"
    lastfile = False
    lastamoutoflines = False
    xaa_filenames = []
    xaa_steps = []
    xaa_lines = []

    # make here an array
    idlast=0
    for idx,c in enumerate(ascii_lowercase):
        for idy,d in enumerate(ascii_lowercase):
            nstepscurr = nsteps
            #print "(stepsremain - nsteps):",(stepsremain - nsteps)
            if (stepsremain - nsteps) > 0:
                stepsremain = stepsremain - nsteps
                exittrue = False
            else:
                nstepscurr = stepsremain
                stepsremain = 0
                exittrue = True
            splitanzahl = splitanzahl+1
            lines = nstepscurr*atoms
            print idlast,"x"+c+d,stepswritten,"\t",nstepscurr,"\t",stepsremain,splitanzahl,lines,exittrue
            xaa_filenames.append("x"+c+d)
            xaa_steps.append(nstepscurr)
            xaa_lines.append(lines)
            lastfile = "x"+c+d
            lastamoutoflines = lines
            idlast = idlast + 1
            if exittrue: break
        if exittrue: break

    #for i in np.arange(20)+1:
    #    atoms = 4*i**3
    #    print  i,atoms,linesokformemory / atoms, int(round(linesokformemory / atoms/10,0)*10)

    if splitanzahl == 0:
        if os.path.isfile(filenamepos) and not os.path.isfile('xaa'):
            os.rename(filenamepos, "xaa")
        if os.path.isfile(filenamepos) and os.path.isfile('xaa'):
            print "move one to other if those are equal"
            lammps_check_if_trj_lammps_and_trj_lammpsnew_equal_also_check_head(filenamepos,'xaa',rm=filenamepos)
        if not os.path.isfile(filenamepos) and os.path.isfile('xaa'):
            pass # this is also ok
        if not os.path.isfile(filenamepos) and not os.path.isfile('xaa'):
            sys.exit("where are the positions?")
        xaa_filenames = [ filenamepos ]
        xaa_filenames = [ 'xaa' ]
        xaa_steps = [ stepswritten ]
        xaa_lines = [ lineswritten ]

    print "splitanzahl:",splitanzahl
    print "xaa_filenames:",xaa_filenames
    print "xaa_steps:",xaa_steps
    print "xaa_lines:",xaa_lines
    print
    print
    print
    print "######################################################################################"
    print "# split"
    print "######################################################################################"
    if splitanzahl > 0:
        print "splitanzhal > 0, in fact it is:",splitanzahl
        if os.path.isfile("xaa"):
            print "xaa file already exists, splitting presumable already started/done!"
        else:
            print "splitting by:",str(int(nsteps*atoms)),"since xaa does not exist! --> splitting"
            print "hmmm ...unexpected... this should have been done by the cluster"

        #    checkscreen = os.popen('echo $STY').read().rstrip()
        #    if checkscreen == "":
        #        onscreen = False
        #        print "onscreen:",onscreen
        #        #sys.exit("ERROR: YOU WANT TO SPLIT A POSSIBLY LARGE FILE and should do this in a screen session!")
        #    else:
        #        onscreen = True
        #        print "onscreen:",onscreen
        #    os.popen('split -l '+str(int(nsteps*atoms))+" "+filenamepos).read()

    print "lastfile:",lastfile
    print "lastamoutoflines:",lastamoutoflines
    if splitanzahl > 0 and os.path.isfile("xaa") and os.path.isfile(lastfile) and os.path.isfile(filenamepos):
        lammps_check_if_trj_lammps_and_trj_lammpsnew_equal(filenamepos,lastfile,rm=filenamepos)

    print "######################################################################################"
    print "now check if last file .."+xaa_filenames[-1]+".. has the correct number of lines:"
    print "######################################################################################"

    if os.path.isfile(xaa_filenames[-1]+"_wcl"):
        print xaa_filenames[-1]+"_wcl","exists, therefore correct number of lines!"

    if not os.path.isfile(xaa_filenames[-1]+"_wcl") and os.path.isfile(xaa_filenames[-1]):
        print "it seems that,",xaa_filenames[-1]+"_wcl","does not exist and ",xaa_filenames[-1],"exists!"
        wclast = int(os.popen('wc -l '+xaa_filenames[-1] + "| awk '{print $1}'").read().rstrip())
        print "wclast:",wclast,type(wclast)
        print "xaa_lines[-1]:",xaa_lines[-1],type(xaa_lines[-1])
        if wclast != xaa_lines[-1]:
            target = open(xaa_filenames[-1]+"_wclerror", 'w')
            target.write("wc -l gives "+str(wclast)+" lines\n")
            target.close()
            sys.exit("wclast not last xaa_lines")
        else:
            print "writing that number of lines is ok",xaa_filenames[-1]+"_wcl"
            target = open(xaa_filenames[-1]+"_wcl", 'w')
            target.write("wc -l gives "+str(wclast)+" lines\n")
            target.close()
            print "######################################################################################"
            print "seems to have correct number of lines!"
            print "######################################################################################"

    if not os.path.isfile(xaa_filenames[-1]+"_wcl") and not os.path.isfile(xaa_filenames[-1]):
        print "######################################################################################"
        print "the last file, ",lastfile," does not exist! I will therefore exit here"
        print "######################################################################################"
        sys.exit()

    print "######################################################################################"
    print "now check if xaa_ folder exist, if not, create them and send to cluster or not"
    print "######################################################################################"
    xaa_folder= glob.glob("xa?_")
    print "len(xaa_folder):",len(xaa_folder)
    if len(xaa_folder) == 0 and check_for_all_space_fft_and_rm_xaa_files_folders(sc,args) == False:
        print "##################################################################################"
        print "create xaa_ folders for c++ skript (even one folder should be submitted to que)"
        print "##################################################################################"
        print "all xaa_filenames:",xaa_filenames
        avail_xaa_filenames = []
        avail_xaa_steps = []
        for idxx,op in enumerate(xaa_filenames):
            #print "xaa_file:",op
            if os.path.isfile(op) == True:
                avail_xaa_filenames.append(xaa_filenames[idxx])
                avail_xaa_steps.append(xaa_steps[idxx])
        print "avail_xaa_filenames:",avail_xaa_filenames,len(avail_xaa_filenames)
        if len(avail_xaa_filenames) > 0:
            #print "I am here!!!"
            #sys.exit()
            submit_c_skript_parallel(avail_xaa_filenames,avail_xaa_steps,sc,alat,atoms,executable="/cmmc/u/aglen/sascha/source/test",args=args)
        else:
            print "##################################################################################"
            print "no xa files available!"
            print "##################################################################################"
    else:
        print "xaa folder exist already."
    return xaa_filenames,xaa_steps

def submit_c_skript_parallel(xaa_filenames,xaa_steps,sc,alat,atoms,executable="/cmmc/u/aglen/sascha/source/test",qsub="/u/aglen/submit.lammps.lifetimejob.c_sascha.cmfe40.sh",args=False):
    print
    print "##############################################################################"
    print "now applying C++ skript (will do this even for one xaa folder)"
    print "this should also create a joblist to submit s++ skript to cluster or directly"
    print "do so if necessary"
    print "this should be started in parallel in any case! -> submit to cluster"
    print "##############################################################################"
    qpoints_all = get_all_qpoints(args.qvec,args)

    hier = os.getcwd()
    for idx,filename in enumerate(xaa_filenames):
        print idx,filename
        os.chdir(hier)
        lammps_create_c_SETTINGS_file(infile=filename,atoms=atoms,steps=xaa_steps[idx],scale=alat*sc,qpoints_all=qpoints_all,executable=executable)
        os.chdir(filename+"_")
        print "cwd:",os.getcwd()
        print "here execture qsub"
        print os.popen('qsub '+qsub).read().rstrip()   # this print statement submits the job
    return

def lammps_create_c_SETTINGS_file(filename="SETTINGS",infile=False,atoms=False,steps=False,scale=False,qpoints_all=False,executable=False):
    ''' infile: xaa'''
    if type(infile) == bool:
        sys.exit("specify infile")
    if type(atoms) == bool:
        sys.exit("specify atoms")
    if type(scale) == bool:
        sys.exit("specify scale")
    if type(infile) == bool:
        sys.exit("specify infile")
    if type(qpoints_all) == bool:
        sys.exit("specify qpoints_all")
    if type(executable) == bool:
        sys.exit("specify executable")

    if not os.path.isfile(infile): sys.exit("file "+infile+" does not exist!")
    if not os.path.isfile(executable): sys.exit("file "+executable+" does not exist!")

    if not os.path.exists(infile+"_"): os.makedirs(infile+"_")

    settingsfile = infile+"_"+"/SETTINGS"
    if os.path.isfile(settingsfile):
        os.remove(settingsfile)

    target = open(settingsfile, 'w')
    target.write("NATOMS = "+str(atoms)+"\n")
    target.write("NSTEPS = "+str(steps)+"\n")
    for idx,i in enumerate(qpoints_all):
        target.write("Q"+str(idx+1)+" = "+str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
    target.write("SCALE = "+str(scale)+"\n")
    target.write("INFILE = ../"+infile+"\n")
    target.close()
    import shutil
    shutil.copy2(executable, infile+"_/test")

#############
# UNUSED
#############
def get_space_fft_from_positions_fast_read_at_once_UNUSED(filename, qpoints_all = False, scale_data = 1., N=False, usestruct=False,args=args):
    ''' dont use engine='python' in read csv which is 6 times slower,
        - this module is aimed for speed, therefore loads the positions as whole into memory
        - chunking however should also work with c engine
        - when checked, cmmc001 had 56GB of free ram
        - 10x10x10sc (4000atoms) with 10^6steps (every 40 wirtten) did not fit into memory
          also it should, I guess it was close ;)

        then, go qpoint by qpoint and do the calculation(s)
        when
    '''
    print "#################################### try_lammpsnew.out #########################"
    #filename = 'trj_lammpsnew.outfirst3'
    #filename = 'trj_lammpsnew.outfirst300'     # 0.86 vs 2.2 sec
    #filename = 'trj_lammpsnew.outfirst3000'     # 3.7 vs 22.0 sec (head -12000000) = 12*10^6 lines of 10*10^8
                                                # size in memory: 0.288 GB
    #filename = 'trj_lammpsnew.outfirst30000'   # 36.7sec; size in memory: 2.88GB (1/10 of whole file)
    #filename = 'trj_lammpsnew.out'             # lines: 1000004000; 23GB; so should be 12GB in memory;
    filename = 'trj_lammpsnew.out'

    import time
    start1 = time.time()
    start = time.time()
    #df = read csv(filename, sep=' ', header=None,engine='c',) # 3.7 sec
    df = read_csv(filename, sep=' ', header=None,engine='c',chunksize=4000*2) # 3.7 sec
    for chunk,pos in enumerate(df):
        print "chunk:",chunk
        xx = pos/float(scale_data)
        print "pos  :",xx
        print xx.as_matrix()
        sys.exit()
    #df = read csv(filename, sep=' ', header=None,engine='c', dtype={'0': np.float32, '1': np.float32, '2': np.float32}) # 3.7 sec
    print sys.getsizeof(df)*1e-09,"GB", "df == pandas read in"#,df.nbytes
    print df[0].dtype # float64
    end = time.time()
    print "#1) read in data :",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC

    start = time.time()
    df = df/float(scale_data)
    print sys.getsizeof(df)*1e-09,"GB", "df == pandas"
    df = df.as_matrix()
    print sys.getsizeof(df)*1e-09,"GB", "df == nupy",df.nbytes*1e-09,"in GB"
    z = np.copy(df)
    print sys.getsizeof(z)*1e-09,"GB", "df == numpy copy",z.nbytes*1e-09,"in GB"
    #z = z.astype(complex64)
    z = np.complex64(z)
    print sys.getsizeof(z)*1e-09,"GB", "df == numpy copy",z.nbytes*1e-09,"in GB"



    qs_all = equivalent_qs([1,1,1])
    qs_all = equivalent_qs([1,0,0])
    print "00000000000:",df.shape
    xxx = df*qs_all[0]
    print sys.getsizeof(xxx)*1e-09,"GB"
    dfqs = np.exp(np.sum(df*qs_all[0],axis=1)).reshape(-1, N**3*args.usestruct).sum(axis=1) # the reshape is quick, it is probably the np.exp which takes time
    print sys.getsizeof(dfqs)*1e-09,"GB"
    print "------>",dfqs.dtype
    #aaa[4000:8000].sum()
    #aaa[8000:12000].sum()

    end = time.time()
    print "#2) np.exp(np.sum()):",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    start = time.time()
    #dfqs = dfqs.reshape(-1, N**3*args.usestruct).sum(axis=1)
    end = time.time()
    print "#3) reshape.sum:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    #aaa = np.sum(np.exp(np.sum(xxx*qs_all[0],axis=1)))
    print
    print dfqs[:10]
    print dfqs.shape
    end1 = time.time()
    print "#5) total:",bcolors.HEADER+  str((end1 - start1)) +" sec"+ bcolors.ENDC
    return dfqs

# needs read csv
def get_space_fft_from_hd5f_UNUSED(filename='dump_h5md.h5',qpoint=[1,0,0],scale_data = False):
    print "#################################### hd5f #########################"
    import time
    start1 = time.time()

    filename = 'dump_h5md.h5'
    start = time.time()
    f = h5py.File(filename, 'r')
    positions = f['particles/all/position']
    end = time.time()
    print "#1) hdf5 read:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    start = time.time()
    print positions['value'].shape
    positions = positions['value'][:]/float(scale_data)
    print "----->",positions.shape
    end = time.time()
    print "#2) tonumpy:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    print "positions.shape:",positions.shape
    print sys.getsizeof(positions)*1e-09,"GB", "df == numpy copy",positions.nbytes*1e-09,"in GB"


    start = time.time()
    out = np.zeros(positions.shape[0])
    end = time.time()
    print "#3) init empty out:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    print "sh:",out.shape
    qs_all = equivalent_qs(qpoint)
    print "qs_all[]:",qs_all[0]
    start = time.time()
    for step,pos in enumerate(positions):
        #xxx = positions['value'][step]/float(scale_data)
        #out[step] = np.sum(np.exp(np.sum(xxx*qs_all[0],axis=1)))
        out[step] = np.sum(np.exp(np.sum(pos*qs_all[0],axis=1)))
    end = time.time()
    print "#4) exp:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    print out.shape
    print out[:10]
    end1 = time.time()
    print "#5) total:",bcolors.HEADER+  str((end1 - start1)) +" sec"+ bcolors.ENDC



    return


##########################################################################################
# functions related to creation of the powerspectrum #####################################
##########################################################################################

def get_lifetimes_for_different_timesteps_parallel(mdstepstocalc,qpoint,space_fft,idxprint=998,idxprintmax=998,args=False):
    ''' text '''
    qpointstr = qpointstring(qpoint)
    finaldata_all_good = np.zeros((len(mdstepstocalc),6))
    finaldata_all_min = np.zeros((len(mdstepstocalc),6))
    filenameout_all = ["ps"+qpointstring(qpoint)+"_"+str(mdstepstocalc[x])+".dat" for x in range(len(mdstepstocalc))]


    jobsa = []
    serial = False
    #print printblue("serial: "+str(serial))

    #if serial:
    #print "args.seriell>>:",args.seriell
    if args.seriell: # seriell !!!!!!!!!!!!!!!!!!!
        print "seriell !!!!!!!!"
        #sys.exit()
        start = time.time()
        for idx,mdsteps in enumerate(mdstepstocalc):
            #print "filenameout:",filenameout
                finaldata_good,finaldata_min,filenameout = get_lifetime_for_particular_mdsteps(mdsteps,mdstepstocalc,space_fft,qpoint,args=args)
                finaldata_all_good[idx] = finaldata_good
                finaldata_all_min[idx] = finaldata_min
                filenameout_all[idx] = filenameout
        end = time.time()
        print printblue("SERIAL GET LIFETIME OF MDSTEPS: "+str((end-start))+" sec.")


    if not args.seriell: # parallel  !!!!!!!!!!!!!!!!!
        #print "parallel !!!!!!! in mdstepstocalc"
        #sys.exit()

        def currentfunc_helper(args):
            return get_lifetime_for_particular_mdsteps(*args)

        # print printblue('not serial!')
        start = time.time()
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        jobs = []
        for idx,i in enumerate(mdstepstocalc):
            arguments = (i,mdstepstocalc,space_fft,qpoint,args,return_dict)
            p = multiprocessing.Process(target=currentfunc_helper, args=[arguments])
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
        #print "#################33"
        #print return_dict.values()
        #print "##99"
        #print return_dict.keys()
        #print "##88"
        #print "##y8"

        d = return_dict
        #for k, v in d.items():
        for k, v in enumerate(d.items()):
            #print "k:",k,'v:',v, "v[0]:",v[0],"v[1]:",v[1][0]
            finaldata_all_good[k]=v[1][0] # good
            finaldata_all_min[k]=v[1][1] # good
            #filenameout_all[k]=v[1][2] # good
        #for k, v in d.items():
        #    print k, v
        #
        finaldata_all_good = finaldata_all_good[np.argsort(finaldata_all_good[:, 5])]
        finaldata_all_min  = finaldata_all_min[np.argsort(finaldata_all_min[:, 5])]
        end = time.time()
        #print printblue("PARALLEL GET LIFETIME OF MDSTEPS: "+str((end-start))+" sec.")


    if args.verbose:
        print "RESULTING: finaldata_all_good:",type(finaldata_all_good)
        print finaldata_all_good
        print "RESULTING: finaldata_all_min :"
        print finaldata_all_min
        print "finaldata_all_good[:,0] # sigma"
        print "finaldata_all_good[:,1] # freq"
        print "finaldata_all_good[:,2] # linewidth ------^"
        print "finaldata_all_good[:,3] # linewidth min -----------------^"
        print "finaldata_all_good[:,4] # linewidth max -----------------------------------^"
        #print "filenameout_all:"
        #print filenameout_all
    return finaldata_all_good, finaldata_all_min, filenameout_all

def space_fft_to_powerspectrum(        qpoint,       execute_timeinversion = True, idxprint =0, idxprintmax = 0, args = False):
    ''' calculates the power_spectrum for a certain qpoint (from space_fft)
    - space_fft contins for a single q point all (symmetry) equivalent q points
    - output of this skript is the power_spectrum for the corresponding q-point
    - qpoint is just for saving the power_spectrum
    '''
    qpointstr = qpointstring(qpoint)
    appendlast = True
    if args.verbose:
        print "qpoint:",qpoint
        print "execute_timeinversion:",execute_timeinversion
        print "idxprint:",idxprint
        print "idxprintmax:",idxprintmax
    dt = args.dt
    mdstepstocalcin = args.calclast
    check_qpoint_for_crossing = check_qpoints_for_crossing(qpoint, args.supercell, args.structure)
    if args.verbose:
        print "mdstepstocalin:",mdstepstocalcin
        print "check_qpoint_for_crossing:",check_qpoint_for_crossing
        print "dt:",dt
        print "----------------------------------- all args ----------------------------"
        print "args 99:",args
        print "----------------------------------- all args ----------------------------"
    loadfile = "space_fft_"+qpointstring(qpoint)+".npy"
    print printblue("np.load("+loadfile+"); executre Timeinversion is "+str(execute_timeinversion)+"; args.seriell is "+str(args.seriell)+"; crossing "+str(check_qpoint_for_crossing)+";")
    if os.path.isfile(loadfile) != True:
        sys.exit("ERROR 44: :"+loadfile+" does not exist! use either -fftpy option (slow but implemented) or use the -fftc optein (fast but currently only working with lammps jobs) to get it!")
    space_fft = np.load(loadfile)
    #print "np.load(space_fft) done!",type(space_fft)
    if args.verbose: print "execute_timeinversion:",execute_timeinversion
    if execute_timeinversion != True:
        print_warning("execute_timeinversion is False!!!!")

    if dt == False:
        sys.exit("ERROR: need dt!")

    def which_mdsteps_to_calc(space_fft, mdstepstocalcin = False):
        ''' just determines which mdsteps you want to run over
        '''

        #####################################################
        # MDSTEPSTOCALC Manual
        #####################################################
        mdstepstocalc_all = np.array([1000,5000,9000,10000,20000,50000,90000,100000,120000,140000,160000,180000,200000,500000,\
                900000,1000000,120000,124000,126000,128000,2000000,5000000,9000000,10000000,20000000,50000000,90000000,100000000])
        #mdstepstocalc_all = np.array([10000])
        #mdstepstocalc_all = np.array([200000])
        #mdstepstocalc_all = np.array([1000])
        #mdstepstocalc_all = np.array([5000])
        #mdstepstocalc_all = np.array([200])
        #mdstepstocalc_all = np.array([100,200])


        #####################################################
        # MDSTEPSTOCALC last order of magnitude to check for error
        #####################################################
        mdstepstocalc_all = np.arange(11)*0
        #print mdstepstocalc_all,type(mdstepstocalc_all)
        steps_tocheck_md_error = int(space_fft.shape[1]/10)
        #print steps_tocheck_md_error,type(steps_tocheck_md_error)
        #print "steps_tocheck_md_error",steps_tocheck_md_error,space_fft.shape[1]
        for idx,i in enumerate((np.arange(10)+1)[::-1]):
            add = int(steps_tocheck_md_error*10 - i*steps_tocheck_md_error)
            #print "add:",add,type(add)
            mdstepstocalc_all[i]=int(add)
            #np.append(mdstepstocalc_all,add)
        mdstepstocalc_all = np.trim_zeros(mdstepstocalc_all)
        #mdstepstocalc_all = np.array([5000])
        #mdstepstocalc_all = np.array([10000])
        #print mdstepstocalc_all.shape,mdstepstocalc_all
        #print lifetimestocheck.shape,lifetimestocheck
        #print "kk:",mdstepstocalc_all

        if appendlast:
            # das /10*10 macht aus 20000{1,9} 200000  aber aus 200010 macht es 200010
            mdstepstocalc_all = np.sort(np.append(mdstepstocalc_all,[space_fft.shape[1]/10*10]))

        #####################################################
        # MDSTEPSTOCALC schoenheitskorrekturen der zahlen
        #####################################################
        #print "mdstepstocalc:",mdstepstocalc
        mdstepstocalc_all = np.sort(np.unique(mdstepstocalc_all[np.where(mdstepstocalc_all <= space_fft.shape[1])[0]]))
        #print "mdstepstocalc_all:",mdstepstocalc_all
        #print "mdstepstocalc:",mdstepstocalc,type(mdstepstocalc)
        #if type(mdstepstocalc) == bool:
        mdstepstocalc = mdstepstocalc_all
        #print "mdstepstocalc:",mdstepstocalc,type(mdstepstocalc)
        ##sys.exit()

        ##print "mdstepstocalc_________________:",mdstepstocalc,type(mdstepstocalc)
        #if type(mdstepstocalc) != bool:
        #    if type(mdstepstocalc) == list:
        #        mdstepstocalc = np.array(mdstepstocalc)
        #        #print "hier 1"
        #    elif type(mdstepstocalc) == int:
        #        # check if we want more in mdstepstocalc than we can offer in space_fft:
        #        if mdstepstocalc > mdstepstocalc_all.max():
        #            mdstepstocalc = np.array([mdstepstocalc_all.max()])
        #        else:
        #            mdstepstocalc = np.array([mdstepstocalc])
        #        #print "hier 2"
        #    elif type(mdstepstocalc) == np.ndarray:
        #        #print "hier 3"
        #        pass
        #else:
        #    sys.exit("error: unknow type for mdstepstocalc")
        #    #print "mdstepstocalc_all:",mdstepstocalc_all

        #print "2 mdstepstocalc_all:",mdstepstocalc_all
        mdstepsmax=mdstepstocalc_all[-1]
        #print "2 mdstepstocalc:",mdstepstocalc
        #print mdstepstocalc.min(),-1.*(int(math.log10(mdstepstocalc.min()))-1)
        ordmag = int(math.log10(mdstepstocalc.min()))
        #print mdstepstocalc.min(),10**ordmag
        factchange=10**ordmag
        factchange=10**ordmag

        mdstepstocalc = np.unique(mdstepstocalc/factchange*factchange) # remove elements wich are too close
        #print "1:",mdstepstocalc
        mdstepstocalc[-1]=mdstepsmax
        mdstepstocalc = mdstepstocalc[np.where(mdstepstocalc > 999)[0]]
        mdstepstocalc_all = mdstepstocalc_all[np.where(mdstepstocalc_all > 999)[0]]
        if mdstepstocalc_all[0] >= 20000:
            mdstepstocalc_all = np.insert(mdstepstocalc_all,0,10000)
        if mdstepstocalc_all[0] >= 10000:
            mdstepstocalc_all = np.insert(mdstepstocalc_all,0,5000)
        if mdstepstocalc_all[0] >= 5000:
            mdstepstocalc_all = np.insert(mdstepstocalc_all,0,1000)
        #print "--> mdstepstocalc:", mdstepstocalc# ,mdstepstocalc[np.where(mdstepstocalc > 999)[0]]
        #print "type:",type(mdstepstocalcin)
        if type(mdstepstocalcin) == bool:
            if mdstepstocalcin == True:
                mdstepstocalc = np.array([mdstepstocalc[-1]])
        elif type(mdstepstocalcin) != bool:
            if mdstepstocalcin == 'last':
                mdstepstocalc = np.array([mdstepstocalc[-1]])
            elif type(mdstepstocalcin) == int:
                mdstepstocalc = np.array([mdstepstocalcin])
            else:
                sys.exit("ERROR 777 check here")
        #print "--> mdstepstocalc:", mdstepstocalc# ,mdstepstocalc[np.where(mdstepstocalc > 999)[0]]


        return mdstepstocalc

    mdstepstocalc = which_mdsteps_to_calc(space_fft=space_fft,mdstepstocalcin = mdstepstocalcin)


    if args.verbose:
        print "-----> mdstepstocalc",mdstepstocalc
    if type(args.mdsteps) != bool:
        mdstepstocalc = np.array(args.mdsteps)
    if args.verbose:
        print "-----> mdstepstocalc",mdstepstocalc

    ##################################################################################
    ##################################################################################
    # get error due to MD convergenz (this takes most of the time)
    ##################################################################################
    ##################################################################################
    #####################################################
    # schleife ueber mdsteps
    # it would be best to start with the longest, write this, and then move to the shorter runs
    #####################################################
    finaldata_good_all, finaldata_min_all, filenameout_all = get_lifetimes_for_different_timesteps_parallel(mdstepstocalc,qpoint,space_fft,idxprint=998,idxprintmax=998,args=args)


    finaldata_good, finaldata_min, filenameout = finaldata_good_all[-1], finaldata_min_all[-1], filenameout_all[-1]
    sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin = finaldata_min
    sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood = finaldata_good


    ##################################################################################
    ##################################################################################
    # get error due to MD convergenz DONE
    ##################################################################################
    ##################################################################################



    #print "filenameout:",filenameout
    #############################################
    # get error due to MDlength
    # mdstepstocald have to be in increasing order
    #############################################
    def get_lt_md_error(lifetimestocheck,mdstepstocalc):
        idxmax = np.where(lifetimestocheck == lifetimestocheck.max())[0].max()

        lifetimestocheckerrormd = np.zeros(len(mdstepstocalc))
        for ind,yy in enumerate(lifetimestocheck):
            #print ind,y,"   ",lifetimestocheck[:ind+1][-6:]
            considerlt = lifetimestocheck[:ind+1][-5:]
            lifetimestocheckerrormd[ind] = considerlt.max() - considerlt.min()

        idxmax = np.where(lifetimestocheckerrormd == lifetimestocheckerrormd.max())[0].max()
        for ind,yy in enumerate(lifetimestocheckerrormd):
            if ind <= idxmax:
                lifetimestocheckerrormd[ind] = lifetimestocheckerrormd[idxmax]
        return lifetimestocheckerrormd

    lttake = ltgood
    lttake = ltmin
    error_md_lt_all = get_lt_md_error(finaldata_good_all[:,2],finaldata_good_all[:,1])
    error_md_lt_all_min = get_lt_md_error(finaldata_min_all[:,2],finaldata_min_all[:,1])
    #print "error_md_lt_all       :",error_md_lt_all
    error_md_lt = np.array([error_md_lt_all[-1],error_md_lt_all_min[-1]]).min()
    error_md_lt = error_md_lt_all[-1]

    ########################################
    # error lifetimes
    ########################################
    def print_error_lt(error_md_lt,error_lt_up,error_pos,ltgood,checkmaxmin,k="UP",mg="good"):
        print
        print "------------------------------"
        print "-->",mg
        print "error_md_lt:",error_md_lt
        print "error_lt_up:",error_lt_up," == ",ltgood," - ",checkmaxmin,"(diff between MIN and GOOD)"
        print "------------------------------"
        print "SUM ERROR LT "+k+":",error_pos
        print "------------------------------"
        print "------------------------------"
        print
        print

    checkmax = np.array([ltgoodmax,ltminmax]).max()
    error_lt_up = abs(lttake - checkmax)
    error_pos =  error_lt_up + error_md_lt
    if args.verbose:
        print_error_lt(error_md_lt,error_lt_up,error_pos,lttake,checkmax,k="UP",mg="FROM GOOD")

    checkmin = np.array([ltmin,ltgoodmin,ltminmin]).min()
    error_lt_down =  abs(lttake- checkmin)
    error_neg =  error_lt_down+ error_md_lt
    if error_neg > lttake: error_neg = lttake  # dont let the error go under 0

    if args.verbose:
        print_error_lt(error_md_lt,error_lt_down,error_neg,lttake,checkmin,k="DOWN",mg="FROM MIN")


    if len(mdstepstocalc) <= 1:
        print "now not doing error due to MD convergence"
        return
    ########################################
    # error due to md
    ########################################
    lifetimestocheck = finaldata_good_all[:,2]
    freqstocheck = finaldata_good_all[:,1]
    error_md_lt = abs(lifetimestocheck[:-6].max() - lifetimestocheck[:-6].min())
    error_md_fq = abs(freqstocheck.max() - freqstocheck.min())

    error_fq = abs(fqgood -fqmin) + error_md_fq
    #print error_gaussian_lt
    #print aerror_gaussian_fq
    #print "mdsteps:",mdsteps

    #############################################
    # get error due to MDlength
    #       - ps*lifetimes_md_convergence.dat
    #       - ps*lifetimes_md_convergencepsec.dat
    #############################################

    if args.verbose:
        print "mdstepstocalc         :",mdstepstocalc
        print "lifetimestocheck      :",[ round(elem, 2) for elem in lifetimestocheck]
        print "error_md_lt_all       :",error_md_lt_all
        print "error_md_lt_all_min   :",error_md_lt_all_min
    print printblue("WRITING5: "+filenameout+".lifetimes_md_convergencepsec.dat")
    #np.savetxt(filenameout+".lifetimes_md_convergence.dat",np.array(np.transpose([mdstepstocalc, lifetimestocheck, error_md_lt_all,error_md_lt_all])),fmt='%.5f')
    ftt = factor_mdsteps_to_time = 10**-15*dt/(10**-12)  # in ps
    np.savetxt(filenameout+".lifetimes_md_convergencepsec.dat",np.array(np.transpose([mdstepstocalc*ftt, lifetimestocheck, error_md_lt_all,error_md_lt_all])),fmt='%.5f')
    np.savetxt(filenameout+".freqs_md_convergencepsec.dat",np.array(np.transpose([mdstepstocalc*ftt, finaldata_good_all[:,1] ])),fmt='%.5f')

    #############################################
    # get final error for lifetimesgood
    # this will only work if loop is done from lowest to highest md step
    #############################################
    #if mdsteps == mdstepstocalc.max():
    if len(mdstepstocalc) >= 1: # it is assumed that sigmagood is for the last step
        filenameoutsmooth = filenameout+".smooth_good"+str(sigmagood)+".dat"

        outlifetimes = np.array([[mdstepstocalc.max(), lttake,error_pos,error_neg,sigmagood,999999,9,9,error_md_lt]])
        outfreq      = np.array([[mdstepstocalc.max(), fqgood,error_fq ,error_fq ,sigmagood,999999,9,9,error_md_fq]])
        #print 'outlifetimes:',outlifetimes
        print printblue("WRITING7: "+filenameout+".lifetimesgood.dat")
        if args.verbose:
            print "                 mdsteps             linewidth       error_pos         error_neg"
            print "outlifetimes:",outlifetimes
        np.savetxt(filenameout+".lifetimesgood.dat",outlifetimes,fmt='%.5f')
        np.savetxt(filenameout+".freqsgood.dat",outfreq,fmt='%.5f')



    return mdstepstocalc, lifetimestocheck

def space_fft_to_powerspectrum_ser_par(qpoints_all, execute_timeinversion = True, args = False):
    if not args.seriell:
        num_cores = multiprocessing.cpu_count()
        print "num_cores 888:",num_cores
        jobs = []
        for idx,i in enumerate(qpoints_all):
            def space_fft_to_powerspectrum_helper(args):
                space_fft_to_powerspectrum(*args)
            arguments = (i,execute_timeinversion,idx,len(qpoints_all),args)
            p = multiprocessing.Process(target=space_fft_to_powerspectrum_helper, args=[arguments])
            jobs.append(p)
            p.start()

        for job in jobs: # just to wait until job is done
            job.join()
    #seriell = True
    if args.seriell:
        for idx,i in enumerate(qpoints_all):
            space_fft_to_powerspectrum(i,execute_timeinversion,idx,len(qpoints_all),args)

    return

def get_lifetime_for_particular_mdsteps(mdsteps,mdstepstocalc,space_fft,qpoint,args=False,return_dict=False):
    qpointstr = qpointstring(qpoint)
    power_spectrum=np.zeros(mdsteps)
    power_spectrumb=np.empty(mdsteps,dtype=np.complex64)
    ''' at this stage the space_fft is already loaded '''

    ##########################################################################
    # calculate power spectrum (needs to stay all the time with both peaks until end for fitting)
    ##########################################################################
    # needs: space_fft, mdsteps, mdstepstocalc
    #print "WWWWWWWWWW this is started for every mdstep"
    #print "kkk",np.arange(len(space_fft))
    #for i in np.arange(len(space_fft)): # laeuft ueber alle symmetrieequivalten qs (6 stueck)
    #    print i
    #sys.exit()
    serial = True
    if serial:  # serial
        #print "start:"
        #print printblue('fftserial !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        #startp = time.time()
        for i in np.arange(len(space_fft)): # laeuft ueber alle symmetrieequivalten qs (6 stueck oder 3 stucek i == [ 0, 1, 2])
            #startt = time.time()
            ### MAIN
            #@ amplitude = (space_fft[i][:mdsteps])[:mdsteps]
            #@ #print 'in ampl.shape:',ampl.shape,ampl  # real and imag
            #@ ##############################################################################
            #@ ######### old definition MAIN (WORKING)
            #@ ##############################################################################
            #@ fourier_transformierte = fft(amplitude)
            #@ #print "old ft.shape:",ft.shape,ft       # real and imag
            #@ a =np.abs(fourier_transformierte)**2.                       # only real

            a =np.abs(fft((space_fft[i][:mdsteps])[:mdsteps]))**2.  # put into one command to save (hopefully) memory
            # timeinversion
            #aa = np.roll(a[::-1],1)
            #a += aa
            ##############################################################################
            ######### new definition MAIN
            ##############################################################################
            #@ampl = (space_fft[i][:mdsteps])[:mdsteps]
            #@amplka=np.empty(len(ampl),dtype=np.complex64)
            #@amplka.real=ampl.real
            #@amplka.imag=-ampl.imag

            #@intensity = np.abs(ampl)**2.            # only real
            #@np.savetxt("i"+str(i),intensity)
            #@intensity = ampl*amplka
            #@#print "new intensity.shape:",intensity.shape,intensity
            #@b = fft(intensity)                      # real and imag
            #@bb = fft(b)                      # real and imag
            #@np.savetxt("bb.real"+str(i),bb.real)
            #@power_spectrumb+=b

            #@#print "a.shape:",a.shape,a
            #@#print "b.shape:",b.shape,b
            #@#print

            # WRONG a =np.abs(fft(((space_fft[i][:mdsteps])[:mdsteps])**2))  # put into one command to save (hopefully) memory
            # DOES NOT WORK a =fft(np.abs((space_fft[i][:mdsteps])[:mdsteps])**2)  # put into one command to save (hopefully) memory
            #b =fft(np.abs(((space_fft[i][:mdsteps])[:mdsteps]))**2.)
            #print "n1a:",a.shape,type(a),a[:2]
            #print "n1b:",b.shape,type(b),b[:2]
            print_single_qpoints = False  # donw for dominique for bcc
            if print_single_qpoints == True and mdsteps == mdstepstocalc.max():
                tmpfilename = "ps"+qpointstring(qpoint)+"__"+str(i)+"__"+str(mdsteps)+".dat"
                parnpsavetxt(tmpfilename,a)

            power_spectrum+=a
        #endp = time.time()
        #print "fft1 ",str((endp-startp)),"sec."

    if not serial: # parallel
        ## DONT! FOR SMALL JOBS THIS IS A BIT FASTER approx 15% but for larger jobs this makes
        ## everything slower by approx 17.5%,  maybe due to using more memory?
        sys.exit()
        ## DONT! FOR SMALL JOBS THIS IS A BIT FASTER approx 15% but for larger jobs this makes
        ## everything slower,  maybe due to using more memory?
        #print printblue('fftparallel !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        startp = time.time()
        def currentfuncc_helper(args):
            return thefunc(*args)
        def thefunc(i,space_fft,mdsteps,return_dictt):
            calc = np.abs(fft((space_fft[i][:mdsteps])[:mdsteps]))**2.
            return_dictt[i] = calc
            return calc

        manager = multiprocessing.Manager()
        return_dictt = manager.dict()
        jobs = []
        for i in np.arange(len(space_fft)): # laeuft ueber alle symmetrieequivalten qs (6 stueck)
            arguments = (i,space_fft,mdsteps,return_dictt)
            p = multiprocessing.Process(target=currentfuncc_helper, args=[arguments])
            jobs.append(p)
            p.start()
        for job in jobs:
            job.join()
        #print "#################33"
        X=np.zeros(mdsteps)
        d = return_dictt
        for k, v in d.items():
            #print k, v
            power_spectrum+=v
        #print "don#################33"
        #print d
        #for k,v in d: X[k] = v
        #print type(d.values())
        #power_spectrum=X
        return_dictt = None
        #endp = time.time()
        #print "fft1 ",str((endp-startp)),"sec."
        #print



    #
    #np.savetxt("psb_pos.real",power_spectrumb.real)

    #power_spectrumb += np.roll(power_spectrumb[::-1],1)
    #print power_spectrumb
    #np.savetxt("ps_pos",power_spectrum)

    execute_timeinversion = True
    if execute_timeinversion:
        #a = np.roll(a[::-1],1)
        # this can only be done as done in the first definition
        #np.savetxt("ps_neg",np.roll(power_spectrum[::-1],1))
        power_spectrum += np.roll(power_spectrum[::-1],1)
        print_single_qpoints = False # for dominique
        if print_single_qpoints == True and mdsteps == mdstepstocalc.max():
            tmpfilename = "ps"+qpointstring(qpoint)+"__"+str(i)+"-__"+str(mdsteps)+".dat"
            parnpsavetxt(tmpfilename,a)
        #power_spectrum+=a
    #end = time.time()
    #print "fft ",mdsteps,i,str((end-startt)),"sec."
    #np.savetxt("ps_pos_neg",power_spectrum)
    #np.savetxt("psb_pos_neg.real",power_spectrumb.real)
    #np.savetxt("psb_pos_neg.imag",power_spectrumb.imag)
    #sys.exit()

    ##########################################################################
    # write out power_spectrum (simple list, only y without x values;
    # length of mdsteps/2 which changes in this loop)
    ##########################################################################
    #power_spectrum[0:10] = 0;  # sometimes very high values , 10 is too much for md with 100 steps
    power_spectrum[:2] = 0;  # power_spectrum is just a list of values (only y without x)
    power_spectrum[-2:] = 0;  # power_spectrum is just a list of values (only y without x)

    ##########################################################################
    # make powerspectrum symmetric (necessery although execute_timeinversion is used!)
    ##########################################################################
    power_spectrum=power_spectrum+power_spectrum[::-1]
    ps=power_spectrum
    factor_psout = 2  # 2 is a save bet since power spectrum is doubled
    filenameout = "ps"+qpointstring(qpoint)+"_"+str(mdsteps)+".dat"

    ##########################################################################
    # WRITE out power_spectrum
    ##########################################################################
    if args.write_full_ps and mdsteps == mdstepstocalc.max():
        startp = time.time()
        print printblue("WRITING0: "+filenameout).rstrip('\n')
        #np.savetxt(filenameout,ps[:len(ps)/factor_psout])  # saves only half the powerspektrum (left peak)
        parnpsavetxt(filenameout,ps[:len(ps)/factor_psout])
        #p = multiprocessing.Process(target=np.savetxt, args=[filenameout,ps[:len(ps)/factor_psout]])
        #p.start()
        #print printblue("WRITING0: "+filenameout+" DONE!")
        print printblue(" DONE!")
        endp = time.time()
        print "WRITING0 ",str((endp-startp)),"sec."

    #sys.exit("ka1")
    ##########################################################################
    # change ps to ps 2d which is smaller and quicker to calculate
    # this is extremely important for crossing of branches
    # xy, ps (need to stay all the time with both peaks until end for fitting)
    #
    # !!!!!!!!!! smooth or not makes at 6 6 20 qp:
    # 0.6  1.6387  0.4432  0.9972  0.016  999999999  0.0  0.554  0.4432 FALSE
    # vs
    # 0.6  1.5973  0.0809  0.5969  0.015  999999999  0.0  0.516  0.0809  TRUE
    #
    # 1 1 20:
    # 0.1  0.1327  0.0015  0.0028  0.00032  999999999  0.0  0.0013  0.0015 FALSE
    # vs
    # 0.1  0.1453  0.0251  0.0265  0.0001  999999999  0.0  0.0013  0.0251  TRUE
    #
    # 10 10 20:
    # 1.0  1.404  0.0126  0.0133  0.0011  999999999  0.0  0.0007  0.0126 TRUE good 2.76439393939
    # vs
    # 1.0  1.4953  0.0271  0.0277  0.0013  999999999  0.0  0.0007  0.0271 FALSE good 2.84533333333
    #
    # 9 9 20:
    # 0.9  1.8893  0.0252  0.0265  0.001  999999999  0.0  0.0013  0.0252 good 2.89 (FALSE)
    # vs
    # 0.9  1.728  0.0156  0.0156  0.00043  999999999  0.0  0.0  0.0156  (True)
    #
    # 8 8 20:
    # 0.8  1.8253  0.0655  0.0668  0.00067  999999999  0.0  0.0013  0.0655 TRUE
    # vs
    # 0.8  1.5727  0.2932  0.2932  0.00095  999999999  0.0  0.0  0.293 FALSE
    ##########################################################################
    smooth=True   # this seems closer to Experiment
    smooth=False  # this seems more accurate from Theoretical standpoint
    #print printblue("smooth is "+str(smooth))
    #np.savetxt("psin",ps)
    #np.savetxt("psin-",ps[::-1])

    ps = power_spectrum
    #############################################################################
    # HIER KANN ps zwar was verschoben sein, aber xy ist genau richtig!
    #############################################################################
    xy,x,ps = powerspectrum_to_powerspectrum_2d_sparse(ps,smooth=smooth,mindestens=5000)  # y = ps, xy should in principally also be scaled, but curr unused


    #np.savetxt("psout",ps)
    #np.savetxt("psout-",ps[::-1])
    #np.savetxt("xy",xy)
    #np.savetxt("xy-",xy[::-1])
    #np.savetxt("x",x)
    #np.savetxt("ps",ps)
    #np.savetxt("xycur.dat",xy)
    #sys.exit()

    ##########################################################################
    # WRITE pssparse if necessary
    ##########################################################################
    writeps_smooth_in = False
    #if args.write_smooth_in: writeps_smooth_in = True
    #if args.write_full_ps and mdsteps == mdstepstocalc.max(): writeps_smooth_in = True
    check_qpoint_for_crossing = check_qpoints_for_crossing(qpoint, args.supercell, args.structure)
    if check_qpoint_for_crossing == True and mdsteps == mdstepstocalc.max(): writeps_smooth_in = True
    if (args.write_smooth_in or args.write_full_ps) == True and mdsteps == mdstepstocalc.max():
        print printblue("WRITING1: "+filenameout+".smooth_in ... (this is written at the last step when mdsteps == mdstepstocalc.max()) "+qpointstring(qpoint))
        parnpsavetxt(filenameout+".smooth_in",xy[:xy.shape[0]/2])

    ps=ps/ps.max()*.9; # full with peaks left and right
    ##########################################################################
    # CROSSING
    # WRITE out power_spectrum_cut
    # (needs to stay all the time with both peaks until end for fitting)
    ##########################################################################
    #np.savetxt("psin",ps)
    #print "1 len(x):",len(x),x
    #print "1 len(ps):",len(ps),ps
    if check_qpoint_for_crossing == True:
        if mdsteps == mdstepstocalc.max():
            checkwriteps = mdsteps
        else:
            checkwriteps = False
        # zurzeit macht get_max_min_max_coords_idx_of_full_ps_when_bandcrossing aus einer ungerade anzahl an mdsteps eine gearade anzahl.
        if mdsteps == mdstepstocalc.max():
            print printblue("WRITING1: "+filenameout+".smooth_in_ ... (this is written at the last step when mdsteps == mdstepstocalc.max()) "+qpointstring(qpoint))
            parnpsavetxt(filenameout+".smooth_in",ps[:xy.shape[0]/2])
        ps = get_max_min_max_coords_idx_of_full_ps_when_bandcrossing(xy,x,ps,qpointstr,checkwriteps,verbose=False)
        #if (args.write_smooth_in or args.write_full_ps) == True and mdsteps == mdstepstocalc.max():
        if mdsteps == mdstepstocalc.max():
            print printblue("WRITING1: "+filenameout+".smooth_in_cut ... (this is written at the last step when mdsteps == mdstepstocalc.max()) "+qpointstring(qpoint))
            parnpsavetxt(filenameout+".smooth_in_cut",ps[:xy.shape[0]/2])
        #if args.verbose:
        #    print "1 ps.max()",ps[:len(ps)/factor_psout].max()
        # if args.write_full_ps and mdsteps == mdstepstocalc.max():
        #if mdsteps == mdstepstocalc.max():  # better write in this cases
        #    print "qp:",qpoint,type(qpoint)
        #    filenameoutcut = "ps"+qpointstr+"_"+str(mdsteps)+".dat.cut"
        #    print printblue("WRITING1: "+filenameoutcut)
        #    np.savetxt(filenameoutcut,ps[:len(ps)/factor_psout])  # saves only half the powerspektrum (left peak)

    #print "2 len(x):",len(x),x
    #print "2 len(ps):",len(ps),ps
    #np.savetxt("psout"+str(ps)
    #np.savetxt("psout-",ps[::-1])

    #sys.exit()

    # set second part to 0
    #print "lenkk",len(x),len(ps)
    #np.savetxt("xy"+str(mdsteps),np.transpose((x,ps)))

    ##########################################################################
    # cut second part of spectrum out
    # ! dont do that ! this sets minimum background to wrong value
    # this was unexpected since get_goodsmoothing should only consider first part of spectrum
    ##########################################################################
    # DONT DO THIS ps[len(ps)/factor_psout:] = 0  DONT DO THIS
    #print "ps.max():",ps.max()


    ##########################################################################
    # write out power_spectrum_smoothed (simple list, only y without x values
    # length of mdsteps/2 which changes in this loop)
    # out1 = smoothing_power_spectrum(ps,smoothingfact)  # out1 is as the powerspectrum have both peaks and has exactly as many as many points
    ##########################################################################
    if True and mdsteps > 999:
        #print "makeps!!!!!!!!!!!!!!!!!!!"
        #np.savetxt("xycur.dat",ps)
        #ps=np.loadtxt("xyold.dat")
        #print "mdsteps:",mdsteps
        #if args.verbose:
            #print "###########################################################################################################################"
            #print "################################## create psin.dat",qpointstring(qpoint)
            #np.savetxt("psin"+qpointstring(qpoint)+".dat",ps)  # aus ps kann man eigntlich schon sehr gut ablesen
        #start= time.time()
        # sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin,minsmoothfunc = get_goodsmoothing_for_ps(ps,args.dt,mdsteps,allowed_func = 1, allowed_der = False,args=args,stringadd=['min',str(mdsteps)])
        #
        ##############################################################################################
        # in general the error due to min <--> vs <--> good is small
        # especially for short mds the error might howevery be very large wor allowed_der = False
        # therefore try first with allowed_der = 1
        # if this should not work well -> try allowed_der = 0 for sgimamin
        ##############################################################################################
        #sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin,minsmoothfunc = get_goodsmoothing_for_ps(ps,args.dt,mdsteps,allowed_func = 1, allowed_der = 1 ,args=args,stringadd=['min',str(mdsteps)])
        sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin,minsmoothfunc = get_goodsmoothing_for_ps(ps,args.dt,mdsteps,allowed_func = 1, allowed_der = 0 ,args=args,stringadd=['min',str(mdsteps)])
        #xo = x[:len(x)/2]
        #m = minsmoothfunc[:len(minsmoothfunc)/2]
        #np.savetxt("7720_"+str(mdsteps)+"_min.dat",np.transpose([xo,m]))
        #print "kk",sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin
        #if args.verbose:
            #print "###########################################################################################################################"
            #print printblue("    -----> min:"),sigmamin,[ round(elem, 2) for elem in [fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin]]
        sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood,goodsmoothfunc = get_goodsmoothing_for_ps(ps,args.dt,mdsteps,allowed_func = 0, allowed_der = 0,args=args,stringadd=['good',str(mdsteps)])
        #g = goodsmoothfunc[:len(goodsmoothfunc)/2]
        #np.savetxt("7720_"+str(mdsteps)+"_good.dat",np.transpose([xo,g]))
        #print "kl",sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood
        #end = time.time()
        #print "find lw ",str((end-start)),"sec."
        #if args.verbose:
        #    print "###########################################################################################################################"
            #print printblue("    -----> good:"),sigmagood, [ round(elem, 2) for elem in [fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood]]
        #print  printblue("good "+str(round(ltgood,3))+" min:"+str(round(ltmin,3))), "white?"
        #sys.exit()
        if check_qpoint_for_crossing == True:
            #sigmagood,sigmamin = sigmamin,sigmagood
            #fqgood, fqmin = fqmin, fqgood  # this is an inbetween, a thorough solution is necessary using michales skript
            #ltgood, ltmin = ltmin, ltgood
            sigmagood = sigmamin
            fqgood = fqmin
            ltgood = ltmin
            #print "ltgoodmin,ltgoodmax",ltgoodmin,ltgoodmax
            #print "ltminmin,ltminmax: ",ltminmin,ltminmax
            ltgoodmin,ltgoodmax = ltminmin,ltminmax
            #sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood,goodsmoothfunc = sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin,minsmoothfunc
            #if args.verbose:
            #    print "################################## FINAL--good--3",sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax#,xmaxwritegood,goodsmoothfunc
            #    print "################################## FINAL--min---3",sigmamin,   fqmin, ltmin, ltminmin, ltminmax#, xmaxwritemin, minsmoothfunc
        #print sigmamin,sigmagood,ltmin,ltgood
        #out1out = out1[:len(out1)/factor_psout]
        #############################################################
        # write here the smoothigfunc, further down the lifetimes
        #############################################################
        finaldata_min  = [sigmamin, fqmin,ltmin,ltminmin,ltminmax,int(mdsteps)]
        finaldata_good = [sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,int(mdsteps)]
        if mdsteps == mdstepstocalc.max():
            filenameoutsmooth = filenameout+".smooth_good"+str(sigmagood)+".dat"
            filenameoutmin = filenameout+".smooth_min"+str(sigmamin)+".dat"
            #print printblue("WRITING3:"+filenameoutsmooth)

            #np.savetxt(filenameoutsmooth,goodsmoothfunc[:xmaxwritegood])   # old way for 1d
            #np.savetxt(filenameoutmin,minsmoothfunc[:xmaxwritemin])        # old way for 1d
            # adapt to 2d which is smaller
            print printblue("WRITING3: "+filenameoutsmooth+" and "+filenameoutmin)
            #print "x:",len(x),"||",x
            #print "goodsmoothfunc:",len(goodsmoothfunc),"||",goodsmoothfunc
            #print "minsmoothfunc:",len(minsmoothfunc),"||",minsmoothfunc
            parnpsavetxt(filenameoutsmooth,np.transpose(np.array([x,goodsmoothfunc]))[:xmaxwritegood],fmt='%.1f %.5f')
            parnpsavetxt(filenameoutmin,   np.transpose(np.array([x,minsmoothfunc]))[:xmaxwritemin],fmt='%.1f %.5f')

    if type(return_dict) != bool:
        return_dict[mdsteps] = [finaldata_good,finaldata_min,filenameout]
    return finaldata_good,finaldata_min,filenameout

def smoothing_power_spectrum(power_spectrum,sigma):
    ''' smoothes a gaussian by sigma '''
    steps = power_spectrum.shape[0]                 # 1000
    nu = np.arange(0,steps)/float(steps)            # [ 0., 0.001, ..., 0.999 ]
    kern = np.exp(-(np.mod(nu+0.5,1)-0.5)**2/(2*sigma**2))  # Gaussian distribution in 1-D
    kern = kern/np.sum(kern)
    ps = np.real(ifft(fft(power_spectrum)*fft(kern)))
    ps=ps/np.max(ps)*.9;
    return ps

def get_uncertainty_smoothed_powerspectrum(power_spectrum_smooth,dt,sigma = "",args=False,stringadd='',verbose=False,color='blue'):
    ''' - dt is the sampled timestep in femtoseconds
        - power_spectrum_smooth is the full power spectrum having both peaks (left and right)
        - sigma is only for plotting to screen
    '''
    verbose = args.verbose
    data = power_spectrum_smooth
    #np.savetxt("ka"+str(sigm)+".dat",data)  #dies schein schon oversmoothed zu sein

    xvalues=data.shape[0]
    #print "xvalues all:",xvalues
    datahalf = data[0:xvalues/2]   # only take one peak
    #np.savetxt("kb.dat",datahalf)
    #print "datahalf:",datahalf

    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]


    # x index of absolute maximum
    x_ind_max = np.where(datahalf == datahalf.max())[0][0]  # 963  (9.63 THz)
    y_val_max = datahalf[x_ind_max]
    #if args.verbose:
    #    print "x_ind_max 77:",x_ind_max

    if x_ind_max == 0:
        return 0 ,0, 0 ,xvalues, 0,0,0,99999999999999999999999999

    def printit(x_ind_max,y_val_max,x_ind_leftdata_min,y_val_leftdata_min,y_val_leftdata_max,y_lifetime_max,y_lifetime_min):
        print "############################################################"
        print "# analyse the hight of the background"
        print "#                          .---------------------- y_val_max:",y_val_max
        print "#                         .|."
        print "#                        . | ."
        print "#                       .  |  ."
        print "#                     .    |   ."
        print "#                   .      |    ."
        print "#                 .        |     ."
        print "#               .          |      ."
        print "#.---------- . ------------|------ .------------- y_val_leftdata_max:",y_val_leftdata_max
        print "#  .      .                |        ."
        print "#     .------------------------------- . -------- y_val_leftdata_min:",y_val_leftdata_min,"(therefore linewidth measured at y = y_lifetime_max:",y_lifetime_max," y_lifetime_min:",y_lifetime_min,")"
        print "#     |                    |               .  .   ."
        print "############################################################"
        print "#     x--------------------|----------------------- x_ind_leftdata_min:",x_ind_leftdata_min
        print "#"
        print "#                          x                        x_ind_max:",x_ind_max
        print "#x                                                  x_val_leftdata_max:"
        print "############################################################"
        return

    y_val_leftdata_min = datahalf[:x_ind_max].min()
    x_ind_leftdata_min = np.where(datahalf[:x_ind_max] == y_val_leftdata_min)[0][0]

    if x_ind_leftdata_min == 0:
        y_val_leftdata_max = datahalf[0]
    else:
        y_val_leftdata_max = datahalf[:x_ind_leftdata_min].max()

    #############################################################
    ## maximal linewidths if background not substracted
    ## NOT substracting the background is not something the experimentalists can do!!!
    ## THEREFORE: get another max .... TODO
    #############################################################
    y_lifetime_max              = datahalf.max()/2.
    y_lifetime_max              = datahalf.max()/2.+y_val_leftdata_min/2  # Therefore changed definition
    y_ind_max_over2_left_max    = find_nearest(datahalf[0:x_ind_max], y_lifetime_max)
    y_ind_max_over2_right_max   = find_nearest(datahalf[x_ind_max:-1], y_lifetime_max)
    x_ind_left_over2_max        = np.where(datahalf[0:x_ind_max] == y_ind_max_over2_left_max)[0][0]
    x_ind_right_over2_max       = np.where(datahalf == y_ind_max_over2_right_max)[0][0]

    #############################################################
    ## minimal linewidths if maximal background is substracted
    #############################################################
    y_lifetime_min = datahalf.max()/2.+y_val_leftdata_max/2.
    y_ind_max_over2_left_min  = find_nearest(datahalf[0:x_ind_max], y_lifetime_min)
    y_ind_max_over2_right_min = find_nearest(datahalf[x_ind_max:-1], y_lifetime_min)
    x_ind_left_over2_min = np.where(datahalf[0:x_ind_max] == y_ind_max_over2_left_min)[0][0]
    x_ind_right_over2_min = np.where(datahalf == y_ind_max_over2_right_min)[0][0]

    if args.verbose > 3:
        printit(x_ind_max,y_val_max,x_ind_leftdata_min,y_val_leftdata_min,y_val_leftdata_max, y_lifetime_max,y_lifetime_min)


    #############################################################
    ## "normal/intermediate"  linewidths if background is substracted
    #############################################################
    y_lifetime = datahalf.max()/2.+y_val_leftdata_min/2.
    #print "datahalf.max()/2     :",datahalf.max()/2.
    #print "y_val_leftdata_min   :",y_val_leftdata_min
    #print "y_lifetime           :",y_lifetime
    #sys.exit()

    y_lifetime_third = datahalf.max()/3.
        #return freq,lifetime, lifetimemin,xvalues, r+l,ll+lr+rl+rr, x_ind_right_over1000

    #############################################################
    # dis ist nur schoen (und korrekt) bei einem relativ glatten gaussian
    #############################################################
    # x index of left
    #print "x_ind_max:",x_ind_max
    #print "datahalf[0:x_ind_max]",datahalf[0:x_ind_max]
    #print "datahalf.max()/2.:",datahalf.max()/2.
    y_ind_max_over2_left  = find_nearest(datahalf[0:x_ind_max], y_lifetime)
    y_ind_max_over3_left  = find_nearest(datahalf[0:x_ind_max], y_lifetime_third)
    #print "y_ind_max_over2_left",y_ind_max_over2_left
    #print "y_ind_max_over3_left",y_ind_max_over3_left
    x_ind_left_over2 = np.where(datahalf[0:x_ind_max] == y_ind_max_over2_left)[0][0]
    x_ind_left_over3 = np.where(datahalf[0:x_ind_max] == y_ind_max_over3_left)[0][0]

    y_ind_max_over2_right = find_nearest(datahalf[x_ind_max:-1], y_lifetime)
    y_ind_max_over3_right = find_nearest(datahalf[x_ind_max:-1], y_lifetime_third)

    y_ind_max_over2_right = find_nearest(datahalf[x_ind_max:-1], y_lifetime)
    y_ind_max_over3_right = find_nearest(datahalf[x_ind_max:-1], y_lifetime_third)


    #print "y_ind_max_over2_right",y_ind_max_over2_right
    x_ind_right_over2 = np.where(datahalf == y_ind_max_over2_right)[0][0]
    x_ind_right_over3 = np.where(datahalf == y_ind_max_over3_right)[0][0]
    #sys.exit()
    #print "x_ind_left_over2",x_ind_left_over2,"y:",y_ind_max_over2_left
    #print "x_ind_right_over2",x_ind_right_over2, "y:",y_ind_max_over2_right


    l = len(np.where(np.diff(datahalf[x_ind_left_over2:x_ind_max])<0)[0])    # using x_ind_left_over2 insted of x_ind_left_over3 makes a TA2 realistic errorbar
    #print l
    #print "bbb"
    r = len(np.where(np.diff(datahalf[x_ind_max:x_ind_right_over2])>0)[0])   # using x_ind_left_over2 insted of x_ind_left_over3 makes a TA2 realistic errorbar
    #print "r+l:",r+l
    #sys.exit()

    #############################################################
    # checke die ableitung des peaks
    #############################################################
    ###############################
    # erstmal die linke haelfte
    ###############################
    #print x_ind_left_over3,x_ind_max,x_ind_right_over3
    ll = lr = 0.
    #print "x_ind_max:",x_ind_max
    if args.verbose > 3:
        print "x_ind_left_over2,y_ind_max_over2_left:",x_ind_left_over2,y_ind_max_over2_left
        print "x_ind_left_over3,y_ind_max_over3_left:",x_ind_left_over3,y_ind_max_over3_left

    #if (x_ind_max - x_ind_left_over3) > 1:
    #    datatake = datahalf[x_ind_left_over3:x_ind_max]
    if (x_ind_max - x_ind_left_over2) > 1:
        datatake = datahalf[x_ind_left_over2:x_ind_max]
        ableitung1 = np.diff(datatake)  # 95:[97] funktioniert k --> ableitung ist nur ein wert
        #print ableitung1
        ableitung1_ind_max = np.where(ableitung1 == ableitung1.max())[0][0]
        ableitung1l = ableitung1[:ableitung1_ind_max]
        ableitung1r = ableitung1[ableitung1_ind_max:]
        ll = len(np.where(np.diff(ableitung1l)<0)[0])
        lr = len(np.where(np.diff(ableitung1r)>0)[0])

        #np.savetxt("ka.dat",datahalf)
        #np.savetxt("kb.dat",datatake)
        #np.savetxt("kc.dat",ableitung1)
        #np.savetxt("kd.dat",ableitung1l)
        #np.savetxt("ke.dat",ableitung1r)
        #sys.exit()
    rl = rr = 0.
    #if (x_ind_right_over3-x_ind_max) > 1:
    #    datatake = datahalf[x_ind_max:x_ind_right_over3]
    if (x_ind_right_over2-x_ind_max) > 1:
        datatake = datahalf[x_ind_max:x_ind_right_over2]
        ableitung1 = np.diff(datatake)  # 95:[97] funktioniert k --> ableitung ist nur ein wert
        #np.savetxt("ka.dat",datatake)
        #np.savetxt("kb.dat",ableitung1)
        #print ableitung1
        ableitung1_ind_max = np.where(ableitung1 == ableitung1.min())[0][0]
        ableitung1l = ableitung1[:ableitung1_ind_max]
        ableitung1r = ableitung1[ableitung1_ind_max:]

        #np.savetxt("kc.dat",ableitung1l)
        #np.savetxt("kd.dat",ableitung1r)


        rl = len(np.where(np.diff(ableitung1l)>0)[0])
        rr = len(np.where(np.diff(ableitung1r)<0)[0])
        #print rl,rr
    #rl = rr = 0.
    #print ll,rr

    #############################################################
    # dis ist zuverlaessiger bei viel zu dicht gesampleten
    #############################################################
    #print
    #print "bbb:",y_lifetime
    #if len(np.where(datahalf>y_lifetime)[0]) == 0:
    #    np.savetxt("kkk.dat",datahalf)
    #print "ccc:",len(datahalf)
    #print "kkk:",len(np.where(datahalf>y_lifetime)[0])
    #if len(np.where(datahalf>y_lifetime)[0]) == 0:
    #    np.savetxt("kkk.dat",datahalf)
    #    print "now 0"

    #print "aaa:",np.where(datahalf>y_lifetime)[0][0]
    absmin = np.where(datahalf>y_lifetime)[0][0]     # nur die spitze des gaussians (linke seite)
    absmax = np.where(datahalf>y_lifetime)[0][-1]    # nur die spitze des gaussians (rechte seite)

    if len(np.where(datahalf>y_lifetime)[0]) == 0:
        np.savetxt("kkk.dat",datahalf)
        sys.exit("now 0")

    #mult = dt*10 / float(xvalues)  ## Wrong!  see example 6x6x6sc 6_0_0 qpoint between dt={10,20}
    mult = 1000./(args.dt*float(xvalues))
    freq = x_ind_max * mult
    lifetime     = (x_ind_right_over2 - x_ind_left_over2) * mult
    lifetimemax  = (x_ind_right_over2_max - x_ind_left_over2_max) * mult
    lifetimemin0 = (x_ind_right_over2_min - x_ind_left_over2_min) * mult
    lifetimemin1 = (absmax - absmin) * mult          # top of gaussian
    lifetimemin = np.array([lifetimemin0, lifetimemin1]).min()
    #print lifetimemin,lifetime, lifetimemax
    #vor_freq = x_ind_max/float(xvalues)
    #vor_lifetime = (x_ind_right_over2 - x_ind_left_over2)/float(xvalues)
    #vor_lifetimemin = (absmax - absmin)/float(xvalues)
    ##print "vor_freq:",vor_freq
    ###print "vor_lifetime:",vor_lifetime
    ###print "vor_lifetimemin:",vor_lifetimemin
    #print "Smoothed function:"
    #if args.verbose and verbose:
    if args.verbose > 2:
        if stringadd[0] == 'min':
            stringadd[0] = 'min '
        #ljust adds zeros (0) to the end
        #ka = printblue(stringadd)
        if color == 'blue':
            print "Freq "+printblue(stringadd[0])+" "+stringadd[1].ljust(10)+":",str(round(freq,2)).ljust(4,'0'),"THz; width: ",str(round(lifetime,2)).ljust(4,'0'),"THz", "||", x_ind_left_over2,x_ind_max,x_ind_right_over2,"||", "sigma:",str(sigma),"r+l",r+l,"rr,ll",ll+lr+rl+rr,"ll:",ll,"rr:",rr,"lr:",lr,"rl:",rl
        if color == 'red':
            print printred("Freq ")+printred(stringadd[0])+" "+stringadd[1].ljust(10)+":",str(round(freq,2)).ljust(4,'0'),printred("THz; width: "),printred(str(round(lifetime,2)).ljust(4,'0')),"THz", "||", x_ind_left_over2,x_ind_max,x_ind_right_over2,"||", "sigma:",str(sigma),"r+l",r+l,"rr,ll",ll+lr+rl+rr,"ll:",ll,"rr:",rr,"lr:",lr,"rl:",rl
            #"lifetime min:",lifetimemin,\
    #print "indizes:",x_ind_left_over2,"<-",x_ind_max,"->",x_ind_right_over2 #,"from" ,xvalues,"xvalues"

    # x cutoff when writing files to harddisk
    y_ind_max_over1000_right = find_nearest(datahalf[x_ind_max:-1], datahalf.max()/1000.)
    x_ind_right_over1000 = np.where(datahalf == y_ind_max_over1000_right)[0][0]
    #print "xxx:",x_ind_right_over1000
    #print "--->>",freq,lifetime, lifetimemin, lifetimemax ,xvalues, r+l,ll+lr+rl+rr, x_ind_right_over1000
    return freq,lifetime, lifetimemin, lifetimemax ,xvalues, r+l,ll+lr+rl+rr, x_ind_right_over1000

def find_nearest(array,value,returnindex=False):
    ''' find value in array
    '''
    idx = (np.abs(array-value)).argmin()
    if returnindex:
        return idx
    return array[idx]

def get_max_min_max_coords_idx_of_full_ps_when_bandcrossing(xy,x,ps,qpointstr,checkwriteps,verbose=False):
    ''' expecially for the 900K runs we also want to find the minimum of the smoothed curve
        a) find min of smoothed curve
        b) make averages of smoothed curve and find the true minimum around the minimum (target would be to have say 100 values in the range of +-1% around minimum
    '''
    if verbose:
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    psout=ps
    #np.savetxt("ps",ps)
    #print "33:",ps.shape,x.shape,xy.shape
    ps=ps[:len(ps)/2]
    x=x[:len(x)/2]
    #print "x:",x.min(),x.max()
    xy=xy[:len(xy)/2]
    #print "44:",ps.shape,x.shape,xy.shape
    #np.savetxt("aa2.dat",ps)
    #np.savetxt("aa3.dat",pss)
    #np.savetxt("aa4.dat",xy)
    def get_max_min_max_of_function(ps,verbose=False,above_y_set=1):
        #np.savetxt("ps",ps)
        verbose = False
        xmax = np.where(ps==ps.max())[0][0]
        if verbose:
            print
            print "INNNN xmax:",xmax
        ymax = ps[xmax]
        if verbose:
            print "xmax,ymax:",xmax,ymax
        # now go to left and right of smoothfunc
        yminleft  = ps[:xmax].min()
        yminright = ps[xmax:].min()
        yminright = ps[-1]
        if verbose:
            print "yminleft :",yminleft
            print "yminright:",yminright
        xminleft  = np.where(ps[:xmax]==yminleft )[0][0]
        #xminrigth = np.where(ps[xmax:]==yminright)[0][0]+xmax
        xminrigth = len(ps)

        if verbose:
            print "xminleft :",xminleft
            print "xminrigth:",xminrigth
            print "now make function which takes all values above"
        above = np.array([0.8,.7,.6,.5,.4,.3,0.2,0.1])

        if verbose:
            print "above_y_set:",above_y_set
        if above_y_set == 1:
            above = (np.arange(9)/10.)[np.arange(9)/10.>np.array([yminleft,yminright]).min()][::-1]
        else:
            above = (np.arange(18)/20.)[np.arange(18)/20.>np.array([yminleft,yminright]).min()][::-1]
        def get_teiler_daten(ps,mindestens):
            lenps = len(ps)
            teiler = 1
            for i in np.arange(16)+1:
                tryteiler = 1*10**i
                #print i,teiler,lenps/teiler,(lenps/teiler)>=5000
                if ((lenps/tryteiler)>=mindestens) == True:
                    #print i,tryteiler,"Ja"
                    teiler = tryteiler
                else:
                    #print i,tryteiler,"NEIN"
                    break
            #print "teiler:",teiler
            N = teiler
            return N

        N = get_teiler_daten(ps,mindestens=500)
        if verbose:
            print "N:",N,len(ps)

        #print "NNNin:",N,len(ps)
        if N > 1:
            ps = mysmooth(ps,N+1)  # mysmooth(y,N*3+1)  ---> dieser ist also immer ungerade!
        #print "NNNout:",N,len(ps)

        if verbose:
            print "N:",N,len(ps)
        #np.savetxt("ps2",ps2)
        if verbose:
            print "above:",above
            print len(above)/2
            #la = len(above)/2
            #von=int(la*1./3.)
            #bis=int(la+la*1./3.)
            #print above[von:bis]
            #print above[:von][::-1]
            #print "get middle 30 percent:",above[2la:2./3.,la+la

        for above_y in above:
            above_y_out = above_y
            allv = np.where(ps[xminleft:xminrigth]>above_y)[0]
            if len(allv) <=1:
                continue
            diff = np.diff(allv)
            diffmax = diff.max()
            if verbose:
                print "above_y:",above_y,"len allv:",len(allv),allv,diff,"len diff:",len(diff),"--------> diffmax:",diffmax
            if diffmax > 1:
                #print above_y,"allv:",allv,len(allv)
                firsthoeckerleft  = xminleft + allv[0]
                secondhoeckerright= xminleft + allv[-1]
                #print "firsthoeckerleft  :",firsthoeckerleft
                #print "secondhoeckerright:",secondhoeckerright
                xmin = np.where(ps[firsthoeckerleft:secondhoeckerright]==ps[firsthoeckerleft:secondhoeckerright].min())[0][0]
                xmin = xmin+firsthoeckerleft
                #if verbose:
                #    print above_y,"----------------------------> xmin:",xmin
                #break
                #print "xmin:",xmin
                #print "diffmax:",diffmax
                #print "allv[0]:",allv[0]
                #print "diffmax idx:",np.where(diff==diffmax)[0][0]
                #print "firsthoeckerleft  :",firsthoeckerleft
                #print "secondhoeckerright:",secondhoeckerright

                xmaxleft = np.where(ps[:xmin]==ps[:xmin].max())[0][0]
                y_xmaxleft = ps[xmaxleft]
                #print "CROSSING --> xmaxleft:",xmaxleft,y_xmaxleft
                y_xmin = ps[xmin]
                #print "CROSSING --> xmin:",xmin,y_xmin
                xmaxright = np.where(ps[xmin:]==ps[xmin:].max())[0][0]+xmin
                y_xmaxright = ps[xmaxright]
                #print "CROSSING --> xmaxright:",xmaxright,y_xmaxright
                if verbose:
                    print "ka1:",(xmin-xmaxleft)
                    print "ka2:",xminrigth
                    print "xmin:",xmin
                    print "xmaxleft:",xmaxleft
                    print "xminrigth:",xminrigth
                    print "(xmin-xmaxleft)/float(xminrigth)"
                aaa = (xmin-xmaxleft)/float(xminrigth)
                bbb = (xmaxright-xmaxleft)/float(xminrigth)
                ccc = (xmin-xmaxleft)/float(len(ps))
                ddd = (xmaxright-xmin)/float(len(ps))
                ccc = np.array([ccc,ddd]).min()
                if verbose:
                    print above_y,"----------------------------> leftpeak,min,rightpeak (X):",xmaxleft ,xmin,xmaxright,"|a>0.002?:",aaa,"|b",bbb,"c>0.01?:",ccc
                    print above_y,"----------------------------> leftpeak,min,rightpeak (Y):",y_xmaxleft,y_xmin,y_xmaxright,"|"
                if aaa > 0.002 and ccc > 0.002:  # a sollte gleich c sein
                    # a gut ########################################################
                    # a ist 0.0299818009697 bei 8_8_20 (t2_t2_0)    # 300K
                    # a ist 0.0876429458276 bei 9_9_20              # 300K
                    # a ist 0.0044179368235 bei 4_4_8 (t2)          # 900K (DFT GGA)
                    #
                    # a ungut #######################################################
                    # a is 0.0068 (== gut)
                    # ein zu kleines a ist 0.000721783607492
                    # ein zu kleines a ist 0.000792169219702
                    # ein zu kleines a ist 0.011 bei 9_9_20 900K /check_4_temperatuer_effects/900K_LDA_4.08

                    # c gut @ 900K 4_4_8    0.0031
                    # c gut @ 300K 10_10_20 0.114 (smothed 0.105)
                    # c gut @ 300K  8_8_20  0.049 (smoothed 0.038)
                    #                                                   0.00066 kann aber auch gut sein....(bei wenigen daten)
                    # c ungut @ 900K  10_10_20  0.055(=gut)  (smoothied 0.00069=ungut)
                    break

                #if aaa < 0.02 or ccc < 0.01:
                #    # checke erst ob zwischen max1 und max2 evtl. ein anderer wert hoecher ist als
                #HIER HIER HIER
                #@@@if above_y == 0.1:
                #@@@    # checke erst ob zwischen max1 und max2 evtl. ein anderer wert hoecher ist als
                #@@@    ps[:xmaxleft


        return xmaxleft ,xmin,xmaxright,xminleft,xminrigth,yminleft,len(ps),y_xmaxleft,y_xmin,y_xmaxright,above_y

    #########################################
    # eigentliches powerspektrum
    #########################################
    #xmaxleft1,xmin1,xmaxright1,xminleft1,xminrigth1,yminleft1,lenps1 = get_max_min_max_of_function(ps,verbose=False)
    #print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&> leftpeak,min,rightpeak (X):",xmaxleft1,xmin1,xmaxright1,"|a",(xmin1-xmaxleft1)/float(xminrigth1),"|b",(xmaxright1-xmaxleft1)/float(xminrigth1)

    #########################################
    # gesmoothtes powerspektrum
    # lets use this minimum in general since it is better
    #########################################

    #print "x:",x.min(),x.max()
    #print "------------------"
    #np.savetxt("pss",ps)
    #print "------------------"
    try:
        xmaxleft,xmin,xmaxright,xminleft,xminrigth,yminleft,lenps,y_xmaxleft,y_xmin,y_xmaxright,above_y = get_max_min_max_of_function(ps,verbose=False)
    except UnboundLocalError:
        #print "HHHHHHHHHHHHHHHHHHHHHHHHH"
        try:
            xmaxleft,xmin,xmaxright,xminleft,xminrigth,yminleft,lenps,y_xmaxleft,y_xmin,y_xmaxright,above_y = get_max_min_max_of_function(ps,verbose=False,above_y_set=2)
        except UnboundLocalError:
            #print "YYYYYYYYYYYYYYYYYYYYY"
            print_warning("IT SEMMS THERE IS ONLY ONE MAXIMUM for ps lenghts of "+str(len(psout)))
            return psout/psout.max()*.9

    #if verbose:
    xmaxleft2 = x[xmaxleft]
    xmin2 = x[xmin]
    xmaxright2 = x[xmaxright]
    #if verbose:
    if False:
        #np.savetxt("ps448_"+str(len(ps)),ps)
        print "left, min, right:",xmaxleft2,xmin2,xmaxright2,"--> above_y:",above_y,"--> len",len(ps),"--> totallen:",len(ps)*2.
        print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&> leftpeak,min,rightpeak (X):",xmaxleft2,xmin2,xmaxright2,"|a",(xmin-xmaxleft)/float(xminrigth),"|b",(xmaxright-xmaxleft)/float(xminrigth)
        print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&> leftpeak,min,rightpeak (Y):",y_xmaxleft,y_xmin,y_xmaxright,"<<<- diese Y sind egal da nur die X'es gebraucht werden und das ps zurechtzuschneiden"
        np.savetxt("pstest",ps)

    #print "x:",x.min(),x.max()
    #np.savetxt("ps2",ps2)
    #np.savetxt("ps3",mysmooth(ps2,101))
    a = ps[:xmaxleft][2:]  # first 2 values are 0!!
    xminfueruntergrund = np.where(a==a.min())[0][0]
    y = a[xminfueruntergrund]
    if False:
        #print "a:",a
        print "xminfueruntergrund:",xminfueruntergrund,"y:",y


    # change ps to psout!!! to keep the length!
    #ps[xminfueruntergrund:xmin] = y
    psout[0:xmin] = y
    psout[-xmin:] = y
    #print "22;",x.shape,ps.shape
    def xyplot(x,y):
        return np.transpose(np.array([x,y]))
    if verbose:
        print "psmax :",psout.max()
    #print "psmax2:",ps[:len(ps)/2.max()
    #ps=ps/ps[:len(ps)/2].max()*.9; # full with peaks left and right
    psout=psout/psout.max()*.9; # full with peaks left and right

    if checkwriteps:
        if verbose:
            print "che",checkwriteps
        #np.savetxt("ps"+qpointstr+"_"+str(checkwriteps)+".dat.cut",xyplot(x,ps))
    #pssave=pssave/pssave[:len(pssave)/2].max()*.9; # full with peaks left and right
    #ps=ps/ps[:len(ps)/2].max()*.9; # full with peaks left and right
    # funktion spiegeln
    #out = np.concatenate([psout,ps[::-1]])
    return psout

def get_goodsmoothing_for_ps_with_crossing(ps):
    min_smoothing = 1./len(ps)
    order_magnitude = -1.*(int(math.log10(min_smoothing))-1)
    #print order_magnitude
    order_magnitude = np.arange(1,order_magnitude+1)[::-1]
    #print "ord:",order_magnitude
    for ido,i in enumerate(order_magnitude):
        #print "------>",ido,i,10**(-i)
        order_magnitude[ido] = 10**(-i)
    loop_smoothing = order_magnitude
    #print "2:",loop_smoothing

    #print "STARTING loop 0 -> (logarithmic loop) 1e-07 .. 1e-01" #,loop_smoothing
    working = []
    minimum = []
    x_max1list = []
    x_max2list = []
    minimumalt = []
    for smoothingfact in loop_smoothing[::-1]:   # [::-1] heisst einfach nur die loopsmoothing andersrum lesen, alle werte bleiben erhaleten
    #for smoothingfact in loop_smoothing:
        print "smoothingfact:",smoothingfact
        if len(working) >= 2: continue   # to orders of magnitude are enough, the third one often fails
        out = smoothing_power_spectrum(ps,smoothingfact)  # ps and out1 are the full power spectrum (left and right peak)
        #np.savetxt("smooth"+str(smoothingfact)+".dat",out)
        if out[0] > 0.5 : continue
        working.append(True)

        # debug line crossing
        #np.diff(out)
        #outsave=np.diff(out)/np.max(np.diff(out))*.9;
        #np.savetxt("pssmoothed"+str(smoothingfact)+"diff.dat",out)

        outhalfright = out[:len(out)/2]
        outhalfleft = np.copy(outhalfright)
        #np.savetxt("newpssmoothed"+str(smoothingfact)+".dat",out1)
        x_max1 = np.where(outhalfright == outhalfright.max())[0][0]    # das ist aber fuer den [ 8 8 20 ] schon der zweite maximalwert
        #print "X_MAXXXXX1:",x_max1


        outhalfright[0:x_max1] = 0
        outhalfleft[x_max1:] = 0

        #np.savetxt("newpssmoothed"+str(smoothingfact)+".dat",outhalfright)
        out1diffright = np.diff(outhalfright)
        out1diffright[x_max1-1] = 0
        #np.savetxt("diff"+str(smoothingfact)+".dat",out1diffright)
        out1diffmaxright = np.where(out1diffright == out1diffright.max())[0][0]

        #np.savetxt("pssmoothed"+str(smoothingfact)+".dat",outhalfleft)
        out1diffleft = np.diff(outhalfleft)
        out1diffleft[x_max1-1] = 0
        #np.savetxt("pssmootheddiff"+str(smoothingfact)+".dat",out1diffleft)
        out1diffmaxleft = np.where(out1diffleft== out1diffleft.max())[0][0]

        out1diffleft[0:out1diffmaxleft] = 0
        out1diffmaxleftright = np.where(out1diffleft == out1diffleft.min())[0][0]
        out1diffleft[0:out1diffmaxleftright] = 0
        out1diffleft = out1diffleft*10000
        lmmm = np.where(out1diffleft == out1diffleft.min())[0][0]
        out1diffleft[0:lmmm] = -1
        lrrr = np.where(out1diffleft == out1diffleft.max())[0][0]
        out1diffleft[lrrr:] = 1

        lminimum = find_nearest(out1diffleft, 0.,returnindex=True)
        minimumalt.append(lminimum)
        #print "lll;",lminimum," this is important for [ 8 8 20 ]

        #out1diffmaxleftrightright = np.where(out1diffleft== out1diffleft.max())[0][0]
        #out1diffleft[0:out1diffmaxleftrightright] = 1

        #np.savetxt("pssmootheddiff"+str(smoothingfact)+".dat",out1diffleft)


        #print "outdiffmaxright:",out1diffmaxright,"out1diffmaxleft:",out1diffmaxleft

        tmp = np.copy(outhalfright)
        tmp[0:out1diffmaxright] = 0
        x_max2right = np.where(tmp == tmp.max())[0][0]

        tmp = np.copy(outhalfright)
        tmp[0:x_max1] = 1
        #np.savetxt("outoutka"+str(smoothingfact)+".dat",tmp)
        tmp[x_max2right:] = 1
        #np.savetxt("outoutkb"+str(smoothingfact)+".dat",tmp)
        x_minimumright = np.where(tmp == tmp.min())[0][0]   # [ 10 10 20 ]: 339748, 345544, nehme dass kleinere

        print smoothingfact,"x_min:", x_minimumright,"l_min:", lminimum
        if x_minimumright < x_max1:
            x_minimumright = len(ps)
        print smoothingfact,"x_max1:",x_max1,out1diffmaxright,"x_max2right:",x_max2right,round(float(x_max2right)/float(x_max1),2),"x_minimumright:",x_minimumright

        minimum.append(x_minimumright)
        #print "fuer den transversalen brauchen wir den hoeheren (rechten) gaussian."
        #print "ps muss entweder auf 0 gesetzt werden beim minimum, oder besser, man kopiert den hintergrund bis zum minimum."
        #print " bis zu 1/3 von x_max1 ist alles nur hintergrund, dass kann kopiert werden"
        #template_len = len(template)
        x_max1list.append(x_max1)
        x_max2list.append(x_max2right)

    print "minialt:",minimumalt
    minimum = np.array(minimum).min()
    minimumalt = np.array(minimumalt).min()
    x_max1 = np.array(x_max1list).min()
    x_max2 = np.array(x_max2list).min()
    if x_max1 == x_max2: minimum = minimumalt
    print
    print "MINIMUM:",minimum
    print "X_MAX1 :",x_max1
    print "X_MAX2 :",x_max2
    print

    #minimum = 339748
    #x_max1 = 258439
    template = np.copy(ps[:x_max1/3+2])
    psout = np.copy(ps)
    #np.savetxt("psin",ps)
    print "x_max1:",x_max1
    print "minimum:",minimum
    for x in np.arange(100):
        ab  = x * x_max1/3
        bis = x * x_max1/3 + x_max1/3
        ab_  = len(psout)-ab    # spiegelverkehrt
        bis_ = len(psout)-bis   # spiegelverkehrt
        if bis > minimum:
            bis = minimum
            bis_ = len(psout)-minimum
        print "ia:",x,"ab:",ab,"bis:",bis,"checks:",psout[ab:bis].shape,template.shape,"lenK:|",len(psout[ab:bis])
        psout[ab:bis+1] = template[:len(psout[ab:bis+1])]   # das minus x ist wichtig weil: np.arange(10)[0:4] == array([0, 1, 2, 3]) und np.arange(10)[5:10] ==  array([5, 6, 7, 8, 9])
        print "ib:",x,"ab:",bis_,"bis:",ab_,"lenK:",len(psout[bis_-1:ab_-1]),len(template[:len(psout[ab:bis])][::-1]),len(template[:len(psout[bis_:ab_])][::-1])
        psout[bis_-1:ab_] = template[:len(psout[ab:bis+1])][::-1]
        print "done"
        if bis >= minimum:
            break
    #np.savetxt("psout",psout)
    return psout

def get_goodsmoothing_for_ps(ps,dt,mdsteps,allowed_func = 0, allowed_der = False,args=False,stringadd=''):
    ''' if allowed_der = False: this will not be checked
        if allowed_der = 0 or another integer it is made sure that this is fulfilled
        an dieser stelle ist ps schon gesmoothed um fakrot 10 oder 100 oder mehr
    '''
    if type(args.smoothing) != bool:
        out1 = smoothing_power_spectrum(ps,args.smoothing)  # ps is the full power spectrum (left and right peak)
        f,lt,ltmin,ltmax,xv,rl,rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,args.smoothing,args=args,stringadd=stringadd)
        od = 4
        return args.smoothing, round(f,od),round(lt,od),round(ltmin,od),round(ltmax,od),xcutoff,out1   # out1 is the smoothed function

    #if args.verbose:
    #    print printblue("###############################  "+stringadd+"   ######### mdsteps: "+str(mdsteps))
    if type(allowed_func) != int:
        sys.exit("ERROR:allowed_func has to be an integer")
    if type(allowed_der) != int:
        allowed_der = 9**40
        #print allowed_der

    #######################################################
    # determine smoothingfactors
    #######################################################
    # determine smoothingfactors
    min_smoothing = 1./len(ps)
    #print "lenps:",len(ps)
    order_magnitude = -1.*(int(math.log10(min_smoothing))-1)
    #print order_magnitude
    order_magnitude = np.arange(1,order_magnitude+1)[::-1]
    #print "ord:",order_magnitude
    for ido,i in enumerate(order_magnitude):
        #print "------>",ido,i,10**(-i)
        order_magnitude[ido] = 10**(-i)
    loop_smoothing = order_magnitude
    #print "loop_smoothing:",loop_smoothing

    #if args.verbose: print "STARTING loop 0 -> (logarithmic loop) 1e-07 .. 1e-01" #,loop_smoothing
    bestworking=False
    #np.savetxt("ps.dat",ps)
    new = True
    if new:
        tmparray = np.zeros((len(loop_smoothing),8))
    #print tmparray
    #print
    for idxtmp,smoothingfact in enumerate(loop_smoothing): # do them all and only afterwards take the broadest which is working, otherwise a local minimum might screw up the result
        #print idxtmp,"smooth:",smoothingfact,"from:",loop_smoothing
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps and out1 are the full power spectrum (left and right peak)
        #print "lenout1:",len(out1)
        #if args.write_full_ps:
        #    np.savetxt("out_loop1_"+str(smoothingfact),out1)


        # debug line crossing
        #np.savetxt("pssmoothed"+str(smoothingfact)+".dat",out1)
        #np.diff(out1)
        #outsave=np.diff(out1)/np.max(np.diff(out1))*.9;
        #np.savetxt("pssmoothed"+str(smoothingfact)+"diff.dat",outsave)
        #np.savetxt("kkb.dat",out1)
        f,lt,ltmin,ltmax,xv,rl,rrll,xcutoff = get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)
        if not new:
            if rl <= allowed_func and rrll <= allowed_der:
                bestworking = smoothingfact
                break
        if new:
            tmparray[idxtmp] = f ,lt,ltmin,ltmax,xv,rl,rrll,xcutoff
            #print tmparray
            #print
    if new:
        printtmparray = False
        if printtmparray:
            print tmparray
        rlall = tmparray[:,5]
        rrllall = tmparray[:,6]
        bestworking = False
        for idx in np.arange(len(loop_smoothing))[::-1]:
            if printtmparray:
                print idx,loop_smoothing[idx],rlall[idx],rrllall[idx],bestworking
            if rlall[idx] <= allowed_func and rrllall[idx] <= allowed_der:
                bestworking = loop_smoothing[idx]
                smoothingfact = bestworking
            else:
                break

    #print "(mdsteps)",mdsteps,"bestworking:",bestworking
    #for idxtmp,smoothingfact in enumerate(loop_smoothing[::-1]):
    #    print "idxtmp:",idxtmp,smoothingfact

    #sys.exit()
    #print "rlall",rlall
    #print "rrllall",rrllall
    #sys.exit()
    loop_smoothing = smoothingfact/10.*np.array([1,2,3,4,5,6,7,8,9,10])  # redo 10 for easier coding, 10 will definitively work.
    loop_smoothing_next = smoothingfact/100.*np.array([1,2,3,4,5,6,7,8,9,10])  # redo 10 for easier coding, 10 will definitively work.
    if type(bestworking) == bool:
        print "no bestworking found! maybe band crossing??? check this!!!!!! in this case you would have 2 gaussians and you'd need to get the linewidths of both"
        bestworking = 1.0
    #if args.verbose: print "DONE loop 1 ->","best:",bestworking,"try:",loop_smoothing

    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps is the full power spectrum (left and right peak)
        f,lt,ltmin,ltmax,xv,rl,rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    if args.verbose > 2:
        print "-----------------------------------------------------------------------------------"
    #print "loop_smoothing:",loop_smoothing
    delta = loop_smoothing[0]/10.
    #print "delta:",delta
    #if args.verbose: print "best:",bestworking
    #print np.arange(1,11)[::-1]
    loop_smoothing = bestworking - delta * np.arange(0,11)[::-1]  # include 0 so bestworking will definitively work (easier code)
    loop_smoothing = np.trim_zeros(loop_smoothing)
    #order_magnitude = -1.*(int(math.log10(min_smoothing))-1)
    loop_smoothing = loop_smoothing[loop_smoothing > 0]
    #if args.verbose: print "DONE loop 2 ->","best:",bestworking,"try:",loop_smoothing
    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps is the full power spectrum (left and right peak)
        f,lt,ltmin,ltmax,xv,rl, rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    if args.verbose > 2:
        print "-----------------------------------------------------------------------------------"
    #if args.verbose:
    #    print "DONE loop 3 ->",smoothingfact,"best:",bestworking

    #    print "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
    #    print "bestworking:",bestworking
    sigma = bestworking
    out1 = smoothing_power_spectrum(ps,sigma)  # ps is the full power spectrum (left and right peak)
    fqgood, ltgood, ltmin, ltmax, xv, rl, rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,sigma,args=args,stringadd=stringadd,verbose=True,color='red')
    #print "wriet out1",mdsteps
    od = 4
    return sigma, round(fqgood,od),round(ltgood,od),round(ltmin,od),round(ltmax,od),xcutoff,out1   # out1 is the smoothed function

def powerspectrum_to_powerspectrum_2d_sparse(ps,smooth=True,mindestens=10000):
    ''' print "len of ps:",len(ps)
        print " immer so fitten dass mindestens 5000 werte drin bleiben "
        print " immer so fitten dass das ps (out) eine gerade anzahl an eintraegen hat -> das macht aber mysmooth nicht!"
    '''

    def get_teiler_daten(ps,mindestens):
        lenps = len(ps)
        teiler = 1
        for i in np.arange(16)+1:
            tryteiler = 1*10**i
            #print i,teiler,lenps/teiler,(lenps/teiler)>=5000
            if ((lenps/tryteiler)>=mindestens) == True:
                #print i,tryteiler,"Ja"
                teiler = tryteiler
            else:
                #print i,tryteiler,"NEIN"
                break
        #print "teiler:",teiler
        N = teiler
        return N

    #print "LEN:",len(ps)
    N = get_teiler_daten(ps,mindestens)
    #print "N:",N
    x = np.arange(len(ps))[::N]+N/2
    #print "x:",x[:3],x[-3:],len(x)

    y = np.mean(ps.reshape(-1, N), 1)
    if smooth:
        print "smooth is TRUE:",N # N = z.b. 100
        y = mysmooth(y,N+1)  # mysmooth(y,N*3+1)  ---> dieser ist also immer ungerade!
    #y = y/y.max()*0.9
    #print "x.sh:",x.shape
    #print "y.sh:",y.shape
    xy=np.transpose(np.array([x,y]))
    ps=y
    #print "xy:",xy[-3:]

    #ps=ps/ps[:len(ps)/2].max()*.9; # full with peaks left and right
    #y=xy[:,1];y=y/y[:len(y)/2].max()*.9;xy[:,1]=y;
    #return np.transpose(np.array([x,y])),x,y
    return xy,x,ps

def check_qpoints_for_crossing(qpoint, N, structure):
    '''
    just gives true for fcc t1 for qpoint [[ 10 10 20 ], [ 9 9 20 ], [ 8 8 20 ]]
    '''
    #print "iiiiiiiiiii",qpoint,N
    returnvalue = False
    #crossingfrom=float(8)/float(10) # there seems to be crossing in DFT GGA 300K [3 3 8]
    #crossingfrom=float(7.5)/float(10)
    crossingfrom=1./np.sqrt(2)
    if qpoint[0]==qpoint[1] and qpoint[2]==2*N:  # bedingutn t1
        if structure == 'fcc': # bedingung fcc
            if float(qpoint[0])/float(N)>=crossingfrom:   # nur [ {8,9,10} {8,9,10} 20]
                returnvalue = True
    elif qpoint[0]==qpoint[1] and qpoint[2]==-2*N:  # bedingutn t1
        if structure == 'fcc': # bedingung fcc
            if float(qpoint[0])/float(N)>=crossingfrom:   # nur [ {8,9,10} {8,9,10} 20]
                returnvalue = True
    else:
        returnvalue = False
    return returnvalue

def get_correct_input_variables(args):
    if args.structure == "fcc": args.usestruct = 4;
    if args.structure == "bcc": args.usestruct = 2;
    if os.path.isfile("struct.bcc"):
        args.structure = "bcc"
        args.usestruct = 2
    #print "ar:",args.lammps_infile
    #print glob.glob("infile.*")
    ##if os.path.isfile(glob.glob("infile.*")):
    ##    args.infile = glob.glob("infile.*")[0]
    ##    print 'ka',args.infile

    N = args.supercell
    alat = args.alat
    dt = args.dt

    # in case of lammps job
    get_inputvariables_from_calculation(positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in", inN=args.supercell, inalat=args.alat, indt=args.dt, inqvec=args.qvec, verbose=args.verbose,args=args)
    if args.supercell == args.alat == args.dt == False:
        get_inputvariables_from_vasp_job()


    print "structure     :",args.structure
    print "N (supercell) :",args.supercell
    print "alat          :",args.alat
    print "dt [fs]       :",args.dt
    print "qvec          :",args.qvec
    if type(args.qvec) == bool:
        sys.exit("ERROR: specify qvec!")
    qpoints_all = get_all_qpoints(args.qvec,args,verbose=True)


    if args.create_lammps_inputfile:
        print "commenot out writing of positions"
        os.popen("sed -i 's|^dump dump1 all xyz|#dump dump1 all xyz|' in_file_dynamics.in").read().rstrip()
    return qpoints_all

######################################################################
# in case of lammps job
######################################################################
def check_for_lammps_job():
    ''' if lammps inputfile is found, lammps job is evaluated '''
    print bcolors.FAIL + "check_from_space_fft" +bcolors.ENDC
    check_for_lammps_infile = glob.glob("in_file_dynamics.in")
    #print "check_for_lammps_infile:",check_for_lammps_infile
    if not len(check_for_lammps_infile) == 1:
        print "This does not look like a lammsp run!"
        return
    else:
        print " -> 1 I assume this is a lammps job since I found the in_file_dynamics.in file"
        if not check_from_space_fft(N,args):
            check_from_xaa_folder(N,args)
    print bcolors.FAIL + "check_from_space_fft DONE" +bcolors.ENDC
    return

def check_for_all_ps_and_create_if_necessary(N,args):
    ''' this currently creates only the ps files and and now also checks if those are
    necessary to create '''
    print bcolors.OKGREEN + "check_for_all_ps_and_create_if_necessary" +bcolors.ENDC
    ####################################################################
    # check for ps xxx files
    ####################################################################
    qpoints_all = get_all_qpoints(args.qvec,args)
    #for i in qpoints_all:
    #    print "hier:",i
    all_ps_files_there = True
    for i in qpoints_all:
        ps_file = glob.glob("ps"+qpointstring(i)+"_*.dat.lifetimesgood.dat")
        #print ps_file,len(ps_file)
        if len(ps_file) < 1:
            all_ps_files_there = False
    print "all_ps_files_there:",all_ps_files_there

    if not all_ps_files_there:
        print "need to create ps_files! (or at least one!)"
        space_fft_to_powerspectrum_ser_par(qpoints_all, execute_timeinversion=True, args = args)
    return

def check_for_all_space_fft_and_rm_xaa_files_folders(N,args):
    print bcolors.OKGREEN + "check_for_all_space_fft_and_rm_xaa_files_folders" +bcolors.ENDC
    qpoints_all = get_all_qpoints(args.qvec,args)
    remove_all = True
    #print "qpall:",qpoints_all,type(qpoints_all)
    for i in qpoints_all:
        #print "iii;",i,type(i)
        #print qpointstring(i)
        if not os.path.isfile("space_fft_"+qpointstring(i)+".npy"):
            remove_all  = False
            print " -> some (or all) space_fft_files are missing ..."
            return False

    if remove_all == True:
        print " -> seems all space_fft_xxx folder are availalbe! -> after removing xaa.. -> get ps"
        # get all xaa files:
        import shutil
        # DONT REMOVE THOSE USER HAS TO DECIDE IF HE WANTS THIS
        #xaa_files = glob.glob("x[a-z]?")
        #for i in xaa_files:
        #    print "removing",i
        #    os.remove(i)
        #xaa_files = glob.glob("x[a-z]?_wcl")
        #for i in xaa_files:
        #    print "removing",i
        #    os.remove(i)
        #xaa_folder= glob.glob("x[a-z]?_")
        #for i in xaa_folder:
        #    print "removing",i
        #    shutil.rmtree(i)
    return remove_all

def get_space_fft_from_xaa_xab_folder_parallel(N,args):
    '''
    ########################################
    # make space_fft
    ########################################'''
    print bcolors.OKGREEN + "get_space_fft_from_xaa_xab_folder_parallel" +bcolors.ENDC
    qpoints_all = get_all_qpoints(args.qvec,args)
    num_cores = multiprocessing.cpu_count()
    print "num_cores++:",num_cores
    #print "qpstring:",qpstring
    #print "qpstring:",qpoints_all
    #sys.exit() # unten
    jobs = []
    for i in qpoints_all:
        print "get_space_fft_from_xaa_xab_folder for",i
        p = multiprocessing.Process(target=get_space_fft_from_xaa_xab_folder, args=[qpointstring(i)])
        jobs.append(p)
        p.start()
    for job in jobs: # this is just to wait until job is finished
        job.join()
    return True

def check_for_xaa_folder_and_make_everything_else(N,args):
    ''' at this stage we assume we do this since no space_fft files exist. '''
    print bcolors.OKGREEN + "check_for_xaa_folder_and_make_everything_else" +bcolors.ENDC
    xaa_file = glob.glob("x??_wcl")
    #print "len(xaa_file):",len(xaa_file)
    if len(xaa_file) == 1: # all xaa files exist; are there also the xaa_ folder?
        xaa_files = glob.glob("x??")
        print " -> 4 All ("+str(len(xaa_files))+") xaa files exists."
        xaa_folder = glob.glob("x??_")
        if len(xaa_folder) == 0:
            print " -> 5 no xaa_ folder, goint to create those and submit to cluster"
            xaa_filenames,xaa_steps = lammps_split_lammps_trajectory_in_xaa_files(filenamelammps="trj_lammps.out",filenamepos="trj_lammpsnew.out",positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in",linesokformemory = 800000000,qpoints_all=qpoints_all,args=args,verbose=args.verbose)
        elif len(xaa_folder) > 0 and len(xaa_folder) == len(xaa_files):
            print " -> 6 All ("+str(len(xaa_folder))+") xaa_ folder seem to exist -> (in case you want to calculate additional space_fft_xxx files delete xaa_ folder first) check in lastfolder if all sum_all_xxx.dat written"
            print "xaa_file:",xaa_file
            last_folder = xaa_file[0].split("_")[0]+"_"
            print len(qpoints_all)
            print "last_folder:",last_folder
            check_space_ffts = glob.glob(last_folder+"/sum_all_*")
            print len(check_space_ffts)
            if len(qpoints_all) == len(check_space_ffts):
                print " -> 7 sum_all_new_xxx files from c++ skript seem to be done"
                if get_space_fft_from_xaa_xab_folder_parallel(N,args):
                    return check_from_space_fft(N,args)
            elif len(qpoints_all) < len(check_space_ffts):
                print " -> 7 check if for the chosen qpoint all sum_all_new_xxx files done"
                for i in qpoints_all:
                    if os.path.isfile(last_folder+'/sum_all_new_'+qpointstring(i)+'.dat'):
                        pass
                    else:
                        sys.exit("it seems that "+last_folder+'/sum_all_new_'+qpointstring(i)+'.dat does not exist!')
                if get_space_fft_from_xaa_xab_folder_parallel(N,args):
                    return check_from_space_fft(N,args)
            else:
                print "len(qpoints_all):",len(qpoints_all),"len(check_space_ffts):",len(check_space_ffts)
                print " -> it seems that xaa_ folder are already created but the sum_all_new is not yet written, is this job in the que?"
                sys.exit()
    else:
        print "--> 8 grep and split necessary"

    return

def check_from_space_fft(N,args):
    print bcolors.FAIL + "check_from_space_fft" +bcolors.ENDC
    if check_for_all_space_fft_and_rm_xaa_files_folders(N,args):
        check_for_all_ps_and_create_if_necessary(N,args)
        print_and_save_lifetimesgood(args = args,printtoscreen=False)
        return True
    else:
        return False

def check_from_xaa_folder(N,args):
    print bcolors.FAIL + "check_from_xaa_folder" +bcolors.ENDC
    if not check_for_xaa_folder_and_make_everything_else(N,args):
        if os.path.isfile("trj_lammps.out") or os.path.isfile("trj_lammpsnew.out"):
            print " -> 9 it seems there are no xaa files! but thre is a trj_lammps(new).out file(s) -> split trajectory and submit c++ to cluster for every xaa folder"
            xaa_filenames,xaa_steps = lammps_split_lammps_trajectory_in_xaa_files(filenamelammps="trj_lammps.out",filenamepos="trj_lammpsnew.out",positionsfilename = "positions.*",infilefilename = "in_file_dynamics.in",linesokformemory = 800000000,qpoints_all=qpoints_all,args=args,verbose=args.verbose)
            sys.exit("I am leaving here since jobs were hopefully submitted to the cluster") # since script above submits to cluster
        else:
            print " -> 9 no xaa files and no trj_lammps(new).out file(s), returning to skript"
            return



if __name__ == '__main__':
    if args.folder == False:
        folder_all = [os.getcwd()]
    else:
        folder_all = glob.glob(os.getcwd()+"/"+args.folder)
    for i in folder_all:
        print i
    print

    for idxfolder,folder_in in enumerate(folder_all):
        os.chdir(folder_in)
        print "#"*(len(os.getcwd())+21)
        print "os.getcwd()  :",idxfolder+1,"/",len(folder_all),os.getcwd()
        print "#"*(len(os.getcwd())+21)

        qpoints_all = get_correct_input_variables(args)
        N = args.supercell
        args.atoms = args.usestruct * (args.supercell**3)

        if args.exit:
            sys.exit()

        if type(args.lifetimesaverage) != bool:
            get_lifetimessummary_from_several_folders(base="run_*/")
            sys.exit()

        if args.make_power_spectrum == False:
            check_for_lammps_job()

        if args.space_fft:
            qpstring_all = [];tmp=0;
            for idx,i in enumerate(qpoints_all):
                tmp+=equivalent_qs(i).shape[0]
                qpstring_all.append(qpointstring(i))

            get_space_fft_from_xaa_xab_folder(qpstring=qpstring_all)

        ########################################
        # look for filenames
        ########################################
        scale_by_alat_N = [ 'trj_lammpsnew.out', 'trj_lammpsnew.out_noJumps_small', 'pos', 'POSITIONs' ] # lammps to direct coords
        scale_by_N = [ 'dum', 'posmichael' ] # michaels sim_fcc_morse (to convert to direct)

        filename = False
        scale_data = False
        for i in scale_by_alat_N:
            if os.path.isfile(i) == True:
                filename = i; scale_data = float(args.alat)*args.supercell
        for i in scale_by_N:
            if os.path.isfile(i) == True:
                filename = i; scale_data = N
        print 'filename     :',filename
        print "scale_data   :",scale_data

        if args.make_space_fft_from_py or args.make_space_fft_from_lammpslog or args.make_power_spectrum:
            for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
                qpointstr = qpointstring(qpoint)
                #print bcolors.FAIL + "######################################################################################" + bcolors.ENDC
                #print bcolors.FAIL + "OUTER LOOP OVER ALL QPOINTS: "+str(idx)+" out of "+str(len(qpoints_all))+" qpointstring: "+qpointstr + bcolors.ENDC
                #print bcolors.FAIL + "######################################################################################" + bcolors.ENDC
                equivalent_qs(qpoint)
                ##############################################################################
                # CHECK IF space_fft_xxxx.npy exist, if not, create it!
                ##############################################################################
                execute_timeinversion = True
                space_fft_file1 = 'space_fft_'+qpointstr+'.npy'
                space_fft_file2 = 'sum_all_new_'+qpointstr+'.npy'
                space_fft_file3 = 'log_for_sum_all' # created from lammps.log
                space_fft_file4 = 'log.lammps'      # lammps.log with comments

                if os.path.isfile(space_fft_file1):
                    pass
                    #print "load existing ",space_fft_file1
                    #sum_all = np.load(space_fft_file1)
                elif os.path.isfile(space_fft_file2):
                    #print "load existing ",space_fft_file2
                    #sum_all = np.load(space_fft_file2)
                    execute_timeinversion=False;
                elif os.path.isfile(space_fft_file3) == True:
                    print "execute lammps_log_to_sum on",space_fft_file3
                    sum_all = lammps_log_to_sum(filename = space_fft_file3, qpoint=qpoint)
                elif os.path.isfile(space_fft_file4) == True and args.make_space_fft_from_lammpslog:
                    print "execute lammps_log_to_sum on",space_fft_file4
                    sum_all = get_space_fft_from_lammps_log(filename = space_fft_file4)
                elif type(filename) != bool and args.make_space_fft_from_py:
                    get_space_fft_from_positions(\
                        filename,\
                        qpoints_all=qpoints_all,\
                        lines_of_one_mdstep=N**3*args.usestruct,\
                        mdsteps_per_chunk = 4000, \
                        scale_data=scale_data, chunkitornot = False, \
                        args = args)

        if args.make_power_spectrum:
            space_fft_to_powerspectrum_ser_par(qpoints_all, execute_timeinversion=True, args = args)
            print_and_save_lifetimesgood(args = args,printtoscreen=False)


        if args.lifetimes or args.freqs: print_and_save_lifetimesgood(args = args)
