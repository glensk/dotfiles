#!/usr/bin/env python

#                                                                     N  alat q1 q2 q3 dt[fs]
# run ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py fcc 10 4.04 10 0  0  10
# run ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py fcc 10 4.04 l  0  0  40
# run ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py fcc 10 4.04 t  t  0  40
# run ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py fcc 10 4.04 l  l  l  10
# run ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py fcc 10 4.07 l  0  0  40 lifetimesgood meV
# run ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py fcc 10 4.07 t1 t1 0  40 lifetimesgood
# run ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py fcc 10 4.07 t1 t1 0  40 freqsgood
# run ~/Thermodynamics/python_thermodynamics/lammps_pos_to_sum.py fcc 10 4.07 t1 t1 0  40 freqsgood mev
#
# lammps_pos_to_sum.py bcc 6 3.253 4 4 4 50   # dominiques job

# /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2016.09_phonon_lifetimes_3_nach_elternzeit/tutild_2016_12_02_10sc_300K_4.04_2_hoch_5_steps

# 2*10^5 schritte --> ein qpunkt hat 19MB
# 1*10^6 schritte --> ein qpunkt hat 200MB  --> ein ast geht locker (20 qpunkte sind 4GB)


import sys
import os
import math
import numpy as np
import h5py
from itertools import permutations
from pandas import read_csv,set_option
from scipy.fftpack import fft
import timeit
import glob
set_option('max_columns', 0)

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

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
        qpointstr = str(qpoint1)+"_"+str(qpoint2)+"_"+str(qpoint3)
    #print qpointstr
    #sys.exit()
    return qpointstr

def get_all_qpoints(N1,N2,N3,N):
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
    appropratestring = [ 'l', 't', 't1', 't2', 'all']
    #print "type(N1):",N1,type(N1)
    #print "str(N1) :",N1,int(N1)
    #print "type(N2):",N2,type(N2)
    #print "str(N2) :",N2,int(N2)
    #print "type(N3):",N3,type(N3)
    liste = [N1,N2,N3]
    qpoints_all =[]
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
    #for i in liste:
    #    print "i:",i,type(i)
    #print ",,,,"
    #print N1 == 0
    #print N1 == 0.
    #print N3 == 0
    #print N3 == 0.
    #print N3 != 0
    #print N3 != 0.

    #print "str(N3) :",N3,int(N3)
    #print "--------"
    if type(N1) != str and type(N2) != str and type(N3) != str:
        qpoints_all = [[N1,N2,N3]] # we have an particular qpoint
        return qpoints_all

    # in case however we had a string in there
    N = int(N)  # supercellsize
    qpoints_all = []
    all = False
    if N1 == 'all' or N2 == 'all' or N3 == 'all':         # [ N 0 0 ] (L)
            all = True
    # string != 0 gibt True
    #
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 0 0 ] >>>>>>>>>>>>>>>>>>>>>>
    if (N1 == 'l' and N2 == 0. and N3 == 0.) or all == True:         # [ N 0 0 ] (L)
        for i in np.arange(N)+1:
            qpoints_all.append([i,0,0])
    if (N1 == 't' and N2 == 0. and N3 == 0.) or all == True:       # [ N 0 0 ] (T)
        for i in np.arange(N)+1:
            qpoints_all.append([2*N,i,0])
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 0 0 ] >>>>>>>>>>>>>>>>>>>>>>


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 0 ] >>>>>>>>>>>>>>>>>>>>>>
    if (N1 == 'l' and N2 == 'l' and N3 == 0.) or all == True:      # [ N N 0 ] (L)
        for i in np.arange(N)+1:
            qpoints_all.append([i,i,0])
    if (N1 == 't1' and N2 == 't1' and N3 == 0.) or all == True:    # [ N N 0 ] (T1)
        for i in np.arange(N)+1:
            qpoints_all.append([i,i,2*N])
    if (N1 == 't2' and N2 == 't2' and N3 == 0.) or all == True:    # [ N N 0 ] (T2)
        for i in np.arange(N)+1:
            qpoints_all.append([2*N+i,2*N-i,0])
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 0 ] >>>>>>>>>>>>>>>>>>>>>>


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 1 ] >>>>>>>>>>>>>>>>>>>>>>
    if (N1 == 'l' and N2 == 'l' and N3 == 'l') or all == True:     # [ N N 0 ] (L)
        for i in np.arange(N/2)+1:
            qpoints_all.append([i,i,i])
    if N1 == 't' and N2 == 't' and N3 == 't':     # [ N N N ] (T)
        sys.exit("ERROR: not yet understood!")
        #for i in np.arange(N/2)+1:
        #    qpoints_all.append([N/2+i,N/2-i,0])
    #print "len:",len(qpoints_all)
    if len(qpoints_all) == 0:
        sys.exit("ERROR: this path is not known!")
    return qpoints_all

def print_and_save_lifetimesgood(keyword):
    ''' keyword is
            - lifetimesgood l 0 0
            - lifetimesgood t2 t2 t2 mev
            - lifetimesgood l l l meV
            - freqsgood t1 t1 0 meV
    '''
    #print "xx:",qpoints_all
    #for i in qpoints_all:
    #    print qpointstring(i)
    in1 = sys.argv[4]
    in2 = sys.argv[5]
    in3 = sys.argv[6]
    print "----in{1,2,3}    :",in1,in2,in3
    keyword = sys.argv[8]  # {lifetimes,freqs,lifetimesmin}
    print "----keyword      :",keyword
    keyword2 = ""
    #print "----len(sys.argv):",len(sys.argv)
    if len(sys.argv) > 9:
        keyword2 = sys.argv[9]  # {mev,meV}
    print "----keyword2     :",keyword2



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
        files.append(glob.glob("ps"+qpointstring(i)+"*.dat."+keyword+".dat")[0])
    #print files
    #sys.exit()
    for file in files:
        #print file
        f1 = file.split("_")
        f2 = f1[3].split(".dat")[0]
        #print "f2:",f2,float(f2),int(float(f2))
        mdsteps.append(int(float(f2)))  # keep int(float())

    mdstepsmax = np.sort(np.array(list(set(mdsteps))))[-1]
    print mdstepsmax
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
    qmax = 0;
    for file in files:
        q1 = file.split("ps")[1].split("_")[0]
        if int(q1) > qmax:
            qmax = int(q1)

    f = open(str(keyword)+'.dat', 'wb')

    for file in files:
        q1 = file.split("ps")[1].split("_")[0]
        #print file,q1
        q2 = file.split("ps")[1].split("_")[1]
        q3 = file.split("ps")[1].split("_")[2]
        lastline = np.loadtxt(file)
        #print lastline
        # 1 0 0
        print type(in1)
        print type(in2)
        print type(in3)
        if in1 == "t" and in2 == "0" and in3 == "0":
            out = q2;faktor=2.;rd=0;sc=-1;
        if (in1 == "l"  and in2 == "0" and in3 == "0") or (type(in1) == str  and in2 == "0" and in3 == "0"):
            out = q1;faktor=1.;rd=0;sc=-1.;

        # 1 1 0
        if in1 == "t1" and in2 == "t1" and in3 == "0":
            out = q1;faktor=1.;rd=1.;sc=1;
        if in1 == "t2" and in2 == "t2" and in3 == "0":
            out = q1;faktor=3.;rd=3.;sc=1;
        if (in1 == "l" and in2 == "l" and in3 == "0") or (type(in1) == str and in2 == "l" and in3 == "0"):
            out = q1;faktor=1.;rd=1.;sc=1.;


        # 1 1 1
        if in1 == "t" and in2 == "t" and in3 == "t":
            out = q2;faktor=2.;rd=0.;sc=-1.;
        if in1 == "l" and in2 == "l" and in3 == "l":
            out = q1;faktor=1.;rd=0.;sc=-1.;


        roundto=4
        print  str(round(sc*(rd-(float(out)*faktor/float(qmax))),4)),"\t",round(lastline[1]*THz_to_meV,roundto),"\t",round(lastline[2]*THz_to_meV,roundto),"\t",round(lastline[3]*THz_to_meV,roundto),"\t",round(lastline[4],7)
        ka = str(round(sc*(rd-(float(out)*faktor/float(qmax))),4))+"\t"+str(round(lastline[1]*THz_to_meV,roundto))+"\t"+str(round(lastline[2]*THz_to_meV,roundto))+"\t"+str(round(lastline[3]*THz_to_meV,roundto))+"\t"+str(round(lastline[4],7))+"\n"
        f.write(ka)
    f.close()

    sys.exit()

def equivalent_qs(q):
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
    return equivalents #*1j*2*np.pi

def pos_compress_to_sum_fast_read_at_once(filename, qpoints_all = False, scale_data = 1., N=False, usestruct=False):
    ''' dont use engine='python' in read_csv which is 6 times slower,
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
    #df = read_csv(filename, sep=' ', header=None,engine='python') # 22.0 sec

    import time
    start1 = time.time()
    start = time.time()
        #reader = read_csv(filename, sep=' ', header=None,chunksize=lines_of_one_mdstep*mdsteps_per_chunk,engine='python')
    #df = read_csv(filename, sep=' ', header=None,engine='c',) # 3.7 sec
    df = read_csv(filename, sep=' ', header=None,engine='c',chunksize=4000*2) # 3.7 sec
    for chunk,pos in enumerate(df):
        print "chunk:",chunk
        xx = pos/float(scale_data)
        print "pos  :",xx
        print xx.as_matrix()
        sys.exit()
    #df = read_csv(filename, sep=' ', header=None,engine='c', dtype={'0': np.float32, '1': np.float32, '2': np.float32}) # 3.7 sec
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
    dfqs = np.exp(np.sum(df*qs_all[0],axis=1)).reshape(-1, N**3*usestruct).sum(axis=1) # the reshape is quick, it is probably the np.exp which takes time
    print sys.getsizeof(dfqs)*1e-09,"GB"
    print "------>",dfqs.dtype
    #aaa[4000:8000].sum()
    #aaa[8000:12000].sum()

    end = time.time()
    print "#2) np.exp(np.sum()):",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    start = time.time()
    #dfqs = dfqs.reshape(-1, N**3*usestruct).sum(axis=1)
    end = time.time()
    print "#3) reshape.sum:",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    #aaa = np.sum(np.exp(np.sum(xxx*qs_all[0],axis=1)))
    print
    print dfqs[:10]
    print dfqs.shape
    end1 = time.time()
    print "#5) total:",bcolors.HEADER+  str((end1 - start1)) +" sec"+ bcolors.ENDC
    return dfqs

# read everything into memory at once
def read_hd5f(filename='dump_h5md.h5',qpoint=[1,0,0],scale_data = False):
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

# this should be in c++ to make it faster
def pos_compress_to_sum(filename, qpoints_all = False, lines_of_one_mdstep = False, mdsteps_per_chunk = 1, scale_data = 1,chunkitornot = False):
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
    qpoints_all=np.array([[10,0,0],[9,0,0]])
    if os.path.isfile(filename) != True:
        print "filename:",filename
        sys.exit("ERROR: filenmae does not exist")
    print "filename:",filename
    print "scale_data:",scale_data,type(scale_data)
    if type(lines_of_one_mdstep) == bool:
        sys.exit("ERROR: lines_of_one_mdstep needs to be a number!")


    ############################################################
    # start doint this for al list of qpoints
    # this would only work only when 2*10**8/amountqpoint (say 20) = 1*10**7 is the max number of timesteps
    ############################################################
    qs_all = dict()
    sum_all = dict()
    #for idx,qpoint in enumerate(qpoints_all):
        #doubles[x] = x * 2
        #qs_all[idx] = equivalent_qs(qpoint)
        #print "define sum_all np.zeros: in case this raises error you need more memory or decrease steps!"
        # - allocation of sum_all depends on memory available (np.zeros(...10^12))
        # - in principle this could be sepereated for every q point but would mean that one
        #   would have to read the file multiple times. (or write the sum to disk once more
        #   steps which could also be a good solution)
        #sum_all = np.zeros((len(qs_all),2*10**8),dtype=np.complex_)  # 4*10^8 seems to be max for cluster
        #sum_all[idx] = np.zeros((len(qs_all[idx]),1*10**7),dtype=np.complex_)  # 4*10^8 seems to be max for cluster; mit 1*10^7 kann man schon 40 qpunkte machen;
        #sum_all[idx] = np.zeros((len(qs_all[idx]),1*10**7),dtype=np.complex_)  # 4*10^8 seems to be max for cluster; mit 1*10^7 kann man schon 40 qpunkte machen;
        #sum_all[idx] = np.zeros((len(qs_all[idx]),1*10**5),dtype=np.complex_)  # 4*10^8 seems to be max for cluster; mit 1*10^7 kann man schon 40 qpunkte machen;
        #print "define sum_all np.zeros DONE"
        #print "lines_of_one_mdstep*mdsteps_per_chunk:",lines_of_one_mdstep*mdsteps_per_chunk
    #print "--------------------------------------------------------"
    #print qs_all[2]
    #sys.exit()



    ############################################################
    # start reading in the positions file
    ############################################################
    import time
    start1 = time.time()
    start = time.time()
    if filename == 'posmichael':
        reader = read_csv(filename, sep=' ', header=None,chunksize=lines_of_one_mdstep*mdsteps_per_chunk,engine='python',usecols=[1,2,3],skiprows=5)
    elif filename == 'dum':
        reader = read_csv(filename, sep=' ', header=None,chunksize=lines_of_one_mdstep*mdsteps_per_chunk,engine='python',usecols=[1,2,3])
    elif filename == 'pos' or filename == 'trj_lammpsnew.out' or filename == 'trj_lammpsnew.out_noJumps_small':
        print "start reader using c engine (nope curr python)!"
        #reader = read_csv(filename, sep=' ', header=None,chunksize=lines_of_one_mdstep*mdsteps_per_chunk,engine='python')
        #reader = read_csv(filename, sep=' ', header=None,engine='c',chunksize=4000*2) # 3.7 sec
        chunksize=lines_of_one_mdstep*mdsteps_per_chunk
        #chunksize=4000*2
        #chunksize=4000*10
        #chunksize=4000*100
        #chunksize=lines_of_one_mdstep*mdsteps_per_chunk
        if chunksize < 1:
            sys.exit("ERROR: chunksize needs to be greater 0!")
        print "chunksize:",chunksize
        print "lines_of_one_mdstep:",lines_of_one_mdstep
        print "mdsteps_per_chunk:",mdsteps_per_chunk
        #reader = read_csv(filename, sep=' ', header=None,engine='c',chunksize=chunksize) # 3.7 sec
        reader = read_csv(filename, sep=' ', header=None,engine='python',chunksize=chunksize) # MUCH SLOWER
        print "start reader done"
    else:
        sys.exit('ERROR: dont know this filename:'+str(filename))
    #reader = read_csv(filename, sep=' ', header=None,chunksize=lines_of_one_mdstep*mdsteps_per_chunk)

    #allavail = glob.glob("sum_all_new_*__"+str(mdsteps_per_chunk)+".npy")
    #for i in allavail:
    maxavail = np.zeros(len(qpoints_all))
    print "qpoints_all:",qpoints_all
    for idx,qpoint in enumerate(qpoints_all):
        qpointstr = qpointstring(qpoint)
        print "qpoint   :",qpoint
        print "qpointstr:",qpointstr,"mdsteps_per_chunk:",str(mdsteps_per_chunk)
        qpointfiles = glob.glob("sum_all_new_"+qpointstr+"__*__"+str(mdsteps_per_chunk)+".npy")

        maxavailqpoint = 0
        for j in qpointfiles: # for particular qpoint
            check = int(j.split("sum_all_new_"+qpointstr+"__")[1].split("__"+str(mdsteps_per_chunk)+".npy")[0])
            if check > maxavailqpoint:
                maxavailqpoint=check
        maxavail[idx] = maxavailqpoint
        #print qpointstr, maxavailqpoint,maxavail,int(maxavail.min())
    maxavail = int(maxavail.min())
    print "maxavail:",maxavail
    if maxavail == 0:
        maxavail = -1;  # just so that nothing goes wrong on step 0

    print "schleife"
    nextcheck = mdsteps_per_chunk
    for chunk,pos in enumerate(reader):
        if chunk < maxavail:
            if chunk == mdsteps_per_chunk or chunk == nextcheck:
                nextcheck = chunk + mdsteps_per_chunk
            continue
        #print "schleife in"
        #print "pos:",pos,"chunk:",chunk
        #print "scale_data:",scale_data,type(scale_data)
        #print "---> 1:",np.abs(pos).min().min()
        #print "---> 1:",np.abs(pos).max().max()
        #print "1"
        #xx1 = pos /float(scale_data)
        #print "2"
        #print xx1
        #sys.exit()
        #xx2 = xx1.as_matrix()
        #print "3"
        #print pos
        #print "done"
        xxx = pos.as_matrix()/float(scale_data)  # keep the float even if you have already a float
        #jprint "4"
        #sys.exit()
        #print "---> 2:",np.abs(xxx).min()
        #print "---> 2:",np.abs(xxx).max()
        #print xxx
        #[[  9.99333895e-03   9.89389675e-01   9.87375024e-01]
        # [  2.67247954e-01   1.13737185e-03   2.62146608e-01]
        # [  4.96956298e-01   1.61162300e-03   9.88854562e-01]
        # [  7.43035943e-01   9.90735492e-01   2.38639338e-01]
        # [  9.93566017e-01   2.39245251e-01   2.57573428e-01]
        # [  2.54780154e-01   2.59485206e-01   9.86915382e-01]
        # [  4.96713240e-01   2.37665810e-01   2.60195589e-01]
        # [  7.46941112e-01   2.38187285e-01   1.52050957e-02]
        # [  9.88076796e-01   4.90893524e-01   9.98361337e-01]
        # [  2.53683834e-01   5.14058681e-01   2.37721229e-01]
        # [  5.14139485e-01   5.09907166e-01   1.20846471e-02]
        # [  7.57383239e-01   5.04911866e-01   2.40059339e-01]
        # [  9.97682816e-01   7.51552287e-01   2.65865167e-01]
        # [  2.39897538e-01   7.56097814e-01   9.53609055e-03]
        # [  5.07634308e-01   7.40612759e-01   2.62847510e-01]
        # [  7.46485580e-01   7.52757799e-01   9.86506089e-01]
        # [  1.74690338e-02   9.98332706e-01   4.97840619e-01]
        # [  2.44518059e-01   9.99013407e-01   7.36423367e-01]
        # [  4.90689042e-01   3.26808985e-03   4.88098507e-01]
        # [  7.58332532e-01   9.94710249e-01   7.41830745e-01]
        # [  3.62759200e-04   2.57165840e-01   7.48354678e-01]
        # [  2.58953031e-01   2.58403959e-01   5.09219127e-01]
        # [  4.97115888e-01   2.47647553e-01   7.48960263e-01]
        # [  7.44870575e-01   2.58012054e-01   5.02703635e-01]
        # [  9.95049745e-01   5.02756221e-01   5.05146712e-01]
        # [  2.40317573e-01   5.12110417e-01   7.57837057e-01]
        # [  5.05100850e-01   4.93188253e-01   4.97974579e-01]
        # [  7.38817310e-01   4.87942121e-01   7.59930599e-01]
        # [  9.91432463e-01   7.52273852e-01   7.40093738e-01]
        # [  2.58558923e-01   7.41484113e-01   5.10353905e-01]
        # [  4.97142578e-01   7.50702173e-01   7.46484653e-01]
        # [  7.47051984e-01   7.64739384e-01   4.98861379e-01]]
        #sys.exit()
        #print xxx.shape
        xxx = np.reshape(xxx,((-1,N**3*usestruct,3)))  # for last step it could be not 3 but 2 oder 1 (32,3)
        #[[[  9.99333895e-03   9.89389675e-01   9.87375024e-01]
        #  [  2.67247954e-01   1.13737185e-03   2.62146608e-01]
        #  [  4.96956298e-01   1.61162300e-03   9.88854562e-01]
        #  [  7.43035943e-01   9.90735492e-01   2.38639338e-01]
        #  [  9.93566017e-01   2.39245251e-01   2.57573428e-01]
        #  [  2.54780154e-01   2.59485206e-01   9.86915382e-01]
        #  [  4.96713240e-01   2.37665810e-01   2.60195589e-01]
        #  [  7.46941112e-01   2.38187285e-01   1.52050957e-02]
        #  [  9.88076796e-01   4.90893524e-01   9.98361337e-01]
        #  [  2.53683834e-01   5.14058681e-01   2.37721229e-01]
        #  [  5.14139485e-01   5.09907166e-01   1.20846471e-02]
        #  [  7.57383239e-01   5.04911866e-01   2.40059339e-01]
        #  [  9.97682816e-01   7.51552287e-01   2.65865167e-01]
        #  [  2.39897538e-01   7.56097814e-01   9.53609055e-03]
        #  [  5.07634308e-01   7.40612759e-01   2.62847510e-01]
        #  [  7.46485580e-01   7.52757799e-01   9.86506089e-01]
        #  [  1.74690338e-02   9.98332706e-01   4.97840619e-01]
        #  [  2.44518059e-01   9.99013407e-01   7.36423367e-01]
        #  [  4.90689042e-01   3.26808985e-03   4.88098507e-01]
        #  [  7.58332532e-01   9.94710249e-01   7.41830745e-01]
        #  [  3.62759200e-04   2.57165840e-01   7.48354678e-01]
        #  [  2.58953031e-01   2.58403959e-01   5.09219127e-01]
        #  [  4.97115888e-01   2.47647553e-01   7.48960263e-01]
        #  [  7.44870575e-01   2.58012054e-01   5.02703635e-01]
        #  [  9.95049745e-01   5.02756221e-01   5.05146712e-01]
        #  [  2.40317573e-01   5.12110417e-01   7.57837057e-01]
        #  [  5.05100850e-01   4.93188253e-01   4.97974579e-01]
        #  [  7.38817310e-01   4.87942121e-01   7.59930599e-01]
        #  [  9.91432463e-01   7.52273852e-01   7.40093738e-01]
        #  [  2.58558923e-01   7.41484113e-01   5.10353905e-01]
        #  [  4.97142578e-01   7.50702173e-01   7.46484653e-01]
        #  [  7.47051984e-01   7.64739384e-01   4.98861379e-01]]]
        #xxx = np.reshape(xxx,((mdsteps_per_chunk,N**3*4,3)))  # for last step it could be not 3 but 2 oder 1 (32,3)
        #print "xxx.shape:",xxx.shape[0]
        #print "xxx.shpae",xxx.shape
        #sum[mdstep] = np.sum(np.exp(np.sum(xxx*qs,axis=1)))  # one value per md step
        mdstepbegin = chunk*mdsteps_per_chunk
        mdstepend = chunk*mdsteps_per_chunk+xxx.shape[0]
        #mdstepend = chunk*mdsteps_per_chunk+mdsteps_per_chunk
        #sys.stdout.write('\r'+str(mdstepend))
        #print "---------------------------------------------------------------------------------------------------------------------------------------------------------------------->"
        print "chunk (=index):",chunk,"step:",mdstepbegin,"-",mdstepend,"nextcheck:",nextcheck
        #print "---------------------------------------------------------------------------------------------------------------------------------------------------------------------->"
        #if mdstepend == 3:
        #    sys.exit()
        #if mdstepend in np.arange(1,10001)*10**5:  # alle 10^5 schritte bis 10^9 schritte
        #    print savestuff

        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
            qs_all[idx] = equivalent_qs(qpoint)
            #sum_all[idx] = np.zeros((len(qs_all[idx]),1*10**5),dtype=np.complex_)  # 4*10^8 seems to be max for cluster; mit 1*10^7 kann man schon 40 qpunkte machen;
            sum_all[idx] = np.zeros((len(qs_all[idx]),mdstepend-mdstepbegin),dtype=np.complex_)  # 4*10^8 seems to be max for cluster; mit 1*10^7 kann man schon 40 qpunkte machen;
            for ind_qs,qs in enumerate(qs_all[idx]):  # fuer alle symmetrieequivalenten qpoints
                #print "ind_qs,qs",ind_qs,qs,mdstepbegin,mdstepend

                #print "--> 1:",qs_all[ind_qs]
                # qvector = [2 2 2]--> 1: [ 0.+12.56637061j  0.+12.56637061j  0.+12.56637061j]

                #aa = xxx*qs_all[ind_qs]
                #[[[ 0. +1.25580001e-01j  0. +1.24330373e+01j  0. +1.24077205e+01j]
                #  [ 0. +3.35833683e+00j  0. +1.42926362e-02j  0. +3.29423144e+00j]
                #  [ 0. +6.24493702e+00j  0. +2.02522519e-02j  0. +1.24263129e+01j]
                #  [ 0. +9.33726504e+00j  0. +1.24499494e+01j  0. +2.99883037e+00j]
                #  [ 0. +1.24855188e+01j  0. +3.00644450e+00j  0. +3.23676315e+00j]
                #  [ 0. +3.20166184e+00j  0. +3.26078726e+00j  0. +1.24019445e+01j]
                #  [ 0. +6.24188267e+00j  0. +2.98659665e+00j  0. +3.26971421e+00j]
                #  [ 0. +9.38633884e+00j  0. +2.99314970e+00j  0. +1.91072868e-01j]
                #  [ 0. +1.24165392e+01j  0. +6.16874995e+00j  0. +1.25457786e+01j]
                #  [ 0. +3.18788508e+00j  0. +6.45985191e+00j  0. +2.98729307e+00j]
                #  [ 0. +6.46086731e+00j  0. +6.40768243e+00j  0. +1.51860154e-01j]
                #  [ 0. +9.51755848e+00j  0. +6.34490964e+00j  0. +3.01667463e+00j]
                #  [ 0. +1.25372520e+01j  0. +9.44428457e+00j  0. +3.34096023e+00j]
                #  [ 0. +3.01464137e+00j  0. +9.50140535e+00j  0. +1.19834048e-01j]
                #  [ 0. +6.37912085e+00j  0. +9.30681441e+00j  0. +3.30303923e+00j]
                #  [ 0. +9.38061445e+00j  0. +9.45943348e+00j  0. +1.23968011e+01j]
                #  [ 0. +2.19522354e-01j  0. +1.25454188e+01j  0. +6.25604972e+00j]
                #  [ 0. +3.07270455e+00j  0. +1.25539727e+01j  0. +9.25416896e+00j]
                #  [ 0. +6.16618036e+00j  0. +4.10680283e-02j  0. +6.13362673e+00j]
                #  [ 0. +9.52948765e+00j  0. +1.24998976e+01j  0. +9.32212007e+00j]
                #  [ 0. +4.55856655e-03j  0. +3.23164125e+00j  0. +9.40410223e+00j]
                #  [ 0. +3.25409976e+00j  0. +3.24719991e+00j  0. +6.39903627e+00j]
                #  [ 0. +6.24694248e+00j  0. +3.11203094e+00j  0. +9.41171224e+00j]
                #  [ 0. +9.36031970e+00j  0. +3.24227509e+00j  0. +6.31716019e+00j]
                #  [ 0. +1.25041639e+01j  0. +6.31782100e+00j  0. +6.34786080e+00j]
                #  [ 0. +3.01991969e+00j  0. +6.43536930e+00j  0. +9.52326132e+00j]
                #  [ 0. +6.34728447e+00j  0. +6.19758637e+00j  0. +6.25773311e+00j]
                #  [ 0. +9.28425214e+00j  0. +6.13166153e+00j  0. +9.54956955e+00j]
                #  [ 0. +1.24587078e+01j  0. +9.45335202e+00j  0. +9.30029220e+00j]
                #  [ 0. +3.24914726e+00j  0. +9.31776417e+00j  0. +6.41329632e+00j]
                #  [ 0. +6.24727788e+00j  0. +9.43360172e+00j  0. +9.38060281e+00j]
                #print aa
                #print "--> 2:"
                #bb = np.sum((xxx*qs_all[ind_qs]),axis=2)
                #print "--> 3:", bb
                #--> 3: [[ 0.+24.96633782j  0. +6.66686091j  0.+18.69150219j  0.+24.78604478j
                #   0.+18.72872645j  0.+18.86439356j  0.+12.49819352j  0.+12.57056141j
                #   0.+31.13106774j  0.+12.63503006j  0.+13.0204099j   0.+18.87914274j
                #   0.+25.32249682j  0.+12.63588077j  0.+18.98897449j  0.+31.23684906j
                #   0.+19.02099086j  0.+24.88084623j  0.+12.34087512j  0.+31.35150536j
                #   0.+12.64030205j  0.+12.90033594j  0.+18.77068566j  0.+18.91975498j
                #   0.+25.16984567j  0.+18.9785503j   0.+18.80260396j  0.+24.96548322j
                #   0.+31.212352j    0.+18.98020775j  0.+25.06148242j  0.+25.26660759j]]

                #cc = np.sum(np.exp(np.sum((xxx*qs_all[ind_qs]),axis=2)),axis=1)
                #print "--> 4:",cc
                # --> 4: [ 31.4394438-0.01434658j]

                #print "--> 5:"
                #print "############################ POSITIONS #############################"
                #print "xxx:",xxx
                #print xxx.shape
                #print
                #print "qs_all:",qs_all[ind_qs]
                #print
                #print "a = xxx*qs_all[ind_qs] ------------------------- hat noch ueberall 0 im realteil"
                #a = xxx*qs_all[ind_qs]
                #print a
                #print a.shape # (1, 432, 3) von dominique 432 ist anzahl atome

                #print
                #print "b = np.sum((xxx*qs_all[ind_qs]),axis=2) ------------------------- hat immernoch ueberall 0 im realteil"
                #b = np.sum((xxx*qs_all[ind_qs]),axis=2)
                #print b
                #print b.shape  # (1, 432) von dominique

                #print
                #print "c = np.exp(b)------------------------- hat zahlen im realteil"
                #c= np.exp(np.sum((xxx*qs_all[ind_qs]),axis=2))
                #print c
                #print "c.shape:",c.shape # (1, 432) von dominique
                #print "d-------------------------"
                #d=np.sum(np.exp(np.sum((xxx*qs_all[ind_qs]),axis=2)),axis=1)
                #print d
                #print "d.shape:",d.shape
                #print "--------------------------"
                #print "--------------------------"
                #print "--------------------------"
                #sys.exit()
                #print "iii:",ind_qs
                #if ind_qs == 5:
                #    sys.exit()



                #print np.sum(np.exp(np.sum((xxx*qs_all[ind_qs]),axis=2)),axis=1)
                #[ 31.4394438-0.01434658j]
                #sum_all[idx][ind_qs,mdstepbegin : mdstepend] = np.sum(np.exp(np.sum((xxx*qs_all[idx][ind_qs]),axis=2)),axis=1)
                sum_all[idx][ind_qs] = np.sum(np.exp(np.sum((xxx*qs_all[idx][ind_qs]),axis=2)),axis=1)
            #############################################
            # schleife ueber equivalent qpoints done
            #############################################
        ##################################################
        # schleife ueber all qpoints [100],[200],... done (they are still all in memory)
        ##################################################
        # once all symmetry equivalent qpoints are done:
        # this is done for every qpoint (and every chunk)
        # the chunk goes from 0 to .... as many as chunks ...
        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
            qs_all[idx] = equivalent_qs(qpoint)
            qpointstr = qpointstring(qpoint)
            np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),sum_all[idx])

        # for every 10-100 chunks, put all qpoints together and remove the last ones.
        # e.g. chunk  0-9  put all together to 9
        # e.g. chunk 10-19 put all together to 19
        if chunk == mdsteps_per_chunk or chunk == nextcheck:

            loadall = np.arange(nextcheck-mdsteps_per_chunk,nextcheck+1)
            #print loadall
            for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
                qpointstr = qpointstring(qpoint)
                a = np.load("sum_all_new_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
                for b in loadall[1:]:
                    b = np.load("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                    a = np.concatenate((a, b), axis=1)
                # this will overwrite the last chunkfile
                np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
                for b in loadall[:-1]:  # delete all but last cunkfile
                    os.remove("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")

            nextcheck = chunk + mdsteps_per_chunk
        #if chunk == 10:
        #    sys.exit()
    # get last steps
    if chunk < nextcheck:
        loadall = np.arange(nextcheck-mdsteps_per_chunk,chunk+1)
        print loadall
        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
            qpointstr = qpointstring(qpoint)
            a = np.load("sum_all_new_"+qpointstr+"__"+str(loadall[0])+"__"+str(mdsteps_per_chunk)+".npy")
            for b in loadall[1:]:
                b = np.load("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")
                a = np.concatenate((a, b), axis=1)
            # this will overwrite the last chunkfile
            np.save("sum_all_new_"+qpointstr+"__"+str(chunk)+"__"+str(mdsteps_per_chunk),a)
            for b in loadall[:-1]:  # delete all but last cunkfile
                os.remove("sum_all_new_"+qpointstr+"__"+str(b)+"__"+str(mdsteps_per_chunk)+".npy")


                #print "--> 6: done"
        # an dieser stelle sollte man schreiben der sum_all alle 10^6
        #if mdstepend in np.arange(1,10001)*10**4:  # alle 10^4 schritte bis 10^8 schritte
        #    for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
        #        qpointstr = qpointstring(qpoint)
        #        print "save sum_all[idx] @",mdstepend
        #        np.save("sum_all_new_"+qpointstr,sum_all[idx][:,:mdstepend])
    print "DONE! mdstepend:",mdstepend,"------------------------------------------------"
    end = time.time()
    print "#1) read in data :",bcolors.FAIL + str((end - start)) +" sec"+ bcolors.ENDC
    sys.exit()
    print "ranaming files:"
    for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
        qpointstr = qpointstring(qpoint)
        file = glob.glob("sum_all_new_"+qpointstr+"__*")
        print len(file),file
        if len(file) == 1:
            os.rename(file[0], 'sum_all_new_'+qpointstr+'.npy')

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

def smoothing_power_spectrum(power_spectrum,sigma):
    ''' smoothes a gaussian by sigma '''
    steps = power_spectrum.shape[0]                 # 1000
    nu = np.arange(0,steps)/float(steps)            # [ 0., 0.001, ..., 0.999 ]
    kern = np.exp(-(np.mod(nu+0.5,1)-0.5)**2/(2*sigma**2))  # Gaussian distribution in 1-D
    kern = kern/np.sum(kern)
    from scipy.fftpack import fft,ifft
    ps = power_spectrum_smooth = np.real(ifft(fft(power_spectrum)*fft(kern)))
    ps=ps/np.max(ps)*.9;
    return ps

def get_uncertainty_smoothed_powerspectrum(power_spectrum_smooth,dt,sigma = ""):
    ''' - dt is the sampled timestep in femtoseconds
        - power_spectrum_smooth is the full power spectrum having both peaks (left and right)
        - sigma is only for plotting to screen
    '''
    data = power_spectrum_smooth
    #np.savetxt("ka.dat",data)

    xvalues=data.shape[0]
    #print "xvalues all:",xvalues
    datahalf = data[0:xvalues/2]   # only take one peak
    #np.savetxt("kb.dat",datahalf)
    #print "datahalf:",datahalf
    #print "datahalf.max():",datahalf.max()
    #print "datahalf.max()/2:",datahalf.max()/2,"(this is used in find_nearest)"

    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]


    # x index of absolute maximum
    x_ind_max0 = np.where(datahalf == datahalf.max())[0]  # 963  (9.63 THz)
    x_ind_max = np.where(datahalf == datahalf.max())[0][0]  # 963  (9.63 THz)
    #print x_ind_max0,x_ind_max # [14], 14
    #############################################################
    # dis ist nur schoen bei einem relativ glatten gaussian
    #############################################################
    # x index of left
    y_ind_max_over2_left  = find_nearest(datahalf[0:x_ind_max], datahalf.max()/2.)
    y_ind_max_over3_left  = find_nearest(datahalf[0:x_ind_max], datahalf.max()/3.)
    #print "y_ind_max_over2_left",y_ind_max_over2_left
    x_ind_left_over2 = np.where(datahalf[0:x_ind_max] == y_ind_max_over2_left)[0][0]
    x_ind_left_over3 = np.where(datahalf[0:x_ind_max] == y_ind_max_over3_left)[0][0]

    y_ind_max_over2_right = find_nearest(datahalf[x_ind_max:-1], datahalf.max()/2.)
    y_ind_max_over3_right = find_nearest(datahalf[x_ind_max:-1], datahalf.max()/3.)
    y_ind_max_over1000_right = find_nearest(datahalf[x_ind_max:-1], datahalf.max()/1000.)
    #print "y_ind_max_over2_right",y_ind_max_over2_right
    x_ind_right_over2 = np.where(datahalf == y_ind_max_over2_right)[0][0]
    x_ind_right_over3 = np.where(datahalf == y_ind_max_over3_right)[0][0]
    x_ind_right_over1000 = np.where(datahalf == y_ind_max_over1000_right)[0][0]
    #print "xxx:",x_ind_right_over1000
    #sys.exit()
    #print "x_ind_max0:",x_ind_max0
    #print "x_ind_left_over2",x_ind_left_over2,"y:",y_ind_max_over2_left
    #print "x_ind_right_over2",x_ind_right_over2, "y:",y_ind_max_over2_right


    #datahalf[x_ind_left_over2:x_ind_max0]
    #np.diff(datahalf[x_ind_left_over2:x_ind_max0])*600.
    #print "indizes:",x_ind_left_over2,"<-",x_ind_max0,"->",x_ind_right_over2 ,"from" ,xvalues,"xvalues, DELTA:",x_ind_right-x_ind_left
    l = len(np.where(np.diff(datahalf[x_ind_left_over3:x_ind_max0])<0)[0])
    #print
    #print datahalf[x_ind_max0:x_ind_right]
    #print np.diff(datahalf[x_ind_max0:x_ind_right])
    r = len(np.where(np.diff(datahalf[x_ind_max0:x_ind_right_over3])>0)[0])
    #print "r+l:",r+l
    #sys.exit()

    #############################################################
    # checke die ableitung des peaks
    #############################################################
    ###############################
    # erstmal die linke haelfte
    ###############################
    #ka= datahalf[x_ind_left_over2:x_ind_max0]
    #print
    #kb =np.diff(datahalf[x_ind_left_over2:x_ind_max0])*600.
    #print
    #np.savetxt("ko.dat",datahalf)
    #np.savetxt("ka.dat",ka)
    #np.savetxt("kb.dat",kb)
    #print datahalf
    #jprint x_ind_left_over3
    #print x_ind_left_over3,x_ind_max0,x_ind_right_over3
    ll = lr = 0.
    if (x_ind_max0-x_ind_left_over3)[0] > 1:
        datatake = datahalf[x_ind_left_over3:x_ind_max0]
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
    if (x_ind_right_over3-x_ind_max0)[0] > 1:
        datatake = datahalf[x_ind_max0:x_ind_right_over3]
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
    absmin = np.where(datahalf>datahalf.max()/2.)[0][0]
    absmax = np.where(datahalf>datahalf.max()/2.)[0][-1]


    #mult = dt*10 / float(xvalues)  ## Wrong!  see example 6x6x6sc 6_0_0 qpoint between dt={10,20}
    mult = 1000./(dt*float(xvalues))
    freq = x_ind_max * mult
    lifetime = (x_ind_right_over2 - x_ind_left_over2) * mult
    lifetimemin = (absmax - absmin) * mult
    #vor_freq = x_ind_max/float(xvalues)
    #vor_lifetime = (x_ind_right_over2 - x_ind_left_over2)/float(xvalues)
    #vor_lifetimemin = (absmax - absmin)/float(xvalues)
    ##print "vor_freq:",vor_freq
    ###print "vor_lifetime:",vor_lifetime
    ###print "vor_lifetimemin:",vor_lifetimemin
    #print "Smoothed function:"
    print "Freq:",str(round(freq,3)).ljust(6,'0'),"(THz); lifetime: ",str(round(lifetime,3)).ljust(6,'0'),"(THz)", "||", x_ind_left_over3,x_ind_max0,x_ind_right_over3,"||", "sigma:",str(sigma),"r+l",r+l,"rr,ll",ll+lr+rl+rr
            #"lifetime min:",lifetimemin,\
    #print "indizes:",x_ind_left_over2,"<-",x_ind_max0,"->",x_ind_right_over2 #,"from" ,xvalues,"xvalues"
    return freq,lifetime, lifetimemin,xvalues, r+l,ll+lr+rl+rr, x_ind_right_over1000

def get_goodsmoothing_for_ps(ps,dt,mdsteps,allowed_func = 0, allowed_der = False):
    ''' if allowed_der = False: this will not be checked
        if allowed_der = 0 or another integer it is made sure that this is fulfilled
    '''
    if type(allowed_func) != int:
        sys.exit("ERROR:allowed_func has to be an integer")
    if type(allowed_der) != int:
        allowed_der = 9**40
        #print allowed_der

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

    print "STARTING loop 0 ->",loop_smoothing
    for smoothingfact in loop_smoothing:
        #print "111"
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps and out1 are the full power spectrum (left and right peak)
        #np.savetxt("xxx.dat",ps)
        #np.savetxt("yyy.dat",out1)
        #sys.exit()
        #print "222"
        freq1, lifetime1, lifetimemin1, xv, rl, rrll, xmaxwrite = get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    loop_smoothing = smoothingfact/10.*np.array([1,2,3,4,5,6,7,8,9,10])  # redo 10 for easier coding, 10 will definitively work.
    loop_smoothing_next = smoothingfact/100.*np.array([1,2,3,4,5,6,7,8,9,10])  # redo 10 for easier coding, 10 will definitively work.
    print "DONE loop 1 ->",smoothingfact,"best:",bestworking,"in zehner",loop_smoothing

    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps is the full power spectrum (left and right peak)
        freq1, lifetime1, lifetimemin1, xv, rl, rrll, xmaxwrite = get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    delta = loop_smoothing[0]/10.
    #print "delta:",delta
    print "best:",bestworking
    #print np.arange(1,11)[::-1]
    loop_smoothing = bestworking - delta * np.arange(0,11)[::-1]  # include 0 so bestworking will definitively work (easier code)
    loop_smoothing = np.trim_zeros(loop_smoothing)
    #order_magnitude = -1.*(int(math.log10(min_smoothing))-1)
    print "DONE loop 2 ->",smoothingfact,"best:",bestworking,"check:",loop_smoothing
    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps is the full power spectrum (left and right peak)
        freq1, lifetime1, lifetimemin1, xv, rl, rrll, xmaxwrite = get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    print "DONE loop 3 ->",smoothingfact,"best:",bestworking

    print "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
    print "bestworking:",bestworking
    sigma = bestworking
    out1 = smoothing_power_spectrum(ps,sigma)  # ps is the full power spectrum (left and right peak)
    fqgood, ltgood, ltmin, xv, rl, rrll, xmaxwrite = get_uncertainty_smoothed_powerspectrum(out1,dt,sigma)

    return sigma, fqgood,ltgood,xmaxwrite,out1   # out1 is the smoothed function

def sum_all_to_powerspectrum(sum_all, qpoint, mdstepstocalc = False, dt = False, appendlast=True):
    ''' calculates the power_spectrum for a certain qpoint (from sum_all)
    - sum_all contins for a single q point all (symmetry) equivalent q points
    - output of this skript is the power_spectrum for the corresponding q-point
    - qpoint is just for saving the power_spectrum
    '''
    if dt == False:
        sys.exit("ERROR: need dt!")



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
    steps_tocheck_md_error = int(sum_all.shape[1]/10)
    #print steps_tocheck_md_error,type(steps_tocheck_md_error)
    #print "steps_tocheck_md_error",steps_tocheck_md_error,sum_all.shape[1]
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
        mdstepstocalc_all = np.sort(np.append(mdstepstocalc_all,[sum_all.shape[1]/10*10]))

    #####################################################
    # MDSTEPSTOCALC schoenheitskorrekturen der zahlen
    #####################################################
    #print "mdstepstocalc:",mdstepstocalc
    mdstepstocalc_all = np.sort(np.unique(mdstepstocalc_all[np.where(mdstepstocalc_all <= sum_all.shape[1])[0]]))
    #print "mdstepstocalc_all:",mdstepstocalc_all

    #mdstepstocalc_all = np.array([500000])
    #mdstepstocalc_all = np.array([333330])
    #print "mdstepstocalc:",mdstepstocalc,type(mdstepstocalc)
    if type(mdstepstocalc) == bool:
        mdstepstocalc = mdstepstocalc_all

    #print "mdstepstocalc_________________:",mdstepstocalc,type(mdstepstocalc)
    if type(mdstepstocalc) != bool:
        if type(mdstepstocalc) == list:
            mdstepstocalc = np.array(mdstepstocalc)
            #print "hier 1"
        elif type(mdstepstocalc) == int:
            # check if we want more in mdstepstocalc than we can offer in sum_all:
            if mdstepstocalc > mdstepstocalc_all.max():
                mdstepstocalc = np.array([mdstepstocalc_all.max()])
            else:
                mdstepstocalc = np.array([mdstepstocalc])
            #print "hier 2"
        elif type(mdstepstocalc) == np.ndarray:
            #print "hier 3"
            pass
    else:
        sys.exit("error: unknow type for mdstepstocalc")
        #print "mdstepstocalc_all:",mdstepstocalc_all

    #print "2 mdstepstocalc_all:",mdstepstocalc_all
    #print "2 mdstepstocalc:",mdstepstocalc
    #print mdstepstocalc.min(),-1.*(int(math.log10(mdstepstocalc.min()))-1)
    ordmag = int(math.log10(mdstepstocalc.min()))
    #print mdstepstocalc.min(),10**ordmag

    mdstepstocalc = np.unique(mdstepstocalc/10**ordmag*10**ordmag) # remove elements wich are too close
    mdstepstocalc = mdstepstocalc[np.where(mdstepstocalc > 999)[0]]
    mdstepstocalc_all = mdstepstocalc_all[np.where(mdstepstocalc_all > 999)[0]]
    lifetimestocheck = np.zeros(len(mdstepstocalc_all))
    freqstocheck = np.zeros(len(mdstepstocalc_all))
    print "--> mdstepstocalc:", mdstepstocalc# ,mdstepstocalc[np.where(mdstepstocalc > 999)[0]]
    #sys.exit()
    #####################################################
    # schleife ueber mdsteps
    #####################################################
    for idx,mdsteps in enumerate(mdstepstocalc):
        qpointstr = qpointstring(qpoint)
        filenameout = "ps"+qpointstr+"_"+str(mdsteps)+".dat"
        #print "###################### idx",idx
        print
        print bcolors.FAIL + "idx "+str(idx+1)+" out of "+str(len(mdstepstocalc))+" mdstep:"+str(mdsteps) + bcolors.ENDC
        #ble_all = np.zeros((sum_all.shape[0],mdsteps))
        power_spectrum=np.zeros(mdsteps)
        print
        print os.getcwd()
        print len(sum_all),"mdsteps:",mdsteps," = ",idx+1,"from",len(mdstepstocalc)
        for i in np.arange(len(sum_all)): # laeuft ueber alle symmetrieequivalten qs (6 stueck)
            #print "symmetry equivalent (should all need same time)",i+1,"of",len(sum_all)
            #print "qpoint:",i+1,"of ",len(sum_all),"mdsteps:",mdsteps,"from",mdstepstocalc
            #bla = sum_all[i][::-1][:mdsteps] # turn around data to get rid of equilibration problem
            #print "bla.shape:",bla.shape,i
            #print "."
            #ble=fft(bla[:mdsteps])
            #ble=fft((sum_all[i][::-1][:mdsteps])[:mdsteps])  # put into one command to save (hopefully) memory
            #print ".."
            #print "ble.shape:",ble.shape
            #ble_all[i] = ble
            #power_spectrum+=np.abs(ble)**2.
            #power_spectrum+=np.abs(fft((sum_all[i][::-1][:mdsteps])[:mdsteps]))**2.  # put into one command to save (hopefully) memory

            ######################################################################
            # fuer dominique
            ######################################################################
            #power_spectrum=np.abs(fft((sum_all[i][:mdsteps])[:mdsteps]))**2.  # put into one command to save (hopefully) memory
            #psd=power_spectrum/power_spectrum.max()*.9; # full with peaks left and right
            #np.savetxt(str(i)+"_"+str(mdsteps)+".dat",psd)

            power_spectrum+=np.abs(fft((sum_all[i][:mdsteps])[:mdsteps]))**2.  # put into one command to save (hopefully) memory
            #print ".. 90.000.000 lines are 42%MEM of cmmc"
            #filenameout = "ps"+str(qpoint[0])+str(qpoint[1])+str(qpoint[2])+"_"+str(mdsteps)+"qpoint"+str(i)+".dat"
            #np.savetxt(filenameout,power_spectrum)

        ##########################################################################
        # write out power_spectrum (simple list, only y without x values;
        # length of mdsteps/2 which changes in this loop)
        ##########################################################################
        #power_spectrum[0:10] = 0;  # sometimes very high values , 10 is too much for md with 100 steps
        power_spectrum[0:2] = 0;  # sometimes very high values
        # power_spectrum is just a list of values (only y without x)
        ps=power_spectrum/power_spectrum.max()*.9; # full with peaks left and right
        factor_psout = 2  # 2 is a save bet since power spectrum is doubled
        #np.save("ble_all",ble_all)

        # write out power_spectrum
        if mdsteps == mdstepstocalc.max():
            print "qp:",qpoint,type(qpoint)
            print "write out power_spectrum ...", filenameout
            np.savetxt(filenameout,ps[:len(ps)/factor_psout])  # saves only half the powerspektrum (left peak)


        ##########################################################################
        # write out power_spectrum_smoothed (simple list, only y without x values
        # length of mdsteps/2 which changes in this loop)
        # out1 = smoothing_power_spectrum(ps,smoothingfact)  # out1 is as the powerspectrum have both peaks and has exactly as many as many points
        ##########################################################################
        # what goes in:
        #(ps,dt,  und evtl. mdsteps,
        # at some point also the information: first derivative = 0, second derivative = 0
        #if False:
        if True:
            if mdsteps > 999:
                sigmamin, fqmin,ltmin,xmaxwritemin,minsmoothfunc = get_goodsmoothing_for_ps(ps,dt,mdsteps,allowed_func = 1, allowed_der = False)
                sigmagood, fqgood,ltgood,xmaxwritegood,goodsmoothfunc = get_goodsmoothing_for_ps(ps,dt,mdsteps,allowed_func = 0, allowed_der = 0)
                #print sigmamin,sigmagood,ltmin,ltgood
                #out1out = out1[:len(out1)/factor_psout]
                #############################################################
                # write here the smoothigfunc, further down the lifetimes
                #############################################################
                if mdsteps == mdstepstocalc.max():
                    #print "YYYYYYY:",mdsteps,mdstepstocalc.max(), mdsteps == mdstepstocalc.max()
                    filenameoutsmooth = filenameout+".smooth_good"+str(sigmagood)+".dat"
                    filenameoutmin = filenameout+".smooth_min"+str(sigmamin)+".dat"
                    print "write out power_spectrum_smoothed sigma:",sigmagood,"...",filenameoutsmooth
                    np.savetxt(filenameoutsmooth,goodsmoothfunc[:xmaxwritegood])  # writing takes considerable time for 90.000.000 steps
                    np.savetxt(filenameoutmin,minsmoothfunc[:xmaxwritemin])  # writing takes considerable time for 90.000.000 steps

        lifetimestocheck[idx] = ltgood  # this is only to check error of mdlength
        freqstocheck[idx] = fqgood      # this is only to check error of mdlength
        #print bcolors.FAIL + "----------------------" + bcolors.ENDC
        print 'idx:',idx,'ltgood:',ltgood,"mdsteps:",mdsteps,"ltcheck:",lifetimestocheck




        #@#########################################################################
        #@# write out frequencies and lifetimes (for xmgrace) (currentlyfor every mdstep)
        #@#########################################################################
        #@sigma=0.0005
        #@print "smoothing power_spectrum_smoothed",sigma
        #@out1 = smoothing_power_spectrum(ps,sigma)  # ps is the full power spectrum (left and right peak)
        #@freq1, lifetime1, dummy, xv, rl, rrll, xmaxwrite = get_uncertainty_smoothed_powerspectrum(out1,dt,sigma)
        #@filenameoutsmooth = filenameout+".smooth_"+str(sigma)+".dat"
        #@print "write out power_spectrum_smoothed sigma:",sigma,"...",filenameoutsmooth
        #@np.savetxt(filenameoutsmooth,out1[:xmaxwrite])  # writing takes considerable time for 90.000.000 steps
        #@print "write out power_spectrum_smoothed",sigma,"...done"


        #@### fixed lifetimes to max
        #@sigmafixmax = 0.0011
        #@sigmafixmax = 0.00075
        #@sigma=sigmafixmax  # this should be the maximum of all used lifetimes
        #@print "smoothing power_spectrum_smoothed",sigma
        #@out2 = smoothing_power_spectrum(ps,sigma)
        #@freqfixmax, ltfixmax, dummy, xv, rl, rrll, xmaxwrite = get_uncertainty_smoothed_powerspectrum(out2,dt,sigma)
        #@filenameoutsmooth = filenameout+".smooth_"+str(sigma)+".dat"
        #@print "write out power_spectrum_smoothed sigma:",sigma,"...",filenameoutsmooth
        #@np.savetxt(filenameoutsmooth,out2[:xmaxwrite]) # # writing takes considerable time for 90.000.000 steps
        #@print "write out power_spectrum_smoothed",sigma,"...done"
        #@outlifetimes = np.array([[xv, ltfixmax,ltfixmax/100.,ltfixmax/100.,sigma]])
        #@outfreq = np.array([[xv, freqfixmax,freqfixmax/100.,freqfixmax/100.,sigma]])
        #@np.savetxt(filenameout+".lifetimesfixmax.dat",outlifetimes)
        #@np.savetxt(filenameout+".freqsfixmax.dat",outfreq)

        #@ps = None # Delete variable from memory
        #@out1 = None # Delete variable from memory
        #@out2 = None # Delete variable from memory

        #fq=np.array([freq1, freq2])
        #lt=np.array([lifetime1,lifetime2])

        #outlifetimes = np.array([[xv, lt.mean(), lt.max()-lt.mean(),lt.min()-lt.mean()]])
        #outfreq = np.array([[xv, fq.mean(), fq.max()-fq.mean(),fq.min()-fq.mean()]])

        #np.savetxt(filenameout+".lifetimes.dat",outlifetimes)
        #np.savetxt(filenameout+".freqs.dat",outfreq)

    #############################################
    # get error due to MDlength
    #############################################
    print "lifetimestocheck,sigmagood",lifetimestocheck,sigmagood
    print lifetimestocheck.min()
    print lifetimestocheck.max()
    error_md_lt = abs(lifetimestocheck.max() - lifetimestocheck.min())
    error_md_fq = abs(freqstocheck.max() - freqstocheck.min())
    error_gaussian_lt =  abs(ltgood - ltmin)
    error_gaussian_fq =  abs(fqgood - fqmin)
    error_fq = error_md_fq + error_gaussian_fq
    #print error_gaussian_lt
    #print aerror_gaussian_fq
    #print "mdsteps:",mdsteps
    if mdsteps == mdstepstocalc.max():
        filenameoutsmooth = filenameout+".smooth_good"+str(sigmagood)+".dat"
        #print "&&&&&&write out power_spectrum_smoothed sigma:",sigmagood,"...",filenameoutsmooth
        #print "XXXXXXX:",mdsteps,mdstepstocalc.max(), mdsteps == mdstepstocalc.max()
        outlifetimes = np.array([[mdsteps, ltgood,error_md_lt,error_md_lt+error_gaussian_lt,sigmagood]])
        outfreq = np.array([[mdsteps, fqgood,error_fq,error_fq,sigmagood]])
        np.savetxt(filenameout+".lifetimesgood.dat",outlifetimes)
        np.savetxt(filenameout+".freqsgood.dat",outfreq)
    print mdstepstocalc
    print lifetimestocheck
    np.savetxt(filenameout+".lifetimes_md_convergence.dat",np.array(np.transpose([mdstepstocalc, lifetimestocheck])))
    return mdstepstocalc, lifetimestocheck

def lammps_log_to_sum(filename, columns = 8, qpoint = False, has_temperature_column = True):
    ''' - this is for the [N,0,0] qpoint
        - the result (sum_all_new) can however be evaluated without knowledge of N
          if lammps could finc c_1 c_2 this could be automated (at least the skipfirst)
          can python also find the end (of a possibly huge file?)  make it fixed first
    '''
    print "starting lammps_log_to_sum ..."
    if type(qpoint) == bool:
        sys.exit("ERROR qpoint is bool")
    skipfirst = None
    skiplast = None

    skipfirst = None
    skiplast = 27

    from itertools import islice
    with open(filename) as myfile:
        head = list(islice(myfile, 200))
        for idx,i in enumerate(head):
            if "c_1 c_2" in i:
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



    #reader = read_csv(filename, sep=' ', header=None,chunksize=1,skiprows=181,skipfooter=30 ,engine='python')
    #data = read_csv(filename, sep=' ', header=None,skiprows=skipfirst,skipfooter=30,engine='python')
    print "read data from lammpsfile with frequencies:",filename
    #print "with 108082514 lines this gives a memory Erroro --> chunk!"
    #data = read_csv(filename,skiprows=skipfirst,engine='python',skipfooter=skiplast,header=None,sep=r"\s*").as_matrix()

    ###################################
    # new approach read in chunks
    ###################################
    print "read data from lammpsfiel with frequencies done!"
    #sum_all_new = np.zeros((columns/2,2*10**8),dtype=np.complex_)  # 4*10^8 seems to be max for cluster
    sum_all_new = np.zeros((columns,2*10**8),dtype=np.complex_)  # 4*10^8 seems to be max for cluster
    chunksize=100000
    chunksize=3;
    #reader = read_csv(filename, skiprows=skipfirst,skipfooter=skiplast,sep=r"\s*",header=None,chunksize=chunksize,engine='python')
    #reader = read_csv(filename,skiprows=skipfirst,engine='python',skipfooter=skiplast,header=None,sep=r"\s*",chunksize=chunksize).as_matrix()

    # works
    #reader = read_csv(filename,skiprows=skipfirst,engine='python',skipfooter=skiplast,header=None,sep=r"\s*").as_matrix()

    # works apparently for log_sum_to_all without header anf footer
    #reader = read_csv(filename, sep=r"\s*",header=None,chunksize=chunksize,engine='python')


    #reader = read_csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize,engine='python',error_bad_lines=False)
    #reader = read_csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize,error_bad_lines=False)
    #reader = read_csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize)
    #reader = read_csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,skipfooter=30)

    # reader for log_for_sum_all
    reader = read_csv(filename, sep=r"\s*",header=None,chunksize=chunksize,engine='python')
    #print "00000000000<<<<<<<<<<<<<<<<--------------"
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
    ##sum_all_new = np.zeros((6,data.shape[0]),dtype=np.complex_)
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
    np.save("sum_all_new_"+qpointstr,sum_all_new)

    #np.save("sum_all_new",sum_all_new)
    return sum_all_new

print "how to run this script: (examples)"
print "     lammps_pos_to_sum.py fcc 10 4.04 10  0  0 40 full"
print "     lammps_pos_to_sum.py fcc 10 4.04 10 10 10 40 chunk"
print "     lammps_pos_to_sum.py bcc  6 3.07 6   6  0 50 lifetimesgood"
print "     lammps_pos_to_sum.py bcc  6 3.07 6   6  6 50 lifetimesgood mev"
print
print
print "convert trj_lammps.out -> trj_lammpsnew.out"
print '     grep "^1 " trj_lammps.out | awk \'{print $2,$3,$4}\' > trj_lammpsnew.out'
print
print
print "sys.argv     :",sys.argv,len(sys.argv)
if len(sys.argv) <= 8:
    sys.exit("ERROR: You need more parameters as input")
structure=sys.argv[1]
N=int(sys.argv[2])
alat=float(sys.argv[3])
qpoints_all = get_all_qpoints(sys.argv[4],sys.argv[5],sys.argv[6],N)
print "structure    :",structure
print "N            :",N
print "alat         :",alat
print "qpoints_all  :",qpoints_all[0],equivalent_qs(qpoints_all[0]).shape[0],equivalent_qs(qpoints_all[0]).shape[0]
sum_all=0
for idx,i in enumerate(qpoints_all):
    sum_all+=equivalent_qs(i).shape[0]
    if idx == 0: continue
    print "              ",i,equivalent_qs(i).shape[0],sum_all

    # This here was just to check that the corresponding entries in sum_all_new.npy
    # are equal: ikk.real == jkk.real and ikk.imag=-1*jkk.imag
    # it is therefore perfectly possible to not put those (symmetry equivalent)
    # qpoints in sum_all_new.new
    # ( This was checked in /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_equivalent_qpoints_if_all_necessary/check)
    #
    #for ikk,kk in enumerate(equivalent_qs(i)):
    #    for jkk,gg in enumerate(equivalent_qs(i)):
    #        #print kk,gg,-1*gg
    #        ch = np.unique(kk==-1*gg)
    #        if len(ch) == 1:
    #            if ch[0] == True:
    #                print ikk,jkk
    #        #if kk == -gg:
    #        #    print ikk,jkk


if structure == "fcc":
    usestruct = 4;
elif structure == "bcc":
    usestruct = 2;
else:
    sys.exit("ERROR: structure:"+str(structure)+" not known!")
dt=int(sys.argv[7])
print "dt [fs]      :",dt
print "sys.argv[8]  :",sys.argv[8]

check_can_be = [ 'lifetimesgood', 'chunk', 'full' ]
if sys.argv[8]in check_can_be:
    check = sys.argv[8]
else:
    sys.exit("you must have: lifetimesgood or chunk or full but do not.")

if check == 'lifetimesgood':
    print_and_save_lifetimesgood(sys.argv)
    sys.exit()


scale_by_alat_N = [ 'trj_lammpsnew.out', 'trj_lammpsnew.out_noJumps_small', 'pos' ] # lammps to direct coords
scale_by_N = [ 'dum', 'posmichael' ] # michaels sim_fcc_morse (to convert to direct)

filename = False
scale_data = False
for i in scale_by_alat_N:
    if os.path.isfile(i) == True:
        filename = i; scale_data = float(alat)*N
for i in scale_by_N:
    if os.path.isfile(i) == True:
        filename = i; scale_data = N
print 'filename     :',filename
print "scale_data   :",scale_data
print "os.getcwd()  :",os.getcwd()
print "##################################################################################"
print

#df = pos_compress_to_sum_fast_read_at_once(filename = filename, qpoints_all = [1,0,0], scale_data = scale_data,N=N,usestruct=usestruct)
#df = pos_compress_to_sum(filename = filename, qpoints_all = [[1,0,0]], scale_data = scale_data, lines_of_one_mdstep = N**3*usestruct, mdsteps_per_chunk = 100, chunkitornot = False)
#
#sys.exit()
#df = pos_compress_to_sum_fast_read_at_once(filename = filename, qpoints_all = [1,0,0], scale_data = scale_data,N=N,usestruct=usestruct)
#
#print
#print
#df = read_hd5f(filename='dump_h5md.h5',qpoint=[1,0,0],scale_data = scale_data)
#sys.exit()



#if __name__ == '__main__':
    # sum_all_new.npy is already twice as large as it would need to be
    # tail -n +16 ../lammps_submit_cmfe.sh.o488667  > log_for_TUTILD
for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
    qpointstr = qpointstring(qpoint)
    print "qpoint:",qpoint,"qpointstring:",qpointstr
    equivalent_qs(qpoint)
    ##########################################
    # check if "sum_all_new.npy" exists
    ##########################################
    if os.path.isfile('sum_all_new_'+qpointstr+'.npy'):
        print "load existing sum_all_new_"+qpointstr+".npy"
        sum_all = np.load('sum_all_new_'+qpointstr+'.npy')
    elif os.path.isfile('log_for_sum_all') == True:
        ##########################################
        # if not, check if "log_for_sum_all" exists
        ##########################################
        sum_all = lammps_log_to_sum(filename = 'log_for_sum_all', qpoint=qpoint)
    elif type(filename) != bool:
        #sum_all = pos_compress_to_sum(\
        pos_compress_to_sum(\
            filename,\
            qpoints_all=qpoints_all,\
            lines_of_one_mdstep=N**3*usestruct,\
            mdsteps_per_chunk = 100, \
            scale_data=scale_data, chunkitornot = check)
        #sum_all = sum_all[idx]  # all qpoint have been evaluated (and written) in pos_compress_to_sum, now load the particular qpoint for sum_all_to_powerspectrum
        sum_all = np.load('sum_all_new_'+qpointstr+'.npy')


    ##########################################
    # after sum_all has been written
    ##########################################
    print "sum_all_new_x_x_x.npy found; now creating powerspectrum ---------------->",os.getcwd()
    print "idx:",idx,"out of:",len(qpoints_all),"qpoint:",qpoint
    print bcolors.FAIL + "######################################################################################" + bcolors.ENDC
    print bcolors.FAIL + "idx "+str(idx)+" out of "+str(len(qpoints_all))+" qpoint: "+str(qpointstr)+ bcolors.ENDC
    print bcolors.FAIL + "######################################################################################" + bcolors.ENDC
    sum_all_to_powerspectrum(sum_all, qpoint=qpoint,dt=dt,appendlast=True)
