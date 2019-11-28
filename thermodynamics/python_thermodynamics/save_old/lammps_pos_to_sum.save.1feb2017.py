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


import time

start = time.time()
from scipy.fftpack import fft,ifft
import sys
import os
import math
import numpy as np
#import h5py
from itertools import permutations
from pandas import read_csv,set_option
from scipy.fftpack import fft
#import pyfftw
#import timeit
import glob
import argparse
import textwrap
from argparse import ArgumentDefaultsHelpFormatter

end = time.time()
print "all python module imported in",str((end-start)),"sec."


np.set_printoptions(suppress=True)   # display arrays withou 000000
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


##########################################################################################
# general functions ######################################################################
##########################################################################################
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
    appropratestring = [ 'l', 't', 't1', 't2', 'all', 'lt', 'tnew' ]
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
    l100 = False
    t100 = False
    t100new = False
    l110 = False
    t1_110 = False
    t2_110 = False
    l111 = False
    t111 = False
    if N1 == 'all' or N2 == 'all' or N3 == 'all': all = True            # [ All ]


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

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 0 0 ] >>>>>>>>>>>>>>>>>>>>>>
    if l100 or all == True:                                         # [ N 0 0 ] (L)
        for i in np.arange(N)+1:
            qpoints_all.append([i,0,0])
    if t100 or all == True:                                         # [ N 0 0 ] (T)
        for i in np.arange(N)+1:
            qpoints_all.append([2*N,i,0])
    #if t100new or all == True:                                         # [ N 0 0 ] (T)
    #    for i in np.arange(N)+1:
    #        qpoints_all.append([2*N,i,2*N])
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 0 0 ] >>>>>>>>>>>>>>>>>>>>>>


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 0 ] >>>>>>>>>>>>>>>>>>>>>>
    if l110 == True or all == True:      # [ N N 0 ] (L)
        for i in np.arange(N)+1:
            qpoints_all.append([i,i,0])
    if t2_110 == True or all == True:    # [ N N 0 ] (T2)
        for i in np.arange(N)+1:
            qpoints_all.append([i,i,2*N])
    if t1_110 == True or all == True:    # [ N N 0 ] (T1)
        for i in np.arange(N)+1:
            qpoints_all.append([2*N+i,2*N-i,0])
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 0 ] >>>>>>>>>>>>>>>>>>>>>>


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 1 ] >>>>>>>>>>>>>>>>>>>>>>
    if l111 == True or all == True:     # [ N N 0 ] (L)
        for i in np.arange(N/2)+1:
            qpoints_all.append([i,i,i])
    if t111 == True or all == True:     # [ N N 0 ] (T)
        for i in (np.arange(N/2)+N/2)[::-1]:
            print "i:",i
            qpoints_all.append([i,i,2*N-i])
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> [ 1 1 1 ] >>>>>>>>>>>>>>>>>>>>>>
        #for i in np.arange(N/2)+1:
        #    qpoints_all.append([N/2+i,N/2-i,0])
    #print "len:",len(qpoints_all)
    if len(qpoints_all) == 0:
        sys.exit("ERROR: this path is not known!")
    return qpoints_all

def print_and_save_lifetimesgood(args = False): #alat, in1, in2, in3, dt):
    ''' keyword is
            - lifetimesgood l 0 0
            - lifetimesgood t2 t2 t2 mev
            - lifetimesgood l l l meV
            - freqsgood t1 t1 0 meV
    '''
    print "-------------------------------------------------------"
    print "---------------- def print_and_save_lifetimesgood() ---"
    print "-------------------------------------------------------"
    #print "xx:",qpoints_all
    #for i in qpoints_all:
    #    print qpointstring(i)
    #alat = sys.argv[3]
    #in1 = sys.argv[4]
    #in2 = sys.argv[5]
    #in3 = sys.argv[6]
    #dt = sys.argv[7]
    #qpoints_all = get_all_qpoints(args.qvec[0],args.qvec[1],args.qvec[2],args.supercell)
    #in1, in2, in3 = args.qvec[0],args.qvec[1],args.qvec[2]
    #print "----in{1,2,3}    :",in1,in2,in3
    #print "----dt           :",args.dt
    #keyword = sys.argv[8]  # {lifetimes,freqs,lifetimesmin}
    #print "----keyword      :",keyword

    dispdir_all = [[args.qvec[0],args.qvec[1],args.qvec[2]]]
    if args.qvec[0]=='all' or args.qvec[1]=='all' or args.qvec[2]=='all':
        dispdir_all = [['l','0','0'],['t','0','0'],['l','l','0'],['t1','t1','0'],['t2','t2','0'],['l','l','l']]
    for dispdir in dispdir_all:
        in1, in2, in3 = dispdir[0],dispdir[1],dispdir[2]
        qpoints_all = get_all_qpoints(in1, in2, in3, args.supercell)
        for keyword in [ 'freqsgood', 'lifetimesgood' ]:
            for keyword2 in [ '', 'meV' ]:
                units = 'THz'
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
                        sys.exit("ERROR: no files foud with name: ps"+qpointstring(i)+"*.dat."+keyword+".dat")



                #print files
                #sys.exit()
                for file in files:
                    #print file
                    f1 = file.split("_")
                    f2 = f1[3].split(".dat")[0]
                    #print "f2:",f2,float(f2),int(float(f2))
                    mdsteps.append(int(float(f2)))  # keep int(float())

                mdstepsmax = np.sort(np.array(list(set(mdsteps))))[-1]
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
                    if in1 == "l" and in2 == "l" and in3 == "l":  # rd=0 and sc=-1 cancells out
                        out = q3;faktor=1.;rd=0.;sc=-1.;
                    if in1 == "t" and in2 == "t" and in3 == "t":  # rd=0 and sc=-1 cancells out
                        out = q3;faktor=1.;rd=args.supercell;sc=-2./args.supercell;

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
                    if in1 == "t" and in2 == "t" and in3 == "t":
                        qmax=1

                    #print "out:",out,"qmax:",qmax
                    #sc*(rd-(float(out)*faktor/float(qmax)))

                #print "keyword:-----------",keyword
                #print "keyword2-----------:",keyword2
                print "-------------------------------------------------------"
                f = open(str(keyword)+"_"+str(in1)+"_"+str(in2)+"_"+str(in3)+"_alat"+str(alat)+"_dt"+str(args.dt)+"_"+str(units)+'.dat', 'wb')
                #outlifetimes = np.array([[mdsteps, ltgood,error_pos,error_neg,sigmagood,999999,error_gaussian_lt_oben_ist_max,error_gaussian_lt_unten_ist_min,lifetimestocheckerrormd[-1]]])
                # idx                         0        1       2         3         4       5          5                                  6                              7

                for file in files:
                    q1 = file.split("ps")[1].split("_")[0]
                    #print file,q1
                    q2 = file.split("ps")[1].split("_")[1]
                    q3 = file.split("ps")[1].split("_")[2]
                    lastline = np.loadtxt(file)
                    #print "file    :",file
                    #print "lastline:",lastline


                    out,faktor,rd,sc = make_phonon_disperion_pic(in1,in2,in3)


                    #print "out:",out,"qmax:",qmax,"faktor:",faktor,"rd:",rd,"sc:",sc,"qmax:",qmax

                    #print lastline

                    roundto=4
                    #print  str(round(sc*(rd-(float(out)*faktor/float(qmax))),4)),"\t",round(lastline[1]*THz_to_meV,roundto),"\t",round(lastline[2]*THz_to_meV,roundto),"\t",round(lastline[3]*THz_to_meV,roundto),"\t",round(lastline[4],7)
                    ka = str(round(sc*(rd-(float(out)*faktor/float(qmax))),4))+"\t"+\
                            str(round(lastline[1]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[2]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[3]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[4],7))+"\t"+\
                            str(999999999)+"\t"+\
                            str(round(lastline[6]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[7]*THz_to_meV,roundto))+"\t"+\
                            str(round(lastline[8]*THz_to_meV,roundto))+"\t"+\
                            "\n"
                    ktoscreen = "\t["+str(q1)+" "+str(q2)+" "+str(q3)+"] qmax:"+str(qmax)+"\n"
                    print ka.rstrip()+ktoscreen.rstrip()
                    f.write(ka)
                f.close()

    sys.exit()

def equivalent_qs(q,show_without_ipi=False,show_lammps_inputfile=False,N=False):
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
    if show_lammps_inputfile == True:
        if type(N) == bool:
            sys.exit('ERROR: you need to specify N if show_lammps_inputfile')
        print out.shape  # z.b. 12,3
        outstr = []
        maxvalue = 0
        qs_computes = []
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
            out2='(('+out+')/${alat}*2.0*PI/'+str(N)+')'
            #print idx,sym_qp
            realcom = [ 'r','i' ]
            for idx2,i in enumerate([ 'cos', 'sin' ]):
                variablename = qpointstring(q)+'_'+str(idx+1)+str(realcom[idx2])
                print 'variable '+variablename+' atom '+i+out2
                qs_computes.append(variablename)
        for variablename in qs_computes:
            print 'compute c'+variablename+' all reduce sum v_'+variablename

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


##########################################################################################
# functions to create space_fft ##########################################################
##########################################################################################
def get_space_fft_from_positions(filename, qpoints_all = False, lines_of_one_mdstep = False, mdsteps_per_chunk = 1, scale_data = 1,chunkitornot = False):
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
        print "qpoint   :",qpoint
        print "qpointstr:",qpointstr,"mdsteps_per_chunk:",str(mdsteps_per_chunk)
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

    ############################################################
    ############################################################
    # start iterating over infile
    ############################################################
    ############################################################
    print "starting iterating over infile ..."
    nextcheck = mdsteps_per_chunk
    for chunk,pos in enumerate(reader):
        #print "chunk:",chunk,"pos:",pos
        #sys.exit()
        if chunk < maxavail:
            if chunk == mdsteps_per_chunk or chunk == nextcheck:
                nextcheck = chunk + mdsteps_per_chunk
            continue
        #print pos.as_matrix()
        #print "scale_data:",scale_data,type(scale_data)
        xxx = pos.as_matrix()/float(scale_data)  # keep the float even if you have already a float
        xxx = np.reshape(xxx,((-1,N**3*usestruct,3)))  # for last step it could be not 3 but 2 oder 1 (32,3)
        mdstepbegin = chunk*mdsteps_per_chunk
        mdstepend = chunk*mdsteps_per_chunk+xxx.shape[0]
        print "chunk (=index):",chunk,"step:",mdstepbegin,"-",mdstepend,"nextcheck:",nextcheck

        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
            qs_all[idx] = equivalent_qs(qpoint)
            sum_all[idx] = np.zeros((len(qs_all[idx]),mdstepend-mdstepbegin),dtype=np.complex_)
            for ind_qs,qs in enumerate(qs_all[idx]):  # fuer alle symmetrieequivalenten qpoints
                # calculate via exp when dumping every chunk (amount number of steps)
                #sum_all[idx][ind_qs] = np.sum(np.exp(np.sum((xxx*qs_all[idx][ind_qs]),axis=2)),axis=1)

                # calculate via exp when dumping only at the very end
                # sum_all[idx][ind_qs,mdstepbegin : mdstepend] = np.sum(np.exp(np.sum((xxx*qs_all[idx][ind_qs]),axis=2)),axis=1
                # OUTPUT TO SCREEN WHEN NECESSARY
                #print "sum_all:",chunk,idx,ind_qs,sum_all[idx][ind_qs].shape,"<-- anzahl mdsteps"
                #pass

                # calculate via eulers formula (sin, cos) when dumping every chunk (amount number of steps)
                a = np.sum((xxx*qs_all[idx][ind_qs]),axis=2)
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

def get_space_fft_from_lammps_log(filename):
    ''' this reads directly from lammps.log; comments above and below the output should be no problem '''

    from_to, from_to_name = get_space_fft_prepare_lammps_log(filename)
    print
    for idx,i in enumerate(from_to):
        print idx,i,from_to_name[idx]
    print
    #reader = read_csv(filename, sep=r"\s*",header=None,engine='c')
    #reader = read_csv(filename, sep=r"\s*",header=None,engine='c')
    print "pandas read_csv ... c engine (15 GB log.lammps takes 29.15% of the cmmc's memory)"
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
            sum_all_new = np.zeros((columns/2,lines),dtype=np.complex_)  # 4*10^8 seems to be max for cluster
            for j in np.arange(columns/2):
                sum_all_new[j] = filecontent[2*j] +1j * filecontent[2*j+1]  #out[:,0] - 1j * out[:,1]
                np.save(filename,sum_all_new)
        else:
            print "columns = 1"
            np.save(filename,filecontent)
            #sum_all_new[5,chunk*chunksize:chunk*chunksize+out.shape[0]] = out[:,0] - 1j * out[:,1]
        #np.save(filename,filecontent.astype(np.complex_))


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
        ##sum_all_new = np.zeros((columns,2*10**8),dtype=np.complex_)  # 4*10^8 seems to be max for cluster
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
    #sum_all_new = np.zeros((columns/2,2*10**8),dtype=np.complex_)  # 4*10^8 seems to be max for cluster
    sum_all_new = np.zeros((columns,2*10**8),dtype=np.complex_)  # 4*10^8 seems to be max for cluster
    chunksize=100000
    chunksize=3;


    #reader = read_csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize,error_bad_lines=False)
    #reader = read_csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,chunksize=chunksize)
    #reader = read_csv(filename, sep=r"\s*",header=None,skiprows=skipfirst,skipfooter=30)

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
    #np.save("sum_all_new_"+qpointstr,sum_all_new)
    np.save("space_fft_"+qpointstr,sum_all_new)

    #np.save("sum_all_new",sum_all_new)
    return sum_all_new

def get_space_fft_from_xaa_xab_folder(qpstring):
    ''' qpstring is a list of all qpoints which are to be evaluated
    e.b. [ '1_0_0', '2_0_0', ...] '''
    foldername="xa[a-z]_"
    filenames_vor_qp="sum_all_new_"
    filenames_nach_qp=".dat"

    folder=sorted(glob.glob(foldername))

    files=sorted(glob.glob(folder[0]+"/"+filenames_vor_qp+"*"+filenames_nach_qp))
    #qpstring= []
    ##print "xx:",qpstring
    ##for i in folder:
    ##    print i
    ##for i in files:
    ##    print i

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
    for idx,i in enumerate(qpstring):
        print "in--->>>:",idx,i
    #sys.exit()

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
    for idxj,j in enumerate(qpstring):
        if os.path.isfile('space_fft_'+j+'.npy'):
            print 'space_fft_'+j,"does already exist...continue"
            continue
        else:
            print 'space_fft_'+j+'.npy',"does not exis"


        #print "idxj:",idxj,"j:",j
        #if j != "1_0_0":
        #    continue
        for idxi,i in enumerate(folder):
            #strout="idxj:",idxj,"j:",j,"         idxi:",idxi,"i:",i
            file=i+"/"+filenames_vor_qp+j+filenames_nach_qp
            #print "importing "+file
            data=np.loadtxt(file)
            #print data
            #print "importing "+file+" done!"

            shape = data.shape
            lines=shape[0]
            columns=shape[1]
            theList=range(columns)
            N=2
            groupby = [theList[n:n+N] for n in range(0, len(theList), N)]
            sum_all_curr = np.zeros((columns/2,lines),dtype=np.complex_)
            if idxi==0: #  and idxj==0:
                print "IDXI +++ 0000000000000000000000"
                sum_all_new = np.zeros((columns/2,lines),dtype=np.complex_)
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

#############
# UNUSED
#############
def get_space_fft_from_positions_fast_read_at_once_UNUSED(filename, qpoints_all = False, scale_data = 1., N=False, usestruct=False):
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

    import time
    start1 = time.time()
    start = time.time()
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
def space_fft_to_powerspectrum(space_fft, qpoint, mdstepstocalcin = False, dt = False, appendlast=True,execute_timeinversion=False, print_single_qpoints=False,idxprint=0,idxprintmax=0, check_qpoint_for_crossing = False, args = False):
    ''' calculates the power_spectrum for a certain qpoint (from space_fft)
    - space_fft contins for a single q point all (symmetry) equivalent q points
    - output of this skript is the power_spectrum for the corresponding q-point
    - qpoint is just for saving the power_spectrum
    '''
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
        #print "--> mdstepstocalc:", mdstepstocalc# ,mdstepstocalc[np.where(mdstepstocalc > 999)[0]]
        if type(mdstepstocalcin) != bool:
            if mdstepstocalcin == 'last':
                mdstepstocalc = np.array([mdstepstocalc[-1]])
            elif type(mdstepstocalcin) == int:
                mdstepstocalc = np.array([mdstepstocalcin])
            else:
                sys.exit("ERROR 777 check here")
        #print "--> mdstepstocalc:", mdstepstocalc# ,mdstepstocalc[np.where(mdstepstocalc > 999)[0]]
        return mdstepstocalc


    mdstepstocalc = which_mdsteps_to_calc(space_fft=space_fft,mdstepstocalcin = mdstepstocalcin)

    #print mdstepstocalc
    #addtothis = np.arange(0,mdstepstocalc[0],mdstepstocalc[0]/10)[1:]
    #mdstepstocalc = np.insert(mdstepstocalc,0,addtothis)
    #print
    #print mdstepstocalc

    #mdstepstocalc = mdstepstocalc[np.where(mdstepstocalc > 999)]
    #sys.exit()









    #mdstepstocalc = np.array([ 300000 ])
    lifetimestocheck = np.zeros(len(mdstepstocalc))
    lifetimestocheckerrormd = np.zeros(len(mdstepstocalc))
    freqstocheck = np.zeros(len(mdstepstocalc))


    ##################################################################################
    ##################################################################################
    # get error due to MD convergenz (this takes most of the time)
    ##################################################################################
    ##################################################################################
    #####################################################
    # schleife ueber mdsteps
    # it would be best to start with the longest, write this, and then move to the shorter runs
    #####################################################

    #for idx,mdsteps in enumerate(mdstepstocalc[::-1]):  # if done in this way, then lifetimesgood.dat and freqsgood.dat wont be written
    for idx,mdsteps in enumerate(mdstepstocalc):
        qpointstr = qpointstring(qpoint)
        power_spectrum=np.zeros(mdsteps)
        print bcolors.FAIL + "  space_fft_to_powerspectrum: "+str(idxprint+1)+"/"+str(idxprintmax)+"   Innerloop: "+str(idx+1)+"/"+str(len(mdstepstocalc))+" mdstep:"+str(mdsteps) +" qpoint:"+str(qpoint) +bcolors.ENDC


        ##########################################################################
        # calculate power spectrum
        ##########################################################################
        for i in np.arange(len(space_fft)): # laeuft ueber alle symmetrieequivalten qs (6 stueck)
            a =np.abs(fft((space_fft[i][:mdsteps])[:mdsteps]))**2.  # put into one command to save (hopefully) memory
            if print_single_qpoints == True and mdsteps == mdstepstocalc.max():
                tmpfilename = "ps"+qpointstr+"__"+str(i)+"__"+str(mdsteps)+".dat"
                np.savetxt(tmpfilename,a)

            power_spectrum+=a

            if execute_timeinversion:
                a = np.roll(a[::-1],1)
                if print_single_qpoints == True and mdsteps == mdstepstocalc.max():
                    tmpfilename = "ps"+qpointstr+"__"+str(i)+"-__"+str(mdsteps)+".dat"
                    np.savetxt(tmpfilename,a)
                power_spectrum+=a

        ##########################################################################
        # write out power_spectrum (simple list, only y without x values;
        # length of mdsteps/2 which changes in this loop)
        ##########################################################################
        #power_spectrum[0:10] = 0;  # sometimes very high values , 10 is too much for md with 100 steps
        power_spectrum[0:2] = 0;  # power_spectrum is just a list of values (only y without x)
        ps=power_spectrum/power_spectrum.max()*.9; # full with peaks left and right
        factor_psout = 2  # 2 is a save bet since power spectrum is doubled
        filenameout = "ps"+qpointstr+"_"+str(mdsteps)+".dat"

        ##########################################################################
        # write out power_spectrum
        ##########################################################################
        if args.write_full_ps and mdsteps == mdstepstocalc.max():
            print "qp:",qpoint,type(qpoint)
            print "write out power_spectrum ...", filenameout
            np.savetxt(filenameout,ps[:len(ps)/factor_psout])  # saves only half the powerspektrum (left peak)

        ##########################################################################
        # write out power_spectrum_cut
        ##########################################################################
        if check_qpoint_for_crossing == True:
            ps = get_goodsmoothing_for_ps_with_crossing(ps)   # ist zurzeit nur correct fuer [ 9 9 20] und [ 10 10 20 ]
            print "1 ps.max()",ps.max()
            ps=ps/ps.max()*.9; # full with peaks left and right
            print "2 ps.max()",ps.max()
            # if args.write_full_ps and mdsteps == mdstepstocalc.max():
            if mdsteps == mdstepstocalc.max():  # better write in this cases
                print "qp:",qpoint,type(qpoint)
                print "write out power_spectrum_cut ...", filenameout
                filenameoutcut = "ps"+qpointstr+"_"+str(mdsteps)+"cut.dat"
                np.savetxt(filenameoutcut,ps[:len(ps)/factor_psout])  # saves only half the powerspektrum (left peak)

        ##########################################################################
        # change ps to ps 2d which is smaller and quicker to calculate
        ##########################################################################
        xy,x,ps = powerspectrum_to_powerspectrum_2d_sparse(ps)  # y = ps

        ##########################################################################
        # write out power_spectrum_smoothed (simple list, only y without x values
        # length of mdsteps/2 which changes in this loop)
        # out1 = smoothing_power_spectrum(ps,smoothingfact)  # out1 is as the powerspectrum have both peaks and has exactly as many as many points
        ##########################################################################
        if True and mdsteps > 999:
            #print "mdsteps:",mdsteps
            sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin,minsmoothfunc = get_goodsmoothing_for_ps(ps,dt,mdsteps,allowed_func = 1, allowed_der = False,args=args,stringadd='min')
            print "-------------1",sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin,minsmoothfunc
            sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood,goodsmoothfunc = get_goodsmoothing_for_ps(ps,dt,mdsteps,allowed_func = 0, allowed_der = 0,args=args,stringadd='good')
            print "-------------2",sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood,goodsmoothfunc

            if check_qpoint_for_crossing == True:
                sigmagood, fqgood,ltgood,ltgoodmin,ltgoodmax,xmaxwritegood,goodsmoothfunc = sigmamin, fqmin,ltmin,ltminmin,ltminmax,xmaxwritemin,minsmoothfunc
            #print sigmamin,sigmagood,ltmin,ltgood
            #out1out = out1[:len(out1)/factor_psout]
            #############################################################
            # write here the smoothigfunc, further down the lifetimes
            #############################################################
            if mdsteps == mdstepstocalc.max():
                filenameoutsmooth = filenameout+".smooth_good"+str(sigmagood)+".dat"
                filenameoutmin = filenameout+".smooth_min"+str(sigmamin)+".dat"
                print "  write out power_spectrum_smoothed sigma:",sigmagood,"...",filenameoutsmooth

                #np.savetxt(filenameoutsmooth,goodsmoothfunc[:xmaxwritegood])   # old way for 1d
                #np.savetxt(filenameoutmin,minsmoothfunc[:xmaxwritemin])        # old way for 1d
                # adapt to 2d which is smaller
                np.savetxt(filenameoutsmooth,np.transpose(np.array([x,goodsmoothfunc]))[:xmaxwritegood],fmt='%.1f %.5f')
                np.savetxt(filenameoutmin,   np.transpose(np.array([x,minsmoothfunc]))[:xmaxwritemin],fmt='%.1f %.5f')

        #lifetimestocheck[::-1][idx] = ltgood  # this is only to check error of mdlength  ## if done with [::-1] lifetimesgood.dat and freqsgood.dat wont be written
        #freqstocheck[::-1][idx] = fqgood      # this is only to check error of mdlength  ## if done with [::-1] lifetimesgood.dat and freqsgood.dat wont be written
        lifetimestocheck[idx] = ltgood  # this is only to check error of mdlength
        freqstocheck[idx] = fqgood      # this is only to check error of mdlength
        if args.verbose: print 'idx:',idx,'ltgood:',ltgood,"mdsteps:",mdsteps,"ltcheck:",lifetimestocheck
    ##################################################################################
    ##################################################################################
    # get error due to MD convergenz DONE
    ##################################################################################
    ##################################################################################
    if len(mdstepstocalc) <= 1:
        print "now not doing error due to MD convergence"
        return



    if args.verbose:
        print
        print
        print "get error due to MD convergenz DONE"
        print
        print

    #############################################
    # get error due to MDlength
    # mdstepstocald have to be in increasing order
    #############################################
    if args.verbose:
        print "lifetimestocheck (== lifetime als funktion des mdsteps),sigmagood",lifetimestocheck,sigmagood
        print lifetimestocheck.min()
        print lifetimestocheck.max()
    error_md_lt = abs(lifetimestocheck[:-6].max() - lifetimestocheck[:-6].min())

    error_md_fq = abs(freqstocheck.max() - freqstocheck.min())
    if args.verbose:
        print "ltgood:",ltgood
        print "ltmin:",ltmin

    error_gaussian_lt_unten_ist_min =  abs(ltgood - np.array([ltmin,ltgoodmin,ltminmin]).min())   # (zweiter teil von) fehler nach unten
    error_gaussian_lt_oben_ist_max =  abs(ltgood - np.array([ltgoodmax,ltminmax]).max())   # (zweiter teil von) fehler nach unten
    error_gaussian_fq =  abs(fqgood - fqmin)
    error_fq = error_md_fq + error_gaussian_fq
    #print error_gaussian_lt
    #print aerror_gaussian_fq
    #print "mdsteps:",mdsteps

    #############################################
    # get error due to MDlength
    #       - ps*lifetimes_md_convergence.dat
    #       - ps*lifetimes_md_convergencepsec.dat
    #############################################

    if args.verbose:
        print "mdstepstocalc   :",mdstepstocalc
        print "lifetimestocheck:",lifetimestocheck
    idxmax = np.where(lifetimestocheck == lifetimestocheck.max())[0].max()
    for ind,yy in enumerate(lifetimestocheck):
        #print ind,y,"   ",lifetimestocheck[:ind+1][-6:]
        considerlt = lifetimestocheck[:ind+1][-5:]
        lifetimestocheckerrormd[ind] = considerlt.max() - considerlt.min()

    idxmax = np.where(lifetimestocheckerrormd == lifetimestocheckerrormd.max())[0].max()
    for ind,yy in enumerate(lifetimestocheckerrormd):
        if ind <= idxmax:
            lifetimestocheckerrormd[ind] = lifetimestocheckerrormd[idxmax]
    np.savetxt(filenameout+".lifetimes_md_convergence.dat",np.array(np.transpose([mdstepstocalc, lifetimestocheck, lifetimestocheckerrormd,lifetimestocheckerrormd])),fmt='%.5f')
    ftt = factor_mdsteps_to_time = 10**-15*dt/(10**-12)  # in ps
    np.savetxt(filenameout+".lifetimes_md_convergencepsec.dat",np.array(np.transpose([mdstepstocalc*ftt, lifetimestocheck, lifetimestocheckerrormd,lifetimestocheckerrormd])),fmt='%.5f')

    #############################################
    # get final error for lifetimesgood
    # this will only work if loop is done from lowest to highest md step
    #############################################
    if args.verbose: print "lifetimestocheckerrormd:",lifetimestocheckerrormd
    if mdsteps == mdstepstocalc.max():
        filenameoutsmooth = filenameout+".smooth_good"+str(sigmagood)+".dat"

        error_pos = error_gaussian_lt_oben_ist_max+lifetimestocheckerrormd[-1]
        error_neg = error_gaussian_lt_unten_ist_min+lifetimestocheckerrormd[-1]
        outlifetimes = np.array([[mdsteps, ltgood,error_pos,error_neg,sigmagood,999999,error_gaussian_lt_oben_ist_max,error_gaussian_lt_unten_ist_min,lifetimestocheckerrormd[-1]]])
        outfreq      = np.array([[mdsteps, fqgood,error_fq ,error_fq ,sigmagood,999999,error_gaussian_fq             ,error_gaussian_fq              ,error_md_fq]])
        #print 'outlifetimes:',outlifetimes
        np.savetxt(filenameout+".lifetimesgood.dat",outlifetimes,fmt='%.5f')
        np.savetxt(filenameout+".freqsgood.dat",outfreq,fmt='%.5f')



    return mdstepstocalc, lifetimestocheck

def smoothing_power_spectrum(power_spectrum,sigma):
    ''' smoothes a gaussian by sigma '''
    steps = power_spectrum.shape[0]                 # 1000
    nu = np.arange(0,steps)/float(steps)            # [ 0., 0.001, ..., 0.999 ]
    kern = np.exp(-(np.mod(nu+0.5,1)-0.5)**2/(2*sigma**2))  # Gaussian distribution in 1-D
    kern = kern/np.sum(kern)
    ps = power_spectrum_smooth = np.real(ifft(fft(power_spectrum)*fft(kern)))
    ps=ps/np.max(ps)*.9;
    return ps

def get_uncertainty_smoothed_powerspectrum(power_spectrum_smooth,dt,sigma = "",args=False,stringadd=''):
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

    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]


    # x index of absolute maximum
    x_ind_max = np.where(datahalf == datahalf.max())[0][0]  # 963  (9.63 THz)
    #print "x_ind_max:",x_ind_max

    if x_ind_max == 0:
        return 0 ,0, 0 ,xvalues, 0,0,0,99999999999999999999999999

    ############################################################
    # analyse the hight of the background
    #                          .
    #                         .|.
    #                        . | .
    #                       .  |  .
    #                     .    |   .
    #                   .      |    .
    #                 .        |     .
    #               .          |      .
    #.           .             |       .
    #  .      .                |        .
    #     .                    |           .
    #     |                    |               .  .   .
    ############################################################
    #     x                    |                        y_val_leftdata_min
    #     x                    |                        x_ind_leftdata_min
    #                          x                        x_ind_max
    #x                                                  y_val_leftdata_max
    ############################################################
    y_val_leftdata_min = datahalf[:x_ind_max].min()
    x_ind_leftdata_min = np.where(datahalf[:x_ind_max] == y_val_leftdata_min)[0][0]
    #print "x_ind_leftdata_min:",x_ind_leftdata_min

    if x_ind_leftdata_min == 0:
        y_val_leftdata_max = datahalf[0]
    else:
        y_val_leftdata_max = datahalf[:x_ind_leftdata_min].max()
    #print "y_val_leftdata_min:",y_val_leftdata_min


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
    # dis ist nur schoen bei einem relativ glatten gaussian
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


    l = len(np.where(np.diff(datahalf[x_ind_left_over3:x_ind_max])<0)[0])
    #print l
    #print "bbb"
    r = len(np.where(np.diff(datahalf[x_ind_max:x_ind_right_over3])>0)[0])
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
    #print "x_ind_left_over3:",x_ind_left_over3
    if (x_ind_max - x_ind_left_over3) > 1:
        datatake = datahalf[x_ind_left_over3:x_ind_max]
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
    if (x_ind_right_over3-x_ind_max) > 1:
        datatake = datahalf[x_ind_max:x_ind_right_over3]
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
    absmin = np.where(datahalf>y_lifetime)[0][0]     # nur die spitze des gaussians (linke seite)
    absmax = np.where(datahalf>y_lifetime)[0][-1]    # nur die spitze des gaussians (rechte seite)


    #mult = dt*10 / float(xvalues)  ## Wrong!  see example 6x6x6sc 6_0_0 qpoint between dt={10,20}
    mult = 1000./(dt*float(xvalues))
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
    if args.verbose:
        #ljust adds zeros (0) to the end
        print "Freq "+stringadd+":",str(round(freq,2)).ljust(4,'0'),"(THz); width: ",str(round(lifetime,2)).ljust(4,'0'),"(THz)", "||", x_ind_left_over2,x_ind_max,x_ind_right_over2,"||", "sigma:",str(sigma),"r+l",r+l,"rr,ll",ll+lr+rl+rr
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

    print "STARTING loop 0 -> (logarithmic loop) 1e-07 .. 1e-01" #,loop_smoothing
    working = []
    minimum = []
    x_max1list = []
    x_max2list = []
    minimumalt = []
    for smoothingfact in loop_smoothing[::-1]:
        if len(working) >= 2: continue   # to orders of magnitude are enough, the third one often fails
        out = smoothing_power_spectrum(ps,smoothingfact)  # ps and out1 are the full power spectrum (left and right peak)
        if out[0] > 0.5 : continue
        working.append(True)

        # debug line crossing
        #np.diff(out)
        #outsave=np.diff(out)/np.max(np.diff(out))*.9;
        #np.savetxt("pssmoothed"+str(smoothingfact)+"diff.dat",out)

        #f,lt,ltmin,ltmax,xv,rl,rrll,xcutoff = get_uncertainty_smoothed_powerspectrum(out,dt,smoothingfact)
        outhalfright = out[:len(out)/2]
        outhalfleft = np.copy(outhalfright)
        #np.savetxt("newpssmoothed"+str(smoothingfact)+".dat",out1)
        x_max1 = np.where(outhalfright == outhalfright.max())[0][0]    # das ist aber fuer den [ 8 8 20 ] schon der zweite maximalwert
        print "X_MAXXXXX1:",x_max1


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

    if args.verbose: print "STARTING loop 0 -> (logarithmic loop) 1e-07 .. 1e-01" #,loop_smoothing
    bestworking=False
    #np.savetxt("ps.dat",ps)
    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps and out1 are the full power spectrum (left and right peak)

        # debug line crossing
        #np.savetxt("pssmoothed"+str(smoothingfact)+".dat",out1)
        #np.diff(out1)
        #outsave=np.diff(out1)/np.max(np.diff(out1))*.9;
        #np.savetxt("pssmoothed"+str(smoothingfact)+"diff.dat",outsave)

        f,lt,ltmin,ltmax,xv,rl,rrll,xcutoff = get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)

        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    loop_smoothing = smoothingfact/10.*np.array([1,2,3,4,5,6,7,8,9,10])  # redo 10 for easier coding, 10 will definitively work.
    loop_smoothing_next = smoothingfact/100.*np.array([1,2,3,4,5,6,7,8,9,10])  # redo 10 for easier coding, 10 will definitively work.
    if type(bestworking) == bool:
        print "no bestworking found! maybe band crossing??? check this!!!!!! in this case you would have 2 gaussians and you'd need to get the linewidths of both"
        bestworking = 1.0
    if args.verbose: print "DONE loop 1 ->","best:",bestworking,"try:",loop_smoothing

    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps is the full power spectrum (left and right peak)
        f,lt,ltmin,ltmax,xv,rl,rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    #print "loop_smoothing:",loop_smoothing
    delta = loop_smoothing[0]/10.
    #print "delta:",delta
    if args.verbose: print "best:",bestworking
    #print np.arange(1,11)[::-1]
    loop_smoothing = bestworking - delta * np.arange(0,11)[::-1]  # include 0 so bestworking will definitively work (easier code)
    loop_smoothing = np.trim_zeros(loop_smoothing)
    #order_magnitude = -1.*(int(math.log10(min_smoothing))-1)
    loop_smoothing = loop_smoothing[loop_smoothing > 0]
    if args.verbose: print "DONE loop 2 ->","best:",bestworking,"try:",loop_smoothing
    for smoothingfact in loop_smoothing:
        out1 = smoothing_power_spectrum(ps,smoothingfact)  # ps is the full power spectrum (left and right peak)
        f,lt,ltmin,ltmax,xv,rl, rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,smoothingfact,args=args,stringadd=stringadd)
        if rl <= allowed_func and rrll <= allowed_der:
            bestworking = smoothingfact
            break
    if args.verbose:
        print "DONE loop 3 ->",smoothingfact,"best:",bestworking

        print "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
        print "bestworking:",bestworking
    sigma = bestworking
    out1 = smoothing_power_spectrum(ps,sigma)  # ps is the full power spectrum (left and right peak)
    fqgood, ltgood, ltmin, ltmax, xv, rl, rrll, xcutoff= get_uncertainty_smoothed_powerspectrum(out1,dt,sigma,args=args,stringadd=stringadd)

    return sigma, fqgood,ltgood,ltmin,ltmax,xcutoff,out1   # out1 is the smoothed function

def powerspectrum_to_powerspectrum_2d_sparse(ps):
    ''' print "len of ps:",len(ps)
        print " immer so fitten dass mindestens 5000 werte drin bleiben "
    '''
    def get_teiler_daten(ps):
        lenps = len(ps)
        teiler = 1
        for i in np.arange(16)+1:
            tryteiler = 1*10**i
            #print i,teiler,lenps/teiler,(lenps/teiler)>=5000
            if ((lenps/tryteiler)>=5000) == True:
                #print i,tryteiler,"Ja"
                teiler = tryteiler
            else:
                #print i,tryteiler,"NEIN"
                break
        #print "teiler:",teiler
        N = teiler
        return N
    N = get_teiler_daten(ps)
    x = np.arange(len(ps))[::N]+N/2
    y = np.mean(ps.reshape(-1, N), 1)
    y = y/y.max()*0.9
    #print "x.sh:",x.shape
    #print "y.sh:",y.shape
    return np.transpose(np.array([x,y])),x,y

def check_qpoint_for_crossing(qpoint, N, structure):
    '''
    just gives true for fcc t1 for qpoint [[ 10 10 20 ], [ 9 9 20 ], [ 8 8 20 ]]
    '''
    #print "iiiiiiiiiii",qpoint,N

    if qpoint[0]==qpoint[1] and qpoint[2]==2*N:  # bedingutn t1
        if structure == 'fcc': # bedingung fcc
            if float(qpoint[0])/float(N)>=float(8)/float(10):   # nur [ {8,9,10} {8,9,10} 20]
                return True
    else:
        return False

def help(p = None):
    string = '''
    examples how this scritp was started previously:
    - lammps_pos_to_sum.py fcc 10 4.04 10  0  0 40 full
    - lammps_pos_to_sum.py fcc 10 4.04 10 10 10 40 chunk
    - lammps_pos_to_sum.py bcc  6 3.07 6   6  0 50 lifetimesgood
    - lammps_pos_to_sum.py bcc  6 3.07 6   6  6 50 lifetimesgood

    new way:
    lammps_pos_to_sum.py -N 10 -a 4.13 -dt 40 -q l 0 0
    lammps_pos_to_sum.py -N 3 -a 4.14 -dt 1 -q lt 0 0 -lt
    lammps_pos_to_sum.py -N 3 -a 4.14 -dt 1 -q l l l -ps -fftpy


    how to convert trj_lammps.out -> trj_lammpsnew.out:
             grep "^1 " trj_lammps.out | awk \'{print $2,$3,$4}\' > trj_lammpsnew.out

    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)

    p.add_argument('-s',   '--structure',choices=[ 'fcc', 'bcc' ],
       help='which structure was calculated? Currently only fcc and bcc are supported. Simple cubic should work too.', type=str, default='fcc')
    p.add_argument('-n', '-N','-sc','--supercell' , required=True,
       help='supercell: 2x2x2 -> 2', type=int, default=False)
    p.add_argument('-dt',   '--dt' , required=True,
       help='timestep of your calculation in femto seconds.', type=float, default=False)
    p.add_argument('-a',   '--alat' , required=True,
       help='alat of your cell in angstromr; e.g. 4.13', type=float, default=False)
    p.add_argument('-q',   '--qvec' , required=True, nargs=3, type=str, default=False,
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

    p.add_argument('-space_fft','--space_fft', action='store_true', default=False,
        help='get space_fft_x_x_x.npy from xaa_, xab_, ... folder')
    p.add_argument('-lt','--lifetimesgood', action='store_true', default=False,
        help='print and save summary of the lifetimes and freqs in the specified q-vecotor direction; (both in meV and THz);')
    p.add_argument('-freqs','--freqsgood', action='store_true', default=False,
        help='print and save summary of the lifetimes and freqs in the specified q-vecotor direction; (both in meV and THz);')

    p.add_argument('-showeq','--showeq', action='store_true', default=False,
        help='print equivalent qpoint to screen;')
    p.add_argument('-showlammps','--showlammpsinput', action='store_true', default=False,
        help='print input for lammps inputfile for this paritcular qpoint(s) to screen;')

    p.add_argument('-fftpy','--make_space_fft_from_py', action='store_true', default=False,
        help='create space_fft_x_x_x.npy file using pandas read_csv;')
    p.add_argument('-fftlammpslog','--make_space_fft_from_lammpslog', action='store_true', default=False,
        help='create space_fft_x_x_x.npy lammps.log file using numpy genfromtxt;')
    p.add_argument('-fftc','--make_space_fft_from_c', action='store_true', default=False,
        help='create space_fft_x_x_x.npy file using c++ skript;')
    p.add_argument('-ps','--make_power_spectrum', action='store_true', default=False,
        help='create ps_x_x_x_MDSTEPS.dat file and evaluate error due to MD/smoothing/background;')
    p.add_argument('-calclast','--calclast', action='store_true', default=False,
        help='calculate powerspectrum only for all (maximum number of) steps and not intermediate step values to estimate the MD error;')
    p.add_argument('-write_full_ps','--write_full_ps', action='store_true', default=False,
        help='Write the full powerspectrum to hd. Files can be relatively larg, 15-150 MB')


    p.add_argument('-smoothing',   '--smoothing' , required=False,
       help='wite to hd a particular smoothing ; e.g. -mdsteps 30000 -smoothing 0.001', type=float, default=False)
    p.add_argument('-mdsteps',   '--mdsteps' , required=False,
       help='wite to hd for a particular md length; e.g. -mdsteps 30000 -smoothing 0.001', type=int, default=False)
    p.add_argument('-v','--verbose',
            help='verbose', action='store_true', default=False)
    return p

if __name__ == '__main__':
    p = help()  # this gives the possibility to change some __init__ settings
    args = p.parse_args()
    structure = args.structure
    N = args.supercell
    alat = args.alat
    dt = args.dt

    if structure == "fcc": usestruct = 4;
    if structure == "bcc": usestruct = 2;
    qpoints_all = get_all_qpoints(args.qvec[0],args.qvec[1],args.qvec[2],N)

    print "structure    :",structure
    print "N(supercell) :",N
    print "alat         :",alat
    print "dt [fs]      :",dt
    print "qvec         :",args.qvec
    print "qpoints_all  : [ qpoint ] \tqpstr\teq\tsum                "
    print "----------------------------------------------------"
    qpstring_all = [];tmp=0;
    for idx,i in enumerate(qpoints_all):
        tmp+=equivalent_qs(i).shape[0]
        print "              ",i,'\t',qpointstring(i),'\t',equivalent_qs(i).shape[0],'\t',tmp
        qpstring_all.append(qpointstring(i))
        if args.showeq: print equivalent_qs(i,show_without_ipi=True)
        if args.showlammpsinput: print equivalent_qs(i,show_lammps_inputfile=True,N=N)
    print "----------------------------------------------------"
    print "qpstring_all :",qpstring_all
    print

    if args.space_fft: get_space_fft_from_xaa_xab_folder(qpstring=qpstring_all)
    #if args.lifetimesgood or args.freqsgood: print_and_save_lifetimesgood(alat,args.qvec[0],args.qvec[1],args.qvec[2],dt)
    if args.lifetimesgood or args.freqsgood: print_and_save_lifetimesgood(args = args)

    ########################################
    # look for filenames
    ########################################
    scale_by_alat_N = [ 'trj_lammpsnew.out', 'trj_lammpsnew.out_noJumps_small', 'pos', 'POSITIONs' ] # lammps to direct coords
    scale_by_N = [ 'dum', 'posmichael' ] # michaels sim_fcc_morse (to convert to direct)

    print
    print
    print "FILENAMES:"
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
    print
    print

    if args.make_space_fft_from_py or args.make_space_fft_from_lammpslog or args.make_power_spectrum:
        for idx,qpoint in enumerate(qpoints_all): # fuer alle qpunkte in der liste
            qpointstr = qpointstring(qpoint)
            print bcolors.FAIL + "######################################################################################" + bcolors.ENDC
            print bcolors.FAIL + "OUTER LOOP OVER ALL QPOINTS: "+str(idx)+" out of "+str(len(qpoints_all))+" qpointstring: "+qpointstr + bcolors.ENDC
            print bcolors.FAIL + "######################################################################################" + bcolors.ENDC
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
                print "load existing ",space_fft_file1
                sum_all = np.load(space_fft_file1)
            elif os.path.isfile(space_fft_file2):
                print "load existing ",space_fft_file2
                sum_all = np.load(space_fft_file2)
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
                    lines_of_one_mdstep=N**3*usestruct,\
                    mdsteps_per_chunk = 4000, \
                    scale_data=scale_data, chunkitornot = False)
            ###############################################################################
            # AT LATEST HERE space_fft_xxxx.npy has to exist
            ###############################################################################


            ##########################################
            # LOAD EXISTING space_fft_xxxx.npy
            ##########################################
            file_to_load = 'space_fft_'+qpointstr+'.npy'
            if os.path.isfile(file_to_load) != True: sys.exit("file "+file_to_load+" does not exist! Run this skript with -fftpy or -fftc first!")

            if args.make_power_spectrum:
                sum_all = np.load(file_to_load)

                ##########################################
                # after sum_all has been written
                ##########################################
                print "sum_all_new_x_x_x.npy found; now creating powerspectrum ---------------->",os.getcwd()
                print bcolors.FAIL + "  ##################################################################################################################" + bcolors.ENDC
                print bcolors.FAIL + "  idx "+str(idx)+" out of "+str(len(qpoints_all))+" qpoint: "+str(qpointstr)+" filename:"+str(file_to_load)+ bcolors.ENDC
                print bcolors.FAIL + "  ##################################################################################################################" + bcolors.ENDC

                mdstepstocalcin = False
                if args.calclast:
                    mdstepstocalcin = 'last'
                space_fft_to_powerspectrum(\
                        space_fft=sum_all, \
                        qpoint=qpoint,\
                        mdstepstocalcin=mdstepstocalcin,\
                        dt=dt,\
                        appendlast=True,\
                        execute_timeinversion=execute_timeinversion,\
                        print_single_qpoints=False,\
                        idxprint=idx,idxprintmax=len(qpoints_all),\
                        check_qpoint_for_crossing = check_qpoint_for_crossing(qpoint, N, structure),\
                        args = args\
                        )
