#!/usr/bin/env python

import sys
import os
import numpy as np
from pandas import read_csv
from scipy.fftpack import fft
import timeit

alat=4.13;
N = 2;


base="/Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2015.08_phonon_lifetimes_2_for_munich/calculations/"
filename_folder=base+"allq_Al_"+str(N)+"_4.13_1000000steps_timestep0.001_totaltime10000ps_900K_dump_10/"
filename=filename_folder+"trj_lammpsnew.outo1"
filename=filename_folder+"trj_lammpsnew.out"


filename=sys.argv[1];
N=int(sys.argv[2]);
print "filename:",filename
print "N:",N
#filename="/Users/glensk/tmp_lammps_pos_workon"


def equivalent_qs(q):
    ''' - finds all the equivalent q vectors for a certain q vector
        - q can be a list or a numpy array e.b. [2,0,0] or np.array([2,0,0])
        currently only for [ 1 0 0 ] and not general
    '''
    #print "type(q):",type(q)
    if type(q) == list:
        q = np.array(q)
    if type(q) != np.ndarray:
        sys.exit("q needs to be a list or a numpy array!")
    if len(q) != 3:
        sys.exit("q needs to be of length 3")

    if q[1] != 0:
        sys.exit("currently for [N,0,0] implemented")
    if q[2] != 0:
        sys.exit("currently for [N,0,0] implemented")
    N = q[0]
    equivalents_fixq = np.array([[N,0,0],[0,N,0],[0,0,N],[0,0,-N],[0,-N,0],[-N,0,0]])
    return equivalents_fixq*1j*2*np.pi

def pos_compress_to_sum_equivalent_qs(filename, qpoint = False, lines_of_one_mdstep = False, mdsteps_per_chunk = 1, scale_data = 1):
    ''' reads a file in chunks to save memory and creates the sum (for a certain list of
        qpoints) which can be directly fourier transformed.

    - per chunk, amount of lines given by 'lines_of_one_mdstep' (=atoms_in_supercell) times
      mdsteps_per_chunk

    - example:
        filename = pos_lammps
        qpoint = np.array([1,0,0])
        lines_of_one_mdstep = 32    # == number of atoms; for a 2x2x2 fcc supercell
        mdsteps_per_chunk = 100     # this would read 32*100 lines per chunk
        scale_data = 4.13*2         # 4.13 Angstrom * 2; this is necessary to get DIRECT coords
    '''
    if os.path.isfile(filename) != True:
        print "filename:",filename
        sys.exit("ERROR: filenmae does not exist")
    print "filename:",filename
    print "scale_data:",scale_data,type(scale_data)
    qs_all = equivalent_qs(qpoint)
    print "define sum_all np.zeros: in case this raises error you need more memory or decrease steps!"
    # - allocation of sum_all depends on memory available (np.zeros(...10^12))
    # - in principle this could be sepereated for every q point but would mean that one
    #   would have to read the file multiple times. (or write the sum to disk once more
    #   steps which could also be a good solution)
    sum_all = np.zeros((len(qs_all),2*10**8),dtype=np.complex_)  # 4*10^8 seems to be max for cluster
    print "define sum_all np.zeros DONE"
    #print "lines_of_one_mdstep*mdsteps_per_chunk:",lines_of_one_mdstep*mdsteps_per_chunk
    reader = read_csv(filename, sep=' ', header=None,chunksize=lines_of_one_mdstep*mdsteps_per_chunk,engine='python')
    #reader = read_csv(filename, sep=' ', header=None,chunksize=lines_of_one_mdstep*mdsteps_per_chunk)
    for chunk,pos in enumerate(reader):
        #print "scale_data:",scale_data,type(scale_data)
        xxx = pos.as_matrix()/float(scale_data)  # keep the float even if you have already a float
        xxx = np.reshape(xxx,((-1,N**3*4,3)))  # for last step it could be not 3 but 2 oder 1 (32,3)
        #xxx = np.reshape(xxx,((mdsteps_per_chunk,N**3*4,3)))  # for last step it could be not 3 but 2 oder 1 (32,3)
        #print "xxx.shape:",xxx.shape[0]
        #print "xxx.shpae",xxx.shape
        #sum[mdstep] = np.sum(np.exp(np.sum(xxx*qs,axis=1)))  # one value per md step
        mdstepbegin = chunk*mdsteps_per_chunk
        mdstepend = chunk*mdsteps_per_chunk+xxx.shape[0]
        #mdstepend = chunk*mdsteps_per_chunk+mdsteps_per_chunk
        #sys.stdout.write('\r'+str(mdstepend))
        print "step:",mdstepend
        for ind_qs,qs in enumerate(qs_all):
            #print "ind_qs,qs",ind_qs,qs,mdstepbegin,mdstepend
            #print "--> 1:",qs_all[ind_qs]
            #aa = xxx*qs_all[ind_qs]
            #print "--> 2:"
            #bb = np.sum((xxx*qs_all[ind_qs]),axis=2)
            #print "--> 3:"
            #cc = np.sum(np.exp(np.sum((xxx*qs_all[ind_qs]),axis=2)),axis=1)
            #print "--> 4:",cc.shape
            #print "--> 5:"
            sum_all[ind_qs,mdstepbegin : mdstepend] = np.sum(np.exp(np.sum((xxx*qs_all[ind_qs]),axis=2)),axis=1)
            #print "--> 6: done"
    print "DONE! mdstepend:",mdstepend,"------------------------------------------------"
    sum_all = sum_all[:,:mdstepend]
    np.save("sum_all",sum_all)
    print "sum_all.shape:",sum_all.shape
    #for ind_qs,qs in enumerate(qs_all):
    #    np.savetxt('sum'+str(qpoint[0])+str(qpoint[1])+str(qpoint[2])+str(ind_qs)+'_'+str(mdstepend)+'.out',sum)
    return sum_all

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

def get_uncertainty_smoothed_powerspectrum(power_spectrum_smooth):
    ''' pass '''
    data = power_spectrum_smooth

    xvalues=data.shape[0]
    #print "xvalues:",xvalues
    datahalf = data[0:xvalues/2]
    #print "datahalf:",datahalf
    #print "datahalf.max():",datahalf.max()
    #print "datahalf.max()/2:",datahalf.max()/2

    def find_nearest(array,value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]


    x_ind_max = np.where(datahalf == datahalf.max())[0][0]  # 963  (9.63 THz)

    x_ind_max_over2_left  = find_nearest(datahalf[0:x_ind_max], datahalf.max()/2)
    x_ind_left = np.where(datahalf[0:x_ind_max] == x_ind_max_over2_left)[0][0]

    x_ind_max_over2_right = find_nearest(datahalf[x_ind_max:-1], datahalf.max()/2)
    x_ind_right = np.where(datahalf == x_ind_max_over2_right)[0][0]
    #print "x_ind_max(zuvor):",x_ind_max
    #print "x_ind_max_over2_left :", x_ind_max_over2_left, x_ind_left
    #print "x_ind_max_over2_right:", x_ind_max_over2_right, x_ind_right

    absmin = np.where(datahalf>0.45)[0][0]
    absmax = np.where(datahalf>0.45)[0][-1]


    vor_freq = x_ind_max/float(xvalues)
    vor_lifetime = (x_ind_right - x_ind_left)/float(xvalues)
    vor_lifetimemin = (absmax - absmin)/float(xvalues)

    #print "vor_freq:",vor_freq
    #print "vor_lifetime:",vor_lifetime
    #print "vor_lifetimemin:",vor_lifetimemin

    teiler = vor_freq/9.5   # hier gehen wir von der festen frequenz 9.5 aus
    #print "teiler (in):",teiler
    if teiler < .000005: sys.exit("no teiler oben,"+str(teiler))
    if teiler > .000005 and teiler < .00005:    teiler = .00001
    if teiler > .00005  and teiler < .0005:     teiler = .0001
    if teiler > .0005   and teiler < .005:      teiler = .001
    if teiler > .005    and teiler < .05:       teiler = .01
    if teiler > .05     and teiler < .5:        teiler = .1
    if teiler > .5      and teiler < 5:         teiler = 1.
    if teiler > 5       and teiler < 50:        teiler = 10.
    if teiler > 50      and teiler < 500:       teiler = 100.
    if teiler > 500     and teiler < 5000:      teiler = 1000.
    if teiler > 5000    and teiler < 50000:     teiler = 10000.
    if teiler > 50000   and teiler < 500000:    teiler = 100000.
    if teiler >                      500000: sys.exit("no teiler unten,"+str(teiler))
    #print "teiler (out):",teiler

    freq = vor_freq/teiler
    lifetime = vor_lifetime/teiler
    lifetimemin = vor_lifetimemin/teiler

    #print ""
    #print "Smoothed function:"
    print "Freq:",freq," THz, lifetime: ", lifetime," THz", "lifetime min:",lifetimemin
    return freq,lifetime, lifetimemin,xvalues

def sum_all_to_powerspectrum(sum_all, qpoint, mdstepstocalc = False):
    ''' calculates the power_spectrum for a certain qpoint (from sum_all)
    - sum_all contins for a single q point all (symmetry) equivalent q points
    - output of this skript is the power_spectrum for the corresponding q-point
    - qpoint is just for saving the power_spectrum
    '''
    mdstepstocalc_all = np.array([1000,5000,9000,10000,20000,50000,90000,100000,200000,500000,\
            900000,1000000,2000000,5000000,9000000,10000000,20000000,50000000,90000000,100000000])
    mdstepstocalc_all = np.sort(np.append(mdstepstocalc_all,[sum_all.shape[1]]))

    #print "mdstepstocalc:",mdstepstocalc
    mdstepstocalc_all = np.sort(np.unique(mdstepstocalc_all[np.where(mdstepstocalc_all <= sum_all.shape[1])[0]]))
    print "mdstepstocalc:",mdstepstocalc,type(mdstepstocalc)
    print "mdstepstocalc_all:",mdstepstocalc_all
    if type(mdstepstocalc) == bool:
        mdstepstocalc = mdstepstocalc_all

    print "mdstepstocalc_________________:",mdstepstocalc,type(mdstepstocalc)
    if type(mdstepstocalc) != bool:
        if type(mdstepstocalc) == list:
            mdstepstocalc = np.array(mdstepstocalc)
        elif type(mdstepstocalc) == int:
            mdstepstocalc = np.array([mdstepstocalc])
        elif type(mdstepstocalc) == np.ndarray:
            pass
    else:
        sys.exit("error: unknow type for mdstepstocalc")
        #print "mdstepstocalc_all:",mdstepstocalc_all

    print "mdstepstocalc:",mdstepstocalc
    mdstepstocalc = np.unique(mdstepstocalc/10*10) # remove elements wich are too close
    print "mdstepstocalc:",mdstepstocalc

    for idx,mdsteps in enumerate(mdstepstocalc):
        power_spectrum=np.zeros(mdsteps)
        print len(sum_all),"mdsteps:",mdsteps," = ",idx,"from",len(mdstepstocalc)
        for i in np.arange(len(sum_all)):
            #print "qpoint:",i+1,"of ",len(sum_all),"mdsteps:",mdsteps,"from",mdstepstocalc
            bla = sum_all[i][::-1][:mdsteps] # turn around data to get rid of equilibration problem
            ble=fft(bla[:mdsteps])
            power_spectrum+=np.abs(ble)**2.
        ps=power_spectrum/power_spectrum.max()*.9;


        # write out power_spectrum
        filenameout = "ps"+str(qpoint[0])+str(qpoint[1])+str(qpoint[2])+"_"+str(mdsteps)+".dat"
        np.savetxt(filenameout,ps)

        # write out power_spectrum_smoothed
        sigma=0.0005
        out1 = smoothing_power_spectrum(ps,sigma)
        np.savetxt(filenameout+".smooth_"+str(sigma)+".dat",out1)
        sigma=0.001
        out2 = smoothing_power_spectrum(ps,sigma)
        np.savetxt(filenameout+".smooth_"+str(sigma)+".dat",out2)

        # write out frequencies and lifetimes (for xmgrace)
        freq1, lifetime1, lifetimemin1, xv = get_uncertainty_smoothed_powerspectrum(out1)
        freq2, lifetime2, lifetimemin2, xv = get_uncertainty_smoothed_powerspectrum(out2)

        fq=np.array([freq1, freq2])
        lt=np.array([lifetime1,lifetime2])

        outlifetimes = np.array([[xv, lt.mean(), lt.max()-lt.mean(),lt.min()-lt.mean()]])
        outfreq = np.array([[xv, fq.mean(), fq.max()-fq.mean(),fq.min()-fq.mean()]])

        np.savetxt(filenameout+".lifetimes.dat",outlifetimes)
        np.savetxt(filenameout+".freqs.dat",outfreq)


if os.path.isfile('sum_all.npy'):
    print "load sum_all.npy"
    sum_all = np.load('sum_all.npy')
else:
    sum_all = pos_compress_to_sum_equivalent_qs(\
        filename,\
        qpoint=np.array([N,0,0]),\
        lines_of_one_mdstep=N**3*4,\
        mdsteps_per_chunk = 10000, \
        scale_data=float(alat)*N)

sum_all_to_powerspectrum(sum_all, qpoint=np.array([N,0,0]))
