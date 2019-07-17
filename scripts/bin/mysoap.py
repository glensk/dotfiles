#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys
import quippy

##################################################################################
## SOAP funcions
##################################################################################
def get_rawsoap(at, gs=0.3, co=3.0, cotw=0.5, nmax=9, lmax=6, cw=1.0, nrm=True, zlist=[], central_z=0, x_cmd=""):
    """
    gs gaussian widths
    co cutoff radius
    cotw you don't care
    n radial basis functions
    l spherical harmonics basis
    zlist list of atomic numbers the SOAPS will see
    central_z list of atomic numbers you will center on : IT DEFINES THE NUMBER OF ROWS of the rawsoaps

    """
    if len(zlist)==0:
        zlist = np.unique(np.asarray(at.z))
    nz = len(zlist)
    lspecies = "n_species="+str(nz)+" species_Z={"+(" ".join(np.asarray(zlist,dtype=str)))+"}"
    soapstr=("soap central_reference_all_species=F normalise="+("T" if nrm else "F")+" central_weight="+str(cw)+
               "  covariance_sigma0=0.0 atom_sigma="+str(gs)+" cutoff="+str(co)+" cutoff_transition_width="+str(cotw)+
               " n_max="+str(nmax)+" l_max="+str(lmax)+' '+lspecies+' Z='+str(central_z)+' '+x_cmd)
    desc = quippy.descriptors.Descriptor(soapstr )
    at.set_cutoff(co);
    at.calc_connect();
    return zlist, desc.calc(at)["descriptor"]

def get_soap2_vec(rawsoap, zlist, zglobal=None,nmax=None,lmax=None):
    if nmax is None:
        sys.exit('Please provide nmax!')
    if lmax is None:
        sys.exit('Please provide lmax!')
    isoap = 0
    isqrttwo = 1.0/np.sqrt(2.0)
    njsoap = {}; ipair = {}

    #print('zglobal:',zglobal)
    #print('zlist  :',zlist,'ka')
    if zglobal is None:
        zglobal = zlist
    #print('len(zglobal)',len(zglobal),nmax,len(zglobal),lmax)
    njsoap = np.zeros((len(zglobal),nmax,len(zglobal),nmax,lmax+1))
    zmap = np.zeros(len(zlist),int)
    for s in xrange(len(zlist)):
        zmap[s] = zglobal.index(zlist[s])
    for s1 in xrange(len(zlist)):
        z1 = zmap[s1]
        for n1 in xrange(nmax):
            for s2 in xrange(s1+1):
                z2 = zmap[s2]
                for n2 in xrange(nmax if s2<s1 else n1+1):
                    soap_element = rawsoap[isoap:isoap+lmax+1]
                    if (s1 != s2 or n1 != n2):
                        soap_element *= isqrttwo  # undo the normalization since we will actually sum over all pairs in all directions!
                    njsoap[z1,n1,z2,n2,:] = soap_element
                    njsoap[z2,n2,z1,n1,:] = soap_element
                    isoap+=lmax+1
    return njsoap

def do_fps(x, d=0):
    if d == 0 : d = len(x)
    n = len(x)
    iy = np.zeros(d, int)
    # faster evaluation of Euclidean distance
    n2 = np.sum(x**2,axis=1)
    iy[0] = np.random.randint(0, n)
    dl = n2 + n2[iy[0]] - 2* np.dot(x, x[iy[0]])
    lmin = np.zeros(d)
    for i in range(1,d):
        iy[i] = np.argmax(dl)
        lmin[i-1] = dl[iy[i]]
        if i%1000 ==0: print("max min dist %f" %( dl[iy[i]]))
        nd = n2 + n2[iy[i]] - 2*np.dot(x,x[iy[i]])
        dl = np.minimum(dl, nd)
    return iy, lmin

def get_soaps(kmcxyz = False, nmax = 8, lmax = 6, co = 4, gs = 0.5, zlist = [12,13,14], central_z=23,verbose=False,showtdqm=True,nrm=True ):
    rsoap = []
    #nmax = 8   # 12
    #lmax = 6  # 9
    #co = 4  # cutoff radius, 3.0  ## is this in angstrom?
    #gs = 0.5 # gaussian width, 0.3
    cotw = gs
    det = "c"+str(co)+"-g"+str(gs)+"-nrm-cw"
    #zlist = [12,13,14]  # why is this working, even without specifying the vacancy/Vanadium ....? [12,13,14,23]
    ztot = len(zlist)
    nsoap = ztot**2 * nmax**2 * (1+lmax)
    ntot = len(kmcxyz)
    if verbose:
        print('rsoap ...')
    if showtdqm == True:
        gothrough = tqdm_notebook(kmcxyz)
    else:
        gothrough = kmcxyz

    for at in gothrough:
        #progress(i,len(kmcxyz))
        # andrea used a zentral_z=0
        # the next line could be potentially done parallel...
        zliat, soaps = get_rawsoap(at, nmax=nmax, lmax=lmax, co=co, gs=gs, cotw=cotw, nrm=nrm, cw=1, zlist=zlist, central_z=central_z,
                                   x_cmd="cutoff_dexp=2 cutoff_scale=3.0 cutoff_rate=2.0")
        rsoap.append((zliat, soaps))
    if verbose:
        print('rsoap done ...')

    soap2 = []
    if verbose:
        print('soap2 ...')
    for i in tnrange(ntot):
        #progress(i,ntot)
        izl, soaps = rsoap[i]
        lenv = np.zeros((len(soaps), ztot, nmax, ztot, nmax, 1+lmax))
        for s in range(len(soaps)):
            lenv[s] = get_soap2_vec(soaps[s], izl, zlist,nmax=nmax,lmax=lmax)
        soap2.append(lenv)
    if verbose:
        print('soap2 done ...')
        print('esoap2 ...')
    esoap2 = []
    for s in soap2:
        #progress(i,len(soap2))
        for e in s:
            esoap2.append(e.flatten())
    if verbose:
        print('esoap2 done ...')
    esoap2 = np.asarray(esoap2)
    return rsoap, esoap2


if __name__ == '__main__':
    pass
