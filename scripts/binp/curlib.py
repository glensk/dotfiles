#!/usr/bin/env python
import argparse
import re
import subprocess,os
import numpy as np
import time
import scipy.linalg as salg
import scipy.sparse.linalg as spalg
import sys
import scipy.integrate as spint
import operator
from iolib import *

def cCURSelOrtho(cov, numSym, costs=1):
    """ Apply (deterministic) CUR selection of numSymm rows & columns of the
    given covariance matrix, including an orthogonalization step. Costs can be weighted if desired. """

    evc,evec = spalg.eigs(cov,numSym)
    weights = np.sum(np.square(evec),axis=1)/costs
    sel = np.argmax(weights)
    vsel = cov[sel].copy()
    vsel *= 1.0/np.sqrt(np.dot(vsel,vsel))
    ocov = cov.copy()
    for i in xrange(len(ocov)):
        # Diagonalize the covariance matrix wrt the chosen column
        ocov[i] -= vsel * np.dot(cov[i],vsel)

    return sel, ocov

def CURSelOrtho(cov, numSym, costs=1):

    ocov = cov.copy()
    rsel = np.zeros(numSym, int)

    for i in xrange(numSym):
        rval, ocov = cCURSelOrtho(ocov, numSym-i, costs)
        rsel[i] = rval
    return rsel

# Orthogonal CUR in all its glory, with SVD flavor
def cCURSelSVDCols(M, nsvd, costs=1, verbose=False):
   """ Apply (deterministic) CUR selection of numSymm rows & columns of the
   given feature matrix, including an orthogonalization step. Costs can be weighted if desired. """
   U,S,VT = spalg.svds(M,nsvd)
   weights = np.sum(np.square(VT),axis=0)/costs
   sel = np.argmax(weights)
   vsel = M[:,sel].copy()
   if verbose:
       print "svd ", VT.shape, np.linalg.norm(vsel), np.sum(weights),  weights[sel], weights.max(), weights.min()
   vsel *= 1.0/np.sqrt(np.dot(vsel,vsel))
   OM = M.copy()
   for i in xrange(len(M[0])):
       OM[:,i] -= vsel * np.dot(M[:,i],vsel)
   if verbose:
       print "check: ", np.linalg.norm(OM[:,sel]), np.sqrt((OM**2).sum())

   return sel, OM

def CURSelSVDCols(M, numSym, costs=1, nsvd = 0, verbose=False):
    print "CURSVD", M.shape
    OM = M.copy()
    rsel = np.zeros(numSym, int)

    for i in xrange(numSym):
        if nsvd == 0:
            nsym = numSym-i
        else:
            nsym = nsvd
        tic = time.time()
        rval, OM = cCURSelSVDCols(OM, nsym, costs, verbose)
        toc = time.time()
        sys.stdout.write("Selected %d. It took %.5f s for the last step.    Remaining lines: %d            \n" % (rval, toc-tic, numSym-i) )

        if verbose:
            print "IS IT ZERO???? ", np.sum(np.square(OM)[:,rval])
        rsel[i] = rval
    return rsel


def CURSelSVDColsAll(Ms, numSym, costs):
    print "CURSVD - all elements"

    OM = {}
    rsel = {}

    for el in Ms:
        OM[el] = Ms[el].copy()
        rsel[el] = []

    for i in xrange(numSym):
        tic = time.time()
        weights = {}
        mxweight = {}
        for el in Ms:
            U,S,VT = spalg.svds(OM[el],numSym-i)
            weights[el] = np.sum(np.square(VT),axis=0)/costs[el]
            mxweight[el] = weights[el].max()
        #finds the best element
        elsel = max(mxweight.iteritems(), key=operator.itemgetter(1))[0]
        # then runs orthogonalization for that element
        sel = weights[elsel].argmax()
        vsel = OM[elsel][:,sel].copy()
        vsel *= 1.0/np.sqrt(np.dot(vsel,vsel))
        for k in xrange(len(OM[elsel][0])):
            OM[elsel][:,k] -= vsel * np.dot(OM[elsel][:,k],vsel)
        rsel[elsel].append(sel)
        toc = time.time()
        sys.stdout.write("Selected %d for element %s. It took %.5f s for the last step.    Remaining lines: %d            \n" %(sel, elsel, toc-tic, numSym-i-1))

    return rsel

def DoCUR(M, cols, rows, rcond=1e-15):
    C = M[:,cols]
    R = M[rows,:]
    Ci = np.linalg.pinv(C, rcond=rcond)
    Ri = np.linalg.pinv(R, rcond=rcond)
    if M.shape[0]>M.shape[1]:
        U = np.dot(np.dot(Ci,M),Ri)
    else:
        U = np.dot(Ci,np.dot(M,Ri))
    if U.shape[0]>U.shape[1]:
        RM = np.dot(np.dot(C,U),R)
    else:
        RM = np.dot(C,np.dot(U,R))
    return RM


def pinv_distance(red_cov,cov):
   mX = np.dot(nalg.pinv(red_cov,rcond=1e-15),cov)
   zz = (np.dot(red_cov,mX) + np.dot(red_cov,mX).T)/2.
   return nalg.norm(cov - zz,'fro')

def GetCov(data, scaling):
    if hasattr(scaling, "__len__") and len(scaling) != len(data[0]):
        raise ValueError("Inconsistent size of the scaling vector")

    avg = data.mean(axis=0)
    scov =  (np.dot(data.T,data) - np.outer(avg,avg) * len(data)) * np.outer(scaling,scaling)
    return scov / len(data)

def GetLandmarks(datafile, nsyms, nat_per_frame, nlandmarks):
    #!TODO nsyms should be discarded, like also nat_per_frame. Both should be read within the function, the only arguments passed
    # should be datafile, logfile, and nlandmarks
    """ This function computes the distance between the frames
    First it reduces the amount of data from function.data in a single array for every frame
    Then it calculates the furthest frames from the ones already selected until the requested number of landmarks is hit. """

    avg_sf = dict.fromkeys(nsyms.keys())

    print 'Preprocessing for landmark selection'
    # reads in symmetry functions from the preprocessed data and computes an "average fingerprint" for the atoms of each specie
    totsym = 0
    natoms = None
    for element in nsyms:
        nsym = nsyms[element]
        totsym += nsym

        nframes = len(nat_per_frame[element])
        if nlandmarks > nframes:
            print "You requested {:d} landmarks and there are {:d} frames. You can use the whole input.data or re-run CurSel requesting fewer landmarks".format(nlandmarks, nframes)
            return range(nframes)
        if natoms is None:
            natoms = np.array(nat_per_frame[element])
        else: natoms += nat_per_frame[element]
        xmat = ReadFuncdata(datafile, element)
        avg_sf[element] = np.zeros((len(nat_per_frame[element]),nsym))

        iframe = 0
        past_row = 0
        for nrow in nat_per_frame[element]:
            if nrow == 0: # this element is missing in this frame, so we just add a block of zeros
                avg_sf[element][iframe,:]=np.zeros(nsym)
            else:
                data = xmat[past_row:past_row+nrow]
                past_row += nrow
                avg_sf[element][iframe,:]=np.sum(data,axis=0)
            iframe += 1

    # now collates the descriptors for each element, and normalizes properly so we can use the distance between
    # the compounded fingerprints as an indicator of the overall nature of each frame
    frame_sf = np.zeros((nframes, totsym))

    ksym = 0
    for element in nsyms:
        nsym = nsyms[element]
        for k in range(nframes):
            frame_sf[k,ksym:ksym+nsym] = avg_sf[element][k,:] / natoms[k]
        ksym += nsym

    print "Selecting the %d most important frames" % nlandmarks
    #This computes the distance and outputs the selected frames
    sel_frames = [0]
    mindist = np.zeros(nframes)
    for k in xrange(nframes):
        mindist[k] = np.linalg.norm(frame_sf[k] - frame_sf[sel_frames[0]])
    while len(sel_frames) < nlandmarks-1:
        sel_frames.append(np.argmax(mindist))
        sys.stdout.write("Selecting frame:\t%d" % (len(sel_frames)+1))
        sys.stdout.flush()
        for k in xrange(nframes):
            mindist[k] = min(mindist[k],np.linalg.norm(frame_sf[k] - frame_sf[sel_frames[-1]]) )
    sel_frames.append(np.argmax(mindist))
    sys.stdout.write("\n")

    return sel_frames
