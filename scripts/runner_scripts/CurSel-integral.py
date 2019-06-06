#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import re
import subprocess,os
import numpy as np
import time
import scipy.linalg as salg
import scipy.sparse.linalg as spalg
import sys
import scipy.integrate as spint
from iolib import *
from curlib import *
from sflib import *

def SymSelect(datafile, deffile, nsymsel = 20, prefix="cursel", nlandmarks=0, verbose=False,  allel=False, scalemode="none", density=0.0, nsvd=0):

    t1=time.time()

    nlandmarks = int(nlandmarks)
    # Read Definitions from log files
    asdef = ReadLOG(deffile)
    nsyms = dict.fromkeys(asdef.keys())
    SF_int = dict.fromkeys(asdef.keys())
    rho = {}
    if density == 0.0: # automatic determination of average atom density
        rho,nat_per_frame = ReadDensity("input.data")
        if verbose:
            print nat_per_frame
    else:
        for el in asdef:
            rho[el] = density
    for element in asdef:
        nsyms[element] = len(asdef[element])
        if not scalemode == "none" :
            print "Now evaluating the integrals for each symmetry function"
            try:
                SF_int[element] = np.loadtxt(prefix+"_"+element+".sfintegral")
                print "Loaded integrals file for ", element
            except:
                SF_int[element] = SF_integrate(asdef[element], rho)
                np.savetxt(prefix+"_"+element+".sfintegral", SF_int[element])
    print "Number of symmetry functions for each element: {}".format(nsyms.items())

    sel1 = {}
    errcur = {}
    if allel:
        xmats = {}
        costs = {}

        for element in nsyms:
            xmats[element] = ReadFuncdata(datafile, element, verbose)
            if scalemode == "full":
                xmats[element] *= (1.0/SF_int[element])
                print "Scaling by full integral"
            elif scalemode == "sqrt":
                print "Scaling by sqrt integral"
                xmats[element] *= np.sqrt(1.0/SF_int[element])
            costs[element] = GetCosts(asdef[element], rho)

        sel1 = CURSelSVDColsAll(xmats, nsymsel*len(nsyms), costs)
        if verbose:
            for element in sel1:
                curx = DoCUR(xmats[element], sel1[element], np.asarray(range(len(xmats[element])),int) )
                errcur[element] = np.sqrt(np.sum((curx-xmats[element])**2)/np.sum(xmats[element]**2))
                print "CUR error: ", element,  errcur[element]
    else:
        for element in nsyms:
            xmat = ReadFuncdata(datafile, element)
            if scalemode == "full":
                xmat *= (1.0/SF_int[element])
            elif scalemode == "sqrt":
                print "Scaling by sqrt integral"
                xmat *= np.sqrt(1.0/SF_int[element])
            xmat/=np.max(np.abs(xmat)) #normalizes xmat to avoid very large or very small values

            # compute cost and select the symmetry functions
            costs = GetCosts(asdef[element], rho)

            # speeds up things by removing SF that are effectively zero
            tot = np.sum(xmat**2,axis=0)
            gthanzero = np.where(tot>1e-12)[0]
            print 'Starting the CUR selection for {}'.format(element)
            cursel = CURSelSVDCols(xmat[:,gthanzero],nsymsel,costs[gthanzero], nsvd, verbose)
            sel1[element] = gthanzero[cursel]
            if verbose:
                np.savetxt(element+'.costs', costs)
                print "Computing CUR error"
                curx = DoCUR(xmat, sel1[element], np.asarray(range(len(xmat)),int))
                errcur[element] = np.sqrt(np.sum((curx-xmat)**2)/np.sum(xmat**2))
                print "CUR error: ",  errcur[element]

    fout = open(prefix+".def","w")
    if verbose:
        fout.write("# Command line: %s\n" % (" ".join(sys.argv[:]) ) )
    for element in asdef:
        if verbose:
            fout.write("# CUR relative error for %s: %15.8e \n" % (element, errcur[element]) )
        WriteDEF([asdef[element][i] for i in sel1[element]], fout)

    # if required compute landmarks based on the selected symmetry functions
    if nlandmarks > 0 :
        sel_landmarks = GetLandmarks(datafile, nsyms, nat_per_frame, nlandmarks)
        np.savetxt(prefix+'.landmarks', sel_landmarks, fmt='%d')
    t2=time.time()
    print 'It took {:.5f} seconds to run'.format((t2-t1))

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("datafile", help="Typically the function.data that has been obtained after running mode 1 on RuNNer")
    parser.add_argument("deffile", help="Typically the log_file obtained after running mode 1 on RuNNer. The input.nn cannot be used as the symmetry functions are listed differently")
    parser.add_argument("-n", type=int,default=20, help="Number of symmetry functions to be selected")
    parser.add_argument("--prefix", type=str,default="cursel", help="Prefix for all output filenames")
    parser.add_argument("--landmarks", type=int,default=0, help="Also print out indices of the N most representative structures")
    parser.add_argument("--nsvd", type=int,default=0, help="Runs CUR with a fixed number of SVD components (default:auto)")
    parser.add_argument("-v", default=False, action='store_true', help="Compute the Frobenius norm relative error in the CUR approximated matrix")
    parser.add_argument("--allelements", default=False, action='store_true', help="Chooses SFs across elements")
    parser.add_argument("--rho", type=float,default=0.0,help="Approximate density of the reference gas (to normalize 2 and 3-body functions) in bohr**-3")
    parser.add_argument("--scale", type=str, default="none", help="Scaling of the SF based on the integral value. Permissible values are [none, sqrt, full]")
    args = parser.parse_args()
    SymSelect(args.datafile, args.deffile, nsymsel=args.n, prefix=args.prefix,
           nlandmarks=args.landmarks, verbose=args.v, allel = args.allelements, scalemode=args.scale,
           density=args.rho, nsvd=args.nsvd)
