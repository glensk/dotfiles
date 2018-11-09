#!/usr/bin/python
import argparse
import re
import subprocess,os
import numpy as np
import time
import scipy.linalg as salg
import scipy.sparse.linalg as spalg
import sys


zmap = {"H": 1,"He": 2,"Li": 3,"Be": 4,"B": 5,"C": 6,"N": 7,"O": 8,"F": 9,"Ne": 10,"Na": 11,"Mg": 12,"Al": 13,"Si": 14,"P": 15,"S": 16,"Cl": 17,"Ar": 18,"K": 19,"Ca": 20,"Sc": 21,"Ti": 22,"V": 23,"Cr": 24,"Mn": 25,"Fe": 26, "Co": 27, "Ni": 28}

def run(cmd, logfile):
    """ To quickly run bash commands from the code for prototyping purposes"""
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
    return p

def CURSelOrtho(cov, numSym, costs=1):
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
    if numSym == 1: return [sel]
    sys.stderr.write("%d more columns to select                    \r" %(numSym) )

    # Calling this function in an iterative way is nice and clean but might be a problem when using too many symmetry functions
    return np.concatenate(([sel],CURSelOrtho(ocov, numSym-1,costs)))

# Functions to read and write symmetry function definitions in a format compatible with RuNNer
def ReadDEF(deffile):
    rsym = []
    with open(deffile,'r') as sf:
        for line in sf:
            rsym.append(line.split()[1:])
    return rsym

def WriteDEF(defs, outstr=sys.stdout):
    for d in defs:
        print >>outstr, "symfunction_short ",
        for s in d:
            print >>outstr, s, " ",
        print >>outstr, ""
    return

def Preproc(datafile, element, mode):
    """ Reads chunks of function.data and appends it to preprocessed datafile"""
    if mode=='fast':
        dataf = []
        nspecies = []
        nframe = 0
        counter = 0
        chunk = 5000
        o = open(element+'_preproc.data','wb')
        with open(datafile,'r') as df:
            for line in df:
                if len(line.split())==1:
                    nspecies.append(0)
                    nframe += 1
                if line.split()[0]==str(zmap[element]) and len(line.split()) > 2:
                    tmp = np.fromstring(line, dtype=float, sep=' ')[1:]
                    dataf.append(tmp)
                    nspecies[nframe - 1] += 1
                counter+=1
                if counter == chunk:
                    np.savetxt(o,dataf,fmt='%.10f')
                    dataf = []
                    counter = 0

        np.savetxt(o, dataf, fmt='%.10f')
        np.savetxt(element+'_species_in_frame', nspecies, fmt='%d')
        print "Saved preconditioned datafile for {}".format(element)
        o.close()

    return

def ReadCov(filename, nsym=100, bufsize=10000):
    f = open(filename, "r")
    vvt = np.zeros((nsym,nsym),float)
    vm  = np.zeros((nsym),float)
    ntot = 0
    i=0
    while True:
        chunk = np.fromfile(f, dtype="float",count=(nsym)*bufsize, sep=" ")
        if len(chunk) ==0: break
        nrow = len(chunk)/(nsym)
        data = chunk.reshape((nrow,nsym))
        vm += data.sum(axis=0)
        vvt += np.dot(data.T,data)
        ntot +=nrow
        i+=1

    return (vvt - np.outer(vm, vm/ntot))/ntot

def ReadLOG(logfile):
    reheader = re.compile('short range atomic symmetry functions element *([a-zA-Z]+)')
    sf = open(logfile,'r')
    fundef = {}
    while True:
        line = sf.readline()
        if line == "":
            break
        rehm = reheader.search(line)
        if not rehm is None:
            print "Parsing symmetry functions for atom", rehm.group(1)
            sf.readline() # parse comment
            rsym = []
            line = sf.readline()
            while line[2] != "-":
                rsym.append(line.split()[1:])
                line = sf.readline()
            fundef[rehm.group(1)] = rsym
            print "READ "
    sf.close()
    return fundef

def GetCosts(sym):
    radius = np.asarray([ s[-1] for s in sym ], float)
    # radial syms have a cost scaling with the cube of radius,
    # angular use triplets so it's the sixth power
    stype = np.asarray([ ( 3 if s[1] == '2' else 6)  for s in sym ], float)

    return radius**stype

def GetCov(data, scaling):
    if hasattr(scaling, "__len__") and len(scaling) != len(data[0]):
        raise ValueError("Inconsistent size of the scaling vector")

    avg = data.mean(axis=0)
    scov =  (np.dot(data.T,data) - np.outer(avg,avg) * len(data)) * np.outer(scaling,scaling)
    return scov / len(data)

def GetLandmarks(nsyms, nlandmarks):
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

        framelen = np.loadtxt(element+'_species_in_frame', int) # reads the number of times this element appears in the frames
        nframes = len(framelen)
        if nlandmarks > nframes:
            print "You requested {:d} landmarks and there are {:d} frames. You can use the whole input.data or re-run CurSel requesting fewer landmarks".format(nlandmarks, nframes)
            return range(nframes)
        if natoms is None:
            natoms = framelen.copy()
        else: natoms += framelen
        f = open(element+'_preproc.data', "r") # opens the file containing the symmetry function data for this element
        avg_sf[element] = np.zeros((len(framelen),nsym))

        iframe = 0
        for nrow in framelen:
            if nrow == 0: # this element is missing in this frame, so we just add a block of zeros
                avg_sf[element][iframe,:]=np.zeros(nsym)
            else:
                chunk = np.fromfile(f, dtype="float",count=nsym*nrow, sep=" ")

                if len(chunk) ==0:
                    raise ValueError("End of file while reading element "+element)
                data = chunk.reshape((nrow,nsym))
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
    distances = [0]
    mindist = np.zeros(nframes)
    for k in xrange(nframes):
        mindist[k] = np.linalg.norm(frame_sf[k] - frame_sf[sel_frames[0]])
    while len(sel_frames) < nlandmarks-1:
        sel_frames.append(np.argmax(mindist))
        distances.append(np.max(mindist))
        sys.stdout.write("\rSelecting frame:\t%d" % (len(sel_frames)+1))
        sys.stdout.flush()
        for k in xrange(nframes):
            mindist[k] = min(mindist[k],np.linalg.norm(frame_sf[k] - frame_sf[sel_frames[-1]]) )
    sel_frames.append(np.argmax(mindist))
    distances.append(np.max(mindist))
    sys.stdout.write("\n")
    #print('dist',distances,type(distances))
    np.savetxt('cursel.distances.dat',np.array(distances))

    return sel_frames


def SymSelect(datafile, deffile, threshold = 1e-5, prefix="cursel", nlandmarks=0, mode="fast", restart=None):

    t1=time.time()
    verbose = False #!TODO Add the option

    threshold = float(threshold)
    nlandmarks = int(nlandmarks)
    # Read Definitions from log files
    asdef = ReadLOG(deffile)
    nsyms = dict.fromkeys(asdef.keys())
    for element in asdef:
        nsyms[element] = len(asdef[element])
    print "Number of symmetry functions for each element: {}".format(nsyms.items())

    if restart is None:
        for element in nsyms:
            Preproc(datafile,element,mode)
    	    preproc = element+'_preproc.data'
            dcov = ReadCov(preproc, nsyms[element], bufsize=10000)
            np.savetxt(prefix+'.cov_'+element, dcov)
            restart = prefix
    else:
        print "Restarting from pre-calculated covariance"


    sel1 = {}
    for element in nsyms:
        dcov = np.loadtxt(restart+'.cov_'+element)
        # estimate number of symfuncts based on threshold
        evc = np.linalg.eigvalsh(dcov)[::-1]
        evc /= evc[0]  # scales eigenvectors
        nH = evc[np.where(evc > threshold)]

        # compute cost and select the symmetry functions
        costs = GetCosts(asdef[element])
        if verbose:
            np.savetxt(element+'.costs', costs)
        print 'Starting the CUR selection for {}'.format(element)
        sel1[element] = CURSelOrtho(dcov,len(nH),costs)

    fout = open(prefix+".def","w")
    for element in asdef:
        WriteDEF([asdef[element][i] for i in sel1[element]], fout)

    # if required compute landmarks based on the selected symmetry functions
    if nlandmarks > 0 :
        sel_landmarks = GetLandmarks(nsyms, nlandmarks)
        np.savetxt(prefix+'.landmarks', sel_landmarks, fmt='%d')
    t2=time.time()
    print 'It took {:.5f} seconds to run'.format((t2-t1))

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument("datafile", help="Typically the function.data that has been obtained after running mode 1 on RuNNer")
    parser.add_argument("deffile", help="Typically the log_file obtained after running mode 1 on RuNNer. The input.nn cannot be used as the symmetry functions are listed differently")
    parser.add_argument("-t", type=float,default=1e-5, help="Threshold (approximate) for selecting symmetry functions based on the eigenspectrum of the covariance")
    parser.add_argument("--prefix", type=str,default="cursel", help="Prefix for all output filenames")
    parser.add_argument("--landmarks", type=int,default=0, help="Also print out indices of the N most representative structures")
    parser.add_argument("--mode", type=str,default="fast",help="How to load the symmetry functions data, if function.data can't be loaded on RAM go for -slow-")
    parser.add_argument("--restart", type=str, default=None, help="If this option is enabled then it will use the covariance dumped from a previous run to re-calculate the symmetry functions. The path to the covariance file must be provided WITHOUT the final .cov_[] part.")
    args = parser.parse_args()
    SymSelect(args.datafile, args.deffile, threshold=args.t, prefix=args.prefix, nlandmarks=args.landmarks, mode=args.mode, restart=args.restart)
