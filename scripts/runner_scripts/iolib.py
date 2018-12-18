#!/usr/bin/python
import argparse
import re
import subprocess,os
import numpy as np
import time
import scipy.linalg as salg
import scipy.sparse.linalg as spalg
import sys
import scipy.integrate as spint


zmap = {"H": 1,"He": 2,"Li": 3,"Be": 4,"B": 5,"C": 6,"N": 7,"O": 8,"F": 9,"Ne": 10,"Na": 11,"Mg": 12,"Al": 13,"Si": 14,"P": 15,"S": 16,"Cl": 17,"Ar": 18,"K": 19,"Ca": 20,"Sc": 21,"Ti": 22,"V": 23,"Cr": 24,"Mn": 25,"Fe": 26, "Ni": 28,"Ga": 31, "As":33}

def run(cmd, logfile):
    """ To quickly run bash commands from the code for prototyping purposes"""
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile)
    return p

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
        for s in d[:-3]:
            print >>outstr, s, " ",
        print >>outstr, ""
    return

def ReadFuncdata(funcfile, element="H", verbose=False):
    environs = []
    nat = 0
    nfr = 0
    zel = zmap[element]
    with open(funcfile) as f:
        for line in f:
            ll = line.split()
            if len(ll) == 1: # start new frame
                nat = 0
            elif len(ll) == 4:
                nfr += 1
                if nfr%100 ==0 and verbose: print(nfr)
            else:
                idx = int(ll[0])
                if idx == zel:
                    nat += 1
                    environs.append(np.asarray(ll[1:], float))
    return np.asarray(environs)

def ReadDensity(indata):
    iframes = []
    nfr = 0
    nat = {}
    tvol = 0
    nat_per_frame = {}
    old_nat = {}
    with open(indata) as inf:
        for l in inf:
            if l.strip() == "begin":
                h = np.zeros((3,3),float)
                nlatt = 0
                pass
            elif l.strip() == "end":
                for el in nat:
                    nat_per_frame[el].append(nat[el]-old_nat[el])
                    old_nat[el] = nat[el]
                tvol += np.abs(np.linalg.det(h))
                nfr +=1
            else:
                sl = l.split()
                if sl[0] == "lattice":
                    h[nlatt] = np.asarray(sl[1:],float)
                    nlatt += 1
                elif sl[0] == "atom":
                    if not sl[4] in nat:
                        nat[sl[4]] = 0
                        nat_per_frame[sl[4]] = []
                        old_nat[sl[4]] = 0
                    nat[sl[4]] += 1
    dens = {}
    for el in nat:
        dens[el] = nat[el] / tvol
    return dens, nat_per_frame

def ReadLOG(logfile):
    reheader = re.compile('Short range atomic symmetry functions element *([a-zA-Z]+)')
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
            sf.readline() # parse comment
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
