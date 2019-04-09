#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse
import myutils

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    #p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def todo(args):
    datafile_try = [ 'data.lmp', 'data.runnerformat.lmp' ]
    datafile = False
    for i in datafile_try:
        if os.path.isfile(i):
            datafile = i
    if datafile == False:
        for i in datafile_try:
            print("lammps datafile",i,"does not exist")
        sys.exit()

    N = 10
    with open(datafile) as myfile:
        head = [next(myfile) for x in xrange(N)]
        #print(head)
    for i in head:
        print(i.rstrip())
        if 'xlo' in i:
            print('i')
            a=i.split()[1]
        if 'ylo' in i:
            b=i.split()[1]
        if 'zlo' in i:
            c=i.split()[1]
        if 'xy xz yz' in i:
            d=i.split()[0]
            e=i.split()[1]
            f=i.split()[2]
    print()
    print('-->a',a,b,c,'def',d,e,f)
    print()
    print('-->x',a,"0","0",d,b,"0",d,f,c)
    print()
    print('cp...')
    myutils.cp("simulation.pos_0.xyz","simulation.pos_libatom.xyz")
    print('sed...')
    #myutils.sed("simulation.pos_libatom.xyz",'#.*',"Lattice="+str(a)+" 0 0 "+str(d)+" "+str(b)+" 0 "+str(d)+" "+str(f)+" "+str(c))  # this is too slow

    #if not os.path.isfile(args.inputfile):
    #    sys.exit("inputfile "+args.inputfile+" does not exist")
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    todo(args)


