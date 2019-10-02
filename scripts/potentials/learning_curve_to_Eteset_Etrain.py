#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-i','--inputfile', required=True, type=str,default="learning-curve.out", help="name of the inputfile inputfile")
    return p

def todo(args):
    if not os.path.isfile(args.inputfile):
        sys.exit("inputfile "+args.inputfile+" does not exist")
    lc = np.loadtxt(args.inputfile)

    lc[:,1] = lc[:,1]*1000.*27.211384
    lc[:,2] = lc[:,2]*1000.*27.211384
    lc[:,3] = lc[:,3]*1000.*51.422063
    lc[:,4] = lc[:,4]*1000.*51.422063

    Etrain = lc[:,[0,1]]
    Etest  = lc[:,[0,2]]
    Ftrain = lc[:,[0,3]]
    Ftest  = lc[:,[0,4]]
    np.savetxt("Etrain.dat",Etrain)
    np.savetxt("Etest.dat",Etest)
    np.savetxt("Ftrain.dat",Ftrain)
    np.savetxt("Ftest.dat",Ftest)
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    todo(args)


