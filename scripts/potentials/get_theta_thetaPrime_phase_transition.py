#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse,glob

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    #p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def todo(args):
    #if not os.path.isfile(args.inputfile):
    #    sys.exit("inputfile "+args.inputfile+" does not exist")
    folder = glob.glob("n2p2_alcu_v2dm_*")
    tall=[]
    for f in folder:
        #print(f)
        t = np.loadtxt(f+"/Theta/fah/Fah_surfacefit_upto1000/thermo_3rd/output_0.0001GPa/Gibbs_energy")
        tp = np.loadtxt(f+"/ThetaPrime/fah/Fah_surfacefit_upto1000/thermo_3rd/output_0.0001GPa/Gibbs_energy")
        #print(f,len(t),len(tp))
        t=t[:1000]
        tp=tp[:1000]
        #t[:,1] = np.array(t-tp)
        #print(t)
        #print()
        #print(t-tp)
        #sys.exit()
        a=t[:,1]-tp[:,1]
        T=np.where(a<0)[0][0]
        print(f,'T',T)
        tall.append(T)
    tall = np.array(tall)
    print('mean',tall.mean())
    print('std',tall.std())

    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    todo(args)
    #folder = glob.glob("n2p2_alcu_v2dm_*")
    #for f in folder:
    #    print(f)
    #    t = np.loadtxt(f+"/Theta/fah/Fah_surfacefit_upto1000/thermo_3rd/output_0.0001GPa/Gibbs_energy")
    #    tp = np.loadtxt(f+"/ThetaPrime/fah/Fah_surfacefit_upto1000/thermo_3rd/output_0.0001GPa/Gibbs_energy")
    #    #t[:,1] = np.array(t-tp)
    #    print(t)
    #    print()
    #    print(t[:,1]-tp[:,1])
    #    sys.exit()


