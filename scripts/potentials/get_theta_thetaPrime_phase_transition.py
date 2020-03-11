#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse,glob
import matplotlib.pyplot as plt

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
    Ball=[]
    bdall=[]
    alle = np.zeros((40,8))
    for idx,f in enumerate(folder):
        print(f)
        nr = int(f.split("n2p2_alcu_v2dm_")[1])
        evt = np.loadtxt(f+"/Theta/evinet/Evinet_1")
        Bt = evt[2]
        Bdt = evt[3]
        evtp = np.loadtxt(f+"/ThetaPrime/evinet/Evinet_1")
        Btp = evtp[2]
        Bdtp = evtp[3]
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
        print(f,'T',T-595.7,'dB',Bt-Btp,'dB',Bdt-Bdtp)
        tall.append(T)
        alle[nr-1][0] = nr
        alle[nr-1][0+1] = T
        alle[nr-1][1+1] = Bt
        alle[nr-1][2+1] = Btp
        alle[nr-1][3+1] = Btp-Bt
        alle[nr-1][4+1] = Bdt
        alle[nr-1][5+1] = Bdtp
        alle[nr-1][6+1] = Bdtp-Bdt
    tall = np.array(tall)
    print('mean',tall.mean())
    print('std',tall.std())

    return alle

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    alle= todo(args)
    plt.plot(alle[:,0],alle[:,1])

    plt.plot([0,40],[alle[:,1].mean(),alle[:,1].mean()],'--',label='average transition temperature (K)')
    plt.legend()
    plt.xlabel('Neural Network nr')
    plt.ylabel('Phase Transition (K)')

    # see no correlation between transition temperateru and B or Ber or differences of those.
    #import pandas as pd
    #df = pd.DataFrame(data=alle,  columns=["T", "Bt","Btp","Bdt","Bdtp",'Bdiff','Bddiff'])

