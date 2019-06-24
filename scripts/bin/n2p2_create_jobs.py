#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse,subprocess
import myutils as my

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-s','--seeds',default=[12345, 8765, 987654, 654911])
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def todo(args):
    hier=os.getcwd()
    if not os.path.isfile(args.inputfile):
        sys.exit("inputfile "+args.inputfile+" does not exist")
    if not os.path.isfile("input.nn"):
        sys.exit("inputfile input.nn does not exist")
    for s in args.seeds:
        #for pl in ["p","t"]:
        for pl in ["p"]:
            for c in [21]:
                os.chdir(hier)
                folder=str(s)+"_"+pl+pl+"l_"+str(c)+"cores_repeated"
                #print('s',s,pl,pl,"l",'cores',c)
                print(folder)
                if os.path.isdir(folder):
                    print(folder+' exists!')
                    continue
                my.mkdir(folder)
                my.cp("input.nn",folder+"/input.nn")
                my.cp(args.inputfile,folder+"/input.data")
                my.sed(folder+"/input.nn",'random_seed.*','random_seed '+str(s))
                my.sed(folder+"/input.nn",'global_activation_short.*','global_activation_short '+str(pl)+" "+str(pl)+" l")
                with my.cd(folder):
                    folder=os.getcwd()
                    my.mkdir(folder)
                    print('strat scaling.data')
                    my.n2p2_get_scaling_and_function_data(cores=c,submitdebug=True,submit_to_que=False,submitmin=2,interactive=True)
                    print('pwd',os.getcwd())
                    os.chdir("get_scaling")
                    subprocess.call(["./submit_n2p2_scaling.sh"],shell=True)
                    print('pwd',os.getcwd())
                    os.chdir(folder) # necessary since n2p2_get_scaling_and_function_data changed the folder, now hopefully not anymore
                    my.n2p2_make_training(cores=c,debugque=False)
                sys.exit()

    my.create_READMEtxt(os.getcwd())
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()

    todo(args)


