#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse
import myutils as my

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-c', '--cores' , required=False, help='cores to use for the job', type=int, default=False)
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    p.add_argument('-d','--debugque', help='use the debug que (maxtime of 1 hour)', action='count', default=False)
    return p


if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    if args.cores == False:
        args.cores = 21
        myhost = os.environ['myhost']
        if myhost in [ "fidis", "helvetios" ]:
            args.cores = 21
        elif myhost == 'daint':
            args.cores = 21
    my.n2p2_make_training(cores=args.cores,days=0,hours=4,minutes=0,submit_to_que=True,submit_to_debug_que=False)
    my.create_READMEtxt(os.getcwd())


