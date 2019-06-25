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
    p.add_argument('-c', '--cores' , required=False, help='cores to use for the job', type=int, default=21)
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    p.add_argument('-d','--debugque', help='use the debug que (maxtime of 1 hour)', action='count', default=False)
    return p


if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    my.n2p2_make_training(cores=args.cores)
    my.create_READMEtxt(os.getcwd())


