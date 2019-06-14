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
    p.add_argument('-c', '--cores' , required=False, help='cores to use for the job', type=int, default=28)
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
<<<<<<< HEAD
    p.add_argument('-n','--normalque', help='use normal que instead of debug que', action='store_true', default=False)
=======
    p.add_argument('-d','--debug', help='use debug que', action='count', default=True)
>>>>>>> 1d634a66223c61d8fe8a4fbee536b2246c2fe485
    return p

def todo(args):
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
<<<<<<< HEAD
    my.print_args(args)
    print('debug?',not args.normalque)
    my.n2p2_get_scaling_and_function_data(cores=args.cores,debug=(not args.normalque))
=======
    my.n2p2_get_scaling_and_function_data(cores=args.cores,debug=args.debug)
>>>>>>> 1d634a66223c61d8fe8a4fbee536b2246c2fe485
    my.create_READMEtxt(os.getcwd())


