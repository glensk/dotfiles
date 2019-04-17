#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    #p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def todo(args):
    if not os.path.isfile('scaling.data'):
        sys.exit('scaling.data does not exist')
    if not os.path.isfile('input.data'):
        sys.exit('input.data does not exist')
    if not os.path.isfile('input.nn'):
        sys.exit('input.nn does not exist')

    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    todo(args)


