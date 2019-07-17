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
    #p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    #ddtodo(args)
    my.n2p2_get_scaling_and_function_data()
    my.create_READMEtxt(os.getcwd())



