#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('inputfile', help="name of the inputfile inputfile i.e. simulation.pos_0.xyz")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def todo(args):
    if not os.path.isfile(args.inputfile):
        sys.exit("inputfile "+args.inputfile+" does not exist")
    # read first two lines of inputfile
    from itertools import islice
    with open(args.inputfile) as myfile:
        head = list(islice(myfile, 2))
    print(head)
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    todo(args)


