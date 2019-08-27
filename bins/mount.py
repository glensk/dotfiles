#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse
from subprocess import call

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    #p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def mount(args):
    call(["sshfs glensk@fidis.epfl.ch:/home/glensk/ /home/glensk"],shell=True)
    print('check /home/glensk/')
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    mount(args)


