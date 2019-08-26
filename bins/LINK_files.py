#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse
from myhost import check_for_known_hosts
import subprocess

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

scratch={}
scratch["cosmopc"] = "/local/scratch"
scratch["datin"]   = "/scratch/snx3000/aglensk"

def LINK(args):
    myhost = check_for_known_hosts()
    print('myhost',myhost)
    #if myhost == 'cosmopc'
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    LINK(args)
    os.environ["DEBUSSY1"] = "1"
    subprocess.call(['export DEBUSSY="3"'],shell=True)


