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
    p.add_argument('nr', type=str,default=False, help="number of the potential inputfile e.g. 18 for weights.012.000018.out")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def n2p2_make_potential_folder_from_nr(argsnr):
    nr_ = str(args.nr).zfill(6)
    if not os.path.isfile('scaling.data'):
        sys.exit('scaling.data does not exist')
    if not os.path.isfile('input.data'):
        sys.exit('input.data does not exist')
    if not os.path.isfile('input.nn'):
        sys.exit('input.nn does not exist')

    checkfor = [ '012', '013', '014' ]
    for i in checkfor:
        weights = 'weights.'+i+'.'+nr_+'.out'
        if not os.path.isfile(weights):
            sys.exit(weights+" does not exist")
        else:
            print(weights,'exist')

    folder = "potential_"+str(args.nr)+'/'
    if os.path.isdir(folder):
        sys.exit(folder+' does already exist!')

    my.mkdir(folder)
    print('cp scaling.dat')
    my.cp('scaling.data',folder+'/scaling.data')
    print('cp input.dat')
    my.cp('input.data',folder+'/input.data')
    print('cp input.nn')
    my.cp('input.nn',folder+'/input.nn')
    for i in checkfor:
        weights = 'weights.'+i+'.'+nr_+'.out'
        print('cp',weights)
        my.cp(weights,folder+'/weights.'+i+'.out')
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    n2p2_make_potential_folder_from_nr(argsnr=args.nr)


