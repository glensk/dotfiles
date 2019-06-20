#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse,subprocess,time
import myutils as my
from ase.io import read as ase_read
from ase.io import write as ase_write

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def convert_datafile_to_extxyz(args):
    #if not os.path.isfile(args.inputfile):
    #    sys.exit("inputfile "+args.inputfile+" does not exist")
    fn = args.inputfile
    bn = fn.replace(".data","")

    if False:
        print('read input_23.extxyz (104MB)')
        start_time = time.time()
        ak = ase_read('input_23.extxyz',":") #,format="extxyz")
        t = time.time() - start_time
        print('done in',t,'sec')
        print()

        print('read input_01.extxyz (11MB)')
        start_time = time.time()
        ak = ase_read('input_01.extxyz',":") #,format="extxyz")
        t = time.time() - start_time
        print('done in',t,'sec')
        print()

    if False:
        fnn = bn+".extxyz.gz"
        print('read',fnn)
        #from ase.io import filetype as ase_ft
        #print(ase_ft(fnn))
        #sys.exit()
        start_time = time.time()
        an = ase_read(fnn,":",parallel=False,format="extxyz")
        t = time.time() - start_time
        print('done in',t,'sec')
        print()
        sys.exit()


    # make sure that to be created file does not already exist
    if os.path.isfile(bn+".extxyz"):
        sys.exit("target file "+bn+".extxyz does already exist!")
    if os.path.isfile(bn+".extxyz.bz2"):
        sys.exit("target file "+bn+".extxyz.bz2 does already exist!")

    # make sure that it ends with .data
    if fn[-5:] != ".data":
        sys.exit("file needs to end with .data")

    print('filename fn',fn)
    print('basename bn',bn)
    bns = bn.split("_")
    print('basenamesplit bns',bns)
    if len(bns) != 2:
        sys.exit("not lenght of 2")
    isint = my.is_int(bns[1])
    if isint != True:
        sys.exit("bns[1] is not an integer")
    print()
    print('reading',fn,' ....')
    a = ase_read(args.inputfile,":",format="runner")

    print('writing',bn+'.extxyz',' ....')
    ase_write(bn+".extxyz",a,format='extxyz')
    size = os.path.getsize(bn+".extxyz")/(1024*1024)
    print("size:",size,"MB")

    # this is necessary if the extxyz is larger than 100MB
    if size > 99:
        print('writing',bn+'.extxyz.gz',' ....')
        ase_write(bn+".extxyz.gz",a,format='extxyz')
        os.remove(bn+".extxyz")

    print('tarring',fn,' ...')
    subprocess.call("tar.bzip2_file_or_folder_and_remove_originals.sh "+fn,shell=True)
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    convert_datafile_to_extxyz(args)


