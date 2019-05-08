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
    p.add_argument('-f','--folder', type=str,default=".", help="folder to evalute")
    p.add_argument('-nr', type=str,default=False, help="number of the potential inputfile e.g. 18 for weights.012.000018.out")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def n2p2_get_best_test_nr_from_learning_curve(folder):
    if not os.path.isdir(folder):
        sys.exit(folder+' does not exist! (7)')
    folder = os.path.abspath(folder)
    print('folder',folder)

    # in case this is a n2p2 job:
    if os.path.isfile(folder+"/learning-curve.out") and not os.path.isfile(folder+"/optweights.012.out"):
        a = np.loadtxt(folder+"/learning-curve.out")
        best_testset = np.argmin(a[:,2])
        #print('best_testset',best_testset)
    elif os.path.isfile(folder+"/learning-curve.out"):
        #a = np.loadtxt(folder+"/learning_curve_test.dat")
        a = np.loadtxt(folder+"/learning-curve.out")
        #print(a)
        #print()
        #print(a[:,1])
        best_testset = np.argmin(a[:,1])
    else:
        sys.exit(folder+"/learning-curve.out does not exist (8)")

    print('best testset :',best_testset)
    #sys.exit()
    return best_testset

def n2p2_make_potential_folder_from_nr(argsnr):
    nr_ = str(args.nr).zfill(6)
    if not os.path.isfile('scaling.data'):
        sys.exit('scaling.data does not exist (1)')
    #print('scaling.data : exists')
    if not os.path.isfile('input.data'):
        sys.exit('input.data does not exist (2)')
    #print('input.data   : exists')
    if not os.path.isfile('input.nn'):
        sys.exit('input.nn does not exist (3)')
    #print('input.nn     : exists')


    checkfor = [ '012', '013', '014' ]
    def get_weightsfile(i,nr_):
        weights1 = 'weights.'+i+'.'+nr_+'.out'
        weights2 = '_weights/weights.'+i+'.'+nr_+'.out'
        weights3 = 'optweights.'+i+'.out'
        if os.path.isfile(weights1):
            weights = weights1
            typ = 'n2p2'
        elif os.path.isfile(weights2):
            weights = weights1
            typ = 'n2p2'
        elif os.path.isfile(weights3):
            weights = weights3
            typ = 'runner'
        else:
            sys.exit(weights1+" does not exist")
            sys.exit(weights2+" does not exist")
            sys.exit(weights3+" does not exist")
            sys.exit("weights files not found (4)")
        print(weights,'exist')
        return weights,typ

    for i in checkfor:
        get_weightsfile(i,nr_)

    folder = "potential_"+str(args.nr)+'/'
    folder = "potential/"
    if os.path.isdir(folder):
        sys.exit(folder+' does already exist!')

    my.mkdir(folder)
    print('cp scaling.dat')
    my.cp('scaling.data',folder+'/scaling.data')
    print('cp input.dat')
    my.cp('input.data',folder+'/input.data')
    print('cp input.nn')
    my.cp('input.nn',folder+'/input.nn')
    if os.path.isfile('learning-curve.out'):
        print('cp learning-curve.out')
        my.cp('learning-curve.out',folder+'/learning-curve.out')

    for i in checkfor:
        weights,typ = get_weightsfile(i,nr_)
        print('cp',weights)
        if typ == 'n2p2':
            my.cp(weights,folder+'/weights.'+i+'.data')
            my.cp(weights,folder+'/weights.'+i+'.'+nr_+'.out')
        elif typ == 'runner':
            my.cp(weights,folder+'/weights.'+i+'.data')
            my.cp(weights,folder+'/optweights.'+i+'.out')
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    if args.nr == False:
        args.nr = n2p2_get_best_test_nr_from_learning_curve(args.folder)
    n2p2_make_potential_folder_from_nr(argsnr=args.nr)


