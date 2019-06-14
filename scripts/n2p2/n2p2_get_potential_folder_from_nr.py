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
    p.add_argument('-l','--get_last_epoch', help='get potential from last epoch', action='count', default=False)
    p.add_argument('-b','--get_best_epoch', help='get potential with lowest testerror', action='count', default=False)
    p.add_argument('-tmp','--tmp', help='make potential_tmp instead of potential', action='count', default=False)
    p.add_argument('--best_testsete',default=False,help=argparse.SUPPRESS, type=int)
    return p

def n2p2_get_best_test_nr_from_learning_curve(args):
    folder = args.folder

    if not os.path.isdir(folder):
        sys.exit(folder+' does not exist! (7)')
    folder = os.path.abspath(folder)
    print('folder',folder)
    print('get_last_epoch',args.get_last_epoch)
    n2p2_or_runner = my.inputnn_runner_or_n2p2(folder+'/input.nn')
    learning_curve_file = my.n2p2_runner_get_learning_curve_filename(folder+"/input.nn")

    if not os.path.isfile(learning_curve_file):
        sys.exit(learning_curve_file+" does not exist!")

    print('learning_curve_file',learning_curve_file)

    learning_curve = lc = my.n2p2_runner_get_learning_curve(folder+'/input.nn')
    #print(lc)
    length = len(lc)-1
    number_of_points_to_save = 20
    geomspace = np.geomspace(1,length,num=number_of_points_to_save,endpoint=True,dtype=int)  # in geomspace this does not make sure the endpoint is really there
    print('g1',geomspace)
    geomspace = np.append(geomspace, length)
    print('g2',geomspace)
    print('folder',folder)
    best_testset = np.argmin(lc[:,2])
    args.get_best_epoch = args.best_testsete = best_testset

    # include best teststet
    geomspace = np.append(geomspace, best_testset)
    print('g3',geomspace)
    geomspace = np.unique(geomspace)
    print('g4',geomspace)
    #sys.exit()

    if args.get_last_epoch:
        last_epoch = len(learning_curve) - 1
        #print('lenxx',last_epoch)
        #print(lc[last_epoch])
        print('last epoch:',last_epoch)
        return [last_epoch]
    #print('learning_curve')
    #print(learning_curve)
    best_testset = np.argmin(lc[:,2])

    #print('best testset :',best_testset)
    #sys.exit()
    #return [best_testset]
    return np.sort(geomspace)

def n2p2_make_potential_folder_from_nr(argsnr):
    if not os.path.isfile('scaling.data'):
        sys.exit('scaling.data does not exist (1)')
    #print('scaling.data : exists')
    if not os.path.isfile('input.data'):
        sys.exit('input.data does not exist (2)')
    #print('input.data   : exists')
    if not os.path.isfile('input.nn'):
        sys.exit('input.nn does not exist (3)')
    #print('input.nn     : exists')
    typ = my.inputnn_runner_or_n2p2('input.nn')


    def get_weightsfile(i,epoch,verbose=False):
        weights1 = 'weights.'+i+'.'+epoch+'.out'
        weights2 = '_weights/weights.'+i+'.'+epoch+'.out'
        weights3 = 'optweights.'+i+'.out'
        if os.path.isfile(weights1):
            weights = weights1
        elif os.path.isfile(weights2):
            weights = weights2
        elif os.path.isfile(weights3):
            weights = weights3
        else:
            sys.exit(weights1+" does not exist")
            sys.exit(weights2+" does not exist")
            sys.exit(weights3+" does not exist")
            sys.exit("weights files not found (4)")
        #print(weights,'exist')
        return weights

    checkfor = [ '012', '013', '014' ]
    for epoch in args.nr:
        epoch000 = str(epoch).zfill(6)
        for i in checkfor:
            weights = get_weightsfile(i,epoch000)


    folder = "potential_"+str(args.nr)+'/'
    folder = "potential/"
    if args.tmp:
        folder = "potential_tmp/"
    if args.get_last_epoch:
        folder = "potential_last/"
    #if os.path.isdir(folder):
    #    sys.exit(folder+' does already exist!')

    my.mkdir(folder)
    print('cp scaling.dat')
    my.cp('scaling.data',folder+'/scaling.data')
    print('cp input.dat')
    my.cp('input.data',folder+'/input.data')
    print('cp input.nn')
    my.cp('input.nn',folder+'/input.nn')
    if os.path.isfile('submit_n2p2_train.sh'):
        print('cp submitfile')
        my.cp('submit_n2p2_train.sh',folder+'/submit_n2p2_train.sh')
    if typ == "n2p2" and os.path.isfile('learning-curve.out'):
        print('cp learning-curve.out')
        my.cp('learning-curve.out',folder+'/learning-curve.out')
    if typ == "runner" and os.path.isfile('learning-curve-runner.out'):
        print('cp learning-curve-runner.out')
        my.cp('learning-curve.out',folder+'/learning-curve-runner.out')
    if os.path.isfile('logfile_mode2'):
        print('cp logfile_mode2')
        my.cp('logfile_mode2',folder+'/log.fit')
    if typ == 'runner' and os.path.isfile('log.fit'):
        print('cp log.fit')
        my.cp('log.fit',folder+'/log.fit')

    np.savetxt(folder+'/best_testsete',np.array([int(args.best_testsete)]),fmt='%i')

    for epoch in args.nr:
        epoch000 = str(epoch).zfill(6)
        for i in checkfor:
            weights = get_weightsfile(i,epoch000)
            #print('cp',weights)
            if typ == 'n2p2':
                my.cp(weights,folder+'/weights.'+i+'.'+epoch000+'.out')
                if epoch == args.get_best_epoch:
                    print('epoch',epoch)
                    print('args.get_best_epoch',args.get_best_epoch)
                    print('best',i,"(i == checkfor)")
                    src = 'weights.'+i+'.'+epoch000+'.out'
                    dest = 'weights.'+i+'.data'
                    print('src',src)
                    print('dest',dest)
                    pwd = os.getcwd()
                    os.chdir(folder)
                    if os.path.islink(dest):
                        os.unlink(dest)
                    os.symlink(src, dest)
                    os.chdir(pwd)
            elif typ == 'runner':
                my.cp(weights,folder+'/weights.'+i+'.data')
                my.cp(weights,folder+'/optweights.'+i+'.out')
    os.chdir(folder)
    my.create_READMEtxt(directory=os.getcwd(),add="# pwd: "+os.getcwd())
    return os.getcwd()

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    if args.nr == False or args.get_best_epoch == False:
        args.nr = n2p2_get_best_test_nr_from_learning_curve(args)
    print('args.nr',args.nr)
    print('args.get_best_epoch',args.get_best_epoch)
    folder = n2p2_make_potential_folder_from_nr(argsnr=args.nr)
    print('folder',folder)
    with my.cd(folder):
        print('getEnergies_byLammps.py -p . -ea # to get all the c44')
        import subprocess
        subprocess.call("getEnergies_byLammps.py -p . -ea",shell=True)
        print()
        print()
        print()
        print('getEnergies_byLammps.py -p . --testkmc_b # to get all the c44')
        subprocess.call("getEnergies_byLammps.py -p . --testkmc_b",shell=True)
        print('getEnergies_byLammps.py -p . --testkmc_l # to get all the c44')
        subprocess.call("getEnergies_byLammps.py -p . --testkmc_l",shell=True)


