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
    p.add_argument('-f','--inputfolder', required=True, type=str,default=False, help="path to the folder to restart from")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def todo(args):
    if not os.path.isdir(args.inputfolder):
        sys.exit("inputfile "+args.inputfolder+" does not exist")
    folder = os.path.abspath(args.inputfolder)
    print('folder', folder)

    my.cp(folder+"/input.nn",".")
    print('changing input.nn to use old weights ...')
    my.sed("input.nn",'#use_old_weights_short','use_old_weights_short')
    my.cp(folder+"/input.data",".")
    my.cp(folder+"/scaling.data",".")
    print('cp function.data ... may take some time')
    my.cp(folder+"/function.data",".")
    my.cp(folder+"/testing.data",".")
    my.cp(folder+"/testforces.data",".")
    my.cp(folder+"/trainforces.data",".")
    my.cp(folder+"/optweights.014.out","weights.014.data")
    my.cp(folder+"/optweights.013.out","weights.013.data")
    my.cp(folder+"/optweights.012.out","weights.012.data")
    my.cp(folder+"/trainstruct.data",".")
    my.cp(folder+"/teststruct.data",".")
    print('now run: ../RuNNer.serial.cosmopc.natascha.x| tee -a logfile_mode2')

    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    todo(args)


