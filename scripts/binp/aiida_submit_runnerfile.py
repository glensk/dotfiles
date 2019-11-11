#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse
from ase.io import read
import aiida_load_runner_dataset_as_aiida_group as aiida_read
#import aiida_launch_workflow_alalloy as aiida_lauch_job
from aiida_utils import create_READMEtxt
from subprocess import call

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-cm','--calc_method', choices=["scf", "relax", "vc-relax"], required=True, help='The type of calculation to perform')
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def aiida_submit_job(args):
    if not os.path.isfile(args.inputfile):
        sys.exit("inputfile "+args.inputfile+" does not exist")

    # read in frame just to varify that everythin is ok.
    frames = read(args.inputfile,":",format="runner")
    for i in frames:
        print(i)
    groupname = os.path.basename(args.inputfile)
    print('groupname:',groupname)


    # aiida_load_runner_dataset_as_aiida_group
    dataset_path = args.inputfile
    group_name   = os.path.basename(args.inputfile)

    print('aaa dataset_path',dataset_path)
    print('aaa group_name',group_name)
    print()

    # callback is necessary since launch is a decorated click funtion and could not be called otherwise
    aiida_read.launch.callback(
            dataset_path=dataset_path,
            group_name=group_name,
            group_description="",
            parse_comments_path=False,
            parse_comments_structure=False,
            supress_readme = False
            )

    # launch the job
    print("launching job ...")
    # callback is necessary since launch is a decorated click funtion and could not be called otherwise
    bpn = "9b370584-3f56-471c-a724-dbaadf022ec5"
    print('submitting to --code_node 114134') #113998') #114027')
    command = ["aiida_launch_workflow_alalloy.py",
            "--code_node","114134", #113998", #114027", #"1" .. "113998",  # this would be the new daint "code_node"
            "--structure_group_name",group_name,
            "--workchain_group_name",group_name+"_calc",
            "--base_parameter_node",bpn,
            "--pseudo_familyname", "SSSP_v1.1_eff",
            "--kptper_recipang", "80",
            "--nume2bnd_ratio", "0.75",
            "--calc_method", args.calc_method,
            ]
    call(command)

    # create README
    create_READMEtxt(os.getcwd())
    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    aiida_submit_job(args)
