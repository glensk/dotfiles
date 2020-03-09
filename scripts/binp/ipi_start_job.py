#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import os,sys,argparse
import myutils as my
from subprocess import check_output,call
from datetime import datetime as datetime   # datetime.datetime.now()
import time
import fah as fah_

def help(p = None):
    string = '''
    myutils.py -ex abc  # executs function abc
    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-ef','--execute_function', required=False, type=str,default='', help="function to run")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p


def ipi_start_job(inputfile="input.xml",sleep=12):
    ''' sleep of 6s is not enough for ~700 atoms '''
    if not os.path.isfile('in.lmp'):
        sys.exit('in.lmp does not exist! Exit')
    if not os.path.isfile(inputfile):
        sys.exit(inputfile+' does not exist! Exit')
    host = os.environ['myhost']
    print('host',host)
    executable = os.environ["HOME"]+"/Dropbox/Albert/scripts/dotfiles/scripts/executables/lmp_"+host
    print('lammps executable',executable)
    if not os.path.isfile(executable):
        sys.exit(executable+' does not exist! Exit')

    python = "python"
    if host == 'mac': python = "/usr/local/bin/python"
    print(python+" $HOME/sources/ipi/bin/i-pi "+inputfile+" &")
    adress_str = my.grep(inputfile,"md_ff_*")[0].split()[1]
    print('ad',adress_str)
    filecheck = '/tmp/ipi_'+adress_str #+'='
    print('fc',filecheck)
    if os.path.exists(filecheck):  # here os.path.isfile did not work
        #print('file exists ')
        os.remove(filecheck)
    else:
        print('file does not exist and this is good;')
    call([python+" $HOME/sources/ipi/bin/i-pi "+inputfile+" &"],shell=True)
    print('now sleep for ',sleep,"(sec)")
    time.sleep(sleep)
    time_now = datetime.now()
    print('sleep done; now sart lammps',time_now)
    call([executable+" < in.lmp"],shell=True)
    time_now = datetime.now()
    print('job done',time_now)
    return

def ipi_start_job_fah_and_evaluate():
    ipi_start_job()

    # evaluate the job
    fah_.get_dudl_from_ipi_job()
    # eveything else is fast
    return

if __name__ == "__main__":
    p = help()
    args = p.parse_args()
    if args.verbose:
        print_args(args)
    if args.execute_function:
        function = eval(args.execute_function)
        function()
        hier = os.path.abspath(os.getcwd())
        readmepath = my.create_READMEtxt(hier)

