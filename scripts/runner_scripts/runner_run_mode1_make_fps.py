#!/usr/bin/env python

import os,sys

import click
import glob
from shutil import copyfile
import numpy as np
from subprocess import check_output,run

#from subprocess import call,Popen,check_output,PIPE,run
from subprocess import call,run,check_output
import myutils as my


CONTEXT_SETTINGS = my.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)

def make_fps():
    '''
    - This scipt uses CurSel.py
    - CurSel.py will choose the most important symmetry functions

    '''
    ######################################
    # check weather input.data exists
    ######################################
    my.check_isfile_or_isfiles(\
            ["function.data","logfile_mode1.1","input.data"],
            ["function.data","logfile_mode1.1","input.data"])


    ######################################
    # run CurSel.py
    ######################################
    structures = check_output(["grep -c begin input.data"],shell=True).decode(sys.stdout.encoding).strip()
    print('len',structures)
    print("CurSel.py runs always about 500 seconds, irrespective of lenght 20-2500 landmarks, therefore, do it for all the files")
    call(["CurSel.py","-t","1e-3","--landmarks",str(structures),"function.data","logfile_mode1.1"]) # this is interactive, CurSel.py output is written to screen!


    ######################################
    # write README
    ######################################
    my.create_READMEtxt(os.getcwd())
    return


if __name__ == "__main__":
    make_fps()
