#!/usr/bin/env python

import os,sys
python_version = sys.version_info[0]
if python_version < 3:
    sys.exit('Your python environment uses a python < 3; Exit;')

import click
import glob
from shutil import copyfile
import numpy as np

#from subprocess import call,Popen,check_output,PIPE,run
from subprocess import call,run,check_output


# from scripts folder
import kmc_createjob

# show default values in click
orig_init = click.core.Option.__init__
def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True
click.core.Option.__init__ = new_init


# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)

def make_fps():
    '''
    This scipt ....
    '''
    ######################################
    # check weather input.data exists
    ######################################
    kmc_createjob.check_isfile(\
            ["function.data","logfile_mode1.1","input.data"],
            ["function.data","logfile_mode1.1","input.data"])

    ######################################
    # write README
    ######################################
    kmc_createjob.create_READMEtxt(os.getcwd())

    return


if __name__ == "__main__":
    make_fps()