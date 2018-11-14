#!/usr/bin/env python
import click
import glob
import os,sys
from subprocess import call

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

@click.argument('folder')
@click.option('-ff','--filename_find',required=False,default="simulation.pos_0.xyz",help="Filename or extension of the searched file")
@click.option('-fo','--filename_out',required=False,default="input.data.all",help="Filename of the output file")

def gather_xyz(folder,filename_find,filename_out):
    '''
    This scipt looks for all files ($filename_find; default: simulation.pos_0.xyz) in
    $folder and converts content to RuNNer readable fileformat $filename_out
    (default: input.data.all).
    '''
    print('folder:',folder)
    print()
    print('files:')
    files = findfiles(folder,filename_find)

    if os.path.isfile(filename_out):
        sys.exit('The file '+filename_out+' does already exist! Exit.')

    f = open(filename_out,"wb")
    for i in files:
        print(i)
        call(["xyz2runner.sh",i],stdout=f)

    kmc_createjob.create_READMEtxt(os.getcwd())
    return filename_out

def findfiles(directory,extension_or_filename):
    listout=[]
    for filename in glob.iglob(directory+'**/*'+extension_or_filename, recursive=True):
        #print(filename)
        listout.append(filename)
    return listout

if __name__ == "__main__":
    gather_xyz()
