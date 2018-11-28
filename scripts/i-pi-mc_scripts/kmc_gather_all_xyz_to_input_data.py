#!/usr/bin/env python
import click
import os,sys
from subprocess import call
import myutils as my

# get help also with -h
CONTEXT_SETTINGS = my.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)

@click.argument('folder')
@click.option('-ff','--filename_find',required=False,default="simulation.pos_0.xyz",help="Filename or extension of the searched file")
@click.option('-fo','--filename_out',required=False,default="input.data",help="Filename of the output file")

def gather_xyz(folder,filename_find,filename_out):
    '''
    This scipt looks for all files ($filename_find; default: simulation.pos_0.xyz) in
    $folder and converts content to RuNNer readable fileformat $filename_out
    (default: input.data).
    '''
    print('folder:',folder)
    print()
    print('files:')
    files = my.findfiles(folder,extension_or_filename=filename_find)

    if os.path.isfile(filename_out):
        sys.exit('The file '+filename_out+' does already exist! Exit.')

    f = open(filename_out,"wb")
    for i in files:
        print(i)
        call(["xyz2runner.sh",i],stdout=f)

    my.create_READMEtxt(os.getcwd(),add="# created input.data for n2p2/runner")
    return filename_out


if __name__ == "__main__":
    gather_xyz()
