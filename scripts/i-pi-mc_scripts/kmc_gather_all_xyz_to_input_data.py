#!/usr/bin/env python
from __future__ import print_function
import click
import numpy as np
import os,sys
from copy import deepcopy
from subprocess import call
import myutils as my
from ase.io import read,write

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
    e.g. \n
    kmc_gather_all_xyz_to_input_data.py `pwd`
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

    print('reading input.data ...')
    frames = read("input.data",":",format="runner")
    print('input.data has',len(frames),"frames.")
    print('removing duplicates ...')
    write("input.data",my.ase_get_unique_frames(frames),format="runner")

    remark0="# created input.data for n2p2/runner which contains only the unique structures;"
    remark1="# Now you need to do an nnp-scaling to get function.data on this input.data run for fps! (fps_considering_oldstruct.py);"
    remark2="# e.g. fps_considering_oldstruct.py -d1 /home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/n2p2_v1ag/ -d2 kmc_gaterh_all_xyz"
    my.mkdir("kmc_gaterh_all_xyz")
    my.create_READMEtxt(os.getcwd()+"/kmc_gaterh_all_xyz",add=[remark0,remark1,remark2])
    my.create_READMEtxt(os.getcwd(),add=[remark0,remark1,remark2])
    os.rename("input.data","kmc_gaterh_all_xyz/input.data")
    print()
    print(remark0)
    print(remark1)
    print(remark2)
    return filename_out


if __name__ == "__main__":
    gather_xyz()
