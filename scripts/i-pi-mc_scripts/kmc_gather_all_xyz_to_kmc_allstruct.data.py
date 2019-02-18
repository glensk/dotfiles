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
@click.option('-fo','--filename_out',required=False,default="kmc_allstruct.data",help="Filename of the output file")

def gather_xyz(folder,filename_find,filename_out):
    '''
    This scipt looks for all files ($filename_find; default: simulation.pos_0.xyz) in
    $folder and converts content to RuNNer readable fileformat $filename_out
    (default: kmc_allstruct.data).

    e.g. kmc_gather_all_xyz_to_kmc_allstruct.data.py `pwd`
    e.g. kmc_gather_all_xyz_to_kmc_allstruct.data.py .
    '''
    frames = read("kmc_allstruct.data",":",format="runner")
    print('readin')
    write("kmc_allstruct.data.new",frames,format="runner",setenergy_eV=0,setforces_ase_units="str_is_settozero")
    #write(filename_out+'tmp',i,format="runner",append=True)
    print('wrote')
    sys.exit()
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

    print('reading',filename_out,'...')
    frames = read(filename_out,":",format="runner")
    print(filename_out,'has',len(frames),"frames.")
    print('removing duplicates ...')
    framesuniq = my.ase_get_unique_frames(frames)
    print()
    print('removing duplicates done.')
    print('len to write',len(framesuniq))
    print('writing',filename_out)
    length = len(framesuniq)
    for idx,i in enumerate(framesuniq):
        print(i.get_positions())
    sys.exit()
    for idx,i in enumerate(framesuniq):
        my.progress(idx, length , status='')
        write(filename_out+'tmp',i,format="runner",append=True)

    print('written tmp')
    os.rename(filename_out+'tmp',filename_out)
    print('written',filename_out)


    remark0="# created "+filename_out+" for n2p2/runner which contains only the unique structures;"
    remark1="# Now you need to do an nnp-scaling to get function.data on this "+filename_out+" run for fps! (fps_considering_oldstruct.py);"
    remark2="# e.g. fps_considering_oldstruct.py -d1 /home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/n2p2_v1ag/ -d2 kmc_gaterh_all_xyz"
    my.mkdir("kmc_gaterh_all_xyz")
    my.create_READMEtxt(os.getcwd()+"/kmc_gaterh_all_xyz",add=[remark0,remark1,remark2])
    my.create_READMEtxt(os.getcwd(),add=[remark0,remark1,remark2])
    os.rename(filename_out,"kmc_gaterh_all_xyz/"+filename_out)
    print()
    print(remark0)
    print(remark1)
    print(remark2)
    #return filename_out  # not good when run from shell??
    return


if __name__ == "__main__":
    gather_xyz()
