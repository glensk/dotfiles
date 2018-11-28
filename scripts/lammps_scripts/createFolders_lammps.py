#!/usr/bin/env python
from __future__ import print_function
import os,sys
from shutil import copyfile
import click
import myutils as my

CONTEXT_SETTINGS = my.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)

@click.argument('foldername')
@click.option('--calc_type','-t', type=click.Choice(['static','geopt','ti']),default='static')


def createFolder_lammps(foldername,calc_type):
    ''' calc_type any of [static|geopt|ti] '''
    if os.path.isdir('foldername'):
        sys.exit(foldername+' does already exist;Exit.')
    my.mkdir(foldername)
    scripts = my.scripts()
    file_inlmp = scripts + "/i-pi-mc_scripts/in.lmp"
    copyfile(file_inlmp, foldername+"/in.lmp")
    return


if __name__ == "__main__":
    createFolder_lammps()


