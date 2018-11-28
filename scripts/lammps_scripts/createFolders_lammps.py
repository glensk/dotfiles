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
@click.option('--pot','-p', type=click.Choice(['eam','runner','n2p2']),default='eam')

def createFolder_lammps(foldername,calc_type,pot):
    ''' calc_type any of [static|geopt|ti] '''
    if os.path.isdir('foldername'):
        sys.exit(foldername+' does already exist;Exit.')
    my.mkdir(foldername)
    scripts = my.scripts()
    file_inlmp = scripts + "/i-pi-mc_scripts/in.lmp"

    # make in.lmp
    get_inlmp(foldername,"../input.data.lmp.runner",pot)

    return

def show_eam_pot():
    scripts = my.scripts()
    scripts + '/potentials/lammps_eam/'
    return


def get_inlmp(folder,positions_file,pot):
    with open(folder+'/in.lmp', "w") as f:
        # stuff which is always added
        f.write("\n")
        f.write("clear\n")
        f.write('variable dump_file string "trj_lammps"\n')
        f.write("\n")
        f.write("units metal\n")
        f.write("boundary  p p p\n")
        f.write("atom_modify sort 0 0.0\n")
        f.write("\n")
        f.write("read_data "+str(positions_file)+"\n")
        f.write("\n")

        # potential
        f.write("## potential\n")
        if pot == 'eam':
            f.write("pair_style eam/alloy\n")
            f.write("pair_coeff * * Al-LEA.eam.alloy Al\n")
        f.write("\n")

        # run
        f.write("## run\n")
        f.write("dump dump_all all custom 50 trj_lammps id type x y z fx fy fz\n")
        f.write("thermo_style custom step cpu pe etotal vol\n")
        f.write("thermo 10\n")
        f.write("\n")
        f.write("min_style cg\n")
        f.write("#minimize 1.0e-9 1.0e-10 1000 1000\n")
        f.write("\n")
        f.write("run 0\n")
    return



if __name__ == "__main__":
    createFolder_lammps()
