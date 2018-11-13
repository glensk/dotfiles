#!/usr/bin/env python
import os,sys,random,massedit,ase
import socket
import datetime
from ase.lattice.cubic import FaceCenteredCubic
from shutil import copyfile
import numpy as np
from subprocess import call
from subprocess import check_output
import click

#@ # show default values in click
#@ orig_init = click.core.Option.__init__
#@ def new_init(self, *args, **kwargs):
#@     orig_init(self, *args, **kwargs)
#@     self.show_default = True
#@ click.core.Option.__init__ = new_init
#@
#@
#@ # get help also with -h
#@ CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
#@ @click.command(context_settings=CONTEXT_SETTINGS)


#from kmc_createjob import commands as group1
#
#@click.group()
#def entry_point():
#    pass
#
#entry_point.add_command(group1.command_group)

#import convert_fileformats
#import kmc_createjob as ka
#
#
#@click.group()
#def cli():
#    pass
#
#@click.command()
#def install():
#    ''' install ipi '''
#    pass
#
#@click.command()
#def kmc():
#    ''' kinetic-monte-carlo '''
#    pass
#
#@click.group()
#def entry_point():
#    pass
#
#cli.add_command(install)
#cli.add_command(kmc)
#cli.add_command(ka.main())
#
#if __name__ == '__main__':
#    cli()

from command_cloudflare import cloudflare
from command_uptimerobot import uptimerobot

cli = click.CommandCollection(sources=[cloudflare, uptimerobot])

if __name__ == '__main__':
    cli()

# input variables (get those with click)
# setting KMC
#@click.option('-ncell',required=True, prompt=True, type=int,
#        help="supercell size of primitive cell")
#@click.option('-nsi'  ,required=True, prompt=True, type=int, help="number of Si atoms")
#@click.option('-nmg'  ,required=True, prompt=True, type=int, help="number of Mg atoms")
#@click.option('-nvac' ,required=True, prompt=True, type=int, help="number of vacancies")
#@click.option('-a0'   ,default = 4.057, type=float, help="fcc lattice constant for Al")
#@click.option('-temp' ,default = 300, type=int, help="KMC temperature")
#@click.option('-nseeds',type=int, default=3, help="number of different seeds")
#@click.option('-seednumber',default=False,multiple=True, type=int,help="define seed number manually (can be defined multiple times)")
#@click.option('-nsteps',type=int, default=200000, help="number of KMC steps to make")
#@click.option('-runnercutoff',type=float, default=10., help="runner cutoff distance ~10Angstrom")
#
#
## environment variables
#@click.option('-scripts', envvar='scripts',help='environment variable $scripts (can alse be set here)')
#@click.option('-nn_pot',type=str, default="v2dg", help="foldername in $scripts containing the neural network potential (can be separately set here for different path)")
#@click.option('-i_pi_mc', envvar='i_pi_mc',help='path to i-pi-mc (or environment variable $i_pi_mc)')
#@click.option('-lmp_exec', envvar='lmp_exec',help='path to lammps executable (or environment variable $lmp_exec)')
#@click.option('-submit/-no-submit', default=False)
#@click.option('-submitdebug/-no-submitdebug', default=False)
