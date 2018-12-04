#!/usr/bin/env python

from __future__ import print_function
import sys,os
import click
import numpy as np
import myutils as my
from ase.io import read
from ast import literal_eval

CONTEXT_SETTINGS = my.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)

@click.option('--infile','-i',required=True,type=str,help='input files containing structures that will be imported by ase')
@click.option('--format_in','-fi',type=str,default='runner',help='ase format for reading files')
@click.option('--pot','-p',type=click.Choice(my.mypot()),required=True,default='n2p2_v1ag')
@click.option('--structures_idx','-idx',default=':',help='which structures to calculate, use ":" for all structues (default), ":3" for structures [0,1,2] etc. (python notation)')
@click.option('--verbose','-v',count=True)


def get_energies(infile,format_in,pot,verbose,structures_idx):
    ''' this is a script which computes for a given set of structures the energies
    for a given potential.
    '''
    scripts = my.scripts()
    my.check_isfile_or_isfiles([infile],verbose=False)
    atoms = read(infile,index=":",format=format_in)
    print('number of structures in total:',len(atoms))
    print('structures_idx               :',structures_idx)
    structures_to_calc = my.string_to_index_an_array(range(len(atoms)),structures_idx)
    print('structures_to_calc           :',structures_to_calc)

    #sys.exit('get uuid of structure and save structure energy somewhere (cache)')
    #sys.exit('find out weather particular structure in test or trainset')
    #sys.exit('make geopt for structures, if different write to runner file')
    for i in structures_to_calc:
        #print('i',i)
        ene_mev_pa_in = my.ase_enepot_mev_pa(atoms[i])
        #print('r',i,mev_pa_in)
        ene,ene_mev_pa = my.ase_calculate_ene_from_pot(atoms[i],pot,verbose=False)
        if verbose:
            #print('structure',i,'/',len(atoms),'ene fromfile',ene_mev_pa_in,"meV/atom")
            #print('         ',i,'/',len(atoms),'ene frompot ',ene_mev_pa,"meV/atom")
            print(np.abs(ene_mev_pa_in-ene_mev_pa),'meV/atom (difference)')

    my.create_READMEtxt(os.getcwd())
    return

if __name__ == "__main__":
    get_energies()
