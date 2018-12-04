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

    ##############################################
    # read in the structures
    ##############################################
    atoms = read(infile,index=":",format=format_in)

    print('number of structures in total:',len(atoms))
    print('structures_idx               :',structures_idx)
    structures_to_calc = my.string_to_index_an_array(range(len(atoms)),structures_idx)
    print('structures_to_calc           :',len(structures_to_calc))
    print('verbose                      :',verbose)

    ene_DFT  = np.empty(len(structures_to_calc));ene_DFT[:]  = np.nan
    ene_pot  = np.empty(len(structures_to_calc));ene_pot[:]  = np.nan
    ene_std  = np.empty(len(structures_to_calc));ene_std[:]  = np.nan
    ene_ste  = np.empty(len(structures_to_calc));ene_ste[:]  = np.nan
    ene_mean = np.empty(len(structures_to_calc));ene_mean[:] = np.nan

    #sys.exit('get uuid of structure and save structure energy somewhere (cache)')
    #sys.exit('find out weather particular structure in test or trainset')
    #sys.exit('make geopt for structures, if different write to runner file')
    # make this parallel at some point
    lmpcmd, atom_types = my.pot_to_ase_lmp_cmd(pot)
    for idx,i in enumerate(structures_to_calc):
        ene_mev_pa_in = my.ase_enepot_mev_pa(atoms[i])
        ene_DFT[idx] = ene_mev_pa_in
        ene,ene_mev_pa = my.ase_calculate_ene_from_pot(atoms[i],lmpcmd=lmpcmd,atom_types=atom_types,verbose=False)
        ene_pot[idx] = ene_mev_pa

        if idx == 0:
            ene_std[idx] = 0.
            ene_ste[idx] = 0.
            ene_mean[idx] = 0.
        else:
            ene_std[idx] = np.std(ene_DFT[:idx+1]-ene_pot[:idx+1])
            ene_ste[idx] = ene_std[idx]/np.sqrt(idx)
            ene_mean[idx] = np.mean(np.abs(ene_DFT[:idx+1]-ene_pot[:idx+1]))
        #print('---> idx',idx,'std',ene_std[idx])
        printed = False
        if idx in range(0,len(structures_to_calc),50):
            print(i,'/',len(structures_to_calc),np.abs(ene_mev_pa_in-ene_mev_pa),'meV/atom (difference)')
            printed = True

        if verbose > 1 and printed == False:
            #print('structure',i,'/',len(atoms),'ene fromfile',ene_mev_pa_in,"meV/atom")
            #print('         ',i,'/',len(atoms),'ene frompot ',ene_mev_pa,"meV/atom")
            print(i,'/',len(structures_to_calc),np.abs(ene_mev_pa_in-ene_mev_pa),'meV/atom (difference)')

    my.create_READMEtxt(os.getcwd())
    np.savetxt("ene_DFT.npy",ene_DFT)
    np.savetxt("ene_pot.npy",ene_pot)
    np.savetxt("ene_std.npy",ene_std)
    np.savetxt("ene_ste.npy",ene_ste)
    np.savetxt("ene_mean.npy",ene_mean)
    return

if __name__ == "__main__":
    get_energies()
