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
@click.option('--units','-u',type=click.Choice(['eV','meV_pa','hartree','hartree_pa']),default='hartree_pa',help='In which units should the output be given')
@click.option('--geopt/--no-geopt','-g',default=False,help='make a geometry optimization of the atoms.')
@click.option('--verbose','-v',count=True)


def get_energies(infile,format_in,pot,verbose,structures_idx,units,geopt):
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
    print()
    print('pot                          :',pot)
    print()
    print('verbose                      :',verbose)
    print('units                        :',units)
    print('geopt                        :',geopt)
    print()

    ene_DFT  = np.empty(len(structures_to_calc));ene_DFT[:]  = np.nan
    ene_pot  = np.empty(len(structures_to_calc));ene_pot[:]  = np.nan
    ene_diff = np.empty(len(structures_to_calc));ene_diff[:]  = np.nan
    ene_diff_abs = np.empty(len(structures_to_calc));ene_diff_abs[:]  = np.nan
    ene_std  = np.empty(len(structures_to_calc));ene_std[:]  = np.nan
    ene_ste  = np.empty(len(structures_to_calc));ene_ste[:]  = np.nan
    ene_mean = np.empty(len(structures_to_calc));ene_mean[:] = np.nan

    #sys.exit('get uuid of structure and save structure energy somewhere (cache)')
    #sys.exit('find out weather particular structure in test or trainset')
    #sys.exit('make geopt for structures, if different write to runner file')
    # make this parallel at some point

    ##############################################
    # get the potential
    ##############################################
    lmpcmd, atom_types = my.pot_to_ase_lmp_cmd(pot,verbose=verbose)
    lmpcmdgeopt, atom_types = my.pot_to_ase_lmp_cmd(pot,geopt=True,verbose=verbose)

    for idx,i in enumerate(structures_to_calc):
        ene_DFT[idx] = my.ase_enepot(atoms[i],units=units)
        ene_pot[idx] = my.ase_calculate_ene_from_pot(atoms[i],lmpcmd=lmpcmd,atom_types=atom_types,units=units, verbose=False)
        print('ENE',ene_pot[idx])
        ene_diff[idx] = ene_DFT[idx]-ene_pot[idx]
        ene_diff_abs[idx] = np.abs(ene_DFT[idx]-ene_pot[idx])

        if idx == 0:
            ene_std[idx] = 0.
            ene_ste[idx] = 0.
            ene_mean[idx] = ene_diff[idx]
        else:
            ene_std[idx] = np.std(ene_DFT[:idx+1]-ene_pot[:idx+1])
            ene_ste[idx] = ene_std[idx]/np.sqrt(idx)
            ene_mean[idx] = np.mean(np.abs(ene_DFT[:idx+1]-ene_pot[:idx+1]))

        printed = False

        if idx in range(0,len(structures_to_calc),50):
            print("%5.0f %5.0f / %6.0f %16.7f %s" % (i,idx,len(structures_to_calc),ene_diff_abs[idx],units))
            printed = True

        if verbose > 1 and printed == False:
            print("%5.0f %5.0f / %6.0f %16.7f %s" % (i,idx,len(structures_to_calc),ene_diff_abs[idx],units))

        if geopt:
            print('goeop',lmpcmdgeopt)
            print(atoms[i].get_positions())
            print('-----------------------------------')
            ene_pot_geopt = my.ase_calculate_ene_from_pot(atoms[i],lmpcmd=lmpcmdgeopt,atom_types=atom_types,units=units, verbose=False)
            print('A------------------AAAAA-----------')
            print(idx,'e diff ',ene_pot[idx],ene_pot_geopt,units)
            print(atoms[i].get_positions())
            print('-----------------------------------')


    my.create_READMEtxt(os.getcwd())
    np.savetxt("ene_DFT.npy",ene_DFT,header=units)
    np.savetxt("ene_pot.npy",ene_pot)
    np.savetxt("ene_std.npy",ene_std)
    np.savetxt("ene_ste.npy",ene_ste)
    np.savetxt("ene_mean.npy",ene_mean)
    ene_all = np.transpose([range(len(ene_DFT)),ene_DFT,ene_pot,ene_diff_abs,ene_mean,ene_std,ene_ste])
    print('en',len(ene_all))
    np.savetxt("ene_all.npy",ene_all,header=units+"\n"+"DFT\t\t"+pot+"\t|diff|\t\t<|diff|>",fmt=' '.join(['%i'] + ['%.10e']*6))
    return

if __name__ == "__main__":
    get_energies()
