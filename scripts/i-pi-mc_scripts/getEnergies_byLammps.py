#!/usr/bin/env python

from __future__ import print_function
import sys,os,copy
import click
import numpy as np
import myutils as my
from myutils import ase_calculate_ene #as ace
from ase.io import read as ase_read
from ase.io import write as ase_write

CONTEXT_SETTINGS = my.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)

@click.option('--infile','-i',required=True,type=str,help='input files containing structures that will be imported by ase')
@click.option('--format_in','-fi',type=str,default='runner',help='ase format for reading files')
@click.option('--pot','-p',type=click.Choice(my.pot_all()),required=True,default='n2p2_v1ag')
@click.option('--structures_idx','-idx',default=':',help='which structures to calculate, use ":" for all structues (default), ":3" for structures [0,1,2] etc. (python notation)')
@click.option('--units','-u',type=click.Choice(['eV','meV_pa','hartree','hartree_pa']),default='hartree_pa',help='In which units should the output be given')
@click.option('--geopt/--no-geopt','-g',default=False,help='make a geometry optimization of the atoms.')
@click.option('--ase/--no-ase','-a',default=True,help='Do the calculations by the ase interface to lammps.')
@click.option('--lmp/--no-lmp','-l',default=False,help='Do the calculations externally by lammps and not through ase interface.')
@click.option('--ipi/--no-ipi','-ipi',default=False,help='Do the calculations externally by ipi-lammps and not through ase interface.')
@click.option('--test/--no-test','-t',default=False,help='Assess formation energies of particular test structures.')
@click.option('--write_runner/--no-write_runner','-wr',required=False,default=False,help='default: runner.out')
@click.option('--verbose','-v',count=True)


def get_energies(infile,format_in,pot,verbose,structures_idx,units,geopt,test,ase,lmp,ipi,write_runner):
    ''' this is a script which computes for a given set of structures the energies
    for a given potential.
    '''

    #### get the potential
    ace = ase_calculate_ene(pot,units=units,geopt=geopt,verbose=verbose)
    units = ace.units

    ### when want to assess some formation energies
    if test:
        test_sisivac(ace)
        sys.exit('test done! Exit')

    ### read in the structures
    my.check_isfile_or_isfiles([infile],verbose=verbose)
    atoms = ase_read(infile,index=":",format=format_in)

    ### print stuff to screen
    print('number of structures in total:',len(atoms))
    print('structures_idx               :',structures_idx)
    structures_to_calc = my.string_to_index_an_array(range(len(atoms)),structures_idx)
    print('structures_to_calc           :',len(structures_to_calc))
    print()
    print('pot                          :',pot)
    print()
    print('units                        :',units)
    print('geopt                        :',geopt)
    print()
    print('verbose                      :',verbose)
    print('lmp                          :',lmp)
    print('ipi                          :',ipi)
    if write_runner:
        write_runner = 'runner.out'
    print('write_runner                 :',write_runner)
    print()

    ana_mg_conz      = np.empty(len(structures_to_calc));ana_mg_conz[:]  = np.nan
    ana_si_conz      = np.empty(len(structures_to_calc));ana_si_conz[:]  = np.nan
    ana_al_conz      = np.empty(len(structures_to_calc));ana_al_conz[:]  = np.nan
    ana_atoms        = np.empty(len(structures_to_calc));ana_atoms[:]  = np.nan

    ene_DFT      = np.empty(len(structures_to_calc));ene_DFT[:]  = np.nan
    ene_pot      = np.empty(len(structures_to_calc));ene_pot[:]  = np.nan
    ene_pot_ase  = np.empty(len(structures_to_calc));ene_pot_ase[:]  = np.nan
    ene_pot_lmp  = np.empty(len(structures_to_calc));ene_pot_lmp[:]  = np.nan
    ene_pot_ipi  = np.empty(len(structures_to_calc));ene_pot_ipi[:]  = np.nan

    ene_diff     = np.empty(len(structures_to_calc));ene_diff[:]  = np.nan
    ene_diff_abs = np.empty(len(structures_to_calc));ene_diff_abs[:]  = np.nan
    ene_std      = np.empty(len(structures_to_calc));ene_std[:]  = np.nan
    ene_ste      = np.empty(len(structures_to_calc));ene_ste[:]  = np.nan
    ene_mean     = np.empty(len(structures_to_calc));ene_mean[:] = np.nan
    ene_diff_lam_ase  = np.empty(len(structures_to_calc));ene_DFT[:]  = np.nan


    #sys.exit('get uuid of structure and save structure energy somewhere (cache)')
    #sys.exit('find out weather particular structure in test or trainset')
    # make this parallel at some point

    ##############################################
    # loop over structures
    ##############################################
    printevery = 50
    if len(structures_to_calc) < 50:
        printevery = 1
    be_very_verbose = 999


    print(',structures_to_calc[:3]:',structures_to_calc[:3])
    for idx,i in enumerate(structures_to_calc):
        d = my.ase_get_chemical_symbols_to_conz(atoms[i])
        #print(d)
        #print(d["Mg"])
        #print(d["Si"])
        ana_mg_conz[idx] = d["Mg"]
        ana_si_conz[idx] = d["Si"]
        ana_al_conz[idx] = d["Al"]
        ana_atoms[idx]   = atoms[i].get_number_of_atoms()

        ### ene from ipi
        if ipi == True:
            atoms_tmp = copy.deepcopy(atoms[i])
            ene_pot_ipi[idx] = my.ipi_ext_calc(atoms_tmp,ace)

        ### ene from DFT
        ene_DFT[idx] = my.ase_enepot(atoms[i],units=ace.units)
        if verbose > be_very_verbose:
            my.show_ase_atoms_content(atoms[i],showfirst=3,comment = "STAT2")
        if verbose > 0:
            print('ene_DFT[idx]     :',ene_DFT[idx],units)

        ### ene from ase (without writing lammps files)
        if ase == True:
            atoms_tmp = copy.deepcopy(atoms[i])  # for other instances, since atoms change when geoopt
            ace.pot_to_ase_lmp_cmd()
            ene_pot_ase[idx] = ace.ene(atoms_tmp)
            if verbose > 0:
                print('ene_pot_ase[idx] :',ene_pot_ase[idx],units)
            if lmp == False:
                ene_pot[idx] = copy.deepcopy(ene_pot_ase[idx])

        ### ene from lammps (by writing lammps files)
        if lmp == True:
            atoms_tmp = copy.deepcopy(atoms[i])  # for other instances, since atoms change when geoopt
            ene_pot_lmp[idx] = my.lammps_ext_calc(atoms_tmp,ace)
            if verbose:
                print("ENE DFT   :",i,ene_DFT[idx],ace.units)
                print("ENE ase   :",i,ene_pot_ase[idx],ace.units)
                print("ENE lammps:",i,ene_pot_lmp[idx],ace.units)
                print("--------------------------------------")
            ene_diff_lam_ase[idx] = ene_pot_ase[idx] - ene_pot_lmp[idx]
            if verbose:
                print("DIFF      :",i,ene_diff_lam_ase[idx],ace.units)
                print("--------------------------------------")
                print("--------------------------------------")
                print("MAKE SURE THAT YOU HAVE the correct species in the lammps run!")
            ene_pot[idx] = ene_pot_lmp[idx]

        ### write runner output
        if write_runner:
            ase_write("out.runner",atoms[i],format='runner',append=True)

        ene_diff[idx] = ene_DFT[idx]-ene_pot[idx]
        ene_diff_abs[idx] = np.abs(ene_DFT[idx]-ene_pot[idx])
        ene_mean[idx] = ene_diff[:idx+1].mean()

        if idx == 0:
            ene_std[idx] = 0.
            ene_ste[idx] = 0.
            ene_mean[idx] = ene_diff[idx]
        else:
            ene_std[idx] = np.std(ene_DFT[:idx+1]-ene_pot[:idx+1])
            ene_ste[idx] = ene_std[idx]/np.sqrt(idx)
            ene_mean[idx] = np.mean(np.abs(ene_DFT[:idx+1]-ene_pot[:idx+1]))

        printed = False

        if idx in range(0,len(structures_to_calc),printevery):
            print("%5.0f %5.0f / %6.0f %16.7f =DFT-ref (%s) %4.0f" % (i,idx,len(structures_to_calc),ene_diff_abs[idx],ace.units,atoms[i].get_number_of_atoms()))
            printed = True

        if verbose > 0 and printed == False:
            print("%5.0f %5.0f / %6.0f %16.7f =DFT-ref (%s)" % (i,idx,len(structures_to_calc),ene_diff_abs[idx],ace.units))

    np.savetxt("ene_diff_lam_ase.dat",ene_diff_lam_ase,header=ace.units)

    if write_runner:
        print('our.runner written')



    np.savetxt("ene_DFT.npy",ene_DFT,header=units)
    np.savetxt("ene_pot.npy",ene_pot,header=units)
    np.savetxt("ene_diff.npy",ene_diff,header=units)
    np.savetxt("ene_diff_abs.npy",ene_diff_abs,header=units)
    np.savetxt("ene_std.npy",ene_std,header=units)

    ene_all = np.transpose([range(len(ene_DFT)),ene_DFT,ene_pot,ene_diff_abs,ene_std])
    np.savetxt("ene_all.npy",ene_all,header=units+"\n"+"DFT\t\t"+pot+"\t|diff|\t\t<|diff|>",fmt=' '.join(['%i'] + ['%.10e']*(ene_all.shape[1]-1)))

    print()
    print(ene_diff_abs)
    print(ana_mg_conz)
    analyze = np.transpose([range(len(ene_DFT)),ene_diff_abs,ana_mg_conz,ana_si_conz,ana_al_conz,ana_atoms])
    print('-analyze')
    print(analyze)
    #np.savetxt("analyze.npy",analyze,header=units,fmt=' '.join(['%i'] + ['%.2e']*(analyze.shape[1]-1)))
    #np.savetxt("analyze.npy",analyze,header=units,fmt='%f')
    np.savetxt("analyze.npy",analyze,header=" i diff Mg   Si   Al  atoms",fmt=' '.join(['%4.0f'] +['%5.2f']*(analyze.shape[1]-2)+['%4.0f']))


    my.create_READMEtxt(os.getcwd())
    return

def test_sisivac(ace):
    scripts = my.scripts()
    sisivac = scripts+'/tests/si-si-vac/'
    ff = '/aiida.in.final.runner'
    e_al108 =       ace.ene(ase_read(sisivac+'/al108'+ff))
    e_al107vac1 =   ace.ene(ase_read(sisivac+'/al107vac1'+ff))
    e_al107si1 =    ace.ene(ase_read(sisivac+'/al107si1'+ff))
    e_al106si2 =    ace.ene(ase_read(sisivac+'/al106si2'+ff))
    e_al105si2va1 = ace.ene(ase_read(sisivac+'/al105si2va1'+ff))
    if ace.verbose:
        print('e_al108',e_al108,ace.units)
        print('e_al107vac1',e_al107vac1,ace.units)
    e_ss_al = e_al108/108.
    e_ss_va = e_al107vac1 - e_al108/108. * 107.
    e_ss_one_si = e_al107si1 - e_al108/108. * 107.
    e_si_si_vac_complex = e_al105si2va1 - 2.*e_ss_one_si - e_ss_va - 105.*e_ss_al

    print('vacancy_formation (unrelaxed)',round(e_ss_va,4),"",ace.units,"DFT 0.69 eV")
    print('si-si-vac-complex (unrelaxed)',round(e_si_si_vac_complex,4),ace.units,"DFT -0.074eV")



if __name__ == "__main__":
    get_energies()
