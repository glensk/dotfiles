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
@click.option('--consider_concentration_al','-cal',default=-1.,type=float,help='only consider structures with particular concentration of element, e.g. -cc Al 1.0')
@click.option('--write_runner/--no-write_runner','-wr',required=False,default=False,help='default: runner.out')
@click.option('--verbose','-v',count=True)


def get_energies(infile,format_in,pot,verbose,structures_idx,units,geopt,test,ase,lmp,ipi,write_runner,consider_concentration_al):
    ''' this is a script which computes for a given set of structures the energies
    for a given potential.
    getEnergies_byLammps.py -p n2p2_v1ag --units meV_pa -i input.data -idx 4850:
    getEnergies_byLammps.py -p n2p2_v1ag --units meV_pa -i input.data -idx :4850
    getEnergies_byLammps.py -p n2p2_v1ag --units hartree -i simulation.pos_0.xyz -fi ipi

    '''

    #### get the potential
    ace = ase_calculate_ene(pot,units=units,geopt=geopt,verbose=verbose)
    ace.pot_to_ase_lmp_cmd()  # just to have lmpcmd defined in case we do test_sisivac
    units = ace.units

    ### when want to assess some formation energies
    if test:
        test_sisivac(ace)
        sys.exit('test done! Exit')

    ### read in the structures
    my.check_isfile_or_isfiles([infile],verbose=verbose)
    atoms = ase_read(infile,index=structures_idx,format=format_in)

    ### print stuff to screen
    print('structures_idx               :',structures_idx)

    #structures_to_calc = my.string_to_index_an_array(range(len(atoms)),structures_idx)
    if type(atoms) == list:
        structures_to_calc = len(atoms)
    else:
        structures_to_calc = 1

    print('structures_to_calc           :',structures_to_calc)
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

    ana_mg_conz      = np.empty(structures_to_calc);ana_mg_conz[:]  = np.nan
    ana_si_conz      = np.empty(structures_to_calc);ana_si_conz[:]  = np.nan
    ana_al_conz      = np.empty(structures_to_calc);ana_al_conz[:]  = np.nan
    ana_atoms        = np.empty(structures_to_calc);ana_atoms[:]  = np.nan
    ana_vol          = np.empty(structures_to_calc);ana_vol[:]  = np.nan
    ana_vol_pa       = np.empty(structures_to_calc);ana_vol_pa[:]  = np.nan
    ana_dist_min     = np.empty(structures_to_calc);ana_dist_min[:]  = np.nan
    ene_DFT          = np.empty(structures_to_calc);ene_DFT[:]  = np.nan
    ene_DFT_wo_atomic= np.empty(structures_to_calc);ene_DFT_wo_atomic[:]  = np.nan
    ene_pot          = np.empty(structures_to_calc);ene_pot[:]  = np.nan
    ene_pot_ase      = np.empty(structures_to_calc);ene_pot_ase[:]  = np.nan
    ene_pot_ase_geop = np.empty(structures_to_calc);ene_pot_ase_geop[:]  = np.nan
    ene_pot_lmp      = np.empty(structures_to_calc);ene_pot_lmp[:]  = np.nan
    ene_pot_ipi      = np.empty(structures_to_calc);ene_pot_ipi[:]  = np.nan
    ene_diff         = np.empty(structures_to_calc);ene_diff[:]  = np.nan
    ene_diff_abs     = np.empty(structures_to_calc);ene_diff_abs[:]  = np.nan
    ene_std          = np.empty(structures_to_calc);ene_std[:]  = np.nan
    ene_ste          = np.empty(structures_to_calc);ene_ste[:]  = np.nan
    ene_mean         = np.empty(structures_to_calc);ene_mean[:] = np.nan
    ene_diff_lam_ase = np.empty(structures_to_calc);ene_DFT[:]  = np.nan
    for_DFTmax       = np.empty(structures_to_calc);for_DFTmax[:]  = np.nan
    #sys.exit('get uuid of structure and save structure energy somewhere (cache)')
    #sys.exit('find out weather particular structure in test or trainset')
    # make this parallel at some point

    ##############################################
    # loop over structures
    ##############################################
    printevery = 50
    if structures_to_calc < 50 or verbose > 0:
        printevery = 1
    be_very_verbose = 999

    ##############################################
    # get atomic energies to substract from total energies
    ##############################################
    if ace.units.split("_")[0] == 'hartree':
        conv = 1.
    elif ace.units.split("_")[0] == 'ev':
        conv = 27.211384    # hartree to ev
    elif ace.units.split("_")[0] == 'mev':
        conv = 27211.384    # hartree to ev
    else:
        print('ace units',ace.units.split("_"))
        sys.exit('ace units not known')
    atom_energy_Mg = ace.mypot.atom_energy["Mg"]*conv
    atom_energy_Si = ace.mypot.atom_energy["Si"]*conv
    atom_energy_Al = ace.mypot.atom_energy["Al"]*conv
    #print("atom_energy_Mg",atom_energy_Mg,ace.units,"one_atom",ace.mypot.atom_energy["Mg"],"hartee_pa")
    #print("atom_energy_Si",atom_energy_Si,ace.units,"one_atom",ace.mypot.atom_energy["Si"],"hartree_pa")
    #print("atom_energy_Al",atom_energy_Al,ace.units,"one_atom",ace.mypot.atom_energy["Al"],"hartree_pa")


    calc_DFT = True
    calc_pot = True
    calc_analysis = True
    calc_ipi = False
    printhead(structures_to_calc,ace.units)
    goover = [1]
    if ace.geopt == True : goover = [1,2]
    if type(atoms) != list:
        atoms = [atoms]
    for idx,i in enumerate(range(structures_to_calc)):
        ### when check for particular concentration
        d = my.ase_get_chemical_symbols_to_conz(atoms[i])
        if consider_concentration_al >= 0:
            if d["Al"] != consider_concentration_al:
                continue


        added=""
        #print(atoms[idx].get_positions())
        if calc_DFT: ### ene from DFT
            ene_DFT[idx] = my.ase_enepot(atoms[i],units=ace.units)

            if verbose > be_very_verbose:
                my.show_ase_atoms_content(atoms[i],showfirst=3,comment = "STAT2")
            if verbose > 1:
                print('ene_DFT[idx]     :',ene_DFT[idx],units)

        if calc_analysis: ### analysis stuff
            ana_vol[idx] = atoms[i].get_volume()
            ana_vol_pa[idx] = atoms[i].get_volume()/atoms[i].get_number_of_atoms()
            ana_dist_min[idx] = np.sort(atoms[i].get_all_distances(mic=True))[:,1:].min()
            n = my.ase_get_chemical_symbols_to_number_of_species(atoms[i])
            ana_mg_conz[idx] = d["Mg"]
            ana_si_conz[idx] = d["Si"]
            ana_al_conz[idx] = d["Al"]
            at_mg = n["Mg"]
            at_si = n["Si"]
            at_al = n["Al"]
            ana_atoms[idx]   = atoms[i].get_number_of_atoms()

            #### analysisstuff when for whole structure (if per atom then it is necessary to use concentration (d) instead of n["Al"]
            if len(ace.units.split("_")) == 1: # per structure
                ene_DFT_wo_atomic[idx] = ene_DFT[idx] - n["Mg"]*atom_energy_Mg - n["Si"]*atom_energy_Si - n["Al"]*atom_energy_Al
            elif len(ace.units.split("_")) == 2: # per atom
                ene_DFT_wo_atomic[idx] = ene_DFT[idx] - d["Mg"]*atom_energy_Mg - d["Si"]*atom_energy_Si - d["Al"]*atom_energy_Al
            else:
                sys.exit("either per atom or per structure")

            ### analysis forces
            for_DFTmax[idx] = np.abs(atoms[i].get_forces()).max()

        if ipi == True:  ### ene from ipi
            atoms_tmp = copy.deepcopy(atoms[i])
            ene_pot_ipi[idx] = my.ipi_ext_calc(atoms_tmp,ace)

        if ase == True:  ### ene from ase (without writing lammps files)
            atoms_tmp = copy.deepcopy(atoms[i])  # for other instances, since atoms change when geoopt
            if ace.geopt == False:
                ace.pot_to_ase_lmp_cmd()
                ene_pot_ase[idx] = ace.ene(atoms_tmp)

            if ace.geopt == True:
                ace.geopt = False
                ace.pot_to_ase_lmp_cmd()
                ene_pot_ase[idx] = ace.ene(atoms_tmp)
                #print(atoms[i].get_positions())

                ace.geopt = True
                ace.pot_to_ase_lmp_cmd()
                ene_pot_ase_geop[idx] = ace.ene(atoms_tmp)
                #print(ace.atomsin.get_positions())
                #print(ace.atoms.get_positions())
                if idx == 0 and os.path.isfile("out_relaxed.runner"):
                    my.rm("out_relaxed.runner")
                if abs(ene_pot_ase[idx]-ene_pot_ase_geop[idx]) > 0.01: # meV_pa
                    #print('adding idx',idx)
                    ase_write("out_relaxed.runner",ace.atoms,format='runner',append=True)
                    added="*"
                #print('ee',ene_pot_ase[idx],ene_pot_ase_geop[idx])

            if verbose > 1:
                print('xx ene_pot_ase[idx] :',ene_pot_ase[idx],units)

            if lmp == False:
                ene_pot[idx] = copy.deepcopy(ene_pot_ase[idx])

        if lmp == True:  ### ene from lammps (by writing lammps files)
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
            x = ene_DFT[:idx+1]-ene_pot[:idx+1]
            x = x[~np.isnan(x)]
            ene_std[idx] = np.std(x)
            ene_ste[idx] = ene_std[idx]/np.sqrt(idx)
            ene_mean[idx] = np.mean(np.abs(x))

        printed = False

        def printhere():
            show = 4
            # enough for meV_pa
            if ace.units == "mev_pa":
                show = 3

            fmt_one = '%16.'+str(show)+'f'
            fmt=' '.join(['%16.'+str(show)+'f']*4)
            fmt=' '.join([fmt_one]*4)
            fmt=' '.join([fmt_one]*5)
            ka1="%5.0f %5.0f / %6.0f %16.2f =DFT-ref (%s) [%4.0f] %16.2f %16.2f %16.2f %16.2f"
            ka2="%5.0f %5.0f / %6.0f %16.2f =DFT-ref (%s) [%4.0f] "+fmt
            ka3="%5.0f %5.0f / %6.0f "+fmt_one+" =DFT-ref (%s) [%4.0f] "+fmt+" "+added
            #print("%5.0f %5.0f / %6.0f %16.2f =DFT-ref (%s) [%4.0f] %16.2f %16.2f %16.2f %16.2f" % (i,idx,structures_to_calc,ene_diff_abs[idx],ace.units,atoms[i].get_number_of_atoms(),ene_DFT[idx],ene_pot[idx],ene_DFT_wo_atomic[idx],for_DFTmax[idx]))
            print(ka3 % (i,idx,structures_to_calc,ene_diff_abs[idx],ace.units,atoms[i].get_number_of_atoms(),ene_DFT[idx],ene_pot[idx],ene_DFT_wo_atomic[idx],for_DFTmax[idx],ene_pot_ase[idx]-ene_pot_ase_geop[idx]))
            return


        if idx in range(0,structures_to_calc,printevery):
            printhere()
            printed = True

        if verbose > 0 and printed == False:
            printhere()

    np.savetxt("ene_diff_lam_ase.dat",ene_diff_lam_ase,header=ace.units)

    if write_runner:
        print('our.runner written')

    def mysavetxt(what,name,units,save=True):
        what = what[~np.isnan(what)]
        if save:
            np.savetxt(name,what,header=units)
        return what

    mysavetxt(ene_DFT,"ene_DFT.npy",units,save=True)
    mysavetxt(ene_pot,"ene_pot.npy",units,save=True)
    mysavetxt(ene_diff,"ene_diff.npy",units,save=True)
    mysavetxt(ene_diff_abs,"ene_diff_abs.npy",units,save=True)
    mysavetxt(ene_std,"ene_std.npy",units,save=True)
    ene_DFT_wo_atomic = mysavetxt(ene_DFT_wo_atomic,"",units,save=False)
    ene_DFTmax = mysavetxt(ene_DFTmax,"",units,save=False)
    ana_mg_conz = mysavetxt(ana_mg_conz,"",units,save=False)
    ana_si_conz = mysavetxt(ana_si_conz,"",units,save=False)
    ana_al_conz = mysavetxt(ana_al_conz,"",units,save=False)
    ana_atoms = mysavetxt(ana_atoms,"",units,save=False)
    ana_vol = mysavetxt(ana_vol ,"",units,save=False)
    ana_vol_pa = mysavetxt(ana_vol_pa ,"",units,save=False)
    ana_dist_min = mysavetxt(ana_dist_min,"",units,save=False)

    np.savetxt("ene_DFT.npy",ene_DFT,header=units)
    np.savetxt("ene_pot.npy",ene_pot,header=units)
    np.savetxt("ene_diff.npy",ene_diff,header=units)
    np.savetxt("ene_diff_abs.npy",ene_diff_abs,header=units)
    np.savetxt("ene_std.npy",ene_std,header=units)

    ene_all = np.transpose([range(len(ene_DFT)),ene_DFT,ene_pot,ene_diff_abs,ene_std])
    np.savetxt("ene_all.npy",ene_all,header=units+"\n"+"DFT\t\t"+pot+"\t|diff|\t\t<|diff|>",fmt=' '.join(['%i'] + ['%.10e']*(ene_all.shape[1]-1)))

    ### write analyze.csv
    analyze = np.transpose([
        np.arange(len(ene_DFT)),   # i
        ene_diff_abs,          # diff
        ene_DFT_wo_atomic,     # E_wo
        for_DFTmax,            # for
        ana_mg_conz,ana_si_conz,ana_al_conz,ana_atoms,ana_vol,ana_vol_pa,ana_dist_min])
    analyze_len = analyze.shape[1] - 1
    np.savetxt("analyze.csv",analyze ,delimiter=',',header=" i   diff  E_wo    for_max  Mg_c   Si_c   Al_c  atoms   vol  vol_pa dist_min") # ,fmt=' '.join(['%4.0f'] +['%6.2f']*analyze_len))


    my.create_READMEtxt(os.getcwd())
    return


def printhead(structures_to_calc,ace_units):
    print('structures_to_calc[:3]:',range(structures_to_calc)[:3],'...',range(structures_to_calc)[-3:])
    print()
    print('    i   idx /    from       diff      =DFT-ref ('+ace_units+') [atms]  ene_DFT('+ace_units+')  ene_pot('+ace_units+')   ene_wo_atomic       for_DFT_max            E-E_geopt')
    print('--------------------------------------------------------------------------------------------------------------------------------------------------------------')
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

    print('vacancy_formation (unrelaxed) NN: ',round(e_ss_va,4),"",ace.units,"DFT: 0.69 eV")
    print('si-si-vac-complex (unrelaxed) NN:',round(e_si_si_vac_complex,4),"",ace.units,"DFT: -0.074eV")



if __name__ == "__main__":
    get_energies()
