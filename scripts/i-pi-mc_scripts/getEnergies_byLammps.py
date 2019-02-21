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

@click.option('--infile','-i',required=False,type=str,help='input files containing structures that will be imported by ase')
@click.option('--format_in','-fi',type=str,default='runner',help='ase format for reading files')
@click.option('--pot','-p',type=click.Choice(my.pot_all()),required=True,default='n2p2_v1ag')
@click.option('--structures_idx','-idx',default=':',help='which structures to calculate, use ":" for all structues (default), ":3" for structures [0,1,2] etc. (python notation)')
@click.option('--units','-u',type=click.Choice(['eV','meV_pa','eV_pa','hartree','hartree_pa']),default='hartree_pa',help='In which units should the output be given')
@click.option('--geopt/--no-geopt','-g',default=False,help='make a geometry optimization of the atoms.')
@click.option('--ase/--no-ase','-a',default=True,help='Do the calculations by the ase interface to lammps.')
@click.option('--lmp/--no-lmp','-l',default=False,help='Do the calculations externally by lammps and not through ase interface.')
@click.option('--ipi/--no-ipi','-ipi',default=False,help='Do the calculations externally by ipi-lammps and not through ase interface.')
@click.option('--test/--no-test','-t',default=False,help='Assess formation energies of particular test structures.')

@click.option('--pick_concentration_al','-pcal',default=-1.,type=float,help='only consider structures with particular concentration of element, e.g. -pcal 1.0')
@click.option('--pick_atoms_al','-paal',default=-1.,type=float,help='only consider structures with particular number of al atoms, e.g. -paal 106 (e.v. 106 of 108)')
@click.option('--pick_number_of_atoms','-pnat',default=-1.,type=float,help='only consider structures with particular number of atoms, e.g. -pnat 107')
@click.option('--pick_forcesmax','-pfm',default=-1.,type=float,help='only consider structures with particular max force, e.g. -pfm 0')
@click.option('--pick_cellshape','-pcs',default=-1.,type=float,help='only consider structures with particular cellshape, e.g. -pfm 0')

@click.option('--write_runner/--no-write_runner','-wr',required=False,default=False,help='default: runner.out')
@click.option('--verbose','-v',count=True)


def get_energies(infile,format_in,pot,verbose,structures_idx,units,geopt,test,ase,lmp,ipi,write_runner,
        pick_concentration_al,pick_atoms_al,pick_number_of_atoms,pick_forcesmax,pick_cellshape):
    ''' this is a script which computes for a given set of structures the energies
    for a given potential.
    getEnergies_byLammps.py -p n2p2_v1ag --units meV_pa -i input.data -idx 4850:
    getEnergies_byLammps.py -p n2p2_v1ag --units meV_pa -i input.data -idx :4850
    getEnergies_byLammps.py -p n2p2_v1ag --units hartree -i simulation.pos_0.xyz -fi ipi

    '''
    ### when want to assess some formation energies
    if test:
        test_formation_energies(pot,geopt,verbose)
        my.create_READMEtxt(os.getcwd())
        sys.exit('test done! Exit')

    ### check infile
    if not infile:
        sys.exit("Error: Missing option \"--infile\" / \"-i\".")

    #### get the potential
    ace = ase_calculate_ene(pot,units=units,geopt=geopt,verbose=verbose)
    ace.pot_to_ase_lmp_cmd()  # just to have lmpcmd defined in case we do test_formation_energies
    units = ace.units


    ### read in the structures
    my.check_isfile_or_isfiles([infile],verbose=verbose)
    frames = ase_read(infile,index=structures_idx,format=format_in)

    ### print stuff to screen
    print('structures_idx               :',structures_idx)

    #structures_to_calc = my.string_to_index_an_array(range(len(frames)),structures_idx)
    if type(frames) == list:
        structures_to_calc = len(frames)
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
    ana_VOL_diff_norm= np.empty(structures_to_calc);ana_VOL_diff_norm[:]  = np.nan
    ana_vol_pa       = np.empty(structures_to_calc);ana_vol_pa[:]  = np.nan
    ana_dist_min     = np.empty(structures_to_calc);ana_dist_min[:]  = np.nan
    ene_DFT          = np.empty(structures_to_calc);ene_DFT[:]  = np.nan
    ene_DFT_atomic   = np.empty(structures_to_calc);ene_DFT_atomic[:]  = np.nan
    ene_DFT_wo_atomic= np.empty(structures_to_calc);ene_DFT_wo_atomic[:]  = np.nan
    ene_pot          = np.empty(structures_to_calc);ene_pot[:]  = np.nan
    ene_pot_wo_atomic= np.empty(structures_to_calc);ene_pot_wo_atomic[:]  = np.nan
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
    if type(frames) != list:
        frames = [frames]
    for idx,i in enumerate(range(structures_to_calc)):
        ana_atoms_ = frames[i].get_number_of_atoms()
        for_DFTmax_ = np.abs(frames[i].get_forces()).max()
        d = my.ase_get_chemical_symbols_to_conz(frames[i])
        n = my.ase_get_chemical_symbols_to_number_of_species(frames[i])
        cell = frames[i].get_cell()
        pos = frames[i].get_positions()
        cellshape = "?"
        if cell[0,1] == cell[0,2] == cell[1,0] == cell[1,2] == cell[2,0] == cell[2,1] == 0:
            cellshape="R"
            if cell[0,0] == cell[1,1] == cell[2,2]:
                cellshape="Q"
                #if ana_atoms_ == 4
                #if pos[0,0] == pos[0,1] == pos[

        ### when check for particular picks
        if pick_number_of_atoms >= 0 and ana_atoms_ != pick_number_of_atoms:
            continue
        if pick_concentration_al >= 0 and d["Al"] != pick_concentration_al:
            continue
        if pick_atoms_al >= 0 and n["Al"] != pick_atoms_al:
            continue
        if pick_forcesmax ==0 and for_DFTmax_ > 0.001:
            continue
        #if pick_cellshape >= 0 and cellshape not in ["Q", "?"]:
        if pick_cellshape >= 0 and cellshape not in ["Q"]:
            continue
        #print('idx',idx,'n_al',n["Al"],"c_al",d["Al"],'consider_atoms_al cnat_al',consider_atoms_al,'cnat consider_number_of_atoms',consider_number_of_atoms)

        if calc_analysis: ### analysis stuff
            ana_atoms[idx] = ana_atoms_
            for_DFTmax[idx] = for_DFTmax_
            ana_vol[idx] = frames[i].get_volume()
            ana_vol_pa[idx] = frames[i].get_volume()/frames[i].get_number_of_atoms()
            VOL_norm = n["Al"]*16.5+n["Mg"]*22.85+n["Si"]*20.5
            VOL_diff = ana_vol[idx] - VOL_norm
            ana_VOL_diff_norm[idx] = VOL_diff/VOL_norm
            #print('ana',ana_atoms_)
            #print('ka',frames[i].get_all_distances(mic=True))
            #print('kb',np.sort(frames[i].get_all_distances(mic=True))[:,1:])
            #print('kc',np.sort(frames[i].get_all_distances(mic=True))[:,1:].min())
            #print('kd')
            if ana_atoms[idx] > 1:
                ana_dist_min[idx] = np.sort(frames[i].get_all_distances(mic=True))[:,1:].min()
            else:
                #print('frames[i].cell',frames[i].cell)
                #sys.exit()
                ana_dist_min[idx] = -1
            ana_mg_conz[idx] = d["Mg"]
            ana_si_conz[idx] = d["Si"]
            ana_al_conz[idx] = d["Al"]
            at_mg = n["Mg"]
            at_si = n["Si"]
            at_al = n["Al"]





        added=""
        #print(frames[idx].get_positions())
        if calc_DFT: ### ene from DFT
            ene_DFT[idx] = my.ase_enepot(frames[i],units=ace.units)

            if verbose > be_very_verbose:
                my.show_ase_atoms_content(frames[i],showfirst=3,comment = "STAT2")
            if verbose > 1:
                print('ene_DFT[idx]     :',ene_DFT[idx],units)


        if len(ace.units.split("_")) == 1: # per structure
            ene_DFT_atomic[idx] = n["Mg"]*atom_energy_Mg + n["Si"]*atom_energy_Si + n["Al"]*atom_energy_Al
        elif len(ace.units.split("_")) == 2: # per atom
            ene_DFT_atomic[idx] = d["Mg"]*atom_energy_Mg + d["Si"]*atom_energy_Si + d["Al"]*atom_energy_Al
        else:
            sys.exit("either per atom or per structure")

        ene_DFT_wo_atomic[idx] = ene_DFT[idx] - ene_DFT_atomic[idx]

        if ipi == True:  ### ene from ipi
            atoms_tmp = copy.deepcopy(frames[i])
            ene_pot_ipi[idx] = my.ipi_ext_calc(atoms_tmp,ace)

        if ase == True:  ### ene from ase (without writing lammps files)
            atoms_tmp = copy.deepcopy(frames[i])  # for other instances, since frames change when geoopt
            if ace.geopt == False:
                ace.pot_to_ase_lmp_cmd()
                ene_pot_ase[idx] = ace.ene(atoms_tmp)
                if pot == "runner_v2dg" and False: # only in case we load DFT energies from new DFT calcs
                    n = my.ase_get_chemical_symbols_to_number_of_species(atoms[i])
                    ### ene_Al, ene_Mg, ene_Si are per atom, since those are later multi-
                    ### -lied by the number of atoms, this can only work if energies are
                    ### not calculated by _pa
                    ene_Al = 468.845752582        # runner - n2p2: -2.4092354 - -19.6286626 = 17.2194272 hartree == 468.56444 eV
                    ene_Mg = -1247.77831679       # runner - n2p2: -62.6068620 - -16.7493346 =  -45.8575274  hartree == -1247.8468 eV
                    ene_Si = -0.0161424965998     # runner_v2dg-n2p2_v1ag atomic energy: -5.5597835 - -5.5274864 = -.0322971 hartree == -0.87884879

                    ### correction in eV
                    correction = ene_Al*n["Al"] + ene_Si*n["Si"] + ene_Mg*n["Mg"]
                    len_units_split = len(ace.units.split("_"))
                    #print("units_split",units_split)
                    if len_units_split == 2:
                        correction = correction / frames[i].get_number_of_atoms()
                        if ace.units.split("_")[0][0] == 'm':  # meV mhartree
                            correction = correction * 1000.

                    ene_pot_ase[idx] = ene_pot_ase[idx] - correction
                stress = atoms_tmp.get_stress()
                #print('stress',stress)

            if ace.geopt == True:
                ace.geopt = False
                ace.pot_to_ase_lmp_cmd()
                ene_pot_ase[idx] = ace.ene(atoms_tmp)
                #print(frames[i].get_positions())

                ace.geopt = True
                ace.pot_to_ase_lmp_cmd()
                ene_pot_ase_geop[idx] = ace.ene(atoms_tmp)
                #print(ace.atomsin.get_positions())
                #print(ace.frames.get_positions())
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
            atoms_tmp = copy.deepcopy(frames[i])  # for other instances, since atoms change when geoopt
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
            ase_write("out.runner",frames[i],format='runner',append=True)
        ene_pot_wo_atomic[idx] = ene_pot[idx] - ene_DFT_atomic[idx]
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
            fmt_one = '%10.'+str(show)+'f'
            fmt_after_atms=' '.join([fmt_one]*8)   # add here if a new entry
            ka3="%5.0f %5.0f / %6.0f "+cellshape+" "+fmt_one+" [%4.0f %4.0f %4.0f %4.0f] "+fmt_after_atms+" "+added
            print(ka3 % (i,idx,structures_to_calc,ene_diff_abs[idx],frames[i].get_number_of_atoms(),n["Si"],n["Mg"],n["Al"],ene_DFT[idx],ene_pot[idx],ene_DFT_wo_atomic[idx],for_DFTmax[idx],ene_pot_ase[idx]-ene_pot_ase_geop[idx],ana_vol_pa[idx],ana_dist_min[idx],ana_VOL_diff_norm[idx]))
            return


        if idx in range(0,structures_to_calc,printevery):
            printhere()
            printed = True

        if verbose > 0 and printed == False:
            printhere()

    np.savetxt("ene_diff_lam_ase.dat",ene_diff_lam_ase,header=ace.units)

    if write_runner:
        print('our.runner written')

    def mysavetxt(what,name,units,save=False):
        whatout = what[~np.isnan(what)]
        if save:
            np.savetxt(name,whatout,header=units)

        #print("in",what.shape,"out",whatout.shape,name)
        return whatout

    ene_DFT         = mysavetxt(ene_DFT,"ene_DFT.npy",units,save=True)
    ene_pot         = mysavetxt(ene_pot,"ene_pot.npy",units,save=True)
    ene_diff        = mysavetxt(ene_diff,"ene_diff.npy",units,save=True)
    ene_diff_abs    = mysavetxt(ene_diff_abs,"ene_diff_abs.npy",units,save=True)
    ene_std         = mysavetxt(ene_std,"ene_std.npy",units,save=True)

    ene_DFT_wo_atomic = mysavetxt(ene_DFT_wo_atomic,"ene_DFT_wo_atomic",units,save=True)
    ene_pot_wo_atomic = mysavetxt(ene_pot_wo_atomic,"ene_pot_wo_atomic",units,save=True)
    for_DFTmax      = mysavetxt(for_DFTmax,"for_DFTmax",units)
    ana_mg_conz     = mysavetxt(ana_mg_conz,"ana_mg_conz",units)
    ana_si_conz     = mysavetxt(ana_si_conz,"ana_si_conz",units)
    ana_al_conz     = mysavetxt(ana_al_conz,"ana_al_conz",units)
    ana_atoms       = mysavetxt(ana_atoms,"ana_atoms",units)
    ana_vol         = mysavetxt(ana_vol ,"ana_vol",units)
    ana_vol_pa      = mysavetxt(ana_vol_pa ,"ana_vol_pa",units)
    ana_VOL_diff_norm = mysavetxt(ana_VOL_diff_norm,"ana_VOL_diff_norm",units)
    ana_dist_min    = mysavetxt(ana_dist_min,"ana_dist_min",units)
    if len(ene_pot) != 0:
        np.savetxt("ene_DFT.npy",ene_DFT,header=units)
        np.savetxt("ene_pot.npy",ene_pot,header=units)
        np.savetxt("ene_diff.npy",ene_diff,header=units)
        np.savetxt("ene_diff_abs.npy",ene_diff_abs,header=units)
        np.savetxt("ene_std.npy",ene_std,header=units)

        ene_all = np.transpose([range(len(ene_DFT)),ene_DFT,ene_pot,ene_diff_abs,ene_std])
        ### write analyze.csv
        try:
            np.savetxt("ene_all.npy",ene_all,header=units+"\n"+"DFT\t\t"+pot+"\t|diff|\t\t<|diff|>",fmt=' '.join(['%i'] + ['%.10e']*(ene_all.shape[1]-1)))
        except IndexError:
            print('len',len(ene_DFT))
            print(ene_DFT.shape)
            print(ene_diff_abs.shape)
            print(ene_DFT_wo_atomic.shape)
            print(ene_pot_wo_atomic.shape)
            print(for_DFTmax.shape)
            print(ana_mg_conz.shape)
            print(ana_si_conz.shape)
            print(ana_al_conz.shape)
            print(ana_atoms.shape)
            print(ana_vol.shape)
            print(ana_vol_pa.shape)
            print(ana_dist_min.shape)
            print('cant save ene_all.npy')

        analyze = np.transpose([
            np.arange(len(ene_DFT)),   # i
            ene_diff_abs,          # diff
            ene_DFT_wo_atomic,     # E_wo
            for_DFTmax,            # for
            ana_mg_conz,
            ana_si_conz,
            ana_al_conz,
            ana_atoms,
            ana_vol,
            ana_vol_pa,
            ana_dist_min,
            ana_VOL_diff_norm])
        #print('a',analyze.shape)
        #print('a',analyze.shape[0])
        #analyze_len = analyze.shape[1] - 1
        analyze_len = analyze.shape[0] - 1
        #np.savetxt("analyze.csv",analyze ,delimiter=',',header=" i   diff  E_wo    for_max  Mg_c   Si_c   Al_c  atoms   vol  vol_pa dist_min") # ,fmt=' '.join(['%4.0f'] +['%6.2f']*analyze_len))
        np.savetxt("analyze.csv",analyze ,delimiter=',') # ,fmt=' '.join(['%4.0f'] +['%6.2f']*analyze_len))


    my.create_READMEtxt(os.getcwd())
    return


def printhead(structures_to_calc,ace_units):
    print('structures_to_calc[:3]:',range(structures_to_calc)[:3],'...',range(structures_to_calc)[-3:])
    print()
    print('#                         ('+ace_units+')                        ('+ace_units+')     ('+ace_units+')    ')
    print('#                         (DFT-ref)                                            ene_wo_atomic    forces    (if geopt)  Vol per')
    print('#   i   idx /    from       diff  [atms   Si   Mg   Al]   ene_DFT     ene_pot                  DFTmax     E-E_geopt   atom')
    print('--------------------------------------------------------------------------------------------------------------------------------------------------------------')
    return

def test_formation_energies(pot,geopt,verbose):
    ace = ase_calculate_ene(pot,units='eV',geopt=geopt,verbose=verbose)
    print('changed units to eV!')
    ace.pot_to_ase_lmp_cmd()  # just to have lmpcmd defined in case we do test_formation_energies
    units = ace.units

    scripts = my.scripts()
    tests = scripts+'/tests/'

    if True: # si-si-vac
        sisivac = tests+'/Al-Mg-Si/si-si-vac/'
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
        print()
        print('vacancy_formation (unrelaxed) NN: ',round(e_ss_va,4),"",ace.units,"DFT: 0.69 eV")
        print('si-si-vac-complex (unrelaxed) NN:',round(e_si_si_vac_complex,4),"",ace.units,"DFT: -0.074eV")
        print()

    if True: # bprime (beta')
        print("########################################################################")
        print("# beta' (Mg9Si5)")
        print("########################################################################")
        bprime = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/BetaPrime_structures_relax.input.data'
        print('bprime input data',bprime)
        frames = ase_read(bprime,index=":",format="runner")
        #print('len',len(frames))
        for idx,i in enumerate(frames):
            ana_atoms_ = frames[idx].get_number_of_atoms()
            d = my.ase_get_chemical_symbols_to_conz(frames[idx])
            n = my.ase_get_chemical_symbols_to_number_of_species(frames[idx])
            print(idx,'nat',ana_atoms_,"nAl:",n["Al"],"nSi:",n["Si"],'nMg:',n["Mg"])

        eDFT_pure_al   = my.ase_enepot(frames[3],units=ace.units) # 108 atoms
        eDFT_dilute_mg = my.ase_enepot(frames[1],units=ace.units) # 108 atoms
        eDFT_dilute_si = my.ase_enepot(frames[2],units=ace.units) # 108 atoms
        eDFT_precip    = my.ase_enepot(frames[0],units=ace.units) # 28 atoms, 10Si, 18Mg
        print('eDFT_pure_al  :',eDFT_pure_al  ,"(eV)",eDFT_pure_al  *0.036749325,"(hartree)")
        print('eDFT_dilute_mg:',eDFT_dilute_mg,"(eV)",eDFT_dilute_mg*0.036749325,"(hartree)")
        print('eDFT_dilute_si:',eDFT_dilute_si,"(eV)",eDFT_dilute_si*0.036749325,"(hartree)")
        print('eDFT_precip   :',eDFT_precip   ,"(eV)",eDFT_precip   *0.036749325,"(hartree)")
        print()
        print('vol  pure  /at:',frames[3].get_volume()/frames[3].get_number_of_atoms())
        print('vol  dilMg /at:',frames[1].get_volume()/frames[1].get_number_of_atoms())
        print('vol  dilSi /at:',frames[2].get_volume()/frames[2].get_number_of_atoms())
        print('vol  precip/at:',frames[0].get_volume()/frames[0].get_number_of_atoms())

        e_pure_al   = ace.ene(frames[3]) # 108 atoms
        e_dilute_mg = ace.ene(frames[1]) # 108 atoms
        e_dilute_si = ace.ene(frames[2]) # 108 atoms
        e_precip    = ace.ene(frames[0]) # 28 atoms, 10Si, 18Mg

        print('e_pure_al  :',e_pure_al  ,"(eV)",e_pure_al  *0.036749325,"(hartree)")
        print('e_dilute_mg:',e_dilute_mg,"(eV)",e_dilute_mg*0.036749325,"(hartree)")
        print('e_dilute_si:',e_dilute_si,"(eV)",e_dilute_si*0.036749325,"(hartree)")
        print('e_precip   :',e_precip   ,"(eV)",e_precip   *0.036749325,"(hartree)")


        dilute_si_f =  e_dilute_si - 107*e_pure_al/108
        dilute_mg_f =  e_dilute_mg - 107*e_pure_al/108
        heat_precip = (e_precip - 18*dilute_mg_f - 10*dilute_si_f)/28

        dilute_si_f_DFT =  eDFT_dilute_si - 107*eDFT_pure_al/108
        dilute_mg_f_DFT =  eDFT_dilute_mg - 107*eDFT_pure_al/108
        heat_precip_DFT = (eDFT_precip - 18*dilute_mg_f_DFT - 10*dilute_si_f_DFT)/28

        print("dilute_si_f pot:",round(dilute_si_f,3)," DFT:",round(dilute_si_f_DFT,3))
        print("dilute_mg_f pot:",round(dilute_mg_f,3)," DFT:",round(dilute_mg_f_DFT,3))
        print("heat_precip pot:",round(heat_precip,3)," DFT:",round(heat_precip_DFT,3))
        print()

        print('bulk')
        ace.murn(frames[3]) # 108 atoms
        print('precip')
        ace.murn(frames[0])

    if True: # beta''
        print("########################################################################")
        print("# beta'' (Mg5Si6)")
        print("########################################################################")
        bprimeprime = tests+'/Al-Mg-Si/beta_prime/BetaPrime_structures_relax.input.data'
        print('bprime input data',bprime)
        frames = ase_read(bprime,index=":",format="runner")

    return



if __name__ == "__main__":
    get_energies()
