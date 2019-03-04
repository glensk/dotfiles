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

def print_compare_ene_vs_DFT(text,pot_ene,DFT_ene,eos=False,f300=False):
    if type(eos) != bool:
        eos = np.round(eos,1)[1:][:2]
    else: eos = ""

    if type(f300) != bool:
        if type(f300) != str:
            f300 = np.round(f300,1)
    else: f300 = ""

    diff = "-"
    if np.abs(DFT_ene) > 0
        try:
            diff = round(np.abs(np.abs(pot_ene/DFT_ene)-1),2)
        except ZeroDivisionError, RuntimeWarning:
            diff = "-"
    print(text.ljust(35)+":",
            str(round(pot_ene,3)).ljust(10),
            "(eV) DFT:",
            str(round(DFT_ene,4)).ljust(12),
            str(diff).ljust(6),eos,f300)
    return

def test_si_si_vac(ace):
    scripts = my.scripts()
    tests = scripts+'/tests/'
    sisivac = tests+'/Al-Mg-Si/si-si-vac/'
    ff = '/aiida.in.final.runner'
    e_al108 =       ace.ene(ase_read(sisivac+'/al108'+ff))
    e_al107vac1 =   ace.ene(ase_read(sisivac+'/al107vac1'+ff))
    e_al107si1 =    ace.ene(ase_read(sisivac+'/al107si1'+ff))
    e_al106si2 =    ace.ene(ase_read(sisivac+'/al106si2'+ff))
    e_al105si2va1 = ace.ene(ase_read(sisivac+'/al105si2va1'+ff))
    e_ss_al = e_al108/108.
    e_ss_va = e_al107vac1 - e_al108/108. * 107.
    e_ss_one_si = e_al107si1 - e_al108/108. * 107.
    e_si_si_vac_complex = e_al105si2va1 - 2.*e_ss_one_si - e_ss_va - 105.*e_ss_al

    #print_compare_ene_vs_DFT('vacancy_formation (unrelaxed) NN: ',e_ss_va,
    print_compare_ene_vs_DFT("vacancy_formation (unrelaxed, oDFT)",e_ss_va,0.69)
    print_compare_ene_vs_DFT("si-si-vac-complex (unrelaxed, oDFT)",e_si_si_vac_complex,0.074)
    return

def get_dilute_si_mg_f(ace):
    scripts = my.scripts()
    tests = scripts+'/tests/'
    bprime = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/BetaPrime_structures_relax.input.data'

    frames = ase_read(bprime,index=":",format="runner")

    struct_pure_al = frames[3]
    struct_dilute_mg = frames[1]
    struct_dilute_si = frames[2]

    # @ 0K (eV/cell)
    e_pure_al   = ace.ene(struct_pure_al)   # 108 atoms
    e_dilute_mg = ace.ene(struct_dilute_mg) # 108 atoms
    e_dilute_si = ace.ene(struct_dilute_si) # 108 atoms

    # (eV/defect)
    dilute_si_f =  e_dilute_si - 107*e_pure_al/108    # eV per defect
    dilute_mg_f =  e_dilute_mg - 107*e_pure_al/108    # eV per defect

    # @ 300K ( for 1 atom, therefore *108 (to get per cell), in meV, therefore /1000)
    f300_pure_al   = ace.get_fh(struct_pure_al)*108/1000.    # for 1 atom, therefore *108
    f300_dilute_mg = ace.get_fh(struct_dilute_mg)*108/1000.  # for 1 atom, therefore *108
    f300_dilute_si = ace.get_fh(struct_dilute_si)*108/1000.  # for 1 atom, therefore *108

    dilute_si_f_300 =  f300_dilute_si - 107*f300_pure_al/108
    dilute_mg_f_300 =  f300_dilute_mg - 107*f300_pure_al/108

    return dilute_si_f, dilute_mg_f, dilute_si_f_300, dilute_mg_f_300

def test_betaprime_mg9si5(ace):
    scripts = my.scripts()
    tests = scripts+'/tests/'

    bprime = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/BetaPrime_structures_relax.input.data'
    bprime2 = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/aiida_exported_group_BetaPrime_structures_relax.input.data'
    relax_unrelax = [ bprime, bprime+'.v1',bprime+'.v3']
    relax_unrelax = [ bprime+'.v3']

    bprime_stable = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/mg9si5_stable_phonons.runner'
    struct_mg9si5_stable = ase_read(bprime_stable,format="runner")

    for bprime in relax_unrelax:
        #print('bprime input data',bprime)
        frames = ase_read(bprime,index=":",format="runner")

        struct_pure_al = frames[3]
        struct_dilute_mg = frames[1]
        struct_dilute_si = frames[2]
        struct_mg9si5    = frames[0]

        eDFT_pure_al   = my.ase_enepot(struct_pure_al  ,units=ace.units) # 108 atoms
        eDFT_dilute_mg = my.ase_enepot(struct_dilute_mg,units=ace.units) # 108 atoms
        eDFT_dilute_si = my.ase_enepot(struct_dilute_si,units=ace.units) # 108 atoms
        eDFT_precip    = my.ase_enepot(struct_mg9si5   ,units=ace.units) # 28 atoms, 10Si, 18Mg

        e_pure_al   = ace.ene(struct_pure_al) # 108 atoms
        e_dilute_mg = ace.ene(struct_dilute_mg) # 108 atoms
        e_dilute_si = ace.ene(struct_dilute_si) # 108 atoms
        e_precip    = ace.ene(struct_mg9si5) # 28 atoms, 10Si, 18Mg
        e_precip_stable = ace.ene(struct_mg9si5_stable) # 28 atoms, 10Si, 18Mg

        dilute_si_f =  e_dilute_si - 107*e_pure_al/108
        dilute_mg_f =  e_dilute_mg - 107*e_pure_al/108
        heat_precip = (e_precip - 18*dilute_mg_f - 10*dilute_si_f)/28
        heat_precip_stable = (e_precip_stable - 18*dilute_mg_f - 10*dilute_si_f)/28

        dilute_si_f_DFT =  eDFT_dilute_si - 107*eDFT_pure_al/108
        dilute_mg_f_DFT =  eDFT_dilute_mg - 107*eDFT_pure_al/108
        heat_precip_DFT = (eDFT_precip - 18*dilute_mg_f_DFT - 10*dilute_si_f_DFT)/28

        ###### pure al
        vinet_al = ace.get_murn(struct_pure_al) # 28 atoms, 10Si, 18Mg
        f300 = ace.get_fh(struct_pure_al)
        print_compare_ene_vs_DFT("Al",e_pure_al,eDFT_pure_al,vinet_al,f300)

        ###### dilute solution
        print_compare_ene_vs_DFT("one Si in Al",dilute_si_f,dilute_si_f_DFT)
        print_compare_ene_vs_DFT("one Mg in Al",dilute_mg_f,dilute_mg_f_DFT)

        ##### beta prime ... seems is not stable
        vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
        f300 = ace.get_fh(struct_mg9si5,debug=False)
        print_compare_ene_vs_DFT("beta prime Mg9Si5",heat_precip,heat_precip_DFT,vinet_mg9si5,f300)

        ###### beta prime after global minimization with minima hopping
        vinet_mg9si5_stable = ace.get_murn(struct_mg9si5_stable) # 28 atoms, 10Si, 18Mg
        f300 = ace.get_fh(struct_mg9si5_stable,debug=False)
        print_compare_ene_vs_DFT("beta prime Mg9Si5 stable",heat_precip_stable,0,vinet_mg9si5_stable,f300)


        ###### beta prime after global minimization relax everything
        e_dilute_mg_r    = ace.ene(struct_dilute_mg,atomrelax=True) # 108 atoms
        e_dilute_si_r    = ace.ene(struct_dilute_si,atomrelax=True) # 108 atoms

        dilute_si_f_r =  e_dilute_si_r - 107*e_pure_al/108
        dilute_mg_f_r =  e_dilute_mg_r - 107*e_pure_al/108
        heat_precip_stable_allrelax = (e_precip_stable - 18*dilute_mg_f_r - 10*dilute_si_f_r)/28
        print_compare_ene_vs_DFT("beta prime Mg9Si5 stable allrelax",heat_precip_stable_allrelax,0,vinet_mg9si5_stable,f300)
        return dilute_mg_f, dilute_si_f

def test_betaprime_mg9si5_find_global_min(ace,find_global_minimum=True):
    scripts = my.scripts()
    tests = scripts+'/tests/'
    bprime = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/BetaPrime_structures_relax.input.data'
    bprime2 = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/aiida_exported_group_BetaPrime_structures_relax.input.data'
    relax_unrelax = [ bprime, bprime+'.v1',bprime+'.v3']
    relax_unrelax = [ bprime+'.v3']

    bprime_stable = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/mg9si5_stable_phonons.runner'

    for bprime in relax_unrelax:
        #print('bprime input data',bprime)
        frames = ase_read(bprime,index=":",format="runner")

        struct_pure_al = frames[3]
        struct_dilute_mg = frames[1]
        struct_dilute_si = frames[2]
        struct_mg9si5    = frames[0]




        eDFT_pure_al   = my.ase_enepot(struct_pure_al  ,units=ace.units) # 108 atoms
        eDFT_dilute_mg = my.ase_enepot(struct_dilute_mg,units=ace.units) # 108 atoms
        eDFT_dilute_si = my.ase_enepot(struct_dilute_si,units=ace.units) # 108 atoms
        eDFT_precip    = my.ase_enepot(struct_mg9si5   ,units=ace.units) # 28 atoms, 10Si, 18Mg

        e_pure_al   = ace.ene(struct_pure_al) # 108 atoms
        e_dilute_mg = ace.ene(struct_dilute_mg) # 108 atoms
        e_dilute_si = ace.ene(struct_dilute_si) # 108 atoms
        e_precip    = ace.ene(struct_mg9si5) # 28 atoms, 10Si, 18Mg

        dilute_si_f =  e_dilute_si - 107*e_pure_al/108
        dilute_mg_f =  e_dilute_mg - 107*e_pure_al/108
        heat_precip = (e_precip - 18*dilute_mg_f - 10*dilute_si_f)/28
        heat_precip_stable = (e_precip_stable - 18*dilute_mg_f - 10*dilute_si_f)/28

        dilute_si_f_DFT =  eDFT_dilute_si - 107*eDFT_pure_al/108
        dilute_mg_f_DFT =  eDFT_dilute_mg - 107*eDFT_pure_al/108
        heat_precip_DFT = (eDFT_precip - 18*dilute_mg_f_DFT - 10*dilute_si_f_DFT)/28


        vinet_al = ace.get_murn(struct_pure_al) # 28 atoms, 10Si, 18Mg
        f300 = ace.get_fh(struct_pure_al)
        print_compare_ene_vs_DFT("Al",e_pure_al,eDFT_pure_al,vinet_al,f300)


        print_compare_ene_vs_DFT("one Si in Al",dilute_si_f,dilute_si_f_DFT)
        print_compare_ene_vs_DFT("one Mg in Al",dilute_mg_f,dilute_mg_f_DFT)
        vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
        f300 = ace.get_fh(struct_mg9si5,debug=False)
        print_compare_ene_vs_DFT("beta prime Mg9Si5",heat_precip,heat_precip_DFT,vinet_mg9si5,f300)
        if find_global_minimum:
            print("#############################")
            print("FIND GLOBAL MINIMUM")
            print("#############################")



            vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
            f300 = ace.get_fh(struct_mg9si5,debug=False)
            print_compare_ene_vs_DFT("beta prime Mg9Si5",e_precip,0,vinet_mg9si5,f300)


            e_pure_al        = ace.ene(struct_pure_al) # 108 atoms
            e_dilute_mg_r    = ace.ene(struct_dilute_mg,atomrelax=True) # 108 atoms
            e_dilute_si_r    = ace.ene(struct_dilute_si,atomrelax=True) # 108 atoms

            e_precip_relaxed = ace.ene(struct_mg9si5   ,atomrelax=True)
            print('v2',my.ase_vpa(struct_mg9si5),struct_mg9si5.get_potential_energy(),my.ase_mepa(struct_mg9si5))
            print('stress',struct_mg9si5.get_stress())
            print()

            ### try to find global minimum
            print('------- 3')
            e = ace.ene(struct_mg9si5   ,atomrelax=True,cellrelax = True)
            print('v3',my.ase_vpa(struct_mg9si5),struct_mg9si5.get_potential_energy(),my.ase_mepa(struct_mg9si5))
            vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
            f300 = ace.get_fh(struct_mg9si5,debug=False)
            print_compare_ene_vs_DFT("beta prime Mg9Si5 3",e,0,vinet_mg9si5,f300)
            print('stress',struct_mg9si5.get_stress())
            print()

            print('------- 4')
            e = ace.ene(struct_mg9si5   ,atomrelax=True,minimizer='mh')
            print('v4',my.ase_vpa(struct_mg9si5),struct_mg9si5.get_potential_energy(),my.ase_mepa(struct_mg9si5))
            vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
            f300 = ace.get_fh(struct_mg9si5,debug=False)
            print_compare_ene_vs_DFT("beta prime Mg9Si5 4",e,0,vinet_mg9si5,f300)
            print('stress',struct_mg9si5.get_stress())
            print()


            print('------- 5')
            e = ace.ene(struct_mg9si5   ,atomrelax=True,cellrelax=True)
            print('v5',my.ase_vpa(struct_mg9si5),struct_mg9si5.get_potential_energy(),my.ase_mepa(struct_mg9si5))
            vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
            f300 = ace.get_fh(struct_mg9si5,debug=False)
            print_compare_ene_vs_DFT("beta prime Mg9Si5 5",e,0,vinet_mg9si5,f300)
            print('stress',struct_mg9si5.get_stress())
            print()
            ase_write("mg9si5_stable_phonons.runner",struct_mg9si5,format='runner')






        heat_precip_relaxed = (e_precip_relaxed - 18*dilute_mg_f - 10*dilute_si_f)/28
        dilute_si_f_r =  e_dilute_si_r - 107*e_pure_al/108
        dilute_mg_f_r =  e_dilute_mg_r - 107*e_pure_al/108
        heat_precip_rall = (e_precip_relaxed - 18*dilute_mg_f_r - 10*dilute_si_f_r)/28





        print('v3',my.ase_vpa(struct_mg9si5))
        print('fm',my.ase_fmax(struct_mg9si5))
        print('st',struct_mg9si5.get_stress())
        print()
        print()
        print()
        print()
        e_precip_relaxed = ace.ene(struct_mg9si5   ,atomrelax=True,cellrelax=True)
        f300 = ace.get_fh(struct_mg9si5,debug=False)
        print_compare_ene_vs_DFT("beta prime Mg9Si5 ??",heat_precip,heat_precip_DFT,vinet_mg9si5,f300)
        print('v4',my.ase_vpa(struct_mg9si5))
        print('fm',my.ase_fmax(struct_mg9si5))
        print('st',struct_mg9si5.get_stress())
        print()
        print()
        print()
        print()
        e_precip_relaxed = ace.ene(struct_mg9si5   ,atomrelax=True,cellrelax=True)
        f300 = ace.get_fh(struct_mg9si5,debug=False)
        print_compare_ene_vs_DFT("beta prime Mg9Si5 aa",heat_precip,heat_precip_DFT,vinet_mg9si5,f300)
        print('v5',my.ase_vpa(struct_mg9si5))
        print('fm',my.ase_fmax(struct_mg9si5))
        print('st',struct_mg9si5.get_stress())
        print()
        #print_compare_ene_vs_DFT("beta prime Mg9Si5 r",heat_precip_rall,heat_precip_DFT)
        #print_compare_ene_vs_DFT("beta prime Mg9Si5 r",heat_precip_rall,heat_precip_DFT,[ace,struct_mg9si5])
        #print_compare_ene_vs_DFT("beta prime Mg9Si5 r",heat_precip_rall,heat_precip_DFT,[ace,struct_mg9si5])
        return dilute_mg_f, dilute_si_f


def test_Mg2Si(ace):
    dilute_si_f, dilute_mg_f, dilute_si_f_300, dilute_mg_f_300 = get_dilute_si_mg_f(ace)
    print("dilute_si_f:",dilute_si_f)
    print("dilute_mg_f:",dilute_mg_f)
    print("dilute_si_f_300:",dilute_si_f_300)
    print("dilute_mg_f_300:",dilute_mg_f_300)

    # @ T=0K
    scripts = my.scripts()
    tests = scripts+'/tests/'
    path = tests+'/Al-Mg-Si/Mg2Si/POSCAR'
    frame = ase_read(path,format="vasp")
    e_mg2si = ace.ene(frame,atomrelax=True,cellrelax=True)
    heat_precip_T0K         = (e_mg2si -         8*dilute_mg_f - 4*dilute_si_f)/12
    vinet_mg2si = ace.get_murn(frame)
    f300_mg2si = ace.get_fh(frame)
    print_compare_ene_vs_DFT("Mg2Si @0K",heat_precip_T0K,0.0,vinet_mg2si,f300_mg2si)


    # @ T300K
    f300_mg2si_cell = f300_mg2si*frame.get_number_of_atoms()/1000.
    heat_precip_300_corr = (f300_mg2si_cell -         8*dilute_mg_f_300 - 4*dilute_si_f_300)/12
    print_compare_ene_vs_DFT("Mg2Si @300K",heat_precip_T0K+heat_precip_300_corr,0.0,vinet_mg2si,f300_mg2si)
    return

def test_betadoubleprime_mg5si6(ace,dilute_mg_f, dilute_si_f):
    scripts = my.scripts()
    tests = scripts+'/tests/'
    bprimeprime = tests+'/Al-Mg-Si/Mg5Si6_beta_doubleprime/POSCAR'
    frame = ase_read(bprimeprime,format="vasp")
    e_mg5si6 = ace.ene(frame) # 108 atoms
    e_mg5si6_relaxed = ace.ene(frame,atomrelax=True) # 28 atoms, 10Si, 18Mg

    vinet_mg5si6 = ace.get_murn(frame) # 28 atoms, 10Si, 18Mg
    f300 = ace.get_fh(frame)

    heat_precip         = (e_mg5si6 -         10*dilute_mg_f - 12*dilute_si_f)/22
    heat_precip_relaxed = (e_mg5si6_relaxed - 10*dilute_mg_f - 12*dilute_si_f)/22
    print_compare_ene_vs_DFT("beta double prime Mg5Si6",heat_precip,0.0,vinet_mg5si6,f300)
    return

def test_formation_energies(pot,geopt,verbose):
    ace = ase_calculate_ene(pot,units='eV',geopt=geopt,verbose=verbose)
    ace.pot_to_ase_lmp_cmd()  # just to have lmpcmd defined in case ...
                              # ... in case we do test_formation_energies
    #@ test_si_si_vac(ace)
    #@ dilute_mg_f, dilute_si_f = test_betaprime_mg9si5(ace)
    #@ #test_betaprime_mg9si5_find_global_min(ace)
    #@ test_betadoubleprime_mg5si6(ace,dilute_mg_f, dilute_si_f)
    test_Mg2Si(ace)
    return



if __name__ == "__main__":
    get_energies()
