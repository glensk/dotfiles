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
@click.option('--pot','-p',type=click.Choice(my.pot_all()),required=True,default=my.get_latest_n2p2_pot())
@click.option('--potpath','-pp',type=str,required=False,default=False,help="In case --pot is set to setpath use --potpath Folder to point to the Folder containing the n2p2/runner potential")
@click.option('--structures_idx','-idx',default=':',help='which structures to calculate, use ":" for all structues (default), ":3" for structures [0,1,2] etc. (python notation)')
@click.option('--units','-u',type=click.Choice(['eV','meV_pa','eV_pa','hartree','hartree_pa']),default='hartree_pa',help='In which units should the output be given')
@click.option('--geopt/--no-geopt','-g',default=False,help='make a geometry optimization of the atoms.')
@click.option('--elastic/--no-elastic','-e',default=False,help='calculate elastic constants (externally by lammps).')
@click.option('--ase/--no-ase','-a',default=True,help='Do the calculations by the ase interface to lammps.')
@click.option('--lmp/--no-lmp','-l',default=False,help='Do the calculations externally by lammps and not through ase interface.')
@click.option('--ipi/--no-ipi','-ipi',default=False,help='Do the calculations externally by ipi-lammps and not through ase interface.')

@click.option('--test/--no-test','-t',default=False,help='Assess formation energies of particular test structures.')
@click.option('--teste/--no-teste','-te',default=False,help='Assess elastic constants.')
@click.option('--test3/--no-test3','-t3',default=False,help='test3')

@click.option('--pick_concentration_al','-pcal',default=-1.,type=float,help='only consider structures with particular concentration of element, e.g. -pcal 1.0')
@click.option('--pick_atoms_al','-paal',default=-1.,type=float,help='only consider structures with particular number of al atoms, e.g. -paal 106 (e.v. 106 of 108)')
@click.option('--pick_number_of_atoms','-pnat',default=-1.,type=float,help='only consider structures with particular number of atoms, e.g. -pnat 107')
@click.option('--pick_forcesmax','-pfm',default=-1.,type=float,help='only consider structures with particular max force, e.g. -pfm 0')
@click.option('--pick_cellshape','-pcs',default=-1.,type=float,help='only consider structures with particular cellshape, e.g. -pfm 0')

@click.option('--write_runner/--no-write_runner','-wr',required=False,default=False,help='default: runner.out')
@click.option('--verbose','-v',count=True)


def get_energies(infile,format_in,pot,potpath,verbose,structures_idx,units,geopt,elastic,test,teste,test3,ase,lmp,ipi,write_runner,
        pick_concentration_al,pick_atoms_al,pick_number_of_atoms,pick_forcesmax,pick_cellshape):
    ''' this is a script which computes for a given set of structures the energies
    for a given potential.
    getEnergies_byLammps.py -p n2p2_v1ag --units meV_pa -i input.data -idx 4850:
    getEnergies_byLammps.py -p n2p2_v1ag --units meV_pa -i input.data -idx :4850
    getEnergies_byLammps.py -p n2p2_v1ag --units hartree -i simulation.pos_0.xyz -fi ipi

    '''
    print('infile           :',infile)

    ### check if ase runneer format is known
    my.ase_get_known_formats(show=False,add_missing_formats=False, copy_formats=True,verbose=False)

    ### check if lammps is working with ase
    if 'LD_LIBRARY_PATH' not in os.environ:
        os.environ['LD_LIBRARY_PATH'] = os.environ['HOME']+'/sources/lammps/src'
    print('LD_LIBRARY_PATH  :',os.environ['LD_LIBRARY_PATH'])
    from lammps import lammps
    lammps()

    ### get ace object for the chosen potential
    if test or teste: units='eV'
    ace = ase_calculate_ene(pot=pot,
            potpath=potpath,
            units=units,
            geopt=geopt,
            elastic=elastic,
            verbose=verbose)
    ace.pot_to_ase_lmp_cmd()  # just to have lmpcmd defined in case ...
    units = ace.units
    print('ace.pot.pot      :',ace.pot.pot)
    print('ace.pot.fullpath :',ace.pot.potpath)
    ace.pot.print_variables(print_nontheless=True,text=">>")

    ### when want to assess some formation energies
    if test:
        print('>> test')
        test_formation_energies(ace)
        my.create_READMEtxt(os.getcwd())
        sys.exit('test done! Exit')

    if teste or elastic:
        print('>> elastic')
        ace.elastic = True
        test2_elastic(ace)
        my.create_READMEtxt(os.getcwd())
        sys.exit('teste done! Exit')

    if test3:
        print('>> test3')
        test3_do(ace)
        sys.exit('test3 done! Exit')


    ### check infile
    if not infile:
        sys.exit("Error: Missing option \"--infile\" / \"-i\".")

    ### read in the structures
    print('reading infile ...',infile)
    if infile == 'POSCAR': format_in = "vasp"

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
    if structures_to_calc == 0:
        print()
        print("ERROR!")
        sys.exit('0 strucutres to calculate? Something went wrong when importing the infile \"'+infile+'\" (mayby you need to change the \"format_in\" of the file? currently \"'+format_in+'\"')

    # show positions of first structure?
    show_positions = False
    if show_positions:
        print(frames[0].positions)
        print()
        print(frames[0].cell)

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
    print('--pick_concentration_al      :',pick_concentration_al)
    print('--pick_atoms_al              :',pick_atoms_al)
    print('--pick_number_of_atoms       :',pick_number_of_atoms)
    print('--pick_forcesmax             :',pick_forcesmax)
    print('--pick_cellshape             :',pick_cellshape)
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
    atom_energy_Mg = ace.pot.atom_energy["Mg"]*conv
    atom_energy_Si = ace.pot.atom_energy["Si"]*conv
    atom_energy_Al = ace.pot.atom_energy["Al"]*conv
    #print("atom_energy_Mg",atom_energy_Mg,ace.units,"one_atom",ace.pot.atom_energy["Mg"],"hartee_pa")
    #print("atom_energy_Si",atom_energy_Si,ace.units,"one_atom",ace.pot.atom_energy["Si"],"hartree_pa")
    #print("atom_energy_Al",atom_energy_Al,ace.units,"one_atom",ace.pot.atom_energy["Al"],"hartree_pa")


    calc_DFT = True
    calc_pot = True
    calc_analysis = True
    calc_ipi = False
    printhead(structures_to_calc,ace.units)
    goover = [1]
    if ace.geopt == True : goover = [1,2]
    if type(frames) != list:
        frames = [frames]
    max_at = 0
    max_at_id = 0
    max_at_orig = 0
    min_at = 100000
    min_at_id = 0
    min_at_orig = 0
    for idx,i in enumerate(range(structures_to_calc)):
        #print('iii',i)
        #print('kkk',frames[i].get_forces())
        ana_atoms_ = frames[i].get_number_of_atoms()
        #print('ama_atoms_',ana_atoms_)
        if verbose > 2:
            print('i',i,'ana_atoms',ana_atoms_)
        for_DFTmax_ = np.abs(frames[i].get_forces()).max()
        #print('kk',for_DFTmax_)
        d = my.ase_get_chemical_symbols_to_conz(frames[i])
        if verbose > 2:
            print('i',i,'d',d)
        n = my.ase_get_chemical_symbols_to_number_of_species(frames[i])
        if verbose > 2:
            print('i',i,'n',n)
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
        if pick_forcesmax >=0 and for_DFTmax_ > pick_forcesmax:
            continue
        #if pick_cellshape >= 0 and cellshape not in ["Q", "?"]:
        if pick_cellshape >= 0 and cellshape not in ["Q"]:
            continue
        #print('idx',idx,'n_al',n["Al"],"c_al",d["Al"],'consider_atoms_al cnat_al',consider_atoms_al,'cnat consider_number_of_atoms',consider_number_of_atoms)
        if cellshape == "?":
            print('frames[i].cell')
            print(frames[i].cell)
        #if cellshape == "Q":
        #    print('frames[i].positions')
        #    print(frames[i].positions)

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

            if verbose > 2: #be_very_verbose:
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

            elif ace.geopt == True:
                ace.geopt = False
                ace.pot_to_ase_lmp_cmd()
                ene_pot_ase[idx] = ace.ene(atoms_tmp)
                #print('b',atoms_tmp.get_positions()[:3])
                #print('ene',ene_pot_ase[idx])

                ace.geopt = True
                ace.pot_to_ase_lmp_cmd()
                ene_pot_ase_geop[idx] = ace.ene(atoms_tmp)
                #print('a',atoms_tmp.get_positions()[:3])
                #print('ene',ene_pot_ase_geop[idx])
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

        ### write runner output
        if write_runner:
            nat = frames[i].get_number_of_atoms()
            listrepeat = [ [2,3,4], [4,10,3], [11,40,2] ]
            listrepeat = [ [2,4,2] ]          # memorize_symfunc_results on
            listrepeat = [ [2,4,3], [5,10,2]] # memorize_symfunc_results off
            listrepeat = [ [2,3,4], [4,12,3], [13,36,2]] # memorize_symfunc_results off, minatoms 44
            listrepeat = [ [2,3,4], [4,12,3], [13,70,2]] # memorize_symfunc_results off, minatoms 108
            #listrepeat = [ [2,2,5], [3,4,4], [5,12,3], [13,70,2]] # memorize_symfunc_results off, minatoms 108
            repeat = 0
            if nat == 1: repeat = 5         # *125  at = 125
            for rp in listrepeat:
                if rp[0] <= nat <= rp[1]: repeat = rp[2]
            #if 2 <= nat <= 3: repeat = 4    # *64   at = 124 - 256
            #if 4 <= nat <= 10: repeat = 3    # *27   at = 135 - 270
            #if 11 <= nat <= 40: repeat = 2   # *8    at = 99 -
            #if nat > 40: reppeat = False
            if repeat > 0:
                atoms_sc,forces_sc = my.ase_repeat_structure(frames[i],repeat)
            #print('ene',ene_DFT[idx],ace.units)
                ene_sc_ev = my.convert_energy(ene_DFT[idx]*repeat**3.,ace.units,"ev",frames[i])
            else:
                atoms_sc = copy.deepcopy(frames[i])
                forces_sc = atoms_sc.get_forces()
                ene_sc_ev = my.convert_energy(ene_DFT[idx],ace.units,"ev",frames[i])
            ase_write("out_runner_repeated.runner",atoms_sc,format='runner',append=True,setforces_ase_units=forces_sc,setenergy_eV=ene_sc_ev)

            # statistics
            nat_new = atoms_sc.get_number_of_atoms()
            if nat_new > max_at:
                max_at = nat_new
                max_at_orig = nat
                max_at_id = idx
            if nat_new <= min_at:
                min_at = nat_new
                min_at_orig = nat
                min_at_id = idx
                print(my.printred("min_at "+str(min_at)+" min_at_orig "+str(min_at_orig)+' (min id '+str(min_at_id)))

                # id: 4, 27,161, 188, 243,267,326,362,414,441,468,534,537 (atoms 58)
            print('nat',nat,'(repeat '+str(repeat)+') -> nat',nat*(3**repeat),"min_at",min_at,"min_at_orig",min_at_orig,'(min id '+str(min_at_id)+") max_at",max_at,"max_at_orig",max_at_orig,'(max id'+str(max_at_id)+")")

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
        print('out.runner written')

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
    print('#                         ('+ace_units+')                        ('+ace_units+')    ('+ace_units+')                                  ('+ace_units+')')
    print('#                         (DFT-ref)                                            ene_wo_atomic  forces   (if geopt)  Vol per')
    print('#   i   idx /    from       diff  [atms   Si   Mg   Al]   ene_DFT     ene_pot                 DFTmax    E-E_geopt   atom')
    print('--------------------------------------------------------------------------------------------------------------------------------------------------------------')
    return

def print_compare_ene_vs_DFT(text,pot_ene,DFT_ene="--",eos=False,f300=False,check=""):
    units="eV"
    conversion = 1
    if False:
        ev_to_kjpmol = 96.48534
        units = "kJ/mol"
        conversion = ev_to_kjpmol
    try:
        f300 = "%.3f" % np.round(f300,3)
    except:
        f300 = "----"

    try:
        eos = np.round(eos,2)[1:][:2]  # roundding to 3 is already beyond what can be well fitted
    except:
        eos = "----"

    ### DFT_ene_out / diff
    if DFT_ene == 0 or type(DFT_ene) == str:
        DFT_ene_out = "----"
        diff        = "----"
    else:
        DFT_ene_out = "%.4f" % round(DFT_ene*conversion,4)
        diff        = "%.2f" % round(np.abs(np.abs(pot_ene/DFT_ene)-1)*conversion,2)

    pot_ene_out = "%.4f" % round(pot_ene*conversion,4)
    #print('p',pot_ene_out,conversion,type(pot_ene_out),type(conversion))
    #print('p',pot_ene_out,conversion,"-->",pot_ene_out*conversion)
    print(text.ljust(45)+":",
            str(pot_ene_out).ljust(9),
            units.ljust(3)+" DFT:",
            str(DFT_ene_out).ljust(9),
            "diff:",
            str(diff).ljust(5),
            "eos:",
            str(eos).ljust(15),
            "df_harm:",
            f300,
            check)
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

    #p0rint_compare_ene_vs_DFT('vacancy_formation (unrelaxed) NN: ',e_ss_va,
    print_compare_ene_vs_DFT("vacancy_formation (unrelaxed, oDFT)",e_ss_va,0.69)
    print_compare_ene_vs_DFT("si-si-vac-complex (unrelaxed, oDFT)",e_si_si_vac_complex,0.074)
    print()
    return

def show_energy_diff_DFT_NN(struct,eDFT,e,text,units="eV"):
    ''' show_energy_diff_DFT_NN(struct_pure_al,ace.eDFT_pure_al,ace.e_pure_al,"pure al",units="eV") '''
    nat = struct.get_number_of_atoms()
    #eDFT = my.ase_enepot(struct_pure_al  ,units="eV")
    #e    = ace.ene(struct_pure_al)
    print(text.ljust(20," "),str(eDFT).ljust(15),"eV",str(e).ljust(15),"eV","diff",str(eDFT-e).ljust(17),"eV/cell","diff mev_pa",str((eDFT-e)/nat*1000).ljust(15),"meV_pa")
    return


def get_dilute_formation_energy(text="dilute formation energy supercell",sc="all",nsi=1,nmg=0,nvac=0,e_si_diamond_pa=0,ace=False,t2="",verbose=False):
    vpa = ace.al_fcc_vol_pa
    if sc == "all":
        sc_check = range(2,6)
    else:
        sc_check = range(sc,sc+1)
    if nsi == 1:  solute_element = "Si"
    if nmg == 1:  solute_element = "Mg"
    if nvac == 1: solute_element = "Vac"

    for i in sc_check:
        sc = i
        nat = 4*(i**3)

        ############## get the dilute frame Al+{Si,Mg} (T=0K ene_al_xx)
        frame_path = ace.savefolder+"frame_solute_Al"+str(nat-1)+solute_element+"1.runner"
        if os.path.isfile(frame_path):
            if ace.verbose:
                print('read frame path:',frame_path)
            frame_al_xx = ase_read(frame_path,format="runner")
            ene_al_xx = ace.ene(frame_al_xx) # this is already relaxed
        else:
            frame_al_xx = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=sc,nsi=nsi,nmg=nmg,nvac=nvac,a0=(vpa*4.)**(1./3.),cubic=True)
            ene_al_xx = ace.ene(frame_al_xx,atomrelax=True) # this takes a bit of time
            ase_write(frame_path,frame_al_xx,format="runner")

        ############## get the bulk al frame (T=0K ene_bulk)
        frame_path = ace.savefolder+"frame_bulk_Al"+str(nat)+".runner"
        if os.path.isfile(frame_path):
            if ace.verbose:
                print('read frame bulk:',filename)
            frame_bulk = ase_read(frame_path,format="runner")
            ene_bulk = ace.ene(frame_bulk) # this is already relaxed
        else:
            frame_bulk  = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=sc,nsi=0  ,nmg=0  ,nvac=0   ,a0=(vpa*4.)**(1./3.),cubic=True)
            ene_bulk  = ace.ene(frame_bulk)



        ############## get the dilute_formation @T=0K && ace.eform_dilute_xx_
        ############## get the dilute_formation @T=0K && ace.eform_dilute_xx_
        ############## get the dilute_formation @T=0K && ace.eform_dilute_xx_
        dilute_formation      = ene_al_xx  - (sc**3.*4. -1.) *ene_bulk/(sc**3.*4.)
        ace.eform_dilute_al_                     = ene_bulk/(sc**3.*4.)
        if solute_element == "Si":
            ace.eform_dilute_si_              = dilute_formation
        if solute_element == "Mg":
            ace.eform_dilute_mg_              = dilute_formation


        ############## get the temperature dependent stuff
        ############## get the temperature dependent stuff
        ############## get the temperature dependent stuff
        ############## get the temperature dependent stuff
        free_ene_path_al = ace.savefolder+"free_ene_formation_dilute_al_"
        free_ene_path_si = ace.savefolder+"free_ene_formation_dilute_si_"
        free_ene_path_mg = ace.savefolder+"free_ene_formation_dilute_mg_"
        sc_str = str(sc)+"x"+str(sc)+"x"+str(sc)+"sc"

        ### ace.free_ene_formation_dilute_al_
        ### ace.free_ene_formation_dilute_al_
        if os.path.isfile(free_ene_path_al) and os.path.isfile(free_ene_path_si) and os.path.isfile(free_ene_path_mg):
            print('loading dilute cells ...')
            ace.free_ene_formation_dilute_al_ = np.loadtxt(free_ene_path_al)
            ace.free_ene_formation_dilute_si_ = np.loadtxt(free_ene_path_si)
            ace.free_ene_formation_dilute_mg_ = np.loadtxt(free_ene_path_mg)
            ace.free_ene_formation_dilute_al_noshift = np.loadtxt(free_ene_path_al+"noshift")
            ace.free_ene_formation_dilute_si_noshift = np.loadtxt(free_ene_path_si+"noshift")
            ace.free_ene_formation_dilute_mg_noshift = np.loadtxt(free_ene_path_mg+"noshift")
            return

        ### in case the stuff can not be loaded...
        ### in case the stuff can not be loaded...
        ### in case the stuff can not be loaded...
        print('could not be loaded... '+solute_element)

        ### bulk free energy
        ### bulk free energy
        filename = ace.savefolder+"/h_"+sc_str+"_pure_al"
        if ace.verbose:
            print('tryread harmonic path dilute frame:',filename)
        free_ene_al_bulk = ace.get_fh(frame_bulk,debug=False,try_readfile=filename)
        ace.free_ene_formation_dilute_al_        = free_ene_al_bulk.ene_atom_only_ev_T0shifted  # save this
        ace.free_ene_formation_dilute_al_noshift = free_ene_al_bulk.ene_atom_only_ev # currently not used
        print('saving bulk al')
        np.savetxt(free_ene_path_al,ace.free_ene_formation_dilute_al_)
        np.savetxt(free_ene_path_al+"noshift",ace.free_ene_formation_dilute_al_noshift)


        ### dilute cell free energy
        ### dilute cell free energy
        filename = ace.savefolder+"/h_"+sc_str+"_al"+str(nat-1)+solute_element+"1"
        if ace.verbose:
            print('read fram path free_ene:',filename)
        print('could not be loaded ...')
        free_ene_al_xx = ace.get_fh(frame_al_xx,debug=False,try_readfile=filename)

        free_dilute_formation = free_ene_al_xx.ene_cell_only_ev_T0shifted - (sc**3.*4. -1.) *free_ene_al_bulk.ene_atom_only_ev_T0shifted
        free_dilute_formation_noshift = free_ene_al_xx.ene_cell_only_ev - (sc**3.*4. -1.) *free_ene_al_bulk.ene_atom_only_ev


        ###################################### GET DILUTE FORMATION
        ###################################### GET DILUTE FORMATION
        ###################################### GET DILUTE FORMATION
        if solute_element == "Si":
            ace.free_ene_formation_dilute_si_ = free_dilute_formation
            ace.free_ene_formation_dilute_si_noshift = free_dilute_formation_noshift
            print('saving dilute si')
            np.savetxt(free_ene_path_si,ace.free_ene_formation_dilute_si_)
            np.savetxt(free_ene_path_si+"noshift",ace.free_ene_formation_dilute_si_noshift)

        if solute_element == "Mg":
            ace.free_ene_formation_dilute_mg_ = free_dilute_formation
            ace.free_ene_formation_dilute_mg_noshift = free_dilute_formation_noshift
            np.savetxt(free_ene_path_mg,ace.free_ene_formation_dilute_mg_)
            np.savetxt(free_ene_path_mg+"noshift",ace.free_ene_formation_dilute_si_noshift)

        print('dilute_formation T0K',text,"sc:",sc,str(round(dilute_formation,3)).ljust(8),"eV; formation energy:", round(dilute_formation - e_si_diamond_pa,3),"eV",t2)
        print('dilute_formation T0K',text,"sc:",sc,str(round(free_dilute_formation[0],3)).ljust(8),"eV; formation energy:", round(free_dilute_formation[0] - e_si_diamond_pa,3),"eV",t2)



    ###################################
    # this seems to be in quality as well converging with supercell size as when volume of the defect is not relaxed
    ###################################
    #for i in range(2,7):
    #    sc = i
    #    frame_al_xx = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=sc,nsi=nsi,nmg=nmg,nvac=nvac,a0=(vpa*4.)**(1./3.),cubic=True)
    #    vp = ace.get_murn(frame_al_xx,verbose=False,return_minimum_volume_frame=True)
    #    frame_bulk  = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=sc,nsi=0  ,nmg=0  ,nvac=0   ,a0=(vpa*4.)**(1./3.),cubic=True)
    #    ene_al_xx = ace.ene(frame_al_xx,atomrelax=True)
    #    ene_bulk  = ace.ene(frame_bulk)
    #    dilute_formation = ene_al_xx - (sc**3.*4. -1.) *ene_bulk/(sc**3.*4.)
    #    print('dilute formation energy supercell (vol def = relaxed)',sc,dilute_formation, dilute_formation - e_si_diamond_pa)
    #return frame_bulk, frame_al_xx
    return

def get_al_fcc_equilibrium(ace):
    filename_frame = ace.savefolder+"frame_al_fcc.runner"
    if os.path.isfile(filename_frame):
        frame_al = ase_read(filename_frame)
    else:
        frame_al = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=1,nsi=0,nmg=0,nvac=0,a0=4.045,cubic=True,create_fake_vacancy=False,whichcell="fcc")
        ace.ase_relax_cellshape_and_volume_only(frame_al,verbose=False)
    ace.al_fcc = frame_al
    ace.al_fcc_ene_pa = my.ase_vpa(ace.al_fcc)
    ace.al_fcc_vol_pa = frame_al.get_volume()/frame_al.get_number_of_atoms()
    if ace.verbose:
        print("NN Al vpa @T=0K",ace.al_fcc_vol_pa,"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    #print('e_ and f_ al should all be done consistently from the NN!!!')
    if not os.path.isfile(filename_frame):
        ase_write(ace.savefolder+"frame_al_fcc.runner",frame_al,format='runner')
    ace.get_elastic_external(atomsin=ace.al_fcc,verbose=False,text="Al_fcc bulk 4at",get_all_constants="C44")

    #frame_al = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=3,nsi=0,nmg=0,nvac=0,a0=4.045,cubic=True,create_fake_vacancy=False,whichcell="fcc")
    #print('nat',frame_al.get_number_of_atoms())
    #ace.ase_relax_cellshape_and_volume_only(frame_al,verbose=False)
    #ace.get_elastic_external(atomsin=frame_al,verbose=False,text="Al_fcc bulk 108at")
    return

def get_mg_hcp_equilibrium(ace):
    filename_frame = ace.savefolder+"frame_mg_hcp.runner"
    if os.path.isfile(filename_frame):
        frame_mg = ase_read(filename_frame)
    else:
        frame_mg = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=1,nsi=0,nmg=0,nvac=0,a0=0,cubic=False,create_fake_vacancy=False,whichcell="hcp")
        ace.ase_relax_cellshape_and_volume_only(frame_mg,verbose=False)
    ace.mg_hcp = frame_mg
    ace.mg_hcp_ene_pa = ace.ene(frame_mg)/frame_mg.get_number_of_atoms()
    ace.mg_hcp_vol_pa = frame_mg.get_volume()/frame_mg.get_number_of_atoms()
    if not os.path.isfile(filename_frame):
        ase_write(ace.savefolder+"frame_mg_hcp.runner",frame_mg,format='runner')
    #ace.get_elastic_external(atomsin=ace.mg_hcp,verbose=False,text="Mg_hcp bulk")
    return

def get_si_dc_equilibrium(ace):
    filename_frame = ace.savefolder+"frame_si_dc.runner"
    if os.path.isfile(filename_frame):
        frame_si = ase_read(filename_frame)
    else:
        frame_si = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=1,nsi=0,nmg=0,nvac=0,a0=0,cubic=False,create_fake_vacancy=False,whichcell="dc")
        ace.ase_relax_cellshape_and_volume_only(frame_si,verbose=False)
    ace.si_dc = frame_si
    ace.si_dc_ene_pa = ace.ene(frame_si)/frame_si.get_number_of_atoms()
    ace.si_dc_vol_pa = frame_si.get_volume()/frame_si.get_number_of_atoms()
    ase_write(ace.savefolder+"frame_si_dc.runner",frame_si,format='runner')
    #ace.get_elastic_external(atomsin=ace.si_dc,verbose=False,text="Si_dc bulk")

    #frame_si = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=4,nsi=0,nmg=0,nvac=0,a0=0,cubic=False,create_fake_vacancy=False,whichcell="dc")
    #print('at',frame_si.get_number_of_atoms(),my.ase_vpa(frame_si))
    #ace.ase_relax_cellshape_and_volume_only(frame_si,verbose=False)
    #print('at',frame_si.get_number_of_atoms(),my.ase_vpa(frame_si))
    #filename = ace.savefolder+"/h_4x4x4sc_pure_si_dc"
    #free_ene_pure_si = ace.get_fh(frame_si,debug=False,try_readfile=filename,return_units="ev_cell")
    #free_ene_pure_si = ace.get_fh(frame_si,debug=False,try_readfile=filename,return_units="mev_pa",shiftT0_ev_pa=0)
    #print('fe',free_ene_pure_si[:5])
    #print('fe',free_ene_pure_si[-5:])
    #print()
    #print('fe',free_ene_pure_si[:5])
    #print('fe',free_ene_pure_si[-5:])
    #sys.exit()
    return

def get_basic_NN_energies_ace(ace):
    print('######## get_basic_NN_energies_ace #############')
    scripts = my.scripts()
    tests = scripts+'/tests/'
    ace.savefolder = tests+'/Al-Mg-Si/save_'+ace.pot.pot+"/"
    if not os.path.isdir(ace.savefolder):
        my.mkdir(ace.savefolder)

    #ase_write('pos_si.lmp',frame_si,format='lammps-runner')
    get_al_fcc_equilibrium(ace)
    get_mg_hcp_equilibrium(ace)
    get_si_dc_equilibrium(ace)
    ace.free_ene_formation_dilute_mg_ = False
    ace.free_ene_formation_dilute_si_ = False
    #ene_pot_lmp = my.lammps_ext_calc(ace.al_fcc,ace)
    #ase_write('pos_al.lmp',ace.al_fcc,format='lammps-runner')
    #get_vpa(
    if ace.verbose:
        print("NN 1 Al vpa @T=0K",my.ase_vpa(ace.al_fcc),"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    #ace.check_frame('bbb',ace.al_fcc)
    #ace.ase_relax_cellshape_and_volume_only(ace.al_fcc,verbose=False)
    #ace.check_frame('aaa',ace.al_fcc)
    if ace.verbose:
        print("NN 2 Al vpa @T=0K",my.ase_vpa(ace.al_fcc),"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    vinet = ace.get_murn(ace.al_fcc,verbose=False,return_minimum_volume_frame = False, atomrelax=False,write_energies=False)
    if ace.verbose:
        print("NN 3 Al vpa @T=0K",my.ase_vpa(ace.al_fcc),"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    #ace.ase_relax_cellshape_and_volume_only(ace.al_fcc,verbose=False)
    #ace.check_frame('ccc',ace.al_fcc)
    #sys.exit()
    if ace.verbose:
        print('vinet',vinet)
        print("NN 4 Al vpa @T=0K",my.ase_vpa(ace.al_fcc)) #,"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    #sys.exit('ace fcc al si dc')
        print('before too lon 12')
    get_dilute_formation_energy(text="NN dilute formation energy Si ",sc=4,nsi=1,nmg=0,e_si_diamond_pa=ace.si_dc_ene_pa,ace=ace,t2="Kobayashi 0.375 eV")
    #sys.exit('too lon 12')
    get_dilute_formation_energy(text="NN dilute formation energy Mg ",sc=4,nsi=0,nmg=1,e_si_diamond_pa=ace.mg_hcp_ene_pa,ace=ace,t2="Kobayashi 0.090 eV")
    get_dilute_formation_energy(text="NN dilute formation energy Vac",sc=4,nsi=0,nmg=0,nvac=1,e_si_diamond_pa=0.,ace=ace,t2="Kobayashi 0.654 eV")
    print()
    return

def get_dilute_si_mg_f(ace):
    print("######## get_dilute_si_mg_f #############")

    scripts = my.scripts()
    tests = scripts+'/tests/'
    bprime = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/dilute_structures_and_beta_prime.input.data'


    frames = ase_read(bprime,index=":",format="runner")

    struct_bprime    = frames[0]
    struct_pure_al   = frames[3]
    struct_dilute_mg = frames[1]
    struct_dilute_si = frames[2]

    # @ 0K (eV/cell)
    ace.eDFT_bprim      = my.ase_enepot(struct_bprime   ,units=ace.units)
    ace.eDFT_pure_al    = my.ase_enepot(struct_pure_al  ,units=ace.units) # 108 atoms
    ace.eDFT_dilute_mg  = my.ase_enepot(struct_dilute_mg,units=ace.units) # 108 atoms
    ace.eDFT_dilute_si  = my.ase_enepot(struct_dilute_si,units=ace.units) # 108 atoms
    ace.e_bprime        = ace.ene(struct_bprime)
    ace.e_pure_al       = ace.ene(struct_pure_al)   # 108 atoms
    ace.e_dilute_mg     = ace.ene(struct_dilute_mg) # 108 atoms
    ace.e_dilute_si     = ace.ene(struct_dilute_si) # 108 atoms


    print("         ","DFT(eV)          NN(eV)                 eV all atoms")
    show_energy_diff_DFT_NN(struct_pure_al,ace.eDFT_pure_al,ace.e_pure_al,"pure al",units="eV")
    show_energy_diff_DFT_NN(struct_dilute_mg,ace.eDFT_dilute_mg,ace.e_dilute_mg,"dilute mg",units="eV")
    show_energy_diff_DFT_NN(struct_dilute_si,ace.eDFT_dilute_si,ace.e_dilute_si,"dilute si",units="eV")
    show_energy_diff_DFT_NN(struct_bprime,ace.eDFT_bprim,ace.e_bprime,"bprime",units="eV")
    print()

    # (eV/defect)
    ace.eform_dilute_al     =  ace.e_pure_al/108.
    ace.eform_dilute_si     =  ace.e_dilute_si    - 107*ace.e_pure_al/108       # eV per defect
    ace.eform_dilute_mg     =  ace.e_dilute_mg    - 107*ace.e_pure_al/108       # eV per defect
    ace.fDFT_dilute_al      =  ace.eDFT_pure_al/108    # eV per defect
    ace.fDFT_dilute_si      =  ace.eDFT_dilute_si - 107*ace.eDFT_pure_al/108    # eV per defect
    ace.fDFT_dilute_mg      =  ace.eDFT_dilute_mg - 107*ace.eDFT_pure_al/108    # eV per defect

    print_compare_ene_vs_DFT('formation dilute si (1Si in bulk Al) 3x3x3',ace.eform_dilute_si,ace.fDFT_dilute_si, '-',"-")
    print_compare_ene_vs_DFT('formation dilute mg (1Mg in bulk Al) 3x3x3',ace.eform_dilute_mg,ace.fDFT_dilute_mg, '-',"-")
    print()
    print_compare_ene_vs_DFT('formation dilute si (1Si in bulk Al) 3x3x3',ace.eform_dilute_si-ace.si_dc_ene_pa ,"-", '-',"-")
    print_compare_ene_vs_DFT('formation dilute mg (1Mg in bulk Al) 3x3x3',ace.eform_dilute_mg-ace.mg_hcp_ene_pa,"-", '-',"-")
    print()


    #sys.exit('dda')

    # @ 300K ( for 1 atom, therefore *108 (to get per cell), in meV, therefore /1000)
    filename = ace.savefolder+"/free_ene_pure_al_3x3x3sc_fixedstruct_108at"
    free_ene_al_bulk = ace.get_fh(struct_pure_al,debug=False,try_readfile=filename) #,T0shift_ev_atom = ene_bulk/nat) #return_units="ev_cell")
    ace.free_ene_formation_dilute_al     = free_ene_al_bulk.ene_atom_only_ev_T0shifted

    ## dilute mg
    filename = ace.savefolder+"/free_ene_dilute_mg_3x3x3sc_fixedstruct_108at"
    free_ene_dilute_mg = ace.get_fh(struct_dilute_mg,debug=False,try_readfile=filename) #,return_units="ev_cell")


    ## dilute si
    filename = ace.savefolder+"/free_ene_dilute_si_3x3x3sc_fixedstruct_108at"
    free_ene_dilute_si = ace.get_fh(struct_dilute_si,debug=False,try_readfile=filename) #,return_units="ev_cell")


    ### get free_ene_formation_dilute_xx
    ### get free_ene_formation_dilute_xx
    #ace.free_ene_al = ace.free_ene_pure_al/108.
    #ace.free_ene_formation_dilute_si =  ace.free_ene_dilute_si - 107.*ace.free_ene_pure_al/108.
    #ace.free_ene_formation_dilute_mg =  ace.free_ene_dilute_mg - 107.*ace.free_ene_pure_al/108.
    ace.free_ene_formation_dilute_al = free_ene_al_bulk.ene_atom_only_ev_T0shifted
    ace.free_ene_formation_dilute_si = free_ene_dilute_si.ene_cell_only_ev_T0shifted - 107.*free_ene_al_bulk.ene_atom_only_ev_T0shifted
    ace.free_ene_formation_dilute_mg = free_ene_dilute_mg.ene_cell_only_ev_T0shifted - 107.*free_ene_al_bulk.ene_atom_only_ev_T0shifted
    return


def get_formation_energy(ace,frame,text,atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=False,try_harmonic_readfile=False,debug=False):
    ''' Bill sais that T=443 is the relevant temperature '''
    d = my.ase_get_chemical_symbols_to_number_of_species(frame)
    conz1 = d["Mg"]
    conz2 = (d["Mg"]+d["Si"])
    conz = np.float(d["Mg"])/np.float((d["Mg"]+d["Si"]))
    #print('mg',d["Mg"],'si',d["Si"])
    nat = frame.get_number_of_atoms()

    heat_precip_T0K_DFT = "-"
    eDFT = ""
    if DFT_ene != False:
        eDFT   = my.ase_enepot(frame  ,units=ace.units)
        e      = ace.ene(frame)
        ediff_ev  = e - eDFT
        ediff_mev_pa = ediff_ev*1000./nat
        print('e - eDFT:',round(ediff_mev_pa,3),"meV/pa")
        heat_precip_T0K_DFT = (eDFT - d["Mg"]*ace.fDFT_dilute_mg - d["Si"]*ace.fDFT_dilute_si - d["Al"]*ace.fDFT_dilute_al)/nat
        if ace.verbose:
            show_energy_diff_DFT_NN(frame,eDFT,e,text,units="eV")
            print("DFT energy precipitate (eV)",eDFT)
            print("DFT energy eform_dilute_mg (eV)",ace.fDFT_dilute_mg,"times",d["Mg"])
            print("DFT energy eform_dilute_si (eV)",ace.fDFT_dilute_si,"times",d["Si"])
            print("divide everything by       ",nat,"to get to the formation energy of",heat_precip_T0K_DFT)
        if "@DFT" in text:
            file = ace.savefolder+"summary_formations_DFT_T0.dat"
            if os.path.isfile(file) and ace.written_summary[0] == False:
                os.remove(file)
                ace.written_summary[0] = True
            f=open(file, "a+")
            f.write(str(conz)+"   "+str(heat_precip_T0K_DFT)+" "+text.replace(" ", "_")+"\n")
            f.close()

    # @ T=0K
    if atomrelax: ace.ase_relax_atomic_positions_only(frame)
    if volumerelax: ace.ase_relax_cellshape_and_volume_only(frame)
    if atomrelax: ace.ase_relax_atomic_positions_only(frame)
    if volumerelax: ace.ase_relax_cellshape_and_volume_only(frame)
    check = ace.check_frame('',frame=frame,verbose=False)
    #print('11 atomrelax',atomrelax,'volumerelax',volumerelax,"check",check)
    vinet = ace.get_murn(frame,verbose=False,return_minimum_volume_frame = volumerelax, atomrelax=atomrelax,write_energies=False)
    e = ace.ene(frame) #,atomrelax=atomrelax,cellrelax=cellrelax)
    #print('d mg:',d["Mg"],'d si:',d["Si"],'d al:',d["Al"])
    #heat_precip_T0K         = (e - d["Mg"]*ace.eform_dilute_mg - d["Si"]*ace.eform_dilute_si - d["Al"]*ace.eform_dilute_al)/nat
    if ace.verbose:
        print('e                   ',e)
        print('ace.eform_dilute_mg_',ace.eform_dilute_mg_)
        print('ace.eform_dilute_si_',ace.eform_dilute_si_)
    heat_precip_T0K         = (e - d["Mg"]*ace.eform_dilute_mg_ - d["Si"]*ace.eform_dilute_si_ - d["Al"]*ace.eform_dilute_al_)/nat

    check = ace.check_frame('',frame=frame,verbose=False)
    #print('22 atomrelax',atomrelax,'volumerelax',volumerelax,"check",check)
    print_compare_ene_vs_DFT(text+" @0K",heat_precip_T0K,heat_precip_T0K_DFT,vinet,"-",check=check)

    if ace.verbose:
        print("NN energy precipitate (eV)",e)
        print("NN energy eform_dilute_mg (eV)",ace.eform_dilute_mg,"times",d["Mg"])
        print("NN energy eform_dilute_si (eV)",ace.eform_dilute_si,"times",d["Si"])
        print("divide everything by       ",nat,"to get to the formation energy of",heat_precip_T0K_DFT)

    if True:
        # @ ace.atTemp K
        #print('vol',my.ase_vpa(frame))
        if ace.verbose:
            print('try_harmonic_readfile:',try_harmonic_readfile)
        free_ene = ace.get_fh(frame,try_readfile=try_harmonic_readfile,atomrelax=atomrelax,debug=debug) #,return_units="ev_cell")
        if type(free_ene) == bool:
            print(text+" @"+str(ace.atTemp)+"K            : --> NEGATIVE EIGENVALUES, structure not in ground state")
            return

        if ace.verbose:
            #print('e   free.ene_cell   ',free_ene.ene_cell)
            #print('e   free.ene_atom   ',free_ene.ene_atom)

            print('e   free            ',free_ene.ene_cell_only_ev_T0shifted)
            print('aceeform_dilute_mg_',ace.free_ene_formation_dilute_mg_)
            print('aceeform_dilute_si_',ace.free_ene_formation_dilute_si_)

        heat_precip_T = (free_ene.ene_cell_only_ev_T0shifted - d["Mg"]*ace.free_ene_formation_dilute_mg_ - d["Si"]*ace.free_ene_formation_dilute_si_ - d["Al"]*ace.free_ene_formation_dilute_al_)/nat

        if ace.verbose:
            print('heat precip',heat_precip_T)
        #free_dilute_formation = free_ene_al_xx.ene_cell_only_ev_T0shifted - (sc**3.*4. -1.) *free_ene_al_bulk.ene_atom_only_ev_T0shifted
        #print('heat',heat_precip_T[0])
        #sys.exit('kka')

        #np.savetxt("free_ene_formation_dilute_mg.dat",ace.free_ene_formation_dilute_mg_)
        #np.savetxt("free_ene_formation_dilute_mg.dat",ace.free_ene_formation_dilute_mg_)
        #np.savetxt("free_ene_formation_dilute_si.dat",ace.free_ene_formation_dilute_si_)
        #np.savetxt("free_ene_formation_dilute_mg_only.dat",ace.free_ene_formation_dilute_mg_noshift)
        #np.savetxt("free_ene_formation_dilute_al_only.dat",ace.free_ene_formation_dilute_al_noshift)
        #np.savetxt("free_ene_formation_dilute_si_only.dat",ace.free_ene_formation_dilute_si_noshift)
        #np.savetxt("free_ene_al.dat",ace.free_ene_formation_dilute_al_)


        print_compare_ene_vs_DFT(text+" @0K zpv",heat_precip_T[0],"","",heat_precip_T[0])
        print_compare_ene_vs_DFT(text+" @"+str(ace.atTemp)+"K",heat_precip_T[ace.atTemp-1],"","",heat_precip_T[0])

        if "@NN" in text:
            #if ace.verbose:
            #    print("free_ene_cell_"+text[:6]+".dat")
            #np.savetxt("free_ene_cell_"+text[:6]+".dat",free_ene.ene_cell_only_ev_T0shifted)
            #np.savetxt("free_ene_atom_"+text[:6]+".dat",free_ene.ene_atom_only_ev_T0shifted)
            #np.savetxt("free_ene_atom_"+text[:6]+"_only.dat",free_ene.ene_atom_only_ev)
            file = ace.savefolder+"summary_formations_NN_T0.dat"
            if os.path.isfile(file) and ace.written_summary[1] == False:
                os.remove(file)
                ace.written_summary[1] = True
            f=open(file, "a+")
            f.write(str(conz)+"   "+str(heat_precip_T0K)+" "+text.replace(" ", "_")+"\n")
            f.close()


            file = ace.savefolder+"summary_formations_NN_T"+str(ace.atTemp)+".dat"
            if os.path.isfile(file) and ace.written_summary[2] == False:
                os.remove(file)
                ace.written_summary[2] = True
            f=open(file, "a+")
            print("write",file,conz,text)
            f.write(str(conz)+"   "+str(heat_precip_T[ace.atTemp-1])+" "+text.replace(" ", "_")+"\n")
            f.close()

    return



def test_betaprime_mg9si5_find_global_min(ace,eform_dilute_si, eform_dilute_mg, f_dilute_si_300, f_dilute_mg_300, find_global_minimum=True):
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

        e_precip    = ace.ene(struct_mg9si5) # 28 atoms, 10Si, 18Mg


        dilute_si_f_DFT =  eDFT_dilute_si - 107*eDFT_pure_al/108
        f_dilute_mg_DFT =  eDFT_dilute_mg - 107*eDFT_pure_al/108
        heat_precip_DFT = (eDFT_precip    - 18*f_dilute_mg_DFT - 10*dilute_si_f_DFT)/28


        if find_global_minimum:
            print("#############################")
            print("FIND GLOBAL MINIMUM")
            print("#############################")


            ### unrelaxed
            vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
            f300 = ace.get_fh(struct_mg9si5,debug=False)
            heat_precip = (e_precip - 18*eform_dilute_mg - 10*eform_dilute_si)/28
            print_compare_ene_vs_DFT("beta prime Mg9Si5 (unrelaxed)",heat_precip,0,vinet_mg9si5,f300)


            ### atomrelax
            e_dilute_mg_r    = ace.ene(struct_dilute_mg,atomrelax=True) # 108 atoms
            e_dilute_si_r    = ace.ene(struct_dilute_si,atomrelax=True) # 108 atoms
            e                = ace.ene(struct_mg9si5   ,atomrelax=True)

            heat_precip = (e - 18*eform_dilute_mg - 10*eform_dilute_si)/28

            print_compare_ene_vs_DFT("beta prime Mg9Si5 (atomrelax)",heat_precip,0,vinet_mg9si5,f300)
            ase_write("mg9si5_stable_phonons_v2_atomrelax.runner",struct_mg9si5,format='runner')
            print('stress',struct_mg9si5.get_stress())
            print()

            ### try to find global minimum
            print('------- 3')
            e = ace.ene(struct_mg9si5   ,atomrelax=True,cellrelax = True)
            vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
            f300 = ace.get_fh(struct_mg9si5,debug=False)
            heat_precip = (e - 18*eform_dilute_mg - 10*eform_dilute_si)/28
            print_compare_ene_vs_DFT("beta prime Mg9Si5 (atomrelax+cellrelax)",heat_precip,0,vinet_mg9si5,f300)
            ase_write("mg9si5_stable_phonons_v2_atomrelax_cellrelax.runner",struct_mg9si5,format='runner')
            print('stress',struct_mg9si5.get_stress())
            print()

            print('------- 4')
            e = ace.ene(struct_mg9si5   ,atomrelax=True,minimizer='mh')
            print('v4',my.ase_vpa(struct_mg9si5),struct_mg9si5.get_potential_energy(),my.ase_mepa(struct_mg9si5))
            vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
            f300 = ace.get_fh(struct_mg9si5,debug=False)
            heat_precip = (e - 18*eform_dilute_mg - 10*eform_dilute_si)/28.
            print_compare_ene_vs_DFT("beta prime Mg9Si5 (minima hopping)",heat_precip,0,vinet_mg9si5,f300)
            print('stress',struct_mg9si5.get_stress())
            print()


            print('------- 5')
            e = ace.ene(struct_mg9si5   ,atomrelax=True,cellrelax=True)
            print('v5',my.ase_vpa(struct_mg9si5),struct_mg9si5.get_potential_energy(),my.ase_mepa(struct_mg9si5))
            vinet_mg9si5 = ace.get_murn(struct_mg9si5) # 28 atoms, 10Si, 18Mg
            f300 = ace.get_fh(struct_mg9si5,debug=False)
            heat_precip = (e - 18*eform_dilute_mg - 10*eform_dilute_si)/28.
            print_compare_ene_vs_DFT("beta prime Mg9Si5 (atomrelax+cellrelax)",heat_precip,0,vinet_mg9si5,f300)
            ase_write("mg9si5_stable_phonons_v2_atomrelax_cellrelax_minimahopping.runner",struct_mg9si5,format='runner')
            print('stress',struct_mg9si5.get_stress())
            print()
            ase_write("mg9si5_stable_phonons_v2.runner",struct_mg9si5,format='runner')
        return

def test_Mg9Si5(ace):
    print("######## test_Mg9Si5 (high Mg conz) #############")
    print("######## test_Mg9Si5 (high Mg conz) #############")
    print("######## test_Mg9Si5 (high Mg conz) #############")

    #########################################
    # @DFT positions
    #########################################
    path = my.scripts()+'/tests/Al-Mg-Si/Mg9Si5_beta_prime/exported_from_aiida/aiida_exported_group_BetaPrime_vc-relaxed__only_relaxed.input.data'
    ### /home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/tests/Al-Mg-Si/Mg9Si5_beta_prime/exported_from_aiida/aiida_exported_group_BetaPrime_vc-relaxed__only_relaxed.input.data

    ### has current energy of -359.771677826279586 hartree (with ase units eV_to_Hartree = 0.03674932247495664) has 22.523 mev_pa diff
    ### has old     energy of -359.771702546166239 hartree (with old eV_to_Hartree = 0.036749325)               has 22.547 mev_pa diff
    ### qe original energy of -9789.88600597 eV (this is what aiida reports)
    ### qe original energy of
    if ace.verbose:
        print('path frame Mg9Si5:',path)
    frame = ase_read(path,format="runner")

    print_compare_ene_vs_DFT("Mg9Si5 (@DFT fully relaxed) @0K Vissers GGA",-0.335,DFT_ene="-",eos=False,f300=False)


    try_read = ace.savefolder+"h_Mg9Si5_at_DFT_relaxed"
    if ace.verbose:
        print('try_read hessematrix',try_read)
    get_formation_energy(ace,frame,"Mg9Si5 (@DFT fully relaxed)",atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=True,try_harmonic_readfile=try_read)

    #########################################
    # relaxed with NN
    #########################################
    relaxed_frame = ace.savefolder+"frame_Mg9Si5_at_NN_relaxed_mh.runner"
    if os.path.isfile(relaxed_frame):
        #print('reading relaxed Mg9Si5')
        frame = ase_read(relaxed_frame)
    else:
        print("################### relaxing with NN by minima hopping ######")
        e = ace.ene(frame   ,atomrelax=True)
        ase_write(relaxed_frame+".local_minimization.runner",frame)
        print('e',e)
        e = ace.ene(frame   ,atomrelax=True,minimizer='mh')
        print('e',e)
        ase_write(relaxed_frame,frame)

    try_read = ace.savefolder+"h_Mg9Si5_at_NN_relaxed"
    if ace.verbose:
        print('try_read free ene Mg9Si5 relaxed:',try_read)
    get_formation_energy(ace,frame,"Mg9Si5 (@NN  fully relaxed)" ,atomrelax=True,cellrelax=True ,volumerelax=True,try_harmonic_readfile=try_read,debug=False)
    #print(np.round(frame.get_positions(),2))
    #print()
    #print(np.round(frame.get_cell(),2))
    return

def test_Mg9Si5_pos(ace):
    path = my.scripts()+'/tests/Al-Mg-Si/Mg9Si5_beta_prime/BetaPrime_structures_relax.input.data.v3'
    frame = ase_read(path,'0',format="runner")
    get_formation_energy(ace,frame,"Mg9Si5 (poslaxed @DFT)",atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=True)
    return


def test_Mg2Si(ace):
    path = my.scripts()+'/tests/Al-Mg-Si/Mg2Si/POSCAR'
    frame = ase_read(path,format="vasp")
    try_read = ace.savefolder+"h_Mg2Si1_at_NN_relaxed"
    #print('try read',try_read)
    get_formation_energy(ace,frame,"Mg2Si (fully relaxed)",atomrelax=True,cellrelax=True,volumerelax=True,try_harmonic_readfile=try_read)
    print('mg2si stress:',frame.get_stress())
    #print('mg2si stress:',frame.get_isotropic_pressure(frame.get_stress()))
    #print('vinet',ace.get_murn(frame),my.ase_vpa(frame))
    print()
    return

def test_beta2_bulk(ace):
    print("######## test_Mg5Si6 (low Mg conz) #############")
    print("######## test_Mg5Si6 (low Mg conz) #############")
    print("######## test_Mg5Si6 (low Mg conz) #############")
    doit = [ "Mg5Si6", "Mg5Al2Si4", "Mg4Al3Si4" ]
    #doit = [ "Mg5Si6" ]
    for i in doit:
        #print()
        #path = my.scripts()+'/tests/Al-Mg-Si/Beta2-bulk/'+i+"/POSCAR'
        #frame = ase_read(path,format="vasp")
        #get_formation_energy(ace,frame,"Mg5Si6 (fully relaxed)",atomrelax=True,cellrelax=True,volumerelax=True)

        print()
        #path = my.scripts()+'/tests/Al-Mg-Si/Beta2-bulk/'+i+'/'
        #frame = ase_read(path+"runner.data",format="runner")

        path = my.scripts()+'/tests/Al-Mg-Si/Beta2-bulk/'+i+'/aiida_exported_group_NN_relaxed_'+i+"_n2p2_v2ag_calc__only_relaxed.input.data"
        frame = ase_read(path,format="runner")

        #get_formation_energy(ace,frame,i+" (unrelaxed    )",atomrelax=False,cellrelax=False,volumerelax=False)
        try_read = ace.savefolder+"h_"+i+"_at_DFT_relaxed"
        get_formation_energy(ace,frame,i+" (@DFT fully relaxed)",atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=True,try_harmonic_readfile=try_read)
        try_read = ace.savefolder+"h_"+i+"_at_NN_relaxed"
        get_formation_energy(ace,frame,i+" (@NN  fully relaxed)",atomrelax=True,cellrelax=True,volumerelax=True,DFT_ene=False,try_harmonic_readfile=try_read)
        ase_write(path+"NN_relaxed_"+i+"_"+ace.pot.pot+".runner",frame,format='runner')

    return


def load_diluete_pure_values():
    scripts = my.scripts()
    filename = scripts+'/tests/Al-Mg-Si/get_dilute_si_mg_f.'+pot+".dat"


def test_formation_energies(ace):
    ace.written_summary = [False,False,False]
    ace.atTemp = 443

    get_basic_NN_energies_ace(ace)
    get_dilute_si_mg_f(ace)

    file = ace.savefolder+"summary_formationsT0.dat"
    if os.path.isfile(file): os.remove(file)
    file = ace.savefolder+"summary_formationsT"+str(ace.atTemp)+".dat"
    if os.path.isfile(file): os.remove(file)
    print()
    print("########### test_formation_energies #########################")
    print()
    #test_si_si_vac(ace)
    #test_Mg2Si(ace)
    test_Mg9Si5(ace)

    ##test_Mg9Si5_pos(ace)
    test_beta2_bulk(ace)
    ##test_betaprime_mg9si5_find_global_min(ace,eform_dilute_si, eform_dilute_mg, f_dilute_si_300, f_dilute_mg_300)
    return




def test2_elastic(ace):
    #get_basic_NN_energies_ace(ace)
    #get_al_fcc_equilibrium(ace)
    ace.elastic_relax = True
    frame_al = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=1,nsi=0,nmg=0,nvac=0,a0=4.045,cubic=True,create_fake_vacancy=False,whichcell="fcc")
    ace.get_elastic_external(atomsin=frame_al,verbose=ace.verbose,text="Al_fcc bulk 4at",get_all_constants=True)
    print('ace.c44:',ace.c44,type(ace.c44))
    print('aa',ace.pot.potpath)
    np.savetxt(ace.pot.potpath+"/elastic.dat",np.array([np.float(ace.c44)]))
    #print('ace.elastic_constants:',ace.elastic_constants)
    #print('cc')
    sys.exit()

    path = my.scripts()+'/tests/Al-Mg-Si/SimpleAlDeformations/SimpleAlDeformations_scf.runner'
    print('path',path)
    frames = ase_read(path,":",format="runner")
    if type(frames) == list:
        structures_to_calc = len(frames)
    else:
        structures_to_calc = 1

    for idx,i in enumerate(range(structures_to_calc)):
        frame = frames[idx]
        ace.get_elastic_external(atomsin=frame,verbose=False,text="Al_fcc simple deformed 4at")

    return

def test3_do(ace):
    path = "/home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/tests//Al-Mg-Si/save_n2p2_v3ag/frame_solute_Al255Si1.runner"
    #print('path frame',path)
    frame = ase_read(path,format="runner")
    try_readfile = "/home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/tests//Al-Mg-Si/save_n2p2_v3ag//h_4x4x4sc_al255Si1"
    #print('fqh try_readfile',try_readfile)
    free_ene = ace.get_fh(frame,debug=False,try_readfile=try_readfile)
    #print(free_ene.ene_atom)
    #print()
    print('############################ Mg9Si5 ################3')
    path = my.scripts()+'/tests/Al-Mg-Si/Mg9Si5_beta_prime/exported_from_aiida/aiida_exported_group_BetaPrime_vc-relaxed__only_relaxed.input.data'
    print('path frame Mg9Si5',path)
    frame = ase_read(path,format="runner")

    try_readfile = "/home/glensk/Dropbox/Albert/scripts/dotfiles/scripts/tests//Al-Mg-Si/save_n2p2_v3ag/h_Mg9Si5_at_DFT_relaxed"
    print('fqh try_readfile',try_readfile)
    free_ene = ace.get_fh(frame,debug=False,try_readfile=try_readfile)
    #print(free_ene.ene_atom)
    print('done kk')

if __name__ == "__main__":
    get_energies()
