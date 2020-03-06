#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import timeit
start = timeit.timeit()
import sys,os,copy,argparse,shutil
from subprocess import check_output,call
#import click
import glob,time
import numpy as np
print('time0, for standardimports',timeit.timeit() - start); start = timeit.timeit()

# check those imports
import myutils as my
print('time1, for myutils',timeit.timeit() - start); start = timeit.timeit()

import fah
from myutils import ase_calculate_ene #as ace
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase import units as aseunits
print('time1, for ase_calc...',timeit.timeit() - start); start = timeit.timeit()

def help(p = None):
    string = '''
    Examples for using this script:
    -------------------------------
    % getEnergies_byLammps.py --units meV_pa -p Ni-Al-Mishin-2009.eam.alloy -thermo -v

    % getEnergies_byLammps.py -p n2p2_v1ag --units meV_pa -i input.data -idx 4850:
    % getEnergies_byLammps.py -p n2p2_v1ag --units meV_pa -i input.data -idx :4850
    % getEnergies_byLammps.py -p n2p2_v1ag --units hartree -i simulation.pos_0.xyz -fi ipi
    % getEnergies_byLammps.py -i $dotfiles/scripts/potentials/runner_v3ag_5000/input.data -p runner_v3ag_5000 -pfm 0.0001 -pc44
    % getEnergies_byLammps.py -i $dotfiles/scripts/potentials/runner_v3ag_5000/input.data -p runner_v3ag_5000 -pc44
    % getEnergies_byLammps.py -p . -e
    % getEnergies_byLammps.py -p runner_v3ag_5000_46489_2 --units meV_pa -i ../data.ipi -fi ipi -ru
    % getEnergies_byLammps.py -p runner_v3ag_5000_46489_2 --units meV_pa -i ../data.runnerformat.lmp -fi lammpsrunner
    % getEnergies_byLammps.py -p runner_v3ag_4998_3 -sys fcc -sys_ele Al -evinet
    % getEnergies_byLammps.py -p runner_v3ag_4998_3 -sys fcc -sys_ele Al -fqh
    % getEnergies_byLammps.py -p runner_v3ag_4998_3 -sys fcc_vac -sys_ele Al -fqh
    % getEnergies_byLammps.py -p runner_v3ag_4998_3 -sys fcc_selfinterstitial -sys_ele Al -fqh
    % getEnergies_byLammps.py -p n2p2_v5ag_ppl_98765_21cores -sys fcc_selfinterstitial -sys_ele Al -sys_ncell 2

    %
    % getEnergies_byLammps.py -p . -sys fcc -sys_ele Al -sys_ncell 1 -e -evinet -v # calculates c44 for al (TABLE I. Kobayashi)
    % getEnergies_byLammps.py -p . --formation_energies beta2 -v --units eV  # get formation energies
    % getEnergies_byLammps.py  -p Al_zhou.eam.alloy  -sys fcc -sys_ele Al -sys_ncell 1 -e -evinet

    % for paper:
        evinet (a,B,B')   : getEnergies_byLammps.py -p . -sys fcc -sys_ele Al -evinet
        elastic constants : getEnergies_byLammps.py -p . -sys fcc -sys_ele Al -sys_ncell 1 -ea
        heats of solute/vac formation: getEnergies_byLammps.py -p . --formation_energies beta2 -v --units eV
        formation energies AlXMgYSiZ : getEnergies_byLammps.py -p . --formation_energies beta2 -v --units eV
        antisite form energies       : getEnergies_byLammps.py -p .. -v --units eV --formation_energies beta2antisites
        interstitials                : getEnergies_byLammps.py -p .. -i ../../aiida_get_structures_new/aiida_exported_group_Al6xxxDB_structures_calc__all_steps.input.data -v --pick_atoms_al 33 --units meV_pa --interstitial_form

        notes on antisites:
        readpath /home/glensk/Dropbox/Albert/scripts/dotfiles//scripts/potentials/aiida_get_structures_new/aiida_exported_group_out_antisites_repMg4Al3Si4_NEW.runner_calc__all_steps.input.data   # DFT antisites
        savepath /scratch/glensk/2020_retrain_almgsi/zeroth20/random_seed_16560/epoch_False/RELAXED_fully_aiida_exported_group_out_antisites_repMg4Al3Si4_NEW.runner_calc__all_steps.input.data   # relaxed antisites




    # creating of structures for quantum-espresso
    # 1) getEnergies_byLammps.py -p runner_v3ag_4998_3 -i /Users/glensk/Dropbox/Albert/scripts/dotfiles//scripts/potentials/aiida_get_structures_new/OQMD-Distorted_scf --units hartree -we
    # 2) createFolders_quantumespresso_from_inputfile.py out.espresso-in_*
    # 3) for every foler sbatch submit.sh
    # 4) getEnergies_byLammps.py -p runner_v3ag_4998_3 -i aiida.out -fi espresso-out --units hartree -wr


    # theta/theta' phase transition:
    # getEnergies_byLammps.py --units meV_pa -p n2p2_alcu_v2dm_11 -i POSCAR_ThetaPrime -fi vasp -thermo --fqh_atoms_max 4 --fah_atoms_max 4

    # normal test to calc energies of structures
    # getEnergies_byLammps.py -i $nninp -p n2p2_v4ag_ppl_987654_28cores -vv

    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-i', '--inputfile', required=False, type=str,default=False, help="input files containing structures that will be imported by ase")
    p.add_argument('-sys','--sys', choices=[ "fcc", "bcc", "dc", 'hcp', 'fcc_vac', 'fcc_dilute', 'fcc_selfinterstitial' ], default=False, help="crystal_structure; can be set in case no inputfile is specified, define the system manually")
    p.add_argument('-sys_ele','--sys_ele', nargs='*', default=False, help="in case no inputfile is given, define the element manually")
    p.add_argument('-sys_ncell','--sys_ncell', type=int, default=1, help="in case no inutfile is given, define how often the primitive/conventional cell is repeated")
    p.add_argument('--sys_cell','-sys_cell',choices = ['cubic','primitive','orthorhombic'],default='cubic',help='Which crystal structure to use')

    p.add_argument('-efm','--exectue_function_myutils', required=False, type=str,default='', help="function to run from myutils.")
    p.add_argument('-ef','--execute_function', required=False, type=str,default='', help="function to run from this file.")
    p.add_argument('-fi','--format_in', required=False, type=str,default='runner', help="ase format for reading files")
    p.add_argument('-ru','--remove_unknown_elements_from_structures'  ,action='store_true',help='Reove unknown elements (at respective atom) from the input structures')
    p.add_argument('--potpath','-p',   required=False, type=str, default="n2p2_v4ag_ppl_987654_21cores", help="In case --pot is set to setpath use --potpath Folder to point to the Folder containing the n2p2/runner potential")
    p.add_argument('--show_availabel_pots','-sp',action='store_true',help='show available potentials from $potentials')
    p.add_argument('--potepoch','-pe',  required=False, type=int, default=False, help="use particular epoch of the potential")
    p.add_argument('--structures_idx','-idx',default=':',help='which structures to calculate, use ":" for all structues (default), ":3" for structures [0,1,2] etc. (python notation)')
    p.add_argument('--units','-u',choices = ['eV','meV_pa','eV_pa','hartree','hartree_pa'],default='hartree_pa',help='In which units should the output be given')
    p.add_argument('--interstitial_form','-interstitial_form' ,action='store_true',help='write interstitials_NN_vs_DFT.txt',default=False)
    p.add_argument('--geopt','-g'               ,action='store_true',help='make a geometry optimization of the atoms.')
    p.add_argument('--elastic','-e'             ,action='store_true',help='calculate elastic constants with given potential for Al (externally by lammps).')
    p.add_argument('--thermo','-thermo'         ,action='store_true',help='calculate free energy surface.')
    p.add_argument('--evinet','-evinet'         ,action='store_true',help='calculate evinet.')
    p.add_argument('--fqh','-fqh'               ,action='store_true',help='calculate fqh.')
    p.add_argument('--fah','-fah'               ,action='store_true',help='calculate fah.')
    p.add_argument('--fqh_atoms_max','-fqh_atoms_max'           ,default=200,type=int,help='maximum number of atoms for supercells used for fqh')
    p.add_argument('--fqh_supercell','-fqh_supercell'           ,default=None,type=int,help='amoun of repetitions of supercell for fqh calculations')
    p.add_argument('--fah_atoms_max','-fah_atoms_max'           ,default=200,type=int,help='maximum number of atoms for supercells used for fah')
    p.add_argument('--fah_supercell','-fah_supercell'           ,default=None,type=int,help='amoun of repetitions of supercell for fah calculations')
    p.add_argument('--elastic_all',     '-ea'   ,action='store_true',help='calculate elastic constants for every epoch for Al (externally by lammps).')
    p.add_argument('--elastic_from_ene','-ee'   ,action='store_true',help='calculate elastic constants for Al from energies (inernally by ase).')
    p.add_argument('--ase'    ,'-a'             ,action='store_true',default=True,help='Do the calculations by the ase interface to lammps.')
    p.add_argument('--lmp'    ,'-lmp'           ,action='store_true',help='Do the calculations externally by lammps and not through ase interface.')
    p.add_argument('--ipi'    ,'-ipi'           ,action='store_true',help='Do the calculations externally by ipi-lammps and not through ase interface.')

    p.add_argument('--test_this_script','-t'    ,action='store_true',help='check if getEnergies_byLammps.py is working correctly.')
    #p.add_argument('--formation_energies','-fe',  action='store_true',help='Assess formation energies of particular test structures.')
    p.add_argument('--formation_energies','-fe',  choices=[ 'inputfile', 'beta2','beta2antisites' ],help='Assess formation energies of particular test structures.')
    p.add_argument('--test3'  ,'-t3'            ,action='store_true',help='test3')
    p.add_argument('--testkmc'  ,'-kmc'         ,action='store_true',help='test accuracy of kmc structures')
    p.add_argument('--testkmc_b','-kmcb'        ,action='store_true',help='test accuracy of kmc structures for epoch with best_test energies')
    p.add_argument('--testkmc_l','-kmcl'        ,action='store_true',help='test accuracy of kmc structures for last epoch')
    p.add_argument('--testkmc_a','-kmca'        ,action='store_true',help='test accuracy of kmc structures for all epochs')
    p.add_argument('--testaccuracy_kmc_approx','-kaka'        ,action='store_true',help='test accuracy of fomation energy in kmc wehn substituting outer shells by Al')
    p.add_argument('--check_testdata', '-ctest' ,action='store_true',help='test accuracy of test.data')
    p.add_argument('--check_traindata','-ctrain',action='store_true',help='test accuracy of train.data')
    p.add_argument('--check_inputdata','-cinput',action='store_true',help='test accuracy of input.data structures used')
    p.add_argument('--check_kmc57data','-ckmc'  ,action='store_true',help='test accuracy of kmc57.data')
    p.add_argument('--check_outliers','-co'     ,action='store_true',help='test for outliers')

    p.add_argument('--analyze_kmc_number_1NN_2NN_ext','-akmc_ext',action='store_true',help='make simulation.pos_0.xyz.1NN.al_mg_si_vac_0.dat files')
    p.add_argument('--analyze_kmc_number_1NN_2NN_ipi','-akmc_ipi',action='store_true',help='make analysis from KMC_analyze')
    p.add_argument('--analyze_kmc_number_1NN_2NN_post','-akmc_post',action='store_true',help='make analysis from KMC_analyze')

    p.add_argument('--pick_concentration_al','-pcal',default=-1.,type=float,help='only consider structures with particular concentration of element, e.g. -pcal 1.0')
    p.add_argument('--pick_concentration_mg','-pcmg',default=-1.,type=float,help='only consider structures with particular concentration of element, e.g. -pcmg 1.0')
    p.add_argument('--pick_concentration_si','-pcsi',default=-1.,type=float,help='only consider structures with particular concentration of element, e.g. -pcsi 1.0')
    p.add_argument('--pick_atoms_al','-paal',default=-1.,type=float,help='only consider structures with particular number of al atoms, e.g. -paal 106 (e.v. 106 of 108)')
    p.add_argument('--pick_atoms_mg','-pamg',default=-1.,type=float,help='only consider structures with particular number of mg atoms, e.g. -pamg 106 (e.v. 106 of 108)')
    p.add_argument('--pick_atoms_si','-pasi',default=-1.,type=float,help='only consider structures with particular number of si atoms, e.g. -pasi 106 (e.v. 106 of 108)')
    p.add_argument('--pick_number_of_atoms','-pnat',default=-1.,type=float,help='only consider structures with particular number of atoms, e.g. -pnat 107')
    p.add_argument('--pick_forcesmax','-pfm',default=-1.,type=float,help='only consider structures with particular max force, e.g. -pfm 0')
    p.add_argument('--pick_cellshape','-pcs',default=-1.,type=float,help='only consider structures with particular cellshape, e.g. -pfm 0')
    p.add_argument('--pick_c44','-pc44'         ,action='store_true',default=False,required=False,help='only consider structures which are candidates for c44 calculations')
    p.add_argument('--pick_amount_1NN','-pa_1NN',action='store_true',default=False,required=False,help='detrmine the amount of Si,Mg,Al in 1NN shell around vacancy')
    p.add_argument('--pick_uuid','-pick_uuid'     ,default=-1.,type=str, nargs='*',required=False,help='detrmine uuis that will be calculated.')

    p.add_argument('--write_runner','-wr'    ,  action='store_true',help='write runnerfile from calculated structures using chosen pot; default filename out: runner.out')
    p.add_argument('--write_runner_ene_forces','-wref'    ,  action='store_true',help='write runnerfile from calculated structures using chosen pot; default filename out: runnerFORCES.out')
    p.add_argument('--write_DFT_ene_forces','-wrDFT'    ,  action='store_true',help='write runnerfile from calculated DFT structures; default filename out: DFT.runner')
    p.add_argument('--write_espresso_job','-we'    ,  action='store_true',help='write inputjobs for quantum espresso to be started on daint, instead of aiida.')
    p.add_argument('--write_runner_DFT','-wrd', action='store_true',help='write runnerfile with selected input structures (e.g. DFT); default filename out: runner_DFT.out')
    p.add_argument('--write_runner_repeated','-wrr',  action='store_true',help='default: runner_repeated.out')
    p.add_argument('--write_forces','-wf',  action='store_true',help='write forces.out of particular strucuture')
    p.add_argument('--write_forcesx','-wfx',action='store_true',help='write forcesx.out of particular strucuture')
    p.add_argument('--write_analysis','-wa',action='store_true',help='write ene_{DFT,pot}... default: False')
    p.add_argument('--write_analysis_full','-waf',action='store_true',help='write ene_{DFT,pot}... default: False')
    p.add_argument('-poe','--print_only_energies',   help='print only the energies and no additional ifno', action='count', default=False)
    p.add_argument('-ghe','--get_harmonic_energy',   help='get the harmonic energy', action='count', default=False)
    p.add_argument('-d','--debug',   help='verbose', action='count', default=False)
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    p.add_argument('-uuid','--uuid', help='show uuid of structure', action='count', default=False)
    return p

print('time3, after help  ...',timeit.timeit() - start); start = timeit.timeit()

def get_energies(args):
    ''' this is a script which computes for a given set of structures the energies
    for a given potential.

    '''
    start = timeit.timeit()
    ##############################################################
    ### check if ase runner/quippy/lammpps-data formats are known
    ### use -v (verbose) option to see known formats
    ##############################################################
    ase_formats = my.ase_get_known_formats_class(verbose=args.verbose)
    ase_formats.check_if_default_formats_known(copy_and_adapt_formatspy_anyhow=False)
    print('time00 after ase_formats_check...',timeit.timeit() - start); start = timeit.timeit()


    #frames = ase_read(args.inputfile,format=args.format_in,index=":")
    #ase_write("DFT0.runner",frames,format='runner',append=True)

    #dudl = fah.get_dudl_from_file_with_energies_lambda_0_1('../simulation.ti',number_of_atoms=32)
    #dudlav = get_dudlav_from_dudl(dudl)
    #save_dudl_to_file('../',dudl,filename="dudl"):
    #if os.path.isfile("../simulation.ti"):
    #    aaatest = np.loadtxt("../simulation.ti")
    #else:
    #    aaatest = np.zeros((1000,3))

    #if args.get_harmonic_energy:
    #    # needs to determine u and h
    #    # u can be obtained from normal positions probably.
    #    import hesse
    #    hessefile = '../HesseMatrix_4.13'
    #    if not os.path.isfile(hessefile):
    #        hessefile = '../hessian.data'
    #    hesse_matrix = hesse.read_Hessematrix(hessefile)

    #    if os.path.isfile('../init.xyz'):
    #        pos0 = ase_read('../init.xyz',format='ipi')

    #    elif os.path.isfile('../OUTCAR'):
    #        print('reading outcar')
    #        #pos0 = ase_read('../OUTCAR',format='vasp-out') # this is as index=0
    #        #pos0 = ase_read('../OUTCAR',format='vasp-out',index=":") # works
    #        pos0 = ase_read('../OUTCAR',format='vasp-out',index=0) # works
    #        ase_write('POSCAR',pos0,format='vasp')
    #        ase_write('data.lmp',pos0,format='lammps-data')
    #        ang_to_bohr = 1.8897261
    #        np.savetxt("x_reference.data",pos0.positions.flatten()*ang_to_bohr,newline=" ",fmt="%3.15f")
    #        ase_write("init.xyz",pos0,format='ipi')
    #    print('pos0 typ',type(pos0))
    #    print('pos0 len',len(pos0))
    #    print(pos0.positions)
    #    print(pos0.cell)


    if args.show_availabel_pots:
        my.list_pot_all()
        sys.exit()

    ##################################
    # create the README
    ##################################
    hier = os.path.abspath(os.getcwd())
    readmepath = my.create_READMEtxt(hier)

    allepochs = [False]
    format_in = args.format_in

    units = args.units
    debug = args.debug
    verbose = args.verbose

    geopt = args.geopt
    ase = args.ase
    lmp = args.lmp
    ipi = args.ipi

    if args.analyze_kmc_number_1NN_2NN_ipi:
        ''' -akmc_ipi '''
        file = "KMC_analyze_head"
        file = "KMC_analyze"
        if not os.path.isfile(file):
            sys.exit(str(file)+" does not exist!")

        a = np.loadtxt(file)
        np.set_printoptions(suppress=False)
        #print(a)
        cdf = a[:,4]
        nrand2 = a[:,5]
        nrandn = np.ones(len(cdf))/2.

        al = a[:,6]
        mg = a[:,7]
        si = a[:,8]
        dt = -1.0/cdf*np.log(1.0 - nrand2)*2.418884326509e-17
        dtn = -1.0/cdf*np.log(1.0 - nrandn)*2.418884326509e-17
        #dt2 = np.zeros(len(cdf))
        #for i in np.arange(len(cdf)):
        #    dt2[i] = -1.0/cdf[i]*np.log(1.0 - nrand2[i])*2.418884326509e-17
        #    print('i',i,'cdf',cdf[i],'nrand2',nrand2[i],'dt2',dt2[i])
        if args.verbose:
            print("cdf",cdf)
            print("nrand2",nrand2)
            print("nrandn",nrandn)
            print("al",al)
            print("mg",mg)
            print("si",si)
            print('dt ',dt)
            print('dtn',dtn)
            #print('dt2',dt2)
        print('len(si)',len(si))
        print('len(mg)',len(mg))
        mode = 'valid'
        mode = 'full'
        mode = 'same'
        N = 400
        si_av = np.convolve(si, np.ones((N,))/N, mode=mode)
        mg_av = np.convolve(mg, np.ones((N,))/N, mode=mode)
        simg_av = np.convolve(si+mg, np.ones((N,))/N, mode=mode)
        dtn_av = np.convolve(dtn, np.ones((N,))/N, mode=mode)
        print('len(si_av)',len(si_av))
        print('len(mg_av)',len(mg_av))

        out = np.zeros((len(cdf),5))
        out[:,0] = np.arange(len(cdf))
        out[:,1] = dtn_av
        out[:,2] = si_av
        out[:,3] = mg_av
        out[:,4] = simg_av
        np.savetxt(file+"_akmc_ipi.dat",out)
        print('written',file+"_akmc_ipi.dat")
        sys.exit()

    if args.analyze_kmc_number_1NN_2NN_ext:
        ''' -akmc_ext '''
        filename = "simulation.pos_0.xyz"
        #filename = '../sim.xyz'
        #filename = 'sim.xyz'
        # 4.23 is too little for Si6Mg6V1.2_ step 797
        # 4.5  is too much   for Si6Mg6V1.4  step 2372
        #
        # 3.5 is too much    for Si6Mg6V1.4 step 2242

        my.count_amount_1NN_around_vacancies(filename,cutoffa=3.4,cutoffb=4.4,skin=0.1,format='ipi')
        sys.exit()

    if args.analyze_kmc_number_1NN_2NN_post:
        a = np.loadtxt("KMC_analyze_akmc_ext")
        b = np.loadtxt("KMC_analyze_np_analyze.dat")
        print('a',a.shape)
        print('b',b.shape)
        sys.exit()

        fig, ax = plt.subplots(1, 1)
        #mesh = ax.scatter(projs[:,0],projs[:,1],c=range(projs.shape[0]),s=5,alpha=1,cmap="jet")
        mesh = ax.scatter(a[:,2],a[:,3],c=range(projs.shape[0]),s=a[:,4],alpha=1,cmap="jet")
        fig.colorbar(mesh)
        plt.savefig(r'corr_residency_time.pdf', transparent=True,dpi=600,bbox_inches='tight', pad_inches=0)


        al1 = a[:,1]
        mg1 = a[:,2]
        si1 = a[:,3]
        al2 = a[:,4]
        mg2 = a[:,5]
        si2 = a[:,6]

        sys.exit()


    if args.pick_c44:
        if args.units != "meV_pa":
            sys.exit('for pick_c44 units have to be in meV_pa')
        args.verbose = True


    hostname = my.hostname()
    if lmp == True:
        ase = False



    ##############################################################
    ### check if lammps is working with ase
    ##############################################################
    if 'LD_LIBRARY_PATH' not in os.environ:
        os.environ['LD_LIBRARY_PATH'] = os.environ['HOME']+'/sources/lammps/src'
    print('getEne(0.1) LD_LIBRARY_PATH            :',os.environ['LD_LIBRARY_PATH'])
    from lammps import lammps
    print('getEne(0.2) imported lammps')
    lammps()
    print('getEne(0.3) lammps')
    os.remove("log.lammps")
    print('getEne(0.4) check if lammps woring with ase done')


    ##############################################################
    ### POTENTIAL
    ### get ace object for the chosen potential (first general)
    ### just to check if the elements (structure vs pot) are the right ones
    ##############################################################
    if args.formation_energies or args.elastic: units='eV'
    print('getEne(p1) args.potepoch             :',args.potepoch)
    ace = ase_calculate_ene(potpath_in=args.potpath,
            use_epoch=args.potepoch,
            units=args.units,
            geopt=geopt,
            elastic=args.elastic,
            fqh_atoms_max = args.fqh_atoms_max,
            fqh_supercell = args.fqh_supercell,
            fah_atoms_max = args.fah_atoms_max,
            fah_supercell = args.fah_supercell,
            verbose=verbose)
    ace.lammps = args.lmp
    ### get the potential
    print('getEne(p2)                           : ace.pot_get_and_ase_lmp_cmd()')
    ace.pot_get_and_ase_lmp_cmd()  # just to have lmpcmd defined in case ...
    if args.verbose > 2:
        print('getEne(p3)                           : pot defined')
    print('ace.lmpcmd',ace.lmpcmd)
    units = ace.units
    ace.pot.print_variables_mypot(print_nontheless=True,text="getEne(P):")
    if args.verbose > 2:
        print('getEne(p4)                           : pot defined')

    ##############################################################
    ### again show args
    ##############################################################
    if args.verbose:
        my.print_args(args)
    if args.verbose > 2:
        print('getEne(p5)                           : pot defined')


    ###################################################################
    ### test accuracy formation energy when substituting Atoms with Al
    ###################################################################
    if args.testaccuracy_kmc_approx:
        my.analyze_accuracy_of_filling_cell_with_Al(ace)
        sys.exit("kaka done")
    if args.verbose > 2:
        print('getEne(p6)                           : pot defined')


    ###################################################################
    ### get the input structures
    ###################################################################
    print('getEne(0) args.inputfile             :',args.inputfile)
    if args.test_this_script or args.testkmc or args.testkmc_b or args.testkmc_l or args.testkmc_a:
        args.inputfile = os.environ["dotfiles"]+"/scripts/potentials/aiida_get_structures_new/aiida_exported_group_KMC57.data"
        units = "meV_pa"
        verbose = True
        if args.test_this_script:
            args.structures_idx = ":3"

    print('392')
    if args.inputfile == False and args.sys != False and args.sys_ele != False:
        if args.sys in [ 'fcc', 'bcc', 'hcp', 'dc', 'fcc_selfinterstitial' ] and len(args.sys_ele) == 1:
            frames = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=args.sys_ncell,nsi=0,nmg=0,nvac=0,matrix_element=args.sys_ele[0],a0=False,cell=args.sys_cell,create_fake_vacancy=False,crystal_structure=args.sys)
        elif args.sys in [ 'fcc_dilute' ] and len(args.sys_ele) == 1:
            sys.exit('please specify 2 elements using -sys_ele e.g. -sys_ele Al Si')
        elif args.sys in [ 'fcc_dilute' ] and len(args.sys_ele) == 2:
            frames = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=args.sys_ncell,nsi=0,nmg=0,nvac=0,matrix_element=args.sys_ele,a0=False,cell=args.sys_cell,create_fake_vacancy=False,crystal_structure=args.sys,verbose=args.verbose)
        if args.verbose:
            print('created structure from given input!')
            print('cell')
            print(frames.cell)
            for idx,i in enumerate(frames.positions):
                print('before relax:',idx,frames.get_chemical_symbols()[idx],frames.positions[idx])
            print('relaxing the structure ...')
        print('406')
        ace.ase_relax_atomic_positions_only(frames,fmax=0.0001,verbose=args.verbose)
        if args.verbose:
            for idx,i in enumerate(frames.positions):
                print('after relax:',idx,frames.get_chemical_symbols()[idx],frames.positions[idx],frames.get_forces()[idx])

    print('410')
    if args.inputfile != False:
        my.check_isfile_or_isfiles([args.inputfile],verbose=args.verbose)
        print('########################################################')
        print('XXX args.inputfile:',args.inputfile)
        print('XXX args.format_in:',args.format_in)
        print('########################################################')
        frames = ase_read(args.inputfile,format=args.format_in,index=":")
        #ase_write("DFT1.runner",frames,format='runner',append=True)
        print('type(frames)',type(frames))
        print('len(frames)',len(frames))
        DFT_FORCES_EV_ANG = []
        DFT_ENERGIES_EV_CELL = []
        print('getting DFT data')
        for idx,i in enumerate(frames):
            #print(frames[idx].get_forces())
            try:
                forces_thisstruct = frames[idx].get_forces()
            except RuntimeError:
                forces_thisstruct = False
            DFT_FORCES_EV_ANG.append(forces_thisstruct)
            try:
                ene_thisstruct = frames[idx].get_potential_energy()
            except RuntimeError:
                ene_thisstruct = False
            DFT_ENERGIES_EV_CELL.append(ene_thisstruct)
        print('getting DFT data DONE!')
        POT_FORCES_EV_ANG = np.copy(DFT_FORCES_EV_ANG)
        POT_ENERGIES_EV_CELL = np.copy(DFT_ENERGIES_EV_CELL)
        if False:  # this prints all the forces from all the structures
            print('DFT_FORCES_EV_ANG')
            print(DFT_FORCES_EV_ANG)
        #print('max',np.abs(DFT_FORCES_EV_ANG).max())
        #np.savetxt('maxforce.txt',np.abs(DFT_FORCES_EV_ANG).max(),fmt='%2.4f')
        if args.write_DFT_ene_forces:
            ase_write("DFT.runner",frames,format='runner',append=True)
        #sys.exit()

    print('430')
    if args.thermo or args.evinet or args.fqh or args.fah:
        # Al evinet: -537461.993661476416 16.579676546844 79.185555426019 2.526937653000
        # Si evinet: -155368.146177827992 20.424692587262 89.078642917550 4.110638394285
        # for one si in al matrix (with 108 atoms):
        # per atom:
        # -533920.334868975333 16.528710460497 79.180381312792 2.610462793516
        # -533920.334868975333*108-(107*-537461.993661476416)-(-155368.146177827992) == 405.302 meV ! that is correct!

        #if args.thermo or args.evinet or args.fqh or args.fah:
        #    if os.path.isdir('evinet'): sys.exit("Exit: Folder evinet exists already!")
        #if False and not os.path.isdir('fah'):

        if False and os.path.isdir('fah'):
            print("##############################")
            print("# extracting anharmonic Fah_surface... #")
            print("##############################")
            hier = os.getcwd()
            print('os.getcwd()',os.getcwd())
            os.chdir(hier+"/fah")

            fvta = glob.glob(os.getcwd()+"/*_*K")
            print('fvta',fvta)
            for fvt in fvta:
                os.chdir(fvt)
                print('os.getcwd() fvt',os.getcwd())
                flsa = glob.glob(os.getcwd()+"/lambda*_*")
                for fls in flsa:
                    os.chdir(fls)
                    print('os.getcwd() fls',os.getcwd())

                    if False:  # read dudl
                        dudl = fah.get_dudl_from_file_with_energies_lambda_0_1('simulation.ti',number_of_atoms=32)
                        dudlav = fah.get_dudlav_from_dudl(dudl)
                        print('dudav[-1]',dudlav[-1])
                    if True:  # run anharmonic
                        if os.path.isfile(os.getcwd()+'/simulation.ti'):
                            print('simulation.ti exists in:')
                            print(os.getcwd()+'/simulation.ti')
                            sys.exit()
                        print('do 1')
                        my.ipi_start_job(inputfile="input.xml")
            sys.exit()



        if args.thermo: args.evinet = True
        if args.thermo: args.fqh = True
        if args.thermo: args.fah = True
        if args.inputfile == 'POSCAR': args.format_in = "vasp"
        if args.verbose > 2:
            print('getEne(p7) reading structure ... ')
        if type(frames) == list and len(frames) == 1:
            print('len(frames)',len(frames))
            frames = frames[0]

        print('args.format_in',args.format_in)
        print('frames.cell',frames.cell)
        print('frames.positions',frames.positions)
        print('frames.get_atomic_numbers()',frames.get_atomic_numbers())
        print('frames.get_atomic_numbers()',list(set(frames.get_atomic_numbers())))
        print('os.getcwd() BEFORE GET CALCULATOR',os.getcwd())
        ace.get_calculator(frames)
        print('forces YYY',frames.get_forces())
        ase_structure_relaxed = my.get_thermo(ace,frames,relax_cellshape_and_volume=True,evinet=args.evinet,fqh=args.fqh,fah=args.fah)
        print('nat ase_structure_relaxed',ase_structure_relaxed.get_number_of_atoms())

        #ace.get_elastic_external(atomsin=ase_structure_relaxed,verbose=ace.verbose,text="structure",get_all_constants=True)
        #print('ace.c44:',ace.c44,type(ace.c44))
        print('done evinet')
        #sys.exit("Thermodynamic properties done up to Fqh. Fah jobs created")
    print('501')
    if args.verbose > 2:
        print('getEne(p7)                           : pot defined')

    ############
    ### testkmc
    ############
    if args.testkmc or args.testkmc_b or args.testkmc_l or args.testkmc_a:
        if args.testkmc_b:
            allepochs = args.potepoch = ace.pot.use_epoch = [ace.pot.potepoch_bestteste]
            print('args.potepoch',args.potepoch)
            print('ace.pot.use_epoch',ace.pot.use_epoch)
        if args.testkmc_l:
            allepochs = args.potepoch = ace.pot.use_epoch = [ace.pot.potepoch_all[-1]]
            print('args.potepoch',args.potepoch)
            print('ace.pot.use_epoch',ace.pot.use_epoch)
        if args.testkmc_a:
            args.potepoch = ace.pot.use_epoch = ace.pot.potepoch_all[-1]
            allepochs = ace.pot.potepoch_all
            print('args.potepoch',args.potepoch)
            print('ace.pot.use_epoch',ace.pot.use_epoch)

        if args.potepoch == False:
            sys.exit('Error: need to specify a particular epoch for kmctest')

        kmc_folder = ace.pot.potpath+"/kmc"
        if not os.path.isdir(kmc_folder):
            my.mkdir(kmc_folder)

    ############
    ### formation energies
    ############
    if args.formation_energies:
        formation_energies(ace,args)
        sys.exit('formation_energies done! Exit')

    if args.exectue_function_myutils:
        print('args.inputfile',args.inputfile)
        print('args.exectue_function_myutils',args.exectue_function_myutils)
        function = eval('my.'+args.exectue_function_myutils)
        #my.get_Mg5Si6_and_other_antisites(ace)
        function(ace)
        sys.exit()

    if args.execute_function:
        print('args.inputfile',args.inputfile)
        print('args.execute_function',args.execute_function)
        function = eval(args.execute_function)
        #my.get_Mg5Si6_and_other_antisites(ace)
        function(ace)
        sys.exit()


    ############
    ### elastic / elastic_all
    ############
    if args.elastic:
        get_elastic_constants_al_ext(ace,frames)
        sys.exit('get_elastic_constants_al_ext done! Exit')

    if args.elastic_all:  # this does not need to initialize potential
        file_elastic_c44_all = ace.pot.potpath+"/elastic_c44_all.dat"
        if os.path.isfile(file_elastic_c44_all):
            elastic_all = np.loadtxt(file_elastic_c44_all)
        else:
            elastic_all = np.empty((0,2), float)

        writeanew = False
        goover = ace.pot.potepoch_all
        goover = np.arange(1,goover[-1]+1)
        if args.potepoch != False:
            goover = np.array([args.potepoch])
        goover = np.array([ace.pot.use_epoch])
        print('goover',goover)
        count = 0
        for epoch in goover: #ace.pot.potepoch_all: #[100]: #np.arange(1,100): #[7]: #ace.pot.potepoch_all:
            print('ELASTIC_ALL',elastic_all)
            if epoch in elastic_all[:,0]:
                print('epoch',str(epoch).ljust(5),'from file')
            else:
                count += 1
                print('epoch',str(epoch).ljust(5),'not already in, count('+str(count)+")")
                ace = ase_calculate_ene(
                        potpath_in=args.potpath,
                        use_epoch=epoch,
                        units=args.units,
                        geopt=geopt,
                        elastic=args.elastic,
                        verbose=verbose)
                ace.lammps = args.lmp
                #print('now getting the lmp cmd')
                ace.pot_get_and_ase_lmp_cmd()  # need to update lmp_cmd when changing the potential
                c44 = get_elastic_constants_al_ext(ace,frames)
                writeanew = True
                #elastic_all.append([epoch,np.float(c44)])
                #print('epoch',str(epoch).ljust(5),c44)
                elastic_all= np.append(elastic_all, np.array([[epoch,np.float(c44)]]), axis=0)
            if count == 20:
                print()
                print('saving ...')
                print()
                elastic_all = elastic_all[elastic_all[:,0].argsort()]
                np.savetxt(file_elastic_c44_all,elastic_all)
                count = 0

        if True: #writeanew:
            elastic_all = elastic_all[elastic_all[:,0].argsort()]
            np.savetxt(file_elastic_c44_all,elastic_all)
        sys.exit('elastic_all done! Exit')

    if args.elastic_from_ene:
        get_elastic_constants_al_from_ene(ace)
        sys.exit('get_elastic_constants_al_from_ene done! Exit')

    if args.test3:
        print('>> args.test3')
        test3_do(ace)
        sys.exit('args.test3 done! Exit')

    if args.check_inputdata or args.check_traindata or args.check_testdata or args.check_kmc57data:
        units = ace.units = "meV_pa"
        args.write_analysis = True
        if args.inputfile == False and args.check_inputdata:
            args.inputfile = ace.pot.inputdata
            ext = "input"
        if args.inputfile == False and args.check_traindata:
            args.inputfile = ace.pot.traindata
            ext = "train"
        if args.inputfile == False and args.check_testdata:
            args.inputfile = ace.pot.testdata
            ext = "test"
        if args.inputfile == False and args.check_kmc57data:
            args.inputfile = os.environ["dotfiles"]+"/scripts/potentials/aiida_get_structures_new/aiida_exported_group_KMC57.data"
            ext = "kmc57"
        print('pa',ace.pot.potepoch_bestteste)
        print('pb',ace.pot.potepoch_all,type(ace.pot.potepoch_all))
        if args.potepoch == False:
            allepochs = [ace.pot.potepoch_bestteste]
        else:
            allepochs = [args.potepoch]

        #allepochs = ace.pot.potepoch_all

        #print('allepochs 0',allepochs)
        #if type(ace.pot.potepoch_all) == bool:
        #    pass
        #elif ace.pot.potepoch_all == []:
        #    pass
        #elif type(ace.pot.potepoch_all) == list:
        #    allepochs = ace.pot.potepoch_all

        #print('allepochs 1',allepochs)
        #if type(ace.pot.potepoch_bestteste) == bool:
        #    pass
        #elif ace.pot.potepoch_bestteste == []:
        #    pass
        #elif type(ace.pot.potepoch_bestteste) == int:
        #    allepochs = np.append(allepochs,np.array([ace.pot.potepoch_bestteste]))
        #print('allepochs 2',allepochs)


        #if len(ace.pot.potepoch_bestteste) == 0 and type(ace.pot.potepoch_all) == bool:
        #    allepochs = []
        #elif len(ace.pot.potepoch_bestteste) != 0 and type(ace.pot.potepoch_all) == bool:
        #    allepochs = [ace.pot.potepoch_bestteste]
        #elif len(ace.pot.potepoch_bestteste) != 0 and type(ace.pot.potepoch_all) != bool:
        #    allepochs = [ace.pot.potepoch_bestteste,ace.pot.potepoch_all[-1]]

        print('allepochs',allepochs)

    ###############################
    ### check args.inputfile
    ###############################
    if not args.inputfile:
        sys.exit("Error: Missing option \"--infile\" / \"-i\".")
    if args.inputfile == 'POSCAR': args.format_in = "vasp"
    my.check_isfile_or_isfiles([args.inputfile],verbose=verbose)

    ##################################################################################
    # go over every chosen potential
    ##################################################################################
    if args.potepoch != False:
        allepochs = [args.potepoch]
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    print('@@  allepochs:',allepochs)
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    for use_epoch in allepochs:
        os.chdir(hier)
        print()
        print()
        print("Epoch:",use_epoch)

        #######################################################################
        ### read in the structures / input.data
        #######################################################################
        print('getEne (23): reading args.inputfile ...  :',args.inputfile)
        if args.inputfile[-6:] == 'extxyz': args.format_in = "extxyz"
        if args.inputfile[-6:] == 'gz': sys.exit("this need to be unzipped first, otherwise takes too long with ase")
        print('YYY getEne (23): args.format_in ...          :',args.format_in)
        print('YYY getEne (23): args.structures_idx ...     :',args.structures_idx)
        #print('try with guessing')
        frames = ase_read(args.inputfile,index=args.structures_idx,format=args.format_in)
        #print('fc')
        #print(frames[0].cell)
        #sys.exit()

        try:
            a = frames[0].cell     # this is just minimal check weather
            v = frames[0].get_volume()  # minimal check,
            #print('cell',a)
            #print(frames[0].get_chemical_symbols())
        except (IndexError,ValueError): # was probably a different filetype
            #print('try with fix ')
            frames = ase_read(args.inputfile,index=args.structures_idx,format=args.format_in)

        try:
            a = frames[0].cell     # this is just minimal check weather
            v = frames[0].get_volume()  # minimal check, (wont work for ipi files)
        except (IndexError,ValueError): # was probably a different filetype
            #print('try with ipi')
            frames = ase_read(args.inputfile,index=args.structures_idx,format='ipi')

        a = frames[0].cell          # this is just minimal check
        v = frames[0].get_volume()  # minimal check
        #print('-->fp vol',frames[0].get_volume())
        #print('-->fp cel',frames[0].cell)

        #sys.exit()


        if type(frames) == list: structures_to_calc = len(frames)
        else: structures_to_calc = 1
        print('getEne (24): structures_to_calc          :',structures_to_calc)

        if structures_to_calc == 0: sys.exit('ERROR! 0 strucutres to calculate? Something \
                went wrong when importing the args.inputfile \"'+args.inputfile+'\" (mayby \
                you need to change the \"args.format_in\" of the file? \
                currently \"'+args.format_in+'\"')

        #################################################################################
        ### set the corresponding potential
        #################################################################################
        if use_epoch == False:
            pass
        else:
            ###############################
            ## define the actual pot
            ###############################
            # the ase_calculate_ene(...) changes the used epoch since mypot(...) is called
            ace = ase_calculate_ene(potpath_in=args.potpath,
                            use_epoch=use_epoch,
                            units=args.units,
                            geopt=geopt,
                            elastic=args.elastic,
                            verbose=args.verbose)
            ace.lammps = args.lmp
            ace.pot_get_and_ase_lmp_cmd()  # need to update lmp_cmd when changing the potential
            if args.testkmc or args.testkmc_b or args.testkmc_l or args.testkmc_a:
                kmc_file = kmc_folder+"/ene_std_epoch_"+str(use_epoch)+".dat"
                if os.path.isfile(kmc_file):
                    print(kmc_file+" does already exist!")
                    print('I will continue with other epoch...')
                    continue

            if args.check_inputdata or args.check_traindata or args.check_testdata or args.check_kmc57data:
                folder = ace.pot.potpath+"/assess_"+ext+"_"+str(use_epoch)
                print('folder',folder)
                if os.path.isdir(folder):
                    print("Exit, folder "+folder+" does already exist!")
                    continue
                else:
                    my.mkdir(folder)
                    os.chdir(folder)
                ace.pot.print_variables_mypot(print_nontheless=True,text="getEne(check_xxxdata):")





        # show positions of first structure?
        show_positions = False
        if show_positions:
            for _idx,_ppb in enumerate(frames[0].positions):
                print(_idx,frames[0].get_chemical_symbols()[_idx],frames[0].positions[_idx])
            print()
            print(frames[0].cell)
            sys.exit('Exit when show_ositions of first structure')

        print()
        print('ase                          :',args.ase)
        print('lmp                          :',args.lmp)
        print()
        print('args.units                   :',args.units)
        print('ace.units                    :',ace.units)

        print('geopt                        :',args.geopt)
        print()
        print('verbose                      :',args.verbose)
        print('lmp                          :',args.lmp)
        print('ipi                          :',args.ipi)
        if args.write_runner:
            args.write_runner = 'runner.out'
        print('write_runner                 :',args.write_runner)
        print('write_runner_repeated        :',args.write_runner_repeated)
        print('--pick_concentration_al      :',args.pick_concentration_al)
        print('--pick_concentration_mg      :',args.pick_concentration_mg)
        print('--pick_concentration_si      :',args.pick_concentration_si)
        print('--pick_atoms_al              :',args.pick_atoms_al)
        print('--pick_atoms_mg              :',args.pick_atoms_mg)
        print('--pick_atoms_si              :',args.pick_atoms_si)
        print('--pick_number_of_atoms       :',args.pick_number_of_atoms)
        print('--pick_forcesmax             :',args.pick_forcesmax)
        print('--pick_cellshape             :',args.pick_cellshape)
        print('--pick_c44                   :',args.pick_c44)
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
        ene_DFT_eV_cell  = np.empty(structures_to_calc);ene_DFT_eV_cell[:]  = np.nan
        ene_pot_eV_cell  = np.empty(structures_to_calc);ene_pot_eV_cell[:]  = np.nan
        ene_har_eV_cell  = np.empty(structures_to_calc);ene_har_eV_cell[:]  = np.nan
        ene_DFT_m0        = np.empty(structures_to_calc);ene_DFT_m0[:]  = np.nan
        ene_pot_m0     = np.empty(structures_to_calc);ene_pot_m0[:]  = np.nan
        ene_har_m0     = np.empty(structures_to_calc);ene_har_m0[:]  = np.nan
        ene_pot_ase_dudl = np.empty(structures_to_calc);ene_pot_ase_dudl[:]  = np.nan
        ene_pot_lmp_dudl = np.empty(structures_to_calc);ene_pot_lmp_dudl[:]  = np.nan
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


        dudl_pot_to_DFT      = np.empty(structures_to_calc);dudl_pot_to_DFT[:]  = np.nan
        dudl_har_to_DFT     = np.empty(structures_to_calc);dudl_har_to_DFT[:]  = np.nan
        dudl_har_to_pot     = np.empty(structures_to_calc);dudl_har_to_pot[:]  = np.nan
        dudlav_pot_to_DFT      = np.empty(structures_to_calc);dudlav_pot_to_DFT[:]  = np.nan
        dudlav_har_to_DFT     = np.empty(structures_to_calc);dudlav_har_to_DFT[:]  = np.nan
        dudlav_har_to_pot     = np.empty(structures_to_calc);dudlav_har_to_pot[:]  = np.nan

        for_DFTmax       = np.empty(structures_to_calc);for_DFTmax[:]  = np.nan
        uuids            = np.empty(structures_to_calc,dtype='|S1');uuids[:]  = ''
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
        if ace.units.split("_")[0].lower() == 'hartree':
            conv = 1.
        elif ace.units.split("_")[0].lower() == 'ev':
            conv = aseunits.Hartree #27.211384    # hartree to ev
        elif ace.units.split("_")[0].lower() == 'mev':
            conv = aseunits.Hartree*1000. #27211.384    # hartree to ev
        else:
            print('ace units',ace.units.split("_"),ace.units.split("_")[0].lower())
            sys.exit('ace units not known')

        #if True: #verbose > 1:
        #    print('kk',ace.pot.atom_energy)
        #    for atom_species in ace.pot.atom_energy:
        #        atom_energy = ace.pot.atom_energy[atom_species]
        #        print(atom_species, atom_energy)
        #    if "Mg" in ace.pot.atom_energy:
        #        print('yes mg exists')
        #sys.exit('77')
        #atom_energy_Mg = ace.pot.atom_energy["Mg"]*conv
        #atom_energy_Si = ace.pot.atom_energy["Si"]*conv
        #atom_energy_Al = ace.pot.atom_energy["Al"]*conv
        #if True: #verbose > 1:
        #    print("atom_energy_Mg",atom_energy_Mg,ace.units,"one_atom",ace.pot.atom_energy["Mg"],"hartee_pa")
        #    print("atom_energy_Si",atom_energy_Si,ace.units,"one_atom",ace.pot.atom_energy["Si"],"hartree_pa")
        #    print("atom_energy_Al",atom_energy_Al,ace.units,"one_atom",ace.pot.atom_energy["Al"],"hartree_pa")
        #sys.exit('77')

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
        print('range(structures_to_calc)',range(structures_to_calc))

        ################################################################
        # begin main loop
        ################################################################
        for idx,i in enumerate(range(structures_to_calc)):
            uuid = False
            try:
                all_comment = frames[i].info['comment']
                all_comment_split = all_comment.split()
                all_comment_split.remove('comment')
                all_comment_split.remove('uuid:')
                if len(all_comment_split) == 1:
                    uuid = all_comment_split[0]
                else:
                    print('776 Exit',all_comment)
                    print('776 Exit',all_comment_split)
                    sys.exit()
            except:
                all_comment = ""
            if verbose > 1:
                print('idxuu (11)',idx,'uuid:',uuid)
            if debug:
                print('iii',i)
            ana_atoms_ = frames[i].get_number_of_atoms()
            if debug:
                print('ama_atoms_',ana_atoms_)
            if verbose > 2:
                print('i',i,'ana_atoms',ana_atoms_,'uuid',uuid)
            try:
                DFT_forces = frames[i].get_forces()
                #print('DFT_forces')
                #print(DFT_forces)
                #sys.exit()
                for_DFTmax_ = np.abs(DFT_forces).max()
            except RuntimeError:
                for_DFTmax_ = -9999
            if debug:
                print('kk',for_DFTmax_)
                print('frame',frames[i].positions)
                print('spec',frames[i].get_chemical_symbols())
                print('pe',ace.pot.elements)
            d = my.ase_get_chemical_symbols_to_conz(frames[i],known_elements_by_pot=ace.pot.elements)
            if debug:
                print('d454:',d)
            n = my.ase_get_chemical_symbols_to_number_of_species(frames[i],known_elements_by_pot=ace.pot.elements)
            if debug:
                print('n454:',n)
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
            if args.pick_number_of_atoms >= 0 and ana_atoms_ != args.pick_number_of_atoms:
                continue
            if args.pick_concentration_al >= 0 and d["Al"] != args.pick_concentration_al:
                continue
            #print('conz Al',d["Al"],'atoms Al',n["Al"])
            #print('conz Si',d["Si"],'atoms Si',n["Si"])
            #print('conz Mg',d["Mg"],'atoms Mg',n["Mg"])

            look_for_particular_structurs = False
            if look_for_particular_structurs:
                found = False
                #if d["Mg"] > 0 and d["Si"] > 0 and d["Mg"]/d["Si"] == 5./6. and d["Al"] == 0:  # look for Mg5Si6
                #if d["Mg"] > 0 and d["Si"] > 0 and d["Al"] > 0 and d["Mg"]/d["Si"] == 5./4. and d["Mg"]/d["Al"] == 5./2.:  # look for Mg5Al2Si4
                if d["Mg"] > 0 and d["Si"] > 0 and d["Al"] > 0 and d["Mg"] == d["Si"] and d["Mg"]/d["Al"] == 4./3.:  # look for Mg4Al3Si4
                    found = True
                if found == False:
                    continue


            if args.pick_concentration_si >= 0 and d["Si"] != args.pick_concentration_si:
                continue
            if args.pick_concentration_mg >= 0 and d["Mg"] != args.pick_concentration_mg:
                continue
            if type(args.pick_uuid) == list and uuid not in args.pick_uuid:
                continue
            if args.pick_atoms_al >= 0 and n["Al"] != args.pick_atoms_al:
                continue
            if args.pick_atoms_mg >= 0 and n["Mg"] != args.pick_atoms_mg:
                continue
            if args.pick_atoms_si >= 0 and n["Si"] != args.pick_atoms_si:
                continue
            if args.pick_forcesmax >=0 and for_DFTmax_ > args.pick_forcesmax:
                continue
            if args.pick_cellshape >= 0 and cellshape not in ["Q"]:
                continue
            if args.pick_c44 == True:
                if n["Al"] != 4:
                    continue
                if n["Mg"] > 0:
                    continue
                if n["Si"] > 0:
                    continue
                if cellshape in [ "R", "Q" ]:
                    continue
                #if cellshape != "?":
                #print('nc')
                #continue
            #print('idx',idx,'n_al',n["Al"],"c_al",d["Al"],'consider_atoms_al cnat_al',consider_atoms_al,'cnat consider_number_of_atoms',consider_number_of_atoms)
            #print('idx',idx,'wow')
            if args.pick_c44 or debug: #cellshape == "?":
                #print('XX frames[i].cell')
                #print(frames[i].cell)
                strain = frames[i].cell[0,1]/frames[i].cell[0,0]
                print('XX volume',frames[i].get_volume(),'strain',round(strain,5))
            #if cellshape == "Q":
            #    print('frames[i].positions')
            #    print(frames[i].positions)

            #print('DFT_forces')
            #print(DFT_forces)
            if False:
                print('cell')
                print(cell)
                print('pos')
                print(pos)
            #sys.exit()
            if args.write_runner_DFT:
                try:
                    addedstr
                except NameError:
                    addedstr = 0

                ase_write("out_selected_uuid.runner",frames[i],format='runner',append=True)
                addedstr += 1
                print('added +1',addedstr,uuid,idx,i)
                args.pick_uuid.remove(uuid)
                try:
                    args.pick_uuid.remove(uuid)
                except ValueError:
                    pass
                try:
                    args.pick_uuid.remove(uuid)
                except ValueError:
                    pass
                if idx == 175: #range(structures_to_calc)[-1]:
                    for i in args.pick_uuid:
                        sys.stdout.write(i+" ")
                    sys.exit()
                continue
            if calc_analysis: ### analysis stuff
                ana_atoms[idx] = ana_atoms_
                for_DFTmax[idx] = for_DFTmax_
                ana_vol[idx] = frames[i].get_volume()
                ana_vol_pa[idx] = frames[i].get_volume()/frames[i].get_number_of_atoms()
                #print('analys')
                #import my_atom
                ##print('kk',my_atom.get_eqvol("Al"))
                VOL_norm = 0
                for ele in n:
                    nat__ = n[ele]
                    #print('nat__',nat__,'ele',ele)

                    rv__ = ace.pot.reference_volumes[ele]
                    #print('ele',ele,'nat',nat,'rv',ace.pot.reference_volumes[ele])
                    VOL_norm += nat__ * rv__

                #sys.exit('22')
                #VOL_norm = n["Al"]*16.5+n["Mg"]*22.85+n["Si"]*20.5
                #print('ana_vol',ana_vol[idx])
                #print('VOL_norm',VOL_norm)
                VOL_diff = ana_vol[idx] - VOL_norm
                #print('VOL_diff',VOL_diff)
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

                try:
                    ana_mg_conz[idx] = d["Mg"]
                    at_mg = n["Mg"]
                except KeyError:
                    ana_mg_conz[idx] = 0
                    at_mg = 0

                try:
                    ana_si_conz[idx] = d["Si"]
                    at_si = n["Si"]
                except KeyError:
                    ana_si_conz[idx] = 0
                    at_si = 0

                try:
                    ana_al_conz[idx] = d["Al"]
                    at_al = n["Al"]
                except KeyError:
                    ana_al_conz[idx] = 0
                    at_al = 0



            added=""
            if debug:
                print("GG")
                print(frames[idx].get_positions())

            if calc_DFT: ### ene from DFT
                if debug:
                    print("DD before ene_DFT")
                f_tmp = copy.deepcopy(frames[i])
                f_atoms = f_tmp.get_number_of_atoms()
                #print()
                #print('iiiiiiiiidx',idx)
                ene_DFT[idx],ene_DFT_eV_cell[idx] = my.ase_enepot(frames[i],units=args.units,verbose=True)
                #print('DFT_FORCES_EV_ANG')
                #print(DFT_FORCES_EV_ANG[idx])
                #sys.exit('23k4')
                #print('ene_DFTxyz        ',idx,ene_DFT[idx])
                # this is only of interest for thermodynamic integration.
                if verbose > 2: #be_very_verbose:
                    print('ene_DFT        ',idx,ene_DFT[idx])
                    print('ene_DFT_eV_cell[idx]',idx,ene_DFT_eV_cell[idx])
                    print('ene_DFT_eV_cell[0]',idx,ene_DFT_eV_cell[0])
                #print('idxxx',idx,f_atoms)
                if f_atoms > 1:
                    ene_DFT_m0[idx] = (ene_DFT_eV_cell[idx] - ene_DFT_eV_cell[0])/(f_atoms-1.)*1000.
                else:
                    ene_DFT_m0[idx] = (ene_DFT_eV_cell[idx] - ene_DFT_eV_cell[0])/(f_atoms)*1000.
                #sys.exit()
                if verbose > 2: #be_very_verbose:
                    my.show_ase_atoms_content(frames[i],showfirst=3,comment = "STAT2")
                if verbose > 1:
                    print('ene_DFT[idx]     :',ene_DFT[idx],units)

            #print('dd',d)
            #print('n',n)
            #print("ace.pot.atom_energy",ace.pot.atom_energy)
            if len(ace.units.split("_")) == 1: # per structure
                ene_DFT_atomic[idx] = my.get_atomc_energy_from_dicts(n,ace.pot.atom_energy)
            elif len(ace.units.split("_")) == 2: # per atom
                ene_DFT_atomic[idx] = my.get_atomc_energy_from_dicts(d,ace.pot.atom_energy)
            else:
                sys.exit("either per atom or per structure")


            if debug:
                print("BB before ene_DFT_wo_atomic")
            ene_DFT_wo_atomic[idx] = ene_DFT[idx] - ene_DFT_atomic[idx]
            if debug:
                print("AAaa ipi")

            if ipi == True:  ### ene from ipi
                atoms_tmp = copy.deepcopy(frames[i])
                ene_pot_ipi[idx] = my.ipi_ext_calc(atoms_tmp,ace)

            if debug:
                print("AAbb ase")
            if ase == True:  ### ene from ase (without writing lammps files)
                atoms_tmp = copy.deepcopy(frames[i])  # for other instances, since frames change when geoopt
                if ace.geopt == False:
                    if debug:
                        print("AAA")
                    ace.pot_get_and_ase_lmp_cmd()
                    if debug:
                        print("BBB before: ene_pot_ase")
                        #kk = atoms_tmp.get_potential_energy()
                        #print('kk',kk)
                    # check if elements in the structrues are know by the pot
                    if args.remove_unknown_elements_from_structures:
                        my.ase_remove_unknown_species_from_structure(atoms_tmp,ace)
                    my.ase_check_chemical_symbols_agree(atoms_tmp,ace)
                    #print('--3')
                    #for idy,cs in enumerate(atoms_tmp.get_chemical_symbols()):
                    #    print('idy',idy,atoms_tmp.get_chemical_symbols()[idy])
                    #print('--4')
                    f_tmp = copy.deepcopy(atoms_tmp)
                    f_atoms = f_tmp.get_number_of_atoms()
                    ene_pot_ase[idx],ene_pot_eV_cell[idx] = ace.ene(atoms_tmp,debug=debug,return_both=True)
                    POT_ENERGIES_EV_CELL[idx] = ene_pot_eV_cell[idx]
                    #print('att')
                    #print(atoms_tmp.positions)
                    forces_thisstruct = atoms_tmp.get_forces()

                    if args.write_runner_ene_forces:
                        ase_write("runnerFORCES.runner",atoms_tmp,format='runner',append=True)
                    #print('fff')
                    #print(forces_thisstruct)
                    #print('kkk')
                    #print(POT_FORCES_EV_ANG)
                    #print('idx',idx)
                    try:
                        POT_FORCES_EV_ANG[idx] = forces_thisstruct
                    except ValueError:
                        pass
                    #print('compare enes')
                    #print(DFT_ENERGIES_EV_CELL[idx],POT_ENERGIES_EV_CELL[idx])
                    #print('compare forces')
                    #print([DFT_FORCES_EV_ANG[idx]-POT_FORCES_EV_ANG[idx]])
                    #print('78uj00')
                    if f_atoms > 1:
                        ene_pot_m0[idx] = (ene_pot_eV_cell[idx] - ene_pot_eV_cell[0])/(f_atoms-1.)*1000.
                    else:
                        ene_pot_m0[idx] = (ene_pot_eV_cell[idx] - ene_pot_eV_cell[0])/(f_atoms)*1000.
                    dudl_pot_to_DFT[idx]   = ene_DFT_m0[idx] - ene_pot_m0[idx]
                    #print('vv0',ene_pot_eV_cell)
                    #print('vv1',ene_DFT_m0) # first element is 0
                    #print('vv2',ene_pot_m0) # first element is 0
                    #print('vv',dudl_pot_to_DFT)
                    #sys.exit()
                    if idx == 0:
                        dudlav_pot_to_DFT[idx] = 0
                    else:
                        dudlav_pot_to_DFT[idx] = np.mean(dudl_pot_to_DFT[:idx])
                    #print('ene_pot_m0   ',ene_pot_m0[idx])
                    if args.get_harmonic_energy:
                        #print('hm')
                        #print(hesse_matrix)
                        #print('cell')
                        #print(atoms_tmp.cell)
                        ##u = crystal.rcar-crystal0.rcar
                        ###print "uorig:",u
                        #print('pos')
                        #print(atoms_tmp.positions)
                        u = hesse.dpos(atoms_tmp.positions,pos0.positions,atoms_tmp.cell)
                        #print('u')
                        #print(u)

                        ene_har_m0[idx] = hesse.qh_energy_atom(dpos = u, h = hesse_matrix)
                        ene_har_eV_cell[idx] = hesse.qh_energy_cell(dpos = u, h = hesse_matrix)
                        dudl_har_to_DFT[idx]   = ene_DFT_m0[idx] - ene_har_m0[idx]
                        dudl_har_to_pot[idx]   = ene_pot_m0[idx] - ene_har_m0[idx]

                        dudlav_har_to_DFT[idx] = np.mean(dudl_har_to_DFT[:idx])
                        dudlav_har_to_pot[idx] = np.mean(dudl_har_to_pot[:idx])



                    #print('ae',ene_pot_ase[idx])
                    if args.write_forces or args.write_forcesx:
                        if args.write_forcesx:
                            np.savetxt("forcesx.dat",atoms_tmp.get_forces()[:,0])
                        if args.write_forces:
                            np.savetxt("forces.dat",atoms_tmp.get_forces())
                        #ase_write("POSCAR_out",atoms_tmp,format="vasp")
                        #ase_write("POSCAR_out40.runner",atoms_tmp,format="runner")
                        #sys.exit()
                    if debug:
                        print("BBB after: ene_pot_ase")
                    if debug:
                        print("CCC")
                    ###### old atom energies (giulio/daniele)
                    # print('at al',atom_energy_Al,"meV",atom_energy_Al/conv,"hartree") # at al -65558.634493 meV -2.4092354 hartree
                    # print('at mg',atom_energy_Mg,"meV",atom_energy_Mg/conv,"hartree") # at mg -1703619.48966 meV -62.606862 hartree
                    # print('at si',atom_energy_Si,"meV",atom_energy_Si/conv,"hartree") # at si -151289.41503 meV -5.5597835 hartree

                    ###### new atom energies
                    # print('at al',atom_energy_Al,"meV",atom_energy_Al/conv,"hartree") # at al -534123.115151 meV -19.6286626 hartree
                    # print('at mg',atom_energy_Mg,"meV",atom_energy_Mg/conv,"hartree") # at mg -455772.609452 meV -16.7493346 hartree
                    # print('at si',atom_energy_Si,"meV",atom_energy_Si/conv,"hartree") # at si -150410.566175 meV -5.5274864 hartree

                    #if pot == "runner_v2dg" and False: # only in case we load DFT energies from new DFT calcs
                    if "Mg" in ace.pot.atom_energy and "Si" in ace.pot.atom_energy and "Al" in ace.pot.atom_energy and ace.pot.atom_energy["Mg"]/conv < -17.0:
                    #if atom_energy_Mg/conv < -17.0:  # we did load the "old" energies
                        atom_energy_Mg = ace.pot.atom_energy["Mg"]*conv
                        atom_energy_Si = ace.pot.atom_energy["Si"]*conv
                        atom_energy_Al = ace.pot.atom_energy["Al"]*conv
                        #n = my.ase_get_chemical_symbols_to_number_of_species(atoms_tmp[i])
                        ### ene_Al, ene_Mg, ene_Si are per atom, since those are later multi-
                        ### -lied by the number of atoms, this can only work if energies are
                        ### not calculated by _pa
                        new_Si = -5.5274864*conv    # == -150410.566175 (meV)
                        new_Mg = -16.7493346*conv    # == -455772.609452 (meV)
                        new_Al = -19.6286626*conv    # == -534123.115151 (meV)

                        ene_Si = (atom_energy_Si - new_Si )/1000.  # (-151289.41503 - -150410.566175) / 1000. = -0.87884879 eV
                        ene_Mg = (atom_energy_Mg - new_Mg )/1000.
                        ene_Al = (atom_energy_Al - new_Al )/1000.

                        #print('atomSi',atom_energy_Si)  # -151289.41503
                        #print("new_Si",new_Si)

                        #print("ene_Si",ene_Si)  # -0.878848855568 (eV)
                        #print("ene_Mg",ene_Mg)  # -1247.8468802 (eV)
                        #print("ene_Al",ene_Al)  # 468.564480658 (eV)

                        #ene_Al = 468.845752582        # runner - n2p2: -2.4092354 - -19.6286626 = 17.2194272 hartree == 468.56444 eV
                        #ene_Mg = -1247.77831679       # runner - n2p2: -62.6068620 - -16.7493346 =  -45.8575274  hartree == -1247.8468 eV
                        #ene_Si = -0.0161424965998     # runner_v2dg-n2p2_v1ag atomic energy: -5.5597835 - -5.5274864 = -.0322971 hartree == -0.87884879 eV

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
                    ace.pot_get_and_ase_lmp_cmd()
                    ene_pot_ase[idx] = ace.ene(atoms_tmp)
                    #print('b',atoms_tmp.get_positions()[:3])
                    #print('ene',ene_pot_ase[idx])

                    ace.geopt = True
                    ace.pot_get_and_ase_lmp_cmd()
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

            if debug:
                print("AAcc lmp")
            if lmp == True:  ### ene from lammps (by writing lammps files)
                atoms_tmp = copy.deepcopy(frames[i])  # for other instances, since atoms change when geoopt
                ene_pot_lmp_eV_cell,ene_pot_lmp[idx] = my.lammps_ext_calc(atoms_tmp,ace)
                if verbose > 2:
                    print("ENE DFT   :",i,ene_DFT[idx],ace.units)
                    print("ENE ase   :",i,ene_pot_ase[idx],ace.units)
                    print("ENE lammps:",i,ene_pot_lmp[idx],ace.units)
                    print("--------------------------------------")
                ene_diff_lam_ase[idx] = ene_pot_ase[idx] - ene_pot_lmp[idx]
                if verbose > 2:
                    print("DIFF      :",i,ene_diff_lam_ase[idx],ace.units)
                    print("--------------------------------------")
                    print("--------------------------------------")
                    print("MAKE SURE THAT YOU HAVE the correct species in the lammps run!")
                ene_pot[idx] = ene_pot_lmp[idx]

            ### write runner output
            if args.write_runner:
                ase_write("out.runner",frames[i],format='runner',append=True)
            if args.write_espresso_job:
                ase_write("out.espresso-in_"+str(idx),frames[i],format='espresso-in')
                # execute:
                # /Users/glensk/Dropbox/Albert/scripts/dotfiles/scripts/qe-aiida/aiida_submitskripts/createFolders_quantumespresso_from_inputfile.py out.espresso-in -vv
                #scripts = os.environ['SCRIPTS']
                #execu=scripts+'/qe-aiida/aiida_submitskripts/createFolders_quantumespresso_from_inputfile.py out.espresso-in_'+str(idx)



            ### write runner output
            if args.write_runner_repeated:
                ace.get_calculator(frames[i])
                listrepeat = [[1,1,5],[2,3,4], [4,12,3], [13,70,2]]
                listrepeat = [[1,1,5],[2,3,3], [4,12,3], [13,70,3]]
                listrepeat = [[1,1,5],[2,3,4], [4,12,4], [13,70,4]]
                listrepeat = [[1,1,5],[2,3,5], [4,12,5], [13,70,5]]
                atoms_sc, forces_sc, ene_sc_ev, repeat = my.ase_repeat_structure_using_listrepeat(frames[i],ene_DFT[idx],ace.units,listrepeat=listrepeat)
                ase_write(args.inputfile+".repeated",atoms_sc,format='runner',append=True,setforces_ase_units=forces_sc,setenergy_eV=ene_sc_ev)
                ase_write(args.inputfile+".repeated_aiida",atoms_sc,format='espresso-in')

                # statistics
                nat_new = atoms_sc.get_number_of_atoms()
                nat = "ka"
                if nat_new > max_at:
                    max_at = nat_new
                    max_at_orig = nat
                    max_at_id = idx
                if nat_new <= min_at:
                    min_at = nat_new
                    min_at_orig = nat
                    min_at_id = idx
                    print(my.printred("min_at "+str(min_at)+" min_at_orig "+str(min_at_orig)+' (min id '+str(min_at_id)),"this is just an output")

                    # id: 4, 27,161, 188, 243,267,326,362,414,441,468,534,537 (atoms 58)
                print('nat',nat,'(repeat '+str(repeat)+') -> nat',nat*(3**repeat),"min_at",min_at,"min_at_orig",min_at_orig,'(min id '+str(min_at_id)+") max_at",max_at,"max_at_orig",max_at_orig,'(max id'+str(max_at_id)+")")

            ene_pot_wo_atomic[idx] = ene_pot[idx] - ene_DFT_atomic[idx]
            if False: #args.verbose > 1:
                print('9q3nav ene_potxx',ene_pot)
                print('9q3nav ene_DFT_atomicxx',ene_DFT_atomic)
                print('9q3nav ene_pot_wo_atomicxx',ene_pot_wo_atomic)

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
		#print('ace.units',ace.units)
                # enough for meV_pa
                if ace.units == "mev_pa":
                    show = 3
                elif ace.units.lower() == "hartree":
                    show = 6
                elif ace.units.lower() == "ev":
                    show = 6

                fmt_one = '%16.'+str(show)+'f'
                fmt_one = '%13.'+str(show)+'f'
                fmt_after_atms=' '.join([fmt_one]*8)   # add here if a new entry
                ka3="%5.0f %5.0f / %6.0f "+cellshape+" "+fmt_one+" [%4.0f %4.0f %4.0f %4.0f] "+fmt_after_atms+" "+added

                def srl(nr,rnd,l):
                    return str(np.round(nr,rnd)).ljust(l)
                def sr(nr,rnd):
                    return str(np.round(nr,rnd))

                if args.print_only_energies:
                    #__alat_half = crystal.cellvec[0,0]/2
                    #__alat      = crystal.cellvec[0,0]
                    ## e.g.: u_ = (np.array([[4.06,4.06,4.06],[4.06,4.07,4.08],[-0.1,-0.3,3]])+2.025)%4.05-2.025
                    #u_ = (crystal.rcal + __alat_half)%__alat-alat_half

                    if True: #ace.units == 'hartree':
                        aatest = aaatest[i][2]
                        print(str(i).ljust(5),
                                #'ene_DFT:',srl(ene_DFT[idx],3,9),
                                'ene_pot',srl(ene_pot[idx],3,9),
                                'ene_pot_eV_cell',srl(ene_pot_eV_cell[idx],3,9),
                                'ene_har_eV_cell',srl(ene_har_eV_cell[idx]-78.86671025,3,9),
                                #'aatest:',srl(aatest,3,9),
                                ##'ene_DFT_m0' ,srl(ene_DFT_m0[idx],3,9),
                                'ene_har_m0:',srl(ene_har_m0[idx],2,9),
                                'ene_pot_m0' ,srl(ene_pot_m0[idx],3,9),

                                #'dudl_pot_to_DFT',sr(dudl_pot_to_DFT[idx],2),
                                #'dudlav_pot_to_DFT',sr(dudlav_pot_to_DFT[idx],2),

                                #'dudl_har_to_DFT',  sr(dudl_har_to_DFT[idx],2),
                                #'dudlav_har_to_DFT',sr(dudlav_har_to_DFT[idx],2),

                                'dudl_har_to_pot',  sr(dudl_har_to_pot[idx],2),
                                'dudlav_har_to_pot',sr(dudlav_har_to_pot[idx],2),
                                )

                    else:
                        print(str(i).ljust(5),ene_diff_abs[idx])
                else:
                    #ka = (ene_DFT[idx]+537462.536751423380)*33/1000
                    #if ka > 7.5:
                    #    print('idxuu (77)',idx,'uuid:',uuid)
                    if True:
                        if d["Al"] == 1: # == concentration_al
                            ene_pot_wo_atomic[idx] = (ene_pot[idx]+537462.536751423380)*n["Al"]/1000
                        else:
                            ene_pot_wo_atomic[idx] = (ene_DFT[idx]-ene_pot[idx])

                    # interstitials
                    if args.interstitial_form:
                        f= open("interstitials_NN_vs_DFT.txt","a+")
                        #f.write("%3.6f %3.6f\n" % (1,2))
                        f.write("%3.6f %3.6f\n" % ((ene_DFT[idx]+537462.536751423380)*33/1000,(ene_pot[idx]+537462.536751423380)*33/1000))
                        f.close()

                    print(ka3 % (
                    i,
                    idx,
                    structures_to_calc,
                    ene_diff_abs[idx],
                    frames[i].get_number_of_atoms(),
                    at_si,at_mg,at_al,
                    #(ene_DFT[idx]+537462.536751423380)*33/1000,
                    #(ene_pot[idx]+537462.536751423380)*33/1000,


                    #(ene_DFT[idx]+537463)*33/1000,
                    #(ene_pot[idx]+537463)*33/1000,
                    #(ene_DFT[idx]+537460)*33/1000,
                    #(ene_pot[idx]+537460)*33/1000,
                    #(ene_DFT[idx]+537462.536751423380)*33/1000,
                    #(ene_pot[idx]+3579.989020373443)*33/1000, # zhou
                    #(ene_pot[idx]+2380.290192865874)*33/1000,
                    ene_DFT[idx],
                    ene_pot[idx],
                    ene_pot_wo_atomic[idx],
                    for_DFTmax[idx],
                    ene_pot_ase[idx]-ene_pot_ase_geop[idx],
                    ana_vol_pa[idx],
                    ana_dist_min[idx],
                    ana_VOL_diff_norm[idx]),"ABC")
                    #print('POT_FORCES_EV_ANG')
                    #print(POT_FORCES_EV_ANG)
                    #print('sss',
                    if ene_diff_abs[idx] > 60 and False:
                        print('uuid:',uuid)

                if idx >= 99:
                    #np.savetxt("e_pot_eV_cell",ene_pot_eV_cell)
                    #np.savetxt("e_har_eV_cell",ene_har_eV_cell)
                    #np.savetxt("e_pot_hartree_cell",ene_pot_eV_cell*0.036749322)
                    #np.savetxt("e_har_hartree_cell",ene_har_eV_cell*0.036749322-2.89829817)
                    dudl = fah.get_dudl_per_atom_min1_from_energies_eV_cell(
                            energies_lambda_0_eV_cell=ene_har_eV_cell,
                            energies_lambda_1_eV_cell=ene_pot_eV_cell,
                            number_of_atoms=32,
                            align_lambda_0_first_to_lambda_1_yes_no=True)
                    #print('dudl')
                    #print(dudl)
                    if len(dudl) > 0:
                        dudlav = fah.get_dudlav_from_dudl(dudl)
                        print('dudl')
                        print(dudl)
                        print('dudlav')
                        print(dudlav)
                        fah.get_dudl_from_file_with_energies_lambda_0_1(filepath,number_of_atoms)
                return
            if verbose > 1 or args.uuid:
                print('idxuu (88)',idx,'uuid:',uuid)
            if idx in range(0,structures_to_calc,printevery):
                printhere()
                printed = True

            if verbose > 0 and printed == False:
                printhere()

        if args.write_analysis_full:
            np.savetxt("ene_diff_lam_ase.dat",ene_diff_lam_ase,header=ace.units)

        if args.write_analysis:
            print('get RMSE FORCES ...')
            #print('range',range(structures_to_calc))
            #print('list',[range(structures_to_calc)])
            indexes = range(structures_to_calc)
            sum_atoms = 0
            aa = []
            bb = []
            #print("POT_FORCES_EV_ANG")
            #print(POT_FORCES_EV_ANG)
            for idx,x in enumerate(indexes):
                aa = aa + list(np.abs(DFT_FORCES_EV_ANG[x]-POT_FORCES_EV_ANG[x]).flatten())
                bb = bb + list(np.abs(DFT_ENERGIES_EV_CELL[x]-POT_ENERGIES_EV_CELL[x]).flatten())
                #print(aa)
                #print('--')
                #print('aa.mean()',np.array(aa).mean())
                #mean = np.array(aa).mean()
                #print('aa.mean()',np.abs(np.array(aa)))
                #print('aa.rmse()',np.sqrt(np.mean((np.abs(np.array(aa)))**2)))
            rmse = np.sqrt(np.mean(np.array(aa)**2.))
            std = np.array(aa).std()
            rmse_e = np.sqrt(np.mean(np.array(bb)**2.))
            #print('mean',mean)
            print('rmse for',rmse)
            print('std  for',std)
            print('rmse_e',rmse_e)
            np.savetxt("RMSE_FORCES.dat",np.array([rmse]),fmt='%.5f')
            #np.savetxt("for_POT.dat",DFT_FORCES_EV_ANG[range(structures_to_calc)])
            #np.savetxt("for_diff_DFT_POT.dat",DFT_FORCES_EV_ANG[range(structures_to_calc)])

        if args.write_runner:
            print('out.runner written')

        def mysavetxt(what,name,units,save=False):
            whatout = what[~np.isnan(what)]
            if save:
                np.savetxt(name,whatout,header=units)

            #print("in",what.shape,"out",whatout.shape,name)
            return whatout

        if os.path.isfile('log.lammps'):
            os.remove('log.lammps')

        if args.testkmc or args.testkmc_b or args.testkmc_l:
            ene_std         = mysavetxt(ene_std,kmc_file,units,save=True)


        if args.write_analysis_full or args.write_analysis:
            ene_diff        = mysavetxt(ene_diff,"ene_diff.npy",units,save=True)
            ene_diff_abs    = mysavetxt(ene_diff_abs,"ene_diff_abs.npy",units,save=True)
            ene_std         = mysavetxt(ene_std,"ene_std_rmse.npy",units,save=True)
            ene_DFT_wo_atomic = mysavetxt(ene_DFT_wo_atomic,"ene_DFT_wo_atomic",units,save=True)
            ene_DFT_wo_atomic_form = mysavetxt(ene_DFT_wo_atomic-17736.263712621,"ene_DFT_wo_atomic_form",units,save=True)
            ene_pot_wo_atomic = mysavetxt(ene_pot_wo_atomic,"ene_pot_wo_atomic",units,save=True)
            ene_pot_wo_atomic_form = mysavetxt(ene_pot_wo_atomic-17736.263712621,"ene_pot_wo_atomic_form",units,save=True)

            if args.write_analysis_full:
                ene_DFT         = mysavetxt(ene_DFT,"ene_DFT.npy",units,save=True)
                ene_pot         = mysavetxt(ene_pot,"ene_pot.npy",units,save=True)

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
                    np.savetxt("ene_std_rmse.npy",ene_std,header=units)

                    ene_all = np.transpose([range(len(ene_DFT)),ene_DFT,ene_pot,ene_diff_abs,ene_std])
                    ### write analyze.csv
                    try:
                        np.savetxt("ene_all.npy",ene_all,header=units+"\n"+"DFT\t\t"+ase.pot+"\t|diff|\t\t<|diff|>",fmt=' '.join(['%i'] + ['%.10e']*(ene_all.shape[1]-1)))
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


    if args.check_inputdata or args.check_traindata or args.check_testdata or args.check_kmc57data:
        if os.path.isdir(ace.pot.pot_tmpdir):
            shutil.rmtree(ace.pot.pot_tmpdir)
    return


def printhead(structures_to_calc,ace_units):
    #print('structures_to_calc[:3]:',range(structures_to_calc)[:3],'...',range(structures_to_calc)[-3:])
    print()
    print('#                         ('+ace_units+')                                ('+ace_units+')       ('+ace_units+')                  ('+ace_units+')     (eV/Ang)')
    print('#                         (DFT-ref)                                                    ene_pot_wo_..      forces      (if geopt)     Vol per atom    min at dist   VOL_diff_norm')
    print('#   i   idx /    from       diff  [atms   Si   Mg   Al]     ene_DFT       ene_pot      .._wo_atomic       DFTmax       E-E_geopt        (Ang^3)      (Angstrom)     (Ang^3)')
    print('--------------------------------------------------------------------------------------------------------------------------------------------------------------')
    return

def print_compare_ene_vs_DFT(text,pot_ene,DFT_ene="--",eos=False,f300=False,check="",verbose=False,formula_unit=1.):
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
        DFT_ene_out = "%.4f" % round(DFT_ene*formula_unit*conversion,4)
        diff        = "%.2f" % round(np.abs(np.abs(pot_ene/DFT_ene)-1)*conversion,2)

    pot_ene_out = "%.4f" % round(pot_ene*formula_unit*conversion,4)
    #print('p',pot_ene_out,conversion,type(pot_ene_out),type(conversion))
    #print('p',pot_ene_out,conversion,"-->",pot_ene_out*conversion)
    print(text.ljust(45)+": NN_ene_form_eV_per_atom:",
            str(pot_ene_out).ljust(9),
            units.ljust(3)+" DFT_ene_form_eV_per_atom:",
            str(DFT_ene_out).ljust(9))
    if verbose:
        print(
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
    e_al108             = ace.ene(ase_read(sisivac+'/al108'+ff))
    e_al107vac1         = ace.ene(ase_read(sisivac+'/al107vac1'+ff))
    e_al107si1          = ace.ene(ase_read(sisivac+'/al107si1'+ff))
    e_al106si2          = ace.ene(ase_read(sisivac+'/al106si2'+ff))
    e_al105si2va1       = ace.ene(ase_read(sisivac+'/al105si2va1'+ff))
    e_ss_al             = e_al108/108.
    e_ss_va             = e_al107vac1 - e_al108/108. * 107.
    e_ss_one_si         = e_al107si1 - e_al108/108. * 107.
    e_si_si_vac_complex = e_al105si2va1 - 2.*e_ss_one_si - e_ss_va - 105.*e_ss_al

    #p0rint_compare_ene_vs_DFT('vacancy_formation (unrelaxed) NN: ',e_ss_va,
    print_compare_ene_vs_DFT("vacancy_formation (unrelaxed, oDFT)",e_ss_va,0.69)
    print_compare_ene_vs_DFT("si-si-vac-complex (unrelaxed, oDFT)",e_si_si_vac_complex,0.074)
    print()
    return

def show_energy_diff_DFT_NN(struct,eDFT,e,text,units="eV"):
    ''' show_energy_diff_DFT_NN(struct_pure_al,ace.eDFT_pure_al,ace.e_pure_al,"pure al",units="eV") '''
    nat = struct.get_number_of_atoms()
    print(text.ljust(20," "),str(eDFT).ljust(15),"eV",str(e).ljust(15),"eV","diff",str(eDFT-e).ljust(17),"eV/cell","diff mev_pa",str((eDFT-e)/nat*1000).ljust(15),"meV_pa")
    return


def get_dilute_formation_energy(matrix_element="Al",text="dilute formation energy supercell",sc="all",nsi=1,nmg=0,nvac=0,e_si_diamond_pa=0,ace=False,t2="",verbose=False):
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

        print('############## get the dilute frame Al+{Si,Mg,Vac} (T=0K ene_al_xx)')
        frame_path = ace.savefolder+"frame_solute_"+matrix_element+str(nat-1)+solute_element+"1.runner"
        print('fp',frame_path)
        if os.path.isfile(frame_path):
            if ace.verbose:
                print('read frame path:',frame_path)
            frame_al_xx = ase_read(frame_path,format="runner")
        else:
            frame_al_xx = my.get_ase_atoms_object_kmc_al_si_mg_vac(matrix_element=matrix_element,ncell=sc,nsi=nsi,nmg=nmg,nvac=nvac,a0=(vpa*4.)**(1./3.),cell='cubic')
            # here, we want to work in he constant pressure approach
            # for this we need to relax atomic positions && volume
            ace.ase_relax_cellshape_volume_positions(frame_al_xx,verbose=False)
            ase_write(frame_path,frame_al_xx,format="runner")



        ############## get the dilute_formation @T=0K && ace.eform_dilute_xx_
        ############## get the dilute_formation @T=0K && ace.eform_dilute_xx_
        ############## get the dilute_formation @T=0K && ace.eform_dilute_xx_
        ene_al_xx             = ace.ene(frame_al_xx) # this is already relaxed
        dilute_formation      = ene_al_xx  - (sc**3.*4. -1.) * ace.al_fcc_ene_pa
        ace.E_SS_Al           = ace.al_fcc_ene_pa
        #print('XX ace.E_SS_Al',ace.E_SS_Al)
        #sys.exit('5555555')
        if solute_element == "Si":
            ace.E_SS_Si              = dilute_formation
            ace.struct_dilute_si_NN  = frame_al_xx
            print('XX ace.E_SS_Si            :',ace.E_SS_Si)
            print('XX ace.struct_dilute_si_NN:',ace.struct_dilute_si_NN,'in sc:',sc)
            print('stress',ace.stress(frame_al_xx))
        if solute_element == "Mg":
            ace.E_SS_Mg              = dilute_formation
            ace.struct_dilute_mg_NN  = frame_al_xx
            print('XX ace.E_SS_Mg:',ace.E_SS_Mg)
            print('XX ace.struct_dilute_mg_NN:',ace.struct_dilute_mg_NN,'in sc:',sc)
            print('stress',ace.stress(frame_al_xx))
        if solute_element == "Vac":
            ace.E_SS_vac             = dilute_formation
            ace.struct_vac_NN        = frame_al_xx
            print('XX ace.E_SS_vac           :',ace.E_SS_vac)
           # print('XX ace.struct_vac_NN      :',ace.struct_dilute_vac_NN,'in sc:',sc)


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
        #print('could not be loaded... '+solute_element)

        ### bulk free energy
        ### bulk free energy
        if False: # if to calculate T>0K
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
    print('# get_al_fcc_equilibrium ...')
    frame_al = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=1,nsi=0,nmg=0,nvac=0,a0=4.045,cell='cubic', create_fake_vacancy=False,matrix_element="Al",crystal_structure="fcc")
    ace.ase_relax_cellshape_and_volume_only(frame_al,verbose=False)
    ace.al_fcc = frame_al
    ace.al_fcc_ene_pa = ace.ene(frame_al)/frame_al.get_number_of_atoms()
    ace.al_fcc_vol_pa = frame_al.get_volume()/frame_al.get_number_of_atoms()
    if ace.verbose:
        print("NN Al vpa @T=0K",ace.al_fcc_vol_pa,"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    #print('e_ and f_ al should all be done consistently from the NN!!!')
    if False:
        # I dont think I need this
        ace.get_elastic_external(atomsin=ace.al_fcc,verbose=False,text="Al_fcc bulk 4at",get_all_constants=False)
    return

def get_mg_hcp_equilibrium(ace):
    print('# get_mg_hcp_equilibrium ...')
    frame_mg = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=1,nsi=0,nmg=0,nvac=0,a0=0,cell='cubic',create_fake_vacancy=False,matrix_element="Mg",crystal_structure="hcp")
    ace.ase_relax_cellshape_and_volume_only(frame_mg,verbose=False)
    ace.mg_hcp        = frame_mg
    ace.mg_hcp_ene_pa = ace.ene(frame_mg)/frame_mg.get_number_of_atoms()
    ace.mg_hcp_vol_pa = frame_mg.get_volume()/frame_mg.get_number_of_atoms()
    # seems also to make the c/a to relax ... good
    #if ace.verbose:
    #    print("NN Mg vpa @T=0K",ace.mg_hcp_vol_pa)
    #    print('frame_mg.cell',frame_mg.cell)
    return

def get_si_dc_equilibrium(ace):
    print('# get_si_dc_equilibrium(ace)...')
    frame_si = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=1,nsi=0,nmg=0,nvac=0,a0=0,cell='cubic',create_fake_vacancy=False,matrix_element="Si",crystal_structure="dc",)
    ace.ase_relax_cellshape_and_volume_only(frame_si,verbose=False)
    ace.si_dc = frame_si
    #print('frame_si.c',frame_si.positions)
    #print('frame_si.p',frame_si.cell)
    #print('at',frame_si.get_chemical_symbols())
    #print('vol',frame_si.get_volume()/frame_si.get_number_of_atoms())
    ace.si_dc_ene_pa = ace.ene(frame_si)/frame_si.get_number_of_atoms()
    ace.si_dc_vol_pa = frame_si.get_volume()/frame_si.get_number_of_atoms()
    return

def get_basic_NN_energies_ace(ace):
    print('######## get_basic_NN_energies_ace #############')
    scripts = os.environ['scripts']
    tests = scripts+'/tests/'
    ace.savefolder = tests+'/Al-Mg-Si/save_'+ace.pot.pot+"/"
    if not os.path.isdir(ace.savefolder):
        my.mkdir(ace.savefolder)

    #ase_write('pos_si.lmp',frame_si,format='lammpsrunner')
    get_al_fcc_equilibrium(ace)
    get_mg_hcp_equilibrium(ace)
    get_si_dc_equilibrium(ace)
    ace.free_ene_formation_dilute_mg_ = False
    ace.free_ene_formation_dilute_si_ = False
    #ene_pot_lmp = my.lammps_ext_calc(ace.al_fcc,ace)
    #ase_write('pos_al.lmp',ace.al_fcc,format='lammpsrunner')
    #get_vpa(
    if ace.verbose:
        print("NN 1 Al vpa @T=0K",my.ase_vpa(ace.al_fcc),"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    #ace.check_frame('bbb',ace.al_fcc)
    #ace.ase_relax_cellshape_and_volume_only(ace.al_fcc,verbose=False)
    #ace.check_frame('aaa',ace.al_fcc)
    if ace.verbose:
        print("NN 2 Al vpa @T=0K",my.ase_vpa(ace.al_fcc),"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    ace.get_murn(ace.al_fcc,verbose=False,return_minimum_volume_frame = False, atomrelax=False)
    if ace.verbose:
        print("NN 3 Al vpa @T=0K",my.ase_vpa(ace.al_fcc),"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    #ace.ase_relax_cellshape_and_volume_only(ace.al_fcc,verbose=False)
    #ace.check_frame('ccc',ace.al_fcc)
    #sys.exit()
    if ace.verbose:
        print("NN 4 Al vpa @T=0K",my.ase_vpa(ace.al_fcc)) #,"(==alat)",(ace.al_fcc_vol_pa*4.)**(1./3.))
    #sys.exit('ace fcc al si dc')
        print('before too lon 12')
    sc = 3;
    get_dilute_formation_energy(text="NN dilute formation energy Si ",sc=sc,nsi=1,nmg=0,e_si_diamond_pa=ace.si_dc_ene_pa,ace=ace,t2="Kobayashi 0.375 eV")
    #sys.exit('too lon 12')
    get_dilute_formation_energy(text="NN dilute formation energy Mg ",sc=sc,nsi=0,nmg=1,e_si_diamond_pa=ace.mg_hcp_ene_pa,ace=ace,t2="Kobayashi 0.090 eV")
    get_dilute_formation_energy(text="NN dilute formation energy Vac",sc=sc,nsi=0,nmg=0,nvac=1,e_si_diamond_pa=0.,ace=ace,t2="Kobayashi 0.654 eV")
    print()
    return

def get_dilute_si_mg_f(ace):
    print("######## get_dilute_si_mg_f #############")

    scripts = os.environ['scripts']
    tests = scripts+'/tests/'
    bprime = tests+'/Al-Mg-Si/Mg9Si5_beta_prime/dilute_structures_and_beta_prime.input.data'


    frames = ase_read(bprime,index=":",format="runner")

    struct_bprime    = frames[0]
    struct_pure_al   = frames[3]
    struct_dilute_mg = frames[1]
    struct_dilute_si = frames[2]

    # @ 0K (eV/cell)
    ace.eDFT_bprim,kb_      = my.ase_enepot(struct_bprime   ,units=ace.units)
    ace.eDFT_pure_al,kb_    = my.ase_enepot(struct_pure_al  ,units=ace.units) # 108 atoms
    ace.eDFT_dilute_mg,kb_  = my.ase_enepot(struct_dilute_mg,units=ace.units) # 108 atoms
    ace.eDFT_dilute_si,kb_  = my.ase_enepot(struct_dilute_si,units=ace.units) # 108 atoms
    sc = 3

    if False:
        print('struct_dilute_si DFT')
        ka = struct_dilute_si
        for idx,i in enumerate(ka.positions):
            print(idx,ka.get_chemical_symbols()[idx],ka.positions[idx])

        print('... get struct dilute si for the NN')
    #struct_dilute_si_NN = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=sc,nsi=0,nmg=0,nvac=0,matrix_element=["Al","Si"], a0=False,cubic=True,create_fake_vacancy=False,crystal_structure="fcc_dilute",verbose=False)
    #ace.ase_relax_atomic_positions_only(struct_dilute_si_NN,fmax=0.0001,verbose=False)
    #if False:
    #    print('... get struct dilute si has ene',my.ase_enepot(struct_dilute_si_NN,units=ace.units))
    #    print()
    #    print('... get struct dilute mg')



    #struct_dilute_mg_NN = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=sc,nsi=0,nmg=0,nvac=0,matrix_element=["Al","Mg"], a0=False,cubic=True,create_fake_vacancy=False,crystal_structure="fcc_dilute",verbose=False)
    #ace.ase_relax_atomic_positions_only(struct_dilute_mg_NN,fmax=0.0001,verbose=False)
    #if False:
    #    print('... get struct dilute mg has ene',my.ase_enepot(struct_dilute_mg_NN,units=ace.units))

    #    print('... done')

    #ace.e_bprime        = ace.ene(struct_bprime)
    ace.e_pure_al       = ace.al_fcc_ene_pa*108   # 108 atoms
    ace.e_dilute_mg     = ace.ene(ace.struct_dilute_mg_NN) # 108 atoms
    ace.e_dilute_si     = ace.ene(ace.struct_dilute_si_NN) # 108 atoms
    #if False:
    #    ka = struct_dilute_si_NN
    #    print('struct_dilute_si_NN')
    #    for idx,i in enumerate(ka.positions):
    #        print(idx,ka.get_chemical_symbols()[idx],ka.positions[idx])
    #nat = struct_bprime.get_number_of_atoms()
    #print('ace.eDFT_bprim  (',nat,'atoms)', ace.eDFT_bprim/nat,"(per atom)")
    nat = struct_pure_al.get_number_of_atoms()
    print('ace.eDFT_pure_al (',nat,'atoms)', ace.eDFT_pure_al/nat,"(per atom)")
    nat = struct_dilute_mg.get_number_of_atoms()
    print('ace.eDFT_dilute_mg  (',nat,'atoms)', ace.eDFT_dilute_mg/nat,"(per atom)")
    nat = struct_dilute_si.get_number_of_atoms()
    print('ace.eDFT_dilute_si (',nat,'atoms)',ace.eDFT_dilute_si/nat,"(per atom)")
    print()
    print('ace.e_pure_al   (',108,'atoms)',ace.al_fcc_ene_pa,"(per atom)")
    print('ace.e_dilute_si (',108,'atoms)',ace.e_dilute_si,"(per cell)")
    print('stress disule_si',ace.stress(ace.struct_dilute_si_NN))
    print('stress disule_mg',ace.stress(ace.struct_dilute_mg_NN))


    print()

    print('ace.E_SS_Al :',ace.E_SS_Al,"energy of one Al")
    print('ace.E_SS_Mg :',ace.E_SS_Mg,"one Mg is left over")
    print('ace.E_SS_Si :',ace.E_SS_Si,"one Si is left over")
    print('ace.E_SS_vac:',ace.E_SS_vac,"vacancy formation energy")

    print("                         ","DFT(eV)          NN(eV)                 eV all atoms")
    print("                         ","per cell         per celper cell        eV all atoms")
    show_energy_diff_DFT_NN(struct_pure_al,ace.eDFT_pure_al,ace.e_pure_al,"pure al",units="eV")
    show_energy_diff_DFT_NN(struct_dilute_mg,ace.eDFT_dilute_mg,ace.e_dilute_mg,"dilute mg",units="eV")
    show_energy_diff_DFT_NN(struct_dilute_si,ace.eDFT_dilute_si,ace.e_dilute_si,"dilute si",units="eV")
    #show_energy_diff_DFT_NN(struct_bprime,ace.eDFT_bprim,ace.e_bprime,"bprime",units="eV")
    print()
    # (eV/defect)
    ace.eform_dilute_al     =  ace.al_fcc_ene_pa
    ace.eform_dilute_si     =  ace.e_dilute_si    - 107*ace.al_fcc_ene_pa # eV per defect
    print('ace.e_dilute_si          (eV):',ace.e_dilute_si,'or atoms:',struct_dilute_si.get_number_of_atoms())
    print('ace.e_pure_al (per atom) (eV):',ace.al_fcc_ene_pa)
    print('ace.eform_dilute_si     =  ace.e_dilute_si    - 107*ace.e_pure_al/108 =',ace.eform_dilute_si,'eV')
    ace.eform_dilute_mg     =  ace.e_dilute_mg    - 107*ace.al_fcc_ene_pa       # eV per defect

    # DFT
    ace.fDFT_dilute_al      =  ace.eDFT_pure_al/108    # eV per defect
    ace.fDFT_dilute_si      =  ace.eDFT_dilute_si - 107*ace.eDFT_pure_al/108    # eV per defect
    ace.fDFT_dilute_mg      =  ace.eDFT_dilute_mg - 107*ace.eDFT_pure_al/108    # eV per defect

    print('solid solution energies, this is an energy for one atom of the specific kind')
    print('ace.eform_dilute_{si,mg},ace.fDFT_dilute_{si,mg}')
    print('in Junge2017 (== doi: 10.1016/j.actamat.2017.08.017):')
    print('E^{SS}_{X}')
    print_compare_ene_vs_DFT('formation dilute si (1Si in bulk Al) 3x3x3',ace.eform_dilute_si,ace.fDFT_dilute_si, '-',"-")
    print_compare_ene_vs_DFT('formation dilute mg (1Mg in bulk Al) 3x3x3',ace.eform_dilute_mg,ace.fDFT_dilute_mg, '-',"-")
    print_compare_ene_vs_DFT('formation dilute al (       bulk Al) 3x3x3',ace.eform_dilute_al,ace.fDFT_dilute_al, '-',"-")
    print(',ace.eform_dilute_si',ace.eform_dilute_si)
    print()
    print('formation energies of the dilute systems (all atoms are acounted for)')
    print('ace.eform_dilute_{si,mg}-ace.{si,mg}_dc_ene_pa')
    print('this is not a quantity that is saved in any variable')
    # dilute_formation = ace.eform_dilute_{si,mg}_ = ene_al_xx  - (sc**3.*4. -1.) * ace.al_fcc_ene_pa
    print('ace.eform_dilute_si',ace.eform_dilute_si)
    print('ace.si_dc_ene_pa   ',ace.si_dc_ene_pa)
    print_compare_ene_vs_DFT('formation_energy dilute si (1Si in bulk Al) 3x3x3',ace.eform_dilute_si-ace.si_dc_ene_pa ,"-", '-',"-")
    print_compare_ene_vs_DFT('formation_energy dilute mg (1Mg in bulk Al) 3x3x3',ace.eform_dilute_mg-ace.mg_hcp_ene_pa,"-", '-',"-")
    print_compare_ene_vs_DFT('formation_energy dilute si (1Si in bulk Al) 4x4x4',ace.E_SS_Si-ace.si_dc_ene_pa ,"-", '-',"-")
    print_compare_ene_vs_DFT('formation_energy dilute mg (1Mg in bulk Al) 4x4x4',ace.E_SS_Mg-ace.mg_hcp_ene_pa,"-", '-',"-")
    np.savetxt("heat_of_formation_si_4x4x4.dat",np.array([ace.E_SS_Si-ace.si_dc_ene_pa]),fmt='%1.4f')
    np.savetxt("heat_of_formation_mg_4x4x4.dat",np.array([ace.E_SS_Mg-ace.mg_hcp_ene_pa]),fmt='%1.4f')
    np.savetxt("heat_of_formation_vac_3x3x3.dat",np.array([ace.E_SS_vac]),fmt='%1.4f')
    print()

    #print('ace.si_dc_ene_pa',ace.si_dc_ene_pa)

    if False:
        ka =  ace.e_dilute_si    - 107*ace.al_fcc_ene_pa  - ace.si_dc_ene_pa     # eV per defect
        print('ka',ka,'eV')
        # for one si in al: from vasp PBE energies
        #  ene 31al1si/32          ene1al                 ene1si
        # -3784.479489982940*32-(-3745.245421433458*31)-(-5424.95)
        # from NN:
        #  ene 107al1si
        # -57663.3933854  - (-58045.8813009/108*107) - (ace.si_dc_ene_pa)
        kb = ace.e_dilute_si - (107*ace.al_fcc_ene_pa) - ace.si_dc_ene_pa
        print('kb',kb)
        print('ace.e_dilute_si         (eV)',ace.e_dilute_si,"?? -57663.3933854 eV (in3x3x3cell)")
        print('ace.si_dc_ene_pa        (eV)',ace.si_dc_ene_pa,"?? -155.368143627 eV")
        print("ace.e_pure_al/108       (eV)",ace.e_pure_al/108,"?? -537461.99 eV")
        print("(107*ace.e_pure_al/108) (eV)",(107*ace.e_pure_al/108),"??, -57508.4338719 eV")
    #
    #   -533920.334868975333*108  -(107*-537461.993661476416) -(-155368.146177827992) = 405 meV ! correct formation energy one si in al
    # = -57663396.1658493         -57508433.321778
    # NN on DFT pos:
    #   -57663.3933854 (eV)       -57508.419437 (eV)
    #   -57663.3933854  - -57508.419437 --155.36 = 0.386
    # @ 300K ( for 1 atom, therefore *108 (to get per cell), in meV, therefore /1000)
    if False: # if T>0K
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


def get_formation_energy(ace,frame,text,atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=False,try_harmonic_readfile=False,debug=False,formula_unit=1.):
    ''' Bill sais that T=443 is the relevant temperature '''
    energy_NN_cell_unrelaxed_or_DFTrelaxed = ace.ene(frame.copy())  # needs a copy here, otherwise the DFT energy is evaluated and not the NN energy
    #print('energy_NN_cell_unrelaxed_or_DFTrelaxed (1)',energy_NN_cell_unrelaxed_or_DFTrelaxed)

    d = my.ase_get_chemical_symbols_to_number_of_species(frame,known_elements_by_pot=ace.pot.elements)
    #print('d',d)
    conz1 = d["Mg"]
    conz2 = (d["Mg"]+d["Si"])
    conz = np.float(d["Mg"])/np.float((d["Mg"]+d["Si"]))
    nat = frame.get_number_of_atoms()
    print('d[Mg]:',d["Mg"],'d[Si]:',d["Si"],'d[Al]:',d["Al"],'number_of_atoms=nat:',nat,'formula_unit',formula_unit)
    print('ace.fDFT_dilute_mg',ace.fDFT_dilute_mg,"ace.E_SS_Mg",ace.E_SS_Mg)
    print('ace.fDFT_dilute_si',ace.fDFT_dilute_si,"ace.E_SS_Si",ace.E_SS_Si)
    print('ace.fDFT_dilute_al',ace.fDFT_dilute_al,"ace.E_SS_Al",ace.E_SS_Al)
    print('energy_NN_cell_unrelaxed_or_DFTrelaxed',energy_NN_cell_unrelaxed_or_DFTrelaxed)

    heat_precip_T0K_DFT = "-"
    eDFT = ""
    if DFT_ene != False:
        eDFT,kb_   		= my.ase_enepot(frame  ,units=ace.units)
        print('eDFT                                  ',eDFT)
        e           = ace.ene(frame)
        #print('energy_DFT',eDFT)
        ediff_ev  	= energy_NN_cell_unrelaxed_or_DFTrelaxed - eDFT
        ediff_mev_pa 	= ediff_ev*1000./nat
        #print("##",text,'e - eDFT:',round(ediff_mev_pa,3),"meV/pa")
        heat_precip_T0K_DFT = (eDFT - d["Mg"]*ace.fDFT_dilute_mg - d["Si"]*ace.fDFT_dilute_si - d["Al"]*ace.fDFT_dilute_al)/nat
        if ace.verbose:
            show_energy_diff_DFT_NN(frame,eDFT,e,text,units="eV")
            print(my.printred("DFT")+" energy precipitate (eV)",eDFT,"containing:",d)
            print(my.printred("DFT")+" energy eform_dilute_mg (eV)",ace.fDFT_dilute_mg,"times",d["Mg"])
            print(my.printred("DFT")+" energy eform_dilute_si (eV)",ace.fDFT_dilute_si,"times",d["Si"])
            print(my.printred("DFT")+" energy eform_dilute_al (eV)",ace.fDFT_dilute_al,"times",d["Al"])
            print(my.printred("DFT: divide everything by       "+str(nat)+" to get to the DFT formation energy of "+str(heat_precip_T0K_DFT)+" eV"))

        # write DFT@DFT
        if ace.verbose:
            print('appending to',ace.written_summary[0])
        f=open(ace.written_summary[0], "a+")
        f.write(str(conz)+"   "+str(heat_precip_T0K_DFT)+" "+text.replace(" ", "_")+"\n")
        f.close()

        if formula_unit != False:
            eform = heat_precip_T0K_DFT*formula_unit
            #print('eform DFT per formula unit:',np.round(eform,3))
            f=open(ace.written_summary[0]+"_per_formula_unit.dat", "a+")
            print('wirttein to',ace.written_summary[0]+"_per_formula_unit.dat")
            f.write(str(conz)+"   "+str(eform)+" "+text.replace(" ", "_")+"\n")
            f.close()
    print('bef el')
    get_elastic_constants_al_ext(ace,frame)
    print('aft el')

    # write NN@DFT
    print('-->>',energy_NN_cell_unrelaxed_or_DFTrelaxed - d["Mg"]*ace.E_SS_Mg - d["Si"]*ace.E_SS_Si - d["Al"]*ace.E_SS_Al)
    heat_precip_T0K_NN_unrelaxed  = (energy_NN_cell_unrelaxed_or_DFTrelaxed - d["Mg"]*ace.E_SS_Mg - d["Si"]*ace.E_SS_Si - d["Al"]*ace.E_SS_Al)/nat

    if ace.verbose:
        show_energy_diff_DFT_NN(frame,eDFT,e,text,units="eV")
        print(my.printred("NN@DFT")+" energy precipitate (eV)",energy_NN_cell_unrelaxed_or_DFTrelaxed,"containing:",d)
        print(my.printred("NN@DFT")+" energy eform_dilute_mg (eV)",ace.E_SS_Mg,"times",d["Mg"])
        print(my.printred("NN@DFT")+" energy eform_dilute_si (eV)",ace.E_SS_Si,"times",d["Si"])
        print(my.printred("NN@DFT")+" energy eform_dilute_al (eV)",ace.E_SS_Al,"times",d["Al"])
        print(my.printred("NN@DFT: divide everything by       "+str(nat)+" to get to the NN@DFT formation energy of "+str(heat_precip_T0K_NN_unrelaxed)+" eV"))

    print_compare_ene_vs_DFT(my.printred("@fixed   pos: NN "+text+" @0K"),heat_precip_T0K_NN_unrelaxed,heat_precip_T0K_DFT,False,"-",check="",verbose=ace.verbose,formula_unit=formula_unit)
    #print('heat_precip_T0K_NN_unrelaxed',heat_precip_T0K_NN_unrelaxed)
    if ace.verbose:
        print('appending to',ace.written_summary[1])
    f=open(ace.written_summary[1], "a+")
    f.write(str(conz)+"   "+str(heat_precip_T0K_NN_unrelaxed)+" "+text.replace(" ", "_")+"\n")
    f.close()

    if formula_unit != False:
        eform = heat_precip_T0K_NN_unrelaxed*formula_unit
        #print('eform NN per formula unit:',np.round(eform,3))
        f=open(ace.written_summary[1]+"_per_formula_unit.dat", "a+")
        print('wirttein to',ace.written_summary[1]+"_per_formula_unit.dat")
        f.write(str(conz)+"   "+str(eform)+" "+text.replace(" ", "_")+"\n")
        f.close()

    # @ T=0K relax
    if True:
        if atomrelax:   ace.ase_relax_atomic_positions_only(frame)
        if volumerelax: ace.ase_relax_cellshape_and_volume_only(frame)
        if atomrelax:   ace.ase_relax_atomic_positions_only(frame)
        if volumerelax: ace.ase_relax_cellshape_and_volume_only(frame)
        check = ace.check_frame('',frame=frame,verbose=False)
        #print('11 atomrelax',atomrelax,'volumerelax',volumerelax,"check",check)
        ace.get_murn(frame,verbose=False,return_minimum_volume_frame = volumerelax, atomrelax=atomrelax,printminimal=False)
        energy_NN_cell_NNrelaxed = ace.ene(frame) #,atomrelax=atomrelax,cellrelax=cellrelax)
        #print('d mg:',d["Mg"],'d si:',d["Si"],'d al:',d["Al"])
        #heat_precip_T0K         = (e - d["Mg"]*ace.eform_dilute_mg - d["Si"]*ace.eform_dilute_si - d["Al"]*ace.eform_dilute_al)/nat
        if ace.verbose:
            print('energy_NN_cell_NNrelaxed',energy_NN_cell_NNrelaxed)
            print(my.printred('NN ace.E_SS_Mg'),ace.E_SS_Mg)
            print(my.printred('NN ace.E_SS_Si'),ace.E_SS_Si)
        heat_precip_T0K         = (energy_NN_cell_NNrelaxed - d["Mg"]*ace.E_SS_Mg - d["Si"]*ace.E_SS_Si - d["Al"]*ace.E_SS_Al)/nat

        if ace.verbose:
            show_energy_diff_DFT_NN(frame,eDFT,e,text,units="eV")
            print(my.printred("NN@NNrelaxed")+" energy precipitate (eV)",energy_NN_cell_NNrelaxed,"containing:",d)
            print(my.printred("NN@NNrelaxed")+" energy eform_dilute_mg (eV)",ace.E_SS_Mg,"times",d["Mg"])
            print(my.printred("NN@NNrelaxed")+" energy eform_dilute_si (eV)",ace.E_SS_Si,"times",d["Si"])
            print(my.printred("NN@NNrelaxed")+" energy eform_dilute_al (eV)",ace.E_SS_Al,"times",d["Al"])
            print(my.printred("NN@NNrelaxed: divide everything by       "+str(nat)+" to get to the NN@NNrelaxed formation energy of "+str(heat_precip_T0K)+" eV"))


        print_compare_ene_vs_DFT(my.printred("@fixed   pos: NN "+text+" @0K"),heat_precip_T0K,False,False,"-",check="",verbose=ace.verbose,formula_unit=formula_unit)
        if ace.verbose:
            print('appending to',ace.written_summary[2])
        f=open(ace.written_summary[2], "a+")
        f.write(str(conz)+"   "+str(heat_precip_T0K)+" "+text.replace(" ", "_")+"\n")
        f.close()
        if formula_unit != False:
            eform = heat_precip_T0K*formula_unit
            f=open(ace.written_summary[2]+"_per_formula_unit.dat", "a+")
            print('wirttein to',ace.written_summary[2]+"_per_formula_unit.dat")
            f.write(str(conz)+"   "+str(eform)+" "+text.replace(" ", "_")+"\n")
            f.close()

        check = ace.check_frame('',frame=frame,verbose=False)
        #print('22 atomrelax',atomrelax,'volumerelax',volumerelax,"check",check)
        #print_compare_ene_vs_DFT(my.printred("@relaxed pos: NN "+text+" @0K"),heat_precip_T0K,heat_precip_T0K_DFT,vinet,"-",check=check)
        #if ace.verbose:
        #    print("NN energy structure/precipitate (eV)",e)
        #    print("NN energy eform_dilute_mg (eV)",ace.eform_dilute_mg,"times",d["Mg"])
        #    print("NN energy eform_dilute_si (eV)",ace.eform_dilute_si,"times",d["Si"])
        #    print("divide everything by       ",nat,"to get to the formation energy @T=0K:",heat_precip_T0K_DFT)


        if False:  # if T>0K
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




        eDFT_pure_al,kb_   = my.ase_enepot(struct_pure_al  ,units=ace.units) # 108 atoms
        eDFT_dilute_mg,kb_ = my.ase_enepot(struct_dilute_mg,units=ace.units) # 108 atoms
        eDFT_dilute_si,kb_ = my.ase_enepot(struct_dilute_si,units=ace.units) # 108 atoms
        eDFT_precip,kb_    = my.ase_enepot(struct_mg9si5   ,units=ace.units) # 28 atoms, 10Si, 18Mg

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

    print_compare_ene_vs_DFT(my.printred("Mg9Si5 (DFT@DFT fully relaxed) @0K Vissers GGA"),-0.335,DFT_ene="-",eos=False,f300=False)


    try_read = ace.savefolder+"h_Mg9Si5_at_DFT_relaxed"
    if ace.verbose:
        print('try_read hessematrix',try_read)
    get_formation_energy(ace,frame,"Mg9Si5 (DFT@DFT fully relaxed)",atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=True,try_harmonic_readfile=try_read)

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
    get_formation_energy(ace,frame,"Mg9Si5 (NN@NN  fully relaxed)" ,atomrelax=True,cellrelax=True ,volumerelax=True,try_harmonic_readfile=try_read,debug=False)
    #print(np.round(frame.get_positions(),2))
    #print()
    #print(np.round(frame.get_cell(),2))
    return

#def test_Mg9Si5_pos(ace):
#    path = my.scripts()+'/tests/Al-Mg-Si/Mg9Si5_beta_prime/BetaPrime_structures_relax.input.data.v3'
#    frame = ase_read(path,'0',format="runner")
#    get_formation_energy(ace,frame,"Mg9Si5 (DFT@DFT)",atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=True)
#    return


def test_Mg2Si(ace):
    path = my.scripts()+'/tests/Al-Mg-Si/Mg2Si/POSCAR'
    frame = ase_read(path,format="vasp")
    try_read = ace.savefolder+"h_Mg2Si1_at_NN_relaxed"
    #print('try read',try_read)
    get_formation_energy(ace,frame,"Mg2Si (NN@NN)",atomrelax=True,cellrelax=True,volumerelax=True,try_harmonic_readfile=try_read)
    print('mg2si stress:',frame.get_stress())
    #print('mg2si stress:',frame.get_isotropic_pressure(frame.get_stress()))
    #print('vinet',ace.get_murn(frame),my.ase_vpa(frame))
    print()
    return





def test_antisites(ace):
    print("######## antisites #############")
    print('taken from:')
    print('verdi group list -A')
    print('  83 out_antisites_repMg4Al3Si4_NEW.runner       albert.glensk.gmail.com')
    print('  84 out_antisites_repMg4Al3Si4_NEW.runner_calc  albert.glensk@gmail.com')
    print('  85 out_antisites_repMg5Al2Si4_NEW.runner       albert.glensk@gmail.com')
    print('  86 out_antisites_repMg5Al2Si4_NEW.runner_calc  albert.glensk@gmail.com')
    print('  87 out_antisites_repMg5Si6_NEW.runner          albert.glensk@gmail.com')
    print('  88 out_antisites_repMg5Si6_NEW.runner_calc     albert.glensk@gmail.com')

    #phase =  "Mg4Al3Si4" # done
    #phase =  "Mg5Al2Si4" # done
    #phase = "Mg5Si6" # done


    potpath_using = ace.pot.potpath+"/epoch_"+str(ace.pot.potepoch_using)

    doit = [ "Mg4Al3Si4", "Mg5Al2Si4", "Mg5Si6" ]
    doit = [ "Mg4Al3Si4"] #, "Mg5Al2Si4", "Mg5Si6" ]
    for phase in doit:
        ######################################
        ## get the antisites
        ######################################
        print('------- antisites (read or calculate) ----------')
        read_antisites = "aiida_exported_group_out_antisites_rep"+phase+"_NEW.runner_calc__all_steps.input.data"

        # read_antisites = xxx  (and relax)
        #antisite_relaxed_by_nn, antisite_initial = my.ase_relax_structure_fully_and_save_to_pot(ace,read_antisites)
        antisite_initial = ase_read(readpath,index=whichstruct,format="runner")
        antisite_relaxed = ase_read(readpath,index=whichstruct,format="runner")
        print('type(antisite_relaxed_by_nn:',type(antisite_relaxed_by_nn))
        print('type(antisite_initial_by_nn:',type(antisite_initial))

        ######################################
        ## get the original phase
        ######################################
        print('------- original phase -----')
        read_orig_phase = "aiida_exported_group_NN_relaxed_"+phase+"_n2p2_v2ag_calc__all_steps.input.data"
        # read_orig_phase = xxx

        read_orig_phase_fullpath = False
        potentialpath = os.environ['potentials']
        print('potentialpath',potentialpath)
        print('potpath_using',potpath_using)
        if True:
            read_orig_phase_fullpath = potentialpath+"/n2p2_v4ag_ppl_987654_21cores/epoch_1383/export/aiida_exported_group_RELAXED_fully_aiida_exported_group_NN_relaxed_Mg5Al2Si4_n2p2_v2ag_calc__all_steps.input.data_calc__all_steps.input.data"
            print('read_orig_phase_fullpath',read_orig_phase_fullpath)
            print('os.path.ispath(read_orig_phase_fullpath)',os.path.isfile(read_orig_phase_fullpath))
            print()
            #read_orig_phase_fullpath = potpath_using+"/export/aiida_exported_group_out_reference_rep"+phase+".runner_calc__all_steps.input.data"
            #print('read_orig_phase_fullpath',read_orig_phase_fullpath)
            #print('os.path.ispath(read_orig_phase_fullpath)',os.path.isfile(read_orig_phase_fullpath))

        origphase_relaxed,origphase_initial = my.ase_relax_structure_fully_and_save_to_pot(ace,read=read_orig_phase,read_fullpath=read_orig_phase_fullpath,whichstruct=-1)
        print('type(origphase_relaxed)',type(origphase_relaxed))
        print('type(origphase_initial)',type(origphase_initial))
        # origphase_relaxe is the DFT relaxed one!
        # origphase_initial is just the first structure, possibly relaxed by a different NN.

        # print('type(origphase_relaxed)',type(origphase_relaxed))
        print('len(origphase_relaxed)',len(origphase_relaxed))
        # print('type(origphase_initial)',type(origphase_initial))
        print('len(origphase_initial)',len(origphase_initial))
        # print('origphase_relaxed',origphase_relaxed)
        # only the origphase initial will (in this case) have the DFT
        origphase_relaxed,origphase_initial = origphase_relaxed[0],origphase_initial[0]
        origphase_initial_ene_dft,kb_ = my.ase_enepot(origphase_initial,units='eV')
        max_dft_forces = np.abs(origphase_initial.get_forces()).max()
        max_dft_forces = np.round(max_dft_forces,5)
        print('(77) origphase_initial_ene_dft:',origphase_initial_ene_dft,'max_dft_forces',max_dft_forces)

        origphase_initial_ene_nn      = ace.ene(origphase_initial) #.copy())
        max_nn_forces = np.abs(origphase_initial.get_forces()).max()
        print('(88) origphase_initial_ene_nn :',origphase_initial_ene_nn ,'max_nn_forces ',max_nn_forces)

        origphase_relaxed_ene_nn      = ace.ene(origphase_relaxed) #.copy())
        max_nn_forces = np.abs(origphase_relaxed.get_forces()).max()
        print('(99) origphase_relaxed_ene_nn :',origphase_relaxed_ene_nn ,'max_nn_forces ',max_nn_forces)

        if False:
            # DEL!!origphase_relaxed_nn,origphase_relaxed_dft = origphase_relaxed[0],origphase_initial[0]
            # this was really minimized by dft
            orig_relaxed_dft_ene_dft,kb_ = my.ase_enepot(origphase_relaxed_dft,units='eV')
            max_dft_forces_o = np.abs(origphase_relaxed_dft.get_forces()).max()
            print('orig_relaxed_dft_ene_dft:',orig_relaxed_dft_ene_dft,"max_dft_forces_o",max_dft_forces_o)
            #print(origphase_relaxed_dft.get_forces())

            orig_relaxed_dft_ene_nn = ace.ene(origphase_relaxed_dft)
            max_nn_forces_o = np.abs(origphase_relaxed_dft.get_forces()).max()
            print('orig_relaxed_dft_ene_nn :',orig_relaxed_dft_ene_nn,"max_nn_forces_o",max_nn_forces_o)
            #print('fnn')
            #print(origphase_relaxed_dft.get_forces())

            orig_relaxed_nn_ene_nn  = ace.ene(origphase_relaxed_nn)
            max_nn_forces0_o = np.abs(origphase_relaxed_nn.get_forces()).max()
            print('orig_relaxed_nn_ene_nn  :',orig_relaxed_nn_ene_nn,"max_nn_forces0_o",max_nn_forces0_o)
            #print('fnn')
            #print(origphase_relaxed_nn.get_forces())

        d = my.ase_get_chemical_symbols_to_number_of_species(origphase_initial,known_elements_by_pot=["Al","Mg","Si"])
        print('11origphase --->',"Al:",d["Al"],"Mg:",d["Mg"],"Si:",d["Si"],'totl number of atoms:',d["Al"]+d["Mg"]+d["Si"])
        #print()
        #print()
        out = []
        for idx,i in enumerate(antisite_relaxed_by_nn):
            f = i.get_number_of_atoms()/origphase_initial.get_number_of_atoms()
            print()
            print()
            antisite_initiali = antisite_initial[idx]
            antisite_relaxedi_by_nn = antisite_relaxed_by_nn[idx]

            antisite_initiali_ene_dft,kb_ = my.ase_enepot(antisite_initiali,units='eV')
            max_dft_forces = np.abs(antisite_initiali.get_forces()).max()
            max_dft_forces = np.round(max_dft_forces,5)
            max_dft_forces2 = np.round(max_dft_forces,2)
            print('(22)antisite_initiali_ene_dft:',antisite_initiali_ene_dft,'max_dft_forces',max_dft_forces)

            antisite_initiali_ene_nn  = ace.ene(antisite_initiali) #.copy())
            max_nn_forces = np.abs(antisite_initiali.get_forces()).max()
            #max_nn_forces = np.round(max_nn_forces,5)
            print('(33)antisite_initiali_ene_nn :',antisite_initiali_ene_nn ,'max_nn_forces ',max_nn_forces)

            #check_forces = 0.00011
            #if max_nn_forces > check_forces:
            #    sys.exit('max forces > '+str(check_forces))

            antisite_relaxedi_ene_nn  = ace.ene(antisite_relaxedi_by_nn.copy())
            max_nn_forces0 = np.abs(antisite_relaxedi_by_nn.get_forces()).max()
            #max_nn_forces0 = np.round(max_nn_forces0,5)
            print('(44)antisite_relaxedi_ene_nn :',antisite_relaxedi_ene_nn ,'max_nn_forces0',max_nn_forces0)

            antisite_nat = i.get_number_of_atoms()
            #print('ace.E_SS_Al',ace.E_SS_Al)
            #print('ace.E_SS_Mg',ace.E_SS_Mg)
            #print('ace.E_SS_Si',ace.E_SS_Si)
            #print('f',f)


            dd = my.ase_get_chemical_symbols_to_number_of_species(i,known_elements_by_pot=["Al","Mg","Si"])
            al_antisite = dd["Al"]
            mg_antisite = dd["Mg"]
            si_antisite = dd["Si"]
            al_orig = d["Al"]*f
            mg_orig = d["Mg"]*f
            si_orig = d["Si"]*f
            print('f',f)
            print('xxx al_orig',al_orig,al_antisite)
            print('xxx mg_orig',mg_orig,mg_antisite)
            print('xxx si_orig',si_orig,si_antisite)
            print('xxx sum_orig',al_orig+mg_orig+si_orig)
            print('xxx sum_anti',al_antisite+mg_antisite+si_antisite)
            sys.eixt()

            al_d = al_antisite - al_orig
            mg_d = mg_antisite - mg_orig
            si_d = si_antisite - si_orig
            print('antisite phase     :',"Al:",dd["Al"],"Mg:",dd["Mg"],"Si:",dd["Si"],'antisite_relaxed_ene_nn   :') #,antisite_relaxedi_ene_nn)
            print('antisite phase     :',"Al:",dd["Al"],"Mg:",dd["Mg"],"Si:",dd["Si"],'antisite_ene_nn_diff      :') #,antisite_relaxedi_ene_nn-orig_relaxed_nn_ene_nn*f)
            print('al_d',al_d,'al_orig',al_orig,-al_d*ace.E_SS_Al)
            print('mg_d',mg_d,'mg_orig',mg_orig,-mg_d*ace.E_SS_Mg)
            print('si_d',si_d,'si_orig',si_orig,-si_d*ace.E_SS_Si)
            if al_d == -1 : text1 = "Al->"
            if mg_d == -1 : text1 = "Mg->"
            if si_d == -1 : text1 = "Si->"
            if al_d == 1 : text2 = "Al"
            if mg_d == 1 : text2 = "Mg"
            if si_d == 1 : text2 = "Si"
            text2=text2+"_"+str(max_dft_forces2)


            ## The errors in the E_SS_{Al,Mg,Si} (between DFT and NN) do not matter at all here, since those are way too small to influence the result since only two Al/Mg/Si
            #Eform_relaxed_nn = antisite_relaxedi_ene_nn - orig_relaxed_nn_ene_nn*f - al_d*ace.E_SS_Al - mg_d*ace.E_SS_Mg - si_d*ace.E_SS_Si
            Eform_relaxed_nn = 0

            #Eform_initial_nn = antisite_initiali_ene_nn - orig_relaxed_nn_ene_nn*f - al_d*ace.E_SS_Al - mg_d*ace.E_SS_Mg - si_d*ace.E_SS_Si
            Eform_initial_nn = antisite_initiali_ene_nn - origphase_initial_ene_nn*f - al_d*ace.E_SS_Al - mg_d*ace.E_SS_Mg - si_d*ace.E_SS_Si

            #Eform_relaxed_dft = antisite_relaxedi_ene_dft - orig_relaxed_dft_ene_dft*f - al_d*ace.E_SS_Al - mg_d*ace.E_SS_Mg - si_d*ace.E_SS_Si
            Eform_relaxed_dft = 0 # does not exist yet
            #Eform_initial_dft = antisite_initiali_ene_dft - orig_relaxed_dft_ene_dft*f - al_d*ace.E_SS_Al - mg_d*ace.E_SS_Mg - si_d*ace.E_SS_Si
            Eform_initial_dft = antisite_initiali_ene_dft - origphase_initial_ene_dft*f - al_d*ace.E_SS_Al - mg_d*ace.E_SS_Mg - si_d*ace.E_SS_Si
            #print("Eform_relaxed_nn :",Eform_relaxed_nn)
            #print("Eform_initial_nn :",Eform_initial_nn)
            #print("Eform_relaxed_dft:",Eform_relaxed_dft)
            #print("Eform_initial_dft:",Eform_initial_dft)
            print('al_d',al_d,mg_d,si_d)
            #print('from to:',text1,text2,
            #        Eform_relaxed_nn,Eform_initial_nn,
            #        Eform_relaxed_dft,Eform_initial_dft)
            #out.append([idx,Eform_relaxed_nn,Eform_initial_nn,Eform_relaxed_dft,Eform_initial_dft,text1+text2])
            out.append([idx,Eform_relaxed_nn,Eform_initial_nn,Eform_relaxed_dft,Eform_initial_dft,text1+text2])
            #print('out',out)
            #if idx == 3:
            #    sys.exit()
        myarray = np.asarray(out)
        print(myarray)
        np.savetxt(phase+"_formation.dat",myarray,fmt='%s')
        print('done-----------2')

    sys.exit('88888866666')

    doit = [ "Mg5Si6", "Mg5Al2Si4", "Mg4Al3Si4" ]
    #doit = [ "Mg5Si6" ]
    #doit = [ "Mg5Al2Si4" ]
    doit = [ "Mg4Al3Si4" ]
    pwd = os.getcwd()
    print('pwd',pwd)
    for i in doit:
        # where to save the formation energies
        #wsp = ace.savefolder+"summary_formations_antisites_"+i+"_"

        wsp = pwd+"/summary_formations_antisites_"+i+"_"
        ace.written_summary = [ wsp+"DFT_T0.dat", wsp+"NNatDFT_T0.dat",wsp+"NN_T0.dat" ]
        for wsp_ in ace.written_summary:
            if os.path.isfile(wsp_): os.remove(wsp_)
            print('saving to',wsp_)

        # get the formation energies
        print('antisites based on',i)
        path = os.environ['potentials']+"/aiida_get_structures_new/aiida_exported_group_out_antisites"+i+".runner_calc__all_steps.input.data"
        path = os.environ['potentials']+"/aiida_get_structures_new/aiida_exported_group_out_antisites_rep"+i+"_NEW.runner_calc__all_steps.input.data"
        frames = ase_read(path,":",format="runner")  # frames with DFT energies.
        #print('frames[0].nat',frames[0].get_number_of_atoms())
        #print('frames[0].nat',frames[0].get_number_of_atoms()/2)
        formula_unit = 11. #frames[0].get_number_of_atoms()/2
        for frame in frames:
            #hessematrix_try = ace.savefolder+"h_"+i+"_at_DFT_relaxed"
            #get_formation_energy(ace,frame,i+" (DFT@DFT fully relaxed)",atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=True,try_harmonic_readfile=hessematrix_try)
            hessematrix_try = ace.savefolder+"h_"+i+"_at_NN_relaxed"
            get_formation_energy(ace,frame,i+" (NN@NN  fully relaxed)", atomrelax=True, cellrelax=True, volumerelax=True, DFT_ene=True,try_harmonic_readfile=hessematrix_try,formula_unit=formula_unit)
            #ase_write(path+"NN_relaxed_"+i+"_"+ace.pot.pot+".runner",frame,format='runner')
    return

def test_inputfile_formation_energy(ace,args):
    # where to save the formation energies
    wsp = os.getcwd()+"/summary_formations_from_inpufile_eV_peratom_"
    ace.written_summary = [ wsp+"DFT_T0.dat", wsp+"NNatDFT_T0.dat",wsp+"NN_T0.dat" ]
    for wsp_ in ace.written_summary:
        if os.path.isfile(wsp_): os.remove(wsp_)
        print('saving to',wsp_)

    # get the formation energies
    frames = ase_read(args.inputfile,":",format="runner")
    for idx,frame in enumerate(frames):
        #hessematrix_try = ace.savefolder+"h_"+i+"_at_DFT_relaxed"
        #get_formation_energy(ace,frame,i+" (DFT@DFT fully relaxed)",atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=True,try_harmonic_readfile=hessematrix_try)
        #hessematrix_try = ace.savefolder+"h_"+i+"_at_NN_relaxed"
        get_formation_energy(ace,frame,str(idx)+" (NN@NN  fully relaxed)", atomrelax=True, cellrelax=True, volumerelax=True, DFT_ene=True,try_harmonic_readfile=False)
        #ase_write(path+"NN_relaxed_"+i+"_"+ace.pot.pot+".runner",frame,format='runner')
    return

def test_beta2_bulk(ace):
    print("######## test_beta2_bulk #############")
    #wsp = ace.savefolder+"summary_formations_beta_"
    wsp = os.getcwd()+"/summary_formations_eV_peratom_beta2_"
    ace.written_summary = [ wsp+"DFT_T0.dat", wsp+"NNatDFT_T0.dat",wsp+"NN_T0.dat" ]
    print('summary will be written to:')
    for i in ace.written_summary:
        if os.path.isfile(i): os.remove(i)
        print(i)
    print()
    # "Mg9Si5" == beta' (beta prime)
    # "Mg5Si6", "Mg5Al2Si4", "Mg4Al3Si4" are the three beta'' (beta double prime)
    doit = [ "Mg9Si5", "Mg5Si6", "Mg5Al2Si4", "Mg4Al3Si4" ]
    #doit = [ "Mg5Si6", "Mg5Al2Si4", "Mg4Al3Si4" ]
    #doit = [ "Mg5Si6" ]
    #doit = [ "Mg9Si5" ]

    for i in doit:
        print("########",i,"############")
        if i in [ "Mg5Si6", "Mg5Al2Si4", "Mg4Al3Si4" ]:
            formula_unit = 11
            path = os.environ['potentials']+"/aiida_get_structures_new/aiida_exported_group_NN_relaxed_"+i+"_n2p2_v2ag_calc__all_steps.input.data"
        elif i in [ "Mg9Si5" ]:
            formula_unit = 14
            path = os.environ['potentials']+"/aiida_get_structures_new/aiida_exported_group_BetaPrime_vc-relaxed__all_steps.input.data"

        print('import from path:',path)
        frame = ase_read(path,":",format="runner")
        #print('type(frame)',type(frame))
        frame = frame[-1] # DFT_relax
        if True:
            #hessematrix_try = ace.savefolder+"h_"+i+"_at_DFT_relaxed"
            #get_formation_energy(ace,frame,i+" (DFT@DFT fully relaxed)",atomrelax=False,cellrelax=False,volumerelax=False,DFT_ene=True,try_harmonic_readfile=hessematrix_try)
            hessematrix_try = ace.savefolder+"h_"+i+"_at_NN_relaxed"
            get_formation_energy(ace,frame,i+" (@DFT  fully relaxed)", atomrelax=True, cellrelax=True, volumerelax=True, DFT_ene=True,try_harmonic_readfile=hessematrix_try,formula_unit=formula_unit)
            ase_write(path+"NN_relaxed_"+i+"_"+ace.pot.pot+".runner",frame,format='runner')
    return


def load_diluete_pure_values():
    scripts = my.scripts()
    filename = scripts+'/tests/Al-Mg-Si/get_dilute_si_mg_f.'+pot+".dat"



def formation_energies(ace,args):
    #print('args.inputfile',args.inputfile)
    #my.get_Mg5Si6_and_other_antisites(ace)
    #sys.exit()
    print("########### formation_energies   #########################")
    ace.atTemp = 443
    print("########### get_basic_NN_energies_ace #########################")
    get_basic_NN_energies_ace(ace)


    print("########### get_dilute_si_mg_f #########################")
    get_dilute_si_mg_f(ace)

    file = ace.savefolder+"summary_formationsT0.dat"
    if os.path.isfile(file): os.remove(file)
    file = ace.savefolder+"summary_formationsT"+str(ace.atTemp)+".dat"
    if os.path.isfile(file): os.remove(file)


    if args.formation_energies == "inputfile":
        test_inputfile_formation_energy(ace,args)
    if args.formation_energies == "beta2":
        test_beta2_bulk(ace)
    if args.formation_energies == "beta2antisites":
        test_antisites(ace)

    #test_si_si_vac(ace)
    #test_Mg2Si(ace)
    #test_Mg9Si5(ace)

    ##test_Mg9Si5_pos(ace)
    ##test_betaprime_mg9si5_find_global_min(ace,eform_dilute_si, eform_dilute_mg, f_dilute_si_300, f_dilute_mg_300)
    return




def get_elastic_constants_al_ext(ace,frames):
    ace.elastic = True
    #get_basic_NN_energies_ace(ace)
    #get_al_fcc_equilibrium(ace)
    if ace.pot.c44_al:
        print('ace.c44:',ace.pot.c44_al,type(ace.pot.c44_al))
    else:
        ace.elastic_relax = True
        frame_al = frames
        ace.get_elastic_external(atomsin=frame_al,verbose=ace.verbose,text="Al_fcc bulk 4at",get_all_constants=True)
        print('ace.c44:',ace.c44,type(ace.c44))
        print('ace.elastic_constants:',ace.elastic_constants)
        file1 = open("elastic_constants_Al_fcc.txt","w")
        file1.write(ace.elastic_constants.decode("utf-8") )
        file1.close()
        #filename = ace.pot.potpath+"/elastic_"+str(ace.pot.potepoch_bestteste)+".dat"
        #if not os.path.isfile(filename):
        #    np.savetxt(filename,np.array([np.float(ace.c44)]))
    return ace.c44

def get_elastic_constants_al_from_ene(ace):
    ace.elastic_relax = True
    frame_al = my.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=1,nsi=0,nmg=0,nvac=0,a0=4.045,cell='cubic',create_fake_vacancy=False,crystal_structure="fcc")
    ace.get_elastic(frame_al,verbose=ace.verbose)
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
    start = timeit.timeit()
    p = help()
    print('__main__ after p = help',timeit.timeit()-start); start = timeit.timeit()
    args = p.parse_args()
    print('__main__ after parsing help',timeit.timeit()-start)
    if args.verbose:
        my.print_args(args)
    get_energies(args)
