#!/usr/bin/env python
from __future__ import print_function
import os,sys,random,argparse
import socket
import numpy as np
from shutil import copyfile

# from scripts folder
import convert_fileformats
import myutils as mu

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--ncell",'-ncell'   ,required=True, type=int,help="supercell size of primitive cell")
    p.add_argument("--nsi",'-nsi'       ,required=True, type=int, help="number of Si atoms")
    p.add_argument("--nmg",'-nmg'       ,required=True, type=int, help="number of Mg atoms")
    p.add_argument("--nvac",'-nvac'     ,required=True, type=int, help="number of vacancies")
    p.add_argument("--a0",'-a0'         ,default = 4.057, type=float, help="fcc lattice constant for Al.")
    p.add_argument("--temp",'-temp'     ,default = 300, type=int, help="KMC temperature")
    p.add_argument("--nseeds",'-nseeds' ,type=int, default=3, help="number of different seeds")
    p.add_argument("--seeds",'-seeds',default=False,nargs='*', type=int,help="define seed(s) manually (otherwise a random number generator; can be used to define multiple seeds)")
    p.add_argument('-nsteps',type=int, default=50000, help="number of KMC steps to make")
    p.add_argument('-fa','--foldername_append',type=str, default="", help="append to foldername")
    p.add_argument('--pot','-p',       required=False, choices=mu.pot_all(), default=mu.get_latest_n2p2_pot())
    p.add_argument('-cubic', '--cubic', action='store_true',default=False,help="Use cubic(==cconventional) fcc cell instead of primitive one")
    p.add_argument('-submit', '--submit', action='store_true')
    p.add_argument('-submitdebug','--submitdebug',action='store_true')
    p.add_argument('-st','--submittime_hours', type=int,default=71,help="slurm time for the job")
    p.add_argument('-n','--nodes', type=int,default=1,help="how many nodes to use?")
    p.add_argument('-t','--test', action='store_true')
    p.add_argument('-tf','--testfiles', action='store_true')
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p


def createjob(args):
    """
    This is an script to create KMC jobs quickly (and submits them if -submit/-submitdebug).

    e.g.

    createFolder_kmc.py -temp 1000 -ncell 5 -nsi 3 -nmg 3 -nvac 1 -submit

    to evaluate energies vs. realtime use:\n
            %kmc_show_time_xmgrace.sh seed*/KMC_AL6XXX
    """
    ncell = args.ncell
    nmg                 = args.nmg
    nsi                 = args.nsi
    nvac                = args.nvac
    a0                  = args.a0
    temp                = args.temp
    nseeds              = args.nseeds
    seeds          = args.seeds
    nsteps              = args.nsteps
    foldername_append   = args.foldername_append
    pot                 = args.pot
    submit              = args.submit
    submitdebug         = args.submitdebug
    submittime_hours    = args.submittime_hours
    test                = args.test
    testfiles           = args.testfiles
    nodes               = args.nodes
    verbose             = args.verbose


    ### check if ase runner/quippy/lammpps-data formats are known
    ase_formats = mu.ase_get_known_formats_class(verbose=True)
    ase_formats.check_if_default_formats_known(copy_and_adapt_formatspy_anyhow=False)

    # definex ffsocket inet/unix
    if nodes == 1:
        ffsocket = "unix"
    elif nodes > 1:
        ffsocket = "inet"
    else:
        sys.exit("Number of nodes has to be positive!")


    # define ntasks, neval
    lmp_par = 2          # = OMP_NUM_THREADS
    ntasks = cores = nodes * 28
    ipi_inst = 4         # for sure best on fidis
    neval  = ipi_inst*2  # was alwasy better, for ompi and impi

    ##### get the seed(s).
    if type(seeds) == bool:
        seeds = random.sample(range(1, 999999), nseeds)
    print('seeds',seeds)
    if test == True:
        nseeds = 1
        seeds = [1]
        print('seeds',seeds)
    nseeds = len(seeds)

    ##### a few checks
    scripts = mu.scripts()
    mypot = mu.mypot(pot)
    if submit is True or submitdebug is True:
        mu.check_for_known_hosts()


    ##### here only chck if the potential can be set up. (in.lmp)
    ##### the same command is then executed for every kmc folder
    ace = mu.ase_calculate_ene(pot=pot,
            potpath=False,
            units='eV',geopt=False,kmc=True,verbose=verbose)
    ace.pot_get_and_ase_lmp_cmd(kmc=True,temp=temp,nsteps=nsteps,ffsocket=ffsocket)

    ##### if test
    if test == True:
        nsteps = 50

    file_ipi_input_runner   = scripts + "/i-pi-mc_scripts/input-runner.xml"


    ####################################
    # get directory
    ####################################
    if verbose:
        print("get directory")
    pcsi = nsi/ncell**3.*100
    pcmg = nmg/ncell**3.*100
    pcvac = nvac/ncell**3.*100
    directory = str(ncell)+"x"+str(ncell)+"x"+str(ncell)+"_"+pot+"_"+\
                str(temp)+"K_"+\
                str(nvac)+"Vac_"+str(nmg)+"Mg_"+str(nsi)+"Si__"+\
                str(round(pcvac,3))+"pctVac_"+str(round(pcmg,3))+"pctMg_"+str(round(pcsi,3))+"pctSi"
    if foldername_append != "":
        directory = directory+"_"+foldername_append

    ###############################################
    # make the structure
    ###############################################
    atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0,create_fake_vacancy = True,cubic=args.cubic)
    atomsc = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0,cubic=args.cubic)

    # make the atomic structure
    # this was to play ... not necessary now?
    if False:
        nndist = a0/np.sqrt(2.)

        from ase.io import read as ase_read
        from ase.io import write as ase_write

        ###############################################
        # get the amount of 1NN in a relly large cell
        ###############################################
        atomsc_fakevac_i = ase_read('dataxx.extxyz3',index=":",format='extxyz') # works, cell ist not changed
        #atomsc_fakevac_i = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=10,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
        #nn = mu.ase_get_neighborlist(atomsc_fakevac_i,atomnr=0,cutoff=3.,skin=0.1)
        #print("nn",nn,'len',len(nn))
        #nn = mu.ase_get_neighborlist(atomsc_fakevac_i,atomnr=0,cutoff=8.5,skin=0.1)
        #print("nn",nn,'len',len(nn))
        #sys.exit()

        print(len(atomsc_fakevac_i),type(atomsc_fakevac_i))

        for idx,i in enumerate(atomsc_fakevac_i):
            print('aa',atomsc_fakevac_i[idx].positions[0])
            #print('aa',i.positions[0])
        print('ipi')
        atomsc_fakevac_i = ase_read('dataxx.ipi2',index=":",format='ipi') # works, cell ist not changed
        print(len(atomsc_fakevac_i),type(atomsc_fakevac_i))
        for idx,i in enumerate(atomsc_fakevac_i):
            print('aa',atomsc_fakevac_i[idx].positions[0])
            #print('aa',i.positions[0])
        print('quippy')
        atomsc_fakevac_i = ase_read('dataxx.quippy.xyz2',index=":",format='quippy') # works, cell ist not changed



        filename = '../sim.xyz'
        filename = '../simulation.pos_0.xyz'
        mu.count_amount_1NN_around_vacancies(filename,cutoffa=nndist,cutoffb=a0,skin=0.1,format='ipi')
        sys.exit()

        def mysave_quippy_xyz(atomsc_fakevac,text=False):
            if type(text) == bool:
                sys.exit('define text')
            atomsc_fakevac.write('data.quippy.xyz',format='quippy',append=True)
            #atomsc_fakevac.write('data.xyz',format="extxyz",append=True)
            atomsc_fakevac.write('data'+text+'.quippy.xyz',format='quippy',append=True)
            #atomsc_fakevac.write('data'+text+'.xyz',format="extxyz",append=True)
            return

        # create Al with single vacancy
        atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
        NN_1_indices, NN_2_indices = mu.ase_get_neighborlist_1NN_2NN(atomsc_fakevac,atomnr=0,cutoffa=nndist,cutoffb=a0,skin=0.1)
        #print('from ....',(atomsc_fakevac.positions)[0])
        #for i in NN_1_indices:
        #    print((atomsc_fakevac.positions)[i])
        print('NN_1_indices (orig  ):',NN_1_indices)
        print('NN_2_indices (orig  ):',NN_2_indices)
        #sys.exit()
        atomsc_fakevac.write('dataxx.quippy.xyz',format='quippy',append=True)
        atomsc_fakevac.write('dataxx.poscar',format='vasp',append=True)
        atomsc_fakevac.write('dataxx.ipi',format='ipi',append=True) # works, currently so implemented that it canges cell
        atomsc_fakevac.write('dataxx.xyz',format='xyz',append=True)
        atomsc_fakevac.write('dataxx.extxyz',format='extxyz',append=True)
        atomsc_fakevac.write('dataxx.lammps-data',format='lammps-data',append=True)
        atomsc_fakevac.write('dataxx.lammps-runner',format='lammps-runner',append=True)

        atomsc_fakevac_a = ase_read('dataxx.extxyz',format='extxyz') # works, cell ist not changed
        atomsc_fakevac_a.write('dataxx.extxyz2',format='extxyz',append=True) # works, cell is not changed

        atomsc_fakevac_b = ase_read('dataxx.xyz',format='xyz') # not working # but this should work
        atomsc_fakevac_b.write('dataxx.xyz2',format='xyz',append=True) # this is working

        atomsc_fakevac_c = ase_read('dataxx.ipi',format='ipi')  # works, currently so implemented that it canges cell
        #print('ipi cell',atomsc_fakevac_c.get_cell())

        atomsc_fakevac_c.write('dataxx.ipi2',format='ipi',append=True) # works, just writes the cell it gests.
        atomsc_fakevac_c.write('dataxx.ipi2_poscar',format='vasp',append=True) # works, just writes the cell it gests.
        NN_1_indices, NN_2_indices = mu.ase_get_neighborlist_1NN_2NN(atomsc_fakevac_c,atomnr=0,cutoffa=nndist,cutoffb=a0,skin=0.1)
        print('NN_1_indices (ipi   ):',NN_1_indices)
        print('NN_2_indices (ipi   ):',NN_2_indices)
        #print('from ....',(atomsc_fakevac_c.positions)[0])
        #for i in NN_1_indices:
        #    print((atomsc_fakevac_c.positions)[i])

        atomsc_fakevac_cc = ase_read('dataxx.ipi2_poscar',format='vasp')  # works, currently so implemented that it canges cell
        atomsc_fakevac_cc.write('dataxx.ipi2_poscar2',format='vasp',append=True)
        atomsc_fakevac_cc.write('dataxx.ipi2_poscar2_ipi',format='ipi',append=True) # works, just writes the cell it gests.
        #print('ipi cell2 (ext):',atomsc_fakevac_cc.get_cell())
        #print()
        #print('now quippy')
        atomsc_fakevac_d = ase_read('dataxx.quippy.xyz',format='quippy')
        #print('quippy cell (ext)',atomsc_fakevac_d.get_cell())
        atomsc_fakevac_d.write('dataxx.quippy.xyz2',format='quippy',append=True)
        atomsc_fakevac_d.write('dataxx.quippy.xyz2_extxyz',format='extxyz',append=True)
        NN_1_indices, NN_2_indices = mu.ase_get_neighborlist_1NN_2NN(atomsc_fakevac_d,atomnr=0,cutoffa=nndist,cutoffb=a0,skin=0.1)
        print('NN_1_indices (quippy):',NN_1_indices)
        print('NN_2_indices (quippy):',NN_2_indices)
        #print('from ....',(atomsc_fakevac_d.positions)[0])
        #for i in NN_1_indices:
        #    print((atomsc_fakevac_d.positions)[i])
        path = "/home/glensk/kmc/run_michele/Si6Mg6V1.1_/simulation.pos_libatom_2struct.xyz"
        atomsc_fakevac_e = ase_read(path,format='quippy')

        NN_1_indices, NN_2_indices = mu.ase_get_neighborlist_1NN_2NN(atomsc_fakevac_e,atomnr=0,cutoffa=nndist,cutoffb=a0,skin=0.1)
        print('NN_1_indices (kmc   ):',NN_1_indices)
        print('NN_2_indices (kmc   ):',NN_2_indices)
        sys.exit()

        NN_1_indices = mu.ase_get_neighborlist(atomsc_fakevac,atomnr=0,cutoff=nndist,skin=0.1)
        NN_1_2_indices_tmp = mu.ase_get_neighborlist(atomsc_fakevac,atomnr=0,cutoff=a0,skin=0.1)
        print('NN_1_indices  :',NN_1_indices)
        NN_2_indices = np.sort(np.array(mu.diff(NN_1_2_indices_tmp,NN_1_indices)))
        print('NN_2_indices  :',NN_2_indices)
        NN_1_2_indices = np.concatenate((NN_1_indices, NN_2_indices ))
        print('NN_1_2_indices:',NN_1_2_indices)


        # fill only 1NN (with one species)
        for i in [ 'Mg', 'Si' ]:
            atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
            mysave_quippy_xyz(atomsc_fakevac,text="1NN")
            for ii in NN_1_indices:
                atomsc_fakevac[ii].symbol = i
                mysave_quippy_xyz(atomsc_fakevac,text="1NN")

        # fill only 2NN (with one species)
        for i in [ 'Mg', 'Si' ]:
            atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
            mysave_quippy_xyz(atomsc_fakevac,text="2NN")
            for ii in NN_2_indices:
                atomsc_fakevac[ii].symbol = i
                mysave_quippy_xyz(atomsc_fakevac,text="2NN")

        # fill 1NN and 2NN (with one species)
        for i in [ 'Mg', 'Si' ]:
            atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
            mysave_quippy_xyz(atomsc_fakevac,text="1and2NN")
            for ii in NN_1_2_indices:
                atomsc_fakevac[ii].symbol = i
                mysave_quippy_xyz(atomsc_fakevac,text="1and2NN")

        # dif compositions in 1NN shell
        filling = [ 2,4,6,8,10]
        for fi in filling:
            atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
            mysave_quippy_xyz(atomsc_fakevac,text="1NN_diffcomp")
            for idx,ii in enumerate(NN_1_indices):
                if idx < fi: ch = "Mg"
                else: ch = "Si"
                atomsc_fakevac[ii].symbol = ch
            mysave_quippy_xyz(atomsc_fakevac,text="1NN_diffcomp")


        sys.exit()

        #mu.ase_get_known_formats(show=True, add_missing_formats=False, copy_formats=False, verbose=False,show_formatspy=True)
        for i in [ 'Mg', 'Si' ]:
            for ii in [ 0,1,2,3,4,5]:
                atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering=i+'_'+str(ii))


        sys.exit()


    # show the input variables
    print('--------------------------- check the input --------------------------------')
    print('JOBS (nseeds) ',nseeds,'(defined by -nseeds / or -seeds)')
    print('seeds         ',seeds)
    print('nsteps        ',nsteps)
    print()
    print('ncell         ',ncell,"(",atomsc.get_number_of_atoms(),"atoms )")
    print('nsi           ',nsi,  "(",pcsi,"at%)")
    print('nmg           ',nmg,"(",pcmg,"at%)")
    print('nvac          ',nvac,"(",pcvac,"at%)")
    print('a0            ',a0,"angstrom")
    print('temp          ',temp,"K")
    print()
    print('mypot.pot     ',mypot.pot)
    print('mypot.potpath ',mypot.potpath)
    print()
    print('directory     ',directory)
    print('submit        ',submit)
    print('submitdebug   ',submitdebug)
    print()
    print('nodes         ',nodes)
    print('ffsocket      ',ffsocket)
    #print('python ver    ',sys.version_info[0])
    #print()
    print('--------------------------- check the input --------------------------------')
    if submit == True or submitdebug == True:
        mu.get_from_prompt_Yy_orexit("Are the ine input variables ok? [y]es: ")

    # make the directory
    if os.path.isdir(directory):
        mu.get_from_prompt_Yy_orexit("This main directory exists already, shall I add jobs? [y]es: ")
    mu.mkdir(directory)

    # create README.md
    IPI_COMMAND = os.environ["IPI_COMMAND"]
    LAMMPS_COMMAND = os.environ["LAMMPS_COMMAND"]
    mu.create_READMEtxt(directory,add=["# to start manually (1): python "+IPI_COMMAND+" input-runner.xml","# to start manually (2):"+LAMMPS_COMMAND+" < in.lmp"])

    for seed in seeds:

        # make jobdirectory
        jobdir = directory+'/seed'+str(seed)
        print('jobdir',jobdir)
        if os.path.exists(jobdir):
            sys.exit("jobdirectory "+str(jobdir)+" already exists!")
        mu.mkdir(jobdir)

        # get data.lmp and data.ipi
        atomsc.write(jobdir+'/data.runnerformat.lmp',format='lammps-runner')
        atomsc_fakevac.write(jobdir+'/data.ipi',format='ipi')
        atomsc_fakevac.write(jobdir+'/data.extxyz',format='extxyz')
        #atomsc_fakevac.write(jobdir+'/data_fakevac.ipi',format='ipi')

        if testfiles == True:
            atomsc.write(jobdir+'/data.lmp',format='lammps-data')
            atomsc.write(jobdir+'/data.POSCAR',format='vasp')
            atomsc.write(jobdir+'/data.xyz',format='xyz')
            atomsc.write(jobdir+'/data.extxyz',format='extxyz')
            atomsc.write(jobdir+'/data.espresso-in',format='espresso-in')

        # create in.lmp
        ace = mu.ase_calculate_ene(pot=pot,potpath=mypot.potpath,
                units='eV',geopt=False,kmc=True,verbose=verbose)
        address = socket.gethostname()+"_"+os.path.basename(jobdir)
        print('address',address)
        ace.pot_get_and_ase_lmp_cmd(kmc=True,temp=temp,nsteps=nsteps,ffsocket=ffsocket,address=address)
        mu.lammps_write_inputfile(folder=jobdir,filename='in.lmp',positions='data.runnerformat.lmp',ace=ace)

        # create input-runner.xml (should be made without copying)
        mu.create_ipi_kmc_inputfile(jobdir,filename="input-runner.xml",nsteps=nsteps,stride=100,seed=seed,a0=a0,ncell=ncell,nsi=nsi,nmg=nmg,nvac=nvac,neval=neval,temp=temp,nodes=nodes,address=address,testrun=test,cubic=args.cubic)

        # create submit-ipi-kmc.sh (should be made without copying)
        mu.create_submitskript_ipi_kmc(jobdir+"/submit-ipi-kmc.sh",nodes,ntasks,
                lmp_par=lmp_par,
                ipi_inst=ipi_inst,
                ffsocket=ffsocket,
                submittime_hours=submittime_hours,
                SBATCH=True)

        # create osubmit-ipi-kmc.sh (should be made without copying)
        mu.create_submitskript_ipi_kmc(jobdir+"/osubmit-ipi-kmc.sh",nodes,ntasks,
                lmp_par=lmp_par,
                ipi_inst=ipi_inst,
                ffsocket=ffsocket,
                submittime_hours=submittime_hours,
                SBATCH=False)

        # submit the job (execute either this or submit-ipi-kmc.sh_all3, not both)
        #mu.submitjob(submit=submit,submitdebug=submitdebug,jobdir=jobdir,submitskript="submit-ipi-kmc.sh")

    # get submit-ipi-kmc.sh_all3 (should be made without copying)
    if nseeds == 3:
        mu.create_submitskript_ipi_kmc(directory+"/submit-ipi-kmc.sh_all3",nodes,ntasks,
                lmp_par=lmp_par,
                ipi_inst=ipi_inst,
                ffsocket=ffsocket,
                submittime_hours=submittime_hours,
                SBATCH=True,
                LOOPFOLDER=True)

        # submit the job (execute either this or submit-ipi-kmc.sh_all3, not both)
        #mu.submitjob(submit=submit,submitdebug=submitdebug,jobdir=directory,submitskript="submit-ipi-kmc.sh_all3")
        if submit == True:
            mu.submitjob(submit_to_que=True,submit_to_debug_que=False,jobdir=directory,submitskript="submit-ipi-kmc.sh_all3")


    print('done')
    return


if __name__ == "__main__":
    p = help()
    args = p.parse_args()
    if args.verbose:
        mu.print_args(args)
    createjob(args)
