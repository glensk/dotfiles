#!/usr/bin/env python
from __future__ import print_function
import os,sys,random
import socket
import numpy as np
from shutil import copyfile
import click

# from scripts folder
import convert_fileformats
import myutils as mu

CONTEXT_SETTINGS = mu.get_click_defaults()
@click.command(context_settings=CONTEXT_SETTINGS)


# input variables (get those with click)
# setting KMC
@click.option('-ncell',required=True, prompt=True, type=int,
        help="supercell size of primitive cell")
@click.option('-nsi'  ,required=True, prompt=True, type=int, help="number of Si atoms")
@click.option('-nmg'  ,required=True, prompt=True, type=int, help="number of Mg atoms")
@click.option('-nvac' ,required=True, prompt=True, type=int, help="number of vacancies")
@click.option('-a0'   ,default = 4.057, type=float, help="fcc lattice constant for Al.")
@click.option('-temp' ,default = 300, type=int, help="KMC temperature")
@click.option('-nseeds',type=int, default=3, help="number of different seeds")
@click.option('-seednumber',default=False,multiple=True, type=int,help="define seed number manually (otherwise a random number generator; can be used to define multiple seeds)")
@click.option('-nsteps',type=int, default=50000, help="number of KMC steps to make")
@click.option('-runnercutoff',type=float, default=10., help="runner cutoff distance ~10Angstrom")
@click.option('-fa','--foldername_append',type=str, default="", help="append to foldername")

# environment variables
@click.option('--pot','-p',type=click.Choice(mu.pot_all()),required=True,default=mu.get_latest_n2p2_pot(),help="potential from $scripts/potentials folder.")

@click.option('-submit/-no-submit', default=False)
@click.option('-submitdebug/-no-submitdebug', default=False)
@click.option('-st','--submittime_hours', type=int,default=71,help="slurm time for the job")
@click.option('-n','--nodes', type=int,default=1,help="how many nodes to use?")
@click.option('-t','--test/--no-test', default=False)
@click.option('-tf','--testfiles/--no-testfiles', default=False)
@click.option('--verbose','-v',count=True)

def createjob(
        ncell,
        nmg,
        nsi,
        nvac,
        a0,
        temp,
        nseeds,
        seednumber,
        nsteps,
        runnercutoff,
        foldername_append,
        pot,
        submit,
        submitdebug,
        submittime_hours,
        test,
        testfiles,
        nodes,
        verbose):
    """
    This is an script to create KMC jobs quickly (and submits them if -submit/-submitdebug).

    e.g.

    createFolder_kmc.py -temp 1000 -ncell 5 -nsi 3 -nmg 3 -nvac 1 -submit

    to evaluate energies vs. realtime use:\n
            %kmc_show_time_xmgrace.sh seed*/KMC_AL6XXX
    """
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
    seeds = random.sample(range(1, 999999), nseeds)
    if test == True:
        nseeds = 1
        seeds = [1]
        print()
    seednumber = list(seednumber)
    if len(seednumber) is not 0:
        seeds = seednumber
        nseeds = len(seednumber)

    ##### a few checks
    scripts = mu.scripts()
    mypot = mu.mypot(pot)
    if submit is True or submitdebug is True:
        if socket.gethostname() != "fidis":    # use != and not: is not
            print('you are NOT on fidis')
            sys.exit('submit or submitdebug is True but you are no fidis! Exit.')



    ##### here only chck if the potential can be set up. (in.lmp)
    ##### the same command is then executed for every kmc folder
    ace = mu.ase_calculate_ene(pot=pot,
            potpath=False,
            units='eV',geopt=False,kmc=True,verbose=verbose)
    mu.ase_calculate_ene.pot_to_ase_lmp_cmd(ace,kmc=True,temp=temp,nsteps=nsteps,ffsocket=ffsocket)

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
                str(round(pcvac,3))+"pctVac_"+str(round(pcmg,3))+"pctMg_"+str(round(pcsi,3))+"pctSi_"+\
                str(runnercutoff)+"rcut"
    if foldername_append != "":
        directory = directory+"_"+foldername_append


    # make the atomic structure
    if True:
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
        mu.count_amount_1NN_around_vacancies(filename,format='ipi')
        sys.exit()


        filename = "al_mg_si.dat"
        if os.path.isfile(filename):
            al_mg_si = np.loadtxt(filename)
        else:
            al_mg_si = np.zeros((len(atomsc_fakevac_i),3))

        for idx in np.arange(go_over):
            vac_idx = ([atom.index for atom in atomsc_fakevac_i[idx] if atom.symbol == 'V'])
            #print('vac_idx',vac_idx)
            for vac in vac_idx:
                #print('aa',atomsc_fakevac_i[idx].positions[vac])
                NN_1_indices, NN_2_indices = mu.ase_get_neighborlist_1NN_2NN(atomsc_fakevac_i[idx],atomnr=vac,cutoffa=nndist,cutoffb=a0,skin=0.1)
                NN_1_sym = [atom.symbol for atom in atomsc_fakevac_i[idx] if atom.index in NN_1_indices]
                NN_1_al = NN_1_sym.count("Al")
                NN_1_mg = NN_1_sym.count("Mg")
                NN_1_si = NN_1_sym.count("Si")
                print(idx,'NN_1_al',NN_1_al,"NN_1_mg",NN_1_mg,'NN_1_si',NN_1_si)
                al_mg_si.append([idx,NN_1_al,NN_1_mg,NN_1_si])
                print(al_mg_si)
                #if idx in np.arange(go_over)[::2]:
                #    np.savetxt(

                #sym = ([atom.symbol for atom in atomsc_fakevac_i[idx] if atom.== 'V'])

                #print('NN_1_indices',NN_1_indices)
                #print('NN_2_indices',NN_2_indices)
        sys.exit()

        def mysave(atomsc_fakevac,text=False):
            if type(text) == bool:
                sys.exit('define text')
            atomsc_fakevac.write('data.quippy.xyz',format='quippy',append=True)
            #atomsc_fakevac.write('data.xyz',format="extxyz",append=True)
            atomsc_fakevac.write('data'+text+'.quippy.xyz',format='quippy',append=True)
            #atomsc_fakevac.write('data'+text+'.xyz',format="extxyz",append=True)
            return

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
            mysave(atomsc_fakevac,text="1NN")
            for ii in NN_1_indices:
                atomsc_fakevac[ii].symbol = i
                mysave(atomsc_fakevac,text="1NN")

        # fill only 2NN (with one species)
        for i in [ 'Mg', 'Si' ]:
            atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
            mysave(atomsc_fakevac,text="2NN")
            for ii in NN_2_indices:
                atomsc_fakevac[ii].symbol = i
                mysave(atomsc_fakevac,text="2NN")

        # fill 1NN and 2NN (with one species)
        for i in [ 'Mg', 'Si' ]:
            atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
            mysave(atomsc_fakevac,text="1and2NN")
            for ii in NN_1_2_indices:
                atomsc_fakevac[ii].symbol = i
                mysave(atomsc_fakevac,text="1and2NN")

        # dif compositions in 1NN shell
        filling = [ 2,4,6,8,10]
        for fi in filling:
            atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering="XX_0")
            mysave(atomsc_fakevac,text="1NN_diffcomp")
            for idx,ii in enumerate(NN_1_indices):
                if idx < fi: ch = "Mg"
                else: ch = "Si"
                atomsc_fakevac[ii].symbol = ch
            mysave(atomsc_fakevac,text="1NN_diffcomp")


        sys.exit()

        #mu.ase_get_known_formats(show=True, add_missing_formats=False, copy_formats=False, verbose=False,show_formatspy=True)
        for i in [ 'Mg', 'Si' ]:
            for ii in [ 0,1,2,3,4,5]:
                atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell=5,nsi=0,nmg=0,nvac=1,a0=a0,cubic=False,create_fake_vacancy = True,normal_ordering=i+'_'+str(ii))


        sys.exit()


    # show the input variables
    print('--------------------------- check the input --------------------------------')
    #print('seednumber   ',seednumber,type(seednumber))
    print('JOBS (nseeds) ',nseeds,'(defined by -nseeds / or -seednumber)')
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
        mu.ase_calculate_ene.pot_to_ase_lmp_cmd(ace,kmc=True,temp=temp,nsteps=nsteps,ffsocket=ffsocket,address=address)
        mu.lammps_write_inputfile(folder=jobdir,filename='in.lmp',positions='data.runnerformat.lmp',ace=ace)


        # get input-runner.xml (should be made without copying)
        stride = 100
        verbosity = "low"
        if test == True:
            nsteps = 5
            stride = 1
            verbosity = "low"   # with high every socket connect is reported ... which is too much info

        copyfile(file_ipi_input_runner, jobdir+"/input-runner.xml")
        mu.sed(jobdir+"/input-runner.xml",'<total_steps>.*</total_steps>','<total_steps> '+str(nsteps)+' </total_steps>')
        mu.sed(jobdir+"/input-runner.xml",'<simulation verbosity=.*','<simulation verbosity="'+verbosity+'">')
        mu.sed(jobdir+"/input-runner.xml",'<trajectory filename.*','<trajectory filename="pos" stride="'+str(stride)+'" cell_units="angstrom" format="xyz" bead="0"> positions{angstrom} </trajectory>')
        mu.sed(jobdir+"/input-runner.xml",'<seed>.*</seed>','<seed> '+str(seed)+' </seed>')
        mu.sed(jobdir+"/input-runner.xml",'<a0 units="angstrom">.*</a0>','<a0 units="angstrom"> '+str(a0)+' </a0>')
        mu.sed(jobdir+"/input-runner.xml",'<ncell>.*</ncell>','<ncell> '+str(ncell)+' </ncell>')
        mu.sed(jobdir+"/input-runner.xml",'<nsi>.*</nsi>','<nsi> '+str(nsi)+' </nsi>')
        mu.sed(jobdir+"/input-runner.xml",'<nmg>.*</nmg>','<nmg> '+str(nmg)+' </nmg>')
        mu.sed(jobdir+"/input-runner.xml",'<nvac>.*</nvac>','<nvac> '+str(nvac)+' </nvac>')
        mu.sed(jobdir+"/input-runner.xml",'<neval>.*</neval>','<neval> '+str(neval)+' </neval>')
        mu.sed(jobdir+"/input-runner.xml",'<temperature units="kelvin">.*','<temperature units="kelvin">'+str(temp)+'</temperature>')


        #mu.sed(jobdir+"/input-runner.xml",'<file mode="xyz" units="angstrom">.*</file>','<file mode="xyz" units="angstrom"> data.ipi </file>')
        activelist = str(range(ncell**3 - nvac))
        activelist = activelist.replace(" ", "")
        mu.sed(jobdir+"/input-runner.xml",'<!-- <activelist> .*','<activelist> '+activelist+' </activelist>')
        atom_x_list = []
        for nvac_idx in range(ncell**3 - nvac,ncell**3):
            atom_x_list.append('atom_x{angstrom}('+str(nvac_idx)+')')
        insert = ", ".join(atom_x_list)
        mu.sed(jobdir+"/input-runner.xml",'atom_x{.*',insert+" ] </properties>")


        mu.sed(jobdir+"/input-runner.xml",'<file mode="xyz" units="angstrom">.*</file>','<file mode="xyz" units="angstrom"> data.ipi </file>')

        mu.sed(jobdir+"/input-runner.xml",'<ffsocket.*','<ffsocket name="lmpserial" mode="'+str(ffsocket)+'">')
        addressline = '<address> '+address+' </address>'
        if ffsocket == "unix":
            mu.sed(jobdir+"/input-runner.xml",'<address.*',addressline+' <latency> 1e-3 </latency>')
        if ffsocket == "inet":
            mu.sed(jobdir+"/input-runner.xml",'<address.*',addressline+' <port> 12345 </port>')



        # get submit-ipi-kmc.sh (should be made without copying)
        mu.create_submitskript_ipi_kmc(jobdir+"/submit-ipi-kmc.sh",nodes,ntasks,
                lmp_par=lmp_par,
                ipi_inst=ipi_inst,
                ffsocket=ffsocket,
                submittime_hours=submittime_hours,
                SBATCH=True)

        # get osubmit-ipi-kmc.sh (should be made without copying)
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
        mu.submitjob(submit=submit,submitdebug=submitdebug,jobdir=directory,submitskript="submit-ipi-kmc.sh_all3")


    print('done')
    return


if __name__ == "__main__":
    createjob()
