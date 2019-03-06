#!/usr/bin/env python
from __future__ import print_function
import os,sys,random
import socket
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
@click.option('--pot','-p',type=click.Choice(mu.pot_all()),required=True,default=my.get_latest_n2p2_pot(),help="potential from $scripts/potentials folder.")

@click.option('-submit/-no-submit', default=False)
@click.option('-submitdebug/-no-submitdebug', default=False)
@click.option('-st','--submittime_hours', type=int,default=71,help="slurm time for the job")
@click.option('-n','--nodes', type=int,default=1,help="how many nodes to use?")
@click.option('-t','--test/--no-test', default=False)
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
        nodes,
        verbose):
    """
    This is an script to create KMC jobs quickly (and submits them if -submit/-submitdebug).

    e.g.

    createFolder_kmc.py -temp 1000 -ncell 5 -nsi 3 -nmg 3 -nvac 1 -submit

    to evaluate energies vs. realtime use:\n
            %kmc_show_time_xmgrace.sh seed*/KMC_AL6XXX
    """
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
    ace = mu.ase_calculate_ene(pot,units='eV',geopt=False,kmc=True,verbose=verbose)
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
    atomsc_fakevac = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0,create_fake_vacancy = True)
    atomsc = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0)


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
    print('mypot.fullpath',mypot.fullpath)
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
    mu.get_from_prompt_Yy_orexit("Are the ine input variables ok? [y]es: ")

    # make the directory
    if os.path.isdir(directory):
        mu.get_from_prompt_Yy_orexit("This main directory exists already, shall I add jobs? [y]es: ")
    mu.mkdir(directory)

    # create README.md
    mu.create_READMEtxt(directory)

    for seed in seeds:

        # make jobdirectory
        jobdir = directory+'/seed'+str(seed)
        print('jobdir',jobdir)
        if os.path.exists(jobdir):
            sys.exit("jobdirectory "+str(jobdir)+" already exists!")
        mu.mkdir(jobdir)

        # get data.lmp and data.ipi
        atomsc.write(jobdir+'/data.lmp.runner',format='lammps-runner')
        atomsc.write(jobdir+'/data.ipi',format='ipi')
        atomsc_fakevac.write(jobdir+'/data_fakevac.ipi',format='ipi')

        if test == True:
            atomsc.write(jobdir+'/data.lmp',format='lammps-data')
            atomsc.write(jobdir+'/data.POSCAR',format='vasp')
            atomsc.write(jobdir+'/data.xyz',format='xyz')
            atomsc.write(jobdir+'/data.extxyz',format='extxyz')
            atomsc.write(jobdir+'/data.espresso-in',format='espresso-in')


        # create in.lmp
        ace = mu.ase_calculate_ene(pot,units='eV',geopt=False,kmc=True,verbose=verbose)
        mu.ase_calculate_ene.pot_to_ase_lmp_cmd(ace,kmc=True,temp=temp,nsteps=nsteps,ffsocket=ffsocket)
        mu.lammps_write_inputfile(folder=jobdir,filename='in.lmp',positions='data.lmp.runner',ace=ace)


        # get input-runner.xml (should be made without copying)
        stride = 100
        verbosity = "low"
        if test == True:
            nsteps = 5
            stride = 1
            verbosity = "high"

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
        mu.sed(jobdir+"/input-runner.xml",'<!-- <activelist> .*','<activelist> '+str(range(ncell**3 - nvac))+' </activelist>')
        mu.sed(jobdir+"/input-runner.xml",'<file mode="xyz" units="angstrom">.*</file>','<file mode="xyz" units="angstrom"> data_fakevac.ipi </file>')

        mu.sed(jobdir+"/input-runner.xml",'<ffsocket.*','<ffsocket name="lmpserial" mode="'+str(ffsocket)+'">')
        addressline = '<address> '+socket.gethostname()+' </address>'
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
