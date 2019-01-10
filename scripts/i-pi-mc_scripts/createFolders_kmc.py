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

# environment variables
@click.option('--pot','-p',type=click.Choice(mu.pot_all()),required=True,default='n2p2_v1ag',help="potential from $scripts/potentials folder.")

@click.option('-submit/-no-submit', default=False)
@click.option('-submitdebug/-no-submitdebug', default=False)
@click.option('-st','--submittime_hours', type=int,default=71,help="slurm time for the job")
@click.option('--ffsocket',type=click.Choice(["unix","inet"]),default='inet',help='ipi fftsocket "unix" or "inet"')
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
        pot,
        submit,
        submitdebug,
        submittime_hours,
        test,
        ffsocket,
        verbose):
    """
    This is an script to create KMC jobs quickly (and submit).

    e.g.

    createFolder_kmc.py -temp 1000 -ncell 5 -nsi 3 -nmg 3 -nvac 1 -submit
    """
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


    file_inlmp              = scripts + "/i-pi-mc_scripts/in.lmp"
    file_submit             = scripts + "/i-pi-mc_scripts/submit-ipi-kmc.sh"
    file_ipi_input_runner   = scripts + "/i-pi-mc_scripts/input-runner.xml"
    mu.check_isfile_or_isfiles([file_inlmp,file_submit],["file_inlmp","file_submit"])
    LAMMPS_COMMAND = mu.test_and_return_environment_var_path('LAMMPS_COMMAND')
    IPI_COMMAND    = mu.test_and_return_environment_var_path('IPI_COMMAND')


    ####################################
    # get directoryname
    ####################################
    if verbose:
        print("get directoryname")
    pcsi = nsi/ncell**3.*100
    pcmg = nmg/ncell**3.*100
    pcvac = nvac/ncell**3.*100
    directory = "ipi_kmc_"+\
                str(ncell)+"x"+str(ncell)+"x"+str(ncell)+"_"+pot+"_"+\
                str(temp)+"K_"+\
                str(nvac)+"Vac_"+str(nmg)+"Mg_"+str(nsi)+"Si__"+\
                str(round(pcvac,3))+"pctVac_"+str(round(pcmg,3))+"pctMg_"+str(round(pcsi,3))+"pctSi_"+\
                str(runnercutoff)+"rcut"


    # make the atomic structure
    atomsc = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0)


    # show the input variables
    print('--------------------------- check the input --------------------------------')
    #print('seednumber   ',seednumber,type(seednumber))
    print('JOBS:         ',nseeds,'!! defined by -nseeds / or -seednumber')
    print()
    print('ncell         ',ncell,"(",atomsc.get_number_of_atoms(),"atoms )")
    print('nsi           ',nsi,  "(",pcsi,"%)")
    print('nmg           ',nmg,"(",pcmg,"%)")
    print('nvac          ',nvac,"(",pcvac,"%)")
    print('a0            ',a0)
    print('temp          ',temp)
    print()
    print('mypot.pot     ',mypot.pot)
    print('mypot.fullpath',mypot.fullpath)
    print()
    print('nseeds        ',nseeds,'seeds',seeds)
    print('nsteps        ',nsteps)
    print()
    print('directory     ',directory)
    print('submit        ',submit)
    print('submitdebug   ',submitdebug)
    #print('python ver    ',sys.version_info[0])
    #print()
    #print('LAMMPS_COMMAND',LAMMPS_COMMAND)
    #print('IPI_COMMAND   ',IPI_COMMAND)
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

        # get data.lmp
        convert_fileformats.save_ase_object_as_ipi_format(atomsc,jobdir+'/data.ipi')

        atomsc.write(jobdir+'/data.lmp.runner',format='lammps-runner')

        #if test == True:
        #    atomsc.write(jobdir+'/data.lmp',format='lammps-data')
        #    atomsc.write(jobdir+'/data.POSCAR',format='vasp')
        #    atomsc.write(jobdir+'/data.xyz',format='xyz')
        #    atomsc.write(jobdir+'/data.extxyz',format='extxyz')
        #    atomsc.write(jobdir+'/data.espresso-in',format='espresso-in')


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
        mu.sed(jobdir+"/input-runner.xml",'<file mode="xyz" units="angstrom">.*</file>','<file mode="xyz" units="angstrom"> data.ipi </file>')

        mu.sed(jobdir+"/input-runner.xml",'<ffsocket.*','<ffsocket name="lmpserial" mode="'+str(ffsocket)+'">')
        addressline = '<address> '+socket.gethostname()+' </address>'
        if ffsocket == "unix":
            mu.sed(jobdir+"/input-runner.xml",'<address.*',addressline)
        if ffsocket == "inet":
            mu.sed(jobdir+"/input-runner.xml",'<address.*',addressline+' <port> 12345 </port>')



        # get submit-ipi-kmc.sh (should be made without copying)
        mu.create_submitskript_ipi_kmc(jobdir+"/submit-ipi-kmc.sh",nodes,ntasks,IPI_COMMAND,LAMMPS_COMMAND,lmp_par,ipi_inst,ffsocket,submittime_hours=submittime_hours)

        # submit the job
        mu.submitjob(submit=submit,submitdebug=submitdebug,jobdir=jobdir,submitskript="submit-ipi-kmc.sh")

    print('done')
    return


if __name__ == "__main__":
    # submitoptions
    if False:
        nodes=2
        ipi_inst = 4
        lmp_par = 14

    if False:
        nodes=3
        ipi_inst = 7
        lmp_par = 12

    # currently on fidis with parallel n2p2 only one node works using unix
    if False:
        nodes=1
        ipi_inst = 1
        lmp_par = 28

    if True:
        nodes=2
        ipi_inst = 4
        lmp_par = 14

    ntasks = cores = nodes * 28
    neval  = ipi_inst*2


    createjob()
