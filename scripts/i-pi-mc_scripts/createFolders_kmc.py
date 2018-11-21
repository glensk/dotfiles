#!/usr/bin/env python
import os,sys,random
#if sys.version_info[0] < 3:
#    raise Exception("Must be using Python 3")

import socket
from shutil import copyfile
import click

# from scripts folder
import convert_fileformats
import myutils as mu

# show default values in click
orig_init = click.core.Option.__init__
def new_init(self, *args, **kwargs):
    orig_init(self, *args, **kwargs)
    self.show_default = True
click.core.Option.__init__ = new_init


# get help also with -h
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
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
@click.option('-seednumber',default=False,multiple=True, type=int,help="define seed number manually (can be defined multiple times)")
@click.option('-nsteps',type=int, default=200000, help="number of KMC steps to make")
@click.option('-runnercutoff',type=float, default=10., help="runner cutoff distance ~10Angstrom")


# environment variables
@click.option('-scripts', envvar='scripts',help='environment variable $scripts (can alse be set here)')
@click.option('-nn_pot',type=str, default="v2dg", help="foldername in $scripts containing the neural network potential (can be separately set here for different path)")
@click.option('-ipi_mc', envvar='ipi_mc',help='path to i-pi-mc (or environment variable $ipi_mc)')
@click.option('-lmp_exec', envvar='lmp_exec',help='path to lammps executable (or environment variable $lmp_exec)')
@click.option('-submit/-no-submit', default=False)
@click.option('-submitdebug/-no-submitdebug', default=False)
@click.option('-t','--test/--no-test', default=False)

def createjob(
        ncell,
        nmg,
        nsi,
        nvac,
        a0,
        temp,
        scripts,
        nn_pot,
        nseeds,
        seednumber,
        nsteps,
        runnercutoff,
        ipi_mc,
        lmp_exec,
        submit,
        submitdebug,
        test):
    """
    This is an script to create KMC jobs quickly (and submit).

    e.g.

    createFolder_kmc.py -temp 1000 -ncell 5 -nsi 3 -nmg 3 -nvac 1 -submit
    """
    verbose = True
    ####################################
    # a few checks
    ####################################
    if verbose:
        print("a few checks")
    if submit is True or submitdebug is True and socket.gethostname() == "fidis":
        if socket.gethostname() != "fidis":    # use != and not: is not
            print('you are NOT on fidis')
            sys.exit('submit or submitdebug is True but you are no fidis! Exit.')

    nn_pot_dir              = scripts + "/potentials/runner_lammps_nn/" + nn_pot
    file_inlmp              = scripts + "/i-pi-mc_scripts/in.lmp"
    file_submit             = scripts + "/i-pi-mc_scripts/submit-ipi-kmc.sh"
    file_ipi_input_runner   = scripts + "/i-pi-mc_scripts/input-runner.xml"
    mu.check_isdir_or_isdirs([nn_pot_dir,scripts])
    mu.check_isfile_or_isfiles([file_inlmp,file_submit],["file_inlmp","file_submit"])
    if not test:
        mu.check_isfile_or_isfiles([ipi_mc,lmp_exec],["ipi_mc","lmp_exec"],envvar=True)


    ####################################
    # get directoryname
    ####################################
    if verbose:
        print("get directoryname")
    pcsi = nsi/ncell**3.*100
    pcmg = nmg/ncell**3.*100
    pcvac = nvac/ncell**3.*100
    directory = "ipi_kmc_"+\
                str(ncell)+"x"+str(ncell)+"x"+str(ncell)+"_"+nn_pot+"_"+\
                str(temp)+"K_"+\
                str(nvac)+"Vac_"+str(nmg)+"Mg_"+str(nsi)+"Si__"+\
                str(round(pcvac,3))+"pctVac_"+str(round(pcmg,3))+"pctMg_"+str(round(pcsi,3))+"pctSi_"+\
                str(runnercutoff)+"rcut"
    seeds = random.sample(range(1, 999999), nseeds)

    if test == True:
        nseeds = 1
        seeds = [1]
        print()
    seednumber = list(seednumber)
    if len(seednumber) is not 0:
        seeds = seednumber
        nseeds = len(seednumber)

    # make the atomic structure
    atomsc = mu.get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0)


    # show the input variables
    print('--------------------------- check the input --------------------------------')
    #print('seednumber   ',seednumber,type(seednumber))
    print('JOBS:        ',nseeds,'!! defined by -nseeds / or -seednumber')
    print()
    print('ncell        ',ncell,"(",atomsc.get_number_of_atoms(),"atoms )")
    print('nsi          ',nsi,  "(",pcsi,"%)")
    print('nmg          ',nmg,"(",pcmg,"%)")
    print('nvac         ',nvac,"(",pcvac,"%)")
    print('a0           ',a0)
    print('temp         ',temp)
    print()
    print('nn_pot       ',nn_pot)
    print('nn_pot_dir   ',nn_pot_dir)
    print()
    print('nseeds       ',nseeds,'seeds',seeds)
    print('nsteps       ',nsteps)
    print()
    print('directory    ',directory)
    print('submit       ',submit)
    print('submitdebug  ',submitdebug)
    print('pytohon      ',sys.version_info[0])
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
        convert_fileformats.save_ase_object_as_lmp_runner(atomsc,jobdir+'/data.lmp.runner')
        if test == True:
            convert_fileformats.save_ase_object_as_lmp(atomsc,jobdir+'/data.lmp')
            convert_fileformats.save_ase_object_in_ase_format(atomsc,jobdir+'/data.qe','espresso-in')
            convert_fileformats.save_ase_object_in_ase_format(atomsc,jobdir+'/data.POSCAR','vasp')
            convert_fileformats.save_ase_object_in_ase_format(atomsc,jobdir+'/data.xyz','xyz')
            convert_fileformats.save_ase_object_in_ase_format(atomsc,jobdir+'/data.extxyz','extxyz')
            # will only work once giulio send me the new runner.py file
            #convert_fileformats.save_ase_object_in_ase_format(atomsc,jobdir+'/data.runner','runner')

        # get and adapt in.lmp
        copyfile(file_inlmp, jobdir+"/in.lmp")
        mu.sed(jobdir+"/in.lmp",'variable runnerDir.*','variable runnerDir string "'+nn_pot_dir+'"')
        mu.sed(jobdir+"/in.lmp",'variable runnerCutoff.*','variable runnerCutoff equal '+str(runnercutoff))
        mu.sed(jobdir+"/in.lmp",'variable nameStartCfg .*','variable nameStartCfg string "data.lmp.runner"')
        mu.sed(jobdir+"/in.lmp",'variable initTemp.*','variable initTemp equal '+str(temp))
        mu.sed(jobdir+"/in.lmp",'variable startTemp.*','variable startTemp equal '+str(temp))
        mu.sed(jobdir+"/in.lmp",'variable stopTemp.*','variable stopTemp equal '+str(temp))
        mu.sed(jobdir+"/in.lmp",'variable numSteps.*','variable numSteps equal '+str(nsteps))


        # get submit-ipi-kmc.sh (could be made without copying)
        copyfile(file_ipi_input_runner, jobdir+"/input-runner.xml")
        mu.sed(jobdir+"/input-runner.xml",'<total_steps>.*</total_steps>','<total_steps> '+str(nsteps)+' </total_steps>')
        mu.sed(jobdir+"/input-runner.xml",'<seed>.*</seed>','<seed> '+str(seed)+' </seed>')
        mu.sed(jobdir+"/input-runner.xml",'<a0 units="angstrom">.*</a0>','<a0 units="angstrom"> '+str(a0)+' </a0>')
        mu.sed(jobdir+"/input-runner.xml",'<ncell>.*</ncell>','<ncell> '+str(ncell)+' </ncell>')
        mu.sed(jobdir+"/input-runner.xml",'<nsi>.*</nsi>','<nsi> '+str(nsi)+' </nsi>')
        mu.sed(jobdir+"/input-runner.xml",'<nmg>.*</nmg>','<nmg> '+str(nmg)+' </nmg>')
        mu.sed(jobdir+"/input-runner.xml",'<nvac>.*</nvac>','<nvac> '+str(nvac)+' </nvac>')
        mu.sed(jobdir+"/input-runner.xml",'<neval>.*</neval>','<neval> '+str(neval)+' </neval>')
        mu.sed(jobdir+"/input-runner.xml",'<temperature units="kelvin">.*','<temperature units="kelvin">'+str(temp)+'</temperature>')
        mu.sed(jobdir+"/input-runner.xml",'<file mode="xyz" units="angstrom">.*</file>','<file mode="xyz" units="angstrom"> '+str("data")+'.ipi </file>')

        # get submit-ipi-kmc.sh (could be made without copying)
        copyfile(file_submit, jobdir+"/submit-ipi-kmc.sh")
        mu.sed(jobdir+"/submit-ipi-kmc.sh",'#SBATCH --nodes=.*','#SBATCH --nodes='+str(nodes))
        mu.sed(jobdir+"/submit-ipi-kmc.sh",'#SBATCH --ntasks.*','#SBATCH --ntasks '+str(ntasks))
        mu.sed(jobdir+"/submit-ipi-kmc.sh",'--exclusive -n .* --mem','--exclusive -n '+str(lmp_par)+' --mem')
        mu.sed(jobdir+"/submit-ipi-kmc.sh",'--mem=4G .* < in.lmp','--mem=4G '+str(lmp_exec)+' < in.lmp')
        mu.sed(jobdir+"/submit-ipi-kmc.sh",'for i in `seq.*','for i in `seq '+str(ipi_inst)+'`')
        mu.sed(jobdir+"/submit-ipi-kmc.sh",'^python .* input-runner','python '+str(ipi_mc)+' input-runner')

        mu.submitjob(submit=submit,submitdebug=submitdebug,jobdir=jobdir,submitskript="submit-ipi-kmc.sh")

    print('done')
    return


if __name__ == "__main__":
    # submitoptions
    if True:
        nodes=2
        ipi_inst = 4
        lmp_par = 14

    if False:
        nodes=3
        ipi_inst = 7
        lmp_par = 12

    ntasks = cores = nodes * 28
    neval  = ipi_inst*2


    createjob()