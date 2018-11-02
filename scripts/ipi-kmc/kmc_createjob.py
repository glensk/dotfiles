#!/usr/bin/env python
import os,sys,random,massedit,ase
from ase.lattice.cubic import FaceCenteredCubic
from shutil import copyfile
import numpy as np
from subprocess import call
import click

# from scripts folder
import convert_fileformats

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
@click.option('-a0'   ,default = 4.057, type=float, help="fcc lattice constant for Al")
@click.option('-temp' ,default = 300, type=int, help="KMC temperature")
@click.option('-nseeds',type=int, default=3, help="number of different seeds")
@click.option('-seednumber',default=False,multiple=True, type=int,help="define seed number manually (can be defined multiple times)")
@click.option('-nsteps',type=int, default=200000, help="number of KMC steps to make")
@click.option('-runnercutoff',type=float, default=10., help="runner cutoff distance ~10Angstrom")


# environment variables
@click.option('-scripts', envvar='scripts',help='path to scripts folder (or environment variable $scripts)')
@click.option('-nn_pot',type=str, default="v2dg", help="foldername for neural network potential")
@click.option('-i_pi_mc', envvar='i_pi_mc',help='path to i-pi-mc (or environment variable $i_pi_mc)')
@click.option('-lmp', envvar='lmp',help='path to lammps executable (or environment variable $lmp)')
@click.option('-submit/-no-submit', default=False)



def main(ncell, nmg, nsi,nvac,a0,temp,scripts,nn_pot,nseeds,seednumber,nsteps,runnercutoff,i_pi_mc,lmp,submit):
    """This is an script to submit KMC jobs quickly."""

    nn_pot_dir = scripts + "pot_nn/" + nn_pot
    file_inlmp = scripts + "ipi-kmc/in.lmp"
    file_submit = scripts + "ipi-kmc/submit-ipi-kmc.sh"
    file_ipi_input_runner = scripts + "ipi-kmc/input-runner.xml"
    check_isdir([nn_pot_dir,scripts])
    check_isfile(\
		[file_inlmp,file_submit,i_pi_mc,lmp],
		["file_inlmp","file_submit","i_pi_mc","lmp"])
    pcsi = nsi/ncell**3.*100
    pcmg = nmg/ncell**3.*100
    pcvac = nvac/ncell**3.*100
    directory = str(ncell)+"x"+str(ncell)+"x"+str(ncell)+"_"+nn_pot+"_"+\
                str(temp)+"K_"+\
                str(nvac)+"Vac_"+str(nmg)+"Mg_"+str(nsi)+"Si__"+\
                str(pcvac)+"pcVac_"+str(pcmg)+"pcMg_"+str(pcsi)+"pcSi_"+\
                str(runnercutoff)+"rcut"
    seeds = random.sample(range(1, 999999), nseeds)
    seednumber = list(seednumber)
    if len(seednumber) is not 0:
        seeds = seednumber
        nseeds = len(seednumber)


    print('--------------------------- check the input --------------------------------')
    #print('seednumber   ',seednumber,type(seednumber))
    print('JOBS:        ',nseeds,'!! defined by nseeds')
    print()
    print('ncell        ',ncell,"(",ncell**3,"atoms )")
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
    print('--------------------------- check the input --------------------------------')
    if os.path.isdir(directory):
        check_prompt("This directory exists already, shall I add jobs? [y]es: ")

    check_prompt("Are the ine input variables ok? [y]es: ")

    atomsc = get_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0)

    for seed in seeds:
        jobdir = directory+'/seed'+str(seed)
        print('jobdir',jobdir)
        if os.path.exists(jobdir):
            sys.exit("jobdirectory "+str(jobdir)+" already exists!")
        mkdir(jobdir)

        # get data.lmp
        convert_fileformats.save_ase_object_as_lmp_runner(atomsc,jobdir+'/data.lmp')
        convert_fileformats.save_ase_object_as_ipi_format(atomsc,jobdir+'/data.ipi')

        # get and adapt in.lmp
        copyfile(file_inlmp, jobdir+"/in.lmp")
        sed(jobdir+"/in.lmp",'variable runnerDir       string ".*','variable runnerDir       string "'+nn_pot_dir+'"')
        sed(jobdir+"/in.lmp",'variable initTemp.*','variable initTemp equal '+str(temp))
        sed(jobdir+"/in.lmp",'variable startTemp.*','variable startTemp equal '+str(temp))
        sed(jobdir+"/in.lmp",'variable stopTemp.*','variable stopTemp equal '+str(temp))
        sed(jobdir+"/in.lmp",'variable numSteps.*','variable numSteps equal '+str(nsteps))
        sed(jobdir+"/in.lmp",'variable runnerCutoff.*','variable runnerCutoff equal '+str(runnercutoff))


        # get submit-ipi-kmc.sh
        copyfile(file_ipi_input_runner, jobdir+"/input-runner.xml")
        sed(jobdir+"/input-runner.xml",'<total_steps>.*</total_steps>','<total_steps> '+str(nsteps)+' </total_steps>')
        sed(jobdir+"/input-runner.xml",'<seed>.*</seed>','<seed> '+str(seed)+' </seed>')
        sed(jobdir+"/input-runner.xml",'<a0 units="angstrom">.*</a0>','<a0 units="angstrom"> '+str(a0)+' </a0>')
        sed(jobdir+"/input-runner.xml",'<ncell>.*</ncell>','<ncell> '+str(ncell)+' </ncell>')
        sed(jobdir+"/input-runner.xml",'<nsi>.*</nsi>','<nsi> '+str(nsi)+' </nsi>')
        sed(jobdir+"/input-runner.xml",'<nmg>.*</nmg>','<nmg> '+str(nmg)+' </nmg>')
        sed(jobdir+"/input-runner.xml",'<nvac>.*</nvac>','<nvac> '+str(nvac)+' </nvac>')
        sed(jobdir+"/input-runner.xml",'<neval>.*</neval>','<neval> '+str(neval)+' </neval>')
        sed(jobdir+"/input-runner.xml",'<temperature units="kelvin">.*','<temperature units="kelvin">'+str(temp)+'</temperature>')
        sed(jobdir+"/input-runner.xml",'<file mode="xyz" units="angstrom">.*</file>','<file mode="xyz" units="angstrom"> '+str("data")+'.ipi </file>')

        # get submit-ipi-kmc.sh
        copyfile(file_submit, jobdir+"/submit-ipi-kmc.sh")
        sed(jobdir+"/submit-ipi-kmc.sh",'#SBATCH --nodes=.*','#SBATCH --nodes='+str(nodes))
        sed(jobdir+"/submit-ipi-kmc.sh",'#SBATCH --ntasks.*','#SBATCH --ntasks '+str(ntasks))
        sed(jobdir+"/submit-ipi-kmc.sh",'--exclusive -n .* --mem','--exclusive -n '+str(lmp_par)+' --mem')
        sed(jobdir+"/submit-ipi-kmc.sh",'--mem=4G .* < in.lmp','--mem=4G '+str(lmp)+' < in.lmp')
        sed(jobdir+"/submit-ipi-kmc.sh",'for i in `seq.*','for i in `seq '+str(ipi_inst)+'`')
        sed(jobdir+"/submit-ipi-kmc.sh",'^python .* input-runner','python '+str(i_pi_mc)+' input-runner')

        #with open(jobdir+"/submit-ipi-kmc.sh", 'r') as fin:
        #    print(fin.read())
        #sys.exit()

        if submit is True:
            cwd = os.getcwd()
            os.chdir(jobdir)
            call(["sbatch","submit-ipi-kmc.sh"])
            os.chdir(cwd)





def sed(file,str_find,str_replace):
    massedit.edit_files([file], ["re.sub('"+str_find+"', '"+str_replace+"', line)"],dry_run=False)


def get_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0):
    atom = ase.build.bulk('Al',crystalstructure='fcc',a=a0)
    atomsc = atom.repeat(ncell)
    number_of_atoms = atomsc.get_number_of_atoms()
    nal = number_of_atoms - nsi - nmg
    for i in np.arange(nsi):
        atomsc[i].symbol = 'Si'
    for i in np.arange(nsi,nmg+nsi):
        atomsc[i].symbol = 'Mg'
    for i in np.arange(nvac):
        del atomsc[-1]
    number_of_atoms = atomsc.get_number_of_atoms()
    nal = number_of_atoms - nsi - nmg
    return atomsc


def check_isdir(path):
    if type(path) is str:
        if not os.path.isdir(path):
            sys.exit('missing directory '+path)

    if type(path) is list:
        for i in path:
            if not os.path.isdir(i):
                sys.exit('missing directory '+i)
    return

def check_isfile(path,pathnames):
    if type(path) is str:
        if not os.path.isfile(path):
            sys.exit('missing file '+path)
    elif type(path) is list:
        for idx,i in enumerate(path):
            #print('ii',i,idx,pathnames[idx])
            if i is None:
                sys.exit('(environment) variable "'+pathnames[idx]+'" is not defined!')
            if not os.path.isfile(i):
                sys.exit('missing file '+i)
    else:
        sys.exit('unknown type file type path '+str(path))
    return

def check_prompt(check):
    #isok = raw_input(check)  # Python 2
    isok = input(check) # Python 3
    if isok == "":
        sys.exit()
    if isok[0] not in [ 'Y', 'y' ]:
        sys.exit('Exist since not Y or y as first letter!')
    return

#def get_positions_files(scriptsipi,ncell,nmg,nsi,nvac,a0):
#    sc=ncell
#    alat=a0
#    nva=nvac
#    number_of_atoms = sc**3 - nva
#    nal = number_of_atoms - nmg - nsi
#    path=scriptsipi+'/kmc_positions/'
#    filename = "al"+str(sc)+"x"+str(sc)+"x"+str(sc)+"_alat"+str(alat)+"_"+str(nal)+"al_"+str(nsi)+"si_"+str(nmg)+"mg_"+str(nva)+"va_"+str(number_of_atoms)+"atoms"
#    print('filename',filename)
#    xyz = path+filename+'.xyz'
#    ipi = path+filename+'.ipi'
#    lmp = path+filename+'.xyz.lmp'
#    if os.path.exists(xyz):
#        return xyz,ipi,lmp,filename
#    else:
#        sys.exit('xyz files '+xyz+" not found!")
#

def mkdir(directory):
    if not os.path.exists(directory):
            os.makedirs(directory)


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



    main()

#    sys.exit()
# sys.exit()
#
#
#
# # kmc_make_xyz_primitive_cell.py  -sc 6 -nmg 6 -nsi 6 -nvac 2
# # kmc_make_xyz_primitive_cell.py  -sc 8 -nmg 6 -nsi 6 -nvac 2
# # kmc_make_xyz_primitive_cell.py  -sc 10 -nmg 6 -nsi 6 -nvac 2
# # python ./kmc_make_xyz_file_to_lammpsinput.py al6x6x6_alat4.057_202al_6si_6mg_2va_214atoms.xyz
# ######################### stop editing #############################
#
#
#
# # import massedit in python 2.7
# python2 = True
# if python2:
#     masseditpath = os.environ['scripts'] + 'ipi-kmc/massedit.py'
#     if not os.path.isfile(masseditpath):
#         sys.exit(masseditpath+' is missing')
#     import imp
#     massedit = imp.load_source('massedit', masseditpath)
#
#
#
#
#
#
#
#
#
#
# directory = str(ncell)+"x"+str(ncell)+"x"+str(ncell)+"_"+str(nvac)+"Vac_"+str(nmg)+"Mg_"+str(nsi)+"Si"+add_to_name
# mkdir(directory)
# seeds=np.arange(nseeds)+seedplus
# xyz,ipi,lmp,positionsname = get_positions_files(scriptsipi,ncell,nmg,nsi,nvac,a0)
# print('xyz',xyz)
# print('ipi',ipi)
# print('lmp',lmp)
# print('positionsname',positionsname)
# print()
