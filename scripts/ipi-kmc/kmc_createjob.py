#!/usr/bin/env python
import os,sys
from shutil import copyfile
import numpy as np
from subprocess import call
scripts = os.getenv("HOME")+'/Dropbox/Albert/scripts/dotfiles/scripts/ipi-kmc/'

# import massedit in python 2.7
python2 = True
if python2:
    import imp
    masseditpath = scripts+'/massedit.py'
    print('me',masseditpath)
    massedit = imp.load_source('massedit', masseditpath)



ncell=8
nmg=6
nsi=6
nvac=1
a0=4.057

nseeds=4
nseeds=2 # for testing
seedplus=32346
steps=200000
steps=200 # for testing

if True:
    nodes=2
    ipi_inst = 4
    lmp_par = 14

if False:
    nodes=3
    ipi_inst = 7
    lmp_par = 12


######################### stop editing #############################
ntasks = cores = nodes * 28
neval = ipi_inst*2


def mkdir(directory):
    if not os.path.exists(directory):
            os.makedirs(directory)

def sed(file,str_find,str_replace):
    massedit.edit_files([file], ["re.sub('"+str_find+"', '"+str_replace+"', line)"],dry_run=False)

def get_positions_files(scripts,ncell,nmg,nsi,nvac,a0):
    sc=ncell
    alat=a0
    nva=nvac
    number_of_atoms = sc**3 - nva
    nal = number_of_atoms - nmg - nsi
    path=scripts+'/kmc_positions/'
    filename = "al"+str(sc)+"x"+str(sc)+"x"+str(sc)+"_alat"+str(alat)+"_"+str(nal)+"al_"+str(nsi)+"si_"+str(nmg)+"mg_"+str(nva)+"va_"+str(number_of_atoms)+"atoms"
    print('filename',filename)
    xyz = path+filename+'.xyz'
    ipi = path+filename+'.ipi'
    lmp = path+filename+'.xyz.lmp'
    if os.path.exists(xyz):
        return xyz,ipi,lmp,filename
    else:
        sys.exit('xyz files '+xyz+" not found!")




#def makejobs(ncell,nvac,nmg,nsi):
directory = str(ncell)+"x"+str(ncell)+"x"+str(ncell)+"_"+str(nvac)+"Vac_"+str(nmg)+"Mg_"+str(nsi)+"Si"
mkdir(directory)
seeds=np.arange(nseeds)+seedplus
xyz,ipi,lmp,positionsname = get_positions_files(scripts,ncell,nmg,nsi,nvac,a0)
print('xyz',xyz)
print('ipi',ipi)
print('lmp',lmp)
print('positionsname',positionsname)
for seed in seeds:
    jobdir = directory+'/seed'+str(seed)
    print('jobdir',jobdir)
    if os.path.exists(jobdir):
        sys.exit("jobdirectory "+str(jobdir)+" already exists!")
    mkdir(jobdir)
    copyfile(scripts+"input-runner.xml", jobdir+"/input-runner.xml")
    #copyfile(xyz, jobdir+"/"+positionsname+".xyz")
    copyfile(ipi, jobdir+"/"+positionsname+".ipi")
    copyfile(lmp, jobdir+"/data.lmp")
    copyfile(scripts+"in.lmp", jobdir+"/in.lmp")
    copyfile(scripts+"submit.sh", jobdir+"/submit.sh")
    sed(jobdir+"/input-runner.xml",'<total_steps>.*</total_steps>','<total_steps> '+str(steps)+' </total_steps>')
    sed(jobdir+"/input-runner.xml",'<seed>.*</seed>','<seed> '+str(seed)+' </seed>')
    sed(jobdir+"/input-runner.xml",'<a0 units="angstrom">.*</a0>','<a0 units="angstrom"> '+str(a0)+' </a0>')
    sed(jobdir+"/input-runner.xml",'<ncell>.*</ncell>','<ncell> '+str(ncell)+' </ncell>')
    sed(jobdir+"/input-runner.xml",'<nsi>.*</nsi>','<nsi> '+str(nsi)+' </nsi>')
    sed(jobdir+"/input-runner.xml",'<nmg>.*</nmg>','<nmg> '+str(nmg)+' </nmg>')
    sed(jobdir+"/input-runner.xml",'<nvac>.*</nvac>','<nvac> '+str(nvac)+' </nvac>')
    sed(jobdir+"/input-runner.xml",'<file mode="xyz" units="angstrom">.*</file>','<file mode="xyz" units="angstrom"> '+str(positionsname)+'.ipi </file>')
    sed(jobdir+"/input-runner.xml",'<neval>.*</neval>','<neval> '+str(neval)+' </neval>')

    sed(jobdir+"/submit.sh",'#SBATCH --nodes=.*','#SBATCH --nodes='+str(nodes))
    sed(jobdir+"/submit.sh",'#SBATCH --ntasks.*','#SBATCH --ntasks '+str(ntasks))
    sed(jobdir+"/submit.sh",'--exclusive -n .* --mem','--exclusive -n '+str(lmp_par)+' --mem')
    sed(jobdir+"/submit.sh",'for i in `seq.*','for i in `seq '+str(ipi_inst)+'`')

    cwd = os.getcwd()
    os.chdir(jobdir)
    call(["sbatch","submit.sh"])
    os.chdir(cwd)
