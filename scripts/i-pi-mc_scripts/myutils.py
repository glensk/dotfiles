#!/usr/bin/env python
import os,sys
import numpy as np
from subprocess import check_output
import datetime
from ase.build import bulk as ase_build_bulk


# from scripts folder
import massedit


def create_READMEtxt(directory,add=False):
    ''' wiretes a README.txt file '''
    # get sha
    hier = os.getcwd()
    os.chdir(os.environ['scripts'])
    sha = check_output(["git","rev-parse","master"]).decode('utf-8')
    os.chdir(hier)

    # get time
    time_now = datetime.datetime.now()

    # name of RADME
    filepath = directory+'/README_'+time_now.strftime("%Y-%m-%d_%H:%M")+'.txt'

    # write README.txt
    strout=os.path.basename(sys.argv[0])+" "+" ".join(sys.argv[1:])
    with open(filepath, "w") as text_file:
        text_file.write("# using https://github.com/glensk/dotfiles/trunk/scripts")
        #print("# using https://github.com/glensk/dotfiles/trunk/scripts",file=text_file)
        text_file.write("# to download it: svn checkout https://github.com/glensk/dotfiles/trunk/scripts")
        #print("# to download it: svn checkout https://github.com/glensk/dotfiles/trunk/scripts",file=text_file)
        #print("# used sha: "+sha, file=text_file)
        #print(f"{strout}", file=text_file)
        if add:
            text_file.write(add)
            #print(add, file=text_file)


    print()
    print('written ',filepath)
    return

def submitjob(submit=False,submitdebug=False,jobdir=False,submitskript=False):
    if submit is True or submitdebug is True:
        check_isdir(jobdir)
        cwd = os.getcwd()
        os.chdir(jobdir)
        if submitdebug is True:  # this works on fidis even with 2 nodes!
            call(["sbatch","-p","debug","-t","01:00:00",submitskript])
        if submit is True:
            call(["sbatch",submitskript])
        os.chdir(cwd)
        return

def sed(file,str_find,str_replace):
    massedit.edit_files([file], ["re.sub('"+str_find+"', '"+str_replace+"', line)"],dry_run=False)
    return

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def check_isdir_or_isdirs(path):
    if type(path) is str:
        if not os.path.isdir(path):
            sys.exit('missing directory '+path)

    if type(path) is list:
        for i in path:
            if not os.path.isdir(i):
                sys.exit('missing directory '+i)
    return


def check_isfile_or_isfiles(path,pathnames,envvar=False):
    add=""
    if envvar == True: add = "(set as environment variable)"
    if type(path) is str:
        if not os.path.isfile(path):
            sys.exit('missing file '+path+add)
    elif type(path) is list:
        for idx,i in enumerate(path):
            #print('ii',i,idx,pathnames[idx])
            if i is None:
                sys.exit(pathnames[idx]+" "+add+' is not defined!')
            if not os.path.isfile(i):
                sys.exit('missing file '+i)
    else:
        sys.exit('unknown type file type path '+str(path))
    return


def get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0):
    atom = ase_build_bulk('Al',crystalstructure='fcc',a=a0)
    atomsc = atom.repeat(ncell)
    number_of_atoms = atomsc.get_number_of_atoms()
    nal = number_of_atoms - nsi - nmg

    for i in np.arange(nmg):
        atomsc[i].symbol = 'Mg'
    for i in np.arange(nmg,nmg+nsi):
        atomsc[i].symbol = 'Si'
    for i in np.arange(nvac):
        del atomsc[-1]
    number_of_atoms = atomsc.get_number_of_atoms()
    nal = number_of_atoms - nsi - nmg
    #ase.io.write('kalmp',atomsc,format='lammps-dump')
    return atomsc

def exit_if_not_python3(exit=False):
    python_version = sys.version_info[0]
    if python_version < 3 and exit:
        sys.exit('Your python environment uses a python < 3; Exit;')

def get_prompt_irrespective_of_python_version(text):
    python_version = sys.version_info[0]
    if python_version < 3:
        fromprompt = raw_input(text)
    else:
        fromprompt = input(text)
    return fromprompt

    ion_or_filename=rint(python_version,type(python_version))
    return


def get_from_prompt_Yy_orexit(text):
    getfromprompt = mu.get_prompt_irrespective_of_python_version(text)
    if getfromprompt == "":
        sys.exit()
    if getfromprompt[0] not in [ 'Y', 'y' ]:
        sys.exit('Exist since not Y or y as first letter!')
    return


def findfiles(directory=False,begin="",contains="",extension_or_filename=""):
    listout=[]
    # for python3
    #for filename in glob.glob(directory+'/**/*'+extension_or_filename, recursive=True):
    #    #print(filename)
    #    listout.append(filename)

    # for python2 and python 3
    import fnmatch
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, begin+'*'+extension_or_filename):
            listout.append(os.path.join(root, filename))
    return listout




if __name__ == "__main__":
    pass
