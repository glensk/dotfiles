#!/usr/bin/env python
from __future__ import print_function
import os,sys
import numpy as np
from subprocess import check_output,call
import datetime
from ase.build import bulk as ase_build_bulk
from socket import gethostname


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
        text_file.write("# using https://github.com/glensk/dotfiles/trunk/scripts")
        text_file.write("# to download it: svn checkout https://github.com/glensk/dotfiles/trunk/scripts")
        text_file.write("# to download it: svn checkout https://github.com/glensk/dotfiles/trunk/scripts")
        text_file.write("# used sha: "+sha)
        text_file.write(strout)
        if add:
            if type(add) == str:
                text_file.write(add)
            elif type(add) == list:
                for i in add:
                    text_file.write(i)


    print()
    print('written ',filepath)
    return

def submitjob(submit=False,submitdebug=False,jobdir=False,submitskript=False):
    if submit is True or submitdebug is True:
        check_isdir_or_isdirs(jobdir)
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


def get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0,cubic=False):
    """Creating bulk systems.

        Crystal structure and lattice constant(s) will be guessed if not
        provided.

        name: str
            Chemical symbol or symbols as in 'MgO' or 'NaCl'.
        crystalstructure: str
            Must be one of sc, fcc, bcc, hcp, diamond, zincblende,
            rocksalt, cesiumchloride, fluorite or wurtzite.
        a: float
            Lattice constant.
        c: float
            Lattice constant.
        covera: float
            c/a ratio used for hcp.  Default is ideal ratio: sqrt(8/3).
        u: float
            Internal coordinate for Wurtzite structure.
        orthorhombic: bool
            Construct orthorhombic unit cell instead of primitive cell
            which is the default.
        cubic: bool
            Construct cubic unit cell if possible.
    """
    atom = ase_build_bulk('Al',crystalstructure='fcc',a=a0,cubic=False)
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
    getfromprompt = get_prompt_irrespective_of_python_version(text)
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

def get_kmesh_size_daniel(ase_structure, kmesh_l):
    reci_cell = ase_structure.get_reciprocal_cell()
    #print("ase_structure.get_cell()",ase_structure.get_cell())
    #print("ase_structure.get_reciprocal_cell()",ase_structure.get_reciprocal_cell())
    kmesh = [np.ceil(kmesh_l * np.linalg.norm(reci_cell[i]))
             for i in range(len(reci_cell))]
    return kmesh

def qe_parse_numelectrons_upfpath(upfpath):
    ''' parses the number of electrons from a quantum espresso potential '''
    for line in open(upfpath):
        if "valence" in line.lower() and "z" in line.lower():
            if len(line.split("=")) == 2:
                num_e = int(float((line.split("=")[-1].strip().strip('"'))))
            elif len(line.split()) == 3:
                num_e = int(float(line.split()[0]))
            else:
                raise Exception("Could not parse {}".format(upfpath))
    return num_e

def qe_full_get_path_to_potential(element,path_to_pseudos):
    potential = findfiles(path_to_pseudos,begin=element,extension_or_filename="UPF")
    if len(potential) != 1:
        print('found potential',potential)
        sys.exit('found more than one or none potential; Exit.')
    return potential


def qe_get_numelectrons(structure_ase, path_to_pseudos):
    element_nume_dict = {}
    for element in structure_ase.get_chemical_symbols():
        if element in element_nume_dict:
            continue
        if element not in element_nume_dict:
            potential = findfiles(path_to_pseudos,begin=element,extension_or_filename="UPF")
            if len(potential) != 1:
                print('found potential',potential)
                sys.exit('found more than one or none potential; Exit.')
            #pot = os.path.basename(potential[0])
            num_elec = qe_parse_numelectrons_upfpath(potential[0])
            element_nume_dict[element] = num_elec
    number_electrons = 0
    for element in structure_ase.get_chemical_symbols():
        number_electrons+= element_nume_dict[element]
    return number_electrons


def get_inputfile_runner(template,filename,
        symfun_old_delete=True,symfun_file=False,
        test_fraction=0,runner_mode=1,number_of_elements=3,
        elements="Al Mg Si",test_input_data=True):
    ''' this creates the runner inputfile
    '''
    if test_input_data:
        len = file_len_linecount(test_input_data)
        #print('len',len)
        if len <= 3:sys.exit('file '+test_input_data+' seems too short! Exit;')

    # read in the runner.in template
    f = open(template,"r")
    lines = f.readlines()
    f.close()

    if symfun_old_delete == True:
        listdelete = []

        # delete the old symmetry functioins
        for idx,line in enumerate(lines):
            #print()
            #print('idx',idx,line)
            #print('idk',idx,line[:18])
            if line[:18] == "symfunction_short ":
                listdelete.append(idx)
            if line[:24] == "# symfunctions for type ":
                listdelete.append(idx)

        for i in np.array(listdelete)[::-1]:
            del lines[i]

        # insert the new symmetry functions
        if symfun_file != False:
            s = open(symfun_file,"r")
            sym = s.readlines()
            s.close()
            for idj,symline in enumerate(np.array(sym)[::-1]):
                #print('sl',symline)
                lines.insert(listdelete[0],symline)

        # set other options
        print('test_fraction        :',test_fraction)
        print('runner_mode          :',runner_mode)
        print('number_of_elements   :',number_of_elements)
        print('elements             :',elements)
        for idx,line in enumerate(lines):
            if line[:14] == "test_fraction ":
                lines[idx] = "test_fraction "+str(test_fraction)+"\n"
            if line[:12] == "runner_mode ":
                lines[idx] = "runner_mode "+str(runner_mode)+"\n"
            if line[:19] == "number_of_elements ":
                lines[idx] = "number_of_elements "+str(number_of_elements)+"\n"
            if line[:9] == "elements ": lines[idx] = "elements "+str(elements)+"\n"


        # write the file
        f = open(filename,"w")
        f.writelines(lines)
        f.close()
        print('written '+filename)
        return

def file_len_linecount(fname):
    i = 0
    with open(fname) as f:
        for i, l in enumerate(f,1):
            pass
    return i


def scripts():
    ''' return environment variable scripts '''
    scripts = os.environ['scripts']
    if not os.path.isdir('scripts'):
        sys.exit('$scripts variable is not defined or is not an existing folder')
    return scripts

def runner_exec():
    ''' return environment variable runner_exec (RuNNer executable)'''
    runner_exec = os.environ['runner_exec']
    if not os.path.isfile('runner_exec'):
        sys.exit('$runner_exec variable is not defined or is not an existing file')
    return runner_exec

def hostname():
    hostname = getgethostname()
    return hostname

if __name__ == "__main__":
    pass
