#!/usr/bin/env python
from __future__ import print_function
import os,sys,re
import click
import numpy as np
import glob #,pathlib
from tqdm import tqdm_notebook, tnrange
from copy import deepcopy
from socket import gethostname
import shutil
from subprocess import check_output,call
from datetime import datetime as datetime   # datetime.datetime.now()
import ase
from ase.build import bulk as ase_build_bulk
from ase.constraints import StrainFilter
from ase.neighborlist import NeighborList, neighbor_list, NewPrimitiveNeighborList
import hesse
try:
    import quippy
except ImportError:
    pass

try:  # not in aiida ase
    from ase.constraints import ExpCellFilter
except ImportError:
    pass

from ase.spacegroup import crystal
from ase.constraints import StrainFilter
try:
    from ase.constraints import ExpCellFilter
except ImportError:
    pass

try:
    from ase.calculators.lammpslib import LAMMPSlib
except ImportError:
    print("ERROR when importing LAMMPSlib ... possibly you have to change your (conda/aiida) environment")

from feos import eos
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.optimize import BFGS
from ase.optimize import LBFGS
from ase.optimize import FIRE
import my_atom
try:
    from ase.optimize import GPMin
except ImportError:
    pass
from ase.optimize.basin import BasinHopping
from ase.optimize.minimahopping import MinimaHopping
import time
import myutils as my
from ase import units as aseunits

start_time = time.time()

def printoutcolor(red,var,ENDC):
    if len(var) == 1:
        return red + str(var[0]) + ENDC
    else:
        return red + str(var) + ENDC



def printred(*var):
    red = '\033[31m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)

def printgreen(*var):
    red = '\033[32m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)

def printorange(*var):
    red = '\033[33m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)

def printblue(*var):
    red = '\033[34m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)

def printpink(*var):
    red = '\033[35m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)

def printturquise(*var):
    red = '\033[36m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)

def printwhite(*var): # stechend
    red = '\033[37m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)

def printgrey(*var):
    red = '\033[38m'
    ENDC = '\033[0m'
    return printoutcolor(red,var,ENDC)

def grep(filepath,string):
    out = []
    file = open(filepath, "r")

    for line in file:
         if re.search(string, line):
            line_ = line.rstrip()
            #print('found:'+line_+":",type(line_))
            out.append(line_)
            #print("foundo",out)
    #print("out",out)
    return out

def create_READMEtxt(directory=False,add=False):
    ''' wiretes a README.txt file '''
    if directory == False:
        directory = os.getcwd()
    # get sha
    pwd = os.getcwd()
    os.chdir(os.environ['scripts'])
    sha = check_output(["git","rev-parse","master"]).decode('utf-8')
    os.chdir(pwd)

    # get time
    time_now = datetime.now()

    # name of RADME
    filepath = directory+'/README_'+time_now.strftime("%Y-%m-%d_%H:%M:%S")+'.txt'

    # write README.txt
    strout=os.path.basename(sys.argv[0])+" "+" ".join(sys.argv[1:])
    with open(filepath, "w") as text_file:
        text_file.write("# using https://github.com/glensk/dotfiles/trunk/scripts\n")
        text_file.write("# to download it: svn checkout https://github.com/glensk/dotfiles/trunk/scripts\n")
        text_file.write("# used sha: "+sha) #+"\n")
        text_file.write("# execution time: "+str(time.time() - start_time)+" seconds.")
        text_file.write("\n")
        if add:
            print('add')
            print(type(add))
            if type(add) == str:
                text_file.write(add+"\n")
            elif type(add) == list:
                for i in add:
                    text_file.write(i+"\n")
        text_file.write("\n")
        text_file.write(strout+"\n")

    print()
    print('written ',filepath)
    print()
    return

def n2p2_get_scaling_and_function_data(submitdebug=True,cores=21,debug=True):
    if not os.path.isfile("input.data"):
        sys.exit("Need input.data file")
    if not os.path.isfile("input.nn"):
        sys.exit("Need input.nn file")

    #submitfile = scripts()+"/n2p2/submit_scaling_debug.sh"
    #if not os.path.isfile(submitfile):
    #    sys.exit("Need "+submitfile+" file!")

    folder="get_scaling"
    if os.path.isdir(folder):
        sys.exit(folder+" already exists!")


    mkdir(folder)
    create_READMEtxt()
    os.chdir(folder)
    cp("../input.data")
    cp("../input.nn")
    #cp(submitfile)
    submitskript = n2p2_write_submit_skript(directory=False,cores=cores,nodes=1,debugque=submitdebug,job="scaling")

    command = ["sbatch"]
    if submitdebug == True:
        command = command + ["-p","debug","-t","01:00:00"]
    command = command + [submitskript]
    call(["sbatch",submitskript])
    #submitjob(submitdebug=submitdebug,jobdir=os.getcwd(),submitskript=submitskript,cores=cores)
    create_READMEtxt(add="submitdebug = "+str(submitdebug))
    return

def n2p2_write_submit_skript(directory=False,nodes=1,cores=28,debugque=False,job=False):
    ''' wiretes a submit_n2p2_{get_scaling,training}.sh file '''
    if job not in ["scaling","train"]:
        sys.exit('job has to be one of nnp-XXX jobs as "train, scaling, ..."')
    hostname = gethostname()
    known_hosts =  ['fidis','helvetios']
    if hostname in known_hosts:
        pass
    else:
        print("known hosts:",known_hosts)
        sys.exi(hostname+" is not in the list of known hosts!")

    if directory == False:
        directory = os.getcwd()

    # name of file
    filepath = directory+'/submit_n2p2_'+job+'.sh'

    # write file
    with open(filepath, "w") as text_file:
        text_file.write("#!/bin/bash\n")
        text_file.write("#SBATCH --job-name=NNP-mpi\n")
        text_file.write("#SBATCH --output=_scheduler-stdout.txt\n")
        text_file.write("#SBATCH --error=_scheduler-stderr.txt\n")
        text_file.write("#SBATCH --nodes="+str(nodes)+"\n")
        text_file.write("#SBATCH --ntasks 28\n")
        if debugque == True:
            text_file.write("#SBATCH --time=00-01:00:00\n")
        else:
            text_file.write("#SBATCH --time=00-72:00:00\n")
        if hostname == 'fidis':
            text_file.write("#SBATCH --constraint=E5v4\n")  # means to only use the fidis nodes
            # to use the Gacrux/Skylake nodes: #SBATCH --constraint=s6g1
        if hostname == 'helvetios':
            print('you can use up to --mem=183G or more')
        text_file.write("#SBATCH --mem=100G\n")
        text_file.write("\n")
        text_file.write("set +e\n")
        text_file.write("# it is necessary to have all the modules which are used when compiling\n")
        text_file.write('export LD_LIBRARY_PATH=""\n')
        text_file.write("module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen\n")
        text_file.write("export LD_LIBRARY_PATH=$HOME/sources/n2p2/lib:${LD_LIBRARY_PATH}\n")
        text_file.write("#echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH\n")
        text_file.write("\n")
        text_file.write("touch time.out\n")
        text_file.write("date +%s >> time.out\n")
        text_file.write("\n")
        if job == 'scaling':
            text_file.write("srun -n "+str(cores)+" $HOME/sources/n2p2/bin/nnp-scaling 1\n")
        elif job == 'train':
            text_file.write("srun -n "+str(cores)+" $HOME/sources/n2p2/bin/nnp-train\n")
        text_file.write("date +%s >> time.out\n")
        text_file.write("cat time.out | xargs | awk '{print $2-$1-10}' > time.sec\n")
        text_file.write("$dotfiles/scripts/n2p2/n2p2_tarfolder_for_scale_train.sh\n")
        if job == 'scaling':
            text_file.write('[ "`pwd | grep -o "/get_scaling$"`" == \'/get_scaling\' ] && echo creating link && ln -s `pwd`/function.data ../function.data')
        text_file.write("\n")
        text_file.write("exit 0\n")

    print()
    print('written ',filepath)
    print()
    return filepath


def n2p2_make_training(cores=21,debugque=False):
    if not os.path.isfile("input.data"):
        sys.exit("Need input.data file")
    if not os.path.isfile("input.nn"):
        sys.exit("Need input.nn file")
    if not os.path.isfile("function.data"):
        sys.exit("Need function.data file")

    #submitfile = scripts()+"/n2p2/submit_training.sh"
    #if not os.path.isfile(submitfile):
    #    sys.exit("Need "+submitfile+" file!")

    if not os.path.isfile("scaling.data"):
        if not os.path.isfile("get_scaling/scaling.data"):
            sys.exit("get_scaling/scaling.data missing")
        else:
            cp("get_scaling/scaling.data","scaling.data")

    #cp(submitfile)
    submitskript = n2p2_write_submit_skript(directory=False,cores=cores,nodes=1,debugque=debugque,job="train")
    command = ["sbatch"]
    if debugque == True:
        command = command + ["-p","debug","-t","01:00:00"]
    command = command + [submitskript]
    call(["sbatch",submitskript])
    #submitjob(submitdebug=submitdebug,jobdir=os.getcwd(),submitskript=submitskript,cores=cores)

    #submitjob(submitdebug=False,submit=True,jobdir=os.getcwd(),submitskript="submit_training.sh",cores=cores)
    create_READMEtxt()
    return

#def submitjob(submit=False,submitdebug=False,jobdir=False,submitskript=False,cores=21):
#    if jobdir == False:
#        jobdir = os.getcwd()
#    if submit is True or submitdebug is True:
#        check_isdir_or_isdirs(jobdir)
#        cwd = os.getcwd()
#        os.chdir(jobdir)
#        if submitdebug is True:  # this works on fidis even with 2 nodes!
#            sed(submitskript,"srun -n.*","srun -n "+str(cores)+" $HOME/sources/n2p2/bin/nnp-train")
#            call(["sbatch","-p","debug","-t","01:00:00",submitskript])
#        if submit is True:
#            sed(submitskript,"srun -n.*","srun -n "+str(cores)+" $HOME/sources/n2p2/bin/nnp-train")
#            call(["sbatch",submitskript])
#        os.chdir(cwd)
#        return

def sed(file,str_find,str_replace):
    # from scripts folder
    import massedit
    massedit.edit_files([file], ["re.sub('"+str_find+"', '"+str_replace+"', line)"],dry_run=False)
    return

def get_absdir_from_relative_filepath(filepath):
    print('fp',filepath)
    return


def cp(src=False,dest=False):
    if dest==False:
        dest = os.getcwd()
    if dest==".":
        dest = os.getcwd()
    src_is = False
    dest_is = False
    #s = pathlib.Path(src)
    #d = pathlib.Path(dest)
    #print('src',src,"file?",s.is_file())
    #print('src',src,"dir? ",s.is_dir())
    #print('dest',dest,"file?",d.is_file())
    #print('dest',dest,"dir ?",d.is_dir())
    if os.path.isfile(src):
        #print("src is file")
        src_is = "file"
    if os.path.isdir(src):
        #print("src is dir")
        src_is = "dir"
    if os.path.isfile(dest):
        #print("dest is file")
        dest_is = "file"
    if os.path.isdir(dest):
        #print("dest is dir")
        dest_is = "dir"
    if src_is == "file" and dest_is == False:
        # possibly because dest does not exist
        shutil.copyfile(src,dest)
    elif src_is == "file" and dest_is == "file":
        shutil.copyfile(src,dest)
    elif src_is == "file" and dest_is == "dir":
        basename = os.path.basename(src)
        #print("basename:",os.path.basename(src))
        shutil.copyfile(src,dest+"/"+basename)
    else:
        print("source is:",src_is,":",src)
        print("dest   is:",dest_is,":",dest)
        basename = os.path.basename(src)
        print("basename:",os.path.basename(src))
        to = dest+'/'+basename
        print('to',to)
        def copyDirectory(src, dest):
            try:
                shutil.copytree(src, dest)
            # Directories are the same
            except shutil.Error as e:
                print('Directory not copied. Error: %s' % e)
            # Any error saying that the directory doesn't exist
            except OSError as e:
                print('Directory not copied. Error: %s' % e)
        print('copy...')
        copyDirectory(src, to)
    return

def rm(src):
    os.remove(src)

def rm_if_exists(src):
    if os.path.isfile(src):
        os.remove(src)

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

class cd:
    """
    Context manager for changing the current working directory
    use:
    import myutils as my
    with my.cd("~/Library"):
        subprocess.call("ls")

    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def check_isdir_or_isdirs(path):
    if type(path) is str:
        if not os.path.isdir(path):
            sys.exit('missing directory '+path)

    if type(path) is list:
        for i in path:
            if not os.path.isdir(i):
                sys.exit('missing directory '+i)
    return


def check_isfile_or_isfiles(path,pathnames=False,envvar=False,verbose=False):
    add=""
    if envvar == True: add = "(set as environment variable)"
    if type(path) is str:
        if not os.path.isfile(path):
            sys.exit('missing file '+path+add)
    elif type(path) is list:
        for idx,i in enumerate(path):
            if verbose:
                print('idx',idx,path[idx],'path',i)
            if i is None and type(pathnames) != bool:
                sys.exit(pathnames[idx]+" "+add+' is not defined!')
            if not os.path.isfile(i):
                sys.exit('missing file '+i)
    else:
        sys.exit('unknown type file type path '+str(path))
    return

def get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0,cubic=False,create_fake_vacancy=False,whichcell="fcc",normal_ordering=True):
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
    if whichcell == "fcc":
        atom = ase_build_bulk('Al',crystalstructure='fcc',a=a0,cubic=cubic)
    elif whichcell == "hcp":
        a = 3.21
        c = 5.21
        atom = crystal('Mg', [(1./3., 2./3., 3./4.)], spacegroup=194, cellpar=[a, a, c, 90, 90, 120])
    elif whichcell == "dc":
        a = 5.47
        atom = crystal('Si', [(0,0,0)], spacegroup=227, cellpar=[a, a, a, 90, 90, 90])
    else:
        sys.exti("whichcell has to be in fcc or hcp")

    atomsc = atom.repeat(ncell)
    number_of_atoms = atomsc.get_number_of_atoms()
    nal = number_of_atoms - nsi - nmg

    #for i in np.arange(nmg):
    #    atomsc[i].symbol = 'Mg'
    #for i in np.arange(nmg,nmg+nsi):
    #    atomsc[i].symbol = 'Si'

    # This is the order which ipi kmc expects
    if type(normal_ordering) == bool:
        for i in np.arange(nsi):
            atomsc[i].symbol = 'Si'
        for i in np.arange(nsi,nmg+nsi):
            atomsc[i].symbol = 'Mg'

        if create_fake_vacancy == False:
            for i in np.arange(nvac):
                del atomsc[-1]
        elif create_fake_vacancy == True:
            startsubst = -1
            for i in np.arange(nvac):
                #print('startsubst',startsubst)
                atomsc[startsubst].symbol = 'V'
                startsubst -= 1
        else:
            sys.exit("create_fake_vacancy has to be True or False")
    elif type(normal_ordering) == str:
        normal_ordering_element = normal_ordering.split("_")[0]
        normal_ordering_pos     = normal_ordering.split("_")[1]

        for i in np.arange(nvac):
            atomsc[i].symbol = 'V'
        #print('normal_ordering_element',normal_ordering_element)
        #symb = 'Si'
        #symb = 'Mg'
        symb = normal_ordering_element
        if normal_ordering_pos == '0':
            return atomsc


        if normal_ordering_pos == '1':
            atomsc[1].symbol = symb
        if normal_ordering_pos == '2':
            atomsc[5].symbol = symb
        if normal_ordering_pos == '3':
            atomsc[25].symbol = symb
        if normal_ordering_pos == '4':
            atomsc[1].symbol = symb
            atomsc[5].symbol = symb
        if normal_ordering_pos == '5':
            atomsc[1].symbol = symb
            atomsc[5].symbol = symb
            atomsc[25].symbol = symb
    else:
        sys.exit('random_ordering has to be of type bool or str!')
    #number_of_atoms = atomsc.get_number_of_atoms()
    #nal = number_of_atoms - nsi - nmg
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

def progress(count, total, status=''):
    ''' use:
    for i in range(DB2len):
        my.progress(i,DB2len)
    '''
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ... %s\r' % (bar, percents, '%', status))
    sys.stdout.flush()  # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)
    return

def diff(first, second):
    second = set(second)
    return [item for item in first if item not in second]

def ase_get_neighborlist(frame,atomnr=0,cutoff=3.,skin=0.1):
    NN = NeighborList([cutoff/2.]*frame.get_number_of_atoms(),skin=skin,self_interaction=False,bothways=True,primitive=NewPrimitiveNeighborList)
    #print('NN')
    NN.update(frame)
    #print(NN.get_connectivity_matrix())
    #for idx,i in enumerate(NN.get_connectivity_matrix()):
    #    print('idx',idx,i)
    NN_indices, offsets = NN.get_neighbors(atomnr)
    #print('NN idx:',np.sort(NN_indices))
    #sys.exit()
    return np.sort(NN_indices)

def print_args(args):
    keys = args.__dict__.keys()
    values = args.__dict__.values()
    print('########################################## argparse values (begin) #########')
    for idx,i in enumerate(keys):
        print("#",str(keys[idx]).ljust(25),str(values[idx]).ljust(30),str(type(values[idx])).ljust(15),"#")
    print('########################################## argparse values (end  ) #########')
    return

def count_amount_1NN_around_vacancies(filename,cutoffa=3.,cutoffb=4.5,skin=0.1,format='ipi',vac_symbol="V",save_every = 10,filename_save="KMC_analyze_akmc_ext"):
    print()
    print("########### count_amount_1NN_around_vacancies(...) #######################")
    print('reading',os.path.abspath(filename),'...')
    frames = ase_read(filename,index=":",format=format)
    print('reading',os.path.abspath(filename),'done.')

    structures = len(frames)
    structures = structures - 1 # just to make the lenght equal between
    print('structures',structures)

    all_vac_idx = ([atom.index for atom in frames[0] if atom.symbol == vac_symbol])
    print('all_vac_idx',all_vac_idx)

    if cutoffa == False:
        print(frames[0].cell)
        print(frames[0].get_number_of_atoms())
        #a0 = frames[0].cell[0,0]/5.*sqrt(2.)
        nndist = cutoffa = 2.95
    if cutoffb == False:
        a0 = cutoffb = 4.29

    if cutoffa == False:sys.exit('cutoffa')
    if cutoffb == False:sys.exit('cutoffb')
    print('cutoffa',cutoffa)
    print('cutoffb',cutoffb)

    filename_analyze_all = []
    al_mg_si_all = []
    for vac_nr,vac_idx in enumerate(all_vac_idx):
        #filename_analyze = filename +  ".1NN.al_mg_si_vac_"+str(vac_nr)+".dat"
        filename_analyze = filename_save+"_"+str(vac_nr)
        filename_analyze_all.append(filename_analyze)
        print('filename_analyze',filename_analyze)

        if os.path.isfile(filename_analyze):
            al_mg_si = np.loadtxt(filename_analyze)
            al_mg_si_all.append(al_mg_si)
        else:
            al_mg_si = np.zeros((structures,7))
            al_mg_si_all.append(al_mg_si)
    print()
    print('filename_analyze_all',filename_analyze_all)
    print()

    def test_anz_nn(al_mg_si_all,vac_nr,step,exit=False):
        #if al_mg_si_all[vac_nr][step][0] == 0:
        #    return True
        anz_1NN = np.sum(al_mg_si_all[vac_nr][step][1:4])
        anz_2NN = np.sum(al_mg_si_all[vac_nr][step][4:7])
        #print(al_mg_si[step][1:4])
        #print('sum',np.sum(al_mg_si[step][1:4]))
        testanz = True
        do_continue = True
        printed = False
        if testanz:
            if anz_1NN != 12:
                #print('step:',str(step).ljust(6),"not 12 1NN but "+str(anz_1NN),'exit',exit)
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_a',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'exit',exit,"not 12 1NN but "+str(anz_1NN),'-> REDO')
                printed = True
                do_continue = False
                if exit == True: sys.exit("ERROR")
            if anz_2NN != 6 and printed == False:
                #print('step:',str(step).ljust(6),"not  6 2NN but "+str(anz_2NN),'exit',exit)
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_b',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'exit',exit,"not  6 2NN but "+str(anz_2NN),'-> REDO')
                printed = True
                do_continue = False
                if exit == True: sys.exit("ERROR")
        if printed == False:
            print('step:',str(step).ljust(6),'vac_nr',vac_nr,'already known_c',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'continue',do_continue)
        return do_continue

    for step in np.arange(structures):
        all_vac_idx = ([atom.index for atom in frames[step] if atom.symbol == vac_symbol])
        #print('step',step,'all_vac_idx',all_vac_idx)
        for vac_nr,vac_idx in enumerate(all_vac_idx):
            if al_mg_si_all[vac_nr][step,0] != 0:
                do_continue = test_anz_nn(al_mg_si_all,vac_nr,step,exit=False)
                if do_continue == True:
                    continue

            NN_1_indices, NN_2_indices = ase_get_neighborlist_1NN_2NN(frames[step],atomnr=vac_idx,cutoffa=cutoffa,cutoffb=cutoffb,skin=skin)
            NN_1_sym = [atom.symbol for atom in frames[step] if atom.index in NN_1_indices]
            NN_2_sym = [atom.symbol for atom in frames[step] if atom.index in NN_2_indices]
            NN_1_al = NN_1_sym.count("Al")
            NN_1_mg = NN_1_sym.count("Mg")
            NN_1_si = NN_1_sym.count("Si")
            NN_2_al = NN_2_sym.count("Al")
            NN_2_mg = NN_2_sym.count("Mg")
            NN_2_si = NN_2_sym.count("Si")
            al_mg_si_all[vac_nr][step,1] = NN_1_al
            al_mg_si_all[vac_nr][step,2] = NN_1_mg
            al_mg_si_all[vac_nr][step,3] = NN_1_si
            al_mg_si_all[vac_nr][step,4] = NN_2_al
            al_mg_si_all[vac_nr][step,5] = NN_2_mg
            al_mg_si_all[vac_nr][step,6] = NN_2_si
            al_mg_si_all[vac_nr][step,0] = step
            anz_1NN = np.sum(al_mg_si_all[vac_nr][step][1:4])
            anz_2NN = np.sum(al_mg_si_all[vac_nr][step][4:7])
            #str_1NN = str(NN_1_al).ljust(3)+str(NN_1_mg).ljust(3)+str(NN_1_si).ljust(3)
            #str_2NN = str(NN_2_al).ljust(3)+str(NN_2_mg).ljust(3)+str(NN_2_si).ljust(3)
            #print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, "||",str(NN_1_al).ljust(3),NN_1_mg,NN_1_si,"||",NN_2_al,NN_2_mg,NN_2_si)
            #print('-------',filename_analyze_all[vac_nr])
            #print(al_mg_si_all[vac_nr])
            if anz_1NN != 12 or anz_2NN != 6:
                print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'ERROR!')
                sys.exit("ERROR see above")
            print('step:',str(step).ljust(6),'vac_nr',vac_nr,'NEW/REDO       ',anz_1NN,anz_2NN, al_mg_si_all[vac_nr][step],'OK')
            #do_continue = test_anz_nn(al_mg_si_all,vac_nr,step,exit=True)

            if step > 0 and step in np.arange(structures)[::save_every]:
                np.savetxt(filename_analyze_all[vac_nr],al_mg_si_all[vac_nr],fmt='%i')
                print('saving',os.path.abspath(filename_analyze_all[vac_nr]),'at step',step)

    # save everything in the very end
    np.savetxt(filename_analyze_all[vac_nr],al_mg_si_all[vac_nr],fmt='%i')
    print('saving (very end)',os.path.abspath(filename_analyze_all[vac_nr]),'at step',step)
    return


def ase_get_neighborlist_1NN_2NN(frame,atomnr=0,cutoffa=3.,cutoffb=4.5,skin=0.1):
    NN_1_indices       = ase_get_neighborlist(frame,atomnr=atomnr,cutoff=cutoffa,skin=skin)
    NN_1_2_indices_tmp = ase_get_neighborlist(frame,atomnr=atomnr,cutoff=cutoffb,skin=skin)
    NN_2_indices       = np.sort(np.array(diff(NN_1_2_indices_tmp,NN_1_indices)))
    return NN_1_indices, NN_2_indices

def ase_get_unique_frames(frames):
    '''
    this function only takes care of exactly same frames;
    structures which are close by another function will be necessary;
    '''
    framesout = deepcopy(frames)
    length = len(frames)
    #for idx,midx in enumerate(tdqm(range(len(frames))[::-1])):
    for idx,midx in enumerate(range(len(frames))[::-1]):
        progress(idx, length , status='')
        isin=False
        if frames[midx] in frames[:midx]:
            isin=True
            del framesout[midx]
        #print(idx,midx,frames[midx].positions[0,0],isin)
    print('returning framesout ... (unique)')
    return framesout


def ase_enepot(atoms,units='eV',verbose=False):
    ''' units: eV, eV_pa, hartree, hartree_pa
        check before if calculator is attached
    '''

    #print('now in ene')
    #print('ac',atoms.cell)
    try:
        # in the case of "DFT"                 , it just retrieves the energy  -> get_stress() can NOT be obtained.
        # in the case of of an ace calculations, it calculates the energy      -> get_stress() CAN     be obtained.
        #print('before')
        ene = atoms.get_potential_energy()
        #print('ene:',ene,atoms.info)
        #uuid = atoms.info["comment"]
        #print('uuid',uuid)
        #stress = atoms.get_stress()
        #print('stress:',stress)
    except: # RuntimeError:
        #print("had runtime error, e.g. ther cant be an energy in the POSCAR")
        ene = 0.
        #stress = False
    if verbose > 1:
        print('ene eV',ene,"(not per atom)")
    units_split = units.split("_")
    #print('us',units_split,units_split[1])
    if units_split[0].lower() == 'ev':
        pass
    elif units_split[0].lower() == 'mev':
        ene = ene*1000.
    elif units_split[0] == "hartree" or units_split[0] == "Hartree":
        ene = ene/aseunits.Hartree

    if len(units_split) == 2:
        if units_split[1] == 'pa':
            ene = ene/atoms.get_number_of_atoms()
        else:
            sys.exit("energy can not have this units (ending must be pa, eV_pa or hartree_pa)")

    return ene

def remove_duplicates_in_numpy_xy_array_and_sort(array,roundto=15):
    disp_vs_force = array
    data = disp_vs_force.copy()  # shifted to simulate a distance @ lattice constants @Tmelt
    sorted_idx = np.lexsort(data.T)
    sorted_data =  data[sorted_idx,:]
    row_mask = np.append([True],np.any(np.diff(sorted_data,axis=0),1))
    disp_vs_force = sorted_data[row_mask]
    disp_vs_force = disp_vs_force[disp_vs_force[:,0].argsort()]
    zero_wh = np.where(disp_vs_force[:,1] == 0)[0]
    #print('out')
    #print(disp_vs_force)
    #print('zw',zero_wh,'len',len(zero_wh))
    if len(zero_wh) > 1:
        #print('lzg',round(disp_vs_force[zero_wh[0]][0],15))
        #print('lzg',round(disp_vs_force[zero_wh[1]][0],15))
        if round(disp_vs_force[zero_wh[0]][0],roundto) == round(disp_vs_force[zero_wh[1]][0],roundto):
            #print("jodel")
            #del disp_vs_force[zero_wh[0]]
            disp_vs_force = np.delete(disp_vs_force,zero_wh[1],0)
    #print('out2')
    #print(disp_vs_force)
    #sys.exit()
    return disp_vs_force

def ase_get_chemical_symbols_to_number_of_species(atoms):
    symbols = atoms.get_chemical_symbols()
    #numat = atoms.get_number_of_atoms()
    #print('symbols',symbols)
    #print('numat',numat)

    uniquesym = set(atoms.get_chemical_symbols())
    d = {}
    for i in uniquesym:
        #print(i,symbols.count(i),numat)
        d[i] = symbols.count(i)

    def dcheck(element):
        if element in d.keys():
            #print(element+" exists")
            pass
        else:
            #print(element+" does not exist")
            d[element] = 0

    dcheck("Mg")
    dcheck("Si")
    dcheck("Al")
    return d

def ase_get_chemical_symbols_to_conz(atoms):
    symbols = atoms.get_chemical_symbols()
    numat = atoms.get_number_of_atoms()
    #print('symbols',symbols)
    #print('numat',numat)

    uniquesym = set(atoms.get_chemical_symbols())
    # print("uniquesym",uniquesym) # --> set(['Mg', 'Al'])
    d = {}
    for i in uniquesym:
        #print(i,symbols.count(i),numat)
        #print('ii',i,'---',float(symbols.count(i))/float(numat))
        d[i] = float(symbols.count(i))/float(numat)
    #print('kk',d)
    def dcheck(element):
        if element in d.keys():
            #print(element+" exists")
            pass
        else:
            #print(element+" does not exist")
            d[element] = 0.0
    #print('dd',d)
    dcheck("Mg")
    dcheck("Si")
    dcheck("Al")
    return d

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
    if not os.path.isdir(scripts):
        print('scripts:',scripts)
        sys.exit('scripts variable is not defined or is not an existing folder')
    return scripts

def test_and_return_environment_var_path(var,path=False,exit=True):
    variable = os.environ[var]
    if path == False:
        if not os.path.isfile(variable):
            message='The variable '+str(var)+' is not defined or is not an existing file'
            if exit == True:
                sys.exit(message)
            else:
                print(message)
    else:
        if not os.path.isdir(variable):
            message = 'directory '+str(var)+' is not defined or does not exist'
            if exit == True:
                sys.exit(message)
            else:
                print(message)
    return variable


def runner_exec(test=False):
    ''' return environment variable runner_exec (RuNNer executable)'''
    runner_exec = os.environ['runner_exec']
    if test == False and not os.path.isfile(runner_exec):
        sys.exit('runner_exec variable is not defined or is not an existing file')
    return runner_exec

def hostname():
    hostname = gethostname()
    return hostname

def get_click_defaults():
    # show default values in click
    orig_init = click.core.Option.__init__
    def new_init(self, *args, **kwargs):
        orig_init(self, *args, **kwargs)
        self.show_default = True
    click.core.Option.__init__ = new_init
    CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'],token_normalize_func=str.lower)
    return CONTEXT_SETTINGS


class mypot( object ):
    ''' return a list of available potentials '''
    def __init__(self,pot=False,potpath=False,verbose=False):
        self.pot        = pot         # n2p2_v1ag
        self.potpath_in = potpath
        self.potpath    = False
        self.pottype    = False       # n2p2/runner
        self.potDONE    = False         # n2p2_v1ag
        self.potlib     = False         # n2p2_v1ag
        self.potcutoff  = False         # n2p2_v1ag

        self.pot_all    = False   # n2p2_v1ag

        self.trigger_set_path_manual = ["setpath","potpath","pp", ".", "..", "../" ]
        self.verbose    = verbose

        self.elements    = False  # list e.g. ['Al', 'Mg', 'Si']
        self.atom_energy = False
        return


    def get_pot_all(self):
        if type(self.pot_all) == bool:
            scripts = os.environ['scripts']
            allpot_fullpath = glob.glob(scripts+'/potentials/*')
            pot_all_tmp = []
            #print('self.pot_all in :',self.pot_all)
            for i in allpot_fullpath:
                pot_all_tmp.append(os.path.basename(i))
            self.pot_all = pot_all_tmp + self.trigger_set_path_manual
            #print('self.pot_all out:',self.pot_all)
        return

    def get_potpath(self):
        scripts = os.environ['scripts']
        allpot_fullpath = glob.glob(scripts+'/potentials/*')
        self.pot_all = self.trigger_set_path_manual
        for i in allpot_fullpath:
            #print(i)
            if type(self.pot) != bool:
                if self.pot == os.path.basename(i):
                    self.potpath = i
                    break
        return

        #print(pot_all)
        #sys.exit()

    def get_elements_and_atomic_energies(self):
        # get from input.nn the atomic energies
        if type(self.pot) != bool and type(self.elements) == bool and type(self.atom_energy) == bool:
            self.pottype = self.pot.split("_")[0]
            if self.pottype == "runner" or self.pottype == "n2p2":
                inputnn = self.potpath+"/input.nn"

            if os.path.isfile(inputnn):
                self.elements, self.atom_energy = inputnn_get_atomic_symbols_and_atom_energy(inputnn)

        return

    def print_variables_mypot(self,text="",print_nontheless=False):
        if self.verbose > 1 or print_nontheless:
            #print("calss mypot      ",text)
            print(text,"self.pot         ",self.pot)
            print(text,"self.potpath     ",self.potpath)
            print(text,"self.potpath_in  ",self.potpath_in)
            print(text,"self.potlib      ",self.potlib)
            print(text,"self.potcutoff   ",self.potcutoff)
            print(text,"self.elements    ",self.elements)
            print(text,"self.atom_energy ",self.atom_energy)
            print(text,"self.pottype     ",self.pottype)
            print(text,"self.potDONE     ",self.potDONE)
            print(text,"self.verbose     ",self.verbose)
            print(text,"self.pot_all     ",self.pot_all)
            print()

    def get(self):
        self.print_variables_mypot('get potential: in')

        ##########################################
        # get potential from path
        ##########################################
        if self.potpath == False and self.potpath_in == False and self.pot in [ ".." , "../", "." ]:
            self.potpath_in = os.path.abspath(self.pot)

        if self.potpath_in != False and self.potpath == False:
            #print('ka')
            #if self.potpath_in == ".":
            #    self.potpath_in =
            if not os.path.isdir(self.potpath_in):
                sys.exit(self.potpath_in+" does not exist! (1)")
            checkfiles = [ "input.nn", "scaling.data", "weights.012.data", "weights.013.data", "weights.014.data" ]
            for i in checkfiles:
                if not os.path.isfile(self.potpath_in+"/"+i):
                    sys.exit(self.potpath_in+"/"+i+" does not exist! (2)")
                #else:
                #    print(self.potpath_in+"/"+i)
            #if self.verbose: print('aa',self.potpath_in)
            self.potpath = os.path.abspath(self.potpath_in)

            with open(self.potpath_in+"/input.nn") as fp:
                for i, line in enumerate(fp):
                    if "NNP" in line:
                        self.pottype = 'n2p2'
                    if "RuNNer" in line:
                        self.pottype = 'runner'
                    elif i > 3:
                        break
            self.pot = self.pottype+"_frompath"
        else:
            ##########################################
            # get potential from string
            ##########################################
            self.get_potpath()
        self.get_elements_and_atomic_energies()
        self.print_variables_mypot('get potential: out')


        ### check if self.pottype can be computed on this host!
        add = ' Your lammps version does not seem to work with '+self.pottype+"!"
        if self.pottype == "runner":
            self.potlib = os.environ["LAMMPSPATH"]+"/src/USER-RUNNER"
            self.potcutoff = 14.937658735
        if self.pottype == "n2p2":
            self.potlib = os.environ["LAMMPSPATH"]+"/src/USER-NNP"
            self.potcutoff = 11.0

        if self.pottype in [ "runner", "n2p2" ] and os.path.isdir(self.potlib) == False:
            sys.exit("ERROR: "+self.potlib+" not found!"+add)
        self.potDONE = True
        return

def pot_all():
    all = mypot()
    all.get_pot_all()
    return all.pot_all

def show_ase_atoms_content(atoms,showfirst=10,comment = ""):
    print()
    print("#####################################################")
    print("# "+comment+" BEGIN #######")
    print("#####################################################")
    print('## atoms')
    print(atoms)
    print('## atoms.get_number_of_atoms()',atoms.get_number_of_atoms())
    #print(atoms.positions)
    print(atoms.get_positions()[:showfirst])
    print('## elements get_chemical_symbols()')
    print(atoms.get_chemical_symbols()) #[:showfirst])
    print(list(set(atoms.get_chemical_symbols())))
    print('## atoms.cell')
    print(atoms.cell)
    print('## aa.get_cell_lengths_and_angles()')
    print(atoms.get_cell_lengths_and_angles())
    print('##atom.numbers',atoms.numbers)
    print("#####################################################")
    print("# "+comment+" END #######")
    print("#####################################################")
    print()
    return

def create_submitskript_ipi_kmc(filepath,nodes,ntasks,lmp_par=False,ipi_inst=False,ffsocket=False,submittime_hours=71,SBATCH=True,LOOPFOLDER=False):
    ''' time is in min
        this should be a class so that it is not necessary to shuffle
        variables back and forth.
    '''

    def check(variable,command_name_str,typehere):
        if type(variable) != typehere:
            print("ERROR In create_submitskript_ipi_kmc")
            print("PROBLEM variable name:",command_name_str,"=",variable,"type(variable)",type(variable),"but should be",str(typehere))
            sys.exit()

    IPI_COMMAND    = test_and_return_environment_var_path('IPI_COMMAND')
    LAMMPS_COMMAND = test_and_return_environment_var_path('LAMMPS_COMMAND')
    N2P2_PATH = test_and_return_environment_var_path('N2P2_PATH',path=True)

    check(IPI_COMMAND,"IPI_COMMAND",str)
    check(LAMMPS_COMMAND,"LAMMPS_COMMAND",str)
    check(N2P2_PATH,"N2P2_PATH",str)
    check(nodes,"nodes",int)
    check(ntasks,"ntasks",int)
    check(lmp_par,"lmp_par",int)
    check(ipi_inst,"ipi_inst",int)
    check(ffsocket,"ffsocket",str)

    text1 = [
    "#!/bin/bash",
    ""]

    text2 = [
    "#SBATCH --job-name=NNP-mpi",
    "#SBATCH --get-user-env",
    "#SBATCH --output=_scheduler-stdout.txt",
    "#SBATCH --error=_scheduler-stderr.txt",
    "#SBATCH --nodes="+str(nodes),
    "#SBATCH --ntasks "+str(ntasks),
    "#SBATCH --time=00-"+str(submittime_hours)+":00:00",
    "#SBATCH --constraint=E5v4",
    ""]

    text3 = [
    "set +e",
    "export LD_LIBRARY_PATH=",
    "#source $MODULESHOME/init/bash    # necessary for zsh or other init shells",
    "module load intel intel-mpi intel-mkl fftw python/2.7.14",
    "export LD_LIBRARY_PATH="+N2P2_PATH+"/lib:${LD_LIBRARY_PATH} # necessary for n2p2",
    "export OMP_NUM_THREADS="+str(lmp_par),  # THIS LETS THE JOBS BE KILLED!
    ""]

    text4 = [
    "touch time.out",
    'date +%s >> time.out',
    "",
    "# sets up the internet/unix socket for connections both for i-PI and on the lammps side",
    'seed=`grep seed input-runner.xml | awk \'{print $3}\'`',
    'hostname=`hostname`',
    'seed_hostname=$hostname\_$seed',
    'sed -i \'s/<ffsocket.*/<ffsocket name="lmpserial" mode="'+ffsocket+'">/\' input-runner.xml',
    #'sed -i \'s/address>.*<.addr/address>\'$(hostname)\'<\/addr/\' input-runner.xml',
    #'sed -i \'s/all ipi [^ ]*/all ipi \'$(hostname)\'/\' in.lmp',
    'sed -i \'s/address>.*<.addr/address> \'\"$seed_hostname\"\' <\/addr/\' input-runner.xml',
    'sed -i \'s/all ipi [^ ]*/all ipi \'\"$seed_hostname\"\'/\' in.lmp',
    '',
    'rm -f /tmp/ipi_*',
    '',
    #'# runs i-PI on the calculating node, so you can use UNIX sockets',
    #'ssh $SLURM_JOB_NODELIST -C \" cd $SLURM_SUBMIT_DIR; nohup python '+IPI_COMMAND+' input-runner.xml &> log.i-pi &\"',
    '',
    'python '+IPI_COMMAND+' input-runner.xml &> log.i-pi &',
    #'',
    'sleep 10',
    '',
    'for i in `seq '+str(ipi_inst)+'`',
    'do',
    #'      srun --hint=nomultithread --exclusive -n '+str(lmp_par)+' --mem=4G '+LAMMPS_COMMAND+' < in.lmp > log.lmp$i  &',
    #'      srun --hint=nomultithread --exclusive --mem=4G '+LAMMPS_COMMAND+' < in.lmp > log.lmp$i  &',
    #'      srun -n '+str(lmp_par)+' --mem=4G '+LAMMPS_COMMAND+' < in.lmp > log.lmp$i  &',
    '      # INTERACTIVE:',
    '      # OMP_NUM_THREADS='+str(lmp_par)+" "+LAMMPS_COMMAND+' < in.lmp > log.lmp$i  &',
    '      srun --exclusive --ntasks=1 --cpus-per-task='+str(lmp_par)+' --mem=4G '+LAMMPS_COMMAND+' < in.lmp > log.lmp$i  &',
    'done',
    '',
    'wait',
    'date +%s >> time.out',
    'cat time.out | xargs | awk \'{print $2-$1-10}\' > tmptime',
    'mv tmptime time.out',
    'exit 0',
    ]

    text4o = [
    'folder_=`ls -1d seed*`',
    'hier=`pwd`',
    'for folder in $folder_;do',
    '    echo folder $folder',
    '    cd $hier',
    '    cd $folder',
    '    ./osubmit-ipi-kmc.sh &',
    'done',
    'wait'
    ]

    if LOOPFOLDER == True:
        text4 = text4o

    if SBATCH == True:
        text = text1 + text2 + text3 + text4
    else:
        text = text1 + text3 + text4


    f = open(filepath,'w')
    for i in text:
        f.write(i+"\n")

    call(['chmod', '0755', filepath])

    return


class ase_calculate_ene( object ):
    '''
    - should be evaluated just once to initialize and NOT FORE EVERY CALCULATION
      (due to the fact that mypot would be called unnecessarily)
    ase_calculate_ene (ace) class which holds lammps commands to be executed
    if only pot is defined, static calculation.
    '''
    def __init__(self,
            pot,
            potpath,
            units=False,
            geopt=False,
            kmc=False,
            verbose=False,
            temp=False,
            elastic=False,
            ):

        #self.pot = pot
        self.potpath = potpath
        self.mypot = False
        self.units = units.lower()
        self.geopt = geopt          # so far only for ene object.
        self.elastic = elastic
        self.elastic_relax = True
        self.nsteps = 0
        self.verbose = verbose
        self.atoms = False          # ase atoms object (frame)
        self.pot = mypot(pot,self.potpath,verbose=self.verbose)

        #####################
        # for the calculator
        #####################
        self.calculator = "lammps"
        self.lmpcmd     = False         # in case we run through ase (also needs lmpcmd) or external lammps
        self.atom_types = False     # for nn pot
        self.keep_alive = True

        #print('init')
        # case of MD or KMC
        self.kmc = kmc
        self.temp = temp

        #self.eos = [ False, False, False, False] # e0_meV/pa, v0_ang^3/pa, B0, B0der]
        return


    def print_variables_ase(self,text=""):
        if self.verbose > 1:
            tt = 'ase_calculate_ene.'
            print()
            print(text,tt+'pot.pot      (1) :',self.pot.pot)    # : n2p2_v2ag
            print(text,tt+'units        (1) :',self.units)      # : ev
            print(text,tt+'geopt        (1) :',self.geopt)      # : False
            print(text,tt+'elastic      (1) :',self.elastic)    # : False
            print(text,tt+'nsteps       (1) :',self.nsteps)     # : 0
            print(text,tt+'lmpcmd       (1) :',self.lmpcmd)     # : False
            print(text,tt+'kmc          (1) :',self.kmc)        # : False
            print(text,tt+'temp         (1) :',self.temp)       # : False
            print(text,tt+'verbose      (1) :',self.verbose)       # : False
            print()
        return


    def lammps_command_potential_n2p2(self):
        #units_giulio_ene = "0.0367493254"
        ase_units_ene    = "0.03674932247495664" # 1./ase.units.Hartree

        #units_giulio_bohr = "1.8897261328"
        ase_units_bohr    = "1.8897261258369282" # 1./ase.units.Bohr
        command = [
        # showewsum 1 showew yes resetew no maxew 1000000
        'variable nnpDir string \"'+self.pot.potpath+'\"',
        "pair_style nnp dir ${nnpDir} showew no resetew yes maxew 100000000 cflength "+ase_units_bohr+" cfenergy "+ase_units_ene,
        "pair_coeff * * "+str(self.pot.potcutoff),
        #"#write_data ./pos.data # would this be the final struct?"
        ]
        return command

    def lammps_command_potential_runner(self):
        command = [
        # comment
        #"# thermo 1 # for geopt",
        'variable nnpDir string \"'+self.pot.potpath+'\"',
        "pair_style runner dir ${nnpDir} showewsum 1 showew yes resetew no maxew 1000000",
        #"# pair_coeff * * 7.937658735",
        #"pair_coeff * *  14.937658735"
        "pair_coeff * * "+str(self.pot.potcutoff)
        ]
        return command

    def lammps_command_masses(self):
        command = [
                "mass 1 24.305",
                "mass 2 26.9815385",
                "mass 3 28.0855",
                ]
        return command

    def pot_to_ase_lmp_cmd(self,kmc=False,temp=False,nsteps=0,ffsocket='inet',address=False):
        ''' geoopt (geometry optimization) is added / or not in
            lammps_write_inputfile(); here only the potential is set.
            ffsocket: ipi ffsocket [ "unix" or "inet" ]
        '''
        if self.pot.potDONE == False:
            self.pot.get()

        self.kmc = kmc
        self.temp = temp
        self.nsteps = nsteps
        self.ffsocket = ffsocket
        if self.ffsocket not in [ "unix", "inet" ]:
            print('ffsocket:',ffsocket)
            sys.exit('ffsocket has to be "unix" or "inet"; Exit!')



        if self.verbose > 2:
            tt = 'pot_to_ase_lmp_cmd_A '
            print(tt+'pot.pot    :',self.pot.pot)        # : n2p2_v2ag
            print(tt+'pot.potpath:',self.pot.potpath)   # :
            print(tt+'pot.pottype:',self.pot.pottype)   # :
            print()

        #sys.exit()
        # this depends only on the potential which is already defined
        # so should be easy to make this general.
        self.lmpcmd = [ "########## lmpcmd.begin #############" ]
        self.lmpcmd = self.lmpcmd + self.lammps_command_masses()

        if self.pot.pottype == "n2p2":
            # showewsum 1 showew yes resetew no maxew 1000000
            self.lmpcmd = self.lmpcmd + self.lammps_command_potential_n2p2()
            self.atom_types = {'Mg':1,'Al':2,'Si':3}

        elif self.pot.pottype == "runner":
            self.lmpcmd = self.lmpcmd + self.lammps_command_potential_runner()
            self.atom_types = {'Mg':1,'Al':2,'Si':3}
        else:
            sys.exit('pot '+str(self.pot.pot)+' not found! (X)')

        if self.kmc:
            if self.ffsocket == "unix": add = "unix"
            if self.ffsocket == "inet": add = ""
            if address == False:
                address = gethostname()
            self.lmpcmd = self.lmpcmd + [
                "",
                "timestep 0.001   # timestep (ps)",
                "velocity all create "+str(self.temp)+" 4928459",  # create initial velocities 4928459 is random seed for velocity initialization"
                "thermo 1   # screen output interval (timesteps)",
                "fix 1 all ipi "+str(address)+" 12345 "+str(add),
                ]
                # with n2p2 in the parallel version, inet is not working
                # "fix 1 all ipi fidis 12345",     # for fidis job
                # "fix 1 all ipi mac 77776 unix",  # for mac job

        self.lmpcmd = self.lmpcmd + [ "########## lmpcmd.end  #############" ]
        if self.verbose > 1:
            print('HERE THE lmpcmd I got',self.lmpcmd)
        self.print_variables_ase("pot_to_ase_lmp_cmd_FIN")
        return

    def define_wrapped_self_atoms(self,atoms=False):
        if atoms == False:
            atoms = self.atoms
        else:
            atoms = atoms
        atoms.wrap()
        #self.atomsin = deepcopy(atoms)
        self.atoms = atoms
        return self.atoms

    def ene_allfix(self,atoms=False):
        atoms = self.define_wrapped_self_atoms(atoms)
        keep_alive = False
        asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        atoms.set_calculator(asecalcLAMMPS)
        ene = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        return ene

    def ene_new(self,atoms=False,
            atomrelax=False,
            volumerelax=False,
            cellshaperelax=False,
            print_minimization_to_screen=False,minimizer="LGBFGS"):
        ''' atoms is an ase object
            if don_change_atomsobject is chosen,
        '''
        ## now the atoms object is not changed
        #atoms = atomsin.copy()
        atoms = self.define_wrapped_self_atoms(atoms)

        if atomrelax == False and self.geopt == False:
            pass
        if atomrelax == True and self.geopt == True:
            pass
        if atomrelax == True and self.geopt == False:
            pass  # realx wins, since locally set explicitely
        if atomrelax == False and self.geopt == True:
            atomrelax = True
        #print('svb',self.verbose)
        if self.verbose > 2:
            print()
            print("#####################################################")
            print('##--lmpcmd:',self.lmpcmd)
            print('##--atom_types',self.atom_types)
            print('##--geopt',self.geopt)
            print('##--atomrelax',atomrelax)
            print('##--minimizer',minimizer)
            print("#####################################################")
            show_ase_atoms_content(atoms,showfirst=10,comment = "START ASE INTERNAL CALUCLATION !!!")
            print()

        ################################
        ## set up the calculator
        ## (by attaching the calculator to atoms)
        ## this is valid for relaxations and static calcs
        ################################
        keep_alive = False
        if atomrelax == False: keep_alive = False
        if atomrelax == True:  keep_alive = True
        asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        atoms.set_calculator(asecalcLAMMPS)

        ### attach to atoms to relax the cell
        constraint = False
        if volumerelax == False and atomrelax == False and cellshaperelax == False:
            constraint = False
        elif volumerelax == True and atomrelax == False and cellshaperelax == False:
            print('before murn')
            vinet = self.get_murn(atoms,verbose=False,return_minimum_volume_frame=True)
            print('after murn')
            #constraint = ExpCellFilter(atoms, hydrostatic_strain=True) does not work
        elif atomrelax == False and cellshaperelax == True:
            constraint = StrainFilter(atoms)  # this relaxes the cell shape & the volume while keeping atomic positions fixed
        else:
            sys.exit('still to define')

        #if cellrelax == True and atomrelax == False:
        #    constraint = StrainFilter(atoms)  # this relaxes the cell shape & the volume while keeping atomic positions fixed
        #    ## in this case it does not work out
        #    sys.exit('This gives a segmentation fault (coredump) when cellrelax == True and atomrelax == True ... in this case do only one, then the other one')
        #elif cellrelax == True and atomrelax == True:
        #    ## when doint both this is recommended
        #    constraint = ExpCellFilter(atoms)  # Modify the supercell and the atom positions.

        ## atomrelax = False and cellrelax = False works
        ## atomrelax = True  and cellrelax = False works
        ## atomrelax = True  and cellrelax = True  works
        ## atomrelax = False and cellrelax = True  NOPE
        if atomrelax == True or cellshaperelax == True or volumerelax == True:
            #if atomrelax == False and cellrelax == True:
            #    print('NOOOOOOOOOOOOOOOWWWWWWWWWWWWWWWWW'*3)
            ################################################
            ### with geometry optimization
            ################################################
            # in case of a verbose run:
            #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, log_file='./xlolg.lammps.log',tmp_dir="./",keep_alive=True,atom_types=self.atom_types)

            # in case  of a non verbose run
            #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
            #atoms.set_calculator(asecalcLAMMPS)
            #from ase.io.trajectory import Trajectory
            #traj = Trajectory('ka', mode='w',atoms=atoms)
            #opt = BFGS(atoms,trajectory="ni.traj")
            #opt.run(steps=20)
            minimizer_choices = [ 'BFGS', 'LGBFGS', 'FIRE', 'GPMin', 'bh', 'mh' ]
            if minimizer not in minimizer_choices:
                print("your minimizer",minimizer)
                print("available:",minimizer_choices)
                sys.exit("choose one of the proper minimizer choices")
            if minimizer == 'mh' and atomrelax == True and cellrelax == True:
                sys.exit("use minimahopping (mh) with either with atomrelax or with cellrelax; or use both but with other minimizer;  but not with both at a time, rather run those sequentially o")

            if print_minimization_to_screen:
                print("AAA nat:",atoms.get_number_of_atoms())
                print("AAA pos:",atoms.get_positions()[:4])
                print("AAA for:",atoms.get_forces()[:4])
                print("AAA fmx:",abs(atoms.get_forces()).max())
                print("AAA vol:",atoms.get_volume())
                print("AAA vpa:",atoms.get_volume()/atoms.get_number_of_atoms())

            ## the syntax apparently vareis a bit depending
            ## if constraint or not
            use = atoms
            if type(constraint) != bool:
                use = constraint

            if print_minimization_to_screen:
                print("minimizer:",minimizer)
                logfile="-" # output to screen
            else:
                logfile="tmp" # output to file and not to screen


            #if print_minimization_to_screen:
            #print('BBB logfile',logfile)

            if minimizer == 'BFGS':
                opt1 = BFGS(use,logfile=logfile) #,trajectory="ni.traj")
            elif minimizer == 'LGBFGS':
                opt1 = LBFGS(use,logfile=logfile) #,trajectory="test.traj")
            elif minimizer == 'GPMin':
                opt1 = GPMin(use,logfile=logfile) #,trajectory="test.traj")
            elif minimizer == 'FIRE':
                opt1 = FIRE(use,logfile=logfile) #,trajectory="ni.traj")
            elif minimizer == 'bh':
                kB = 1.38064852e-23
                kB = 1.6021765e-19
                opt1 = BasinHopping(atoms=use, # the system to optimize
                  temperature=1*kB, # 'temperature' to overcome barriers
                  dr=0.5,      # maximal stepwidth
                  optimizer=LBFGS, # optimizer to find local minima
                  fmax=0.1,      # maximal force for the optimizer
                  logfile=logfile)
            elif minimizer == 'mh':
                if os.path.isfile("tmp"):
                    os.remove("tmp")
                opt1 = MinimaHopping(atoms=use,logfile=logfile)
                opt1(totalsteps=10)
            #print('startrun....')
            maxsteps = 200
            if minimizer == 'bh':
                maxsteps = 3

            ######################################################
            ## MINIMIZE
            ######################################################
            if minimizer not in ['mh','bh']: # in all cases but
                #opt1.run(steps=maxsteps,fmax=0.005)
                #opt1.run(steps=maxsteps,fmax=0.0001)
                #print('start')
                opt1.run(fmax=0.0001)
                #print('maxsteps                ',maxsteps,type(maxsteps))
                #print('opt1.get_number_of_steps',opt1.get_number_of_steps(),type(opt1.get_number_of_steps()))

                if maxsteps == opt1.get_number_of_steps():
                    print('DID NOT CONVErGE IN '+str(maxsteps)+' number of minimizer steps!')
                    #print('---- cell -----')
                    #print(atoms.get_cell())
                    #print('---- positions -----')
                    #print(atoms.get_positions())
                    return np.nan

        if print_minimization_to_screen:
            print('UUU atomrelax:',atomrelax)
            print('UUU cellrelax:',cellrelax)
            print("UUU nat:",atoms.get_number_of_atoms())
            print("UUU pos[:4]:",atoms.get_positions()[:4])
            print("UUU for[:4]:",atoms.get_forces()[:4])
            print("UUU fmx:",abs(atoms.get_forces()).max())
            print("UUU vol:",atoms.get_volume())
            print("UUU vpa:",atoms.get_volume()/atoms.get_number_of_atoms())

        if self.verbose > 1:
            print('ZZ done2')
            print('ZZ self.units',self.units)
        ######################################################
        # calculate the energy
        ######################################################
        #print('atxxx',atoms)
        ene = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        if print_minimization_to_screen:
            print('atoms')
            print(atoms.get_positions()[:3])
        if self.verbose > 1:
            print('ZZ ene:',ene,self.units)
        if not print_minimization_to_screen and os.path.isfile("tmp"):
            os.remove("tmp")
        #sys.exit()
        #ene = atoms.get_total_energy()
        #if self.verbose:
        #    print('ene',ene)
        #return ene,ene/atoms.get_number_of_atoms()*1000.
        if self.verbose > 2:
            show_ase_atoms_content(atoms,showfirst=10,comment="FINISHED ASE INTERNAL CALUCLATION")
            print()
            print()

        #print('forces out',atomrelax,cellrelax)
        #print(atoms.get_forces()[:3])
        return ene

    def ene(self,atoms=False,
            atomrelax=False,
            cellrelax=False,
            cellshaperelax=False,
            cellvolumerelax=False,
            print_minimization_to_screen=False,
            minimizer="LGBFGS",
            debug=False):
        ''' atoms is an ase object
            if don_change_atomsobject is chosen,
        '''
        if debug:
            print_minimization_to_screen=True
            print('777 degub is on for calculation of ene')
        #unique_elements = list(set(atoms.get_chemical_symbols()))
        ## now the atoms object is not changed
        #atoms = atomsin.copy()
        atoms = self.define_wrapped_self_atoms(atoms)

        if atomrelax == False and self.geopt == False:
            pass
        if atomrelax == True and self.geopt == True:
            pass
        if atomrelax == True and self.geopt == False:
            pass  # realx wins, since locally set explicitely
        if atomrelax == False and self.geopt == True:
            atomrelax = True
        #print('svb',self.verbose)
        if self.verbose > 2 or debug:
            print()
            print("#####################################################")
            print('##--lmpcmd:',self.lmpcmd)
            print('##--atom_types',self.atom_types)
            print('##--geopt',self.geopt)
            print('##--atomrelax',atomrelax)
            print('##--minimizer',minimizer)
            print("#####################################################")
            show_ase_atoms_content(atoms,showfirst=10,comment = "START ASE INTERNAL CALUCLATION !!!")
            print()

        ################################
        ## set up the calculator
        ## (by attaching the calculator to atoms)
        ## this is valid for relaxations and static calcs
        ################################
        if atomrelax == False: keep_alive = False
        if atomrelax == True:  keep_alive = True
        self.keep_alive = keep_alive
        if debug:
            print('ATTACHING CALCULATOR!!')
        self.get_calculator(atoms)

        ### attach to atoms to relax the cell
        if debug:
            print('attaching constraints')
        constraint = False
        if cellrelax == True: # and atomrelax == False:
            constraint = StrainFilter(atoms)  # this relaxes the cell shape & the volume while keeping atomic positions fixed
            ## in this case it does not work out
            #sys.exit('This gives a segmentation fault (coredump) when cellrelax == True and atomrelax == True ... in this case do only one, then the other one')
        elif cellrelax == True and atomrelax == True:
            ## when doint both this is recommended
            constraint = ExpCellFilter(atoms)  # Modify the supercell and the atom positions.

        ## atomrelax = False and cellrelax = False works
        ## atomrelax = True  and cellrelax = False works
        ## atomrelax = True  and cellrelax = True  works
        ## atomrelax = False and cellrelax = True  NOPE
        if atomrelax == True or cellrelax == True:
            if atomrelax == False and cellrelax == True:
                print('NOOOOOOOOOOOOOOOWWWWWWWWWWWWWWWWW'*3)
            ################################################
            ### with geometry optimization
            ################################################
            # in case of a verbose run:
            #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, log_file='./xlolg.lammps.log',tmp_dir="./",keep_alive=True,atom_types=self.atom_types)

            # in case  of a non verbose run
            #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
            #atoms.set_calculator(asecalcLAMMPS)
            #from ase.io.trajectory import Trajectory
            #traj = Trajectory('ka', mode='w',atoms=atoms)
            #opt = BFGS(atoms,trajectory="ni.traj")
            #opt.run(steps=20)
            minimizer_choices = [ 'BFGS', 'LGBFGS', 'FIRE', 'GPMin', 'bh', 'mh' ]
            if minimizer not in minimizer_choices:
                print("your minimizer",minimizer)
                print("available:",minimizer_choices)
                sys.exit("choose one of the proper minimizer choices")
            if minimizer == 'mh' and atomrelax == True and cellrelax == True:
                sys.exit("use minimahopping (mh) with either with atomrelax or with cellrelax; or use both but with other minimizer;  but not with both at a time, rather run those sequentially o")

            if print_minimization_to_screen:
                print("AAA nat:",atoms.get_number_of_atoms())
                print("AAA pos:",atoms.get_positions()[:4])
                print("AAA for:",atoms.get_forces()[:4])
                print("AAA fmx:",abs(atoms.get_forces()).max())
                print("AAA vol:",atoms.get_volume())
                print("AAA vpa:",atoms.get_volume()/atoms.get_number_of_atoms())

            ## the syntax apparently vareis a bit depending
            ## if constraint or not
            atoms_or_constraint = atoms
            if type(constraint) != bool:
                atoms_or_constraint = constraint

            if print_minimization_to_screen:
                print("minimizer:",minimizer)
                logfile="-" # output to screen
            else:
                logfile="tmp" # output to file and not to screen


            #if print_minimization_to_screen:
            #print('BBB logfile',logfile)

            if minimizer == 'BFGS':
                opt1 = BFGS(atoms_or_constraint,logfile=logfile) #,trajectory="ni.traj")
            elif minimizer == 'LGBFGS':
                opt1 = LBFGS(atoms_or_constraint,logfile=logfile) #,trajectory="test.traj")
            elif minimizer == 'GPMin':
                opt1 = GPMin(atoms_or_constraint,logfile=logfile) #,trajectory="test.traj")
            elif minimizer == 'FIRE':
                opt1 = FIRE(atoms_or_constraint,logfile=logfile) #,trajectory="ni.traj")
            elif minimizer == 'bh':
                kB = 1.38064852e-23
                kB = 1.6021765e-19
                opt1 = BasinHopping(atoms=atoms_or_constraint, # the system to optimize
                  temperature=1*kB, # 'temperature' to overcome barriers
                  dr=0.5,      # maximal stepwidth
                  optimizer=LBFGS, # optimizer to find local minima
                  fmax=0.1,      # maximal force for the optimizer
                  logfile=logfile)
            elif minimizer == 'mh':
                if os.path.isfile("tmp"):
                    os.remove("tmp")
                opt1 = MinimaHopping(atoms=atoms_or_constraint,logfile=logfile)
                opt1(totalsteps=10)
            #print('startrun....')
            maxsteps = 200
            if minimizer == 'bh':
                maxsteps = 3

            ######################################################
            ## MINIMIZE
            ######################################################
            if minimizer not in ['mh','bh']: # in all cases but
                #opt1.run(steps=maxsteps,fmax=0.005)
                #opt1.run(steps=maxsteps,fmax=0.0001)
                #print('start')
                opt1.run(fmax=0.0001)
                #print('maxsteps                ',maxsteps,type(maxsteps))
                #print('opt1.get_number_of_steps',opt1.get_number_of_steps(),type(opt1.get_number_of_steps()))

                if maxsteps == opt1.get_number_of_steps():
                    print('DID NOT CONVErGE IN '+str(maxsteps)+' number of minimizer steps!')
                    #print('---- cell -----')
                    #print(atoms.get_cell())
                    #print('---- positions -----')
                    #print(atoms.get_positions())
                    return np.nan

        if print_minimization_to_screen:
            print('UUA atomrelax:',atomrelax)
            print('UUA cellrelax:',cellrelax)
            print("UUA nat:",atoms.get_number_of_atoms())
            print("UUA pos[:4]:",atoms.get_positions()[:4])
            #print('for',atoms.get_forces.__module__)
            #print('for',atoms.get_forces.__globals__)
            print("UUA for:",atoms.get_forces()[:4])
            print("UUA fmx:",abs(atoms.get_forces()).max())
            print("UUA vol:",atoms.get_volume())
            print("UUA vpa:",atoms.get_volume()/atoms.get_number_of_atoms())

        if self.verbose > 1:
            print('ZZ done243')

        if self.verbose > 1:
            print('ZZ done236')
            print('ZZ self.units',self.units)
        ######################################################
        # calculate the energy
        ######################################################
        #print('atxxx',atoms)
        ene = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        #print('jo')
        if print_minimization_to_screen:
            print('atoms')
            print(atoms.get_positions()[:3])
        if self.verbose > 1:
            print('ZZ ene:',ene,self.units)
        if not print_minimization_to_screen and os.path.isfile("tmp"):
            os.remove("tmp")
        #sys.exit()
        #ene = atoms.get_total_energy()
        #if self.verbose:
        #    print('ene',ene)
        #return ene,ene/atoms.get_number_of_atoms()*1000.
        if self.verbose > 2:
            show_ase_atoms_content(atoms,showfirst=10,comment="FINISHED ASE INTERNAL CALUCLATION")
            print()
            print()

        #print('forces out',atomrelax,cellrelax)
        #print(atoms.get_forces()[:3])
        return ene

    def stress(self,atoms=False):
        atoms = self.define_wrapped_self_atoms(atoms)
        try:
            stress = atoms.get_stress()
        except:
            stress = False
        #print('stress',stress)
        return stress

    def get_v0(self,atomsin=False):
        ''' the function will never change the atomsobject '''
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        atoms_murn = atomsin.copy()
        atoms_murn.wrap()

        keep_alive = False
        atomrelax = False
        if atomrelax == False: keep_alive = False
        if atomrelax == True:  keep_alive = True
        asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        atoms_murn.set_calculator(asecalcLAMMPS)

        ### relax the atoms_murn first to the equilibrium
        self.ene(atoms_murn,cellrelax=True,atomrelax=True)
        return atoms_murn.get_volume()

    def get_v0_pa(self,atomsin=False):
        ''' the function will never change the atomsobject '''
        return self.get_v0(atomsin=atomsin)/atomsin.get_number_of_atoms()

    def get_elastic_external(self,atomsin=False,verbose=False,text=False,get_all_constants=False):
        ''' the function will never change the atomsobject '''
        print('######## get_elastic_external #############')
        #print('HOME             :',os.environ['HOME'])
        #print('LD_LIBRARY_PATH  :',os.environ['LD_LIBRARY_PATH'])
        #print('hier1')
        #print(os.environ['HOME'])
        #print('hier2')
        #print(os.environ['LD_LIBRARY_PATH'])
        #print('hier2.2')
        #print(os.environ['PYTHONPATH'])
        #from lammps import lammps
        #print('hier3')
        #lmp = lammps()
        #print('hier4')
        #sys.exit('hier')
        #print("LMP",LAMMPS_COMMAND)

        #print('1')
        #if self.elastic != True:
        #    return
        #print('2')
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        frame = atomsin.copy()
        frame.wrap()
        #print('3')

        self.get_calculator(frame)  # to be able to calculate stress
        #print('4x')
        #print('stress original frame:',frame.get_stress())
        #print('stress original frame:',frame.get_stress())
        #print('5x')
        if verbose:
            print('frame cell::',frame.get_cell())
            print('verbose   ::',verbose)
        if self.elastic_relax == True:
            self.ase_relax_cellshape_and_volume_only(frame,verbose=verbose)
        print('stress relaxed frame :',frame.get_stress())
        print('volume relaxed frame :',frame.get_volume())
        print('ene    relaxed frame :',frame.get_potential_energy())
        print('frame cell',frame.get_cell())
        #ase_write('pos.runner',frame,format='runner')
        print('-------------- lammps_ext_calc -----------')
        ene_pot_lmp = lammps_ext_calc(frame,self,get_elastic_constants=get_all_constants)
        #print('ene_pot_lmp...kk',ene_pot_lmp)
        #sys.exit('88')

        if verbose:
            if type(text) != bool:
                text = printred(text)
                #if get_all_constants == True:
                ene_pot_lmp = ene_pot_lmp.replace('Elastic Constant ', 'Elastic Constant '+text+" ")
                #else:
                #    ene_pot_lmp = ene_pot_lmp.replace('Elastic Constant ', 'Elastic Constant '+text+" ")
                #ene_pot_lmp = ene_pot_lmp.replace('C44all =', printred('C44all ='))
            print("relaxed? "+str(self.elastic_relax)+";",ene_pot_lmp)
        return

    def get_elastic(self,atomsin=False,verbose=False):
        ''' the function will never change the atomsobject '''

        if atomsin == False:
            sys.exit('need to define atoms in this case XX')

        atoms_h = atomsin.copy()
        atoms_h.wrap()

        keep_alive = False
        atomrelax = False
        if atomrelax == False: keep_alive = False
        if atomrelax == True:  keep_alive = True
        #asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        #atoms_h.set_calculator(asecalcLAMMPS)  # wird in parcalc.py gesetzt
        #/home/glensk/miniconda2/lib/python2.7/site-packages/parcalc/parcalc.py


        #### load the elastic stuff
        #from parcalc import ClusterVasp, ParCalculate
        #from elastic import get_pressure, BMEOS, get_strain
        #from elastic import get_elementary_deformations, scan_volumes
        #from elastic import get_BM_EOS, get_elastic_tensor

        print('############## from elastic ################')
        if False:
            print('ene   ',self.ene(atoms_h))
            print('stress1',atoms_h.get_stress())
        self.ase_relax_cellshape_and_volume_only(atoms_h,verbose=False)
        if False:
            print('stress2',atoms_h.get_stress())
            print('cell',atoms_h.get_cell())

        cell_ref = (atoms_h.copy()).get_cell()
        atoms_work = atoms_h.copy()
        self.ase_relax_cellshape_and_volume_only(atoms_work,verbose=False)

        for i in [0.98,0.99,1.00,1.01, 1.02, 1.03]:
            cell_work = cell_ref.copy()
            #print('cell_ref',atoms_h.get_cell())
            cell_work[0,0] = cell_ref[0,0]*i
            #print('cell_ref',cell_ref)
            atoms_work.set_cell(cell_work,scale_atoms=True)
            #print('cell_ref?',atoms_h.get_cell())
            if False:
                print('strain',i,'stress?',atoms_work.get_stress())
            #print('cell_ref',cell_ref)



        ################################################################
        # from elastic
        # http://wolf.ifj.edu.pl/elastic/lib-usage.html
        ################################################################
        try:
            from elastic.elastic import get_cart_deformed_cell, get_lattice_type, get_elementary_deformations
        except ImportError:
            return
        from elastic import get_pressure, BMEOS, get_strain
        from elastic import get_BM_EOS, get_elastic_tensor
        from parcalc import ParCalculate

        sym = get_lattice_type(atoms_h)
        print('sym',sym)
        # Create 10 deformation points on the a axis
        systems = []
        ss=[]
        for d in np.linspace(-0.2,0.2,11):
            # get_cart_deformed_cell:
            # The axis is specified as follows: 0,1,2 = x,y,z ;
            # sheers: 3,4,5 = yz, xz, xy.
            # d: The size of the deformation is in percent and degrees, respectively.
            struct = get_cart_deformed_cell(atoms_h, axis=0, size=d)
            if verbose:
                print()
            strc = struct.get_cell()
            if verbose:
                print('d      ',d)
                print('struct :',strc[0],strc[1],strc[2])
            strca = atoms_h.get_cell()
            if verbose:
                print('atoms_h:',strca[0],strca[1],strca[2])
            stress = struct.get_stress()
            strain = get_strain(struct, atoms_h)
            pressure = get_pressure(stress)
            if verbose:
                print("stress :",stress)
                print("strain :",strain)
                print('pressure:',pressure)
            ss.append([strain, stress])
            systems.append(struct)

        def myparcalc():
            systems_all = get_elementary_deformations(atoms_h, n=5, d=0.33)
            if type(systems_all) != type([]) :
                sysl=[systems_all]
                #print('11')
            else:
                sysl=systems_all
                #print('22')

            res = []
            for n,s in enumerate(sysl):
                s.get_potential_energy()
                res.append([n,s])
            return [r for ns,s in enumerate(sysl) for nr,r in res if nr==ns]
        res = myparcalc()
        Cij, Bij = get_elastic_tensor(atoms_h, systems=res)
        print("Cij (GPa):", Cij/aseunits.GPa)

        ss=np.array(ss)
        lo=min(ss[:,0,0])
        hi=max(ss[:,0,0])
        mi=(lo+hi)/2
        wi=(hi-lo)/2
        xa=np.linspace(mi-1.1*wi,mi+1.1*wi, 50)

        # Now fit the polynomials to the data to get elastic constants
        # C11 component
        f=np.polyfit(ss[:,0,0],ss[:,1,0],3)
        c11=f[-2]/aseunits.GPa
        #print('ffff')
        #print(f)
        #print()
        #print(f[-2])

        # C12 component
        f=np.polyfit(ss[:,0,0],ss[:,1,1],3)
        c12=f[-2]/aseunits.GPa

        #np.savetxt('c11.dat',np.transpose([ss[:,0,0],ss[:,1,0]]))
        print('C11 = %.3f GPa, C12 = %.3f GPa => K= %.3f GPa' % (
                    c11, c12, (c11+2*c12)/3))

        ################################################################
        # daniels manual way
        ################################################################
        from daniel_strainstuff import _gen_strainmatrix, _apply_strain
        from daniel_strainstuff import _2lammpslattice
        from daniel_lmprun import _find_compliance_viaenergy
        from daniel_lmprun import run_fcc, _print_compliance_components
        from lammps import lammps
        lmp = lammps()

        print('############## daniels way  ################')
        for d in np.linspace(-0.2,0.2,5):
            struct = get_cart_deformed_cell(atoms_h, axis=3, size=d)
            if verbose:
                print()
            strc = struct.get_cell()
            if verbose:
                print('d      ',d)
                print('struct :',strc[0],strc[1],strc[2])
            strca = atoms_h.get_cell()
            if verbose:
                print('atoms_h:',strca[0],strca[1],strca[2])
            stress = struct.get_stress()
            strain = get_strain(struct, atoms_h)
            pressure = get_pressure(stress)
            if verbose:
                print("stress :",stress)
                print("strain :",strain)

        def my_get_cart_deformed_cell(base_cryst, size=1,verbose=False,vol=False):
            from ase.atoms import Atoms
            cryst = Atoms(base_cryst)
            uc = base_cryst.get_cell()
            if vol != False:
                uc = base_cryst.get_cell()*vol
            s = size/100.0
            L = np.diag(np.ones(3))
            #L = L * 0.9997
            if verbose:
                print(L)
                print()
            if False: # not volume conserving
                L[1, 2] += s/2.
                L[2, 1] += s/2.
            if True: # volume conserving
                L[0, 1] += s/2.
                L[1, 0] += s/2.
                L[2, 2] += (s**2.)/(4.-s**2.)
            if verbose:
                print(L)
                print()
            uc = np.dot(uc, L)
            cryst.set_cell(uc, scale_atoms=True)
            return cryst


        print()
        print("########### now only one deformed cell ###########")
        print('atoms_h.get_cell()     :')
        print(atoms_h.get_cell())
        print('atoms_h.get_stress()   :')
        print(atoms_h.get_stress())
        volfact = 1.0111  # C44 36.86712322917455  vol (68.43265399591313) = 17.108163499  d = 0.557 ang^3
        volfact = 1.0051  # C44 39.770796753113444 vol (67.22160401031846) = 16.805401003  d = 0.254 ang^3
        volfact = 1.0011  # C44 41.41759291379788  vol (66.42222758795769) = 16.605556897  d = 0.055 ang^3
        volfact = 1.0001  # C44 41.793101000472575 vol (66.2233786205119)  = 16.555844655  d = 0.005 ang^3
        volfact = 1.0000  # C44 41.82985611457602  vol (66.20351557966632) = 16.550878895  d = 0.000 ang^3
        volfact = 0.995477787394 # C44 43.3418676  vol (65.30941200550778) = 16.550878895  d = 0.000 ang^3
        volfact = 0.             # C44 42.3672012  vol (65.90412920697601) = 16.476032301  d = 0.075 ang^3

        # DFT                                      vol (65.905655194)      = 16.476413798491  (first converged)
        # DFT                                      vol (65.904129207)      = 16.476032301744  (beset converge)
        V0DFT = 16.476413798491 #  first converged
        V0DFT = 16.476032301744 #  beset converge
        a0DFT = (4.*V0DFT)**(1./3.)
        volfact = ((V0DFT*4.)/66.20351557966632)
        print('volfact',volfact)

        #cryst.set_cell(np.diag(np.ones(3))*a0, scale_atoms=True)
        #cell = cryst.get_cell()
        #cryst = my_get_cart_deformed_cell(cryst, size=0.2)


        ### This is to be at volume of DFT! (JUST IF YOU WANT TO CHECK HWO LARGE THE
        ### ERROR WOULD BE IF WE USE THIS (DFT) VOLUME
        if False:
            if volfact >= 1.0:
                atoms_h.set_cell(atoms_h.get_cell()*volfact, scale_atoms=True)
            else:
                atoms_h.set_cell(np.diag(np.ones(3))*a0DFT, scale_atoms=True)

        e0 = atoms_h.get_potential_energy()
        V0 = atoms_h.get_volume()
        print('atoms_h.get_cell()     :')
        print(atoms_h.get_cell())
        print('atoms_h.get_potential():',e0)
        print('atoms_h.V0',V0)
        print()
        print('volum0',V0,'   s 0               ene0: e0',e0)
        print('----------------------------------------------------------------------------------------')
        #for d in np.linspace(-0.1,0.1,4):
        points = 10
        ene_vs_strain = np.zeros((points,2))
        ene_vs_strain_wo = np.zeros((points,2))
        for idx,d in enumerate(np.linspace(-10.,10.,points)):
            sd = 0.2
            sd = d
            s = sd/100.

            cryst = my_get_cart_deformed_cell(atoms_h, size=sd,vol=False)
            cell = cryst.get_cell()
            scheck =strain= cell[0,1]/cell[0,0]*2.
            print('straincheck',s,scheck)
            ene_vs_strain[idx,0] = s
            ene_vs_strain[idx,0] = scheck
            ene_vs_strain_wo[idx,0] = scheck
            stress = cryst.get_stress()
            enecryst_eV_cell = cryst.get_potential_energy()  # eV for whole cell
            ene_vs_strain[idx,1] = (enecryst_eV_cell-e0)*1000./cryst.get_number_of_atoms()
            ene_vs_strain_wo[idx,1] = (enecryst_eV_cell)*1000./cryst.get_number_of_atoms()
            vol = cryst.get_volume()

            if True:
                if False:
                    print()
                    print('cryst.get_cell()     :')
                    print(cryst.get_cell())
                    print('cryst.get_stress()   :')
                    print(cryst.get_stress())
                if False:
                    print('cryst.energy:',enecryst_eV_cell)
                    print('cryst.get_volume()')
            C44 = (enecryst_eV_cell - e0)/vol*(2./(strain**2.))
            ase_write("out_c_check_vol_cons_widerange_DFTV0.runner",cryst,format='runner',append=True)
            #print(d,'de',de2)
            atb = ANGSTROM_TO_BOHRRADIUS = 1./aseunits.Bohr
            print("volume",str(vol).ljust(20),'s',str(round(s,5)).ljust(10),'ene',enecryst_eV_cell,'c44:',str(round(C44/aseunits.GPa,2)).ljust(5),cell[0]*atb,cell[2]*atb) #stress)
            print("volume",str(vol).ljust(20),'s',str(round(s,5)).ljust(10),'ene',enecryst_eV_cell,'c44:',str(round(C44/aseunits.GPa,2)).ljust(5),cell[0],cell[2]) #stress)
        np.savetxt("elastic_ene.dat",np.array([C44]))
        np.savetxt("ene_vs_strain_NN.dat",ene_vs_strain)
        np.savetxt("ene_vs_strain_NN_wo.dat",ene_vs_strain_wo)
        sys.exit()
        print()
        cryst = my_get_cart_deformed_cell(atoms_h, size=0.2,vol=False)
        print(cryst.get_cell())
        print('stress                ',stress)
        print('stress/2              ',stress/2.)
        print('stress/aeunits.GPa2/2.',stress/aseunits.GPa/2.)
        print('stress/aeunits.GPa2   ',stress/aseunits.GPa)
        print('stress',1000*stress/aseunits.GPa/2.)
        print('st C44',1000*stress[3]/aseunits.GPa/2.)
        print('st C44',1000*stress[5]/aseunits.GPa/2.)
        print()
        print()
        print()
        print("########### get murn structures ###########")
        sys.exit()
        print('linsp',np.linspace(-0.03,0.03,9))
        for d in np.linspace(-0.03,0.03,9):
            cryst.set_cell(atoms_h.get_cell()*(1.+d), scale_atoms=True)
            #print('ah--',atoms_h.get_cell()/4.)
            #print('volh',atoms_h.get_volume()/4.)
            print('volc',cryst.get_volume()/4.,cryst.get_potential_energy())
            #ase_write("out_murn.runner",cryst,format='runner',append=True)

        V0 = 16.476413798491 #  first converged
        V0 = 16.476032301744 #  beset converge
        a0 = (4.*V0)**(1./3.)
        print("V0",V0)
        print("a0",a0)
        cryst.set_cell(np.diag(np.ones(3))*a0, scale_atoms=True)
        print(np.diag(np.ones(3)))
        print(np.diag(np.ones(3))*a0)
        cell = cryst.get_cell()
        print('-->cryst.cell:',cryst.get_cell())
        print()
        #for d in [0,0.01,-0.01]: #np.linspace(-0.01,0.01,3):
        for d in np.linspace(-0.01,0.01,3):
            print()
            cryst.set_cell(cell*(1.+d), scale_atoms=True)
            print('cell0',cryst.get_cell())
            print('cell1',cell)
            print('-->d',d,'    volc',cryst.get_volume()/4.,cryst.get_potential_energy())
            #ase_write("out_c_check.runner",cryst,format='runner',append=True)
        print()
        print("##########")
        cryst.set_cell(cell, scale_atoms=True)
        cryst = my_get_cart_deformed_cell(cryst, size=0.2)
        print('cell0 ',cryst.get_cell())
        print('cell0v',cryst.get_volume()/4.)
        ase_write("out_c_check_vol_cons.runner",cryst,format='runner',append=True)
        return


    def get_calculator(self,atoms):
        if self.verbose > 1:
            for i in self.lmpcmd:
                print("lmpcmds    :",i)
            print("atom_types :",self.atom_types)
            print("keep_alive :",self.keep_alive)
        if self.calculator == "lammps":
            asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=self.keep_alive)
            atoms.set_calculator(asecalcLAMMPS)
        return

    def ase_relax_atomic_positions_only(self,atoms,fmax=0.0001,verbose=False):
        ''' The strain filter is for optimizing the unit cell while keeping scaled positions fixed. '''
        self.keep_alive = True
        self.get_calculator(atoms)

        if verbose:
            print('relax atomic positions; stress:',atoms.get_stress(),"volume per atom:",ase_vpa(atoms))

        logfile="-" # output to screen
        logfile="tmp" # output to file and not to screen
        opt = LBFGS(atoms,logfile=logfile)
        opt.run(fmax=fmax)
        if os.path.isfile("tmp"):
            os.remove("tmp")
        if verbose:
            print('relax atomic positions; stress:',atoms.get_stress(),"volume per atom:",ase_vpa(atoms))
        return

    def ase_relax_cellshape_and_volume_only(self,atoms,verbose=False):
        ''' The strain filter is for optimizing the unit cell while keeping scaled positions fixed. '''
        self.keep_alive = True
        #print('1')
        self.get_calculator(atoms)
        #print('2')

        if verbose: self.check_frame('ase_relax_cellshape_and_volume_only in',frame=atoms,verbose=verbose)
        #print('3')

        sf = StrainFilter(atoms)
        logfile="-" # output to screen
        logfile="tmp" # output to file and not to screen
        opt = BFGS(sf,logfile=logfile)
        opt.run(0.005)
        if os.path.isfile("tmp"):
            os.remove("tmp")
        if verbose: self.check_frame('ase_relax_cellshape_and_volume_only out',frame=atoms)
        return

    def get_murn(self,atomsin=False,verbose=False,
            return_minimum_volume_frame=False,
            return_frame_with_volume_per_atom=False,
            atomrelax=False,
            write_energies=False,
            get_to_minvol_first=True):
        ''' the murn will never change the atomsobject
        return_frame_with_volume_per_atom : volume can be specified and he frame scaled
        '''
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        if return_minimum_volume_frame == True or type(return_frame_with_volume_per_atom) != bool:
            atoms_murn = atomsin
        else:
            atoms_murn = atomsin.copy()

        atoms_murn.wrap()

        # probably not necessary since this is set up when necessary
        self.keep_alive = True
        self.get_calculator(atoms_murn)

        ### relax the atoms_murn first to the equilibrium
        if verbose: self.check_frame('get_murn 1 in',frame=atoms_murn)
        self.ase_relax_cellshape_and_volume_only(atoms_murn,verbose=verbose)
        if verbose: self.check_frame('get_murn 2 atfer cellshape relax',frame=atoms_murn)

        if atomrelax:
            self.ase_relax_atomic_positions_only(atoms_murn,fmax=0.0001,verbose=False)
            if verbose: self.check_frame('get_murn 2 atfer atomrelax only',frame=atoms_murn)

        dvol_rel=[0.97,0.975,0.98,0.985,0.99,0.995,0.998,1.0,1.002,1.005,1.01,1.015,1.02,1.025,1.03]
        #dvol_rel = np.arange(0.97,1.03,0.001)
        #dvol_rel = np.arange(0.995,1.005,0.0003)
        vol_pa = np.zeros(len(dvol_rel))
        ene_pa = np.zeros(len(dvol_rel))


        cell_ref = atoms_murn.get_cell()
        nat = atoms_murn.get_number_of_atoms()

        atoms_murn_loop = atoms_murn.copy()

        for idx,i in enumerate(dvol_rel):
            if verbose > 2:
                print()
                print('000 idx:',idx,'i:',i)
            atoms_murn_loop.set_cell(cell_ref*i,scale_atoms=True)
            if verbose > 2:
                print('111 cell',atoms_murn_loop.get_cell())
            #print('111 cell',atoms_murn_loop.get_cell())
            #print('111 cell',atoms_murn_loop.get_cell()[0,0]/atoms_murn_loop.get_cell()[1,1])
            vol=atoms_murn_loop.get_volume()
            if verbose > 2:
                print('222 vol',atoms_murn_loop.get_volume())
            #ene = self.ene(atoms_murn_loop)                         # works
            #ene = self.ene_new(atoms_murn_loop)                         # works
            ene = self.ene_allfix(atoms_murn_loop)                       # works
            #print('ams3',atoms_murn_loop.get_stress(),ase_vpa(atoms_murn_loop))
            if verbose > 2:
                stress = atoms_murn_loop.get_stress()[:3]
                cell = atoms_murn_loop.get_cell()
                pos = atoms_murn_loop.get_positions()[1]/cell[0,0]

                if cell[0,1] == cell[0,2] == cell[1,0] == cell[1,2] == cell[2,0] == cell[2,1] == 0:
                    if cell[0,0] == cell[1,1] == cell[2,2]:
                        if round(stress[0],6) == round(stress[1],6) == round(stress[2],6):
                            print('murn cell++  ',round(cell[0,0],6),pos)
                        else:
                            print('murn cell--  ',cell[0,0],stress)
                    else:
                        print('murn cell---  ',cell[0,0],cell[1,1],cell[2,2])
            #ene = ase_enepot(atoms_murn_loop) #,units=ace.units)    # core dump
            #ene = atoms_murn_loop.get_potential_energy()            # core dump
            #ene = ase_enepot(atoms_murn_loop,units=self.units,verbose=self.verbose)  # core dump
            if verbose > 2:
                print('333 ene',ene) #,ene2)
            vol_pa[idx] = vol/nat
            ene_pa[idx] = ene/nat
            if verbose > 2:
                print('idx:',str(idx).ljust(3),'i:',str(i).ljust(10),'vol:',str(vol).ljust(10),'ene:',ene)
            if verbose:
                stress = atoms_murn_loop.get_stress()[:3]
                print('i',str(i).ljust(5),'vol/nat',str(round(vol/nat,7)).ljust(10),'ene/nat',str(ene/nat).ljust(19),stress)

        if verbose: self.check_frame('get_murn 3 atfer loop           ',frame=atoms_murn)
        if verbose > 1:
            stress = atoms_murn.get_stress()[:3]
            cell = atoms_murn.get_cell()
            pos = atoms_murn.get_positions()[1]/cell[0,0]
            print("---1-->>",pos)
        if write_energies:
            print('we1')
            if type(write_energies) == bool:
                print('we2')
                write_energies = "energies.dat"
            np.savetxt(write_energies,np.transpose([vol_pa,ene_pa]))
        if verbose > 1:
            print('loop done')
        vinet = eos()
        data=np.transpose([vol_pa,ene_pa])
        vinet.fit_to_energy_vs_volume_data(datax=vol_pa,datay=ene_pa)
        #self.eos = vinet.parameters
        if verbose > 1:
            print('pars',vinet.parameters)
        if verbose: self.check_frame('get_murn 4 before min vol ret   ',frame=atoms_murn)
        #if return_minimum_volume_frame == True or type(return_frame_with_volume_per_atom) != bool:
        #    print('adapt atoms in')
        #    atomsin = atoms_murn.copy()
        #        # old
        #        #volume_in  = ase_vpa(atomsin)
        #        #volume_out = vinet.parameters[1]
        #        #if type(return_frame_with_volume_per_atom) != bool:
        #        #    volume_out = return_frame_with_volume_per_atom
        #        #volume_scale = (volume_out/volume_in)**(1./3.)
        #        ##print('volume_in',volume_in)
        #        ##print('volume_out',volume_out)
        #        ##print('scale',volume_scale)
        #        #atomsin.set_cell(atomsin.get_cell()*volume_scale,scale_atoms=True)

        if verbose: self.check_frame('get_murn 5 after  min vol ret   ',frame=atoms_murn)



        if verbose > 2:
            stress = atoms_murn.get_stress()[:3]
            cell = atoms_murn.get_cell()
            pos = atoms_murn.get_positions()[1]/cell[0,0]
            print("---2-->>",pos)
        return vinet.parameters

    def get_fh(self,atomsin=False,disp=0.03,debug=False,try_readfile=False,atomrelax=True):
        ''' the function will never change the atomsobject
        atomrelax: in most cases it is desirable to relax the atomic positions
        in few cases we might be tempted to assess the free energy for particular positions (e.g. to check weather the DFT equilibrium position is the stable position of a NN)
        '''
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        atoms_h = atomsin.copy()
        atoms_h.wrap()

        T0shift_ev_atom = self.ene(atoms_h)/atoms_h.get_number_of_atoms()
        #print('--> T0shift_ev_atom',T0shift_ev_atom)
        #sys.exit('T0shift_ev_atom')

        #if return_units == "mev_pa":
        #    return_mult = 1.
        #elif return_units == "ev_cell":
        #    return_mult = atomsin.get_number_of_atoms()/1000.
        #else:
        #    raise NameError('return_units conversion not yet set for '+return_units)

        if try_readfile:
            if os.path.isfile(try_readfile+"_hessematrix"):
                if debug: print('tryread 1x .....',try_readfile+"_hessematrix")
                hessematrix = hesse.read_Hessematrix(try_readfile+"_hessematrix")
                if debug: print('hesse 2x ...',hessematrix)
                if debug: print('hesse 3x ...',atoms_h.get_chemical_symbols())
                hes = hesse.hesseclass(listin=atoms_h.get_chemical_symbols(),H=hessematrix,show_negative_eigenvalues = False, Tmax=1000, T0shift_ev_atom = T0shift_ev_atom)
                if hes.has_negative_eigenvalues == True:
                    print("Negative Eigenvalues",hes.freqs)

                if debug: print('done 1x ...')

                #free_ene      = hes.ene_atom
                #print('fe pa',free_ene[:3])
                #print('fe pa',free_ene[-3:])
                #free_ene_cell = hes.ene_cell
                #print()
                #print('fe cell',free_ene_cell[:3])
                #print('fe cell',free_ene_cell[-3:])
                #print()
                #print('fe cell ev',hes.ene_cell_only_ev)
                #print()
                #print('ka',hes.ene_cell_only_ev_T0shifted[:3])
                #print('ka',hes.ene_cell_only_ev_T0shifted[-3:])
                #sys.exit()
                #get = np.loadtxt(try_readfile)
                #return np.transpose([get[:,0],get[:,1]*return_mult])
                return hes


        self.keep_alive = False
        self.get_calculator(atoms_h)

        if debug:
            print("###########################################")
            #`print("forces harmonic 3:",atoms_h.get_forces()[:3])
            #print("###########################################")
            print('in atoms_h str',atoms_h.get_stress())
            maxforce = np.abs(atoms_h.get_forces()).max()
            print('in maxforce in',maxforce)

        #ene = self.ene(atoms_h,cellrelax=True,atomrelax=True,print_minimization_to_screen=debug)
        if atomrelax == True:
            self.ase_relax_atomic_positions_only(atoms_h,fmax=0.0001,verbose=False)
            print('harmonic cell stress',atoms_h.get_stress())
        maxforce = np.abs(atoms_h.get_forces()).max()
        if debug:
            print("###########################################")
            #print("forces harmonic 4:",atoms_h.get_forces()[:3])
            #print("###########################################")
            print('out atoms_h str',atoms_h.get_stress())
            print('out maxforce out',maxforce)

        if atomrelax == True and maxforce > 0.0001:
            print("forces harmonic 4:",atoms_h.get_forces()[:3])
            print('maxforce',maxforce)
            sys.exit('maxforce is too large')

        if maxforce > 0.0001:
            # from this it needs to be deduced that NEGATIVE EIGENVALUES
            # calculating the hesse makes a segmentation fault with ace lammps
            return False
        if debug:
            print("###########################################")
            print('atoms_h.get_cell() before',atoms_h.get_cell())

        nat = atoms_h.get_number_of_atoms()
        if nat < 20:
            atoms_h *= (2,2,2)
        if debug:
            nat = atoms_h.get_number_of_atoms()
            print("###########################################")
            print('!!!!!!! nat:',nat)
            print("forces harmonic 2:",atoms_h.get_forces()[:3])
            print("stress:",atoms_h.get_stress())
            print("forces max harmonic 2:",abs(atoms_h.get_forces()).max())
            print('atoms_h.get_cell() after mult',atoms_h.get_cell())
            print("###########################################")
        pos0 = atoms_h.get_positions()
        hessematrix=np.zeros((pos0.shape[0]*3,pos0.shape[0]*3))

        ### schleife ueber alle atome, 1..32

        for iidx,i in enumerate(pos0): # loop over all atoms
            progress(iidx,len(pos0),status=try_readfile)
            #print(iidx,"/",pos0.shape[0]) # loop over xyz 1..3
            for jidx,j in enumerate(i):
                pos1 = np.copy(pos0)
                pos1[iidx,jidx] = pos1[iidx,jidx]+disp

                atoms_h.set_positions(pos1)
                fah = atoms_h.get_forces()
                hessematrix[iidx*3+jidx] = fah.reshape((1,pos0.shape[0]*        3))/(-disp)

        #np.savetxt("HesseMatrix.dat",hessematrix/97.173617,fmt="%.13f")
        if debug:
            print("get free energy ...")
            print('get_chemical_symbols()',atoms_h.get_chemical_symbols())
            print()

        hes = hesse.hesseclass(listin=atoms_h.get_chemical_symbols(),H=hessematrix,show_negative_eigenvalues = True, Tmax=1000, T0shift_ev_atom = T0shift_ev_atom)
        print('eigenvalues',hes.freqs)
        #try:
        #    #free_ene = (hes.ene_atom[300]-hes.ene_atom[0])[1]
        #    free_ene      = hes.ene_atom
        #    free_ene_cell = hes.ene_cell
        #except IndexError:
        #    free_ene = "UNSTABLE"
        #if debug and type(free_ene) != str:
        #    hes.write_ene_atom()
        if try_readfile:
            if True: #hes.has_negative_eigenvalues == False:
                # write in any case, if it has negative eigenvalues so be it
                hes.write_hessematrix(try_readfile+"_hessematrix")
                if hes.has_negative_eigenvalues == False:
                    # only write for ground state
                    hes.write_ene_atom(try_readfile+"_per_atom")
                    hes.write_ene_cell(try_readfile+"_per_cell")
            #np.savetxt(try_readfile,free_ene)

        #print('k T)shift   ',T0shift_ev_atom)
        #print('hes.ene_cell',hes.ene_cell)
        #print('hes.ene_atom',hes.ene_atom)
        #get = np.loadtxt(try_readfile)
        #return np.transpose([get[:,0],get[:,1]*return_mult])
        #return free_ene*return_mult
        return hes


    def check_frame(self,text,frame=False,verbose=True,setupcalc=True):
        #print('a')
        if type(text) != str:
            raise TypeError("need a text!")
        #print('b')

        if setupcalc == True:
            #print('c',verbose)
            self.keep_alive = True
            self.get_calculator(frame)
        #print('d',verbose)

        if verbose:
            print('check_frame::')
            print('check_frame::',frame.get_stress())
        check_stress_max = round(abs(frame.get_stress()).max(),5)
        check_vpa = round(ase_vpa(frame),2)  # only 2 digits can be nicely fitted
        check_force_max = round(abs(frame.get_forces()).max(),5)
        if verbose:
            print(text.ljust(42),'fm,sm,curr vol ',[check_force_max,check_stress_max,check_vpa]) # forces max , stress max
        return ['fm,sm,curr vol ',check_force_max,check_stress_max,check_vpa]



    def submit_aiida(self,atomsin=False):
        ## o) move this from the ase calass to something separate
        ## a) create inputxxx.data
        ## b) follow kmc_submit_inputdata_to_aiida.sh
        return

def ase_vpa(atoms):
    return atoms.get_volume()/atoms.get_number_of_atoms()
def ase_epa(atoms):
    return atoms.get_potential_energy()/atoms.get_number_of_atoms()
def ase_mepa(atoms):
    return atoms.get_potential_energy()/atoms.get_number_of_atoms()*1000.
def ase_fmax(atoms):
    return abs(atoms.get_forces()).max()


def string_to_index_an_array(array,string):
    ''' array will be indexed according to string e.g.
        array["string"] where "string" = ":7"
        array[":7"] = [ 0,1,2,3,4,5,6,7]
    '''
    if string == ":":                  # ":"
        return array
    elif string[0] == ':':             # begins with : e.g. :7
        #print('t1',string[1:])
        try:
            if type(int(string[1:])) == int:
                return array[:int(string[1:])]
        except ValueError:
            pass
    elif string[-1] == ':':             # ends with : e.g. 7:
        #print('t2',string)
        #print('t3',string[:-1])
        try:
            if type(int(string[:-1])) == int:
                return array[int(string[:-1]):]
        except ValueError:
            pass
    else:
        try:                           # 7
            a = int(string)
            return [array[int(string)]]
        except ValueError:
            pass
        return False


def q():
    host = hostname()
    #print('host',host)
    if host != 'fidis':
        return [],[],[]
    out=check_output(['q'])
    debug=False
    out2=out.split('\n')
    id=[]
    stat=[]
    cores=[]
    runtime=[]
    path=[]
    for idx,i in enumerate(out2):
        #if debug:
        #    print('i.split:',i.split(" "))
        str_list = filter(None, i.split(" "))
        if debug:
            print('str_list:',str_list,'len:',len(str_list))
        #if idx > 0 and len(str_list) > 2:
        if len(str_list) > 2 and str_list[0] != "JOBID":
            if debug:
                print('xxx str_list:',str_list,'len:',len(str_list))
            id.append(int(str_list[0]))
            stat.append(str_list[1])
            cores.append(int(str_list[2]))
            runtime.append(str_list[3])
            path.append(str_list[4])

    return id,stat,path

def lammps_write_inputfile_from_command(folder,filename='in.lmp',command=False):
    ''' command is holding the lammps commands '''
    mkdir(folder)
    f = open(folder+'/'+filename,'w')
    for i in command:
        f.write(i+"\n")
    f.close()
    return

def lammps_write_inputfile(folder,filename='in.lmp',positions=False,ace=False):
    ''' ace is the ace object holding the lammps commands '''
    mkdir(folder)
    f = open(folder+'/'+filename,'w')

    f.write("clear\n")
    f.write("units metal\n")
    f.write("boundary p p p\n")
    f.write("atom_style atomic\n")

    # necessary for structure 344, (343 shourd be 476846.3656628 probably eV)
    # if this structure would be wired, it would make sense to check the box tilt large for all other structures.
    f.write("box tilt large\n")
    f.write("read_data \"" +str(positions)+"\"\n")
    f.write("\n")

    # masses
    # variable nnpDir string
    # pair_style nnp dir
    # pair_coeff
    # neighbor
    for i in ace.lmpcmd:
        f.write(i+"\n")

    f.write("\n")

    if ace.geopt:
        ###############################################################
        # works but has the types wrong, can be read in by ase,
        # atomtypes (however) need to be changed manually
        ###############################################################
        f.write("dump dump_all all custom 1 %s id type x y z vx vy vz fx fy fz\n")

        ###############################################################
        # works and writes instaed of he types the actual element
        # e.g. Mg which can not be read in by ase
        ###############################################################
        # works and writes instaed of he types the actual element
        #f.write("dump dump_all all custom 1 %s id element x y z vx vy vz fx fy fz\n")
        #f.write("dump_modify dump_all element Mg Al Si\n")


        ###############################################################
        # try with edoardos comment
        ###############################################################
        #stride = 1  # print all
        #stride = 1000  # print first and last
        #f.write("dump dump_all all xyz "+str(stride)+" traj.out\n")  # dumps with edoardos pach
        #f.write("dump_modify dump_all element Mg Al Si\n")


        #f.write("min_style cg\n")
        f.write("min_style fire\n")
        f.write("minimize 1.0e-9 1.0e-10 1000 1000\n")

    if ace.kmc:
        ace.nsteps = 66000000   # this needs to by any very high number
    f.write("run "+str(ace.nsteps)+"\n")
    return

def ipi_write_inputfile(folder=False,filename='input.xml',positions='init.xyz',ace=False):
    ipi_cmd =[
    "<simulation mode='static' verbosity='high'>",
    "  <output prefix='simulation'>",
    "    <properties stride='1' filename='out'>  [ step, potential ] </properties>",
    "    <trajectory filename='pos' stride='1'> positions </trajectory>",
    "  </output>",
    "  <total_steps> 1000 </total_steps>",
    "  <prng>",
    "    <seed> 32342 </seed>",
    "  </prng>",
    "  <ffsocket name='lammps' mode='unix' pbc='true'>",
    "    <address> geop </address>",
    "  </ffsocket>",
    "  <system>",
    "    <initialize nbeads='1'>",
    "      <file mode='xyz'> init.xyz </file>",
    "    </initialize>",
    "    <forces>",
    "      <force forcefield='lammps'> </force>",
    "    </forces>",
    "    <motion mode='minimize'>",
    "      <optimizer mode='lbfgs'>",
    "        <tolerances>",
    "          <energy> 1e-6 </energy>",
    "          <force> 1e-6 </force>",
    "          <position> 1e-6 </position>",
    "        </tolerances>",
    "      </optimizer>",
    "    </motion>",
    "  </system>",
    "</simulation>",
    ]
    mkdir(folder)
    f = open(folder+'/'+filename,'w')
    for i in ipi_cmd:
        f.write(i+"\n")
    return

def ipi_ext_calc(atoms,ace):
    ''' atoms is an ase atoms object which can hold several frames, or just one'''
    ### mkdir tmpdir
    tmpdir = os.environ['HOME']+"/._tmp_ipi/"
    mkdir(tmpdir)
    ipi_write_inputfile(folder=tmpdir,filename='input.xml',positions='init.xyz',ace=ace)
    ase_write(folder+'/init.xyz',atoms,format='ipi')
    return



def lammps_ext_elastic_init_mod(ace,positions='pos.lmp'):
    command = [
    "# NOTE: This script can be modified for different atomic structures,",
    "# units, etc. See in.elastic for more info.",
    "#",
    "",
    "# Define the finite deformation size. Try several values of this",
    "# variable to verify that results do not depend on it.",
    "variable up equal 1.0e-8",
    "",
    "# Define the amount of random jiggle for atoms",
    "# This prevents atoms from staying on saddle points",
    "variable atomjiggle equal 1.0e-5",
    "",
    "# Uncomment one of these blocks, depending on what units",
    "# you are using in LAMMPS and for output",
    "",
    "# metal units, elastic constants in eV/A^3",
    "#units    metal",
    "#variable cfac equal 6.2414e-7",
    "#variable cunits string eV/A^3",
    "",
    "# metal units, elastic constants in GPa",
    "units    metal",
    "variable cfac equal 1.0e-4",
    "variable cunits string GPa",
    "",
    "# real units, elastic constants in GPa",
    "#units    real",
    "#variable cfac equal 1.01325e-4",
    "#variable cunits string GPa",
    "",
    "# Define minimization parameters",
    "variable etol equal 0.0",
    "variable ftol equal 1.0e-10" ]
    if ace.elastic_relax == True:
        command = command + [
    "variable maxiter equal 100 # do the relaxation of atomic positions" ]
    else:
        command = command + [
    "variable maxiter equal 0.0 # dont do relaxation of atomic positions"]

    command = command + [
    "variable maxeval equal 1000",
    "variable dmax equal 1.0e-2",
    "",
    "boundary  p p p",
    "box tilt large",
    "",
    "# generate the box and atom positions using a diamond lattice",
    'read_data "'+positions+'"',
    "",
    "# Need to set mass to something, just to satisfy LAMMPS",
    ] + ace.lammps_command_masses()
    return command


def lammps_ext_elastic_potential_mod(ace):
    command = [
    "# NOTE: This script can be modified for different pair styles",
    "# See in.elastic for more info.",
    "",
    "# Choose potential"]

    # change the potential
    if ace.pot.pottype == "n2p2":
        command = command + ace.lammps_command_potential_n2p2()
    elif ace.pot.pottype == "runner":
        command = command + ace.lammps_command_potential_runner()
    else: sys.exit('potential now known yet')

    command = command + [
    "",
    "# Setup neighbor style",
    "neighbor 1.0 nsq",
    "neigh_modify once no every 1 delay 0 check yes",
    "",
    "# Setup minimization style",
    "min_style       cg",
    "min_modify       dmax ${dmax} line quadratic",
    "",
    "# Setup output",
    "thermo    1",
    "thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol",
    "thermo_modify norm no"
    ]
    return command

def lammps_ext_calc(atoms,ace,get_elastic_constants=False):
    ''' atoms is an ase atoms object which can hold several frames, or just one'''

    ###############################################################
    # find LAMMPS executable
    ###############################################################
    #if ace.pot.pottype == "runner":
    #    LAMMPS_COMMAND = my.scripts()+'/executables/lmp_fidis_par_runner'
    #    #print("LAMMPS_COMMAND 1",LAMMPS_COMMAND)
    #else:
    #    #print("LAMMPS_COMMAND 2",LAMMPS_COMMAND)
    if ace.verbose:
        print("get_elastic_constants",get_elastic_constants)

    LAMMPS_COMMAND = os.environ['LAMMPS_COMMAND']
    if ace.verbose:
        print("LAMMPS_COMMAND",LAMMPS_COMMAND)
    #sys.exit()

    ###############################################################
    # make tmpdir (/home/glensk/._tmp_lammps/)
    ###############################################################
    tmpdir = os.environ['HOME']+"/._tmp_lammps/"
    mkdir(tmpdir)
    folder = tmpdir

    # delete all previous file in folder
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            #elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

    print('doing the calculation in lammps, externally, in folder',tmpdir)

    ###############################################################
    # write input structure (pos.lmp)
    ###############################################################
    if ace.verbose > 2:
        show_ase_atoms_content(atoms,showfirst=10,comment="START LAMMPS EXTERNALLY")
    atoms.set_calculator(None)
    atoms.write(tmpdir+'pos.lmp',format='lammps-runner')
    if ace.verbose:
        print('written ',tmpdir+'/pos.lmp')

    ###############################################################
    # write in.lmp (for mormal energy calculation)
    ###############################################################
    if not get_elastic_constants:
        execute_file = 'in.lmp'
        lammps_write_inputfile(folder=tmpdir,filename=execute_file,positions='pos.lmp',ace=ace)
        if ace.verbose:
            print('written ',tmpdir+'/'+execute_file)
    if ace.verbose > 1:
        print("written lammsp inputfile to ",tmpdir)

    ### calculate with lammps (trigger externally)
    ene = False
    if os.path.isfile(tmpdir+'log.lammps'):
        os.remove(tmpdir+"log.lammps")


    ###############################################################
    # IF get elastic constants
    ###############################################################
    if get_elastic_constants:
        execute_file = 'in.elastic'
        if ace.verbose:
            print('written ',tmpdir+'/potential.mod')
        lammps_write_inputfile_from_command(folder=tmpdir,filename='potential.mod',command=lammps_ext_elastic_potential_mod(ace))
        my.cp(scripts()+'/lammps_scripts/elastic/in.elastic',tmpdir+'/in.elastic')
        my.cp(scripts()+'/lammps_scripts/elastic/displace.mod',tmpdir+'/displace.mod')
        lammps_write_inputfile_from_command(folder=tmpdir,filename='init.mod',command=lammps_ext_elastic_init_mod(ace,positions='pos.lmp'))
        if ace.verbose:
            print('written ',tmpdir+'/'+execute_file)

    ###############################################################
    # cd to folder and run lammps
    ###############################################################
    with cd(tmpdir):  # this cd's savely into folder
        # RUN LAMMPS
        # without SHELL no LD LIBRARY PATH
        #print('pwd',os.getcwd())
        if ace.verbose:
            print(LAMMPS_COMMAND+" < "+execute_file+" > /dev/null")
        call([LAMMPS_COMMAND+" < "+execute_file+" > /dev/null"],shell=True)

        ###############################################################
        # extract energy and forces
        ###############################################################
        if not get_elastic_constants:
            ene = check_output(["tail -300 log.lammps | grep -A 1 \"Step Temp E_pai\" | tail -1 | awk '{print $3}'"],shell=True).strip()
            ene=float(ene)
            #print('ene',ene,'lammps in eV')
            #print('ace units',ace.units)
            print('ace.units.lower()',ace.units.lower())
            if ace.units.lower() == 'ev':
                pass
            elif ace.units.lower() == 'hartree':
                #ene = ene*0.036749325
                ene = ene*(1./aseunits.Hartree)
            elif ace.units.lower() == 'hartree_pa':
                ene = ene*(1./aseunits.Hartree)/atoms.get_number_of_atoms()
            elif ace.units.lower() == 'mev_pa':
                ene = ene*1000./atoms.get_number_of_atoms()
            else:
                sys.exit('units '+ace.units+' unknown! Exit!')
            #print('ene out',ene,ace.units)
        else:
            ace.elastic_constants = elastic_constants = check_output(["tail -300 log.lammps | grep \"^Elastic Constant\""],shell=True).strip()
            #ace.elastic_constants_ = elastic_constants = check_output(["tail -300 log.lammps | grep \"^Elastic Constant\""],shell=True)
            co = check_output(["tail -300 log.lammps | grep \"^Elastic Constant C44\" | sed 's|.*= ||' | sed 's|GPa.*||'"],shell=True)
            #print('co')
            #print(co)

            co1 = co.strip()
            #print('co1')
            #print(co1)
            co11 = str(co1.decode("utf-8"))
            #print('co11')
            #print(co11)

            #co2 = co1.split(" ")
            #print('co2')
            #print(co2)
            #ace.c44 = co2[4]
            ace.c44 = co11
            #print('aa c44',float(ace.c44))
            #print('aa el',ace.elastic_constants_)
            #for ij in ace.elastic_constants_:
            #    print('ij',ij)
            #sys.exit('ace44')

            #if get_elastic_constants == True:
            #    elastic_constants = check_output(["tail -300 log.lammps | grep \"^Elastic Constant\""],shell=True).strip()
            #else:
            #    elastic_constants = check_output(["tail -300 log.lammps | grep \"^Elastic Constant "+get_elastic_constants+"\""],shell=True).strip()
            #print('ec',elastic_constants.split(" "))
            #sys.exit('ec')
            ene = ace.elastic_constants

    if ace.verbose > 2:
        show_ase_atoms_content(atoms,showfirst=10,comment="FINISHED LAMMPS EXTERNALLY")
    return ene


def ase_repeat_structure(atoms,repeat):
    atomsc = atoms.repeat(repeat)
    cell_o = atoms.get_cell()
    cell_n = cell_o * repeat


    pos_o = atoms.get_positions()
    for_o = atoms.get_forces()
    nat_o = atoms.get_number_of_atoms()
    #print()
    #print(pos_o)
    #print()
    #print('sh old',pos_o.shape)
    #print('nat_o',nat_o)
    nat_n = nat_o*repeat**3
    #print('nat_n',nat_n)
    pos_n = np.zeros((nat_n,3))
    for_n = np.zeros((nat_n,3))
    #print('sh new',pos_n.shape)
    ## create translation matrix
    N1 = N2 = N3 = repeat
    var1 = np.mgrid[ 0:N1, 0:N2, 0:N3 ]
    ntrans1 = var1[0,:,:].flatten(1)
    ntrans2 = var1[1,:,:].flatten(1)
    ntrans3 = var1[2,:,:].flatten(1)
    transmat = np.zeros( [ len(ntrans1),3] )
    transmat[:,0] = ntrans1
    transmat[:,1] = ntrans2
    transmat[:,2] = ntrans3
    #print('tr',transmat)
    sclattpos = np.dot( transmat, cell_o )

    for ij, sclatt in enumerate(sclattpos):
        #print('-> ij',ij,sclatt)
        for idx,i in enumerate(pos_o):
            #print('ij',ij,'i',i,type(ij),type(i))
            #print('-> pos_o',i,'--->',sclatt+i)
            pos_n[ij*nat_o+idx] = sclatt + pos_o[idx]
            for_n[ij*nat_o+idx] = for_o[idx]
            #print('i',i,pos_n[ij+idx])

    #print('pos c ??')
    atomsc.set_positions(pos_n)
    #print(atomsc.get_positions()[:5])
    #sys.exit()
    return atomsc,for_n


def convert_energy(ene,units_in_,units_out_,frame,verbose=False):
    known = [ "hartree_pa", "ev_pa", "mev_pa", 'ev', "hartree" ]
    units_in  = units_in_.lower()
    units_out = units_out_.lower()
    if verbose:
        print('units_in :',units_in)
        print('units_out:',units_out)
    if units_in not in known or units_out not in known:
        print('units_in :',units_in)
        print('units_out:',units_out)
        print('known    :',known)
        sys.exit("units_in or units_out not in known")

    units_in_split  = units_in.split("_")
    units_out_split = units_out.split("_")
    #print('us in ',units_in_split,units_in_split[0])
    #print('us out',units_out_split,units_out_split[0])
    #print('ene in',ene)

    if units_in_split[0] == units_out_split[0]:
        pass
    elif units_in_split[0] == 'ev' and units_out_split[0] == 'hartree':
        if verbose:
            print('111',aseunits.Hartree)
        ene = ene / aseunits.Hartree
    elif units_in_split[0] == 'hartree' and units_out_split[0] == 'ev':
        ene = ene * aseunits.Hartree
        if verbose:
            print('2222',aseunits.Hartree)
    else: sys.exit('conversion from to now yet implemented')
    if verbose:
        print('ene (1)',ene)

    if len(units_in_split) == 2 and len(units_out_split) == 1:
        ene = ene * frame.get_number_of_atoms()
    elif len(units_in_split) == 1 and len(units_out_split) == 2:
        ene = ene / frame.get_number_of_atoms()
    if verbose:
        print('ene (2)',ene)
    return ene


def get_latest_n2p2_pot():
    checkfor = scripts()+"/potentials/n2p2_v*ag/"
    found = glob.glob(checkfor)
    ver = []
    for i in found:
        try:
            name2 = i.split("/potentials/n2p2_v")[-1]
            name3 = name2.split("ag")[0]
            name4 = int(name3)
            ver.append(name4)
            #print(i,name2,name3,name4,"ver",ver)
        except:
            pass
    ver = np.sort(np.array(ver))
    potout = "n2p2_v"+str(ver[-1])+"ag"
    #print('ver',ver,"->",potout)
    return potout


class ase_get_known_formats_class():
    """
    helptext
    """
    def __init__(self,verbose=False):
        self.formatspy              = os.path.dirname(ase.io.__file__)+"/formats.py"
        self.all_known_formats      = []
        self.all_known_formats_ase  = False
        self.my_formats_shall       = [ 'runner',   'lammps-runner', 'lammps-data'   ,'ipi'   , 'quippy'    ]
        self.my_formats_filenames   = [ "runner.py","lammpsrunner.py","lammpsdata.py","ipi.py", "quippy.py" ]
        self.my_formats_is          = []
        self.verbose                = verbose
        self.needcopy_since_missing = False
        return

    def get_all_known_formats(self):
        if len(self.all_known_formats) == 0:
            self.all_known_formats_ase = ase.io.formats.all_formats
            for i in self.all_known_formats_ase:
                self.all_known_formats.append(i)
        return

    def check_if_format_in_know_formats(self,typ):
        if typ in self.all_known_formats:
            if self.verbose: print(">> formats.py knows", typ)
            return True
        else:
            self.verbose = True
            if self.verbose: print(">> ERROR, formats.py does not know",typ)
            return False

    def copy(self):
        scripts = my.scripts()
        from_ = scripts+"/runner_scripts/ase_fileformat_for_"
        to = os.path.dirname(ase.io.__file__)+"/"
        for ff in self.my_formats_filenames:
            #print('copying ',from_+ff,'to',to+ff)
            if self.verbose:
                print('copying ',from_+ff)
                print('                    to',to+ff)
            shutil.copyfile(from_+ff,to+ff)

    def adapt_formatspy(self,writeformatspy = False):
        # check if formatspy exist
        if not os.path.isfile(self.formatspy):
            print('formatspy',self.formatspy)
            sys.exit('did not find '+str(self.formatspy))

        f = open(self.formatspy, "r")
        contents = f.readlines()
        f.close()
        insert=0
        insert2=0
        for idx,i in enumerate(contents):
            #print('i',idx,i)
            #print("|"+i[:20]+"|")
            if i[:20] == "    'abinit': ('ABIN":
                insert = idx
            if i[:30] == "    'lammps-data': 'lammpsdata":
                insert2 = idx

        for fo in [ 'runner', 'ipi', 'quippy' ]:
            if fo in self.all_known_formats_ase:
                if self.verbose:
                    print(fo.ljust(14)+'format are already added in formats.py (of ase).')
            else:
                print(fo.ljust(14)+'format NOT KNOWN in formats.py (of ase, WILL BE ADDED).')
                contents.insert(insert, "    '"+fo+"': ('"+fo+" input file', '+F'),\n")
                writeformatspy = True


        if 'lammps-runner' in self.all_known_formats_ase:
            if self.verbose:
                print('lammps-runner format are already added in formats.py (of ase).')
        else:
            contents.insert(insert, "    'lammps-runner': ('LAMMPS data input file for n2p2 or runner', '1F'),\n")
            contents.insert(insert2,"    'lammps-runner': 'lammpsrunner',\n")
            writeformatspy = True

        if writeformatspy == True:
            print('! adapting formatspy ... ! (seems to be save in any case)')
            #print('insert',insert)

            f = open(self.formatspy, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()
        else:
            if self.verbose:
                print('everything was already in',self.formatspy)
        return

    def check_if_default_formats_known(self,copy_and_adapt_formatspy_anyhow=False):
        self.get_all_known_formats()
        if self.verbose: print(">> formats.py", self.formatspy)
        for i in self.my_formats_shall:
            if self.check_if_format_in_know_formats(i) == False: self.needcopy_since_missing = True

        if self.needcopy_since_missing == True or copy_and_adapt_formatspy_anyhow == True: # only in this case copy stuff and add to formats.py
            self.verbose = True
            self.copy()
            self.adapt_formatspy(writeformatspy = copy_and_adapt_formatspy_anyhow)
        return

def is_upper_triangular(arr, atol=1e-8):
    """test for upper triangular matrix based on numpy"""
    # must be (n x n) matrix
    assert len(arr.shape)==2
    assert arr.shape[0] == arr.shape[1]
    return np.allclose(np.tril(arr, k=-1), 0., atol=atol)

def convert_cell(cell,pos):
    """
    Convert a parallelepipedal (forming right hand basis)
    to lower triangular matrix LAMMPS can accept. This
    function transposes cell matrix so the bases are column vectors
    """
    cell = np.matrix.transpose(cell)
    if not is_upper_triangular(cell):
        # rotate bases into triangular matrix
        tri_mat = np.zeros((3, 3))
        A = cell[:, 0]
        B = cell[:, 1]
        C = cell[:, 2]
        tri_mat[0, 0] = np.linalg.norm(A)
        Ahat = A / np.linalg.norm(A)
        AxBhat = np.cross(A, B) / np.linalg.norm(np.cross(A, B))
        tri_mat[0, 1] = np.dot(B, Ahat)
        tri_mat[1, 1] = np.linalg.norm(np.cross(Ahat, B))
        tri_mat[0, 2] = np.dot(C, Ahat)
        tri_mat[1, 2] = np.dot(C, np.cross(AxBhat, Ahat))
        tri_mat[2, 2] = np.linalg.norm(np.dot(C, AxBhat))

        # create and save the transformation for coordinates
        volume = np.linalg.det(cell)
        trans = np.array([np.cross(B, C), np.cross(C, A), np.cross(A, B)])
        trans /= volume
        coord_transform = np.dot(tri_mat, trans)
        pos = np.dot(coord_transform, pos.transpose())
        pos = pos.transpose()

        return tri_mat, pos
    else:
        return cell, pos


#def ase_get_known_formats(show=False, add_missing_formats=False, copy_formats=False, verbose=False,show_formatspy=False):
#    ''' adds formats runner and lammps-runner to ase '''
#
#    ### get the known formats
#    known_formats = []
#    x = ase.io.formats.all_formats
#    for i in x:
#        known_formats.append(i)
#
#    ### show the known formats
#    if show_known_formats:
#        import pprint
#        pp = pprint.PrettyPrinter(indent=4)
#        pp.pprint(x)
#
#    ### get formatspy
#    formatspy = os.path.dirname(ase.io.__file__)+"/formats.py"
#    if verbose or show_formatspy:
#        print('>> formatspy        :',formatspy)
#
#    ### check if formats are known by ase
#    missing = [ "runner.py","lammpsrunner.py", "lammpsdata.py", "ipi.py", "quippy.py" ]
#    def checkformats(typ,verbose):
#        if typ in known_formats:
#            if verbose: print(">> formats.py knows", typ)
#            return True
#        else:
#            if verbose: print(">> ERROR, formats.py does not know",typ)
#            return False
#
#    for i in [ 'runner', 'lammps-runner', 'lammps-data' ]:
#        formats_known = checkformats(i,verbose)
#        if formats_known == False: add_missing_formats = True
#
#
#    ### copies the missing format files
#    if copy_formats or add_missing_formats:
#        if verbose:
#            print('cc',ase.io.__file__)
#        scripts = my.scripts()
#        from_ = scripts+"/runner_scripts/ase_fileformat_for_"
#        to = os.path.dirname(ase.io.__file__)+"/"
#        for ff in missing:
#            #print('copying ',from_+ff,'to',to+ff)
#            if verbose: print('copying ',ff,'to',to+ff)
#            shutil.copyfile(from_+ff,to+ff)
#
#
#    ### check if necessary files for formats are known
#    if add_missing_formats:  # copies the missing format files
#        if verbose:
#            print('adapting ase formats.py .... ')
#
#
#        if not os.path.isfile(formatspy):
#            print('formatspy',formatspy)
#            sys.exit('did not find '+str(formatspy))
#
#
#        f = open(formatspy, "r")
#        contents = f.readlines()
#        f.close()
#        insert=0
#        insert2=0
#        for idx,i in enumerate(contents):
#            #print('i',idx,i)
#            #print("|"+i[:20]+"|")
#            if i[:20] == "    'abinit': ('ABIN":
#                insert = idx
#            if i[:30] == "    'lammps-data': 'lammpsdata":
#                insert2 = idx
#
#        writeformatspy = False
#        if 'runner' in x:
#            if verbose:
#                print('runner        format are already added in formats.py (of ase).')
#        else:
#            contents.insert(insert, "    'runner': ('Runner input file', '+F'),\n")
#            writeformatspy = True
#
#        if 'ipi' in x:
#            if verbose:
#                print('ipi           format are already added in formats.py (of ase).')
#        else:
#            contents.insert(insert, "    'ipi': ('ipi input file', '+F'),\n")
#            writeformatspy = True
#
#        if 'lammps-runner' in x:
#            if verbose:
#                print('lammps-runner format are already added in formats.py (of ase).')
#        else:
#            contents.insert(insert, "    'lammps-runner': ('LAMMPS data input file for n2p2 or runner', '1F'),\n")
#            contents.insert(insert2,"    'lammps-runner': 'lammpsrunner',\n")
#            writeformatspy = True
#
#        if writeformatspy == True:
#            print('now changing formatspy')
#            print('insert',insert)
#
#            f = open(formatspy, "w")
#            contents = "".join(contents)
#            f.write(contents)
#            f.close()
#        else:
#            if verbose:
#                print('everything was already in formats.py')
#
#    return known_formats
def inputnn_runner_or_n2p2(file):
    rn = grep(file,"runner_mode")
    #print('type)',type(rn),rn)
    if type(rn) == list and len(rn) == 0:
        return "n2p2"
    else:
        return "runner"

def inputdata_get_nuber_of_structures(inputdata):
    #print('inputdata',inputdata)
    inputdatanr = inputdata.replace("input.data","input.data_nr")
    #print('inputdatanr',inputdatanr)
    if os.path.isfile(inputdatanr):
        nr = np.loadtxt(inputdatanr)
        nr = int(nr)
        #print('from loaded file',nr)
        return nr
    else:
        if os.path.isfile(inputdata):
            nr = check_output(["grep","-c","begin",inputdata]).decode('utf-8')
            nr = int(nr)
            #print('from grep nr',nr,type(nr),int(nr))
            np.savetxt(inputdatanr,np.array([nr]))
        else:
            nr = 0
    return nr

def inputnn_get_potential_number_from_weightsfile(inputnn):
    weights = inputnn.replace('input.nn', 'weights.013.*.out')
    print('weights',weights)
    f = glob.glob(weights)
    print('f',f)
    if len(f) == 1:
        taken = f[0].split("weights.013.")[1].split(".out")[0]
        print('t',taken)
        leading_removed = [s.lstrip("0") for s in [taken]][0]
        # Remove leading
        print('lr',leading_removed)
        return int(leading_removed)
    return False

def inputnn_get_testfraction(file):
    test_fraction = np.float(grep(file,"test_fraction")[0].split()[1])
    return test_fraction

def inputnn_get_random_seed(file):
    random_seed = np.float(grep(file,"random_seed")[0].split()[1])
    return random_seed

def inputnn_get_nodes_short(file,as_string=False):
    nn = (grep(file,"global_nodes_short")[0]).split()
    #print(nn.index("#"))
    nna_str = nn[1:nn.index("#")]
    if as_string == True:
        return "_".join(nna_str)
    nna_int = list(map(int, nna_str))
    #print('nna_str:',nna_int)
    return np.array(nna_int)

def inputnn_get_activation_short(file):
    nn = (grep(file,"global_activation_short")[0]).split()
    nna = nn[1:nn.index("#")]
    nnb = "_".join(nna)
    return nnb

def inputnn_get_trainfraction(file):
    test_fraction = inputnn_get_testfraction(file)
    return test_fraction - 1.

def inputnn_get_atomic_symbols_and_atom_energy(inputnn,verbose=False):
    elements = []
    if os.path.isfile(inputnn):
        ##### get elements
        lines = grep(inputnn,"^elements")
        if len(lines) == 1:
            line = lines[0]
            line_elements_ = line.split()[1:]
            elements = []
            for i in line_elements_:
                #print(i)
                if i in my_atom.atomic_symbols:
                    #print("yo",i)
                    elements.append(i)
                else:
                    break
        if verbose:
            print('elements ++',elements)

        ##### get atomic energies
        lines = grep(inputnn,"^atom_energy")
        ele_list = []
        ene_list = []
        for i in lines:
            if i.split()[1] in elements:
                ele = i.split()[1]
                ene = float(i.split()[2])
                #print('lines',i.split(),"--------->>",ele,ene,type(ene))
                ele_list.append(ele)
                ene_list.append(ene)
                #print("ele_list",ele_list)
        if verbose:
            print('1 ele_final:',ele_list,len(ele_list))
            print('1 ene_final',ene_list,len(ene_list))

        if len(ele_list) == len(ene_list) == 0:
            ele_list = elements
        else:
            elements = ele_list

        if verbose:
            print('2 ele_final:',ele_list)
            print('2 ene_final',ene_list)
        if len(ene_list) == 0 and len(ele_list) != 0:
            ene_list = list(np.zeros(len(ele_list)))
        if verbose:
            print('3 ele_final:',ele_list)
            print('3 ene_final',ene_list)

        d = {}
        if len(ele_list) == len(ene_list):
            d = {}
            for idx,i in enumerate(ele_list):
                d[i] = ene_list[idx]
            elements = ele_list
            atom_energy = d
        if verbose:
            print("elements,",elements)
            print("atom_energy",atom_energy)
        return elements, atom_energy

def n2p2_runner_get_learning_curve(filename,only_get_filename=False,verbose=False):
    ''' filename is path to log.fit (runner) or learning-curve.out '''
    if verbose:
        print()
        print('filename in:',filename)
    basename = os.path.basename(filename)  # "learning-curve.out"
    folder = os.path.abspath(filename.replace(basename,''))
    type = inputnn_runner_or_n2p2(folder+'/input.nn')
    tryname = [ "logfiele_mode2", "log.fit", "logfile_mode2" ]
    changefilename = False

    if verbose:
        print('filename mid',filename)
        print('basename mid',basename)
        print('type     mid',type)

    if os.path.isfile(folder+'/optweights.012.out'): # and basename == "learning-curve.out":
        changefilename = True
    if os.path.isfile(folder+'/tmpweights.012.out'): # and basename == "learning-curve.out":
        changefilename = True
    if verbose:
        print('changefilename',changefilename)

    if changefilename == True:
        for i in tryname:
            filename = folder+"/"+i
            if os.path.isfile(filename):
                type = 'runner'
                break
    if not os.path.isfile(filename):
        sys.exit(filename+" does not exist!")

    basename = os.path.basename(filename)
    if verbose:
        print('filename out',filename)
        print('basename out',basename)
        print('type     out',type)
    if only_get_filename == True:
        return filename

    if basename == "learning-curve.out": # n2p2
        lc = np.loadtxt(filename) #+'/learning-curve.out')
        lc[:,1] = lc[:,1]*1000.*27.211384
        lc[:,2] = lc[:,2]*1000.*27.211384
        lc[:,3] = lc[:,3]*1000.*51.422063
        lc[:,4] = lc[:,4]*1000.*51.422063
    elif basename in tryname:          # runner
        f = open(filename, "r")
        contents = f.readlines()
        f.close()
        ene = []
        force = []
        all = []
        for idx,ii in enumerate(contents):
            if ii[:7] == " ENERGY":
                lst = ii.split()[1:4]
                #print('lst',lst)
                if lst[2] == 'NaN':
                    lst = [ '0','0','0']
                eneone = [float(iii) for iii in lst]
                ene.append(eneone)
                allone = [0,0,0,0,0]
                allone[0] = eneone[0]
                allone[1] = eneone[1]*1000.
                allone[2] = eneone[2]*1000.
            if ii[:7] == " FORCES":
                lst = ii.split()[1:4]
                forceone = [float(iii) for iii in lst]
                force.append(forceone)
                allone[3] = forceone[1]*1000.
                allone[4] = forceone[2]*1000.
                all.append(allone)
        ene = np.asarray(ene)
        force = np.asarray(force)
        all = np.asarray(all)
        lc = all
        #print(ene[:3])
        #print(force[:3])
        #print(all[:3])
        #print('lc',lc)
        #sys.exit()

        if len(lc.shape) == 1:
            lc = np.array([lc])

    if lc.shape == (1,0):
        lc = np.zeros((1,5))
    if type == 'runner':
        np.savetxt("learning-curve-runner.out",lc)
    #print('lc')
    #print(lc)
    #sys.exit()
    return lc


#######################################################################################
## SOAP funcions
#######################################################################################
def get_rawsoap(at, gs=0.3, co=3.0, cotw=0.5, nmax=9, lmax=6, cw=1.0, nrm=True, zlist=[], central_z=0, x_cmd=""):
    """
    gs gaussian widths
    co cutoff radius
    cotw you don't care
    n radial basis functions
    l spherical harmonics basis
    zlist list of atomic numbers the SOAPS will see
    central_z list of atomic numbers you will center on : IT DEFINES THE NUMBER OF ROWS of the rawsoaps

    """
    if len(zlist)==0:
        zlist = np.unique(np.asarray(at.z))
    nz = len(zlist)
    lspecies = "n_species="+str(nz)+" species_Z={"+(" ".join(np.asarray(zlist,dtype=str)))+"}"
    soapstr=("soap central_reference_all_species=F normalise="+("T" if nrm else "F")+" central_weight="+str(cw)+
               "  covariance_sigma0=0.0 atom_sigma="+str(gs)+" cutoff="+str(co)+" cutoff_transition_width="+str(cotw)+
               " n_max="+str(nmax)+" l_max="+str(lmax)+' '+lspecies+' Z='+str(central_z)+' '+x_cmd)
    desc = quippy.descriptors.Descriptor(soapstr )
    at.set_cutoff(co);
    at.calc_connect();
    return zlist, desc.calc(at)["descriptor"]

def get_soap2_vec(rawsoap, zlist, zglobal=None,nmax=False,lmax=False):
    isoap = 0
    isqrttwo = 1.0/np.sqrt(2.0)
    njsoap = {}; ipair = {}

    if zglobal is None:
        zglobal = zlist
    njsoap = np.zeros((len(zglobal),nmax,len(zglobal),nmax,lmax+1))
    zmap = np.zeros(len(zlist),int)
    for s in xrange(len(zlist)):
        zmap[s] = zglobal.index(zlist[s])
    for s1 in xrange(len(zlist)):
        z1 = zmap[s1]
        for n1 in xrange(nmax):
            for s2 in xrange(s1+1):
                z2 = zmap[s2]
                for n2 in xrange(nmax if s2<s1 else n1+1):
                    soap_element = rawsoap[isoap:isoap+lmax+1]
                    if (s1 != s2 or n1 != n2):
                        soap_element *= isqrttwo  # undo the normalization since we will actually sum over all pairs in all directions!
                    njsoap[z1,n1,z2,n2,:] = soap_element
                    njsoap[z2,n2,z1,n1,:] = soap_element
                    isoap+=lmax+1
    return njsoap

def do_fps(x, d=0):
    if d == 0 : d = len(x)
    n = len(x)
    iy = np.zeros(d, int)
    # faster evaluation of Euclidean distance
    n2 = np.sum(x**2,axis=1)
    iy[0] = np.random.randint(0, n)
    dl = n2 + n2[iy[0]] - 2* np.dot(x, x[iy[0]])
    lmin = np.zeros(d)
    for i in range(1,d):
        iy[i] = np.argmax(dl)
        lmin[i-1] = dl[iy[i]]
        if i%1000 ==0: print("max min dist %f" %( dl[iy[i]]))
        nd = n2 + n2[iy[i]] - 2*np.dot(x,x[iy[i]])
        dl = np.minimum(dl, nd)
    return iy, lmin

def get_soaps(kmcxyz = False, nmax = 8, lmax = 6, co = 4, gs = 0.5, zlist = [12,13,14], central_z=23,verbose=False,showtdqm=True ):
    rsoap = []
    #nmax = 8   # 12
    #lmax = 6  # 9
    #co = 4  # cutoff radius, 3.0  ## is this in angstrom?
    gs = 0.5 # gaussian width, 0.3
    cotw = gs
    det = "c"+str(co)+"-g"+str(gs)+"-nrm-cw"
    #zlist = [12,13,14]  # why is this working, even without specifying the vacancy/Vanadium ....? [12,13,14,23]
    ztot = len(zlist)
    nsoap = ztot**2 * nmax**2 * (1+lmax)
    ntot = len(kmcxyz)
    if verbose:
        print('rsoap ...')
    if showtdqm == True:
        gothrough = tqdm_notebook(kmcxyz)
    else:
        gothrough = kmcxyz

    for at in gothrough:
        #progress(i,len(kmcxyz))
        # andrea used a zentral_z=0
        # the next line could be potentially done parallel...
        zliat, soaps = get_rawsoap(at, nmax=nmax, lmax=lmax, co=co, gs=gs, cotw=cotw, nrm=False, cw=1, zlist=zlist, central_z=central_z,
                                   x_cmd="cutoff_dexp=2 cutoff_scale=3.0 cutoff_rate=2.0")
        rsoap.append((zliat, soaps))
    if verbose:
        print('rsoap done ...')

    soap2 = []
    if verbose:
        print('soap2 ...')
    for i in tnrange(ntot):
        #progress(i,ntot)
        izl, soaps = rsoap[i]
        lenv = np.zeros((len(soaps), ztot, nmax, ztot, nmax, 1+lmax))
        for s in range(len(soaps)):
            lenv[s] = get_soap2_vec(soaps[s], izl, zlist,nmax=nmax,lmax=lmax)
        soap2.append(lenv)
    if verbose:
        print('soap2 done ...')
        print('esoap2 ...')
    esoap2 = []
    for s in soap2:
        #progress(i,len(soap2))
        for e in s:
            esoap2.append(e.flatten())
    if verbose:
        print('esoap2 done ...')
    esoap2 = np.asarray(esoap2)
    return rsoap, esoap2


if __name__ == "__main__":
    pass
