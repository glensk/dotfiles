#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import os,sys,re,fnmatch
import filecmp
#import click
import numpy as np
import glob #,pathlib

try:
    from tqdm import tqdm_notebook, tnrange
except ImportError:
    pass

from copy import deepcopy
from socket import gethostname
import shutil
from subprocess import check_output,call
from datetime import datetime as datetime   # datetime.datetime.now()

#try:  # not in aiida ase
#    from ase.constraints import ExpCellFilter
#except ImportError:
#    pass

#from ase.spacegroup import crystal
#from ase.constraints import StrainFilter
#try:
#    from ase.constraints import ExpCellFilter
#except ImportError:
#    pass

import hesse
from feos import eos
import my_atom
#try:
#    from ase.optimize import GPMin
#except ImportError:
#    pass
import time
#import myutils as my

start_time = time.time()

##################################################################################
## genereal helper funcions
##################################################################################
def ptest():
    print("yes works 1")
    return

def printoutcolor(red,var,ENDC):
    if len(var) == 1:
        return red + str(var[0]) + ENDC
    else:
        return red + str(var) + ENDC

def is_int(x):
    try:
        a = float(x)
        b = int(a)
    except ValueError:
        return False
    else:
        return a == b

def isfiledir(file,exit=False,check_extension=False):
    if os.path.isfile(file):
        return file
    if check_extension != False:
        if type(check_extension) == str:
            check_extension = [check_extension]
        for ext in check_extension:
            if os.path.isfile(file+ext):
                return file+ext
            if os.path.isfile(file+"."+ext):
                return file+"."+ext

    if exit == True:
        sys.exit(file+" does not exist! (97)")
    return False

def read_lastline(file):
    with open(file, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        return f.readline().decode()

def printnormal(*var):
    ENDC = '\033[0m'
    return printoutcolor(ENDC,var,ENDC)

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

def print_args(args):
    ''' prints arguments from argparse '''
    keys = args.__dict__.keys()
    values = args.__dict__.values()
    print('########################################## argparse values (begin) #########')
    for idx,i in enumerate(keys):
        print("#",str(keys[idx]).ljust(35),str(values[idx]).ljust(20),str(type(values[idx])).ljust(15),"#")
    print('########################################## argparse values (end  ) #########')
    return

def check_for_known_hosts(exit=False):
    hostname = gethostname()   # fidis, helvetios, daint105, h332, g037, f221
    known_hosts =  ['fidis','helvetios', "daint", 'mac', 'cosmopc' ]
    for i in known_hosts:
        if i in hostname: return i
    if hostname[0] == 'h' and is_int(hostname[1:]) == True:
        return 'helvetios'
    if hostname[0] in ['g','f'] and is_int(hostname[1:]) == True:
        return 'fidis'

    # if not known host
    if exit == True:
        print("known hosts:",known_hosts)
        sys.exit(hostname+" is not in the list of known hosts!")
    else:
        return False

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
    getfromprompt = get_prompt_irrespective_of_python_version(text+" ")
    if getfromprompt == "":
        sys.exit()
    if getfromprompt[0] not in [ 'Y', 'y' ]:
        sys.exit('Exist since not Y or y as first letter!')
    return

def get_from_prompt_True_or_False(text):
    getfromprompt = get_prompt_irrespective_of_python_version(text+" ")
    #print('getfromprompt',getfromprompt,getfromprompt[0])
    if getfromprompt == "":
        return False
    if getfromprompt[0] in [ 'Y', 'y' ]:
        return True
    return False

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
                message = "ERROR "+message
                sys.exit(message)
            else:
                print(message)
    else:
        if not os.path.isdir(variable):
            message = 'directory '+str(var)+' is not defined or does not exist (21)'
            if exit == True:
                message = "ERROR "+message
                sys.exit(message)
            else:
                print(message)
    return variable

def hostname():
    hostname = gethostname()
    return hostname

#def get_click_defaults():
#    # show default values in click
#    orig_init = click.core.Option.__init__
#    def new_init(self, *args, **kwargs):
#        orig_init(self, *args, **kwargs)
#        self.show_default = True
#    click.core.Option.__init__ = new_init
#    CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'],token_normalize_func=str.lower)
#    return CONTEXT_SETTINGS

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

def q(verbose=False):
    id=[]
    stat=[]
    cores=[]
    runtime=[]
    path=[]

    host = check_for_known_hosts(exit=False)
    if host == False:
        return id,stat,path
    out=check_output(['q'])
    debug=False
    out2=out.split('\n')
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

    if verbose:
        print('----- currently in the que ----------------------')
        for idx,i in enumerate(id):
            print(idx,id[idx],stat[idx],path[idx])
        print('----- currently in the que done -----------------')
    return id,stat,path

def find_files(folder,searchpattern,maxdepth=False):
    ''' before :
    #import subprocess
    #files = check_output(["find "+folder+" -name \""+searchpattern+"\""],shell=True, stderr=subprocess.STDOUT)
    '''
    import utils_rename
    if maxdepth == False:
        files = utils_rename.run2(command='find '+folder+' -name "'+searchpattern+'"').split()
    else:
        files = utils_rename.run2(command='find '+folder+' -maxdepth '+str(maxdepth)+' -name "'+searchpattern+'"').split()
    #print('--------1',type(files))
    #for i in files:
    #    print(i)
    #print('--------1')
    #sys.exit()
    return files

def glob_recursively(folder,searchpattern):
    ''' searchpater e.g. "*.c"  '''
    matches = []
    for root, dirnames, filenames in os.walk(folder):
        #print('root',root)
        #print('dirnames',dirnames)
        #print('filenames',filenames)
        for filename in fnmatch.filter(filenames, searchpattern):
            matches.append(os.path.join(root, filename))
    return matches

##################################################################################
## quantum espresso funcions
##################################################################################
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

##################################################################################
## functions related to potential
##################################################################################
class mypot( object ):
    ''' return a list of available potentials '''
    def __init__(self,pot=False,potpath=False,use_different_epoch=False,verbose=False):
        self.pot                        = pot         # n2p2_v1ag
        self.potpath_in                 = potpath
        self.use_different_epoch        = use_different_epoch
        self.potpath                    = False        # this is the source where the potential recides
        self.potpath_work               = False        # this is ususally the self.potpath but for cases where different epoch is used ->
        self.pottype                    = False        # n2p2/runner
        self.potepoch_all               = False
        self.assessed_epochs            = []
        self.assessed_test              = []
        self.assessed_train             = []
        self.assessed_kmc57             = []
        self.assessed_c44               = []
        self.assessed_input             = []
        self.potepoch_bestteste         = False
        self.potepoch_bestteste_checked = False
        self.c44_al_file                = False
        self.c44_al                     = False
        self.potDONE                    = False         # n2p2_v1ag
        self.potlib                     = False         # n2p2_v1ag
        self.potcutoff                  = False         # n2p2_v1ag
        self.learning_curve_file        = False
        timestr                         = str(int(time.time()))
        self.lammps_tmpdir              = os.environ['HOME']+"/._tmp_lammps_"+str(gethostname())+"/" #+timestr+"/"
        self.pot_tmpdir                 = os.environ['HOME']+"/._tmp_pot_"   +str(gethostname())+"/" #+timestr+"/"
        self.inputnn                    = False         # path to input.nn
        self.inputdata                  = False         # path to input.data
        self.pot_all                    = False   # n2p2_v1ag
        self.trigger_set_path_manual    = ["setpath","potpath","pp", ".", "..", "../" ]
        self.verbose                    = verbose
        self.elements                   = False  # list e.g. ['Al', 'Mg', 'Si']
        self.atom_energy                = False
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
                self.elements, self.atom_energy = inputnn_get_atomic_symbols_and_atom_energy_dict(inputnn)

        return

    def print_variables_mypot(self,text="",print_nontheless=False):
        if self.verbose > 1 or print_nontheless:
            print(text,"self.pot                    ",self.pot)
            print(text,"self.potpath (source = src) ",self.potpath)
            print(text,"self.potpath_work (typ src) ",self.potpath_work)
            print(text,"self.potpath_in             ",self.potpath_in)
            print(text,"self.inputnn                ",self.inputnn)
            print(text,"self.inputdata              ",self.inputdata)
            print(text,"self.learning_curve_file    ",self.learning_curve_file)
            print(text,"self.potepoch_all           ",self.potepoch_all)
            print(text,"self.assessed_epochs        ",self.assessed_epochs)
            print(text,"self.assessed_test          ",self.assessed_test)
            print(text,"self.assessed_train         ",self.assessed_train)
            print(text,"self.assessed_kmc57         ",self.assessed_kmc57)
            print(text,"self.assessed_c44           ",self.assessed_c44)
            print(text,"self.use_different_epoch    ",self.use_different_epoch)
            print(text,"self.potepoch_bestteste     ",self.potepoch_bestteste,"(It was checked that this is the linked potential)")
            print(text,"self.potepoch_bestteste_chk?",self.potepoch_bestteste_checked)
            print(text,"self.c44_al_file            ",self.c44_al_file)
            if type(self.c44_al) != bool:
                print(text,"self.c44_al                 ",np.round(self.c44_al,2),"GPa")
            else:
                print(text,"self.c44_al                 ",self.c44_al)

            print(text,"self.pot_tmpdir             ",self.pot_tmpdir)
            print(text,"self.lammps_tmpdir          ",self.lammps_tmpdir)
            print(text,"self.potlib                 ",self.potlib)
            print(text,"self.potcutoff (Angstrom)   ",self.potcutoff)
            print(text,"self.elements               ",self.elements)
            print(text,"self.atom_energy            ",self.atom_energy)
            print(text,"self.pottype                ",self.pottype)
            print(text,"self.potDONE                ",self.potDONE)
            print(text,"self.verbose                ",self.verbose)
            print(text,"self.pot_all                ",self.pot_all)
            print()

    def get_my_assessments(self):
        self.test_b = self.test_l = self.train_b = self.train_l = self.kmc57_b = self.kmc57_l = False
        if len(self.potepoch_all) == 0:
            return
        #print('folder',self.potpath)
        #print('all',self.potepoch_all)
        epoch_last = self.potepoch_all[-1]
        #print('epoch_last',epoch_last)
        #self.print_variables_mypot('ka',print_nontheless=True)
        #print('epoch_best',self.potepoch_bestteste)

        ext_ = [ 'test', 'train', 'kmc57', 'input' ]
        # it may be that best == last or that only one epoch was shown in which case also
        # best == last
        epochs_ = [ str(self.potepoch_bestteste), str(epoch_last) ]

        #print('epochs_',epochs_)
        self.test_b = self.test_l = self.train_b = self.train_l = self.kmc57_b = self.kmc57_l = False
        allepochs = []
        for ext in ext_:
            filex = glob.glob(self.potpath+"/assess_"+ext+"_*")
            for i in filex:
                #print('i99',i)
                idx = i.replace(self.potpath+"/assess_"+ext+"_","")
                allepochs.append(int(idx))
        self.assessed_epochs = np.sort(list(set(allepochs)))
        self.assessed_test  = [None] * len(self.assessed_epochs)
        self.assessed_train = [None] * len(self.assessed_epochs)
        self.assessed_kmc57 = [None] * len(self.assessed_epochs)
        self.assessed_input = [None] * len(self.assessed_epochs)
        self.assessed_c44   = [None] * len(self.assessed_epochs)
        self.assessed_test_b  = "-"
        self.assessed_train_b = "-"
        self.assessed_kmc57_b = "-"
        self.assessed_input_b = "-"
        self.assessed_test_l  = "-"
        self.assessed_train_l = "-"
        self.assessed_kmc57_l = "-"
        self.assessed_input_l = "-"
        #print('allep',self.assessed_epochs)
        for ext in ext_:
            for eidx, epoch in enumerate(self.assessed_epochs):
                #print('ext',ext,eidx,epoch)
                file = self.potpath+"/assess_"+ext+"_"+str(epoch)+'/ene_std.npy'
                if ext == 'test':
                    if os.path.isfile(file): self.assessed_test[eidx] = float(my.read_lastline(file))
                    else: self.assessed_test[eidx] = "-"
                    if epoch == self.potepoch_bestteste: self.assessed_test_b = self.assessed_test[eidx]
                    if epoch == self.assessed_epochs[-1]: self.assessed_test_l = self.assessed_test[eidx]
                if ext == 'train':
                    if os.path.isfile(file): self.assessed_train[eidx] = float(my.read_lastline(file))
                    else: self.assessed_train[eidx] = "-"
                    if epoch == self.potepoch_bestteste: self.assessed_train_b = self.assessed_train[eidx]
                    if epoch == self.assessed_epochs[-1]: self.assessed_train_l = self.assessed_train[eidx]
                if ext == 'kmc57':
                    if os.path.isfile(file): self.assessed_kmc57[eidx] = float(my.read_lastline(file))
                    else: self.assessed_kmc57[eidx] = "-"
                    if epoch == self.potepoch_bestteste: self.assessed_kmc57_b = self.assessed_kmc57[eidx]
                    if epoch == self.assessed_epochs[-1]: self.assessed_kmc57_l = self.assessed_kmc57[eidx]
                if ext == 'input':
                    if os.path.isfile(file): self.assessed_input[eidx] = float(my.read_lastline(file))
                    else: self.assessed_input[eidx] = "-"
                    if epoch == self.potepoch_bestteste: self.assessed_input_b = self.assessed_input[eidx]
                    if epoch == self.assessed_epochs[-1]: self.assessed_input_l = self.assessed_input[eidx]

        if os.path.isfile(self.potpath+"/elastic_c44_all.dat"):
            elastic_c44_all = np.loadtxt(self.potpath+"/elastic_c44_all.dat")
            #print(elastic_c44_all)
            #print(len(elastic_c44_all))
            for eidx, epoch in enumerate(self.assessed_epochs):
                if len(elastic_c44_all) >= epoch:
                    #print(eidx,epoch,elastic_c44_all[epoch-1])
                    self.assessed_c44[eidx] = elastic_c44_all[epoch-1][1]
                else:
                    self.assessed_c44[eidx] = -1
            #sys.exit()
        #print('self.assessed_test',self.assessed_test)
        #print('self.assessed_train',self.assessed_train)
        #sys.exit()
                    #if ext == 'test' and eidx == 0: self.test_b = float(my.read_lastline(file))
                    #if ext == 'test' and eidx == 1: self.test_l = float(my.read_lastline(file))
                    #if ext == 'train' and eidx == 0: self.train_b = float(my.read_lastline(file))
                    #if ext == 'train' and eidx == 1: self.train_l = float(my.read_lastline(file))
                    #if ext == 'kmc57' and eidx == 0: self.kmc57_b = float(my.read_lastline(file))
                    #if ext == 'kmc57' and eidx == 1: self.kmc57_l = float(my.read_lastline(file))
        return


    def get_my_assessments_check_outliers(self,specific_epoch = False,greater = 60,verbose=False):
        ''' greater 60 sets outlier when diff (NN-DFT) > 60meV/atom '''
        has_outliers = False
        has_outliers_ = "  "
        struct_idx = []
        diffs = []
        file = False
        checked_epochs = []
        for ext in [ 'test' ]:  # only check for outliers in test!
            checked_epochs  = self.assessed_epochs
            if specific_epoch != False:
                checked_epochs  = [specific_epoch]
            for epoch in checked_epochs :
                file = self.potpath+"/assess_"+ext+"_"+str(epoch)+'/ene_diff_abs.npy'
                if not os.path.isfile(file):
                    has_outliers = True
                    has_outliers_ = "??O"
                    # getEnergies_byLammps.py -p . -ctest -pe 1244
                    # getEnergies_byLammps.py -p . -ctest -pe epoch
                    return has_outliers,has_outliers_,checked_epochs,struct_idx,diffs
                else:
                    greater = 60 # meV/atom
                    if ext == 'test':
                        c = np.loadtxt(file)
                        #print(c)
                        struct_idx = np.where(c > greater)[0]
                        if len(struct_idx) > 0:
                            has_outliers = True
                            has_outliers_ = "!!O"
                            diffs = np.round(np.array(c[struct_idx]),1)
                            if verbose:
                                print('vOv diff greater >',greater,'in epoch',epoch,'| struct_idx:',struct_idx,'| diff',diffs,'file',file)
                            for gf in struct_idx:
                                frame = ase_read(self.testdata,index=gf,format='runner')
                                if verbose:
                                    print('frame',gf,frame.info['comment'],'number_of_atoms:',frame.get_number_of_atoms())
                                #print('now check also in original input.data if this struct is there')
                            #sys.exit()
        #print('...',self.test_b,self.test_l,self.train_b,self.train_l,self.kmc57_b,self.kmc57_l)
        #return has_outliers,checked_epochs,struct_idx,diffs
        return has_outliers,has_outliers_,checked_epochs,struct_idx,np.clip(diffs,0,999.9)

    def get(self,exit=True,showerrors=True):
        ''' exit is True if for instance weithts files do not exist '''
        self.potepoch_bestteste_checked = False
        self.print_variables_mypot('PPk get potential: in')

        ##########################################
        # get potential from path
        ##########################################
        if self.potpath == False and self.potpath_in == False and self.pot in [ ".." , "../", "." ]:
            self.potpath_in = os.path.abspath(self.pot)

        if self.potpath_in != False and self.potpath == False:
            if not os.path.isdir(self.potpath_in):
                sys.exit(self.potpath_in+" does not exist! (1)")

            if not os.path.isfile(self.potpath_in+"/input.nn"):
                print("PPa could not find "+self.potpath_in+"/input.nn")
                self.potpath_in = "/".join(self.potpath_in.split("/")[:-1])
                if not os.path.isfile(self.potpath_in+"/input.nn"):
                    print("PPb could not find "+self.potpath_in+"/input.nn")
                    self.potpath_in = "/".join(self.potpath_in.split("/")[:-1])
                    if not os.path.isfile(self.potpath_in+"/input.nn"):
                        sys.exit("PPc could not find "+self.potpath_in+"/input.nn")

            if exit == True:
                checkfiles = [ "input.nn", "scaling.data", "weights.012.data", "weights.013.data", "weights.014.data" ]
                for i in checkfiles:
                    if not os.path.isfile(self.potpath_in+"/"+i):
                        sys.exit(self.potpath_in+"/"+i+" does not exist! (2)")

            self.potpath = os.path.abspath(self.potpath_in)
            self.inputnn     = isfiledir(self.potpath_in+"/input.nn",exit=True)
            self.scalingdata = isfiledir(self.potpath_in+"/scaling.data",exit=True)

            self.inputdata = isfiledir(self.potpath_in+"/input.data",check_extension=[".extxyz"])
            self.testdata  = isfiledir(self.potpath_in+"/test.data",check_extension=[".extxyz"])
            self.traindata = isfiledir(self.potpath_in+"/train.data",check_extension=[".extxyz"])


            self.pottype = inputnn_runner_or_n2p2(self.inputnn)
            self.learning_curve_file = n2p2_runner_get_learning_curve_filename(self.inputnn)
            try:
                self.potepoch_all = inputnn_get_potential_number_from_all_weightsfiles(self.inputnn)
            except NameError:
                self.potepoch_all = []
            self.potepoch_bestteste = n2p2_runner_get_bestteste_idx(self.inputnn)

            #### check if current weights.xxx.data files are the ones from weights.xxx. self.potepoch_bestteste
            if exit == True:
                file1 = self.potpath_in+"/weights.012.data"
                epstr = str(self.potepoch_bestteste).zfill(6)
                file2 = self.potpath_in+"/weights.012."+epstr+".out"
                if not filecmp.cmp(file1, file1):
                    sys.exit("PPd File "+file1+" is not "+file2)
                else:
                    self.potepoch_bestteste_checked = True

            self.pot = self.pottype+"_frompath"
        else:
            ##########################################
            # get potential from string
            ##########################################
            self.get_potpath()
        self.get_elements_and_atomic_energies()


        ### check if self.pottype can be computed on this host!
        add = 'PPe Your lammps version does not seem to work with '+self.pottype+"!"
        if self.pottype == "runner":
            self.potlib = os.environ["LAMMPSPATH"]+"/src/USER-RUNNER"
            self.potcutoff = 14.937658735
        if self.pottype == "n2p2":
            self.potlib = os.environ["LAMMPSPATH"]+"/src/USER-NNP"
            self.potcutoff = 10.6  # the minimum cutoff of the SF was 20bohrradius = 10.58Angstrom


        if self.pottype in [ "runner", "n2p2" ] and os.path.isdir(self.potlib) == False and showerrors == True:
            #sys.exit("ERROR: "+self.potlib+" not found!"+add)
            print("ERROR: "+self.potlib+" not found!"+add)

        if self.use_different_epoch == False:
            self.potpath_work = self.potpath
            #print('se false')
        elif int(self.use_different_epoch) == self.potepoch_bestteste:
            self.potpath_work = self.potpath
        else:  # use_different_epoch
            if self.verbose:
                print("PPf copy files ...")
            #print('se true')
            self.potpath_work = self.pot_tmpdir
            #print('self.potpath_work (1)',self.potpath_work)
            epstr = str(self.use_different_epoch).zfill(6)
            if not os.path.isdir(self.pot_tmpdir):
                mkdir(self.pot_tmpdir)
            f12 = self.potpath+"/weights.012."+epstr+".out"
            f13 = self.potpath+"/weights.013."+epstr+".out"
            f14 = self.potpath+"/weights.014."+epstr+".out"

            f12a = self.potpath+"/../_weights/weights.012."+epstr+".out"
            f13a = self.potpath+"/../_weights/weights.013."+epstr+".out"
            f14a = self.potpath+"/../_weights/weights.014."+epstr+".out"
            f12b = self.potpath+"/../weights.012."+epstr+".out"
            f13b = self.potpath+"/../weights.013."+epstr+".out"
            f14b = self.potpath+"/../weights.014."+epstr+".out"
            if not os.path.isfile(f12):
                #print("f12 dne",f12)
                f12 = f12a
                #print("f12 new",f12)
            if not os.path.isfile(f13): f13 = f13a
            if not os.path.isfile(f14): f14 = f14a
            if not os.path.isfile(f12): f12 = f12b
            if not os.path.isfile(f13): f13 = f13b
            if not os.path.isfile(f14): f14 = f14b
            for ff in [f12,f13,f14]:
                if not os.path.isfile(ff):
                    sys.exit(ff+" does not exist! (65)")
            my.cp(f12,self.pot_tmpdir+"/weights.012.data")
            my.cp(f13,self.pot_tmpdir+"/weights.013.data")
            my.cp(f14,self.pot_tmpdir+"/weights.014.data")
            my.cp(self.scalingdata,self.pot_tmpdir)
            my.cp(self.inputnn,self.pot_tmpdir)
            print('... copying potential',epstr,'to',self.pot_tmpdir)
            #print('self.potpath_work (2)',self.potpath_work)
        self.potDONE = True
        #print('self.potpath_work (3)',self.potpath_work)
        self.print_variables_mypot('get potential: out')
        return

def pot_all():
    all = mypot()
    all.get_pot_all()
    return all.pot_all

##################################################################################
## ipi related functions
##################################################################################
def create_ipi_kmc_inputfile(jobdir,filename="input-runner.xml", nsteps=False,stride=100,seed=12345,a0=4.057,ncell=4,nsi=3,nmg=3,nvac=1,neval=8,temp=300,verbosity="low",checkpoint_restart_stride=10,nodes=1,address="hostname_jobdir",testrun=False):
    '''
        nodes is only necessary to define ffsockets
        for calculations on one node, use ffsocket = "unix"; when > 1 noes: ffsocket = "inet"
        for verbosity = "high" every socket connect is reported ... which seems too much info
    '''
    if testrun == True:
        nsteps = 3
        strie = 1
        verbosity = 'low'

    if nodes == 1:
        ffsocket = "unix"
        addressline = '        <address> '+str(address)+' </address> <latency> 1e-3 </latency>'
    elif nodes > 1:
        ffsocket = "inet"
        addressline = '        <address> '+str(address)+' </address> <port> 12345 </port>'
    else:
        sys.exit("Number of nodes has to be positive!")

    activelist = str(range(ncell**3 - nvac))
    activelist = activelist.replace(" ", "")
    atom_x_list = []
    for nvac_idx in range(ncell**3 - nvac,ncell**3):
        atom_x_list.append('atom_x{angstrom}('+str(nvac_idx)+')')
    insert = ", ".join(atom_x_list)

    ipi_cmd = [
    '<simulation verbosity="'+str(verbosity)+'">',
    '    <output prefix="simulation">',
    '        <properties stride="1" filename="out"> [ step, time{picosecond}, potential] </properties>',
    '        <properties stride="1" filename="vac"> [ step, time{picosecond}, '+str(insert)+' ] </properties>',
    '        <trajectory filename="pos" stride="'+str(stride)+'" cell_units="angstrom" format="xyz" bead="0"> positions{angstrom} </trajectory>',
    '        <checkpoint stride="'+str(checkpoint_restart_stride)+'" filename="checkpoint.restart"/>',
    '    </output>',
    '    <total_steps>'+str(nsteps)+'</total_steps>',
    '    <prng> <seed>'+str(seed)+'</seed> </prng>',
    '    <ffsocket name="lmpserial" mode="'+str(ffsocket)+'">',
    addressline,
    #'        <!-- <activelist> '+activelist+' </activelist> -->',
    '        <activelist> '+activelist+' </activelist>',
    '    </ffsocket>',
    '    <total_time> 258800 </total_time>',
    '    <system>',
    '        <initialize nbeads="1">',
    #'            <file mode="xyz" units="angstrom"> al10x10x10_alat4.057_988al_6si_5mg_1va_999atoms.ipi </file>',
    '            <file mode="xyz" units="angstrom"> data.ipi </file>',
    '            <velocities mode="thermal" units="kelvin"> 0 </velocities>',
    '        </initialize>',
    '        <forces>',
    '            <force forcefield="lmpserial"> </force>',
    '        </forces>',
    '        <motion mode="al-kmc">',
    '            <al6xxx_kmc>',
    '               <geop mode="lbfgs">',
    '                    <ls_options> <iter> 3 </iter> </ls_options>',
    '               </geop>',
    '               <a0 units="angstrom"> '+str(a0)+' </a0>',
    '               <nstep> 5 </nstep>',
    '               <ncell> '+str(ncell)+' </ncell>',
    '               <nsi> '+str(nsi)+' </nsi>',
    '               <nmg> '+str(nmg)+' </nmg>',
    '               <nvac> '+str(nvac)+' </nvac>',
    '               <neval> '+str(neval)+' </neval>',
    '               <diffusion_barrier_al units="electronvolt"> 0.52  </diffusion_barrier_al>',
    '               <diffusion_prefactor_al units="terahertz"> 16.6 </diffusion_prefactor_al>',
    '                                                    <!-- from Mantina2009 dHm && v*  -->',
    '               <ecache_file> KMC_ECACHE </ecache_file>',
    '               <qcache_file> KMC_QCACHE </qcache_file>',
    '            </al6xxx_kmc>',
    '        </motion>',
    '        <ensemble>',
    '            <temperature units="kelvin">'+str(temp)+'</temperature>',
    '        </ensemble>',
    '    </system>',
    '</simulation>'
    ]

    if not os.path.isdir(jobdir):
        mkdir(jobdir)
    f = open(jobdir+'/'+filename,'w')
    for i in ipi_cmd:
        f.write(i+"\n")
    return
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

    #IPI_COMMAND    = test_and_return_environment_var_path('IPI_COMMAND')
    IPI_COMMAND     = test_and_return_environment_var_path('IPI_COMMAND_PLAY')
    LAMMPS_COMMAND  = get_LAMMPS_executable(exit=True)
    N2P2_PATH       = test_and_return_environment_var_path('N2P2_PATH',path=True)

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
    "#SBATCH --time=00-"+str(submittime_hours)+":00:00"]

    hostname = gethostname()
    if hostname == 'fidis' or hostname[0] == 'f':
        text2 = text2 + ["#SBATCH --constraint=E5v4"]

    text2 = text2 + [""]

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
    'seed=`grep "<seed>" input-runner.xml | sed "s|.*<seed>||" | sed "s|</seed>.*||" | awk \'{print $1}\'`',
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
    '',
    #'sleep 10',
    '# depending on the size fo the cell, this is necessary to initialize ipi',
    'for i in `seq 100`;do',
    '   [ ! -e "KMC_AL6XXX" ] && echo $i sleep 10 && sleep 100',
    '   [ -e "KMC_AL6XXX" ] && echo $i yes && break',
    'done',
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
    'echo "before wait text4" `pwd`',
    'wait',
    'echo "after wait text4" `pwd`',
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
    'echo "before wait text4o" `pwd`',
    'wait'
    'echo "after wait text4o" `pwd`',
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

##################################################################################
## crystal structure / atom cell  functions material scienc related
##################################################################################
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

def submitjob(submit_to_que=True,submit_to_debug_que=False,jobdir=False,submitskript=False):
    if jobdir == False:
        jobdir = os.getcwd()
    hier = os.getcwd()

    if submit_to_que or submit_to_debug_que:
        if submitskript == False:
            sys.exit("Error. please provite the submitskript")
        os.chdir(jobdir)
        command = ["sbatch"]
        if submit_to_debug_que == True:
            command = command + ["-p","debug","-t","01:00:00"]
        command = command + [submitskript]
        print('sbatch '+submitskript)
        call(["sbatch",submitskript])

    os.chdir(hier)
    return


##################################################################################
## lammps functions
##################################################################################
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

def ipi_write_static_inputfile(folder=False,filename='input.xml',positions='init.xyz',ace=False):
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
    ipi_write_static_inputfile(folder=tmpdir,filename='input.xml',positions='init.xyz',ace=ace)
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

def get_LAMMPS_executable(exit=True,verbose=False):
    SCR = os.environ["scripts"]
    onhost = gethostname()
    if verbose:
        print('SCR',SCR)
        print('onhost', onhost)
    LAMMPS_COMMAND = SCR+"/executables/lmp_"+onhost
    if verbose:
        print("LAMMPS_COMMAND",LAMMPS_COMMAND)
    if os.path.isfile(LAMMPS_COMMAND):
        return LAMMPS_COMMAND
    else:
        if onhost[0] in ['f','g','h']:
            if verbose:
                print('onhost[0]',onhost[0])
            checkint = onhost[1:]
            if is_int(checkint):
                return  SCR+"/executables/lmp_fidis"  # fidis and helvetios have
    if exit == True:
        sys.exit("ERROR: could not fine LAMMPS_executable")
    return

def lammps_ext_calc(atoms,ace,get_elastic_constants=False):
    ''' atoms is an ase atoms object which can hold several frames, or just one'''

    ###############################################################
    # find LAMMPS executable
    ###############################################################
    if ace.verbose:
        print("get_elastic_constants",get_elastic_constants)

    ###############################################################
    # make ace.pot.lammps_tmpdir (/home/glensk/._tmp_lammps/)
    ###############################################################
    #ace.pot.lammps_tmpdir = os.environ['HOME']+"/._tmp_lammps/"
    if not os.path.isdir(ace.pot.lammps_tmpdir):
        mkdir(ace.pot.lammps_tmpdir)
    #folder = tmpdir

    # delete all previous file in ace.pot.lammps_tmpdir
    for the_file in os.listdir(ace.pot.lammps_tmpdir):
        file_path = os.path.join(ace.pot.lammps_tmpdir, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            #elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

    print('doing the calculation in lammps, externally, in ace.pot.lammps_tmpdir',ace.pot.lammps_tmpdir)

    ###############################################################
    # write input structure (pos.lmp)
    ###############################################################
    if ace.verbose > 2:
        show_ase_atoms_content(atoms,showfirst=10,comment="START LAMMPS EXTERNALLY")
    atoms.set_calculator(None)
    atoms.write(ace.pot.lammps_tmpdir+'pos.lmp',format='lammps-runner')
    if ace.verbose:
        print('written ',ace.pot.lammps_tmpdir+'/pos.lmp')

    ###############################################################
    # write in.lmp (for mormal energy calculation)
    ###############################################################
    if not get_elastic_constants:
        execute_file = 'in.lmp'
        lammps_write_inputfile(folder=ace.pot.lammps_tmpdir,filename=execute_file,positions='pos.lmp',ace=ace)
        if ace.verbose:
            print('written ',ace.pot.lammps_tmpdir+'/'+execute_file)
    if ace.verbose > 1:
        print("written lammsp inputfile to ",ace.pot.lammps_tmpdir)

    ### calculate with lammps (trigger externally)
    ene = False
    if os.path.isfile(ace.pot.lammps_tmpdir+'log.lammps'):
        os.remove(ace.pot.lammps_tmpdir+"log.lammps")


    ###############################################################
    # IF get elastic constants
    ###############################################################
    if get_elastic_constants:
        execute_file = 'in.elastic'
        if ace.verbose:
            print('written ',ace.pot.lammps_tmpdir+'/potential.mod')
        lammps_write_inputfile_from_command(folder=ace.pot.lammps_tmpdir,filename='potential.mod',command=lammps_ext_elastic_potential_mod(ace))
        my.cp(scripts()+'/lammps_scripts/elastic/in.elastic',ace.pot.lammps_tmpdir+'/in.elastic')
        my.cp(scripts()+'/lammps_scripts/elastic/displace.mod',ace.pot.lammps_tmpdir+'/displace.mod')
        lammps_write_inputfile_from_command(folder=ace.pot.lammps_tmpdir,filename='init.mod',command=lammps_ext_elastic_init_mod(ace,positions='pos.lmp'))
        if ace.verbose:
            print('written ',ace.pot.lammps_tmpdir+'/'+execute_file)

    ###############################################################
    # cd to folder and run lammps
    ###############################################################
    with cd(ace.pot.lammps_tmpdir):  # this cd's savely into folder
        # RUN LAMMPS
        # without SHELL no LD LIBRARY PATH
        #print('pwd',os.getcwd())
        if ace.verbose:
            print(ace.LAMMPS_COMMAND+" < "+execute_file+" > /dev/null")
        call([ace.LAMMPS_COMMAND+" < "+execute_file+" > /dev/null"],shell=True)

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



##################################################################################
## n2p2/runner functions
##################################################################################
def inputnn_runner_or_n2p2(file):
    with open(file) as fp:
        for i, line in enumerate(fp):
            if "NNP" in line:
                return 'n2p2'
            if "RuNNer" in line:
                return 'runner'
            elif i > 3:
                break

    ## old way
    # rn = grep(file,"runner_mode")
    # #print('type)',type(rn),rn)
    # if type(rn) == list and len(rn) == 0:
    #     return "n2p2"
    # else:
    #     return "runner"
    return False

def inputdata_get_nuber_of_structures(inputnn):
    inputdatanr = inputnn.replace("input.nn","input.data_nr")
    inputdata   = inputnn.replace("input.nn","input.data")
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

def inputnn_get_potential_number_from_all_weightsfiles(inputnn):
    weights = inputnn.replace('input.nn', 'weights.013.*.out')
    #print('weights',weights)
    f = glob.glob(weights)
    #print('f',f)
    #sys.exit()
    if len(f) == 1:
        taken = f[0].split("weights.013.")[1].split(".out")[0]
        #print('t',taken)
        leading_removed = [s.lstrip("0") for s in [taken]][0]
        # Remove leading
        if leading_removed == "":
            leading_removed = 0
        #print('lr',leading_removed)
        return [int(leading_removed)]
    elif len(f) > 1:
        all = []
        for idx in np.arange(len(f)):
            taken = f[idx].split("weights.013.")[1].split(".out")[0]
            #print('t',taken)
            leading_removed = [s.lstrip("0") for s in [taken]][0]
            # Remove leading
            if leading_removed == "":
                leading_removed = 0
            all.append(int(leading_removed))
        return np.sort(np.asarray(all))
    else: # len(f) == 0:
        raise NameError('no weights files found')
    return

def inputnn_get_testfraction(file):
    test_fraction = np.float(grep(file,"test_fraction")[0].split()[1])
    return test_fraction

def inputnn_get_random_seed(file):
    random_seed = np.float(grep(file,"random_seed")[0].split()[1])
    return int(random_seed)

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
    #print('nn',nn)
    nni = nn[1:]
    #print('nni',nni)
    nna = []
    known = [ 't', 'p', 'l', 's']
    for i in nni:
        #print('i',i)
        if i in known:
            nna.append(i)
    #print('nna',nna)
    #try:
    #    nna = nni[nni.index("#")]
    #except ValueError:
    #    nna = nni
    #print('nna',nna,len(nna))
    nnb = "_".join(nna)
    #print('nnb',nnb,len(nnb))
    if len(nnb) != 5:
        sys.exit()
    return nnb

def inputnn_get_trainfraction(file):
    test_fraction = inputnn_get_testfraction(file)
    return test_fraction - 1.

def inputnn_get_atomic_symbols_and_atom_energy_dict(inputnn,verbose=False):
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

def inputnn_get_atomic_symbols_and_atom_energy_list(inputnn,verbose=False):
    pot_elements, pot_atom_energy = my.inputnn_get_atomic_symbols_and_atom_energy_dict(inputnn)
    try:
        mg = pot_atom_energy["Mg"]*-1
    except KeyError:
        mg = 0
    try:
        al = pot_atom_energy["Al"]*-1
    except KeyError:
        al = 0
    try:
        si = pot_atom_energy["Si"]*-1
    except KeyError:
        si = 0

    return ["Al","Mg","Si"], [al,mg,si]

def n2p2_runner_get_learning_curve_filename(inputnn):
    runner_n2p2 = inputnn_runner_or_n2p2(inputnn)
    folder = os.path.abspath(inputnn.replace('input.nn',''))
    if runner_n2p2 == 'n2p2':
        return folder+'/learning-curve.out'
    elif runner_n2p2 == 'runner':
        tryname = [ "logfiele_mode2", "log.fit", "logfile_mode2" ]
        for i in tryname:
            filename = folder+"/"+i
            if os.path.isfile(filename):
                return filename
    return False

def n2p2_check_SF_inputnn(inputnn):
    ''' inputnn can also be a file containing the symmetry functions only'''
    #inputnn = "../n2p2_v3ag_5000_new_2424_new_atomene_new_SF/get_scaling/cursel_64.def"
    out = grep(inputnn,"symfunction_short.*Si.*3.*Al.*Si") # echo $out | tail -1 | wc -w
    for line in out:
        if len(line.split()) != 9:
            print(line.split())
            print(len(line.split()))
            sys.exit("ERROR "+inputnn+" has 10 entries in a 3body line, should have 9! Exit!")
    return

def n2p2_runner_get_learning_curve(inputnn,only_get_filename=False,verbose=False):
    ''' filename is path to log.fit (runner) or learning-curve.out '''
    filename = n2p2_runner_get_learning_curve_filename(inputnn)
    if verbose:
        print('learning_curve_filename:',filename)
    basename = os.path.basename(filename)  # "learning-curve.out"
    folder = os.path.abspath(filename.replace(basename,''))
    if False: #verbose:
        print('learning_curve folder  :',folder)
        print('learning_curve basename:',basename)
    n2p2_runner = type = inputnn_runner_or_n2p2(folder+'/input.nn')
    if False: #verbose:
        print('nn',n2p2_runner)
    if not os.path.isfile(filename):
        sys.exit(filename+" does not exist! (32)")

    finished = False
    if n2p2_runner == "n2p2": # basename == "learning-curve.out": # n2p2
        lc = np.loadtxt(filename) #+'/learning-curve.out')
        #print('lc')
        #print(lc)
        #print(lc.shape)
        #print(len(lc.shape))
        if len(lc.shape) == 1:
            lc = np.array([lc])
        lc[:,1] = lc[:,1]*1000.*27.211384
        lc[:,2] = lc[:,2]*1000.*27.211384
        lc[:,3] = lc[:,3]*1000.*51.422063
        lc[:,4] = lc[:,4]*1000.*51.422063
    elif n2p2_runner == 'runner': #basename in tryname:          # runner
        runner_lc = folder+"/learning-curve-runner.out"
        if os.path.isfile(runner_lc):
            lc = np.loadtxt(runner_lc)
        else: # make runner_lc
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
                if ii[:48] == " Best short range fit has been obtained in epoch":
                    finished = True
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

    #print('aa',lc.shape,len(lc.shape))
    if len(lc.shape) == 1:
        #print('bb',lc)
        lc = np.array([lc])
    if lc.shape == (1,0):
        lc = np.zeros((1,5))
    if n2p2_runner == 'runner' and finished:
        np.savetxt(folder+"/learning-curve-runner.out",lc)
    #else:
    #    #os.remove(folder+"/learning-curve-runner.out")
    #    print('finished?',finished,folder+'/learning-curve-runner.out')
    #print('lc')
    #print(lc)
    #sys.exit()
    #print('ll',lc.shape,len(lc.shape))
    return lc

def n2p2_runner_job_finished(inputnn):
    #print('inputnn',inputnn)
    pp = os.environ["potentials"]
    #print('pp',pp,len(pp))
    if inputnn[:len(pp)] == pp:
        # job is already in potentials folder
        return True
        #print('same')
    else:
        #print('ns',pp)
        #print('ns',inputnn[:len(pp)])
        return False
    return False

def n2p2_runner_get_bestteste_idx(inputnn):
    #print('inputnn',inputnn)
    job_finished = n2p2_runner_job_finished(inputnn)
    best_testsete_file = inputnn.replace('input.nn', 'best_testsete')
    #print('best_testsete_file',best_testsete_file)
    if os.path.isfile(best_testsete_file):
        best_testsete = int(np.loadtxt(best_testsete_file))
        #print('from best_testsete_file',best_testsete_file)
        return best_testsete
    else:
        #print('from lc')
        learning_curve = lc = n2p2_runner_get_learning_curve(inputnn)
        best_testsete = np.argmin(lc[:,2])
        if job_finished:
            # np.savetxt(best_testsete_file,aa)
            np.savetxt(best_testsete_file,np.array([int(best_testsete)]),fmt='%i')
        #print(lc)
        #print('mmm',best_testsete)
        return best_testsete
        #save best_sestset_file
    return False

def n2p2_get_scaling_and_function_data(cores=28,days=0,hours=0,minutes=5,submit_to_que=False,submit_to_debug_que=True,interactive=False):
    ''' interactive == False -> que
        interactive == True  -> execture directly
        3 minutes is not enough for repeated structures (input.data ~ 130MB) but 4 minutes is enough
    '''
    if not os.path.isfile("input.data"):
        sys.exit("Need input.data file")
    if not os.path.isfile("input.nn"):
        sys.exit("Need input.nn file")
    n2p2_check_SF_inputnn("input.nn")

    hier=os.getcwd()
    if interactive == submit_to_que == True or interactive == submit_to_debug_que == True:
        sys.exit('either submit_to_{debug}_que or interacive')

    folder="get_scaling"
    if os.path.isdir(folder):
        sys.exit(folder+" already exists!")

    mkdir(folder)
    hier=os.getcwd()
    create_READMEtxt()
    os.chdir(folder)
    cp("../input.data")
    cp("../input.nn")
    submitskript = n2p2_write_submit_skript(directory=False,cores=cores,nodes=1,job="scaling",days=days,hours=hours,minutes=minutes,interactive=interactive)
    submitjob(submit_to_que=submit_to_que,submit_to_debug_que=submit_to_debug_que,jobdir=False,submitskript=submitskript)
    #submitjob(submit=True,submit_to_que=submit_to_que,submit_to_debug_que=submit_to_debug_que,jobdir=False,submitskript=submitskript)
    create_READMEtxt(add="submit_to_debug_que= "+str(submit_to_debug_que))
    os.chdir(hier)
    return

def n2p2_write_submit_skript(directory=False,nodes=1,cores=28,job=False,interactive=False,days=0,hours=72,minutes=0,seconds=0):
    ''' wiretes a submit_n2p2_{get_scaling,training}.sh file '''
    if job not in ["scaling","train"]:
        sys.exit('job has to be one of nnp-XXX jobs as "train, scaling, ..."')
    myhost = check_for_known_hosts(exit=True)
    cores_per_node = 28
    if myhost in [ "fidis", "helvetios"]:
        cores_per_node = 28
    elif myhost == 'daint':
        cores_per_node = 36


    if directory == False:
        directory = os.getcwd()

    # name of file
    filepath = directory+'/submit_n2p2_'+job+'.sh'

    # write file
    with open(filepath, "w") as text_file:
        text_file.write("#!/bin/bash\n")

        if interactive == False:
            text_file.write("#SBATCH --job-name=NNP-mpi\n")
            text_file.write("#SBATCH --output=_scheduler-stdout.txt\n")
            text_file.write("#SBATCH --error=_scheduler-stderr.txt\n")
            text_file.write("#SBATCH --nodes="+str(nodes)+"\n")
            text_file.write("#SBATCH --ntasks "+str(cores_per_node*nodes)+"\n")
            days_   = str(days).zfill(2)
            hours_  = str(hours).zfill(2)
            min_    = str(minutes).zfill(2)
            sec_    = str(seconds).zfill(2)
            text_file.write("#SBATCH --time="+str(days_)+"-"+str(hours_)+":"+str(min_)+":00\n")
            if myhost == 'fidis':
                text_file.write("#SBATCH --constraint=E5v4\n")  # means to only use the fidis nodes
                # to use the Gacrux/Skylake nodes: #SBATCH --constraint=s6g1
            if myhost == 'helvetios':
                print('you can use up to --mem=183G or more')
            if myhost == 'daint':
                text_file.write("#SBATCH --constraint=mc\n")
            text_file.write("#SBATCH --mem=100G\n")


        text_file.write("\n")
        text_file.write("set +e\n")
        text_file.write("# it is necessary to have all the modules which are used when compiling\n")
        text_file.write('export LD_LIBRARY_PATH=""\n')
        if myhost in ["helvetios", "fidis"]:
            if job == 'scaling':
                text_file.write("module load intel intel-mpi intel-mkl fftw gsl eigen\n")
            else:
                text_file.write("module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen\n")
        if myhost == "daint":
            text_file.write("module load daint-mc && module switch PrgEnv-cray PrgEnv-intel && module unload cray-libsci && module load GSL/2.5-CrayIntel-18.08 cray-python/2.7.15.1 cray-fftw")

        text_file.write("module list")
        text_file.write("export LD_LIBRARY_PATH=$HOME/sources/n2p2/lib:${LD_LIBRARY_PATH}\n")
        text_file.write("#echo LD_LIBRARY_PATH: $LD_LIBRARY_PATH\n")
        text_file.write("\n")
        text_file.write("touch time.out\n")
        text_file.write("date +%s >> time.out\n")
        text_file.write("\n")
        if job == 'scaling':
            if interactive == False:
                text_file.write("srun -n "+str(cores)+" $HOME/sources/n2p2/bin/nnp-scaling 1\n")
            else:
                text_file.write("$HOME/sources/n2p2/bin/nnp-scaling 1\n")
        elif job == 'train':
            text_file.write("srun -n "+str(cores)+" $HOME/sources/n2p2/bin/nnp-train\n")
            #text_file.write("strigger --set --jobid=$partA_ID --time --offset=-1200 --program=$dotfiles/scripts/bin/n2p2_get_potential_folder_from_nr.py\n")
            #text_file.write("strigger --set --jobid=$partA_ID --time --offset=-1200 --program=$dotfiles/scripts/bin/n2p2_get_potential_folder_from_nr.py\n")
        text_file.write("date +%s >> time.out\n")
        text_file.write("cat time.out | xargs | awk '{print $2-$1-10}' > time.sec\n")
        text_file.write("$dotfiles/scripts/n2p2/n2p2_tarfolder_for_scale_train.sh\n")
        if job == 'train':
            text_file.write("$dotfiles/scripts/bin/n2p2_get_potential_folder_from_nr.py\n")  # makes the potential folder
        if job == 'scaling':
            text_file.write('[ "`pwd | grep -o "/get_scaling$"`" == \'/get_scaling\' ] && echo creating link && ln -s `pwd`/function.data ../function.data')
        text_file.write("\n")
        text_file.write("exit 0\n")

    print('written (39)',filepath)
    call(["chmod u+x "+filepath],shell=True)
    return filepath

def n2p2_make_training(cores=21,days=7,hours=0,minutes=0,submit_to_que=True,submit_to_debug_que=False):
    if not os.path.isfile("input.data"):
        sys.exit("Need input.data file")
    if not os.path.isfile("input.nn"):
        sys.exit("Need input.nn file")
    n2p2_check_SF_inputnn("input.nn")

    if not os.path.isfile("scaling.data"):
        if not os.path.isfile("get_scaling/scaling.data"):
            sys.exit("get_scaling/scaling.data missing")
        else:
            cp("get_scaling/scaling.data","scaling.data")

    submitskript = n2p2_write_submit_skript(directory=False,cores=cores,nodes=1,days=days,hours=hours,minutes=minutes,job="train")
    submitjob(submit_to_que=True,submit_to_debug_que=False,jobdir=False,submitskript=submitskript)
    create_READMEtxt()
    return

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

def runner_exec(test=False):
    ''' return environment variable runner_exec (RuNNer executable)'''
    runner_exec = os.environ['runner_exec']
    if test == False and not os.path.isfile(runner_exec):
        sys.exit('runner_exec variable is not defined or is not an existing file')
    return runner_exec


if __name__ == "__main__":
    pass
    #n2p2_check_SF_inputnn(inputnn="cursel_64.def")
