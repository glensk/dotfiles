#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import os,sys,re,fnmatch
import filecmp
#import click
import numpy as np
import glob,random #,pathlib
#from my_atom import #atom as my_atom
import my_atom
from itertools import islice

import fah as fah_
import argparse
from copy import deepcopy
from socket import gethostname
import shutil
from tempfile import mkstemp
from subprocess import check_output,call
from datetime import datetime as datetime   # datetime.datetime.now()

import ase
from ase import Atoms
from ase import Atom
from ase.build import bulk as ase_build_bulk
from ase.neighborlist import NeighborList, neighbor_list, NewPrimitiveNeighborList
from ase.spacegroup import crystal
from ase.constraints import StrainFilter, ExpCellFilter, UnitCellFilter
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.optimize import BFGS
from ase.optimize import LBFGS
from ase.optimize import FIRE
from ase.optimize import GPMin
from ase.optimize.basin import BasinHopping
from ase.optimize.minimahopping import MinimaHopping
from ase import units as aseunits
from ase.calculators.singlepoint import SinglePointCalculator
try:
    from ase.calculators.lammpslib import LAMMPSlib
except ImportError:
    print("ERROR when importing LAMMPSlib ... possibly you have to change your (conda/aiida) environment")

from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms



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
#try:
#    from ase.optimize import GPMin
#except ImportError:
#    pass
import time
#import myutils as my

start_time = time.time()

def help(p = None):
    string = '''
    myutils.py -ex abc  # executs function abc
    '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-ef','--execute_function', required=False, type=str,default='', help="function to run")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

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

def read_firstlines(filepath,lines):
    # grep first four lines
    with open(filepath) as myfile:
        head = list(islice(myfile, lines))
    #print('head')
    return head


def printnormal(*var):
    ENDC = '\033[0m'
    return printoutcolor(ENDC,var,ENDC)

def printred(*var):
    ''' print(my.printred("min_at "+str(min_at)+" min_at_orig "+str(min_at_orig)+' (min id '+str(min_at_id)),"this is just an output") '''
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
    #print('keys',keys)
    #print('values',values)
    print('########################################## argparse values (begin) #########')
    for idx,i in enumerate(keys):
        #print('aaaaaaaa',idx)
        #print('kkkkkkkk',list(keys)[idx])
        #print('kkkkkkkk',list(values)[idx])
        #sys.exit()
        print("#",str(list(keys)[idx]).ljust(35),str(list(values)[idx]).ljust(20),str(type(list(values)[idx])).ljust(15),"#")
    print('########################################## argparse values (end  ) #########')
    return

#def sed(file,str_find,str_replace):
#    # from scripts folder
#    import massedit
#    massedit.edit_files([file], ["re.sub('"+str_find+"', '"+str_replace+"', line)"],dry_run=False)
#    return


def sed(source,pattern, replace, dest=None, count=0):
    """Reads a source file and writes the destination file.

    In each line, replaces pattern with replace.

    Args:
        pattern (str): pattern to match (can be re.pattern)
        replace (str): replacement str
        source  (str): input filename
        count (int): number of occurrences to replace
        dest (str):   destination filename, if not given, source will be over written.
    """

    fin = open(source, 'r')
    num_replaced = count

    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = mkstemp()
        fout = open(name, 'w')
    idx=0
    for line in fin:
        #print('line',idx,":"+line+":")
        out = re.sub(pattern, replace, line)
        #print('out ',idx,out)
        #print('pat',pattern)
        #print('rep',replace)
        #print()
        fout.write(out)

        #if idx > 760:
        #    sys.exit('hha445')
        if out != line:
            num_replaced += 1
        if count and num_replaced > count:
            break
    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E

    fin.close()
    fout.close()

    if not dest:
        shutil.move(name, source)
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
                print('idx',idx,path[idx],'path (i)',i)
            if type(i) == bool:
                sys.exit('path to file is not defined! (got a boolean instead)')
            if i is False or i is None and type(pathnames) != bool:
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

def progress_orig(count, total, suffix=''):
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()  # As suggested by Rom Ruben

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

#def scripts():
#    ''' return environment variable scripts '''
#    scripts = os.environ['scripts']
#    if not os.path.isdir(scripts):
#        print('scripts:',scripts)
#        sys.exit('scripts variable is not defined or is not an existing folder')
#    return scripts

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
    host = os.environ["myhost"]
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

def get_atomc_energy_from_dicts(dict_atomic_species_structure,dict_atom_energies):
    ''' dict_atomic_species: {'Ca': 1.0, 'Si': 0.0, 'Mg': 0.0, 'Al': 0.0}
                             {'Ca': 4, 'Si': 0, 'Mg': 0, 'Al': 0}
        dict_atom_energies:  {'Mg': -16.7493351, 'Si': -5.5274864, 'Al': -19.62}
    '''
    out = 0
    for atom_species in dict_atomic_species_structure:
        atom_energy = dict_atom_energies[atom_species]
        en = dict_atomic_species_structure[atom_species]*atom_energy
        out += en
        #print('atom_species',atom_species,'atom_energy',atom_energy,'en',en,'out',out)
    return out

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
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class mypot( object ):
    ''' return a list of available potentials '''
    def __init__(self,potpath_in=False,use_epoch=False,verbose=False):
        self.potpath_in                 = potpath_in
        self.pot                        = False       # n2p2_v1ag
        self.use_epoch                  = use_epoch
        self.potpath                    = False        # this is the source where the potential recides
        self.potpath_work               = False        # this is ususally the self.potpath but for cases where different epoch is used ->
        self.pottype                    = False        # n2p2/runner
        self.pottype_all                = []        # [ 'n2p2','runner','eam',.. ]
        self.potepoch_all               = False
        self.assessed_epochs            = []
        self.assessed_test              = []
        self.assessed_train             = []
        self.assessed_kmc57             = []
        self.assessed_c44               = []
        self.assessed_input             = []
        self.potepoch_bestteste         = False
        self.potepoch_linked            = False
        self.potepoch_using             = False
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
        self.scalingdata                = False         # path to scaling.data
        self.trigger_set_path_manual    = ["setpath","potpath","pp", ".", "..", "../" ]
        self.verbose                    = verbose
        self.elements                   = False  # list e.g. ['Al', 'Mg', 'Si']
        self.reference_volumes          = False
        self.atom_energy                = False
        self.atom_masses                = False # necessary for n2p2 / runner
        self.atom_masses_str            = False # necessary for n2p2 / runner
        self.atom_types                 = None # needs to be defines for runner/n2p2
        self.atom_numbers               = None # needs to be defines for runner/n2p2
                                        # self.atom_types = {'Mg':1,'Al':2,'Si':3}
        self.lmpcmd                     = False
        self.pot_all_dict               = {}
        self.get_pot_all()
        return

    def get_pot_all(self):
        if len(self.pot_all_dict) == 0:
            scripts = os.environ['SCRIPTS']
            home = os.environ['HOME']
            p1 = scripts+'/potentials/'
            p2 = home+'/sources/lammps/potentials/'
            check = [ [p1+'runner*',    'runner'],
                      [p1+'n2p2*',      'n2p2'],
                      [p1+'*eam.alloy*','eam'],
                      [p2+'*eam.alloy*','eam'],
                      [p2+'*.adp',      'adp'],
                    ]
            allpot_fullpath = []
            #allpot_type = []
            for idx,i in enumerate(check):
                allpaths = glob.glob(check[idx][0])
                pottype = check[idx][1]
                if pottype not in self.pottype_all:
                    self.pottype_all += [pottype]
                for jdx,potpath in enumerate(allpaths):
                    pot = os.path.basename(potpath)
                    self.pot_all_dict[pot] = [potpath,pottype,pot]
                    #print('pot',pot)
                    #self.pot_all.append(pot)

            #self.pot_all = self.pot_all + self.trigger_set_path_manual

            #if len(self.pot_all) == 0:
            #    sys.exit('no potentials found; ERROR (92)')
        return


    def get_pottype_elements_and_atomic_energies(self):
        ''' help '''
        # get from input.nn the atomic energies
        if type(self.pot) != bool and type(self.elements) == bool and type(self.atom_energy) == bool:
            #print('self.pot (0):',self.pot,self.pot[:6])
            if self.pot[:6] == 'runner': self.pottype = 'runner'
            if self.pot[:4] == 'n2p2': self.pottype = 'n2p2'
            #print('self.pottype (0)',self.pottype)
            #if 'runner' in self.pot: self.pottype('runner')
            #if 'n2p2' in self.pot: self.pottype('n2p2')
            if self.pottype in self.trigger_set_path_manual or self.pot in self.trigger_set_path_manual:
                self.pottype = inputnn_runner_or_n2p2(self.potpath+"/input.nn")

            if self.pottype == "runner" or self.pottype == "n2p2":
                inputnn = self.potpath+"/input.nn"
                if os.path.isfile(inputnn):
                    self.elements, self.atom_energy = inputnn_get_atomic_symbols_and_atom_energy_dict(inputnn)


                #print('fin (3)',self.atom_types)
                #print('fin (M)',self.atom_masses_str)
                #sys.exit('887')
            elif 'eam' in self.pot and 'alloy' in self.pot: # eam-alloy
                self.pottype = 'eam-alloy'
                # grep first four lines
                with open(self.potpath) as myfile:
                    head = list(islice(myfile, 4))
                #print('head')
                for idx,i in enumerate(head):
                    if idx == 3:
                        self.elements = i.split()[1:]
                    #print(idx,i)
                #print('ele',elements.split()[1:])
                #self.elements = elements.split()[1:]

                self.atom_energy = {}
                for idx,i in enumerate(self.elements):
                    self.atom_energy[i] = 0 #self.elements[idx]
            elif '.adp' in self.pot: # adp potential
                self.pottype = 'adp'
                #print('self.pot',self.pot)
                #print('self.potpath',self.potpath)
                with open(self.potpath) as myfile:
                    head = list(islice(myfile, 4))
                #print('head')
                for idx,i in enumerate(head):
                    #print('idx',idx,i)
                    if idx == 3:
                        self.elements = i.split()[1:]
                    #print(idx,i)
                #print('elemetns_',elements_.split()[1:])
                #sys.exit('error 77 define self.elemetns self.atom_energy')
                self.atom_energy = {}
                for idx,i in enumerate(self.elements):
                    self.atom_energy[i] = 0 #self.elements[idx]
            if self.pottype == False:
                list_pot_all()
                print()
                print('self.pot    :',self.pot)
                print('self.pottype:',self.pottype)
                print('self.potpath:',self.potpath)
                sys.exit("self.pot unknown (Error 91)")
        #print('self.elements',self.elements)
        #print('self.pot',self.pot)
        #print('self.potpath',self.potpath)
        #print('self.elements',self.elements)
        #print('self.atom_energy',self.atom_energy)
        self.atom_numbers = {}
        self.atom_masses = {}
        for i in self.elements:
            self.atom_numbers[i]         = my_atom.atom([i]).number[0]
            self.atom_masses[i]         = my_atom.atom([i]).mass[0]
            #self.reference_volumes[i]   = my_atom.atom([i]).reference_volume[0]
        #print('fin (M)',self.atom_masses)
        list_atom_nr = []
        for i  in self.atom_numbers:
            list_atom_nr.append(self.atom_numbers[i])
            #print(i,self.atom_numbers[i])
        #print('fin (2):',np.sort(np.array(list_atom_nr)))
        self.atom_types = {}
        self.atom_masses_str = []
        for idx,i in enumerate(np.sort(np.array(list_atom_nr))):
            for j in self.atom_numbers:
                #print(idx,i,"||",j,self.atom_numbers[j])
                if i == self.atom_numbers[j]:
                    self.atom_types[j] = idx+1
                    self.atom_masses_str = self.atom_masses_str+ [ "mass "+str(idx+1)+" "+str(self.atom_masses[j])]
                    #print('-->',j,self.atom_types[j],idx+1)
                    #print('-->xx',self.atom_types)
                    break

        self.reference_volumes = {}
        self.atomic_numbers = {}
        for i in self.elements:
            #self.atom_masses[i]         = my_atom.atom([i]).mass[0]
            #print('iii',i)
            #print()
            #print('all',my_atom.atom([i]))
            #print()
            #print('mass',my_atom.atom([i]).mass[0])
            #print()
            #print('rv',(my_atom.atom([i]).reference_volume)[0])
            rv  = (my_atom.atom([i]).reference_volume)[0]
            an  = (my_atom.atom([i]).number)[0]
            self.reference_volumes[i] = rv
            self.atomic_numbers[i] = an
            #print('rv',my_atom.atom([i]).reference_volume[0])
        #print('rv',self.reference_volumes)
        #sys.exit('123')
        return

    def print_variables_mypot(self,text="",print_nontheless=False,exit=False):
        if self.verbose > 1 or print_nontheless:
            print(text,"self.potpath_in             ",self.potpath_in)
            print(text,"self.pot                    ",self.pot)
            print(text,"self.potpath (source = src) ",self.potpath)
            print(text,"self.pottype                ",self.pottype)
            print(text,"self.potpath_work (typ src) ",self.potpath_work)
            print(text,"self.inputnn                ",self.inputnn)
            print(text,"self.inputdata              ",self.inputdata)
            print(text,"self.learning_curve_file    ",self.learning_curve_file)
            print(text,"self.potepoch_all           ",self.potepoch_all)
            print(text,"self.assessed_epochs        ",self.assessed_epochs)
            print(text,"self.assessed_test          ",self.assessed_test)
            print(text,"self.assessed_train         ",self.assessed_train)
            print(text,"self.assessed_kmc57         ",self.assessed_kmc57)
            print(text,"self.assessed_c44           ",self.assessed_c44)
            print(text,"self.use_epoch (from input) ",self.use_epoch)
            print(text,"self.potepoch_bestteste     ",self.potepoch_bestteste)
            print(text,"self.potepoch_linked        ",self.potepoch_linked)
            print(text,"self.potepoch_using (finaly)",self.potepoch_using)
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
            print(text,"self.reference_volumes      ",self.reference_volumes)
            print(text,"self.atom_types             ",self.atom_types)
            print(text,"self.atom_numbers           ",self.atom_numbers)
            print(text,"self.atom_masses            ",self.atom_masses)
            print(text,"self.atom_masses_str        ",self.atom_masses_str)
            print(text,"self.potDONE                ",self.potDONE)
            print(text,"self.verbose                ",self.verbose)
            print(text,"len(self.pot_all_dict)      ",len(self.pot_all_dict))
        if exit == True:
            sys.exit('778899')
        return

    def get_my_assessments(self):
        self.test_b = self.test_l = self.train_b = self.train_l = self.kmc57_b = self.kmc57_l = False
        if len(self.potepoch_all) == 0:
            return
        #print('folder',self.potpath)
        #print('all',self.potepoch_all)
        epoch_last = self.potepoch_all[-1]
        #print('epoch_last',epoch_last)
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
                idx = i.replace(self.potpath+"/assess_"+ext+"_","")
                #print('i99',i,idx)
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
                    if os.path.isfile(file): self.assessed_test[eidx] = float(read_lastline(file))
                    else: self.assessed_test[eidx] = "-"
                    if epoch == self.potepoch_bestteste: self.assessed_test_b = self.assessed_test[eidx]
                    if epoch == self.assessed_epochs[-1]: self.assessed_test_l = self.assessed_test[eidx]
                if ext == 'train':
                    if os.path.isfile(file): self.assessed_train[eidx] = float(read_lastline(file))
                    else: self.assessed_train[eidx] = "-"
                    if epoch == self.potepoch_bestteste: self.assessed_train_b = self.assessed_train[eidx]
                    if epoch == self.assessed_epochs[-1]: self.assessed_train_l = self.assessed_train[eidx]
                if ext == 'kmc57':
                    if os.path.isfile(file): self.assessed_kmc57[eidx] = float(read_lastline(file))
                    else: self.assessed_kmc57[eidx] = "-"
                    if epoch == self.potepoch_bestteste: self.assessed_kmc57_b = self.assessed_kmc57[eidx]
                    if epoch == self.assessed_epochs[-1]: self.assessed_kmc57_l = self.assessed_kmc57[eidx]
                if ext == 'input':
                    if os.path.isfile(file): self.assessed_input[eidx] = float(read_lastline(file))
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
        #self.potepoch_bestteste_checked = False
        self.print_variables_mypot('PPk get potential: in')

        #############################################
        # get potential from path  1. get potpath_in
        #############################################
        if self.verbose > 2:
            print('PPk in0000')
            print('PPk self.potpath',self.potpath)
            print('PPk self.potpath_in',self.potpath)

        ##### get first list of know potentials
        if len(self.pot_all_dict) == 0:
            self.get_pot_all()

        try:
            # here a proper name of a potential was specified
            self.potpath = self.pot_all_dict[self.potpath_in][0]
            self.pottype = self.pot_all_dict[self.potpath_in][1]
            self.pot     = self.pot_all_dict[self.potpath_in][2]
        except:
            # here only a relative path or something like this has been specified as ../ or . or ./
            #print('seems like a relative path ....')
            self.potpath = os.path.abspath(self.potpath_in)
            self.pot     = os.path.basename(self.potpath_in)
            #print('self.potpath',self.potpath)
            #print('self.pot    ',self.pot)
        self.get_pottype_elements_and_atomic_energies()  # defines self.pottype if not defined


        if self.potpath == False:
            self.print_variables_mypot(text="self.potpath not found! (1)",print_nontheless=False,exit=True)

        if self.verbose:
            print('PPk self.potpath_in',self.potpath_in)
            print('PPk self.potpath   ',self.potpath)
            print('PPk self.pottype   ',self.pottype)
            print('PPk self.pot       ',self.pot)

        if self.pottype in [ "runner", "n2p2" ]:
            checkfiles = [ "input.nn", "scaling.data"]
            # add to checkfiles: "weights.012.data", "weights.013.data", "weights.014.data"
            if False:
                for i in self.atomic_numbers:
                    filename = "weights."+(str(self.atomic_numbers[i]).zfill(3))+".data"
                    checkfiles += [filename]

            for i in checkfiles:
                if not os.path.isfile(self.potpath+"/"+i):
                    print('checked for:',self.potpath+"/"+i)
                    sys.exit(self.potpath+"/"+i+" does not exist! (3)")

            self.inputnn     = isfiledir(self.potpath+"/input.nn",exit=True)
            self.scalingdata = isfiledir(self.potpath+"/scaling.data",exit=True)
            self.inputdata   = isfiledir(self.potpath+"/input.data",check_extension=[".extxyz"])
            self.testdata    = isfiledir(self.potpath+"/test.data",check_extension=[".extxyz"])
            self.traindata   = isfiledir(self.potpath+"/train.data",check_extension=[".extxyz"])
            #self.pottype    = inputnn_runner_or_n2p2(self.inputnn)
            self.learning_curve_file = n2p2_runner_get_learning_curve_filename(self.inputnn)
            try:
                self.potepoch_all = inputnn_get_potential_number_from_all_weightsfiles(self.inputnn)
            except NameError:
                self.potepoch_all = []
            self.potepoch_bestteste = n2p2_runner_get_bestteste_idx(self.inputnn)

            self.potepoch_linked = False
            weightfiles = glob.glob(self.potpath+'/weights.*.data')
            for i in weightfiles:
                if os.path.islink(i):
                    a = os.path.basename(os.path.realpath(i))
                    self.potepoch_linked = int(a[12:18])


            ##### check if current weights.xxx.data files are the ones from weights.xxx. self.potepoch_bestteste
            #if False: # dont do this test, accept any potential that is linke, even if this one is not the best (forced c44 agreement)
            #    if exit == True:
            #        file1 = self.potpath+"/weights.012.data"
            #        epstr = str(self.potepoch_bestteste).zfill(6)
            #        file2 = self.potpath+"/weights.012."+epstr+".out"
            #        print('file1',file1)
            #        print('file2',file2)
            #        print('epstr',epstr)
            #        filecompare = filecmp.cmp(file1, file2)
            #        print('filecompare',filecompare)
            #        sys.exit()
            #        if not filecmp.cmp(file1, file2):
            #            sys.exit("PPd File "+file1+" is not "+file2)
            #        else:
            #            self.potepoch_bestteste_checked = True


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


            ###################################
            # check if we use the right epoch
            ###################################
            # if no epoch to use specified, use whatever is available
            #print('useepo',self.use_epoch,type(self.use_epoch))
            #print('linked',self.potepoch_linked,type(self.potepoch_linked))
            if self.potepoch_linked == False and self.use_epoch == False:
                self.potpath_work = self.use_epoch = self.potepoch_all[-1]
                self.potpath_work = self.pot_tmpdir
            print('self.use_epoch',self.use_epoch)
            print('self.potpath_work',self.potpath_work)
            print('self.potepoch_linked',self.potepoch_linked)
            print('self.potepoch_using',self.potepoch_using)
            print('self.potpath',self.potpath)
            print('self.potepoch_linked',self.potepoch_linked)
            if self.use_epoch == False:
                self.potpath_work   = self.potpath
                self.potepoch_using  = self.potepoch_linked
            elif self.use_epoch != False and self.potepoch_linked == self.use_epoch:
                # some epoch is specified but it is already the linked one!
                self.potpath_work   = self.potpath
                self.potepoch_using  = self.potepoch_linked
            else:
                # some epoch is specified but it is not the linked one
                self.potpath_work = self.pot_tmpdir
                print('We need to work in a different folder since other epoch is linked!')
                epstr = str(self.use_epoch).zfill(6)
                print("PPf copy files ...",epstr)
                if not os.path.isdir(self.pot_tmpdir):
                    mkdir(self.pot_tmpdir)
                f  = []
                fa = []
                fb = []
                wfn = []
                for i in self.atomic_numbers:
                    #filename = "weights."+(str(self.atomic_numbers[i]).zfill(3))+".data"
                    atnrstr_ = str(self.atomic_numbers[i]).zfill(3)
                    weightsfile = "weights."+atnrstr_+"."+epstr+".out"
                    wfn += ["weights."+atnrstr_+".data"]
                    f_  = self.potpath+"/"+weightsfile
                    fa_ = self.potpath+"/../_weights/"+weightsfile
                    fb_ = self.potpath+"/../"+weightsfile

                    print('f_',f_)
                    if not os.path.isfile(f_): f_ = fa_
                    if not os.path.isfile(f_): f_ = fb_
                    if not os.path.isfile(f_): sys.exit(f_+" does not exist! (65)")

                    f  += [f_]
                    fa += [fa_]
                    fb += [fb_]

                #f12 = self.potpath+"/weights.012."+epstr+".out"
                #f13 = self.potpath+"/weights.013."+epstr+".out"
                #f14 = self.potpath+"/weights.014."+epstr+".out"

                #f12a = self.potpath+"/../_weights/weights.012."+epstr+".out"
                #f13a = self.potpath+"/../_weights/weights.013."+epstr+".out"
                #f14a = self.potpath+"/../_weights/weights.014."+epstr+".out"

                #f12b = self.potpath+"/../weights.012."+epstr+".out"
                #f13b = self.potpath+"/../weights.013."+epstr+".out"
                #f14b = self.potpath+"/../weights.014."+epstr+".out"
                #if not os.path.isfile(f[0]): f[0] = fa[0]
                #if not os.path.isfile(f[1]): f[1] = fa[1]
                #if not os.path.isfile(f[2]): f[2] = fa[2]
                #if not os.path.isfile(f[0]): f[0] = fb[0]
                #if not os.path.isfile(f[1]): f[1] = fb[1]
                #if not os.path.isfile(f[2]): f[2] = fb[2]
                #if not os.path.isfile(f[0]): sys.exit(f[0]+" does not exist! (65)")
                #if not os.path.isfile(f[1]): sys.exit(f[1]+" does not exist! (66)")
                #if not os.path.isfile(f[2]): sys.exit(f[2]+" does not exist! (67)")


                #for ff in [f12,f13,f14]:
                #    if not os.path.isfile(ff):
                #        sys.exit(ff+" does not exist! (65)")
                for idx,i in enumerate(f):
                    if self.verbose:
                        print('copying',i,'to',self.pot_tmpdir+'/'+wfn[idx])
                    cp(i,self.pot_tmpdir+"/"+wfn[idx]) #weights.012.data")
                #print('copying f13',f[1],'to',self.pot_tmpdir)
                #print('copying f14',f[2],'to',self.pot_tmpdir)
                #print('copying input.nn',self.inputnn)
                #print('copying scalingdata',self.scalingdata)
                #cp(f12,self.pot_tmpdir+"/weights.012.data")
                #cp(f13,self.pot_tmpdir+"/weights.013.data")
                #cp(f14,self.pot_tmpdir+"/weights.014.data")
                print("self.pot_tmpdir",self.pot_tmpdir)
                cp(self.scalingdata,self.pot_tmpdir)
                cp(self.inputnn,self.pot_tmpdir)
                if self.verbose:
                    print('... copying potential',epstr,'to',self.pot_tmpdir)
                #print('self.potpath_work (2)',self.potpath_work)
        self.potDONE = True
        #print('self.potpath_work (3)',self.potpath_work)
        self.print_variables_mypot('get potential: out')
        return

def pot_all():
    mp = mypot()
    mp.get_pot_all()
    return mp.pot_all
    #return mp.pot_all_dict_potpath
    ##return mp.pot_all

def list_pot_all():
    mp = mypot()
    mp.get_pot_all()
    for i in mp.pottype_all:
        print('---------'+"-"*len(i))
        print('pottype:',i)
        print('---------'+"-"*len(i))
        my_dict =  mp.pot_all_dict
        for item in my_dict:
            if my_dict[item][1] == i:
                print(" {} , potpath: {}".format(item,my_dict[item][0]))
    return

##################################################################################
## ipi related functions
##################################################################################


def create_ipi_kmc_inputfile(jobdir,filename="input-runner.xml", nsteps=False,stride=100,seed=12345,a0=4.057,ncell=4,nsi=3,nmg=3,nvac=1,neval=8,temp=300,verbosity="low",checkpoint_restart_stride=10,nodes=1,address="hostname_jobdir",testrun=False,cubic=False):
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
    if cubic == True:
        activelist = str(range(4*ncell**3 - nvac))
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

def submitskript_1():
    text1 = ["#!/bin/bash",""]
    return text1

def submitskript_2(nodes,ntasks,submittime_hours):
    text2 = [
    "#SBATCH --job-name=NNP-mpi",
    "#SBATCH --get-user-env",
    "#SBATCH --output=_scheduler-stdout.txt",
    "#SBATCH --error=_scheduler-stderr.txt",
    "#SBATCH --nodes="+str(nodes),
    "#SBATCH --ntasks "+str(ntasks),
    "#SBATCH --time=00-"+str(submittime_hours)+":00:00"]
    return text2

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

    text1 = submitskript_1()
    test2 = submitskript_2(nodes,ntasks,submittime_hours)

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

def cell_is_cubic(cell):
    if cell[0,1] == cell[0,2] == cell[1,0] == cell[1,2] == cell[2,0] == cell[2,1] == 0:
        if cell[0,0] == cell[1,1] == cell[2,2]:
            return True
    return False

def create_READMEtxt(directory=False,add=False):
    ''' wiretes a README.txt file '''
    if directory == False:
        directory = os.getcwd()
    # get sha
    pwd = os.getcwd()
    os.chdir(os.environ['SCRIPTS'])
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
    return filepath

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
            command = ["sbatch","-p","debug","-t","01:00:00",submitskript]
        else:
            command = ["sbatch",                             submitskript]
        print('command',command)
        call(command)

    os.chdir(hier)
    return


def get_number_of_atoms_as_function_of_cutoff():
    '''
    %myutils.py
    making cell
    starting loop
    2 0
    3 12
    4 18
    5 42
    6 54
    7 86
    8 140
    9 200
    10 320  # 10.583544 Angstrom (==20 bohrradius) is current max cutoff of NN potential
    10.6 320  # 10.583544 Angstrom (==20 bohrradius) is current max cutoff of NN potential
    --------------------------------------------------------------------------------------
    10.8 368  # 10.583544 Angstrom (==20 bohrradius) is current max cutoff of NN potential
    11 368  # 10.583544 Angstrom (==20 bohrradius) is current max cutoff of NN potential
    12 458
    13 602
    14 766

    --> maxcutoff is 10.6 Angstrom which needs to be considered
    '''
    a0 = 4
    matrix_element = "Al"
    ncell = 2
    vacidx = 0

    print('making cell')
    atoms = ase_build_bulk(matrix_element,crystalstructure='fcc',a=a0,cubic=True)
    atoms = atoms.repeat(ncell)
    atoms[vacidx].symbol = "V"
    atomsc_fakevac = atoms.copy()
    print('starting loop')
    for cutoff in [10.8]: #np.arange(2,15):
        NN = ase_get_neighborlist(atomsc_fakevac,atomnr=vacidx,cutoff=cutoff,skin=0.1)
        print(cutoff,len(NN))
    return

def get_NN1_from_positions(atoms,alat,atomnr=0,verbose=False):
    NN1         = ase_get_neighborlist(atoms,atomnr=atomnr,cutoff=0.8*alat,skin=0.1)
    return NN1

def get_NN1_idx_for_every_atom(atoms,crystal_structure='fcc'):
    ''' crystal_strucure = fcc,bcc,hcp ...
    '''
    if crystal_structure == 'fcc':
        amount_first_NN        = 12
        amount_first_NN_unuque = 6
    elif crystal_structure == 'bcc':
        amount_first_NN = 8
        amount_first_NN_unuque = 4
    else:
        print('crystal_structure',crystal_structure)
        sys.exit('so far only for fcc and bcc')

    sc,alat = try_to_get_alat_and_sc_from_cell_and_positions(atoms.cell,atoms.positions)
    #print('sc',sc,'alat',alat)

    ## get first NN for every atom
    NN1_       = np.zeros((atoms.get_number_of_atoms(),amount_first_NN))
    #NN1_unique = np.zeros((atoms.get_number_of_atoms(),amount_first_NN_unuque))
    for i in np.arange(atoms.get_number_of_atoms()):
        lookat = i
        ka = get_NN1_from_positions(atoms,alat*1.1,atomnr=lookat,verbose=False)
        #print('ka',ka,'len',len(ka))
        NN1_[i] = get_NN1_from_positions(atoms,alat*1.1,atomnr=lookat,verbose=False)
    NN1_ = NN1_.astype(int)
    #print('NN1_')
    #print(NN1_)
    #get_NN1_idx_unique_for_every_atom(NN1_)
    return NN1_

def get_NN1_idx_unique_for_every_atom(NN1):
    print("NN1shape--->",NN1.shape)
    nat  = NN1.shape[0]
    neigh = NN1.shape[1]
    if ((neigh % 2) == 0) == False:
        print("NN1shape--->",NN1.shape)
        sys.exit('odd number of neighbors, exi')
    print('nat  ',nat  )
    print('neigh',neigh)
    print('neigh',neigh/2)

    def make_NN1uniqe_one_free_spot(NN1uniq,verbose=True):
        # if line is full
        #   check if corresponding atoms have fewer lines and if so swich
        possible_all = np.zeros(len(NN1uniq))
        for idx,i in enumerate(NN1uniq):
            possible = np.where(NN1uniq[idx]==-1)[0]
            possible_all[idx] = int(len(possible))

        for idx,idss in enumerate(NN1uniq):
            print(idx,NN1uniq[idx],'possible',possible_all[idx],type(possible_all[idx]))

        # check if all equal length:
        if np.all(possible_all == possible_all[0]):
            print('all equal')
            return NN1uniq

        # if different lenght (if diff == 1, skip it)
        #print('poss',possible_all)
        #print('diff',np.diff(possible_all))
        #diff = np.abs(np.diff(possible_all)).max()
        minlength = possible_all.min()
        maxlength = possible_all.max()
        if np.abs(minlength-maxlength) <= 2:
            print('all close to equal (diff 2)')
            return NN1uniq

        at1_weg_all = np.where(possible_all == minlength)[0]
        print('at1_weg_all',at1_weg_all)
        for at1_weg in at1_weg_all:
            at2_check_all = NN1uniq[at1_weg]
            print('at1_weg (a)',at1_weg,'at2_check_all',at2_check_all)
            at2_check_all = at2_check_all[np.where(at2_check_all >= 0)]
            at2_check_possibilities = possible_all[at2_check_all]
            at2_take_check = np.where(at2_check_possibilities == at2_check_possibilities.max())[0][0]
            print('at1_weg (b)',at1_weg,'at2_check_all',at2_check_all,at2_check_possibilities,'take',at2_take_check)
            NN1_unique = add_to_line_which_atom(NN1uniq,line_nr=at2_take_check,atom=at1_weg)
            #for at2_check in at2_check_all:
            #    at = np.where(possible_all == minlength)[0]

        sys.exit('not close')
        for idx,i in enumerate(NN1uniq):
            if len(possible) == 0: # line is full
                otheratoms = NN1uniq[idx] # == i
                if verbose:
                    print('#############################')
                    for ddd,idss in enumerate(NN1uniq):
                        print(ddd,NN1uniq[ddd])
                    print('#############################')
                    print('full',idx,'otheratoms',otheratoms)
                    sys.exit()
                shortest = 0
                for o in otheratoms:
                    possible_o = np.where(NN1uniq[o]==-1)[0]
                    print('posible_o',possible_o,len(possible_o))
                    #if len(possible_o) > shortest:
                    #    shortest =
        for ddd,idss in enumerate(NN1uniq):
            print(ddd,NN1uniq[ddd])
        sys.exit()
        return NN1uniq


    def add_to_line_which_atom(NN1uniq,line_nr,atom):
        ''' line is the line/atom where to add
            atom is the corresponding atom to add in this line
        '''
        for ddd,idss in enumerate(NN1uniq):
            print('before',ddd,NN1uniq[ddd])
        line = NN1uniq[line_nr]
        print('line_nr',line_nr,'line(a)',line)
        number_in_line = np.where(line==-1)[0]
        print('number_in_line',number_in_line)
        NN1uniq[line,nnext] = atom
        line = NN1uniq[line_nr]
        print('line(b)',line)
        for ddd,idss in enumerate(NN1uniq):
            print('after',ddd,NN1uniq[ddd])
        sys.exit()
        return NN1uniq

    def get_atom_with_most_empty_slots(at1,at2_all,NN1uniq):
        ''' this needs to be done before it is decided which atoms are used '''
        print('at1',at1)
        print('at2_all',at2_all)
        at_all = np.concatenate((np.array([at1]),at2_all))
        print('at_all',at_all)
        line_possible = np.array([-1,-1])
        for idx,i in enumerate(at_all):
            print('idx',idx,'i',i)
            print('NN1uniq[i]',NN1uniq[i])
            possible = np.where(NN1uniq[i]==-1)[0]
            if len(possible) >= line_possible[1]:
                line_possible = np.array([i,len(possible)])
        atnext = line_possible[0]
        nnext  = np.where(NN1uniq[atnext]==-1)[0]
        if line_possible[1] == 0 or len(nnext) == 0:
            print('line_possible:',line_possible,"(if 0 -> problem)")
            for ddd,idss in enumerate(NN1uniq):
                print(ddd,NN1uniq[ddd])
            sys.exit('neither nnextat1 nor nnextat2')
        # here, atnext can either be from at1 or at2_all
        return atnext,nnext[0]

    NN1uniq = (np.ones((nat,(neigh/2)))*-1).astype(int)
    liste_in = []
    for at1 in np.arange(nat):
        at2idx = 0
        at2_all = NN1[at1]
        # -----> at1 (a) 12 at2_all [ 4  5  6  7 16 18 20 22 24 25 28 29]
        print('-----> at1 (a)',at1,'at2_all',at2_all)
        for at2idx,at2 in enumerate(at2_all):
            at2_all_rest = at2_all[at2idx:]
            ww = [at1,at2]
            ww2 = [at2,at1]
            if ww in liste_in or ww2 in liste_in:
                print('have already',ww,ww2)
            else:

                NN1uniq = make_NN1uniqe_one_free_spot(NN1uniq,verbose=True)
                atnext,nnext = get_atom_with_most_empty_slots(at1,at2_all_rest,NN1uniq)
                print('need to append (a)',at1,at2_all,'-> atnext',atnext,nnext)
                if at1 == 20:
                    sys.exit('kk')

                if atnext == at1:
                    ww = [at1,at2]
                    ww2 = [at2,at1]
                    atom_ = at2
                else:  # one of the at2_all
                    ww = [atnext,at1]
                    ww2= [at1,atnext]
                    atom_ = at1

                if ww in liste_in or ww2 in liste_in:
                    print('have already',ww,ww2)
                else:
                #if ww not in liste_in and ww2 not in liste_in:
                    liste_in.append(ww)
                    liste_in.append(ww2)
                    print('appending      (b)','-> atnext',atnext,nnext)
                    NN1uniq[atnext,nnext] = atom_
                #nnextat1 = np.where(NN1uniq[at1]==-1)[0]
                #nnextat2 = np.where(NN1uniq[at2]==-1)[0]
                #print('need to append',ww,ww2,nnextat1,nnextat2,'-> atnext',atnext)
                #if len(nnextat2) >= len(nnextat1):
                #    if len(nnextat2) == len(nnextat1) == 0:
                #        for ddd,idss in enumerate(NN1uniq):
                #            print(ddd,NN1uniq[ddd])
                #        sys.exit('neither nnextat1 nor nnextat2')

                #    NN1uniq[at2,nnextat2[0]] = at1
                #else:
                #    NN1uniq[at1,nnextat1[0]] = at2

                #if len(nnextat1) > 0:
                #    NN1uniq[at1,nnextat1[0]] = at2
                #else:
                #    if len(nnextat2) > 0:
                #        NN1uniq[at2,nnextat2[0]] = at1
                #    else:
                #        for ddd,idss in enumerate(NN1uniq):
                #            print(ddd,NN1uniq[ddd])
                #        sys.exit('neither nnextat1 nor nnextat2')
                #if at2idx < neigh/2.:
                #    nnextat1 = np.where(NN1uniq[at1]==-1)[0]
                #    NN1uniq[at1,nnextat1[0]] = at2
                #elif at2idx >= neigh/2.:
                #    #print('nnext1',NN1uniq[at2])
                #    #print('nnext2',NN1uniq[at2]<0)
                #    nnextat2 = np.where(NN1uniq[at2]==-1)[0]
                #    NN1uniq[at2,nnextat2[0]] = at1
            at2idx += 1
        print('-----> at1 (b)',at1,'at2_all',at2_all)
        for ddd,idss in enumerate(NN1uniq):
            print(ddd,NN1uniq[ddd])
        #if at1 == 3:
        #    sys.exit()

    print('NN1uniq')
    print(NN1uniq)
    sys.exit()
    return


def get_NN1_NN2_NN3_from_positions(atoms,alat,atomnr=0,verbose=False):
    NN1         = ase_get_neighborlist(atoms,atomnr=atomnr,cutoff=0.8*alat,skin=0.1)
    NN1_NN2     = ase_get_neighborlist(atoms,atomnr=atomnr,cutoff=1.1*alat,skin=0.1)
    NN1_NN2_NN3 = ase_get_neighborlist(atoms,atomnr=atomnr,cutoff=1.23*alat,skin=0.1)
    NN2         = np.sort(np.array(diff(NN1_NN2,NN1)))
    NN2_NN3     = np.sort(np.array(diff(NN1_NN2_NN3,NN1)))
    NN3         = np.sort(np.array(diff(NN2_NN3,NN2)))
    if verbose:
        print("NN1",NN1,'len',len(NN1))
        #print("NN2_NN1    :",NN1_NN2)
        #print("NN2_NN1_NN3:",NN1_NN2_NN3)
        #print("NN2_NN3:",NN2_NN3)
        print("NN2",NN2,'len',len(NN2))
        print("NN3",NN3,'len',len(NN3))
    return NN1,NN2,NN3

def get_forces_on_atom_by_considering_all_its_1NNs(params,forcefunc,atoms,alat,atomnr,pos):
    ''' alat  : 4.13 (angstrom)
        atoms : an ase object of the perfect cell.
    '''
    NN1_1st_neigh = get_NN1_from_positions(atoms,alat,atomnr=atomnr,verbose=False)
    #print('NN1_1st_neigh',NN1_1st_neigh)
    F = 0
    for idxnn11,nn11 in enumerate(NN1_1st_neigh):    # for all neighbors of a particular 1NN
        D  = vec = pos[atomnr]-pos[nn11]
        vecnorm = r = np.linalg.norm(vec)
        fnorm = forcefunc(vecnorm,*params)
        forces = vec/vecnorm*fnorm
        #print('a1',atomnr,'a2',nn11,'r',str(np.round(r,3)).ljust(5),'fnorm',str(np.round(fnorm,3)).ljust(3),'p1',pos[atomnr],'p2',pos[nn11],"D",D,'forces',forces)
        F += forces
    #print('FORCES',F)
    return F

def create_al_sphere(a0=4.05,matrix_element="Al",cubic=True,ncell=4,nvac=1,cutoff=4.05,vacidx=0):
    ''' cutoff = a0 # a0 includes up to 2NN '''
    nndist = a0/np.sqrt(2.)
    atoms = ase_build_bulk(matrix_element,crystalstructure='fcc',a=a0,cubic=cubic)
    atoms = atoms.repeat(ncell)
    atoms[vacidx].symbol = "V"
    atomsc_fakevac = atoms.copy()

    #atomsc_fakevac = get_ase_atoms_object_kmc_al_si_mg_vac(ncell=ncell,nsi=0,nmg=0,matrix_element=matrix_element,nvac=nvac,a0=a0,cubic=cubic,create_fake_vacancy = True,normal_ordering=True)
    #NN_1_indices, NN_2_indices = ase_get_neighborlist_1NN_2NN(atomsc_fakevac,atomnr=0,cutoffa=nndist,cutoffb=a0,skin=0.1)
    NN = ase_get_neighborlist(atomsc_fakevac,atomnr=vacidx,cutoff=cutoff,skin=0.1)
    #print('NN_1_indices (orig  ):',NN_1_indices)
    #print('NN_2_indices (orig  ):',NN_2_indices)
    #print('NN_1_2_indices (orig  ):',NN_1_2_indices_tmp,type(NN_1_2_indices_tmp))
    atomsc_sphere = ase_get_atoms(atomsc_fakevac,np.append(NN,[vacidx]))
    return atomsc_sphere

def analyze_accuracy_of_filling_cell_with_Al(ace,ncell=3,nsi=3,nmg=3,nvac=1,a0=4.05,create_fake_vacancy=True):
    atoms = get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0,matrix_element="Al",cell='cubic',create_fake_vacancy=create_fake_vacancy,crystal_structure="fcc",normal_ordering=True)
    ele_list = atoms.get_chemical_symbols()
    print(ele_list,type(ele_list))
    import random
    random.shuffle(ele_list)
    print(ele_list,type(ele_list))
    V_idx = [i for i, x in enumerate(ele_list) if x == "V"]
    Si_idx = [i for i, x in enumerate(ele_list) if x == "Si"]
    Mg_idx = [i for i, x in enumerate(ele_list) if x == "Mg"]
    Al_idx = [i for i, x in enumerate(ele_list) if x == "Al"]
    print('ele_list.index("V")',V_idx)
    print('ele_list.index("Si")',Si_idx)
    ele_list[V_idx[0]], ele_list[Al_idx[0]] = ele_list[Al_idx[0]], ele_list[V_idx[0]]
    print(ele_list,type(ele_list))
    atoms.set_chemical_symbols(ele_list)
    if False:
        for idx,i in enumerate(atoms):
            print('idx',idx,atoms.get_chemical_symbols()[idx],atoms.positions[idx])

    def get_ene_geop(atoms_tmp):
        sys.exit('this was only done for mg, al, si and is not genereal')
        atom_types = {'Mg':1,'Al':2,'Si':3}
        #atom_types = {'Mg':1,'Si':3,'Al':2}  # same same
        lmpcmd = ['mass 1 24.305', 'mass 2 26.9815385', 'mass 3 28.0855', 'variable nnpDir string "/Users/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/runner_v3ag_5000_46489_2"', 'pair_style runner dir ${nnpDir} showewsum 1 showew yes resetew no maxew 1000000', 'pair_coeff * * 14.937658735']
        from ase.calculators.lammpslib import LAMMPSlib
        asecalcLAMMPS = LAMMPSlib(lmpcmds=lmpcmd, atom_types=atom_types,keep_alive=True)
        atoms_tmp.set_calculator(asecalcLAMMPS)
        #ene = atoms_tmp.get_potential_energy()
        #print('ene sta2',ene*0.036749322,"Hartree")

        from ase.optimize import GPMin
        opt = GPMin(atoms_tmp,logfile=None)
        fmax=0.03
        opt.run(fmax=fmax)
        ene = atoms_tmp.get_potential_energy()
        #    ace.geopt = True
        #    ace.pot_get_and_ase_lmp_cmd()
        #    ene = ace.ene(atoms_tmp,minimizer="GPMin")
        #print('ene geop',ene*0.036749322,"Hartree")
        return ene*0.036749322

    def get_ene_geop_from_filename(filename):
        normal = ase_read(filename)
        atoms_tmp = normal.copy()  # for other instances, since frames change when geoopt
        ace.pot_get_and_ase_lmp_cmd()
        ase_remove_unknown_species_from_structure(atoms_tmp,ace)
        ene = get_ene_geop(atoms_tmp)
        return ene


    ace.pot.print_variables_mypot(print_nontheless=True,text="myutils 2304")

    #for step in np.arange(10):
    for step in np.arange(1):
        ref_n = glob.glob('pos_step_'+str(step)+'_normal.extxyz')
        #ref_f = glob.glob('pos_step_'+str(step)+'_filled.extxyz')
        if ref_n == []:
            continue
        print(step,ref_n)
        swap = glob.glob('pos_step_'+str(step)+'_swap*_normal.extxyz')
        if swap == []:
            continue
        print(swap)
        sneigh = []
        for i in swap:
            print(i.split("_"))
            svac = int(i.split("_")[4])
            sneigh.append(int(i.split("_")[5]))
        #print('sneigh',sneigh)
        sneigh = np.sort(np.array(sneigh))
        print('sneigh',sneigh)


        ene_ref_n = get_ene_geop_from_filename(ref_n[0])
        #ene_ref_f = get_ene_geop_from_filename(ref_f[0])
        print('ene_ref_n geop ref ',ene_ref_n,"Hartree")
        #print('ene_ref_f geop ref ',ene_ref_f,"Hartree")
        hartree_to_mev = 27211.386
        diff_mev_n = np.zeros((12,2))
        diff_mev_f = np.zeros((12,2))
        for idx,i in enumerate(sneigh):
            #print('idx',idx,'i',i,len(sneigh))
            filename = 'pos_step_'+str(step)+'_swap_'+str(svac)+"_"+str(i)+'_normal.extxyz'
            ene_n = get_ene_geop_from_filename(filename)
            print('idx',idx,'svac',svac,'sneigh',i,'diff normal ref',ene_ref_n,'fin:',ene_n,'diff:',(ene_n-ene_ref_n)*hartree_to_mev/32)

            filename = 'pos_step_'+str(step)+'_filled_'+str(svac)+"_"+str(i)+'_normal.extxyz'
            ene_f_ref = get_ene_geop_from_filename(filename)
            filename = 'pos_step_'+str(step)+'_swap_'+str(svac)+"_"+str(i)+'_filled.extxyz'
            ene_f = get_ene_geop_from_filename(filename)
            #print('diff filled',(ene_f-ene_f_ref)*hartree_to_mev/32)
            print('idx',idx,'svac',svac,'sneigh',i,'diff filled ref',ene_f_ref,'fin:',ene_f,'diff:',(ene_f-ene_f_ref)*hartree_to_mev/32)

            diff_mev_n[idx][0] = step*12+idx
            diff_mev_n[idx][1] = (ene_n-ene_ref_n)*hartree_to_mev/32
            diff_mev_f[idx][0] = step*12+idx
            diff_mev_f[idx][1] = (ene_f-ene_f_ref)*hartree_to_mev/32

        np.savetxt("nrates_step_"+str(step)+"_filled",diff_mev_f)
        np.savetxt("nrates_step_"+str(step)+"_normal",diff_mev_n)



    return

def create_al_structures_for_analysis_SOAP():
    ''' create kmc structures for analyzing nearest neighbors with soap '''
    #############################################################
    ### create ase object with sphere of atoms around the vacancy
    #############################################################
    atomsc_sphere = create_al_sphere(a0=4.05,matrix_element="Al",cubic=True,ncell=4,nvac=1,cutoff=4.05,vacidx=0)

    f1 = "/Users/glensk/tmp/dataxx.quippy.xyz"
    f1n = "/Users/glensk/tmp/dataxx.name"
    f1_ = "/Users/glensk/tmp/dataxx_.quippy.xyz"
    f1n_ = "/Users/glensk/tmp/dataxx_.name"

    f11 = "/Users/glensk/tmp/dataxx11.quippy.xyz"
    f11n = "/Users/glensk/tmp/dataxx11.name"
    f111 = "/Users/glensk/tmp/dataxx111.quippy.xyz"
    f111n = "/Users/glensk/tmp/dataxx111.name"
    for i in [f1,f1n,f1_,f1n_,f11,f11n,f111,f111n]:
        if os.path.isfile(i):
            os.remove(i)

    def qwrite(f1,f1n,atomsc_sphere,text):
        atomsc_sphere.write(f1,format='quippy',append=True)
        f = open(f1n, "a")
        f.write(text+"\n")
        f.close()
        return

    #############################################################
    ### create fcc {Al,Mg,Si} with single vacancy
    #############################################################
    for matrix_element in ["Al", "Mg", "Si"]:
        atomsc_sphere = create_al_sphere(a0=4.05,matrix_element=matrix_element,cubic=True,ncell=4,nvac=1,cutoff=4.05,vacidx=0)
        ase_showpos(atomsc_sphere)
        qwrite(f1,f1n,atomsc_sphere,matrix_element+" pure")

        if matrix_element == "Al":
            for change in [1,2]: #,3]:  this symmetry seems to work
                for c in ["Mg", "Si"]:
                    frame = atomsc_sphere.copy()
                    frame[change].symbol = c
                    qwrite(f1,f1n,frame,"1NN_"+c+"_on_"+str(change))

        if matrix_element == "Al":
            for c in ["Mg", "Si"]:
                frame = atomsc_sphere.copy()
                frame[1].symbol = c
                frame[3].symbol = c
                #frame.write('dataxx.quippy.xyz',format='quippy',append=True)
                qwrite(f1,f1n,frame,"1NN_"+c+"on1__1NN_"+c+"_on3")
                qwrite(f1_,f1n_,frame,"1NN_"+c+"on1__1NN_"+c+"_on3")

        if matrix_element == "Al":
            for c in ["Mg", "Si"]:
                frame = atomsc_sphere.copy()
                frame[18].symbol = c
                frame[3].symbol = c
                #frame.write('dataxx.quippy.xyz',format='quippy',append=True)
                qwrite(f1,f1n,frame,"1NN_"+c+"on18__1NN_"+c+"_on3")
                qwrite(f1_,f1n_,frame,"1NN_"+c+"on18__1NN_"+c+"_on3")

        if matrix_element == "Al":
            frame = atomsc_sphere.copy()
            frame[1].symbol = "Mg"
            frame[3].symbol = "Si"
            #frame.write('dataxx.quippy.xyz',format='quippy',append=True)
            qwrite(f1,f1n,frame,"1NN_Mg__1NN_Si")

        if matrix_element == "Al":
            for change in [13,8]:  # beides 2NN
                for c in ["Mg", "Si"]:
                    frame = atomsc_sphere.copy()
                    frame[3].symbol = c            # 1NN
                    frame[change].symbol = c       # 2NN
                    qwrite(f1,f1n,frame,"1NN_"+c+"on3__2NN_"+c+"_on"+str(change))

        if matrix_element == "Al":
            for c in ["Mg", "Si"]:
                frame = atomsc_sphere.copy()
                frame[1].symbol = c
                frame[2].symbol = c
                frame[3].symbol = c
                qwrite(f1,f1n,frame,"1NN_"+c+"__1NN_"+c+"__1NN_"+c)
                qwrite(f11,f11n,frame,"1NN_"+c+"__1NN_"+c+"__1NN_"+c)
                qwrite(f111,f111n,frame,"1NN_"+c+"__1NN_"+c+"__1NN_"+c)

        if matrix_element == "Al":
            frame = atomsc_sphere.copy()
            frame[1].symbol = "Si"
            frame[2].symbol = "Si"
            frame[3].symbol = "Mg"
            #frame.write('dataxx.quippy.xyz',format='quippy',append=True)
            qwrite(f111,f111n,frame,"1NN_Si__1NN_Si__1NN_Mg")
    return



##################################################################################
## parametrize foces
##################################################################################
def get_michaels_paramerization(pos_all,force_all,NN1,alat,atoms,parametrize_only_idx=False,rcut=0.88,function=False,save_parametrization=False):
    ''' NN1 are the NN1idx
        atoms : an ase object of the perfect cell.

         1. get a function that for a given ATOM, gives all its 1NN (NN=a,b,c,d,...) (or its 2NN, 3NN, ...), basically, D for particular neighbors.
         2. write down the equation to solve: FORCE_VASP(ATOM)__x = SUM_over_all_neighbors[v_o_v__NN__x*Forcefunction(r)]
                                              FORCE_VASP(ATOM)__y = SUM_over_all_neighbors[v_o_v__NN__y*Forcefunction(r)]
                                              FORCE_VASP(ATOM)__z = SUM_over_all_neighbors[v_o_v__NN__z*Forcefunction(r)]
         3. This gives, for every displacement (=7), for every NN (=12), three (=x,y,z) equations.
         was ich jetzt loesen muss ist ein gleichungssystem von vielen vektoren D and die VASP force, mache , basically, D for particular neighbors.

    '''
    if save_parametrization != False:
        abcd_name = "poly_abcd_rcut"+str(rcut)+".dat"
        cd_name   = "poly_cd_rcut"+str(rcut)+".dat"
        parametrizationfile_abcd = save_parametrization+"/disp_fit.parameters."+abcd_name
        parametrizationfile_cd   = save_parametrization+"/disp_fit.parameters."+cd_name
        if os.path.isfile(parametrizationfile_abcd) and os.path.isfile(parametrizationfile_cd):
            params_abcd = np.loadtxt(parametrizationfile_abcd)
            params_cd   = np.loadtxt(parametrizationfile_cd)
            return params_abcd,params_cd

    a = []
    b = []

    a_ = []
    b_ = []

    a__ = []
    b__ = []
    rmin =  10000000000000000
    rmax = -10000000000000000
    for disp_idx in np.arange(len(pos_all)):    # for all displacements
    #disp_idx=7
        parametrize_over_atoms = NN1
        if parametrize_only_idx != False:
            parametrize_over_atoms = parametrize_only_idx

        for idxnn1,nn1 in enumerate(parametrize_over_atoms): # for all first NN
            Fv = force_all[disp_idx,nn1]
            for xyz in [0,1,2]:         # for xyz
                Fvx = Fv[xyz]

                NN1_1st_neigh = get_NN1_from_positions(atoms,alat,atomnr=nn1,verbose=False)
                ##print('idxnn1',idxnn1,nn1,"NN1_1st_neigh:",NN1_1st_neigh)
                #print("Fv",Fv,"===")
                Dx = np.zeros(12)
                R  = np.zeros(12)
                Rcut  = np.zeros(12)
                for idxnn11,nn11 in enumerate(NN1_1st_neigh):    # for all neighbors of a particular 1NN
                    D  = pos_all[disp_idx,nn1]-pos_all[disp_idx,nn11]
                    #if parametrize_only_idx != False:
                    #    #print("D",D)
                    #    #print('i will take this D!')
                    #    print("D",D,'nn1',nn1,'nn11',nn11)
                    Drel = D/alat               # given for every particle
                    vecnorm = r = np.linalg.norm(Drel)   # given for every particle
                    if r > rmax:
                        rmax = r
                    if r < rmin:
                        rmin = r
                    v_o_v = Drel/vecnorm
                    R[idxnn11] = r
                    Rcut[idxnn11] = rcut
                    Dx[idxnn11] = v_o_v[xyz]
                    if function != False:
                        force_add = function(R)*Dx[idxnn11]

                #print('idxnn11',idxnn11,'p1',pos_all[disp_idx,nn1],'p2',pos_all[disp_idx,nn11],'D',D,"Drel",Drel,'v_o_v',v_o_v,'r',r)
                # this is the model without fixing anything:
                # Forces: a + b*r**(-1) + c*r**(-2) + d*r**(-3)
                # Forces * xyz : (a + b*r**(-1) + c*r**(-2) + d*r**(-3) #+ e*r**(-4)) * x  # where x = xyz = v_o_v
                # for every vasp force:
                # Force_vasp_x(r=R) = SUM_over_all_neigh[ax + bx/R + cx/R^2 + dx/R^3]  ## R is given; x = Dx
                # ( a + b/(r_1)^1 + c/(r_1)^2 + d/(r_1)^3 ) * x_1 +
                # ( a + b/(r_2)^1 + c/(r_2)^2 + d/(r_2)^3 ) * x_2
                #  ==
                #  a ( x_1          + x_2          + x_3 ... )  # = np.sum(Dx)
                #  b ( x_1/(r_1)    + x_2/(r_2)    + ...        # = np.sum(Dx/R)
                #  c ( x_1/(r_1)^2  + x_2/(r_2)^2  + ...        # = np.sum(Dx/R**2.)
                #  d ( x_1/(r_1)^3  + x_2/(r_2)^3  + ...        # = np.sum(Dx/R**3.)
                #  --> a*np.sum(Dx) + b*np.sum(Dx/R) + c*np.sum(Dx/R**2.) + d*np.sum(Dx/R**3) = Fvx  # **1)

                #############################
                # allgemein ohne constraints **1)
                #############################
                a_add =  [np.array([ np.sum(Dx), np.sum(Dx/R), np.sum(Dx/R**2.), np.sum(Dx/R**3.) ])]
                a += a_add
                b += [Fv[xyz]]

                #############################
                # Forces(rcut=0.88) = 0
                #############################
                a_add_ =  [np.array([ np.sum(Dx/R), np.sum(Dx/R**2.), np.sum(Dx/R**3.) ])]
                a_ += a_add_
                no_a_corr_F =  np.sum(Dx/rcut) + np.sum(Dx/rcut**2.)  + np.sum(Dx/rcut**3.)
                b_ += [Fv[xyz]+no_a_corr_F]


                # with constraints:
                # V[r_] := aa + dd/r^3 + cc/r^2 + bb/r
                # V'[r]  = -((3 dd)/r^4) - (2 cc)/r^3 - bb/r^2
                # V''[r] = (12 dd)/r^5 + (6 cc)/r^4 + (2 bb)/r^3
                #
                # Solve[V[rcut] == 0, aa]       --> aa -> (-dd - cc rcut - bb rcut^2)/rcut^3
                # Vnoa[r_] := (-dd - cc rcut - bb rcut^2)/rcut^3 + dd/r^3 + cc/r^2 + bb/r
                # Solve[Vnoa'[rcut] == 0, bb]   -->  bb -> (-3 dd - 2 cc rcut)/rcut^2
                #
                #
                #     a = - b*rcut**(-1) - c*rcut**(-2) - d*rcut**(-3) # a(b,c,d)   -- > && a * np.sum(Dx)
                a_add_ =  [np.array([ np.sum(Dx/R), np.sum(Dx/R**2.), np.sum(Dx/R**3.) ])]
                a_ += a_add_
                no_a_corr_F =  np.sum(Dx/rcut) + np.sum(Dx/rcut**2.)  + np.sum(Dx/rcut**3.)
                b_ += [Fv[xyz]+no_a_corr_F]

                # with constraints:
                #     a = - b*rcut**(-1) - c*rcut**(-2) - d*rcut**(-3) # a(b,c,d)   -- > && a * np.sum(Dx)
                #     b = - 3.*d/rcut**2. - 2.*c/rcut  # from first der b(d,c)
                a_add__ =  [np.array([ np.sum(Dx/R**2.)+np.sum(Dx/(rcut**2.))-np.sum((Dx*2.)/(R*rcut)), np.sum(Dx/R**3.)+np.sum((2.*Dx)/(rcut**3.))-np.sum(3*Dx/(R*rcut**2.)) ])]
                a__ += a_add__

                #print('xyz',xyz,a_add)
            #if parametrize_only_along_110 == True:
            #    print('Dx',Dx)
            #    print('R ',R)
            #print('Dx/R ',Dx/R)
            #print('a',np.sum(Dx)) #,np.sum(Dy),np.sum(Dz))
            #print('b',np.sum(Dx/R)) #,np.sum(Dy/R),np.sum(Dz/R))
            #print('c',np.sum(Dx/R**2.)) #,np.sum(Dy/R**2.),np.sum(Dz/R**2.))
            #print('d',np.sum(Dx/R**3.)) #,np.sum(Dy/R**3.),np.sum(Dz/R**3.))
    a = np.array(a)
    b = np.array(b)
    a_ = np.array(a_)
    b_ = np.array(b_)
    a__ = np.array(a__)
    b__ = np.array(b__)
    #print('b',b,'len(b)',len(b))
    #print('a',a,'len(a)',len(a))
    params  = np.linalg.lstsq(a, b,rcond=None)[0]
    #params_ = np.linalg.lstsq(a_, b_,rcond=None)[0]
    params__ = np.linalg.lstsq(a__, b,rcond=None)[0]
    #print('params      ',params)
    #aa = - params_[0]*rcut**(-1) - params_[1]*rcut**(-2) - params_[2]*rcut**(-3)
    #print('params_ (i) ',params_)
    #params_ = [aa,params_[0],params_[1],params_[2]]
    #print('params_ (ii)',params_)


    bb = -3*params__[1]*rcut**(-2) -2.*params__[0]/rcut
    aa = - bb*rcut**(-1) - params__[0]*rcut**(-2) - params__[1]*rcut**(-3)
    #print('params__(i) ',params__)
    params_abcd = [params[0],params[1],params[2],params[3],rcut,rmin,rmax]
    params_cd   = [aa,bb,params__[0],params__[1],rcut,rmin,rmax]
    #print('params__(ii)',params__)
    #for i in params:
    #    print(i)
    print('params_abcd:',params_abcd)
    print('params_cd  :',params_cd)
    print('params_abcdx',params_abcd[:4])
    print('params_cd  x',params_cd[:4])
    if save_parametrization != False:
        mmodel_abcd = hesse.Michael_poly_der
        print('rmin',rmin)
        print('rmax',rmax)
        if rmax < 0.9:
            rmax = 1.2
        if rmin > 0.5:
            rmin = 0.5
        xred_dense = np.arange(rmin,rmax,0.001)
        np.savetxt(save_parametrization+"/disp_vs_forces_"+abcd_name,np.array([xred_dense,-1*mmodel_abcd(xred_dense,*(params_abcd[:4]))]).T)
        np.savetxt(save_parametrization+"/disp_vs_forces_"+cd_name,np.array([xred_dense,-1*mmodel_abcd(xred_dense,*(params_cd[:4]))]).T)
        print('saved',save_parametrization+"/disp_vs_forces_"+abcd_name)
        print('saved',save_parametrization+"/disp_vs_forces_"+cd_name)
        np.savetxt(parametrizationfile_abcd,params_abcd)
        np.savetxt(parametrizationfile_cd,params_cd)
        print('saved',parametrizationfile_abcd)
        print('saved',parametrizationfile_cd)


    return params_abcd,params_cd

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
    if ace.ipi:
        # write fix 1 all ipi md_ff 77776 unix
        f.write('fix 1 all ipi md_ff 77776 unix\n')
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
    ] + ace.pot.atom_masses_str #lammps_command_masses()
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
    elif ace.pot.pottype  in [ 'eam', 'eam-alloy' ]:
        command = command + ace.lammps_command_potential_eam_alloy(pair_style="eam/alloy")
    elif ace.pot.pottype == 'adp':
        command = command + ace.lammps_command_potential_eam_alloy(pair_style="apd")
    else:
        print('ace.pot.pottype:',ace.pot.pottype)
        sys.exit('potential now known yet')

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
    SCR = os.environ["SCRIPTS"]
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
    if ace.verbose > 2:
        print("lammps_ext_calc (1): get_elastic_constants",get_elastic_constants)

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

    if ace.verbose > 2:
        print('lammps_ext_calc (2) doing the calculation in lammps, externally, in ace.pot.lammps_tmpdir',ace.pot.lammps_tmpdir)

    ###############################################################
    # write input structure (pos.lmp)
    ###############################################################
    if ace.verbose > 2:
        show_ase_atoms_content(atoms,showfirst=10,comment="START LAMMPS EXTERNALLY")
    atoms.set_calculator(None)
    if ace.verbose > 2:
        print('ace.pot.pottype',ace.pot.pottype)

    if ace.pot.pottype in [ 'runner' , 'n2p2' ]:
        atoms.write(ace.pot.lammps_tmpdir+'pos.lmp',format='lammps-runner',pot=ace.pot)
    elif ace.pot.pottype in [ 'eam', 'eam-alloy' ]:
        atoms.write(ace.pot.lammps_tmpdir+'pos.lmp',format='lammps-data')
    elif ace.pot.pottype in [ 'adp' ]:
        atoms.write(ace.pot.lammps_tmpdir+'pos.lmp',format='lammps-data')
    else:
        sys.exit('44323 Error: dont know how to write this lammps file')
    if ace.verbose > 2:
        print('lammps_ext_calc (3) written ',ace.pot.lammps_tmpdir+'/pos.lmp')
    ###############################################################
    # write in.lmp (for mormal energy calculation)
    ###############################################################
    if not get_elastic_constants:
        lammps_in_file = 'in.lmp'
        lammps_write_inputfile(folder=ace.pot.lammps_tmpdir,filename=lammps_in_file,positions='pos.lmp',ace=ace)
        if ace.verbose > 2:
            print('written (77) ',ace.pot.lammps_tmpdir+'/'+lammps_in_file)
    if ace.verbose > 2:
        print("lammps_ext_calc (4) written lammsp inputfile to ",ace.pot.lammps_tmpdir)

    ### calculate with lammps (trigger externally)
    ene = False
    if os.path.isfile(ace.pot.lammps_tmpdir+'log.lammps'):
        os.remove(ace.pot.lammps_tmpdir+"log.lammps")


    ###############################################################
    # IF get elastic constants
    ###############################################################
    if get_elastic_constants:
        lammps_in_file = 'in.elastic'
        if ace.verbose > 2:
            print('written ',ace.pot.lammps_tmpdir+'/potential.mod')
        lammps_write_inputfile_from_command(folder=ace.pot.lammps_tmpdir,filename='potential.mod',command=lammps_ext_elastic_potential_mod(ace))
        cp(os.environ['scripts']+'/lammps_scripts/elastic/in.elastic',ace.pot.lammps_tmpdir+'/in.elastic')
        cp(os.environ['scripts']+'/lammps_scripts/elastic/displace.mod',ace.pot.lammps_tmpdir+'/displace.mod')
        lammps_write_inputfile_from_command(folder=ace.pot.lammps_tmpdir,filename='init.mod',command=lammps_ext_elastic_init_mod(ace,positions='pos.lmp'))
        if ace.verbose > 2:
            print('written ',ace.pot.lammps_tmpdir+'/'+lammps_in_file)

    ###############################################################
    # cd to folder and run lammps
    ###############################################################
    if ace.verbose > 2:
        print('lammps_ext_calc (5) cd ..')
    with cd(ace.pot.lammps_tmpdir):  # this cd's savely into folder
        # RUN LAMMPS
        # without SHELL no LD LIBRARY PATH
        #print('pwd',os.getcwd())
        if ace.verbose > 2:
            print('lammps_ext_calc (6)',ace.LAMMPS_COMMAND+" < "+lammps_in_file+" > /dev/null")
        call([ace.LAMMPS_COMMAND+" < "+lammps_in_file+" > /dev/null"],shell=True)

        ###############################################################
        # extract energy and forces
        ###############################################################
        if not get_elastic_constants:
            if ace.verbose > 2:
                print('lammps_ext_calc (7) pwd',os.getcwd())
            ene = check_output(["tail -300 log.lammps | grep -A 1 \"Step Temp E_pai\" | tail -1 | awk '{print $3}'"],shell=True).strip()
            if ace.verbose > 2:
                print('lammps_ext_calc (8) ene from grep:',ene,'(eV/cell)')
            ene=float(ene)
            ene_ev_cell = deepcopy(ene)
            #print('ene',ene,'lammps in eV')
            #print('ace units',ace.units)
            #print('ace.units.lower()',ace.units.lower())
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
            # get_elastic_constants == True
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
            ene_ev_cell = deepcopy(ene)

    if ace.verbose > 2:
        print('lammps_ext_calc (9) ene final:',ene,"chosen units")
    if ace.verbose > 2:
        show_ase_atoms_content(atoms,showfirst=10,comment="LAMMPS EXTERNALLY")
    #print('eneneee',ene_ev_cell,ene)
    #sys.exit()
    return ene_ev_cell, ene

##################################################################################
## phonons / phonopy
##################################################################################

def ase_to_phonopy(ase_structure):
    phonopy_structure    = PhonopyAtoms(
        symbols          = ase_structure.get_chemical_symbols(),
        scaled_positions = ase_structure.get_scaled_positions(),
        cell             = ase_structure.get_cell()
        )
    return phonopy_structure

def phonopy_to_ase(phonopy_structure):
    phonopy_structure    = Atoms(
        symbols          = phonopy_structure.get_chemical_symbols(),
        scaled_positions = phonopy_structure.get_scaled_positions(),
        cell             = phonopy_structure.get_cell(),
        pbc              = True)
    return phonopy_structure

def phonopy_pre_process(cell, supercell_matrix=None):
    if supercell_matrix is None:
        smat = [[2,0,0], [0,2,0], [0,0,2]],
    else:
        smat = supercell_matrix
    from phonopy.units import AbinitToTHz
    phonon = Phonopy(cell,
                     smat,
                     primitive_matrix='auto')
                     #primitive_matrix=[[0, 0.5, 0.5],
                     #                  [0.5, 0, 0.5],
                     #                  [0.5, 0.5, 0]],
                     #factor=AbinitToTHz)
    phonon.generate_displacements(distance=0.03)
    print("[Phonopy] Atomic displacements:")
    disps = phonon.get_displacements()
    for d in disps:
        print("[Phonopy] %d %s" % (d[0], d[1:]))
    return phonon


def phonopy_post_process(phonon, set_of_forces):
    phonon.produce_force_constants(forces=set_of_forces)
    print('')
    print("[Phonopy] Phonon frequencies at Gamma:")
    for i, freq in enumerate(phonon.get_frequencies((0, 0, 0))):
        print("[Phonopy] %3d: %10.5f THz" %  (i + 1, freq)) # THz

    # DOS
    phonon.set_mesh([21, 21, 21])
    phonon.set_total_DOS(tetrahedron_method=True)
    print('')
    print("[Phonopy] Phonon DOS:")
    for omega, dos in np.array(phonon.get_total_DOS()).T:
        print("%15.7f%15.7f" % (omega, dos))

    ## band structure
    if False:
        phonon.run_mesh([20, 20, 20])
        phonon.run_total_dos()
        phonon.plot_total_dos().show()

    if False:
        phonon.run_mesh([20, 20, 20], with_eigenvectors=True, is_mesh_symmetry=False)
        phonon.run_projected_dos()
        phonon.plot_projected_dos().show()

    if True:
        phonon.run_mesh([20, 20, 20])
        phonon.run_thermal_properties(t_step=1,
                              t_max=1000,
                              t_min=0)
        tp_dict = phonon.get_thermal_properties_dict()
        temperatures = tp_dict['temperatures']
        free_energy = tp_dict['free_energy']
        entropy = tp_dict['entropy']
        heat_capacity = tp_dict['heat_capacity']

        for t, F, S, cv in zip(temperatures, free_energy, entropy, heat_capacity):
            print(("%12.3f " + "%15.7f" * 3) % ( t, F, S, cv ))

    return

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
    #print('inputnn',inputnn)
    inputdata   = inputnn.replace("input.nn","input.data")
    inputdatanr = inputnn.replace("input.nn","input.data_nr")
    #print('inputdatanr',inputdatanr)
    nr = False
    if os.path.isfile(inputdatanr):
        nr = np.loadtxt(inputdatanr)
        nr = int(nr)
        #print('from loaded file',nr)
        return nr
    if nr == False:
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
                #print('iii',i)
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
    pot_elements, pot_atom_energy = inputnn_get_atomic_symbols_and_atom_energy_dict(inputnn)
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
    verbose=True
    if verbose > 2:
        print('myutils.py (7) inputnn',inputnn)
    filename = n2p2_runner_get_learning_curve_filename(inputnn)
    if verbose > 2:
        print('myutils.py (7) learning_curve_filename:',filename)
    if filename == False:
        return False
    basename = os.path.basename(filename)  # "learning-curve.out"
    folder = os.path.abspath(filename.replace(basename,''))
    if False: #verbose:
        print('learning_curve folder  :',folder)
        print('learning_curve basename:',basename)
    n2p2_runner = type = inputnn_runner_or_n2p2(folder+'/input.nn')
    if False: #verbose:
        print('nn',n2p2_runner)
    if not os.path.isfile(filename):
        #sys.exit(filename+" does not exist! (32)")
        print(filename+" does not exist! (32)")
        return False

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

def n2p2_runner_get_bestteste_idx(inputnn,verbose=False):
    if verbose:
        print('inputnn',inputnn)
    job_finished = n2p2_runner_job_finished(inputnn)
    if verbose:
        print('job_finished',job_finished)
    best_testsete_file = inputnn.replace('input.nn', 'best_testsete')
    if verbose:
        print('best_testsete_file',best_testsete_file)
    if os.path.isfile(best_testsete_file):
        best_testsete = int(np.loadtxt(best_testsete_file))
        #print('from best_testsete_file',best_testsete_file)
        return best_testsete
    else:
        #print('from lc')
        learning_curve = lc = n2p2_runner_get_learning_curve(inputnn)
        if verbose:
            print('lc',lc)
        if type(lc) == bool:
            return False
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
    readmepath = create_READMEtxt()
    os.chdir(folder)
    cp("../input.data")
    cp("../input.nn")
    submitskript = n2p2_write_submit_skript(directory=False,cores=cores,nodes=1,job="scaling",days=days,hours=hours,minutes=minutes,interactive=interactive)
    submitjob(submit_to_que=submit_to_que,submit_to_debug_que=submit_to_debug_que,jobdir=False,submitskript=submitskript)
    #submitjob(submit=True,submit_to_que=submit_to_que,submit_to_debug_que=submit_to_debug_que,jobdir=False,submitskript=submitskript)
    readmepath = create_READMEtxt(add="submit_to_debug_que = "+str(submit_to_debug_que))
    os.chdir(hier)
    return

def n2p2_write_submit_skript(directory=False,nodes=1,cores=28,job=False,interactive=False,days=0,hours=72,minutes=0,seconds=0):
    ''' wiretes a submit_n2p2_{get_scaling,training}.sh file '''
    if job not in ["scaling","train"]:
        sys.exit('job has to be one of nnp-XXX jobs as "train, scaling, ..."')
    myhost = os.environ["myhost"]
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

        text_file.write("module list\n")
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
    readmepath = create_READMEtxt()
    return

def get_latest_n2p2_pot():
    checkfor = scripts()+"/potentials/n2p2_v*ag/"
    #checkfor = scripts()+"/potentials/n2p2_v4ag_ppl_987654_21cores/"
    print('checkfor',checkfor)
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



##################################################################################
## getting positions
##################################################################################
def folder_get_pos_forces_cell(i):
    ''' i is the folder '''
    create_dofor_POSITIONs = False
    create_dofor_u_OUTCAR = False
    if not os.path.isfile(i+'/pos'):
        with cd(i):
            call(["OUTCAR_positions-last-ARRAY.sh > pos"],shell=True)
    if not os.path.isfile(i+'/forces'):
        with cd(i):
            call(["OUTCAR_forces-last-ARRAY.sh > forces"],shell=True)

    if not os.path.isfile(i+'/POSITIONs'):
        with cd(i):
            #print('pwd',os.getcwd())
            call(["extractPOSITIONS.sh"],shell=True)
    pos_forces  = np.loadtxt(i+'/POSITIONs')
    pos         = pos_forces[:,[0,1,2]]
    forces      = pos_forces[:,[3,4,5]]
    if not os.path.isfile(i+'/u_OUTCAR'):
        with cd(i):
            call(["OUTCAR_ene-potential_energy.sh > u_OUTCAR"],shell=True)

    if not os.path.isfile(i+'/u_OUTCAR_without_first_substracted'):
        with cd(i):
            call(["OUTCAR_ene-potential_energy_without_first_substracted.sh > u_OUTCAR_without_first_substracted"],shell=True)

    #@if create_dofor_POSITIONs == True:
    #@    with cd(i):
    #@        call(["cat POSITIONs >> ../POSITIONs"],shell=True)
    #@if create_dofor_u_OUTCAR == True:
    #@    with cd(i):
    #@        call(["cat u_OUTCAR >> ../u_OUTCAR"],shell=True)
    #alat = float(i.split("vasp4/")[1].split("Ang")[0])
    if os.path.isfile(i+'/cell'):
        cell = np.loadtxt(i+'/cell')
        #print('cell1',cell)
    elif os.path.isfile(i+'/POSCAR'):
        cell = POSCAR_get_cell(i+'/POSCAR')
    else:
        sys.exit('could not determine cell dimensions')
    #print(i)
    #print('cell',cell)
    #sys.exit()
    return pos,forces,cell

def check_vec_in_array(test,array):
    #for idx,x in enumerate(array):
    #    print('idx',idx,'x',x,'test',test)
    #print()
    #return any(np.array_equal(x, test) for x in array)
    return any(np.allclose(x, test) for x in array)

def try_to_get_from_positions_only__nat_cryst_sc_cell(positions,forces,or_folder_to_positions=False,verbose=False):
    '''
    nat,cryst,sc,alat,cell,pos0 = try_to_get_from_positions_only__nat_cryst_sc_cell(pos,forces,verbose=False)

    pos: (all) positions from POSITIONs
    forces: (all) forces from POSITIONs
    %head POSITIONs
    0.00000 0.00000 0.00000 0.000000 -0.000000 -0.000000
    0.00000 0.00000 4.13000 0.000000 -0.000000 0.000000
    0.00000 4.13000 0.00000 0.000000 -0.000000 -0.000000
    0.00000 4.13000 4.13000 -0.000000 -0.000000 -0.000000
    4.13000 0.00000 0.00000 -0.000000 -0.000000 0.000000
    4.13000 0.00000 4.13000 -0.000000 0.000000 -0.000000
    4.13000 4.13000 0.00000 -0.000000 0.000000 -0.000000
    4.13000 4.13000 4.13000 0.000000 -0.000000 0.000000
    0.00000 2.06500 2.06500 -0.000000 0.000000 -0.000000
    0.00000 2.06500 6.19500 0.000000 -0.000000 0.000000
    '''
    verbose = True
    if or_folder_to_positions != False:
        if not os.path.isfile(or_folder_to_positions):
            sys.exit('ftp does not exist 85;Exit')
        else:
            positions_forces = np.loadtxt(or_folder_to_positions)
            positions = positions_forces[:,[0,1,2]]
            forces = positions_forces[:,[3,4,5]]

    #########################
    # get nat (e.g. 32 atoms)
    #########################
    f = forces
    p = positions

    if verbose:
        print('positions')
        print(p)
        print('forces')
        print(f)
    nat = 1
    for i in np.arange(len(forces)):
        nat = i+1
        fa = np.abs(f[i]).max()
        if verbose:
            print('fa',i,fa)
        if fa > 1e-4:
            nat = i
            break
        if i == 950:
            print('hit 950 entries')
            sys.exit()
    if verbose:
        print('nat:',nat)
    if nat < 4:
        sys.exit('did not find the number of atoms: suggest to take first step, and scan for next occurenct of positions of first step +small epsilon (due to MD step)')
    # first check fcc up to sc = 6:
    #########################
    # get cryst (fcc or bcc)
    # get sc (N=2)
    #########################
    cryst = False
    sc = False
    for i in np.arange(1,7):
        ifcc = ((i**3.)*4).astype(int)
        ibcc = ((i**3.)*2).astype(int)
        #print(i,ifcc,'fcc') # case of fcc
        #print(i,ibcc,'bcc') # case of fcc
        if ifcc == nat:
            cryst = 'fcc'
            sc = i
        if ibcc == nat:
            cryst = 'fcc'
            sc = i
    if verbose:
        print('cryst',cryst)
        print('sc',sc)
    #########################
    # get alat (e.g. 4.05 angstrom)
    #########################
    pos0 = p[:nat]
    if verbose:
        print('pos0')
        print(pos0)
        print('pos0/sc')
        print(pos0/sc)
    alat = 0
    if pos0[1,0] == 0 == pos0[1,1]:
        #print('y1')
        if pos0[1,2] > 2.0:
            alat = pos0[1,2]
            #print('y2',alat)
    if verbose:
        print('alat',alat)
    if alat == 0:
        sys.exit('did not find alat')


    #########################
    # get cell (e.g. [[4.05,0,0],[0,4.05,0],[0,0,4.05]])
    #########################
    cell = np.eye(3, dtype=int)*sc*alat
    if verbose:
        print('cell')
        print(cell)
    return nat, cryst, sc, alat, cell, pos0

def try_to_get_alat_and_sc_from_cell_and_positions(cell,positions):
    ''' this should work for all cubic cells (fcc, bcc, dc) '''
    sc = 0
    alat = 0

    c = cell
    celltype = False
    if c[0,0] == c[1,1] == c[2,2] and c[0,1] == c[0,2] == c[1,0] == c[1,2] == c[2,0] == c[2,1] == 0:
        celltype = 'cubic'
        sc = 1
        scalat = c[0,0]
        #print('c[0,0] == scalat',scalat)
        for i in np.arange(1,10):
            alltrue = False
            #print('i (=sc) =',i,'scalat/i',scalat/i)
            for j in np.arange(1,i):
                pos = np.array([0,0,(scalat/i)*j])
                chk = check_vec_in_array(pos, positions)
                #print('pos?',pos,chk)
                if chk == True:
                    alltrue = True
                elif chk == False:
                    alltrue = False
                    break
            #print('i (=sc) =',i,'scalat/i',scalat/i,"alltue?",alltrue)
            if alltrue == True:
                sc = i
                alat = c[0,0]/sc
            #print('i (=sc) =',i,'scalat/i',scalat/i,"alltue?",alltrue,"--> sc:",sc,"alat",alat)
    #print('sc --------->',sc)
    #print('alat ------->',alat)

    if alat > 6 or alat < 3:
        print('alat',alat)
        print('cell')
        print(cell)
        print('celltype:',celltype)
        print()
        print('c[0,0] == c[1,1] == c[2,2]:',c[0,0] == c[1,1] == c[2,2])
        print('c[0,1] == c[0,2] == c[1,0] == c[1,2] == c[2,0] == c[2,1] == 0:',c[0,1] == c[0,2] == c[1,0] == c[1,2] == c[2,0] == c[2,1] == 0)
        sys.exit("this does not seem to be a correct alat")
    return sc,alat

def center_positions_around_0(pos,sc,alat):
    pos_ = np.copy(pos)
    # this brings the displaced atom into the middle of the cell.
    for iidx,i in enumerate(pos):  # every atomic positions
        #print(iidx,i)
        for jdx,j in enumerate(i):
            #print(jdx,j,'--',pos_[iidx,jdx],'sc',self.sc,'alat',self.alat,'--',self.sc*self.alat/2)
            if pos_[iidx,jdx] > sc*alat/2.:
                pos_[iidx,jdx] = pos_[iidx,jdx]-sc*alat
    return pos_

def POSCAR_get_cell(path_to_POSCAR):
    head = read_firstlines(path_to_POSCAR,7)
    wc3 = []
    for idx,i in enumerate(head):
        #print(i)
        #print(i.split())
        line = i.split()
        wc = len(line)
        #print('idx',idx,'line',line,'wc',wc)
        if wc == 3: wc3 += [idx]
    #print('wc3',wc3)
    for idx,i in enumerate(wc3):
        linenr = wc3[idx]
        if wc3[idx] == linenr and wc3[idx+1] == linenr+1 and wc3[idx+2] == linenr+2:
            #print('aha',linenr,'idx',idx)
            cell = head[wc3[idx]:wc3[idx]+3]
            break
    #print('ok',cell)
    cell_ = np.zeros((3,3))
    for idx,j in enumerate(cell):
        #print('-->',[float(i) for i in j.split()])
        cell_[idx] = [float(i) for i in j.split()]
    cell = cell_
    #frame = ase_read(i+'/POSCAR',format='vasp')
    #cell = frame.cell
    #print('cell2',cell)
    return cell

##################################################################################
## ase funcions
##################################################################################
def load_POSITIONs_to_ase_frames(folder,matrix_element="Al"):
    ''' to make ase object ouf ot this I would need the cell '''
    posxyz = folder+"/POSITIONs.extxyz"
    if os.path.isfile(posxyz):
        print('reading frames',posxyz)
        frames = ase_read(posxyz,index=":")
        print('reading frames done')
        #print('type',type(frames))
        #print(frames[0].cell)
        #print(frames[0].info)
        return frames

    print('reading POSITIONs ...')
    pos_forces = np.loadtxt(folder+"/POSITIONs")
    print('done ...')
    #print('reading outcar...')
    #pos = ase_read(folder+'/OUTCAR',format='vasp-out',index=0) # works
    #pos = ase_read(folder+'/OUTCAR',format='vasp-out',index=":") # works
    #print('pos.type',type(pos))
    #print('len(pos)',len(pos))
    #print('done ...')
    #pos     = pos_forces[:,[0,1,2]].reshape((-1,nat,3))
    #forces  = pos_forces[:,[3,4,5]].reshape((-1,nat,3))
    pos     = pos_forces[:,[0,1,2]].reshape((-1,3))
    forces  = pos_forces[:,[3,4,5]].reshape((-1,3))
    #print(pos[0])
    #print()
    #print(pos[1])
    #print(np.abs(pos[:]).max())
    #print(np.abs(pos[0]).max())
    #print(np.abs(pos[1]).max())
    #print(np.abs(pos[2]).max())
    #print(np.abs(pos[3]).max())
    #print(np.abs(pos[1]).max())
    #print(forces[0])
    nat,cryst,sc,alat,cell,pos0 = try_to_get_from_positions_only__nat_cryst_sc_cell(pos,forces,verbose=False)

    ## reshape pos forces
    pos     = pos_forces[:,[0,1,2]].reshape((-1,nat,3))
    forces  = pos_forces[:,[3,4,5]].reshape((-1,nat,3))
    #print(pos.shape)
    steps = (pos.shape[0])
    #print('steps',steps)
    if os.path.isfile(folder+"/u_OUTCAR"):
        energies = np.loadtxt(folder+"/u_OUTCAR")*(nat-1.)/1000.
    elif os.path.isfile(folder+"/u_OUTCAR_without_first_substracted"):
        energies = np.loadtxt(folder+"/u_OUTCAR_without_first_substracted")*(nat-1.)/1000.
    else:
        energies = np.zeros(steps)

    #print('pos0')
    #print(pos0)
    bulk = get_ase_atoms_object_kmc_al_si_mg_vac(sc,a0=alat,matrix_element=matrix_element,cell='cubic',create_fake_vacancy=False,crystal_structure=cryst,normal_ordering=True,pos0=pos0,verbose=False)
    #print('bulk')
    #print(bulk.positions[:3])
    #print(bulk.cell)
    #print(bulk.get_pbc())
    #print(type(bulk))

    # reassign all positions/forces
    #print('setting pos forces')
    #print('ene3',energies[:3])
    frames = []
    for idx in np.arange(steps):
        bulk_new = bulk.copy()  # needs to be a new object
        #print(idx,'bn',type(bulk_new),bulk_new.cell)
        #print(idx,'bp',type(bulk_new),bulk_new.positions[:2])
        #print(idx,'b?',type(bulk_new),pos[idx][:2])
        bulk_new.set_positions(pos[idx])
        #if idx < 4:
        #    print('ene',idx,energies[idx])
        #print('bni',bulk_new.info)
        #bulk_new.info = {'comment': 'ka'}
        bulk_new.info['cryst'] = cryst
        bulk_new.info['alat'] = alat
        bulk_new.info['sc'] = sc
        calc = SinglePointCalculator(bulk_new, energy=energies[idx], forces=forces[idx])
        bulk_new.set_calculator(calc)
        #bulk_new.set_forces(forces[idx])
        #bulk_new.forces = forces[idx]
        #print(idx,'bx',type(bulk_new),bulk_new.cell)
        #frames.append(bulk.set_positions(pos[idx]))
        frames.append(bulk_new)
        #if idx == 4:
        #    sys.exit()
    verbose = False
    if verbose:
        print('done; len(frames)',len(frames),frames[0])
        print('done; type(frames)',type(frames),type(frames[0]))
        print()
    #for idx in np.arange(3):
    #    print(idx,'pp',pos[idx][:2])
    #    frames[idx].set_positions(pos[idx])
    #    #frames[idx].forces = pos[idx]
    if verbose:
        print()
        print('f[x]: out')
        print('f[0]:',frames[0].positions[:1],frames[0].get_forces()[:1])
        print('f[1]:',frames[1].positions[:1],frames[1].get_forces()[:1])
        print('f[2]:',frames[2].positions[:1],frames[2].get_forces()[:1])
        print()
        print('ref: pos')
        print('ref pos[0]:',pos[0][:1],forces[0][:1])
        print('ref pos[1]:',pos[1][:1],forces[1][:1])
        print('ref pos[2]:',pos[2][:1],forces[2][:1])
    print('write POSITIONs.extxyz',posxyz)
    ase_write(posxyz,frames,format='extxyz')
    print('write POSITIONs.extxyz done')
    return frames

def get_all_1NN_dist_vec_from_POSITIONs(folder,matrix_element="Al",verbose=False):
    '''
    dimensions:
        NN1_distances = np.zeros((steps,nat,NN1_neighbors,3))
    '''
    frames = load_POSITIONs_to_ase_frames(folder,matrix_element=matrix_element)
    cryst = frames[0].info['cryst']

    print('cryst',cryst)
    import h5py
    if os.path.isfile(folder+'/POSITIONs.distances.h5'):
        print('found',folder+'/POSITIONs.distances.h5')
        h5f = h5py.File(folder+'/POSITIONs.distances.h5','r')
        NN1_distances = h5f['NN1_distances'][:]
        NN1 = h5f['NN1'][:]
        h5f.close()
    else:
        print('creaeing h5f ...')
        steps = len(frames)
        print('steps',steps)

        ## get the ase object
        #if False:
        #    bulk = poscar = ase_read(folder+'/POSCAR',format='vasp')
        #    print(poscar.positions)

        ## get the list of first NN
        NN1 = get_NN1_idx_for_every_atom(frames[0],crystal_structure=cryst)
        if verbose:
            print('NN1.shape',NN1.shape)
            print(NN1)
        NN1_neighbors = NN1.shape[1]
        nat = NN1.shape[0]
        if verbose:
            print('nat          :',nat)
            print('NN1_neighbors:',NN1_neighbors)


        NN1_distances = np.zeros((steps,nat,NN1_neighbors,3))
        # get distance vectors of 1NN for every step
        for step,frame in enumerate(frames):
            for atom in np.arange(nat):
                dall = frame.get_distances(atom,NN1[atom],mic=True,vector=True)
                #print('step',step,'atom',atom,'atompos',frame.positions[atom],'dall',dall.shape)
                NN1_distances[step,atom] = dall
        print('saving POSITIONs.distances.h5')
        h5f = h5py.File(folder+'/POSITIONs.distances.h5', 'w')
        h5f.create_dataset('NN1_distances', data=NN1_distances)
        h5f.create_dataset('NN1', data=NN1)
        h5f.close()
        print('NN1_distances done')

    # it would be nice to have only the unique distances
    # meaning, every interaction between atoms only once (and not twice)
    # this would ~half size of the NN1_distacnes
    return NN1_distances,NN1,frames


def get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi=0,nmg=0,nvac=0,a0=False,matrix_element="Al",cell='cubic',create_fake_vacancy=False,crystal_structure="fcc",normal_ordering=True,pos0=False,verbose=False):
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
        cell: {cubic, primitive, orthorhombic}
            Construct cubic unit cell
    """
    ###################################################
    # define the cell: cubic, primitive, orthorhombic
    ###################################################
    if cell not in [ 'cubic', 'primitive', 'orthorhombic' ]:
        sys.exit('cell has to be one of: cubic, primitive, orthorhombic')
    cubic = orthorhombic = primitive = False
    if cell == 'cubic': cubic = True
    if cell == 'orthorhombic': orthorhombic = True

    ##########################################################################
    # define the crystal structure: sc, fcc, bcc, hcp, diamond, zincblende,
    # rocksalt, cesiumchloride, fluorite or wurtzite.
    ##########################################################################
    #print('crystal_structure:',crystal_structure)
    defect = False
    defect_element = False
    if "_" in crystal_structure:
        #print('yes _',crystal_structure.split("_"))
        defect    = crystal_structure.split("_")[1]  # dilute, vac, selfinterstitial ...
        crystal_structure = crystal_structure.split("_")[0]
        #print('me',matrix_element)
        if defect == 'dilute' and type(matrix_element) != list and len(matrix_element) != 2:
            print('matrix_element:',matrix_element)
            print('type(matrix_element):',type(matrix_element))
            print('len(matrix_element):',len(matrix_element))
            sys.exit('please provide two matrix elements, first = matrix, second is the solute atom type')
        elif defect in [ 'vac', 'selfinterstitial' ] and type(matrix_element) != str:
            print('matrix_element:',matrix_element)
            print('type(matrix_element):',type(matrix_element))
            print('len(matrix_element):',len(matrix_element))
            sys.exit('Error, myutils: please provide only one matrix element')

        if defect == 'dilute':
            defect_element = matrix_element[1]
            matrix_element = matrix_element[0]

    if verbose:
        print('defect',defect)
        print('crystal_structure',crystal_structure)
        print('matrix_element',matrix_element)
        print('defect_element',defect_element)

    if a0 == False:
        state = my_atom.atom([matrix_element]).reference_state[0]
        #print('a0',a0)
        #print('state',state)
        #print('state',state['symmetry'])
        #print('state',state['a'])
        a0 = state['a']
    if crystal_structure == "fcc":
        atom = ase_build_bulk(matrix_element,crystalstructure='fcc',a=a0,cubic=cubic,orthorhombic=orthorhombic)
    elif crystal_structure == "hcp":
        a = 3.21
        c = 5.21
        atom = crystal(matrix_element, [(1./3., 2./3., 3./4.)], spacegroup=194, cellpar=[a, a, c, 90, 90, 120])
    elif crystal_structure == "dc":
        #a = 5.47
        atom = crystal(matrix_element, [(0,0,0)], spacegroup=227, cellpar=[a0, a0, a0, 90, 90, 90])
    elif crystal_structure == "bcc":
        atom = ase_build_bulk(matrix_element,crystalstructure='bcc',a=a0,cubic=cubic,orthorhombic=orthorhombic)

    else:
        sys.exti("crystal_structure: "+str(crystal_structure)+" has to be in fcc/hcp/dc/bcc")

    atomsc = atom.repeat(ncell)
    number_of_atoms = atomsc.get_number_of_atoms()
    nal = number_of_atoms - nsi - nmg

    print('kk')
    print(atomsc.positions)
    if defect != False:
        if defect == "dilute":
            atomsc[0].symbol = defect_element
        elif defect == 'vac':
            del atomsc[0]
        elif defect == 'selfinterstitial':
            atomsc.append(Atom(matrix_element,(a0/2., 0, 0)))

    #print('kki')
    #print(atomsc.positions)
    #sys.exit()

    #for i in np.arange(nmg):
    #    atomsc[i].symbol = 'Mg'
    #for i in np.arange(nmg,nmg+nsi):
    #    atomsc[i].symbol = 'Si'

    # This is the order which ipi kmc expects
    if type(normal_ordering) == bool:
        for i in np.arange(nsi):            atomsc[i].symbol = 'Si'
        for i in np.arange(nsi,nmg+nsi):    atomsc[i].symbol = 'Mg'
        for i in np.arange(1,nvac+1)*-1:    atomsc[i].symbol = 'V'  # for the case: create_fake_vacancy == True

        if create_fake_vacancy != False and create_fake_vacancy != True:
            sys.exit('create_fake_vacancy has to be False or True')

        if create_fake_vacancy == False:
            for i in np.arange(nvac):
                del atomsc[-1]
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
    if type(pos0) != bool:
        for i in np.arange(len(pos0)):
            atomsc[i].position = pos0[i]
    return atomsc

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
        this gives the DFT energy if no calculator is attached
        and the NN/eam energy if a corresponding (lammps) calculator is attached
    '''

    #print('now in ene')
    #print('ac',atoms.cell)
    try:
        # in the case of "DFT"                 , it just retrieves the energy  -> get_stress() can NOT be obtained.
        # in the case of of an ace calculations, it calculates the energy      -> get_stress() CAN     be obtained.
        #print('before')
        ene = atoms.get_potential_energy()
        if False:
            print('ene1:',ene,'atoms.info',atoms.info)
            forces = atoms.get_forces()
            print('pos')
            print(atoms.get_positions())
            print('forces --> this are the NN/EAM forces')
            print(forces)
        #sys.exit()
        #uuid = atoms.info["comment"]
        #print('uuid',uuid)
        #stress = atoms.get_stress()
        #print('stress:',stress)
    except: # RuntimeError:
        print("had runtime error, e.g. ther cant be an energy in the POSCAR")
        ene = 0.
        #stress = False
    if verbose > 1:
        print('ene eV',ene,"(not per atom)")
    units_split = units.split("_")
    einheit = units_split[0].lower()
    #print('us1  :',units_split) #,units_split[1])
    #print('uslen:',len(units_split)) #,units_split[1])
    #print('einheit1:',einheit)
    if len(units_split) > 2:
        print('us',units_split) #,units_split[1])
        sys.exit('len of units should not be > 2!')
    #print('einheit',einheit)
    #print('ene1 (eV)',ene)
    if einheit == 'ev':
        ene_units = ene
    elif einheit == 'mev':
        ene_units = ene*1000.
    elif einheit == "hartree":
        ene_units = ene/aseunits.Hartree
    #print('ene_units (in einheit)',einheit,ene_units)

    if len(units_split) == 2:
        if units_split[1] == 'pa':
            ene_units = ene_units/atoms.get_number_of_atoms()
        else:
            sys.exit("energy can not have this units (ending must be pa, eV_pa or hartree_pa)")
    #print('ene_units (in einheit)',einheit,ene_units)
    return ene_units, ene

def ase_get_chemical_symbols_to_number_of_species(atoms,known_elements_by_pot=[]):
    symbols = atoms.get_chemical_symbols()
    #numat = atoms.get_number_of_atoms()
    #print('symbols',symbols)
    #print('numat',numat)

    uniquesym = set(atoms.get_chemical_symbols())
    d = {}
    for i in uniquesym:
        #print(i,symbols.count(i),numat)
        d[i] = int(symbols.count(i))

    def dcheck(element):
        if element in d.keys():
            #print(element+" exists")
            pass
        else:
            d[element] = 0

    for i in known_elements_by_pot:
        dcheck(i)
    #dcheck("Mg")
    #dcheck("Si")
    #dcheck("Al")
    return d

def ase_get_chemical_symbols_to_conz(atoms,known_elements_by_pot=[]):
    symbols = atoms.get_chemical_symbols()
    numat = atoms.get_number_of_atoms()
    #print('symbols',symbols)
    #print('numat',numat)

    uniquesym = set(atoms.get_chemical_symbols())
    #print("uniquesym",uniquesym) # --> set(['Mg', 'Al'])
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
            d[element] = 0
    #print('dd',d)
    # here, only check for pot.elements
    for i in known_elements_by_pot:
        dcheck(i)
    #dcheck("Mg")
    #dcheck("Si")
    #dcheck("Al")
    return d

def ase_check_chemical_symbols_agree(atoms,ace):
    #print('aa',ace.pot.elements)
    unique_elements_structure = list(set(atoms.get_chemical_symbols()))
    for element in unique_elements_structure:
        #print(element)
        if element not in ace.pot.elements:
            sys.exit('Given structure contains the element '+element+' which is not in the ace potential '+ace.pot.pot)
    return

def ase_remove_unknown_species_from_structure(atoms,ace):
    chemical_symbols = atoms.get_chemical_symbols()
    #print('cs',chemical_symbols)
    unique_elements_structure = list(set(atoms.get_chemical_symbols()))
    #print('aa',unique_elements_structure)
    for idx in range(len(chemical_symbols))[::-1]:
        #print('idx',idx,chemical_symbols[idx])
        if chemical_symbols[idx] not in ace.pot.elements:
            del atoms[idx]
    #print('--1')
    #for idx,cs in enumerate(atoms.get_chemical_symbols()):
    #    print('idx',idx,atoms.get_chemical_symbols()[idx])
    #print('--2')
    return


class ase_calculate_ene( object ):
    '''
    - should be evaluated just once to initialize and NOT FORE EVERY CALCULATION
      (due to the fact that mypot would be called unnecessarily)
    ase_calculate_ene (ace) class which holds lammps commands to be executed
    if only pot is defined, static calculation.
    '''
    def __init__(self,
            potpath_in,
            use_epoch=False,
            units=False,
            geopt=False,
            kmc=False,
            verbose=False,
            temp=False,
            elastic=False,
            fqh_atoms_max=None,
            fqh_supercell=None,
            fah_atoms_max=None,
            fah_supercell=None
            ):

        self.fqh_atoms_max  = fqh_atoms_max
        self.fqh_supercell  = fqh_supercell
        self.fah_atoms_max  = fah_atoms_max
        self.fah_supercell  = fah_supercell
        #self.pot = pot
        self.LAMMPS_COMMAND = False
        self.mypot          = False
        self.units          = units.lower()
        self.geopt          = geopt          # so far only for ene object.
        self.elastic        = elastic
        self.elastic_relax  = True
        self.nsteps         = 0
        self.verbose        = verbose
        self.atoms          = False          # ase atoms object (frame)
        self.lammps         = False          # do the calculation by lammps externally and not by ace.
        if self.verbose:
            print('>> ase_calculate_ene: initializing mypot .... to self.pot')
        self.pot            = mypot(potpath_in,use_epoch = use_epoch,verbose = self.verbose)

        #####################
        # for the calculator
        #####################
        self.calculator = "lammps"
        self.lmpcmd     = False         # in case we run through ase (also needs lmpcmd) or external lammps
        self.atom_types = None          # None for eam-alloy, dict for runner/n2p2
        self.keep_alive = True

        #print('init')
        # case of MD or KMC
        self.kmc  = kmc
        self.ipi = False
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
            print(text,tt+'LAMMPS_COMMAND1) :',self.LAMMPS_COMMAND)    # : n2p2_v2ag
            print()
        return

    def lammps_command_potential_n2p2(self):
        #units_giulio_ene = "0.0367493254"
        ase_units_ene    = "0.03674932247495664" # 1./ase.units.Hartree
        #units_giulio_bohr = "1.8897261328"
        ase_units_bohr    = "1.8897261258369282" # 1./ase.units.Bohr

        command = [
        # showewsum 1 showew yes resetew no maxew 1000000
        'variable nnpDir string \"'+self.pot.potpath_work+'\"',
        "pair_style nnp dir ${nnpDir} showew no resetew yes maxew 100000000 cflength "+ase_units_bohr+" cfenergy "+ase_units_ene,
        "pair_coeff * * "+str(self.pot.potcutoff),
        #"#write_data ./pos.data # would this be the final struct?"
        ]
        return command

    def lammps_command_potential_runner(self):
        command = [
        # comment
        #"# thermo 1 # for geopt",
        'variable nnpDir string \"'+self.pot.potpath_work+'\"',
        "pair_style runner dir ${nnpDir} showewsum 1 showew yes resetew no maxew 1000000",
        #"# pair_coeff * * 7.937658735",
        #"pair_coeff * *  14.937658735"
        "pair_coeff * * "+str(self.pot.potcutoff)
        ]
        return command

    def lammps_command_potential_eam_alloy(self,pair_style="eam/alloy"):
        ''' pair_style eam/alloy or adp '''
        ## This should acutally be the elements of the structure ....
        out = ""
        for idx,i in enumerate(self.pot.elements):
            #print('self.pot.elements['+str(idx)+']:',i)
            out = out +" "+i

        command = [
        "pair_style "+pair_style,
        "pair_coeff * * "+self.pot.potpath+out
        ]

        # bad hack
        #command = [
        #"pair_style eam/alloy",
        #"pair_coeff * * /Users/glensk/Dropbox/Albert/scripts/dotfiles/scripts/potentials/Ni-Al-Mishin-2009.eam.alloy Ni",
        #   "#dump coord2 all xyz 1 ./trajNiH.xyz"
        #]
        return command

    def pot_get_and_ase_lmp_cmd(self,kmc=False,temp=False,nsteps=0,ffsocket='inet',address=False,ipi=False):
        ''' geoopt (geometry optimization) is added / or not in
            lammps_write_inputfile(); here only the potential is set.
            ffsocket: ipi ffsocket [ "unix" or "inet" ]
        '''
        if self.verbose > 1:
            print('PPh potDONE:',self.pot.potDONE)

        ###########################################################
        # here the magic happens where the potential is defined
        ###########################################################
        if self.pot.potDONE == False:
            if self.verbose:
                print("PPi self.pot.get()")
            self.pot.get()

        self.kmc = kmc
        self.temp = temp
        self.nsteps = nsteps
        self.ffsocket = ffsocket
        if self.ffsocket not in [ "unix", "inet" ]:
            print('ffsocket:',ffsocket)
            sys.exit('ffsocket has to be "unix" or "inet"; Exit!')



        if self.verbose > 2:
            tt = 'PPj pot_get_and_ase_lmp_cmd_A '
            print(tt+'pot.pot           :',self.pot.pot)
            print(tt+'pot.potpath       :',self.pot.potpath)
            print(tt+'pot.potpath_work  :',self.pot.potpath_work)
            print(tt+'pot.pottype       :',self.pot.pottype)
            print()

        #sys.exit()
        # this depends only on the potential which is already defined
        # so should be easy to make this general.
        self.lmpcmd = [ "### lmpcmd.begin ###" ]
        if self.pot.pottype == "n2p2" or self.pot.pottype == "runner":
            self.lmpcmd = self.lmpcmd + self.pot.atom_masses_str #self.lammps_command_masses()
        else:
            self.lmpcmd = self.lmpcmd

        if self.pot.pottype == "n2p2":
            # showewsum 1 showew yes resetew no maxew 1000000
            self.lmpcmd = self.lmpcmd + self.lammps_command_potential_n2p2()
            #self.atom_types = {'Mg':1,'Al':2,'Si':3}
            self.atom_types = self.pot.atom_types

        elif self.pot.pottype == "runner":
            self.lmpcmd = self.lmpcmd + self.lammps_command_potential_runner()
            #self.atom_types = {'Mg':1,'Al':2,'Si':3}
            self.atom_types = self.pot.atom_types

        elif self.pot.pottype == 'eam-alloy':
            self.lmpcmd = self.lmpcmd + self.lammps_command_potential_eam_alloy()
            self.atom_types = None
            #for idx,i in enumerate(self.pot.elements):
            #    #print('i',i)
            #    self.atom_types[i] = idx
        elif self.pot.pottype == 'adp':
            self.atom_types = None
            self.lmpcmd = self.lmpcmd + self.lammps_command_potential_eam_alloy(pair_style='adp')
        else:
            print('self.pot.pottype:',self.pot.pottype)
            print('self.pot.pot    :',self.pot.pot)
            sys.exit('pot '+str(self.pot.pot)+' not found! (X)')
        #print('self.lm (1) ',self.lmpcmd)
        #elements = np.unique(atoms.get_chemical_symbols())

        #self.lmpcmd = self.lmpcmd + [ 'ka ka test\n' ]
        #print('self.lm (2) ',self.lmpcmd)
        if self.kmc:
            self.lmpcmd = self.lmpcmd + [
                "",
                "timestep 0.001   # timestep (ps)",
                "velocity all create "+str(self.temp)+" 4928459",  # create initial velocities 4928459 is random seed for velocity initialization"
                "thermo 1   # screen output interval (timesteps)"
                ]
                # with n2p2 in the parallel version, inet is not working
                # "fix 1 all ipi fidis 12345",     # for fidis job
                # "fix 1 all ipi mac 77776 unix",  # for mac job

        if self.ffsocket == "unix": add = "unix"
        if self.ffsocket == "inet": add = ""
        if address == False:
            address = gethostname()

        if self.kmc or ipi:
            self.lmpcmd = self.lmpcmd + [ "", "## in case run by ipi: ",
            "## fix 1 all ipi "+str(address)+" 12345 "+str(add), "" ]

        self.lmpcmd = self.lmpcmd + [ "### lmpcmd.end  ####" ]
        if self.verbose > 1:
            print('HERE THE lmpcmd I got',self.lmpcmd)
        self.print_variables_ase("pot_get_and_ase_lmp_cmd_FIN")
        self.LAMMPS_COMMAND = get_LAMMPS_executable(exit=True) #,verbose=self.verbose)
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
        ene,ene_ev_tot = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        return ene

    def ene_new(self,atoms=False,
            atomrelax=False,
            volumerelax=False,
            cellshaperelax=False,
            print_minimization_to_screen=False,minimizer="LGBFGS",
            return_both=False):
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
            self.get_murn(atoms,verbose=False,return_minimum_volume_frame=True)
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
        ene,ene_ev_tot = ase_enepot(atoms,units=self.units,verbose=self.verbose)
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
        if return_both == True:
            return ene,ene_ev_tot
        else:
            return ene

    def ene(self,atoms=False,
            atomrelax=False,
            cellrelax=False,
            cellshaperelax=False,
            cellvolumerelax=False,
            print_minimization_to_screen=False,
            minimizer="LGBFGS",
            debug=False,
            return_both=False):
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
                #logfile=None
                #logfile="-"


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
            print('ZZ done236')
            print('ZZ self.units',self.units)
        ######################################################
        # calculate the energy
        ######################################################
        #print('atxxx',atoms)
        ene1 = atoms.get_potential_energy()
        if self.verbose > 1:
            print('ZZ done237',ene1)
        ene,ene_ev_tot = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        if self.verbose > 1:
            print('ZZ done238')
        #print('self.units:',self.units,'ene pot',ene,'ene_ev_tot',ene_ev_tot)
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
        if return_both == True:
            return ene,ene_ev_tot
        else:
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
        #print('######## get_elastic_external #############')
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
            print('frame nat ::',frame.get_number_of_atoms())
            print('verbose   ::',verbose)
        if self.elastic_relax == True:
            self.ase_relax_cellshape_and_volume_only(frame,verbose=verbose)
        #print('stress !relaxed! frame :',frame.get_stress())
        #print('volume !relaxed! frame :',frame.get_volume())
        #print('ene    !relaxed! frame :',frame.get_potential_energy())
        print('### get_elastic_external: !RELAXED! stress(max),vol,ene',abs(frame.get_stress()).max(),frame.get_volume(),frame.get_potential_energy())
        #print('frame cell',frame.get_cell())
        #ase_write('pos.runner',frame,format='runner')
        #print('-------------- lammps_ext_calc -----------')
        ene_ev_cell, ene_pot_lmp = lammps_ext_calc(frame,self,get_elastic_constants=get_all_constants)
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
        ''' here it would be good to have the pot values... '''
        #print('hhhhhhhhhhhh',self.pot.pottype)
        #print('hhhhhhhhhhhh in ',self.lmpcmd)
        #if self.pot.pottype in [ 'eam', 'eam-alloy' ]:
        #    elements = np.unique(atoms.get_chemical_symbols())
        #    if len(elements) != 1:
        #        sys.exit('Error: you have more than one element in you simulation box. I am not sure how to define the eam-alloy potential. using e.g. Al Ni will give different results to Ni Al; in the docs they propose to define all the atoms which does not work in my case')
        #    else:
        #        out = " "+elements[0]

        #    self.atom_types = None
        #    ### This should acutally be the elements of the structure ....
        #    #out = ""
        #    #for idx,i in enumerate(atoms.get_chemical_symbols()):
        #    #    print('self.pot.elements['+str(idx)+']:',i)
        #    #    out = out +" "+i
        #    #    out = " Al"
        #    #    out = " Ni"
        #    #    out = " Al Ni"
        #    #    out = " Ni Al"

        #    self.lmpcmd = [
        #    "pair_style eam/alloy",
        #    "pair_coeff * * "+self.pot.potpath+out
        #    ]
        #print('hhhhhhhhhhhh out ',self.lmpcmd)
        #sys.exit()
        #print('hhhhhhhhhhhh----',self.lmpcmd)

        if self.verbose > 1:
            for i in self.lmpcmd:
                print("lmpcmds    :",i)
            print("get_calculator: atom_types :",self.atom_types)
            print("get_calculator: keep_alive :",self.keep_alive)
            print("structure/atoms types      :",atoms.get_chemical_symbols())
        if self.calculator == "lammps":
            if self.verbose > 1:
                print('get_calculator (Y1): lmpcmds = self.lmpcmd    :',self.lmpcmd)
                print('get_calculator (Y2): atom_types=self.atom_types:',self.atom_types)
            asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=self.keep_alive)
            if self.verbose > 1:
                print('now asecalcLAMMPS defined')
            atoms.set_calculator(asecalcLAMMPS)
            if self.verbose > 1:
                print('now with calculator')
        return

    def ase_relax_atomic_positions_only(self,atoms,fmax=0.0001,verbose=False,output_to_screen=False):
        ''' The strain filter is for optimizing the unit cell while keeping scaled positions fixed. '''
        self.keep_alive = True
        self.get_calculator(atoms)
        #print('sssss')
        if False: #verbose:
            print('1: relax atomic positions; stress:',atoms.get_stress(),"volume per atom:",ase_vpa(atoms))

        logfile="-" # output to screen
        logfile="tmpxxase" # output to file and not to screen
        if output_to_screen == True:
            logfile = '-'
        #opt = LBFGS(atoms,logfile=logfile)
        opt = FIRE(atoms,logfile=logfile,dt = 0.1)
        #opt = BFGS(atoms,logfile=logfile)
        opt.run(fmax=fmax)
        if os.path.isfile(logfile):
            os.remove(logfile)
        if verbose:
            print('2: relax atomic positions; stress:',atoms.get_stress(),"volume per atom:",ase_vpa(atoms))
        return


    def ase_relax_volume_only_from_murn_roughly(self,atoms,verbose=False):
        atoms.wrap()
        # probably not necessary since this is set up when necessary
        self.keep_alive = True
        self.get_calculator(atoms)

        alat_rel = np.arange(0.9865,1.015,0.0025)  # seems most symmetric for cu
        alat_rel = np.array([0.9865,1-0.007,1,1+0.007,1.015]) # 4 values enough for quick estimate
        cell_ref = atoms.get_cell()
        nat = atoms.get_number_of_atoms()  # in case stuff is repeated
        cell_ref_vol_pa = atoms.get_volume()/nat
        vol_pa = np.zeros(len(alat_rel))
        ene_pa = np.zeros(len(alat_rel))
        for idx,i in enumerate(alat_rel):
            atoms_murn_loop = atoms.copy()  # in case cell is repeated
            atoms_murn_loop.set_cell(cell_ref*i,scale_atoms=True)
            vol=atoms_murn_loop.get_volume()
            vol_per_atom = vol/atoms_murn_loop.get_number_of_atoms()
            self.get_calculator(atoms_murn_loop)
            percent = np.round(vol_per_atom/cell_ref_vol_pa,2)
            ene = atoms_murn_loop.get_potential_energy()
            if self.lammps == True:
                ene, ene_pa = lammps_ext_calc(atoms_murn_loop,self)
            print('idx!!',idx,'i:',np.round(i,2),np.round(vol,2),'vpa',np.round(vol_per_atom,2),'ref',cell_ref_vol_pa,'v/vo',percent,'ene_pa',ene/nat,'ene',ene)
            vol_pa[idx] = vol/nat
            ene_pa[idx] = ene/nat

        from feos import eos
        vinet = eos()
        data=np.transpose([vol_pa,ene_pa])
        vinet.fit_to_energy_vs_volume_data(datax=vol_pa,datay=ene_pa)
        print('pars',vinet.parameters)
        if True: # scale volume
            print('adapt atoms in')
            atomsin = atoms.copy()
            # old
            volume_in  = ase_vpa(atomsin)
            volume_out = vinet.parameters[1]
            #if type(return_frame_with_volume_per_atom) != bool:
            #    volume_out = return_frame_with_volume_per_atom
            volume_scale = (volume_out/volume_in)**(1./3.)
            print('volume_in   (calc)',volume_in)
            print('volume_out  (calc)',volume_out)
            ##print('scale',volume_scale)
            atoms.set_cell(atomsin.get_cell()*volume_scale,scale_atoms=True)
            print('atomsout vol:',atoms.get_volume()/atoms.get_number_of_atoms())
            if volume_in/volume_out > 1.05 or volume_in/volume_out < 0.96:
                self.ase_relax_volume_only_from_murn_roughly(atoms,verbose=self.verbose)
        return

    def ase_relax_volume_only(self,atoms,verbose=False):
        self.keep_alive = True
        self.get_calculator(atoms)
        #if False: # this works but is quite slow! better do manually
        sf = UnitCellFilter(atoms,hydrostatic_strain=True)
        logfile="-" # output to screen
        #logfile="tmp" # output to file and not to screen
        opt = FIRE(sf,logfile=logfile)
        #from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
        #opt = SciPyFminCG(sf,logfile=logfile)
        opt.run()
        if os.path.isfile("tmp"):
            os.remove("tmp")
        if verbose > 1: self.check_frame('ase_relax_volume_only out',frame=atoms)
        return

    def ase_relax_cellshape_and_volume_only(self,atoms,verbose=False):
        ''' The strain filter is for optimizing the unit cell while keeping scaled positions fixed. '''
        self.keep_alive = True
        #print('1aaaaaaaaaaaaa')
        #print('self.lmpcmd: ddxx',self.lmpcmd)
        self.get_calculator(atoms)
        #print('2bbbbbbbbbbbbbb')
        #print(atoms.positions)
        #print(atoms.get_chemical_symbols())
        #enRef = atoms.get_potential_energy()
        #print('3bbbbbbbbbbbbbb enRef',enRef)
        if verbose > 1: self.check_frame('ase_relax_cellshape_and_volume_only in',frame=atoms,verbose=verbose)
        #print('3cccccccccccccccc')

        sf = StrainFilter(atoms)
        logfile="-" # output to screen
        logfile="tmp" # output to file and not to screen
        opt = BFGS(sf,logfile=logfile)
        opt.run(0.005)
        if os.path.isfile("tmp"):
            os.remove("tmp")
        if verbose > 1: self.check_frame('ase_relax_cellshape_and_volume_only out',frame=atoms)
        return

    def ase_relax_cellshape_volume_positions(self,atoms,verbose=False):
        print('relaxing cellshape, volume, positions .... this may take a while.')
        self.get_calculator(atoms)
        delta_pos_max = 0.00001
        stress_max = 0.00001
        #delta_pos_max = 0.1
        #stress_max = 0.1

        i = 0
        stressmax = 10
        deltaposmax = 10
        while stressmax > stress_max or deltaposmax > delta_pos_max:
            i+=1
            pos_in = atoms.copy().positions
            #print('px0 p',atoms.positions[1])
            self.ase_relax_atomic_positions_only(atoms,fmax=0.0001,verbose=False,output_to_screen=False)
            #print('px1 p',atoms.positions[1])
            self.ase_relax_cellshape_and_volume_only(atoms,verbose=False)
            #print('px2 p',atoms.positions[1])
            pos_out = atoms.positions

            if verbose:
                print(i,'px?in  p',pos_in[1])
                print(i,'px?out p',pos_out[1])


            deltaposmax = np.abs(pos_in - pos_out).max()
            if verbose:
                print(i,'deltaposmax',deltaposmax)
            stressmax = np.abs(atoms.get_stress()).max()
            if verbose:
                print(i,i, 'stressmax',stressmax)
        return



    def get_murn(self,atomsin=False,verbose=False,
            return_minimum_volume_frame=False,
            return_frame_with_volume_per_atom=False,
            atomrelax=False,
            write_evinet=False,
            write_fqh=False,
            get_to_minvol_first=True,
            printminimal=True):
        ''' the murn will never change the atomsobject
        return_frame_with_volume_per_atom : volume can be specified and he frame scaled
        '''
        if write_evinet and os.path.isdir('evinet'):
            sys.exit('evinet folder already exists')
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

        if write_evinet == False and write_fqh == False:
            return

        ### only get here if either write_evinet or write_fqh
        # range for eos!
        # the range for volume expansion until melting has to be larger!
        if write_evinet == True:
            alat_rel=[0.97,0.975,0.98,0.985,0.99,0.995,0.998,1.0,1.002,1.005,1.01,1.015,1.02,1.025,1.03]
            alat_rel = np.arange(0.97,1.03,0.005)
            alat_rel = np.arange(0.985,1.015,0.0025)
            alat_rel = np.arange(0.99,1.015,0.0025)  # otherwise repulsive part too strong
            alat_rel = np.arange(0.9865,1.015,0.0025)  # seems most symmetric for cu
        if write_evinet == False and write_fqh != False:
            print('getting volume range for Fqh, V0 (T=0K) - V0*1.1 (T=Tmelt) seems good for Cu and Al')
            alat_rel = np.arange(0.965,1.036,0.005)
            print('alat_rel',alat_rel)

        vol_pa = np.zeros(len(alat_rel))
        ene_pa = np.zeros(len(alat_rel))

        cell_ref = atoms_murn.get_cell()
        cell_ref_vol_pa = atoms_murn.get_volume()/atoms_murn.get_number_of_atoms()
        nat = atoms_murn.get_number_of_atoms()

        atoms_murn_loop = atoms_murn.copy()
        if verbose:
            print('               vol [Ang^2/at] ene [??]:')
        print('nat:',nat)
        print('vol cell:',atoms_murn.get_volume())
        print('vol cell/pa:',atoms_murn.get_volume()/nat)

        for idx,i in enumerate(alat_rel):
            if verbose > 2:
                print('000 idx:',idx,'i:',i)
            atoms_murn_loop = atoms_murn.copy()  # in case cell is repeated
            atoms_murn_loop.set_cell(cell_ref*i,scale_atoms=True)

            print('atoms_murn_loop.cell (1)',atoms_murn_loop.cell)
            print('atoms_murn_loop.pos  (1)',atoms_murn_loop.positions)
            print('atoms_murn_loop nat  (1)',atoms_murn_loop.get_number_of_atoms())
            if False: # this all seems to work
                atoms_murn_loop = ase_repeat_cell(atoms_murn_loop,2,2,2)
                atoms_murn_loop.write(self.pot.lammps_tmpdir+'pos.lmp',format='lammps-data')
                ase_write(self.pot.lammps_tmpdir+'pos.lmp2',atoms_murn_loop,format='lammps-data')
                ase_write(self.pot.lammps_tmpdir+'pos.POSCAR',atoms_murn_loop,format='vasp')
                ase_write(self.pot.lammps_tmpdir+'pos.extxyz',atoms_murn_loop,format='extxyz')
                ang_to_bohr = 1.8897261
                np.savetxt(self.pot.lammps_tmpdir+"/x_reference.data",atoms_murn_loop.positions.flatten()*ang_to_bohr,newline=" ",fmt="%3.15f")
                ase_write(self.pot.lammps_tmpdir+"/pos.ipi.xyz",atoms_murn_loop,format='ipi')
            print('atoms_murn_loop.cell (2)',atoms_murn_loop.cell)
            print('atoms_murn_loop nat  (2)',atoms_murn_loop.get_number_of_atoms())
            #print('atoms_murn_loop.pos  (2)',atoms_murn_loop.positions)
            nat = atoms_murn_loop.get_number_of_atoms()  # in case stuff is repeated

            if verbose > 2:
                print('111 cell',atoms_murn_loop.get_cell())
            #print('111 cell',atoms_murn_loop.get_cell())
            #print('111 cell',atoms_murn_loop.get_cell()[0,0]/atoms_murn_loop.get_cell()[1,1])
            vol=atoms_murn_loop.get_volume()
            vol_per_atom = vol/atoms_murn_loop.get_number_of_atoms()
            percent = np.round(vol_per_atom/cell_ref_vol_pa,2)
            if verbose > 2:
                print('222 vol',atoms_murn_loop.get_volume())
            #ene = self.ene(atoms_murn_loop)                         # works
            #ene = self.ene_new(atoms_murn_loop)                         # works

            ########################################
            # for vacancies, self interstitials we need a geometry relaxation here
            # for highly symmetric structures this will not affect the result
            ########################################
            #print('volume',atoms_murn_loop.get_volume())
            #print('p1    :',atoms_murn_loop.positions[0])
            self.ase_relax_atomic_positions_only(atoms_murn_loop,fmax=0.0001,verbose=False,output_to_screen=False)
            #print('volume',atoms_murn_loop.get_volume())
            #print('p1    :',atoms_murn_loop.positions[0])
            #sys.exit()
            ene = self.ene_allfix(atoms_murn_loop)                       # works
            ene_ev = atoms_murn_loop.get_potential_energy()
            ene = atoms_murn_loop.get_potential_energy()
            print('ene_ev ase',ene_ev)
            print('ene    ase',ene)
            print('self.lammps',self.lammps)
            if self.lammps == True:
                ene_ev_cell, ene = lammps_ext_calc(atoms_murn_loop,self)
                ene_ev = ene_ev_cell
                ene = ene_ev_cell
                print('ene_ev_cell (lmpammps external):',ene_ev_cell)
                #print('ene    lmp',ene)
            print('idx!!',idx,np.round(i,2),np.round(vol,2),'vpa',np.round(vol_per_atom,2),'ref',cell_ref_vol_pa,'v/vo',percent,'ene_pa',ene/nat)
            #if printminimal == True:
            if verbose:
                print('ene ase tot (eV):',str(round(ene_ev,7)).ljust(10),'nat:'+str(nat),"vol/pa:",str(round(vol/nat,7)).ljust(10),'(eV/at):',str(ene/nat).ljust(19)) #self.units)
            #print('ams3',atoms_murn_loop.get_stress(),ase_vpa(atoms_murn_loop))
            if verbose > 2:
                stress = atoms_murn_loop.get_stress()[:3]
                cell = atoms_murn_loop.get_cell()
                pos = atoms_murn_loop.get_positions()[1]/cell[0,0]

                if cell_is_cubic == True:
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
            if write_fqh and alat_rel[idx] >= 0.999:
                #print('write_fqh: idx:',idx,'i:',i,'alat_rel:',alat_rel[idx])
                print('now phonopy (not yet)')
                #self.get_fh_phonopy(atomsin=atoms_murn_loop)
                #sys.exit('phonopy 77 done')
                self.get_fh_onevolume(atomsin=atoms_murn_loop,disp=0.03,debug=False,try_readfile=False,atomrelax=True,write_fqh=write_fqh)
            if verbose > 2:
                print('idx:',str(idx).ljust(3),'i:',str(i).ljust(10),'vol:',str(vol).ljust(10),'ene:',ene)
            if verbose:
                stress = 0
                if self.lammps == True:  # external calculation
                    stress = 0
                else:
                    stress = atoms_murn_loop.get_stress()[:3]
                print('i',str(i).ljust(5),'vol/nat',str(round(vol/nat,7)).ljust(10),'ene/nat',str(ene/nat).ljust(19),stress)


        if verbose: self.check_frame('get_murn 3 atfer loop           ',frame=atoms_murn)
        if verbose > 1:
            stress = atoms_murn.get_stress()[:3]
            cell = atoms_murn.get_cell()
            pos = atoms_murn.get_positions()[1]/cell[0,0]
            print("---1-->>",pos)
        if write_evinet:
            hier = os.getcwd()
            #print('hier1',os.getcwd())
            if not os.path.isdir('evinet'):
                os.mkdir("evinet")
            os.chdir("evinet")
            print('vol_pa',vol_pa)
            print('ene_pa',ene_pa)

            write_energies_atom = "energy_eV_per_atom_vs_ang3_per_atom.dat"
            write_energies_cell = "energy_eV_per_cell_vs_ang3_per_cell.dat"
            np.savetxt(write_energies_atom,np.transpose([vol_pa,ene_pa]))
            np.savetxt(write_energies_cell,np.transpose([vol_pa*nat,ene_pa*nat]))
            os.chdir(hier)
        if verbose > 1:
            print('loop done')

        ########################
        # feos
        ########################
        from feos import eos
        vinet = eos()
        data=np.transpose([vol_pa,ene_pa])
        if printminimal == True:
            print('vol_pa',vol_pa)
            print('ene_pa',ene_pa)
        for idx,i in enumerate(vol_pa):
            print(idx,vol_pa[idx],ene_pa[idx])

        vinet.fit_to_energy_vs_volume_data(datax=vol_pa,datay=ene_pa)
        #self.eos = vinet.parameters
        if True: #verbose > 1:
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
        #if type(write_evinet) == bool:
        if write_evinet == True: # this shoul never (over) write the Evinet file in case write_evinet == False
            hier = os.getcwd()
            os.chdir("evinet")
            vinet.write_data()
            print('vinet.parameters:',vinet.parameters)
            os.chdir(hier)

        return

    def get_fh_phonopy(self,atomsin=False):
        ''' the function will never change the atomsobject '''
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        atoms_h = atomsin.copy()
        atoms_h.wrap()
        nat = atoms_h.get_number_of_atoms()
        print('nat',nat)
        if nat < 20:
            atoms_h *= (2,2,2)
        nat = atoms_h.get_number_of_atoms()
        print('nat',nat)

        self.keep_alive = False
        self.get_calculator(atoms_h)

        phonopy_atoms = ase_to_phonopy(atoms_h)
        phonon = phonopy_pre_process(phonopy_atoms, supercell_matrix=np.eye(3, dtype='intc'))
        supercells = phonon.get_supercells_with_displacements()

        # Force calculations by calculator
        set_of_forces = []
        for scell in supercells:
            cell = phonopy_to_ase(scell)
            #cell = Atoms(symbols=scell.get_chemical_symbols(),
            #             scaled_positions=scell.get_scaled_positions(),
            #             cell=scell.get_cell(),
            #             pbc=True)
            #cell.set_calculator(cell)
            self.get_calculator(atoms_h)
            forces = atoms_h.get_forces()
            drift_force = forces.sum(axis=0)
            print(("[Phonopy] Drift force:" + "%11.5f" * 3) % tuple(drift_force))
            # Simple translational invariance
            for force in forces:
                force -= drift_force / forces.shape[0]
            set_of_forces.append(forces)
        phonopy_post_process(phonon, set_of_forces)

        return

    def get_fh_onevolume(self,atomsin=False,disp=0.03,debug=False,try_readfile=False,atomrelax=True,write_fqh=False):
        ''' the function will never change the atomsobject
        atomrelax: in most cases it is desirable to relax the atomic positions
        in few cases we might be tempted to assess the free energy for particular positions (e.g. to check weather the DFT equilibrium position is the stable position of a NN)
        '''
        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        atoms_h = atomsin.copy()
        atoms_h.wrap()
        #sys.exit('79')
        T0shift_ev_atom = self.ene(atoms_h)/atoms_h.get_number_of_atoms()
        #sys.exit('80')
        #print('--> T0shift_ev_atom',T0shift_ev_atom)
        #sys.exit('T0shift_ev_atom')

        #if return_units == "mev_pa":
        #    return_mult = 1.
        #elif return_units == "ev_cell":
        #    return_mult = atomsin.get_number_of_atoms()/1000.
        #else:
        #    raise NameError('return_units conversion not yet set for '+return_units)

        #sys.exit('81')
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


        #sys.exit('80')
        #print('xx80.1')
        self.keep_alive = False
        self.get_calculator(atoms_h)
        #sys.exit('81')
        #print('xx80.2')

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
            #print('harmonic cell stress',atoms_h.get_stress())
        maxforce = np.abs(atoms_h.get_forces()).max()
        #print('xx80.3 maxforce',maxforce)
        if debug:
            print("###########################################")
            #print("forces harmonic 4:",atoms_h.get_forces()[:3])
            #print("###########################################")
            print('out atoms_h str',atoms_h.get_stress())
            print('out maxforce out',maxforce)

        #print('xx80.4')
        if atomrelax == True and maxforce > 0.0001:
            print("forces harmonic 4:",atoms_h.get_forces()[:3])
            print('maxforce',maxforce)
            sys.exit('maxforce is too large')

        #print('xx80.5')
        if maxforce > 0.0001:
            # from this it needs to be deduced that NEGATIVE EIGENVALUES
            # calculating the hesse makes a segmentation fault with ace lammps
            return False
        if debug:
            print("###########################################")
            print('atoms_h.get_cell() before',atoms_h.get_cell())

        #print('xx80.6')
        #################################
        # decide abou supercell size
        # a) defaut cellsize, or take the defined one
        #################################
        nat = atoms_h.get_number_of_atoms()
        print('self.fqh_atoms_max (1):',self.fqh_atoms_max)
        print('self.supercell     (1):',self.fqh_supercell)
        print('number of atoms    (1):',nat)
        if self.fqh_supercell is not None:
            rep = self.fqh_supercell
            print('rep from fqh_supercell',rep)
        elif self.fqh_atoms_max is not None:
            rep = 1
            for i in np.arange(2,11):
                at_ = nat*i**3
                print('i',i,'at_',at_,'rep',rep)
                if at_ <= self.fqh_atoms_max:
                    rep = i
                elif at_ > self.fqh_atoms_max:
                    break
            print('rep from fqh_atoms_max',rep)
        else:
            rep = 1
            print('rep not defined (therefore: rep = 1)')


        atoms_h *= (rep,rep,rep)
        nat = atoms_h.get_number_of_atoms()
        print('rep',rep,'nat',nat)
        #ase_write("POSCAR_ka",atoms_h,format='vasp')
        #sys.exit()



        if debug:
            nat = atoms_h.get_number_of_atoms()
            print("###########################################")
            print('!!!!!!! nat:',nat)
            print("forces harmonic 2:",atoms_h.get_forces()[:3])
            print("stress:",atoms_h.get_stress())
            print("forces max harmonic 2:",abs(atoms_h.get_forces()).max())
            print('atoms_h.get_cell() after mult',atoms_h.get_cell())
            print("###########################################")
        #print('xx80.8')
        pos0 = atoms_h.get_positions()
        #print('xx80.9')
        hessematrix=np.zeros((pos0.shape[0]*3,pos0.shape[0]*3))
        print('xx81.0',nat)

        ### schleife ueber alle atome, 1..32

        for iidx,i in enumerate(pos0): # loop over all atoms
            progress(iidx,len(pos0),status=try_readfile) #try_readfile)
            #print(iidx,"/",pos0.shape[0]) # loop over xyz 1..3
            for jidx,j in enumerate(i):
                pos1 = np.copy(pos0)
                pos1[iidx,jidx] = pos1[iidx,jidx]+disp
                #print('iidx (a)',iidx,"/",pos0.shape[0],'jidx',jidx)
                atoms_h.set_positions(pos1)
                #print('iidx (b)',iidx,"/",pos0.shape[0],'jidx',jidx)
                #print(atoms_h.cell)
                #print(atoms_h.positions)
                self.keep_alive = False
                self.get_calculator(atoms_h)
                #print('iidx (c)',iidx,"/",pos0.shape[0],'jidx',jidx)
                fah = atoms_h.get_forces()
                #print('iidx (d)',iidx,"/",pos0.shape[0],'jidx',jidx)
                hessematrix[iidx*3+jidx] = fah.reshape((1,pos0.shape[0]*        3))/(-disp)

        print('xx81.1',nat)
        #np.savetxt("HesseMatrix.dat",hessematrix/97.173617,fmt="%.13f")
        if debug:
            print("get free energy ...")
            print('get_chemical_symbols()',atoms_h.get_chemical_symbols())
            print()

        Tmax_elements = list(set(atoms_h.get_atomic_numbers()))
        Tmax = my_atom.get_meltingpoint_rounded(Tmax_elements,maximum=True)
        hes = hesse.hesseclass(listin=atoms_h.get_chemical_symbols(),H=hessematrix,show_negative_eigenvalues = True, Tmax=Tmax, T0shift_ev_atom = T0shift_ev_atom)
        print('eigenvalues',hes.freqs)
        #try:
        #    #free_ene = (hes.ene_atom[300]-hes.ene_atom[0])[1]
        #    free_ene      = hes.ene_atom
        #    free_ene_cell = hes.ene_cell
        #except IndexError:
        #    free_ene = "UNSTABLE"
        #if debug and type(free_ene) != str:
        #    hes.write_ene_atom()
        #if try_readfile:
        if write_fqh:
            #try_readfile = "Ni"
            vol = atoms_h.get_volume()/atoms_h.get_number_of_atoms()
            volstr = "_"+str(vol)
            folder = "Fqh_"+str(rep)+"x"+str(rep)+"x"+str(rep)
            if True: #hes.has_negative_eigenvalues == False:
                # write in any case, if it has negative eigenvalues so be it
                #hes.write_hessematrix(try_readfile+"_hessematrix"+volstr)
                print('os.getcwd xx',os.getcwd())
                print('write_fqh',write_fqh)
                self.fqh_folder = "fqh"
                print('self.fqh_folder',self.fqh_folder)  # folder+fqh
                if not os.path.isdir(self.fqh_folder):
                    os.makedirs(self.fqh_folder)
                # here we are now in the fqh folder
                hes.write_hessematrix(self.fqh_folder+"/Hessematrix"+volstr)
                if hes.has_negative_eigenvalues == False:
                    # only write for ground state
                    subfolder = self.fqh_folder+"/Fqh_"+str(nat)+"at_cell_per_"
                    hes.write_ene_atom(subfolder+"atom"+volstr)
                    hes.write_ene_cell(subfolder+"cell"+volstr)
                    ase_write(self.fqh_folder+"/positions"+volstr+".extxyz",atoms_h,format='extxyz')
                    ase_write(self.fqh_folder+"/positions"+volstr+".POSCAR",atoms_h,format='vasp')
            #self.fqh_files.append([volstr,])
            #np.savetxt(try_readfile,free_ene)

        #print('k T)shift   ',T0shift_ev_atom)
        #print('hes.ene_cell',hes.ene_cell)
        #print('hes.ene_atom',hes.ene_atom)
        #get = np.loadtxt(try_readfile)
        #return np.transpose([get[:,0],get[:,1]*return_mult])
        #return free_ene*return_mult
        return hes

    def check_frame(self,text,frame=False,verbose=True,setupcalc=True):
        if verbose:
            print('YYX',text)
        if type(text) != str:
            raise TypeError("need a text!")
        #print('b')


        if setupcalc == True:
            #print('c',verbose)
            self.keep_alive = True
            self.get_calculator(frame)
        #print('d',verbose)

        if verbose:
            print('check_frame: (YY)')
            print('check_frame:',frame.get_stress())
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

def ase_repeat_cell(atoms,r1,r2,r3):
    ''' this should never change the original ase object'''
    atoms_sc = atoms.copy()
    atoms_sc *= (r1,r2,r3)
    return atoms_sc

def ase_vpa(atoms):
    return atoms.get_volume()/atoms.get_number_of_atoms()
def ase_epa(atoms):
    return atoms.get_potential_energy()/atoms.get_number_of_atoms()
def ase_mepa(atoms):
    return atoms.get_potential_energy()/atoms.get_number_of_atoms()*1000.
def ase_fmax(atoms):
    return abs(atoms.get_forces()).max()

def ase_relax_structure_fully_and_save_to_pot(ace,read=False,read_fullpath=False,whichstruct=":"):
    if read != False:
        pathbase = os.environ['potentials']+"/aiida_get_structures_new/"
        readpath = pathbase + read
    #print('ace',ace.pot.potpath)
    #sys.exit('77')
    if read_fullpath != False:
        read = os.path.basename(read_fullpath)
        readpath = read_fullpath

    if not os.path.isfile(readpath):
        sys.exit('readpath '+readpath+" does NOT exist (55); Exit")

    #print('read from',readpath)
    #print('whichstruct',whichstruct)
    print('readpath',readpath,'whichstruct:',whichstruct)
    frames = ase_read(readpath,index=whichstruct,format="runner")  # frames with DFT energies.
    #print('type(frames)',type(frames))
    #print('frames.cell',frames.cell)
    if type(frames) != list:
        frames = [frames]
    #print('whichstruct',whichstruct)
    #print('frame len',len(frames))
    #print('ace.pot',ace.pot.pot)
    #print('ace.pot',ace.pot.potpath)
    #print('ace.pot',ace.pot.use_epoch)
    savefolder = ace.pot.potpath+"/epoch_"+str(ace.pot.potepoch_using)
    if not os.path.isdir(savefolder):
        os.mkdir(savefolder)
    savepath = savefolder+"/RELAXED_fully_"+read

    ######################################
    # just read in and return if it exists
    ######################################
    if os.path.isfile(savepath):
        print('savepath (exists)',savepath)
        frames_out = ase_read(savepath,":",format="runner")  # frames with DFT energies.
        if len(frames) == len(frames_out):
            return frames_out, frames
        else:
            print('savepath exists but has different name than readpath')
            sys.exit('savepath '+savepath+" does already exist (55); Exit")

    print('savepath (will be created)',savepath)
    ######################################
    # relax if necessary
    ######################################
    for idx,atoms in enumerate(frames):
        #print('idx',idx,'out of',len(frames))
        #print('aaa')
        #print(atoms.get_cell())
        #print(atoms.positions)
        #print('c')
        ace.ase_relax_cellshape_volume_positions(atoms,verbose=True)
        ase_write(savepath,atoms,format='runner',append=True)
    frames_out = ase_read(savepath,":",format="runner")  # frames with DFT energies.
    print('savepath (created)',savepath)
    return frames_out, frames


def get_thermo(ace,atoms,relax_cellshape_and_volume=True,evinet=True,fqh=False,fah=False):
    ''' only one of evinet, fqh, fah can be true at one time '''

    #################################
    ## show what we have to deal with
    #################################
    print('atoms.cell (angstrom) before relaxing cellshape and volume:')
    print(atoms.cell)
    print('atoms.volume/nat before relaxation:',atoms.get_volume()/atoms.get_number_of_atoms())
    print('atoms.positions before relaxing cellshape and volume (first 4) in angstrom:')
    #for idx,i in enumerate(atoms.positions):
    for idx,i in enumerate(np.arange(4)):
        print(idx,atoms.get_chemical_symbols()[idx],atoms.positions[idx])
    print()
    ace.ase_relax_volume_only_from_murn_roughly(atoms,verbose=ace.verbose) # necessary for bcc primitive cells which might have a bad starting guess
    ace.ase_relax_volume_only(atoms,verbose=ace.verbose)
    ace.ase_relax_cellshape_and_volume_only(atoms,verbose=ace.verbose)
    print('atoms.cell after relaxing cellshape and volume:')
    print(atoms.cell)
    print('atoms.volume/nat after relaxation:',atoms.get_volume()/atoms.get_number_of_atoms())

    print('atoms.positions after relaxing cellshape and volume (first 4) in angstrom:')
    #for idx,i in enumerate(atoms.positions):
    for idx,i in enumerate(np.arange(4)):
        print(idx,atoms.get_chemical_symbols()[idx],atoms.positions[idx])
    print('atoms.volume:',atoms.get_volume())
    print('atoms.volume/nat:',atoms.get_volume()/atoms.get_number_of_atoms())
    print('a 125 atom bcc cell will relax to fcc if this is thre ground state!, this is great since we than have a 125 atom fcc cell (maybe also possible to do this with smaller cells?)')
    ase_write("POSCAR_eq",atoms,format='vasp')
    ase_write("POSCAR_eq_dir",atoms,format='vasp',direct=True)
    print('atoms.positions')
    print(atoms.positions)
    print()
    print(atoms.get_forces())
    print()


    print("####################################################################")
    print("# GETTING Evinet ... first without Fqh to see if EVinet is working #")
    print("# make this in anycase to relaxe the system fully ...              #")
    print("####################################################################")
    ace.get_murn(atoms,verbose=ace.verbose,return_minimum_volume_frame=True,write_evinet=evinet,write_fqh=False,atomrelax=True)

    # now evinet needs to be set to False! otherwise the evinet/Evinet_1 file is overwritten with fqh volume range (only pos volumes) which is wrong!
    evinet = False
    #print('atoms vol',atoms.get_volume())

    if fqh:
        print("#############################")
        print("# NOW GETTING Fqh_files .. .#")
        print("#############################")
        if False:
            print("#############################")
            print("# fqh from phonopy.. .#")
            print("#############################")
            phonopy_atoms = ase_to_phonopy(atoms)
            phonon = phonopy_pre_process(phonopy_atoms, supercell_matrix=np.eye(3, dtype='intc'))
            supercells = phonon.get_supercells_with_displacements()

            # Force calculations by calculator
            set_of_forces = []
            for scell in supercells:
                cell = phonopy_to_ase(scell)
                #cell = Atoms(symbols=scell.get_chemical_symbols(),
                #             scaled_positions=scell.get_scaled_positions(),
                #             cell=scell.get_cell(),
                #             pbc=True)
                cell.set_calculator(cell)
                forces = cell.get_forces()
                drift_force = forces.sum(axis=0)
                print(("[Phonopy] Drift force:" + "%11.5f" * 3) % tuple(drift_force))
                # Simple translational invariance
                for force in forces:
                    force -= drift_force / forces.shape[0]
                set_of_forces.append(forces)
            sys.exit()

        print("#############################")
        print("# fqh from manual displacements (in get_murn() )... #")
        print("#############################")
        print('oscwd',os.getcwd())
        ace.get_murn(atoms,verbose=ace.verbose,return_minimum_volume_frame=True,write_evinet=False,write_fqh=fqh)
        print('ace.fqh_folder',ace.fqh_folder)
        print('osggg11 tmp3?',os.getcwd())
        import fah
        with cd("fqh"):
            fah.make_fqh_thermos_again()
        #    call(["fqh.py -i Fqh_*at_cell_per_atom* -wqh1 -wqh2 -wqh3"],shell=True)
        #print('osggg222 tmp3',os.getcwd())
        #hier = os.getcwd()
        #for i in ["1st","2nd","3rd"]:
        #    os.chdir(hier)
        #    os.mkdir("fqh/thermo_"+str(i))
        #    with cd("fqh/thermo_"+str(i)):
        #        print('now fqh/termo_'+str(i),os.getcwd())
        #        call(["cp ../../evinet/EVinet_1 EVinet"],shell=True)
        #        call(["cp ../Fqh_*at_cell_per_atom_Surface_"+i+"_order__* Fqh"],shell=True)
        #        print('now gethermodynamcs.sh',os.getcwd())
        #        call(["$dotfiles/thermodynamics/getThermodynamics.sh"],shell=True)



    if fah:
        print("##############################")
        print("# creating anharmonic jobs ... #")
        print("##############################")
        fah_.fah_create_jobs_and_joblist_from_fqh(ace)
        fah_.fah_submit_ipi_ti_job()
    return atoms


def ipi_submit_missing_jobs_to_que():
    return

def ipi_start_job(inputfile="input.xml",sleep=6):
    if not os.path.isfile('in.lmp'):
        sys.exit('in.lmp does not exist! Exit')
    if not os.path.isfile(inputfile):
        sys.exit(inputfile+' does not exist! Exit')
    host = os.environ['myhost']
    print('host',host)
    executable = os.environ["HOME"]+"/Dropbox/Albert/scripts/dotfiles/scripts/executables/lmp_"+host
    print('lammps executable',executable)
    if not os.path.isfile(executable):
        sys.exit(executable+' does not exist! Exit')

    python = "python"
    if host == 'mac': python = "/usr/local/bin/python"
    print(python+" $HOME/sources/ipi/bin/i-pi "+inputfile+" &")
    adress_str = grep(inputfile,"md_ff_*")[0].split()[1]
    print('ad',adress_str)
    filecheck = '/tmp/ipi_'+adress_str #+'='
    print('fc',filecheck)
    if os.path.exists(filecheck):  # here os.path.isfile did not work
        #print('file exists ')
        os.remove(filecheck)
    else:
        print('file does not exist and this is good;')
    call([python+" $HOME/sources/ipi/bin/i-pi "+inputfile+" &"],shell=True)
    print('now sleep for ',sleep,"(sec)")
    time.sleep(sleep)
    time_now = datetime.now()
    print('sleep done; now sart lammps',time_now)
    call([executable+" < in.lmp"],shell=True)
    time_now = datetime.now()
    print('job done',time_now)
    return

def load_hessefiles_volumes_positionsfiles_from_fqh_folder():
    folder = os.getcwd()+"/fqh"
    paths = glob.glob(folder+'/Hessematrix_*')
    out = []
    for path in paths:
        volstr = path.split("_")[-1]
        print("->",path,volstr)
        pos = folder+'/positions_'+volstr+'.extxyz'
        if not os.path.isfile(pos):
            print('pos',pos)
            sys.exit('not found 77')
        out.append([path,float(volstr),folder+'/positions_'+volstr+'.extxyz'])
    if len(out) == 0:
        print('no Hessematrix_* files found in',folder)
    return out

def get_Mg5Si6_and_other_antisites(ace):
    '''
    for Mg5Si6:
    replaces Mg -> Si and Si -> Mg  (22 configurations in Mg5Si6)
    also replace Mg -> Al? and Si -> Al? (since there is no Al)

    for Mg4Al3Si4:
    Mg -> Al or Si
    Si -> Al or Mg
    Al -> Si or Mg
    '''
    doit = [ "Mg5Si6", "Mg5Al2Si4", "Mg4Al3Si4" ]
    for iname in doit:
        #iname = "Mg5Si6"
        path = scripts()+'/tests/Al-Mg-Si/Beta2-bulk/'+iname+'/aiida_exported_group_NN_relaxed_'+iname+"_n2p2_v2ag_calc__only_relaxed.input.data"
        print('path',path)
        frame = ase_read(path,format="runner")
        print('fram.cell()')
        print(frame.cell)
        for idx,i in enumerate(frame.positions):
            print('idx',idx,i,frame.get_chemical_symbols()[idx])

        print()
        print()
        frame_rep = frame.repeat([1,3,2])
        print('fram_rep.cell()')
        print(frame_rep.cell)
        for idx,i in enumerate(frame_rep.positions):
            print('idx',idx,i,frame_rep.get_chemical_symbols()[idx])

        ace.ase_relax_atomic_positions_only(frame_rep,verbose=True,output_to_screen=True)
        ase_write("out_reference_rep"+iname+".runner",frame_rep,format='runner',append=True)
        ase_write("out_reference_rep.runner",frame_rep,format='runner',append=True)
        #chk = ase_get_chemical_symbols_to_number_of_species(frame,known_elements_by_pot=["Al","Mg","Si"])
        #print('chk',chk)

        if False:
            frame_ = frame.copy()
            frame_rep_ = frame_rep.copy()

            struct_written = 0
            for idx,i in enumerate(frame_.positions):  # go through every position
                if idx == 11:
                    print('now done')
                    break
                print('iname',iname,'idx',idx)
                for replace_idx in [0,1]:
                    frame_ = frame.copy()                  # get original frame
                    frame_rep_ = frame_rep.copy()                  # get original frame

                    curr_ele = frame_.get_chemical_symbols()[idx]
                    curr_ele_rep = frame_rep.get_chemical_symbols()[idx]
                    if curr_ele != curr_ele_rep:
                        sys.exit('something went wrong in repetition')
                    if curr_ele == "Mg": replace = [ "Al", "Si"]
                    if curr_ele == "Si": replace = [ "Al", "Mg"]
                    if curr_ele == "Al": replace = [ "Si", "Mg"]
                    replace_now = replace[replace_idx]
                    frame_[idx].symbol = replace_now
                    frame_rep_[idx].symbol = replace_now

                    #for jdx,j in enumerate(frame_.positions):
                    #    print(idx,'-->jdx',jdx,j,frame_.get_chemical_symbols()[jdx])
                    #print()

                    #if idx == 3:
                    #    sys.exit()
                    #ace.ase_relax_atomic_positions_only(frame_)
                    #ace.ase_relax_cellshape_and_volume_only(frame_)
                    print('iname',iname,'idx',idx,'replace_idx',replace_idx,'relax pos')
                    ace.ase_relax_atomic_positions_only(frame_rep_,verbose=True,output_to_screen=True)
                    print('iname',iname,'idx',idx,'replace_idx',replace_idx,'relax cellshape')
                    #ace.ase_relax_cellshape_and_volume_only(frame_rep_,verbose=True)
                    #ase_write("out_antisites"+iname+".runner",frame_,format='runner',append=True)
                    ase_write("out_antisites_rep"+iname+".runner",frame_rep_,format='runner',append=True)
                    struct_written += 1
    return

def ase_repeat_structure(atoms,repeat):
    ''' only repeats cubically; e.g. repeat = 2 repeats the cell 2x2x2'''
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

def ase_repeat_structure_using_listrepeat(atoms,ene_DFT,givenunits,listrepeat=[[1,1,5], [2,3,4], [4,12,3], [13,70,2]]):
    ''' listrepeat = [ [2,3,4], [4,12,3], [13,70,2]] means that
    structures having 1     atoms will be repeated 5^3 times == 5x5x5 = 125 times (125 atoms)
    structures having 2-3   atoms will be repeated 4^3 times == 4x4x4 = 64 times  (128-192 atoms)
    structures having 4-12  atoms will be repeated 3^3 times == 3x3x3 = 27 times  (108-324 atoms)
    structures having 13-70 atoms will be repeated 2^3 times == 2x2x2 = 8 times   (104-560 atoms)
    you can specify different listrepeats to make this more fine grained
    ene_DFT is the DFT energy of the cell
    '''
    nat = atoms.get_number_of_atoms()
    repeat = 0
    if nat == 1: repeat = 5         # *125  at = 125
    for rp in listrepeat:
        if rp[0] <= nat <= rp[1]: repeat = rp[2]
    #if 2 <= nat <= 3: repeat = 4    # *64   at = 124 - 256
    #if 4 <= nat <= 10: repeat = 3    # *27   at = 135 - 270
    #if 11 <= nat <= 40: repeat = 2   # *8    at = 99 -
    #if nat > 40: reppeat = False
    if repeat > 0:
        atoms_sc,forces_sc = my.ase_repeat_structure(atoms,repeat)
        ene_sc_ev = my.convert_energy(ene_DFT*repeat**3.,givenunits,"ev",nat)
    else:
        atoms_sc = copy.deepcopy(atoms)
        forces_sc = atoms_sc.get_forces()
        ene_sc_ev = my.convert_energy(ene_DFT,givenunits,"ev",nat)
    return atoms_sc, forces_sc, ene_sc_ev

class ase_get_known_formats_class():
    """
    helptext
    """
    def __init__(self,verbose = False):
        self.formatspy              = os.path.dirname(ase.io.__file__)+"/formats.py"
        self.vasp_out               = os.path.dirname(ase.io.__file__)+"/vasp.py"
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

    def copy_if_necessary(self):
        #scripts = os.environ['scripts']
        scripts = os.environ['SCRIPTS']  # needs to be capitalized in python 3.x
        from_ = scripts+"/runner_scripts/ase_fileformat_for_"
        to = os.path.dirname(ase.io.__file__)+"/"
        import filecmp
        for ff in self.my_formats_filenames:
            #print('copying ',from_+ff,'to (7666)',to+ff)
            if os.path.isfile(from_+ff) and os.path.isfile(to+ff):
                seemsequal = filecmp.cmp(from_+ff,to+ff)
            else:
                seemsequal = False
            #print('seemsequal',seemsequal)
            #if self.verbose:
            #    print('copying ',from_+ff)
            #    print('                    to',to+ff)
            #    print('seem equal?',seemsequal)
            if seemsequal == False:
                print('copying ',from_+ff,'to (7666)',to+ff)
                shutil.copyfile(from_+ff,to+ff)

    def adapt_formatspy(self,writeformatspy = False):
        # check if formatspy exist
        if not os.path.isfile(self.formatspy):
            print('formatspy',self.formatspy)
            sys.exit('did not find '+str(self.formatspy))

        f = open(self.formatspy, "r")
        contents = f.readlines()
        contents_before = contents.copy()
        f.close()
        insert1=0
        insert2=0
        format_new = False
        for idx,i in enumerate(contents):
            #print('i',idx,i)
            #print("|"+i[:20]+"|")
            if i[:20] == "    'abinit': ('ABIN":
                insert1 = idx
            if i[:30] == "    'lammps-data': 'lammpsdata":
                insert2 = idx
        print('insert1a',insert1)
        print('insert2a',insert2)
        if insert1 == 0 or insert1 == 0:
            format_new = True
            for idx,i in enumerate(contents):
                if i[:10] == "F('abinit'":
                    insert1 = idx
                if i[:15] == "F('lammps-data'":
                    insert2 = idx

        print('insert1b',insert1)
        print('insert2b',insert2)
        if insert1 == 0 or insert1 == 0:
            print('insert1',insert1)
            print('insert2',insert2)
            sys.exit('did not find inserts 98')
        print('insert1c',insert1)
        print('insert2c',insert2)
        for fo in [ 'runner', 'ipi', 'quippy' ]:
            if fo in self.all_known_formats_ase:
                if self.verbose:
                    print(fo.ljust(14)+'format are already added in formats.py (of ase).')
            else:
                print(fo.ljust(14)+'format NOT KNOWN in formats.py (of ase, WILL BE ADDED).')
                if format_new == False:
                    contents.insert(insert1, "    '"+fo+"': ('"+fo+" input file', '+F'),\n")
                else:
                    contents.insert(insert1, "F('"+fo+"', '"+fo+" input file', '+F'),\n")
                writeformatspy = True


        if 'lammps-runner' in self.all_known_formats_ase:
            if self.verbose:
                print('lammps-runner format are already added in formats.py (of ase).')
        else:
            if format_new == False:
                contents.insert(insert1,"    'lammps-runner': ('LAMMPS data input file for n2p2 or runner', '1F'),\n")
                contents.insert(insert2,"    'lammps-runner': 'lammpsrunner',\n")
            else:
                contents.insert(insert1,"F('lammps-runner', 'LAMMPS data input file for n2p2 or runner', '1F', module='lammpsrunner'),\n")
            writeformatspy = True

        if writeformatspy == True:
            print('! adapting formatspy ... ! (seems to be save in any case)')
            #print('insert',insert)

            f = open(self.formatspy, "w")
            contents = "".join(contents)
            f.write(contents)
            f.close()
            sys.exit('! since adapting formatspy need to exit here!')
        else:
            if self.verbose:
                print('everything was already in',self.formatspy)
        return

    def check_if_default_formats_known(self,copy_and_adapt_formatspy_anyhow=False):
        self.get_all_known_formats()
        #print('needcopy_since_missing',self.needcopy_since_missing)
        if self.verbose: print(">> formats.py", self.formatspy)
        for idx,i in enumerate(self.my_formats_shall):
            if self.check_if_format_in_know_formats(i) == False: self.needcopy_since_missing = True
            #filename = self.my_formats_filenames[idx]

        #import filecmp
        #for i in self.my_formats_filenames: # = [ "runner.py","lammpsrunner.py","lammpsdata.py","ipi.py", "quippy.py" ]
        #    # make diff
        #    f1 =
        #    f2 =
        #    print(filecmp.cmp(f1,f2))
        #print('needcopy_since_missing',self.needcopy_since_missing)
        self.copy_if_necessary()
        #sys.exit('88')

        if self.needcopy_since_missing == True or copy_and_adapt_formatspy_anyhow == True: # only in this case copy stuff and add to formats.py
        #    self.verbose = True
            self.copy_if_necessary()
            self.adapt_formatspy(writeformatspy = copy_and_adapt_formatspy_anyhow)

        #################################################
        # check if vasp-out format writes (rather reads) free energies -> makes force_consistent=False --> force_consistent=True
        vaspout = grep(self.vasp_out,'force_consistent=')
        #print('vo',vaspout)
        #print('vo len',len(vaspout))
        if len(vaspout) > 0: # old version
            vo = vaspout[0].split("force_consistent=")
            #print('vo',vo)
            #print('vo[1]',vo[1])
            if vo[1][:5] == "False":
                print('make vasp.out write free energies by default (not sigma->0 energies)')
                sed(self.vasp_out,"force_consistent=False","force_consistent=True")
            if self.verbose: print(">> vasp.py adapted to write free energies by default.")
        else:
            #print('ss',self.vasp_out)
            sed(self.vasp_out,"energy = free_energy \+ de","energy = free_energy #+ de # since I always want the free energy")
            if self.verbose: print(">> vasp.py adapted to write free energies by default.")
        return

def convert_energy(ene,units_in_,units_out_,nat,verbose=False):
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
        ene = ene * nat
    elif len(units_in_split) == 1 and len(units_out_split) == 2:
        ene = ene / nat
    if verbose:
        print('ene (2)',ene)
    return ene

def show_ase_atoms_content(atoms,showfirst=10,comment = ""):
    print()
    print("#####################################################")
    print("# "+comment+" BEGIN #######")
    print("#####################################################")
    print('## atoms')
    print(atoms)
    print('## atoms.get_number_of_atoms()',atoms.get_number_of_atoms())
    print('atoms.get_positions[fistsome]')
    print(atoms.get_positions()[:showfirst])
    print('## elements get_chemical_symbols()')
    print(atoms.get_chemical_symbols()) #[:showfirst])
    print(list(set(atoms.get_chemical_symbols())))
    print('## atoms.cell (in angstrom)')
    print(atoms.cell)
    print('## aa.get_cell_lengths_and_angles()')
    print(atoms.get_cell_lengths_and_angles())
    print('##atom.numbers',atoms.numbers)
    print("#####################################################")
    print("# "+comment+" END #######")
    print("#####################################################")
    print()
    return

def ase_get_atoms(frame,keep_atoms):
    ''' the function will never change the atomsobject '''
    #print('frame')
    #for idx,i in enumerate(frame.get_positions()):
    #    print(idx,i,idx in keep_atoms)
    #print()

    frame_out = frame.copy()
    for idx,midx in enumerate(range(len(frame_out))[::-1]):
        #print(idx,midx,frame_out.get_positions()[midx],midx in keep_atoms)
        if midx in keep_atoms:
            pass
        else:
            del frame_out[midx]
    #print()
    #for idx,i in enumerate(frame_out.get_positions()):
    #    print(idx,i)
    #print('keep_atoms')
    #print(keep_atoms)
    return frame_out

def ase_showpos(atomsc_sphere):
    for idx,i in enumerate(atomsc_sphere.get_positions()):
        print(str(idx).ljust(5),atomsc_sphere.get_chemical_symbols()[idx].ljust(4),i)
    print()
    return

def fakefun():
    print('dd')
    return
############################################################
# math function
############################################################
def crosscorr(datax, datay, lag=0, wrap=False):
    """ Lag-N cross correlation.
    Shifted data filled with NaNs

    Parameters
    ----------
    lag : int, default 0
    datax, datay : pandas.Series objects of equal length
    Returns
    ----------
    crosscorr : float
    """
    if wrap:
        shiftedy = datay.shift(lag)
        shiftedy.iloc[:lag] = datay.iloc[-lag:].values
        return datax.corr(shiftedy)
    else:
        return datax.corr(datay.shift(lag))


if __name__ == "__main__":
    p = help()
    args = p.parse_args()
    if args.verbose:
        print_args(args)
    if args.execute_function:
        function = eval(args.execute_function)
        function()
        hier = os.path.abspath(os.getcwd())
        readmepath = create_READMEtxt(hier)
    #n2p2_check_SF_inputnn(inputnn="cursel_64.def")
    #get_number_of_atoms_as_function_of_cutoff()
    #create_al_structures_for_analysis_SOAP()
    #get_Mg5Si6_and_other_antisites()
