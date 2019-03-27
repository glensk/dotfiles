#!/usr/bin/env python
from __future__ import print_function
import os,sys,re
import click
import numpy as np
import glob #,pathlib
from copy import deepcopy
from socket import gethostname
from shutil import copyfile
from subprocess import check_output,call
from datetime import datetime as datetime   # datetime.datetime.now()
from ase.build import bulk as ase_build_bulk
from ase.constraints import StrainFilter
try:  # not in aiida ase
    from ase.constraints import ExpCellFilter
except ImportError:
    pass
from ase.spacegroup import crystal
from ase.constraints import StrainFilter
from ase.constraints import ExpCellFilter
try:
    from ase.calculators.lammpslib import LAMMPSlib
except ImportError:
    pass

from feos import eos
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.optimize import BFGS
from ase.optimize import LBFGS
from ase.optimize import FIRE
try:
    from ase.optimize import GPMin
except ImportError:
    pass
from ase.optimize.basin import BasinHopping
from ase.optimize.minimahopping import MinimaHopping
import shutil
import time
import myutils as my


start_time = time.time()

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

def n2p2_get_scaling_and_function_data():
    if not os.path.isfile("input.data"):
        sys.exit("Need input.data file")
    if not os.path.isfile("input.nn"):
        sys.exit("Need input.nn file")
    submitfile = scripts()+"/n2p2/submit_scaling_debug.sh"
    if not os.path.isfile(submitfile):
        sys.exit("Need "+submitfile+" file!")

    folder="get_scaling"
    if os.path.isdir(folder):
        sys.exit(folder+" already exists!")

    mkdir(folder)
    os.chdir(folder)
    cp("../input.data")
    cp("../input.nn")
    cp(submitfile)

    submitjob(submitdebug=True,jobdir=os.getcwd(),submitskript="submit_scaling_debug.sh")
    create_READMEtxt()
    return



def submitjob(submit=False,submitdebug=False,jobdir=False,submitskript=False):
    if jobdir == False:
        jobdir = os.getcwd()
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
        copyfile(src,dest)
    elif src_is == "file" and dest_is == "file":
        copyfile(src,dest)
    elif src_is == "file" and dest_is == "dir":
        basename = os.path.basename(src)
        #print("basename:",os.path.basename(src))
        copyfile(src,dest+"/"+basename)
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
    with cd("~/Library"):
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

def get_ase_atoms_object_kmc_al_si_mg_vac(ncell,nsi,nmg,nvac,a0,cubic=False,create_fake_vacancy=False,whichcell="fcc"):
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

def progress(count, total, status=''):
    ''' use:
    for i in range(DB2len):
        my.progress(i,DB2len)
    '''
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ...%s\r' % (bar, percents, '%', status))
    sys.stdout.flush()  # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)
    return



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
    ''' units: eV, eV_pa, hartree, hartree_pa '''
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
    except RuntimeError:
        print("had runtime error")
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
        ene = ene*0.036749325

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
    def __init__(self,pot=False):
        self.elements = False  # list e.g. ['Al', 'Mg', 'Si']
        self.pot = pot         # n2p2_v1ag
        self.pot_all = False   # n2p2_v1ag
        self.fullpath = False  #
        return


    def get_pothpath_and_pot_all(self):
        # get the self.pot && self.pot_all
        scripts = os.environ['scripts']
        allpot_fullpath = glob.glob(scripts+'/potentials/*')
        self.pot_all = []
        onepot_fullpath = False
        onepot_basename = False
        for i in allpot_fullpath:
            #print(i)
            self.pot_all.append(os.path.basename(i))
            if type(self.pot) != bool:
                if self.pot == os.path.basename(i):
                    self.fullpath = i
        return

        #print(pot_all)
        #sys.exit()

    def get_elements_and_atomic_energies(self):
        # get from input.nn the atomic energies
        if type(self.pot) != bool:
            self.pottype = self.pot.split("_")[0]
            if self.pottype == "runner" or self.pottype == "n2p2":
                inputnn = self.fullpath+"/input.nn"
                if os.path.isfile(inputnn):
                    ##### get elements
                    lines = grep(inputnn,"^elements")
                    if len(lines) == 1:
                        line = lines[0]
                        line_elements_ = line.split()[1:]
                        self.elements = []
                        import my_atom
                        for i in line_elements_:
                            #print(i)
                            if i in my_atom.atomic_symbols:
                                #print("yo",i)
                                self.elements.append(i)
                            else:
                                break
                        #print('elements ++',self.elements)

                    ##### get atomic energies
                    lines = grep(inputnn,"^atom_energy")
                    ele_list = []
                    ene_list = []
                    d = {}
                    for i in lines:
                        if i.split()[1] in self.elements:
                            ele = i.split()[1]
                            ene = float(i.split()[2])
                            #print('lines',i.split(),"--------->>",ele,ene,type(ene))
                            ele_list.append(ele)
                            ene_list.append(ene)
                            #print("ele_list",ele_list)
                    #print
                    self.elements = ele_list
                    #print('ele_final:',ele_list)
                    #print('ene_final',ene_list)
                    if len(ele_list) == len(ene_list):
                        d = {}
                        for idx,i in enumerate(ele_list):
                            d[i] = ene_list[idx]
                        self.elements = ele_list
                        self.atom_energy = d
        return

    def get(self):
        self.get_pothpath_and_pot_all()
        self.get_elements_and_atomic_energies()
        return

def pot_all():
    all = mypot()
    all.get()
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
            units=False,
            geopt=False,
            kmc=False,
            verbose=False,
            temp=False,
            ):

        self.pot = pot
        self.units = units.lower()
        self.geopt = geopt          # so far only for ene object.
        self.nsteps = 0
        self.verbose = verbose
        self.atoms = False          # ase atoms object (frame)

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
        self.mypot = mypot(self.pot)
        self.mypot.get()
        #self.eos = [ False, False, False, False] # e0_meV/pa, v0_ang^3/pa, B0, B0der]

        if self.verbose:
            print('ase_calculate_ene, self.pot           :',self.pot)
            print('ase_calculate_ene, self.units         :',self.units)
            print('ase_calculate_ene, self.geopt         :',self.geopt)
            print('ase_calculate_ene, self.nsteps        :',self.nsteps)
            print('ase_calculate_ene, self.lmpcmd        :',self.lmpcmd)
            print('ase_calculate_ene, self.kmc           :',self.kmc)
            print('ase_calculate_ene, self.temp          :',self.temp)
            print('ase_calculate_ene, self.mypot.pot     :',self.mypot.pot)
            #print('ase_calculate_ene, self.mypot.pot_all :',self.mypot.pot_all)
            print('ase_calculate_ene, self.mypot.elements:',self.mypot.elements)
        return

    def pot_to_ase_lmp_cmd(self,kmc=False,temp=False,nsteps=0,ffsocket='inet'):
        ''' geoopt (geometry optimization) is added / or not in
            lammps_write_inputfile(); here only the potential is set.
            ffsocket: ipi ffsocket [ "unix" or "inet" ]
        '''
        self.kmc = kmc
        self.temp = temp
        self.nsteps = nsteps
        self.ffsocket = ffsocket
        if self.ffsocket not in [ "unix", "inet" ]:
            print('ffsocket:',ffsocket)
            sys.exit('ffsocket has to be "unix" or "inet"; Exit!')

        pot = mypot(self.pot)
        pot.get()

        basename = pot.pot
        fullpath = pot.fullpath

        if self.verbose > 1:
            print('...get_potential--basename:',basename)
            print('...get_potential--fullpath:',fullpath)
            print('...pot',self.pot.split("_"))

        # this depends only on the potential which is already defined
        # so should be easy to make this general.
        self.lmpcmd = [
                "mass 1 24.305",
                "mass 2 26.9815385",
                "mass 3 28.0855",
                'variable nnpDir string \"'+fullpath+'\"'
                ]

        #print('spp',self.pot.split("_")[0])
        #sys.exit()
        if self.pot.split("_")[0] == "n2p2":
            # showewsum 1 showew yes resetew no maxew 1000000
            self.lmpcmd = self.lmpcmd + [
                "pair_style nnp dir ${nnpDir} showew no resetew yes maxew 100000000 cflength 1.8897261328 cfenergy 0.0367493254",
                "pair_coeff * * 11.0",
                "#write_data ./pos.data # would this be the final struct?",
            ]
            self.atom_types = {'Mg':1,'Al':2,'Si':3}

        elif self.pot.split("_")[0] == "runner":
            self.lmpcmd = self.lmpcmd + [
                "# thermo 1 # for geopt",
                "pair_style runner dir ${nnpDir} showewsum 1 showew yes resetew no maxew 1000000",
                "pair_coeff * * 7.937658735",
            ]
            self.atom_types = {'Mg':1,'Al':2,'Si':3}

        else:
            sys.exit('pot '+str(self.pot)+' not found!')

        if self.kmc:
            if self.ffsocket == "unix": add = "unix"
            if self.ffsocket == "inet": add = ""
            self.lmpcmd = self.lmpcmd + [
                "",
                "timestep 0.001   # timestep (ps)",
                "velocity all create "+str(self.temp)+" 4928459",  # create initial velocities 4928459 is random seed for velocity initialization"
                "thermo 1   # screen output interval (timesteps)",
                "fix 1 all ipi "+str(gethostname())+" 12345 "+str(add),
                ]
                # with n2p2 in the parallel version, inet is not working
                # "fix 1 all ipi fidis 12345",     # for fidis job
                # "fix 1 all ipi mac 77776 unix",  # for mac job

        if self.verbose > 1:
            print('...lmpcmd',self.lmpcmd)
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
        if self.verbose > 1:
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
            print("UUU pos:",atoms.get_positions()[:4])
            print("UUU for:",atoms.get_forces()[:4])
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
        if self.verbose > 1:
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
            minimizer="LGBFGS"):
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
        if self.verbose > 1:
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
        self.get_calculator(atoms)

        ### attach to atoms to relax the cell
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
            print('UUU atomrelax:',atomrelax)
            print('UUU cellrelax:',cellrelax)
            print("UUU nat:",atoms.get_number_of_atoms())
            print("UUU pos:",atoms.get_positions()[:4])
            print("UUU for:",atoms.get_forces()[:4])
            print("UUU fmx:",abs(atoms.get_forces()).max())
            print("UUU vol:",atoms.get_volume())
            print("UUU vpa:",atoms.get_volume()/atoms.get_number_of_atoms())
        if print_minimization_to_screen:
            print('XXX atomrelax:',atomrelax)
            print('XXX cellrelax:',cellrelax)
            print("XXX nat:",atoms_or_constraint.get_number_of_atoms())
            print("XXX pos:",atoms_or_constraint.get_positions()[:4])
            print("XXX for:",atoms_or_constraint.get_forces()[:4])
            print("XXX fmx:",abs(atoms_or_constraint.get_forces()).max())
            print("XXX vol:",atoms_or_constraint.get_volume())
            print("XXX vpa:",atoms_or_constraint.get_volume()/atoms.get_number_of_atoms())

        if self.verbose > 1:
            print('ZZ done2')

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
        if self.verbose > 1:
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
        asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=keep_alive)
        #atoms_h.set_calculator(asecalcLAMMPS)  # wird in parcalc.py gesetzt
        #/home/glensk/miniconda2/lib/python2.7/site-packages/parcalc/parcalc.py


        #### load the elastic stuff
        from parcalc import ClusterVasp, ParCalculate
        from elastic import get_pressure, BMEOS, get_strain
        from elastic import get_elementary_deformations, scan_volumes
        from elastic import get_BM_EOS, get_elastic_tensor

        print('ene   ',self.ene(atoms_h))
        print('stress1',atoms_h.get_stress())
        self.ase_relax_cellshape_and_volume_only(atoms_h,verbose=False)
        print('stress2',atoms_h.get_stress())
        cryst = atoms_h.copy()
        #print('stress2',cryst.get_isotropic_pressure(cryst.get_stress()))

        sys.exit()


        # Create elementary deformations (systems are ase frames)
        #print('sys',elastic.__file__)
        systems = get_elementary_deformations(atoms_h, n=5, d=0.33)
        print(systems)
        print()
        print(systems[0])

        # Run the stress calculations on deformed cells
        res = ParCalculate(systems, asecalcLAMMPS)

        ## Elastic tensor by internal routine
        #Cij, Bij = get_elastic_tensor(atoms_h, systems=res)
        #print("Cij (GPa):", Cij/units.GPa)
        return

    def get_calculator(self,atoms):
        if self.calculator == "lammps":
            asecalcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=self.keep_alive)
            atoms.set_calculator(asecalcLAMMPS)
        return

    def ase_relax_atomic_positions_only(self,atoms,fmax=0.001,verbose=False):
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
        self.get_calculator(atoms)

        if verbose: self.check_frame('ase_relax_cellshape_and_volume_only in',frame=atoms)

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
            self.ase_relax_atomic_positions_only(atoms_murn,fmax=0.0002,verbose=False)
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

    def get_fh(self,atomsin=False,disp=0.03,debug=False,try_readfile=False,return_units="mev_pa"):
        ''' the function will never change the atomsobject '''
        if return_units == "mev_pa":
            return_mult = 1.
        elif return_units == "ev_cell":
            return_mult = atomsin.get_number_of_atoms()/1000.
        else:
            raise NameError('return_units unknown')

        if try_readfile:
            if os.path.isfile(try_readfile):
                return np.loadtxt(try_readfile)*return_mult

        if atomsin == False:
            sys.exit('need to define atoms in this case XX')
        atoms_h = atomsin.copy()
        atoms_h.wrap()

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
        self.ase_relax_atomic_positions_only(atoms_h,fmax=0.001,verbose=False)
        maxforce = np.abs(atoms_h.get_forces()).max()
        if debug:
            print("###########################################")
            #print("forces harmonic 4:",atoms_h.get_forces()[:3])
            #print("###########################################")
            print('out atoms_h str',atoms_h.get_stress())
            print('out maxforce out',maxforce)

        if maxforce > 0.001:
            print("forces harmonic 4:",atoms_h.get_forces()[:3])
            print('maxforce',maxforce)
            sys.exit('maxforce is too large')

        nat = atoms_h.get_number_of_atoms()
        if nat < 20:
            atoms_h *= (2,2,2)
        if debug:
            nat = atoms_h.get_number_of_atoms()
            print("###########################################")
            print('!!!!!!!nat',nat)
            print("forces harmonic 2:",atoms_h.get_forces()[:3])
            print("###########################################")
        pos0 = atoms_h.get_positions()
        hessematrix=np.zeros((pos0.shape[0]*3,pos0.shape[0]*3))

        ### schleife ueber alle atome, 1..32

        for iidx,i in enumerate(pos0): # loop over all atoms
            progress(iidx,len(pos0))
            #print(iidx,"/",pos0.shape[0]) # loop over xyz 1..3
            for jidx,j in enumerate(i):
                pos1 = np.copy(pos0)
                pos1[iidx,jidx] = pos1[iidx,jidx]+disp

                atoms_h.set_positions(pos1)
                fah = atoms_h.get_forces()
                hessematrix[iidx*3+jidx] = fah.reshape((1,pos0.shape[0]*        3))/(-disp)
        #np.savetxt("HesseMatrix.dat",hessematrix/97.173617,fmt="%.13f")
        import hesse as h
        if debug:
            print("!!!!!!!!!!!!!!!!! get free energy")
        #hes = h.hesseclass(listin=['al'],H=hessematrix)
        if debug:
            print('get_chemical_symbols()')
            print(atoms_h.get_chemical_symbols())
        hes = h.hesseclass(listin=atoms_h.get_chemical_symbols(),H=hessematrix,show_negative_eigenvalues = False, Tmax=1000)
        try:
            #free_ene = (hes.ene_atom[300]-hes.ene_atom[0])[1]
            free_ene = hes.ene_atom
        except IndexError:
            free_ene = "UNSTABLE"
        if debug and type(free_ene) != str:
            hes.write_ene_atom()
        if try_readfile:
            np.savetxt(try_readfile,free_ene)
        return free_ene*return_mult


    def check_frame(self,text,frame=False,verbose=True,setupcalc=True):
        if type(text) != str:
            raise TypeError("need a text!")

        if setupcalc == True:
            self.keep_alive = True
            self.get_calculator(frame)

        #print('fgs',frame.get_stress())
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
    basename,fullpath = mypot(getbasename=ace.pot)
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

def lammps_ext_calc(atoms,ace):
    ''' atoms is an ase atoms object which can hold several frames, or just one'''

    ### mkdir tmpdir
    tmpdir = os.environ['HOME']+"/._tmp_lammps/"
    mkdir(tmpdir)
    print('doint the calculation in lammps, externally, in folder',tmpdir)

    ### write input structure
    if ace.verbose > 1:
        show_ase_atoms_content(atoms,showfirst=10,comment="START LAMMPS EXTERNALLY")
    atoms.set_calculator(None)
    atoms.write(tmpdir+'pos.lmp',format='lammps-runner')
    #sys.exit('pos now written or not')

    ### write inputfile  # geopt is here added or not
    lammps_write_inputfile(folder=tmpdir,filename='in.lmp',positions='pos.lmp',ace=ace)
    if ace.verbose > 1:
        print("written lammsp inputfile to ",tmpdir)

    ### calculate with lammps (trigger externally)
    ene = False
    if os.path.isfile(tmpdir+'log.lammps'):
        os.remove(tmpdir+"log.lammps")


    #ka=ace.pot.split("_")[0]
    #print('aa:'+ka+":",type(ka))
    if ace.pot.split("_")[0] == "runner":
        LAMMPS_COMMAND = my.scripts()+'/executables/lmp_fidis_par_runner'
        #print("LAMMPS_COMMAND 1",LAMMPS_COMMAND)
    else:
        #print("LAMMPS_COMMAND 2",LAMMPS_COMMAND)
        LAMMPS_COMMAND = os.environ['LAMMPS_COMMAND']

    #print("LAMMPS_COMMAND",LAMMPS_COMMAND)
    #sys.exit()

    with cd(tmpdir):  # this cd's savely into folder
        # RUN LAMMPS
        # without SHELL no LD LIBRARY PATH
        #print('pwd',os.getcwd())
        call([LAMMPS_COMMAND+" < in.lmp > /dev/null"],shell=True)

        #print('pwd2',os.getcwd())
        #sys.exit()
        ### extract energy and forces
        ene = check_output(["tail -300 log.lammps | grep -A 1 \"Step Temp E_pai\" | tail -1 | awk '{print $3}'"],shell=True).strip()
        ene=float(ene)
        #print('ene',ene,'lammps in eV')
        #print('ace units',ace.units)
        if ace.units.lower() == 'ev':
            pass
        elif ace.units.lower() == 'hartree':
            ene = ene*0.036749325
        elif ace.units.lower() == 'hartree_pa':
            ene = ene*0.036749325/atoms.get_number_of_atoms()
        elif ace.units.lower() == 'mev_pa':
            ene = ene*1000./atoms.get_number_of_atoms()
        else:
            sys.exit('units '+ace.units+' unknown! Exit!')
        #print('ene out',ene,ace.units)

    if ace.verbose > 1:
        show_ase_atoms_content(atoms,showfirst=10,comment="FINISHED LAMMPS EXTERNALLY")
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

if __name__ == "__main__":
    pass
