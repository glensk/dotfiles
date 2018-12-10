#!/usr/bin/env python
from __future__ import print_function
import os,sys
import click
import numpy as np
import glob
from socket import gethostname
from shutil import copyfile
from subprocess import check_output,call
from datetime import datetime as datetime   # datetime.datetime.now()
from ase.build import bulk as ase_build_bulk
from ase.calculators.lammpslib import LAMMPSlib



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
    time_now = datetime.now()

    # name of RADME
    filepath = directory+'/README_'+time_now.strftime("%Y-%m-%d_%H:%M")+'.txt'

    # write README.txt
    strout=os.path.basename(sys.argv[0])+" "+" ".join(sys.argv[1:])
    with open(filepath, "w") as text_file:
        text_file.write("# using https://github.com/glensk/dotfiles/trunk/scripts\n")
        text_file.write("# to download it: svn checkout https://github.com/glensk/dotfiles/trunk/scripts\n")
        text_file.write("# used sha: "+sha+"\n")
        text_file.write(strout+"\n")
        if add:
            if type(add) == str:
                text_file.write(add+"\n")
            elif type(add) == list:
                for i in add:
                    text_file.write(i+"\n")


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

def cp(src,dest):
    copyfile(src,dest)

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


def ase_enepot(atoms,units='eV'):
    ''' units: eV, eV_pa, hartree, hartree_pa '''
    ene = atoms.get_potential_energy()

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
        sys.exit('$scripts variable is not defined or is not an existing folder')
    return scripts

def runner_exec(test=False):
    ''' return environment variable runner_exec (RuNNer executable)'''
    runner_exec = os.environ['runner_exec']
    if test == False and not os.path.isfile(runner_exec):
        sys.exit('$runner_exec variable is not defined or is not an existing file')
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

def mypot(getbasename=False):
    ''' return a list of available potentials '''
    scripts = os.environ['scripts']
    allpot_fullpath = glob.glob(scripts+'/potentials/*')
    allpot_basenames = []
    onepot_fullpath = False
    onepot_basename = False
    for i in allpot_fullpath:
        #print(i)
        allpot_basenames.append(os.path.basename(i))
        if type(getbasename) != bool:
            if getbasename == os.path.basename(i):
                onepot_fullpath = i
                onepot_basename = getbasename

    #print(allpot_basenames)
    #sys.exit()
    if getbasename == False:
        return allpot_basenames
    else:
        return onepot_basename,onepot_fullpath

def ase_calculate_ene_from_pot(atoms,lmpcmd=False,atom_types=False,geopt=False,units='eV',verbose=False):
    ''' atoms is an ase object '''
    atoms.wrap()

    if verbose:
        print('XXpot',pot)
        showfirst = 10
        print('XXatomsXX',atoms)
        print('XX atoms')
        #print(atoms.positions)
        print(atoms.get_positions()[:showfirst])
        print('XX elements get_chemical_symbols()')
        print(atoms.get_chemical_symbols()) #[:showfirst])
        print('XX atoms.cell')
        print(atoms.cell)
        print('XX aa.get_cell_lengths_and_angles()')
        print(atoms.get_cell_lengths_and_angles())
        print('XXatom.numbers',atoms.numbers)
        print('XX-------------------------------------------YY')
        print('YY--lmpcmd:',lmpcmd)
        print('YY--atom_types',atom_types)
    if geopt == False:
        calcLAMMPS = LAMMPSlib(lmpcmds=lmpcmd, atom_types=atom_types)
        atoms.set_calculator(calcLAMMPS)
    else:
        calcLAMMPS = LAMMPSlib(lmpcmds=lmpcmd, log_file='./xlolg.lammps.log',tmp_dir="./",keep_alive=True,atom_types=atom_types)
        atoms.set_calculator(calcLAMMPS)
        #from ase.io.trajectory import Trajectory
        #traj = Trajectory('ka', mode='w',atoms=atoms)
        from ase.optimize import BFGS
        from ase.optimize import LBFGS
        from ase.optimize import FIRE
        #opt = BFGS(atoms,trajectory="ni.traj")
        #opt.run(steps=20)
        opt1 = LBFGS(atoms,trajectory="ni.traj")
        opt1.run(steps=20)
        #opt2 = FIRE(atoms,trajectory="ni.traj")
        #opt2.run(steps=20)
        #calcLAMMPS.attach(traj)
        #calcLAMMPS.run()
    if verbose:
        print('done2')
        print('units',units)
    ene = ase_enepot(atoms,units=units)
    if verbose:
        print('ene',ene)
    #sys.exit()
    #ene = atoms.get_total_energy()
    #if verbose:
    #    print('ene',ene)
    #return ene,ene/atoms.get_number_of_atoms()*1000.
    return ene

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
    else:
        try:                           # 7
            a = int(string)
            return [array[int(string)]]
        except ValueError:
            pass
        return False


def pot_to_ase_lmp_cmd(pot,geopt=False,verbose=False):
    basename,fullpath = mypot(getbasename=pot)
    if verbose:
        print('...get_potential--basename:',basename)
        print('...get_potential--fullpath:',fullpath)
        print('...pot',pot.split("_"))
    if pot.split("_")[0] == "n2p2":
        lmpcmd = [
            "mass 1 24.305",
            "mass 2 26.9815385",
            "mass 3 28.0855",
            "variable nnpDir string "+fullpath,
            "pair_style nnp dir ${nnpDir} showew no resetew yes maxew 1000000  cflength 1.8897261328 cfenergy 0.0367493254",
            "pair_coeff * * 17.0",
            "neighbor 0.4 bin",
            "#thermo 1 # only for geopt",
        ]
        att = {'Mg':1,'Al':2,"Si":3}

    elif pot.split("_")[0] == "runner":
        lmpcmd = [
            "mass 1 24.305",
            "mass 2 26.9815385",
            "mass 3 28.0855",
            "variable runnerDir       string "+fullpath,
            "variable runnerCutoff    equal  10.0",
            "# thermo 1 # for geopt",
            "pair_style runner dir ${runnerDir} showew no resetew yes maxew 1000000",
            "pair_coeff * * ${runnerCutoff}",
        ]
        att = {'Mg':1,'Al':2,"Si":3}

    else:
        sys.exit('pot '+str(pot)+' not found!')

    if geopt:
        lmpcmd.append("min_style cg")
        lmpcmd.append("minimize 1.0e-9 1.0e-10 1000 1000")
    if verbose:
        print('...lmpcmd',lmpcmd)
    return lmpcmd, att



if __name__ == "__main__":
    pass
