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
from ase.io import read as ase_read
from ase.io import write as ase_write
from ase.optimize import BFGS
from ase.optimize import LBFGS
from ase.optimize import FIRE
from ase.optimize.basin import BasinHopping
from ase.optimize.minimahopping import MinimaHopping
import shutil
import time

# from scripts folder
import massedit

start_time = time.time()

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
        text_file.write("# used sha: "+sha) #+"\n")
        text_file.write("# execution time: "+str(time.time() - start_time)+" seconds.")
        text_file.write("\n")
        if add:
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

class cd:
    """Context manager for changing the current working directory"""
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


def ase_enepot(atoms,units='eV',verbose=False):
    ''' units: eV, eV_pa, hartree, hartree_pa '''
    ene = atoms.get_potential_energy()
    if verbose:
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

def ase_get_chemical_symbols_to_conz(atoms):
    symbols = atoms.get_chemical_symbols()
    num = atoms.get_number_of_atoms()
    #print('symbols',symbols)
    #print('num',num)

    uniquesym = set(atoms.get_chemical_symbols())
    d = {}
    for i in uniquesym:
        #print(i,symbols.count(i),num)
        d[i] = float(symbols.count(i))/float(num)

    def dcheck(element):
        if element in d.keys():
            #print(element+" exists")
            pass
        else:
            #print(element+" does not exist")
            d[element] = 0.0

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
        sys.exit('$scripts variable is not defined or is not an existing folder')
    return scripts

def test_and_return_environment_var_path(var):
    variable = os.environ[var]
    if not os.path.isfile(variable):
        sys.exit('variable '+str(var)+' is not defined or is not an existing file')
    return variable


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

def mypotold(getbasename=False):
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

class mypot( object ):
    ''' return a list of available potentials '''
    def __init__(self,pot=False):
        scripts = os.environ['scripts']
        allpot_fullpath = glob.glob(scripts+'/potentials/*')
        pot_all = []
        onepot_fullpath = False
        onepot_basename = False
        for i in allpot_fullpath:
            #print(i)
            pot_all.append(os.path.basename(i))
            if type(pot) != bool:
                if pot == os.path.basename(i):
                    onepot_fullpath = i
                    onepot_basename = pot

        #print(pot_all)
        #sys.exit()
        self.pot = onepot_basename
        self.all = pot_all
        self.fullpath = onepot_fullpath
        #if pot == False:
        #    return pot_all
        #else:
        #    return onepot_basename,onepot_fullpath
        return

def pot_all():
    all = mypot()
    return all.all

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

def create_submitskript_ipi_kmc(filepath,nodes,ntasks,IPI_COMMAND=False,LAMMPS_COMMAND=False,lmp_par=False,ipi_inst=False,ffsocket=False,submittime_hours=71):
    ''' time is in min '''

    def check(variable,command_name_str,typehere):
        if type(variable) != typehere:
            print("ERROR In create_submitskript_ipi_kmc")
            print("PROBLEM",command_name_str,variable)
            sys.exit(command_name_str+" is not a "+str(typehere))

    check(IPI_COMMAND,"IPI_COMMAND",str)
    check(LAMMPS_COMMAND,"LAMMPS_COMMAND",str)
    check(nodes,"nodes",int)
    check(ntasks,"ntasks",int)
    check(lmp_par,"lmp_par",int)
    check(ipi_inst,"ipi_inst",int)
    check(ffsocket,"ffsocket",str)

    text = [
    "#!/bin/bash",
    "",
    "#SBATCH --job-name=NNP-mpi",
    "#SBATCH --get-user-env",
    "#SBATCH --output=_scheduler-stdout.txt",
    "#SBATCH --error=_scheduler-stderr.txt",
    "#SBATCH --nodes="+str(nodes),
    "#SBATCH --ntasks "+str(ntasks),
    "#SBATCH --time=00-"+str(submittime_hours)+":00:00",
    "#SBATCH --constraint=E5v4",
    "",
    "set +e",
    "source $MODULESHOME/init/bash    # necessary for zsh or other init shells",
    "module load intel intel-mpi intel-mkl fftw python/2.7.14",
    "#export OMP_NUM_THREADS=1",  # THIS LETS THE JOBS BE KILLED!
    "touch time.out",
    'date +%s >> time.out',
    "",
    "# sets up the internet/unix socket for connections both for i-PI and on the lammps side",
    'sed -i \'s/<ffsocket.*/<ffsocket name="lmpserial" mode="'+ffsocket+'">/\' input-runner.xml',
    'sed -i \'s/address>.*<.addr/address>\'$(hostname)\'<\/addr/\' input-runner.xml',
    'sed -i \'s/all ipi [^ ]*/all ipi \'$(hostname)\'/\' in.lmp',
    '',
    'rm -f /tmp/ipi_'+gethostname(),
    '',
    'python '+IPI_COMMAND+' input-runner.xml &> log.i-pi &',
    '',
    'sleep 10',
    '',
    'for i in `seq '+str(ipi_inst)+'`',
    'do',
    '      srun --hint=nomultithread --exclusive -n '+str(lmp_par)+' --mem=4G '+LAMMPS_COMMAND+' < in.lmp > log.lmp$i  &',
    #'      srun -n '+str(lmp_par)+' --mem=4G '+LAMMPS_COMMAND+' < in.lmp > log.lmp$i  &',
    'done',
    '',
    'wait',
    'date +%s >> time.out',
    'exit 0',
    ]

    f = open(filepath,'w')
    for i in text:
        f.write(i+"\n")

    return


class ase_calculate_ene( object ):
    '''
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
        self.units = units
        self.geopt = geopt   # so far only for ene object.
        self.nsteps = 0
        self.verbose = verbose
        self.atoms = False # ase atoms object (frame)

        self.lmpcmd = False # in case we run through lammps or ase+lammps
        self.atom_types = False    # for nn pot

        # case of MD or KMC
        self.kmc = kmc
        self.temp = temp

        self.mypot = mypot(self.pot)
        return

    def pot_to_ase_lmp_cmd(self,kmc=False,temp=False,nsteps=False,ffsocket='inet'):
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

        basename = pot.pot
        fullpath = pot.fullpath

        if self.verbose > 1:
            print('...get_potential--basename:',basename)
            print('...get_potential--fullpath:',fullpath)
            print('...pot',self.pot.split("_"))

        self.lmpcmd = [
                "mass 1 24.305",
                "mass 2 26.9815385",
                "mass 3 28.0855",
                "variable nnpDir string \""+fullpath+"\""
                ]

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
                "pair_style runner dir ${nnpDir} showew no resetew yes maxew 1000000",
                "pair_coeff * * 10.0",
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


    def ene(self,atoms=False):
        ''' atoms is an ase object '''
        if atoms == False:
            atoms = self.atoms
        else:
            atoms = atoms

        atoms.wrap()

        if self.verbose > 1:
            show_ase_atoms_content(atoms,showfirst=10,comment = "START ASE INTERNAL CALUCLATION !!!")
            print()
            print("#####################################################")
            print('##--lmpcmd:',self.lmpcmd)
            print('##--atom_types',self.atom_types)
            print('##--geopt',self.geopt)
            print("#####################################################")
            print()


        if self.geopt == False:
            calcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types)
            atoms.set_calculator(calcLAMMPS)
        else:
            # in case of a verbose run:
            #calcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, log_file='./xlolg.lammps.log',tmp_dir="./",keep_alive=True,atom_types=self.atom_types)

            # in case  of a non verbose run
            calcLAMMPS = LAMMPSlib(lmpcmds=self.lmpcmd, atom_types=self.atom_types,keep_alive=True)
            atoms.set_calculator(calcLAMMPS)
            #from ase.io.trajectory import Trajectory
            #traj = Trajectory('ka', mode='w',atoms=atoms)
            #opt = BFGS(atoms,trajectory="ni.traj")
            #opt.run(steps=20)
            minimizer_choices = [ 'BFGS', 'LGBFGS', 'FIRE', 'bh' ]
            minimizer = 'FIRE'
            #minimizer = 'BFGS'
            #minimizer = 'LGBFGS'
            #minimizer = 'bh'
            #minimizer = 'mh'
            print('startminimize....')
            if minimizer == 'BFGS':
                opt1 = BFGS(atoms) #,trajectory="ni.traj")
            elif minimizer == 'LGBFGS':
                opt1 = LBFGS(atoms) #,trajectory="ni.traj")
            elif minimizer == 'FIRE':
                opt1 = FIRE(atoms) #,trajectory="ni.traj")
            elif minimizer == 'bh':
                kB = 1.38064852e-23
                kB = 1.6021765e-19
                opt1 = BasinHopping(atoms=atoms,         # the system to optimize
                  temperature=1*kB, # 'temperature' to overcome barriers
                  dr=0.5,               # maximal stepwidth
                  optimizer=LBFGS,      # optimizer to find local minima
                  fmax=0.1,             # maximal force for the optimizer
                  )
            elif minimizer == 'mh':
                opt1 = MinimaHopping(atoms=atoms)
                opt1(totalsteps=10)
            print('startrun....')
            maxsteps = 200
            if minimizer == 'bh':
                maxsteps = 3
            if minimizer != 'mh': # in all cases but
                opt1.run(steps=maxsteps,fmax=0.005)
                #print('maxsteps                ',maxsteps,type(maxsteps))
                #print('opt1.get_number_of_steps',opt1.get_number_of_steps(),type(opt1.get_number_of_steps()))
                if maxsteps == opt1.get_number_of_steps():
                    print('DID NOT CONVErGE IN '+str(maxsteps)+' number of minimizer steps!')
                    return np.nan
            #print('ene struct 42: without geoopt: -6852.254 eV lmp and ace')
            #print('ene struct 42: without geoopt: -456816.99 meV_pa lmp and ace')
            #print()
            #print('ene struct 42: ENE lammps: -6853.9673 eV 44 steps')
            #print('ene struct 42: ENE    ace: -6852.44299271 eV on 15 steps BFGS')

            #opt1.replay_trajectory('history.traj')
            #print('done....')
            #opt2 = FIRE(atoms,trajectory="ni.traj")
            #opt2.run(steps=20)
            #calcLAMMPS.attach(traj)
            #calcLAMMPS.run()
        if self.verbose > 1:
            print('ZZ done2')
            print('ZZ self.units',self.units)
        ene = ase_enepot(atoms,units=self.units,verbose=self.verbose)
        if self.verbose > 1:
            print('ZZ ene',ene,self.units)
        #sys.exit()
        #ene = atoms.get_total_energy()
        #if self.verbose:
        #    print('ene',ene)
        #return ene,ene/atoms.get_number_of_atoms()*1000.
        if self.verbose > 1:
            show_ase_atoms_content(atoms,showfirst=10,comment="FINISHED ASE INTERNAL CALUCLATION")
            print()
            print()
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
    out2=out.split('\n')
    id=[]
    stat=[]
    cores=[]
    runtime=[]
    path=[]
    for idx,i in enumerate(out2):
        #print(i.split(" "))
        str_list = filter(None, i.split(" "))
        if idx > 0 and len(str_list) > 2:
            #print(str_list)
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

    ### write input structure
    if ace.verbose > 1:
        show_ase_atoms_content(atoms,showfirst=10,comment="START LAMMPS EXTERNALLY")
    atoms.write(tmpdir+'pos.lmp',format='lammps-runner')
    #sys.exit('pos now written or not')

    ### write inputfile  # geopt is here added or not
    lammps_write_inputfile(folder=tmpdir,filename='in.lmp',positions='pos.lmp',ace=ace)

    ### calculate with lammps (trigger externally)
    ene = False
    if os.path.isfile(tmpdir+'log.lammps'):
        os.remove(tmpdir+"log.lammps")

    LAMMPS_COMMAND = os.environ['LAMMPS_COMMAND']

    with cd(tmpdir):  # this cd's savely into folder
        # without SHELL no LD LIBRARY PATH
        call([LAMMPS_COMMAND+" < in.lmp > /dev/null"],shell=True)

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

if __name__ == "__main__":
    pass
