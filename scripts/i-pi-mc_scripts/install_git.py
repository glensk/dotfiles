#!/usr/bin/env python
from __future__ import print_function
import argparse

import sys,os,socket
import myutils as my
import subprocess

myhostname = my.check_for_known_hosts()
#known = ["ipi","ipi_cosmo","n2p2","lammps_runner", "lammps_n2p2","lbzip","lbzip2","atomsk", "vmd", "aiida-alloy" ]
known = ["ipi","ipi_cosmo","eigen", "n2p2","lammps", "lammps_runner", "lammps_n2p2","lbzip","lbzip2","atomsk", "vmd", "aiida-alloy", 'units', "cosmo_tools", "cosmo-tools", 'mlip','miniconda2', 'miniconda3', 'notes', 'ncdu', 'n2p2_edo' ]
# git clone https://github.com/glensk/i-pi.git
# create pull request
# i-pi/tools/py/mux-positions.py
address = {};branch={};install_folder_different={}

install_folder_different["notes"] = True

address["eigen"]        = "https://github.com/eigenteam/eigen-git-mirror.git"

address["ipi_old"]      = "https://github.com/ceriottm/i-pi-mc"
branch['ipi_old']       = "kmc-al6xxx"

address["ipi_cosmo"]    = "https://github.com/cosmo-epfl/i-pi.git";
branch['ipi_cosmo']     = "feat/kmc"   #"feat/kmc-al6xxx"

address["ipi"]          = "https://github.com/glensk/i-pi.git";
branch['ipi']           = "feat/kmc"
branch['ipi']           = False # get all the branches and not just feat/kmc

address["aiida-alloy"]  = "https://gitlab.com/daniel.marchand/aiida-alloy.git"
branch['aiida-alloy']   = False

address["lammps"]        = "https://github.com/lammps/lammps.git";               branch["lammps"]        = False   # this is the preferred way
address["lammps_n2p2"]   = "https://github.com/lammps/lammps.git";               branch["lammps_n2p2"]   = False
#address["lammps_runner"] = "https://github.com/cosmo-epfl/lammps.git";           branch["lammps_runner"] = False
#address["lammps_runner"] = "https://github.com/glensk/lammps.git";               branch["lammps_runner"] = False
address["lammps_runner"] = "https://github.com/lammps/lammps.git";               branch["lammps_runner"] = False

address["n2p2"]          = "https://github.com/CompPhysVienna/n2p2.git";
branch['n2p2']           = 'develop' # this branch is necessary to get scaling.data (or function.data... one of both)

address["units"]         = "http://ftp.gnu.org/gnu/units/units-2.18.tar.gz"
branch["units"]          = False

address["cosmo_tools"]   = "https://github.com/cosmo-epfl/cosmo-tools.git"
branch["cosmo_tools"]    = False
address["cosmo-tools"]   = address["cosmo_tools"]
branch["cosmo-tools"]    = branch["cosmo_tools"]

address["mlip"]         = "http://gitlab.skoltech.ru/shapeev/mlip.git"
branch["mlip"]          = False

def help(p = None ,known=known):
    string='''e.g. install_git.py -i atomsk'''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-i','--install', choices=known, required=True,
            help='choose what to install')

    ### from now on we save everything in the source folder! (not in the dropbox)
    p.add_argument('-sf','--sources_folder', action='store_true', default=os.environ.get('HOME')+"/sources/",
            help='The target directory for installation. Default: $HOME/sources/')
    p.add_argument('-if','--install_folder', type=str, default=False,
            help='The target foldername for installation. Default: if False: the --install name')
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

p = help()
args = p.parse_args()
args.verbose = True

def install_(args,known):
    # make sources folder if it does not exist ($HOME/sources)
    if not os.path.isdir(args.sources_folder):
        os.makedirs(args.sources_folder)

    install        = args.install
    print('>> dotfiles folder',os.environ.get('dotfiles'))
    #with my.cd(os.environ.get('dotfiles')):
    #    print('>> pwd (ssh key should was added: https://github.com/settings/keys) so that it should not be necessary to use git passwords.')
    if args.verbose: print(">> pwd:",os.getcwd())

        # git config --global credential.helper store   ; before git push makes it work without password
        #subprocess.call(["git","config","--global","credential.helper","store"])
        #get_my_address = "https://github.com/glensk/dotfiles.git"  # git config --get remote.origin.url
        #subprocess.call(["git","pull"],shell=True)
        #subprocess.call(["git","push",get_my_address],shell=True)

    # check if the program is known
    if args.install not in known:
        sys.exit("Not known how to install "+args.install+"; Exit!")
    if args.verbose:
        print('>> args.install              : \"'+args.install+'\" (is a known program)')


    # get the address
    try:
        args.git = address[install]
    except KeyError:
        args.git = False
    try:
        args.branch = branch[install]
    except KeyError:
        args.branch = False


    if args.verbose:
        print(">> args.git (address)        :",args.git)
        print(">> args.branch               :",args.branch)

    if args.verbose:
        print(">> args.sources_folder       :",args.sources_folder)
        print(">> args.install_folder (1)   :",args.install_folder)


    # get the install folder
    try:
        args.install_folder = install_folder_different[install]
    except KeyError:
        args.install_folder = False
    if args.verbose:
        print(">> args.install_folder (2)   :",args.install_folder,type(args.install_folder))

    if args.install_folder == False:
        args.install_folder = args.sources_folder+args.install
    if args.verbose:
        print(">> args.install_folder (3)   :",args.install_folder)
        print(">> args.install_folder==True :",args.install_folder==True)
        print(">> args.install_folder==False:",args.install_folder==False)

    if args.install_folder == False:
        sys.exit("args.install_folder == False")
    elif args.install_folder == True:
        pass
    else:
        if os.path.isdir(args.install_folder):
            sys.exit('args.install_folder '+args.install_folder+' does already exist; Exit')

    print("cd "+args.sources_folder)
    with my.cd(args.sources_folder):
        #if args.install   in ['ipi']                    : git_clone(args,specify_depth = False,checkout="feat/kmc")
        if args.install   in ['atomsk']                 : install_atomsk(args)
        elif args.install in ['miniconda','miniconda2'] : install_miniconda(args)
        elif args.install in ['lbzip','lbzip2']         : install_lbzip(args)
        elif args.install in ['n2p2']                   : install_n2p2(args)
        elif args.install in ['n2p2_edo']               : install_n2p2_edo(args)
        elif args.install in ['vmd']                    : install_vmd(args)
        elif args.install in ['units']                  : install_units(args)
        elif args.install in ['notes']                  : install_notes(args)
        elif args.install in ['lammps','lammps_n2p2','lammps_runner']    : install_lammps(args)
        elif args.install in ['ncdu']: install_ncdu(args)
        else: git_clone(args)  # ipi_cosmo, aiia-alloy, mlip, ....

        # not working yet
        if args.install == 'xmgrace': install_xmgrace(args)

    print("DONE")
    return


def git_clone(args,specify_depth = True,checkout=False):
    do = ["git","clone"]
    if specify_depth == True:
        do = do + ["--depth","1"]
    if args.branch != False:
        do = do + ["-b",args.branch]
    do = do + [args.git,args.install_folder]
    print('do: ',do)
    subprocess.call(do)
    os.chdir(args.install_folder)
    print('pwd:',os.getcwd())
    subprocess.call(["git","branch"])
    if checkout != False:
        subprocess.call(["git","checkout",checkout])
    return

def install_ncdu(args):
    subprocess.call(["wget https://dev.yorhel.nl/download/ncdu-linux-i486-1.14.tar.gz && untargz ncdu-linux-i486-1.14.tar.gz"],shell=True)
    print('ka')
    return

def install_units(args):
    print('pwd',os.getcwd())
    subprocess.call(["wget","http://ftp.gnu.org/gnu/units/units-2.18.tar.gz"])
    subprocess.call(["tar","-xvf","units-2.18.tar.gz"])
    print('HOME:',home)
    with my.cd("units-2.18"):
        print('configure ....')
        subprocess.call(['./configure','--prefix='+os.environ["HOME"]+'/.local'])
        print('make ....')
        subprocess.call(['make'])
        print('make install ....')
        subprocess.call(['make','install'])
    return

def install_notes(args):
    ''' raw file from https://github.com/pimterry/notes/blob/master/notes '''
    # cd dotfiles/
    os.chdir(os.environ["dotfiles"]+'/aliases')
    subprocess.call(["wget https://raw.githubusercontent.com/pimterry/notes/master/notes && chmod +x notes && sed -i 's|^NOTES_EXT=.*|NOTES_EXT=\"txt\"|' notes"],shell=True)
    return

def install_lbzip(args):
    ''' works on fidis
        works on mac (also not necessary, since alredy installed)
    '''
    print('pwd',os.getcwd())
    print("wget http://archive.lbzip2.org/lbzip2-2.5.tar.gz")
    subprocess.call(["wget","http://archive.lbzip2.org/lbzip2-2.5.tar.gz"])
    subprocess.call(["tar","-xvf","lbzip2-2.5.tar.gz"])
    home = os.environ["HOME"]
    print('HOME:',home)
    with my.cd("lbzip2-2.5/"):
        print('configure ....')
        subprocess.call(['./configure','--prefix='+home+'/.local'])
        print('make ....')
        subprocess.call(['make'])
        print('make install ....')
        subprocess.call(['make','install'])
    return

def install_lammps(args):
    '''
        fidis   : works for runner & n2p2
        mac     : works for runner
        cosmopc : works for runner
    '''
    if myhostname == 'fidis':
        serialfidis = 'fidis'
        ser_or_par = "par"
    elif myhostname == 'mac':
        serialfidis = 'serial'
        ser_or_par  = "ser"
    else:
        serialfidis = 'serial'
        ser_or_par = "ser"



    #if args.install in ["lammps_runner"]:
    #    git_clone(args,specify_depth = False,checkout="runner-lammps")  # like this it is 405 MB; do without depth or runner-lammps branch wont be there;
    #    extension = [ "runner" ]
    if args.install in [ 'lammps', "lammps_n2p2" ,"lammps_runner"]:  # this is thre preferred way and tries to install both, lammps and runner
        git_clone(args,specify_depth = True)

        n2p2_folder=args.sources_folder+"/n2p2"
        if os.path.isdir(n2p2_folder):
            #extension = [ "n2p2" ]
            extension = [ "n2p2", "runner" ] # this is thre preferred way and tries to install both, lammps and runner (e.g. on fidis)
        else:
            extension = [ "runner" ]

        if hostname == "mac":
            extension = [ "runner" ]

        #if not os.path.isdir(n2p2_folder): sys.exit("please downlaod is enough? or need to install? n2p2 first")
    print('extension:',extension)

    os.chdir(args.install_folder)
    if "n2p2" in extension:
        # ln -s $n2p2_folder lib/nnp
        # cp -r $n2p2_folder/src/interface/LAMMPS/src/USER-NNP src
        os.symlink(n2p2_folder, "lib/nnp")
        my.cp(n2p2_folder+'/src/interface/LAMMPS/src/USER-NNP','src/')

    print('os.getcwd() (1)',os.getcwd())
    os.chdir("src")
    print('os.getcwd() (2)',os.getcwd())
    subprocess.call(["make", "clean-all",])


    list=["yes-CLASS2","yes-KSPACE","yes-MANYBODY","yes-MISC","yes-MOLECULE","yes-REPLICA","yes-RIGID","yes-USER-MISC" ]
    if "runner" in extension:
        list = list + [ "yes-USER-RUNNER" ]
        my.cp(my.scripts()+"/lammps_scripts/src_runner/pair_runner.cpp", args.install_folder+'/src/pair_runner.cpp')
        my.cp(my.scripts()+"/lammps_scripts/src_runner/pair_runner.h", args.install_folder+'/src/pair_runner.h')
        my.cp(my.scripts()+"/lammps_scripts/src_runner/USER-RUNNER", args.install_folder+'/src/USER-RUNNER')
    if "n2p2"   in extension:  list = list + [ "yes-user-nnp" ]
    for i in list:
        subprocess.call(["make", i])

    if args.install == "lammps_runner":
        my.sed("pair_runner.h","^#define MAXNEIGH.*","#define MAXNEIGH 500")

    checkdir = []
    if "runner" in extension: checkdir = checkdir + ['USER-RUNNER']
    if "n2p2" in extension: checkdir = checkdir + ['USER-NNP']
    for i in checkdir:
        if not os.path.isdir(i):
            sys.exit(checkdir+" does not exist; Exit")

    ## copy the pythonfiel (adapted) to make
    my.cp(my.scripts()+'/lammps_makefiles/lammps.py',args.install_folder+'/python/lammps.py')  # this makes sure later on that liblammps.so is also found on mac (which was working when python;from lammps import lammps;lmp = lammps() but not in a script due to issues with LD_LIBRARY_PATH which was not recognized (even if set in python)

    ## compile serial or parallel
    if hostname == 'fidis':
        if not os.path.isdir(os.getcwd()+'/MAKE/MINE'):
            my.cp(my.scripts()+'/lammps_makefiles/fidis_deneb_2018-10-31/MINE',os.getcwd()+'/MAKE')

        print("module load ... && make fidis")
        bash_command("source $MODULESHOME/init/bash && module purge && module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen && module list && make fidis",os.getcwd())
        print()


    elif hostname == 'mac':
        print("####################################################################")
        print("# for mac you want to have the true gcc working: conda install gcc #")
        print("# libnnp and nnp-predict can be compiled to make n2p2 predict energies from input.data")
        print("# libnnpif is not working yet.")
        print("#")
        print("# in n2p2: cd src; make libnnp COMP=gnu # did not work in following steps in the end")
        print("#")
        print("# %g++ --version  --> needs to point to ---> g++ (GCC) 4.8.5")
        print("# in n2p2: cd src; cd libnnp; make static COMP=gnu # !!! THIS worked!!!! (with g++)")
        print("# in n2p2: cd src/application/nnp-predict; make static COMP=gnu # !!! THIS worked!!!! (with g++)")
        print("#")
        print("# /Users/glensk/sources/n2p2/src/application/nnp-predict/nnp-predict input.data     # !!! wow, worked on mac!!!")
        print("#")
        print("# %find ~/sources/n2p2 -name \"InterfaceLammps.h\"")
        print("#      /Users/glensk/sources/n2p2/include/InterfaceLammps.h")
        print("#      /Users/glensk/sources/n2p2/src/libnnpif/InterfaceLammps.h")
        print("#")
        print("####################################################################")
        subprocess.call(["make", serialfidis])  # serialfidis = "serial"
        print()
    else:
        subprocess.call(["make", serialfidis])   # serialfidis = "serial"
        #sys.exit("hostname "+hostname+" not set up yet")

    print()
    print("************ make done ************")
    print()
    print("************ copy executable ************")
    print("************ copy executable ************")
    print("************ copy executable ************")
    print("************ copy executable ************")
    print("************ copy executable ************")
    print("************ copy executable ************")
    print("************ copy executable ************")
    print("************ copy executable ************")
    #### copy the executable
    os.chdir(args.install_folder+"/src")
    executable = 'lmp_'+serialfidis
    if not os.path.isfile(executable):
        sys.exit(executable +" does not exist, .... was not created; Exit")
    print('copy ',executable," to",my.scripts()+"/executables/"+executable+"_"+ser_or_par) #_"+extension)
    my.cp(executable,my.scripts()+"/executables/lmp_"+hostname+"_"+ser_or_par) #_"+extension)
    print()

    ##### now get the lammps libraries for python (to be able to use getEnergies_byLammps.py
    #subprocess.call(["make", 'mpi-stubs'])
    #subprocess.call(["make", 'g++_serial','mode=shlib'])
    print()
    print("************ make mode=shlib xxxx ************")
    print("**** this creates the library liblammps.so in the lammps/scr foler which is necessary for ase lammps calcs")
    print("**** make chmod u+x $dotfiles/scripts/executables/lmp_fidis_par ")
    os.chdir(args.install_folder+"/src")
    print()
    if hostname == 'fidis':
        bash_command("source $MODULESHOME/init/bash && module purge && module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen && module list && make mode=shlib fidis",os.getcwd())
    else:
        bash_command("make mode=shlib "+serialfidis,os.getcwd())
    print()
    os.chdir(args.install_folder+"/python")
    print()
    #print("************ install.py ************")  # is not necessary anymore since I added to $PYTHONPATH and LD_LIBRARY_PATH manually
    #print('pwd:',os.getcwd())
    #subprocess.call(["chmod", 'u+x','install.py'])
    #subprocess.call(['./install.py'])
    print('copied ',executable," to",my.scripts()+"/executables/"+executable+"_"+ser_or_par) #_"+extension)
    return

def install_n2p2_edo(args):
    # git clone --depth 1 --branch cmake https://github.com/edoardob90/n2p2.git n2p2_edo
    return

def install_n2p2(args):
    '''
    see https://compphysvienna.github.io/n2p2/ for more details
    '''
    subprocess.call(["git","clone","--depth","1","-b","develop","https://github.com/CompPhysVienna/n2p2.git",args.install_folder])
    os.chdir(args.install_folder)
    subprocess.call(["git","branch"])
    os.chdir("src")
    print('pwd aa:',os.getcwd())
    #my.cp("makefile","makefile.back")
    #my.cp("makefile.intel","makefile.intel.back")
    #my.cp("libnnptrain/makefile","libnnptrain/makefile.back")

    if myhostname in ["fidis", "helvetios"]:
        COMP="intel"
        GSL_ROOT = "GSL_ROOT"
        PROJECT_CC = "icpc # or icc"
        PROJECT_MPICC = "mpiicpc # mpic++ does not work on fidis/helvetios"
        EIGENPATH = "${EIGEN_ROOT}/include/eigen3"
        includelibs = "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"

    elif myhostname == "daint":
        # edo actually excluded the ipo in makefie.intel
        COMP="intel"
        GSL_ROOT = "EBROOTGSL"
        PROJECT_CC = "CC"
        PROJECT_MPICC = "CC"
        EIGENPATH = "${EIGEN_ROOT}"
        includelibs = "-mkl=parallel -liomp5 -lpthread -lm -ldl"
        # cd $MKLROOT/tools
        # ./mkl_link_tool -libs -opts -c intel_c
        # ./mkl_link_tool -libs -opts -env -c intel_c --os=lnx -a intel64  -> shows that Compiler option: -mkl=parallel is necessary && Linking line:  -liomp5 -lpthread -lm -ldl
        # export EIGEN_ROOT="/users/aglensk/sources/eigen/" && module load daint-mc && module switch PrgEnv-cray PrgEnv-intel && module unload cray-libsci && module load GSL/2.5-CrayIntel-18.08 cray-python/2.7.15.1 cray-fftw

    elif myhostname == 'mac':
        COMP = "intel"  # makes problems on mac
        COMP = "gnu"
        # on mac you can always try
        # cd /Users/glensk/sources/n2p2/src/libnnp
        # %make or %make COMP=gnu -j4 shared or make COMP=gnu  # seem all to be fairly similar
        # with COMP=gnu this works for the standard clang g++ but not for the conda gcc
        # GSL & Eigen are only necessary for training see https://compphysvienna.github.io/n2p2/ (section code structure)
        # the compilatin of libnnpif fails ... could try with intel suite to compile (would also need mpic++)
        #
        # now, first installed all latest intel compilers such that icc, icpc and ifort are version 19.0
        # need: Intel@ Parallel Studio XE Composer Edition for C++ macOS
        # need: Intel@ Parallel Studio XE Composer Edition for Fortran macOS
        #
        # then install openmpi with the intel compilers
        # see: https://software.intel.com/en-us/articles/performance-tools-for-software-developers-building-open-mpi-with-the-intel-compilers
        # for openmpi: ./configure --prefix=/usr/local CC=icc CXX=icpc F77=ifort FC=ifort CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64  (use /usr/local since icc,icpc,ifort are also from /usr/local
        # for openmpi:  make
        # if this fails -> download the whole intes suite for c++ (and not just parts)
        #
        # for openmpi:  make all install  -> this did not work out
        #  in makefile gnu can try CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64
        # with openmpi ithen should have something like mpic++ (gnu)
        #
        # in makefile gnu can try CFLAGS=-m64 CXXFLAGS=-m64 FFLAGS=-m64 FCFLAGS=-m64
    else:
        COMP="gnu" # cosmopc15
        # on mac now EIGEN_ROOT and GSL_ROOT are defined in $scripts/dotfiles/scripts/source_to_add_to_path.sh
        #GLS = "/Users/glensk/miniconda2/pkgs/gsl-2.4-ha2d443c_1005/include/gsl"
        #EIGEN = /Users/glensk/miniconda2/

    # makefile
    my.sed("makefile","^COMP=.*","COMP="+COMP)
    my.sed("makefile","^PROJECT_DIR.*","PROJECT_DIR=./")
    my.sed("makefile","^LIB=libnnp.so libnnpif.so libnnptrain.so pynnp.so","LIB=libnnp.so libnnpif.so libnnptrain.so") # remove pynnp.so
    # makefile.intel
    my.sed("makefile.intel","^PROJECT_GSL=.*","PROJECT_GSL=${"+GSL_ROOT+"}/include")
    my.sed("makefile.intel","^PROJECT_EIGEN=.*","PROJECT_EIGEN="+EIGENPATH)
    my.sed("makefile.intel","^PROJECT_LDFLAGS_BLAS=.*","PROJECT_LDFLAGS_BLAS=-L${"+GSL_ROOT+"}/lib "+includelibs)
    my.sed("makefile.intel","^PROJECT_CC=.*","PROJECT_CC="+PROJECT_CC)
    my.sed("makefile.intel","^PROJECT_MPICC=.*","PROJECT_MPICC="+PROJECT_MPICC)


    if myhostname == 'mac':
        my.sed("makefile.gnu","-fopenmp","#-fopenmp")  # try also with the acutal paths
        # with this libnnp compiles, current probs with libnnpif
        # if in makefiel PROJECT_OPTIONS+= -DNOMPI this defines weather wether to compie with/without MPI and used mpic++ (with) and g++ (without)
    if False:
        ## conda activate basegcc !!!!!!!!!!!!
        my.sed("makefile.gnu","^PROJECT_GSL=.*","PROJECT_GSL=./")  # try also with the acutal paths
        my.sed("makefile.gnu","^PROJECT_EIGEN=.*","PROJECT_EIGEN=./")  # try also with the actual paths
        # on mac: conda install gcc
        # on mac: conda install -c conda-forge gsl   # has now the libgsl... in ~/miniconda2/lib
        # on mac: conda install -c omnia eigen3
        # or
        # on mac: conda install -c conda-forge eigen

    ##############################################################################
    # do on every machine
    ##############################################################################
    f = open("makefile.intel", "r");contents = f.readlines();f.close();insert=0
    for idx,i in enumerate(contents):
        if i[:14] == "PROJECT_EIGEN=": insert = idx
        if i[:14] == "MKL_INCLUDE=${": insert = False
    print('insert bb',insert,type(insert))
    if type(insert) != bool:
        contents.insert(insert, "MKL_INCLUDE=${MKLROOT}/include\n")
        my.rm_if_exists("makefile.intel.new")
        f = open("makefile.intel.new", "w"); contents = "".join(contents);f.write(contents);f.close()
        my.cp("makefile.intel.new","makefile.intel")
    # libnnptrain/makefile
    my.sed("libnnptrain/makefile","^INCLUDES=.*","INCLUDES=-I./ -I${PROJECT_INCLUDE}/ -I${PROJECT_GSL} -I${PROJECT_EIGEN} -I${MKL_INCLUDE}")
    my.sed("libnnptrain/makefile","^PROJECT_DIR.*","PROJECT_DIR=../..")
    # libnnp/makefile  --> this had to be changed in the latest version of the development branch

    makefiles_to_change = [ "libnnp","libnnpif", "libnnptrain"]
    for lib in makefiles_to_change:
        my.sed(lib+"/makefile","^PROJECT_DIR.*","PROJECT_DIR=../..")
        #my.sed(lib+"/makefile","^COMP=.*","COMP="+COMP)

    # module load on fidis
    print("*****cc",os.getcwd())
    if myhostname in ['fidis',"helvetios"]:
        bash_command("module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen && module list && make libnnpif-shared && make",os.getcwd())
    if myhostname == 'daint':
        bash_command('export EIGEN_ROOT="/users/aglensk/sources/eigen/" && module load daint-mc && module switch PrgEnv-cray PrgEnv-intel && module unload cray-libsci && module load GSL/2.5-CrayIntel-18.08 cray-python/2.7.15.1 cray-fftw && module list && make libnnpif-shared && make',os.getcwd())
    return

def install_xmgrace(args):
    subprocess.call(["git","clone","--depth","1","https://github.com/fxcoudert/xmgrace",args.install_folder])
    with my.cd(args.install_folder):
        subprocess.call(["./configure"])  # this once complains about missing: configure: error: M*tif has not been found

def install_miniconda(args):
    with my.cd(os.environ['HOME']):
        #if os.path.isdir(args.install):
        #    sys.exit(args.install+" does already exist!")
        #print('-----------')
        #print(os.getcwd())
        #print()
        #print()
        #print('1. donwload miniconda')
        #print('------------------------------')
        #if args.install == "miniconda2":
        #    version = "2"
        #elif args.install == "miniconda3":
        #    version = "3"
        #else:
        #    sys.exit('either miniconda2 or miniconda3')

        #subprocess.call(["wget https://repo.anaconda.com/miniconda/Miniconda"+version+"-latest-Linux-x86_64.sh -O ~/miniconda"+version+".sh"],shell=True)

        #print('2. install it (in silent mode)')
        #print('------------------------------')
        #subprocess.call(["bash ~/miniconda"+version+".sh -b -p"],shell=True)

        print('3. install modules')
        print('------------------------------')
        subprocess.call(["conda config --add channels intel"],shell=True)
        subprocess.call(["conda install -c -y conda-forge ase"],shell=True)
    sys.exit()


def install_atomsk(args):
    print("atoms was in the end installed with conda!")
    print("by: conda install -y -c conda-forge atomsk")
    sys.exit()
    subprocess.call(["git","clone","--depth","1","https://github.com/pierrehirel/atomsk",args.install_folder])
    os.chdir(args.install_folder)
    os.chdir("src")
    bash_command("make atomsk",os.getcwd())
    return

def install_vmd(args):
    vmdfolder="vmd-1.9.3"

    print("### check if already installed")
    if os.path.exists(vmdfolder):
        sys.exit(vmdfolder+" does already exist in the source folder ! Exit!")

    print("### make sure VMDINSTALLNAME is defined")
    VMDINSTALLNAME = os.environ["VMDINSTALLNAME"]
    print("### make sure VMDINSTALLBINDIR is defined")
    VMDINSTALLBINDIR = os.environ["VMDINSTALLBINDIR"]
    print("### make sure VMDINSTALLLIBRARYDIR is defined")
    VMDINSTALLLIBRARYDIR = os.environ["VMDINSTALLLIBRARYDIR"]


    print("### copy to sources")
    home = os.environ["HOME"]
    vmd = home+"/Dropbox/Albert/scripts/dotfiles/sources/"+vmdfolder+".tar.bzip2"
    if not os.path.exists(vmd):
        sys.exit("vmd is not saved @ "+vmd)
    my.cp(vmd,".")

    print("### untar")
    import tarfile
    tar = tarfile.open("vmd-1.9.3.tar.bzip2", "r:bz2")
    tar.extractall()
    tar.close()

    print("### configure")
    os.chdir(vmdfolder)
    subprocess.call(["./configure"])
    #subprocess.call(['./configure','--prefix='+home+'/.local'])

    ## This would only be necessary if VMD variables were not set
    #vmdpath = home+"/sources/"+vmdfolder+"/"
    #binhome = vmdpath+"/bin"
    #libhome = vmdpath+"/lib/lib"  # make lib/lib to distinguish
    #bash_command("export VMDINSTALLNAME=\"vmd\" && export VMDINSTALLBINDIR=\""+binhome+"\" && export VMDINSTALLLIBRARYDIR=\""+libhome+"/$VMDINSTALLNAME\" && ./configure",os.getcwd())

    print("###",os.getcwd())
    print("### make install")
    os.chdir("src")
    print("###",os.getcwd())
    subprocess.call(["make", "install"])
    return


def bash_command(bashCommand,cwd=False):
    ''' help
    '''
    #print("HIER:",cwd)
    #print("HIER:",os.getcwd())
    #print("BASH:",bashCommand)
    #subprocess.Popen(['/bin/bash', '-c', cmd])
    #process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE,cwd=cwd)
    #output, error = process.communicate()
    subprocess.call(bashCommand,shell=True)
    return

if __name__ == "__main__":
    install_(args,known)
#[ ! -e "$ipi_folder_name" ] && echo "git clone https://github.com/ceriottm/i-pi-mc" && git clone --depth 1 -b kmc-al6xxx https://github.com/ceriottm/i-pi-mc $ipi_folder_name
#cd $ipi_folder_name
#git checkout kmc-al6xxx

