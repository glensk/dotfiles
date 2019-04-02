#!/usr/bin/env python
from __future__ import print_function
import argparse

import sys,os
import myutils as my
import subprocess
import myutils

known = ["ipi","ipi_cosmo","n2p2","lammps_runner", "lammps_n2p2","lbzip","lbzip2","atomsk", "vmd", "aiida-alloy" ]
# git clone https://github.com/glensk/i-pi.git
# create pull request
# i-pi/tools/py/mux-positions.py
address = {};branch={}
address["ipi_old"]      = "https://github.com/ceriottm/i-pi-mc"
branch['ipi_old']       = "kmc-al6xxx"

address["ipi_cosmo"]    = "https://github.com/cosmo-epfl/i-pi.git";
branch['ipi_cosmo']     = "feat/kmc"   #"feat/kmc-al6xxx"

address["ipi"]          = "https://github.com/glensk/i-pi.git";
branch['ipi']           = "feat/kmc"
branch['ipi']           = False # get all the branches and not just feat/kmc

address["aiida-alloy"]  = "https://gitlab.com/daniel.marchand/aiida-alloy.git"
branch['aiida-alloy']   = False

address["lammps_n2p2"]   = "https://github.com/lammps/lammps.git";               branch["lammps_n2p2"]  = False
address["lammps_runner"] = "https://github.com/cosmo-epfl/lammps.git";           branch["lammps_runner"] = False
address["lammps_runner"] = "https://github.com/glensk/lammps.git";               branch["lammps_runner"] = False


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
    return p

p = help()
args = p.parse_args()

def install_(args,known):
    install        = args.install

    with my.cd(os.environ.get('dotfiles')):
        print('pwd (ssh key should was added: https://github.com/settings/keys)',os.getcwd())

        # git config --global credential.helper store   ; before git push makes it work without password
        #subprocess.call(["git","config","--global","credential.helper","store"])
        #get_my_address = "https://github.com/glensk/dotfiles.git"  # git config --get remote.origin.url
        #subprocess.call(["git","pull"],shell=True)
        #subprocess.call(["git","push",get_my_address],shell=True)

    # check if the program is known
    if args.install not in known:
        sys.exit("Not known how to install "+args.install+"; Exit!")

    # get the address
    args.git = address[install]
    args.branch = branch[install]

    # mkdir the install folder
    if not os.path.isdir(args.sources_folder):
        my.mkdir(args.sources_folder)

    if args.install_folder == False:
        args.install_folder = args.sources_folder+args.install

    print('args.install       :',args.install)
    print('args.sources_folder:',args.sources_folder)
    print('args.install_folder:',args.install_folder)
    print()
    if os.path.isdir(args.install_folder):
        sys.exit('args.install_folder '+args.install_folder+' does already exist; Exit')

    print("cd "+args.sources_folder)
    with my.cd(args.sources_folder):
        if args.install in ['ipi']              : git_clone(args,specify_depth = False,checkout="feat/kmc")
        if args.install in ['ipi_cosmo']        : git_clone(args)
        if args.install in ['aiida-alloy']      : git_clone(args)
        if args.install in ['atomsk']           : install_atomsk(args)
        if args.install in ['lbzip','lbzip2']   : install_lbzip(args)
        if args.install in ['n2p2']             : install_n2p2(args)
        if args.install in ['vmd']              : install_vmd(args)
        if args.install in ['lammps_runner']    : install_lammps(args)
        if args.install in ['lammps_n2p2']      : install_lammps(args)

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
    subprocess.call(do)
    os.chdir(args.install_folder)
    print('pwd:',os.getcwd())
    subprocess.call(["git","branch"])
    if checkout != False:
        subprocess.call(["git","checkout",checkout])
    return

def install_lbzip(args):
    ''' works on fidis
        works on mac (also not necessary, since alredy installed)
    '''
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
    ''' lammps_runner works on fidis && mac
        lammps_n2p2   works on fidis
    '''
    if args.install == "lammps_runner":
        git_clone(args,specify_depth = False,checkout="runner-lammps")  # like this it is 405 MB; do without depth or runner-lammps branch wont be there;
        extension = "runner"
    elif args.install == "lammps_n2p2":
        git_clone(args,specify_depth = True)
        extension = "n2p2"
        n2p2_folder=args.sources_folder+"/n2p2"
        if not os.path.isdir(n2p2_folder): sys.exit("please downlaod is enough? or need to install? n2p2 first")

    os.chdir(args.install_folder)
    if extension == "n2p2":
        # ln -s $n2p2_folder lib/nnp
        # cp -r $n2p2_folder/src/interface/LAMMPS/src/USER-NNP src
        os.symlink(n2p2_folder, "lib/nnp")
        my.cp(n2p2_folder+'/src/interface/LAMMPS/src/USER-NNP','src/')

    print('os.getcwd() (1)',os.getcwd())
    os.chdir("src")
    print('os.getcwd() (2)',os.getcwd())
    subprocess.call(["make", "clean-all",])


    list=["yes-CLASS2","yes-KSPACE","yes-MANYBODY","yes-MISC","yes-MOLECULE","yes-REPLICA","yes-RIGID","yes-USER-MISC" ]
    if args.install == "lammps_runner": list = list + [ "yes-USER-RUNNER" ]
    if args.install == "lammps_n2p2":   list = list + [ "yes-user-nnp" ]
    for i in list:
        subprocess.call(["make", i])

    if args.install == "lammps_runner":
        my.sed("pair_runner.h","^#define MAXNEIGH.*","#define MAXNEIGH 500")

    if extension == "runner": checkdir = 'USER-RUNNER'
    if extension == "n2p2": checkdir = 'USER-NNP'
    if not os.path.isdir(checkdir):
        sys.exit(checkdir+" does not exist; Exit")

    import socket
    hostname = socket.gethostname()
    print('hostname',hostname)
    if hostname == 'fidis':
        serialfidis = 'fidis'
    elif hostname == 'mac':
        serialfidis = 'serial'

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
        subprocess.call(["make", "serial"])
        print()
    else:
        sys.exit("hostname "+hostname+" not set up yet")

    print()
    print("************ make done ************")
    print()
    print("************ copy executable ************")
    #### copy the executable
    os.chdir(args.install_folder+"/src")
    executable = 'lmp_'+serialfidis
    if not os.path.isfile(executable):
        sys.exit(executable +" does not exist, .... was not created; Exit")
    print('copy ',executable," to",my.scripts()+"/executables/"+executable+"_par_"+extension)
    my.cp(executable,my.scripts()+"/executables/"+executable+"_par_"+extension)
    print()

    ##### now get the lammps libraries for python (to be able to use getEnergies_byLammps.py
    #subprocess.call(["make", 'mpi-stubs'])
    #subprocess.call(["make", 'g++_serial','mode=shlib'])
    print()
    print("************ make mode=shlib xxxx ************")
    os.chdir(args.install_folder+"/src")
    subprocess.call(["make", 'mode=shlib',serialfidis])  # serialfidis can be fidis,serial,mpi
    os.chdir(args.install_folder+"/python")
    print()
    #print("************ install.py ************")  # is not necessary anymore since I added to $PYTHONPATH and LD_LIBRARY_PATH manually
    #print('pwd:',os.getcwd())
    #subprocess.call(["chmod", 'u+x','install.py'])
    #subprocess.call(['./install.py'])
    return

def install_n2p2(args):
    subprocess.call(["git","clone","--depth","1","-b","develop","https://github.com/CompPhysVienna/n2p2.git",args.install_folder])
    os.chdir(args.install_folder)
    subprocess.call(["git","branch"])
    os.chdir("src")
    print('pwd aa:',os.getcwd())
    #my.cp("makefile","makefile.back")
    #my.cp("makefile.intel","makefile.intel.back")
    #my.cp("libnnptrain/makefile","libnnptrain/makefile.back")

    hostname = socket.gethostname()

    if hostname == 'fidis':
        COMP="intel"
    if hostname == 'mac':
        COMP="gnu"
        #GLS = "/Users/glensk/miniconda2/pkgs/gsl-2.4-ha2d443c_1005/include/gsl"
        #EIGEN = /Users/glensk/miniconda2/

    # makefile
    my.sed("makefile","^COMP=.*","COMP="+COMP)
    my.sed("makefile","^PROJECT_DIR.*","PROJECT_DIR=./")
    my.sed("makefile","^LIB=libnnp.so libnnpif.so libnnptrain.so pynnp.so","LIB=libnnp.so libnnpif.so libnnptrain.so") # remove pynnp.so
    # makefile.intel
    my.sed("makefile.intel","^PROJECT_GSL=.*","PROJECT_GSL=${GSL_ROOT}/include")
    my.sed("makefile.intel","^PROJECT_EIGEN=.*","PROJECT_EIGEN=${EIGEN_ROOT}/include/eigen3")
    my.sed("makefile.intel","^PROJECT_LDFLAGS_BLAS=.*","PROJECT_LDFLAGS_BLAS=-L${GSL_ROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
    if hostname == 'mac':
        my.sed("makefile.gnu","^PROJECT_GSL=.*","PROJECT_GSL=./")  # try also with the acutal paths
        my.sed("makefile.gnu","^PROJECT_EIGEN=.*","PROJECT_EIGEN=./")  # try also with the actual paths
        # on mac: conda install gcc
        # on mac: conda install -c conda-forge gsl   # has now the libgsl... in ~/miniconda2/lib
        # on mac: conda install -c omnia eigen3
        # or
        # on mac: conda install -c conda-forge eigen
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
        my.sed(lib+"/makefile","^COMP=.*","COMP=intel")

    # module load on fidis
    print("cc",os.getcwd())
    # @fidis
    bash_command("module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen && module list && make libnnpif-shared && make",os.getcwd())
    # @daint
    #bash_command("module load daint-mc intel craype cray-mpich intel-mkl fftw python/2.7.14 gsl eigen && module list && make libnnpif-shared && make",os.getcwd())

    #subprocess.call(["module","load","intel"],shell=True)
    #subprocess.call(["module","load","intel-mpi"],shell=True)
    #subprocess.call(["module","load","intel-mkl"],shell=True)
    #subprocess.call(["module","load","fftw"],shell=True)
    #subprocess.call(["module","load","python/2.7.14"],shell=True)
    #subprocess.call(["module","load","gsl"],shell=True)
    #subprocess.call(["module","load","eigen"],shell=True)
    return

def install_xmgrace(args):
    subprocess.call(["git","clone","--depth","1","https://github.com/fxcoudert/xmgrace",args.install_folder])
    with my.cd(args.install_folder):
        subprocess.call(["./configure"])  # this once complains about missing: configure: error: M*tif has not been found

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
    myutils.cp(vmd,".")

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

