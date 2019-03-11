#!/usr/bin/env python
from __future__ import print_function
import argparse

import sys,os
import myutils as my
import subprocess
import myutils

known = ["ipi","n2p2","lammps","lbzip","lbzip2","atomsk", "vmd", "aiida-alloy" ]
# git clone https://github.com/glensk/i-pi.git
# create pull request
# i-pi/tools/py/mux-positions.py
address = {};branch={}
address["ipi_old_DONT_USE"] = "https://github.com/ceriottm/i-pi-mc";                branch['ipi_old']       = "kmc-al6xxx"
address["ipi_cosmo"]        = "https://github.com/cosmo-epfl/i-pi.git";             branch['ipi_cosmo']     = "feat/kmc"   #"feat/kmc-al6xxx"
address["ipi"]              = "https://github.com/glensk/i-pi.git";                 branch['ipi']           = "feat/kmc"
address["aiida-alloy"]      = "https://gitlab.com/daniel.marchand/aiida-alloy.git"; branch['aiida-alloy']   = False
address["lammps"]           = "https://github.com/lammps/lammps.git";               branch["lammps"]        = False


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
        print('pwd',os.getcwd())
        # git config --global credential.helper store   ; before git push makes it work without password
        subprocess.call(["git","config","--global","credential.helper","store"])
        get_my_address = "https://github.com/glensk/dotfiles.git"  # git config --get remote.origin.url
        subprocess.call(["git","pull"],shell=True)
        subprocess.call(["git","push",get_my_address],shell=True)

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
        args.install_folder = args.install

    print("cd "+args.sources_folder)
    with my.cd(args.sources_folder):
        if args.install in ['ipi']              : git_clone(args)
        if args.install in ['aiida-alloy']      : git_clone(args)
        if args.install in ['atomsk']           : install_atomsk(args)
        if args.install in ['lbzip','lbzip2']   : install_lbzip(args)
        if args.install in ['n2p2']             : install_n2p2(args)
        if args.install in ['vmd']              : install_vmd(args)
        if args.install in ['lammps']           : install_lammps(args)

        # not working yet
        if args.install == 'xmgrace': install_xmgrace(args)

    print("DONE")
    return

def git_clone(args):
    print('pwd:',os.getcwd())
    if args.branch == False:
        print("git clone --depth 1 "+args.git+" "+args.install_folder)
        subprocess.call(["git","clone","--depth","1",args.git,args.install_folder])
    else:
        print("git clone --depth 1 -b "+args.branch+" "+args.git+" "+args.install_folder)
        subprocess.call(["git","clone","--depth","1","-b",args.branch,args.git,args.install_folder])
    os.chdir(args.install_folder)
    print('pwd:',os.getcwd())
    subprocess.call(["git","branch"])
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
    sys.exit("not yet finished")
    subprocess.call(["git","clone","--depth","1","https://github.com/lammps/lammps.git",args.install_folder])
    os.chdir(args.install_folder)
    os.chdir("src")
    subprocess.call(["make", "clean-all",])
    list=["yes-CLASS2","yes-KSPACE","yes-MANYBODY","yes-MISC","yes-MOLECULE","yes-REPLICA","yes-RIGID","yes-USER_MISC"]
    for i in list:
        subprocess.call(["make", i])


def install_n2p2(args):
    subprocess.call(["git","clone","--depth","1","-b","develop","https://github.com/CompPhysVienna/n2p2.git",args.install_folder])
    os.chdir(args.install_folder)
    subprocess.call(["git","branch"])
    os.chdir("src")
    print('pwd aa:',os.getcwd())
    #my.cp("makefile","makefile.back")
    #my.cp("makefile.intel","makefile.intel.back")
    #my.cp("libnnptrain/makefile","libnnptrain/makefile.back")

    # makefile
    my.sed("makefile","^COMP=.*","COMP=intel")
    my.sed("makefile","^PROJECT_DIR.*","PROJECT_DIR=./")
    my.sed("makefile","^LIB=libnnp.so libnnpif.so libnnptrain.so pynnp.so","LIB=libnnp.so libnnpif.so libnnptrain.so") # remove pynnp.so
    # makefile.intel
    my.sed("makefile.intel","^PROJECT_GSL=.*","PROJECT_GSL=${GSL_ROOT}/include")
    my.sed("makefile.intel","^PROJECT_EIGEN=.*","PROJECT_EIGEN=${EIGEN_ROOT}/include/eigen3")
    my.sed("makefile.intel","^PROJECT_LDFLAGS_BLAS=.*","PROJECT_LDFLAGS_BLAS=-L${GSL_ROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
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

