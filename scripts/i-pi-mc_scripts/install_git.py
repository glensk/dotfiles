#!/usr/bin/env python
from __future__ import print_function
import argparse

import sys,os
import myutils as my
import subprocess

known = ["ipi","n2p2","lbzip"]

def help(p = None ,known=known):
    string='''This is the Help'''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('-i','--install', choices=known, required=True,
            help='choose what to install')
    p.add_argument('-if','--install_folder', action='store_true', default=os.environ.get('HOME')+"/sources/",
            help='The target folder for installation.')
    return p

p = help()
args = p.parse_args()

#import click
#CONTEXT_SETTINGS = my.get_click_defaults()
#@click.command(context_settings=CONTEXT_SETTINGS)
#@click.option('--install','-i',type=click.Choice(known),required=True,help="choose what to install")
#@click.option('--install_folder','-if',type=str,default=os.environ.get('HOME')+"/sources/",required=False,help="the target folder for installation.")


def install_(args,known):
    install_folder = args.install_folder
    install = args.install
    if not os.path.isdir(install_folder):
        my.mkdir(install_folder)
    if install not in known:
        sys.exit("Not known how to install "+install+"; Exit!")
    install_folder_prog = install_folder+"/"+install
    #if os.path.isdir(install_folder_prog):
    #    sys.exit(install_folder_prog+" does already exist; Exit!")

    print("cd "+install_folder)
    with my.cd(install_folder):

        if install == 'ipi':
            print("git clone --depth 1 -b kmc-al6xxx https://github.com/ceriottm/i-pi-mc "+install)
            subprocess.call(["git","clone","--depth","1","-b","kmc-al6xxx","https://github.com/ceriottm/i-pi-mc",install])
            print(os.getcwd())
            with my.cd(install):
                print(os.getcwd())
                subprocess.call(["git","checkout","kmc-al6xxx"])
                print("git branch")
                subprocess.call(["git","branch"])

        if install == 'xmgrace':  # currently not working, see below
            subprocess.call(["git","clone","--depth","1","https://github.com/fxcoudert/xmgrace",install])
            with my.cd(install):
                subprocess.call(["./configure"])  # this once complains about missing: configure: error: M*tif has not been found

        if install == 'n2p2':
            install_n2p2(install)
            #os.chdir(install+"/src")
            #print("cc",os.getcwd())
            #bash_command("module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen && module list && make libnnpif-shared && make",os.getcwd())

    print("DONE")
    return


def install_n2p2(install):
    subprocess.call(["git","clone","--depth","1","-b","develop","https://github.com/CompPhysVienna/n2p2.git",install])
    os.chdir(install)
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
    bash_command("module load intel intel-mpi intel-mkl fftw python/2.7.14 gsl eigen && module list && make libnnpif-shared && make",os.getcwd())

    #subprocess.call(["module","load","intel"],shell=True)
    #subprocess.call(["module","load","intel-mpi"],shell=True)
    #subprocess.call(["module","load","intel-mkl"],shell=True)
    #subprocess.call(["module","load","fftw"],shell=True)
    #subprocess.call(["module","load","python/2.7.14"],shell=True)
    #subprocess.call(["module","load","gsl"],shell=True)
    #subprocess.call(["module","load","eigen"],shell=True)
    return

def bash_command(bashCommand,cwd=False):
    print("HIER:",cwd)
    print("HIER:",os.getcwd())
    print("BASH:",bashCommand)
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

