#!/usr/bin/env python
from __future__ import print_function
import argparse

import sys,os
import myutils as my
#import click
from subprocess import check_output,call

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
print('args.if',args.install_folder)
print('args.i',args.install)

#CONTEXT_SETTINGS = my.get_click_defaults()
#@click.command(context_settings=CONTEXT_SETTINGS)
#@click.option('--install','-i',type=click.Choice(known),required=True,help="choose what to install")
#@click.option('--install_folder','-if',type=str,default=os.environ.get('HOME')+"/sources/",required=False,help="the target folder for installation.")


def install_(args,known):
    install = args.install
    install_folder = args.install_folder
    if not os.path.isdir(install_folder):
        my.mkdir(install_folder)
    if install not in known:
        sys.exit("Not known how to install "+install+"; Exit!")
    install_folder_prog = install_folder+"/"+install
    if os.path.isdir(install_folder_prog):
        sys.exit(install_folder_prog+" does already exist; Exit!")

    print("cd "+install_folder)
    my.cd(install_folder)

    if install == 'ipi':
        print("git clone --depth 1 -b kmc-al6xxx https://github.com/ceriottm/i-pi-mc "+install)
        call(["git","clone","--depth","1","-b","kmc-al6xxx","https://github.com/ceriottm/i-pi-mc",install])
        print(os.getcwd())
        with my.cd(install_folder_prog):
            print(os.getcwd())
            call(["git","checkout","kmc-al6xxx"])
            print("git branch")
            call(["git","branch"])

    if install == 'xmgrace':  # currently not working, see below
        call(["git","clone","--depth","1","https://github.com/fxcoudert/xmgrace",folder])
        with my.cd(install_folder_prog):
            call(["./configure"])  # this once complains about missing: configure: error: M*tif has not been found

    print("DONE")
    return


if __name__ == "__main__":
    install_(args,known)
#[ ! -e "$ipi_folder_name" ] && echo "git clone https://github.com/ceriottm/i-pi-mc" && git clone --depth 1 -b kmc-al6xxx https://github.com/ceriottm/i-pi-mc $ipi_folder_name
#cd $ipi_folder_name
#git checkout kmc-al6xxx

