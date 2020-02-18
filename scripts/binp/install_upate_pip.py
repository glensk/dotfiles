#!/usr/bin/env python
 # -*- coding: utf-8 -*-
from __future__ import print_function
import os,sys,argparse
from subprocess import check_output,call

def help(p = None):
    string = ''' helptext '''
    p = argparse.ArgumentParser(description=string,
            formatter_class=argparse.RawTextHelpFormatter)
    #p.add_argument('-i','--inputfile', required=True, type=str,default=False, help="name of the inputfile inputfile")
    p.add_argument('-v','--verbose', help='verbose', action='count', default=False)
    return p

def todo(args):
    #if not os.path.isfile(args.inputfile):
    #   sys.exit("inputfile "+args.inputfile+" does not exist")
    packages = [ 'numpy', 'pandas', 'phonopy', 'cython', 'pyspglib', 'ase' , 'lmfit', 'intel-numpy', 'tqdm', 'jupyter', 'jupyter_contrib_nbextensions', 'plotly', 'virtualenv', 'pathlib', 'iterm2', 'yfinance', 'quandl', 'yahoo_finance','termcolor', 'BeautifulSoup4', 'html5lib','pandas_datareader' ] # 'aiida-core', #, 'aiida-quantumespresso==3.0.0a5' ]
    # while installing aiida-core I got following errors:
    # ERROR: notebook 6.0.2 has requirement pyzmq>=17, but you'll have pyzmq 16.0.4 which is incompatible.
    # ERROR: notebook 6.0.2 has requirement tornado>=5.0, but you'll have tornado 4.5.3 which is incompatible.
    # ERROR: jupyter-console 6.0.0 has requirement prompt-toolkit<2.1.0,>=2.0.0, but you'll have prompt-toolkit 1.0.18 which is incompatible.
    # brew install pyenv --> https://opensource.com/article/19/5/python-3-default-mac
    for idx,i in enumerate(packages):
        print()
        print()
        print()
        print("##########################################################")
        print("#",idx+1,'package out of',len(packages),i)
        print("##########################################################")
        #pip install --upgrade --user phonopy
        call(['python','-m','pip', 'install', '--upgrade', '--user',i])

    return

if __name__ == '__main__':
    p = help()
    args = p.parse_args()
    todo(args)


