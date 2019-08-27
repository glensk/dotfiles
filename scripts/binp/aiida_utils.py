#!/usr/bin/env python
from __future__ import print_function
import os,sys
from subprocess import check_output
from datetime import datetime as datetime


def create_READMEtxt(directory=False,add=False):
    ''' wiretes a README.txt file '''
    if directory == False:
        directory = os.getcwd()

    # get sha
    curr_folder = os.path.dirname(os.path.abspath(__file__))
    pwd = os.getcwd()
    os.chdir(curr_folder)
    sha = check_output(["git","rev-parse","master"]).decode('utf-8')
    os.chdir(pwd)

    # get time
    time_now = datetime.now()

    # name of RADME
    filepath = directory+'/README_'+time_now.strftime("%Y-%m-%d_%H:%M")+'.txt'

    # write README.txt
    strout=os.path.basename(sys.argv[0])+" "+" ".join(sys.argv[1:])
    with open(filepath, "w") as text_file:
        text_file.write("# using https://gitlab.com/daniel.marchand/aiida-alloy\n")
        text_file.write("# used sha: "+sha) #+"\n")
        text_file.write("\n")
        text_file.write(strout+"\n")

    print('written ',filepath)
    return
