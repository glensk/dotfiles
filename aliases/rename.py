#!/usr/bin/env python
from __future__ import print_function
import argparse,glob,re,shutil

def rename(args):
    ######## find files
    print("renaming ...")
    for i in args.files:
        #new_name = re.sub(args.replace, r'\1_\2\3', i)
        #print(i,"to",new_name)
        #l = re.compile("(?<!^)\s+(?=[A-Z])(?!.\s)").split(s)
        l = re.compile(args.replace).split(i)
        for idxj,j in enumerate(l):
            if j == "":
                l[idxj] = args.withexp
        #print(i,"to",i.split(args.replace))
        newname = ''.join(l)
        print(i,"to",newname)
        if args.doit:
            shutil.move(i, newname)
    return

if __name__ == '__main__':
    def help(p = None):
        string = '''
        e.g. rename.py -f simulation.pos_* -r .extxyz$ -w .xyz     # shows all changes (.extexy to xyz)
        e.g. rename.py -f simulation.pos_* -r .extxyz$ -w .xyz -d  # ranemes all .extexy to xyz

        Nothing will every happen until you specity the option: -d / --doit
        '''
        p = argparse.ArgumentParser(description=string,
                formatter_class=argparse.RawTextHelpFormatter)
                #ArgumentDefaultsHelpFormatter)
        p.add_argument('-f',   '--files', required=True, help='regular expression used to search files', type=str, nargs='+')
        p.add_argument('-r','--replace',  required=True, help='regular expression to be replaced e.g. \".extxyz$\"', type=str)
        p.add_argument('-w','--withexp',  required=True, help='replace with this expression/string e.g. \"xyz\"', type=str)
        p.add_argument('-d','--doit',     required=False, help="If you really want to renae the files/folders, specify this option",
                default=False,action='store_true')
        p.add_argument('-v','--verbose',help='verbose', action='count', default=False)
        return p
    p = help()
    args = p.parse_args()
    if args.verbose:
        print("------------------- args ------------------------------------------------")
        print(args)
        print("-------------------------------------------------------------------------")
    rename(args)

