#!/usr/bin/env python
from __future__ import print_function
import argparse,glob,re,shutil,sys

def rename(args):
    ######## find files
    #print("renaming ...")
    print("found following files/folders ...")
    #print('args.files:',args.files)
    for i in args.files:
        print(i)
    print()
    print('found files/folders',len(args.files))
    for i in args.files:
        #new_name = re.sub(args.replace, r'\1_\2\3', i)
        print('i:',i)
        #l = re.compile("(?<!^)\s+(?=[A-Z])(?!.\s)").split(s)
        print('args.replace:'+str(args.replace)+":")
        print('args.withexp:'+str(args.withexp)+":")
        if args.withexp != False:
            print('--> args.replace:'+str(args.replace)+":")
            print('--> re.compile(args.replace):',re.compile(args.replace))
            if type(args.replace) == bool:
                sys.exit('ERROR: please use -r to define what to replace')
            l = re.compile(args.replace).split(i)
            print('--> l:',l)
            for idxj,j in enumerate(l):
                if j == "":
                    l[idxj] = args.withexp
        elif args.append != False:
            print('args.append',args.append)
            l = i+args.append
        else:
            sys.exit('ERROR: please use -w to define what to replace with')
        #print(i,"to",i.split(args.replace))
        newname = ''.join(l)
        print("-->> ## newname:",newname)
        if args.doit:
            shutil.move(i, newname)
    return

if __name__ == '__main__':
    def help(p = None):
        string = '''
        tipp: you could also use mmv! (try nots o mmv.commands)
        caution: use -f without quotes!
        e.g. rename.py -f simulation.pos_* -r .extxyz$ -w .xyz     # shows all changes (.extexy to xyz)
        e.g. rename.py -f simulation.pos_* -r .extxyz$ -w .xyz -d  # ranemes all .extexy to xyz
        e.g. rename.py -f * -a .md # justa append .md to every file
        e.g. rename.py -f * -r .md.txt -w .txt -v -d  # change .md.txt to .txt



        Nothing will every happen until you specity the option: -d / --doit
        '''
        p = argparse.ArgumentParser(description=string,
                formatter_class=argparse.RawTextHelpFormatter)
                #ArgumentDefaultsHelpFormatter)
        p.add_argument('-f','--files',   required=True,  help='regular expression used to search files', type=str, nargs='+')
        p.add_argument('-r','--replace',  required=False, default=False,help='regular expression to be replaced e.g. \".extxyz$\"', type=str)
        p.add_argument('-a','--append',   required=False, default=False,help='append \"string\" to filenema', type=str)
        p.add_argument('-w','--withexp',  required=False, default=False,help='replace with this expression/string e.g. \"xyz\"', type=str)
        p.add_argument('-d','--doit',     required=False, default=False,help="If you really want to renae the files/folders, specify this option",
                action='store_true')
        p.add_argument('-v','--verbose',help='verbose', action='count', default=False)
        return p
    p = help()
    args = p.parse_args()
    if args.verbose:
        print("------------------- args ------------------------------------------------")
        print(args)
        print("-------------------------------------------------------------------------")
    rename(args)

