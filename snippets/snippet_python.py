#!/usr/bin/env python
from __future__ import print_function
import argparse



if __name__ == '__main__':
    def help(p = None):
        string = '''
        Helpstring
        '''
        p = argparse.ArgumentParser(description=string,
                formatter_class=argparse.RawTextHelpFormatter)
                #ArgumentDefaultsHelpFormatter)
        p.add_argument('-v','--verbose',help='verbose', action='count', default=False)
        return p
    p = help()
    args = p.parse_args()

