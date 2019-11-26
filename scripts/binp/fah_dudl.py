#!/usr/bin/env python

import argparse
import sys
import numpy as np
from argparse import ArgumentDefaultsHelpFormatter
import pot_energy_forces
import os
import utils
import matplotlib.pyplot as plt
import timeit

class dudl( object ):
    '''Class to evaluate dudl file'''
    def __init__(self, filename = False):
        '''
        '''
        self.verbose            = False
        self._filenamein        = filename
        self.filename           = False
        self.dudl               = False
        self.correlation_length = 15

        self.load_dudl()
        self.u                  = self.dudl[:,[0,4]]
        self.uref               = self.dudl[:,[0,5]]
        self.dudl               = np.array([self.dudl[:,0],self.u[:,1]-self.uref[:,1]]).transpose()
        self.dudl               = self.dudl[::self.correlation_length]
        for i in range(1,len(self.dudl)):
            #self.dudl = utils.append_row_to_2d_array(self.dudl,[i,self.dudl[:i,1].std()])
            self.dudl = utils.append_row_to_2d_array(self.dudl,[self.dudl[i,0],self.dudl[:i,1].std()])
        self.dudl       = False
        for i in range(1,len(self.dudl)):
            #self.dudl = utils.append_row_to_2d_array(self.dudl,[i,self.dudl[:i,1].std()/np.sqrt(i)])
            self.dudl = utils.append_row_to_2d_array(self.dudl,[self.dudl[i,0],self.dudl[:i,1].std()/np.sqrt(i)])


    def load_dudl(self):
        ''' load in dudl file '''
        if type(self.dudl) != bool:
            return

        if type(self._filenamein) == bool:
            sys.exit("you have to provide the path + filename for the dudl file")


        if os.path.isfile(self._filenamein) != True:
            if os.path.isdir(self._filenamein) != True:
                print "self._filenamein:",self._filenamein
                sys.exit("path to dUdL file is not correct")

        if os.path.isfile(self._filenamein) != True:
            if os.path.isdir(self._filenamein) == True:
                possible_filenames = [ 'dUdL', 'dudl' ]
                for i in possible_filenames:
                    if os.path.isfile(self._filenamein+"/"+i) == True:
                        self.filename =  self._filenamein+"/"+i

        if os.path.isfile(self._filenamein) == True:
            self.filename = os.path.abspath(self._filenamein)

        self.filename = os.path.abspath(self.filename)

        if os.path.isfile(self.filename) != True:
            print "self.filename:",self.filename
            sys.exit("filename not found")

        self.dudl = np.loadtxt(self.filename)
        return

    def plt(self, show ):
        ''' shows corresponding plot '''
        plt.plot(show[:,0],show[:,1])
        plt.show()
        return()


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='''help string''')
    p.add_argument('-f', '--filename',
            help='define path to dudl file', default=False)
    p.add_argument('-v', '--verbose',action='count',
            help='verbosity level: v:verbose, vv:even more verbose', default=False)
    args = p.parse_args()
    print args
    if args.filename:
        dudl = dudl(args.filename)
