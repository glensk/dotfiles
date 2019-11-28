#!/bin/python

import sys, os
from P import utils

path      = os.path.dirname(sys.argv[0])
script    = os.path.basename(sys.argv[0])
optionsIn = sys.argv[1:]

usage   = script+'  [Options]'

needed  = [ '' ]

options = [ '-p       create parameters.dat and exit',
            '-I       create INCAR template and exit'  ]

details = [ 'script creates ...' ]

if '-h'    in optionsIn: utils.printHelpAndExit(usage,needed,options)
if '-help' in optionsIn: utils.printDetailsAndExit(script,details)

utils.checkOptions(optionsIn,options)

