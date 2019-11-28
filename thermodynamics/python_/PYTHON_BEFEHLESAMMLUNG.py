#!/usr/bin/env python

print "hallo"*3
import os
print os.system('ls')
## hallohallohallo
## INCAR  KPOINTS  POTCAR
## 0


##########################################################################
import numpy as np
a = np.arange(10000001)
b = np.arange(10000001)
c = a + b
print a
print b
print c
## [       0        1        2 ...,  9999998  9999999 10000000]
## [       0        1        2 ...,  9999998  9999999 10000000]
## [       0        2        4 ..., 19999996 19999998 20000000]


##########################################################################
import math
print math.sin(45)
## 0.850903524534

##########################################################################
import subprocess
def run(befehl):
    import subprocess
    return subprocess.check_output(befehl, shell=True, stderr=subprocess.STDOUT)
stoich = run('frompath_stoich.sh')
##########################################################################

print type.("hall")
## <type 'str'>
