#!/usr/bin/env python
import os,sys
print(os.name)
print()
print(sys.executable)
print()
print(sys.version)
print()
print(os.__file__)
#for i in os.environ['PYTHONPATH']:
#    print i
print(os.environ['PYTHONPATH'])
print()
print('test1','HOME' in os.environ)
print('test2','LD_LIBRARY_PATH' in os.environ)
print()
if 'LD_LIBRARY_PATH' not in os.environ:
    os.environ['LD_LIBRARY_PATH'] = os.environ['HOME']+'/sources/lammps/src'
print(os.environ['LD_LIBRARY_PATH'])

#for i in os.environ:
#    print i
#print(os.environ['LD_LIBRARY_PATH'])
