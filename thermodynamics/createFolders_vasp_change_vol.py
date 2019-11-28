#!/usr/bin/env python
import numpy as np
import sys
import os
import shutil
import utils

vol_alats = np.arange(19.5*2,21.5*2,0.3)


for i in [ 'KPOINTS' , "POSCAR", "POTCAR", "INCAR" ]:
    if os.path.isfile(i) != True:
        sys.exit(i+' does not exist')

hier = os.getcwd()
for ind,i in enumerate(vol_alats):
    print i
    os.chdir(hier)
    folder = str(i)+"_angvol"
    if os.path.isdir(folder):
        print folder,"does exist"
        continue
    else:
        os.makedirs(folder)
        shutil.copyfile("KPOINTS",folder+'/KPOINTS')
        shutil.copyfile("INCAR",folder+'/INCAR')
        shutil.copyfile("POTCAR",folder+'/POTCAR')
        shutil.copyfile("POSCAR",folder+'/POSCAR')
        utils.run2('sed -i \'s|xxxALATxxx|'+str(i)+'|\' '+folder+'/POSCAR')



