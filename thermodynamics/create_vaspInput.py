#!/bin/python

#-----set parameters and paths------------------------
defKP='8'
defType='fcc'
defElem='Al'
potType='PAW-GGA-PBE_vasp5.2'
supportedTypes=['fcc', 'bcc']
#-----------------------------------------------------

import sys, os, shutil
from P import utils

path      = os.path.dirname(sys.argv[0])
script    = os.path.basename(sys.argv[0])
optionsIn = sys.argv[1:]

usage   = script+'  [Options]'

needed  = [ '' ]

options = [ 
            '-e ELEM  use ELEM as element (default: '+defElem+')',
            '-a ALAT  use ALAT as lattice constant (default: adjust to ELEM)',
            '-s STR   use STR as atomic structure (default: '+defType+')',
            '-k KP    use KPxKPxKP for KPOINTS file (default: '+defKP+')',
            '-t       include thermodynamic integration block in INCAR',
            '-m       include MD block in INCAR',
            '-r       include relaxation block in INCAR',
            '-f       force overwriting previous input files'
          ]

details = [ 
            'script creates basic VASP input (INCAR,KPOINTS,POSCAR,POTCAR)',
            '',
            'the default values are set to a primitive fcc Al (PBE) cell with',
            'low convergence parameters to perform quick test runs',
            '',
            'the element (effectively the POTCAR) can be changed with "-e ELEM"',
            'where ELEM must correspond to one of the folder names in:',
            os.path.join(path,'vasp_potentials','PAW-GGA-PBE_vasp5.2'),
            'Example: '+script+' -e Al',
            '',
            'alternatively ELEM can correspond to the path to a POTCAR file',
            'you can for instance use one of the POTCARs provided in:',
            os.path.join(path,'vasp_potentials','PAW-{GGA-PBE,LDA}{,_vasp5.2}','ELEM','POTCAR')
          ]

if '-h'    in optionsIn: utils.printHelpAndExit(usage,needed,options)
if '-help' in optionsIn: utils.printDetailsAndExit(script,details)

utils.checkOptions(optionsIn,options)

# check -f option for overwriting
if '-f' not in optionsIn:
  if os.path.isfile('POTCAR')  or os.path.isfile('INCAR') or os.path.isfile('KPOINTS') or os.path.isfile('POSCAR'):
     utils.error("previous vasp input exists; use -f to force overwriting")

# check -e option for element
if '-e' in optionsIn:
  elem = utils.getValue('-e',optionsIn)
  if elem.find('POTCAR') > -1:
    tmpElem = elem.split('/')[-2]
    potType = elem.split('/')[-3]
    pot = elem
    elem = tmpElem
  else:
    pot = os.path.join(path,'vasp_potentials',potType,elem,'POTCAR')
else:
  elem = defElem
  pot = os.path.join(path,'vasp_potentials',potType,elem,'POTCAR')

# check and copy POTCAR
utils.check(pot)
shutil.copyfile(pot,'./POTCAR')

# read ENMAX from POTCAR
f = open(pot,'r')
l = f.readline()
while 'ENMAX' not in l: l = f.readline()
encut = l.split(';')[0].split('=')[1]
encut = str(int(round(float(encut))))
f.close()

INCAR =   [
           ' ISTART  = 0         ! 0 = start from scratch,    1 = read WAVECAR',
           ' ICHARG  = 2         ! 2 = atomic superpositions, 0 = from WAVECAR',
           ' ISMEAR  = 1         ! 1 = Methfessel-Paxton,    -1 = Fermi',
           ' SIGMA   = 0.1',
           ' ENCUT   = '+encut,
           ' PREC    = Accurate ! or Normal or Medium or Low',
           ' ADDGRID = TRUE',
           ' ALGO    = NORMAL',
           ' EDIFF   = 1E-5',
           ' LWAVE   = FALSE     ! write WAVECAR ?',
           ' LCHARG  = FALSE     ! write CHGCAR ?'
          ]

RELAX =   [
           ' IBRION  = 1        ! 0 = MD, 1 = QuasiNewton, 2 = conj. gradient',
           ' NSW     = 50',
           ' EDIFFG  = 1E-3     ! a negative value will use forces as criterion'
          ]

MD =      [
           ' IBRION  = 0        ! 0 = MD, 1 = QuasiNewton, 2 = conj. gradient',
           ' SMASS   = 0        ! if 0 Nose mass corresponding to period of 40 time steps will be chosen',
           ' POTIM   = 5        ! times step in femto seconds',
           ' TEBEG   = 1000     ! temperature in K',
           ' NSW     = 1000',
           ' ISYM    = 0'
          ]
           
TI =      [
           ' IBRION   = 0        ! 0 = MD, 1 = QuasiNewton, 2 = conj. gradient',
           ' POTIM    = 5        ! times step in femto seconds',
           ' TEBEG    = 1000     ! temperature in K',
           ' NSW      = 1000',
           ' ISYM     = 0',
           ' SEED     = -1',
           ' LAMBDA   = 0.5',
           ' GAMMA_LD = 0.01',
           ' REF_TYPE = hesse',
           ' PRE_EQ_N = 10000',
           ' DUDL_OFFSET = 20'
          ]

if '-r' in optionsIn and '-m' in optionsIn: error("options -r and -m cannot be used simultaneously")
if '-r' in optionsIn and '-t' in optionsIn: error("options -r and -t cannot be used simultaneously")
if '-m' in optionsIn and '-t' in optionsIn: error("options -m and -t cannot be used simultaneously")

# write INCAR file
f = open('INCAR','w')
for i in INCAR: f.write(str(i)+'\n')
f.write('\n')
if '-r' in optionsIn:
  for i in RELAX: f.write(str(i)+'\n')
if '-m' in optionsIn:
  for i in MD:    f.write(str(i)+'\n')
if '-t' in optionsIn:
  for i in TI:    f.write(str(i)+'\n')
f.close()

# check -k option for k-points
if '-k' in optionsIn:
  kp = utils.getValue('-k',optionsIn,typ='int')
  kp = ' '+str(kp)+' '+str(kp)+' '+str(kp)
else:
  kp = ' '+defKP+' '+defKP+' '+defKP

KPOINTS = [
           'K-Points',
           ' 0',
           'Monkhorst Pack',
           kp,
           ' 0 0 0'
          ]

# write KPOINTS file
f = open('KPOINTS','w')
for i in KPOINTS: f.write(str(i)+'\n')
f.close()

# check -s option for strType
if '-s' in optionsIn:
  strType = utils.getValue('-s',optionsIn,supported=supportedTypes)
else:
  strType = defType

if strType == 'fcc':
  a1='0 0.5 0.5'
  a2='0.5 0 0.5'
  a3='0.5 0.5 0'

if strType == 'bcc':
  a1='-0.5 0.5 0.5'
  a2='0.5 -0.5 0.5'
  a3='0.5 0.5 -0.5'

# check -a option for aLat
if '-a' in optionsIn:
  aLat = utils.getValue('-a',optionsIn,typ='numeric')
else:
  # extract aLat from calculated equilibrium volumes
  f = open(os.path.join(path,'utilities','volumes.dat'),'r')
  l = f.readline()
  while l and l.split()[-1] != elem and l.split()[-1].split('_')[0] != elem:
    l = f.readline()
  if not l:
    utils.warning('element '+elem+' not found in volumes.dat; using 4.0 for aLat')
    aLat = '4.0'
  else:
    if strType == 'bcc': aLat = l.split()[3]
    if strType == 'fcc': aLat = l.split()[4]
  f.close()
  
POSCAR =  [
           'Comment',
           aLat,
           a1,
           a2,
           a3,
           1,
           'Cartesian',
           '0 0 0'
          ]

# write POSCAR
f = open('POSCAR','w')
for i in POSCAR: f.write(str(i)+'\n')
f.close()

# printout
print '\n vasp input created for:   '+strType+' '+elem+'   '+potType+'   '+aLat+' Ang   '+encut+' eV   '+kp+' kp'

