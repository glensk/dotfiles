from
https://aiida-core.readthedocs.io/en/latest/install/details/database.html

ssh cosmopc
su local  (115 urgences)
sudo apt install rabbitmq-server
sudo apt install postgresql postgresql-server-dev-all postgresql-client
sudo su - postgres
psql
CREATE USER aiida WITH PASSWORD 'aiidadb';   
CREATE DATABASE aiida_al6xxx OWNER aiida ENCODING 'UTF8' LC_COLLATE='en_US.UTF-8' LC_CTYPE='en_US.UTF-8' TEMPLATE=template0; # create the database
CREATE DATABASE aiida_al6xxx_test OWNER aiida ENCODING 'UTF8' LC_COLLATE='en_US.UTF-8' LC_CTYPE='en_US.UTF-8' TEMPLATE=template0;
\l # to see db's
GRANT ALL PRIVILEGES ON DATABASE aiida_al6xxx to aiida;
GRANT ALL PRIVILEGES ON DATABASE aiida_al6xxx_test to aiida;
\q # to quit
psql -h localhost -d aiida_al6xxx -U aiida -W  # pw: aiidadb

# rest can now be done as user
ssh cosmopc
pip install --user virtualenv
virtualenv aiida
source aiida/bin/activate
cd aiida
mkdir .aiida
vi bin/activate # add at last line: export AIIDA_PATH="${VIRTUAL_ENV}/.aiida"
deactivate
source aiida/bin/activate

cd $HOME
mcd code/aiida
#git clone git@github.com:aiidateam/aiida_core aiida-core
git clone https://github.com/aiidateam/aiida_core aiida-core
cd aiida-core
pip install -e .'[atomic_tools]'

# read https://aiida-core.readthedocs.io/en/latest/working_with_aiida/index.html#command-line-interface (concept section)


verdi profile list
verdi setup  --help
verdi setup al6xxx

##########################################################################################
# setup the databases
##########################################################################################
%verdi setup al6xxx
Email Address (identifies your data when sharing): albert.glensk@gmail.com
First name: Albert
Last name: Glensk
Institution: EPFL
Setting up profile al6xxx.
AiiDA backend (available: django, sqlalchemy - sqlalchemy is in beta mode): sqlalchemy
Default user email: albert.glensk@gmail.com
Database engine (available: postgresql_psycopg2 - mysql is deprecated): postgresql_psycopg2
PostgreSQL host: localhost
PostgreSQL port: 5432
AiiDA Database name: aiida_al6xxx
AiiDA Database user: aiida
AiiDA Database password: aiidadb
AiiDA repository directory: /home/glensk/aiida/.aiida/repository-al6xxx/
The repository /home/glensk/aiida/.aiida/repository-al6xxx/ will be created.
Executing now a migrate command...

Info: enter "?" for help
First name: Albert
Last name: Glensk
Institution: EPFL
Password [***]:  (aiidadb)
Repeat for confirmation:
>> User Albert Glensk (albert.glensk@gmail.com) created. <<
Setup finished.

verdi profile list
verdi profile show

vi ~/aiida/bin/activate
 80 # This will set the aiida instance
 81 export AIIDA_PATH="${VIRTUAL_ENV}/.aiida"
 82
 83 # this will enable autocompletion for verdi
 84 eval "$(_VERDI_COMPLETE=source verdi)"


%verdi setup al6xxx_test
Email Address (identifies your data when sharing): albert.glensk@gmail.com
First name: Albert
Last name: Glensk
Institution: EPFL
Setting up profile al6xxx_test.
AiiDA backend (available: django, sqlalchemy - sqlalchemy is in beta mode): sqlalchemy
Default user email: albert.glensk@gmail.com
Database engine (available: postgresql_psycopg2 - mysql is deprecated): postgresql_psycopg2
PostgreSQL host: localhost
PostgreSQL port: 5432
AiiDA Database name: aiida_al6xxx_test
AiiDA Database user: aiida
AiiDA Database password: aiidadb
AiiDA repository directory: /home/glensk/aiida/.aiida/repository-al6xxx_test/
The repository /home/glensk/aiida/.aiida/repository-al6xxx_test/ will be created.
Executing now a migrate command...
...for SQLAlchemy backend
It is time to perform your first SQLAlchemy migration.
Migrating to the last version
Warning: Invalidating all the hashes of all the nodes. Please run verdi rehash
Database was created successfully
Loading new environment...
Installing default AiiDA user...
Starting user configuration for albert.glensk@gmail.com...
Info: enter "?" for help
First name: Albert
Last name: Glensk
Institution: EPFL
Password [***]:
Repeat for confirmation:
>> User Albert Glensk (albert.glensk@gmail.com) created. <<
Setup finished.

# to make backup of database: psql -U aiida -h localhost -p 5432 -d aiida_al6xxx > aida.psql
# the repository needs also to be saved



##########################################################################################
# setup the profiles
##########################################################################################
verdi profile show shows aiidadb_repository_uri  file:///home/glensk/aiida/.aiida/repository-al6xxx/
as the repository
verdi profile list
verdi profile setdefault PROFILENAME # to change profile


##########################################################################################
# setup the computers
##########################################################################################

###########################
# setup the local computer
###########################
%verdi computer setup
Info: enter "?" for help
Computer label: cosmopc18            # now this is localhost
Hostname: cosmopc18                  # now this is localhost
Description []: cosmopc18
Enable the computer? [True]:
Transport plugin: local
Scheduler plugin: ?
Info: Scheduler type.
Select one of:
  direct
  lsf
  pbspro
  sge
  slurm
  torque
Scheduler plugin: direct
Shebang line (first line of each script, starting with #!) [#!/bin/bash]:
Work directory on the computer [/scratch/{username}/aiida/]: /local/scratch/glensk/aiida_scratch
Mpirun command [mpirun -np {tot_num_mpiprocs}]:
Default number of CPUs per machine: 1

# :wq to sve the vim file


###########################
# configure the local computer
###########################
%verdi computer configure local cosmopc18
Info: enter "?" for help
Connection cooldown time (sec) [2.0]: 0
Info: Configuring computer cosmopc18 for user albert.glensk@gmail.com.
Success: cosmopc18 successfully configured for albert.glensk@gmail.com

%verdi computer configure show cosmopc18
* safe_interval  0

verdi computer test cosmopc18


###########################
# setup computer daint 
###########################
(aiida) 11:56:09 glensk@lammm@cosmopc18 /home/glensk/code/aiida/aiida-core
%verdi computer setup
Info: enter "?" for help
Computer label: daint
Hostname: daint.cscs.ch
Description []: daint
Enable the computer? [True]:
Transport plugin: ssh
/s|u\-
Scheduler plugin: Shebang line (first line of each script, starting with #!) [#!/bin/bash]:
Work directory on the computer [/scratch/{username}/aiida/]: /scratch/snx3000/aglensk/aiida_scratch      # is apparently created by the setup
Mpirun command [mpirun -np {tot_num_mpiprocs}]: srun
Default number of CPUs per machine: 18


# on cosmopc: %ssh-keygen
Generating public/private rsa key pair.
Enter file in which to save the key (/home/glensk/.ssh/id_rsa):
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Your identification has been saved in /home/glensk/.ssh/id_r



###########################
# configure computer daint 
###########################
%verdi computer configure ssh daint
Info: enter "?" for help
User name [aglensk]:
port Nr [22]:
Look for keys [True]:
SSH key file [/home/glensk/.ssh/id_rsa]:
Connection timeout in s [60]:
Allow ssh agent [True]:
SSH proxy command [ssh aglensk@ela.cscs.ch /usr/bin/netcat daint.cscs.ch 22]:
Compress file transfers [True]:
GSS auth [False]:
GSS kex [False]:
GSS deleg_creds [False]:
GSS host [daint.cscs.ch]:
Load system host keys [True]:
Key policy [RejectPolicy]:
Connection cooldown time (sec) [30]:
Info: Configuring computer daint for user albert.glensk@gmail.com.
Success: daint successfully configured for albert.glensk@gmail.com


verdi computer test daint



cd /home/glensk/code/aiida
git clone https://github.com/aiidateam/aiida-quantumespresso
cd aiida-quantumespresso
pip install -e .
%verdi calculation plugins

verdi work plugins

launch_calculation_pw --help

vim /home/glensk/code/aiida/aiida-quantumespresso/aiida_quantumespresso/cli/calculations/pw/base.py # have a look at the click stuff !

verdi computer list

verdi shell
from aiida.orm lkmport load_computer
computer = load_computer(1)
computer = load_computer(2)
    computer.uuid
computer = load_computer(computer.uuid)
computer = load_computer('daint')


%vim aiida_quantumespresso/cli/calculations/pw/base.py
(aiida) 10:32:41 glensk@lammm@cosmopc18 /home/glensk/code/aiida/aiida-quantumespresso
%verdi code setup
Info: enter "?" for help
Label: pw-v6.3
Description []: Quantum ESPRESSO pw.x v6.3
Default calculation input plugin: ?
Info: Default calculation plugin to use for this code.
Select one of:
  calculation
  function
  inline
  job
  quantumespresso.cp
  quantumespresso.dos
  quantumespresso.matdyn
  quantumespresso.namelists
  quantumespresso.neb
  quantumespresso.ph
  quantumespresso.pp
  quantumespresso.projwfc
  quantumespresso.pw
  quantumespresso.pw2wannier90
  quantumespresso.pwimmigrant
  quantumespresso.q2r
  simpleplugins.arithmetic.add
  simpleplugins.templatereplacer
  work
Default calculation input plugin: quantumespresso.pw?
Error: entry point 'quantumespresso.pw?' is not valid for any of the allowed entry point groups: aiida.calculations
Info: Default calculation plugin to use for this code.
Select one of:
  calculation
  function
  inline
  job
  quantumespresso.cp
  quantumespresso.dos
  quantumespresso.matdyn
  quantumespresso.namelists
  quantumespresso.neb
  quantumespresso.ph
  quantumespresso.pp
  quantumespresso.projwfc
  quantumespresso.pw
  quantumespresso.pw2wannier90
  quantumespresso.pwimmigrant
  quantumespresso.q2r
  simpleplugins.arithmetic.add
  simpleplugins.templatereplacer
  work
Default calculation input plugin: quantumespresso.pw
Installed on target computer? [True]: ?
Info: Whether the code is installed on the target computer or should be copied each time from a local path.
Installed on target computer? [True]:
Computer: daint
Remote absolute path: /users/dmarchan/Install/software/QuantumESPRESSO/6.3-backports-20181003-CrayIntel-18.08/bin/pw.x

verdi code list

11:03:40 glensk@lammm@cosmopc18 /local/scratch/glensk
%source ~/aiida/bin/activate
(aiida) 11:06:36 glensk@lammm@cosmopc18 /local/scratch/glensk
%verdi data upf uploadfamily --help
Usage: verdi data upf uploadfamily [OPTIONS] FOLDER GROUP_NAME
                                   GROUP_DESCRIPTION

  Upload a new pseudopotential family.

  Returns the numbers of files found and the number of nodes uploaded.

  Call without parameters to get some help.

Options:
  --stop-if-existing  Interrupt pseudos import if a pseudo was already present
                      in the AiiDA database
  -h, --help          Show this message and exit.
(aiida) 11:07:02 glensk@lammm@cosmopc18 /local/scratch/glensk
%mv pseudo_SSSPefV1.1_alalloy ~/Dropbox/Albert/scripts/
(aiida) 11:08:49 glensk@lammm@cosmopc18 /local/scratch/glensk
%verdi data upf uploadfamily pseudo_SSSPefV1.1_alalloy 'SSSP_v1.1_eff_custom' 'SSSP v1.1 efficiency with custom Cu pseudo'




verdi process list
verdi process list -a

iverdi calculation show 106
verdi profile show
verdi calculation inputcat 106 aiida.in
verdi calculation outputls
launch_calculation_pw -c pw-v6.3@daint -p SSSP_v1.1_eff_custom -s 87 -d
verdi shell
verdi work plugins # shows workflows

launch_calculation_pw -c pw-v6.3@daint -p SSSP_v1.1_eff_custom -s 87 -d
launch_workflow_pw_base -c pw-v6.3@daint -p SSSP_v1.1_eff_custom -s 87 -d

 1704  xyz2runner_mode1.sh
 1705  wc -l symfun.output
 1707  dropbox_.py status
 1708  cd Downloads/
 1709  which verdi
 1710  which aiida
 1711  /scratch/glensk/raw_data/test_3_at_1000K_for_fps
 1712  cd runner_scratch/
 1713  aiida
 1715  cd co
 1716  cd code
 1717  cd ai
 1718  cd aiida
 1719  git clone https://github.com/aiidateam/aiida-quantumespresso
 1720  cd aiida-quantumespresso
 1721  pip install -e .
 1722  verdi calculation plugins
 1724  vim aiida_quantumespresso/cli/calculations/pw/base.py
 1726  git branch
 1727  -
 1728  mv /local/scratch/glensk/runner_scratch ~/
 1729  ls ~/
 1730  untargz pseudo_SSSPefV1.1_alalloy.tar.gz
 1731  source ~/aiida/bin/activate
 1733  verdi data upf uploadfamily pseudo_SSSPefV1.1_alalloy 'SSSP_v1.1_eff_custom' 'SSSP v1.1 efficiency with custom Cu pseudo'
 1734  verdi data upf listfamilies
 1737  verdi data structure list
 1738  launch_calculation_pw --help
 1740  ssh daint -c 'sbatch /scratch/snx3000/aglensk/aiida_scratch/fb/10/56f1-8185-44b2-a583-2d369f31f201/_aiidasubmit.sh'
 1741  ssh --help
 1742  which ssh
 1743  ssh daint 'sbatch /scratch/snx3000/aglensk/aiida_scratch/fb/10/56f1-8185-44b2-a583-2d369f31f201/_aiidasubmit.sh'
 1744  verdi calculation gotocomputer 91
 1745  cd ../aiida-core
 1746  verdi devel setproperty logging.aiida_loglevel DEBUG
 1750  verdi calculation logshow 91
 1751  launch_calculation_pw -c pw-v6.3@daint -p SSSP_v1.1_eff_custom -s 87
 1754  verdi calculation list
 1755  /ssh daint "cd '/scratch/snx3000/aglensk/aiida_scratch/dc/04/71a5-40da-4226-96af-5eae749aacb1' && sbatch '_aiidasubmit.sh'"
 1756  verdi calculation gotocomputer 101
 1757  unalias ssh
 1758  ssh daint "cd '/scratch/snx3000/aglensk/aiida_scratch/dc/04/71a5-40da-4226-96af-5eae749aacb1' && sbatch '_aiidasubmit.sh'"
 1759  /local/scratch/glensk/runner_scratch
 1760  rm pseudo_SSSPefV1.1_alalloy.tar.gz
 1761  iterm_title runner mode 1
 1762  cd Dropbox
 1763  ..
 1764  cd glensk
 1765  diff input.nn $dotfiles/scripts/runner_scripts/input.nn
 1766  vimdiff input.nn $dotfiles/scripts/runner_scripts/input.nn
 1767  goo sed.commands
 1768  cat ~/scripts/dotfiles/commands/sed.commands
 1769  ls ~/google_drive
 1770  ls ~/google_drive/
 1771  ls ~/google_drive/KMC
 1772  ls ~/google_drive/KMC/
 1773  ls /local/scratch/glensk/google_drive_from_dropbox_to_delete/
 1774  ls /local/scratch/glensk/google_drive_from_dropbox_to_delete/KMC
 1775  cd
 1776  rm *
 1777  rm -rf save
 1778  ls /local/scratch/glensk
 1779  mv runner_scratch /local/scratch/glensk
 1780  cd /local/scratch/glensk/runner_scratch/
 1781  rm debug.out gmon.out
 1783  cp -r runner_scratch/ runner_scratch_tmp
 1784  cp -r runner_scratch/ runner_scratch_loesch
 1785  screen -r
 1786  RuNNer.serial.x > logfile_mode1.1
 1787  rm logfile_mode1.1
 1788  ./RuNNer.serial.cosmopc.natascha.x > logfile_mode1.1
 1789  tail -f logfile_mode1.1
 1790  ps aux | grep RuNNer
 1791  grep nelem input.nn
 1792  grep number_of_elements input.nn
 1793  grep "number_of_elements" input.nn
 1794  cp input.nn ../runner_scratch/
 1795  cp input.nn ../runner_scratch_tmp
 1796  ./RuNNer.serial.cosmopc.natascha.x
 1797  rm -rf runner_scratch_loesch/
 1798  cd runner_scratch
 1799  vi input.nn
 1800  grep begin input.data.all | wc -l
 1801  cd google_drive
 1802  cd KMC
 1803  cd results
 1804  cd 20171114_Al6xxxDB_NN_training_runner/
 1805  cd training-set
 1806  cd RuNNer
 1807  grep -c begin input.nn
 1808  grep -c begin training_input.data
 1809  mv runner_scratch runner_scratch_new_data
 1810  cp -r runner_scratch_new_data runner_scratch_old_data
 1811  /local/scratch/glensk
 1812  cd runner_scratch_new_data
 1813  cp -r runner_scratch_new_data_all
 1814  cp -r /runner_scratch_old_data runner_scratch_new_data_all
 1815  cp -r runner_scratch_old_data/ runner_scratch_new_data_all
 1816  cd runner_scratch_new_data_all
 1817  grep input.data input.nn
 1818  mv input.data.all input.data
 1819  cd runner_scratch_old_data/
 1820  rm input.data*
 1821  cd ..
 1822  cd 20181003_Al6xxxDB_recalculate_formation_energies/
 1823  cd 20181017_Al6xxxDB_NN_training_runner_5000inputstructures_from_giulio/
 1824  vi fps-5000.data
 1825  cp /local/scratch/glensk/google_drive/KMC/results/20181017_Al6xxxDB_NN_training_runner_5000inputstructures_from_giulio/fps-5000.data input.data
 1826  ./RuNNer.serial.cosmopc.natascha.x > logfile_mode1.1&
 1827  gd
 1828  ecgd
 1829  dot
 1830  gs
 1831  cd scratch
 1832  cd /local/scratch/glensk
 1833  ssh fidis
 1834  zsh
 1835  hostname
 1836  echo cosmopc18 | sed 's/\..*$//'
 1837  echo cmmd018 | sed 's/\..*$//'
 1838  echo cmmd018 | sed 's/\..*$//' | sed -r 's/[0-9]{1,10}$//'
 1839  bash
 1840  top
 1841  c
 1842  which c
 1843  cd runner_scratch_new_data/
 1844  grep -c begin input.data.all
 1846  runner_scratch_new_data_all/
 1847  du -sh *
 1848  grep -c begin input.data
 1849  python CurSel.py -t 1e-3 --landmarks `grep -c begin input.data` function.data logfile_mode1.1
 1850  iterm_title fps kmc
 1851  CurSel.py -h
 1852  ll -rtl
 1853  du -sh cursel*
 1854  wc -l cursel*
 1855  CurSel.py -t 1e-3 --landmarks `grep -c begin input.data` function.data logfile_mode1.1
 1856  mcd tmp
 1857  cp ../* .
 1858  CurSel.py -t 1e-3 --landmarks 200 function.data logfile_mode1.1
 1860  cd ../
 1861  frame_selector.py -h
 1862  cp cursel.landmarks cursel.landmarks.all
 1863  vi cursel.landmarks
 1864  frame_selector.py input.data precomp cursel.landmarks
 1865  ls -rtl
 1866  vi input_selected.data
 1867  grep -c input_selected.data
 1868  grep -c begin input_selected.data
 1869  gu
 1870  vi $scripts/runner_scripts/xyz2runner_mode1.sh
 1871  exit
 1873  htop
 1874  cd code/aiida/
 1875  cd aiida-core
 1876  vim aiida/transport/plugins/local.py
 1878  ssh daintr
 1880  ssh ela
 1892  verdi process list -a -p2
 1901  vim aiida/scheduler/plugins/slurm.py
 1905  vim aiida/scheduler/__init__.py
 1906  verdi daemon restart --reset
 1907  verdi calculation gotocomputer 106
 1909  \ssh daint
 1911  ssh daiont
 1912  ssh daint
 1914  verdi process list
 1915  daint
 1916  verdi devel setproperty logging.aiida_loglevel INFO
 1917  verdi daemon restart
 1918  launch_calculation_pw -c pw-v6.3@daint -p SSSP_v1.1_eff_custom -s 87 -d
 1919  pwd
 1920  git status
 1921  git diff .
 1922  git checkout -- .
 1924  verdi node delete 91 96 101
 1925  verdi process list -a
 1926  verdi calculation show 106
 1928  verdi calculation inputls 106
 1929  verdi calculation inputcat 106 aiida.in
 1930  verdi calculation inputcat 106 aiida.in | less
 1931  verdi calculation inputcat 106 aiida.in | vi
 1932  verdi calculation inputcat 106 aiida.in | vi -
 1933  verdi calculation outputls 106
 1934  verdi calculation outputcat 106 aiida.out | less
 1935  verdi data parameter show 110
 1936  verdi calculation
 1937  verdi calculation logshow 106
 1938  history
 1939  history -100
 1943* verdi work plugins
 1944* verdi daemon logshow
 1945* launch_workflow_pw_base -c pw-v6.3@daint -p SSSP_v1.1_eff_custom -s 87 -d
 1947  history -100 | grep launch
 1948* verdi work show 127
 1949* verdi work report 127
 1950* verdi work status 127
 1952  source aiida/bin/activate
 1953  verdu work status 138
 1954  verdi work status 138
 1957* verdi daemon stop
 1958* verdi profile show
 1959* rm -rf /home/glensk/aiida/.aiida/repository-al6xxx/repository
 1960* sudo su - postgres
 1962  cat ~/scripts/dotfiles/commands/cosmo_cosmopc_setup
 1964* ssh cosmopc
 1965* su local
 1973* verdi computer setup
 1976* verdi code setup
 1978* verdi export create --help
 1979* verdi export create -Y localhost daint -C pw-v6.3 -- base_computers_code.aiida
 1980* verdi export create -Y localhost daint -X pw-v6.3 -- base_computers_code.aiida
 1982* ls -al
 1983* ls
 1985* verdi profile setdefault al6xxx
 1988* verdi import base_computers_code.aiida
 1989* verdi computer list
 1990* verdi computer list -A
 1991* verdi computer list -a
 1992* verdi computer configure localhost
 1993* verdi computer configure local localhost
 1994* verdi computer configure ssh daint
 1995* verdi computer test localhost
 1996* verdi computer test daint
 1997* verdi code list
 1999* verdi profile setdefault al6xxx_test
 2000* ll
 2001  /s
 2002  s
 2003  untargz pseudo_SSSPefV1.1_alalloy
 2004* verdi data upf uploadfamily --help
 2005* verdi data upf uploadfamily /local/scratch/glensk/pseudo_SSSPefV1.1_alalloy SSSP_v1.1_eff 'SSSP v1.1 efficiency PBE with updated Cu pseudo'
 2006* verdi shell
 2007* verdi daemon status
 2008* verdi daemon start
 2009* launch_workflow_pw_base -c pw-v6.3@daint -p SSSP_v1.1_eff -s 87 -d
 2011* verdi profile list
 2012* verdi process list -a -p1
 2013  history -200
 2014  history -200 > ~/scripts/dotfiles/commands/cosmo_aiida_cosmopc