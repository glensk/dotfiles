###########################################################################
###########################################################################
for vasp 4.6 (cmmd)
###########################################################################
###########################################################################
old: 
$setenv | grep MKL
MKLROOT=/opt/intel/Compiler/11.1/073/mkl

- 1st compile lib_i4
    rm -f vasp *.o libdmy.a
    module load intelcompiler/11.1 openmpi/1.6/intelcompiler/11.1
    make

- little script to install automatically
tar -xvf lib_i4.tar.gz
installfolder=sources_tdi
source=$installfolder.tar.gz
makefile=`pwd`/makefile.tdi/makefile.parallel.cmmd002.mpie.de
[ ! -e "$makefile" ] && echo makefile $makefile does not exist && exit 
rm -rf $installfolder
tar -xvf $source
cd $installfolder
cp $makefile makefile
make


what might improve the performance on cmmd or mac is the avx flag in the makefile



###########################################################################
###########################################################################
for vasp 5.2.12 (cmmd)
###########################################################################
###########################################################################
cd vasp.5.lib
- in in diesem folder die .o dateien loschen da auc linpack.o oder so gebraucht wird! (ande o files aber nicht)
- if fftw/3.3.2 module is really necessary I am not sure
in vasp.5.lib einfach make machen

module load intelcompiler/11.1 openmpi/1.6/intelcompiler/11.1 fftw/3.3.2/intelcompiler/11.1
make


sp.5.lib.tar.gz
installfolder=sources_5.2.12
source=$installfolder.tar.gz
makefile=`pwd`/makefile.`hostname`
[ ! -e "$makefile" ] && echo makefile $makefile does not exist && exit
rm -rf $installfolder
tar -xvf $source
cd $installfolder
cp $makefile .
make






what might improve the performance on cmmd or mac is the avx flag in the makefile


vi src/lib/getshmem.c



###########################################################################
###########################################################################
for vasp 5.4.1  (for cmmd)
###########################################################################
###########################################################################
$ifort -V
Intel(R) Fortran Intel(R) 64 Compiler Professional for applications running on Intel(R) 64, Version 11.1    Build 20100806 Package ID: l_cprof_p_11.1.073
Copyright (C) 1985-2010 Intel Corporation.  All rights reserved.
--->>>>> is this evtl tooo old??




###########################################################################
###########################################################################
for vasp 5.4.1  (for mac)
###########################################################################
###########################################################################
a) make sure you have installed the icc (intel c compiler): icc -V
b) create fresh folder
c) cp vasp.5.4.1.tar.gz to this folder and uncompress vasp and cd to vasp.5.4.1/ folder
d) cp patch.5.4.1.08072015 .
e) apply the patch to vasp 5.4.1:
patch -p1 < patch.5.4.1.08072015
neu hinzugefuegt f) cp /Users/glensk/Dropbox/scripts/vasp_sources/vasp.5.4.1/makefile.include.mac.works makefile.include
                - dann musste ich noch aus dem build folder was 
                - dann zeilen im code ausgetauscht (getshmem.c siehe unten)

f) install mpich with following steps:
- first install mpich (NOT WITH BREW!, at least not for vasp compilation)
  downloaded from https://www.mpich.org/downloads/
- untargz mpich
- cd mpich
- export CC=icc       # intel c compiler
- export CXX=icc
- export FC=ifort     # intel fortran compiler
- export F77=ifort
export CC=icc;export CXX=icc;export FC=ifort;export F77=ifort
# - ./configure   # or better
- mkdir mac_conf
- ./configure 
- ./configure --prefix=/Users/glensk/Dropbox/scripts/vasp_sources/vasp.5.4.1/mpich-3.1.4/mac_conf  # did not write anything to mac_conf
- ./configure --prefix=/Users/glensk/local 
- make          # 
- sudo make install  # dont sudo make install! breaks lammps with brew! 
                # creates mpirun in mac_conf/bin/mpirun

source /opt/intel/composer_xe_2015.1.108/mkl/bin/mklvars.sh intel64  # sets DYLD

% mpif90 -V  
mpifort for MPICH version 3.1.4
ifort version 15.0.1

% mpicc -V
Intel(R) C Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 15.0.1.108 Build 20141022
Copyright (C) 1985-2014 Intel Corporation.  All rights reserved.

%ifort -V
Intel(R) Fortran Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 15.0.1.108 Build 20141022
Copyright (C) 1985-2014 Intel Corporation.  All rights reserved.


%env | grep MKL

# folgende 4 lines (ging nicht)
# cd /Users/glensk/Dropbox/scripts/vasp_sources/vasp.5.4.1
# cp /opt/intel/composer_xe_2015.1.108/mkl/interfaces/fftw3xf .
# make libintel64

d) compile the fftw3xf
cd /opt/intel/composer_xe_2015.1.108/mkl/interfaces/fftw3xf
make libintel64           # or at make lib64
sudo make libintel64      # or at 

e) change getshmem.c file
then cd vasp folder  (cd ..)
vi src/lib/getshmem.c
# now set in line 6:
define SHM_NORESERVE 0
#

f) compile vasp!  (what might improve the performance on cmmd or mac is the avx flag in the makefile)

g) set LD_LIBRARY_PATH to the folder where mpich was compiled to (z.B. /Users/glensk/local/lib)

h) otool -L ~/scripts/vasp_sources/vasp.5.4.1/bin/vasp_std

###########################################################################
for vasp 5.4.1  (for cmmc)
###########################################################################
patch -p1 < patch.5.4.1.08072015
module load impi/4.1.3
module load intel/14.0
module load mkl/11.1
- in tcsh set:
setenv LD_LIBRARY_PATH "$HOME/lib:/afs/@cell/common/soft/intel/ics2013/14.0/compiler/lib/intel64:/afs/@cell/common/soft/intel/ics2013/14.0/mkl/lib/intel64/"
- ankit had to create the bin folder by hand

################## --> OUTPUT of versions
mpiifort -V     (and)   ifort -V
Intel(R) Fortran Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 14.0.4.211 Build 20140805
Copyright (C) 1985-2014 Intel Corporation.  All rights reserved.

GNU ld (GNU Binutils; SUSE Linux Enterprise 11) 2.23.1
/afs/ipp-garching.mpg.de/common/soft/intel/ics2013/composer_xe_2013_sp1.4.211/compiler/lib/intel64/for_main.o: In function `main':
/export/users/nbtester/efi2linux_nightly/branch-14_0/20140806_000000/libdev/frtl/src/libfor/for_main.c:(.text+0x42): undefined reference to `MAIN__'


###########################################################################
###########################################################################
stefan used:
###########################################################################
###########################################################################
mpif90 -V
mpicc -V
ifort -V


setenv
setenv | grep MKL
cd /opt/intel/Compiler/11.1/073/mkl/interfaces/
cd fftw3xf/
make
make libintel64
make intel
cd ../
cp -R fftw3xf/ ~/
cd ~/fftw3xf/
make lib64

mpichversion
icpc -V


###########################################################################
to start vasp:
###########################################################################
mpirun -np 20 /u/system/SLES11/soft/vasp5.4/5.4.1/intel-14.0/impi-4.1/bin/vasp
mpirun -np 20 ~/Thermodynamics/vasp_Langevin/vasp_5.3.5/vasp-stable-cmmc-24.06.2016
