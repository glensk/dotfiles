#!/bin/sh
echo ""
echo "n2p2 should be compiled before compiling lammps!"
echo ""

installfolder=$HOME  # folder where n2p2 will be installed
[ "$USER" ==  "glensk" ] && installfolder=$HOME/Dropbox/Albert/git

#### load modules
#conda deactivate
source $MODULESHOME/init/bash
module purge
module load intel
module load intel-mpi
module load intel-mkl
module load fftw
module load python/2.7.14
module load gsl
module load eigen

#### clone n2p2 to installfolder
cd $installfolder
echo "#########################"
echo installfolder: `pwd`
echo "#########################"
[ ! -e "n2p2" ] && git clone https://github.com/CompPhysVienna/n2p2.git

cd n2p2
git checkout develop
git branch
cd src

# adapt makefile
cp makefile makefile.back
sed -i 's|^COMP=.*|COMP=intel|' makefile
sed -i 's|^PROJECT_DIR.*|PROJECT_DIR=./|' makefile
sed -i 's|^LIB=libnnp.so libnnpif.so libnnptrain.so pynnp.so|LIB=libnnp.so libnnpif.so libnnptrain.so|' makefile  # remove pynnp.so

# adapt makefile.intel
cp makefile.intel makefile.intel.back
sed -i 's|^PROJECT_GSL=.*|PROJECT_GSL=${GSL_ROOT}/include|' makefile.intel
sed -i 's|^PROJECT_EIGEN=.*|PROJECT_EIGEN=${EIGEN_ROOT}/include/eigen3|' makefile.intel
line_eigen=`grep -n PROJECT_EIGEN makefile.intel | sed 's|:.*||'`
sed -i ''"$line_eigen"'i\MKL_INCLUDE=${MKLROOT}/include' makefile.intel
sed -i 's|^PROJECT_LDFLAGS_BLAS=.*|PROJECT_LDFLAGS_BLAS=-L${GSL_ROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl|' makefile.intel

# adapt libnnptrain/makefile
cp libnnptrain/makefile libnnptrain/makefile.back 
sed -i 's|^INCLUDES=.*|INCLUDES=-I./ -I${PROJECT_INCLUDE}/ -I${PROJECT_GSL} -I${PROJECT_EIGEN} -I${MKL_INCLUDE}|' libnnptrain/makefile

# was added recently
make libnnpif-shared

# make
make

# on cosmopc:
# in my case /usr/include/eigen3/Eigen/src/Core/util/MKL_support.h has in line 56 a preprocessor command #if defined EIGEN_USE_MKL
# the EIGEN_USE_MKL is actually not acitvaed in the makefiel but is somehow, during the compilation set (which makes the compilation crash)
