#!/bin/sh

##installfolder=$HOME  # folder where n2p2 will be installed
## on fidis do: intel/18.0.2   2) intel-mpi/2018.2.199   3) gsl/2.4   4) eigen/3.3.4
#module load intel
#exit
#module load intel-mpi
#module load intel-mkl
#module load gsl
#module load eigen
# gives intel/18.0.2   2) intel-mpi/2018.2.199   3) gsl/2.4   4) eigen/3.3.4
#
lammpsfolder="lammps_source_cosmo"
#####################################################################
# get n2p2 
#####################################################################
cd $HOME #$installfolder
[ "$USER" ==  "glensk" ] && cd $HOME/Dropbox/Albert/git
echo "#########################"
pwd
echo "#########################"
[ ! -e "n2p2" ] && git clone https://github.com/CompPhysVienna/n2p2.git
cd n2p2
git checkout develop
git branch
cd src

# adapt makefile
sed -i 's|^COMP=.*|COMP=intel|' makefile
sed -i 's|^PROJECT_DIR.*|PROJECT_DIR=./|' makefile
sed -i 's|^LIB=libnnp.so libnnpif.so libnnptrain.so pynnp.so|LIB=libnnp.so libnnpif.so libnnptrain.so|' makefile  # remove pynnp.so

# adapt makefile.intel
sed -i 's|^PROJECT_GSL=.*|PROJECT_GSL=${GSL_ROOT}/include|' makefile.intel
sed -i 's|^PROJECT_EIGEN=.*|PROJECT_EIGEN=${EIGEN_ROOT}/include/eigen3|' makefile.intel
line_eigen=`grep -n PROJECT_EIGEN makefile.intel | sed 's|:.*||'`
sed -i ''"$line_eigen"'i\MKL_INCLUDE=${MKLROOT}/include' makefile.intel
sed -i 's|^PROJECT_LDFLAGS_BLAS=.*|PROJECT_LDFLAGS_BLAS=-L${GSL_ROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl|' makefile.intel

# adapt libnnptrain/makefile
sed -i 's|^INCLUDES=.*|INCLUDES=-I./ -I${PROJECT_INCLUDE}/ -I${PROJECT_GSL} -I${PROJECT_EIGEN} -I${MKL_INCLUDE}|' libnnptrain/makefile

# make
make

# on cosmopc:
# in my case /usr/include/eigen3/Eigen/src/Core/util/MKL_support.h has in line 56 a preprocessor command #if defined EIGEN_USE_MKL
# the EIGEN_USE_MKL is actually not acitvaed in the makefiel but is somehow, during the compilation set (which makes the compilation crash)

#####################################################################
# now the part to make lammps work with n2p2
#####################################################################
exit
# from https://compphysvienna.github.io/n2p2/a01203.html
#pwd
#make libnnpif-shared

cd $HOME/sources/$lammpsfolder
ln -s $HOME/Dropbox/Albert/git/n2p2 lib/nnp
cp -r $HOME/Dropbox/Albert/git/n2p2/src/interface/LAMMPS/src/USER-NNP $HOME/sources/$lammpsfolder/src

cd $HOME/sources/$lammpsfolder/src
make yes-user-nnp
#make mpi
make fidis #(on fidis)

