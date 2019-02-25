ist die integerbreite in den bibliotheken: blas lapac fft bleich der integerbreite in vasp?
otool -L extract_forceconstants :to check for links, like ldd:


############################################
cp $MKLROOT/include/fftw/fftw3.f .
und nun in garching:
cp ../vasp_4.6/fftw3.f .
###############################################

change:
fftmpiw.F and make INTEGER :: planx, plany, planz   to INTEGER(8) :: planx, plany, planz





garching:  # DAS IST NUN FALSCH
module load intel mkl impi   # -> laed echo $LD_LIBRARY_PATH/   -> /afs/@cell/common/soft/intel/mkl/lib/intel64//    # DAS IST NUN FALSCH
module load vasp/4.6.28      # -> laed echo $LD_LIBRARY_PATH/   -> /afs/@cell/common/soft/intel/mkl/lib/intel64//    # DAS IST NUN FALSCH


# DAS IST RICHTIG
garching now:
module load intel mkl impi:  macht aus dem LD_LIBRARY_PATH
/afs/@cell/common/soft/intel/ics2013/14.0/compiler/lib/intel64/
