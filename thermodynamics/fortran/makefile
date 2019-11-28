FLAGS=-free  #ifort
#FLAGS=ifort
#FLAGS=-ffree-form -ffixed-line-length-none #gfortran

FC=ifort    # for mac
FCG=ifort   # for mac

#FC=gfortran  # for cmpc (at cmpc getGibbsEnergyOfFormation will not be compiled correctly)
#FCG=gfortran # for cmpc (at cmpc getGibbsEnergyOfFormation will not be compiled correctly)


FLAGS=-free

inc=utils.o constants.o free_energy_surface.o

lapack=-L//afs/@cell/common/soft/intel/ics2017.4/17.0/linux/mkl/lib/intel64/ \
	-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
lapack=/usr/lib/liblapack.dylib ~/BLAS/blas_MACOS.a # gfortran on apple

all=toReduced toCartesian dos dos2d getFel getThermodynamics \
    getJump getStacking norm dotProduct removeMappingInMD debye \
    get_d_from_POSITIONs get_dF_from_POSITIONs getGibbsEnergyOfFormation

all: $(all)

clean:
	rm -f *.o *.mod *.x

removeMappingInMD:   # for cmpc01
	$(FCG) invert.f90 $@.f90 -o $@.x

toReduced: $(inc)
	$(FC) $(FLAGS) invert.f90 $@.f90 -o $@.x

toCartesian: $(inc)
	$(FC) $(FLAGS) $@.f90 -o $@.x

dotProduct: $(inc)
	$(FC) $(FLAGS) $@.f90 -o $@.x

debye: $(inc)
	$(FC) $(FLAGS) spline.f90 $@.f90 -o $@.x

dos: $(inc)
	$(FC) $(FLAGS) quicksort.f90 $@.f90 -o $@.x

dos2d: $(inc)
	$(FC) $(FLAGS) quicksort.f90 $@.f90 -o $@.x

getFel: $(inc)
	$(FC) $(FLAGS) quicksort.f90 $@.f90 -o $@.x

getEinstein: $(inc)
	$(FC) $(FLAGS) quicksort.f90 $@.f90 -o $@.x

getThermodynamics: $(inc)
	$(FC) $(FLAGS) $(inc) $@.f90 -o $@.x

getGibbsEnergyOfFormation: $(inc)
	$(FC) $(FLAGS) $(inc) $(lapack) $@.f90 -o $@.x

getJump:
	$(FC) $(FLAGS) invert.f90 $@.f90 -o $@.x

getStacking:
	$(FC) $(FLAGS) invert.f90 $@.f90 -o $@.x

norm:
	$(FC) $(FLAGS) $@.f90 -o $@.x

get_d_from_POSITIONs:
	$(FC) $(FLAGS) $@.f90 -o $@.x

get_dF_from_POSITIONs:
	$(FC) $(FLAGS) $@.f90 -o $@.x

%.o: %.f90
	$(FC) $(FLAGS) -c $<

