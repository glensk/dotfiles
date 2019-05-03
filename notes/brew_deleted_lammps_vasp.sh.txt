#!/bin/sh


# Unexpected header files:
rm    /usr/local/include/mpi.h
rm    /usr/local/include/mpicxx.h
rm    /usr/local/include/mpif.h
rm    /usr/local/include/mpio.h
rm    /usr/local/include/mpiof.h
rm    /usr/local/include/opa_config.h
rm    /usr/local/include/opa_primitives.h
rm    /usr/local/include/opa_queue.h
rm    /usr/local/include/opa_util.h
rm    /usr/local/include/primitives/opa_by_lock.h
rm    /usr/local/include/primitives/opa_emulated.h
rm    /usr/local/include/primitives/opa_gcc_ia64.h
rm    /usr/local/include/primitives/opa_gcc_intel_32_64.h
rm    /usr/local/include/primitives/opa_gcc_intel_32_64_barrier.h
rm    /usr/local/include/primitives/opa_gcc_intel_32_64_ops.h
rm    /usr/local/include/primitives/opa_gcc_intel_32_64_p3.h
rm    /usr/local/include/primitives/opa_gcc_intel_32_64_p3barrier.h
rm    /usr/local/include/primitives/opa_gcc_intrinsics.h
rm    /usr/local/include/primitives/opa_gcc_ppc.h
rm    /usr/local/include/primitives/opa_gcc_sicortex.h
rm    /usr/local/include/primitives/opa_nt_intrinsics.h
rm    /usr/local/include/primitives/opa_sun_atomic_ops.h
rm    /usr/local/include/primitives/opa_unsafe.h
rm    /usr/local/include/ublio.h



# Unexpected dylibs:
rm     /usr/local/lib/libmpi.12.dylib
rm     /usr/local/lib/libmpicxx.12.dylib
rm     /usr/local/lib/libmpifort.12.dylib
rm     /usr/local/lib/libpmpi.12.dylib

# Unexpected .la files:
rm    /usr/local/lib/libmpi.la
rm    /usr/local/lib/libmpicxx.la
rm    /usr/local/lib/libmpifort.la
rm    /usr/local/lib/libpmpi.la


#Unexpecte d.pc files:
rm    /usr/local/lib/pkgconfig/mpich.pc
rm    /usr/local/lib/pkgconfig/openpa.pc


#Unexpected static libraries:
rm    /usr/local/lib/libmpi.a
rm    /usr/local/lib/libmpicxx.a
rm    /usr/local/lib/libmpifort.a
rm    /usr/local/lib/libpmpi.a



#Warning: Broken symlinks were found. Remove them with `brew prune`:
brew prune  /usr/local/lib/libfmpich.dylib
brew prune  /usr/local/lib/libmpi.dylib
brew prune  /usr/local/lib/libmpich.dylib
brew prune  /usr/local/lib/libmpichcxx.dylib
brew prune  /usr/local/lib/libmpichf90.dylib
brew prune  /usr/local/lib/libmpicxx.dylib
brew prune  /usr/local/lib/libmpifort.dylib
brew prune  /usr/local/lib/libmpl.dylib
brew prune  /usr/local/lib/libopa.dylib
brew prune  /usr/local/lib/libpmpi.dylib
