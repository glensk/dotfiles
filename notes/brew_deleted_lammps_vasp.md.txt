11:07:37 glensk@mac /Users/glensk
%brew doctor
Please note that these warnings are just used to help the Homebrew maintainers
with debugging if you file an issue. If everything you use Homebrew for is
working fine: please don't worry and just ignore them. Thanks!

Warning: "config" scripts exist outside your system or Homebrew directories.
`./configure` scripts often look for *-config scripts to determine if
software packages are installed, and what additional flags to use when
compiling and linking.

Having additional scripts in your path can confuse software installed via
Homebrew if the config script overrides a system or Homebrew provided
script of the same name. We found the following "config" scripts:

  /Users/glensk/anaconda/bin/curl-config
  /Users/glensk/anaconda/bin/freetype-config
  /Users/glensk/anaconda/bin/libdynd-config
  /Users/glensk/anaconda/bin/libpng-config
  /Users/glensk/anaconda/bin/libpng16-config
  /Users/glensk/anaconda/bin/llvm-config
  /Users/glensk/anaconda/bin/python-config
  /Users/glensk/anaconda/bin/python2-config
  /Users/glensk/anaconda/bin/python2.7-config
  /Users/glensk/anaconda/bin/xml2-config
  /Users/glensk/anaconda/bin/xslt-config

Warning: Putting non-prefixed findutils in your path can cause python builds to fail.

Warning: Your default Python does not recognize the Homebrew site-packages
directory as a special site-packages directory, which means that .pth
files will not be followed. This means you will not be able to import
some modules after installing them with Homebrew, like wxpython. To fix
this for the current user, you can run:

  mkdir -p /Users/glensk/.local/lib/python2.7/site-packages
  echo 'import site; site.addsitedir("/usr/local/lib/python2.7/site-packages")' >> /Users/glensk/.local/lib/python2.7/site-packages/homebrew.pth

Warning: Unbrewed dylibs were found in /usr/local/lib.
If you didn't put them there on purpose they could cause problems when
building Homebrew formulae, and may need to be deleted.

Unexpected dylibs:
    /usr/local/lib/libmpi.12.dylib
    /usr/local/lib/libmpicxx.12.dylib
    /usr/local/lib/libmpifort.12.dylib
    /usr/local/lib/libpmpi.12.dylib

Warning: Unbrewed header files were found in /usr/local/include.
If you didn't put them there on purpose they could cause problems when
building Homebrew formulae, and may need to be deleted.

Unexpected header files:
    /usr/local/include/mpi.h
    /usr/local/include/mpicxx.h
    /usr/local/include/mpif.h
    /usr/local/include/mpio.h
    /usr/local/include/mpiof.h
    /usr/local/include/opa_config.h
    /usr/local/include/opa_primitives.h
    /usr/local/include/opa_queue.h
    /usr/local/include/opa_util.h
    /usr/local/include/primitives/opa_by_lock.h
    /usr/local/include/primitives/opa_emulated.h
    /usr/local/include/primitives/opa_gcc_ia64.h
    /usr/local/include/primitives/opa_gcc_intel_32_64.h
    /usr/local/include/primitives/opa_gcc_intel_32_64_barrier.h
    /usr/local/include/primitives/opa_gcc_intel_32_64_ops.h
    /usr/local/include/primitives/opa_gcc_intel_32_64_p3.h
    /usr/local/include/primitives/opa_gcc_intel_32_64_p3barrier.h
    /usr/local/include/primitives/opa_gcc_intrinsics.h
    /usr/local/include/primitives/opa_gcc_ppc.h
    /usr/local/include/primitives/opa_gcc_sicortex.h
    /usr/local/include/primitives/opa_nt_intrinsics.h
    /usr/local/include/primitives/opa_sun_atomic_ops.h
    /usr/local/include/primitives/opa_unsafe.h
    /usr/local/include/ublio.h

Warning: Unbrewed .la files were found in /usr/local/lib.
If you didn't put them there on purpose they could cause problems when
building Homebrew formulae, and may need to be deleted.

Unexpected .la files:
    /usr/local/lib/libmpi.la
    /usr/local/lib/libmpicxx.la
    /usr/local/lib/libmpifort.la
    /usr/local/lib/libpmpi.la

Warning: Unbrewed .pc files were found in /usr/local/lib/pkgconfig.
If you didn't put them there on purpose they could cause problems when
building Homebrew formulae, and may need to be deleted.

Unexpected .pc files:
    /usr/local/lib/pkgconfig/mpich.pc
    /usr/local/lib/pkgconfig/openpa.pc

Warning: Unbrewed static libraries were found in /usr/local/lib.
If you didn't put them there on purpose they could cause problems when
building Homebrew formulae, and may need to be deleted.

Unexpected static libraries:
    /usr/local/lib/libmpi.a
    /usr/local/lib/libmpicxx.a
    /usr/local/lib/libmpifort.a
    /usr/local/lib/libpmpi.a

Warning: You have unlinked kegs in your Cellar
Leaving kegs unlinked can lead to build-trouble and cause brews that depend on
those kegs to fail to run properly once built. Run `brew link` on these:

    lesstif
    open-mpi
