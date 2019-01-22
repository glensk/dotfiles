#!/bin/sh
#########################################################################
# define system dependent (and system independent) path variables
# e.g. $PATH / $PYTHONPATH / $LD_LIBRARY_PATH
# all prior defined environment variables are known
#########################################################################


#########################################################################
# on every system: ( $PATH / $PYTHONPATH / $LD_LIBRARY_PATH / ...)
#########################################################################
#anacondaold="/usr/local/bin"   # use this to install macvim @mac otherwise interference with homebrew python
#anaconda="$HOME/anaconda/bin"  # don't load anaconda on cmmc since the module is available
svnctags="/usr/local/bin"   # svn / ctags / ifort / icc  (ifort = /opt/intel/bin/ifort)
homebrew="/usr/local/sbin"
bin="/usr/bin"              # on mac it is there anyhow.
bin6="$HOME/.local/bin"


dropboxbin="$dotfiles/bin"       # Dropbox binaries
bin2="$dotfiles/sources_bin"
bin3="$dotfiles/bin/stefan"
bin5="$dotfiles/bin/phonon_lifetimes"
bin6="$dotfiles/aliases"  

export PATH="$PATH:$dropboxbin:$svnctags:$homebrew:$bin:$bin2:$bin3:$bin4:$bin5:$bin6"
# ifrot / icc / mpif90 / mpicc / (mpirun) are installed both in /usr/local/bin and ~/local/bin !


#########################################################################
# on mac: ( $PATH / $PYTHONPATH / $LD_LIBRARY_PATH / ... )
#########################################################################
if [ "$onhost" = "mac" ];then
    # $PATH (on mac)
    sed="/usr/local/opt/gnu-sed/libexec/gnubin"
    
    # !! too much time !!! brew="$(brew --prefix coreutils)/libexec/gnubin" 
    brew="/usr/local/opt/coreutils/libexec/gnubin"   # without brew ls options -- are not working
    phonopy="$HOME/scripts/phonons/phonopy_at_$host"
    #tdep="$HOME/Dropbox/Albert/scripts/phonons/tdep-devel/bin"
    #phonopybin="$phonopy/bin"
    #phonopybin="$HOME/.local/bin"
    #phonopybin2="$HOME/.miniconda2/bin"
    #phonopybin2=""
    #vasp="$HOME/local/bin"   # dont confuse this path with /usr/local/bin where also icc and ifort are installed!
    #vimctags="$HOME/Dropbox/scripts/dotfiles/vim/ctags-5.8/installfolder/bin"
    #sphinx="$HOME/Dropbox/scripts/phonons/sphinx/bin"


    # (on mac)
    pygraceplot="$HOME/scripts/python/graceplot"
    pylammps1="/usr/local/Cellar/lammps/2014.02.12/lib/python2.7/site-packages"
    #pyphonopy1="$phonopy/lib/python"     # may be lib64 instead
    #pyphonopy2="$phonopy/lib64/python"   # may be lib64 instead
    pyphonopy1="$HOME/.local/lib/python2.7"     # may be lib64 instead
    pyphonopy2=""   # may be lib64 instead  no path needed now when using phono3py


    # $LD_LIBRARY_PATH add on mac
    ldvasp="$HOME/local/lib"
    ldphonopylapack="/usr/local/opt/lapack/lib"
    ldphonopyopenblas="/usr/local/opt/openblas/lib"
    #setenv LD_LIBRARY_PATH '/afs/@cell/common/soft/intel/mkl/lib/intel64/'
    #setenv LD_LIBRARY_PATH "$HOME/lib:/afs/@cell/common/soft/intel/ics2013/14.0/compiler/lib/intel64:/afs/@cell/common/soft/intel/ics2013/14.0/mkl/lib/intel64/"
    
    export PATH="$sed:$brew:$phonopy:$PATH"   # brew necessary to load coreutils
    export LD_LIBRARY_PATH="$ldvasp:$ldphonopylapack:$ldphonopyopenblas:$LD_LIBRARY_PATH"
    export C_INCLUDE_PATH="/Users/glensk/.miniconda3/include:$C_INCLUDE_PATH"
fi
