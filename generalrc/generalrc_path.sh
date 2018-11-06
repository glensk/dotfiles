#!/bin/sh
#########################################################################
# define system dependent (and system independent) path variables
# e.g. $PATH / $PYTHONPATH / $LD_LIBRARY_PATH
# all prior defined environment variables are known
#########################################################################


#########################################################################
# on every system: ( $PATH / $PYTHONPATH / $LD_LIBRARY_PATH / ...)
#########################################################################
dropboxbin="$dotfiles/bin"       # Dropbox binaries
#anacondaold="/usr/local/bin"   # use this to install macvim @mac otherwise interference with homebrew python
#anaconda="$HOME/anaconda/bin"  # don't load anaconda on cmmc since the module is available
svnctags="/usr/local/bin"   # svn / ctags / ifort / icc  (ifort = /opt/intel/bin/ifort)
homebrew="/usr/local/sbin"
bin="/usr/bin"
bin2="$dotfiles/sources_bin"
bin3="$dotfiles/bin/stefan"
bin5="$dotfiles/bin/phonon_lifetimes"
bin6="$HOME/.local/bin"
ipi="$HOME/google_drive/scripts_epfl/i-pi-mc/bin"
kmc="$HOME/google_drive/Glensk/KMC-Tutorial/scripts"
kmc2="$dotfiles/scripts/ipi-kmc"
runnertools="$dotfiles/scripts/runner_scripts"
aliases="$dotfiles/aliases"
aiida="$dotfiles/scripts/qe-aiida/aiida_solutejobs_scripts"
aiida2="$dotfiles/scripts/qe-aiida/aiida_analyze"
lammps1="$dotfiles/scripts/lammps_executables"
lammps2="$dotfiles/scripts/lammps_scripts"

#PATH="$dropboxbin:bin3:$PATH"  !! not like this! $PATH first !!! 
PATH="$PATH:$aiida:$aiida2:$aliases:$lammps1:$lammps2:$kmc:$kmc2:$dropboxbin:$runnertools:$svnctags:$homebrew:$bin:$bin2:$bin3:$bin4:$bin5:$bin6"
# ifrot / icc / mpif90 / mpicc / (mpirun) are installed both in /usr/local/bin and ~/local/bin !


#pythermodynamics=
#########################################################################
# on mac: ( $PATH / $PYTHONPATH / $LD_LIBRARY_PATH / ... )
#########################################################################
if [ "$currenthost" = "onmac" ];then
    # $PATH (on mac)
    sed="/usr/local/opt/gnu-sed/libexec/gnubin"
    # !! too much time !!! brew="$(brew --prefix coreutils)/libexec/gnubin" 
    #anaconda="$HOME/anaconda2/bin"  # don't load anaconda on cmmc since the module is available
    #anaconda="$HOME/miniconda3/bin"  # don't load anaconda on cmmc since the module is available
    #anaconda=""  # don't load anaconda on cmmc since the module is available
    brew="/usr/local/opt/coreutils/libexec/gnubin" 
    #tdep="$HOME/Dropbox/Albert/scripts/phonons/tdep-devel/bin"
    phonopy="$HOME/scripts/phonons/phonopy_at_$host"
    #phonopybin="$phonopy/bin"
    #phonopybin="$HOME/.local/bin"
    #phonopybin2="$HOME/.miniconda2/bin"
    #phonopybin2=""
    vasp="$HOME/local/bin"   # dont confuse this path with /usr/local/bin where also icc and ifort are installed!
    #vimctags="$HOME/Dropbox/scripts/dotfiles/vim/ctags-5.8/installfolder/bin"
    #sphinx="$HOME/Dropbox/scripts/phonons/sphinx/bin"


    # (on mac)
    pygraceplot="$HOME/scripts/python/graceplot"
    pylammps1="/usr/local/Cellar/lammps/2014.02.12/lib/python2.7/site-packages"
    pylammps2="/usr/local/lib/python2.7/site-packages"  # this is too general an might import libs from non anaconda python 
    pylammps2=""                                        # this is too general an might import libs from non anaconda python 
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
    
    #[ "$host" = "$mylaptop" ] && \
    PATH="$sed:$brew:$phonopy:$vasp:$PATH" #&& \
    #PATH="$sed:$brew:$tdep:$phonopy:$vasp:$PATH" #&& \
    #PATH="$sed:$brew:$anaconda:$tdep:$phonopy:$phonopybin:$phonopybin2:$vasp:$PATH" #&& \
    #PYTHONPATH="$pygraceplot:$pylammps1:$pylammps2:$pyphonopy1:$PYTHONPATH" #&& \
    # for Miniconda2 use:
    #PYTHONPATH="/Users/glensk/.miniconda2"
    LD_LIBRARY_PATH="$ldvasp:$ldphonopylapack:$ldphonopyopenblas:$LD_LIBRARY_PATH"
fi


#########################################################################
# on cmmd: ( $PATH / $PYTHONPATH / $LD_LIBRARY_PATH / ... )
# cmmd does not exist anymore
#########################################################################
#if [ "$currenthost" = "oncmmd" ];then
#    
#    xmgracebin="/data/grabowski/xmgrace/bin/"
#    add1="/lib64"
#    PATH="$xmgracebin:$add1:$PATH" 
#
#    argparse="${HOME}/Dropbox/scripts/python/add_cmmd/argparse-1.2.1"
#    PYTHONPATH="$argparse:$PYTHONPATH" 
#
#    # ldgibbsform ist fuer getGibbsEnergyOfFormation.sh
#    ldgibbsform="/opt/intel/Compiler/11.1/073/mkl/lib/em64t/"   
#    ld1="/home/grabowski/libs/"
#    ld2="/home/glensk/"
#    ldxmgrace="/opt/intel/lib/intel64"
#    LD_LIBRARY_PATH="$ldgibbsform:$ld1:$ld2:$ldxmgrace:$LD_LIBRARY_PATH"
#fi
    
#########################################################################
# on cmmd001:  module load openmpi/1.5/intelcomp
#########################################################################


#########################################################################
# remove trailing / beginning ":" = (0.10s) CAN BE DROPPED IF EVERYTHING IS CORRECT! 
#########################################################################
#!! NOT NECESSARY !! PATH=`echo $PATH | sed 's|:$||' | sed 's|^:||'` 
#!! NOT NECESSARY !! PYTHONPATH=`echo $PYTHONPATH | sed 's|:$||' | sed 's|^:||'` 
#!! NOT NECESSARY !! LD_LIBRARY_PATH=`echo $LD_LIBRARY_PATH | sed 's|:$||' | sed 's|^:||'` 



#########################################################################
# remove duplicates from $PATH / $PYTHONPATH / $LD_LIBRARY_PATH = (0.10s)
# we let this be for to save the 0.10 seconds
#########################################################################
#echo echo $PATH
#PATH=`perl -e 'print join ":", grep {!$h{$_}++} split ":", $ENV{PATH}'`
#PATH=`perl -e 'print join ":", grep {!$h{$_}++} split ":", ${PATH}'`
#PATH="$(printf "%s" "${PATH}" | /usr/bin/awk -v RS=: -v ORS=: '!($0 in a) {a[$0]; print}')"
#PYTHONPATH=`perl -e 'print join ":", grep {!$h{$_}++} split ":", $ENV{PYTHONPATH}'`
#LD_LIBRARY_PATH=`perl -e 'print join ":", grep {!$h{$_}++} split ":", $ENV{LD_LIBRARY_PATH}'`

#PATH="$(echo $PATH | perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, scalar <>))')"
#PYTHONPATH="$(echo $PYTHONPATH | perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, scalar <>))')"
#LD_LIBRARY_PATH="$(echo $LD_LIBRARY_PATH | perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, scalar <>))')"

#PATH=$(echo "$PATH" | awk -v RS=':' -v ORS=":" '!a[$1]++')  DONT!! INCLUDES : at the end

#########################################################################
# export $PATH / $PYTHONPATH / $LD_LIBRARY_PATH   (dies sollte VOR! module load geschehen!)
#########################################################################
if [ "$currentshell" != "tcsh" ];then
    #echo 'export PATH="'"$PATH"'"; export PYTHONPATH="'"$PYTHONPATH"'"; export LD_LIBRARY_PATH="'"$LD_LIBRARY_PATH"'"'
    echo 'export PATH="'"$PATH"'"; export PYTHONPATH="'"$PYTHONPATH"'"; export LD_LIBRARY_PATH="'"$LD_LIBRARY_PATH"'"; export C_INCLUDE_PATH="/Users/glensk/.miniconda2/include"'
else
    tcshpath=`echo $PATH | sed 's|:| |g'`
    tcshpythonpath=`echo $PYTHONPATH | sed 's|:| |g'`
    echo 'set path = ( '"$tcshpath"' ); set PYTHONPATH = ( '"$tcshpythonpath"' );' # export C_INCLUDE_PATH="/Users/glensk/.miniconda2/include";'
fi

#export C_INCLUDE_PATH="/Users/glensk/.miniconda2/include"

#########################################################################
# unused / oldstuff 
#########################################################################
######### pycharm
#setenv JAVA_HOME /usr/local/java/jdk1.6.0_26
#alias pycharm '$HOME/Dropbox/scripts/unix_tools/pycharm.sh'  # this is just an alias

########## vmd (to show atomic structures)
#setenv VMDINSTALLNAME vmd
#setenv VMDINSTALLBINDIR /home/glensk/scripts/unix_tools/vmd/bin
#setenv VMDINSTALLLIBRARYDIR /home/glensk/scripts/unix_tools/vmd/lib

########## for phon dario alfe
#alias phon '$HOME/Dropbox/scripts/phonons/phon_dario_alfe/install_folder_phon/src/phon'

########## for xmgrace
#set path = ($HOME/Dropbox/scripts/mac_tools/apps/xmgrace/xmgrace_bin $path)

########## for ase
#setenv ASE_TAGS https://svn.fysik.dtu.dk/projects/ase/tags/
#setenv PYTHONPATH ${HOME}/ase-3.7.1:${PYTHONPATH}
#setenv PATH ${HOME}/ase/tools:${PATH}

########## for gracelot in ipython
#setenv PYTHONPATH "${PYTHONPATH}:$HOME/scripts/python/graceplot"

########## for phonopy (works everywhere)
#setenv phonpath $HOME/scripts/phonons/atushi_togo_phonopy/install_folder_phonopy/phonopy-1.7.2


#### for trunk olle Hellman
#setenv olletrunk $HOME/scripts/codes_to_try/olle_trunk
#setenv olletrunkshared $HOME/scripts/codes_to_try/olle_trunk_Shared
#setenv DYLD_LIBRARY_PATH /opt/intel/composer_xe_2013.3.171/mkl/lib:/opt/intel/composer_xe_2013.3.171/compiler/lib/
#set path = ($olletrunk/bin $path) # only workin at cmmd002
##set path = ($olletrunk/extract_forceconstants $path)

########## for VASP
#source /opt/intel/composer_xe_2015.1.108/mkl/bin/mklvars.sh intel64  # sets DYLD library path
# The next line is for vasp but seems now unnecessary
#export LD_LIBRARY_PATH=/opt/intel/composer_xe_2015.1.108/compiler/lib/intel64:/opt/intel/composer_xe_2015.1.108/mkl/lib/
#export LD_LIBRARY_PATH=/opt/intel/composer_xe_2015.1.108/mkl/lib/
#export MKLROOT=/opt/intel/composer_xe_2015.1.108/mkl
#export LD_LIBRARY_PATH="$HOME/local/lib:$LD_LIBRARY_PATH"
#export PATH="$HOME/local/bin:$PATH"
