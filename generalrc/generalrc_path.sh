#!/bin/sh
#########################################################################
# define system dependent (and system independent) path variables
# e.g. $PATH / $PYTHONPATH / $LD_LIBRARY_PATH
# all prior defined environment variables are known
#########################################################################


#########################################################################
# on every system: ( $PATH / $PYTHONPATH / $LD_LIBRARY_PATH / ...)
#########################################################################
dropboxbin="$HOME/Dropbox/scripts/dotfiles/bin"       # Dropbox binaries
anaconda="$HOME/anaconda/bin"
svnctags="/usr/local/bin"
homebrew="/usr/local/sbin"
bin="/usr/bin"
PATH="$dropboxbin:$anaconda:$svnctags:$homebrew:$bin:$PATH"


#pythermodynamics=
#PYTHONPATH="$PYTHONPATH:"
#########################################################################
# on mac: ( $PATH / $PYTHONPATH / $LD_LIBRARY_PATH / ... )
#########################################################################
if [ "$currenthost" = "onmac" ];then
    # $PATH (on mac)
    sed="/usr/local/opt/gnu-sed/libexec/gnubin"
    # !! too much time !!! brew="$(brew --prefix coreutils)/libexec/gnubin" 
    brew="/usr/local/opt/coreutils/libexec/gnubin" 
    tdep="$HOME/Dropbox/scripts/phonons/tdep-devel/bin"
    phonopy="$HOME/scripts/phonons/phonopy_at_$host"
    phonopybin="$phonopy/bin"
    vasp="$HOME/local/bin"
    #vimctags="$HOME/Dropbox/scripts/dotfiles/vim/ctags-5.8/installfolder/bin"
    #sphinx="$HOME/Dropbox/scripts/phonons/sphinx/bin"


    # $PYTHONPATH (on mac)
    pygraceplot="$HOME/scripts/python/graceplot"
    pylammps1="/usr/local/Cellar/lammps/2014.02.12/lib/python2.7/site-packages"
    pylammps2="/usr/local/lib/python2.7/site-packages"
    pyphonopy1="$phonopy/lib/python"     # may be lib64 instead
    pyphonopy2="$phonopy/lib64/python"   # may be lib64 instead


    # $LD_LIBRARY_PATH add on mac
    ldvasp="$HOME/local/lib"
    #setenv LD_LIBRARY_PATH '/afs/@cell/common/soft/intel/mkl/lib/intel64/'
    
    #[ "$host" = "$mylaptop" ] && \
    PATH="$sed:$brew:$tdep:$phonopy:$phonopybin:$vasp:$PATH" #&& \
    PYTHONPATH="$pygraceplot:$pylammps1:$pylammps2:$pyphonopy1:$PYTHONPATH" #&& \
    LD_LIBRARY_PATH="$ldvasp:$LD_LIBRARY_PATH"
fi


#########################################################################
# on cmmd: ( $PATH / $PYTHONPATH / $LD_LIBRARY_PATH / ... )
#########################################################################
if [ "$currenthost" = "oncmmd" ];then
    
    xmgracebin="/data/grabowski/xmgrace/bin/"
    add1="/lib64"
    PATH="$xmgracebin:$add1:$PATH" 

    argparse="${HOME}/Dropbox/scripts/python/add_cmmd/argparse-1.2.1"
    PYTHONPATH="$argparse:$PYTHONPATH" 

    # ldgibbsform ist fuer getGibbsEnergyOfFormation.sh
    ldgibbsform="/opt/intel/Compiler/11.1/073/mkl/lib/em64t/"   
    ld1="/home/grabowski/libs/"
    ld2="/home/glensk/"
    ldxmgrace="/opt/intel/lib/intel64"
    LD_LIBRARY_PATH="$ldgibbsform:$ld1:$ld2:$ldxmgrace:$LD_LIBRARY_PATH"
fi
    
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
# export $PATH / $PYTHONPATH / $LD_LIBRARY_PATH
#########################################################################

if [ "$currentshell" != "tcsh" ];then
    echo 'export PATH="'"$PATH"'"; export PYTHONPATH="'"$PYTHONPATH"'"; export LD_LIBRARY_PATH="'"$LD_LIBRARY_PATH"'"'
else
    tcshpath=`echo $PATH | sed 's|:| |g'`
    tcshpythonpath=`echo $PYTHONPATH | sed 's|:| |g'`
    echo 'set path = ( '"$tcshpath"' ); set PYTHONPATH = ( '"$tcshpythonpath"' );'
fi


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
