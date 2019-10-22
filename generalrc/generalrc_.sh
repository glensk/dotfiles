#!/bin/sh
##################################################################################
# this script defines environment variables independent of the used shell
# (tcsh;bash;zsh). Although it uses "setenv" (not any more) which is usually used in the tcsh 
# shell, all set variables will be loaded by either shell (bash or tcsh or zsh) 
# since the setenv function is defined for the bash shell in the first loaded 
# $generalrc/generalrc_necessary_bash.sh script.
##################################################################################

##################################################################################
# this has to be the first line  (loads senenv for zsh/bash); defines: tab-color; mkalias
##################################################################################

##################################################################################
# set global variables: currentshell, host, scripts, dotfiles
##################################################################################
#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (1) : $gett"
myprompttime="black"
[ "$BASH_VERSION" != "" ] && currentshell="bash" && myprompttime="red"
[ "$ZSH_VERSION" != "" ]  && currentshell="zsh"  && myprompttime="magenta"
export generalrc="$HOME/Dropbox/Albert/scripts/dotfiles/generalrc"
source $generalrc/generalrc_necessary_bash.sh  # loads setenv for bash/zsh (not any more)

echo "myhost generalrc_.sh in1:$myhost:"
#######################
# this is crucial, even if it is loaded a second time!
#######################
myhost=`myhost_get`
export myhost=$myhost
echo "myhost generalrc_.sh in2:$myhost:"
#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (0) : $gett before s0"

host=`hostname`   # 0.001s
export host=$host
#onhost=`echo $host | sed 's/\..*$//' | sed -r 's/[0-9]{1,10}$//'` # mac, cmpc, cosmopc, cmmc, daint, fidis
# currenly hostname:h258 myhost: helvetios onhost:UNKNOWN
#echo "myhost:$myhost"
##echo "onhost:$onhost"
#[ "$myhost" != "UNKNOWN" ] && [ "$onhost" = "UNKNOWN" ] && onhost=$myhost
#export onhost=$onhost
#echo "22myhost:$myhost"
#echo "22onhost:$onhost"

export dotfiles="$HOME/Dropbox/Albert/scripts/dotfiles/";
export potentials="$HOME/Dropbox/Albert/scripts/dotfiles/scripts/potentials";
export MYVIMRC="$dotfiles/nvim/init.vim"
export MYVIM="$HOME/sources/nvim/bin/nvim"
#jecho "copy1:$copy1:"
#echo "copy2:$copy2:"
[ ! -e "$HOME/.local/binp" ] && $dotfiles/bins/LINK_files.sh
[ ! -e "$HOME/.local/bins" ] && $dotfiles/bins/LINK_files.sh
[ ! -d "$HOME/sources" ] && mkdir $HOME/sources 


#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (1) : $gett time setenv diverses (change setenv to export)"

##################################################################################
# COSMOSTUFF: PATH, PYTHONPATH, LD_LIBRARY_PATH, ESPRESSO_PSEUDO, IPI_COMMAND, LAMMPS_COMMAND, scripts,
##################################################################################
source $dotfiles/scripts/source_to_add_to_path.sh
#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (1) : $gett time source_to_add_to_path.sh"

##################################################################################
# HOST dependent variables (myshell{=zsh,bash,tcsh}, module load, promptcolor, whichalias ...);  PATH due to module load
##################################################################################
source $generalrc/generalrc_hostdependent.sh
#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (3) : $gett time generalrc_hostdependent"

##################################################################################
# PATH, PYTHONPATH, LD_LIBRARY_PATH, C_INCLUDE_PATH (PYTHONPATH should not be set)
##################################################################################
source $generalrc/generalrc_path.sh $myhost

##############################################
# ALIASES & PROMPT & tabcolor
##############################################
source $generalrc/aliases.sh   # shellscript containing aliases
source $generalrc/generalrc_prompt_$currentshell.sh
    
#limit coredumpsize 0    # Disable core dumps # limit command is not know in bash

# this would be for tcsh which will not be enabled since I curretnly dont use tcsh
#source $generalrc/generalrc_alias_.sh $whichalias $currentshell

tab-color $mypromptpath

##############################################
# conda anaconda virtualenv (takes most of the time when loading)
##############################################
#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (4) : $gett before conda/aiida activate"

# on mac currently base, aiida, intelpy, python2 (12GB) (anaconda 2GB)
# the conda activate step takes all the time (not the source)
##############################################
# lets see weather we can survive without conda
# on mac: install stuff with easy_install / pip / brew
# on fidis: might be necessary to install click somehow but for now lets try without;
# fidis: pip install --install-option="--prefix=$HOME/.local" click
##############################################
# s3 goo zsh_set
# s1 CONDA autojump AFTER_CONDA
# s1 thermodynamics autojump generalrc
# conda is necessary for python lammps, for jupyter (mac) ... better load it.
case $myhost in
#    #cosmopc) source $HOME/aiida/bin/activate; ;;
#    # mac)       source $HOME/miniconda2/etc/profile.d/conda.sh && conda activate; ;;
#    # on mac: 
#    # pip install --upgrade --user phonopy
#    # pip install --upgrade --user ase
#    # pip install --upgrade --user lmfit
#    # pip install --upgrade --user intel-numpy    # to make numpy faster
#    # pip install --upgrade --user tqdm
#    # pip install --upgrade --user jupyter  # necessary to open ipynb notebooks
#    # pip install --upgrade --user jupyter_contrib_nbextensions   # notebooks table of contents
#    # pip install --upgrade --user sklearn  # necessary for andreas soap stuff
#    # pip install --upgrade --user plotly_express # to make nice scatterplots
#
#    fidis)     source $HOME/miniconda2/etc/profile.d/conda.sh && conda activate; ;;
#    helvetios) source $HOME/miniconda2/etc/profile.d/conda.sh && conda activate; ;;
#    # on fidis/helvetios:
#    # conda deactivate
#    # module load intel python/2.7.16
#    # pip install --upgrade --user 'ase<3.18.0'  # since since version 3.18.0 ase supports only python 3
#    # pip install --upgrade --user ipython
    fidis) moduel load intel python/2.7.16;;
    helvetios) moduel load intel python/2.7.16;;
esac
#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (5) : $gett CONDA"

##############################################
# set Thermodynamics stuff
# I currently dont use the Thermodynamics scritps
##############################################
export thermodynamics="$HOME/Thermodynamics"
export userme="glensk"
#export convcrit=0.5
if [ -e "$thermodynamics/utilities/" ];then
source $thermodynamics/utilities/bashrc_add
fi
#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (6) : $gett thermodynamics"

##############################################
# completion for {notes} commands  -> this needs to go to zshrc! (autoload/compdef)
##############################################
fpath=($dotfiles/completions_fpath $fpath)

##############################################
# other completitions 
##############################################
#if [ "$currentshell" = "zsh" ];then
##function goo() { $dotfiles/aliases/goo "$1" 
##}  # this nees to be a in second line
##autoload _goo    # do not forget BEFORE the next cmd! 
##compdef _goo goo # binds the completion function to a command
#fi

#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (6) : $gett goo"


##############################################
# general variables
##############################################
#export PAGER=most   # dont! makes problems with %git branch fatal: cannot run most: No such file or directory
export EDITOR=$MYVIM   
export SVN_EDITOR=$MYVIM        # for svn (Thermodynamics folder)
export GIT_EDITOR=$MYVIM 
export LESS="-R"
#export LC_ALL=C   # necessary for perl git svn
#export LANG C     # necessary for perl
export LANG=en_US.UTF-8     # necessary for powerline in zsh
export LC_ALL=en_US.UTF-8
export PAGER=less 
export BROWSER=chrome
export BROWSER=open # is necessary for opening jupyter notebook files
export GREP_COLOR=31    # red; some greps have colorized ouput. enable...
export GREPCOLOR=31     # dito here GREP_COLOR=1;32  # green
export NOTES_DIRECTORY=$dotfiles/notes  # /Users/glensk/Dropbox/Albert/scripts/dotfiles/aliases/notes
export NOTES_EXT="txt" # /Users/glensk/Dropbox/Albert/scripts/dotfiles/aliases/notes
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (6) : $gett setenv"

##############################################
# aiida-alloy (daniel marchand)
##############################################

##############################################
# autojump 
# takes a bit too long to laod ... and i dont use it
##############################################
if [ ! -e "$HOME/.local/autojump" ];then
    mkdir -p $HOME/.local/
    cd $HOME/.local/
    echo "installing autojump"
    git clone git://github.com/wting/autojump.git
    cd autojump 
    ./install.py #or ./uninstall.py
fi
#case "$currentshell" in 
#    Alternative: /Users/glensk/.local/zsh
#    zsh) source ~/.autojump/share/autojump/autojump.zsh;;
#    bash) [[ -s $HOME/.autojump/etc/profile.d/autojump.sh ]] && source $HOME/.autojump/etc/profile.d/autojump.sh;;
#esac
#[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (7) : $gett autojump"

##############################################
# shell dependent settings; 
# defines colors for ls; 
# bindkeys for history-search-bakcward ...
# bash/zsh completion
##############################################
source $dotfiles/$currentshell/$currentshell\_set
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (7) : $gett zsh_set"

##############################################
[ "$printloadstat" = "true" ] && echo " myhost $myhost"

[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (6) : $gett generalrc_AFTER_CONDA"

########################## make sure that nvim is istalled
#if [ "`command -v curl`" = "" ] && pip install --upgrade --user curl
[ ! -e "$MYVIM" ] && echo "MYVIM is not installed!" && install_git.py -i nvim

########################## make sure that nvim has xclip/pbcopy
copy1=`command -v pbcopy`
copy2=`command -v xclip`
if [ "$copy1" = "" ] && [ "$copy2" = "" ]; then 
    echo "######################################################"
    echo "no local clipboard pbcopy (mac) / xclip (linux) / xsel" 
    echo "######################################################"
    xclip_local_install.sh
fi

