#!/bin/sh
##################################################################################
# this script defines environment variables independent of the used shell
# (tcsh;bash;zsh). Although it uses "setenv" which is usually used in the tcsh 
# shell, all set variables will be loaded by either shell (bash or tcsh or zsh) 
# since the setenv function is defined for the bash shell in the first loaded 
# $generalrc/generalrc_necessary_bash.sh script.
##################################################################################

##################################################################################
# this has to be the first line  (loads senenv for zsh/bash); defines: tab-color; module; mkalias
##################################################################################

##################################################################################
# set global variables: currentshell, host, onhost, scripts, dotfiles
##################################################################################
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (1) : $gett"
[ "$BASH_VERSION" != "" ] && currentshell="bash"
[ "$ZSH_VERSION" != "" ] && currentshell="zsh"
export generalrc="$HOME/Dropbox/Albert/scripts/dotfiles/generalrc"
source $generalrc/generalrc_necessary_bash.sh  # loads setenv for bash/zsh
host=`hostname`   # 0.001s
onhost=`echo $host | sed 's/\..*$//' | sed -r 's/[0-9]{1,10}$//'` # mac, cmpc, cosmopc, cmmc, daint, fidis
setenv host $host;
setenv onhost $onhost;
setenv dotfiles "$HOME/Dropbox/Albert/scripts/dotfiles/";

##################################################################################
# COSMOSTUFF: PATH, PYTHONPATH, LD_LIBRARY_PATH, ESPRESSO_PSEUDO, IPI_COMMAND, LAMMPS_COMMAND, scripts,
##################################################################################
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (0) : $gett before s0"
source $dotfiles/scripts/source_to_add_to_path.sh
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (1) : $gett time s1"

##################################################################################
# HOST dependent variables (myshell{=zsh,bash,tcsh}, module load, promptcolor, whichalias ...);  PATH due to module load
##################################################################################
source $generalrc/generalrc_hostdependent.sh
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (3) : $gett time s3"

##################################################################################
# PATH, PYTHONPATH, LD_LIBRARY_PATH, C_INCLUDE_PATH (PYTHONPATH should not be set)
##################################################################################
source $generalrc/generalrc_path.sh $onhost

##############################################
# ALIASES & PROMPT & tabcolor
##############################################
if [ "$onhost" != "UNKNOWN" ];then
    [ ! -e "$generalrc/generalrc_alias_$onhost" ] && $generalrc/generalrc_alias_renew.sh
    # this works only for bash and zsh
    source $generalrc/generalrc_alias_$onhost
fi
source $generalrc/generalrc_alias_$currentshell.sh
source $generalrc/generalrc_prompt_$currentshell.sh
    
#limit coredumpsize 0    # Disable core dumps # limit command is not know in bash

# this would be for tcsh which will not be enabled since I curretnly dont use tcsh
#source $generalrc/generalrc_alias_.sh $whichalias $currentshell

tab-color $mypromptpath

##############################################
# conda anaconda virtualenv (takes most of the time when loading)
##############################################
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (4) : $gett before conda/aiida activate"

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
case $onhost in
    mac) source $HOME/miniconda2/etc/profile.d/conda.sh && conda activate; ;;
    cosmopc) source $HOME/aiida/bin/activate; ;;
    fidis) source $HOME/miniconda2/etc/profile.d/conda.sh && conda activate; ;;
esac
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (5) : $gett CONDA"

##############################################
# set Thermodynamics stuff
##############################################
setenv thermodynamics "$HOME/Thermodynamics"
setenv userme "glensk"
setenv convcrit 0.5
if [ -e "$thermodynamics/utilities/" ];then
#[ "$currentshell" = "tcsh" ]  && source $thermodynamics/utilities/tcshrc_add
#[ "$currentshell" != "tcsh" ] && source $thermodynamics/utilities/bashrc_add
source $thermodynamics/utilities/bashrc_add
fi
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (6) : $gett thermodynamics"

##############################################
# completion goo for commands
##############################################
fpath=($dotfiles/completions_fpath $fpath)
function goo() { $dotfiles/aliases/goo $1 }
autoload _goo # do not forget BEFORE the next cmd! 
compdef _goo goo # binds the completion function to a command

[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (6) : $gett goo"


##############################################
# general variables
##############################################
#setenv PAGER most   # dont! makes problems with %git branch fatal: cannot run most: No such file or directory
setenv EDITOR vim   
setenv SVN_EDITOR vim        # for svn (Thermodynamics folder)
setenv LESS "-R"
setenv LC_ALL C   # necessary for perl git svn
setenv LANG C     # necessary for perl
setenv GIT_EDITOR vim
setenv PAGER less 
setenv BROWSER chrome
setenv BROWSER open # is necessary for opening jupyter notebook files
setenv LC_ALL en_US.UTF-8
setenv GREP_COLOR 31    # red; some greps have colorized ouput. enable...
setenv GREPCOLOR 31     # dito here GREP_COLOR=1;32  # green
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (6) : $gett setenv"

##############################################
# autojump 
##############################################
case "$currentshell" in 
    zsh) source ~/.autojump/share/autojump/autojump.zsh;;
    bash) [[ -s $HOME/.autojump/etc/profile.d/autojump.sh ]] && source $HOME/.autojump/etc/profile.d/autojump.sh;;
esac
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (7) : $gett autojump"

##############################################
# shell dependent settings; defines colors for ls; bindkeys for history-search-bakcward ...
##############################################
source $dotfiles/$currentshell/$currentshell\_set
[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (7) : $gett zsh_set"

##############################################
# set host variables  
##############################################
[ "$printloadstat" = "true" ] && echo " onhost $onhost"

[ "$gettime" = "true" ] && gett=`gt $gett` && echo "general (6) : $gett generalrc_AFTER_CONDA"
