#!/bin/sh
###################################################################################
# this script defines environment variables independent of the used shell
# (tcsh;bash;zsh). Although it uses "setenv" which is usually used in the tcsh 
# shell, all set variables will be loaded by either shell (bash or tcsh or zsh) 
# since the setenv function is defined for the bash shell in the first loaded 
# $generalrc/generalrc_necessary_bash.sh script.
###################################################################################

###################################################################################
# this has to be the first line  (loads senenv for zsh/bash); defines: tab-color; module; mkalias
###################################################################################

###################################################################################
# set global variables: currentshell, host, onhost, scripts, dotfiles
###################################################################################
case "$0" in bash) export currentshell="bash";;tcsh) setenv currentshell "tcsh";;*) export currentshell="zsh";esac
export generalrc="$HOME/Dropbox/Albert/scripts/dotfiles/generalrc"
source $generalrc/generalrc_necessary_bash.sh  # loads setenv for bash/zsh
host=`hostname`   # 0.001s
onhost=`echo $host | sed 's/\..*$//' | sed -r 's/[0-9]{1,10}$//'` # mac, cmpc, cosmopc, cmmc, daint, fidis
setenv host $host;
setenv onhost $onhost;
setenv dotfiles "$HOME/Dropbox/Albert/scripts/dotfiles/";
setenv scripts  "$HOME/Dropbox/Albert/scripts/dotfiles/scripts/";
source $scripts/source_to_add_to_path.sh

###################################################################################
# HOST dependent variables (module load, (on cmmc with module load) ...)
###################################################################################
eval `$generalrc/generalrc_hostdependent.sh`   
            # myshell to {zsh,bash, tcsh}; promptcolor; whichalias = {alias or mkalias}

###################################################################################
# PATH, PYTHONPATH
# PYTHONPATH should not be set (according to anaconda/miniconda)
###################################################################################
eval `$generalrc/generalrc_path.sh $onhost` 

#eval "$(_kmc_createjob_py_COMPLETE=source kmc_createjob.py)"
#@# time without aliases          : 0.38
#@# time with generalrc_alias_.sh : 0.62
#@# time with    alias            : 0.38 (in zsh: alias > alias; add as 1st row alias)

##############################################
# ALIASES & PROMPT & tabcolor
##############################################
[ ! -e "$generalrc/generalrc_alias_$onhost" ] && $generalrc/generalrc_alias_renew.sh
source $generalrc/generalrc_alias_$currentshell.sh
source $generalrc/generalrc_prompt_$currentshell.sh

# this works only for bash and zsh
source $generalrc/generalrc_alias_$onhost
#limit coredumpsize 0    # Disable core dumps # limit command is not know in bash

# this would be for tcsh which will not be enabled since I curretnly dont use tcsh
#source $generalrc/generalrc_alias_.sh $whichalias $currentshell

tab-color $mypromptpath

##############################################
# conda anaconda virtualenv (takes most of the time when loading)
##############################################
[ "$onhost" = "mac" ] && source $HOME/miniconda3/etc/profile.d/conda.sh && conda activate intelpy 


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

##############################################
# completion goo for commands
##############################################
export PATH="$dotfiles/commands/:$PATH"  # if this is loaded before: complete -W ... $tocomplete the commands are executed which is not intended
#tocomplete=`ls -1d $dotfiles/commands/* | sed 's|.*commands/||g'`
#complete -W "n/*/`echo $tocomplete `/" goo


setenv EDITOR vim   
setenv SVN_EDITOR vim        # for svn (Thermodynamics folder)
setenv LESS "-R"
setenv LC_ALL C   # necessary for perl git svn
setenv LANG C     # necessary for perl
    
setenv GIT_EDITOR vim
#setenv PAGER most   # dont! makes problems with %git branch fatal: cannot run most: No such file or directory
setenv PAGER less 
setenv BROWSER chrome
setenv BROWSER open # is necessary for opening jupyter notebook files
setenv LC_ALL en_US.UTF-8
setenv GREP_COLOR 31    # red; some greps have colorized ouput. enable...
setenv GREPCOLOR 31     # dito here GREP_COLOR=1;32  # green


##############################################
# autojump 
##############################################
case "$currentshell" in 
    zsh) source ~/.autojump/share/autojump/autojump.zsh;;
    bash) [[ -s $HOME/.autojump/etc/profile.d/autojump.sh ]] && source $HOME/.autojump/etc/profile.d/autojump.sh;;
esac

##############################################
# shell dependent settings; defines colors for ls; bindkeys for history-search-bakcward ...
##############################################
source $dotfiles/$currentshell/$currentshell\_set

##############################################
# set host variables  
##############################################
[ "$printloadstat" = "true" ] && \
    echo " onhost $onhost"
