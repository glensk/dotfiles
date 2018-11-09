#!/bin/sh
###################################################################################
# this script defines environment variables independent of the used shell
# (tcsh;bash;zsh). Although it uses "setenv" which is usually used in the tcsh 
# shell, all set variables will be loaded by either shell (bash or tcsh or zsh) 
# since the setenv function is defined for the bash shell in the first loaded 
# $generalrc/generalrc_necessary_bash.sh script.
###################################################################################

###################################################################################
# this has to be the first line  (loads for zsh and bash); defines: tab-color; module; mkalias
###################################################################################
[ "$currentshell" != "tcsh" ] && source $generalrc/generalrc_necessary_bash.sh

host=`hostname`   # 0.001s
onhost=`echo $host | sed 's/\..*$//' | sed -r 's/[0-9]{1,10}$//'` # mac, cmpc, cosmopc, cmmc, daint, fidis
setenv host $host;
setenv onhost $onhost;

###################################################################################
# HOST dependent variables (module load, anaconda (on cmmc with module load) ...)
###################################################################################
eval `$generalrc/generalrc_hostdependent.sh`   
            # myshell to {zsh,bash, tcsh}; promptcolor; whichalias = {alias or mkalias}

###################################################################################
# PATH, PYTHONPATH  (anaconda on mac)
# PYTHONPATH should not be set (according to anaconda/miniconda)
###################################################################################
eval `$generalrc/generalrc_path.sh $onhost` 

#export C_INCLUDE_PATH="/Users/glensk/.miniconda2/include"
#@# time without aliases          : 0.38
#@# time with generalrc_alias_.sh : 0.62
#@# time with    alias            : 0.38 (in zsh: alias > alias; add as 1st row alias)

##############################################
# ALIASES & PROMPT & tabcolor
##############################################
[ ! -e "$generalrc/generalrc_alias_$onhost" ] && $generalrc/generalrc_alias_renew.sh
source $generalrc/generalrc_alias_$myshell.sh
source $generalrc/generalrc_prompt_$myshell.sh

case $currentshell in 
tcsh)
    source $generalrc/generalrc_alias_.sh $whichalias $currentshell
 ;;
*)  # for zsh and bash
    source $generalrc/generalrc_alias_$onhost
  ;;
esac

[ "$onhost" = "mac" ] && tab-color magenta
[ "$ohhost" = "daint" ] && tab-color blue

##############################################
# set host unspecific variables  
##############################################
setenv thermodynamics "$HOME/Thermodynamics"
setenv userme "glensk"
setenv i_pi_mc "$HOME/Dropbox/Albert/scripts/i-pi-mc/bin/i-pi"
setenv scripts "$HOME/Dropbox/Albert/scripts/dotfiles/scripts/"
setenv convcrit 0.5

# completion goo for commands
export PATH="$dotfiles/commands/:$PATH"  # if this is loaded before: complete -W ... $tocomplete the commands are executed which is not intended
tocomplete=`ls -1d $dotfiles/commands/* | sed 's|.*commands/||g'`
complete -W "n/*/`echo $tocomplete `/" goo

[ "$host" = "fidis" ] && setenv lmp "$scripts/lammps_executables/lmp_fidis_fidis_runner_2018_10_31"
[ "$host" = "mac" ]   && setenv lmp "$scripts/lammps_executables/lmp_serial_mac_runner_2018_10_30"

if [ -e "$thermodynamics/utilities/" ];then
[ "$currentshell" = "tcsh" ]  && source $thermodynamics/utilities/tcshrc_add
[ "$currentshell" != "tcsh" ] && source $thermodynamics/utilities/bashrc_add
fi

[ "$currentshell" != "bash" ] && limit coredumpsize 0    # Disable core dumps
setenv EDITOR vim   
setenv SVN_EDITOR vim        # for svn (Thermodynamics folder)
setenv LESS "-R"
setenv LC_ALL C   # necessary for perl git svn
setenv LANG C     # necessary for perl
    
setenv GIT_EDITOR vim
setenv PAGER most
setenv BROWSER chrome
setenv BROWSER open # is necessary for opening jupyter notebook files
setenv LC_ALL en_US.UTF-8
setenv GREP_COLOR 31    # red; some greps have colorized ouput. enable...
setenv GREPCOLOR 31     # dito here GREP_COLOR=1;32  # green

##############################################
# set host variables  
##############################################
[ "$printloadstat" = "true" ] && \
    echo " onhost $onhost"
