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

###################################################################################
# HOST dependent variables (prompt, module load, anaconda (on cmmc with module load) ...)
###################################################################################
eval `$generalrc/generalrc_hostdependent.sh`   
            # defines:  
            # --------
            # onmac, oncmmc, oncmdft} to true or false
            # currenthost to {onmac, oncmmc, oncmdft}
            # myshell to {zsh,bash, tcsh}
            # prompt, promptcolor
            # printloadstats
            # echo ONMAC $onmac
            # echo WHICHALIAS $whichalias  
            # echo CURRENTHOST $currenthost


###################################################################################
# PATH, PYTHONPATH  (anaconda on mac)
# PYTHONPATH should not be set (according to anaconda/miniconda)
###################################################################################
eval `$generalrc/generalrc_path.sh $currenthost` 

#export C_INCLUDE_PATH="/Users/glensk/.miniconda2/include"
#@# time without aliases          : 0.38
#@# time with generalrc_alias_.sh : 0.62
#@# time with    alias            : 0.38 (in zsh: alias > alias; add as 1st row alias)

##############################################
# ALIASES
##############################################
[ ! -e "$generalrc/generalrc_alias_$currenthost" ] && $generalrc/generalrc_alias_renew.sh

# aliases which worc everywhere
[ "$currentshell" != "tcsh" ] && source $generalrc/generalrc_alias_$currenthost
[ "$currentshell" = "tcsh" ] && source $generalrc/generalrc_alias_.sh $whichalias $currentshell

# shell SPECIFIC aliases
[ "$currentshell"  = "zsh" ] && source $generalrc/generalrc_alias_zsh.sh
[ "$currentshell"  = "bash" ] && source $generalrc/generalrc_alias_bash.sh
[ "$currentshell"  = "tcsh" ] && source $generalrc/generalrc_alias_tcsh.sh

##############################################
# PROMPT & tabcolor
##############################################
[ "$currentshell" = "tcsh" ] && source $generalrc/generalrc_prompt_tcsh.sh
[ "$currentshell" = "zsh" ] && source $generalrc/generalrc_prompt_zsh.sh
[ "$currentshell" = "bash" ] && source $generalrc/generalrc_prompt_bash.sh
[ "$onmac" = "true" ] && tab-color magenta
[ "$ondaint" = "true" ] && tab-color blue

#echo
#echo MYSHELL $myshell;
#echo CURRENTSHELL $currentshell;
#echo MYPROMPTTIME $myprompttime;
#eval `$generalrc/generalrc_prompt.sh`
#[ "$printloadstat" = "true" ] && echo "... generalrc_prompt.sh done"
#echo MYPROMPTTIME $myprompttime;
#echo 

#[ "`hostname`" = "cmmc002" ] && load_local_anaconda
#[ "`hostname`" = "cmmc001" ] && load_local_anaconda
##############################################
# set host unspecific variables  
##############################################
setenv thermodynamics "$HOME/Thermodynamics"
setenv userme "glensk"
setenv i_pi_mc "$HOME/Dropbox/Albert/scripts/i-pi-mc/bin/i-pi"
setenv scripts "$HOME/Dropbox/Albert/scripts/dotfiles/scripts/"
setenv convcrit 0.5
[ "$host" = "fidis" ] && setenv lmp "$scripts/lammps_executables/lmp_fidis_fidis_runner_2018_10_31"
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
    echo " +onmac $onmac"   && echo " +oncmmc $oncmmc" && \
    echo " +oncmmd $oncmmd" && echo " +oncmpc $oncmpc"

