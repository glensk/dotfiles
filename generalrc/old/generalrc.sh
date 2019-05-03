#!/bin/sh

[ "$currentshell" != "tcsh" ] && source $generlrc/generalrc_necessary_bash.sh
###################################################################################
# this script defines environment variables independent of the used shell
# (tcsh;bash;zsh). Although it uses "setenv" --- which is usually used in the tcsh 
# shell --- all set variables will be loaded by either shell (bash or tcsh or zsh) 
# if the setenv function is defined in the bash shell (previously to sourcing this script).
###################################################################################
#[ "$currentshell" -ne "tcsh" ] && exe() { echo "\$ $@" ; "$@" ; }
setenv host `hostname`
setenv myhost cmpc34               # cmpc25=wolfgang; cmpc07=biswanath;
setenv myhostmpie $myhost.mpie.de   # cmpc25=wolfgang; cmpc07=biswanath;
setenv mylaptop mac
setenv submithost cmmd002.mpie.de
setenv printloadstat true

##############################################
# determine the current host; from generalrc_hostdependent.sh:
#   on{mac,cmmc,cmmc,cmpc} == true/false
#   myshell == zsh / tcsh / bash
#   whichalias == alias / mkalias
##############################################
setenv onmac false;setenv oncmmc false; setenv oncmmd false; setenv oncmpc false
echo 22
eval `$generalrc/generalrc_hostdependent.sh`
eval `$generalrc/generalrc_path.sh` 
source $generalrc/generalrc_alias.sh $whichalias
#eval `$HOME/Dropbox/scripts/dotfiles/generalrc_prompt.sh` # --> onXXX true

##############################################
# set host unspecific variables  
##############################################
[ "$currentshell" != "bash" ] && limit coredumpsize 0    # Disable core dumps
setenv EDITOR vim   
setenv LC_ALL C   # necessary for perl git svn
setenv LANG C     # necessary for perl
    
setenv GIT_EDITOR vim
setenv PAGER most
setenv BROWSER chrome
setenv LC_ALL en_US.UTF-8
setenv GREP_COLOR 31    # red; some greps have colorized ouput. enable...
setenv GREPCOLOR 31     # dito here GREP_COLOR=1;32  # green

##############################################
# set host variables  
##############################################
#eval `$HOME/Dropbox/scripts/dotfiles/generalrc_module_load.sh`  # has to go seperately

[ "$printloadstat" = "true" ] && \
    echo " +onmac $onmac"   && echo " +oncmmc $oncmmc" && \
    echo " +oncmmd $oncmmd" && echo " +oncmpc $oncmpc"

##############################################
# set host specific variables  
# aliase dont work yet; 
##############################################

# for local setting see tcsh/tcsh.local!!!!!!!

#@ [ "$onmac" = "true" ] && \
#@     #alias trash="rmtrash" && \   
#@     #alias   del="rmtrash" && \
#@     #alias rm="echo Use 'del', or the full path i.e. '/bin/rm'" && \
#@     #alias ctagsnew='/usr/local/Cellar/ctags/5.8/bin/ctags -R .'
#alias wiki='web_search wiki'
#alias google='web_search google'

#aliasgeneral () {
#    alias $1 '"$2"'
#}
#

