#!/bin/bash
#################################################################
# this script returns host dependent environment variables 
# it does know all prior defined environment variables           
# it will not know the here defined environment variables 
# all here defined environment variables are defined by setenv  
# which will (in the end) work for bash / tcsh / zsh                         
#################################################################

#################################################################
# determine the current host;
#   on{mac,cmmc,cmmc,cmpc} == true/false
#   myshell == zsh / tcsh / bash
#   whichalias == alias / mkalias  (alias for tcsh) (mkalias for bash/zsh)
#   currenthost == on{mac,cmmc,cmmc,cmpc}
#################################################################

host=`hostname`   # 0.001s
echo "setenv host $host;
      setenv myhost cmpc34;
      setenv myhostmpie cmpc34.mpie.de;
      setenv mylaptop mac;
      setenv submithost cmmd002.mpie.de;
      setenv printloadstat false;"
# cmpc25=wolfgang; cmpc07=biswanath;

# Prompt colors: use 
# time:             zsh: magenta;    bash: orange;      tcsh: green
# host& path:       mac: magenta;    cmpc: red;       cmmd: blue;    cmmc:turquoise
if [ "$host" = "mac" ];then
    onxxx="setenv onmac true"
    currenthost="setenv currenthost onmac"
    myshell="setenv myshell zsh"
    myprompttime="setenv myprompttime magenta"  # yellow green blue
    myprompthostuser="setenv myprompthostuser magenta"
    mypromptpath="setenv mypromptpath magenta"
    module="setenv GZIP -9"   # an empty string ("") is not possible, better use dummy
    add="setenv ka kb"      # an empty string ("") is not possible, better use dummy
else 
    if [ "`echo $host | grep -o cmpc`" = "cmpc" ];then
    onxxx="setenv oncmpc true"
    currenthost="setenv currenthost oncmpc"
    myshell="setenv myshell zsh"
    myprompttime="setenv myprompttime magenta"
    myprompthostuser="setenv myprompthostuser red"
    mypromptpath="setenv mypromptpath red"
    module="module load mathematica1002/10.0.2"
    # christoph changed in /etc/profile.d/gxhivemgr.sh: if test x"$PS1" = x || echo
    add="limit stacksize unlimited;setenv TMPDIR /scratch/$USER;setenv PRINTER cmcopy1166"
else
    if [ "`echo $host | grep -o cmmc`" = "cmmc" ];then
    onxxx="setenv oncmmc true"
    currenthost="setenv currenthost oncmmc"
    myshell="setenv myshell zsh"
    myprompttime="setenv myprompttime cyan"
    myprompthostuser="setenv myprompthostuser cyan"
    mypromptpath="setenv mypromptpath cyan"
    module="module load git intel impi mkl grace mathematica/9.0.1 vasp/5.3.5" 
    add="limit stacksize unlimited;setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1;setenv submithost cmmc002"
else
    if [ "`echo $host | grep -o cmmd`" = "cmmd" ];then
    onxxx="setenv oncmmd true"
    currenthost="setenv currenthost oncmmd"
    myshell="setenv myshell zsh"
    myprompttime="setenv myprompttime blue"
    myprompthostuser="setenv myprompthostuser blue"
    mypromptpath="setenv mypromptpath blue"
    module="module load sge vasp/parallel/4.6-tdi sphinx/serial/2.0.4"
    add="setenv TMPDIR /scratch/$USER;limit stacksize unlimited;setenv GRACE_HOME /data/grabowski/xmgrace/grace-5.1.22/"
fi;fi;fi;fi

###################################################
# output defined variables
# currentshell:  (echo "echo JOJOJOJO $currentshell")
#   zsh : zsh setenv onmac true
#   tcsh: tcsh setenv onmac true
#   bash: bash setenv onmac true ; setenv myshell zsh ; ; ;
###################################################

whichalias="setenv whichalias mkalias"
[ "$currentshell" = "tcsh" ] && whichalias="setenv whichalias alias"


#[ "$currentshell" = "zsh" ] && echo  "echo  ...SHELL  zsh ALIAS $whichalias PRINT $printloadstat HOST $host;"
#[ "$currentshell" = "tcsh" ] && echo "echo  ...SHELL tcsh ALIAS $whichalias PRINT $printloadstat HOST $host;"
#[ "$currentshell" = "bash" ] && echo "echo  ...SHELL bash ALIAS $whichalias PRINT $printloadstat HOST $host;"

echo "
    setenv onmac false; 
    setenv oncmmc false; 
    setenv oncmmd false; 
    setenv oncmpc false;
    $myshell; $onxxx; $myprompttime; $myprompthostuser; $mypromptpath; $module; 
    $add; $currenthost; $whichalias"

