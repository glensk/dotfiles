#!/bin/sh
#####################################################
# this script shall return host dependent settings  #
#####################################################
host=`hostname`
mylaptop=mac
if [ "$host" = "$mylaptop" ];then
    onxxx="setenv onmac true"
    myshell="setenv myshell zsh"
    module=""
else 
    if [ "`echo $host | grep -o cmpc`" = "cmpc" ];then
    onxxx="setenv oncmpc true"
    myshell="setenv myshell zsh"
    module="module load mathematica1002/10.0.2"
    # christoph changed in /etc/profile.d/gxhivemgr.sh: if test x"$PS1" = x || echo
    add="limit stacksize unlimited;setenv TMPDIR /scratch/$USER"
else
    if [ "`echo $host | grep -o cmmc`" = "cmmc" ];then
    onxxx="setenv oncmmc true"
    myshell="setenv myshell tcsh"
    module="module load git intel impi mkl grace mathematica/9.0.1 vasp/5.3.5" 
    add="limit stacksize unlimited;setenv vaspq /u/aglen/vasp/vasp_4.6_lj_morse_alles_copied_FFTWPLANS_works/vasp.v1"
else
    if [ "`echo $host | grep -o cmmd`" = "cmmd" ];then
    onxxx="setenv oncmmd true"
    myshell="setenv myshell tcsh"
    module="module load sge vasp/parallel/4.6-tdi sphinx/serial/2.0.4"
    add="setenv TMPDIR /scratch/$USER;limit stacksize unlimited;setenv GRACE_HOME /data/grabowski/xmgrace/grace-5.1.22/"
fi;fi;fi;fi

whichalias="setenv whichalias mkalias"; 
#[ "$currentshell" = "tcsh" ] && whichalias="setenv whichalias alias"
#[ "$currentshell" = "bash" ] && whichalias="setenv whichalias mkalias"

if [ "$currentshell" = "bash" ];then
echo "$onxxx \; $myshell \; $module \; $add \; $whichalias"
else
echo $onxxx \; $myshell \; $module \; $add \; $whichalias
fi
