#!/bin/sh

iterm2_get_window_number () {
    #echo $ITERM_SESSION_ID 
    echo "${ITERM_SESSION_ID:1:1}" | awk '{print $1+1}'
}

iterm2_set_window_name () {
    window_nr=`iterm2_get_window_number`
    export ITERM2_WINDOW$window_nr\_NAME="$*"
    rm -f $HOME/.iterm2_windowname_$window_nr && echo "$*" > $HOME/.iterm2_windowname_$window_nr
}

iterm2_get_window_name() {
    window_nr=`iterm2_get_window_number`
    name=`cat $HOME/.iterm2_windowname_$window_nr`
    echo $name
}

iterm2_write_window_name () {
    window_nr=`iterm2_get_window_number`
    [ "$1" != "" ] && iterm2_set_window_name "$*"
    name=`iterm2_get_window_name`
    [ "$name" != "" ] && export AWS_PROFILE="$name"
}

README_timestamp () {
    timestamp=`python -c 'from datetime import datetime;time_now = datetime.now();print(time_now.strftime("%Y-%m-%d_%H:%M:%S"))'`
    echo $timestamp
}


#setenv () {
#    if [ "x$1" = "x" ] ; then
#      echo "$0: environment variable name required" >&2
#    elif [ "x$2" = "x" ] ; then
#      echo "$0: environment variable value required" >&2
#    else
#      export $1=$2
#    fi
#}

#getms() {
#    # for mac: /usr/local/opt/coreutils/libexec/gnubin/date +%s%N
#    $date +%s%N
#}

# setenv wird nur von der tcsh benutzt
#setenv () {
#      export $1=$2
#}

function iterm_title {
        echo -ne "\033]0;"$*"\007"
}



#setenv () { export $1=$2 }

#aliastcshtobash () {
#	if [ "x$2" = "x" ]; then
#		echo alias ${1}="''"
#	elif echo "$2" | egrep -s '(\!|#)' >/dev/null 2>&1; then
#		comm=`echo $2 | sed  's/\\!\*/"$\@"/g
#				      s/\\!:\([1-9]\)/"$\1"/g
#			              s/#/\#/g'`
#		echo $1 \(\) "{" command "$comm"  "; }"
#	else
#		echo alias ${1}=\'`echo "${2}" | sed "s:':'\\\\\\\\'':"`\'
#	fi
#}

#aliastcshtobash () {
#	#if [ "x$2" = "x" ]; then
#	#	echo alias ${1}="''"
#	if echo "$2" | egrep -s '(\!|#)' >/dev/null 2>&1; then
#		comm=`echo $2 | sed  's/\\!\*/"$\@"/g
#				      s/\\!:\([1-9]\)/"$\1"/g
#			              s/#/\#/g'`
#		echo $1 \(\) "{" command "$comm"  "; }"
#	else
#		echo alias ${1}=\'`echo "${2}" | sed "s:':'\\\\\\\\'':"`\'
#	fi
#}

#aliastcshtobash () {
#	#if [ "x$2" = "x" ]; then
#	#	echo alias ${1}="''"
#	#if echo "$2" | egrep -s '(\!|#)' >/dev/null 2>&1; then
#	#	comm=`echo $2 | sed  's/\\!\*/"$\@"/g
#	#			      s/\\!:\([1-9]\)/"$\1"/g
#	#		              s/#/\#/g'`
#	#	echo $1 \(\) "{" command "$comm"  "; }"
#	#else
#		echo alias ${1}=\'`echo "${2}" | sed "s:':'\\\\\\\\'':"`\'
#	#fi
#}

#aliastcshtobash () {
#		echo alias ${1}=\'`echo "${2}" | sed "s:':'\\\\\\\\'':"`\'
#}

#mkalias () {
#    # for bash $2 need to be in ""
#    eval `aliastcshtobash $1 "$2"`
#}

# THis was ment only for garching
# if [ "$ZSH_VERSION" = "5.0.5" ];then
# module () {
#     #eval `/afs/ipp/common/usr/modules.2014/amd64_sles11/Modules/$MODULE_VERSION/bin/modulecmd sh $*`
#     # log in to garching by \ssh aglen@cmmc001.bc.rzg.mpg.de
#     # echo $SHELL -> /bin/bash
#     # echo $MODULEHOME  ## !! and add .../bin/modulecmd
#     eval `/afs/ipp-garching.mpg.de/common/usr/modules.2017/amd64_sles12/Modules/current/bin/modulecmd sh $*`
# 
# }
# fi

tab-color() {
    [ "$1" = "red" ]       && a=255 b=0   c=0           # works
    [ "$1" = "orange" ]    && a=255 b=165 c=0           # works
    [ "$1" = "yellow" ]    && a=255 b=165 c=0           # works
    [ "$1" = "black" ]     && a=0   b=0   c=0           
    [ "$1" = "green" ]     && a=0   b=255 c=0           # works
    [ "$1" = "blue" ]      && a=0   b=0   c=255         # works
    [ "$1" = "cyan" ]      && a=64  b=244 c=208         # works
    [ "$1" = "turquoise" ] && a=64  b=244 c=208         # works
    [ "$1" = "mediumturquoise" ]    && a=72 b=209 c=204  
    [ "$1" = "zomp" ]    && a=57 b=167 c=142  
    [ "$1" = "magenta" ] && a=255 b=0 c=255         # works but you need to write it correctly
    #[ "$1" = "turquoise" ] && a=100 b=100 c=208
echo -ne "\033]6;1;bg;red;brightness;$a\a"
echo -ne "\033]6;1;bg;green;brightness;$b\a"
echo -ne "\033]6;1;bg;blue;brightness;$c\a"
}

load_local_anaconda() {
    echo 'module unload anaconda'
    module unload anaconda
    echo
    echo 'export PATH="/u/aglen/conda-envs/my_root/bin:$PATH"'
    export PATH="/u/aglen/conda-envs/my_root/bin:$PATH"
    echo
    echo 'source activate /u/aglen/conda-envs/my_root'
    source activate /u/aglen/conda-envs/my_root
}

myhost_get() {
    # shouls recognise: mac, helvetios, fidis, cmmd, cosmopc, f209, g308, h211
    host=`hostname`
    case $host in
        "mac") echo "mac" && return;;
        "daint") echo "daint" && return;;
        "helvetios") echo "helvetios" && return;;
        "fidis") echo "fidis" && return;;
        "cmmd") echo "cmmd" && return;;
        "cosmopc") echo "cosmopc" && return;;
    esac
    first_five="${host:0:5}"
    #echo first_five: $first_five
    case $first_five in
        "daint") echo "daint" && return;;
        "fidis") echo "fidis" && return;;
    esac
    
    first_seven="${host:0:7}"
    case $first_seven in
        "cosmopc") echo "cosmopc" && return;;
    esac

    first_letter="${host:0:1}"
    remaining_letter="${host:1:999}"
    re='^[0-9]+$'
    #echo first_letter $first_letter
    #echo remaining_letter $remaining_letter
    if [[ $remaining_letter =~ $re ]];then # remaining is an integer
        case $first_letter in
            "h") echo "helvetios" && return;;
            "f") echo "fidis" && return;;
            "g") echo "fidis" && return;;
        esac
    fi
    echo UNKNOWN 
    
}

conda_activate() {
    echo myhost `myhost`
    if [ "`myhost`" = "fidis" ];then
        echo "module purge"
        module purge
    elif [ "`myhost`" = "daint" ];then
        echo 00
        echo "/store/marvel/mr23/aglensk/miniconda3"
        source /store/marvel/mr23/aglensk/miniconda3/etc/profile.d/conda.sh
        conda activate
    fi
    
    if [ -e "$HOME/miniconda2" ];then
        echo 11
        echo "$HOME/miniconda2"
        source $HOME/miniconda2/etc/profile.d/conda.sh
        conda activate
    elif [ -e "$HOME/miniconda3" ];then
        echo 22
        echo "$HOME/miniconda3"
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate
    elif [ -e "$SCRATCH/miniconda2" ];then
        echo 33
        echo "$SCRATCH/miniconda2"
        source $SCRATCH/miniconda2/etc/profile.d/conda.sh
        conda activate
    elif [ -e "$SCRATCH/miniconda3" ];then
        echo 44
        echo "$SCRATCH/miniconda3"
        echo 55
        echo $SCRATCH
        echo 66
        source $SCRATCH/miniconda3/etc/profile.d/conda.sh
        echo 77
        conda activate
    fi
}

ca() {
    conda_activate
}

virtualenv_activate () {
     ACTIVATE_SCRIPT=$HOME/virtualenvs/$1/bin/activate;
     if [[ ! -f "$ACTIVATE_SCRIPT" ]]; then
         echo $1 is not a valid virtual env environment;
         echo Try one of: $(ls $HOME/virtualenvs/);
         return;
     fi;
     source $ACTIVATE_SCRIPT
}

start_notebook_fidis () {
    ipynbs=`find . -maxdepth 1 -name "*.ipynb"`
    echo ipynbs ":$ipynbs:"
    if [ "$ipynbs" != "" ];then
        echo
        echo 'starting: conda activate gap'
        echo 'starting: jupyter notebook --no-browser --port=8886'
        echo 'to open the notebook on mac: open_notebok_from_fidis'
        echo
        echo
        conda activate gap
        jupyter notebook --no-browser --port=8886
    else
        echo "no ipynb notebooks found"
    fi
}

open_notebook_from_fidis () {
    lsof -ti:8889 | xargs kill -9
    ssh -N -f -L localhost:8889:localhost:8886 glensk@fidis.epfl.ch
    open http://localhost:8889
}

# was before in aliases.sh function
mcd () { 
    mkdir -p $1 
    cd $1 
}

daint_modules () {
    module load daint-mc && module switch PrgEnv-cray PrgEnv-intel && module unload cray-libsci && module load GSL/2.5-CrayIntel-18.08 cray-python/2.7.15.1 cray-fftw
    module list
}

mld () {
    daint_modules
}
