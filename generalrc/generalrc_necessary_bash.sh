#!/bin/sh


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

setenv () {
      export $1=$2
}

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

aliastcshtobash () {
		echo alias ${1}=\'`echo "${2}" | sed "s:':'\\\\\\\\'':"`\'
}

mkalias () {
    # for bash $2 need to be in ""
    eval `aliastcshtobash $1 "$2"`
}

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
    [ "$1" = "red" ] && a=255 b=0 c=0               # works
    [ "$1" = "green" ] && a=0 b=255 c=0             # works
    [ "$1" = "blue" ] && a=0 b=0 c=255              # works
    [ "$1" = "cyan" ] && a=64 b=244 c=208           # works
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
        echo 'open the notebook on mac: open_notebok_from_fidis'
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

