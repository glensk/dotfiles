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

setenv () {
      export $1=$2
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

if [ "$ZSH_VERSION" = "4.3.6" ];then
module () {
    eval `/afs/ipp/common/usr/modules.2014/amd64_sles11/Modules/$MODULE_VERSION/bin/modulecmd sh $*`
}
fi


tab-color() {
    [ "$1" = "red" ] && a=255 b=0 c=0
    [ "$1" = "gren" ] && a=0 b=255 c=0
    [ "$1" = "blue" ] && a=0 b=0 c=255
    [ "$1" = "cyan" ] && a=64 b=244 c=208
    [ "$1" = "magenta" ] && a=255 b=0 c=255
    #[ "$1" = "turquoise" ] && a=100 b=100 c=208
echo -ne "\033]6;1;bg;red;brightness;$a\a"
echo -ne "\033]6;1;bg;green;brightness;$b\a"
echo -ne "\033]6;1;bg;blue;brightness;$c\a"
}

