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

