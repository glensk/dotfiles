#!/bin/sh


#echo $PATH
PATH="$(echo $PATH | perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, scalar <>))')"
PYTHONPATH="$(echo $PYTHONPATH | perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, scalar <>))')"
LD_LIBRARY_PATH="$(echo $LD_LIBRARY_PATH | perl -e 'print join(":", grep { not $seen{$_}++ } split(/:/, scalar <>))')"

#echo
#echo $PATH
export PATH="$PATH"
export PYTHONPATH="$PYTHONPATH"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH"
