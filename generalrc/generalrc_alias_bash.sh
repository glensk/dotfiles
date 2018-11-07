#!/bin/sh

# everything works perfect even mcd k/k/k
# for bash transferability ls; needs the ; at the end
#function cd { builtin cd "$@" ; settitlepath.sh; ls; };
#function cd { builtin cd "$@" ; ls; };  # even this seems to make probs with aiida verdi computer test daint
function mcd () { mkdir -p $1;cd $1; }
alias -- -='cd -' #alias \\-='cd -' 
function .. { cd ..; }
function ../ { cd ..; }
function ../.. { cd ../..; }
function ../../ { cd ../..; }
