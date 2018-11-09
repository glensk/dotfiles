#!/bin/sh
# everything works perfect even mcd k/k/k
# for bash transferability ls; needs the ; at the end
#function cd { builtin cd "$@" ; settitlepath.sh; ls; } # not necassary when preexec() { ODIR="$(pwd)" } and precmd() { [[ "$(pwd)" != $ODIR ]] && ls } are set in zsh_set
function mcd () { mkdir -p $1 && cd $1 }
alias -- -='cd -' #alias \\-='cd -' 
alias "../"="cd .."
