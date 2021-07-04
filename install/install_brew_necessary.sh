#!/bin/sh


######################################################################################
######################################################################################
# only execute rarely
######################################################################################
######################################################################################
# check if necessary stuff installed
brewinstall () {
    if [ "`type $1`" = "$1 not found" ];then
        type $1
        brew install $1
    else 
        echo "                         .... $1 installed"
    fi
}

brewinstall2 () {
    if [ "`type $1`" = "$1 not found" ];then
        type $1
        brew install $2
    else 
        echo "                         .... $1 installed (/opt/homebrew/bin/)"
    fi
}

    
echo "checking what to install .... (2)"
if [ "true" = "true" ];then

    if [ "`type brew`" = "brew not found" ];then
        type brew
        exit
    else
        echo "                         .... brew installed"
    fi
    brewinstall wget
    #brew install wget --with-iri
    brewinstall ncdu
    brewinstall gls
    brewinstall git 
    git config --global credential.helper osxkeychain
    brewinstall lbzip2 # to parallel zip files or try install_git.py -i lbzip2 

    # to use the coreutis, the have (now) to be put in the path
    brewinstall gawk      # gawk
    brewinstall2 gtar gnu-tar   # PATH="/opt/homebrew/opt/gnu-tar/libexec/gnubin:$PATH"
    brewinstall2 gsed gnu-sed   # PATH="/opt/homebrew/opt/gnu-sed/libexec/gnubin:$PATH"
    brewinstall coreutils       # PATH="/opt/homebrew/opt/coreutils/libexec/gnubin:$PATH"
    brewinstall findutils       # PATH="/opt/homebrew/opt/findutils/libexec/gnubin:$PATH"
    brewinstall gnu-indent      # PATH="/opt/homebrew/opt/gnu-indent/libexec/gnubin:$PATH"
    brewinstall gnutls       
    brewinstall grep            # PATH="/opt/homebrew/opt/grep/libexec/gnubin:$PATH"
    brewinstall colordiff
    brewinstall youtube-dl
    #brew install gnuplot --with-aquaterm
fi
