#!/bin/sh

host=`hostname`
if [[ "$host" = "cmmc002" || "$host" == "cmmc001" ]];then
    echo "in",`hostname`
    cd /cmmc/ptmp/aglen/Understand_phonon_lifetimes
    pwd
    sendx `find . -type f \( -name "lifetimesgood_*" -o -name "freqsgood_*" \) | xargs`
    echo "done"
    fi


if [ "$host" = "mac" ];then
    echo "in",`hostname`
    if [ -e "/Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_20_python_plotting" ];then
        cd /Users/glensk/Dropbox/proj/proj_current/__2017.01_phonon_linewidth_al/__2017.01_phonon_lifetimes_4_ab2017/check_20_python_plotting
        echo Password for cmmc002:
        recievex
    fi
fi
