if [ "`hostname`" = "fidis" ];then
    echo module purge
    module purge
    if [ -e "$HOME/miniconda2" ];then
        echo source miniconda2
        source $HOME/miniconda2/etc/profile.d/conda.sh 
        echo conada activate
        conda activate
    fi
fi

