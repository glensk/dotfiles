if [ "`hostname`" = "fidis" ];then
    echo fidis
    echo module purge
    module purge
    if [ -e "$HOME/miniconda2" ];then
        echo source miniconda2
        source $HOME/miniconda2/etc/profile.d/conda.sh 
        echo conada activate
        conda activate
    fi
fi

if [ "`hostname | grep -o cosmopc`" = "cosmopc" ];then
    echo cosmopc
    if [ -e "$HOME/miniconda3" ];then
        echo source miniconda3
        source $HOME/miniconda3/etc/profile.d/conda.sh 
        echo conada activate
        conda activate
    fi
fi
