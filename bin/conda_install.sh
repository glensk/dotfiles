
[ ! -e "$HOME/Downloads" ] && mkdir $HOME/Downloads 
[ ! -e "$HOME/Downloads" ] && echo not folder $HOME/Downloads && exit
cd $HOME/Downloads
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -f 

echo "now source"
source $HOME/miniconda3/etc/profile.d/conda.sh && conda
#source ~/.zshrc
echo "now conda config"
conda config --add channels intel    # this makes probs;
echo "now conda activate"
conda activate

echo "SHOW WHICH CONDA ENVIRONMENT!! cant see it when installing"
conda info --envs
echo "now install stuff CURRENTLY NEED TO ACCEPT EVERYTHIN WHITH IS NOT GOOD"
conda install -y colorama
conda install -y click  # actually comes with ase
conda install -y -c conda-forge ase
conda install -y seaborn
conda install -y -c conda-forge lmfit # for phonon_lifetimes.py kette.py ...
[ "`hostname`" = "mac" ] && conda install -y -c conda-forge jupyter_contrib_nbextensions   # get the notebook extensions for jupyter notebooks
[ "`hostname`" = "mac" ] && jupyter contrib nbextension install --user                  # also necessary to get the notebook extensions working

#conda install -c conda-forge pyfftw=0.10.4  # is this already installed using 
#conda install argcomplete           # to get argcompletion of python scripts in bash/zsh

