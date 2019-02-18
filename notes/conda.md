# to activate conda on mac: 
source $HOME/miniconda2/etc/profile.d/conda.sh && conda activate


conda update conda
conda info --envs           # Conda show available virtual environments
conda env list              # Conda show available virtual environments
conda activate base         # active a certain virtual environment (here base)
conda activate python2

conda deactivate            # deactivated a certain virtual environment
conda install package       # will install package in the current environment
