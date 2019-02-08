
condaversion=2 # currently I dont see a reason to use python3


[ ! -e "$HOME/Downloads" ] && mkdir $HOME/Downloads 
[ ! -e "$HOME/Downloads" ] && echo not folder $HOME/Downloads && exit
cd $HOME/Downloads
Linux_or_MacOSX=`uname -a | awk '{print $1}'`
if [ "$Linux_or_MacOSX" == "Linux" ];then
    Linux_or_MacOSX="Linux"
elif [ "$Linux_or_MacOSX" == "Darwin" ];then 
    Linux_or_MacOSX="MacOSX"
else
    uname -a
    echo "is neither Linux nor Darwin" 
    exit
fi
RED='\033[0;31m'
#RED="\033[31m"
NC="\033[0m" # No Color

echo -e "${RED}###############################################################${NC}"
echo -e "${RED}## try go get everything running currently with python2 #######${NC}"
echo -e "${RED}###############################################################${NC}"
echo -e "${RED}Linux_or_MacOSX: $Linux_or_MacOSX"
echo -e "condaversion   : $condaversion${NC}"

if [ "$condaversion" == "3" ];then
wget https://repo.continuum.io/miniconda/Miniconda3-latest-$Linux_or_MacOSX-x86_64.sh
bash Miniconda3-latest-$Linux_or_MacOSX-x86_64.sh -b -f 
fi

if [ "$condaversion" == "2" ];then
wget https://repo.continuum.io/miniconda/Miniconda2-latest-$Linux_or_MacOSX-x86_64.sh
bash Miniconda2-latest-$Linux_or_MacOSX-x86_64.sh -b -f 
fi


echo "${RED}now source${NC}"
source $HOME/miniconda$condaversion/etc/profile.d/conda.sh && conda
#source ~/.zshrc

echo 
echo "${RED}now conda add chennels intel${NC}"
# actually, my mac/cosmopc and all the clusters are intel based (clusters are all Xeon, cosmopc is i5)
conda config --add channels intel    

echo
echo "${RED}now conda activate${NC}"
conda activate

echo ""
echo "${RED}SHOW WHICH CONDA ENVIRONMENT!! cant see it when installing${NC}"
conda info --envs
echo "now install stuff CURRENTLY NEED TO ACCEPT EVERYTHIN WHITH IS NOT GOOD"
echo "${RED}Snow install stuff CURRENTLY NEED TO ACCEPT EVERYTHIN WHITH IS NOT GOOD${NC}"
echo
echo "${RED}scikit-learn${NC}"
conda install -y scikit-learn # for analyzing lifangs ipynb of correlationis
echo "${RED}colorame${NC}"
conda install -y colorama
echo "${RED}click${NC}"
conda install -y click       # actually comes with ase
echo "${RED}ase${NC}"
conda install -y -c conda-forge ase
echo "${RED}seaborn${NC}"
conda install -y seaborn
echo "${RED}lmfit${NC}"
conda install -y -c conda-forge lmfit # for phonon_lifetimes.py kette.py ...
echo "${RED}jupyter_contrib_nbextension${NC}"
conda install -y -c conda-forge jupyter_contrib_nbextensions   # get the notebook extensions for jupyter notebooks
echo "${RED}jupyter contrib nbextension install --user${NC}"
jupyter contrib nbextension install --user # also necessary to get the notebook extensions working
echo "${RED}nglview${NC}"  # to look at structures
conda install -c omnia nglview
conda install -y -c conda-forge atomsk  # checking for unwrapping atomic structures

conda install -y -c numba # tried for faster SF_integrate with CurSel
#conda install -c conda-forge pyfftw=0.10.4  # is this already installed using 
#conda install argcomplete           # to get argcompletion of python scripts in bash/zsh

