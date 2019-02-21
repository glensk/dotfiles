1) on mac: (open the tunnel)
shimmer -N -L 8080:localhost:8080 -t glensk@cosmopc15.epfl.ch   "[ -e `th.sh` ] &&
cd `th.sh`; zsh"'  # open this first on mac (maybe it is possible to put this in
the background

2) on cosmopc, open the notebook:
jupyter-notebook SOAPS.ipynb
