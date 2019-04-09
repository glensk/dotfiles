# .bashrc
#echo "Hello moto"
# Source global definitions
if [ -f /etc/bashrc ]; then
  . /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=
export LD_LIBRARY_PATH=/home/lopanits/source/cosmo-lammps/lib/nnp/lib:${LD_LIBRARY_PATH}
#User specific aliases and functions
alias sme='squeue -u $USER -o "%.18i %.9P %.8u %.2t %.10M %.6D %R %.20j"'
alias slist='tsqueue -o "%.8i %4t %.3D %.8g %.12u %.20S %.8Q %.20j"'
alias keepcalm='for i in `seq 1 100`; do sme; sleep 1; done'

