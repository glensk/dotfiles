#!/bin/sh

#function mcd () { mkdir -p $1 && cd $1 }  # does not work properly for bash

####################################################################
## ls alias
## -v sorts numerical
## -H follows symbolic links
## -F shows control chars like / or *
# This following ls command seems to work everywhere  (-v added to sort numbered files by true numerical number)
alias ls='ls -F --color=auto --group-directories-first --show-control-chars --hide="Icon?" -v'
#ls='ls -b -F -H -v --color=tty --show-control-chars --ignore="*.log" --group-directories-first'
#ls -F -H -v --color=tty --hide=Icon? --show-control-chars --ignore="*.log" --group-directories-first $*
#ls 'ls -F -H -v --color=tty'    
#################################################################
# i-pi 
#################################################################
#$1 i-pi-dev-venkat '$HOME/Dropbox/Albert/scripts/i-pi-dev/bin/i-pi'
#$1 i-pi-mc        '$HOME/Dropbox/Albert/scripts/i-pi-mc/bin/i-pi'
#################################################################
# alias navigation   use ' '  instead of ""  since only those work with `` correctly
#################################################################
alias gar='sshgar.sh'
alias finder-cd-to-current-finder-dir='cd `finder-list-current-finder-dir`'
#$1 n 'cd `tn.sh`'  # go to corresponding nas folder
#$1 d 'cd `td.sh`'  # go to corresponding data folder
#alias l='cd `ti_high_0_create_Folders_vasp.sh -l`'  # go to corresponding high folder
#alias h='cd `ti_low_0_create_Folders_vasp.sh -hit`'  # go to corresponding high folder
alias rl='source $HOME/.`echo $currentshell`rc'    # ist zsh unter zsh
alias c='j'    # for autojump

alias aliases='cd $dotfiles/aliases'
alias dl='cd $HOME/Downloads'
alias scripts='cd $HOME/Dropbox/Albert/scripts'
alias dot='cd $dotfiles'
alias pot='cd $dotfiles/scripts/potentials'
alias v='cd $HOME/Dropbox/Albert/v'
alias db='cd $HOME/Dropbox/Albert/'
alias kmc_='cd $HOME/Dropbox/Albert/kmc'
alias td='cd $HOME/Dropbox/Albert/Thermodynamics'

alias proj='cd $HOME/Dropbox/Albert/proj'
alias curr='cd $HOME/Dropbox/Albert/proj/proj_current'
alias paper='cd $HOME/Dropbox/Albert/proj/0000_paper_antraege'
alias prx='cd $HOME/Dropbox/Albert/proj/0000_paper_antraege/2013.12_PRX_Al_Cu_vakanz'
alias prl='cd $HOME/Dropbox/Albert/proj/0000_paper_antraege/2015.02_PRL_origin_anharmonicity'
alias diss='cd $HOME/Dropbox/Albert/proj/0000_paper_antraege/2015.12_Dissertation'
alias ptmp='cd /cmmc/ptmp/aglen/Understand_phonon_lifetimes'

alias vimrc='vi $HOME/.vimrc'
alias vip='vi $HOME/.ipython/profile_default/ipython_config.py'
alias virc='vi $HOME/.`echo $currentshell`rc'
#alias virca='vi $dotfiles/generalrc/generalrc_alias_.sh'
alias virca='vi $dotfiles/generalrc/aliases.sh'
alias vircg='vi $dotfiles/generalrc/generalrc_.sh'

alias -- -='cd -'
alias cd..='cd ..'
alias ..='cd ..'
alias ....='cd ../..'
alias ......='cd ../../..'
alias e='exit'

#$1 sxlp 'hier=`pwd`;cd ~/Thermodynamics/python_thermodynamics/save_old; cp ../lammps_pos_to_sum.py .; sendx lammps_pos_to_sum.py;cd $hier'
#$1 rxlp 'hier=`pwd`;cd ~/Thermodynamics/python_thermodynamics/save_old; recievex; cp lammps_pos_to_sum.py .. ;cd $hier'



##############################
# system
##############################
alias clearpath='source $dotfiles/aliases/clearpath.sh'
alias ra='$dotfiles/generalrc/generalrc_alias_renew.sh'
alias la='ls -la'
alias ll='ls -la'
alias sl='ls'
#$1 rm 'rm -rf'   # too general and maybe too dangerous?
#$1 vi 'vim -v'   # is way too slow, better use plain vi
#$1 vim 'vim -v'  # is way too slow, better use plain vi

alias gu='gitup.sh'
#$1 gd 'gitdown.sh'
alias gd='gitup.sh'
alias ggd='git clone --recursive https://github.com/glensk/dotfiles.git'
alias gurl='git config --get remote.origin.url'
alias gs='git status -u'
alias ss='svn status -u'

alias fd='find . -type d -name'
alias ff='find . -type f -name'
#alias run='./run.sh'
#alias s='./run.sh'
alias pu='pushd `pwd`'
alias po='popd'
alias extract='$dotfiles/aliases/extract.sh'
alias monthly_update='cd $dotfiles/cron;./monthlyupdate.sh'


#$1 du 'echo ______better ncdu______;echo; du'
alias lstarfile='tar -ztvf'
#$1 diff 'colordiff'
alias ep='epstopdf'
alias grep='grep --color=auto -d skip'    # skips directories
#$1 grep'echo use ack instead'
alias bc='bc -l'
#alias cp='cp -a'  # this one did not copy symliks as files but as symlinks -> remove -a
alias xmgrace='DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -param ~/.xmgracerc'
alias x='DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -param ~/.xmgracerc'
alias xll="DISPLAY=:0.0 xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -param $dotfiles/xmgrace/tpl_log_log.par"
alias g='gnuplot.py'
alias gp='gnuplot.py'
alias untargz='tar -xvf'
alias untar='tar xfv'
alias untarbz2='tar -xf'
alias ungz='gzip -d'
alias untgz='tar -xvf'
alias unz='uncompress'
alias unbz2='tar -jxvf'
alias killall='killall -9'

#$1 gnuplot 'gnuplot -persist'
alias i='ipython'
#$1 load_last_saved_settings='$dotfiles/mac_KeyRemap4MacBook/load_last_saved_settings.sh'
alias cpdf='$dotfiles/bin/cpdf-binaries/OSX-Intel/cpdf'
alias df='df -h' # | grep $USER'
alias sp='settitlepath.sh'
#$1 con 'google contacts list --fields name,email,phone_number,address,birthday --title "(?i).*\!*"'
alias pingtest='ping 69.20.11.136'
alias pingtestgoogle='ping 173.194.69.102'   # ping google.com   # this also get you directly to google.com in the browser without DNS\
alias pingmit='ping 18.62.0.96'				# this can be directly opened without DNS in google/firefox
#$1 echor 'echo "\033[1;31m\!*\033[m"'  # problems with zsh
#$1 echog 'echo "\033[1;32m\!*\033[m"'  # problems with zsh    
#$1 echob 'echo "\033[1;34m\!*\033[m"'  # problems with zsh    
#$1 echoB 'echo "\033[1m\!*\033[m"'     # problems with zsh    
#$1 xc 'xclock -d -brief -padding -8 -geometry 1x1-0+999'
#$1 undo 'undo -i'

alias memHogsTop='top -l 1 -o rsize | head -20'
alias memHogsPs='ps wwaxm -o pid,stat,vsize,rss,time,command | head -10'
alias cpu_hogs='ps wwaxr -o pid,stat,%cpu,time,command | head -10'


##############################
# module 
##############################
alias ma='module avail'
alias ml='module load'
alias mli='module list'
alias mu='module unload'


##############################
# que
##############################
alias qls='qls.sh -u'
alias qstatn='qstatn -d'


##################################################################################
# host specific
##################################################################################
##############################
# on mac
##############################
if [ "$onhost" = "mac" ];then
    alias math='/Applications/Mathematica.app/Contents/MacOS/MathKernel' #&& \
    alias mathematica='/Applications/Mathematica.app/Contents/MacOS/Mathematica' #&& \
    #$1 xmgrace '$HOME/scripts/mac_tools/apps/xmgrace/grace-5.1.23_bigbuf/src/xmgrace' && \
    alias top='top -o cpu' #&& \
    alias trash="rmtrash" #&& \
    alias del "rmtrash" #&& \
    #$1 rm 'echo Use del command, or the full path i.e. /bin/rm' && \
    #$1 rm 'del' && \
    alias ctagsnew='/usr/local/Cellar/ctags/5.8/bin/ctags -R .' && \
    alias mvim='open /Applications/MacVim.app' 
    #$1 edit "open -a MacVim.app $1"
    #$1 vi       'vim' 
fi

    #alias sed='$dotfiles/bin/sed/sed-4.2/build/bin/sed-4.2'
    #alias mvim='open /Applications/MacVim.app'
    #alias cat 'colorize'  # colorize seems to have sometimes problems

if [ "$onhost" = "cmmc" ];then
    alias garblazej='ssh bzg@cmmc001.bc.rzg.mpg.de' && \
    alias scpfromcmmc='scp -r aglen@cmmc001.bc.rzg.mpg.de:\!:1 \!:2' && \
    alias scptocmmc='scp -r \!:1 aglen@cmmc001.bc.rzg.mpg.de:\!:2' && \
    alias s4='~/start.4.6.28.tid.quick.sh' && \
    alias v4='~/start.4.6.28.tid.quick.sh' && \
    alias s5='~/start.5.3.5.20cores.sh' && \
    alias s52='~/start.5.3.5.20cores.sh' && \
    alias v5='~/start.5.3.5.20cores.sh' && \
    alias v52='~/start.5.3.5.20cores.sh' && \
    alias tmux='/u/aglen/local/bin/tmux' && \
    alias tar='gtar'
fi
#$1 vi 'gvim -v' && \
#$1 vim 'gvim -v' && \

#[ "$onhost" = "cmpc" ] && \
#    $1 tmux '/home/glensk/local/bin/tmux'

#$1 vi 'gvim -v' && \
#$1 vim 'gvim -v' && \

#[ "$host" = "cmmc001.mpie.de" ] && $1 vi 'vi -v' && $1 vim 'vi -v' # check if cmmc has vim or only vi

##############################
# ssh
##############################
alias ssh="ssh.sh -Y -X -o ServerAliveInterval=1600 -o ServerAliveCountMax=1200"   # -X for loading locas stuff locally (xmgrace)  ServerAliveInterval=1600 is now enough for cmpc  ... Blazej uses ServerAliveCountMax=1200 which seems to be good
alias sshimmer='ssh.sh -Y -X -o ServerAliveInterval=1600 -o ServerAliveCountMax=1200 '

alias pc='sshimmer -t glensk@cosmopc15.epfl.ch   "[ -e `th.sh` ] && cd `th.sh`; zsh"'
alias cosmopc='sshimmer -t glensk@cosmopc15.epfl.ch   "[ -e `th.sh` ] && cd `th.sh`; zsh"'
#$1 cosmopc_tunnel   'sshimmer -N -L 8080:localhost:8080 -t glensk@cosmopc15.epfl.ch   "[ -e `th.sh` ] && cd `th.sh`; zsh"'  # open this first

alias cluster='sshimmer -t glensk@cmmc002.mpie.de     "[ -e `th.sh` ] && cd `th.sh`; zsh"'
alias cmpc='sshimmer -t glensk@cmpc34.mpie.de      "[ -e `th.sh` ] && cd `th.sh`; zsh"' 
alias gate='sshimmer -t aglen@gate.rzg.mpg.de      "[ -e `th.sh` ] && cd `th.sh`; zsh"'
alias gatelinux='sshimmer -t aglen@woese.opt.rzg.mpg.de "[ -e `th.sh` ] && cd `th.sh`; zsh"'
alias gatehydra='sshimmer -t aglen@hydra.rzg.mpg.de     "[ -e `th.sh` ] && cd `th.sh`; zsh"'

alias fidis='ssh fidis'
alias f='ssh fidis'
alias daint='ssh daint'
alias hel='ssh helvetios'
alias h='ssh helvetios'
alias cosmopc='ssh cosmopc'
alias pc='ssh cosmopc'

##############################
# autocorrect
##############################
alias moduel='module'

################################
# Slurm que jobsceduler
# #############################
#$1 mount_fidis 'sshfs -o idmap=user glensk@fidis@epfl.ch:/scratch/glensk ~/fidis'
# removed: -o reconnect

alias mount_fidis_scratch='sshfs glensk@fidis.epfl.ch:/scratch/glensk/ /scratch/glensk -o reconnect -C; echo mounted /scratch/glensk' 
alias umount_fidis_scratch='umount -f /scratch/glensk; echo umounted /scratch/glensk'
alias mount_fidis_home='sshfs glensk@fidis.epfl.ch:/home/glensk/ /home/glensk -o reconnect -C; echo mounted /home/glensk' 
alias umount_fidis_home='umount -f /home/glensk; echo umounted /home/glensk'
alias mount_cosmopc_home='sshfs glensk@cosmopc15.epfl.ch:/home/glensk/ /home/glensk -o reconnect -C; echo mounted /home/glensk' 
alias umount_cosmopc_home='umount -f /home/glensk; echo umounted /home/glensk'
alias mount_pc_scratch='sshfs glensk@cosmopc15.epfl.ch:/local/scratch/glensk/ /local/scratch/glensk -o reconnect -C; echo mounted /local/scratch/glensk/'
alias umount_pc_scratch='umount -f /local/scratch/glensk'
alias mount_pc_home='sshfs glensk@cosmopc15.epfl.ch:/home/glensk /home/glensk -o reconnect -C; echo mounted /home/glensk/'
alias umount_pc_home='umount -f /home/glensk'
alias mount_daint='sshfs aglensk@ela.cscs.ch:/users/aglensk/ ~/daint -o reconnect -C; echo mounted '
alias umount_daint='umount -f /users/aglensk'
alias make_home_accessible='chmod -R ga+r $HOME'

#$1 conda_activate "[ -e "$HOME/miniconda2" ] && source $HOME/miniconda2/etc/profile.d/conda.sh && conda activate"
alias gap="conda activate gap"
alias base="conda activate base"
alias aiida="source $HOME/aiida/bin/activate"
alias s='cd $SCRATCH'
alias aiida_activate="source $HOME/aiida/bin/activate"
alias aiida_deactivate="deactivate"
alias jupyter_port='jupyter notebook --no-browser --port=8080'
