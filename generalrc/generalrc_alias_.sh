#!/bin/sh

##############################################
# ALIASES
##############################################
# how to create a new alias (which works on every host)
# 1. shell: virca  --> add the alias you need to this file here
# 2. to renew aliasfile: ra   (this only creates a new alias file)
# 3. START A NEW SHELL (also you typed ra!) 




#[ "$printloadstat" = "true" ] && echo ... alias before && alias && echo ... alias before done
#echo ------------------- SHELL $currentshell WHICH $whichalias 
#[ "$currentshell" = "tcsh" ] && kkk="$alias"
#[ "$currentshell"  = "bash" ] && unalias -a
#whichalias=$1
#currentshell=$2

####################################################################
## ls alias
## -v sorts numerical
## -H follows symbolic links
## -F shows control chars like / or *
# This following ls command seems to work everywhere  (-v added to sort numbered files by true numerical number)
$1 ls 'ls -F --color=auto --group-directories-first --show-control-chars --hide="Icon?" -v'
#ls='ls -b -F -H -v --color=tty --show-control-chars --ignore="*.log" --group-directories-first'
#ls -F -H -v --color=tty --hide=Icon? --show-control-chars --ignore="*.log" --group-directories-first $*
#ls 'ls -F -H -v --color=tty'    
#################################################################
# i-pi 
#################################################################
$1 i-pi-dev-venkat '$HOME/google_drive/scripts_epfl/i-pi-dev/bin/i-pi'
$1 i-pi-kmc '$HOME/google_drive/scripts_epfl/i-pi-mc/bin/i-pi'
#################################################################
# alias navigation   use ' '  instead of ""  since only those work with `` correctly
#################################################################
$1 gar 'sshgar.sh'
$1 finder-cd-to-current-finder-dir 'cd `finder-list-current-finder-dir`'
$1 n 'cd `tn.sh`'  # go to corresponding nas folder
$1 d 'cd `td.sh`'  # go to corresponding data folder
$1 l 'cd `ti_high_0_create_Folders_vasp.sh -l`'  # go to corresponding high folder
$1 h 'cd `ti_low_0_create_Folders_vasp.sh -hit`'  # go to corresponding high folder
$1 rl 'source $HOME/.`echo $currentshell`rc'    # ist zsh unter zsh

$1 scripts 'cd $HOME/Dropbox/Albert/scripts'
$1 dot     'cd $dotfiles'
$1 v       'cd $HOME/Dropbox/Albert/v'


$1 proj  'cd $HOME/Dropbox/Albert/proj'
$1 curr  'cd $HOME/Dropbox/Albert/proj/proj_current'
$1 paper 'cd $HOME/Dropbox/Albert/proj/0000_paper_antraege'
$1 prx   'cd $HOME/Dropbox/Albert/proj/0000_paper_antraege/2013.12_PRX_Al_Cu_vakanz'
$1 prl   'cd $HOME/Dropbox/Albert/proj/0000_paper_antraege/2015.02_PRL_origin_anharmonicity'
$1 diss  'cd $HOME/Dropbox/Albert/proj/0000_paper_antraege/2015.12_Dissertation'
$1 ptmp  'cd /cmmc/ptmp/aglen/Understand_phonon_lifetimes'

$1 vimrc 'vi $HOME/.vimrc'
$1 vip   'vi $HOME/.ipython/profile_default/ipython_config.py'
$1 virc  'vi $HOME/.`echo $currentshell`rc'
$1 virca 'vi $dotfiles/generalrc/generalrc_alias_.sh'

$1 cd.. 'cd ..'
$1 .. 'cd ..'
$1 .... 'cd ../..'
$1 ...... 'cd ../../..'

#$1 rx 'recievex'
$1 rx '$dotfiles/bin/stefan/recievex'
$1 rxo '$dotfiles/bin/stefan/recievexoxford'
$1 sx 'sendx'
$1 sxlp 'hier=`pwd`;cd ~/Thermodynamics/python_thermodynamics/save_old; cp ../lammps_pos_to_sum.py .; sendx lammps_pos_to_sum.py;cd $hier'
$1 rxlp 'hier=`pwd`;cd ~/Thermodynamics/python_thermodynamics/save_old; recievex; cp lammps_pos_to_sum.py .. ;cd $hier'



##############################
# system
##############################
$1 clearpath 'source $dotfiles/bin/clearpath.sh'
$1 ra '$dotfiles/generalrc/generalrc_alias_renew.sh'
$1 la 'ls -la'
$1 ll 'ls -la'
$1 sl 'ls'
#$1 rm 'rm -rf'   # too general and maybe too dangerous?
#$1 vi 'vim -v'   # is way too slow, better use plain vi
#$1 vim 'vim -v'  # is way too slow, better use plain vi

$1 gu 'gitup.sh'
$1 gd 'gitdown.sh'
$1 ggd 'git clone --recursive https://github.com/glensk/dotfiles.git'
$1 gurl 'git config --get remote.origin.url'
$1 gs 'git status -u'
$1 ss 'svn status -u'

$1 fd 'find . -type d -name'
$1 ff 'find . -type f -name'

$1 run './run.sh'
$1 s './run.sh'
$1 pu 'pushd `pwd`'
$1 po 'popd'

$1 extract '$dotfiles/bin/extract.sh'

$1 monthly_update 'cd $dotfiles/cron;./monthlyupdate.sh'


$1 du 'echo ______better ncdu______;echo; du'
$1 lstarfile 'tar -ztvf'
#$1 diff 'colordiff'
$1 ep 'epstopdf'
$1 grep 'grep --color=auto -d skip'    # skips directories
#$1 grep'echo use ack instead'
$1 bc 'bc -l'
#alias cp='cp -a'  # this one did not copy symliks as files but as symlinks -> remove -a
$1 xmgrace 'xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -param ~/.xmgracerc'
$1 x 'xmgrace -maxpath 1000000 -geom 1100x860 -nosigcatch -param ~/.xmgracerc'
$1 g gnuplotfile
$1 untargz 'tar -xvf'
$1 untar 'tar xfv'
$1 untarbz2 'tar -xf'
$1 ungz 'gzip -d'
$1 untgz 'tar -xvf'
$1 unz 'uncompress'
$1 unbz2 'tar -jxvf'
$1 killall 'killall -9'

$1 gnuplot 'gnuplot -persist'
$1 i 'ipython'
#$1 load_last_saved_settings='$dotfiles/mac_KeyRemap4MacBook/load_last_saved_settings.sh'
$1 cpdf '$dotfiles/bin/cpdf-binaries/OSX-Intel/cpdf'
$1 df 'df -h' # | grep $USER'
$1 sp 'settitlepath.sh'
#$1 con 'google contacts list --fields name,email,phone_number,address,birthday --title "(?i).*\!*"'
$1 pingtest 'ping 69.20.11.136'
$1 pingtestgoogle 'ping 173.194.69.102'   # ping google.com   # this also get you directly to google.com in the browser without DNS\
$1 pingmit 'ping 18.62.0.96'				# this can be directly opened without DNS in google/firefox
#$1 echor 'echo "\033[1;31m\!*\033[m"'  # problems with zsh
#$1 echog 'echo "\033[1;32m\!*\033[m"'  # problems with zsh    
#$1 echob 'echo "\033[1;34m\!*\033[m"'  # problems with zsh    
#$1 echoB 'echo "\033[1m\!*\033[m"'     # problems with zsh    
#$1 xc 'xclock -d -brief -padding -8 -geometry 1x1-0+999'
#$1 undo 'undo -i'

$1 memHogsTop 'top -l 1 -o rsize | head -20'
$1 memHogsPs 'ps wwaxm -o pid,stat,vsize,rss,time,command | head -10'
$1 cpu_hogs 'ps wwaxr -o pid,stat,%cpu,time,command | head -10'


##############################
# module 
##############################
$1 ma 'module avail'
$1 ml 'module load'
$1 mli 'module list'
$1 mu 'module unload'


##############################
# que
##############################
$1 qls 'qls.sh -u'
$1 qstatn 'qstatn -d'


##################################################################################
# host specific
##################################################################################
##############################
# on mac
##############################
[ "$onmac" = "true" ] && \
    $1 math '/Applications/Mathematica.app/Contents/MacOS/MathKernel' && \
    $1 mathematica   '/Applications/Mathematica.app/Contents/MacOS/Mathematica' && \
    #$1 xmgrace '$HOME/scripts/mac_tools/apps/xmgrace/grace-5.1.23_bigbuf/src/xmgrace' && \
    $1 units '$dotfiles/bin/units/units-1.88/units' && \
    $1 top 'top -o cpu' && \
    $1 trash "rmtrash" && \
    $1 del "rmtrash" && \
    #$1 rm 'echo Use del command, or the full path i.e. /bin/rm' && \
    $1 rm 'del' && \
    $1 ctagsnew '/usr/local/Cellar/ctags/5.8/bin/ctags -R .' && \
    $1 mvim     'open /Applications/MacVim.app' 
    #$1 edit "open -a MacVim.app $1"
    #$1 vi       'vim' 


    #alias sed='$dotfiles/bin/sed/sed-4.2/build/bin/sed-4.2'
    #alias mvim='open /Applications/MacVim.app'
    #alias cat 'colorize'  # colorize seems to have sometimes problems

# cmmd does not exist anymore
#[ "$oncmmd" != "true" ] && \
#    $1 qhost 'qhost.sh' && \
#    $1 qstat 'qstat.sh' && \
#    $1 qdel 'qdel.sh'   && \
#    $1 qalter 'qalter.sh' && \
#    $1 qsub 'qsub.sh'   && \
#    $1 qls 'qls.garching -u'  && \
#    $1 qhold 'qhold.sh' && \
#    $1 qrls 'qrls.sh'         

[ "$oncmmc" = "true" ] && \
    $1 garblazej 'ssh bzg@cmmc001.bc.rzg.mpg.de' && \
    $1 scpfromcmmc 'scp -r aglen@cmmc001.bc.rzg.mpg.de:\!:1 \!:2' && \
    $1 scptocmmc 'scp -r \!:1 aglen@cmmc001.bc.rzg.mpg.de:\!:2' && \
    $1 s4 '~/start.4.6.28.tid.quick.sh' && \
    $1 v4 '~/start.4.6.28.tid.quick.sh' && \
    $1 s5 '~/start.5.3.5.20cores.sh' && \
    $1 s52 '~/start.5.3.5.20cores.sh' && \
    $1 v5 '~/start.5.3.5.20cores.sh' && \
    $1 v52 '~/start.5.3.5.20cores.sh' && \
    $1 tmux '/u/aglen/local/bin/tmux' && \
    $1 tar 'gtar'

#$1 vi 'gvim -v' && \
#$1 vim 'gvim -v' && \

[ "$oncmpc" = "true" ] && \
    $1 tmux '/home/glensk/local/bin/tmux'

#$1 vi 'gvim -v' && \
#$1 vim 'gvim -v' && \

[ "$host" = "cmmd001.mpie.de" ] && $1 vi 'vi -v' && $1 vim 'vi -v'
#[ "$host" = "cmmc001.mpie.de" ] && $1 vi 'vi -v' && $1 vim 'vi -v' # check if cmmc has vim or only vi

##############################
# ssh
##############################
$1 ssh 'ssh.sh -Y -X -o ServerAliveInterval=1600 -o ServerAliveCountMax=1200'   # -X for loading locas stuff locally (xmgrace)  ServerAliveInterval=1600 is now enough for cmpc  ... Blazej uses ServerAliveCountMax=1200 which seems to be good
$1 sshimmer 'ssh.sh -Y -X -o ServerAliveInterval=1600 -o ServerAliveCountMax=1200'

$1 cluster 'sshimmer -t glensk@$submithost         "[ -e `th.sh` ] && cd `th.sh`; zsh"'
$1 001 'sshimmer -t glensk@cmmd001.mpie.de         "[ -e `th.sh` ] && cd `th.sh`; zsh"'
$1 002 'sshimmer -t glensk@cmmd002.mpie.de         "[ -e `th.sh` ] && cd `th.sh`; zsh"'
$1 cmpc 'sshimmer -t glensk@$myhost.mpie.de         "[ -e `th.sh` ] && cd `th.sh`; zsh"' # this loads tcsh twice due to the necessary tcsh at the end; once zsh is installed on the cmpc (and cmmd?) the tcsh can be left empty and instead of tcsh one could simply write zsh

$1 gate 'sshimmer -t aglen@gate.rzg.mpg.de      "[ -e `th.sh` ] && cd `th.sh`; zsh"'
$1 gatelinux 'sshimmer -t aglen@woese.opt.rzg.mpg.de "[ -e `th.sh` ] && cd `th.sh`; zsh"'
$1 gatehydra 'sshimmer -t aglen@hydra.rzg.mpg.de     "[ -e `th.sh` ] && cd `th.sh`; zsh"'

#$1 gar 'sshimmer -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh` ]  && cd `th.sh`; zsh"'
#$1 gar1 'sshimmer -t aglen@cmmc001.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh` ] && cd `th.sh`; zsh"'
#$1 gar2 'sshimmer -t aglen@cmmc002.bc.rzg.mpg.de -R 48540:cmcc1.mpie.de:80 "[ -e `th.sh` ] && cd `th.sh`; zsh"'

##############################
# autocorrect
##############################
$1 moduel 'module'

################################
# Slurm que jobsceduler
# #############################
#$1 q 'squeue -u "$USER" --format "%A %t %N %C %M %Z" | awk '{printf "%8s %3s %20s %4s %6s %10s\n",$1,$2,$3,$4,$5,$6}''
#$1 mount_fidis 'sshfs -o idmap=user glensk@fidis@epfl.ch:/scratch/glensk ~/fidis'
$1 mount_fidis 'sshfs glensk@fidis.epfl.ch:/scratch/glensk/ /scratch/glensk -o reconnect -C' 
$1 umount_fidis 'umount -f /scratch/glensk'
$1 mount_daint 'sshfs aglensk@ela.cscs.ch:/users/aglensk/ ~/daint -o reconnect -C'
$1 make_home_accessible 'chmod -R ga+r $HOME'

$1 getconda ". $HOME/miniconda3/etc/profile.d/conda.sh && conda activate"
