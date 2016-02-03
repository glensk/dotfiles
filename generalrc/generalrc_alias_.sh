#!/bin/sh


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
# This following ls command seems to work everywhere
$1 ls 'ls -F --color=auto --group-directories-first --show-control-chars --hide="Icon?"'
#ls='ls -b -F -H -v --color=tty --show-control-chars --ignore="*.log" --group-directories-first'
#ls -F -H -v --color=tty --hide=Icon? --show-control-chars --ignore="*.log" --group-directories-first $*
#ls 'ls -F -H -v --color=tty'    

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

$1 prx 'cd $HOME/Dropbox/proj/0000_paper_antraege/2013.12_PRX_Al_Cu_vakanz'
$1 paper 'cd $HOME/dropbox/proj/0000_paper_antraege/2015.02_PRL_origin_anharmonicity'
$1 pap 'paper'
$1 prl 'paper'
$1 diss 'cd $HOME/Dropbox/proj/0000_paper_antraege/2015.12_Dissertation'

$1 vimrc 'vi $HOME/.vimrc'
$1 vip 'vi $HOME/.ipython/profile_default/ipython_config.py'
$1 virc 'vi $HOME/.`echo $currentshell`rc'
$1 virca 'vi $HOME/Dropbox/scripts/dotfiles/generalrc/generalrc_alias_.sh'

$1 cd.. 'cd ..'
$1 .. 'cd ..'
$1 .... 'cd ../..'
$1 ...... 'cd ../../..'



##############################
# system
##############################
$1 clearpath 'source $HOME/Dropbox/scripts/dotfiles/bin/clearpath.sh'
$1 ra '$HOME/Dropbox/scripts/dotfiles/generalrc/generalrc_alias_renew.sh'
$1 la 'ls -la'
$1 ll 'ls -la'
$1 sl 'ls'
#$1 rm 'rm -rf'   # too general and maybe too dangerous?
$1 vi 'vim -v'
$1 vim 'vim -v'

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

$1 extract '$HOME/Dropbox/scripts/dotfiles/bin/extract.sh'

$1 weekly_update 'cd $HOME/Dropbox/scripts/dotfiles/cron;./weeklyupdate.sh'


$1 du 'echo ______better ncdu______;echo; du'
$1 lstarfile 'tar -ztvf'
$1 diff 'colordiff'
$1 ep 'epstopdf'
$1 grep 'grep --color=auto -d skip'    # skips directories
#$1 grep'echo use ack instead'
$1 bc 'bc -l'
#alias cp='cp -a'  # this one did not copy symliks as files but as symlinks -> remove -a
$1 xmgrace 'xmgrace -geom 1100x860 -nosigcatch -param ~/.xmgracerc'
$1 x 'xmgrace -geom 1100x860 -nosigcatch -param ~/.xmgracerc'
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
#$1 load_last_saved_settings='$HOME/Dropbox/scripts/dotfiles/mac_KeyRemap4MacBook/load_last_saved_settings.sh'
$1 cpdf '$HOME/Dropbox/scripts/dotfiles/bin/cpdf-binaries/OSX-Intel/cpdf'
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
    $1 units '$HOME/scripts/dotfiles/bin/units/units-1.88/units' && \
    $1 top 'top -o cpu' && \
    $1 trash "rmtrash" && \
    $1 del "rmtrash" && \
    $1 rm 'echo Use del command, or the full path i.e. /bin/rm' && \
    $1 ctagsnew '/usr/local/Cellar/ctags/5.8/bin/ctags -R .' && \
    $1 mvim     'open /Applications/MacVim.app'


    #alias sed='$HOME/scripts/dotfiles/bin/sed/sed-4.2/build/bin/sed-4.2'
    #alias mvim='open /Applications/MacVim.app'
    #alias cat 'colorize'  # colorize seems to have sometimes problems


[ "$oncmmd" != "true" ] && \
    $1 qhost 'qhost.sh' && \
    $1 qstat 'qstat.sh' && \
    $1 qdel 'qdel.sh'   && \
    $1 qalter 'qalter.sh' && \
    $1 qsub 'qsub.sh'   && \
    $1 qls 'qls.alexej.cmmc.sh -u'  && \
    $1 qhold 'qhold.sh' && \
    $1 qrls 'qrls.sh'         

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
    $1 vi 'gvim -v' && \
    $1 vim 'gvim -v' && \
    $1 tmux '/u/aglen/local/bin/tmux'

[ "$oncmpc" = "true" ] && \
    $1 vi 'gvim -v' && \
    $1 vim 'gvim -v' && \
    $1 tmux '/home/glensk/local/bin/tmux'

[ "$host" = "cmmd001.mpie.de" ] && $1 vi 'vi -v' && $1 vim 'vi -v'
#[ "$host" = "cmmc001.mpie.de" ] && $1 vi 'vi -v' && $1 vim 'vi -v' # check if cmmc has vim or only vi

##############################
# ssh
##############################
$1 ssh 'ssh.sh -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200000'   # -X for loading locas stuff locally (xmgrace)  ServerAliveInterval=1600 is now enough for cmpc
$1 sshimmer 'ssh.sh -Y -X -o ServerAliveInterval=160 -o ServerAliveCountMax=1200000'

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
