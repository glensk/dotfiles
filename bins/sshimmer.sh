ssh.sh -Y -X -o ServerAliveInterval=1600 -o ServerAliveCountMax=1200 -t $* "[ -e `th.sh` ] && cd `th.sh`; zsh"

