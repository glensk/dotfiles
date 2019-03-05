tcsh -c 'ls'  ## fuehrt das in der tcsh aus
bash -c 'ls'  ## in bash


[ "$host" = "$mylaptop" ] && bash -c 'echo -ne "\033]6;1;bg;green;brightness;255\a"'  # sets default tab color to green, for other colors see ssh.sh
