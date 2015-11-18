#######################################################################
# PROMPT: set only after sourcing oh-my-zsh
#######################################################################
#echo yo $myprompttime
#PROMPT="%{$fg[red]%}%n%{$reset_color%}@%{$fg[blue]%}%m %{$fg_no_bold[yellow]%}%1~ %{$reset_color%}%#"
#RPROMPT="[%{$fg_no_bold[yellow]%}%?%{$reset_color%}]"
local prompttime="%{$fg[$myprompttime]%}%*%{$reset_color%}"
local user_host="%{$fg[$myprompthostuser]%}%n@%m%{$reset_color%}"
local current_dir="%{$fg[$mypromptpath]%}%d%{$reset_color%}"
PROMPT="${prompttime} ${user_host} ${current_dir} ${svn_branch}
%#"

