###################################################################################
#### XXX ### system specific 
###################################################################################
# Set the PS1 prompt (with colors).
# Based on http://www-128.ibm.com/developerworks/linux/library/l-tip-prompt/
# And http://networking.ringofsaturn.com/Unix/Bash-prompts.php .
#\u - Username. The original prompt also has 
#\h - host name.
#\H - Hostname
#\w - Current absolute path. Use \W for current relative path.
#\$ - The prompt character (eg. # for root, $ for regular users).
#\[ and \] - These tags should be placed around color codes so bash knows how to properly place the cursor.
#   \a		an ASCII bell character (07)
#	\d		the date in "Weekday Month Date" format (e.g., "Tue May 26")
#	\D{format}	the format is passed to strftime(3) and the result
#			  is inserted into the prompt string an empty format
#			  results in a locale-specific time representation.
#			  The braces are required
#	\e		an ASCII escape character (033)
#	\h		the hostname up to the first `.'
#	\H		the hostname
#	\j		the number of jobs currently managed by the shell
#	\l		the basename of the shell's terminal device name
#	\n		newline
#	\r		carriage return
#	\s		the name of the shell, the basename of $0 (the portion following
#			  the final slash)
#	\t		the current time in 24-hour HH:MM:SS format
#	\T		the current time in 12-hour HH:MM:SS format
#	\@		the current time in 12-hour am/pm format
#	\A		the current time in 24-hour HH:MM format
#	\u		the username of the current user
#	\v		the version of bash (e.g., 2.00)
#	\V		the release of bash, version + patch level (e.g., 2.00.0)
#	\w		the current working directory, with $HOME abbreviated with a tilde
#	\W		the basename of the current working directory, with $HOME
#			 abbreviated with a tilde
#	\!		the history number of this command
#	\#		the command number of this command
#	\$		if the effective UID is 0, a #, otherwise a $
#	\nnn		the character corresponding to the octal number nnn
#	\\		a backslash
#	\[		begin a sequence of non-printing characters, which could be used
#			  to embed a terminal control sequence into the prompt
#	\]		end a sequence of non-printing characters
#
# Reset
Color_Off='\e[0m'       # Text Reset

# Regular Colors
Black='\e[0;30m'        # Black
Red='\e[0;31m'          # Red
Green='\e[0;32m'        # Green
Yellow='\e[0;33m'       # Yellow
Blue='\e[0;34m'         # Blue
Purple='\e[0;35m'       # Purple
Cyan='\e[0;36m'         # Cyan
White='\e[0;37m'        # White

set     green='{\033[1;32m%}'
set   setback='{\033[0m%}'
[ "$myprompthostuser" = "black" ]       && myprompthostuser='\e[0;30m'
[ "$myprompthostuser" = "red" ]         && myprompthostuser='\e[0;31m'
[ "$myprompthostuser" = "green" ]       && myprompthostuser='\e[0;32m'
[ "$myprompthostuser" = "orange" ]      && myprompthostuser='\e[0;33m'
[ "$myprompthostuser" = "yellow" ]      && myprompthostuser='\e[0;33m'
[ "$myprompthostuser" = "blue" ]        && myprompthostuser='\e[0;34m'
[ "$myprompthostuser" = "magenta" ]     && myprompthostuser='\e[0;35m'
[ "$myprompthostuser" = "turquoise" ]   && myprompthostuser='\e[0;36m'
[ "$myprompthostuser" = "white" ]       && myprompthostuser='\e[0;37m'
set myprompthostuser=33

#PS1="[\j] \A \e[04;34m\u\e[m@\e[04;33m\H $PWD\e[m\n$"
#PS1="$Yellow[\j] \t $Red\u@\H $PWD\e[m\n$"
#PS1="$Yellow\t $myprompthostuser\u@\H $PWD\e[m\n$"
#PS1="$Yellow\t $myprompthostuser\u@\H $PWD\e[m\n$"
#PS1="${RESET}${YELLOW}\u@\h${NORMAL} \`${SELECT}\` ${YELLOW}\w \$(__git_ps1) >${NORMAL} "
#PS1="${RESET}${YELLOW}\t \u@\H \`${SELECT}\` ${YELLOW}PWD >${NORMAL} "
PS1='\[`[ $? = 0 ] && X=2 || X=1; tput setaf $X`\]\t \u@\H\[`tput sgr0`\] $PWD\n\$ '
#\H Host

