###################################################################################
#### XXX ### system specific 
###################################################################################
set     green='{\033[1;32m%}'
set   setback='{\033[0m%}'
[ "$myprompthostuser" = "black" ]       && set myprompthostuser='{\033[1;30m%}'
[ "$myprompthostuser" = "red" ]         && set myprompthostuser='{\033[1;31m%}'
[ "$myprompthostuser" = "green" ]       && set myprompthostuser='{\033[1;32m%}'
[ "$myprompthostuser" = "orange" ]      && set myprompthostuser='{\033[1;33m%}'
[ "$myprompthostuser" = "blue" ]        && set myprompthostuser='{\033[1;34m%}'
[ "$myprompthostuser" = "magenta" ]     && set myprompthostuser='{\033[1;35m%}'
[ "$myprompthostuser" = "turquoise" ]   && set myprompthostuser='{\033[1;36m%}'
[ "$myprompthostuser" = "white" ]       && set myprompthostuser='{\033[1;37m%}'
set prompt='\n%'"$green"'%P%'"$myprompthostuser"' %B%n%b%'"$myprompthostuser"'@%m %B%'"$myprompthostuser"'%/%b%u \n%'"$setback"'%#'
#@ #             %T  Time            %n glensk   %m mac %/ gibt den pfad an
#@ # %T        Time (no seconds)
#@ # %P        Time (including seconds)
#@ # %B (%b)   Start (stop) boldfacing mode.
#@ # %n        $USER 
#@ # %U (%u)   Start (stop) underline mode.
#@ # %m        hostname
#@ # %/        current path
