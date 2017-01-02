(qstat -j 12345) 2>&1  kann gespeichert werden ohne () nicht

echo y | read y    # answers directly
$ echo "cherry apple peach" | tr " " "\n"
cherry
apple
peach 

$ COL_BLUE="\x1b[34;01m"
$ COL_RESET="\x1b[39;49;00m"
$ echo -e $COL_BLUE"Important Message: "$COL_RESET"This is a message"


echo "..." 1>&2         ## displays only on screen, cant be redirected/find in files
echo "..." > /dev/null  ## do not print to screen



echo -en "\r text" ersetzt diese zeile danach

schreibe rot
  echo -e "\033[31m\033[1mERROR\033[0m: $math not available" 1>&2
\033[31m\033[1m   ROT           \033[0m
\033[32m\033[1m   GRUEN         \033[0m
\033[33m\033[1m   GELB          \033[0m
\033[34m\033[1m   BLAU          \033[0m
\033[35m\033[1m   ROSA          \033[0m
\033[38m\033[1m   Schwarz Bold  \033[0m



  ERROR\033[0m: $math not available" 1>&2


both on mac and cmpc(fedora):
start script with #!/bin/bash
echo -e "\n\033[1;32m   only bash and -e: Light Colors\033[0m  \t\t\033[1;4;31m  Dark Colors\033[0m"
echo -e "\033[1;4;31m  RED \033[0m"
echo -e "\033[1;4;32m  GREEN \033[0m"

echo "\033[1;4;31m  RED \033[0m"
echo "\033[1;4;32m  GREEN \033[0m"
