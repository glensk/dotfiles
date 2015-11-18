#!/bin/sh

[ "$1" = "" ] && echo 'provide $1 to be the inputfile' && exit
[ ! -f "$1" ] && echo file $1 not found && exit


pdftk $1 burst


rm doc_data.txt
