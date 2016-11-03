#!/bin/sh

[ ! -e "$1" ] && echo 'please specify the file by as first argument $1' && exit -1
file=$1; string=$(wc -c $file); bite=${string% *}; okay=$(echo "scale=2; $bite/1024" | bc);friend=$(echo -e "$file $okay" "kb"); echo -e "$friend"
