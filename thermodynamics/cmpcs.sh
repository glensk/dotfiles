#!/bin/bash

rm -f cmpc_list; touch cmpc_list

for (( i=1; i<=99; i++ )); do
  s=`echo $i | awk '{printf("%02i",$1)}'`
  ping -c1 cmpc$s 2> _tmp_err > _tmp_cmpc
  c1=`grep "unknown host" _tmp_err`
  c2=`grep "0 received" _tmp_cmpc`
  if [ -z "$c1" -a -z "$c2" -a -z "$c3" ]; then 
    who=`ssh cmpc$s who 2> _tmp_err | awk '{print $1}' | sort -u | xargs`
    c=`grep "Connection refused" _tmp_err`
    if [ -z "$c" ]; then
      last=`ssh cmpc$s last 2> /dev/null | awk '{print $1}' | sort -u | grep -v grabowsk | xargs`
      echo -e "cmpc$s  \033[31m\033[1m$who\033[0m   $last" >> cmpc_list
      tail -n1 cmpc_list; rm _tmp_err _tmp_cmpc
    fi
  fi
done
