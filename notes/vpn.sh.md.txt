#!/bin/sh

idline=`ps aux | grep "sbin/racoon" | grep -v "grep "`
echo "idline: $idline"
id=`echo "$idline" | awk '{print $2}'`
echo "id: $id"
[ "$id" != "" ] && echo "sudo kill -9 $id" && sudo kill -9 $id
[ "$id" == "" ] && echo "starting /usr/sbin/racoon ... now vpn should work again :)" && sudo /usr/sbin/racoon
