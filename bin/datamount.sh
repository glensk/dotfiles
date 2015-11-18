#!/bin/sh

base=`basename $0 | sed 's|mount.sh||'`
volume="/$base/$USER"

mount_immer.sh $base $myhost
