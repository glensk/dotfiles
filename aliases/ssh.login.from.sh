#!/bin/bash

tty=`tty | sed 's|/dev/||'`
echo tty:$tty:
from=`finger $user | grep $tty | awk '{print $NF}' | sed 's/\..*$//' | sed 's/(//g'`
echo from:$from:
