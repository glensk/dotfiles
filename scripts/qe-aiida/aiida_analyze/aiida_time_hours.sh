#!/bin/sh
tac $1 | grep -m 1 time | awk '{print $9/3600}'

