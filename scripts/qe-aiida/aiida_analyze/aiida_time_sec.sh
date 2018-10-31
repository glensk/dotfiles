#!/bin/sh
tac aiida.out | grep -m 1 time | awk '{print $9}'

