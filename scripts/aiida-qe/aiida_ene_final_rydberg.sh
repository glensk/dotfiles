#!/bin/sh
tac aiida.out | grep -m 1 "Final energy" | awk '{print $4}'

