#!/bin/sh

tac $1 | grep -m 1 "Final energy" | awk '{print $4}'

