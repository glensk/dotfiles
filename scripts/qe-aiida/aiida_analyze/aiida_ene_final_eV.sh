#!/bin/sh

tac $1 | grep -m 1 "Final energy" | awk '{printf "%.10f\n", $4*13.605692}'

