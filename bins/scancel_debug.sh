#!/bin/sh
DEBUG_JOBID=$(squeue -u $USER -p debug | tail -n 1 | awk '{print $1}');
scancel $DEBUG_JOBID

