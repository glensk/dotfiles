#!/bin/sh
cd cd ~/proj/0000_Dienstreise_Urlaub_Stundenzettel/Stundenzettel_novaTime
cd ~/Thermodynamics
echo "5630" | ./novaTime.sh -u Glensk -p 70407 -i
open MPIE-Timesheet.xlsx

