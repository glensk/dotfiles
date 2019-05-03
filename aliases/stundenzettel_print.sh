#!/bin/sh
#cd ~/proj/0000_Dienstreise_Urlaub_Stundenzettel/Stundenzettel_novaTime
cd ~/Thermodynamics
echo "5630" | ./novaTime.sh -u Glensk -p 70407 -i   # Personalnummer 5630
open MPIE-Timesheet.xlsx

