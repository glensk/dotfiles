#!/bin/sh
#export GDFONTPATH="/usr/share/fonts/truetype/msttcorefonts/"
gnuplot << EOF
set term png font "arialbd.ttf" 18 size 950,600
set output "Calibration Curve - Full.png"
set title "Calibration Curve - Full"
set key noautotitles
unset mouse
set bmargin 4
set grid xtics ytics
set xlabel "10^5/ADC"
set format x "%3.0f"
set ylabel "Resistance - Ohm"
set format y "%3.0f"
set yrange [0:100]
set datafile separator "\t"
f(x) = m*x + c
fit f(x) "Measurements/Calibration.csv" using 3:1 via m,c
set label 1 sprintf("m = %3.4f",m) at 510,75 font "arialbd,18"
set label 2 sprintf("c = %3.4f",c) at 510,70 font "arialbd,18"
plot    \
 "Measurements/Calibration.csv" \
 using 3:1 with linespoints lt 3 lw 3 pt 3 , \
f(x) lt 4 lw 2 
EOF

