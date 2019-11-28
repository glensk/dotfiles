set label "G" at 0.0,-0.3
set label "X" at 0.5,-0.3
set label "M" at 1.0,-0.3
set label "R" at 1.5,-0.3
set label "G" at 2.366,-0.3

set arrow from 0.5,0.0 to 0.5,10.0 nohead lw 2
set arrow from 1.5,0.0 to 1.5,10.0 nohead lw 2
set arrow from 1.0,0.0 to 1.0,10.0 nohead lw 2

set ylabel "Energy (eV)"
unset xtics
set arrow from 0.0,5.5053 to 2.366,5.5053 nohead lw 2
plot [0:2.366][0:10] 'bands.dat' u 1:2 w l lw 2 notitle,  'bands.dat' u 1:3 w l lw 2 notitle, 'bands.dat' u 1:4 w l lw 2 notitle, 'bands.dat' u 1:5 w l lw 2 notitle, 'bands.dat' u 1:6 w l lw 2 notitle, 'bands.dat' u 1:7 w l lw 2 notitle, 'bands.dat' u 1:8 w l lw 2 notitle, 'bands.dat' u 1:9 w l lw 2 notitle, 'bands.dat' u 1:10 w l lw 2 notitle, 'bands.dat' u 1:11 w l lw 2 notitle, 'bands.dat' u 1:12 w l lw 2 notitle, 'bands.dat' u 1:13 w l lw 2 notitle, 'bands.dat' u 1:14 w l lw 2 notitle, 'bands.dat' u 1:15 w l lw 2 notitle, 'bands.dat' u 1:16 w l lw 2 notitle, 'bands.dat' u 1:17 w l lw 2 notitle, 'bands.dat' u 1:18 w l lw 2 notitle, 'bands.dat' u 1:19 w l lw 2 notitle
