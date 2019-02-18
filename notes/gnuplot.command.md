to test if gnuplot works:
    type in terminal: gnuplot
    then: plot sin(x)
    if this is not showing a graph gnuplot has to be installed differently or you have to login to ssh with -X

in terminal type: gnuplot
    plot 'file.dat' using 1:2

in terminal type: gnuplot
    plot 'file.dat' 1:2

*********************************#gnuplot**************************
plot 'plot.dat' 1:2 with points           means: x=1:y=2
f(x)=a*x**2+b*x+c                        set fit function f(x)=a*x²+b*x+c
fit f(x) 'plot.dat' via a, b, c
plot 'plot.dat' with points, f(x)            ohne 1:2
    set term png (or jpg, eps, ...)
    set out "test.png"
    plot 'plot.dat' with points, f(x)
