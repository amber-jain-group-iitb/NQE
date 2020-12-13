set terminal png size
set output 'pop.png'
set title 'ground state population'

set ylabel 'population'
set xlabel 'time(ps)'
set yrange[0:1.1]
#set xrange[0:2001]

plot 'pops.dat' using ($1/40000):2 w lines
