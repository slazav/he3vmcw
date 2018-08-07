#!/usr/bin/gnuplot

set terminal "x11"
set style data lines
set nokey

#set terminal "fig" metric size 16 6
#set output "out_m.fig"

plot [] \
  "mesh.txt" using 2:1 with linespoints pt 6,\
  "mesh.txt" using 2:($3*300) with linespoints pt 6,\
  "mesh.txt" using 2:($4) with linespoints pt 6,\

pause -1
