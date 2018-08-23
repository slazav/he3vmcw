#!/usr/bin/gnuplot

set terminal "x11"
set style data lines
set nokey


pref="hpd"
num="0"

plot [] \
  pref.".mesh".num.".dat" using 2:1 with linespoints pt 7 ps 0.8,\

pause -1

#set terminal "fig" metric size 16 6
#set output pref.".mesh".num.".fig"
#replot