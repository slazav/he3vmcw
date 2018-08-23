#!/usr/bin/gnuplot

set terminal "x11"
set style data lines


pref="hpd"
num="0"

plot [] \
  pref.".magn.dat" using 1:3 with lines lw 2 title "Mx",\
  pref.".magn.dat" using 1:4 with lines lw 2 title "My",\
  pref.".magn.dat" using 1:5 with lines lw 2 title "Mz",\
0 lc 2 notitle,\
1 lc 3 notitle,\
-1/4.0 lc 3 notitle

pause -1

#set terminal "fig" metric size 16 6
#set output pref.".magn.fig"
#replot