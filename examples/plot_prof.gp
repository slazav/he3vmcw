#!/usr/bin/gnuplot

set terminal "x11"
set style data lines


pref="hpd"
num="0"

plot [] \
  pref.".prof".num.".dat" using 1:2 with lines lw 2 title "Mx",\
  pref.".prof".num.".dat" using 1:3 with lines lw 2 title "My",\
  pref.".prof".num.".dat" using 1:4 with lines lw 2 title "Mz",\
  pref.".prof".num.".dat" using 1:5 with lines lw 2 title "nx",\
  pref.".prof".num.".dat" using 1:6 with lines lw 2 title "ny",\
  pref.".prof".num.".dat" using 1:7 with lines lw 2 title "nz",\
  pref.".prof".num.".dat" using 1:8 with lines lw 2 title "theta",\
0 lc 1 notitle,\

pause -1

#set terminal "fig" metric size 16 6
#set output pref.".prof".num.".fig"
#replot