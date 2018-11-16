#!/usr/bin/gnuplot

set terminal "x11"
set style data lines


pref="deform/hpd_th2a"
num="2"

plot [] \
  pref.".prof".num.".dat" using 0:2 with lines lw 2 title "Mx",\
  pref.".prof".num.".dat" using 0:3 with lines lw 2 title "My",\
  pref.".prof".num.".dat" using 0:4 with lines lw 2 title "Mz",\
  pref.".prof".num.".dat" using 0:5 with lines lw 2 title "nx",\
  pref.".prof".num.".dat" using 0:6 with lines lw 2 title "ny",\
  pref.".prof".num.".dat" using 0:7 with lines lw 2 title "nz",\
  pref.".prof".num.".dat" using 0:(sqrt($5*$5+$6*$6+$7*$7)) with lines lw 2 title "nz",\
0 lc 1 notitle,\

pause -1

#set terminal "fig" metric size 16 6
#set output pref.".prof".num.".fig"
#replot