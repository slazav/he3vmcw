#!/usr/bin/gnuplot

#set terminal "x11"
set style data lines


pref="deform/hpd_th2a"
num="6"


set arrow 2 from -1,0 to 1,0 nohead
set arrow 1 from 0,-1 to 0,1 nohead
set object 2 ellipse at 0,0 size 2*acos(-0.25)/pi, 2*acos(-0.25)/pi
set object 1 ellipse at 0,0 size 2,2

plot [-2:2] [-2:2]\
  pref.".prof".num.".dat" using (sqrt($5**2+$6**2)/pi):($7/pi) with linespoints pt 6 ps 0.5 lw 2 title "(n_p vs nz)*theta",\
  pref.".prof".num.".dat" using (sqrt($2**2+$3**2)):4 with linespoints pt 6 ps 0.5 lw 2 title "m_p vs mz"

#  pref.".prof".num.".dat" using 2:4 with linespoints pt 6 ps 0.5 lw 2 title "(mx-mz)",\
#  pref.".prof".num.".dat" using 3:4 with linespoints pt 6 ps 0.5 lw 2 title "(my-mz)"
#  pref.".prof".num.".dat" using ($6*$8/pi):($7*$8/pi) with linespoints pt 6 ps 0.5 lw 2 title "(ny vs nz)*theta",\
#  pref.".prof".num.".dat" using ($5*$8/pi):($7*$8/pi) with linespoints pt 6 ps 0.5 lw 2 title "(nx vs nz)*theta",\

pause -1

#set terminal "fig" metric size 16 6
#set output pref.".prof".num.".fig"
#replot