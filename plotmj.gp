#!/usr/bin/gnuplot

set terminal "x11"
set style data lines
set nokey

#set terminal "fig" metric size 16 6
#set output "out_m.fig"

set style line 1 lc 1
set style line 2 lc 2
set style line 3 lc 3

set title "Mx, My, Mz" 

cell_len=0.18

plot [] [-1:1.1]\
  "mj0.0000.dat" using 1:3with lines lc 1 lw 1.5,\
  "mj0.0090.dat" using 1:3with lines lc 1,\
  "mj0.0181.dat" using 1:3with lines lc 1,\
  "mj0.0271.dat" using 1:3with lines lc 1,\
  "mj0.0362.dat" using 1:3with lines lc 1,\
  "mj0.0450.dat" using 1:3with lines lc 1,\
  "mj0.0561.dat" using 1:3with lines lc 1,\
  "mj0.0651.dat" using 1:3with lines lc 1,\
  "mj0.0742.dat" using 1:3with lines lc 1,\
  "mj0.0832.dat" using 1:3with lines lc 1,\
  "mj0.0900.dat" using 1:3with lines lc 1,\
sqrt(15/16.0)with lines lc 1,\
  "mj0.0000.dat" using 1:4with lines lc 2 lw 1.5,\
  "mj0.0090.dat" using 1:4with lines lc 2,\
  "mj0.0181.dat" using 1:4with lines lc 2,\
  "mj0.0271.dat" using 1:4with lines lc 2,\
  "mj0.0362.dat" using 1:4with lines lc 2,\
  "mj0.0450.dat" using 1:4with lines lc 2,\
  "mj0.0561.dat" using 1:4with lines lc 2,\
  "mj0.0651.dat" using 1:4with lines lc 2,\
  "mj0.0742.dat" using 1:4with lines lc 2,\
  "mj0.0832.dat" using 1:4with lines lc 2,\
  "mj0.0900.dat" using 1:4with lines lc 2,\
0with lines lc 2,\
  "mj0.0000.dat" using 1:5with lines lc 3 lw 1.5,\
  "mj0.0090.dat" using 1:5with lines lc 3,\
  "mj0.0181.dat" using 1:5with lines lc 3,\
  "mj0.0271.dat" using 1:5with lines lc 3,\
  "mj0.0362.dat" using 1:5with lines lc 3,\
  "mj0.0450.dat" using 1:5with lines lc 3,\
  "mj0.0561.dat" using 1:5with lines lc 3,\
  "mj0.0651.dat" using 1:5with lines lc 3,\
  "mj0.0742.dat" using 1:5with lines lc 3,\
  "mj0.0832.dat" using 1:5with lines lc 3,\
  "mj0.0900.dat" using 1:5with lines lc 3,\
1with lines lc 3,\
-1/4.with lines lc 3,\
0 lc 0

pause -1

set terminal png
set output "m1.png"
replot

set terminal x11
set title "Nx, Ny, Nz" 

plot [] [-1:1.1]\
  "mj0.0000.dat" using 1:6with lines lc 1 lw 1.5,\
  "mj0.0090.dat" using 1:6with lines lc 1,\
  "mj0.0181.dat" using 1:6with lines lc 1,\
  "mj0.0271.dat" using 1:6with lines lc 1,\
  "mj0.0362.dat" using 1:6with lines lc 1,\
  "mj0.0450.dat" using 1:6with lines lc 1,\
  "mj0.0561.dat" using 1:6with lines lc 1,\
  "mj0.0651.dat" using 1:6with lines lc 1,\
  "mj0.0742.dat" using 1:6with lines lc 1,\
  "mj0.0832.dat" using 1:6with lines lc 1,\
  "mj0.0900.dat" using 1:6with lines lc 1,\
  "mj0.0000.dat" using 1:7with lines lc 2 lw 1.5,\
  "mj0.0090.dat" using 1:7with lines lc 2,\
  "mj0.0181.dat" using 1:7with lines lc 2,\
  "mj0.0271.dat" using 1:7with lines lc 2,\
  "mj0.0362.dat" using 1:7with lines lc 2,\
  "mj0.0450.dat" using 1:7with lines lc 2,\
  "mj0.0561.dat" using 1:7with lines lc 2,\
  "mj0.0651.dat" using 1:7with lines lc 2,\
  "mj0.0742.dat" using 1:7with lines lc 2,\
  "mj0.0832.dat" using 1:7with lines lc 2,\
  "mj0.0900.dat" using 1:7with lines lc 2,\
  "mj0.0000.dat" using 1:8with lines lc 3 lw 1.5,\
  "mj0.0090.dat" using 1:8with lines lc 3,\
  "mj0.0181.dat" using 1:8with lines lc 3,\
  "mj0.0271.dat" using 1:8with lines lc 3,\
  "mj0.0362.dat" using 1:8with lines lc 3,\
  "mj0.0450.dat" using 1:8with lines lc 3,\
  "mj0.0561.dat" using 1:8with lines lc 3,\
  "mj0.0651.dat" using 1:8with lines lc 3,\
  "mj0.0742.dat" using 1:8with lines lc 3,\
  "mj0.0832.dat" using 1:8with lines lc 3,\
  "mj0.0900.dat" using 1:8with lines lc 3,\
1with lines lc 3,\
0with lines lc 2

pause -1

set terminal png
set output "n1.png"
replot


set terminal x11
set title "cos(theta)" 

plot [] []\
  "mj0.0000.dat" using 1:(cos($9))with lines lc 1 lw 1.5,\
  "mj0.0090.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0181.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0271.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0362.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0450.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0561.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0651.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0742.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0832.dat" using 1:(cos($9))with lines lc 1,\
  "mj0.0900.dat" using 1:(cos($9))with lines lc 1,\
-1/4.0with lines lc 2

pause -1

set terminal png
set output "t1.png"
replot
