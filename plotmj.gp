#!/usr/bin/gnuplot

set terminal "x11"
set style data lines
set nokey

#set terminal "fig" metric size 16 6
#set output "out_m.fig"

set style line 1 lc 1
set style line 2 lc 0
set style line 3 lc 3

set title "Mx, My, Mz" 

cell_len=0.18

plot [] [-0.4:1.1]\
1 ls 2,\
  "mj0.0000.dat" using 1:3 ls 1 lw 1.5,\
  "mj0.0181.dat" using 1:3 ls 1,\
  "mj0.0360.dat" using 1:3 ls 1,\
  "mj0.0541.dat" using 1:3 ls 1,\
  "mj0.0721.dat" using 1:3 ls 1,\
  "mj0.0900.dat" using 1:3 ls 1,\
  "mj0.1081.dat" using 1:3 ls 1,\
  "mj0.1260.dat" using 1:3 ls 1,\
  "mj0.1441.dat" using 1:3 ls 1,\
  "mj0.1621.dat" using 1:3 ls 1,\
  "mj0.1800.dat" using 1:3 ls 1,\
sqrt(15/16.0) ls 2,\
  "mj0.0000.dat" using 1:4 ls 1 lw 1.5,\
  "mj0.0181.dat" using 1:4 ls 1,\
  "mj0.0360.dat" using 1:4 ls 1,\
  "mj0.0541.dat" using 1:4 ls 1,\
  "mj0.0721.dat" using 1:4 ls 1,\
  "mj0.0900.dat" using 1:4 ls 1,\
  "mj0.1081.dat" using 1:4 ls 1,\
  "mj0.1260.dat" using 1:4 ls 1,\
  "mj0.1441.dat" using 1:4 ls 1,\
  "mj0.1621.dat" using 1:4 ls 1,\
  "mj0.1800.dat" using 1:4 ls 1,\
0 ls 2,\
  "mj0.0000.dat" using 1:5 ls 1 lw 1.5,\
  "mj0.0181.dat" using 1:5 ls 1,\
  "mj0.0360.dat" using 1:5 ls 1,\
  "mj0.0541.dat" using 1:5 ls 1,\
  "mj0.0721.dat" using 1:5 ls 1,\
  "mj0.0900.dat" using 1:5 ls 1,\
  "mj0.1081.dat" using 1:5 ls 1,\
  "mj0.1260.dat" using 1:5 ls 1,\
  "mj0.1441.dat" using 1:5 ls 1,\
  "mj0.1621.dat" using 1:5 ls 1,\
  "mj0.1800.dat" using 1:5 ls 1,\
-1/4. ls 2,\
0 lc 0

pause -1

#set output "out_n.fig"

set title "Nx, Ny, Nz" 

plot [] [-0.2:1.1]\
  "mj0.0000.dat" using 1:6 ls 1 lw 1.5,\
  "mj0.0181.dat" using 1:6 ls 1,\
  "mj0.0360.dat" using 1:6 ls 1,\
  "mj0.0541.dat" using 1:6 ls 1,\
  "mj0.0721.dat" using 1:6 ls 1,\
  "mj0.0900.dat" using 1:6 ls 1,\
  "mj0.1081.dat" using 1:6 ls 1,\
  "mj0.1260.dat" using 1:6 ls 1,\
  "mj0.1441.dat" using 1:6 ls 1,\
  "mj0.1621.dat" using 1:6 ls 1,\
  "mj0.1800.dat" using 1:6 ls 1,\
0 ls 2,\
  "mj0.0000.dat" using 1:7 ls 1 lw 1.5,\
  "mj0.0181.dat" using 1:7 ls 1,\
  "mj0.0360.dat" using 1:7 ls 1,\
  "mj0.0541.dat" using 1:7 ls 1,\
  "mj0.0721.dat" using 1:7 ls 1,\
  "mj0.0900.dat" using 1:7 ls 1,\
  "mj0.1081.dat" using 1:7 ls 1,\
  "mj0.1260.dat" using 1:7 ls 1,\
  "mj0.1441.dat" using 1:7 ls 1,\
  "mj0.1621.dat" using 1:7 ls 1,\
  "mj0.1800.dat" using 1:7 ls 1,\
1 ls 2,\
  "mj0.0000.dat" using 1:8 ls 1 lw 1.5,\
  "mj0.0181.dat" using 1:8 ls 1,\
  "mj0.0360.dat" using 1:8 ls 1,\
  "mj0.0541.dat" using 1:8 ls 1,\
  "mj0.0721.dat" using 1:8 ls 1,\
  "mj0.0900.dat" using 1:8 ls 1,\
  "mj0.1081.dat" using 1:8 ls 1,\
  "mj0.1260.dat" using 1:8 ls 1,\
  "mj0.1441.dat" using 1:8 ls 1,\
  "mj0.1621.dat" using 1:8 ls 1,\
  "mj0.1800.dat" using 1:8 ls 1,\
0. ls 2,\
0 lc 0

pause -1

#set output "out_t.fig"

set title "Theta" 

plot []\
  "mj0.0000.dat" using 1:(cos($9)) ls 1 lw 1.5,\
  "mj0.0181.dat" using 1:(cos($9)) ls 1,\
  "mj0.0360.dat" using 1:(cos($9)) ls 1,\
  "mj0.0541.dat" using 1:(cos($9)) ls 1,\
  "mj0.0721.dat" using 1:(cos($9)) ls 1,\
  "mj0.0900.dat" using 1:(cos($9)) ls 1,\
  "mj0.1081.dat" using 1:(cos($9)) ls 1,\
  "mj0.1260.dat" using 1:(cos($9)) ls 1,\
  "mj0.1441.dat" using 1:(cos($9)) ls 1,\
  "mj0.1621.dat" using 1:(cos($9)) ls 1,\
  "mj0.1800.dat" using 1:(cos($9)) ls 1,\
-1/4. lc 0

pause -1

