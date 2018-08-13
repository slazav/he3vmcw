#!/usr/bin/perl
use strict;
use warnings;

# create vmcw_pars.fh from vmcw_pars.h

# parameters (type, name, default, comment):
my $pars ="
double time      0.0     current time, s (always 0 in the begining)
double tend      0.0     end of current sweep/wait commant, s (always 0 in the begining)
double tstep     5e-3    time step, s (can be changed at any time)
double t11       0.0
double H         284.0   constant part of magnetic field, G
double grad      0.1     gradient of magnetic field, G/cm
double HR0       1e-3    RF field, G
double HR_SWR    0.0     RF field sweep rate, G/s
double LP0       0.0     Larmor position, cm
double LP_SWR    0.0     Larmor position sweep rate, cm/s
double DF0       0.1     Spin diffusion, cm^2/s
double DF_SWR    0.0     Spin diffusion sweep rate (cm^2/s)/s
double TF0       1e-5    Leggett-Takagi relaxation time, s
double TF_SWR    0.0     Leggett-Takagi relaxation time sweep rate s/s
double LF0       1e5     Leggett frequency, Hz
double LF_SWR    0.0     Leggett frequency sweep rate, Hz/s
double CPAR0     5e2     spin-wave velocity (c_par?), cm/s
double CPAR_SWR  0.0     spin-wave velocity sweep rate, (cm/s)/s
double AER_LEN   0.5     aerogel length / cell length
double AER_CNT   0.0     center of aerogel area / cell length
double AER_TRW   6e-3    transition width / cell length
double CELL_LEN  0.18    cell length, cm
double XMESH_K   0.2     mesh step: CELL_LEN/(NPTS-1) / (1+XMESH_K*AER_STEP)
double XMESH_ACC 1e-10   mesh accuracy
double BETA      0.0     initial tipping angle
double PRESS     0.0     pressure, bar
double TTC       0.0     temperature T/Tc
double TTC_ST    0.0     temperature step
double T1C       6e14    T1=T1C*T
int    AER       0       use aerogel
int    IBN       2       type of boundary condition: 1 - open cell, 2 - closed cell
";

# type conversion C->F
my %types = (
  'double' => 'real*8',
  'int'    => 'integer',
);

open CH, "> vmcw_pars.h"   or die "can't open vmcw_pars.h: $!\n";
open FH, "> vmcw_pars.fh"  or die "can't open vmcw_pars.fh: $!\n";
open CC, "> vmcw_pars.cpp" or die "can't open vmcw_pars.cpp: $!\n";

######
print CH qq*// This file is created by $0 script. Do not modify.
#ifndef VMCW_PARS_H
#define VMCW_PARS_H
void set_def_pars(struct pars_t \*p);
extern "C" {
  struct pars_t {
*;

######
print FH qq*! This file is created by $0 script. Do not modify.

        implicit real\*8(A-H,O-Z)
*;

######
print CC qq*// This file is created by $0 script. Do not modify.
#include "vmcw_pars.h"
void
set_def_pars(struct pars_t \*p) {
*;

my @names;
foreach (split /\n/, $pars) {
  next unless /(\S+)\s+(\S+)\s+(\S+)(\s+(.*))?/;
  my ($ctype, $name, $def, $comm) = ($1,$2,$3,$5);
  next unless $ctype && $name;
  my $ftype = $types{$ctype};
  next unless $ftype;

  my $ccomm = $comm? " ///< $comm": '';
  my $fcomm = $comm? " ! $comm": '';
  printf FH "      %8s %12s$fcomm\n", $ftype, $name;
  printf CH "    %8s %12s;$ccomm\n", $ctype, $name;
  printf CC "  p->%-12s = $def;\n", $name;
  push @names, $name;
}

print FH "\n", ' 'x6, "common /PARS/";
my $l=0;
for (my $i=0; $i<=$#names; $i++) {
  my $fl = length($names[$i])+2;
  if ($l==0 || $l+$fl > 68){ print FH "\n     +"; $l=9;}
  print FH " $names[$i]", ($i<$#names?',':'');
  $l+=$fl;
}
print FH "\n";
print CH "  };\n}\n#endif\n";
print CC "}\n";

close CH;
close CC;
close FH;
