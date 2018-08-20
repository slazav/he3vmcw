#!/usr/bin/perl
use strict;
use warnings;

# create vmcw_pars.h and vmcw_pars.cpp

# parameters (type, name, default, comment):
my $pars ="
double t11       1.0     1/t1, s^-1
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
double CELL_LEN  0.4     cell length, cm
double XMESH_K   0.2     mesh step: CELL_LEN/(NPTS-1) / (1+XMESH_K*AER_STEP)
double XMESH_ACC 1e-10   mesh accuracy
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
open CC, "> vmcw_pars.cpp" or die "can't open vmcw_pars.cpp: $!\n";

######
print CH qq*// This file is created by $0 script. Do not modify.
#ifndef VMCW_PARS_H
#define VMCW_PARS_H
void set_def_pars(struct pars_t \*p);
struct pars_t {
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
  printf CH "    %8s %12s;$ccomm\n", $ctype, $name;
  printf CC "  p->%-12s = $def;\n", $name;
  push @names, $name;
}

print CH "};\n#endif\n";
print CC "}\n";

close CH;
close CC;
