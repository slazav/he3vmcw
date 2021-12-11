#!../../vmcw

# Pulsed NMR in field gradient

npts 257
mindt 1e-20
tstep 5e-5
set_t1 0.1
set_leggett_freq 10000

set_field_hz    20
set_field_grad  0
set_field_quad  0.1
set_rf_field 0

start
pnm_start
wait 5e-3
pnm_hline
set_rf_field 3e-3
wait 20e-3
pnm_hline
set_rf_field 0
wait 150e-3

write_profile
