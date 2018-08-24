# Pulsed NMR in field gradient

npts 257
mindt 1e-20
tstep 5e-5
set_t1 0.1
set_leggett_freq 10000

set_field_hz    0
set_field_grad  0.1
set_field_quad  0
set_rf_field 0

start
pnm_start
wait 5
pnm_hline
set_rf_field 5e-3
wait 20
pnm_hline
set_rf_field 0
wait 50

write_profile
