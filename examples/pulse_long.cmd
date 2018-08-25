npts 257
mindt 1e-20
tstep 5e-4
set_t1 0.05
set_leggett_freq 10000

set_field_hz    0
set_field_grad  0
set_field_quad  0
set_rf_field 0

start
pnm_start
wait 0.05
pnm_hline
set_rf_field 1e-3
wait 0.35
pnm_hline
set_rf_field 0
wait 1

