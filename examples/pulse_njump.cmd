npts 257
mindt 1e-20
acc2 18
tstep 5e-4
set_t1 0.05
set_leggett_freq 10000

set_field_hz    0
set_field_grad  0
set_field_quad  0
set_rf_field 0

start
pnm_start
wait 50
pnm_hline
set_rf_field 3e-3
wait 50
pnm_hline
set_rf_field 0
wait 1000

