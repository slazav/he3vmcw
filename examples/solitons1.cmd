# Short RF pulse (20ms) in a large gradient.
# Two n-solitons appear.

npts 257
mindt 1e-20
acc2 16
tstep 5e-5
set_t1 0.1

set_field_hz    0
set_field_grad  0.1
set_field_quad  0
set_rf_field 0

start
pnm_start
wait 5
pnm_hline
set_rf_field 7e-3
wait 20
pnm_hline
set_rf_field 0


write_profile
wait 2
write_profile
wait 30
write_profile
