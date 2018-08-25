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
wait 5e-3
pnm_hline
set_rf_field 7e-3
wait 0.02
pnm_hline
set_rf_field 0


write_profile
wait 2e-3
write_profile
wait 0.03
write_profile

