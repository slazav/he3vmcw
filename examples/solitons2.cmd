# Short RF pulse (20ms) in a large gradient.
# Two n-solitons appear.
# Same as solitons1.txt, but gradient is larger and RF pumping is smaller

npts 257
mindt 1e-20
acc2 18
tstep 5e-5

set_field_hz    0
set_field_grad  0.2
set_field_quad  0
set_rf_field 0

start
pnm_start
wait 5
pnm_hline
set_rf_field 1.8e-2
wait 20
pnm_hline
set_rf_field 0

write_profile
wait 2
write_profile
wait 30
write_profile