# Short RF pulse (20ms) in a large gradient.
# Two n-solitons appear.
# Same as solitons1.txt, but gradient is larger and RF pumping is smaller

npts 129
mindt 1e-20
acc2 18
tstep 5e-5

set_field_hz    0
set_field_grad  0.2
set_field_quad  0
set_rf_field 0

start
pnm_start
wait 5e-3
pnm_hline
set_rf_field 1.8e-2
wait 0.02
pnm_hline
set_rf_field 0

write_profile
wait 2e-3
write_profile
wait 0.03
write_profile
