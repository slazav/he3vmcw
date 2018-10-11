# CW nmr in a field minimum.
# Standing wave with 0,2,4 nodes appear.

npts 257
mindt 1e-20
acc2 20
tstep 5e-2
set_diff 0.1
set_t1 0.01

set_field_grad  0
set_field_quad  0.4
set_rf_field 5e-4

set_field_hz    20
start

# wait with high spin diffusion seting
wait 1
# sweep spin diffusion to a small value
sweep_diff 0.015 1


pnm_start
sweep_field_hz    0 2
# mark larmor position with a horizontal line
pnm_hline
sweep_field_hz    -60 2
write_profile

