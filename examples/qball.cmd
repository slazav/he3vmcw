# CW nmr in a field minimum.
# Standing wave with 0,2,4 nodes appear.

npts 257
mindt 1e-20
acc2 22
tstep 5e-2
set_diff 0.02
set_t1 0.01


set_field_grad  0
set_field_quad  0.2
set_rf_field 5e-4

set_field_hz    20
start
pnm_start vectors.pnm

sweep_field_hz    0 2

# mark larmor position with a horizontal line
pnm_hline
sweep_field_hz    -40 2
write_profile
