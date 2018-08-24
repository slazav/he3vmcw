# HPD in the field minimum

set_rf_field   2e-3    # RF-field
set_field_quad  1
set_field_grad  0.0

acc2 18
npts 256
mindt 1e-20

set_field_hz   40
start                   # start the solver
pnm_start

sweep_field_hz -80 50
