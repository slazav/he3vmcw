#!../../vmcw

# HPD in the field minimum

set_rf_field   3e-3    # RF-field
set_field_quad  1
set_field_grad  0.0

acc2 18
npts 256
mindt 1e-20

set_field_hz   40
start                   # start the solver
pnm_start

sweep_field_hz -150 50

set_rf_field   0    # RF-field
tstep 5e-4
wait 0.2

