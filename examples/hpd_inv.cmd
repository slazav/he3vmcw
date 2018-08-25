acc2 10
icond_type 11
npts 512
tstep 1e-4

CELL_LEN 0.9
f0 1123000

set_t1 1e7
set_tf 7.6e-6
set_ttc_press 0.93 24.8
set_diff 0.01

set_rf_field    1e-3    # RF-field
set_field_grad 15.7e-3  # field gradient, 15.7 G/cm/A
set_rf_prof     0 2

set_field_hz   -30
start                   # start the solver
pnm_start
wait 0.099
