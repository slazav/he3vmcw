acc2 16
icond_type 12
npts 512
tstep 5e-5

CELL_LEN 0.9
f0 1123000

set_t1 1e-3
set_tf 7.6e-6
set_ttc_press 0.93 24.8
set_diff 0.01

set_rf_field    2e-3    # RF-field
set_field_grad 15.7e-3  # field gradient, 15.7 G/cm/A

set_field_hz   -25
start                   # start the solver
pnm_start hpd_2pi.pnm   # draw M and N vectors to this file
wait 0.04
