#! ../../vmcw

# Start from flipped HPD

set_rf_field   2e-3    # RF-field
set_field_cm    0.2    # larmor position
set_field_grad  0.0
set_field_quad -0.05
tstep 1e-2
acc2 16
npts 256
mindt 1e-20
set_tf 1e-7
set_diff 0.015

load_state hpd_quad.state
reset_time
pnm_start
wait 0.02

tstep 1e-4
deform half_turn
reset_time
pnm_hline
wait 0.055
