# save/restore state test

set_rf_field   2e-3    # RF-field
set_field_cm    0.2    # larmor position
set_field_grad  0.11
tstep 1e-2
acc2 16
npts 256
mindt 1e-20
set_tf 1e-7
set_diff 0.015

load_state hpd_grad.state
reset_time
pnm_start
wait 0.02


tstep 1e-4
acc2 14
deform rotation 2
reset_time
pnm_hline
wait 0.03
