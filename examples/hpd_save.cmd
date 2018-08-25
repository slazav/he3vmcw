# save/restore state test

set_rf_field   2e-3    # RF-field
set_field_cm   -0.3    # larmor position
set_field_grad  0.1    # field gradient, G/cm
tstep 1e-2

acc2 16
npts 256
mindt 1e-20

start
pnm_start

sweep_field_cm 0.20 0.1   # sweep field
wait 200
save_state hpd_save.state

