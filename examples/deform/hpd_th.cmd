set_rf_field   2e-3    # RF-field
set_field_cm    0.2    # larmor position
set_field_grad  0
set_field_quad  0
tstep 1e-2
acc2 22
npts 256
mindt 1e-20
set_tf 1e-7
set_diff 0.025

load_state hpd_flat.state
reset_time
pnm_start
wait 0.02

tstep 1e-5
deform th_soliton 0.01
reset_time
pnm_hline
wait 5e-3
write_profile
