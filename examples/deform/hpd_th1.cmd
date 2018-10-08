set_rf_field   2e-3    # RF-field
set_field    -20e-3    # larmor position
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

tstep 1e-6
deform th_soliton1 0.02
reset_time
pnm_hline


wait 2e-6
write_profile

wait 15e-6
write_profile

wait 2e-6
write_profile

wait 2e-6
write_profile

wait 2e-5
write_profile

wait 2e-4
write_profile
