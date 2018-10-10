bcond_type 2

set_rf_field  20e-3    # RF-field
set_field    -20e-3    # larmor position
set_field_grad  0
set_field_quad  0
tstep 1e-5
acc2 22
npts 256
mindt 1e-20
set_tf 0.3e-7
set_diff 0.5

npts 512
load_state hpd_flat.state
deform none
start
reset_time
pnm_start
wait 5e-5

tstep 1e-5
deform th_soliton2 0.02
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

set_cpar 200

wait 2e-4
write_profile

tstep 1e-5
wait 2e-2
write_profile
save_state th2