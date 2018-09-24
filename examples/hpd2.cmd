acc2 16
set_icond_hpd2
npts 512
tstep 1e-5

cell_len 0.9
f0 1123000


set_tf 1.6e-6
set_ttc_press 0.93 24.8
set_diff 0.01

set_rf_field    2e-3    # RF-field
set_field_grad  -0.1  # field gradient, 15.7 G/cm/A

set_field_hz   -120
start                   # start the solver
pnm_start
wait 0.0002
write_profile

wait 0.002
write_profile

wait 0.002
write_profile

wait 0.002
write_profile

wait 0.002
write_profile

wait 0.002
write_profile

wait 0.002
write_profile
wait 0.002
write_profile
wait 0.002
write_profile
wait 0.002
write_profile
wait 0.002
write_profile
wait 0.002
write_profile
wait 0.002
write_profile
wait 0.002
write_profile
wait 0.002
write_profile
wait 0.002
write_profile

wait 0.002
write_profile
