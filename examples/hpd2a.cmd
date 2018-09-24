acc2 16
set_icond_hpd2
npts 512
tstep 1e-4

cell_len 0.9
f0 1123000

set_ttc_press 0.93 24.8

set_diff 0.01

set_rf_field    2e-3    # RF-field
set_field_grad  0  # field gradient, 15.7 G/cm/A

set_field_hz   -120
start                   # start the solver
pnm_start

set_tf 0.5e-5
wait 0.02
write_profile

sweep_tf 1.0e-7 1e-4
wait 0.02
write_profile

sweep_tf 2.0e-7 1e-4
wait 0.02
write_profile

sweep_tf 4.0e-7 1e-4
wait 0.02
write_profile

sweep_tf 8.0e-7 1e-4
wait 0.02
write_profile

