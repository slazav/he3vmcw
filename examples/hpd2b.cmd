acc2 16
set_icond_hpd2
npts 512
tstep 2e-3

cell_len 0.9

set_ttc_press 0.93 24.8

set_diff 0.01

set_rf_field    2e-3    # RF-field
set_field_grad  0.0  # field gradient, 15.7 G/cm/A
set_field_quad  0.1  # field gradient, 15.7 G/cm/A

set_field_hz   -120
start                   # start the solver
pnm_start

set_tf 1e-7
wait 0.02
write_profile

sweep_field_hz   30 100
