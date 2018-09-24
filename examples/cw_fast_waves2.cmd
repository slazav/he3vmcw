acc2 16
npts 512
tstep 5e-4

cell_len 0.9
f0 1123000

set_ttc_press 0.93 24.8

set_diff 0.01

set_rf_field    1.5e-3    # RF-field
set_field_grad  0.0  # field gradient, 15.7 G/cm/A
set_field_quad  0.1  # field gradient, 15.7 G/cm/A

set_field_hz   120

start                   # start the solver
pnm_start
#pnm_legend

sweep_field_hz -200 1000

#save_state hpd
#wait 0.05
#hpd_deform
#wait 0.05
