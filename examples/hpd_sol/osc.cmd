#!../../vmcw

# load soliton from sol.state

set_rf_field   1.4e-3    # RF-field
set_field_hz   100    # larmor position
set_field_grad  0.0    # field gradient, G/cm
set_field_quad -0.1    # field gradient, G/cm
set_tf 1e-7
cell_len 0.9

set_rf_prof 0 2
set_field_hz -20
set_rf_field  8e-3

acc2 18
mindt 1e-20
tstep 1e-4


load_state hz/hz_20
reset_time
pnm_start

wait 3e-4

set_field_hz -50
wait 3e-2

set_field_hz -80
wait 3e-2

set_field_hz -110
wait 3e-2
