#!../../vmcw

npts 400
cell_len 0.9

set_rf_field   0.73e-3     # RF-field
set_field_grad 15.7e-3  # field gradient, 15.7 G/cm/A
set_field_hz   30

tstep 0.005
start                   # start the solver
wait 0.5
tstep 0.001

pnm_start

sweep_field_hz -30 50  # sweep field, 1.5kHz/min


