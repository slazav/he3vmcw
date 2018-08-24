npts 400
CELL_LEN 0.9
f0 1124000

set_rf_field   0.61e-3     # RF-field
set_field_grad 5.7e-3  # field gradient, 15.7 G/cm/A
set_field_quad -25e-3     # field gradient, 15.7 G/cm/A
set_field_hz   30

tstep 0.005
start                   # start the solver
wait 500
tstep 0.001

pnm_start

sweep_field_hz -30 25  # sweep field, 1.5kHz/min

