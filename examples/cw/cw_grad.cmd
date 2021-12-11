#!../../vmcw

# Fast sweep of magnetic field in presence of gradient and RF pumping.
# Ringing at Larmor frequency + some spin waves.

acc2 16
npts 512
tstep 5e-4

cell_len 0.9

set_ttc_press 0.93 24.8
set_diff 0.01

set_rf_field    1.5e-3    # RF-field
set_field_grad  0.1  # field gradient, 15.7 G/cm/A
set_field_quad  0.0  # field gradient, 15.7 G/cm/A

set_field_hz   200

start                   # start the solver
pnm_start
# pnm_legend

sweep_field_hz -200 1000
