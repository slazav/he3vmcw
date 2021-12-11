#!../../vmcw

# Create HPD by sweeping magnetic field with a gradient,
# then switch off the excitation and sobserve decay of HPD.

set_rf_field   1.2e-3    # RF-field
set_field_hz   40    # larmor position
set_field_grad  0.0    # field gradient, G/cm
set_field_quad -0.1    # field gradient, G/cm

acc2 18
npts 256
mindt 1e-20

start
pnm_start

tstep 1e-3
sweep_field_hz -20 100


