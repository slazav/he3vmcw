#!../../vmcw

# Using HPD initial condition.
# Gradient of magnetic field applied,
# Larmor condition is in the middle of the cell,
# HPD on the right side is unstable,
# boundary is formed

set_rf_field   2e-3    # RF-field
set_field_grad  0.1    # field gradient, G/cm
set_field_cm    0    # larmor position

acc2 18
npts 256
mindt 1e-20

set_icond_hpd_simple

tstep 5e-4
start
pnm_start

wait 0.4

