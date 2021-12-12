#!../../vmcw

# Using HPD- initial condition (unstable).

set_rf_field   2e-3    # RF-field
set_field_grad  0.1    # field gradient, G/cm
set_field_cm   +0.3    # larmor position

acc2 18
npts 256
mindt 1e-20

set_tf 1e-6  # 0.1 of default value, to have faster stabilization

set_icond_hpd -1

tstep 5e-4
start
pnm_start

wait 0.2

