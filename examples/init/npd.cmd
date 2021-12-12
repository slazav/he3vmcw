#!../../vmcw

# Using NPD initial condition.
# It works well only far from resonance because spin diffusion
# is not taken into account.

set_rf_field   2e-3    # RF-field
set_field_grad  0.1    # field gradient, G/cm
set_field_cm    0.3    # larmor position

acc2 18
npts 256
mindt 1e-20

set_tf 1e-6  # 0.1 of default value, to have faster stabilization

set_icond_npd

start
pnm_start

tstep 1e-7
wait 1e-7
# initial profile
write_profile

tstep 5e-4
wait 0.2

# equilibrium profile
write_profile

