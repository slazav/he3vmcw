set_rf_field   10e-3    # RF-field
set_field_cm   -0.3    # larmor position
set_field_grad  0.05    # field gradient, G/cm

bcond_type_r 3
tstep 0.02

acc2 20
npts 512
mindt 1e-20

start
pnm_start
pnm_legend

sweep_field_cm 0.1 0.1   # sweep field
sweep_field_cm 0.2 0.02   # sweep field

save_state wall1