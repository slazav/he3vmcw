# save s simple HPD state with flat/grad/quad field profile

set_rf_field   2e-3    # RF-field
set_field_cm   -0.3    # larmor position
set_field_grad  0.1    # field gradient, G/cm
tstep 1e-2
acc2 16
npts 256
mindt 1e-20

set_tf 1e-7
set_diff 0.015


start
pnm_start
sweep_field_cm 0.25 0.2
save_state hpd_grad.state
write_profile
pnm_hline

set_field_grad  0.0
wait 0.1
save_state hpd_flat.state
write_profile
pnm_hline

set_field_quad  -0.1
wait 0.1
save_state hpd_quad.state
write_profile
pnm_hline
