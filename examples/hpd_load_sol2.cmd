# save/restore state test

set_rf_field   2e-3    # RF-field
set_field_cm    0.2    # larmor position
set_field_grad  0.1    # field gradient, G/cm
tstep 1e-2

acc2 20
load_state hpd_save.state
pnm_start

acc2 18
make_npd_soliton 0.002

tstep 2e-6
wait 0.00035
write_profile

#wait 0.001
#write_profile


