
set_rf_field   2e-3    # RF-field
set_field_cm   -0.3    # larmor position
set_field_grad  0.1    # field gradient, G/cm

acc2 18
npts 256
mindt 1e-20

start
pnm_start
pnm_legend

sweep_field_cm -0.1 0.1   # sweep field

pnm_hline
write_profile
sweep_field_cm  0.0 0.1

pnm_hline
write_profile
sweep_field_cm  0.1 0.1

pnm_hline
write_profile
sweep_field_cm 0.2 0.1

write_profile

acc2 20
set_rf_field   0    # RF-field
tstep 5e-4
wait 125