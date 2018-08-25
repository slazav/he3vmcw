# * create HPD
# * artificially create 1-pi soliton in the center
# * soliton is transformed to NPD wall

set_rf_field   2e-3    # RF-field
set_field_cm   -0.3    # larmor position
set_field_grad  0.1    # field gradient, G/cm

acc2 16
npts 256
mindt 1e-20

start
pnm_start
pnm_legend


sweep_field_cm 0.2 0.1
make_pi_soliton -0.02

tstep 1e-4
wait 0.020

write_profile
