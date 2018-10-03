# Create a soliton inside HPD

set_rf_field   1.4e-3    # RF-field
set_field_hz   100    # larmor position
set_field_grad  0.0    # field gradient, G/cm
set_field_quad -0.1    # field gradient, G/cm
set_tf 1e-7
cell_len 0.9
set_rf_prof 0 2
acc2 18
npts 256
mindt 1e-20
tstep 1e-3



start
pnm_start
sweep_field_hz 0 250

set_rf_field  10e-3  # increase the excitation
reset_time
#pnm_hline
tstep 1e-4
sweep_field_hz -2 250

set_rf_field  8e-3
reset_time
#pnm_hline
write_profile
wait 0.004

reset_time
write_profile
wait 0.03

write_profile
save_state sol.state




