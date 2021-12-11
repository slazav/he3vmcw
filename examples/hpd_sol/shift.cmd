#!../../vmcw

# load soliton from sol.state

set_rf_field   1.4e-3    # RF-field
set_field_hz   100    # larmor position
set_field_grad  0.0    # field gradient, G/cm
set_field_quad -0.1    # field gradient, G/cm
set_tf 1e-7
cell_len 0.9

set_rf_prof 0 2
set_field_hz -2
set_rf_field  8e-3

acc2 18
mindt 1e-20
tstep 1e-4


load_state sol.state

#npts 512
#hpd_deform 0
#start

reset_time
pnm_start

sweep_field_hz -6 200
write_profile hz_6.pr
save_state hz_6

sweep_field_hz -8 200
write_profile hz_08.pr
save_state hz_08

sweep_field_hz -10 200
write_profile hz_10.pr
save_state hz_10

sweep_field_hz -12 200
write_profile hz_12.pr
save_state hz_12

sweep_field_hz -14 200
write_profile hz_14.pr
save_state hz_14

sweep_field_hz -16 200
write_profile hz_16.pr
save_state hz_16

sweep_field_hz -18 200
write_profile hz_18.pr
save_state hz_18

sweep_field_hz -20 200
write_profile hz_20.pr
save_state hz_20

sweep_field_hz -22 200
write_profile hz_22.pr
save_state hz_22

sweep_field_hz -24 200
write_profile hz_24.pr
save_state hz_24

sweep_field_hz -26 200
write_profile hz_26.pr
save_state hz_26

sweep_field_hz -28 200
write_profile hz_28.pr
save_state hz_28
