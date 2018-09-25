npts 512
acc2 16
tstep 5e-3
mindt 1e-20
cell_len 0.9
set_freq 1123000
set_ttc_press 0.93 24.8
set_diff 0.01

set_rf_field    5e-3    # RF-field
set_field_grad  0.0157  # field gradient, 15.7 G/cm/A
set_field_quad -0.01    #

#set_field_hz   280
#start                   # start the solver
#pnm_start
#sweep_field_hz -40 25
#save_state hpd

tstep 5e-5
set_field_hz   -40
load_state hpd
pnm_start

wait 0.001
tstep 5e-5
pnm_hline
hpd_deform 0
#reset_time
wait 0.03
