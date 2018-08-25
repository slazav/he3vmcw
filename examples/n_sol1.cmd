acc2 9
set_icond_nsol 0
mindt 1e-10
tstep 0.001

CELL_LEN 0.9
f0 1124000


set_rf_field   1.05e-3     # RF-field
set_field_grad -15.7e-3  # field gradient, 15.7 G/cm/A
set_field_hz   30

set_leggett_freq 50000

npts 267
mindt 1e-10

start                    # start the solver
pnm_start

sweep_field_hz -30 100  # sweep field, 1.5kHz/min

