### vmcw -- 1D Legget-Takagi equation solver

Usage: `vmcw <command file>`

The numerical experiment is done according with the command file. Each
line of the file containes a command for controlling the solver,
modifying 3He parameters etc. Comments are started with `#` symbol,
empty lines are ignored. Command file name without extention is used as a prefix
for some default output files.

### Command list

#### Solver:

* `start` -- (Re)start the solver. Coordinate mesh is saved into
  <prefix>.mesh file (4 columns: point number, coordinate value, "aerogel
  density", "aerogel density" derivative).

* `stop` -- Stop the solver.

* `exit` -- Stop the solver and exit.

* `save state <file>` -- Save solver state to a file. Only the solver
state is saved, not experimental conditions (magnetic fields, 3He
parameters, etc.). If solver is not running, then do nothing and print
a warning message.

* `load_state <file>` -- Load the solver state from a file and (re)start the
solver.

* `wait <time, s>` -- Do calculations for some time period.

* `reset_time` -- Reset time to 0. This command is useful
when the time is much larger then time step and the calculation fails
with `T + DT = T` error.

* `acc <value>` -- Change solver accuracy. If solver is not running,
the value will be used after it starts.

* `acc2 <value>` -- Same, but use log2 value (10 for accuracy 2^(-10)).

* `mindt <value>` -- Change min. time step (recommended 1e-10). It
can be changes at any time, but real change happenes when the new
solver is started.

* `npts <number>` -- Change number of points. It can be changed at any time,
but real change happenes when the new solver is started.

* `tstep <value>` -- Set time step. Can be changed during calculation.

* `cell_len <value, cm>` -- Change cell length. Do not change while
the solver is running.

* `mesh_k`, `aer_len`, `aer_cnt`, `aer_trw` -- Non-uniform mesh density,
"aerogel" length, "aerogel" center, "aerogel" transition width. Do not
change while the solver is running. Used in non-uniform "aerogel"
calculation, not tested for a long time. For normal calculations
`aer_len` should be less or equal then zero.

#### Initial and boundary conditions:

* `bcond_type <value>` -- set boundary condition type (not tested)

Following commands change default boundary condition without restarting the solver.

* `set_icond_uniform [<nz value>]` -- set uniform i.c. with nz=-1 or nz=+1 (default)
* `set_icond_hpd` -- set HPD i.c. with ny=-1 or ny=+1 (default).
* `set_icond_hpd2`
* `set_icond_nsol <width>` -- set i.c witn a simple n-soliton. Width >0 does not work yet
* `set_icond_tsol <width>` -- set i.c witn a simple t-soliton

Following commands update initial condition from current function values,
modify it somehow and restart the solver.

* `make_2pi_soliton` -- Update initial condition from current function values,
make 2-pi hpd soliton, and restart the solver.

* `make_pi_soliton` -- Update initial condition from current function values,
make pi hpd soliton, restart solver.

* `make_npd_soliton` -- Update initial condition from current function values,
 make npd soliton, restart solver

* `hpd_deform <type>` -- Update initial condition from current function values,
make some modification, and restart the solver. Types:
  * 0 - no modifications;
  * 1 - invert theta angle


#### Data output:

* `write_profile [<file>]` -- Write function profiles to the file.
8 columns: coordinate, mx,my,mz, nx,ny,nz, theta. Default filename is
`<prefix>.prof<n>.dat`. If solver is not running then only prints a warning.

* `pnm_start [<file>]` -- Start writing a pnm-file with vectors n and M
orientation. If solver is not started do nothing and print a warning
message. Default filename is `<prefix>.pic.pnm`. A few pnm writers
can be opened. Pnm writers are stopped either then the solver is stopped
or then `pnm_stop` command is executed.

* `pnm_legend [<file>]` -- Start drawing a legend in the pnm file.
If the pnm writer for this file is not running only prints a warning.

* `pnm_hline [<file>]` -- Draw a horizontal line in the pnm file.
If the pnm writer for this file is not running only prints a warning.

* `pnm_stop [<file>]` -- Stop the pnm writer. If the pnm writer for this
file is not running only prints a warning.


#### External conditions:

* `set_freq <value, Hz>` -- Set NMR frequency.

* `set_field <value, G>` -- Set uniform field, in Gauss from larmor.

* `step_field <value, G>` -- Uniform field step [G].

* `set_field_grad <value, G/cm>` -- Set field gradient.

* `set_field_quad <value, G/cm^2>` -- Set field quadratic term.

* `sweep_field <destination, G> <rate, G/s>` -- Sweep uniform field:
destination [G from Larmor], rate [G/s]

* `set_field_hz <value, Hz>` --  Set uniform field in frequency shift
units [Hz from NMR freq].

* `step_field_hz <value, Hz>` -- Uniform field step [Hz].

* `sweep_field_hz <destination, Hz> <rate, Hz/s>` -- Sweep uniform field:
destination [Hz from NMR frequency], rate [Hz/s].

* `set_field_cm <value, cm>` -- Set uniform field in Larmor position
units [cm]. Gradient term is used to convert field to cm. Quadratic term
is not used. Lower Larmor positon means higher field. If the gradient
term is zero do nothing and print a warning message.

* `step_field_cm <value, cm>` -- Uniform field step [cm].

* `sweep_field_cm <destination, cm> <rate, cm/s>` -- Sweep uniform field.

* `set_rf_field <value, G>` -- Set RF field.

* `step_rf_field <value, G>` -- Step RF field.

* `sweep_rf_field <destination, G> <rate, G/s>` -- Sweep RF field.

#### He3 properties:

* `set_t1 <value, s>`   -- Set relaxation time t_1 [s]
* `sweep_t1 <value, s> <rate, s/s>` -- Sweep relaxation time t_1 [s]
* `set_tf <value, s>`   -- Set Leggett-Takagi relaxation time tau_f [s]
* `sweep_tf <value, s> <rate, s/s>` -- Sweep Leggett-Takagi relaxation time tau_f [s]
* `set_diff <value, cm^2/s>` -- Set spin diffusion [cm^2/s]
* `sweep_diff <value, cm^2/s> <rate, (cm^2/s)/s>` -- Sweep spin diffusion [cm^2/s]
* `set_cpar <value, cm/s>` -- Set  spin-wave velocity c_parallel [cm/s]
* `sweep_cpar <value, cm/s> <rate, (cm/s)/s>` -- Sweep spin-wave velocity c_parallel [cm/s]
* `set_leggett_freq <value, Hz>` -- Set Leggett frequency [Hz]
* `sweep_leggett_freq <value, Hz> <rate, Hz/s>` -- Sweep Leggett frequency [Hz]

* `set_ttc_press <T/Tc> <pressure, bar>` -- Set He3 parameters for a
given temperature and pressure. Program should be built with `he3lib` support.

