### vmcw -- 1D Legget-Takagi equation solver

Usage: `vmcw <command file>`

The numerical experiment is done according with the command file. Each
line of the file containes a command for controlling the solver,
modifying 3He parameters etc. Comments are started with `#` symbol,
empty lines are ignored. Command file name without extention is used as a prefix
for some default output files.

### Command list

#### Solver:

* `start` -- (Re)start the solver using a uniform mesh with number of
points specified by `npts` command. If solver was running, use its data
as initial conditions.

* `stop` -- Stop the solver, save its data as
initial conditions.

* `exit` -- Stop the solver and exit.

* `save state <file>` -- Save solver state to a file. Only the solver
state is saved, not experimental conditions (magnetic fields, 3He
parameters, etc.). If solver is not running, then do nothing and print
a warning message.

* `load_state <file>` -- Load solver state from a file. It does not matter
was a solver running previously or not.

* `wait <time, s>` -- Do calculations for some time period.

* `average <duration> <units=s> <npts=20>` -- Average functions over some time
and restart solver using the average value as initial conditions.
Duration can be specified using different units: `s` (secnds), `ms`
(milliseconds), `us` (microseconds), `per_rf` (period of RF pumping,
1/f0), `per_hpd` (period of low-frequency HPD oscillations, depends on
f0, fB, HR0, uses RF field in the cell center if it is non-uniform).

* `reset_time` -- Reset time to 0. This command is useful
when the time is much larger then time step and the calculation fails
with `T + DT = T` error.

* `acc <value>` -- Change solver accuracy. If solver is not running,
the value will be used after it starts. Accuracy can be changed for
the running solver. Default value for accuracy is 2^-20.

* `acc2 <value>` -- Same, but use log2 value (10 for accuracy 2^-10).

* `mindt <value>` -- Change min. time step (seconds). It can be changes
at any time, but real change happenes when the new solver is started.
It's recommended not to change default value (1e-10).

* `npts <number>` -- Change number of points. It is used when a new solver
is started. It is also used for setting non-uniform initial conditions.

* `tstep <value>` -- Set time step. Can be changed during calculation.

* `cell_len <value, cm>` -- Change cell length. Do not change while
the solver is running.

* `adaptive_mesh [<N>] [<density>]` -- Set a non-uniform mesh using
solution from a running solver, restart the solver using the new mesh.
`<N>` - number of points in the mesh (default - use old value).
`<density>` - Non-uniform mesh step: `cell_len/(N-1) / (1+<density>*RMS(U'))`.
0 means that mesh is uniform,  1 means that mesh is twice mere dense if
RMS of function derivatives is 1.

#### Initial and boundary conditions:

* `bcond_type_l <value>`, `bcond_type_r <value>`, `bcond_type <value>` --
  set boundary condition type on left/right/both sides.
  BC types:
  * default - open BC
  * 2 - closed cell, no spin flow through walls
  * 3 - NPD wall, constant functions on the boundary

Following commands change default boundary condition and have an effect
only when the solver is starded next time.

* `set_icond_uniform [<nz value>=1]` -- set uniform i.c. with nz=-1 or nz=+1.
* `set_icond_npd_simple [<nz value>=1]` -- same

* `set_icond_hpd_simple [<ny value>=1]` -- set simple HPD initial
condition: th = acos(-0.25), m=(ny*sin(th),0,cos(th)), n = (0, ny, 0).

* `set_icond_hpd [<ny value>=1]` -- set HPD initial condition with ny=-1
or ny=+1. This function uses 3He parameters, cell size and field profile
to calculate HPD state. Use it after all parameters are set. State with
ny=+1 is stable, ny=-1 is unstable.

* `set_icond_npd [<ny value>=1]` -- set NPD initial condition with ny=-1
or ny=+1. This function uses 3He parameters, cell size and field profile
to calculate HPD state. Use it after all parameters are set. State with
ny=+1 is stable, ny=-1 is unstable. This initial condition
could be bad near resonance because spin diffusion is not taken into account
and sharp resonance peak appeares.


* `set_icond_nsol <width>` -- set i.c witn a simple n-soliton. Width >0 does not work yet
* `set_icond_tsol <width>` -- set i.c witn a simple t-soliton

* `deform <type> [<par>]` -- Update initial condition from current
function values, make some modification, and apply it to the solver. Types:

  * `deform none` - no modifications;
  * `deform half_turn` - rotate alpha_n, alpha_m by pi.
  * `deform rotation <N>` - add constant gradient to alpha_n, alpha_m
    to have N full turns on the cell length.
  * `deform 2pi_soliton <w>` - add a 2-pi soliton in alpha_n, alpha_m
    in the middle of the cell, with width `<w>`.
  * `deform th_soliton <w>` - add a theta soliton in the middle
    of the cell, with width `<w>`.
  * `deform th_soliton1 <w>` - th-soliton type 1 with sharp theta change
  * `deform th_soliton2 <w>` - th-soliton type 2 with sharp theta change
  * `deform th_soliton1a <w>` - th-soliton type 1 with artificial theta change (works only with BC=3)
  * `deform th_soliton2a <w>` - th-soliton type 2 with artificial theta change (works only with BC=3)

After the deformation the solver uses same parameters as before,
pnm_writer is not restarted.

#### Data output:

* `write_profile [<file>] [<N>]` -- Write function profiles to the file.
Meaning of columns depends on `NPDE` setting. For `NPDE=7` it is `x,
mx,my,mz, nx,ny,nz, th, mx',my'...` If filename is not set or set to '-'
then a default filename `<prefix>.prof<n>.dat` is used. If `N` is not set
or `0`, then use solver's mesh. If `N>2`, build a uniform mesh with `N`
points. If solver is not running then only prints a warning.

* `write_mesh [<file>]` -- Write coordinate breakpoints x(i).
Default filename is `<prefix>.mesh<n>.dat`.

* `magn_start [<file>]` -- Start recording components of total
magnetization (what you see in usual NMR experiment) to a file. If
recording is already started then close old file and start new one.
Default filename is `<prefix>.magn<n>.dat`.

* `magn_stop` -- Stop writing magnetization file.

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

* `set_t1 <value, s>`   -- Set relaxation time t_1 [s], or 0 for no relaxation
* `sweep_t1 <value, s> <rate, s/s>` -- Sweep relaxation time t_1 [s]
* `set_tf <value, s>`   -- Set Leggett-Takagi relaxation time tau_f [s], or 0 for no relaxation
* `sweep_tf <value, s> <rate, s/s>` -- Sweep Leggett-Takagi relaxation time tau_f [s]
* `set_diff <value, cm^2/s>` -- Set spin diffusion [cm^2/s]
* `sweep_diff <value, cm^2/s> <rate, (cm^2/s)/s>` -- Sweep spin diffusion [cm^2/s]
* `set_cpar <value, cm/s>` -- Set  spin-wave velocity c_parallel [cm/s]
* `sweep_cpar <value, cm/s> <rate, (cm/s)/s>` -- Sweep spin-wave velocity c_parallel [cm/s]
* `set_leggett_freq <value, Hz>` -- Set Leggett frequency [Hz]
* `sweep_leggett_freq <value, Hz> <rate, Hz/s>` -- Sweep Leggett frequency [Hz]

* `set_ttc_press <T/Tc> <pressure, bar>` -- Set He3 parameters for a
given temperature and pressure (Leggett frequency, spin-wave velocity, spin diffusion, Leggett-takagi relaxation).
Program should be built with `he3lib` support.

* `pin <k> <w>` -- set pinning profile: minimum of Leggett frequency in the center with width xiD*w.
May be useful for stabilizing theta-solitons. 0 (defailt) - no
pinning, 1 - full pinning. `wB = wB*(1- k*exp(-(x/xiD/w)^2))`.
