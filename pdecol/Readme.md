## A C++ wrapper for PDECOL algorythm

`pdecol_solver` class is a c++ wrapper for Algorithm 540R (remark on
alg.540 aka PDECOL) for solving systems of non-linear differencial equations.

References:

*Algorithm 540: PDECOL, General Collocation Software for Partial Differential Equations.
ACM Transactions on Mathematical Software (TOMS), 5, 326-351 (1979),
https://dl.acm.org/citation.cfm?id=355849

* Algorithm 688: EPDCOL: a more efficient PDECOL code.
ACM Transactions on Mathematical Software (TOMS), 17, 153-166 (1991),
https://dl.acm.org/citation.cfm?id=108558

* Remark on algorithm 540.
ACM Transactions on Mathematical Software (TOMS), 18, 343-344 (1992),
https://dl.acm.org/citation.cfm?id=131773

Algorithm 688 is not working in this wrapper yet.

## Usage

User should provide functions for pdecol: `f_`, `uinit_`, `bndry_`, `derivf_`
(see comments in the PDECOL code).

Constructor:
```c++
pdecol_solver( double t0, double mindt, double eps, std::vector<double> & XBKPT,
               int NPDE, int KORD=4, int NCC=2, int MF=22, int verbose=1

```
* t0    -- Initial time.
* mindt -- Min time step (recommended 1e-10).
* eps   -- Accuracy. Can be changed during calculations (see ch_eps() below)
* XBKPT -- strictly increasing array of piecewise polynomial breakpoints
* NPDE  -- number of differential equations
* KORD  -- polynom.order (recommended 4)
* NCC   -- number of cont.cond (recommended 2)
* MF    -- The method flag (11,12,21 or 22). Can be changed during calculations (see ch_mf() below)
* verbose -- verbosity level.


Do the calculation until time `t`:
```c++
int step(double t);
```


