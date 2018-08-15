FFLAGS= -Werror -Wconversion\
  -Wintrinsic-shadow -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -Wno-align-commons

TARGETS=vmcw
LDLIBS=-lhe3 -lgfortran -lm
FFLAGS=-I/usr/include -fno-range-check
CC=g++

all: $(TARGETS)
clean:
	rm -f $(TARGETS) *.o vmcw_pars.h vmcw_pars.cpp

# To share main parameter structure beteen C++ and Fortran
# we create a .fh file from .h
vmcw_pars.h vmcw_pars.cpp: vmcw_pars.pl
	./$<

# Fortran objects
FOBJ=pdecol/pde_dp.o vmcw_func.o he3_funcs.o

# C++ object files
COBJ=vmcw.o pdecol/pdecol_solver.o vmcw_pars.o vmcw_mesh.o pnm_writer.o

$(COBJ): %.o: %.cpp
$(FOBJ): %.o: %.f par.fh
vmcw.o pdecol_solver.o: pdecol/pdecol_solver.h vmcw_pars.h pnm_writer.h

vmcw: $(FOBJ) $(COBJ)

