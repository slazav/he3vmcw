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
	rm -f $(TARGETS) *.o vmcw_pars.fh vmcw_pars.h

# To share main parameter structure beteen C++ and Fortran
# we create a .fh file from .h
vmcw_pars.fh vmcw_pars.h vmcw_pars.cpp: vmcw_pars.pl
	./$<

# Fortran objects
FOBJ=vmcw_f.o libs/pde_dp.o\
     vmcw_mon.o\
     vmcw_cmd.o vmcw_func.o he3_funcs.o

# C++ object files
COBJ=vmcw.o vmcw_pdecol.o vmcw_pars.cpp vmcw_mesh.o

$(COBJ): %.o: %.cpp
$(FOBJ): %.o: %.f par.fh vmcw_pars.fh
vmcw.o vmcw_pdecol.o: vmcw_pdecol.h vmcw_pars.h

vmcw: $(FOBJ) $(COBJ)

