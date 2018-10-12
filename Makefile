FFLAGS= -Werror -Wconversion\
  -Wintrinsic-shadow -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -Wno-align-commons

# posible solvers:
#  - pde_dp -- oroginal one
#  - epde_dp -- does not work yet
SOLVER = pde_dp
# link with he3lib
HE3LIB = 1


TARGETS=vmcw grad_test grad_speed

LDLIBS= -lgfortran -lm
FFLAGS=-I/usr/include -fno-range-check
CPPFLAGS=-std=c++11
CC=g++

ifeq ($(SOLVER), pde_dp)
  CPPFLAGS+=-DSOLVER_PDE_DP
else ifeq ($(SOLVER), epde_dp)
  CPPFLAGS+=-DSOLVER_EPDE_DP
endif

ifeq ($(HE3LIB), 1)
  CPPFLAGS+=-DHE3LIB
  LDLIBS+=-lhe3
endif


all: $(TARGETS)
clean:
	rm -f $(TARGETS) *.o pdecol/*.o

# Fortran objects
FOBJ=pdecol/$(SOLVER).o vmcw_func.o

# C++ object files
COBJ=vmcw.o pdecol/pdecol_solver.o pnm_writer.o
pdecol/pdecol_solver.o: pdecol/pdecol_solver.h

$(COBJ): %.o: %.cpp
$(FOBJ): %.o: %.f
vmcw.o: pdecol/pdecol_solver.h pnm_writer.h


# grad functions and tests
grad.c grad_test.c: grad.h
grad_test: grad_test.o grad.o
grad_speed: grad_speed.o grad.o

vmcw: $(FOBJ) $(COBJ)

