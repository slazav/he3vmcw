FFLAGS= -Werror -Wconversion\
  -Wintrinsic-shadow -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -Wno-align-commons

TARGETS=vmcw vft
LDLIBS=-lhe3 -lgfortran
FFLAGS=-I/usr/include -fno-range-check
CC=g++

all: $(TARGETS)
clean:
	rm -f $(TARGETS) *.o

FOBJ=vmcw_f.o libs/pde_dp.o\
     vmcw_mesh.o vmcw_mon.o\
     vmcw_cmd.o vmcw_func.o he3_funcs.o
COBJ=vmcw.o vmcw_pdecol.o

$(COBJ): %.o: %.cpp

$(FOBJ): %.o: %.f vmcw.fh par.fh
vmcw.o vmcw_pdecol.o: vmcw_pdecol.h

vmcw: $(FOBJ) $(COBJ)

vft: vft.f

