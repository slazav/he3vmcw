# Makefile

# SHELL=/bin/sh

# Fortran compiler, compiler flags, linker flags
# object files and name of executable

FC=gfortran
OBJS=mainlib.o general.o free.o modu.o texture.o gap.o profiles.o nmr.o blas.o tn.o 
OBJ2=mainprog.o
EXE=texture
LIBR=texture.so
FFLAGS=-O
FFLAGX=-fPIC

all:	$(EXE)

$(EXE):		$(OBJ2) $(LIBR)
		$(FC) -o $(EXE) $(OBJ2) $(LIBR) $(LDFLAGS) -Wl,-rpath=$(PWD)

mainprog.o:	mainprog.f90
		$(FC) -c mainprog.f90 $(FFLAGS)

$(LIBR):	$(OBJS)
		$(FC) -shared -o $(LIBR) $(OBJS) $(LDFLAGS) $(FFLAGS)


mainlib.o:	mainlib.f90 modu.o free.o nmr.o general.o texture.o profiles.o 
		$(FC) -c mainlib.f90 $(FFLAGS) $(FFLAGX)


free.o:		free.f90 texture.o general.o modu.o 
		$(FC) -c free.f90 $(FFLAGS) $(FFLAGX)


nmr.o:		nmr.f90 general.o
		$(FC) -c nmr.f90 $(FFLAGS) $(FFLAGX)


modu.o:		modu.f90 general.o
		$(FC) -c modu.f90 $(FFLAGS) $(FFLAGX)

profiles.o:	profiles.f90 general.o modu.o free.o
		$(FC) -c profiles.f90 $(FFLAGS) $(FFLAGX)


general.o:	general.f90 
		$(FC) -c general.f90 $(FFLAGS) $(FFLAGX)

texture.o:	texture.f90 gap.o general.o
		$(FC) -c texture.f90 $(FFLAGS) $(FFLAGX)

gap.o:		gap.f90 general.o
		$(FC) -c gap.f90 $(FFLAGS) $(FFLAGX)

blas.o:		blas.f 
		$(FC) -c -FI -zero -fno-automatic -std=legacy blas.f $(FFLAGS) $(FFLAGX)

tn.o:		tn.f
		$(FC) -c -FI -zero -fno-automatic -std=legacy tn.f $(FFLAGS) $(FFLAGX)
