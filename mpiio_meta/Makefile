HOST := $(shell hostname)

ifeq ($(findstring edison, $(HOST)),edison)
 CC = cc
 FC = ftn
else
 CC = mpicc
 FC = mpif90
endif

F90FLAGS = -O3

OBJF90 = case1.o

all:  case1

case1: $(OBJF90)
	$(FC) $(F90FLAGS)  $(OBJF90) -o $@

.SUFFIXES: .o .F90

.F90.o:
	$(FC) $(F90FLAGS) -c $< -o $@

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@

clean:
	rm -f  case1  *.o *.mod 

spotless:
	rm -f  *.o *.mod case1 *~ timing datafile.mpio
