HOST := $(shell hostname)

ifeq ($(findstring edison, $(HOST)),edison)
 CC = cc
 FC = ftn
else
 CC = mpicc
 FC = mpif90
endif

F90FLAGS = -O3

OBJ = ../NonContig/noncon.o
OBJF90 = case2.o

all:  case2

case2: $(OBJ) $(OBJF90)
	$(FC) $(F90FLAGS) $(OBJ) $(OBJF90) -o $@

.SUFFIXES: .o .F90

.F90.o:
	$(FC) $(F90FLAGS) -c $< -o $@

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@

clean:
	rm -f  ../NonContig/*.o *.o *.mod 

spotless:
	rm -f ../NonContig/*.o *.o *.mod case2 *~ timing datafile.mpio

