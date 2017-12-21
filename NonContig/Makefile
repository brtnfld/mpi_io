CC = mpicc
FC = mpif90
CC = cc
FC = ftn

F90FLAGS = -O3

OBJ = noncon.o
OBJF90 = noncontig.int.o

all:  noncontig.int

noncontig.int: $(OBJ) $(OBJF90)
	$(FC) $(F90FLAGS) $(OBJ) $(OBJF90) -o $@

.SUFFIXES: .o .F90

.F90.o:
	$(FC) $(F90FLAGS) -c $< -o $@

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@

clean:
	rm -f *.o *.mod noncontig.int

spotless:
	rm -f *.o *.mod noncontig.int *~ timing datafile

