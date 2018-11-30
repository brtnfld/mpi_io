HOST := $(shell hostname)

ifeq ($(findstring edison, $(HOST)),edison)
 CC = cc
 FC = ftn
else
 ifeq ($(findstring cetus, $(HOST)),cetus)
   CC = mpixlc
   FC = mpixlf90
 else
   CC = mpicc
   FC = mpif90
 endif
endif

F90FLAGS = -O3 -Wall

OBJF90 = bcast.o 
OBJF90LG = bcastlg.o 

all:  bcast bcastlg

bcast: $(OBJF90)
	$(FC) $(F90FLAGS)  $(OBJF90) -o $@

bcastlg: $(OBJF90LG)
	$(FC) $(F90FLAGS)  $(OBJF90LG) -o $@

.SUFFIXES: .o .F90

.F90.o:
	$(FC) $(F90FLAGS) -c $< -o $@

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@

clean:
	rm -f  bcast bcastlg *.o *.mod 

spotless:
	rm -f  *.o *.mod bcast bcastlg *~ timing datafile.mpio
