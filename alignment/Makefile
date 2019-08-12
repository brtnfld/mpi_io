HDF5_DIR=$(HOME)/packages/hdf5/build/hdf5
HOSTNAME:=$(shell hostname -d)

DEF=
ifneq (,$(findstring summit,$(HOSTNAME)))
   DEF=-DSUMMIT
endif

CC=mpicc $(DEF)
CC2=$(HDF5_DIR)/bin/h5pcc $(DEF)

ifneq (,$(findstring nersc,$(HOSTNAME)))
	CC=cc
	CC2=./h5pcc
endif

all: h5pset_alignment

h5pset_alignment: h5pset_alignment.c
	$(CC2) -DPROC0_MPIIO -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o  h5pset_alignment
