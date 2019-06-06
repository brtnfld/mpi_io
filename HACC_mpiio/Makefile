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

all: Sample_mpio_measure_time Sample_hdf5_measure_time

Sample_mpio_measure_time: Sample_mpio_measure_time.c
	$(CC) -o $@ $^ $(CFLAGS)

Sample_hdf5_measure_time: Sample_hdf5_measure_time.c
	$(CC2) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o Sample_mpio_measure_time Sample_hdf5_measure_time
