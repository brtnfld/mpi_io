HDF5_DIR=$(HOME)/packages/hdf5_msb/build/hdf5

CC=mpicc
CC2=$(HDF5_DIR)/bin/h5pcc

OBJ1=Sample_mpio_measure_time.o
OBJ2=Sample_hdf5_measure_time.o

all: Sample_mpio_measure_time Sample_hdf5_measure_time

Sample_mpio_measure_time: $(OBJ1)
	$(CC) -o $@ $^ $(CFLAGS)

Sample_hdf5_measure_time: $(OBJ2)
	$(CC2) -o $@ $^ $(CFLAGS)
