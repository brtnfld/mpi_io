HDF5_DIR=$(HOME)/packages/hdf5_msb/build/hdf5

CC=mpicc
CC2=$(HDF5_DIR)/bin/h5pcc


all: Sample_mpio_measure_time Sample_hdf5_measure_time

Sample_mpio_measure_time: Sample_mpio_measure_time.c
	$(CC) -o $@ $^ $(CFLAGS)

Sample_hdf5_measure_time: Sample_hdf5_measure_time.c
	$(CC2) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o Sample_mpio_measure_time Sample_hdf5_measure_time
