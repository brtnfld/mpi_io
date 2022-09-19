HDF5_DIR=$(HOME)/packages/hdf5/build/hdf5

CC=mpicc
CFLAGS=-O3

LIBS=-L$(HDF5_DIR)/lib64 -lhdf5 -lm
INC=-I$(HDF5_DIR)/include

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

all: example

example: example.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	-rm -f example example.o
	-rm -f example_cmpd.o
	-rm -f *.h5