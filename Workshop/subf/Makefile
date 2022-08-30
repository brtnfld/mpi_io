HDF5_DIR=$(HOME)/packages/hdf5/build.polaris/hdf5

CC=cc
CFLAGS=-O3

LIBS=-L$(HDF5_DIR)/lib -lhdf5
INC=-I$(HDF5_DIR)/include

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

all: example example_cmpd

example: example.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

example_cmpd: example_cmpd.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	-rm -f example example.o
	-rm -f example_cmpd example_cmpd.o
	-rm -f *.h5