HDF5_DIR=$(HOME)/packages/hdf5/build/hdf5/bin

CC=$(HDF5_DIR)/h5pcc
CFLAGS=-O3

all: example example_cmpd

example: example.o
	$(CC) $(CFLAGS) -o $@ $<

example_cmpd: example_cmpd.o
	$(CC) $(CFLAGS) -o $@ $<

clean:
	-rm -f example example.o
	-rm -f example_cmpd example_cmpd.o
	-rm -f *.h5