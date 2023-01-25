HDF5_DIR=$(HOME)/packages/hdf5/build/hdf5

CC=mpicc
CFLAGS=-O3

HDF5_LIB=$(HDF5_DIR)/lib
ifneq "$(wildcard $(HDF5_DIR)/lib64) )" ""
  # if directory xists:
  HDF5_LIB=$(HDF5_DIR)/lib64
endif

LIBS=-L$(HDF5_LIB) -lhdf5 -lm
INC=-I$(HDF5_DIR)/include

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(INC)

all: example

example: example.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	-rm -f example example.o