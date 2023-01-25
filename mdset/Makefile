HDF5_DIR=$(HOME)/packages/hdf5/build/hdf5
ifneq ($(UNAME),"")
  HDF5_DIR=$(HOME)/packages/hdf5/build.$(UNAME)/hdf5
endif

CC=mpicc
ifeq ($(UNAME),crusher)
   CC=cc
endif

CFLAGS=-O3

HDF5_LIB=$(HDF5_DIR)/lib
ifneq ($(wildcard $(HDF5_DIR)/lib64),)
  # if directory exists:
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
