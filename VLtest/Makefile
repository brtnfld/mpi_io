HDF5_DIR=$(HOME)/packages/hdf5/build/hdf5/bin

CC=$(HDF5_DIR)/h5pcc
CFLAGS=-O3

all: t_pseudo_vlen_bench # t_vlen_bench

t_pseudo_vlen_bench: t_pseudo_vlen_bench.o
	$(CC) $(CFLAGS) -o $@ $<

t_vlen_bench: t_vlen_bench.o
	$(CC) $(CFLAGS) -o $@ $<

clean:
	-rm -f t_pseudo_vlen_bench t_vlen_bench.o
	-rm -f t_pseudo_vlen_bench t_pseudo_vlen_bench.o