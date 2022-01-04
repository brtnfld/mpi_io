CC = h5pcc

SRC = H5GH5L.c 

all: H5GH5L 

H5GH5L: $(SRC)
	$(CC) $(SRC) -o $@

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@

clean:
	rm -f  *.o H5GH5L *.h5 

