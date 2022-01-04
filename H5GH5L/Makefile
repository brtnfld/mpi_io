SRC=$(wildcard *.c)
BIN=$(SRC:.c=)

CC = h5pcc

all: $(BIN)

H5GH5L: $(SRC)
	$(CC) $(SRC) -o $@

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@

clean:
	rm -f  *.o $(BIN) *.h5 

