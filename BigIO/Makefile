HOST := $(shell hostname)

CC = mpicc
FC = mpif90

OBJ = parallel_posix.o 

all: parallel_posix 

parallel_posix: $(OBJ) 
	$(CC) $(OBJ) -o $@

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@

clean:
	rm -f  *.o *.mod 

spotless:
	rm -f *.o *.mod  *~ 

