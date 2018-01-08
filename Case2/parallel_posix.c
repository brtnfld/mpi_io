#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/types.h>
#include <mpi.h>
#include<time.h>

int main() {

  off_t expand_fs;
  FILE *fd;

  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  // Get the rank of the process
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if( myid == 0) {

    fd=fopen("datafile","w");
    fclose(fd);
    
    expand_fs = 33554432;
    
    clock_t tic = clock();
    truncate("datafile", expand_fs);
    clock_t toc = clock();
    
    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
