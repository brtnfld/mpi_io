#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/types.h>
#include <mpi.h>
#include<time.h>

static void raw(MPI_File fh, int rank, int nprocs, int bufsize,  MPI_Offset offset) {

  MPI_Offset LOCoffset;
  int *buf = NULL;
  int k;
  MPI_Status wstatus;
    
  buf = (int *)malloc(bufsize*sizeof(int));
  if (buf == NULL)
    printf("alloc failed \n");

  for(k = 0; k < bufsize; k++ ) {
    buf[k] = k;
  }
  LOCoffset = offset + rank*bufsize*sizeof(int);

  MPI_File_set_view(fh, LOCoffset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, buf, bufsize, MPI_INTEGER, &wstatus);

  free(buf);
  return;

}

int main(int argc, char *argv[] ) {

  int64_t N;
  FILE *fd;
  int rc, GB;
  MPI_File fdp;
  MPI_Offset offset;
  int bufsize;

  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  // Get the rank of the process
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  rc = MPI_File_open( MPI_COMM_WORLD, "datafile", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fdp);

  offset = 0;
  N = 1024*1024*1024/sizeof(int);
  GB = 2; /* SIZE OF FILE, GB */
  bufsize = GB*N/nprocs;

  raw(fdp, myid, nprocs, bufsize, offset);  

  MPI_File_close( &fdp );

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();

  return 0;
}
