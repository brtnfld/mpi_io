#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<sys/types.h>
#include <mpi.h>
#include<time.h>

int main() {

  off_t expand_fs;
  FILE *fd;
  int rc;
  MPI_File fdp;
  MPI_Offset offset, sb_sz;
  int proc_cnt, k;
  int sz_superblock = 2048;
  int superblock[2048];
  MPI_Status *wstatus;

  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  // Get the rank of the process
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  rc = MPI_File_open( MPI_COMM_WORLD, "datafile", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fdp);

  sb_sz = 0;
  proc_cnt = 0;
  offset = 0;

  for(k = 0; k < sz_superblock; ++k ) {
    superblock[k] = k;
  }

  if(myid == proc_cnt){
    MPI_File_write_all(fdp, superblock, sz_superblock, MPI_INTEGER, wstatus);
  } else {
    MPI_File_write_all(fdp, superblock, 0, MPI_INTEGER, wstatus);
  }

  MPI_File_close( &fdp );

  MPI_Barrier(MPI_COMM_WORLD);
  if( myid == 0) {

#if 0
    fd=fopen("datafile","w");
    fclose(fd);
#endif
    
    expand_fs = 134217728;
    
    clock_t tic = clock();
    truncate("datafile", expand_fs);
    clock_t toc = clock();
    
    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}
