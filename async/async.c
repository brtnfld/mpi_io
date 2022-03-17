/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* Simple MPI-IO program to measure the performance
 * of async pattern.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>
#include <math.h>

#define PRINTID printf("Proc %d: ", mpi_rank)

static int PImpi(int pe, int processes, int intervals) {

    double time1 = MPI_Wtime();

    int count = intervals / processes;
    int start = count * pe;
    int end = count * pe + count;

    int i;
    double subtotal, total = 0;
    for (i = start; i < end; ++i) {
        subtotal += pow(-1, i) / (2 * i + 1);
    }

    MPI_Reduce(&subtotal, &total, 1, MPI_DOUBLE, MPI_SUM,
        0, MPI_COMM_WORLD);

    double time2 = MPI_Wtime();

    if (pe == 0) {
        total = total * 4;
        printf("Result:   %.10lf\n", total);
        printf("Time:     %.10lf\n", time2 - time1);
    }

    return 0;
}

int main(int ac, char **av)
{
    int  mpi_size, mpi_rank;
    MPI_File fh;
    char *filename = "./mpitest.data";
    char mpi_err_str[MPI_MAX_ERROR_STRING];
    int  mpi_err_strlen;
    int  mpi_err;
    char expect_val;
    int  i=0 ; 
    int  nerrors = 0;		/* number of errors */
    /* buffer size is the total size for one variable. */
    // This will make the variable size 5GB(5*1024*1024*1024). 
    // Uncomment the following line
    //int64_t  buf_size = 5368709120LL;
    //Uncomment the following line makes the variable size 50GB.
    //int64_t  buf_size = 53687091200LL;
    
    //int64_t  buf_size = 536870912LL;
    int64_t  buf_size; // = 10240LL;
    /* Number of variables, currently is 9 like Generic IO. */
    int num_vars  = 2;
    MPI_Offset  mpi_off = 0;
    MPI_Status  mpi_stat;
    char* writedata = NULL;

    long  buf_size_per_proc = 0;
    double mpiio_stime =0;
    double mpiio_etime = 0;
    double total_time = 0;
    double Max_total_time = 0;
    double Min_total_time = 0;
    double Sum_total_time = 0;
    double rate = 0;
    
    MPI_Request* request;
    MPI_Info info = MPI_INFO_NULL;
    MPI_Init(&ac, &av);

    //    MPI_Info_create(&info);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    buf_size = 1*1024*1024*1024;

    /* Buffer size per process per variable, The smallest unit for one process to write */
    buf_size_per_proc = buf_size/mpi_size;

    if (mpi_rank==0){
	printf("Testing simple C MPI collective IO program with %d processes accessing file %s\n",
	    mpi_size, filename);
        printf("This tests the MPIO different patterns for %d variables with MPIO write.\n", num_vars);
    }

    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, filename,
	    MPI_MODE_RDWR | MPI_MODE_CREATE ,
	    info, &fh))
	    != MPI_SUCCESS){
	    MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	    PRINTID;
	    printf("MPI_File_open failed (%s)\n", mpi_err_str);
	    MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Allocate total amount of data to the buffer
    writedata = malloc(buf_size*num_vars);
    request = malloc(num_vars*sizeof(MPI_Request));

    /* each process writes some data */
  
    if(mpi_rank == 0) 
      printf("Start write\n");

    mpi_off = buf_size_per_proc*mpi_rank;

    for (i=0; i < num_vars; i++) {
      if ((mpi_err = MPI_File_iwrite_at_all(fh, mpi_off, writedata, buf_size_per_proc, MPI_BYTE,&request[i]))
          != MPI_SUCCESS){
        MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
        PRINTID;
        printf("MPI_File_write_at_all offset(%ld), bytes (%d), failed (%s)\n",
               (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
      mpi_off+=buf_size;     
    }

    mpiio_stime = MPI_Wtime();
    PImpi(mpi_rank, mpi_size, 1E8);
    mpiio_etime = MPI_Wtime();
    total_time = mpiio_etime - mpiio_stime;
    if(mpi_rank == 0) {
      printf("Compute time %f seconds. \n", total_time);
    }

    mpiio_stime = MPI_Wtime();
    for (i=0; i < num_vars; i++) {
      MPI_Wait( &request[i], &mpi_stat);
    }
    mpiio_etime = MPI_Wtime();

    total_time = mpiio_etime - mpiio_stime;
     
    MPI_Reduce(&total_time, &Max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &Sum_total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time, &Min_total_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if(mpi_rank == 0) {
      rate = (double)(buf_size*num_vars)/Max_total_time/(1024.*1024.);
      printf("%d Procs Wrote %d variables in %f seconds. \n",mpi_size,num_vars,Max_total_time);
      printf(" Bandwidth is %f MB/s.\n",rate);
      printf("Total IO time for all processes is %f seconds.\n",Sum_total_time);
      printf("Minimum IO time for all processes is %f seconds.\n",Min_total_time);
      rate = (double)(buf_size*num_vars)/(Sum_total_time/mpi_size)/(1024.*1024.);
      printf("Average IO time for all processes is %f seconds.\n",Sum_total_time/mpi_size);
      printf(" Average Bandwidth is %f MB/s.\n",rate);
    }

    free(writedata);
    free(request);
    MPI_File_close(&fh);
    MPI_Finalize();

#if 0
    /* each process reads all data and verify. */
    for (irank=0; irank < mpi_size; irank++){
	mpi_off = irank*DIMSIZE;
	if ((mpi_err = MPI_File_read_at(fh, mpi_off, readdata, DIMSIZE, MPI_BYTE,
		&mpi_stat))
		!= MPI_SUCCESS){
	    MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	    PRINTID;
	    printf("MPI_File_read_at offset(%ld), bytes (%d), failed (%s)\n",
		    (long) mpi_off, (int) DIMSIZE, mpi_err_str);
	    MPI_Abort(MPI_COMM_WORLD, 1);
	};
	for (i=0; i < DIMSIZE; i++){
	    expect_val = irank*DIMSIZE + i;
	    if (readdata[i] != expect_val){
		PRINTID;
		printf("read data[%d:%d] got %d, expect %d\n", irank, i,
			readdata[i], expect_val);
		nerrors++;
	    }
	}
    }
    if (nerrors)
	MPI_Abort(MPI_COMM_WORLD, 1);
    MPI_File_close(&fh);
#endif


    return 0;
}
