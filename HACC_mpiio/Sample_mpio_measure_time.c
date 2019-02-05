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
 * of interleaved pattern vs contiguous pattern.
 * KY
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>
#ifndef MPI_FILE_NULL           /*MPIO may be defined in mpi.h already       */
#   include <mpio.h>
#endif

#define PRINTID printf("Proc %d: ", mpi_rank)

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
    /* The buffer size will be 50 GB */
    //int64_t  buf_size = 53687091200LL;
    // The buffer size will be 5 GB(5*1024*1024*1024).
    int64_t  buf_size = 5368709120LL;
    //For debugging uncomment the following line
    //int64_t  buf_size = 1024LL;
    
    /* Number of variables, currently is 9 like Generic IO. */
    int num_vars  = 9;
    int rest_num =0;
    MPI_Offset  mpi_off = 0;
    MPI_Status  mpi_stat;
    char* writedata = NULL;

    int64_t buf_size_per_proc = 0;
    double mpiio_stime =0;
    double mpiio_etime = 0;
    double total_time = 0;
    double Max_total_time = 0;
    double Min_total_time = 0;
    double Sum_total_time = 0;
    double rate = 0;


    MPI_Init(&ac, &av);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    /* Buffer size per process per variable, The smallest unit for one process to write */
    buf_size_per_proc = buf_size/mpi_size;

    if (mpi_rank==0){
	    printf("Testing simple C MPIO program with %d processes accessing file %s\n",
	    mpi_size, filename);
        printf(" This tests the MPIO different patterns for 9 variables with MPIO write.\n");
        printf("There are four patterns. -c 9 contiguous writes -i 9 interleaved writes -p 3 writes -t 1 write.\n");
    }


    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, filename,
	    MPI_MODE_RDWR | MPI_MODE_CREATE ,
	    MPI_INFO_NULL, &fh))
	    != MPI_SUCCESS){
	    MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
	    PRINTID;
	    printf("MPI_File_open failed (%s)\n", mpi_err_str);
	    MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Allocate total amount of data per process to the buffer
    writedata = malloc(buf_size*num_vars/mpi_size);

    /* each process writes some data */
    if(ac >1) {
    
        if(strcmp(av[1],"-c")==0) {
            if(mpi_rank == 0) 
                printf("coming to contiguous pattern\n"); 
            mpi_off = buf_size_per_proc*mpi_rank;
            mpiio_stime = MPI_Wtime(); 
            for (i=0; i < num_vars; i++) {
                if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, buf_size_per_proc, MPI_BYTE,
  	                &mpi_stat))
  	                != MPI_SUCCESS){
  	                 MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
  	                 PRINTID;
  	                 printf("MPI_File_write_at offset(%ld), bytes (%d), failed (%s)\n",
  		             (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
  	                 MPI_Abort(MPI_COMM_WORLD, 1);
                }
                mpi_off+=buf_size;
      
           }
           mpiio_etime = MPI_Wtime();
           total_time = mpiio_etime - mpiio_stime;
        }
        else if(strcmp(av[1],"-i")==0) {
            if(mpi_rank == 0) 
                printf("Coming to the interleaved pattern.\n");

            // Each process has a contiuous write.
            mpi_off = mpi_rank*buf_size_per_proc*num_vars; 
            mpiio_stime = MPI_Wtime();
            for (i=0; i < num_vars; i++) {
                if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, buf_size_per_proc, MPI_BYTE,
                    &mpi_stat))
                    != MPI_SUCCESS){
                    MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
                    PRINTID;
                    printf("MPI_File_write_at offset(%ld), bytes (%d), failed (%s)\n",
                            (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                mpi_off+=buf_size_per_proc;
            }
            mpiio_etime = MPI_Wtime();
            total_time = mpiio_etime - mpiio_stime;
        }
    
        else if(strcmp(av[1],"-t")==0) {
            buf_size_per_proc = buf_size*num_vars/mpi_size;
            mpi_off = mpi_rank *buf_size_per_proc;
            if(mpi_rank == 0) 
                printf("Coming to one MPI-IO write with the atomic datatype.\n");

            // Each process has a contiuous write.
            mpiio_stime = MPI_Wtime();
            if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, buf_size_per_proc, MPI_BYTE,
                 &mpi_stat))
                 != MPI_SUCCESS){
                MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
                PRINTID;
                printf("MPI_File_write_at offset(%ld), bytes (%d), failed (%s)\n",
                       (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            mpiio_etime = MPI_Wtime();
            total_time = mpiio_etime - mpiio_stime;
        }
  
        else if(strcmp(av[1],"-p")==0) {
            // Leave the 2 variables out, handle the rest variables.
            mpiio_stime = MPI_Wtime();
            if(num_vars >2) {
                buf_size_per_proc = buf_size*(num_vars-2)/mpi_size;
                mpi_off = mpi_rank *buf_size_per_proc;
                if(mpi_rank == 0) 
                    printf("Coming to three  MPI-IO writes with the atomic datatype.\n");
                // Each process has a contiuous write.
                if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, buf_size_per_proc, MPI_BYTE,
                    &mpi_stat))
                    != MPI_SUCCESS){
                    MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
                    PRINTID;
                    printf("MPI_File_write_at offset(%ld), bytes (%d), failed (%s)\n",
                    (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }
       
            // This is just in case the total number of variables is <=2
            if(num_vars <=2) 
                rest_num = num_vars;
            else
                rest_num = 2;
            for (i=0; i < rest_num; i++) {
                buf_size_per_proc = buf_size/mpi_size;
                if(rest_num == 2) 
                    mpi_off = buf_size*(num_vars-2)+mpi_rank*buf_size_per_proc;
                else
                    mpi_off = mpi_rank*buf_size_per_proc;
        		if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, buf_size_per_proc, MPI_BYTE,
          			&mpi_stat))
          			!= MPI_SUCCESS){
           			MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
           			PRINTID;
           			printf("MPI_File_write_at offset(%ld), bytes (%d), failed (%s)\n",
          			(long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
          			MPI_Abort(MPI_COMM_WORLD, 1);
        		}	
       			 mpi_off+=buf_size;
     		}	
    
     		mpiio_etime = MPI_Wtime();
     		total_time = mpiio_etime - mpiio_stime;
    	}
    
  
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
  
  	}
  
  	free(writedata);

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
