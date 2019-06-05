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
#include <sys/stat.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#ifndef MPI_FILE_NULL           /*MPIO may be defined in mpi.h already       */
#   include <mpio.h>
#endif

#define PRINTID printf("Proc %d: ", mpi_rank)

#define CHCK_VAL 0

bool dequal(double a, double b, double epsilon)
{
 return fabs(a-b) < epsilon;
}

int main(int ac, char **av)
{
    int  mpi_size, mpi_rank;
    MPI_File fh;
    char *filename = "./data.mpi";
    char mpi_err_str[MPI_MAX_ERROR_STRING];
    int  mpi_err_strlen;
    int  mpi_err;
    double expect_val;
    int  i=0, j=0 ; 
    int  nerrors = 0;		/* number of errors */
    long int gb;
    /* buffer size is the total size for one variable. */
    /* The buffer size will be 8 GB, 72 GB total (9* 1073741824*8/(1024*1024*1024)). */

    // (multiples of 32)
    //int64_t buf_size = 1073741824LL;
    //int64_t buf_size = 4194304;
    // Summit (multiples of 42)
    //int64_t buf_size = 1213857792LL;
    //For debugging uncomment the following line
    int64_t  buf_size = 32768LL;
    //int64_t  buf_size = 1024;
    double dexpect_val;
    
    /* Number of variables, currently is 9 like Generic IO. */
    // 12
    int num_vars  = 9;
    int rest_num =0;
    MPI_Offset  mpi_off = 0;
    MPI_Status  mpi_stat;
    double* writedata = NULL;
    double* readdata = NULL;

    int64_t buf_size_per_proc = 0;
    double mpiio_stime =0;
    double mpiio_etime = 0;
    double total_time = 0;
    double Max_total_time = 0;
    double Min_total_time = 0;
    double Sum_total_time = 0;
    int    Sum_nerrors = 0;
    double rate_avg=0., rate_min=0., rate_max=0.;
    FILE *pFile;

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
      pFile = fopen("timing.txt", "a");
    }

    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, filename,
                                 MPI_MODE_CREATE|MPI_MODE_WRONLY,
                                 MPI_INFO_NULL, &fh))
        != MPI_SUCCESS){
      MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
      PRINTID;
      printf("MPI_File_open failed (%s)\n", mpi_err_str);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if(ac >1) {
        if(strcmp(av[1],"-t")==0) {
          num_vars = 1;
        }
    }


    // Allocate total amount of data per process to the buffer
    writedata = malloc(buf_size_per_proc*num_vars*sizeof(double));

    i = 0;
    for (j=0; j < buf_size_per_proc*num_vars; j++){
      if( j % (buf_size_per_proc) == 0) {
        i = 0.;
      } else {
        i+=1;
      }
      writedata[j] = (double)(mpi_rank*buf_size_per_proc + i);
      //   printf(" WRITE data[%d:%d] %lf\n", mpi_rank, j, writedata[j]);
    }


    /* each process writes some data */
    if(ac >1) {
        if(strcmp(av[1],"-c")==0) {
            if(mpi_rank == 0) 
                printf("coming to contiguous pattern\n"); 
            mpi_off = buf_size_per_proc*mpi_rank*sizeof(double);
            mpiio_stime = MPI_Wtime(); 
            for (i=0; i < num_vars; i++) {
                if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, buf_size_per_proc*sizeof(double), MPI_BYTE,
  	                &mpi_stat))
  	                != MPI_SUCCESS){
  	                 MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
  	                 PRINTID;
  	                 printf("MPI_File_write_at offset(%ld), bytes (%d), failed (%s)\n",
  		             (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
  	                 MPI_Abort(MPI_COMM_WORLD, 1);
                }
                //   printf(" offest %d %ld \n", mpi_rank, mpi_off);
                mpi_off+=buf_size*sizeof(double);
           }
           mpiio_etime = MPI_Wtime();
           total_time = mpiio_etime - mpiio_stime;
        }
        else if(strcmp(av[1],"-i")==0) {
            if(mpi_rank == 0) 
                printf("Coming to the interleaved pattern.\n");

            // Each process has a contiuous write.
            mpi_off = buf_size_per_proc*mpi_rank*num_vars*sizeof(double); 
            mpiio_stime = MPI_Wtime();
            for (i=0; i < num_vars; i++) {
                if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, buf_size_per_proc*sizeof(double), MPI_BYTE,
                    &mpi_stat))
                    != MPI_SUCCESS){
                    MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
                    PRINTID;
                    printf("MPI_File_write_at offset(%ld), bytes (%d), failed (%s)\n",
                            (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
                mpi_off+=buf_size_per_proc*sizeof(double);
            }
            mpiio_etime = MPI_Wtime();
            total_time = mpiio_etime - mpiio_stime;
        }
    
        else if(strcmp(av[1],"-t")==0) {
            mpi_off = buf_size_per_proc*mpi_rank*num_vars*sizeof(double);
            if(mpi_rank == 0)
                printf("Coming to one MPI-IO write with the atomic datatype.\n");

            // Each process has a contiuous write.
            mpiio_stime = MPI_Wtime();

            if ((mpi_err = MPI_File_write_at(fh, mpi_off, writedata, buf_size_per_proc*sizeof(double), MPI_BYTE,
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
          gb = (long int)(buf_size*sizeof(double)*num_vars)/(1024*1024*1024);
          rate_min = (double)(gb/Max_total_time);
          rate_max = (double)(gb/Min_total_time);
          rate_avg = (double)(gb/(Sum_total_time/mpi_size));
          printf("%d Procs WRITE %d variables, %ld GB \n",mpi_size, num_vars, gb);
          printf("WRITE time is avg, min, max: %lf %f %f s.\n", Sum_total_time/mpi_size, Min_total_time, Max_total_time);
          printf("WRITE Bandwidth is avg, min, max: %f %f %f GB/s.\n", rate_avg, rate_min, rate_max);
          fprintf(pFile, "%s %d %f %f %f\n", av[1], mpi_size, rate_avg, rate_min, rate_max);
        }	
  
    }
  
    free(writedata);

    MPI_File_close(&fh);

#if 0
    if(mpi_rank == 0) {
    FILE *ptr;
    ptr = fopen(filename,"rb");
    double buffer[32];

    fread(buffer,sizeof(buffer),1,ptr); // read 10 bytes to our buffer

    for(i = 0; i<32; i++)
      printf("%lf \n", buffer[i]); // prints a series of bytes

    }
    MPI_Barrier(MPI_COMM_WORLD);
    abort();
#endif
 

    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, filename,
                                 MPI_MODE_RDONLY,
                                 MPI_INFO_NULL, &fh))
        != MPI_SUCCESS){
      MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
      PRINTID;
      printf("MPI_File_open failed (%s)\n", mpi_err_str);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Allocate total amount of data per process to the buffer
    readdata = malloc(buf_size*num_vars/mpi_size*sizeof(double));
  

    if(ac >1) {

      if(strcmp(av[1],"-c")==0 || strcmp(av[1],"-t")==0) {

        if(mpi_rank == 0)
          printf("coming to contiguous pattern\n"); 
        mpi_off = buf_size_per_proc*mpi_rank*sizeof(double);
        mpiio_stime = MPI_Wtime();
        for (i=0; i < num_vars; i++) {
          if ((mpi_err = MPI_File_read_at(fh, mpi_off, readdata, buf_size_per_proc*sizeof(double), MPI_BYTE,
                                           &mpi_stat))
              != MPI_SUCCESS){
            MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
            PRINTID;
            printf("MPI_File_read_at offset(%ld), bytes (%d), failed (%s)\n",
                   (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
          mpi_off+=buf_size*sizeof(double);
#if CHCK_VAL
          for (j=0; j < buf_size_per_proc; j++){

	    dexpect_val = (double)(mpi_rank*buf_size_per_proc + j);
            
            if(!dequal(readdata[j], dexpect_val, 1.0e-6)) {
              PRINTID;
              printf("read data[%d:%d] got %f, expect %f\n", mpi_rank, j,
                     readdata[j], dexpect_val);
              nerrors++;
            }
          }
#endif    
        }
        mpiio_etime = MPI_Wtime();
        total_time = mpiio_etime - mpiio_stime;
      }
      else if(strcmp(av[1],"-i")==0) {
        if(mpi_rank == 0) 
          printf("Coming to the READ interleaved pattern.\n");
        
        // Each process has a contiuous write.
        mpi_off = buf_size_per_proc*mpi_rank*num_vars*sizeof(double); 
        mpiio_stime = MPI_Wtime();
        for (i=0; i < num_vars; i++) {
          if ((mpi_err = MPI_File_read_at(fh, mpi_off, readdata, buf_size_per_proc*sizeof(double), MPI_BYTE,
                                           &mpi_stat))
                    != MPI_SUCCESS){
            MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
            PRINTID;
            printf("MPI_File_read_at offset(%ld), bytes (%d), failed (%s)\n",
                   (long) mpi_off, (int) buf_size_per_proc, mpi_err_str);
            MPI_Abort(MPI_COMM_WORLD, 1);
          }
          mpi_off+=buf_size_per_proc*sizeof(double);
#if CHCK_VAL
          for (j=0; j < buf_size_per_proc; j++){

	    dexpect_val = (double)(mpi_rank*buf_size/mpi_size + j);

            if(!dequal(readdata[j], dexpect_val, 1.0e-6)) {
              PRINTID;
              printf("read data[%d:%d] got %f, expect %f\n", mpi_rank, j,
                     readdata[j], dexpect_val);
              nerrors++;
            }
          }
#endif
        }
        mpiio_etime = MPI_Wtime();
        total_time = mpiio_etime - mpiio_stime;
      }
 
      MPI_Reduce(&total_time, &Max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&total_time, &Sum_total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&total_time, &Min_total_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(&nerrors, &Sum_nerrors, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      if(mpi_rank == 0) {
        gb = (long int)(buf_size*sizeof(double)*num_vars)/(1024*1024*1024);
        rate_min = (double)(gb/Max_total_time);
        rate_max = (double)(gb/Min_total_time);
        rate_avg = (double)(gb/(Sum_total_time/mpi_size));
        printf("%d Procs READ %d variables, %ld GB \n",mpi_size,num_vars, gb);
        printf("READ time is avg, min, max: %lf %f %f s.\n", Sum_total_time/mpi_size, Min_total_time, Max_total_time);
        printf("READ Bandwidth is avg, min, max: %f %f %f GB/s.\n", rate_avg, rate_min, rate_max);
#if CHCK_VAL
        printf("NUMBER OF ERRORS = %d \n",Sum_nerrors);
#endif
        fprintf(pFile, "%s %d %f %f %f\n", av[1], mpi_size, rate_avg, rate_min, rate_max);
        fclose(pFile);
      }

    }

    free(readdata);
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

        i = 0;
        for (j=0; j < buf_size*num_vars/mpi_size; j++){
          writedata[j] = mpi_rank*buf_size/mpi_size + i;
          if(i%(buf_size/mpi_size-1) == 0)
            i = 0;
          else
            i=+1;
        }

	for (i=0; i < DIMSIZE; i++){
          
	    expect_val = irank*DIMSIZE + i;
	    if (readdata[i] != dexpect_val){
		PRINTID;
		printf("read data[%d:%d] got %d, expect %d\n", irank, i,
			readdata[i], dexpect_val);
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
