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
/* HDF5 header file */
#include "hdf5.h"

#define RANK 1

#define PRINTID printf("Proc %d: ", mpi_rank)

int main(int ac, char **av)
{
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
    // OLD int64_t  buf_size = 5368709120LL;
    // 1Gib = 1073741824
    // 9*1073741824 ( 9 GB total, 1GB per var.)
    
    int64_t buf_size = 9663676416LL;
    
    //int64_t buf_size = 36864LL;

    //For debugging uncomment the following line
    //int64_t  buf_size = 1024LL;
    
    /* Number of variables, currently is 9 like Generic IO. */
    int num_vars  = 9;
    int rest_num =0;
    MPI_Offset  mpi_off = 0;
    MPI_Status  mpi_stat;
    char* writedata = NULL;
    double *writedata2;

    int64_t buf_size_per_proc = 0;
    double mpiio_stime =0;
    double mpiio_etime = 0;
    double total_time = 0;
    double Max_total_time = 0;
    double Min_total_time = 0;
    double Sum_total_time = 0;
    double rate = 0;

    // HDF5
    hid_t file_id;              /* File ID */
    hid_t fapl_id;		/* File access property list */
    hid_t dset_id;		/* Dataset ID */
    hid_t file_space_id;	/* File dataspace ID */
    hid_t mem_space_id;		/* Memory dataspace ID */
    hsize_t file_dims[RANK];   	/* Dataset dimemsion sizes */
    hsize_t mem_dims[1];   	/* Memory buffer dimemsion sizes */
    hsize_t file_start[RANK];	/* File dataset selection start coordinates (for hyperslab setting) */
    hsize_t file_count[RANK];	/* File dataset selection count coordinates (for hyperslab setting) */
    hsize_t mem_start[1];	/* Memory buffer selection start coordinates (for hyperslab setting) */
    hsize_t mem_count[1];	/* Memory buffer selection count coordinates (for hyperslab setting) */
    int mpi_size, mpi_rank;	/* MPI variables */
    herr_t ret;         	/* Generic return value */
    size_t ii;
    int hdf5=1;

    const char * DATASETNAME[] = {
      "id",
      "mask",
      "x",
      "y",
      "z",
      "vx",
      "vy",
      "vz",
      "phi"
    };

    typedef struct {
      double id;
      double mask;
      double x;
      double y;
      double z;
      double vx;
      double vy;
      double vz;;
      double phi;
    } hacc_t;
    hacc_t *Hdata;
    hid_t Hmemtype;
    hid_t Hfiletype;

    
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


    if(hdf5) {

      if(strcmp(av[1],"-i")==0) {

          Hmemtype = H5Tcreate (H5T_COMPOUND, sizeof (hacc_t));
          H5Tinsert (Hmemtype, "id",
                     HOFFSET (hacc_t, id), H5T_NATIVE_DOUBLE);
          H5Tinsert (Hmemtype, "mask", 
                     HOFFSET (hacc_t, mask), H5T_NATIVE_DOUBLE);
          H5Tinsert (Hmemtype, "x",
                     HOFFSET (hacc_t, x), H5T_NATIVE_DOUBLE);
          H5Tinsert (Hmemtype, "y",
                     HOFFSET (hacc_t, y), H5T_NATIVE_DOUBLE);
          H5Tinsert (Hmemtype, "z",
                     HOFFSET (hacc_t, z), H5T_NATIVE_DOUBLE);
          H5Tinsert (Hmemtype, "vx",
                     HOFFSET (hacc_t, vx), H5T_NATIVE_DOUBLE);
          H5Tinsert (Hmemtype, "vy",
                     HOFFSET (hacc_t, vy), H5T_NATIVE_DOUBLE);
          H5Tinsert (Hmemtype, "vz",
                     HOFFSET (hacc_t, vz), H5T_NATIVE_DOUBLE);
          H5Tinsert (Hmemtype, "phi",
                     HOFFSET (hacc_t, phi), H5T_NATIVE_DOUBLE);

      }

      if( mpi_rank == 0) {
      
        file_dims[0] = buf_size/num_vars/sizeof(double);

        if(strcmp(av[1],"-c")==0) {

          for (i=0; i < num_vars; i++) {
            
            /* Create the file */
            file_id = H5Fcreate(DATASETNAME[i], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

            file_space_id = H5Screate_simple(1, file_dims, NULL);
            
            /* Create the dataset collectively */
            dset_id = H5Dcreate2(file_id, DATASETNAME[i], H5T_NATIVE_DOUBLE,
                                 file_space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            ret = H5Sclose(file_space_id);
            ret = H5Dclose(dset_id);
            ret = H5Fclose(file_id);
          }
        }
      }

      MPI_Barrier(MPI_COMM_WORLD);
      

    }
    else {

      if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, filename,
                                   MPI_MODE_RDWR | MPI_MODE_CREATE ,
                                   MPI_INFO_NULL, &fh))
          != MPI_SUCCESS){
        MPI_Error_string(mpi_err, mpi_err_str, &mpi_err_strlen);
        PRINTID;
        printf("MPI_File_open failed (%s)\n", mpi_err_str);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

    }

    // Allocate total amount of data per process to the buffer
    if(hdf5) {
      mem_dims[0] = buf_size/num_vars/mpi_size/sizeof(double);

      if(strcmp(av[1],"-i")==0) {

        Hdata = (hacc_t *) malloc ( mem_dims[0] * sizeof (hacc_t));
        for (ii=0; ii < mem_dims[0]; ii++) {
          Hdata[ii].id = (double)(ii+1+(double)mpi_rank/1000.);
          Hdata[ii].mask =  Hdata[ii].id;
          Hdata[ii].x =  Hdata[ii].id;
          Hdata[ii].y =  Hdata[ii].id;
          Hdata[ii].z =  Hdata[ii].id;
          Hdata[ii].vx =  Hdata[ii].id;
          Hdata[ii].vy =  Hdata[ii].id;
          Hdata[ii].vz =  Hdata[ii].id;
          Hdata[ii].phi =  Hdata[ii].id;
        }


      } else {

        writedata2 = malloc( mem_dims[0]*sizeof(hacc_t));
        for (ii=0; ii < mem_dims[0]; ii++) {
          writedata2[ii] = (double)(ii+1+(double)mpi_rank/1000.);
        }

      }
        file_dims[0] = mem_dims[0]*mpi_size;
        file_space_id = H5Screate_simple(1, file_dims, NULL);


    } else {
      writedata = malloc(buf_size*num_vars/mpi_size);
    }
    /* each process writes some data */

    if(ac >1) {
        if(strcmp(av[1],"-c")==0) {
            if(mpi_rank == 0) 
                printf("coming to contiguous pattern\n");

            mpiio_stime = MPI_Wtime();
            if(hdf5) {
               
              for (i=0; i < num_vars; i++) {

                /* Create an HDF5 file access property list */
                fapl_id = H5Pcreate (H5P_FILE_ACCESS);
                
                /* Set file access property list to use the MPI-IO file driver */
                ret = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
                
                /* Create the file collectively */
                file_id = H5Fopen(DATASETNAME[i], H5F_ACC_RDWR, fapl_id);
      
                /* Release file access property list */
                ret = H5Pclose(fapl_id);

                dset_id = H5Dopen(file_id, DATASETNAME[i], H5P_DEFAULT);

                /* Create memory dataspace for write buffer */
                
                mem_space_id = H5Screate_simple(1, mem_dims, NULL);
                
                /* Select column of elements in the file dataset */
                file_start[0] = mpi_rank*mem_dims[0];
                file_count[0] = mem_dims[0];
                ret = H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
                /* Select all elements in the memory buffer */
                mem_start[0] = 0;
                mem_count[0] = mem_dims[0];
                ret = H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET, mem_start, NULL, mem_count, NULL);
                
                /* Write data independently */
                ret = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, H5P_DEFAULT, writedata2);
                /* Close memory dataspace */
                ret = H5Sclose(mem_space_id);
                /* Close dataset collectively */
                ret = H5Dclose(dset_id);
                ret = H5Fclose(file_id);
              }
            } else {
              mpi_off = buf_size_per_proc*mpi_rank;
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
            }

           mpiio_etime = MPI_Wtime();
           total_time = mpiio_etime - mpiio_stime;
        }
        else if(strcmp(av[1],"-i")==0) {

            if(mpi_rank == 0) 
                printf("Coming to the interleaved pattern.\n");
 
            mpiio_stime = MPI_Wtime();
            if(hdf5) {

              dset_id = H5Dopen(file_id, "ALLVAR", H5P_DEFAULT);

              /* Create memory dataspace for write buffer */
                
              mem_space_id = H5Screate_simple(1, mem_dims, NULL);
                
              /* Select column of elements in the file dataset */
              file_start[0] = mpi_rank*mem_dims[0];
              file_count[0] = mem_dims[0];
              ret = H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
              /* Select all elements in the memory buffer */
              mem_start[0] = 0;
              mem_count[0] = mem_dims[0];
              ret = H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET, mem_start, NULL, mem_count, NULL);
                
              /* Write data independently */
              ret = H5Dwrite(dset_id, Hmemtype, mem_space_id, file_space_id, H5P_DEFAULT, Hdata);
              /* Close memory dataspace */
              ret = H5Sclose(mem_space_id);
              /* Close dataset collectively */
              ret = H5Dclose(dset_id);

            } else {
              
              // Each process has a contiuous write.
              mpi_off = mpi_rank*buf_size_per_proc*num_vars;
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
            }


            mpiio_etime = MPI_Wtime();
            total_time = mpiio_etime - mpiio_stime;
            
            ret = H5Tclose (Hmemtype);
            free(Hdata);
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



        if(hdf5) {
          /* Close file dataspace */
          ret = H5Sclose(file_space_id);
          // ret = H5Fclose(file_id);
        } else {
          MPI_File_close(&fh);
        }


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
