/************************************************************

  This example shows how to read and write data to a dataset
  using gzip compression (also called zlib or deflate).  The
  program first checks if gzip compression is available,
  then if it is it writes integers to a dataset using gzip,
  then closes the file.  Next, it reopens the file, reads
  back the data, and outputs the type of compression and the
  maximum value in the dataset to the screen.

  This file is intended for use with HDF5 Library version 1.8

 ************************************************************/

#include "hdf5.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define FILE            "h5ex_d_compress.h5"
#define DATASET         "DS1"

/**** PROGRAM PARAMETERS *****/
#define DIM0            8388608 /* data buffer size 0-dir */
#define DIM1            2048 /* data buffer size 1-dir */
#define CHUNK0          16384 /* chunk size in 0-dir = DIM0/CHUNK0 */
#define CHUNK1          1 /* chunk size in 1-dir = DIM1/CHUNK1 */
/**** PROGRAM PARAMETERS *****/

#define RANK            2

#define DEBUG 0

int
main (int argc, char **argv)
{
    hid_t           file, dset, dcpl;    /* Handles */
     hid_t          filespace, memspace;      /* file and memory dataspace identifiers */
    herr_t          status;
    htri_t          avail;
    hsize_t         chunk[2];
    int             *wdata; /* pointer to data buffer to write */
    int             *rdata; 
    hsize_t         dimsf[2];                 /* dataset dimensions */
    hsize_t         i;
    int             max;
    hsize_t	    count[2];	          /* hyperslab selection parameters */
    hsize_t	    offset[2];
    hid_t	    plist_id;                 /* property list identifier */
    unsigned int    filter_info;

    /*
     * MPI variables
     */
    int mpi_size, mpi_rank;
    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;

    /*
     * Initialize MPI
     */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    /* 
     * Set up file access property list with parallel I/O access
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    H5Pset_fclose_degree(plist_id,H5F_CLOSE_WEAK);

    /*
     * Check if gzip compression is available and can be used for both
     * compression and decompression.  Normally we do not perform error
     * checking in these examples for the sake of clarity, but in this
     * case we will make an exception because this filter is an
     * optional part of the hdf5 library.
     */
    avail = H5Zfilter_avail(H5Z_FILTER_DEFLATE);
    if (!avail) {
        printf ("gzip filter not available.\n");
        return 1;
    }
    status = H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info);
    if ( !(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) ||
                !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED) ) {
        printf ("gzip filter not available for encoding and decoding.\n");
        return 1;
    }

    if(mpi_rank == 0)
      printf("chunk size (MB) = %d\n ", DIM0/CHUNK0*DIM1/CHUNK1*sizeof(int)/1048576);
    
    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
   

    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = DIM0;
    dimsf[1] = DIM1;
    filespace = H5Screate_simple(RANK, dimsf, NULL); 

    /*
     * Create the dataset creation property list, add the gzip
     * compression filter and set the chunk size.
     */
    
    chunk[0] = DIM0/CHUNK0;
    chunk[1] = DIM1/CHUNK1;

    dcpl = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_layout(dcpl, H5D_CHUNKED);
    H5Pset_fill_time(dcpl, H5D_FILL_TIME_NEVER);
    status = H5Pset_deflate (dcpl, 9);
    status = H5Pset_chunk (dcpl, 2, chunk);

    /*
     * Create the dataset and close filespace.
     */
    dset = H5Dcreate (file, DATASET, H5T_STD_I32LE, filespace, H5P_DEFAULT, dcpl,
                H5P_DEFAULT);
    status = H5Sclose (filespace);
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    memspace = H5Screate_simple(RANK, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Initialize data.
     */
    wdata = (int *) malloc(sizeof(int)*count[0]*count[1]);
    for (i=0; i < count[0]*count[1]; i++) {
        wdata[i] = mpi_rank + 1;
    }

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    /*
     * Write the data to the dataset.
     */

    status = H5Dwrite (dset, H5T_NATIVE_INT, memspace, filespace, plist_id, wdata);

    free(wdata);
    H5Pclose(plist_id);

    /*
     * Close and release resources.
     */
    status = H5Pclose (dcpl);
    status = H5Dclose (dset);
    status = H5Fclose (file);
    /*
     * Now we begin the read section of this example.
     */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    /*
     * Open file and dataset using the default properties.
     */
    file = H5Fopen (FILE, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);


    dset = H5Dopen (file, DATASET, H5P_DEFAULT);

    /*
     * Retrieve dataset creation property list.
     */
    dcpl = H5Dget_create_plist (dset);

#if 0 // DEBUG
    /*
     * Retrieve and print the filter type.  Here we only retrieve the
     * first filter because we know that we only added one filter.
     */
    size_t nelmts = 0;
    unsigned int flags;
    H5Z_filter_t filter_type = H5Pget_filter (dcpl, 0, &flags, &nelmts, NULL, 0, NULL,
                &filter_info);

    if(mpi_rank == 0) {
      printf ("Filter type is: ");
      switch (filter_type) {
      case H5Z_FILTER_DEFLATE:
	printf ("H5Z_FILTER_DEFLATE\n");
	break;
      case H5Z_FILTER_SHUFFLE:
	printf ("H5Z_FILTER_SHUFFLE\n");
	break;
      case H5Z_FILTER_FLETCHER32:
	printf ("H5Z_FILTER_FLETCHER32\n");
	break;
      case H5Z_FILTER_SZIP:
	printf ("H5Z_FILTER_SZIP\n");
	break;
      case H5Z_FILTER_NBIT:
	printf ("H5Z_FILTER_NBIT\n");
	break;
      case H5Z_FILTER_SCALEOFFSET:
	printf ("H5Z_FILTER_SCALEOFFSET\n");
      }
    }
#endif
    /*
     * Initialize data.
     */
    rdata = (int *) malloc(sizeof(int)*count[0]*count[1]);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

    /*
     * Read the data using the default properties.
     */
    status = H5Dread (dset, H5T_NATIVE_INT, memspace, filespace, plist_id, rdata);
#if DEBUG
      hsize_t icnt = 0;
      hsize_t j;
      for (i=0; i<count[0]; i++) {
	printf("[ ");
	for (j=0; j<count[1]; j++) {
	  printf("%d ",rdata[icnt]);
	  icnt +=1;
	}
	printf("]\n");
      }
#endif
    /*
     * Find the maximum value in the dataset, to verify that it was
     * read correctly.
     */
    max = rdata[0];
    for (i=1; i < count[0]*count[1]; i++)
      if (max < rdata[i])
	max = rdata[i];
    /*
     * Print the maximum value.
     */
    if(mpi_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, &max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD);
      printf ("Maximum value in %s is: %d\n", DATASET, max);
    } else {
      MPI_Reduce(&max, &max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    status = H5Pclose (plist_id);
    /*
     * Close and release resources.
     */
    status = H5Pclose (dcpl);
    status = H5Dclose (dset);
    status = H5Fclose (file);
    free(rdata); // free allocated memory
    status = H5Sclose (memspace);
    status = H5Sclose (filespace);

    MPI_Finalize();
    
    return 0;
}
