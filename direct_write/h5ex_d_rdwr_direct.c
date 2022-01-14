/************************************************************

  This example shows how to read and write data to a
  dataset.  The program first writes integers to a dataset
  with dataspace dimensions of DIM0xDIM1, then closes the
  file.  Next, it reopens the file, reads back the data, and
  outputs it to the screen.

  This file is intended for use with HDF5 Library version 1.8

 ************************************************************/

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>

#define FILE            "h5ex_d_rdwr.h5"
#define DATASET         "DS1"
#define DIM0            4096
#define DIM1            4096

int
main (void)
{
  hid_t       file, space, dset, dcpl_id;          /* Handles */
    herr_t      status;
    hsize_t     dims[2] = {DIM0, DIM1};
    int         wdata[DIM0][DIM1],          /* Write buffer */
                rdata[DIM0][DIM1],          /* Read buffer */
                i, j;
  
    /*
     * Initialize data.
     */
    for (i=0; i<DIM0; i++)
        for (j=0; j<DIM1; j++)
            wdata[i][j] = i * j - j;

    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple (2, dims, NULL);

    /*
     * Create the dataset.  We will use all default properties for this
     * example.
     */

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

    dset = H5Dcreate (file, DATASET, H5T_STD_I32LE, space, H5P_DEFAULT,
               dcpl_id, H5P_DEFAULT);
    
    off_t offset = (off_t)H5Dget_offset(dset);

    // status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata[0]);

    /*
     * Close and release resources.
     */
    status = H5Dclose (dset);
    status = H5Pclose (dcpl_id);
    status = H5Sclose (space);
    status = H5Fclose (file);

    /*
     * Open file and dataset using the default properties.
     */


    int fd;
    fd=open(FILE, O_WRONLY);

    size_t nbytes = sizeof(int)*dims[0]*dims[1];

    clock_t t;
    t = clock();
    pwrite(fd, wdata, nbytes, offset);
    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
    printf("pwrite took %f seconds to execute \n", time_taken);

    close(fd);

    /*
     * Now we begin the read section of this example.
     */

    /*
     * Open file and dataset using the default properties.
     */
    file = H5Fopen (FILE, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen (file, DATASET, H5P_DEFAULT);

    /*
     * Read the data using the default properties.
     */
    status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                rdata[0]);

    /*
     * Check the data
     */
    for (i=0; i<DIM0; i++)
        for (j=0; j<DIM1; j++)
          if(rdata[i][j] != wdata[i][j]) {
            printf (" error in writing/reading data");
            return -1;
          }

    /*
     * Close and release resources.
     */
    status = H5Dclose (dset);
    status = H5Fclose (file);

    return 0;
}
