/************************************************************

  This example shows how to read and write variable-length
  datatypes to a dataset.  The program first writes two
  variable-length integer arrays to a dataset then closes
  the file.  Next, it reopens the file, reads back the data,
  and outputs it to the screen.

  This file is intended for use with HDF5 Library version 1.8

 ************************************************************/

#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define DEBUG 0
#define FILENAME        "h5ex_t_vlen.h5"
#define DATASET_VL      "DSVL"
#define DATASET         "DS"

int
main (int argc, char *argv[] )
{
    FILE * pFile;

    hid_t       file, filetype, memtype, space, dset;
                                    /* Handles */
    herr_t      status;
    hvl_t       *wdataVL,           /* Array of vlen structures */
                *rdataVL;             /* Pointer to vlen structures */
    hsize_t     dims[1];
    hsize_t     dims2D[1];
    int         *ptr,
                ndims;
    int         *wdata,
                *rdata;
    hsize_t     i, j;
    double w_vl, w, r_vl, r;
    hsize_t DSsize;
    hsize_t NROWS = 2048;
    hsize_t NVL = 4096;
    struct timeval  tic, toc;

    /*
     * Initialize variable-length data.  wdataVL[0] is a countdown of
     * length NVL
     */
    if (argc > 1) {
      NROWS = strtoimax(argv[1], NULL, 10);;
    }

    printf("(NROWS,NVL) = (%ld,%ld)\n",NROWS,NVL);

    dims[0] = NROWS;
    dims2D[0] =NROWS*NVL;

    wdataVL = malloc (NROWS * sizeof (hvl_t));

    for (j=0; j<NROWS; j++) {
      wdataVL[j].len = NVL;
      ptr = (int *) malloc (wdataVL[j].len * sizeof (int));
      for (i=0; i<wdataVL[j].len; i++)
        ptr[i] = wdataVL[j].len - (size_t)(i);       /* n-1 */
      wdataVL[j].p = (void *) ptr;
    }

    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate (FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create variable-length datatype for file and memory.
     */
    filetype = H5Tvlen_create (H5T_STD_I32LE);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple (1, dims, NULL);

    /*
     * Create the dataset and write the variable-length data to it.
     */
    dset = H5Dcreate (file, DATASET_VL, filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (filetype);

    space = H5Screate_simple (1, dims2D, NULL);

    dset = H5Dcreate (file, DATASET, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dclose (dset);
    status = H5Fclose (file);

#if 1
    file = H5Fopen(FILENAME, H5F_ACC_RDWR, H5P_DEFAULT);
    dset = H5Dopen(file, DATASET_VL, H5P_DEFAULT);
    memtype = H5Tvlen_create (H5T_NATIVE_INT);
    space = H5Screate_simple (1, dims, NULL);

    gettimeofday(&tic, NULL);
    status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdataVL);

    /*
     * Close and release resources.  Note the use of H5Dvlen_reclaim
     * removes the need to manually free() the previously malloc'ed
     * data.
     */
    status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, wdataVL);

    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (memtype);
    

    H5Fclose(file);

    gettimeofday(&toc, NULL);
    w_vl = (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);

#endif
 

    wdata = (int *)malloc(NROWS*NVL*sizeof(int));
    for (i = 0; i <  NROWS; i++)
      for (j = 0; j < NVL; j++)
	*(wdata + i*NVL + j) = NVL-j;

    file = H5Fopen(FILENAME, H5F_ACC_RDWR, H5P_DEFAULT);
    dset = H5Dopen(file, DATASET, H5P_DEFAULT);

    DSsize = H5Dget_storage_size(dset);
    /*
     * Write the data to the dataset.
     */
    gettimeofday(&tic, NULL);

    status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &wdata[0]);
    status = H5Dclose (dset);
    free(wdata);
    status = H5Fclose (file);

    gettimeofday(&toc, NULL);

    w = (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);

    /*
     * Now we begin the read section of this example.  Here we assume
     * the dataset has the same name and rank, but can have any size.
     * Therefore we must allocate a new array to read in data using
     * malloc().
     */

    /* Try to clear the cache */

    //  system("dd if=/dev/zero of=bigfile bs=1M count=100000");
    //system("cat bigfile > /dev/null");

    /*system("perl -e '$a = \"A\" x 1_000_000_000; sleep 5' &");*/
    //  system("find / -type f -exec cat {} >>/dev/null ;");
/*     FILE * fp; */
/*     /\* open the file for writing*\/ */
/*     fp = fopen ("scr","wb"); */
/*     fwrite(wdata, sizeof(*wdata),NROWS*NVL,fp); */
/*     fwrite(wdata, sizeof(*wdata),NROWS*NVL,fp); */
/*     fwrite(wdata, sizeof(*wdata),NROWS*NVL,fp); */
/*     fwrite(wdata, sizeof(*wdata),NROWS*NVL,fp); */
/*     system("cp scr scr2"); */
    // fread(wdata, sizeof(*wdata),NROWS*NVL,fp);
    //  printf("%d\n", *(wdata + 1*NVL + 1));
    
/*     fclose(fp); */


#if 1
    /*
     * Open file and dataset.
     */
    file = H5Fopen (FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    dset = H5Dopen (file, DATASET_VL, H5P_DEFAULT);

    /*
     * Get dataspace and allocate memory for array of vlen structures.
     * This does not actually allocate memory for the vlen data, that
     * will be done by the library.
     */
    space = H5Dget_space (dset);
    ndims = H5Sget_simple_extent_dims (space, dims, NULL);
    rdataVL = (hvl_t *) malloc (dims[0] * sizeof (hvl_t));

    /*
     * Create the memory datatype.
     */
    memtype = H5Tvlen_create (H5T_NATIVE_INT);

    /*
     * Read the data.
     */
    gettimeofday(&tic, NULL);
    status = H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdataVL);
    gettimeofday(&toc, NULL);
    r_vl = (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);
#if DEBUG
    printf("Total %ld MB, VL write time = %f seconds\n", DSsize, r_vl);
    /*
     * Output the variable-length data to the screen.
     */
    for (i=0; i<dims[0]; i++) {
        printf ("%s[%u]:\n  {",DATASET_VL,i);
        ptr = rdataVL[i].p;
        for (j=0; j<rdataVL[i].len; j++) {
            printf (" %d", ptr[j]);
            if ( (j+1) < rdataVL[i].len )
                printf (",");
        }
        printf (" }\n");
    }
#endif


    

    /*
     * Close and release resources.  Note we must still free the
     * top-level pointer "rdataVL", as H5Dvlen_reclaim only frees the
     * actual variable-length data, and not the structures themselves.
     */
    status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, rdataVL);
    free (rdataVL);
    
    DSsize = H5Dget_storage_size(dset);
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (memtype);


    dset = H5Dopen (file, DATASET, H5P_DEFAULT);

    /*
     * Get dataspace and allocate memory for array of vlen structures.
     * This does not actually allocate memory for the vlen data, that
     * will be done by the library.
     */
    space = H5Dget_space (dset);
    ndims = H5Sget_simple_extent_dims (space, dims2D, NULL);
    rdata = (int *)malloc(dims2D[0]*sizeof(int));

    /*
     * Read the data.
     */
    gettimeofday(&tic, NULL);
    status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rdata[0]);
    gettimeofday(&toc, NULL);
    r = (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);

    /*
     * Output the variable-length data to the screen.
     */
#if DEBUG
    printf("Total %ld Bytes, VL read time = %f seconds\n",DSsize, r);
    printf (" {\n");
    for (i=0; i<dims2D[0]; i++) {
      printf(" %d ",rdata[i]);
    }
    printf (" }\n");
#endif

    /*
     * Close and release resources.  Note we must still free the
     * top-level pointer "rdataVL", as H5Dvlen_reclaim only frees the
     * actual variable-length data, and not the structures themselves.
     */
    free (rdata);


    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Fclose (file);
#endif

    DSsize = NROWS*NVL*sizeof(int)/1048576;
    printf("Total %ld MB, %f %f %f %f MB/s \n",DSsize, DSsize/w_vl,DSsize/w,DSsize/r_vl,DSsize/r);
    // printf("Total %ld MB, %f %f %f %f MB/s \n",DSsize, w_vl, w, r_vl, r);

    pFile = fopen ("VL_timing.txt", "a");
    fprintf(pFile, "%f %f %f %f \n", w_vl,w,r_vl,r);
    fclose(pFile);
    return 0;
}
