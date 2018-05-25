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
#include <unistd.h>
#include <stdbool.h>

#define DEBUG 0
#define FILENAME        "h5ex_t_vlen.h5"
#define DATASET_VL      "DSVL"
#define DATASET         "DS"
#define DATASET_INDX    "DS_INDX"
//#define core 67108864

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
    int opt, cnt=0;
    double w_vl=0., w=0., r_vl=0., r=0.;
    hsize_t DSsize;
    hsize_t NROWS = 4096;
    hsize_t NVL = 4096;
    struct timeval  tic, toc;
    hid_t   plist_id, fcpl;
    int write=0,read=0,vl=0;

    hsize_t *indx;
    int *nvl_len;
    
    bool vlvl = true;

    while ((opt = getopt(argc, argv, "rwv")) != -1) {
        cnt=cnt+1;
        switch (opt) {
        case 'r': read = 1; break;
        case 'w': write = 1; break;
        case 'v': vl = 1; break;
        default:
            fprintf(stderr, "Usage: %s [-rwv] [vl]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    printf("options: w=%d r=%d vl=%d \n", write, read, vl);
    if (argc > 2) {
      NROWS = strtoimax(argv[2], NULL, 10);
    }


    /*
     * Initialize variable-length data.  wdataVL[0] is a countdown of
     * length NVL
     */

    dims[0] = NROWS;
    if(vlvl)
      dims[0] = 2*NROWS;

    indx = malloc(dims[0]*sizeof(hsize_t));
    nvl_len = malloc(dims[0]*sizeof(hsize_t));

    printf("VL_2D(NROWS,NVL*) = (%ld,%ld)\n", dims[0], NVL);


    if( write==1 && !(vl == 1)) {
    /*
     * Create a new file using the default properties.
     */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef core
      H5Pset_fapl_core(plist_id, core, 1);
#endif
      H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

      fcpl = H5Pcreate(H5P_FILE_CREATE);

/* typedef enum H5F_fspace_strategy_t { */
/*           H5F_FSPACE_STRATEGY_FSM_AGGR = 0, /\* FSM, Aggregators, VFD *\/  */
/*           H5F_FSPACE_STRATEGY_PAGE = 1	   /\* Paged FSM, VFD *\/ */
/*           H5F_FSPACE_STRATEGY_AGGR = 2	   /\* Aggregators, VFD *\/ */
/*           H5F_FSPACE_STRATEGY_NONE = 3,     /\* VFD *\/ */
/*           H5F_FSPACE_STRATEGY_NTYPES      */
/*     } H5F_fspace_strategy_t; */
#if 0
      H5Pset_file_space_strategy(fcpl,H5F_FSPACE_STRATEGY_PAGE,0,(hsize_t)1);
      H5Pset_file_space_page_size(fcpl, (hsize_t)(1048576));
      
      H5Pset_page_buffer_size(plist_id, (size_t)(8*1048576), 0, 0);
#endif

      file = H5Fcreate (FILENAME, H5F_ACC_TRUNC, fcpl, plist_id);

      int nd = dims[0]/2/NVL;
      hsize_t k = 1;
      hsize_t cnt = 0;
      for (j=0; j< dims[0]; j++) {
	cnt = cnt + k;
	nvl_len[j] = k;
	if(j <= dims[0]/2-1) {
	  if( j != dims[0]/2-1) k++;
	} else {
	  k--;
	}
      }

      dims2D[0] = cnt;

      space = H5Screate_simple (1, dims2D, NULL);

      dset   = H5Dcreate (file, DATASET, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose (dset);
      status = H5Sclose (space);

      space = H5Screate_simple (1, dims, NULL);
      dset   = H5Dcreate (file, DATASET_INDX, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose (dset);
      status = H5Sclose (space);

      status = H5Fclose (file);
      
      status = H5Pclose (plist_id);
      status = H5Pclose (fcpl);
    }
 
    if(write) {
      wdata = (int *)malloc(dims2D[0]*sizeof(int));
      int icnt2 = 0;
      for (i = 0; i <  dims[0]; i++){
	NVL = nvl_len[i];
	nvl_len[i] = icnt2;
	for (j = 0; j < NVL; j++) {
	  *(wdata + nvl_len[i] + j) = NVL-j;
	}
	icnt2 = icnt2 + NVL;
      }
      
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef core
      H5Pset_fapl_core(plist_id, core, 1);
#endif
      file = H5Fopen(FILENAME, H5F_ACC_RDWR, plist_id);
      dset = H5Dopen(file, DATASET, H5P_DEFAULT);
      
      DSsize = H5Dget_storage_size(dset);
      /*
       * Write the data to the dataset.
       */
      gettimeofday(&tic, NULL);
      
      status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &wdata[0]);
      status = H5Dclose (dset);
      free(wdata);

      dset = H5Dopen(file, DATASET_INDX, H5P_DEFAULT);
      status = H5Dwrite (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nvl_len[0]);
      status = H5Dclose (dset);
      free(nvl_len);

      status = H5Fclose (file);
      
      gettimeofday(&toc, NULL);
      
      w = (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);
      H5Pclose(plist_id);
    }
    if(read) {
      /*
       * Open file and dataset.
       */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef core
      H5Pset_fapl_core(plist_id, core, 1);
#endif
      file = H5Fopen (FILENAME, H5F_ACC_RDONLY, plist_id);
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

      dset = H5Dopen (file, DATASET_INDX, H5P_DEFAULT);

      /*
       * Get dataspace and allocate memory for array of vlen structures.
       * This does not actually allocate memory for the vlen data, that
       * will be done by the library.
       */
      space = H5Dget_space (dset);
      ndims = H5Sget_simple_extent_dims (space, dims, NULL);
      rdata = (int *)malloc(dims[0]*sizeof(int));
      /*
       * Read the data.
       */
      gettimeofday(&tic, NULL);
      status = H5Dread (dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rdata[0]);
      
      /*
       * Output the variable-length data to the screen.
       */
#if DEBUG
      printf("Total %ld Bytes, VL read time = %f seconds\n",DSsize, r);
      printf (" {\n");
      for (i=0; i<dims[0]; i++) {
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
      gettimeofday(&toc, NULL);
      r = r + (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);
      status = H5Pclose (plist_id);

    }

    DSsize = (dims2D[0]+dims[0])*sizeof(int)/1048576;
    printf("Total %ld MB, %f %f %f %f MB/s \n",DSsize, -99., DSsize/w, -99.,DSsize/r);
    //printf("Total time %ld MB, %f %f %f %f MB/s \n",DSsize, w_vl, w, r_vl, r);

    pFile = fopen ("VL_timing.txt", "a");
    if(write && !vl){
      fprintf(pFile, "%f ", w);
    } else if(read && !vl){
      fprintf(pFile, "%f ", r);
    } else if(write && vl){
      fprintf(pFile, "%f ", w_vl);
    } else if(read && vl){
      fprintf(pFile, "%f \n", r_vl);
    }
    fclose(pFile);
    return 0;
}
