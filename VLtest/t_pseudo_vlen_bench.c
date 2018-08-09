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

    hid_t       file, space, dset;
                                    /* Handles */
    herr_t      status;
    hsize_t     dims[1];
    hsize_t     dims2D[1];
    int         *wdata,
                *rdata;
    hsize_t     *rdata_indx;
    hsize_t     i, j;
    int opt, cnt=0;
    double w=0., r=0.;
    hsize_t DSsize;
    hsize_t NROWS = 4096; //4096;
    hsize_t NVL = 4096; //4096;
    struct timeval  tic, toc;
    hid_t   plist_id, fcpl;
    int write=0,read=0;

    hsize_t *nvl_len;
    hsize_t nelm;
    hsize_t nelm_indx = 0;
    
    bool vlvl = true;

    while ((opt = getopt(argc, argv, "rwv")) != -1) {
        cnt=cnt+1;
        switch (opt) {
        case 'r': read = 1; break;
        case 'w': write = 1; break;
        default:
            fprintf(stderr, "Usage: %s [-rw]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    printf("options: w=%d r=%d \n", write, read);
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

    nelm = NROWS*NVL;

    if( !(nvl_len = malloc(dims[0]*sizeof(hsize_t)) ) ) {
      printf("malloc nvl_len failed \n");
      abort();
    }

    printf("VL_2D(NROWS,NVL*) = (%lld,%lld)\n", dims[0], NVL);

    if(write==1) {
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
      hsize_t k = 0;
      int increase = 1;
      dims2D[0] = 0;
      for (j=0; j < nelm; j++) {
	nelm_indx += 1;
	if(k == NVL)
	  increase = 0;
	if(increase)
	  k += 1;
	else
	  k -= 1;
	if(k < 1 && increase == 0){
	  k = 2;
	  increase = 1;
	}
	if(dims2D[0] == nelm) {
	  nelm_indx -= 1;
	  break;
	} else if(dims2D[0] + k > nelm) {
	  nvl_len[j] = nelm - dims2D[0];
	  break;
	} else {
	  dims2D[0] = dims2D[0] + k;
	  nvl_len[j] = k;
	}
	//	printf("%lld %lld %lld\n", nelm_indx, nvl_len[j], dims2D[0]);
      }

      printf("%lld %lld \n", nelm_indx, dims[0]);
      if(nelm_indx>dims[0]) {
	printf("nv length does not match\n");
	abort();
      }

      dims2D[0] = 0;
      for (j=0; j < nelm_indx ; j++) {
	dims2D[0] += nvl_len[j];
      }
      
      //  printf("nct %lld %lld %lld \n", nelm_indx, NROWS*NVL, dims2D[0]);
      if(dims2D[0] != NROWS*NVL){
	printf("Failed VL sizes %lld %lld\n",dims2D[0],NROWS*NVL);
	abort();
      }

      space = H5Screate_simple (1, dims2D, NULL);

      dset   = H5Dcreate (file, DATASET, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose (dset);
      status = H5Sclose (space);

      space  = H5Screate_simple (1, &nelm_indx, NULL);
      dset   = H5Dcreate (file, DATASET_INDX, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose (dset);
      status = H5Sclose (space);

      status = H5Fclose (file);
      
      status = H5Pclose (plist_id);
      status = H5Pclose (fcpl);
      //   printf("nct %lld \n", nelm_indx);
    }
    // printf("nct %lld \n", nelm_indx);
 
    if(write) {
      wdata = (int *)malloc(dims2D[0]*sizeof(int));
      if(wdata) {
	hsize_t icnt2 = 0;
	for (i = 0; i <  nelm_indx; i++){
	  NVL = nvl_len[i];
	  nvl_len[i] = icnt2;
	  //printf("aa %lld \n", nvl_len[i]);
	  for (j = 0; j < NVL; j++) {
	    *(wdata + nvl_len[i] + j) = NVL-j;
	    // printf("aa %lldd %lldd \n",i,NVL-j);
	  }
	  icnt2 = icnt2 + NVL;
	}
      } else {
	printf("wdata malloc failed \n");
	abort();
      }
      //  printf("here \n");
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
      status = H5Dwrite (dset, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nvl_len[0]);
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
      int ndims = H5Sget_simple_extent_dims (space, dims2D, NULL);

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
#endif

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
      rdata_indx = (hsize_t *)malloc(dims[0]*sizeof(hsize_t));
      /*
       * Read the data.
       */
      gettimeofday(&tic, NULL);
      status = H5Dread (dset, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rdata_indx[0]);
      
      /*
       * Output the variable-length data to the screen.
       */
#if DEBUG
      printf("Total %lld Bytes, VL read time = %f seconds\n",DSsize, r);
      for (i=0; i<dims[0]; i++) {
	hsize_t istart = rdata_indx[i];
	hsize_t inum;
	if(i != dims[0]-1) {
	  inum =  rdata_indx[i+1]-rdata_indx[i];
	} else {
	  inum = dims2D[0] - rdata_indx[i]; 
	}

	printf (" {");
	//	printf("is %ld %ld", istart, iend);
	for (j=istart; j<istart + inum; j++) {
	  printf(" %d ",rdata[j]);
	}
	printf (" }\n");
      }
#endif
      
      /*
       * Close and release resources.  Note we must still free the
       * top-level pointer "rdataVL", as H5Dvlen_reclaim only frees the
       * actual variable-length data, and not the structures themselves.
       */
      free (rdata_indx);
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

    DSsize = (dims2D[0]*sizeof(int)+dims[0]*sizeof(hsize_t))/1048576;
    //printf("Total time %ld MB, %f %f  MB/s \n", DSsize,  w, r);

    pFile = fopen ("VL_timing.txt", "a");
    if(write){
      printf("Total %lld MB,%f MB/s \n",DSsize,DSsize/w);
      fprintf(pFile, "%f ", w);
    } else if(read){
      printf("Total %lld MB,%f MB/s \n",DSsize,DSsize/r);
      fprintf(pFile, "%f ", r);
    }
    fclose(pFile);
    return 0;
}
