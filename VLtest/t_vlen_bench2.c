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
#include <sys/stat.h>

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
    hsize_t     dims[1];
    hsize_t     dims2D;
    int         *ptr;
    hvl_t       *wdataVL,           /* Array of vlen structures */
                *rdataVL;             /* Pointer to vlen structures */
    hsize_t     i, j;
    int opt, cnt=0;
    double w=0., r=0., r_vl=0.;
    hsize_t DSsize;
    hsize_t NROWS =4096; //4096;
    hsize_t NVL = 4096; //4096;
    struct timeval  tic, toc;
    hid_t   plist_id, fcpl;
    int write=0,read=0;
    hsize_t k;

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

    printf("VL_2D(NROWS,NVL*) = (%lld,%lld)\n", dims[0], NVL);

    if(write==1) {
    /*
     * Create a new file using the default properties.
     */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fclose_degree(plist_id, H5F_CLOSE_WEAK);
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

      int increase = 1;
      dims2D = 0;
      k = 0;
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
	if(dims2D == nelm) {
	  nelm_indx -= 1;
	  break;
	} else if(dims2D + k > nelm) {
	  break;
	} else {
	  dims2D = dims2D + k;
	}
	//	printf("%lld %lld %lld\n", nelm_indx, nvl_len[j], dims2D);
      }
      
      if( !(wdataVL = malloc (nelm_indx * sizeof (hvl_t)) ) ) {
	printf("malloc nvl_len failed \n");
	abort();
      }

      increase = 1;
      dims2D = 0;
      k = 0;
      for (j=0; j < nelm_indx; j++) {
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
	if(dims2D == nelm) {
	  break;
	} else if(dims2D + k > nelm) {
	  wdataVL[j].len = nelm - dims2D;
	  break;
	} else {
	  dims2D = dims2D + k;
	  wdataVL[j].len = k;
	}
	//	printf("%lld %lld %lld\n", nelm_indx, nvl_len[j], dims2D);
      }

      // printf("%lld %lld \n", nelm_indx, dims[0]);
      if(nelm_indx>dims[0]) {
	printf("nv length does not match\n");
	abort();
      }

      dims2D = 0;
      for (j=0; j < nelm_indx ; j++) {
	dims2D += wdataVL[j].len;
      }

      for (i = 0; i <  nelm_indx; i++){
	ptr = (int *) malloc (wdataVL[i].len * sizeof (int));
	for (j = 0; j <  wdataVL[i].len; j++) {
	  ptr[j] =  wdataVL[i].len-(size_t)j;
	}
	wdataVL[i].p = (void *) ptr;
      }

      //  printf("nct %lld %lld %lld \n", nelm_indx, NROWS*NVL, dims2D);
      if(dims2D != NROWS*NVL){
	printf("Failed VL sizes %lld %lld\n",dims2D,NROWS*NVL);
	abort();
      }
      filetype = H5Tvlen_create (H5T_STD_I32LE);
      memtype = H5Tvlen_create (H5T_NATIVE_INT);

      space = H5Screate_simple (1, &nelm_indx, NULL);

      dset   = H5Dcreate (file, DATASET_VL, filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      gettimeofday(&tic, NULL);
      status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdataVL);

      status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, wdataVL);

      // DSsize = H5Dget_storage_size(dset);

      status = H5Dclose (dset);
      status = H5Sclose (space);
      status = H5Tclose (filetype);
      status = H5Tclose (memtype);

      status = H5Fclose (file);

      gettimeofday(&toc, NULL);
      w = (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);
      free(wdataVL);
      H5Pclose(plist_id);
      H5Pclose(fcpl);
    }

    if(read) {
      /*
       * Open file and dataset.
       */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fclose_degree(plist_id, H5F_CLOSE_WEAK);
#ifdef core
      H5Pset_fapl_core(plist_id, core, 1);
#endif
      file = H5Fopen (FILENAME, H5F_ACC_RDONLY, plist_id);
      dset = H5Dopen (file, DATASET_VL, H5P_DEFAULT);

      /*
       * Get dataspace and allocate memory for array of vlen structures.
       * This does not actually allocate memory for the vlen data, that
       * will be done by the library.
       */
      space = H5Dget_space (dset);
      H5Sget_simple_extent_dims (space, &dims2D, NULL);
      rdataVL = (hvl_t *) malloc (dims2D * sizeof (hvl_t));

      
      /*
       * Create the memory datatype.
       */
      memtype = H5Tvlen_create (H5T_NATIVE_INT);
      /*
       * Read the data.
       */
      gettimeofday(&tic, NULL);
      status = H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdataVL);
      //   gettimeofday(&toc, NULL);
      //r_vl = (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);
      
#if DEBUG
      printf("Total %lld MB, VL write time = %f seconds\n", DSsize, r_vl);
      /*
       * Output the variable-length data to the screen.
       */
      for (i=0; i<dims2D; i++) {
        printf ("%s[%lld]:\n  {",DATASET_VL,i);
        ptr = rdataVL[i].p;
        for (j=0; j<rdataVL[i].len; j++) {
	  printf (" %d", ptr[j]);
	  if ( (j+1) < rdataVL[i].len )
	    printf (",");
        }
        printf (" }\n");
      }
#endif

      status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, rdataVL);

      //  DSsize = H5Dget_storage_size(dset);
      status = H5Dclose (dset);
      status = H5Sclose (space);
      status = H5Tclose (memtype);

      status = H5Fclose (file);
      gettimeofday(&toc, NULL);
      r = r + (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);
      status = H5Pclose (plist_id);
      free (rdataVL);

    }
    struct stat st;
    stat(FILENAME, &st);
    DSsize = st.st_size/1048576;

    //  DSsize = (dims2D*sizeof(int))/1048576;
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
