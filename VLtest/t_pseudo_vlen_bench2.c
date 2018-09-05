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
#define DATASET2         "DS2"
#define DATASET_INDX    "DS_INDX"
//#define core 67108864

int
main (int argc, char *argv[] )
{
    FILE * pFile;

    hid_t       file, space, dset, dcpl,mem_dataspace;
                                    
    hid_t       space_indx, dset_indx;/* Handles */
    herr_t      status;
    hsize_t     dims[1];
    hsize_t     dims2D;
    int         *wdata,
                *rdata;
    hsize_t     *rdata_indx;
    hsize_t     i, j, k, m, n;
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
    hsize_t maxdims[1];
    hsize_t extdims[1],
      chunk[1] = {16};
    
    bool vlvl = true;

    hsize_t start[1],
      count[1];

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

    if( !(nvl_len = malloc((dims[0])*sizeof(hsize_t)) ) ) {
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
      H5Pset_fclose_degree(plist_id, H5F_CLOSE_WEAK);
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

      hsize_t k = 0;
      int increase = 1;
      int dir = 0;
      dims2D = 0;
      extdims[0]=1;
      int number = 0;
      int forend = 0;
      int jj=0;
      w = 0.;
      for (j=0; j < nelm; j++) {

	count[0] = (j/NVL%2?(NVL-1)-j%NVL:j%NVL+1);
	if(count[0] == 0) {
	  j = j + 1;
	  continue;
	}
	nelm_indx += 1;
	k = count[0];

	start[0] = dims2D;

	if(dims2D == nelm) {
	  nelm_indx -= 1;
	  break;
	} else if(dims2D + k > nelm) {
	  nvl_len[jj] = nelm - dims2D;
	  forend = 1;
	  //  break;
	} else {
	  dims2D = dims2D + k;
	  nvl_len[jj] = k;
	}
	count[0] = nvl_len[jj];
	//	printf("nelm_indx, nvl_len[j], dims2D, %lld %lld %lld\n", nelm_indx, nvl_len[j], dims2D);

	//	dims2D = 0;
	// 	for (m=0; m < nelm_indx ; m++) {
	//	  dims2D += nvl_len[j];
	  //	}


	dcpl = H5Pcreate (H5P_DATASET_CREATE);
	status = H5Pset_chunk (dcpl, 1, chunk);

	if(j == 0) {
	  maxdims[0] = H5S_UNLIMITED;
	  space = H5Screate_simple (1, &nvl_len[jj], maxdims);
	  dset  = H5Dcreate (file, DATASET2, H5T_NATIVE_INT, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
	} else {
	  extdims[0] +=  count[0];
	  //  printf("extdims %lld \n", extdims[0]);
	  dset = H5Dopen (file, DATASET2, H5P_DEFAULT);
	  status = H5Dset_extent (dset, extdims);
	  space = H5Dget_space (dset);
	}
	//	count[0] = nvl_len[j];

	status = H5Sselect_hyperslab (space, H5S_SELECT_SET, start, NULL, count, NULL);

	//	printf("start,count %lld %lld \n",start[0], count[0]);
	wdata = (int *)malloc(nvl_len[jj]*sizeof(int));
	
	if(wdata) {
	    for (m = 0; m < nvl_len[jj]; m++) {
	       *(wdata + m) = nvl_len[jj]-m;
	      //  *(wdata + m) = 99-m;
	      //  printf("aa %lld %lld \n",i,NVL-m);
	    }
	} else {
	  printf("wdata malloc failed \n");
	  abort();
	}
	//	for (m = 0; m < nvl_len[jj]; m++)
	//  printf("%d \n", wdata[m]);

	mem_dataspace = H5Screate_simple (1, count, NULL);

	gettimeofday(&tic, NULL);
	status = H5Dwrite (dset, H5T_NATIVE_INT, mem_dataspace, space, H5P_DEFAULT, &wdata[0]);
	if(status < 0) printf("H5Dwrite FAILED \n");

	status = H5Dclose (dset);
	status = H5Sclose (space);
	status = H5Sclose (mem_dataspace);
	status = H5Pclose (dcpl);
	free(wdata);

	gettimeofday(&toc, NULL);

	w += (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);

	if(forend == 1) break;
	jj +=1;
      }

      //printf("%lld %lld \n", nelm_indx, dims[0]);
 /*      if(nelm_indx>dims[0]) { */
/* 	printf("nv length does not match\n"); */
/* 	abort(); */
/*       } */

/*       printf("nct %lld %lld %lld \n", nelm_indx, NROWS*NVL, dims2D); */

      dims2D = 0;
      for (j=0; j < nelm_indx ; j++) {
	dims2D += nvl_len[j];
      }
 

      //  printf("nct %lld %lld %lld \n", nelm_indx, NROWS*NVL, dims2D);
      if(dims2D != NROWS*NVL){
	printf("Failed VL sizes %lld %lld\n",dims2D,NROWS*NVL);
	abort();
      }

      space_indx  = H5Screate_simple (1, &nelm_indx, NULL);
      dset_indx   = H5Dcreate (file, DATASET_INDX, H5T_NATIVE_INT, space_indx, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      //   DSsize = H5Dget_storage_size(dset);
      /*
       * Write the data to the dataset.
       */
      gettimeofday(&tic, NULL);

      status = H5Dwrite (dset_indx, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nvl_len[0]);
      status = H5Dclose (dset_indx);
      status = H5Sclose (space_indx);
      free(nvl_len);
      status = H5Fclose (file);
      
      gettimeofday(&toc, NULL);
      
      w += (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);
      H5Pclose(plist_id);
      H5Pclose(fcpl);
    } /* end write */

    r = 0.;
    if(read) {
      /*
       * Open file and dataset.
       */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef core
      H5Pset_fapl_core(plist_id, core, 1);
#endif
      file = H5Fopen (FILENAME, H5F_ACC_RDONLY, plist_id);

      dset = H5Dopen (file, DATASET_INDX, H5P_DEFAULT);

      /*
       * Get dataspace and allocate memory for array of vlen structures.
       * This does not actually allocate memory for the vlen data, that
       * will be done by the library.
       */
      space = H5Dget_space (dset);
      int ndims = H5Sget_simple_extent_dims (space, dims, NULL);
      rdata_indx = (hsize_t *)malloc(dims[0]*sizeof(hsize_t));
      /*
       * Read the data.
       */
      gettimeofday(&tic, NULL);
      status = H5Dread (dset, H5T_NATIVE_HSIZE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rdata_indx[0]);
      gettimeofday(&toc, NULL);
      r += (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);
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
	  inum = dims2D - rdata_indx[i]; 
	}

	printf (" {");
	//	printf("is %ld %ld", istart, iend);
	for (j=istart; j<istart + inum; j++) {
	  printf(" %d ",rdata[j]);
	}
	printf (" }\n");
      }
#endif

      status = H5Dclose (dset);
      status = H5Sclose (space);

      dset = H5Dopen (file, DATASET2, H5P_DEFAULT);

      /*
       * Get dataspace and allocate memory for array of vlen structures.
       * This does not actually allocate memory for the vlen data, that
       * will be done by the library.
       */
      space = H5Dget_space (dset);
      ndims = H5Sget_simple_extent_dims (space, &dims2D, NULL);

      
      start[0] = 0;
      for (j=0; j<dims[0] ; j++) {

	rdata = (int *)malloc(rdata_indx[j]*sizeof(int));
	mem_dataspace = H5Screate_simple (1, &rdata_indx[j], NULL);
	count[0] = rdata_indx[j];
	if (j != 0) {
	    start[0] += rdata_indx[j-1];
	}
	//	printf("start,count %lld %lld \n",start[0], count[0]);
	status = H5Sselect_hyperslab (space, H5S_SELECT_SET, start, NULL, count, NULL);

      /*
       * Read the data.
       */
	gettimeofday(&tic, NULL);
	status = H5Dread (dset, H5T_NATIVE_INT, mem_dataspace, space, H5P_DEFAULT, &rdata[0]);
	gettimeofday(&toc, NULL);

	r += (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);

#if DEBUG
	printf (" {");
	for (i=0; i<count[0]; i++) {
	  printf(" %d ",rdata[i]);
	}
	printf (" }\n");
#endif
	status = H5Sclose (mem_dataspace);
	free (rdata);
      }

      /*
       * Output the variable-length data to the screen.
       */
#if DEBUG
      printf("Total %ld Bytes, VL read time = %f seconds\n",DSsize, r);
#endif

      status = H5Dclose (dset);
      status = H5Sclose (space);

      status = H5Fclose (file);
      status = H5Pclose (plist_id);
      
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

    }

    struct stat st;
    stat(FILENAME, &st);
    DSsize = st.st_size/1048576;

    //DSsize = (dims2D*sizeof(int)+dims[0]*sizeof(hsize_t))/1048576;
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
