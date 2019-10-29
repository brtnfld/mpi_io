/************************************************************

  This benchmark if for reading and writing variable-length
  datatypes to a dataset.  The program first writes
  variable-length integer arrays to a dataset then closes
  the file.  Next, it reopens the file, reads back the data, and
  closes the file. The VL is random for each array element.

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
#define FILENAME        "t_vlen_rand.h5"
#define DATASET_VL      "DSVL"
#define MiB 1048576

int
main (int argc, char *argv[] )
{
    FILE * pFile;

    hid_t       file, filetype, memtype, space, dset;
                                    /* Handles */
    herr_t      status;
    hsize_t     dims_r;
    int         *ptr;
    hvl_t       *wdataVL,           /* Array of vlen structures */
                *rdataVL;             /* Pointer to vlen structures */
    hsize_t     i, j;
    int opt;
    double w=0., r=0.;
    hsize_t DSsize;
    hsize_t dims_w =1048576; //4096;
    hsize_t VLmax = 1024; //4096;
    struct timeval  tic, toc;
    hsize_t vl_size;
    hid_t   plist_id, fcpl;
    int write=0, read=0;
    unsigned seed = 5;
    int fsm = 0;
    hsize_t fs_page_size = 4;
    size_t buf_page_size = 16;


    while ((opt = getopt(argc, argv, "rwhs:f:p:b:n:v:")) != -1) {
        switch (opt) {
        case 'r': read = 1; break;
        case 'w': write = 1; break;
        case 's':
          seed = atoi(optarg);
          break;
        case 'n':
          dims_w = atoi(optarg);
          break;
        case 'v':
          VLmax = atoi(optarg);
          break;
        case 'f': 
          fsm = atoi(optarg);
          break;
        case 'p': 
          fs_page_size = atoi(optarg);
          break;
        case 'b':
          buf_page_size = atoi(optarg);
          break;
        case 'h':
          printf("OPTIONS:\n");
          printf("   -r           read file [default - no]\n");
          printf("   -w           write file [default - no]\n");
          printf("   -n <int>     number of array elements [default = %lld]\n", dims_w);
          printf("   -v <int>     maximum variable length [default = %lld]\n", VLmax);
          printf("   -s <seed>    seed for random number generator [default = %d]\n", seed);
          printf("   -f <int>     free space manager: \n");
          printf("                      0 - FSM, Aggregators [default] \n");
          printf("                      1 - Paged FSM\n");
          printf("                      2 - Aggregators (no FSM)\n");
          printf("                      3 - (no FSM)\n");
          printf("   -p <int>      paged buffering (-f 1) option: \n");
          printf("                      int - file space page size (in MiB)  [default = %lld] \n", fs_page_size);
          printf("   -b <int>      paged buffering (-f 1) option: \n");
          printf("                      int - buffer size of the page (in MiB) [default = %ld] \n", buf_page_size);
          printf("   -h            help\n");
          return 0;
        case '?':
          if ( (optopt == 's'))  //  | (optopt == 'f')
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          return 1;
        default:
            fprintf(stderr, "Usage: %s [-rws]\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    printf("\n OPTION SUMMARY:\n  write=%d\n  read=%d\n  nelem=%lld\n  vlmax=%lld\n  seed=%d\n  fsm=%d\n ", 
           write, read, dims_w, VLmax, seed, fsm);

    if(fsm == 1) {
       printf("  file space page size = %lld MiB \n", fs_page_size);
       printf("   buffer size of the page = %ld MiB \n", buf_page_size);
    }

    if(write == 0 && read == 0) return 0;

    printf("\n %s (NROWS, MAX(VL)) = (%lld,%lld)\n\n", DATASET_VL, dims_w, VLmax);

    srand(seed);

    /*
     * Initialize variable-length data.
     */
    if(write == 1) {
      /*
       * Create a new file using the default properties.
       */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fclose_degree(plist_id, H5F_CLOSE_WEAK);

      H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

      fcpl = H5Pcreate(H5P_FILE_CREATE);

      /*
       *typedef enum H5F_fspace_strategy_t {
       *           H5F_FSPACE_STRATEGY_FSM_AGGR = 0, FSM, Aggregators, VFD 
       *           H5F_FSPACE_STRATEGY_PAGE = 1      Paged FSM, VFD 
       *           H5F_FSPACE_STRATEGY_AGGR = 2      Aggregators, VFD 
       *           H5F_FSPACE_STRATEGY_NONE = 3,     VFD 
       *           H5F_FSPACE_STRATEGY_NTYPES
       *     } H5F_fspace_strategy_t; 
       */

      H5Pset_file_space_strategy(fcpl, fsm, 0, (hsize_t)1);

      if(fsm == 1) {
        H5Pset_file_space_page_size(fcpl, fs_page_size*(hsize_t)MiB);
        H5Pset_page_buffer_size(plist_id, (size_t)(buf_page_size*(hsize_t)MiB), 0, 0);
      }

      vl_size = 0;
      for (j=0; j < dims_w; j++) {
	vl_size += (rand() % (VLmax)) + 1;
      }

      if( !(wdataVL = malloc (vl_size * sizeof (hvl_t)) ) ) {
	printf("malloc wdataVL failed \n");
	abort();
      }

      for (j=0; j < dims_w; j++) {
        wdataVL[j].len = (rand() % (VLmax)) + 1;
      }

      for (i = 0; i <  dims_w; i++){
	ptr = (int *) malloc (wdataVL[i].len * sizeof (int));
	for (j = 0; j <  wdataVL[i].len; j++) {
	  ptr[j] =  wdataVL[i].len-(size_t)j;
	}
	wdataVL[i].p = (void *) ptr;
      }

      gettimeofday(&tic, NULL);

      file = H5Fcreate(FILENAME, H5F_ACC_TRUNC, fcpl, plist_id);

      filetype = H5Tvlen_create (H5T_NATIVE_INT);
      memtype = H5Tvlen_create (H5T_NATIVE_INT);

      space = H5Screate_simple (1, &dims_w, NULL);

      dset = H5Dcreate (file, DATASET_VL, filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdataVL);
      if(status != 0) {
        fprintf (stderr, "H5Dwrite Failed. \n");
        return status;
      }

      status = H5Dvlen_reclaim (memtype, space, H5P_DEFAULT, wdataVL);

      status = H5Dclose (dset);
      status = H5Sclose (space);
      status = H5Tclose (filetype);
      status = H5Tclose (memtype);

      status = H5Fclose (file);

      gettimeofday(&toc, NULL);
      w = (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);

      H5Pclose(plist_id);
      H5Pclose(fcpl);

      free(wdataVL);
    }

    if(read == 1) {
      /*
       * Open file and dataset.
       */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fclose_degree(plist_id, H5F_CLOSE_WEAK);

      gettimeofday(&tic, NULL);

      if( (file = H5Fopen (FILENAME, H5F_ACC_RDONLY, plist_id)) < 0 ) {
        printf("error opening file: %s\n", FILENAME);
        return 1;
      }

      status = H5Pclose (plist_id);
      dset = H5Dopen (file, DATASET_VL, H5P_DEFAULT);

      /*
       * Get dataspace and allocate memory for array of vlen structures.
       * This does not actually allocate memory for the vlen data, that
       * will be done by the library.
       */
      space = H5Dget_space (dset);
      H5Sget_simple_extent_dims (space, &dims_r, NULL);
      rdataVL = (hvl_t *) malloc (dims_r * sizeof (hvl_t));

      /*
       * Create the memory datatype.
       */
      memtype = H5Tvlen_create (H5T_NATIVE_INT);
      /*
       * Read the data.
       */
      status = H5Dread (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdataVL);
      if(status != 0) {
        fprintf (stderr, "H5Dread Failed. \n");
        return status;
      }

#if DEBUG
      /*
       * Output the variable-length data to the screen.
       */
      for (i=0; i<dims_r; i++) {
        printf ("%s[%lld]: {",DATASET_VL,i);
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

      status = H5Dclose (dset);
      status = H5Sclose (space);
      status = H5Tclose (memtype);

      status = H5Fclose (file);

      gettimeofday(&toc, NULL);
      r = r + (double) (toc.tv_usec - tic.tv_usec) / 1000000 + (double) (toc.tv_sec - tic.tv_sec);

      free (rdataVL);

    }

    struct stat st;
    stat(FILENAME, &st);
    DSsize = st.st_size/MiB;

    printf("Total time %lld MiB, %f s \n", DSsize,  w+r);

    pFile = fopen ("time_VLrandom.txt", "a");
    if(write){
      printf(" -- Write -- Total %lld MiB, %f MiB/s \n",DSsize,DSsize/w);
      fprintf(pFile, "%f ", w);
    }
    if(read){
      printf(" -- Read  -- Total %lld MiB, %f MiB/s \n",DSsize,DSsize/r);
      fprintf(pFile, "%f ", r);
    }
    fclose(pFile);
    return 0;

}