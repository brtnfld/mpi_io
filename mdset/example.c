/* 
 *  This example writes eight datasets to the HDF5 file by rows.
 */

#include "hdf5.h"
#include "stdlib.h"
#include "string.h"
#include "timer.h"

#define H5FILE_NAME "SDS.h5"
#define DSET_NAME   "data"
#define StringBool(x) ((x) ? "True" : "False")

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

int
main (int argc, char **argv)
{
    /*
     * HDF5 APIs definitions
     */ 	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[1];                 /* dataset dimensions */
    hid_t *dsets_ids;
    hid_t *mem_type_ids;
    hid_t *mem_space_ids;
    hid_t *file_space_ids;
    hid_t dxpl_id;
    const void **ptr;
    void **ptr_r;
    double *wbufs;
    double *rbufs;
    hsize_t	count[1];	          /* hyperslab selection parameters */
    hsize_t	offset[1];
    hid_t	plist_id;                 /* property list identifier */
    size_t      i,j;
    char        dset_name[14];
    double      write_time, read_time;
    int write, read, collective;

    hsize_t ndsets=1000;
    hsize_t nelem = 1024;

    char filename[80];

    FILE *fptr=NULL;

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

    write = 0;
    read  = 0;
    collective = 0;
    for (i = 1; i < argc; i++) {
      if(strcmp(argv[i],"-w")==0) {
        write=1;
      } else if(strcmp(argv[i],"-r")==0) {
        read=1;
      } else if(strcmp(argv[i],"-c")==0) {
        collective=1;
      } else if(strcmp(argv[i],"-d")==0) {
        i++;
        sscanf((argv)[i],"%zu",&ndsets);
      } else if(strcmp(argv[i],"-n")==0) {
        i++;
        sscanf((argv)[i],"%zu",&nelem);
      }
    }
                       
    if( write == 0 && read == 0) {
      write=1;
      read=1;
    }

    if(mpi_rank == 0) {
      printf(MAG);
      printf("SUMMARY\n-------\n");
      printf(" WRITE: %s\n", StringBool(write));
      printf(" READ: %s\n", StringBool(read));
      printf(" COLLECTIVE: %s\n", StringBool(collective));
      printf(" NUMEL: %ld\n",nelem);
      printf(" NUMDSETS: %ld\n\n",ndsets);
      printf(RESET);

      sprintf(filename, "ex_np_%d_nd_%zu_c_%d.dat", mpi_size, ndsets, collective);

      fptr = fopen(filename,"a");
      // fprintf(fptr,"# WRITE (MD,SD) , READ (MD,SD) \n",);

    }
    
    if( nelem%mpi_size != 0 ) {
      if(mpi_rank == 0) printf("The program assumes nelem is divisible by the number of ranks, stopping...\n");
      MPI_Abort(comm,1);
    }
    dsets_ids = (hid_t *) malloc(sizeof(hid_t)*ndsets);
    mem_type_ids = (hid_t *) malloc(sizeof(hid_t)*ndsets);
    mem_space_ids = (hid_t *) malloc(sizeof(hid_t)*ndsets);
    file_space_ids = (hid_t *) malloc(sizeof(hid_t)*ndsets);

    dimsf[0] = nelem;

    /*
     *   _   _   _   _   _
     *  / \ / \ / \ / \ / \
     * ( W | R | I | T | E )
     *  \_/ \_/ \_/ \_/ \_/
     *
     */
    
    if(write ==1) {
      wbufs = (double *) malloc(sizeof(double)*nelem);

      /*
       * Set up file access property list with parallel I/O access
       */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
      hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);
      /* OPTIMIZE */
      H5Pset_all_coll_metadata_ops(plist_id, 1);
      H5Pset_coll_metadata_write(plist_id, 1);
      H5Pset_alignment(plist_id, 0, 16777216);

      H5Pset_fapl_mpio(plist_id, comm, info);
        
      /*
       * Create a new file collectively and release property list identifier.
       */
      file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
      H5Pclose(plist_id);


      count[0] =  dimsf[0]/mpi_size;
      offset[0] = (mpi_rank)*count[0];

      for (i=0; i < count[0]; i++) {
        wbufs[i] = (double)(offset[0]+i);
      }

      ptr = malloc(sizeof(double*) * ndsets);

      size_t j=0;
      for (i=0; i < ndsets; i++) {
        snprintf(dset_name, 14, "%s_%08zu", DSET_NAME, i);

        file_space_ids[i] = H5Screate_simple(1, dimsf, NULL);
        mem_type_ids[i] = H5T_NATIVE_DOUBLE;
        dsets_ids[i] = H5Dcreate(file_id, dset_name, mem_type_ids[i], file_space_ids[i], H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

        mem_space_ids[i] = H5Screate_simple(1, count, NULL);

        H5Sselect_hyperslab(file_space_ids[i], H5S_SELECT_SET, offset, NULL, count, NULL);

        ptr[i] = &wbufs[0];

      }

      dxpl_id = H5Pcreate(H5P_DATASET_XFER);
      if(collective == 1)
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);

      timer_tick(&write_time, comm, 1);
      H5Dwrite_multi(ndsets, dsets_ids, mem_type_ids, mem_space_ids, file_space_ids, dxpl_id, ptr);
      timer_tock(&write_time);
      timer_collectprintstats(write_time, comm, 0, "Time for H5Dwrite_multi to complete:", fptr, mpi_size, ndsets, collective);

      timer_tick(&write_time, comm, 1);
      for (i=0; i < ndsets; i++) {
        H5Dwrite(dsets_ids[i], mem_type_ids[i], mem_space_ids[i], file_space_ids[i], dxpl_id, ptr[i]);
      }
      timer_tock(&write_time);
      timer_collectprintstats(write_time, comm, 0, "Time for H5Dwrite to complete:", fptr, mpi_size, ndsets, collective);


      for (i=0; i < ndsets; i++) {
        H5Sclose(file_space_ids[i]);
        H5Dclose(dsets_ids[i]);
        H5Sclose(mem_space_ids[i]);
      }

      /*
       * Close/release resources.
       */
      H5Pclose(dxpl_id);
      H5Fclose(file_id);

      free(ptr);
      free(wbufs);

    }
    else {
      if(mpi_rank == 0) {
        if(fptr != NULL){
          fprintf(fptr,"%d %zu %d %s ",mpi_size, ndsets, collective, "*");
        }
      }
    }

    /*
     *   _   _   _   _
     *  / \ / \ / \ / \
     * ( R | E | A | D )
     *  \_/ \_/ \_/ \_/
     */

    if( read == 1) {

      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(plist_id, comm, info);

      /* OPTIMIZE */
      H5Pset_all_coll_metadata_ops(plist_id, 1);
      H5Pset_coll_metadata_write(plist_id, 1);

      /*
       * Create a new file collectively and release property list identifier.
       */
      file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, plist_id);
      H5Pclose(plist_id);

      count[0] =  dimsf[0]/mpi_size;
      offset[0] = (mpi_rank)*count[0];

      rbufs = (double *) malloc(sizeof(double)*nelem);

      ptr_r = malloc(sizeof(double*) * ndsets);

      size_t j=0;
      for (i=0; i < ndsets; i++) {
        snprintf(dset_name, 14, "%s_%08zu", DSET_NAME, i);

        file_space_ids[i] = H5Screate_simple(1, dimsf, NULL);
        mem_type_ids[i] = H5T_NATIVE_DOUBLE;
        dsets_ids[i] = H5Dopen(file_id, dset_name, H5P_DEFAULT);

        mem_space_ids[i] = H5Screate_simple(1, count, NULL);

        H5Sselect_hyperslab(file_space_ids[i], H5S_SELECT_SET, offset, NULL, count, NULL);

        ptr_r[i] = &rbufs[0];

      }

      dxpl_id = H5Pcreate(H5P_DATASET_XFER);
      if(collective == 1)
        H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE);

      timer_tick(&read_time, comm, 1);
      H5Dread_multi(ndsets, dsets_ids, mem_type_ids, mem_space_ids, file_space_ids, dxpl_id, ptr_r);
      timer_tock(&read_time);
      timer_collectprintstats(read_time, comm, 0, "Time for H5Dread_multi to complete:", fptr, mpi_size, ndsets, collective);

      
      timer_tick(&read_time, comm, 1);
      for (i=0; i < ndsets; i++) {
        H5Dread(dsets_ids[i], mem_type_ids[i], mem_space_ids[i], file_space_ids[i], dxpl_id, ptr_r[i]);
      }
      timer_tock(&read_time);
      timer_collectprintstats(read_time, comm, 0, "Time for H5Dread to complete:", fptr, mpi_size, ndsets, collective);


      for (i=0; i < ndsets; i++) {
        H5Sclose(file_space_ids[i]);
        H5Dclose(dsets_ids[i]);
        H5Sclose(mem_space_ids[i]);
      }

      H5Pclose(dxpl_id);
      H5Fclose(file_id);

#if 0
      if(mpi_rank == 1) {
        for (i=0; i < ndsets; i++) {
          for (j=0; j < count[0]; j++) {
            printf("%lf\n",*(double *)ptr_r[i]);
            ptr_r[i] = ptr_r[i] + sizeof(double);
          }
        }
      }
#endif

      free(ptr_r);
      free(rbufs);

    } else {
      if(mpi_rank == 0) {
        if(fptr != NULL){
          fprintf(fptr,"%d %zu %d %s ",mpi_size, ndsets, collective, "*");
        }
      }
    }

    if(mpi_rank == 0) {
      if(fptr != NULL){
        fprintf(fptr,"\n");
        fclose(fptr);
      }

    }


    free(dsets_ids);
    free(mem_type_ids);
    free(mem_space_ids);
    free(file_space_ids);

    MPI_Finalize();

    return 0;
}    
