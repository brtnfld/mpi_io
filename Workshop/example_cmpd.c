/* 
 *  This example writes eight datasets to the HDF5 file by rows.
 */

#include "hdf5.h"
#include "stdlib.h"
#include "string.h"
#include "timer.h"

#define H5FILE_NAME "SDS_cmpd.h5"
#define DSET_NAME   "data"
#define NDSETS 9
#define StringBool(x) ((x) ? "True" : "False")

int
main (int argc, char **argv)
{
    /*
     * HDF5 APIs definitions
     */ 	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[1];                 /* dataset dimensions */
    hsize_t	count[1];	          /* hyperslab selection parameters */
    hsize_t	offset[1];
    hid_t	plist_id;                 /* property list identifier */
    int         i;
    char        dset_name[12];
    double      write_time, read_time;
    hid_t       filetype;
    
    hsize_t TotSize = 16777216;

    typedef struct { 
        int data1;
        int data2;
        int data3;
        int data4;
        int data5;
        int data6;
        int data7;
        int data8;
        int data9;
    } data_t;

    data_t *data;
    int write, read, collective;

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
      } else if(strcmp(argv[i],"-n")==0) {
        i++;
        TotSize=atoi(argv[i]);
      } 
    }
                       
    if( write == 0 && read == 0) {
      write=1;
      read=1;
    }

    if(mpi_rank == 0) {
      printf("SUMMARY\n-------\n");
      printf(" WRITE: %s\n", StringBool(write));
      printf(" READ: %s\n", StringBool(read));
      printf(" COLLECTIVE: %s\n", StringBool(collective));
      printf(" NUMEL: %ld\n",TotSize);
    }
    
    if( TotSize%mpi_size != 0 ) {
      if(mpi_rank == 0) printf("The program assumes TotSize is divisible by the number of ranks, stopping...\n");
      MPI_Abort(comm,1);
    }

    dimsf[0] = TotSize;

    filetype = H5Tcreate (H5T_COMPOUND, sizeof(data_t));
    
    size_t offst=0;
    for (i=0; i < NDSETS; i++) {
      snprintf(dset_name, 12, "%s_%06d", DSET_NAME, i);
      H5Tinsert (filetype, dset_name, offst, H5T_NATIVE_INT);
      offst += sizeof(int);
    }
    /*
     *   _   _   _   _   _
     *  / \ / \ / \ / \ / \
     * ( W | R | I | T | E )
     *  \_/ \_/ \_/ \_/ \_/
     *
     */
    if(write ==1) {
      if(mpi_rank == 0) {
        /*
         * Set up file access property list with parallel I/O access
         */
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_EARLY);
        H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);
        
        /*
         * Create a new file collectively and release property list identifier.
         */
        file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);
        
        /*
         * Create the dataspace for the dataset.
         */
        filespace = H5Screate_simple(1, dimsf, NULL);
        
        /*
         * Create the dataset with default properties and close filespace.
         */
        dset_id = H5Dcreate(file_id, DSET_NAME, filetype, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        
        H5Dclose(dset_id);
        H5Sclose(filespace);
        H5Fclose(file_id);
        H5Pclose(dcpl_id);
      }
      
      MPI_Barrier(comm);
      
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(plist_id, comm, info);

      /*
       * Create a new file collectively and release property list identifier.
       */
      file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, plist_id);
      H5Pclose(plist_id);

      /*
       * Each process defines dataset in memory and writes it to the hyperslab
       * in the file.
       */
      count[0]  = dimsf[0]/mpi_size;
      offset[0] = mpi_rank*count[0];
      memspace  = H5Screate_simple(1, count, NULL);

      /*
       * Initialize data buffer
       */
      data = (data_t *) malloc(sizeof(data_t)*count[0]);
      for (i=0; i < count[0]; i++) {
        data[i].data1 = mpi_rank+10;
        data[i].data2 = mpi_rank+10;
        data[i].data3 = mpi_rank+10;
        data[i].data4 = mpi_rank+10;
        data[i].data5 = mpi_rank+10;
        data[i].data6 = mpi_rank+10;
        data[i].data7 = mpi_rank+10;
        data[i].data8 = mpi_rank+10;
        data[i].data9 = mpi_rank+10;
      }

      // printf("proc %d: count NX = %ld \n", mpi_rank, count[0]);
      // printf("proc %d: offset NX = %ld \n", mpi_rank, offset[0]);

      dset_id = H5Dopen(file_id, DSET_NAME, H5P_DEFAULT);
      filespace = H5Dget_space(dset_id);
      
      /*
       * Select hyperslab in the file.
       */
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
      
      timer_tick(&write_time, comm, 1);
      
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      
      if(collective == 1)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      
      H5Dwrite(dset_id, filetype, memspace, filespace, plist_id, data);
      
      H5Pclose(plist_id);
      H5Dclose(dset_id);
      
      timer_tock(&write_time);
      timer_collectprintstats(write_time, comm, 0, "Time for H5Dwrite to complete");
      
      /*
       * Close/release resources.
       */
      //   H5Tclose(filetype);
      H5Sclose(memspace);
      H5Sclose(filespace);
      H5Fclose(file_id);
      free(data);
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
      
      /*
       * Create a new file collectively and release property list identifier.
       */
      file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, plist_id);
      H5Pclose(plist_id);
      /*
       * Each process defines dataset in memory and writes it to the hyperslab
       * in the file.
       */
      count[0] =  dimsf[0]/mpi_size; 
      offset[0] = (mpi_rank)*count[0];
      memspace = H5Screate_simple(1, count, NULL);

      dset_id = H5Dopen(file_id, DSET_NAME, H5P_DEFAULT);
      filespace = H5Dget_space (dset_id);

      data = (data_t *) malloc(sizeof(data_t)*count[0]);

      /*
       * Select hyperslab in the file.
       */
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

      timer_tick(&read_time, comm, 1);

      plist_id = H5Pcreate(H5P_DATASET_XFER);

      if(collective == 1)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

      H5Dread(dset_id, filetype, memspace, filespace, plist_id, data);

      H5Pclose(plist_id);
      H5Dclose(dset_id);

      timer_tock(&read_time);
      timer_collectprintstats(read_time, comm, 0, "Time for H5Dread  to complete");

      /*
       * Close/release resources.
       */
      H5Tclose(filetype);
      H5Sclose(memspace);
      H5Sclose(filespace);
      H5Fclose(file_id);
      free(data);
    }

    MPI_Finalize();

    return 0;
}    
