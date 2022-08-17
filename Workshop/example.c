/* 
 *  This example writes eight datasets to the HDF5 file by rows.
 */

#include "hdf5.h"
#include "stdlib.h"
#include "string.h"
#include "timer.h"

#define H5FILE_NAME "SDS.h5"
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
    int         *data1,*data2,*data3,*data4,*data5,*data6,*data7,*data8,*data9; /* pointer to data buffer to write */
    hsize_t	count[1];	          /* hyperslab selection parameters */
    hsize_t	offset[1];
    hid_t	plist_id;                 /* property list identifier */
    int         i;
    char        dset_name[12];
    double      write_time, read_time;
    int write, read, collective;

    hsize_t TotSize = 16777216;

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
    if(mpi_rank == 0) {
      printf("SUMMARY\n-------\n");
      printf(" WRITE: %s\n", StringBool(write));
      printf(" READ: %s\n", StringBool(read));
      printf(" COLLECTIVE: %s\n", StringBool(collective));
      printf(" NUMEL: %ld\n",TotSize);
    }
                       
    if( write == 0 && read == 0) {
      write=1;
      read=1;
    }
    
    if( TotSize%mpi_size != 0 ) {
      if(mpi_rank == 0) printf("The program assumes TotSize is divisible by the number of ranks, stopping...\n");
      MPI_Abort(comm,1);
    }
    
    dimsf[0] = TotSize;
    /*
     *   _   _   _   _   _
     *  / \ / \ / \ / \ / \
     * ( W | R | I | T | E )
     *  \_/ \_/ \_/ \_/ \_/
     *
     */
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
      
      for (i=0; i < NDSETS; i++) {
        snprintf(dset_name, 12, "%s_%06d", DSET_NAME, i);
        /*
         * Create the dataset with default properties and close filespace.
         */
        dset_id = H5Dcreate(file_id, dset_name, H5T_NATIVE_INT, filespace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        H5Dclose(dset_id);
      }
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
     * Create the dataspace for the dataset.
     */
    filespace = H5Screate_simple(1, dimsf, NULL);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] =  dimsf[0]/mpi_size; 
    offset[0] = (mpi_rank)*count[0];
    memspace = H5Screate_simple(1, count, NULL);

    /*
     * Initialize data buffer
     */
    data1 = (int *) malloc(sizeof(int)*count[0]);
    data2 = (int *) malloc(sizeof(int)*count[0]);
    data3 = (int *) malloc(sizeof(int)*count[0]);
    data4 = (int *) malloc(sizeof(int)*count[0]);
    data5 = (int *) malloc(sizeof(int)*count[0]);
    data6 = (int *) malloc(sizeof(int)*count[0]);
    data7 = (int *) malloc(sizeof(int)*count[0]);
    data8 = (int *) malloc(sizeof(int)*count[0]);
    data9 = (int *) malloc(sizeof(int)*count[0]);

    for (i=0; i < count[0]; i++) {
      data1[i] = mpi_rank+10;
      data2[i] = mpi_rank+10;
      data3[i] = mpi_rank+10;
      data4[i] = mpi_rank+10;
      data5[i] = mpi_rank+10;
      data6[i] = mpi_rank+10;
      data7[i] = mpi_rank+10;
      data8[i] = mpi_rank+10;
      data9[i] = mpi_rank+10;
    }

    //  printf("proc %d: count NX = %ld \n", mpi_rank, count[0]);
    // printf("proc %d: offset NX = %ld \n", mpi_rank, offset[0]);

    /*
     * Select hyperslab in the file.
     */
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    timer_tick(&write_time, comm, 1);

    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 0);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data1);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 1);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data2);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 2);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data3);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 3);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data4);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 4);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data5);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 5);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data6);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 6);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data7);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 7);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data8);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 8);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data9);
    H5Dclose(dset_id);

    timer_tock(&write_time);
    timer_collectprintstats(write_time, comm, 0, "Time for H5Dwrite to complete");

    /*
     * Close/release resources.
     */
    H5Pclose(plist_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);

    /*
     *   _   _   _   _
     *  / \ / \ / \ / \
     * ( R | E | A | D )
     *  \_/ \_/ \_/ \_/
     */

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    /*
     * Create a new file collectively and release property list identifier.
     */
    file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);

    /*
     * Create the dataspace for the dataset.
     */
    filespace = H5Screate_simple(1, dimsf, NULL);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] =  dimsf[0]/mpi_size; 
    offset[0] = (mpi_rank)*count[0];
    memspace = H5Screate_simple(1, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    timer_tick(&read_time, comm, 1);

    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 0);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data1);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 1);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data2);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 2);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data3);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 3);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data4);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 4);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data5);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 5);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data6);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 6);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data7);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 7);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data8);
    H5Dclose(dset_id);
    snprintf(dset_name, 12, "%s_%06d",  DSET_NAME, 8);
    dset_id = H5Dopen(file_id, dset_name, H5P_DEFAULT);
    H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data9);
    H5Dclose(dset_id);

    timer_tock(&read_time);
    timer_collectprintstats(read_time, comm, 0, "Time for H5Dread  to complete");

    free(data1);
    free(data2);
    free(data3);
    free(data4);
    free(data5);
    free(data6);
    free(data7);
    free(data8);
    free(data9);

    /*
     * Close/release resources.
     */
    H5Pclose(plist_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);

    MPI_Finalize();

    return 0;
}    
