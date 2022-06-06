/*  
 *  This example writes data to the HDF5 file by rows.
 *  Number of processes is assumed to be 1 or multiples of 2 (up to 8)
 */
 
#include "hdf5.h"
#include "stdlib.h"
#include <unistd.h>
#include <math.h>

#define H5FILE_NAME     "SDS_row.h5"
#define DATASETNAME     "IntArray" 
#define NX     2048                      /* dataset dimensions */
#define NY     4 
#define RANK   2

static int PImpi(int pe, int processes, unsigned int intervals) {

    double time1 = MPI_Wtime();

    int count = intervals / processes;
    int start = count * pe;
    int end = count * pe + count;

    int i;
    double subtotal, total = 0;
    for (i = start; i < end; ++i) {
        subtotal += pow(-1, i) / (2 * i + 1);
    }

    MPI_Reduce(&subtotal, &total, 1, MPI_DOUBLE, MPI_SUM,
        0, MPI_COMM_WORLD);

    double time2 = MPI_Wtime();

    if (pe == 1) {
        total = total * 4;
        printf("Result:   %.10lf\n", total);
        printf("Time:     %.10lf\n", time2 - time1);
    }

    return 0;
}

int
main (int argc, char **argv)
{
    /*
     * HDF5 APIs definitions
     */ 
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[2];                 /* dataset dimensions */
    int         *data;                    /* pointer to data buffer to write */
    hsize_t     count[2];                 /* hyperslab selection parameters */
    hsize_t     offset[2];
    hid_t       plist_id;                 /* property list identifier */
    int         i;
    herr_t      status;

    /*
     * MPI variables
     */
    int mpi_size, mpi_rank;
    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;

    hid_t es_id_g;
    size_t num_in_progress;
    hbool_t op_failed;

    /*
     * Initialize MPI
     */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);  
 
    /* 
     * Set up file access property list with parallel I/O access
     */
     plist_id = H5Pcreate(H5P_FILE_ACCESS);
     H5Pset_fapl_mpio(plist_id, comm, info);

    /*
     * Create a new file collectively and release property list identifier.
     */
     file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
     H5Pclose(plist_id);
   

     es_id_g = H5EScreate();
     
    /*
     * Create the dataspace for the dataset.
     */
    dimsf[0] = NX;
    dimsf[1] = NY;
    filespace = H5Screate_simple(RANK, dimsf, NULL); 

    /*
     * Create the dataset with default properties and close filespace.
     */
#if 0
    dset_id = H5Dcreate_async(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_id_g);
#else
    dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
    H5Sclose(filespace);

    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    memspace = H5Screate_simple(RANK, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    /*
     * Initialize data buffer 
     */
    data = (int *) malloc(sizeof(int)*count[0]*count[1]);
    for (i=0; i < count[0]*count[1]; i++) {
        data[i] = mpi_rank + 10;
    }

    /*
     * Create property list for collective dataset write.
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    
    status = H5Dwrite_async(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, data, es_id_g);

    /*
     * Close/release resources.
     */
#if 0
    H5Dclose_async(dset_id, es_id_g);
#endif

    if(mpi_rank == 0) {
      H5ESwait(es_id_g, 10E9, &num_in_progress, &op_failed);
    } else {
      PImpi(mpi_rank, mpi_size-1, 1E9);
    }

    double time1 = MPI_Wtime();
    H5ESwait(es_id_g, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    MPI_Barrier(MPI_COMM_WORLD);
    double time2 = MPI_Wtime();
    printf("H5ESwait Time:     %.10lf\n", time2 - time1);

#if 0
    H5ESwait(es_id_g, 0, &num_in_progress, &op_failed);
    printf("Before compute phase %ld \n", num_in_progress);
    MPI_Barrier(MPI_COMM_WORLD);
    PImpi(mpi_rank, mpi_size, 1E9);
    H5ESwait(es_id_g, 1E8, &num_in_progress, &op_failed);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("After compute phase %ld \n", num_in_progress);
    PImpi(mpi_rank, mpi_size, 1E9);
    H5ESget_count(es_id_g, &num_in_progress);
    printf("BB %ld \n", num_in_progress);
    MPI_Barrier(MPI_COMM_WORLD);
    double time1 = MPI_Wtime();
    H5ESwait(es_id_g, H5ES_WAIT_FOREVER, &num_in_progress, &op_failed);
    MPI_Barrier(MPI_COMM_WORLD);
    double time2 = MPI_Wtime();
    printf("H5ESwait Time:     %.10lf\n", time2 - time1);
#endif

    

    printf("C %ld \n", num_in_progress);
    free(data);
    H5Dclose(dset_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Fclose(file_id);
 
    MPI_Finalize();

    return 0;
}     

