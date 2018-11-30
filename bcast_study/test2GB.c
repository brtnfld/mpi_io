/*  
 *  This example writes data to the HDF5 file by rows.
 *  Number of processes is assumed to be 1 or multiples of 2 (up to 8)
 */
 
#include "hdf5.h"
#include "stdlib.h"

#define H5FILE_NAME     "2GB.h5"
#define DATASETNAME 	"IntArray" 
#define NX     134217729  /* dataset dimensions */
#define NY     4 
#define RANK   2

int
main (int argc, char **argv)
{
    /*
     * HDF5 APIs definitions
     */ 	
    hid_t       file_id, dset_id;         /* file and dataset identifiers */
    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t     dimsf[2];                 /* dataset dimensions */
    int         *data=NULL;                    /* pointer to data buffer to write */
    hsize_t	count[2];	          /* hyperslab selection parameters */
    hsize_t	offset[2];
    hid_t	plist_id;                 /* property list identifier */
    int         i;
    herr_t	status;

    /*
     * MPI variables
     */
    int mpi_size, mpi_rank;
    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;
    MPI_Comm newcomm; 
    int color=1;

    /*
     * Initialize MPI
     */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    if(mpi_rank == 0) 
      color=0;
    MPI_Comm_split(comm, color, 0, &newcomm );

    dimsf[0] = NX;
    dimsf[1] = NY;
 
    if(color == 0) {

      /* 
       * Set up file access property list with parallel I/O access
       */
      plist_id = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(plist_id, newcomm, info);
      
      /*
       * Create a new file collectively and release property list identifier.
       */
      file_id = H5Fcreate(H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
      H5Pclose(plist_id);
      
      /*
       * Create the dataspace for the dataset.
       */
      filespace = H5Screate_simple(RANK, dimsf, NULL); 
      
      /*
       * Create the dataset with default properties and close filespace.
       */
      dset_id = H5Dcreate(file_id, DATASETNAME, H5T_NATIVE_INT, filespace,
			  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      printf("dataset size = %ld MiB\n", H5Dget_storage_size( dset_id )/1048576 ); 
      /*
       * Initialize data buffer 
       */
      data = (int *) malloc(sizeof(int)*dimsf[0]*dimsf[1]);
      for (i=0; i < dimsf[0]*dimsf[1]; i++) {
        data[i] = 20;
      }
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

      status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, data);

      free(data);
      H5Dclose(dset_id);
      H5Pclose(plist_id);
      H5Sclose(filespace);
      H5Fclose(file_id);
    }

    MPI_Barrier(comm);

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);

    dset_id = H5Dopen(file_id, DATASETNAME, H5P_DEFAULT);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    if((data = (int *) malloc(sizeof(int)*dimsf[0]*dimsf[1])) == NULL)
      printf("error in allocation \n");

    printf("Start H5Dread \n");
    MPI_Barrier(comm);

    status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, plist_id, data);

    printf("End H5Dread \n");
    MPI_Barrier(comm);

    printf("[%d]: data(NX,NY) = %d\n",mpi_rank, data[NX*NY-1]);

    /*
     * Close/release resources.
     */
    H5Dclose(dset_id);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    free(data);

    MPI_Finalize();

    return 0;
}     
