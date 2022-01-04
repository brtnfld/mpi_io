/*  
 *  This example writes data to the HDF5 file by rows.
 *  Number of processes is assumed to be 1 or multiples of 2 (up to 8)
 */
 
#include "hdf5.h"
#include "stdlib.h"
#include <string.h>

#define H5FILE_NAME     "2GB.h5"
#define DATASETNAME 	"IntArray" 
#define NX     134217727  /* dataset dimensions */
#define NY     4 
#define RANK   2

int
main (int argc, char **argv)
{
    /*
     * HDF5 APIs definitions
     */ 	
  hid_t       file_id, gid, lcpl_id;         /* file and dataset identifiers */
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
    int ngroups=100000;
    int color=1;
    char buffer[1000];
    double t1, t2, t3;
    double Max_total_time3 = 0;
    double Min_total_time3 = 0;
    double Sum_total_time3 = 0;
    double Max_total_time2 = 0;
    double Min_total_time2 = 0;
    double Sum_total_time2 = 0;
    htri_t exists;
    int h5g=0, h5l=0;
    
    /*
     * Initialize MPI
     */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    if(mpi_rank == mpi_size-1) 
      color=0;
    MPI_Comm_split(comm, color, 1, &newcomm );


    if (argc==1){
      if (mpi_rank==0) printf("REQUIRES INPUT ARGUMENT -g or -l \n");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if(strcmp(argv[1],"-g")==0) {
      h5g=1;
    } else {
      h5l=1;
    }

 
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

      lcpl_id = H5Pcreate(H5P_LINK_CREATE);
      H5Pset_create_intermediate_group( lcpl_id, 1);

      for (i=0; i < ngroups; i++) {
        sprintf(buffer, "%s%02d", "link01/link02/link03/link04/",i);
        gid = H5Gcreate2(file_id, buffer, lcpl_id, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(gid);
      }

      H5Fclose(file_id);
      H5Pclose(lcpl_id);
    }

    MPI_Barrier(comm);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);


    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);

    if( h5g == 1) {
      file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, plist_id);
      t1 = MPI_Wtime();
      
      for (i=0; i < ngroups; i++) {
        sprintf(buffer, "%s%02d", "link01/link02/link03/link04/",i);
        gid = H5Gopen2(file_id, buffer, H5P_DEFAULT);
        if( gid > 0) {
          H5Gclose(gid);
        }
      }
      t2 = MPI_Wtime() - t1;
      
      H5Pclose(plist_id);
      H5Fclose(file_id);
    }

    if( h5l == 1) {
      file_id = H5Fopen(H5FILE_NAME, H5F_ACC_RDWR, plist_id);

      t1 = MPI_Wtime();
      for (i=0; i < ngroups; i++) {
        sprintf(buffer, "%s%02d", "link01/link02/link03/link04/",i);
        exists = H5Lexists(file_id, buffer, H5P_DEFAULT);
        if( exists > 0 ) {
            gid = H5Gopen2(file_id, buffer, H5P_DEFAULT);
            if( gid > 0)
              H5Gclose(gid);
        }
      }
      t2 = MPI_Wtime() - t1;

      H5Pclose(plist_id);
      H5Fclose(file_id);
    }

    MPI_Reduce(&t2, &Max_total_time2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t2, &Sum_total_time2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t2, &Min_total_time2, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

#if 0
    MPI_Reduce(&t3, &Max_total_time3, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t3, &Sum_total_time3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t3, &Min_total_time3, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
#endif

    if(mpi_rank == 0) {
      printf("%d Procs WRITE %d groups\n", mpi_size, ngroups);
      printf("h5g=%d, h5l=%d: time is avg, min, max: %lf %f %f s.\n", h5g, h5l, Sum_total_time2/mpi_size, Min_total_time2, Max_total_time2);
    }

    MPI_Finalize();

    return 0;
}     
