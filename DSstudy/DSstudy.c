/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the COPYING file, which can be found at the root of the source code       *
 * distribution tree, or in https://support.hdfgroup.org/ftp/HDF5/releases.  *
 * If you do not have access to either file, you may request a copy from     *
 * help@hdfgroup.org.                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "hdf5.h"
#include "hdf5_hl.h"
#include <mpi.h>
#include <stdlib.h> 
#include <stdio.h>

#define RANK      2
#define DIM_DATA  12
#define DIM1_SIZE 3
#define DIM2_SIZE 4
#define DIM0      0
#define DIM1      1
#define NDSET     20

#define DSET_NAME "data0"
#define DS_1_NAME "Yaxis"
#define DS_2_NAME "Xaxis"


int
main(int argc, char *argv[])
{

    hid_t   fid;                                                          /* file ID */
    hid_t   did;                                                          /* dataset ID */
    hid_t   dsid;                                                         /* DS dataset ID */
    int     rank               = RANK;                                    /* rank of data dataset */
    int     rankds             = 1;                                       /* rank of DS dataset */
    hsize_t dims[RANK]         = {DIM1_SIZE, DIM2_SIZE};                  /* size of data dataset */
    int     buf[DIM_DATA]      = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}; /* data of data dataset */
    hsize_t s1_dim[1]          = {DIM1_SIZE};                             /* size of DS 1 dataset */
    hsize_t s2_dim[1]          = {DIM2_SIZE};                             /* size of DS 2 dataset */
    float   s1_wbuf[DIM1_SIZE] = {10, 20, 30};                            /* data of DS 1 dataset */
    int     s2_wbuf[DIM2_SIZE] = {10, 20, 50, 100};                       /* data of DS 2 dataset */

    int OneRank = 0;
    double e1,et;

    MPI_Init(&argc, &argv);

    if( argc == 2 ) {
      OneRank = atoi(argv[1]);
    }

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(myid == 0) {

      /* create a file using default properties */
      if ((fid = H5Fcreate("ex_ds1.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
        goto out;

      for ( int i=0; i < NDSET; i++) {

        char dset_name[9];
        char ds_1_name[9];
        char ds_2_name[9];

        snprintf(dset_name, 9, "%s_%d",  DSET_NAME, i); 
        snprintf(ds_1_name, 9, "%s_%d",  DS_1_NAME, i);
        snprintf(ds_2_name, 9, "%s_%d",  DS_2_NAME, i);

        /* make a dataset */
        if (H5LTmake_dataset_int(fid, dset_name, rank, dims, buf) < 0)
          goto out;
        
        /* make a DS dataset for the first dimension */
        if (H5LTmake_dataset_float(fid, ds_1_name, rankds, s1_dim, s1_wbuf) < 0)
          goto out;
        
        /* make a DS dataset for the second dimension */
        if (H5LTmake_dataset_int(fid, ds_2_name, rankds, s2_dim, s2_wbuf) < 0)
          goto out;

      }

      /* close file */
      H5Fclose(fid);

    }

    hid_t fapl;
    fapl = H5Pcreate(H5P_FILE_ACCESS);

    MPI_Comm newcomm = MPI_COMM_WORLD;

    if( OneRank ) {
      MPI_Comm_split(MPI_COMM_WORLD, myid, 0, &newcomm);
    }
    H5Pset_fapl_mpio(fapl, newcomm, MPI_INFO_NULL);

    e1 = MPI_Wtime();

    if ( OneRank ) {
      if (myid > 0) goto skip;
    }

    fid = H5Fopen("ex_ds1.h5", H5F_ACC_RDWR, fapl);

    for ( int i=0; i < NDSET; i++) {

      char dset_name[9];
      char ds_1_name[9];
      char ds_2_name[9];
      
      snprintf(dset_name, 9, "%s_%d",  DSET_NAME, i); 
      snprintf(ds_1_name, 9, "%s_%d",  DS_1_NAME, i);
      snprintf(ds_2_name, 9, "%s_%d",  DS_2_NAME, i);
    /*-------------------------------------------------------------------------
     * attach the DS_1_NAME dimension scale to DSET_NAME at dimension 0
     *-------------------------------------------------------------------------
     */

    /* get the dataset id for DSET_NAME */
      if ((did = H5Dopen2(fid, dset_name, H5P_DEFAULT)) < 0)
        goto out;
      
      /* get the DS dataset id */
      if ((dsid = H5Dopen2(fid, ds_1_name, H5P_DEFAULT)) < 0)
        goto out;
      
      /* attach the DS_1_NAME dimension scale to DSET_NAME at dimension index 0 */
      if (H5DSattach_scale(did, dsid, DIM0) < 0)
        goto out;

      htri_t is_scale= H5DSis_scale(dsid);
      if ( is_scale != 1)
        goto out;

      //char dimscale_name[11];
      //if (H5DSget_scale_name(dsid, dimscale_name, sizeof(dimscale_name)) < 0 )
      // goto out;

      // printf("%s \n", dimscale_name);

      /* close DS id */
      if (H5Dclose(dsid) < 0)
        goto out;
      
      /*-------------------------------------------------------------------------
       * attach the DS_2_NAME dimension scale to DSET_NAME
       *-------------------------------------------------------------------------
       */
      
      /* get the DS dataset id */
      if ((dsid = H5Dopen2(fid, ds_2_name, H5P_DEFAULT)) < 0)
        goto out;
      
      /* attach the DS_2_NAME dimension scale to DSET_NAME as the 2nd dimension (index 1)  */
      if (H5DSattach_scale(did, dsid, DIM1) < 0)
        goto out;
      
      /* close DS id */
      if (H5Dclose(dsid) < 0)
        goto out;
      if (H5Dclose(did) < 0)
        goto out;
    }

    /* close file */
    H5Fclose(fid);

 skip:
   H5Pclose(fapl);

   MPI_Barrier(MPI_COMM_WORLD);

   et = MPI_Wtime() - e1;
   if(myid == 0){
     FILE *fp;
     fp = fopen("rate", "w+");
     printf("nprocs = %d time = %lf \n", nprocs, et);
     fprintf(fp, "%lf\n", et);
     fclose(fp);
   }

   return 0;

out:
    printf("Error on return function...Exiting\n");
    return 1;
}
