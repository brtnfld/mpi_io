!
! Number of processes is assumed to be 1 or multiples of 2 (1,2,4,6,8)
!

     PROGRAM DATASET_BY_COL

     USE HDF5 ! This module contains all necessary modules 
        
     IMPLICIT NONE

     INCLUDE 'mpif.h'
     INTEGER, PARAMETER :: NDIM=10
     CHARACTER(LEN=10), PARAMETER :: filename = "sds_col.h5"  ! File name
     CHARACTER(LEN=11), PARAMETER :: filename2 = "sds_col2.h5"  ! File name
     CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

     INTEGER(HID_T) :: file_id       ! File identifier
     INTEGER(HID_T) :: gid
     INTEGER(HID_T) :: dset_id       ! Dataset identifier 
     INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
     INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
     INTEGER(HID_T) :: plist_id,dcpl      ! Property list identifier 

     INTEGER(HSIZE_T), DIMENSION(1) :: dimsf = (/NDIM/) ! Dataset dimensions.
     INTEGER(HSIZE_T), DIMENSION(1) :: dimsfi = (/NDIM/)

     INTEGER(HSIZE_T), DIMENSION(1) :: count  
     INTEGER(HSSIZE_T), DIMENSION(1) :: offset 
     INTEGER, DIMENSION(1:NDIM):: DATA ! Data to write
     INTEGER :: rank = 1 ! Dataset rank 

     INTEGER :: i
     INTEGER :: error, error_n  ! Error flags
     !
     ! MPI definitions and calls.
     !
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, info
     INTEGER :: mpi_size, mpi_rank
     REAL*8, DIMENSION(1:4) :: t, timing
     REAL*8 :: t1

     comm = MPI_COMM_WORLD
     info = MPI_INFO_NULL
     CALL MPI_INIT(mpierror)
     CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
     CALL MPI_COMM_RANK(comm, mpi_rank, mpierror) 
     !
     ! Initialize FORTRAN predefined datatypes
     !
     CALL h5open_f(error) 

     t = 0.d0

     IF(mpi_rank.EQ.0)THEN

        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        
        CALL H5Pset_libver_bounds_f(plist_id, H5F_LIBVER_LATEST_F, H5F_LIBVER_LATEST_F, error)
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
        
        CALL H5Pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
        CALL H5Pset_alloc_time_f(dcpl, H5D_ALLOC_TIME_EARLY_F, error)
        CALL h5pset_fill_time_f(dcpl, H5D_FILL_TIME_NEVER_F, error) 
        !
        ! Create the data space for the  dataset. 
        !
        CALL h5screate_simple_f(rank, dimsf, filespace, error)
        
        t1 = MPI_Wtime()
        CALL h5gcreate_f(file_id, "G1", gid, error)
        t(1) = MPI_Wtime() - t1

        !
        ! Create the dataset with default properties.
        !
        CALL h5dcreate_f(gid, dsetname, H5T_NATIVE_INTEGER, filespace, &
             dset_id, error, dcpl_id=dcpl)
        CALL h5sclose_f(filespace, error)

        CALL h5pclose_f(dcpl, error)
        !
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file. 
        !
        COUNT(1) = dimsf(1)
        offset(1) = 0
        CALL h5screate_simple_f(rank, count, memspace, error) 
        ! 
        ! Select hyperslab in the file.
        !
        CALL h5dget_space_f(dset_id, filespace, error)
        CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)
        ! 
        ! Initialize data buffer with trivial data.
        !
        DO i = 1, NDIM
           DATA(i) = i
        ENDDO
        
        !
        ! Write the dataset
        !
        t1 = MPI_Wtime()
        CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, DATA, dimsf, error, &
             file_space_id = filespace, mem_space_id = memspace)
        t(2) = MPI_Wtime() - t1
        !
        ! Close dataspaces.
        !
        CALL h5sclose_f(filespace, error)
        CALL h5sclose_f(memspace, error)
        
        !
        ! Close the dataset and property list.
        !
        CALL h5dclose_f(dset_id, error)
        CALL h5pclose_f(plist_id, error)
        
        !
        ! Close the file.
        !
        CALL h5gclose_f(gid, error)
        CALL h5fclose_f(file_id, error)
        
     ENDIF

     CALL MPI_Barrier(MPI_COMM_WORLD,error)

     ! 
     ! Setup file access property list with parallel I/O access.
     !
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
     CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)

     !
     ! Create the file collectively.
     ! 
     CALL h5fcreate_f(filename2, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     CALL h5pclose_f(plist_id, error)
     !
     ! Create the data space for the  dataset. 
     !
     CALL h5screate_simple_f(rank, dimsf, filespace, error)

     t1 = MPI_Wtime()
     CALL h5gcreate_f(file_id, "G1", gid, error)
     t(3) = MPI_Wtime() - t1
     !
     ! Create the dataset with default properties.
     !
     CALL h5dcreate_f(gid, dsetname, H5T_NATIVE_INTEGER, filespace, &
                      dset_id, error)
     CALL h5sclose_f(filespace, error)
     !
     ! Each process defines dataset in memory and writes it to the hyperslab
     ! in the file. 
     !
     count(1) = dimsf(1)
     offset(1) = 0
     CALL h5screate_simple_f(rank, count, memspace, error) 
     ! 
     ! Select hyperslab in the file.
     !
     CALL h5dget_space_f(dset_id, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)
     ! 
     ! Initialize data buffer with trivial data.
     !
     DO i = 1, NDIM
        DATA(i) = i ! mpi_rank
     ENDDO
     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     
     t1 = MPI_Wtime()
     !
     ! Write the dataset collectively. 
     !
!     CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsfi, error, &
!                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
     !
     ! Write the dataset independently. 
     !

     CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, DATA, dimsfi, error, &
          file_space_id = filespace, mem_space_id = memspace)
     t(4) = MPI_Wtime() - t1

     !
     ! Close dataspaces.
     !
     CALL h5sclose_f(filespace, error)
     CALL h5sclose_f(memspace, error)

     !
     ! Close the dataset and property list.
     !
     CALL h5dclose_f(dset_id, error)
     CALL h5pclose_f(plist_id, error)
     CALL h5gclose_f(gid,error)
     !
     ! Close the file.
     !
     CALL h5fclose_f(file_id, error)

     !
     ! Close FORTRAN predefined datatypes.
     !
     CALL h5close_f(error)

     CALL MPI_Reduce(t, timing, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, error)


     IF(mpi_rank.EQ.0)THEN
        ! Total time
        WRITE(*,'( (I0,X), 4(e14.6,X))') mpi_size, timing(1), timing(2), timing(3)/mpi_size, timing(4)/mpi_size
     ENDIF

     CALL MPI_FINALIZE(mpierror)

     END PROGRAM DATASET_BY_COL
