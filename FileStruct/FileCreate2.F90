!
! Number of processes is assumed to be 1 or multiples of 2 (1,2,4,6,8)
!


PROGRAM DATASET_BY_COL

  USE HDF5 ! This module contains all necessary modules 
  USE MPI
  IMPLICIT NONE

  CHARACTER(LEN=10), PARAMETER :: filename = "sds_col.h5"  ! File name
  CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

  INTEGER(HID_T) :: file_id       ! File identifier 
  INTEGER(HID_T) :: dset_id, datatype_id       ! Dataset identifier 
  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace, mid      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier
  INTEGER(HID_T) :: lcpl_id, dcpl

  INTEGER(HSIZE_T), DIMENSION(1) :: dimsf = (/10/) ! Dataset dimensions.

  INTEGER(HSIZE_T), DIMENSION(1) :: count  
  INTEGER(HSSIZE_T), DIMENSION(1) :: offset 
  INTEGER, ALLOCATABLE :: data (:,:)  ! Data to write
  INTEGER :: rank = 1 ! Dataset rank 

  REAL*8 :: t0, t1, t2, t3, t4, t5, timing

  INTEGER :: error  ! Error flags
  !
  ! MPI definitions and calls.
  !
  INTEGER :: mpierror       ! MPI error flag
  INTEGER :: comm, info
  INTEGER :: mpi_size, mpi_rank

  INTEGER :: i, j
  CHARACTER(LEN=4) :: id1,id2
  INTEGER :: depth1, depth2
  character(LEN=180) :: dset

  CHARACTER(len=128) :: arg
  CHARACTER(len=1) :: argv
  INTEGER k
  INTEGER PROC0
  TYPE(C_PTR) :: f_ptr

#define DEBUG 0

  depth1 = 128
  depth2 = 32

  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL
  CALL MPI_INIT(mpierror)
  CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
  CALL MPI_COMM_RANK(comm, mpi_rank, mpierror)


  PROC0 = 0
  k = 0
  argv=""
  DO
     CALL get_command_argument(k, arg)
     IF (LEN_TRIM(arg) == 0) EXIT
     argv(1:1) = arg(1:1)
     k = k + 1
  END DO

  IF(argv .EQ. '0')THEN
     PROC0=1
  ENDIF

  !
  ! Initialize FORTRAN predefined datatypes
  !
  CALL h5open_f(error) 

  IF(mpi_rank.EQ.0)THEN
        
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        
     CALL H5Pset_libver_bounds_f(plist_id, H5F_LIBVER_LATEST_F, H5F_LIBVER_LATEST_F, error)
     CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     CALL H5Pclose_f(plist_id, error)
        
     CALL H5Pcreate_f(H5P_LINK_CREATE_F, lcpl_id, error)
     CALL H5Pset_create_inter_group_f(lcpl_id, 1, error)
     !
     ! Create the data space for the  dataset. 
     !
     CALL h5screate_simple_f(rank, dimsf, filespace, error)
        
     CALL H5Pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
     CALL H5Pset_alloc_time_f(dcpl, H5D_ALLOC_TIME_EARLY_F, error)
     CALL h5pset_fill_time_f(dcpl, H5D_FILL_TIME_NEVER_F, error)
        
     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
     DO i = 1, depth1
        WRITE(id1,"(I4.4)") i
        DO j = 1, depth2
           WRITE(id2,"(I4.4)") j
           !
           ! Create the dataset with default properties.
           !
           CALL h5dcreate_f(file_id, "/base_"//id1//"/zone_"//id2//"/"//TRIM(dsetname), H5T_NATIVE_INTEGER, & 
                filespace, dset_id, error,lcpl_id=lcpl_id, dcpl_id=dcpl)
           CALL h5dclose_f(dset_id, error)
        ENDDO
     ENDDO
     
     !
     ! Close dataspaces.
     !
     CALL h5sclose_f(filespace, error)
     !
     ! Close the dataset and property list.
     !
     CALL h5pclose_f(plist_id, error)
     CALL h5pclose_f(dcpl, error)
     CALL h5pclose_f(lcpl_id, error)
     !
     ! Close the file.
     !
     CALL h5fclose_f(file_id, error)
     t2 = MPI_Wtime() - t1
  ENDIF

  CALL MPI_BARRIER(MPI_COMM_WORLD, error)

  t0 = MPI_Wtime()

  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)

  t2 = 0.
  DO i = 1, depth1
     WRITE(id1,"(I4.4)") i
     DO j = 1, depth2
        WRITE(id2,"(I4.4)") j
        !
        ! Create the dataset with default properties.
        !
        dset = "/base_"//id1//"/zone_"//id2//"/"//TRIM(dsetname)
        CALL h5dopen_f(file_id, TRIM(dset), dset_id, error)

        t1 = MPI_Wtime()
        CALL h5dget_type_f(dset_id, datatype_id, error)
        CALL H5Tget_native_type_f(datatype_id, H5T_DIR_ASCEND_F, mid, error)

        
        t2 = t2 + MPI_Wtime() - t1
        CALL h5dclose_f(dset_id, error)
        CALL h5tclose_f(mid, error)
     ENDDO
  ENDDO


  CALL h5fclose_f(file_id,error)

  call MPI_Reduce(t2, timing, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, error)

  t4 = MPI_Wtime()
  IF(mpi_rank.EQ.0)THEN
     ! Total time, H5Fclose
     WRITE(*,'( (I0,X), 2(f14.6,X))') mpi_size, timing/mpi_size
  ENDIF
  !
  ! Close FORTRAN predefined datatypes.
  !
  CALL h5close_f(error)

  CALL MPI_FINALIZE(mpierror)

END PROGRAM DATASET_BY_COL
