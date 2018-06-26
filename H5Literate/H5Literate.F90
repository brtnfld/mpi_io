!
! Number of processes is assumed to be 1 or multiples of 2 (1,2,4,6,8)
!
MODULE liter_cb_mod

  USE HDF5
  USE, INTRINSIC :: ISO_C_BINDING
  IMPLICIT NONE

  ! Custom group iteration callback data 
  TYPE, bind(c) ::  iter_info
     CHARACTER(KIND=C_CHAR), DIMENSION(1:10) :: name !  The name of the object 
     INTEGER(c_int) :: TYPE    !  The TYPE of the object 
     INTEGER(c_int) :: command ! The TYPE of RETURN value 
  END TYPE iter_info

CONTAINS
!***************************************************************
!**
!**  liter_cb(): Custom link iteration callback routine.
!**
!***************************************************************

  INTEGER(KIND=C_INT) FUNCTION liter_cb(group, name, link_info, op_data) bind(C)

    IMPLICIT NONE

    INTEGER(HID_T), VALUE :: group
    CHARACTER(LEN=1), DIMENSION(1:10) :: name
    CHARACTER(LEN=10) :: name2

    INTEGER :: error
    TYPE (H5L_info_t) :: link_info
    TYPE(H5O_info_t) :: infobuf

    TYPE(iter_info) :: op_data
    INTEGER i

    INTEGER, SAVE :: count
    INTEGER, SAVE :: count2

    DO i = 1, 10
       IF(name2(i:i).EQ.'char(0)') EXIT
       name2(i:i) = name(i)
    ENDDO
    !PRINT*,name2
    op_data%type = op_data%type + 1

 !   CALL H5Oget_info_by_name_f (group, TRIM(name2), infobuf, error)

  !  PRINT*,error

   ! IF(infobuf%type .EQ. H5O_TYPE_GROUP_F)THEN
  !     liter_cb = liter_cb + 1
  !  ENDIF

  END FUNCTION liter_cb
END MODULE liter_cb_mod

PROGRAM DATASET_BY_COL

  USE liter_cb_mod
  USE HDF5 ! This module contains all necessary modules 
  USE MPI
  IMPLICIT NONE

  CHARACTER(LEN=10), PARAMETER :: filename = "sds_col.h5"  ! File name
  CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

  INTEGER(HID_T) :: file_id       ! File identifier 
  INTEGER(HID_T) :: dset_id       ! Dataset identifier 
  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier
  INTEGER(HID_T) :: lcpl_id, dcpl

  INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/4,131072/) ! Dataset dimensions.

  INTEGER(HSIZE_T), DIMENSION(2) :: count  
  INTEGER(HSSIZE_T), DIMENSION(2) :: offset 
  INTEGER, ALLOCATABLE :: data (:,:)  ! Data to write
  INTEGER :: rank = 2 ! Dataset rank 

  REAL*8 :: t0, t1, t2, t3, t4, t5

  INTEGER :: error  ! Error flags
  !
  ! MPI definitions and calls.
  !
  INTEGER :: mpierror       ! MPI error flag
  INTEGER :: comm, cinfo
  INTEGER :: mpi_size, mpi_rank

  INTEGER :: i, j
  CHARACTER(LEN=4) :: id1,id2
  INTEGER :: depth1, depth2

  CHARACTER(len=128) :: arg
  CHARACTER(len=1) :: argv
  INTEGER k
  INTEGER PROC0
  TYPE(C_PTR) :: f_ptr

  TYPE(iter_info), TARGET :: info
  INTEGER(hsize_t) :: idx ! Index in the group 

  TYPE(C_FUNPTR) :: f1
  TYPE(C_PTR) :: f2
  INTEGER :: ret_value

#define DEBUG 0

  depth1 = 64
  depth2 = 4

  comm = MPI_COMM_WORLD
  cinfo = MPI_INFO_NULL
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


  IF(PROC0.eq.1)THEN
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
              CALL h5dcreate_f(file_id, "/"//id1//"/"//id2//"/"//TRIM(dsetname), H5T_NATIVE_INTEGER, & 
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
        t1 = MPI_Wtime()
        CALL h5fclose_f(file_id, error)
        t2 = MPI_Wtime() - t1
     ENDIF

  ELSE

     ! 
     ! Setup file access property list with parallel I/O access.
     !
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)

     CALL h5pset_fapl_mpio_f(plist_id, comm, cinfo, error)

     !
     ! Create the file collectively.
     !

     CALL H5Pset_libver_bounds_f(plist_id, H5F_LIBVER_LATEST_F, H5F_LIBVER_LATEST_F, error)

     CALL h5pset_coll_metadata_write_f(plist_id, .TRUE., error)
     CALL h5pset_all_coll_metadata_ops_f(plist_id, .TRUE., error)
     
     CALL H5Fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     CALL H5Pclose_f(plist_id, error)

     CALL H5Pcreate_f(H5P_LINK_CREATE_F, lcpl_id, error)
     CALL H5Pset_create_inter_group_f(lcpl_id, 1, error)

     !
     ! Create the data space for the  dataset. 
     !
     CALL h5screate_simple_f(rank, dimsf, filespace, error)
     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

     CALL H5Pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
     CALL h5pset_fill_time_f(dcpl, H5D_FILL_TIME_NEVER_F, error)

     DO i = 1, depth1
        WRITE(id1,"(I4.4)") i
        DO j = 1, depth2
           WRITE(id2,"(I4.4)") j
           !
           ! Create the dataset with default properties.
           !
           CALL h5dcreate_f(file_id, "/"//id1//"/"//id2//"/"//TRIM(dsetname), H5T_NATIVE_INTEGER, & 
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
     CALL h5pclose_f(lcpl_id, error)
     !
     ! Close the file.
     !
     t1 = MPI_Wtime()
     CALL h5fclose_f(file_id, error)
     CALL MPI_BARRIER( MPI_COMM_WORLD, error)
     t2 = MPI_Wtime() - t1
  ENDIF


  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

  ! begin the iteration
  info%command = 0
  idx = 0
  f1 = C_FUNLOC(liter_cb)
  f2 = C_LOC(info)

  t0 = MPI_Wtime()

  CALL H5Literate_f(file_id, H5_INDEX_NAME_F, H5_ITER_INC_F, idx, f1, f2, ret_value, error)

  CALL MPI_BARRIER( MPI_COMM_WORLD, error)
  t4 = MPI_Wtime()

!  IF(mpi_rank.EQ.0) PRINT*,info%type !ret_value

  call h5fclose_f(file_id,error)

  IF(mpi_rank.EQ.0)THEN
     ! Total time, H5Fclose
     WRITE(*,'( 1(I0,X), 2(f7.4,X))') mpi_size, t4-t0
  ENDIF
  !
  ! Close FORTRAN predefined datatypes.
  !
  CALL h5close_f(error)

  CALL MPI_FINALIZE(mpierror)

END PROGRAM DATASET_BY_COL
