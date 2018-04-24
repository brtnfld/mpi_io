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
  INTEGER(HID_T) :: dset_id       ! Dataset identifier 
  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier
  INTEGER(HID_T) :: lcpl_id

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
  INTEGER :: comm, info
  INTEGER :: mpi_size, mpi_rank

  INTEGER :: i, j
  CHARACTER(LEN=4) :: id1,id2
  INTEGER :: depth1, depth2

  TYPE(H5AC_cache_config_t) :: config

#define DEBUG 0

  depth1 = 32
  depth2 = 128

  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL
  CALL MPI_INIT(mpierror)
  CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
  CALL MPI_COMM_RANK(comm, mpi_rank, mpierror) 
  !
  ! Initialize FORTRAN predefined datatypes
  !
  CALL h5open_f(error) 

  ! 
  ! Setup file access property list with parallel I/O access.
  !
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)

  !
  ! Create the file collectively.
  !
  CALL H5Pset_libver_bounds_f(plist_id, H5F_LIBVER_LATEST_F, H5F_LIBVER_LATEST_F, error)

  ! CALL h5pset_coll_metadata_write_f(plist_id, .TRUE., error)
  ! CALL h5pset_all_coll_metadata_ops_f(plist_id, .TRUE., error)

  config%version = H5AC_CURR_CACHE_CONFIG_VER_F;
  CALL h5pget_mdc_config_f(plist_id,  config, error)
#if DEBUG
  IF(mpi_rank.EQ.0)THEN
     PRINT*,  config%version
     PRINT*,  config%rpt_fcn_enabled
     PRINT*,  config%open_trace_file
     PRINT*,  config%close_trace_file
     PRINT*,  config%trace_file_name
     PRINT*,  config%evictions_enabled
     PRINT*,  config%set_initial_size
     PRINT*,  config%initial_size
     PRINT*,  config%min_clean_fraction
     PRINT*,  config%max_size
     PRINT*,  config%min_size
     PRINT*,  config%epoch_length
     PRINT*,  config%incr_mode ! H5C_cache_incr_mode
     PRINT*,  config%lower_hr_threshold
     PRINT*,  config%increment
     PRINT*,  config%apply_max_increment
     PRINT*,  config%max_increment
     PRINT*,  config%flash_incr_mode ! enum H5C_cache_flash_incr_mode 
     PRINT*,  config%flash_multiple
     PRINT*,  config%flash_threshold
     PRINT*,  config%decr_mode ! enum H5C_cache_decr_mode 
     PRINT*,  config%upper_hr_threshold
     PRINT*,  config%decrement
     PRINT*,  config%apply_max_decrement
     PRINT*,  config%max_decrement
     PRINT*,  config%epochs_before_eviction
     PRINT*,  config%apply_empty_reserve
     PRINT*,  config%empty_reserve
     PRINT*,  config%dirty_bytes_threshold
     PRINT*,  config%metadata_write_strategy
  ENDIF
#endif
  config%metadata_write_strategy = 0

  CALL h5pset_mdc_config_f(plist_id,  config, error)

#if DEBUG
  CALL h5pget_mdc_config_f(plist_id,  config, error)
  IF(mpi_rank.EQ.0)THEN
     PRINT*,'h5pget_mdc_config_f'
     PRINT*,  config%version
     PRINT*,  config%rpt_fcn_enabled
     PRINT*,  config%open_trace_file
     PRINT*,  config%close_trace_file
     PRINT*,  config%trace_file_name
     PRINT*,  config%evictions_enabled
     PRINT*,  config%set_initial_size
     PRINT*,  config%initial_size
     PRINT*,  config%min_clean_fraction
     PRINT*,  config%max_size
     PRINT*,  config%min_size
     PRINT*,  config%epoch_length
     PRINT*,  config%incr_mode ! H5C_cache_incr_mode
     PRINT*,  config%lower_hr_threshold
     PRINT*,  config%increment
     PRINT*,  config%apply_max_increment
     PRINT*,  config%max_increment
     PRINT*,  config%flash_incr_mode ! enum H5C_cache_flash_incr_mode 
     PRINT*,  config%flash_multiple
     PRINT*,  config%flash_threshold
     PRINT*,  config%decr_mode ! enum H5C_cache_decr_mode 
     PRINT*,  config%upper_hr_threshold
     PRINT*,  config%decrement
     PRINT*,  config%apply_max_decrement
     PRINT*,  config%max_decrement
     PRINT*,  config%epochs_before_eviction
     PRINT*,  config%apply_empty_reserve
     PRINT*,  config%empty_reserve
     PRINT*,  config%dirty_bytes_threshold
     PRINT*,  config%metadata_write_strategy
  endif
#endif

  t0 = MPI_Wtime()
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  CALL h5pclose_f(plist_id, error)

  CALL H5Pcreate_f(H5P_LINK_CREATE_F, lcpl_id, error)
  CALL H5Pset_create_inter_group_f(lcpl_id, 1, error)

  !
  ! Create the data space for the  dataset. 
  !
  CALL h5screate_simple_f(rank, dimsf, filespace, error)
  !
  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file. 
  !
  count(1) = dimsf(1)
  count(2) = dimsf(2)/mpi_size
  offset(1) = 0
  offset(2) = mpi_rank * count(2) 
  CALL h5screate_simple_f(rank, count, memspace, error) 
  ! 
  ! Initialize data buffer with trivial data.
  !
  ALLOCATE ( data(count(1),count(2)))
  data = mpi_rank + 10
  !
  ! Create property list for collective dataset write
  !
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  DO i = 1, depth1
     WRITE(id1,"(I4.4)") i
     DO j = 1, depth2
        WRITE(id2,"(I4.4)") j
        !
        ! Create the dataset with default properties.
        !
        CALL h5dcreate_f(file_id, "/"//id1//"/"//id2//"/"//TRIM(dsetname), H5T_NATIVE_INTEGER, & 
             filespace, dset_id, error,lcpl_id=lcpl_id)
        ! 
        ! Select hyperslab in the file.
        !
        CALL h5dget_space_f(dset_id, filespace, error)
        CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)

        !
        ! Write the dataset collectively. 
        !
        CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, DATA, count, error, &
             file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

        CALL h5dclose_f(dset_id, error)
     ENDDO
  ENDDO
  !
  ! Close dataspaces.
  !
  CALL h5sclose_f(filespace, error)
  CALL h5sclose_f(memspace, error)

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

  !
  ! Close FORTRAN predefined datatypes.
  !
  CALL h5close_f(error)
  CALL MPI_BARRIER( MPI_COMM_WORLD, error)
  t3 = MPI_Wtime()
  IF(mpi_rank.EQ.0)THEN
     WRITE(*,'(f7.4)') t3-t0, t2
  ENDIF

  !
  ! Deallocate data buffer.
  !
  DEALLOCATE(data)

  CALL MPI_FINALIZE(mpierror)

END PROGRAM DATASET_BY_COL
