!
! Example how to write and read a 1D array in parallel
!
! Number of processes is assumed to be 1 or multiples of 2 (1,2,4,6,8)
!
! The program can be compiled with "h5pfc" wrapper included as part
! of the HDF5 Library.
!

PROGRAM main
  
  USE ISO_C_BINDING
  USE HDF5 ! This module contains all necessary modules 
        
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  CHARACTER(LEN=10), PARAMETER :: filename = "sds.h5"  ! File name
  CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name
  CHARACTER(LEN=6), PARAMETER :: groupname = "group1"  ! group name
  
  INTEGER, PARAMETER :: rank = 1  ! Dataset rank 
  INTEGER, PARAMETER :: real_kind_15 = KIND(1.d0)

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: group_id      ! Group identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier 
  INTEGER(HID_T) :: dset_id1      ! Dataset identifier 
  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
  INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
  INTEGER(HID_T) :: plist_id      ! Property list identifier 

  INTEGER(HSIZE_T), DIMENSION(1:rank) :: dimsf = (/1024/) ! Dataset dimensions.

  INTEGER(HSIZE_T), DIMENSION(1:rank)  :: count  
  INTEGER(HSSIZE_T), DIMENSION(1:rank) :: offset 
  REAL(real_kind_15), ALLOCATABLE, TARGET :: DATA (:)   ! Data to write
  REAL(real_kind_15), ALLOCATABLE, TARGET :: data_r(:)  ! Data to read

  TYPE(C_PTR) :: f_ptr

  INTEGER :: error, error_n  ! Error flags
  !
  ! MPI definitions and calls.
  !
  INTEGER :: mpierror       ! MPI error flag
  INTEGER :: comm, info
  INTEGER :: mpi_size, mpi_rank
  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL
  CALL MPI_INIT(mpierror)
  CALL MPI_COMM_SIZE(comm, mpi_size, mpierror)
  CALL MPI_COMM_RANK(comm, mpi_rank, mpierror) 
  !
  ! Initialize FORTRAN predefined datatypes
  !
  CALL h5open_f(error) 

  ! __          _______  _____ _______ _____ _   _  _____ 
  ! \ \        / /  __ \|_   _|__   __|_   _| \ | |/ ____|
  !  \ \  /\  / /| |__) | | |    | |    | | |  \| | |  __ 
  !   \ \/  \/ / |  _  /  | |    | |    | | | . ` | | |_ |
  !    \  /\  /  | | \ \ _| |_   | |   _| |_| |\  | |__| |
  !     \/  \/   |_|  \_\_____|  |_|  |_____|_| \_|\_____|

  ! 
  ! Setup file access property list with parallel I/O access.
  !
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
  !
  ! Create the file collectively.
  ! 
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  CALL h5pclose_f(plist_id, error)

  CALL h5gcreate_f(file_id, groupname, group_id, error)
  !
  ! Create the data space for the dataset. 
  !
  CALL h5screate_simple_f(rank, dimsf, filespace, error)
  !
  ! Create the dataset with default properties.
  !
  CALL h5dcreate_f(group_id, dsetname, h5kind_to_type(real_kind_15,H5_REAL_KIND), filespace, dset_id, error)
  CALL h5sclose_f(filespace, error)
  !
  ! Each process defines dataset in memory and writes it to the hyperslab
  ! in the file. 
  !
  count(1) = dimsf(1)/mpi_size
  offset(1) = mpi_rank * count(1) 
  CALL h5screate_simple_f(rank, count, memspace, error) 
  ! 
  ! Select hyperslab in the file.
  !
  CALL h5dget_space_f(dset_id, filespace, error)
  CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)
  ! 
  ! Initialize data buffer with trivial data.
  !
  ALLOCATE ( DATA(count(1)))
  DATA = mpi_rank + 1
  !
  ! Create property list for collective dataset write
  !
  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)   
  !
  ! Write the dataset collectively. 
  !

  f_ptr = C_LOC(data(1))

  CALL h5dwrite_f(dset_id, h5kind_to_type(real_kind_15,H5_REAL_KIND), f_ptr, error, &
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
  ! alternative:
  ! Write the dataset independently. 
  !
  ! CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsfi, error, &
  !                     file_space_id = filespace, mem_space_id = memspace)
  !
  ! Deallocate data buffer.
  !
  DEALLOCATE(DATA)

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

  ! Close the group
  CALL h5gclose_f(group_id,error)
  
  !
  ! Close the file.
  !
  CALL h5fclose_f(file_id, error)

  !  _____  ______          _____ _____ _   _  _____ 
  ! |  __ \|  ____|   /\   |  __ \_   _| \ | |/ ____|
  ! | |__) | |__     /  \  | |  | || | |  \| | |  __ 
  ! |  _  /|  __|   / /\ \ | |  | || | | . ` | | |_ |
  ! | | \ \| |____ / ____ \| |__| || |_| |\  | |__| |
  ! |_|  \_\______/_/    \_\_____/_____|_| \_|\_____|
 
  ! 
  ! Setup file access property list with parallel I/O access.
  !
  CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
  !
  ! Open the file collectively.
  !
  CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
  CALL h5pclose_f(plist_id, error)
  !
  ! Open the group and dataset.
  !
  CALL h5gopen_f(file_id, groupname, group_id, error)
  CALL h5dopen_f(group_id, dsetname, dset_id, error)

  !
  ! Get dataset's dataspace identifier.
  !
  CALL h5dget_space_f(dset_id, filespace, error)

  ! Create the memory data space
  count(1) = dimsf(1)/mpi_size
  offset(1) = mpi_rank * count(1)
  CALL h5screate_simple_f(rank, count, memspace, error) 

  !
  ! Select hyperslab in the dataset.
  !
  count(1) = dimsf(1)/mpi_size
  offset(1) = mpi_rank * count(1) 
  CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, error)

  CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  ! Read collectively
  ! CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  ! Read independently
  CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error) 

  !
  ! Read data from hyperslab in the file into the hyperslab in
  ! memory and display, independently (default)
  !
  ALLOCATE (data_r(1:count(1)))

  f_ptr = C_LOC(data_r(1))
  CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

  IF(mpi_rank.EQ.mpi_size-1)THEN
     WRITE(*,'(A,I0,A,I0,A,*(X,f7.1))')"data(",offset,":",offset+count,") =", data_r
  ENDIF

  !
  ! Close the dataspace for the dataset.
  !
  CALL h5sclose_f(filespace, error)

  !
  ! Close the memoryspace.
  !
  CALL h5sclose_f(memspace, error)

  !
  ! Close the group, dataset and property list.
  !
  CALL h5pclose_f(plist_id, error)
  CALL h5gclose_f(group_id, error)
  CALL h5dclose_f(dset_id, error)

  !
  ! Close the file.
  !
  CALL h5fclose_f(file_id, error)
  ! 
  !
  ! Close FORTRAN predefined datatypes.
  !
  CALL h5close_f(error)
  
  CALL MPI_FINALIZE(mpierror)

END PROGRAM MAIN
