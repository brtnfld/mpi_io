!
! Number of processes is assumed to be 1 or multiples of 2 (1,2,4,6,8)
!

     PROGRAM DATASET_BY_COL

     USE HDF5 ! This module contains all necessary modules 
        
     IMPLICIT NONE
! 2147483648
     include 'mpif.h'
     CHARACTER(LEN=10), PARAMETER :: filename = "sds_col.h5"  ! File name
     CHARACTER(LEN=8), PARAMETER :: dsetname = "IntArray" ! Dataset name

     INTEGER(HID_T) :: file_id, dcpl       ! File identifier 
     INTEGER(HID_T) :: dset_id       ! Dataset identifier 
     INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
     INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
     INTEGER(HID_T) :: plist_id      ! Property list identifier 

!     INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/65536,8192/) ! Dataset dimensions.

     INTEGER(HSIZE_T), DIMENSION(2) :: dimsf = (/65536,8182/) ! Dataset dimensions.

     INTEGER(HSIZE_T), DIMENSION(2) :: count, maxdims
     INTEGER(HSIZE_T), DIMENSION(2) :: chunk
     INTEGER(HSSIZE_T), DIMENSION(2) :: offset
     INTEGER, ALLOCATABLE, TARGET :: DATA (:,:)  ! Data to write
!     INTEGER*8, ALLOCATABLE, TARGET :: DATA(:,:)  ! Data to write
     INTEGER :: rank = 2 ! Dataset rank
     INTEGER :: i, j, k

     INTEGER :: error, error_n  ! Error flags
     !
     ! MPI definitions and calls.
     !
     INTEGER :: mpierror       ! MPI error flag
     INTEGER :: comm, comm2, info
     INTEGER :: mpi_size, mpi_rank
     INTEGER :: color, icode
     TYPE(C_PTR) :: f_ptr
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
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
     CALL h5pclose_f(plist_id, error)
     !
     ! Create the data space for the  dataset. 
     !
     maxdims(1:2) = H5S_UNLIMITED_F

!     CALL h5screate_simple_f(rank, dimsf, filespace, error, maxdims)
     CALL h5screate_simple_f(rank, dimsf, filespace, error)

     !
     ! Create the dataset with default properties.
     !
     chunk(1) = dimsf(1)/2 + 512
     chunk(2) = dimsf(2)

     CALL h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
!     CALL h5pset_chunk_f(dcpl, 2, chunk, error)

     CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, filespace, &
                      dset_id, error, dcpl)
     CALL h5sclose_f(filespace, error)
     CALL h5pclose_f(dcpl, error)
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
     ! Select hyperslab in the file.
     !
     CALL h5dget_space_f(dset_id, filespace, error)
     CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error)
     ! 
     ! Initialize data buffer with trivial data.
     !
     ALLOCATE ( data(count(1),count(2)))
!     data = mpi_rank + 10
     data = 99
     !
     ! Create property list for collective dataset write
     !
     CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
     CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
     
     !
     ! Write the dataset collectively. 
     !
     CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsf, error, &
                     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
     !
     ! Write the dataset independently. 
     !
     !
     ! Deallocate data buffer.
     !
     DEALLOCATE(data)

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
     CALL h5fclose_f(file_id, error)

     PRINT*,"READING"
     CALL MPI_BARRIER(MPI_COMM_WORLD, error)
     
     color = 0
     IF(mpi_rank .GE. 4) color = 1
     CALL MPI_Comm_split(MPI_COMM_WORLD, color, 1, comm2, error)
     comm2=comm

 !    IF(mpi_rank.GE.4)THEN

        CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        CALL h5pset_fapl_mpio_f(plist_id, comm2, info, error)

     !
     ! Create the file collectively.
     ! 
        CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error, access_prp = plist_id)
        CALL h5pclose_f(plist_id, error)

        CALL h5dopen_f(file_id, dsetname, dset_id, error)
        
        
        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        ALLOCATE(DATA(dimsf(1),dimsf(2)))
        f_ptr = C_LOC(DATA)
        CALL H5Dread_f(dset_id, H5T_NATIVE_INTEGER, f_ptr, error, H5S_ALL_F, H5S_ALL_F, plist_id)
!     CALL H5Dread_f(dset_id, H5T_STD_I64LE, f_ptr, error, H5S_ALL_F, H5S_ALL_F, plist_id)


        CALL H5Pclose_f(plist_id, error)

        DO i = 1, dimsf(1)
           DO j = 1, dimsf(2)
              IF(DATA(i,j).NE.99) THEN
                 PRINT*,"Wrong value", DATA(i,j)
                 CALL MPI_ABORT(comm2, icode, error)
              ENDIF      
           ENDDO
        ENDDO
!!$
!!$        DO k = 1, mpi_size
!!$           
!!$           IF(mpi_rank.EQ.k-1)THEN
!!$              DO i = 1, dimsf(1)
!!$                 WRITE(*,'(I0,a)',ADVANCE="no") mpi_rank,":["
!!$                 DO j = 1, dimsf(2)
!!$                    WRITE(*,'(X, I0)',ADVANCE="no") DATA(i,j)
!!$                 ENDDO
!!$                 WRITE(*,'(a)') "]"
!!$              ENDDO
!!$           ENDIF
!!$           CALL MPI_BARRIER(comm2, error)
!!$
!!$        ENDDO

        CALL h5dclose_f(dset_id, error)

        
        CALL h5fclose_f(file_id, error)

        !
        ! Close FORTRAN predefined datatypes.
        !

!     ENDIF

     CALL h5close_f(error)

     CALL MPI_FINALIZE(mpierror)

     END PROGRAM DATASET_BY_COL
