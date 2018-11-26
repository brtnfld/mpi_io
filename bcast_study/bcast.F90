#define DEBUG 0

PROGRAM case1

  USE mpi
  USE iso_c_binding
  use iso_fortran_env
  IMPLICIT NONE

  INTEGER, DIMENSION(mpi_status_size) :: wstatus
  INTEGER :: fh, k, i
  INTEGER :: ierr, rank
  
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset = 0
  DOUBLE PRECISION, DIMENSION(1:2) :: t, tg1, tg2, tg3
  DOUBLE PRECISION :: t1
  CHARACTER(len=32) :: arg
  character(len=128) filename
  INTEGER nprocs
  INTEGER ndim

  INTEGER, DIMENSION(:), ALLOCATABLE :: ndata

  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  IF( COMMAND_ARGUMENT_COUNT () .NE. 0) THEN
     CALL get_command_argument(1, arg)
     READ(arg(1:32),'(I32)') ndim
  ELSE
     arg(1:2) = "16"
     ndim = 16
  ENDIF
  filename = TRIM(arg)//".mpio"

  ALLOCATE(ndata(1:ndim))
  
!
! Create the file on one process
!
  IF(rank.EQ.0)THEN

     CALL MPI_File_open(MPI_COMM_SELF, filename,     &
          IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
          MPI_INFO_NULL, fh, ierr)

     DO k = 1, ndim
        ndata(k) = k
     ENDDO

     CALL MPI_File_write_all(fh, ndata, ndim, MPI_INTEGER, wstatus, ierr)
     CALL MPI_File_close(fh, ierr)

  ENDIF

  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
!
! Read the file
!  Option 0 -- All the processes read the same data
!  Option 1 -- Read all the metadata on one process, then Bcast all the metadata (i.e. using one Bcast).

  t(1:2) = 0.

! All the processes read the same data

  t1 = MPI_Wtime()

  CALL MPI_File_open(MPI_COMM_WORLD, filename,     &
       MPI_MODE_RDWR, MPI_INFO_NULL, fh, ierr)
     
  CALL MPI_File_set_view(fh, offset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr)

  CALL MPI_File_read_all(fh, ndata, ndim, MPI_INTEGER, wstatus, ierr );
#if DEBUG
  DO k = 0, nprocs
     IF(k.EQ.rank)THEN
        WRITE(*,'(I0,A)',ADVANCE="NO") rank,":"
        DO i = 1, ndim
           WRITE(*,'(X,I0)',ADVANCE="NO") ndata(i)
        ENDDO
        WRITE(*,"()")
     ENDIF
     CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
  ENDDO
#endif

  CALL MPI_File_close(fh, ierr)

  t(1) = (MPI_Wtime() - t1)

  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
! read-proc0-and-bcast

  t1 = MPI_Wtime()

  IF(rank.EQ.0)THEN

     CALL MPI_File_open(MPI_COMM_SELF, "datafile.mpio",     &
          MPI_MODE_RDWR, &
          MPI_INFO_NULL, fh, ierr)

     CALL MPI_File_set_view(fh, offset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr)
     CALL MPI_File_read(fh, ndata, ndim, MPI_INTEGER, wstatus, ierr )
     CALL MPI_File_close(fh, ierr)        
  ENDIF

  CALL MPI_BCAST(ndata, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#if DEBUG 
  DO k = 0, nprocs
     IF(k.EQ.rank)THEN
        WRITE(*,'(I0,A)',ADVANCE="NO") rank,":"
        DO i = 1, ndim
           WRITE(*,'(X,I0)',ADVANCE="NO") ndata(i)
        ENDDO
        WRITE(*,"()")
     ENDIF
     CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
  ENDDO
#endif
  t(2) = (MPI_Wtime() - t1);
   
  DEALLOCATE( ndata )

  CALL MPI_reduce(t, tg1, 2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_reduce(t, tg2, 2, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)
  CALL MPI_reduce(t, tg3, 2, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

  IF(rank.EQ.0)THEN
     WRITE(*,'(I0,X,I0,6F14.6)') nprocs, ndim, tg1(1)/DBLE(nprocs), tg2(1), tg3(1), &
          tg1(2)/DBLE(nprocs), tg2(2), tg3(2)
  ENDIF
    
  CALL MPI_Finalize(ierr)

END PROGRAM case1
