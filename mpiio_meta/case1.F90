PROGRAM case1

! Large chunks of data are interspersed with randomly scattered 
! small clusters of small (i.e. < 5 KiB) pieces of data totaling 1-2 MiB. 
! Large data is written in several collective writes.  
! Small data is also written collectively just before file close, with at 
! most one piece of small data per process.  File is truncated 
! to its current size at close via MPI_File_set_size()

  USE mpi
  USE iso_c_binding
  use iso_fortran_env
  IMPLICIT NONE

!  INTEGER(KIND=int64) ,PARAMETER :: N = 524288_int64
  INTEGER, PARAMETER :: N = 524288

  INTEGER, DIMENSION(mpi_status_size) :: wstatus
  INTEGER :: fh,i
  INTEGER :: ierr, rank
  
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset
  INTEGER, PARAMETER :: bufsize = 8
  INTEGER, DIMENSION(1:bufsize) :: buf
  DOUBLE PRECISION, DIMENSION(1:1) :: t
  DOUBLE PRECISION :: t1, t2, t3
  CHARACTER(len=1) :: argv
  LOGICAL :: exist
  CHARACTER(len=128) :: arg
  INTEGER k
  INTEGER nprocs

  INTEGER*4, DIMENSION(1:N) :: ndata
  INTEGER:: proc_cnt
  INTEGER :: incr

  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  k = 0
  argv=""
  DO
     CALL get_command_argument(k, arg)
     IF (LEN_TRIM(arg) == 0) EXIT
     argv(1:1) = arg(1:1)
     k = k + 1
  END DO

  IF(rank.EQ.0)THEN

     CALL MPI_File_open(MPI_COMM_SELF, "datafile.mpio",     &
          IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
          MPI_INFO_NULL, fh, ierr)

     DO k = 1, N
        ndata(k) = k
     ENDDO

     CALL MPI_File_write_all(fh, ndata, N, MPI_INTEGER, wstatus, ierr)
     CALL MPI_File_close(fh, ierr)

  ENDIF

  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

  CALL MPI_File_open(MPI_COMM_WORLD, "datafile.mpio",     &
       MPI_MODE_RDWR, &
       MPI_INFO_NULL, fh, ierr)

  incr = 64
  t(1) = 0.
  offset = 0
  DO k = 1, N, incr
    CALL MPI_File_set_view(fh, offset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr)
    t1 = MPI_Wtime()
    CALL MPI_File_read_all(fh, buf, bufsize, MPI_INT, wstatus, ierr );
    t(1) = t(1) + (MPI_Wtime() - t1);
    offset = offset + incr*sizeof(k)

#if 0
    IF(rank.EQ.1)THEN
       DO i = 1, bufsize
          PRINT*,buf(i)
       ENDDO
    ENDIF
#endif
  ENDDO

  CALL MPI_File_close(fh, ierr)
  
  CALL MPI_reduce(t, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr);
  CALL MPI_reduce(t, t2, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr);
  CALL MPI_reduce(t, t3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr);

  IF(rank.EQ.0)THEN
    WRITE(*,'(3F14.6)') t1/DBLE(nprocs), t2, t3
  ENDIF
    
  CALL MPI_Finalize(ierr)

END PROGRAM case1
