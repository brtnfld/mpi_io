#define DEBUG 0


PROGRAM case1

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
  INTEGER :: fullbufsize
  INTEGER, DIMENSION(1:bufsize) :: buf
  INTEGER, ALLOCATABLE, DIMENSION(:) :: fullbuf
  DOUBLE PRECISION, DIMENSION(1:1) :: t
  DOUBLE PRECISION :: t1, t2, t3
  CHARACTER(len=1) :: argv
  LOGICAL :: exist
  CHARACTER(len=128) :: arg
  INTEGER k
  INTEGER nprocs

  INTEGER*4, DIMENSION(1:N) :: ndata
  INTEGER:: proc_cnt
  INTEGER :: incr, icnt, narg

  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  narg = 0
  argv(1:1)=""
  DO
     CALL get_command_argument(narg, arg)
     IF (LEN_TRIM(arg) == 0) EXIT
     argv(1:1) = arg(1:1)
     narg = narg + 1
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

  incr = 64
  offset = 0
  t(1) = 0.

  IF( narg .EQ. 1)THEN
     argv(1:1)="0"
     CALL MPI_File_open(MPI_COMM_WORLD, "datafile.mpio",     &
          MPI_MODE_RDWR, &
          MPI_INFO_NULL, fh, ierr)
     
     DO k = 1, N, incr
        CALL MPI_File_set_view(fh, offset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr)
        t1 = MPI_Wtime()
        CALL MPI_File_read_all(fh, buf, bufsize, MPI_INTEGER, wstatus, ierr );
        t(1) = t(1) + (MPI_Wtime() - t1);
        offset = offset + incr*sizeof(k)

#if DEBUG
        IF(rank.EQ.1)THEN
           DO i = 1, bufsize
              PRINT*,buf(i)
           ENDDO
        ENDIF
#endif
     ENDDO

     CALL MPI_File_close(fh, ierr)
  ENDIF

! Read on proc 0, and broadcast


! OPTION 1

  IF( argv .EQ. '1')THEN

     fullbufsize = N/incr*bufsize
     ALLOCATE(fullbuf(1:fullbufsize))

     IF(rank.EQ.0)THEN

        CALL MPI_File_open(MPI_COMM_SELF, "datafile.mpio",     &
             MPI_MODE_RDWR, &
             MPI_INFO_NULL, fh, ierr)

        icnt = 1
        DO k = 1, N, incr
           CALL MPI_File_set_view(fh, offset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr)
           t1 = MPI_Wtime()
           CALL MPI_File_read_all(fh, fullbuf(icnt), bufsize, MPI_INTEGER, wstatus, ierr );
           t(1) = t(1) + (MPI_Wtime() - t1);
           offset = offset + incr*sizeof(k)
           icnt = icnt + bufsize
        ENDDO
        CALL MPI_File_close(fh, ierr)        
     ENDIF
     
     t1 = MPI_Wtime()
     CALL MPI_BCAST(fullbuf, fullbufsize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     t(1) = t(1) + (MPI_Wtime() - t1);
#if DEBUG    
     IF(rank.EQ.1)THEN
        DO i = 1, fullbufsize
           PRINT*,fullbuf(i)
        ENDDO
     ENDIF
#endif     
     DEALLOCATE( fullbuf )
  ENDIF

! OPTION 2

  IF( argv .EQ. '2')THEN

     IF(rank.EQ.0)THEN
        CALL MPI_File_open(MPI_COMM_SELF, "datafile.mpio",     &
             MPI_MODE_RDWR, &
             MPI_INFO_NULL, fh, ierr)
     ENDIF

     DO k = 1, N, incr
        t1 = MPI_Wtime()
        IF(rank.EQ.0)THEN
           CALL MPI_File_set_view(fh, offset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr)
           CALL MPI_File_read_all(fh, buf, bufsize, MPI_INTEGER, wstatus, ierr );
        ENDIF
        CALL MPI_BCAST(buf, bufsize, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        t(1) = t(1) + (MPI_Wtime() - t1);
        offset = offset + incr*sizeof(k)
#if DEBUG
        IF(rank.EQ.1)THEN
           DO i = 1, bufsize
              PRINT*,buf(i)
           ENDDO
        ENDIF
#endif
     ENDDO
     
     IF(rank.EQ.0)THEN
        CALL MPI_File_close(fh, ierr)
     ENDIF
  ENDIF

  CALL MPI_reduce(t, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr);
  CALL MPI_reduce(t, t2, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr);
  CALL MPI_reduce(t, t3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr);

  IF(rank.EQ.0)THEN
    WRITE(*,'(I0,X,A,3F14.6)') nprocs, argv(1:1), t1/DBLE(nprocs), t2, t3
  ENDIF
    
  CALL MPI_Finalize(ierr)

END PROGRAM case1
