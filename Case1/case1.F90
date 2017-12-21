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

  INTEGER, PARAMETER :: i64 = SELECTED_INT_KIND(64)

  INTEGER(KIND=int64), PARAMETER :: MegaB = 2097152_int64

! 2**20
!  INTEGER*8, PARAMETER :: N = 1048576
! 2**27
!  INTEGER*8, PARAMETER :: N = 67108864 
! 2**30
!  INTEGER*8, PARAMETER :: N=1073741824
! 2**31
  INTEGER(KIND=int64), PARAMETER :: N=2147483648_int64
! 2**32
!  INTEGER(KIND=int64), PARAMETER :: N=4294967296
!
! 2**35
!  INTEGER(KIND=int64), PARAMETER :: 34359738368
! 2**36	
!  INTEGER(KIND=int64), PARAMETER :: 68719476736
! 2**37	
!  INTEGER(KIND=int64) ,PARAMETER :: 137438953472

! Write 4 large datasets 

  INTEGER, DIMENSION(mpi_status_size) :: wstatus
  INTEGER :: fh,i
  INTEGER, DIMENSION(1:4) :: message
  INTEGER :: ierr, rank, size
  
  INTEGER :: filetype, contig
  INTEGER (KIND=MPI_ADDRESS_KIND) :: extent
  INTEGER(KIND=MPI_OFFSET_KIND) :: disp, offset, expand_fs,sb_sz
  INTEGER, ALLOCATABLE, DIMENSION(:) :: buf
  DOUBLE PRECISION, DIMENSION(1:3) :: t
  DOUBLE PRECISION :: t1, t2, t3
  INTEGER :: bufsize
  CHARACTER(len=1) :: argv
  LOGICAL :: exist
  CHARACTER(len=128) :: arg
  INTEGER k
  INTEGER(KIND=MPI_OFFSET_KIND) f_sz
  INTEGER nranks, nprocs

  INTEGER, PARAMETER :: sz_superblock = 2048 ! 8,192 Bytes
  INTEGER*4, DIMENSION(1:sz_superblock) :: superblock
  INTEGER, DIMENSION(:), ALLOCATABLE :: mb_buf
  INTEGER:: proc_cnt

  CALL MPI_Init(ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)




  IF(rank.EQ.0)  CALL EXECUTE_COMMAND_LINE("rm -f datafile.mpio")

  CALL MPI_File_open(MPI_COMM_WORLD, "datafile.mpio",     &
       IOR(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
       MPI_INFO_NULL, fh, ierr)

  sb_sz = 0
  proc_cnt = 0
  offset = 0

! WRITE OUR PSEUDO SUPERBLOCK COLLECTIVELY,
! ONLY PROC 0 WRITES SOMETHING.
#if 0
  DO k = 1, sz_superblock
     superblock(k) = k
  ENDDO
  IF(rank.EQ.proc_cnt)THEN
     CALL MPI_File_write_all(fh, superblock, sz_superblock, MPI_INTEGER, wstatus, ierr)
  ELSE
     CALL MPI_File_write_all(fh, superblock, 0, MPI_INTEGER, wstatus, ierr)
  ENDIF
  offset = offset + sz_superblock*sizeof(k)

! WRITE THE RAW DATA AFTER THE SUPERBLOCK, ALL PROCESSES CONTRIBUTE
! TO WRITING THE DATA

  IF(MOD(N,nprocs).NE.0)THEN
     IF(rank.EQ.0) PRINT*, "ERROR: SIZE OF DATASET NOT DIVISABLE BY NUMBER OF PROCESSES"
     CALL MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
  ELSE
     bufsize = N/nprocs
  ENDIF

  CALL raw(fh, rank, nprocs, bufsize, offset)
#endif
! A) WRITE META-DATA INTER-DISPERSED BETWEEN RAW DATA WRITES, 4 SETS

  DO k = 1, 4

!    A.1) WRITE META-DATA, WRITE 512 INDIVIUAL WRITES, EACH INDIVIDUAL WRITE IS WRITTEN 
!         BY A SEPERATE PROCESS 


     DO i = 1, 512 ! writes of metadata  
        
        proc_cnt = proc_cnt + 1
        IF(proc_cnt.GT.nprocs-1) proc_cnt = 0
        
        bufsize = MegaB/sizeof(k)/512
        
        ALLOCATE(buf(1:bufsize))
        
        buf(:) = rank
        
        
        CALL MPI_File_set_view(fh, offset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr)
        
        IF(rank.EQ.proc_cnt)THEN
           CALL MPI_File_write_all(fh, buf, bufsize, MPI_INTEGER, wstatus, ierr)
        ELSE
           CALL MPI_File_write_all(fh, buf, 0, MPI_INTEGER, wstatus, ierr)
        ENDIF
        
        offset = offset + bufsize*sizeof(k)
        
        DEALLOCATE(buf)
        
     ENDDO
     
  ENDDO

  ! A.2) write raw data collectively

  bufsize = N/nprocs

  CALL raw(fh, rank, nprocs, bufsize, offset)

! EXPAND THE FILE

!  expand_fs = sizeof(i)*N + 524288 + sb_sz
  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
  t = 0.
  t1 = MPI_Wtime()

! (1) Expand using MPI IO
  CALL MPI_FILE_GET_SIZE(fh, f_sz, ierr) 
  t2 = MPI_Wtime() 
  CALL MPI_File_set_size(fh, f_sz, ierr)

  t(2) = MPI_Wtime() - t2;

  t3 = MPI_Wtime()
  CALL MPI_File_close(fh, ierr)
  t(3) = MPI_Wtime() - t3;
   
  t(1) = MPI_Wtime() - t1;

  CALL MPI_Allreduce(MPI_IN_PLACE, t, 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);

  IF (rank .EQ. (nprocs-1)) THEN
     INQUIRE(file="timing", exist=exist)
     WRITE(*, *) "TOTAL, MPI_File_set_size, MPI_File_close"
     PRINT*,t
     IF (exist) THEN
        OPEN(12, file="timing", status="old", position="append", action="write")
     ELSE
        OPEN(12, file="timing", status="new", action="write")
     END IF
     WRITE(12, *) t(1:3)
     CLOSE(12)
  ENDIF
  CALL MPI_Barrier(MPI_COMM_WORLD, ierr)

    
  CALL MPI_Finalize(ierr)

END PROGRAM case1

SUBROUTINE raw(fh, rank, nprocs, bufsize, offset)

    USE MPI

    IMPLICIT NONE
    INTEGER :: fh, ierr
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset
    INTEGER(KIND=MPI_OFFSET_KIND) :: LOCoffset
    INTEGER, ALLOCATABLE, DIMENSION(:) :: buf
    INTEGER k
    INTEGER, DIMENSION(mpi_status_size) :: wstatus
    INTEGER rank, nprocs
    INTEGER :: bufsize
    
    ALLOCATE(buf(1:bufsize))
    
    DO k = 1, bufsize
       buf(k) = k
    ENDDO
    
    LOCoffset = offset + rank*bufsize*sizeof(k)
    CALL MPI_File_set_view(fh, LOCoffset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr)
    CALL MPI_File_write_all(fh, buf, bufsize, MPI_INTEGER, wstatus, ierr)

    offset = offset + bufsize*nprocs*sizeof(k)

    DEALLOCATE(buf)
  
  END SUBROUTINE raw
