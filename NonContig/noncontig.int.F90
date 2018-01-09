PROGRAM noncontig
    use mpi
    use iso_c_binding
    implicit none

    INTERFACE
       INTEGER FUNCTION ctrunc(a) BIND(C,NAME="ctrunc")
         USE ISO_C_BINDING
         USE mpi
         INTEGER(KIND=MPI_OFFSET_KIND), VALUE :: a
       END FUNCTION ctrunc
    END INTERFACE

    INTEGER, DIMENSION(mpi_status_size) :: wstatus
    INTEGER :: fh,i
    INTEGER, DIMENSION(1:4) :: message
    INTEGER :: ierr, rank, size
    INTEGER, PARAMETER :: N=1073741824 ! 2147483648 
    
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

    CALL MPI_Init(ierr)
    CALL MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    k = 0
    argv=""
    DO
       CALL get_command_argument(k, arg)
       IF (LEN_TRIM(arg) == 0) EXIT
       argv(1:1) = arg(1:1)
       k = k + 1
    END DO
    
    call MPI_File_open(MPI_COMM_WORLD, "datafile",     &
                       ior(MPI_MODE_CREATE,MPI_MODE_WRONLY), &
                       MPI_INFO_NULL, fh, ierr)

! WRITE OUR PSEUDO SUPERBLOCK COLLECTIVELY

    message = rank

    CALL MPI_Type_contiguous(1, MPI_INTEGER, contig, ierr)

    extent = size*4
    CALL MPI_Type_create_resized(contig, 0_MPI_ADDRESS_KIND, extent, filetype, ierr)
    CALL MPI_Type_commit(filetype, ierr)

    disp = rank*sizeof(i)
    CALL MPI_File_set_view(fh, disp, MPI_INTEGER, filetype, "native", MPI_INFO_NULL, ierr) 

    CALL MPI_File_write_all(fh, message, 4, MPI_INTEGER, wstatus, ierr)

! WRITE THE RAW DATA AFTER THE SUPERBLOCK
    bufsize = N/size

    ALLOCATE(buf(1:bufsize))

    buf = rank

    sb_sz = 16*size

    offset = sb_sz + rank*bufsize*sizeof(i)
    CALL MPI_File_set_view(fh, offset, MPI_INTEGER, MPI_INTEGER, "native", MPI_INFO_NULL, ierr);
    CALL MPI_File_write_all(fh, buf, bufsize, MPI_INTEGER, wstatus, ierr);
! EXPAND THE FILE

    expand_fs = 0
!    expand_fs = sizeof(i)*N + 524288 + sb_sz
    CALL MPI_Barrier(MPI_COMM_WORLD, ierr)
    CALL MPI_FILE_GET_SIZE(fh, f_sz, ierr)

    expand_fs = f_sz

    t = 0.
    t1 = MPI_Wtime()

! (1) Expand using MPI IO
    IF(argv .EQ. '1')THEN
    
       t2 = MPI_Wtime()
       CALL MPI_File_set_size(fh, expand_fs, ierr)
       t(2) = MPI_Wtime() - t2;
    ENDIF

    t3 = MPI_Wtime()
    CALL MPI_File_close(fh, ierr)
    t(3) = MPI_Wtime() - t3;
   
 ! or (2) Expand using POSIX
    IF( argv .EQ. '2' .AND. rank .EQ. 0)THEN
       t2 = MPI_Wtime()
       i = ctrunc(expand_fs)
       t(2) = MPI_Wtime() - t2;
    ENDIF

    t(1) = MPI_Wtime() - t1;

    CALL MPI_Allreduce(MPI_IN_PLACE, t, 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr);

    IF (rank .EQ. (size-1)) THEN
       INQUIRE(file="timing", exist=exist)
       WRITE(*, *) "TOTAL,    MPI_File_set_size,    MPI_File_close"
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

    DEALLOCATE(buf)
    
    call MPI_Finalize(ierr)

end program noncontig
