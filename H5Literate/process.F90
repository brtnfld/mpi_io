PROGRAM main

  CHARACTER(LEN=180) :: filename
  CHARACTER(LEN=1) :: string
  INTEGER :: i, j, ii, icnt
  REAL, DIMENSION(100) :: timing
  REAL, DIMENSION(6) :: stat
  REAL :: aw,minw,maxw

  PRINT*,"Input file name"
  READ*, filename

  OPEN(10, file=filename)

  stat(1) = 0.
  stat(2) = HUGE(1.0)
  stat(3) = -HUGE(1.0)
  stat(4) = 0.
  stat(5) = HUGE(1.0)
  stat(6) = -HUGE(1.0)

  icnt = 0
  DO i = 1, 10
     READ(10,*, IOSTAT=iostate) iaux1, iaux2, timing(1)
     IF(iostate.NE.0) EXIT
     icnt = icnt + 1
     IF(icnt.LE.10)THEN
        stat(1) = stat(1) + timing(1) 
        stat(2) = MIN(stat(2),timing(1))
        stat(3) = MAX(stat(3),timing(1))
     ELSE
        stat(4) = stat(4) + timing(1) 
        stat(5) = MIN(stat(5),timing(1))
        stat(6) = MAX(stat(6),timing(1))
     ENDIF
  ENDDO

  OPEN(11, file=TRIM(filename)//".0.dat")
  WRITE(11,*) iaux1, stat(1)/10, stat(2:3)
  close(11)

  OPEN(11, file=TRIM(filename)//".1.dat")
  WRITE(11,*) iaux1,stat(4)/10, stat(5:6)
  close(11)

  

END PROGRAM main
