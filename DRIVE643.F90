PROGRAM drive643
!                                   Driving program to read and
!                                   process the data file
USE Fisher_Exact
IMPLICIT NONE

INTEGER   :: i, ncol, nrow
REAL (dp) :: emin, expect, percnt, pre, prt
REAL (dp), ALLOCATABLE :: table(:,:)

!                                  Read table dimensions
10 READ (*,*) nrow, ncol, expect, percnt, emin
!                                  Terminate on 0, 0
IF (nrow == 0 .AND. ncol == 0) GO TO 9000
!                                  Read and output TABLE
ALLOCATE( table(nrow, ncol) )
WRITE (*,99998)
DO i=1, nrow
  READ (*, *) table(i, 1:ncol)
  WRITE (*,99997) table(i, 1:ncol)
  99997 FORMAT (' ', 10F7.0)
END DO

CALL fexact (nrow, ncol, table, nrow, expect, percnt, emin, prt, pre)

WRITE (*,99999) prt, pre
DEALLOCATE( table )
GO TO 10

99998 FORMAT (/' The contingency table for this problem is:')
99999 FORMAT ('  PRT = ', f8.6, '   PRE = ', f8.6)

9000 STOP
END PROGRAM drive643

