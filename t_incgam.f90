PROGRAM GTST

! Code converted using TO_F90 by Alan Miller
! Date: 2001-07-11  Time: 14:09:39

!     ALGORITHM 654, COLLECTED ALGORITHMS FROM ACM.
!     THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!     VOL. 13, NO. 3, P. 318.

! ----------------------------------------------------------------------
!  SAMPLE PROGRAM EMPLOYING GRATIO AND GAMINV.  GIVEN A AND X,
!  P AND Q ARE COMPUTED BY GRATIO. THEN FOR A, THE INVERSE OF
!  P AND Q, DENOTED BY XN, IS OBTAINED BY GAMINV.  D IS THE
!  RELATIVE DIFFERENCE BETWEEN X AND XN.

!  NO DATA IS READ. THE OUTPUT FOR THE PROGRAM IS WRITTEN ON
!  UNIT 6. THE FIRST STATMENT OF THIS TEXT MAY BE USED TO
!  BEGIN THE PROGRAM FOR THE CDC 6000-7000 SERIES COMPUTERS.

! ----------------------------------------------------------------------

USE Incomplete_Gamma
IMPLICIT NONE
REAL (dp)  :: a, d, p, q, x, x0, xn
INTEGER    :: i, ierr, l

WRITE (6,5000)

a = 0.1_dp
x0 = 0.0_dp
DO  l = 1, 10
  WRITE (6,5200)
  x = 0.1_dp
  DO  i = 1, 10
    CALL gratio(a, x, p, q, 0)
    CALL gaminv(a, xn, x0, p, q, ierr)
    d = ABS((x-xn)/x)
    WRITE (6,5100) a, x, p, q, xn, d, ierr
    x = x + 0.1_dp
  END DO
  a = a + 0.1_dp
END DO
STOP

5000 FORMAT ('    A     X', t20, 'P', t36, 'Q', t52, 'XN', t66, 'D       IERR')
5100 FORMAT (2F6.2, 3G16.6, G12.2, i5)
5200 FORMAT (' ')
END PROGRAM GTST
