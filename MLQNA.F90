MODULE lqna_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."
 
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS
 

SUBROUTINE lqna(n, x, qn, qd)

!    =====================================================
!    Purpose: Compute Legendre functions Qn(x) and Qn'(x)
!    Input :  x  --- Argument of Qn(x) ( -1 ף x ף 1 )
!             n  --- Degree of Qn(x) ( n = 0,1,2,תתת )
!    Output:  QN(n) --- Qn(x)
!             QD(n) --- Qn'(x)
!             ( 1.0D+300 stands for infinity )
!    =====================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: qn(0:n)
REAL (dp), INTENT(OUT)  :: qd(0:n)

REAL (dp)  :: q0, q1, qf
INTEGER    :: k

IF (ABS(x) == 1.0D0) THEN
  DO  k = 0, n
    qn(k) = 1.0D+300
    qd(k) = -1.0D+300
  END DO
ELSE IF (ABS(x) < 1.0D0) THEN
  q0 = 0.5D0 * LOG((1.0D0+x)/(1.0D0-x))
  q1 = x * q0 - 1.0D0
  qn(0) = q0
  qn(1) = q1
  qd(0) = 1.0D0 / (1.0D0 - x*x)
  qd(1) = qn(0) + x * qd(0)
  DO  k = 2, n
    qf = ((2*k-1)*x*q1 - (k-1)*q0) / k
    qn(k) = qf
    qd(k) = (qn(k-1) - x*qf) * k / (1.0D0 - x*x)
    q0 = q1
    q1 = qf
  END DO
END IF
RETURN
END SUBROUTINE lqna
 
END MODULE lqna_func
 
 
 
PROGRAM mlqna
USE lqna_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:44

!    ======================================================
!    Purpose: This program computes the Legendre functions
!             Qn(x) and Qn'(x) using subroutine LQNA
!    Input :  x  --- Argument of Qn(x)  ( -1 ף x ף 1 )
!             n  --- Degree of Qn(x)  ( n = 0,1,תתת )
!    Output:  QN(n) --- Qn(x)
!             QD(n) --- Qn'(x)
!    Example:  x = 0.50
!              n        Qn(x)         Qn'(x)
!             ---------------------------------
!              0      .54930614     1.33333333
!              1     -.72534693     1.21597281
!              2     -.81866327     -.84270745
!              3     -.19865477    -2.87734353
!              4      .44017453    -2.23329085
!              5      .55508089     1.08422720
!    ======================================================

REAL (dp)  :: qn(0:100), qd(0:100), x
INTEGER    :: k, n

WRITE (*,*) '  Please enter Nmax and x '
READ (*,*) n, x
WRITE (*,5100) x
WRITE (*,*)
CALL lqna(n, x, qn, qd)
WRITE (*,*) '  n        Qn(x)         Qn''(X)'
WRITE (*,*) ' ---------------------------------'
DO  k = 0, n
  WRITE (*,5000) k, qn(k), qd(k)
END DO
STOP

5000 FORMAT (' ', i3, 2F15.8)
5100 FORMAT ('   x =', f5.2)
END PROGRAM mlqna
