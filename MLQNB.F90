MODULE lqnb_func
 
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
 

SUBROUTINE lqnb(n, x, qn, qd)

!     ===================================================
!     Purpose: Compute Legendre functions Qn(x) & Qn'(x)
!     Input :  x  --- Argument of Qn(x)
!              n  --- Degree of Qn(x)  ( n = 0,1,2,תתת)
!     Output:  QN(n) --- Qn(x)
!              QD(n) --- Qn'(x)
!     ===================================================

INTEGER, INTENT(IN)      :: n
REAL (dp), INTENT(IN)    :: x
REAL (dp), INTENT(OUT)   :: qn(0:n)
REAL (dp), INTENT(OUT)   :: qd(0:n)

REAL (dp)  :: eps, q0, q1, qc1, qc2, qf, qf0, qf1, qf2, qr, x2
INTEGER    :: j, k, l, nl

eps = 1.0D-14
IF (ABS(x) == 1.0_dp) THEN
  DO  k = 0, n
    qn(k) = 1.0D+300
    qd(k) = 1.0D+300
  END DO
  RETURN
END IF

IF (x <= 1.021_dp) THEN
  x2 = ABS((1.0_dp+x)/(1.0_dp-x))
  q0 = 0.5_dp * LOG(x2)
  q1 = x * q0 - 1.0_dp
  qn(0) = q0
  qn(1) = q1
  qd(0) = 1.0_dp / (1.0_dp - x*x)
  qd(1) = qn(0) + x * qd(0)
  DO  k = 2, n
    qf = ((2*k-1)*x*q1 - (k-1)*q0) / k
    qn(k) = qf
    qd(k) = (qn(k-1) - x*qf) * k / (1.0_dp - x*x)
    q0 = q1
    q1 = qf
  END DO
ELSE
  qc2 = 1.0_dp / x
  DO  j = 1, n
    qc2 = qc2 * j / ((2*j+1)*x)
    IF (j == n-1) qc1 = qc2
  END DO
  DO  l = 0, 1
    nl = n + l
    qf = 1.0_dp
    qr = 1.0_dp
    DO  k = 1, 500
      qr = qr * (0.5_dp*nl+k-1) * (0.5_dp*(nl-1)+k) / ((nl+k-0.5_dp)*k*x*x)
      qf = qf + qr
      IF (ABS(qr/qf) < eps) EXIT
    END DO
    IF (l == 0) THEN
      qn(n-1) = qf * qc1
    ELSE
      qn(n) = qf * qc2
    END IF
  END DO
  qf2 = qn(n)
  qf1 = qn(n-1)
  DO  k = n, 2, -1
    qf0 = ((2*k-1)*x*qf1 - k*qf2) / (k-1)
    qn(k-2) = qf0
    qf2 = qf1
    qf1 = qf0
  END DO
  qd(0) = 1.0_dp / (1.0_dp - x*x)
  DO  k = 1, n
    qd(k) = k * (qn(k-1) - x*qn(k)) / (1.0_dp - x*x)
  END DO
END IF
RETURN
END SUBROUTINE lqnb
 
END MODULE lqnb_func
 
 
 
PROGRAM mlqnb
USE lqnb_func
IMPLICIT NONE


! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:44

!       ===============================================================
!       Purpose: This program computes the Legendre functions Qn(x)
!                and Qn'(x) using subroutine LQNB
!       Input :  x  --- Argument of Qn(x)
!                n  --- Degree of Qn(x)  ( n = 0,1,תתת)
!       Output:  QN(n) --- Qn(x)
!                QD(n) --- Qn'(x)
!       Examples:     x1 = 0.50,    x2 = 2.50

!       n      Qn(x1)        Qn'(x1)       Qn(x2)          Qn'(x2)
!     ----------------------------------------------------------------
!       0     .54930614    1.33333333   .42364893D+00  -.19047619D+00
!       1    -.72534693    1.21597281   .59122325D-01  -.52541546D-01
!       2    -.81866327    -.84270745   .98842555D-02  -.13109214D-01
!       3    -.19865477   -2.87734353   .17695141D-02  -.31202687D-02
!       4     .44017453   -2.23329085   .32843271D-03  -.72261513D-03
!       5     .55508089    1.08422720   .62335892D-04  -.16437427D-03
!       ===============================================================

REAL (dp)  :: qn(0:100), qd(0:100), x
INTEGER    :: k, n

WRITE (*,*) 'Please enter Nmax and x '
READ (*,*) n, x
WRITE (*,5200) x
WRITE (*,*)
WRITE (*,*) '  n          Qn(x)           Qn''(X)'
WRITE (*,*) '--------------------------------------'
CALL lqnb(n, x, qn, qd)
DO  k = 0, n
  IF (ABS(x) < 1.0_dp) THEN
    WRITE (*,5000) k, qn(k), qd(k)
  ELSE
    WRITE (*,5100) k, qn(k), qd(k)
  END IF
END DO
STOP

5000 FORMAT (' ', i3, 2F17.8)
5100 FORMAT (' ', i3, 2g17.8)
5200 FORMAT ('   x =', f5.2)
END PROGRAM mlqnb
