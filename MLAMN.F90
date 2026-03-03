MODULE lamn_func
 
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


SUBROUTINE lamn(n, x, nm, bl, dl)

!    =========================================================
!    Purpose: Compute lambda functions and their derivatives
!    Input:   x --- Argument of lambda function
!             n --- Order of lambda function
!    Output:  BL(n) --- Lambda function of order n
!             DL(n) --- Derivative of lambda function
!             NM --- Highest order computed
!    Routines called:
!             MSTA1 and MSTA2 for computing the start
!             point for backward recurrence
!    =========================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(OUT)    :: nm
REAL (dp), INTENT(OUT)  :: bl(0:n)
REAL (dp), INTENT(OUT)  :: dl(0:n)

REAL (dp)  :: bg, bk, bs, f, f0, f1, r, r0, uk, x2
INTEGER    :: i, k, m

nm = n
IF (ABS(x) < 1.0D-100) THEN
  DO  k = 0, n
    bl(k) = 0.0D0
    dl(k) = 0.0D0
  END DO
  bl(0) = 1.0D0
  dl(1) = 0.5D0
  RETURN
END IF
IF (x <= 12.0D0) THEN
  x2 = x * x
  DO  k = 0, n
    bk = 1.0D0
    r = 1.0D0
    DO  i = 1, 50
      r = -0.25D0 * r * x2 / (i*(i+k))
      bk = bk + r
      IF (ABS(r) < ABS(bk)*1.0D-15) EXIT
    END DO
    bl(k) = bk
    IF (k >= 1) dl(k-1) = -0.5D0 * x / k * bk
  END DO
  uk = 1.0D0
  r = 1.0D0
  DO  i = 1, 50
    r = -0.25D0 * r * x2 / (i*(i+n+1.0D0))
    uk = uk + r
    IF (ABS(r) < ABS(uk)*1.0D-15) EXIT
  END DO
  dl(n) = -0.5D0 * x / (n+1.0D0) * uk
  RETURN
END IF
IF (n == 0) nm = 1
m = msta1(x, 200)
IF (m < nm) THEN
  nm = m
ELSE
  m = msta2(x, nm, 15)
END IF
bs = 0.0D0
f0 = 0.0D0
f1 = 1.0D-100
DO  k = m, 0, -1
  f = 2.0D0 * (k+1.0D0) * f1 / x - f0
  IF (k <= nm) bl(k) = f
  IF (k == 2*INT(k/2)) bs = bs + 2.0D0 * f
  f0 = f1
  f1 = f
END DO
bg = bs - f
DO  k = 0, nm
  bl(k) = bl(k) / bg
END DO
r0 = 1.0D0
DO  k = 1, nm
  r0 = 2.0D0 * r0 * k / x
  bl(k) = r0 * bl(k)
END DO
dl(0) = -0.5D0 * x * bl(1)
DO  k = 1, nm
  dl(k) = 2.0D0 * k / x * (bl(k-1)-bl(k))
END DO
RETURN
END SUBROUTINE lamn


FUNCTION msta1(x, mp) RESULT(fn_val)

!    ===================================================
!    Purpose: Determine the starting point for backward
!             recurrence such that the magnitude of
!             Jn(x) at that point is about 10^(-MP)
!    Input :  x     --- Argument of Jn(x)
!             MP    --- Value of magnitude
!    Output:  MSTA1 --- Starting point
!    ===================================================

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: mp
INTEGER                :: fn_val

REAL (dp)  :: a0, f, f0, f1
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
n0 = INT(1.1*a0) + 1
f0 = envj(n0,a0) - mp
n1 = n0 + 5
f1 = envj(n1,a0) - mp
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dp - f0/f1)
  f = envj(nn,a0) - mp
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn
RETURN
END FUNCTION msta1



FUNCTION msta2(x, n, mp) RESULT(fn_val)

!    ===================================================
!    Purpose: Determine the starting point for backward
!             recurrence such that all Jn(x) has MP
!             significant digits
!    Input :  x  --- Argument of Jn(x)
!             n  --- Order of Jn(x)
!             MP --- Significant digit
!    Output:  MSTA2 --- Starting point
!    ===================================================

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: n
INTEGER, INTENT(IN)    :: mp
INTEGER                :: fn_val

REAL (dp)  :: a0, ejn, f, f0, f1, hmp, obj
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
hmp = 0.5_dp * mp
ejn = envj(n, a0)
IF (ejn <= hmp) THEN
  obj = mp
  n0 = INT(1.1*a0)
ELSE
  obj = hmp + ejn
  n0 = n
END IF
f0 = envj(n0,a0) - obj
n1 = n0 + 5
f1 = envj(n1,a0) - obj
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dp - f0/f1)
  f = envj(nn, a0) - obj
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn + 10
RETURN
END FUNCTION msta2



FUNCTION envj(n, x) RESULT(fn_val)

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

fn_val = 0.5_dp * LOG10(6.28_dp*n) - n * LOG10(1.36_dp*x/n)
RETURN
END FUNCTION envj

END MODULE lamn_func
 
 
 
PROGRAM mlamn
USE lamn_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:43

!    ====================================================
!    Purpose: This program computes the lambda functions
!             and their derivatives using subroutine LAMN
!    Input:   x --- Argument of lambda function
!             n --- Order of lambda function
!                   ( n = 0,1,..., n ó 250 )
!    Output:  BL(n) --- Lambda function of order n
!             DL(n) --- Derivative of lambda function
!    Example: Nmax = 5,  x = 10.00

!              n       lambda(x)        lambda'(x)
!             ---------------------------------------
!              0    -.24593576D+00    -.43472746D-01
!              1     .86945492D-02    -.50926063D-01
!              2     .20370425D-01    -.46703503D-02
!              3     .28022102D-02     .10540929D-01
!              4    -.84327431D-02     .89879627D-02
!              5    -.89879627D-02     .55521954D-03
!    ====================================================

REAL (dp)  :: bl(0:250), dl(0:250), x
INTEGER    :: k, n, nm, ns

WRITE (*,*) '  Please enter n,x = ? '
READ (*,*) n, x
WRITE (*,5000) n, x
IF (n <= 10) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns '
  READ (*,*) ns
END IF
CALL lamn(n, x, nm, bl, dl)
WRITE (*,*)
WRITE (*,*) '  n       lambda(x)        lambda''(X)'
WRITE (*,*) ' ---------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5100) k, bl(k), dl(k)
END DO
STOP

5000 FORMAT (' N =', i4, '      x =', f8.2)
5100 FORMAT (' ', i3, 2g18.8)
END PROGRAM mlamn
