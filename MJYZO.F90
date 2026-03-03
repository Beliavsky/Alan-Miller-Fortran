MODULE jyzo_func
 
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
 

SUBROUTINE jyzo(n, nt, rj0, rj1, ry0, ry1)

!    ======================================================
!    Purpose: Compute the zeros of Bessel functions Jn(x),
!             Yn(x), and their derivatives
!    Input :  n  --- Order of Bessel functions ( n ó 101 )
!             NT --- Number of zeros (roots)
!    Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
!             RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
!             RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
!             RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
!    Routine called: JYNDD for computing Jn(x), Yn(x), and
!                    their first and second derivatives
!    ======================================================

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: nt
REAL (dp), INTENT(OUT)  :: rj0(nt)
REAL (dp), INTENT(OUT)  :: rj1(nt)
REAL (dp), INTENT(OUT)  :: ry0(nt)
REAL (dp), INTENT(OUT)  :: ry1(nt)

REAL (dp)  :: bjn, byn, djn, dyn, fjn, fyn, x, x0
INTEGER    :: l

IF (n <= 20) THEN
  x = 2.82141 + 1.15859 * n
ELSE
  x = n + 1.85576 * n ** 0.33333 + 1.03315 / n ** 0.33333
END IF
l = 0

10 x0 = x
CALL jyndd(n, x, bjn, djn, fjn, byn, dyn, fyn)
x = x - bjn / djn
IF (ABS(x-x0) > 1.0D-9) GO TO 10
l = l + 1
rj0(l) = x
x = x + 3.1416 + (0.0972 + 0.0679*n - 0.000354*n**2) / l
IF (l < nt) GO TO 10

IF (n <= 20) THEN
  x = 0.961587 + 1.07703 * n
ELSE
  x = n + 0.80861 * n ** 0.33333 + 0.07249 / n ** 0.33333
END IF
IF (n == 0) x = 3.8317
l = 0

20 x0 = x
CALL jyndd(n, x, bjn, djn, fjn, byn, dyn, fyn)
x = x - djn / fjn
IF (ABS(x-x0) > 1.0D-9) GO TO 20
l = l + 1
rj1(l) = x
x = x + 3.1416 + (0.4955 + 0.0915*n - 0.000435*n**2) / l
IF (l < nt) GO TO 20

IF (n <= 20) THEN
  x = 1.19477 + 1.08933 * n
ELSE
  x = n + 0.93158 * n ** 0.33333 + 0.26035 / n ** 0.33333
END IF
l = 0

30 x0 = x
CALL jyndd(n, x, bjn, djn, fjn, byn, dyn, fyn)
x = x - byn / dyn
IF (ABS(x-x0) > 1.0D-9) GO TO 30
l = l + 1
ry0(l) = x
x = x + 3.1416 + (0.312 + 0.0852*n - 0.000403*n**2) / l
IF (l < nt) GO TO 30

IF (n <= 20) THEN
  x = 2.67257 + 1.16099 * n
ELSE
  x = n + 1.8211 * n ** 0.33333 + 0.94001 / n ** 0.33333
END IF
l = 0

40 x0 = x
CALL jyndd(n, x, bjn, djn, fjn, byn, dyn, fyn)
x = x - dyn / fyn
IF (ABS(x-x0) > 1.0D-9) GO TO 40
l = l + 1
ry1(l) = x
x = x + 3.1416 + (0.197 + 0.0643*n - 0.000286*n**2) / l
IF (l < nt) GO TO 40

RETURN
END SUBROUTINE jyzo


SUBROUTINE jyndd(n, x, bjn, djn, fjn, byn, dyn, fyn)

!    =========================================================
!    Purpose: Compute Bessel functions Jn(x) and Yn(x), and
!             their first and second derivatives
!    Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 )
!             n   ---  Order of Jn(x) and Yn(x)
!    Output:  BJN ---  Jn(x)
!             DJN ---  Jn'(x)
!             FJN ---  Jn"(x)
!             BYN ---  Yn(x)
!             DYN ---  Yn'(x)
!             FYN ---  Yn"(x)
!    =========================================================

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: bjn
REAL (dp), INTENT(OUT)  :: djn
REAL (dp), INTENT(OUT)  :: fjn
REAL (dp), INTENT(OUT)  :: byn
REAL (dp), INTENT(OUT)  :: dyn
REAL (dp), INTENT(OUT)  :: fyn

REAL (dp)  :: bj(102), bs, by(102), e0, ec, f, f0, f1, s1, su
INTEGER    :: k, m, mt, nt

DO  nt = 1, 900
  mt = INT(0.5*LOG10(6.28*nt) - nt*LOG10(1.36*ABS(x)/nt))
  IF (mt > 20) EXIT
END DO
m = nt
bs = 0.0D0
f0 = 0.0D0
f1 = 1.0D-35
su = 0.0D0
DO  k = m, 0, -1
  f = 2.0D0 * (k+1) * f1 / x - f0
  IF (k <= n+1) bj(k+1) = f
  IF (k == 2*INT(k/2)) THEN
    bs = bs + 2.0D0 * f
    IF (k /= 0) su = su + (-1) ** (k/2) * f / k
  END IF
  f0 = f1
  f1 = f
END DO
DO  k = 0, n + 1
  bj(k+1) = bj(k+1) / (bs-f)
END DO
bjn = bj(n+1)
ec = 0.5772156649015329D0
e0 = 0.3183098861837907D0
s1 = 2.0D0 * e0 * (LOG(x/2.0D0)+ec) * bj(1)
f0 = s1 - 8.0D0 * e0 * su / (bs-f)
f1 = (bj(2)*f0-2.0D0*e0/x) / bj(1)
by(1) = f0
by(2) = f1
DO  k = 2, n + 1
  f = 2.0D0 * (k-1) * f1 / x - f0
  by(k+1) = f
  f0 = f1
  f1 = f
END DO
byn = by(n+1)
djn = -bj(n+2) + n * bj(n+1) / x
dyn = -by(n+2) + n * by(n+1) / x
fjn = (n*n/(x*x) - 1.0D0) * bjn - djn / x
fyn = (n*n/(x*x) - 1.0D0) * byn - dyn / x
RETURN
END SUBROUTINE jyndd
 
END MODULE jyzo_func
 
 
 
PROGRAM mjyzo
USE jyzo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:42

!    ==========================================================
!    Purpose: This program computes the zeros of Bessel
!             functions Jn(x), Yn(x), and their derivatives
!             using subroutine JYZO
!    Input :  n --- Order of Bessel functions ( n ó 100 )
!             NT --- Number of zeros
!    Output:  RJ0(m) --- m-th zero of Jn(x),  m=1,2,...,NT
!             RJ1(m) --- m-th zero of Jn'(x), m=1,2,...,NT
!             RY0(m) --- m-th zero of Yn(x),  m=1,2,...,NT
!             RY1(m) --- m-th zero of Yn'(x), m=1,2,...,NT
!    Example: n = 1, NT =5

!   Zeros of Bessel funcions Jn(x), Yn(x) and their derivatives
!                              ( n = 1 )
!    m       jnm           j'nm          ynm           y'nm
!   -----------------------------------------------------------
!    1     3.8317060     1.8411838     2.1971413     3.6830229
!    2     7.0155867     5.3314428     5.4296810     6.9415000
!    3    10.1734681     8.5363164     8.5960059    10.1234047
!    4    13.3236919    11.7060049    11.7491548    13.2857582
!    5    16.4706301    14.8635886    14.8974421    16.4400580
!    ==========================================================

REAL (dp)  :: rj0(101), rj1(101), ry0(101), ry1(101)
INTEGER    :: m, n, nt

WRITE (*,*) 'Please enter n and NT '
READ (*,*) n, nt
WRITE (*,*)
CALL jyzo(n, nt, rj0, rj1, ry0, ry1)
WRITE (*,5000)
WRITE (*,5100) n
WRITE (*,*) '  m       jnm           j''NM          YNM           y''NM'
WRITE (*,*) ' -----------------------------------------------------------'
DO  m = 1, nt
  WRITE (*,5200) m, rj0(m), rj1(m), ry0(m), ry1(m)
END DO
STOP

5000 FORMAT ('  Zeros of Bessel funcions Jn(x), Yn(x) and their derivatives')
5100 FORMAT (t31, '( n =', i2, ' )')
5200 FORMAT (' ', i3, 4F14.7)
END PROGRAM mjyzo
