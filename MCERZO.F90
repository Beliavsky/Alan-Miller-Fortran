MODULE cerzo_func
 
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
 

SUBROUTINE cerzo(nt, zo)

!    ===============================================================
!    Purpose : Evaluate the complex zeros of error function erf(z)
!              using the modified Newton's iteration method
!    Input :   NT --- Total number of zeros
!    Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT
!    Routine called: CERF for computing erf(z) and erf'(z)
!    ===============================================================

INTEGER, INTENT(IN)        :: nt
COMPLEX (dp), INTENT(OUT)  :: zo(nt)

REAL (dp)     :: pu, pv, px, py, w = 0.0_dp, w0
COMPLEX (dp)  :: z, zd, zf, zfd, zgd, zp, zq, zw
INTEGER       :: i, it, j, nr
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp

DO  nr = 1, nt
  pu = SQRT(pi*(4.0_dp*nr - 0.5_dp))
  pv = pi * SQRT(2.0_dp*nr - 0.25_dp)
  px = 0.5 * pu - 0.5 * LOG(pv) / pu
  py = 0.5 * pu + 0.5 * LOG(pv) / pu
  z = CMPLX(px, py, KIND=dp)
  it = 0

  10 it = it + 1
  CALL cerf(z, zf, zd)
  zp = (1.0_dp, 0.0_dp)
  DO  i = 1, nr - 1
    zp = zp * (z - zo(i))
  END DO
  zfd = zf / zp
  zq = (0.0_dp, 0.0_dp)
  DO  i = 1, nr - 1
    zw = (1.0_dp, 0.0_dp)
    DO  j = 1, nr - 1
      IF (j /= i) THEN
        zw = zw * (z-zo(j))
      END IF
    END DO
    zq = zq + zw
  END DO
  zgd = (zd - zq*zfd) / zp
  z = z - zfd / zgd
  w0 = w
  w = ABS(z)
  IF (it <= 50 .AND. ABS((w-w0)/w) > 1.0D-11) GO TO 10
  zo(nr) = z
END DO
RETURN
END SUBROUTINE cerzo



SUBROUTINE cerf(z, cer, cder)

!    ==========================================================
!    Purpose: Compute complex Error function erf(z) & erf'(z)
!    Input:   z   --- Complex argument of erf(z)
!             x   --- Real part of z
!             y   --- Imaginary part of z
!    Output:  CER --- erf(z)
!             CDER --- erf'(z)
!    ==========================================================

COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cer
COMPLEX (dp), INTENT(OUT)  :: cder

REAL (dp)  :: c0, cs, eps, ei1, ei2, er, er0, er1, er2, eri, ERR, pi, r, ss, &
              w, w1, w2, x, x2, y
INTEGER    :: k, n

eps = 1.0D-12
pi = 3.141592653589793_dp
x = REAL(z, KIND=dp)
y = AIMAG(z)
x2 = x * x
IF (x <= 3.5_dp) THEN
  er = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 100
    r = r * x2 / (k + 0.5_dp)
    er = er + r
    IF (ABS(er-w) <= eps*ABS(er)) EXIT
    w = er
  END DO
  c0 = 2.0_dp / SQRT(pi) * x * EXP(-x2)
  er0 = c0 * er
ELSE
  er = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 12
    r = -r * (k - 0.5_dp) / x2
    er = er + r
  END DO
  c0 = EXP(-x2) / (x*SQRT(pi))
  er0 = 1.0_dp - c0 * er
END IF
IF (y == 0.0_dp) THEN
  ERR = er0
  eri = 0.0_dp
ELSE
  cs = COS(2.0_dp*x*y)
  ss = SIN(2.0_dp*x*y)
  er1 = EXP(-x2) * (1.0_dp-cs) / (2.0_dp*pi*x)
  ei1 = EXP(-x2) * ss / (2.0_dp*pi*x)
  er2 = 0.0_dp
  DO  n = 1, 100
    er2 = er2 + EXP(-.25_dp*n*n) / (n*n + 4.0_dp*x2) * (2.0_dp*x -   &
          2.0_dp*x*COSH(n*y)*cs + n*SINH(n*y)*ss)
    IF (ABS((er2-w1)/er2) < eps) EXIT
    w1 = er2
  END DO
  c0 = 2.0_dp * EXP(-x2) / pi
  ERR = er0 + er1 + c0 * er2
  ei2 = 0.0_dp
  DO  n = 1, 100
    ei2 = ei2 + EXP(-.25_dp*n*n) / (n*n + 4.0_dp*x2) * (2.0_dp*x*COSH(n*y)*ss  &
          + n*SINH(n*y)*cs)
    IF (ABS((ei2-w2)/ei2) < eps) EXIT
    w2 = ei2
  END DO
  eri = ei1 + c0 * ei2
END IF
cer = CMPLX(ERR, eri, KIND=dp)
cder = 2.0_dp / SQRT(pi) * EXP(-z*z)
RETURN
END SUBROUTINE cerf
 
END MODULE cerzo_func
 
 
 
PROGRAM mcerzo
USE cerzo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!       ===============================================================
!       Purpose : This program evaluates the complex zeros of error
!                 function erf(z) using subroutine CERZO
!       Input:    NT --- Total number of zeros
!       Example:  NT = 10

!    n     complex zeros of erf(z)     n     complex zeros of erf(z)
!   -------------------------------------------------------------------
!    1   1.450616163 + i 1.880943000   6   4.158998400 + i 4.435571444
!    2   2.244659274 + i 2.616575141   7   4.516319400 + i 4.780447644
!    3   2.839741047 + i 3.175628100   8   4.847970309 + i 5.101588043
!    4   3.335460735 + i 3.646174376   9   5.158767908 + i 5.403332643
!    5   3.769005567 + i 4.060697234  10   5.452192201 + i 5.688837437
!       ===============================================================

COMPLEX (dp)  :: zo(100)
INTEGER       :: i, nt

WRITE (*,*) 'Please Enter NT '
READ (*,*) nt
WRITE (*,5000) nt
CALL cerzo(nt, zo)
WRITE (*,*) '  *****    Please Wait !    *****'
WRITE (*,*)
WRITE (*,*) '  n        Complex zeros of erf(z)'
WRITE (*,*) '-------------------------------------'
DO  i = 1, nt
  WRITE (*,5100) i, zo(i)
END DO
STOP

5000 FORMAT ('  NT=', i3)
5100 FORMAT (' ', i3, '  ', f13.9, '  +i', f13.9)
END PROGRAM mcerzo
