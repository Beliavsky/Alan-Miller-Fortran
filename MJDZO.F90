MODULE jdzo_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Zhang & Jin's errata applied to the lower limits of 2 arrays (ZO & ZOC).
! Latest revision - 5 January 2002
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Variables x, zoc & p were used before values had been assigned to them.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE jdzo(nt, n, m, p, zo)

!       ===========================================================
!       Purpose: Compute the zeros of Bessel functions Jn(x) and
!                Jn'(x), and arrange them in the order of their
!                magnitudes
!       Input :  NT    --- Number of total zeros ( NT ף 1200 )
!       Output:  ZO(L) --- Value of the L-th zero of Jn(x)
!                          and Jn'(x)
!                N(L)  --- n, order of Jn(x) or Jn'(x) associated
!                          with the L-th zero
!                M(L)  --- m, serial number of the zeros of Jn(x)
!                          or Jn'(x) associated with the L-th zero
!                          ( L is the serial number of all the
!                            zeros of Jn(x) and Jn'(x) )
!                P(L)  --- TM or TE, a code for designating the
!                          zeros of Jn(x)  or Jn'(x).
!                          In the waveguide applications, the zeros
!                          of Jn(x) correspond to TM modes and
!                          those of Jn'(x) correspond to TE modes
!       Routine called:    BJNDD for computing Jn(x), Jn'(x) and
!                          Jn''(x)
!       =============================================================

INTEGER, INTENT(IN)             :: nt
INTEGER, INTENT(OUT)            :: n(1400)
INTEGER, INTENT(OUT)            :: m(1400)
CHARACTER (LEN=1), INTENT(OUT)  :: p(1400)
REAL (dp), INTENT(OUT)          :: zo(0:1400)

CHARACTER (LEN=4)  :: p1(70)
INTEGER    :: n1(70), m1(70)
REAL (dp)  :: zoc(0:70), bj(101), dj(101), fj(101)
INTEGER    :: i, j, k, l, l0, l1, l2, mm, nm
REAL (dp)  :: x, x0, x1, x2, xm

IF (nt < 600) THEN
  xm = -1.0 + 2.248485 * nt ** 0.5 - .0159382 * nt + 3.208775E-4 * nt ** 1.5
  nm = INT(14.5+.05875*nt)
  mm = INT(.02*nt) + 6
ELSE
  xm = 5.0 + 1.445389 * nt ** .5 + .01889876 * nt - 2.147763E-4 * nt ** 1.5
  nm = INT(27.8 + .0327*nt)
  mm = INT(.01088*nt) + 10
END IF
l0 = 0
p = ' '
DO  i = 1, nm
  x1 = .407658 + .4795504 * (i-1) ** .5 + .983618 * (i-1)
  x2 = 1.99535 + .8333883 * (i-1) ** .5 + .984584 * (i-1)
  x = x1
  l1 = 0
  DO  j = 1, mm
    IF (i /= 1 .OR. j /= 1) THEN
      x = x1
      10 CALL bjndd(i, x, bj, dj, fj)
      x0 = x
      x = x - dj(i) / fj(i)
      IF (x1 > xm) GO TO 20
      IF (ABS(x-x0) > 1.0D-10) GO TO 10
    END IF
    l1 = l1 + 1
    n1(l1) = i - 1
    m1(l1) = j
    IF (i == 1) m1(l1) = j - 1
    p1(l1) = 'TE'
    zoc(l1) = x
    IF (i <= 15) THEN
      x1 = x + 3.057 + .0122 * (i-1) + (1.555+.41575*(i-1)) / (j+1 ) ** 2
    ELSE
      x1 = x + 2.918 + .01924 * (i-1) + (6.26+.13205*(i-1)) / (j+1 ) ** 2
    END IF

    20 x = x2
    30 CALL bjndd(i, x, bj, dj, fj)
    x0 = x
    x = x - bj(i) / dj(i)
    IF (x <= xm) THEN
      IF (ABS(x-x0) > 1.0D-10) GO TO 30
      l1 = l1 + 1
      n1(l1) = i - 1
      m1(l1) = j
      p1(l1) = 'TM'
      zoc(l1) = x
      IF (i <= 15) THEN
        x2 = x + 3.11 + .0138 * (i-1) + (.04832 + .2804*(i-1)) / (j+1) ** 2
      ELSE
        x2 = x + 3.001 + .0105 * (i-1) + (11.52 + .48525*(i-1)) / (j+3) ** 2
      END IF
    END IF
  END DO
  l = l0 + l1
  l2 = l
  zoc(0) = 0.0_dp

  50 IF (l0 == 0) THEN
    DO  k = 1, l
      zo(k) = zoc(k)
      n(k) = n1(k)
      m(k) = m1(k)
      p(k) = p1(k)
    END DO
    l1 = 0
  ELSE IF (l0 /= 0) THEN
    IF (zo(l0) >= zoc(l1)) THEN
      zo(l0+l1) = zo(l0)
      n(l0+l1) = n(l0)
      m(l0+l1) = m(l0)
      p(l0+l1) = p(l0)
      l0 = l0 - 1
    ELSE
      zo(l0+l1) = zoc(l1)
      n(l0+l1) = n1(l1)
      m(l0+l1) = m1(l1)
      p(l0+l1) = p1(l1)
      l1 = l1 - 1
    END IF
  END IF
  IF (l1 /= 0) GO TO 50
  l0 = l2
END DO
RETURN
END SUBROUTINE jdzo



SUBROUTINE bjndd(n, x, bj, dj, fj)

!       =====================================================
!       Purpose: Compute Bessel functions Jn(x) and their
!                first and second derivatives ( n= 0,1,תתת )
!       Input:   x ---  Argument of Jn(x)  ( x ע 0 )
!                n ---  Order of Jn(x)
!       Output:  BJ(n+1) ---  Jn(x)
!                DJ(n+1) ---  Jn'(x)
!                FJ(n+1) ---  Jn"(x)
!       =====================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: bj(101)
REAL (dp), INTENT(OUT)     :: dj(101)
REAL (dp), INTENT(OUT)     :: fj(101)

REAL (dp)  :: bs, f, f0, f1
INTEGER    :: k, m, mt, nt

DO  nt = 1, 900
  mt = INT(0.5*LOG10(6.28*nt) - nt*LOG10(1.36*ABS(x)/nt))
  IF (mt > 20) EXIT
END DO
m = nt
bs = 0.0_dp
f0 = 0.0_dp
f1 = 1.0D-35
DO  k = m, 0, -1
  f = 2.0_dp * (k+1) * f1 / x - f0
  IF (k <= n) bj(k+1) = f
  IF (k == 2*INT(k/2)) bs = bs + 2.0_dp * f
  f0 = f1
  f1 = f
END DO
DO  k = 0, n
  bj(k+1) = bj(k+1) / (bs-f)
END DO
dj(1) = -bj(2)
fj(1) = -1.0_dp * bj(1) - dj(1) / x
DO  k = 1, n
  dj(k+1) = bj(k) - k * bj(k+1) / x
  fj(k+1) = (k*k/(x*x) - 1.0_dp) * bj(k+1) - dj(k+1) / x
END DO
RETURN
END SUBROUTINE bjndd
 
END MODULE jdzo_func
 
 
 
PROGRAM mjdzo
USE jdzo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:41

!       =============================================================
!       Purpose: This program computes the zeros of Bessel functions
!                Jn(x) and Jn'(x), and arranges them in the order
!                of their values
!       Input :  NT    --- Number of total zeros ( NT ף 1200 )
!       Output:  ZO(L) --- Value of the L-th zero of Jn(x) and Jn'(x)
!                N(L)  --- n, order of Jn(x) or Jn'(x) associated
!                          with the L-th zero
!                M(L)  --- m, serial number of the zeros of Jn(x)
!                          or Jn'(x) associated with the L-th zero
!                          ( L is the serial number of all the
!                            zeros of Jn(x) and Jn'(x) )
!                P(L)  --- TM or TE, a code for designating the
!                          zeros of Jn(x) or Jn'(x)
!                          In the waveguide applications, the zeros
!                          of Jn(x) correspond to TM modes and those
!                          of Jn'(x) correspond to TE modes.
!       =============================================================

CHARACTER (LEN=4)  :: p(1400)
INTEGER            :: n(1400), m(1400)
REAL (dp)          :: zo(0:1400)
INTEGER            :: j1, j2, k, k0, ks, nt

WRITE (*,*) 'NT=? '
READ (*,*) nt
WRITE (*,5000) nt
WRITE (*,5200)
CALL jdzo(nt, n, m, p, zo)
WRITE (*,*)
ks = nt / 101 + 1
DO  k0 = 1, ks
  WRITE(*,*) ' Table           Zeros of Bessel functions Jn(x) and Jn''(X)'
  WRITE(*,*)
  WRITE(*,*) ' --------------------------------------------------------------------'
  DO  k = 1, 50
    j1 = 100 * (k0-1) + k + 1
    j2 = j1 + 50
    IF (j1 <= nt+1 .AND. j2 <= nt+1) THEN
      WRITE (*,5100) j1-1, p(j1), n(j1), m(j1), zo(j1), j2-1,  &
                     p(j2), n(j2), m(j2), zo(j2)
    ELSE IF (j1 <= nt+1 .AND. j2 > nt+1) THEN
      WRITE (*,5100) j1-1, p(j1), n(j1), m(j1), zo(j1)
    END IF
  END DO
  WRITE (*,*) ' --------------------------------------------------------------------'
  WRITE (*,5300)
END DO
STOP

5000 FORMAT (' Total number of the zeros:', i5)
5100 FORMAT (' ', i4, '   ', a2, i4, ' -', i2, f14.8, '   |  ', i4, '   ',  &
             a2, i4, ' -', i2, f14.8)
5200 FORMAT (t16, '***  Please wait !  The program is running  ***')
5300 FORMAT (/)
END PROGRAM mjdzo
