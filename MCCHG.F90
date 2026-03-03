MODULE cchg_func
 
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
 

SUBROUTINE cchg(a, b, z, chg)

!    ===================================================
!    Purpose: Compute confluent hypergeometric function
!             M(a,b,z) with real parameters a, b and a
!             complex argument z
!    Input :  a --- Parameter
!             b --- Parameter
!             z --- Complex argument
!    Output:  CHG --- M(a,b,z)
!    Routine called: GAMMA for computing gamma function
!    ===================================================

REAL (dp), INTENT(IN OUT)     :: a
REAL (dp), INTENT(IN)         :: b
COMPLEX (dp), INTENT(IN OUT)  :: z
COMPLEX (dp), INTENT(OUT)     :: chg

COMPLEX (dp)  :: cfac, chg1, chg2, chw, ci, cr, cr1, cr2, crg, cs1, cs2,  &
                 cy0, cy1, z0
REAL (dp)     :: a0, a1, ba, g1, g2, g3, phi, x, x0, y
INTEGER       :: i, j, k, la, m, n, nl, ns
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp

ci = (0.0_dp, 1.0_dp)
a0 = a
a1 = a
z0 = z
IF (b == 0.0 .OR. b == -INT(ABS(b))) THEN
  chg = (1.0D+300, 0.0_dp)
ELSE IF (a == 0.0_dp .OR. z == 0.0_dp) THEN
  chg = (1.0_dp, 0.0_dp)
ELSE IF (a == -1.0_dp) THEN
  chg = 1.0_dp - z / b
ELSE IF (a == b) THEN
  chg = EXP(z)
ELSE IF (a-b == 1.0_dp) THEN
  chg = (1.0_dp + z/b) * EXP(z)
ELSE IF (a == 1.0_dp .AND. b == 2.0_dp) THEN
  chg = (EXP(z) - 1.0_dp) / z
ELSE IF (a == INT(a) .AND. a < 0.0_dp) THEN
  m = INT(-a)
  cr = (1.0_dp,0.0_dp)
  chg = (1.0_dp,0.0_dp)
  DO  k = 1, m
    cr = cr * (a+k-1) / k / (b+k-1) * z
    chg = chg + cr
  END DO
ELSE
  x0 = REAL(z)
  IF (x0 < 0.0_dp) THEN
    a = b - a
    a0 = a
    z = -z
  END IF
  IF (a < 2.0_dp) nl = 0
  IF (a >= 2.0_dp) THEN
    nl = 1
    la = INT(a)
    a = a - la - 1.0_dp
  END IF
  DO  n = 0, nl
    IF (a0 >= 2.0_dp) a = a + 1.0_dp
    IF (ABS(z) < 20.0_dp + ABS(b) .OR. a < 0.0_dp) THEN
      chg = (1.0_dp,0.0_dp)
      crg = (1.0_dp,0.0_dp)
      DO  j = 1, 500
        crg = crg * (a+j-1) / (j*(b+j-1)) * z
        chg = chg + crg
        IF (ABS((chg-chw)/chg) < 1.d-15) GO TO 40
        chw = chg
      END DO
    ELSE
      CALL gamma(a, g1)
      CALL gamma(b, g2)
      ba = b - a
      CALL gamma(ba, g3)
      cs1 = (1.0_dp,0.0_dp)
      cs2 = (1.0_dp,0.0_dp)
      cr1 = (1.0_dp,0.0_dp)
      cr2 = (1.0_dp,0.0_dp)
      DO  i = 1, 8
        cr1 = -cr1 * (a+i-1) * (a-b+i) / (z*i)
        cr2 = cr2 * (b-a+i-1) * (i-a) / (z*i)
        cs1 = cs1 + cr1
        cs2 = cs2 + cr2
      END DO
      x = REAL(z)
      y = AIMAG(z)
      IF (x == 0.0 .AND. y >= 0.0) THEN
        phi = 0.5_dp * pi
      ELSE IF (x == 0.0 .AND. y <= 0.0) THEN
        phi = -0.5_dp * pi
      ELSE
        phi = ATAN(y/x)
      END IF
      IF (phi > -0.5*pi .AND. phi < 1.5*pi) ns = 1
      IF (phi > -1.5*pi .AND. phi <= -0.5*pi) ns = -1
      cfac = EXP(ns*ci*pi*a)
      IF (y == 0.0_dp) cfac = COS(pi*a)
      chg1 = g2 / g3 * z ** (-a) * cfac * cs1
      chg2 = g2 / g1 * EXP(z) * z ** (a-b) * cs2
      chg = chg1 + chg2
    END IF
    40 IF (n == 0) cy0 = chg
    IF (n == 1) cy1 = chg
  END DO
  IF (a0 >= 2.0_dp) THEN
    DO  i = 1, la - 1
      chg = ((2.0_dp*a-b+z)*cy1 + (b-a)*cy0) / a
      cy0 = cy1
      cy1 = chg
      a = a + 1.0_dp
    END DO
  END IF
  IF (x0 < 0.0_dp) chg = chg * EXP(-z)
END IF
a = a1
z = z0
RETURN
END SUBROUTINE cchg


SUBROUTINE gamma(x, ga)

!    ==================================================
!    Purpose: Compute gamma function ג(x)
!    Input :  x  --- Argument of ג(x)
!                    ( x is not equal to 0,-1,-2,תתת)
!    Output:  GA --- ג(x)
!    ==================================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ga

REAL (dp), PARAMETER  :: g(26) = (/  1.0_dp, 0.5772156649015329_dp,  &
      -0.6558780715202538_dp, -0.420026350340952D-1, 0.1665386113822915_dp,  &
      -0.421977345555443D-1, -.96219715278770D-2,  &
      .72189432466630D-2, -.11651675918591D-2, -.2152416741149D-3,  &
      .1280502823882D-3, -.201348547807D-4, -.12504934821D-5,  &
      .11330272320D-5, -.2056338417D-6, .61160950D-8,  &
      .50020075D-8, -.11812746D-8, .1043427D-9, .77823D-11,  &
      -.36968D-11, .51D-12, -.206D-13, -.54D-14, .14D-14, .1D-15 /)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: gr, r, z
INTEGER    :: k, m, m1

IF (x == INT(x)) THEN
  IF (x > 0.0_dp) THEN
    ga = 1.0_dp
    m1 = x - 1
    DO  k = 2, m1
      ga = ga * k
    END DO
  ELSE
    ga = 1.0D+300
  END IF
ELSE
  IF (ABS(x) > 1.0_dp) THEN
    z = ABS(x)
    m = z
    r = 1.0_dp
    DO  k = 1, m
      r = r * (z-k)
    END DO
    z = z - m
  ELSE
    z = x
  END IF

  gr = g(26)
  DO  k = 25, 1, -1
    gr = gr * z + g(k)
  END DO
  ga = 1.0_dp / (gr*z)
  IF (ABS(x) > 1.0_dp) THEN
    ga = ga * r
    IF (x < 0.0_dp) ga = -pi / (x*ga*SIN(pi*x))
  END IF
END IF
RETURN
END SUBROUTINE gamma
 
END MODULE cchg_func
 
 
 
PROGRAM mcchg
USE cchg_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!    ===========================================================
!    Purpose: This program computes confluent hypergeometric
!             function M(a,b,z) with real parameters a, b, and
!             a complex argument z using subroutine CCHG
!    Input :  a --- Parameter
!             b --- Parameter
!             z --- Complex argument
!    Output:  CHG --- M(a,b,z)
!    Examples:
!       a      b        z        Re[M(a,b,z)]   Im[M(a,b,z)]
!      -------------------------------------------------------
!      3.3   4.25    10 + 0i    .61677489D+04    0
!      3.3   4.25    25 + 0i    .95781835D+10  -.15738228D-03
!      3.3   4.25     3 -  i    .75828716D+01  -.86815474D+01
!      3.3   4.25    15 +10i   -.58313765D+06  -.48195426D+05
!    ===========================================================

COMPLEX (dp)  :: chg, z
REAL (dp)     :: a, b, x, y

DO
  WRITE (*,*) 'Please enter a, b, x and y (z=x+iy): '
  READ (*,*) a, b, x, y
  WRITE (*,5100) a, b, x, y
  z = CMPLX(x, y, KIND=dp)
  CALL cchg(a, b, z, chg)
  WRITE (*,5000) chg
END DO
STOP

5000 FORMAT(t11, 'M(a,b,z) =', g18.8, ' + i ', g18.8)
5100 FORMAT(' a =', f5.1, ',  b =', f5.1, ',  x =', f5.1, ',  y =', f5.1)
END PROGRAM mcchg
