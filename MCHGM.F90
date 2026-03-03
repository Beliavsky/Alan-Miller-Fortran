MODULE chgm_func
 
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


SUBROUTINE chgm(a, b, x, hg)

!       ===================================================
!       Purpose: Compute confluent hypergeometric function
!                M(a,b,x)
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,... )
!                x  --- Argument
!       Output:  HG --- M(a,b,x)
!       Routine called: GAMMA for computing ג(x)
!       ===================================================

REAL (dp), INTENT(IN OUT)  :: a
REAL (dp), INTENT(IN)      :: b
REAL (dp), INTENT(IN OUT)  :: x
REAL (dp), INTENT(OUT)     :: hg

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)  :: a0, a1, hg1, hg2, r, r1, r2, rg, sum1, sum2, ta, tb, tba,  &
              x0, xg, y0, y1
INTEGER    :: i, j, k, la, m, n, nl

a0 = a
a1 = a
x0 = x
hg = 0.0_dp
IF (b == 0.0_dp .OR. b == -ABS(INT(b))) THEN
  hg = 1.0D+300
ELSE IF (a == 0.0_dp .OR. x == 0.0_dp) THEN
  hg = 1.0_dp
ELSE IF (a == -1.0_dp) THEN
  hg = 1.0_dp - x / b
ELSE IF (a == b) THEN
  hg = EXP(x)
ELSE IF (a-b == 1.0_dp) THEN
  hg = (1.0_dp + x/b) * EXP(x)
ELSE IF (a == 1.0_dp .AND. b == 2.0_dp) THEN
  hg = (EXP(x)-1.0_dp) / x
ELSE IF (a == INT(a) .AND. a < 0.0_dp) THEN
  m = INT(-a)
  r = 1.0_dp
  hg = 1.0_dp
  DO  k = 1, m
    r = r * (a+k-1) / k / (b+k-1) * x
    hg = hg + r
  END DO
END IF
IF (hg /= 0.0_dp) RETURN
IF (x < 0.0_dp) THEN
  a = b - a
  a0 = a
  x = ABS(x)
END IF
IF (a < 2.0_dp) nl = 0
IF (a >= 2.0_dp) THEN
  nl = 1
  la = INT(a)
  a = a - la - 1
END IF
DO  n = 0, nl
  IF (a0 >= 2.0_dp) a = a + 1.0_dp
  IF (x <= 30.0_dp + ABS(b) .OR. a < 0.0_dp) THEN
    hg = 1.0_dp
    rg = 1.0_dp
    DO  j = 1, 500
      rg = rg * (a+j-1) / (j*(b+j-1)) * x
      hg = hg + rg
      IF (ABS(rg/hg) < 1.0D-15) GO TO 40
    END DO
  ELSE
    CALL gamma(a, ta)
    CALL gamma(b, tb)
    xg = b - a
    CALL gamma(xg, tba)
    sum1 = 1.0_dp
    sum2 = 1.0_dp
    r1 = 1.0_dp
    r2 = 1.0_dp
    DO  i = 1, 8
      r1 = -r1 * (a+i-1) * (a-b+i) / (x*i)
      r2 = -r2 * (b-a+i-1) * (a-i) / (x*i)
      sum1 = sum1 + r1
      sum2 = sum2 + r2
    END DO
    hg1 = tb / tba * x ** (-a) * COS(pi*a) * sum1
    hg2 = tb / ta * EXP(x) * x ** (a-b) * sum2
    hg = hg1 + hg2
  END IF
  40 IF (n == 0) y0 = hg
  IF (n == 1) y1 = hg
END DO
IF (a0 >= 2.0_dp) THEN
  DO  i = 1, la - 1
    hg = ((2.0_dp*a-b+x)*y1 + (b-a)*y0) / a
    y0 = y1
    y1 = hg
    a = a + 1.0_dp
  END DO
END IF
IF (x0 < 0.0_dp) hg = hg * EXP(x0)
a = a1
x = x0
RETURN
END SUBROUTINE chgm


SUBROUTINE gamma(x,ga)

!       ==================================================
!       Purpose: Compute gamma function ג(x)
!       Input :  x  --- Argument of ג(x)
!                       ( x is not equal to 0,-1,-2,תתת)
!       Output:  GA --- ג(x)
!       ==================================================


REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ga

REAL (dp), PARAMETER  :: g(26) = (/1.0_dp, 0.5772156649015329_dp,  &
      -0.6558780715202538_dp, -0.420026350340952D-1, 0.1665386113822915_dp,  &
      -0.421977345555443D-1, -.96219715278770D-2, .72189432466630D-2,  &
      -0.11651675918591D-2, -.2152416741149D-3, .1280502823882D-3,  &
      -0.201348547807D-4, -.12504934821D-5, 0.11330272320D-5,  &
      -0.2056338417D-6, .61160950D-8, .50020075D-8, -.11812746D-8,  &
       0.1043427D-9, .77823D-11, -.36968D-11, .51D-12, -.206D-13,   &
       -.54D-14, .14D-14, .1D-15 /)
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
    m = INT(z)
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
 
END MODULE chgm_func
 
 
 
PROGRAM mchgm
USE chgm_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!       ===========================================================
!       Purpose: This program computes the confluent hypergeometric
!                function M(a,b,x) using subroutine CHGM
!       Input  : a  --- Parameter
!                b  --- Parameter ( b <> 0,-1,-2,... )
!                x  --- Argument
!       Output:  HG --- M(a,b,x)
!       Example:
!                   a       b       x          M(a,b,x)
!                 -----------------------------------------
!                  1.5     2.0    20.0     .1208527185D+09
!                  4.5     2.0    20.0     .1103561117D+12
!                 -1.5     2.0    20.0     .1004836854D+05
!                 -4.5     2.0    20.0    -.3936045244D+03
!                  1.5     2.0    50.0     .8231906643D+21
!                  4.5     2.0    50.0     .9310512715D+25
!                 -1.5     2.0    50.0     .2998660728D+16
!                 -4.5     2.0    50.0    -.1806547113D+13
!       ===========================================================

REAL (dp)  :: a, b, x, hg

DO
  WRITE (*,*) 'Please enter a, b and x: '
  READ (*,*) a, b, x
  WRITE (*,*) '   a       b       x          M(a,b,x)'
  WRITE (*,*) ' -----------------------------------------'
  CALL chgm(a, b, x, hg)
  WRITE (*,5000) a, b, x, hg
END DO
STOP

5000 FORMAT (' ', f5.1, '   ', f5.1, '   ', f5.1, g20.10)
END PROGRAM mchgm
