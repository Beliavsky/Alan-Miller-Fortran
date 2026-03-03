MODULE incog_func
 
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


SUBROUTINE incog(a, x, gin, gim, gip)

!    ================================================
!    Purpose: Compute the incomplete gamma function
!             r(a,x), ג(a,x) and P(a,x)
!    Input :  a   --- Parameter ( a ף 170 )
!             x   --- Argument
!    Output:  GIN --- r(a,x)
!             GIM --- ג(a,x)
!             GIP --- P(a,x)
!    Routine called: GAMMA for computing ג(x)
!    ================================================

REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: gin
REAL (dp), INTENT(OUT)  :: gim
REAL (dp), INTENT(OUT)  :: gip

REAL (dp)  :: ga, r, s, t0, xam
INTEGER    :: k

xam = -x + a * LOG(x)
IF (xam > 700.0 .OR. a > 170.0) THEN
  WRITE (*,*) 'a and/or x too large'
  STOP
END IF
IF (x == 0.0) THEN
  gin = 0.0
  CALL gamma(a, ga)
  gim = ga
  gip = 0.0
ELSE IF (x <= 1.0+a) THEN
  s = 1.0D0 / a
  r = s
  DO  k = 1, 60
    r = r * x / (a+k)
    s = s + r
    IF (ABS(r/s) < 1.0D-15) EXIT
  END DO
  gin = EXP(xam) * s
  CALL gamma(a, ga)
  gip = gin / ga
  gim = ga - gin
ELSE IF (x > 1.0+a) THEN
  t0 = 0.0D0
  DO  k = 60, 1, -1
    t0 = (k-a) / (1.0D0 + k/(x+t0))
  END DO
  gim = EXP(xam) / (x+t0)
  CALL gamma(a,ga)
  gin = ga - gim
  gip = 1.0D0 - gim / ga
END IF
RETURN
END SUBROUTINE incog



SUBROUTINE gamma(x, ga)

!    ==================================================
!    Purpose: Compute gamma function ג(x)
!    Input :  x  --- Argument of ג(x)
!                    ( x is not equal to 0,-1,-2,תתת)
!    Output:  GA --- ג(x)
!    ==================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ga

REAL (dp), PARAMETER  :: g(26) = (/ 1.0_dp, 0.5772156649015329_dp,  &
     -0.6558780715202538_dp, -0.420026350340952D-1, 0.1665386113822915_dp,   &
     -0.421977345555443D-1, -.96219715278770D-2, .72189432466630D-2,  &
     -0.11651675918591D-2, -.2152416741149D-3, .1280502823882D-3,  &
     -0.201348547807D-4, -.12504934821D-5, .11330272320D-5, -.2056338417D-6, &
      0.61160950D-8, .50020075D-8, -.11812746D-8, .1043427D-9, .77823D-11,  &
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

END MODULE incog_func
 
 
 
PROGRAM mincog
USE incog_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:40

!    ============================================================
!    Purpose: This program computes the incomplete gamma function
!             r(a,x), ג(a,x) and P(a,x) using subroutine INCOG
!    Input :  a   --- Parameter
!             x   --- Argument
!    Output:  GIN --- r(a,x)
!             GIM --- ג(a,x)
!             GIP --- P(a,x)
!    Example:
!         a     x      r(a,x)         ג(a,x)         P(a,x)
!        -------------------------------------------------------
!        3.0   5.0  .17506960D+01  .24930404D+00  .87534798D+00
!    =============================================================

REAL (dp) :: a, x, gin, gim, gip

WRITE (*,*) 'Plese enter a and x '
READ (*,*) a, x
WRITE (*,*)
WRITE (*,*) '   a     x      r(a,x)         ג(a,x)         P(a,x)'
WRITE (*,*) ' --------------------------------------------------------'
CALL incog(a, x, gin, gim, gip)
WRITE (*,5000) a, x, gin, gim, gip
STOP

5000 FORMAT (' ', f5.1, ' ', f5.1, 3g15.8)
END PROGRAM mincog
