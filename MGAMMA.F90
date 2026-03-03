MODULE gamma_func
 
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
 

SUBROUTINE gamma(x, ga)

!    ==================================================
!    Purpose: Compute the gamma function ג(x)
!    Input :  x  --- Argument of ג(x)
!                    ( x is not equal to 0,-1,-2,תתת )
!    Output:  GA --- ג(x)
!    ==================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ga

REAL (dp), PARAMETER  :: g(26) = (/  &
       1.0_dp, 0.5772156649015329_dp, -0.6558780715202538_dp,  &
      -0.420026350340952D-1, 0.1665386113822915_dp,   &
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
 
END MODULE gamma_func
 
 
 
PROGRAM mgamma
USE gamma_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:39

!    ====================================================
!    Purpose: This program computes the gamma function
!             ג(x) using subroutine GAMMA
!    Examples:
!                x            ג(x)
!             ----------------------------
!               1/3       2.678938534708
!               0.5       1.772453850906
!              -0.5      -3.544907701811
!              -1.5       2.363271801207
!               5.0      24.000000000000
!    ====================================================

REAL (dp)  :: x, ga
INTEGER    :: k

REAL (dp)  :: a(5) = (/ 0.333333333333333333_dp, 0.5_dp, -0.5_dp, -1.5_dp, &
                        5.0_dp /)
WRITE (*,*) '     x            ג(x)'
WRITE (*,*) ' ----------------------------'
DO  k = 1, 5
  x = a(k)
  CALL gamma(x, ga)
  WRITE (*,5000) x, ga
END DO
WRITE (*,*) 'Please enter x:'
READ (*,*) x
CALL gamma(x, ga)
WRITE (*,5000) x, ga
STOP

5000 FORMAT (' ', f8.4, g20.12)
END PROGRAM mgamma
