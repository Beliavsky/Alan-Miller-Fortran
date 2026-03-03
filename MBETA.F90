MODULE beta_func
 
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
 

SUBROUTINE beta(p, q, bt)

!    ==========================================
!    Purpose: Compute the beta function B(p,q)
!    Input :  p  --- Parameter  ( p > 0 )
!             q  --- Parameter  ( q > 0 )
!    Output:  BT --- B(p,q)
!    Routine called: GAMMA for computing ג(x)
!    ==========================================

REAL (dp), INTENT(IN)   :: p
REAL (dp), INTENT(IN)   :: q
REAL (dp), INTENT(OUT)  :: bt

REAL (dp)  :: gp, gpq, gq, ppq

CALL gamma(p, gp)
CALL gamma(q, gq)
ppq = p + q
CALL gamma(ppq, gpq)
bt = gp * gq / gpq
RETURN
END SUBROUTINE beta



SUBROUTINE gamma(x, ga)

!    ==================================================
!    Purpose: Compute gamma function ג(x)
!    Input :  x  --- Argument of ג(x)
!                    ( x is not equal to 0,-1,-2,תתת)
!    Output:  GA --- ג(x)
!    ==================================================

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: ga

INTEGER    :: k, m, m1
REAL (dp)  :: gr, r, z
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp), PARAMETER  :: g(26) = (/ 1.0_dp, 0.5772156649015329_dp,  &
      -0.6558780715202538_dp, -0.420026350340952D-1, 0.1665386113822915_dp,   &
      -0.421977345555443D-1, -0.96219715278770D-2, 0.72189432466630D-2,  &
      -0.11651675918591D-2, -0.2152416741149D-3, 0.1280502823882D-3,  &
      -0.201348547807D-4, -0.12504934821D-5, 0.11330272320D-5,  &
      -0.2056338417D-6, 0.61160950D-8, 0.50020075D-8, -0.11812746D-8,  &
       0.1043427D-9, 0.77823D-11, -0.36968D-11, 0.51D-12, -0.206D-13,  &
      -0.54D-14, 0.14D-14, 0.1D-15 /)

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
 
END MODULE beta_func
 
 
 
PROGRAM mbeta
USE beta_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:34

!    ====================================================
!    Purpose: This program computes the beta function
!             B(p,q) for p > 0 and q > 0 using
!             subroutine BETA
!    Input :  p  --- Parameter  ( p > 0 )
!             q  --- Parameter  ( q > 0 )
!    Output:  BT --- B(p,q)
!    Examples:
!              p       q           B(p,q)
!            ---------------------------------
!             1.5     2.0     .2666666667D+00
!             2.5     2.0     .1142857143D+00
!             1.5     3.0     .1523809524D+00
!    ====================================================

REAL (dp)  :: p, q, bt

DO
  WRITE (*,*) 'Please enter p and q: '
  READ (*,*) p, q
  WRITE (*,*)
  WRITE (*,*) '    p       q           B(p,q)'
  WRITE (*,*) '  ---------------------------------'
  CALL beta(p, q, bt)
  WRITE (*,5000) p, q, bt
END DO
STOP

5000 FORMAT (t3, f5.1, '   ', f5.1, g20.10)
END PROGRAM mbeta
