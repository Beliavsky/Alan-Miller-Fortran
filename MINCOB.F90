MODULE incob_func
 
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
 

SUBROUTINE incob(a, b, x, bix)

!    ========================================================
!    Purpose: Compute the incomplete beta function Ix(a,b)
!    Input :  a --- Parameter
!             b --- Parameter
!             x --- Argument ( 0 ף x ף 1 )
!    Output:  BIX --- Ix(a,b)
!    Routine called: BETA for computing beta function B(p,q)
!    ========================================================

REAL (dp), INTENT(IN)   :: a
REAL (dp), INTENT(IN)   :: b
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(OUT)  :: bix

REAL (dp)  :: bt, dk(51), fk(51), s0, t1, t2, ta, tb
INTEGER    :: k

s0 = (a+1.0_dp) / (a+b+2.0_dp)
CALL beta(a, b, bt)
IF (x <= s0) THEN
  DO  k = 1, 20
    dk(2*k) = k * (b-k) * x / (a+2.0_dp*k-1.0_dp) / (a+2.0_dp*k)
  END DO
  DO  k = 0, 20
    dk(2*k+1) = -(a+k) * (a+b+k) * x / (a+2._dp*k) / (a+2.0*k+1.0)
  END DO
  t1 = 0.0_dp
  DO  k = 20, 1, -1
    t1 = dk(k) / (1.0_dp+t1)
  END DO
  ta = 1.0_dp / (1.0_dp+t1)
  bix = x ** a * (1.0_dp-x) ** b / (a*bt) * ta
ELSE
  DO  k = 1, 20
    fk(2*k) = k * (a-k) * (1.0_dp-x) / (b+2.*k-1.0) / (b+2.0*k)
  END DO
  DO  k = 0, 20
    fk(2*k+1) = -(b+k) * (a+b+k) * (1._dp-x) / (b+2._dp*k) / (b+2._dp *k+1._dp)
  END DO
  t2 = 0.0_dp
  DO  k = 20, 1, -1
    t2 = fk(k) / (1.0_dp+t2)
  END DO
  tb = 1.0_dp / (1.0_dp+t2)
  bix = 1.0_dp - x ** a * (1.0_dp-x) ** b / (b*bt) * tb
END IF
RETURN
END SUBROUTINE incob



SUBROUTINE beta(p, q, bt)

!    ==========================================
!    Purpose: Compute the beta function B(p,q)
!    Input :  p --- Parameter  ( p > 0 )
!             q --- Parameter  ( q > 0 )
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

END MODULE incob_func
 
 
 
PROGRAM mincob
USE incob_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:40

!       =========================================================
!       Purpose: This program computes the incomplete beta
!                function Ix(a,b) using subroutine INCOB
!       Input :  a --- Parameter
!                b --- Parameter
!                x --- Argument ( 0 ף x ף 1 )
!       Output:  BIX --- Ix(a,b)
!       Example:
!                  a       b       x       Ix(a,b)
!                -----------------------------------
!                 1.0     3.0     .25     .57812500
!       =========================================================

REAL (dp)  :: a, b, bix, x

WRITE (*,*) 'Please enter a, b and x ( 0 ף x ף 1 ) '
READ (*,*) a, b, x
WRITE (*,*)
WRITE (*,*) '   a       b       x       Ix(a,b)'
WRITE (*,*) ' -----------------------------------'
CALL incob(a, b, x, bix)
WRITE (*,5000) a, b, x, bix
STOP

5000 FORMAT (' ', f5.1, '   ', f5.1, '   ', f5.2, f14.8)
END PROGRAM mincob
