MODULE enxa_func
 
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
 

SUBROUTINE enxa(n, x, en)

!     ============================================
!     Purpose: Compute exponential integral En(x)
!     Input :  x --- Argument of En(x) ( x ó 20 )
!              n --- Order of En(x)
!     Output:  EN(n) --- En(x)
!     Routine called: E1XB for computing E1(x)
!     ============================================

INTEGER, INTENT(IN)      :: n
REAL (dp), INTENT(IN)    :: x
REAL (dp), INTENT(OUT)   :: en(0:n)

REAL (dp)  :: e1, ek
INTEGER    :: k

en(0) = EXP(-x) / x
CALL e1xb(x, e1)
en(1) = e1
DO  k = 2, n
  ek = (EXP(-x) - x*e1) / (k-1)
  en(k) = ek
  e1 = ek
END DO
RETURN
END SUBROUTINE enxa



SUBROUTINE e1xb(x, e1)

!     ============================================
!     Purpose: Compute exponential integral E1(x)
!     Input :  x  --- Argument of E1(x)
!     Output:  E1 --- E1(x)
!     ============================================

REAL (dp), INTENT(IN)    :: x
REAL (dp), INTENT(OUT)   :: e1

REAL (dp)  :: ga, r, t0, t
INTEGER    :: k, m

IF (x == 0.0) THEN
  e1 = 1.0D+300
ELSE IF (x <= 1.0) THEN
  e1 = 1.0_dp
  r = 1.0_dp
  DO  k = 1, 25
    r = -r * k * x / (k+1) ** 2
    e1 = e1 + r
    IF (ABS(r) <= ABS(e1)*1.0D-15) EXIT
  END DO
  ga = 0.5772156649015328_dp
  e1 = -ga - LOG(x) + x * e1
ELSE
  m = 20 + INT(80.0/x)
  t0 = 0.0_dp
  DO  k = m, 1, -1
    t0 = k / (1.0_dp+k/(x+t0))
  END DO
  t = 1.0_dp / (x+t0)
  e1 = EXP(-x) * t
END IF
RETURN
END SUBROUTINE e1xb
 
END MODULE enxa_func
 
 
 
PROGRAM menxa
USE enxa_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!     =======================================================
!     Purpose: This program computes the exponential integral
!              En(x) using subroutine ENXA
!     Example: x = 10.0
!                 n         En(x)
!               ----------------------
!                 0     .45399930D-05
!                 1     .41569689D-05
!                 2     .38302405D-05
!                 3     .35487626D-05
!                 4     .33041014D-05
!                 5     .30897289D-05
!     =======================================================

REAL (dp)  :: en(0:100), x
INTEGER    :: k, n

WRITE (*,*) 'Please enter n and x '
READ (*,*) n, x
WRITE (*,5000) n, x
WRITE (*,*)
WRITE (*,*) '   n         En(x)'
WRITE (*,*) ' ----------------------'
CALL enxa(n, x, en)
DO  k = 0, n
  WRITE (*,5100) k, en(k)
END DO
STOP

5000 FORMAT (t6, i3, ',   x=', f5.1)
5100 FORMAT ('  ', i3, g18.8)
END PROGRAM menxa
