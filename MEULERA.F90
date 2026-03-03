MODULE eulera_func
 
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


SUBROUTINE eulera(n, en)

!       ==================================
!       Purpose: Compute Euler number En
!       Input :  n --- Serial number
!       Output:  EN(n) --- En
!       ==================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: en(0:n)

REAL (dp)  :: r, s
INTEGER    :: j, k, m

en(0) = 1.0_dp
DO  m = 1, n / 2
  s = 1.0_dp
  DO  k = 1, m - 1
    r = 1.0_dp
    DO  j = 1, 2 * k
      r = r * (2*(m-k)+j) / REAL(j, KIND=dp)
    END DO
    s = s + r * en(2*k)
  END DO
  en(2*m) = -s
END DO
RETURN
END SUBROUTINE eulera
 
END MODULE eulera_func
 
 
 
PROGRAM meulera
USE eulera_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!    =====================================================
!    Purpose: This program computes Euler number En using
!             subroutine EULERA
!    Example: Compute Euler number En for n = 0,2,...,10
!             Computed results:

!                n            En
!              --------------------------
!                0     .100000000000D+01
!                2    -.100000000000D+01
!                4     .500000000000D+01
!                6    -.610000000000D+02
!                8     .138500000000D+04
!               10    -.505210000000D+05
!    =====================================================

REAL (dp)  :: e(0:200)
INTEGER    :: k, n

WRITE (*,*) '  Please enter Nmax '
READ (*,*) n
CALL eulera(n, e)
WRITE (*,*) '   n            En'
WRITE (*,*) ' --------------------------'
DO  k = 0, n, 2
  WRITE (*,5000) k, e(k)
END DO
STOP

5000 FORMAT (t3, i3, g22.12)
END PROGRAM meulera
