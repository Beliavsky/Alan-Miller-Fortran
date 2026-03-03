MODULE enxb_func
 
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


SUBROUTINE enxb(n, x, en)

!       ===============================================
!       Purpose: Compute exponential integral En(x)
!       Input :  x --- Argument of En(x)
!                n --- Order of En(x)  (n = 0,1,2,...)
!       Output:  EN(n) --- En(x)
!       ===============================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: en(0:n)

REAL (dp)  :: ens, ps, r, rp, s, s0, t, t0
INTEGER    :: j, k, l, m

IF (x == 0.0) THEN
  en(0) = 1.0D+300
  en(1) = 1.0D+300
  DO  k = 2, n
    en(k) = 1.0_dp / (k-1)
  END DO
  RETURN
ELSE IF (x <= 1.0) THEN
  en(0) = EXP(-x) / x
  DO  l = 1, n
    rp = 1.0_dp
    DO  j = 1, l - 1
      rp = -rp * x / j
    END DO
    ps = -0.5772156649015328_dp
    DO  m = 1, l - 1
      ps = ps + 1.0_dp / m
    END DO
    ens = rp * (-LOG(x) + ps)
    s = 0.0_dp
    DO  m = 0, 20
      IF (m /= l-1) THEN
        r = 1.0_dp
        DO  j = 1, m
          r = -r * x / j
        END DO
        s = s + r / (m-l+1)
        IF (ABS(s-s0) < ABS(s)*1.0D-15) EXIT
        s0 = s
      END IF
    END DO
    en(l) = ens - s
  END DO
ELSE
  en(0) = EXP(-x) / x
  m = 15 + INT(100.0/x)
  DO  l = 1, n
    t0 = 0.0_dp
    DO  k = m, 1, -1
      t0 = (l+k-1) / (1.0_dp + k/(x+t0))
    END DO
    t = 1.0_dp / (x+t0)
    en(l) = EXP(-x) * t
  END DO
END IF
RETURN
END SUBROUTINE enxb
 
END MODULE enxb_func
 
 
 
PROGRAM menxb
USE enxb_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!       =======================================================
!       Purpose: This program computes the exponential integral
!                En(x) using subroutine ENXB
!       Example: x = 10.0
!                   n         En(x)
!                 ----------------------
!                   0     .45399930D-05
!                   1     .41569689D-05
!                   2     .38302405D-05
!                   3     .35487626D-05
!                   4     .33041014D-05
!                   5     .30897289D-05
!       =======================================================

REAL (dp)  :: en(0:100), x
INTEGER    :: k, n

WRITE (*,*) 'Please enter n and x '
READ (*,*) n, x
WRITE (*,5000) n, x
WRITE (*,*)
WRITE (*,*) '   n         En(x)'
WRITE (*,*) ' ----------------------'
CALL enxb(n, x, en)
DO  k = 0, n
  WRITE (*,5100) k, en(k)
END DO
STOP

5000 FORMAT (t6, i3, ',   x=', f5.1)
5100 FORMAT (t3, i3, g18.8)
END PROGRAM menxb
