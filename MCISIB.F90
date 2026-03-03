MODULE cisib_func
 
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
 

SUBROUTINE cisib(x, ci, si)

!       =============================================
!       Purpose: Compute cosine and sine integrals
!                Si(x) and Ci(x) ( x ò 0 )
!       Input :  x  --- Argument of Ci(x) and Si(x)
!       Output:  CI --- Ci(x)
!                SI --- Si(x)
!       =============================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ci
REAL (dp), INTENT(OUT)     :: si

REAL (dp)  :: fx, gx, x2

x2 = x * x
IF (x == 0.0) THEN
  ci = -1.0D+300
  si = 0.0_dp
ELSE IF (x <= 1.0_dp) THEN
  ci = ((((-3.0D-8*x2 + 3.10D-6)*x2 - 2.3148D-4)*x2 + 1.041667D-2)*x2 -  &
        0.25) * x2 + 0.577215665_dp + LOG(x)
  si = ((((3.1D-7*x2 - 2.834D-5)*x2 + 1.66667D-003)*x2 - 5.555556D-002)* x2 + 1.0) * x
ELSE
  fx = ((((x2 + 38.027264_dp)*x2 + 265.187033_dp)*x2 + 335.67732_dp)*x2 +   &
       38.102495_dp) / ((((x2 + 40.021433_dp)*x2 + 322.624911_dp)*x2 +   &
       570.23628_dp)*x2 + 157.105423_dp)
  gx = ((((x2 + 42.242855_dp)*x2 + 302.757865_dp)*x2 + 352.018498_dp)*x2 +   &
       21.821899_dp) / ((((x2 + 48.196927_dp)*x2 + 482.485984_dp)*x2 +   &
       1114.978885_dp)*x2 + 449.690326_dp) / x
  ci = fx * SIN(x) / x - gx * COS(x) / x
  si = 1.570796327_dp - fx * COS(x) / x - gx * SIN(x) / x
END IF
RETURN
END SUBROUTINE cisib
 
END MODULE cisib_func
 
 
 
PROGRAM mcisib
USE cisib_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:36

!       ========================================================
!       Purpose: This program computes the cosine and sine
!                integrals using subroutine CISIB
!       Input :  x  --- Argument of Ci(x) and Si(x)
!       Output:  CI --- Ci(x)
!                SI --- Si(x)
!       Example:

!                   x        Ci(x)           Si(x)
!                ------------------------------------
!                  0.0    - ì                 0
!                  5.0    -.190030D+00      1.549931
!                 10.0    -.454563D-01      1.658348
!                 20.0     .444201D-01      1.548241
!                 30.0    -.330326D-01      1.566757
!                 40.0     .190201D-01      1.586985
!       ========================================================

REAL (dp)  :: ci, si, x

DO
  WRITE (*,*) 'Please enter x '
  READ (*,*) x
  WRITE (*,*) '   x        Ci(x)           Si(x)'
  WRITE (*,*) '------------------------------------'
  CALL cisib(x, ci, si)
  IF (x /= 0.0_dp) WRITE (*,5000) x, ci, si
  IF (x == 0.0_dp) WRITE (*,5100)
END DO
STOP

5000 FORMAT (' ', f5.1, g16.6, f14.6)
5100 FORMAT ('    .0    - ì', t32, '0')
END PROGRAM mcisib
