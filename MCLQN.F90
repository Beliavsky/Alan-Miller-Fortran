MODULE clqn_func
 
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
 

SUBROUTINE clqn(n, x, y, cqn, cqd)

!       ==================================================
!       Purpose: Compute the Legendre functions Qn(z) and
!                their derivatives Qn'(z) for a complex
!                argument
!       Input :  x --- Real part of z
!                y --- Imaginary part of z
!                n --- Degree of Qn(z), n = 0,1,2,...
!       Output:  CQN(n) --- Qn(z)
!                CQD(n) --- Qn'(z)
!       ==================================================

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x
REAL (dp), INTENT(IN OUT)  :: y
COMPLEX (dp), INTENT(OUT)  :: cqn(0:n)
COMPLEX (dp), INTENT(OUT)  :: cqd(0:n)

COMPLEX (dp)  :: cq0, cq1, cqf0, cqf1, cqf2, z
INTEGER       :: k, km, ls

z = CMPLX(x, y, KIND=dp)
IF (z == 1.0_dp) THEN
  DO  k = 0, n
    cqn(k) = (1.0D+300,0.0_dp)
    cqd(k) = (1.0D+300,0.0_dp)
  END DO
  RETURN
END IF
ls = 1
IF (ABS(z) > 1.0_dp) ls = -1
cq0 = 0.5_dp * LOG(ls*(1.0_dp + z)/(1.0_dp - z))
cq1 = z * cq0 - 1.0_dp
cqn(0) = cq0
cqn(1) = cq1
IF (ABS(z) < 1.0001_dp) THEN
  cqf0 = cq0
  cqf1 = cq1
  DO  k = 2, n
    cqf2 = ((2*k-1)*z*cqf1 - (k-1)*cqf0) / k
    cqn(k) = cqf2
    cqf0 = cqf1
    cqf1 = cqf2
  END DO
ELSE
  IF (ABS(z) > 1.1_dp) THEN
    km = 40 + n
  ELSE
    km = (40+n) * INT(-1.0-1.8*LOG(ABS(z-1.0)))
  END IF
  cqf2 = 0.0_dp
  cqf1 = 1.0_dp
  DO  k = km, 0, -1
    cqf0 = ((2*k+3)*z*cqf1 - (k+2)*cqf2) / (k+1)
    IF (k <= n) cqn(k) = cqf0
    cqf2 = cqf1
    cqf1 = cqf0
  END DO
  DO  k = 0, n
    cqn(k) = cqn(k) * cq0 / cqf0
  END DO
END IF
cqd(0) = (cqn(1) - z*cqn(0)) / (z*z - 1.0_dp)
DO  k = 1, n
  cqd(k) = (k*z*cqn(k) - k*cqn(k-1)) / (z*z - 1.0_dp)
END DO
RETURN
END SUBROUTINE clqn
 
END MODULE clqn_func
 
 
 
PROGRAM mclqn
USE clqn_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:37

!       ==========================================================
!       Purpose: This program computes the Legendre polynomials
!                Qn(z) and Qn'(z) for a complex argument using
!                subroutine CLQN
!       Input :  x --- Real part of z
!                y --- Imaginary part of z
!                n --- Degree of Qn(z), n = 0,1,...
!       Output:  CQN(n) --- Qn(z)
!                CQD(n) --- Qn'(z)
!       Examples:

!       z = 0.5 + 0.5 i
!       n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]
!      -----------------------------------------------------------
!       0   .402359D+00   .553574D+00   .800000D+00   .400000D+00
!       1  -.107561D+01   .477967D+00   .602359D+00   .115357D+01
!       2  -.136636D+01  -.725018D+00  -.242682D+01   .183390D+01
!       3   .182619D+00  -.206146D+01  -.622944D+01  -.247151D+01
!       4   .298834D+01  -.110022D+01  -.114849D+01  -.125963D+02
!       5   .353361D+01   .334847D+01   .206656D+02  -.123735D+02

!       z = 3.0 + 2.0 i
!       n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn'(z)]    Im[Qn'(z)]
!      -----------------------------------------------------------
!       0   .229073D+00  -.160875D+00  -.250000D-01   .750000D-01
!       1   .896860D-02  -.244805D-01   .407268D-02   .141247D-01
!       2  -.736230D-03  -.281865D-02   .190581D-02   .155860D-02
!       3  -.264727D-03  -.227023D-03   .391535D-03   .314880D-04
!       4  -.430648D-04  -.443187D-05   .527190D-04  -.305592D-04
!       5  -.481362D-05   .265297D-05   .395108D-05  -.839883D-05
!       ==========================================================

COMPLEX (dp)  :: cqn(0:100), cqd(0:100)
REAL (dp)     :: x, y
INTEGER       :: k, n

WRITE (*,*) '  Please enter Nmax, x and y (z=x+iy): '
READ (*,*) n, x, y
WRITE (*,5100) x, y
WRITE (*,*)
CALL clqn(n, x, y, cqn, cqd)
WRITE (*,*) '  n    Re[Qn(z)]     Im[Qn(z)]     Re[Qn''(Z)]    Im[Qn''(Z)]'
WRITE (*,*) ' -----------------------------------------------------------'
DO  k = 0, n
  WRITE (*,5000) k, cqn(k), cqd(k)
END DO
STOP

5000 FORMAT (' ', i3, 4g14.6)
5100 FORMAT ('   x =', f5.1, ',  y =', f5.1)
END PROGRAM mclqn
