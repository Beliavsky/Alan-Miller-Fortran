MODULE fcszo_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 13 January 2002
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! In the first DO loop in routine FCSZO, the square root of an integer was
! taken instead of the square root of the double precision value of the integer.
! The value of W was undefined the first time through the iteratiom in
! routine FCSZO.
! The value of WB0 was undefined the first time through the DO k = 1, 80 loop
! in routine CFS.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE fcszo(kf, nt, zo)

!     ===========================================================
!     Purpose: Compute the complex zeros of Fresnel integral C(z)
!              or S(z) using modified Newton's iteration method
!     Input :  KF  --- Function code
!                      KF=1 for C(z) or KF=2 for S(z)
!              NT  --- Total number of zeros
!     Output:  ZO(L) --- L-th zero of C(z) or S(z)
!     Routines called:
!          (1) CFC for computing Fresnel integral C(z)
!          (2) CFS for computing Fresnel integral S(z)
!     ===========================================================

INTEGER, INTENT(IN)        :: kf
INTEGER, INTENT(IN)        :: nt
COMPLEX (dp), INTENT(OUT)  :: zo(nt)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
REAL (dp)     :: psq, px, py, w, w0
COMPLEX (dp)  :: z, zd, zf, zfd, zgd, zp, zq, zw
INTEGER       :: i, it, j, nr

DO  nr = 1, nt
  IF (kf == 1) psq = SQRT(4.0_dp*nr - 1.0_dp)
  IF (kf == 2) psq = SQRT(4.0_dp*nr)
  px = psq - LOG(pi*psq) / (pi*pi*psq**3.0)
  py = LOG(pi*psq) / (pi*psq)
  z = CMPLX(px, py, KIND=dp)
  w = ABS(z)
  IF (kf == 2) THEN
    IF (nr == 2) z = (2.8334, 0.2443)
    IF (nr == 3) z = (3.4674, 0.2185)
    IF (nr == 4) z = (4.0025, 0.2008)
  END IF
  it = 0

  10 it = it + 1
  IF (kf == 1) CALL cfc(z, zf, zd)
  IF (kf == 2) CALL cfs(z, zf, zd)
  zp = (1.0_dp,0.0_dp)
  DO  i = 1, nr - 1
    zp = zp * (z-zo(i))
  END DO
  zfd = zf / zp
  zq = (0.0_dp,0.0_dp)
  DO  i = 1, nr - 1
    zw = (1.0_dp,0.0_dp)
    DO  j = 1, nr - 1
      IF (j /= i) THEN
        zw = zw * (z-zo(j))
      END IF
    END DO
    zq = zq + zw
  END DO
  zgd = (zd - zq*zfd) / zp
  z = z - zfd / zgd
  w0 = w
  w = ABS(z)
  IF (it <= 50 .AND. ABS((w-w0)/w) > 1.0D-12) GO TO 10
  zo(nr) = z
END DO
RETURN
END SUBROUTINE fcszo



SUBROUTINE cfc(z, zf, zd)

!     =========================================================
!     Purpose: Compute complex Fresnel integral C(z) and C'(z)
!     Input :  z --- Argument of C(z)
!     Output:  ZF --- C(z)
!              ZD --- C'(z)
!     =========================================================

COMPLEX (dp), INTENT(IN)    :: z
COMPLEX (dp), INTENT(OUT)   :: zf
COMPLEX (dp), INTENT(OUT)   :: zd

REAL (dp), PARAMETER  :: eps = 1.0D-14, pi = 3.141592653589793_dp
COMPLEX (dp)  :: c, cf, cf0, cf1, cg, cr, z0, zp, zp2
REAL (dp)     :: w0, wa, wa0
INTEGER       :: k, m

w0 = ABS(z)
zp = 0.5_dp * pi * z * z
zp2 = zp * zp
z0 = (0.0_dp,0.0_dp)
IF (z == z0) THEN
  c = z0
ELSE IF (w0 <= 2.5) THEN
  cr = z
  c = cr
  DO  k = 1, 80
    cr = -.5_dp * cr * (4.0_dp*k-3.0_dp) / k / (2.0_dp*k-1.0_dp) / (4.0_dp*k+1.0_dp) * zp2
    c = c + cr
    wa = ABS(c)
    IF (ABS((wa-wa0)/wa) < eps .AND. k > 10) GO TO 50
    wa0 = wa
  END DO
ELSE IF (w0 > 2.5 .AND. w0 < 4.5) THEN
  m = 85
  c = z0
  cf1 = z0
  cf0 = (1.0D-100,0.0_dp)
  DO  k = m, 0, -1
    cf = (2*k+3) * cf0 / zp - cf1
    IF (k == INT(k/2)*2) c = c + cf
    cf1 = cf0
    cf0 = cf
  END DO
  c = SQRT(2.0_dp/(pi*zp)) * SIN(zp) / cf * c
ELSE
  cr = (1.0_dp,0.0_dp)
  cf = (1.0_dp,0.0_dp)
  DO  k = 1, 20
    cr = -.25_dp * cr * (4*k-1) * (4*k-3) / zp2
    cf = cf + cr
  END DO
  cr = 1.0_dp / (pi*z*z)
  cg = cr
  DO  k = 1, 12
    cr = -.25_dp * cr * (4*k+1) * (4*k-1) / zp2
    cg = cg + cr
  END DO
  c = .5_dp + (cf*SIN(zp) - cg*COS(zp)) / (pi*z)
END IF
50 zf = c
zd = COS(0.5*pi*z*z)
RETURN
END SUBROUTINE cfc



SUBROUTINE cfs(z, zf, zd)

!     =========================================================
!     Purpose: Compute complex Fresnel Integral S(z) and S'(z)
!     Input :  z  --- Argument of S(z)
!     Output:  ZF --- S(z)
!              ZD --- S'(z)
!     =========================================================

COMPLEX (dp), INTENT(IN)    :: z
COMPLEX (dp), INTENT(OUT)   :: zf
COMPLEX (dp), INTENT(OUT)   :: zd

REAL (dp), PARAMETER  :: eps = 1.0D-14, pi = 3.141592653589793_dp
COMPLEX (dp)  :: cf, cf0, cf1, cg, cr, s, z0, zp, zp2
REAL (dp)     :: w0, wb, wb0
INTEGER       :: k, m

w0 = ABS(z)
zp = 0.5_dp * pi * z * z
zp2 = zp * zp
z0 = (0.0_dp,0.0_dp)
IF (z == z0) THEN
  s = z0
ELSE IF (w0 <= 2.5) THEN
  s = z * zp / 3.0_dp
  cr = s
  wb0 = ABS(s)
  DO  k = 1, 80
    cr = -.5_dp * cr * (4.0_dp*k-1.0_dp) / k / (2.0_dp*k+1.0_dp) / (4.0_dp*k+3.0_dp) * zp2
    s = s + cr
    wb = ABS(s)
    IF (ABS(wb-wb0) < eps .AND. k > 10) GO TO 50
    wb0 = wb
  END DO
ELSE IF (w0 > 2.5 .AND. w0 < 4.5) THEN
  m = 85
  s = z0
  cf1 = z0
  cf0 = (1.0D-100,0.0_dp)
  DO  k = m, 0, -1
    cf = (2*k+3) * cf0 / zp - cf1
    IF (k /= INT(k/2)*2) s = s + cf
    cf1 = cf0
    cf0 = cf
  END DO
  s = SQRT(2.0_dp/(pi*zp)) * SIN(zp) / cf * s
ELSE
  cr = (1.0_dp,0.0_dp)
  cf = (1.0_dp,0.0_dp)
  DO  k = 1, 20
    cr = -.25_dp * cr * (4*k-1) * (4*k-3) / zp2
    cf = cf + cr
  END DO
  cr = 1.0_dp / (pi*z*z)
  cg = cr
  DO  k = 1, 12
    cr = -.25_dp * cr * (4*k+1) * (4*k-1) / zp2
    cg = cg + cr
  END DO
  s = .5_dp - (cf*COS(zp) + cg*SIN(zp)) / (pi*z)
END IF
50 zf = s
zd = SIN(0.5*pi*z*z)
RETURN
END SUBROUTINE cfs
 
END MODULE fcszo_func
 
 
 
PROGRAM mfcszo
USE fcszo_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:39

!     ===========================================================
!     Purpose : This program computes the complex zeros of the
!               Fresnel integral C(z) or S(z) using subroutine FCSZO
!     Input :   KF  --- Function code
!                       KF=1 for C(z) or KF=2 for S(z)
!               NT  --- Total number of zeros
!     Output:   ZO(L) --- L-th zero of C(z) or S(z)
!     Example:  NT=10

!     n     Complex zeros of C(z)        Complex zeros of S(z)
!   ------------------------------------------------------------
!     1    1.7436675 + i .30573506      2.0092570 + i .28854790
!     2    2.6514596 + i .25290396      2.8334772 + i .24428524
!     3    3.3203593 + i .22395346      3.4675331 + i .21849268
!     4    3.8757345 + i .20474747      4.0025782 + i .20085103
!     5    4.3610635 + i .19066973      4.4741893 + i .18768859
!     6    4.7976077 + i .17970801      4.9006784 + i .17732036
!     7    5.1976532 + i .17081930      5.2929467 + i .16884418
!     8    5.5690602 + i .16339854      5.6581068 + i .16172492
!     9    5.9172173 + i .15706585      6.0011034 + i .15562108
!    10    6.2460098 + i .15156826      6.3255396 + i .15030246
!     ===========================================================

COMPLEX (dp)  :: zo(100)
INTEGER       :: i, kf, nt

WRITE (*,*) 'Please Enter KF and NT '
READ (*,*) kf, nt
WRITE (*,5000) kf, nt
WRITE (*,*) ' *****     Please Wait !     *****'
CALL fcszo(kf, nt, zo)
WRITE (*,*)
IF (kf == 1) WRITE (*,*) '  n        Complex zeros of C(z)'
IF (kf == 2) WRITE (*,*) '  n        Complex zeros of S(z)'
WRITE (*,*) '-----------------------------------'
DO  i = 1, nt
  WRITE (*,5100) i, zo(i)
END DO
STOP

5000 FORMAT ('  KF=', i2, ',     NT=', i3)
5100 FORMAT (' ', i3, f13.8, '  +i', f13.8)
END PROGRAM mfcszo
