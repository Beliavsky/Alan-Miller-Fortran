MODULE cyzo_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 3 January 2002
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Value of W was undefined in routine CYZO.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS
 

SUBROUTINE cyzo(nt, kf, kc, zo, zv)

!       ===========================================================
!       Purpose : Compute the complex zeros of Y0(z), Y1(z) and
!                 Y1'(z), and their associated values at the zeros
!                 using the modified Newton's iteration method
!       Input:    NT --- Total number of zeros/roots
!                 KF --- Function choice code
!                        KF=0 for  Y0(z) & Y1(z0)
!                        KF=1 for  Y1(z) & Y0(z1)
!                        KF=2 for  Y1'(z) & Y1(z1')
!                 KC --- Choice code
!                        KC=0 for complex roots
!                        KC=1 for real roots
!       Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
!                 ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
!                           at the L-th zero
!       Routine called: CY01 for computing Y0(z) and Y1(z), and
!                       their derivatives
!       ===========================================================

INTEGER, INTENT(IN)           :: nt
INTEGER, INTENT(IN)           :: kf
INTEGER, INTENT(IN)           :: kc
COMPLEX (dp), INTENT(IN OUT)  :: zo(nt)
COMPLEX (dp), INTENT(OUT)     :: zv(nt)

COMPLEX (dp)  :: z, zd, zero, zf, zfd, zgd, zp, zq, zw
REAL (dp)     :: h, w, w0, x, y
INTEGER       :: i, it, j, nr

IF (kc == 0) THEN
  x = -2.4D0
  y = 0.54D0
  h = 3.14
ELSE IF (kc == 1) THEN
  x = 0.89
  y = 0.0
  h = -3.14
END IF
IF (kf == 1) x = -0.503
IF (kf == 2) x = 0.577
zero = CMPLX(x, y, KIND=dp)
z = zero
w = z
DO  nr = 1, nt
  IF (nr /= 1) z = zo(nr-1) - h
  it = 0

  10 it = it + 1
  CALL cy01(kf, z, zf, zd)
  zp = (1.0D0,0.0D0)
  DO  i = 1, nr - 1
    zp = zp * (z-zo(i))
  END DO
  zfd = zf / zp
  zq = (0.0D0,0.0D0)
  DO  i = 1, nr - 1
    zw = (1.0D0,0.0D0)
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
DO  i = 1, nt
  z = zo(i)
  IF (kf == 0 .OR. kf == 2) THEN
    CALL cy01(1, z, zf, zd)
    zv(i) = zf
  ELSE IF (kf == 1) THEN
    CALL cy01(0, z, zf, zd)
    zv(i) = zf
  END IF
END DO
RETURN
END SUBROUTINE cyzo



SUBROUTINE cy01(kf, z, zf, zd)

!       ===========================================================
!       Purpose: Compute complex Bessel functions Y0(z), Y1(z)
!                and their derivatives
!       Input :  z  --- Complex argument of Yn(z) ( n=0,1 )
!                KF --- Function choice code
!                    KF=0 for ZF=Y0(z) and ZD=Y0'(z)
!                    KF=1 for ZF=Y1(z) and ZD=Y1'(z)
!                    KF=2 for ZF=Y1'(z) and ZD=Y1''(z)
!       Output:  ZF --- Y0(z) or Y1(z) or Y1'(z)
!                ZD --- Y0'(z) or Y1'(z) or Y1''(z)
!       ===========================================================

INTEGER, INTENT(IN)        :: kf
COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: zf
COMPLEX (dp), INTENT(OUT)  :: zd

COMPLEX (dp)  :: cbj0, cbj1, cby0, cby1, cdy0, cdy1, ci, cp, cp0, cp1,  &
                 cq0, cq1, cr, cs, ct1, ct2, cu, z1, z2
REAL (dp)     :: a0, rp2, w0, w1
INTEGER       :: k, k0
REAL (dp), PARAMETER  :: pi = 3.141592653589793D0, el = 0.5772156649015329D0
REAL (dp), PARAMETER  :: a(12) = (/ -.703125D-01, .112152099609375D+00, -  &
      .5725014209747314D+00, .6074042001273483D+01, -  &
      .1100171402692467D+03, .3038090510922384D+04, -  &
      .1188384262567832D+06, .6252951493434797D+07, -  &
      .4259392165047669D+09, .3646840080706556D+11, -  &
      .3833534661393944D+13, .4854014686852901D+15 /)
REAL (dp), PARAMETER  :: b(12) = (/ .732421875D-01, -.2271080017089844D+00,  &
      .1727727502584457D+01, -.2438052969955606D+02,  &
      .5513358961220206D+03, -.1825775547429318D+05,  &
      .8328593040162893D+06, -.5006958953198893D+08,  &
      .3836255180230433D+10, -.3649010818849833D+12,  &
      .4218971570284096D+14, -.5827244631566907D+16 /)
REAL (dp), PARAMETER  :: a1(12) = (/ .1171875D+00, -.144195556640625D+00,  &
      .6765925884246826D+00, -.6883914268109947D+01,  &
      .1215978918765359D+03, -.3302272294480852D+04,  &
      .1276412726461746D+06, -.6656367718817688D+07,  &
      .4502786003050393D+09, -.3833857520742790D+11,  &
      .4011838599133198D+13, -.5060568503314727D+15 /)
REAL (dp), PARAMETER  :: b1(12) = (/ -.1025390625D+00, .2775764465332031D+00, -  &
      .1993531733751297D+01, .2724882731126854D+02, -  &
      .6038440767050702D+03, .1971837591223663D+05, -  &
      .8902978767070678D+06, .5310411010968522D+08, -  &
      .4043620325107754D+10, .3827011346598605D+12, -  &
      .4406481417852278D+14, .6065091351222699D+16 /)

rp2 = 2.0D0 / pi
ci = (0.0D0,1.0D0)
a0 = ABS(z)
z2 = z * z
z1 = z
IF (a0 == 0.0D0) THEN
  cbj0 = (1.0D0,0.0D0)
  cbj1 = (0.0D0,0.0D0)
  cby0 = -(1.0D300,0.0D0)
  cby1 = -(1.0D300,0.0D0)
  cdy0 = (1.0D300,0.0D0)
  cdy1 = (1.0D300,0.0D0)
  GO TO 130
END IF

IF (REAL(z) < 0.0) z1 = -z
IF (a0 <= 12.0) THEN
  cbj0 = (1.0D0,0.0D0)
  cr = (1.0D0,0.0D0)
  DO  k = 1, 40
    cr = -0.25D0 * cr * z2 / (k*k)
    cbj0 = cbj0 + cr
    IF (ABS(cr) < ABS(cbj0)*1.0D-15) EXIT
  END DO

  cbj1 = (1.0D0,0.0D0)
  cr = (1.0D0,0.0D0)
  DO  k = 1, 40
    cr = -0.25D0 * cr * z2 / (k*(k+1.0D0))
    cbj1 = cbj1 + cr
    IF (ABS(cr) < ABS(cbj1)*1.0D-15) EXIT
  END DO

  cbj1 = 0.5D0 * z1 * cbj1
  w0 = 0.0D0
  cr = (1.0D0,0.0D0)
  cs = (0.0D0,0.0D0)
  DO  k = 1, 40
    w0 = w0 + 1.0D0 / k
    cr = -0.25D0 * cr / (k*k) * z2
    cp = cr * w0
    cs = cs + cp
    IF (ABS(cp) < ABS(cs)*1.0D-15) EXIT
  END DO

  cby0 = rp2 * (LOG(z1/2.0D0)+el) * cbj0 - rp2 * cs
  w1 = 0.0D0
  cr = (1.0D0,0.0D0)
  cs = (1.0D0,0.0D0)
  DO  k = 1, 40
    w1 = w1 + 1.0D0 / k
    cr = -0.25D0 * cr / (k*(k+1)) * z2
    cp = cr * (2.0D0*w1+1.0D0/(k+1.0D0))
    cs = cs + cp
    IF (ABS(cp) < ABS(cs)*1.0D-15) EXIT
  END DO
  cby1 = rp2 * ((LOG(z1/2.0D0)+el)*cbj1-1.0D0/z1-.25D0*z1*cs)
ELSE
  k0 = 12
  IF (a0 >= 35.0) k0 = 10
  IF (a0 >= 50.0) k0 = 8
  ct1 = z1 - .25D0 * pi
  cp0 = (1.0D0,0.0D0)
  DO  k = 1, k0
    cp0 = cp0 + a(k) * z1 ** (-2*k)
  END DO
  cq0 = -0.125D0 / z1
  DO  k = 1, k0
    cq0 = cq0 + b(k) * z1 ** (-2*k-1)
  END DO
  cu = SQRT(rp2/z1)
  cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
  cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
  ct2 = z1 - .75D0 * pi
  cp1 = (1.0D0,0.0D0)
  DO  k = 1, k0
    cp1 = cp1 + a1(k) * z1 ** (-2*k)
  END DO
  cq1 = 0.375D0 / z1
  DO  k = 1, k0
    cq1 = cq1 + b1(k) * z1 ** (-2*k-1)
  END DO
  cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
  cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
END IF
IF (REAL(z) < 0.0) THEN
  IF (AIMAG(z) < 0.0) cby0 = cby0 - 2.0D0 * ci * cbj0
  IF (AIMAG(z) > 0.0) cby0 = cby0 + 2.0D0 * ci * cbj0
  IF (AIMAG(z) < 0.0) cby1 = -(cby1-2.0D0*ci*cbj1)
  IF (AIMAG(z) > 0.0) cby1 = -(cby1+2.0D0*ci*cbj1)
  cbj1 = -cbj1
END IF
cdy0 = -cby1
cdy1 = cby0 - 1.0D0 / z * cby1
130 IF (kf == 0) THEN
  zf = cby0
  zd = cdy0
ELSE IF (kf == 1) THEN
  zf = cby1
  zd = cdy1
ELSE IF (kf == 2) THEN
  zf = cdy1
  zd = -cdy1 / z - (1.0D0-1.0D0/(z*z)) * cby1
END IF
RETURN
END SUBROUTINE cy01
 
END MODULE cyzo_func
 
 
 
PROGRAM mcyzo
USE cyzo_func
IMPLICIT NONE


! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:38

!       ===========================================================
!       Purpose : This program evaluates the complex zeros of
!                 Y0(z), Y0'(z), Y1(z) and Y1'(z), and their
!                 associated values at the zeros using the
!                 modified Newton's iteration method
!       Input:    NT --- Total number of roots/zeros
!                 KF --- Function choice code
!                        KF=0 for  Y0(z) & Y1(z0)
!                        KF=1 for  Y1(z) & Y0(z1)
!                        KF=2 for  Y1'(z) & Y1(z1')
!                 KC --- Choice code
!                        KC=0 for complex roots
!                        KC=1 for real roots
!       Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
!                 ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
!                           at the L-th zero
!       Examples: NT = 5

!   No.      z0, Zeros of Y0(z)                Y1(z0)
!  -----------------------------------------------------------------
!    1   -2.403016632 + i .5398823130   .1007476893 - i .8819677101
!    2   -5.519876702 + i .5471800106  -.0292464182 + i .5871695027
!    3   -8.653672403 + i .5484120673   .0149080637 - i .4694587524
!    4  -11.791512030 + i .5488191184  -.0093736817 + i .4023045429
!    5  -14.930906564 + i .5490008289   .0065788031 - i .3575673214

!   No.      z1, Zeros of Y1(z)                 Y0(z1)
!  -----------------------------------------------------------------
!    1    -.502743273 + i .7862437145  -.4595276847 + i1.3171019361
!    2   -3.833535193 + i .5623565382   .0483019087 - i .6925128842
!    3   -7.015903683 + i .5533930459  -.0201269494 + i .5186425332
!    4  -10.173573834 + i .5512733877   .0116140017 - i .4320329636
!    5  -13.323739307 + i .5504585830  -.0077719300 + i .3779698048

!   No.      z1', Zeros of Y1'(z)                Y1(z1')
!   ----------------------------------------------------------------
!    1     .576785129 + i .9039847922  -.7634970879 + i .5892448647
!    2   -1.940477342 - i .7211859189   .1620640057 + i .9520278864
!    3   -5.333478617 - i .5672196368  -.0317940081 - i .5968536736
!    4   -8.536768577 - i .5560607040   .0154177166 + i .4726011652
!    5  -11.706175219 - i .5528590607  -.0095443768 - i .4037533396
!       ============================================================

COMPLEX (dp)  :: zo(50), zv(50)
INTEGER       :: i, kc, kf, nt

WRITE (*,*) 'Please Enter NT, KF and KC'
WRITE (*,*) '  NT --- Total number of the roots'
WRITE (*,*) '  KF  --- Function choice code'
WRITE (*,*) '          KF=0 for Y0(z) & Y1(z0)'
WRITE (*,*) '          KF=1 for Y1(z) & Y0(z1)'
WRITE (*,*) '          KF=2 for Y1''(Z) & Y1(Z1'')'
WRITE (*,*) '  KC  --- Choice code'
WRITE (*,*) '          KC=0 for complex roots'
WRITE (*,*) '          KC=1 for real roots: '
READ (*,*) nt, kf, kc
WRITE (*,5100) nt, kf, kc
WRITE (*,*)
WRITE (*,5000)
CALL cyzo(nt, kf, kc, zo, zv)
WRITE (*,*)
IF (kf == 0) THEN
  WRITE (*,*) ' No.          z0, Zeros of Y0(z)', '                 Y1(z0)'
ELSE IF (kf == 1) THEN
  WRITE (*,*) ' No.          z1, Zeros of Y1(z)', '                 Y0(z1)'
ELSE IF (kf == 2) THEN
  WRITE (*,*) ' No.        z1'', ZEROS OF Y1''(Z)', '                Y1(z1'')'
END IF
WRITE (*,*) '------------------------------------------------------------------'
DO  i = 1, nt
  WRITE (*,5200) i, zo(i), zv(i)
END DO
STOP

5000 FORMAT (t21, '*****    Please Wait !    *****')
5100 FORMAT ('  NT=', i3, ',  KF=', i3, ',  KC=', i3)
5200 FORMAT (' ', i3, '  ', f15.9, 3F15.10)
END PROGRAM mcyzo
