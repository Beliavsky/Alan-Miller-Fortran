MODULE cjyna_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! Latest revision - 31 December 2001
! Corrections by Alan Miller (amiller @ bigpond.net.au)
! Variable lb0 not initialized.
! Array bounds exceeded in routine CJYNA for arrays CBJ, CDJ, CBY & CDY.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE cjyna(n, z, nm, cbj, cdj, cby, cdy)

!       =======================================================
!       Purpose: Compute Bessel functions Jn(z), Yn(z) and
!                their derivatives for a complex argument
!       Input :  z --- Complex argument of Jn(z) and Yn(z)
!                n --- Order of Jn(z) and Yn(z)
!       Output:  CBJ(n) --- Jn(z)
!                CDJ(n) --- Jn'(z)
!                CBY(n) --- Yn(z)
!                CDY(n) --- Yn'(z)
!                NM --- Highest order computed
!       Routines called:
!            (1) CJY01 to calculate J0(z), J1(z), Y0(z), Y1(z)
!            (2) MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       =======================================================

INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(OUT)       :: nm
COMPLEX (dp), INTENT(OUT)  :: cbj(0:)
COMPLEX (dp), INTENT(OUT)  :: cdj(0:)
COMPLEX (dp), INTENT(OUT)  :: cby(0:)
COMPLEX (dp), INTENT(OUT)  :: cdy(0:)

COMPLEX (dp)  :: cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1,  &
                 cf, cf1, cf2, cg0, cg1, ch0, ch1, ch2, cj0, cj1, cjk,  &
                 cp11, cp12, cp21, cp22, cs, cyk, cyl1, cyl2, cylk
REAL (dp)     :: a0, wa, ya0, ya1, yak
INTEGER       :: k, lb, lb0, m
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp

a0 = ABS(z)
nm = n
IF (a0 < 1.0D-100) THEN
  DO  k = 0, n
    cbj(k) = (0.0_dp,0.0_dp)
    cdj(k) = (0.0_dp,0.0_dp)
    cby(k) = -(1.0D+300,0.0_dp)
    cdy(k) = (1.0D+300,0.0_dp)
  END DO
  cbj(0) = (1.0_dp,0.0_dp)
  cdj(1) = (0.5_dp,0.0_dp)
  RETURN
END IF
CALL cjy01(z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1)
cbj(0) = cbj0
cbj(1) = cbj1
cby(0) = cby0
cby(1) = cby1
cdj(0) = cdj0
cdj(1) = cdj1
cdy(0) = cdy0
cdy(1) = cdy1
IF (n <= 1) RETURN
IF (n < INT(0.25*a0)) THEN
  cj0 = cbj0
  cj1 = cbj1
  DO  k = 2, n
    cjk = 2 * (k-1) / z * cj1 - cj0
    cbj(k) = cjk
    cj0 = cj1
    cj1 = cjk
  END DO
ELSE
  m = msta1(a0, 200)
  IF (m < n) THEN
    nm = m
  ELSE
    m = msta2(a0, n, 15)
  END IF
  cf2 = (0.0_dp,0.0_dp)
  cf1 = (1.0D-100,0.0_dp)
  DO  k = m, 0, -1
    cf = 2 * (k+1) / z * cf1 - cf2
    IF (k <= nm) cbj(k) = cf
    cf2 = cf1
    cf1 = cf
  END DO
  IF (ABS(cbj0) > ABS(cbj1)) THEN
    cs = cbj0 / cf
  ELSE
    cs = cbj1 / cf2
  END IF
  DO  k = 0, nm
    cbj(k) = cs * cbj(k)
  END DO
END IF
DO  k = 2, nm
  cdj(k) = cbj(k-1) - k / z * cbj(k)
END DO
ya0 = ABS(cby0)
lb = 0
cg0 = cby0
cg1 = cby1
DO  k = 2, nm
  cyk = 2 * (k-1) / z * cg1 - cg0
  IF (ABS(cyk) <= 1.0D+290) THEN
    yak = ABS(cyk)
    ya1 = ABS(cg0)
    IF (yak < ya0 .AND. yak < ya1) lb = k
    cby(k) = cyk
    cg0 = cg1
    cg1 = cyk
  END IF
END DO
lb0 = 0
IF (lb > 4 .AND. AIMAG(z) /= 0.0_dp) THEN
  70 IF (lb /= lb0) THEN
    ch2 = (1.0_dp,0.0_dp)
    ch1 = (0.0_dp,0.0_dp)
    lb0 = lb
    DO  k = lb, 1, -1
      ch0 = 2.0_dp * k / z * ch1 - ch2
      ch2 = ch1
      ch1 = ch0
    END DO
    cp12 = ch0
    cp22 = ch2
    ch2 = (0.0_dp,0.0_dp)
    ch1 = (1.0_dp,0.0_dp)
    DO  k = lb, 1, -1
      ch0 = 2 * k / z * ch1 - ch2
      ch2 = ch1
      ch1 = ch0
    END DO
    cp11 = ch0
    cp21 = ch2
    IF (lb == nm) cbj(lb+1) = 2.0_dp * lb / z * cbj(lb) - cbj(lb-1)
    IF (ABS(cbj(0)) > ABS(cbj(1))) THEN
      cby(lb+1) = (cbj(lb+1)*cby0 - 2.0_dp*cp11/(pi*z)) / cbj(0)
      cby(lb) = (cbj(lb)*cby0 + 2.0_dp*cp12/(pi*z)) / cbj(0)
    ELSE
      cby(lb+1) = (cbj(lb+1)*cby1 - 2.0_dp*cp21/(pi*z)) / cbj(1)
      cby(lb) = (cbj(lb)*cby1 + 2.0_dp*cp22/(pi*z)) / cbj(1)
    END IF
    cyl2 = cby(lb+1)
    cyl1 = cby(lb)
    DO  k = lb - 1, 0, -1
      cylk = 2 * (k+1) / z * cyl1 - cyl2
      cby(k) = cylk
      cyl2 = cyl1
      cyl1 = cylk
    END DO
    cyl1 = cby(lb)
    cyl2 = cby(lb+1)
    DO  k = lb + 1, nm - 1
      cylk = 2 * k / z * cyl2 - cyl1
      cby(k+1) = cylk
      cyl1 = cyl2
      cyl2 = cylk
    END DO
    DO  k = 2, nm
      wa = ABS(cby(k))
      IF (wa < ABS(cby(k-1))) lb = k
    END DO
    GO TO 70
  END IF
END IF
DO  k = 2, nm
  cdy(k) = cby(k-1) - k / z * cby(k)
END DO
RETURN
END SUBROUTINE cjyna



SUBROUTINE cjy01(z, cbj0, cdj0, cbj1, cdj1, cby0, cdy0, cby1, cdy1)

!       ===========================================================
!       Purpose: Compute complex Bessel functions J0(z), J1(z)
!                Y0(z), Y1(z), and their derivatives
!       Input :  z --- Complex argument
!       Output:  CBJ0 --- J0(z)
!                CDJ0 --- J0'(z)
!                CBJ1 --- J1(z)
!                CDJ1 --- J1'(z)
!                CBY0 --- Y0(z)
!                CDY0 --- Y0'(z)
!                CBY1 --- Y1(z)
!                CDY1 --- Y1'(z)
!       ===========================================================

COMPLEX (dp), INTENT(IN)   :: z
COMPLEX (dp), INTENT(OUT)  :: cbj0
COMPLEX (dp), INTENT(OUT)  :: cdj0
COMPLEX (dp), INTENT(OUT)  :: cbj1
COMPLEX (dp), INTENT(OUT)  :: cdj1
COMPLEX (dp), INTENT(OUT)  :: cby0
COMPLEX (dp), INTENT(OUT)  :: cdy0
COMPLEX (dp), INTENT(OUT)  :: cby1
COMPLEX (dp), INTENT(OUT)  :: cdy1

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
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, el = 0.5772156649015329_dp
COMPLEX (dp)  :: ci, cp, cp0, cp1, cq0, cq1, cr, cs, ct1, ct2, cu, z1, z2
REAL (dp)     :: a0, rp2, w0, w1
INTEGER       :: k, k0

rp2 = 2.0_dp / pi
ci = (0.0_dp,1.0_dp)
a0 = ABS(z)
z2 = z * z
z1 = z
IF (a0 == 0.0_dp) THEN
  cbj0 = (1.0_dp,0.0_dp)
  cbj1 = (0.0_dp,0.0_dp)
  cdj0 = (0.0_dp,0.0_dp)
  cdj1 = (0.5_dp,0.0_dp)
  cby0 = -(1.0D300,0.0_dp)
  cby1 = -(1.0D300,0.0_dp)
  cdy0 = (1.0D300,0.0_dp)
  cdy1 = (1.0D300,0.0_dp)
  RETURN
END IF
IF (REAL(z) < 0.0) z1 = -z
IF (a0 <= 12.0) THEN
  cbj0 = (1.0_dp,0.0_dp)
  cr = (1.0_dp,0.0_dp)
  DO  k = 1, 40
    cr = -0.25_dp * cr * z2 / (k*k)
    cbj0 = cbj0 + cr
    IF (ABS(cr/cbj0) < 1.0D-15) EXIT
  END DO

  cbj1 = (1.0_dp,0.0_dp)
  cr = (1.0_dp,0.0_dp)
  DO  k = 1, 40
    cr = -0.25_dp * cr * z2 / (k*(k+1))
    cbj1 = cbj1 + cr
    IF (ABS(cr/cbj1) < 1.0D-15) EXIT
  END DO

  cbj1 = 0.5_dp * z1 * cbj1
  w0 = 0.0_dp
  cr = (1.0_dp,0.0_dp)
  cs = (0.0_dp,0.0_dp)
  DO  k = 1, 40
    w0 = w0 + 1.0_dp / k
    cr = -0.25_dp * cr / (k*k) * z2
    cp = cr * w0
    cs = cs + cp
    IF (ABS(cp/cs) < 1.0D-15) EXIT
  END DO
  cby0 = rp2 * (LOG(z1/2.0_dp) + el) * cbj0 - rp2 * cs
  w1 = 0.0_dp
  cr = (1.0_dp,0.0_dp)
  cs = (1.0_dp,0.0_dp)
  DO  k = 1, 40
    w1 = w1 + 1.0_dp / k
    cr = -0.25_dp * cr / (k*(k+1)) * z2
    cp = cr * (2.0_dp*w1 + 1.0_dp/(k+1))
    cs = cs + cp
    IF (ABS(cp/cs) < 1.0D-15) EXIT
  END DO
  cby1 = rp2 * ((LOG(z1/2.0_dp) + el)*cbj1 - 1.0_dp/z1 - .25_dp*z1*cs)
ELSE
  k0 = 12
  IF (a0 >= 35.0) k0 = 10
  IF (a0 >= 50.0) k0 = 8
  ct1 = z1 - 0.25_dp * pi
  cp0 = (1.0_dp,0.0_dp)
  DO  k = 1, k0
    cp0 = cp0 + a(k) * z1 ** (-2*k)
  END DO
  cq0 = -0.125_dp / z1
  DO  k = 1, k0
    cq0 = cq0 + b(k) * z1 ** (-2*k-1)
  END DO
  cu = SQRT(rp2/z1)
  cbj0 = cu * (cp0*COS(ct1) - cq0*SIN(ct1))
  cby0 = cu * (cp0*SIN(ct1) + cq0*COS(ct1))
  ct2 = z1 - 0.75_dp * pi
  cp1 = (1.0_dp,0.0_dp)
  DO  k = 1, k0
    cp1 = cp1 + a1(k) * z1 ** (-2*k)
  END DO
  cq1 = 0.375_dp / z1
  DO  k = 1, k0
    cq1 = cq1 + b1(k) * z1 ** (-2*k-1)
  END DO
  cbj1 = cu * (cp1*COS(ct2) - cq1*SIN(ct2))
  cby1 = cu * (cp1*SIN(ct2) + cq1*COS(ct2))
END IF
IF (REAL(z) < 0.0) THEN
  IF (AIMAG(z) < 0.0) cby0 = cby0 - 2.0_dp * ci * cbj0
  IF (AIMAG(z) > 0.0) cby0 = cby0 + 2.0_dp * ci * cbj0
  IF (AIMAG(z) < 0.0) cby1 = -(cby1 - 2.0_dp*ci*cbj1)
  IF (AIMAG(z) > 0.0) cby1 = -(cby1 + 2.0_dp*ci*cbj1)
  cbj1 = -cbj1
END IF
cdj0 = -cbj1
cdj1 = cbj0 - 1.0_dp / z * cbj1
cdy0 = -cby1
cdy1 = cby0 - 1.0_dp / z * cby1
RETURN
END SUBROUTINE cjy01



FUNCTION msta1(x, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that the magnitude of
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point
!       ===================================================

REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(IN)        :: mp
INTEGER                    :: fn_val

REAL (dp)  :: a0, f, f0, f1
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
n0 = INT(1.1*a0) + 1
f0 = envj(n0,a0) - mp
n1 = n0 + 5
f1 = envj(n1,a0) - mp
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dp - f0/f1)
  f = envj(nn,a0) - mp
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn
RETURN
END FUNCTION msta1



FUNCTION msta2(x, n, mp) RESULT(fn_val)

!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================

REAL (dp), INTENT(IN)      :: x
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: mp
INTEGER                    :: fn_val

REAL (dp)  :: a0, ejn, f, f0, f1, hmp, obj
INTEGER    :: it, n0, n1, nn

a0 = ABS(x)
hmp = 0.5_dp * mp
ejn = envj(n, a0)
IF (ejn <= hmp) THEN
  obj = mp
  n0 = INT(1.1*a0)
ELSE
  obj = hmp + ejn
  n0 = n
END IF
f0 = envj(n0,a0) - obj
n1 = n0 + 5
f1 = envj(n1,a0) - obj
DO  it = 1, 20
  nn = n1 - (n1-n0) / (1.0_dp - f0/f1)
  f = envj(nn, a0) - obj
  IF (ABS(nn-n1) < 1) EXIT
  n0 = n1
  f0 = f1
  n1 = nn
  f1 = f
END DO

fn_val = nn + 10
RETURN
END FUNCTION msta2



FUNCTION envj(n, x) RESULT(fn_val)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x
REAL (dp)                  :: fn_val

fn_val = 0.5_dp * LOG10(6.28_dp*n) - n * LOG10(1.36_dp*x/n)
RETURN
END FUNCTION envj

END MODULE cjyna_func
 
 
 
PROGRAM mcjyna
USE cjyna_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:36

!       ================================================================
!       Purpose: This program computes Bessel functions Jn(z), Yn(z)
!                and their derivatives for a complex argument using
!                subroutine CJYNA
!       Input :  z --- Complex argument of Jn(z) and Yn(z)
!                n --- Order of Jn(z) and Yn(z)
!                      ( n = 0,1,תתת, n ף 250 )
!       Output:  CBJ(n) --- Jn(z)
!                CDJ(n) --- Jn'(z)
!                CBY(n) --- Yn(z)
!                CDY(n) --- Yn'(z)
!       Eaxmple: z = 4.0 + i 2.0

!     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
!    -------------------------------------------------------------------
!     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
!     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
!     2   .93050039D+00  -.77959350D+00  -.72363400D+00  -.72836666D+00
!     3   .93991546D+00   .23042918D+00   .29742236D+00  -.63587637D+00
!     4   .33565567D+00   .49215925D+00   .47452722D+00  -.29035945D-01
!     5  -.91389835D-02   .28850107D+00   .20054412D+00   .19908868D+00

!     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
!   --------------------------------------------------------------------
!     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
!     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
!     2   .79074211D+00   .86842120D+00   .78932897D+00  -.70142425D+00
!     3  -.29934789D+00   .89064431D+00   .70315755D+00   .24423024D+00
!     4  -.61557299D+00   .37996071D+00   .41126221D-01   .34044655D+00
!     5  -.38160033D+00   .20975121D+00  -.33884827D+00  -.20590670D-01

!                z = 20.0 + i 10.0 ,      Nmax = 5

!     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
!   --------------------------------------------------------------------
!     0   .15460268D+04  -.10391216D+04  -.10601232D+04  -.15098284D+04
!     1   .10601232D+04   .15098284D+04   .14734253D+04  -.10783122D+04
!     2  -.14008238D+04   .11175029D+04   .11274890D+04   .13643952D+04
!     3  -.11948548D+04  -.12189620D+04  -.11843035D+04   .11920871D+04
!     4   .96778325D+03  -.12666712D+04  -.12483664D+04  -.93887194D+03
!     5   .13018781D+04   .65878188D+03   .64152944D+03  -.12682398D+04

!     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
!   --------------------------------------------------------------------
!     0   .10391216D+04   .15460268D+04   .15098284D+04  -.10601232D+04
!     1  -.15098284D+04   .10601232D+04   .10783122D+04   .14734253D+04
!     2  -.11175029D+04  -.14008238D+04  -.13643952D+04   .11274890D+04
!     3   .12189620D+04  -.11948548D+04  -.11920871D+04  -.11843035D+04
!     4   .12666712D+04   .96778324D+03   .93887194D+03  -.12483664D+04
!     5  -.65878189D+03   .13018781D+04   .12682398D+04   .64152944D+03
!       ================================================================

COMPLEX (dp)  :: cbj(0:251), cdj(0:251), cby(0:251), cdy(0:251), z
REAL (dp)     :: x, y
INTEGER       :: k, n, nm, ns

WRITE (*,*) '  Please enter n, x,y (z=x+iy) '
READ (*,*) n, x, y
z = CMPLX(x, y, KIND=dp)
WRITE (*,5100) x, y, n
IF (n <= 8) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns'
  READ (*,*) ns
END IF
CALL cjyna(n, z, nm, cbj, cdj, cby, cdy)
WRITE(*,*)
WRITE(*,*) '   n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn''(Z)]      IM[JN''(Z)]'
WRITE(*,*) ' --------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, cbj(k), cdj(k)
END DO
WRITE(*,*)
WRITE(*,*) '   n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn''(Z)]      IM[YN''(Z)]'
WRITE(*,*) ' --------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, cby(k), cdy(k)
END DO
STOP

5000 FORMAT (' ', i4, 4g16.8)
5100 FORMAT ('   z =', f5.1, ' + i ', f5.1, ' ,      Nmax =', i3)
END PROGRAM mcjyna
