MODULE cjyva_func
 
! From the book "Computation of Special Functions"
!      by Shanjie Zhang and Jianming Jin
!   Copyright 1996 by John Wiley & Sons, Inc.
! The authors state:
!   "However, we give permission to the reader who purchases this book
!    to incorporate any of these programs into his or her programs
!    provided that the copyright is acknowledged."

! The output from running the driver program does NOT agree with the
! the output contained in the comments.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
 
CONTAINS


SUBROUTINE cjyva(v, z, vm, cbj, cdj, cby, cdy)

!       ===========================================================
!       Purpose: Compute Bessel functions Jv(z), Yv(z) and their
!                derivatives for a complex argument
!       Input :  z --- Complex argument
!                v --- Order of Jv(z) and Yv(z)
!                      ( v = n+v0, n = 0,1,2,..., 0 ף v0 < 1 )
!       Output:  CBJ(n) --- Jn+v0(z)
!                CDJ(n) --- Jn+v0'(z)
!                CBY(n) --- Yn+v0(z)
!                CDY(n) --- Yn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ===========================================================

REAL (dp), INTENT(IN)      :: v
COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(OUT)     :: vm
COMPLEX (dp), INTENT(OUT)  :: cbj(0:)
COMPLEX (dp), INTENT(OUT)  :: cdj(0:)
COMPLEX (dp), INTENT(OUT)  :: cby(0:)
COMPLEX (dp), INTENT(OUT)  :: cdy(0:)

REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp, rp2 = .63661977236758_dp
COMPLEX (dp)  :: ca0, cb, cck, cec, cf, cf0, cf1, cf2, cfac0, cfac1, cg0, cg1,  &
                 ch0, ch1, ch2, ci, cju0, cju1, cjv0, cjv1, cjvl, &
                 cp11, cp12, cp21, cp22, cpz, cqz, cr, cr0, cr1, crp, crq, &
                 cs, cs0, cs1, csk, cyk, cyl1, cyl2, cylk, cyv0, cyv1, z1, z2, zk
REAL (dp)     :: a0, ca, ga, gb, pv0, pv1, v0, vl, vg, vv, w0, w1, wa,  &
                 ya0, ya1, yak
INTEGER       :: j, k, k0, l, lb, lb0, m, n

ci = (0.0_dp,1.0_dp)
a0 = ABS(z)
z1 = z
z2 = z * z
n = INT(v)
v0 = v - n
pv0 = pi * v0
pv1 = pi * (1.0_dp + v0)
IF (a0 < 1.0D-100) THEN
  DO  k = 0, n
    cbj(k) = (0.0_dp,0.0_dp)
    cdj(k) = (0.0_dp,0.0_dp)
    cby(k) = -(1.0D+300,0.0_dp)
    cdy(k) = (1.0D+300,0.0_dp)
  END DO
  IF (v0 == 0.0) THEN
    cbj(0) = (1.0_dp,0.0_dp)
    cdj(1) = (0.5_dp,0.0_dp)
  ELSE
    cdj(0) = (1.0D+300,0.0_dp)
  END IF
  vm = v
  RETURN
END IF
IF (REAL(z) < 0.0) z1 = -z
IF (a0 <= 12.0) THEN
  DO  l = 0, 1
    vl = v0 + l
    cjvl = (1.0_dp,0.0_dp)
    cr = (1.0_dp,0.0_dp)
    DO  k = 1, 40
      cr = -0.25_dp * cr * z2 / (k*(k+vl))
      cjvl = cjvl + cr
      IF (ABS(cr) < ABS(cjvl)*1.0D-15) EXIT
    END DO

    vg = 1.0_dp + vl
    CALL gamma(vg,ga)
    ca = (0.5_dp*z1) ** vl / ga
    IF (l == 0) cjv0 = cjvl * ca
    IF (l == 1) cjv1 = cjvl * ca
  END DO
ELSE
  k0 = 11
  IF (a0 >= 35.0) k0 = 10
  IF (a0 >= 50.0) k0 = 8
  DO  j = 0, 1
    vv = 4.0_dp * (j+v0) * (j+v0)
    cpz = (1.0_dp,0.0_dp)
    crp = (1.0_dp,0.0_dp)
    DO  k = 1, k0
      crp = -0.78125D-2 * crp * (vv - (4*k-3)**2) * (vv - (4*k-1)**2) / (k*(2*k-1)*z2)
      cpz = cpz + crp
    END DO
    cqz = (1.0_dp,0.0_dp)
    crq = (1.0_dp,0.0_dp)
    DO  k = 1, k0
      crq = -0.78125D-2 * crq * (vv - (4*k-1)**2) * (vv - (4*k+1)**2) / (k*(2*k+1)*z2)
      cqz = cqz + crq
    END DO
    cqz = 0.125_dp * (vv-1.0) * cqz / z1
    zk = z1 - (0.5_dp*(j+v0) + 0.25_dp) * pi
    ca0 = SQRT(rp2/z1)
    cck = COS(zk)
    csk = SIN(zk)
    IF (j == 0) THEN
      cjv0 = ca0 * (cpz*cck - cqz*csk)
      cyv0 = ca0 * (cpz*csk + cqz*cck)
    ELSE IF (j == 1) THEN
      cjv1 = ca0 * (cpz*cck - cqz*csk)
      cyv1 = ca0 * (cpz*csk + cqz*cck)
    END IF
  END DO
END IF
IF (a0 <= 12.0) THEN
  IF (v0 /= 0.0) THEN
    DO  l = 0, 1
      vl = v0 + l
      cjvl = (1.0_dp,0.0_dp)
      cr = (1.0_dp,0.0_dp)
      DO  k = 1, 40
        cr = -0.25_dp * cr * z2 / (k*(k-vl))
        cjvl = cjvl + cr
        IF (ABS(cr) < ABS(cjvl)*1.0D-15) EXIT
      END DO

      vg = 1.0_dp - vl
      CALL gamma(vg, gb)
      cb = (2.0_dp/z1) ** vl / gb
      IF (l == 0) cju0 = cjvl * cb
      IF (l == 1) cju1 = cjvl * cb
    END DO
    cyv0 = (cjv0*COS(pv0) - cju0) / SIN(pv0)
    cyv1 = (cjv1*COS(pv1) - cju1) / SIN(pv1)
  ELSE
    cec = LOG(z1/2.0_dp) + .5772156649015329_dp
    cs0 = (0.0_dp,0.0_dp)
    w0 = 0.0_dp
    cr0 = (1.0_dp,0.0_dp)
    DO  k = 1, 30
      w0 = w0 + 1.0_dp / k
      cr0 = -0.25_dp * cr0 / (k*k) * z2
      cs0 = cs0 + cr0 * w0
    END DO
    cyv0 = rp2 * (cec*cjv0 - cs0)
    cs1 = (1.0_dp,0.0_dp)
    w1 = 0.0_dp
    cr1 = (1.0_dp,0.0_dp)
    DO  k = 1, 30
      w1 = w1 + 1.0_dp / k
      cr1 = -0.25_dp * cr1 / (k*(k+1)) * z2
      cs1 = cs1 + cr1 * (2.0_dp*w1 + 1.0_dp/(k+1))
    END DO
    cyv1 = rp2 * (cec*cjv1 - 1.0_dp/z1 - 0.25_dp*z1*cs1)
  END IF
END IF
IF (REAL(z) < 0.0_dp) THEN
  cfac0 = EXP(pv0*ci)
  cfac1 = EXP(pv1*ci)
  IF (AIMAG(z) < 0.0_dp) THEN
    cyv0 = cfac0 * cyv0 - 2.0_dp * ci * COS(pv0) * cjv0
    cyv1 = cfac1 * cyv1 - 2.0_dp * ci * COS(pv1) * cjv1
    cjv0 = cjv0 / cfac0
    cjv1 = cjv1 / cfac1
  ELSE IF (AIMAG(z) > 0.0_dp) THEN
    cyv0 = cyv0 / cfac0 + 2.0_dp * ci * COS(pv0) * cjv0
    cyv1 = cyv1 / cfac1 + 2.0_dp * ci * COS(pv1) * cjv1
    cjv0 = cfac0 * cjv0
    cjv1 = cfac1 * cjv1
  END IF
END IF
cbj(0) = cjv0
cbj(1) = cjv1
IF (n >= 2 .AND. n <= INT(0.25*a0)) THEN
  cf0 = cjv0
  cf1 = cjv1
  DO  k = 2, n
    cf = 2.0_dp * (k+v0-1.0_dp) / z * cf1 - cf0
    cbj(k) = cf
    cf0 = cf1
    cf1 = cf
  END DO
ELSE IF (n >= 2) THEN
  m = msta1(a0, 200)
  IF (m < n) THEN
    n = m
  ELSE
    m = msta2(a0, n, 15)
  END IF
  cf2 = (0.0_dp,0.0_dp)
  cf1 = (1.0D-100,0.0_dp)
  DO  k = m, 0, -1
    cf = 2.0_dp * (v0 + k + 1) / z * cf1 - cf2
    IF (k <= n) cbj(k) = cf
    cf2 = cf1
    cf1 = cf
  END DO
  IF (ABS(cjv0) > ABS(cjv1)) cs = cjv0 / cf
  IF (ABS(cjv0) <= ABS(cjv1)) cs = cjv1 / cf2
  DO  k = 0, n
    cbj(k) = cs * cbj(k)
  END DO
END IF
cdj(0) = v0 / z * cbj(0) - cbj(1)
DO  k = 1, n
  cdj(k) = -(k+v0) / z * cbj(k) + cbj(k-1)
END DO
cby(0) = cyv0
cby(1) = cyv1
ya0 = ABS(cyv0)
lb = 0
cg0 = cyv0
cg1 = cyv1
DO  k = 2, n
  cyk = 2.0_dp * (v0+k-1) / z * cg1 - cg0
  IF (ABS(cyk) <= 1.0D+290) THEN
    yak = ABS(cyk)
    ya1 = ABS(cg0)
    IF (yak < ya0 .AND. yak < ya1) lb = k
    cby(k) = cyk
    cg0 = cg1
    cg1 = cyk
  END IF
END DO
IF (lb > 4 .AND. AIMAG(z) /= 0.0_dp) THEN
  180   IF (lb /= lb0) THEN
    ch2 = (1.0_dp,0.0_dp)
    ch1 = (0.0_dp,0.0_dp)
    lb0 = lb
    DO  k = lb, 1, -1
      ch0 = 2.0_dp * (k+v0) / z * ch1 - ch2
      ch2 = ch1
      ch1 = ch0
    END DO
    cp12 = ch0
    cp22 = ch2
    ch2 = (0.0_dp,0.0_dp)
    ch1 = (1.0_dp,0.0_dp)
    DO  k = lb, 1, -1
      ch0 = 2.0_dp * (k+v0) / z * ch1 - ch2
      ch2 = ch1
      ch1 = ch0
    END DO
    cp11 = ch0
    cp21 = ch2
    IF (lb == n) cbj(lb+1) = 2.0_dp * (lb+v0) / z * cbj(lb) - cbj(lb-1)
    IF (ABS(cbj(0)) > ABS(cbj(1))) THEN
      cby(lb+1) = (cbj(lb+1)*cyv0 - 2.0_dp*cp11/(pi*z)) / cbj(0)
      cby(lb) = (cbj(lb)*cyv0 + 2.0_dp*cp12/(pi*z)) / cbj(0)
    ELSE
      cby(lb+1) = (cbj(lb+1)*cyv1 - 2.0_dp*cp21/(pi*z)) / cbj(1)
      cby(lb) = (cbj(lb)*cyv1 + 2.0_dp*cp22/(pi*z)) / cbj(1)
    END IF
    cyl2 = cby(lb+1)
    cyl1 = cby(lb)
    DO  k = lb - 1, 0, -1
      cylk = 2.0_dp * (k+v0+1.0_dp) / z * cyl1 - cyl2
      cby(k) = cylk
      cyl2 = cyl1
      cyl1 = cylk
    END DO
    cyl1 = cby(lb)
    cyl2 = cby(lb+1)
    DO  k = lb + 1, n - 1
      cylk = 2.0_dp * (k+v0) / z * cyl2 - cyl1
      cby(k+1) = cylk
      cyl1 = cyl2
      cyl2 = cylk
    END DO
    DO  k = 2, n
      wa = ABS(cby(k))
      IF (wa < ABS(cby(k-1))) lb = k
    END DO
    GO TO 180
  END IF
END IF
cdy(0) = v0 / z * cby(0) - cby(1)
DO  k = 1, n
  cdy(k) = cby(k-1) - (k+v0) / z * cby(k)
END DO
vm = n + v0
RETURN
END SUBROUTINE cjyva



SUBROUTINE gamma(x, ga)

!       ==================================================
!       Purpose: Compute gamma function ג(x)
!       Input :  x  --- Argument of ג(x)
!                       ( x is not equal to 0,-1,-2,תתת)
!       Output:  GA --- ג(x)
!       ==================================================

REAL (dp), INTENT(IN)      :: x
REAL (dp), INTENT(OUT)     :: ga

REAL (dp), PARAMETER  :: g(26) = (/   &
      1.0_dp, 0.5772156649015329_dp, -0.6558780715202538_dp,  &
      -0.420026350340952D-1, 0.1665386113822915_dp, -  &
      .421977345555443D-1, -.96219715278770D-2,  &
      .72189432466630D-2, -.11651675918591D-2, -.2152416741149D-3,  &
      .1280502823882D-3, -.201348547807D-4, -.12504934821D-5,  &
      .11330272320D-5, -.2056338417D-6, .61160950D-8,  &
      .50020075D-8, -.11812746D-8, .1043427D-9, .77823D-11,  &
      -.36968D-11, .51D-12, -.206D-13, -.54D-14, .14D-14, .1D-15 /)
REAL (dp), PARAMETER  :: pi = 3.141592653589793_dp
INTEGER    :: k, m, m1
REAL (dp)  :: gr, r, z

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

END MODULE cjyva_func
 
 
 
PROGRAM mcjyva
USE cjyva_func
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2001-12-25  Time: 11:55:36

!       ===============================================================
!       Purpose: This program computes Bessel functions Jv(z), Yv(z),
!                and their derivatives for a complex argument using
!                subroutine CJYVA
!       Input :  z --- Complex argument
!                v --- Order of Jv(z) and Yv(z)
!                      ( v = n+v0, 0 ף n ף 250, 0 ף v0 < 1 )
!       Output:  CBJ(n) --- Jn+v0(z)
!                CDJ(n) --- Jn+v0'(z)
!                CBY(n) --- Yn+v0(z)
!                CDY(n) --- Yn+v0'(z)
!       Example:
!                v = n +v0,  v0 = 1/3,   z = 4.0 + i 2.0

!     n     Re[Jv(z)]       Im[Jv(z)]      Re[Jv'(z)]      Im[Jv'(z)]
!    ------------------------------------------------------------------
!     0  -.13829878D+01  -.30855145D+00  -.18503756D+00   .13103689D+01
!     1   .82553327D-01  -.12848394D+01  -.12336901D+01   .45079506D-01
!     2   .10843924D+01  -.39871046D+00  -.33046401D+00  -.84574964D+00
!     3   .74348135D+00   .40665987D+00   .45318486D+00  -.42198992D+00
!     4   .17802266D+00   .44526939D+00   .39624497D+00   .97902890D-01
!     5  -.49008598D-01   .21085409D+00   .11784299D+00   .19422044D+00

!     n     Re[Yv(z)]      Im[Yv(z)]       Re[Yv'(z)]      Im[Yv'(z)]
!    ------------------------------------------------------------------
!     0   .34099851D+00  -.13440666D+01  -.13544477D+01  -.15470699D+00
!     1   .13323787D+01   .53735934D-01  -.21467271D-01  -.11807457D+01
!     2   .38393305D+00   .10174248D+01   .91581083D+00  -.33147794D+00
!     3  -.49924295D+00   .71669181D+00   .47786442D+00   .37321597D+00
!     4  -.57179578D+00   .27099289D+00  -.12111686D+00   .23405313D+00
!     5  -.25700924D+00   .24858555D+00  -.43023156D+00  -.13123662D+00
!       ===============================================================

COMPLEX (dp)  :: cbj(0:251), cdj(0:251), cby(0:251), cdy(0:251), z
REAL (dp)     :: v, v0, vm, x, y
INTEGER       :: k, n, nm, ns

WRITE (*,*) '  Please enter v, x and y ( z=x+iy ): '
READ (*,*) v, x, y
z = CMPLX(x, y, KIND=dp)
n = INT(v)
v0 = v - n
WRITE (*,5100) v0, x, y
IF (n <= 8) THEN
  ns = 1
ELSE
  WRITE (*,*) '  Please enter order step Ns: '
  READ (*,*) ns
END IF
CALL cjyva(v, z, vm, cbj, cdj, cby, cdy)
nm = INT(vm)
WRITE(*,*)
WRITE(*,*) '  n       Re[Jv(z)]       Im[Jv(z)]       Re[Jv''(Z)]      IM[JV''(Z)]'
WRITE(*,*) ' ---------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, cbj(k), cdj(k)
END DO
WRITE(*,*)
WRITE(*,*) '  n       Re[Yv(z)]       Im[Yv(z)]       Re[Yv''(Z)]      IM[YV''(Z)]'
WRITE(*,*) ' ---------------------------------------------------------------------'
DO  k = 0, nm, ns
  WRITE (*,5000) k, cby(k), cdy(k)
END DO
STOP

5000 FORMAT (' ', i3, '  ', 4g16.8)
5100 FORMAT (t9, 'v = n+v0,  v0 =', f5.2, ',  z =', f7.2, ' +', f7.2, 'i')
END PROGRAM mcjyva
