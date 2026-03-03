PROGRAM erftst
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------
! FORTRAN 90 program to test ERF and related functions.

!  Method:

!     Accuracy test compare function values against local Taylor's series
!     expansions.  Derivatives for erfc(x) are expressed as repeated
!     integrals of erfc(x).  These are generated from the recurrence
!     relation using a technique due to Gautschi (see references).

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - An environmental inquiry program providing information
!              on the floating-point arithmetic system.
!              Note that the call to MACHAR can be deleted provided the
!              following five parameters are assigned the values indicated

!              IBETA  - The radix of the floating-point system
!              IT     - The number of base-ibeta digits in the
!                       significant of a floating-point number
!              EPS    - The smallest positive floating-point
!                       number such that 1.0+EPS .NE. 1.0
!              XMIN   - The smallest non-vanishing floating-point
!                       integral power of the radix
!              XMAX   - The largest finite floating-point number

!     REN(K) - A function subprogram returning random real
!              numbers uniformly distributed over (0,1)


!  Intrinsic functions required are:

!      ABS, DBLE, LOG, MAX, REAL, SQRT

!  References: "Performance evaluation of programs for the error
!               and complementary error functions", W. J. Cody,
!               submitted for publication.

!              "Evaluation of the repeated integrals of the coerror
!               function", W. Gautschi, TOMS 3, 1977, pp. 240-252.

!  Latest modification: March 12, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!------------------------------------------------------------------
USE ErrorFunction
USE Toms715_Utilities

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, irnd, it, j, k, k1, k2, k3,  &
              machep, maxexp, minexp, n, negep, ngrd, n0, n1
REAL (dp)  :: a, ait, albeta, b, beta, c, c2, del, eps, epscon,  &
              epsneg, ff, f0, r, r6, r7, sc, u, v, w, x, xbig, xc,  &
              xl, xmax, xmin, xn, xn1, z, zz
REAL (dp)  :: r1(500)
!------------------------------------------------------------------
!   C1 = 1/sqrt(pi)
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, ten = 10.0_dp, sixten = 16.0_dp,  &
                         thresh = 0.46875_dp, x99 = -999.0_dp,  &
                         c1 = 5.6418958354775628695D-1
INTEGER, PARAMETER  :: iout = 6

!------------------------------------------------------------------
!  Determine machine parameters and set constants
!------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
c = ABS(ait*albeta) + ten
a = zero
b = thresh
n = 2000
xn = n
n0 = (ibeta/2) * (it+5) / 6 + 4
!-----------------------------------------------------------------
!  Determine largest argument for ERFC test by Newton iteration
!-----------------------------------------------------------------
c2 = LOG(xmin) + LOG(one/c1)
xbig = SQRT(-c2)
10 x = xbig
f0 = x * x
ff = f0 + half / f0 + LOG(x) + c2
f0 = x + x + one / x - one / (x*f0)
xbig = x - ff / f0
IF (ABS(x-xbig)/x > ten*eps) GO TO 10
!-----------------------------------------------------------------
!  Random argument accuracy tests
!-----------------------------------------------------------------
DO  j = 1, 5
  k1 = 0
  k3 = 0
  xc = zero
  r6 = zero
  r7 = zero
  del = (b-a) / xn
  xl = a
  DO  i = 1, n
    x = del * ren() + xl
    IF (j == 1) THEN
!-----------------------------------------------------------------
!  Test erf against double series expansion
!-----------------------------------------------------------------
      f0 = n0
      ff = f0 + f0 + one
      z = x * x
      w = z + z
      u = zero
      v = zero
      DO  k = 1, n0
        u = -z / f0 * (one+u)
        v = w / ff * (one+v)
        f0 = f0 - one
        ff = ff - two
      END DO
      v = c1 * (x+x) * (((u*v + (u+v)) + half) + half)
      u = erf(x)
    ELSE
!-----------------------------------------------------------------
!  Test erfc or scaled erfc against expansion in repeated
!   integrals of the coerror function.
!-----------------------------------------------------------------
      z = x + half
      x = z - half
      r = zero
      IF (x <= one) THEN
        n0 = 499
      ELSE
        n0 = MIN(499,INT(c/(ABS(LOG(z)))))
      END IF
      n1 = n0
      xn1 = n1+1
      DO  k = 1, n0
        r = half / (z+xn1*r)
        r1(n1) = r
        n1 = n1 - 1
        xn1 = xn1 - one
      END DO
      ff = c1 / (z+r1(1))
      IF ((j/2)*2 == j) THEN
        f0 = erfc(z) * EXP(x+half*half)
        u = erfc(x)
      ELSE
        f0 = erfcx(z)
        u = erfcx(x)
      END IF
      sc = f0 / ff
!-----------------------------------------------------------------
!  Scale things to avoid premature underflow
!-----------------------------------------------------------------
      epscon = f0
      ff = sixten * ff / eps
      DO  n1 = 1, n0
        ff = r1(n1) * ff
        r1(n1) = ff * sc
        IF (r1(n1) < epscon) THEN
          k = n1
          EXIT
        END IF
      END DO
      v = r1(k)
      DO  n1 = 1, k - 1
        v = v + r1(k-n1)
      END DO
!-----------------------------------------------------------------
!  Remove scaling here
!-----------------------------------------------------------------
      v = v * eps / sixten + f0
    END IF
!--------------------------------------------------------------------
!  Accumulate results
!--------------------------------------------------------------------
    w = (u-v) / u
    IF (w > zero) THEN
      k1 = k1 + 1
    ELSE IF (w < zero) THEN
      k3 = k3 + 1
    END IF
    w = ABS(w)
    IF (w > r6) THEN
      r6 = w
      xc = x
    END IF
    r7 = r7 + w * w
    xl = xl + del
  END DO
!------------------------------------------------------------------
!  Gather and print statistics for test
!------------------------------------------------------------------
  k2 = n - k3 - k1
  r7 = SQRT(r7/xn)
  IF (j == 1) THEN
    WRITE (iout,5000)
    WRITE (iout,5300) n, a, b
    WRITE (iout,5400) k1
  ELSE IF ((j/2)*2 == j) THEN
    WRITE (iout,5100)
    WRITE (iout,5300) n, a, b
    WRITE (iout,5500) k1
  ELSE
    WRITE (iout,5200)
    WRITE (iout,5300) n, a, b
    WRITE (iout,5600) k1
  END IF
  WRITE (iout,5700) k2, k3
  WRITE (iout,5800) it, ibeta
  IF (r6 /= zero) THEN
    w = LOG(ABS(r6)) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5900) r6, ibeta, w, xc
  w = MAX(ait+w,zero)
  WRITE (iout,6000) ibeta, w
  IF (r7 /= zero) THEN
    w = LOG(ABS(r7)) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,6100) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,6000) ibeta, w
!------------------------------------------------------------------
!  Initialize for next test
!------------------------------------------------------------------
  IF (j == 1) THEN
    a = b
    b = two
  ELSE IF (j == 3) THEN
    a = b
    b = AINT(xbig*sixten) / sixten - half
  ELSE IF (j == 4) THEN
    b = ten + ten
  END IF
END DO
!-----------------------------------------------------------------
!  Special tests
!  First check values for negative arguments.
!-----------------------------------------------------------------
WRITE (iout,6200)
WRITE (iout,6300) ibeta
x = zero
del = -half
DO  i = 1, 10
  u = erf(x)
  a = u + erf(-x)
  IF (a*u /= zero) a = ait + LOG(ABS(a/u)) / albeta
  v = erfc(x)
  b = u + v - one
  IF (b /= zero) b = ait + LOG(ABS(b)) / albeta
  w = erfcx(x)
  c = AINT(x*sixten) / sixten
  r = (x-c) * (x+c)
  c = (EXP(c*c)*EXP(r)*v-w) / w
  IF (c /= zero) c = MAX(zero,ait+LOG(ABS(c))/albeta)
  WRITE (iout,6400) x, a, b, c
  x = x + del
END DO
!-----------------------------------------------------------------
!  Next, test with special arguments
!-----------------------------------------------------------------
WRITE (iout,6500)
z = xmax
zz = erf(z)
WRITE (iout,6600) z, zz
z = zero
zz = erf(z)
WRITE (iout,6600) z, zz
zz = erfc(z)
WRITE (iout,6700) z, zz
z = -xmax
zz = erfc(z)
WRITE (iout,6700) z, zz
!-----------------------------------------------------------------
!  Test of error returns
!-----------------------------------------------------------------
WRITE (iout,6800)
w = xbig
z = w * (one-half*half)
WRITE (iout,6900) z
zz = erfc(z)
WRITE (iout,7400) zz
z = w * (one+ten*eps)
WRITE (iout,7000) z
zz = erfc(z)
WRITE (iout,7400) zz
w = xmax
IF (c1 < xmax*xmin) w = c1 / xmin
z = w * (one-one/sixten)
WRITE (iout,7100) z
zz = erfcx(z)
WRITE (iout,7500) zz
w = -SQRT(LOG(xmax/two))
z = w * (one-one/ten)
WRITE (iout,7200) z
zz = erfcx(z)
WRITE (iout,7500) zz
z = w * (one+ten*eps)
WRITE (iout,7300) z
zz = erfcx(z)
WRITE (iout,7500) zz
WRITE (iout,7600)
STOP
!-----------------------------------------------------------------
5000 FORMAT ('1Test of erf(x) vs double series expansion'//)
5100 FORMAT (///' Test of erfc(x) vs exp(x+1/4) SUM i^n erfc(x+1/2)'//)
5200 FORMAT ('1Test of exp(x*x) erfc(x) vs SUM i^n erfc(x+1/2)'//)
5300 FORMAT (i7,' Random arguments were tested from the interval (',f7  &
    .3,',',f7.3,')'//)
5400 FORMAT ('    ERF(X) was larger',i6,' times,')
5500 FORMAT ('   ERFC(X) was larger',i6,' times,')
5600 FORMAT ('  ERFCX(X) was larger',i6,' times,')
5700 FORMAT (t15, ' agreed', i6,' times, and'/ t11,'was smaller',i6, ' times.'//)
5800 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number'//)
5900 FORMAT (' The maximum relative error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =', e13.6)
6000 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
6100 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
6200 FORMAT ('1Special Tests'//)
6300 FORMAT (t8, 'Estimated loss of base', i3, 'significant digits in'//  &
    t4, 'X     Erf(x)+Erf(-x)   Erf(x)+Erfc(x)-1   ',  &
    'Erfcx(x)-exp(x*x)*erfc(x)'/)
6400 FORMAT (f7.3,3F16.2)
6500 FORMAT (//' Test of special arguments'//)
6600 FORMAT ('   ERF (',e13.6,') = ',e13.6//)
6700 FORMAT ('  ERFC (',e13.6,') = ',e13.6//)
6800 FORMAT (' Test of Error Returns'///)
6900 FORMAT (' ERFC will be called with the argument',e13.6,/  &
    ' This should not underflow'//)
7000 FORMAT (' ERFC will be called with the argument',e13.6,/  &
    ' This may underflow'//)
7100 FORMAT (' ERFCX will be called with the argument',e13.6,/  &
    ' This should not underflow'//)
7200 FORMAT (' ERFCX will be called with the argument',e13.6,/  &
    ' This should not overflow'//)
7300 FORMAT (' ERFCX will be called with the argument',e13.6,/  &
    ' This may overflow'//)
7400 FORMAT (' ERFC returned the value',e13.6///)
7500 FORMAT (' ERFCX returned the value',e13.6///)
7600 FORMAT (' This concludes the tests')
!---------- Last line of ERF test program ----------
END PROGRAM erftst
