PROGRAM anrtst
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------
! FORTRAN 77 program to test ANORM, the normal probability integral.

!  Method:

!     Accuracy test compares function values against local Taylor's
!     series expansions.  Derivatives for anorm(x) are expressed as
!     repeated integrals of erfc(-x/sqrt(2)).  These are generated
!     from the recurrence relation using a technique due to Gautschi
!     (see references).

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - An environmental inquiry program providing information on
!              information on the floating-point arithmetic system.
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
!               TOMS 16, 1990, pp. 29-37.

!              "Evaluation of the repeated integrals of the coerror
!               function", W. Gautschi, TOMS 3, 1977, pp. 240-252.

!  Latest modification: March 15, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!------------------------------------------------------------------
USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, irnd, it, j, k, k1, k2, k3,  &
              machep, maxexp, minexp, n, negep, ngrd, n0, n1
REAL (dp)  :: a, ait, albeta, b, beta, c, c2, del, eps, epscon,  &
              epsneg, ff, f0, r, r6, r7, sc,  &
              u, v, w, x, xbig, xc, xl, xmax, xmin, xn, xn1, y, z
REAL (dp)  :: r1(500)
!------------------------------------------------------------------
!   C1 = 1/sqrt(pi)
!------------------------------------------------------------------
REAL (dp), PARAMETER :: half = 0.5_dp, zero = 0.0_dp, one = 1.0_dp,  &
                        two = 2.0_dp, ten = 10.0_dp, sixten = 16.0_dp,  &
                        thresh = 0.66291_dp, x99 = -999.0_dp,  &
                        c1 = 5.6418958354775628695D-1, root32 = -5.65685_dp, &
                        sqrtwo = 7.0710678118654752440D-1
INTEGER, PARAMETER  :: iout = 6

INTERFACE
  FUNCTION ANORM(ARG) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: arg
    REAL (dp)              :: fn_val
  END FUNCTION ANORM
END INTERFACE
!------------------------------------------------------------------
!  Determine machine parameters and set constants
!------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
c = ABS(ait*albeta) + ten
a = -thresh
b = thresh
n = 2000
xn = n
n0 = (ibeta/2) * (it+5) / 6 + 4
!-----------------------------------------------------------------
!  Determine largest argument for ANORM test by Newton iteration
!-----------------------------------------------------------------
c2 = LOG(xmin/c1/sqrtwo)
xbig = SQRT(-c2)
10 x = xbig
f0 = x * x
ff = f0 * half + one / f0 + LOG(x) + c2
f0 = x + one / x - two / (x*f0-x)
xbig = x - ff / f0
IF (ABS(x-xbig)/x > ten*eps) GO TO 10
!-----------------------------------------------------------------
!  Random argument accuracy tests
!-----------------------------------------------------------------
DO  j = 1, 3
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
!  Test anorm against double series expansion
!-----------------------------------------------------------------
      f0 = n0
      ff = f0 + f0 + one
      w = x * x
      z = w * half
      u = zero
      v = zero
      DO  k = 1, n0
        u = -z / f0 * (one+u)
        v = w / ff * (one+v)
        f0 = f0 - one
        ff = ff - two
      END DO
      v = half + sqrtwo * c1 * x * (((u*v+(u+v))+half)+half)
      u = anorm(x)
    ELSE
!-----------------------------------------------------------------
!  Test anorm against expansion in repeated
!   integrals of the coerror function.
!-----------------------------------------------------------------
      y = x - half
      x = y + half
      u = anorm(x)
      z = -y / SQRT(two)
      r = zero
      IF (z <= one) THEN
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
      f0 = anorm(y) * EXP((-y-half*half)*half)
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
      v = sqrtwo * r1(k)
      DO  n1 = 1, k - 1
        v = (v+r1(k-n1)) * sqrtwo
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
  ELSE
    WRITE (iout,5100)
  END IF
  WRITE (iout,5200) n, a, b
  WRITE (iout,5300) k1
  WRITE (iout,5400) k2, k3
  WRITE (iout,5500) it, ibeta
  IF (r6 /= zero) THEN
    w = LOG(ABS(r6)) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5600) r6, ibeta, w, xc
  w = MAX(ait+w,zero)
  WRITE (iout,5700) ibeta, w
  IF (r7 /= zero) THEN
    w = LOG(ABS(r7)) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5800) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,5700) ibeta, w
!------------------------------------------------------------------
!  Initialize for next test
!------------------------------------------------------------------
  IF (j == 1) THEN
    b = a
    a = root32
  ELSE
    b = a
    a = -AINT(xbig*sixten) / sixten + half
  END IF
END DO
!-----------------------------------------------------------------
!  Special tests
!  First check values for positive arguments.
!-----------------------------------------------------------------
WRITE (iout,5900)
WRITE (iout,6000)
del = ten
DO  i = 1, 10
  x = ren() * del
  u = anorm(-x)
  a = u + anorm(x)
  y = (a-half) - half
  WRITE (iout,6100) x, u, y
END DO
!-----------------------------------------------------------------
!  Test with special arguments
!-----------------------------------------------------------------
WRITE (iout,6200)
z = xmax
y = anorm(z)
WRITE (iout,6300) z, y
z = zero
y = anorm(z)
WRITE (iout,6300) z, y
z = -xmax
y = anorm(z)
WRITE (iout,6300) z, y
!-----------------------------------------------------------------
!  Test of error returns
!-----------------------------------------------------------------
WRITE (iout,6400)
w = xbig
z = -w * (one-half*half)
WRITE (iout,6500) z
y = anorm(z)
WRITE (iout,6700) y
z = -w * (one+ten*eps)
WRITE (iout,6600) z
y = anorm(z)
WRITE (iout,6700) y
WRITE (iout,6800)
STOP
!-----------------------------------------------------------------
5000 FORMAT ('1Test of anorm(x) vs double series expansion'//)
5100 FORMAT (///' Test of anorm(x) vs Taylor series about x-1/2'//)
5200 FORMAT (i7,' Random arguments were tested from the interval (',f7.3,  &
    ',', f7.3, ')'//)
5300 FORMAT ('  ANORM(X) was larger',i6,' times,')
5400 FORMAT (t15, ' agreed',i6,' times, and'/ t11, 'was smaller',i6, ' times.'//)
5500 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number'//)
5600 FORMAT (' The maximum relative error of', e15.4, ' = ', i4, ' **',f7.2  &
    / t5, 'occurred for X =', e13.6)
5700 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5800 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5900 FORMAT ('1Special Tests'//)
6000 FORMAT (t8, 'Check of identity anorm(X) + anorm(-X) = 1.0'// t10, 'X', &
    t22, 'ANORM(-x)   ANORM(x)+ANORM(-x)-1'/)
6100 FORMAT (3(t4, e13.6)/)
6200 FORMAT (//' Test of special arguments'//)
6300 FORMAT (' ANORM (',e13.6,') = ',e13.6//)
6400 FORMAT (' Test of Error Returns'///)
6500 FORMAT (' ANORM will be called with the argument ',e13.6,/  &
    ' The result should not underflow'//)
6600 FORMAT (' ANORM will be called with the argument ',e13.6,/  &
    ' The result may underflow'//)
6700 FORMAT (' ANORM returned the value',e13.6///)
6800 FORMAT (' This concludes the tests')
!---------- Last line of ANORM test program ----------
END PROGRAM anrtst
