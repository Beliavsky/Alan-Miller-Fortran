PROGRAM eitest
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------
! FORTRAN 90 program to test EI, EONE, and EXPEI.

!  Method:

!     Accuracy test compare function values against local Taylor's series
!     expansions.  Derivatives for Ei(x) are generated from the recurrence
!     relation using a technique due to Gautschi (see references).
!     Special argument tests are run with the related functions E1(x) and
!     exp(-x)Ei(x).

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - An environmental inquiry program providing information on
!              the floating-point arithmetic system.  Note that the call to
!              MACHAR can be deleted provided the following five parameters
!              are assigned the values indicated

!              IBETA  - The radix of the floating-point system
!              IT     - The number of base-ibeta digits in the
!                       significant of a floating-point number
!              XMAX   - The largest finite floating-point number

!     REN(K) - A function subprogram returning random real
!              numbers uniformly distributed over (0,1)


!  Intrinsic functions required are:

!      ABS, AINT, DBLE, LOG, MAX, REAL, SQRT

!  References: "The use of Taylor series to test accuracy of
!               function programs", Cody, W. J., and Stoltz, L.,
!               submitted for publication.

!              "Recursive computation of certain derivatives -
!               A study of error propagation", Gautschi, W., and
!               Klein, B. J., Comm. ACM 13 (1970), 7-9.

!              "Remark on Algorithm 282", Gautschi, W., and Klein,
!               B. J., Comm. ACM 13 (1970), 53-54.

!  Latest modification: March 9, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!------------------------------------------------------------------

USE Toms715_Utilities
USE Bessel_EI
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, irnd, it, j, k, k1, k2, k3,  &
              machep, maxexp, minexp, n, negep, ngrd, n1
REAL (dp)  :: a, ait, albeta, b, beta, c1, del, dx, en, eps, epsneg,   &
              r6, r7, sum, two, u, v, w, x, xbig, xc, xden, xl,  &
              xlge, xmax, xmin, xn, xnp1, xnum, y, z
REAL (dp)  :: d(0:25)
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, fourth = 0.25_dp, one = 1.0_dp,  &
                         four = 4.0_dp, six = 6.0_dp, ten = 10.0_dp,  &
                         x0 = 0.3725_dp, x99 = -999.0_dp, p0625 = 0.0625_dp, &
                         fiv12 = 512.0_dp, rem = -7.424779065800051695596D-5
INTEGER, PARAMETER    :: iout = 6
!------------------------------------------------------------------
!  Determine machine parameters and set constants
!------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
dx = -p0625
a = fourth + dx
b = x0 + dx
n = 2000
n1 = 25
xn = n
!-----------------------------------------------------------------
!  Random argument accuracy tests
!-----------------------------------------------------------------
DO  j = 1, 8
  k1 = 0
  k3 = 0
  xc = zero
  r6 = zero
  r7 = zero
  del = (b-a) / xn
  xl = a
  DO  i = 1, n
    y = del * ren() + xl
    x = y - dx
    y = x + dx
!-----------------------------------------------------------------
!  Test Ei against series expansion
!-----------------------------------------------------------------
    v = ei(x)
    z = ei(y)
    sum = zero
    u = x
    CALL dsubn(u, n1, xmax, d)
    en = n1 + one
    sum = d(n1) * dx / en
    DO  k = n1, 1, -1
      en = en - one
      sum = (sum+d(k-1)) * dx / en
    END DO
    u = v + sum
!--------------------------------------------------------------------
!  Accumulate results
!--------------------------------------------------------------------
    w = z - u
    w = w / z
    IF (w > zero) THEN
      k1 = k1 + 1
    ELSE IF (w < zero) THEN
      k3 = k3 + 1
    END IF
    w = ABS(w)
    IF (w > r6) THEN
      r6 = w
      xc = y
    END IF
    r7 = r7 + w * w
    xl = xl + del
  END DO
!------------------------------------------------------------------
!  Gather and print statistics for test
!------------------------------------------------------------------
  k2 = n - k3 - k1
  r7 = SQRT(r7/xn)
  WRITE (iout,5000)
  WRITE (iout,5100) n, a, b
  WRITE (iout,5200) k1
  WRITE (iout,5300) k2, k3
  WRITE (iout,5400) it, ibeta
  IF (r6 /= zero) THEN
    w = LOG(ABS(r6)) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5500) r6, ibeta, w, xc
  w = MAX(ait+w,zero)
  WRITE (iout,5600) ibeta, w
  IF (r7 /= zero) THEN
    w = LOG(ABS(r7)) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5700) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,5600) ibeta, w
!------------------------------------------------------------------
!  Initialize for next test
!------------------------------------------------------------------
  IF (j == 1) THEN
    dx = -dx
    a = x0 + dx
    b = six
  ELSE IF (j <= 4) THEN
    a = b
    b = b + b
  ELSE IF (j == 5) THEN
    a = -fourth
    b = -one
  ELSE IF (j == 6) THEN
    a = b
    b = -four
  ELSE
    a = b
    b = -ten
  END IF
END DO
!-----------------------------------------------------------------
!  Special tests.  First, check accuracy near the zero of Ei(x)
!-----------------------------------------------------------------
WRITE (iout,5800)
x = (four-one) / (four+four)
y = ei(x)
WRITE (iout,5900) x, y
z = ((y-(four+one)/(fiv12))-rem) / y
IF (z /= zero) THEN
  w = LOG(ABS(z)) / albeta
ELSE
  w = x99
END IF
WRITE (iout,6000) z, ibeta, w
w = MAX(ait+w,zero)
WRITE (iout,5600) ibeta, w
!-----------------------------------------------------------------
!  Check near XBIG, the largest argument acceptable to EONE, i.e.,
!    the negative of the smallest argument acceptable to EI.
!    Determine XBIG with Newton iteration on the equation
!                  EONE(x) = XMIN.
!---------------------------------------------------------------------
WRITE (iout,6100)
two = one + one
v = SQRT(eps)
c1 = minexp * LOG(beta)
xn = -c1
40 xnum = -xn - LOG(xn) + LOG(one+one/xn) - c1
xden = -(xn*xn+xn+xn+two) / (xn*(xn+one))
xnp1 = xn - xnum / xden
w = (xn-xnp1) / xnp1
IF (ABS(w) > v) THEN
  xn = xnp1
  GO TO 40
END IF
xbig = xnp1
x = AINT(ten*xbig) / ten
WRITE (iout,6200) x
y = eone(x)
WRITE (iout,6700) y
x = xbig * (one+v)
WRITE (iout,6300) x
y = eone(x)
WRITE (iout,6700) y
!---------------------------------------------------------------------
!  Check near XMAX, the largest argument acceptable to EI.  Determine
!    XLGE with Newton iteration on the equation
!                  EI(x) = XMAX.
!---------------------------------------------------------------------
c1 = maxexp * LOG(beta)
xn = c1
50 xnum = xn - LOG(xn) + LOG(one+one/xn) - c1
xden = (xn*xn-two) / (xn*(xn+one))
xnp1 = xn - xnum / xden
w = (xn-xnp1) / xnp1
IF (ABS(w) > v) THEN
  xn = xnp1
  GO TO 50
END IF
xlge = xnp1
x = AINT(ten*xlge) / ten
WRITE (iout,6400) x
y = ei(x)
WRITE (iout,6800) y
x = xlge * (one+v)
WRITE (iout,6500) x
y = ei(x)
WRITE (iout,6800) y
!---------------------------------------------------------------------
!  Check with XHUGE, the largest acceptable argument for EXPEI
!---------------------------------------------------------------------
IF (xmin*xmax <= one) THEN
  x = xmax
ELSE
  x = one / xmin
END IF
WRITE (iout,6600) x
y = expei(x)
WRITE (iout,6900) y
x = zero
WRITE (iout,6500) x
y = ei(x)
WRITE (iout,6800) y
WRITE (iout,7000)
STOP
!-----------------------------------------------------------------
5000 FORMAT ('1Test of Ei(x) vs series expansion'//)
5100 FORMAT (i7, ' Random arguments were tested from the interval (', f7.3,  &
             ',', f7.3, ')'//)
5200 FORMAT ('     EI(X) was larger',i6,' times,')
5300 FORMAT (t15, ' agreed', i6, ' times, and'/   &
             t11, 'was smaller', i6, ' times.'//)
5400 FORMAT (' There are',i4,' base', i4,  &
             ' significant digits in a floating-point number'//)
5500 FORMAT (' The maximum relative error of', e15.4, ' = ', i4, ' **', f7.2  &
    /t5, 'occurred for X =', e13.6)
5600 FORMAT (' The estimated loss of base', i4, ' significant digits is',  &
    f7.2//)
5700 FORMAT (' The root mean square relative error was', e15.4, ' = ', i4,  &
    ' **', f7.2)
5800 FORMAT (//' Test of special arguments'//)
5900 FORMAT ('   EI (', e13.6, ') = ', e13.6//)
6000 FORMAT (' The relative error is', e15.4, ' = ', i4, ' **', f7.2/)
6100 FORMAT (' Test of Error Returns'///)
6200 FORMAT (' EONE will be called with the argument', e13.6/  &
    ' This should not underflow'//)
6300 FORMAT (' EONE will be called with the argument', e13.6/  &
    ' This should underflow'//)
6400 FORMAT (' EI will be called with the argument', e13.6/  &
    ' This should not overflow'//)
6500 FORMAT (' EI will be called with the argument', e13.6/  &
    ' This should overflow'//)
6600 FORMAT (' EXPEI will be called with the argument', e13.6/  &
    ' This should not underflow'//)
6700 FORMAT (' EONE returned the value', e13.6///)
6800 FORMAT (' EI returned the value', e13.6///)
6900 FORMAT (' EXPEI returned the value', e13.6///)
7000 FORMAT (' This concludes the tests')
!---------- Last line of EI test program ----------
END PROGRAM eitest
