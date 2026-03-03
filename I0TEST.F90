PROGRAM i0test
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------
! FORTRAN 90 program to test BESI0, BESEI0

!  Method:

!     Accuracy tests compare function values against values generated with
!     the multiplication formula for small arguments and values generated
!     from a Taylor's Series Expansion using Amos' Ratio Scheme for initial
!     values for large arguments.

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - An environmental inquiry program providing information on the
!         floating-point arithmetic system.  Note that the call to MACHAR can
!         be deleted provided the following five parameters are assigned
!         the values indicated

!              IBETA  - the radix of the floating-point system
!              IT     - the number of base-IBETA digits in the
!                       significant of a floating-point number
!              EPS    - the smallest positive floating-point
!                       number such that 1.0+EPS .NE. 1.0
!              XMIN   - the smallest non-vanishing normalized
!                       floating-point power of the radix
!              XMAX   - the largest finite floating-point number

!     REN(K) - a function subprogram returning random real
!              numbers uniformly distributed over (0,1)


!  Intrinsic functions required are:

!      ABS, DBLE, INT, LOG, MAX, REAL, SQRT

!  Reference: "Computation of Modified Bessel Functions and
!              Their Ratios," D. E. Amos, Math. of Comp.,
!              Volume 28, Number 24, January, 1974.

!             "Performance evaluation of programs for certain
!              Bessel functions", W. J. Cody and L. Stoltz,
!              ACM Trans. on Math. Software, Vol. 15, 1989, pp 41-48.

!             "Use of Taylor series to test accuracy of function programs,"
!              W. J. Cody and L. Stoltz, submitted for publication.

!  Latest modification: March 12, 1992

!  Authors:  W. J. Cody and L. Stoltz
!            Mathematics and Computer Science Division
!            Argonne National Laboratory
!            Argonne, IL 60439

!------------------------------------------------------------------

USE Bessel_I0
USE Toms715_Utilities

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER :: i, ibeta, iexp, ii, ind, irnd, it, j, j1, j2,  &
           k1, k2, k3, machep, maxexp, mb, mborg, minexp, mb2, n, negep, ngrd
REAL (dp)  :: a, ait, ak, akk, albeta, b, beta, const, d, del, delta, e,   &
              eps, epsneg, f, ovrchk, r6, r7, sum, temp, t1, t2, w, x, xa, &
              xb, xbad, xj1, xl, xlarge, xmax, xmb, xmin, xn, x1, y, z, zz
REAL (dp)  :: u(560), u2(560)
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, sixten = 16.0_dp,  &
                         hund = 100.0_dp, x99 = -999.0_dp,   &
                         xlam = 1.03125_dp, c = 0.9189385332_dp
REAL (dp), PARAMETER  :: ARR(8,6) = RESHAPE(  &
        (/ 0.0D0, 1.0D0, -1.0D0, 1.0D0, -2.0D0, 1.0D0, -3.0D0, 1.0D0,    &
          -999.0D0, -999.0D0, -999.0D0, 3.0D0, -12.0D0, 9.0D0, -51.0D0,  &
          18.0D0, -5040.0D0, 720.0D0, 0.0D0, -999.0D0, -999.0D0,         &
          60.0D0, -360.0D0, 345.0D0, -1320.0D0, 192.0D0, -120.0D0,       &
          24.0D0, 0.0D0, -999.0D0, -999.0D0, 2520.0D0, -96.0D0, 15.0D0,  &
          -33.0D0, 7.0D0, -6.0D0, 2.0D0, 0.0D0, -999.0D0, -4.0D0, 1.0D0, &
          -3.0D0, 1.0D0, -2.0D0, 1.0D0, -1.0D0, 1.0D0 /), (/ 8, 6 /) )
INTEGER, PARAMETER  :: iout = 6
!------------------------------------------------------------------

INTERFACE
  FUNCTION top(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION top

  FUNCTION bot(x) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION bot
END INTERFACE

!------------------------------------------------------------------
!  Determine machine parameters and set constants
!------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
xlarge = 1.0D4
beta = ibeta
ait = it
albeta = LOG(beta)
a = zero
b = two
const = c + LOG(xmax)
delta = xlam - one
f = (xlam-one) * (xlam+one) * half
!-----------------------------------------------------------------
!  Random argument accuracy tests
!-----------------------------------------------------------------
DO  j = 1, 4
!-------------------------------------------------------------------
!  Calculate the number of terms needed for convergence of the series by
!  using Newton's iteration on the asymptotic form of the multiplication
!  theorem
!-------------------------------------------------------------------
  xbad = b
  d = ait * albeta - c + one
  e = LOG(xbad*f) + one
  akk = one
  10 ak = akk
  z = d + e * ak - (ak+half) * LOG(ak+one)
  zz = e - (ak+half) / (ak+one) - LOG(ak+one)
  akk = ak - z / zz
  IF (ABS(ak-akk) > hund*eps*ak) GO TO 10
  mborg = INT(akk) + 1
  n = 2000
  xn = n
  k1 = 0
  k2 = 0
  k3 = 0
  x1 = zero
  r6 = zero
  r7 = zero
  del = (b-a) / xn
  xl = a
  DO  i = 1, n
    x = del * ren() + xl
!------------------------------------------------------------------
!   Carefully purify arguments
!------------------------------------------------------------------
    IF (j == 1) THEN
      y = x / xlam
    ELSE
      y = x - delta
    END IF
    temp = sixten * y
    t1 = temp + y
    t1 = temp + t1
    y = t1 - temp
    y = y - temp
    IF (j == 1) THEN
      x = y * xlam
    ELSE
      x = y + delta
    END IF
!------------------------------------------------------------------
!   Use Amos' Ratio Scheme
!------------------------------------------------------------------
    d = f * y
    mb = mborg + mborg
    mb2 = mb - 1
    xmb = mb2
    temp = (xmb+one+half) * (xmb+one+half)
    u2(mb) = y / (xmb+half+SQRT(temp+y*y))
!------------------------------------------------------------------
!   Generate ratios using recurrence
!------------------------------------------------------------------
    DO  ii = 2, mb
      ovrchk = xmb / (y*half)
      u2(mb2) = one / (ovrchk+u2(mb2+1))
      xmb = xmb - one
      mb2 = mb2 - 1
    END DO
    u(1) = besi0(y)
    IF (j == 1) THEN
!------------------------------------------------------------------
!   Accuracy test is based on the multiplication theorem
!------------------------------------------------------------------
      mb = mb - mborg
      DO  ii = 2, mb
        u(ii) = u(ii-1) * u2(ii-1)
      END DO
!------------------------------------------------------------------
!   Accurate Summation
!------------------------------------------------------------------
      mb = mb - 1
      xmb = mb
      sum = u(mb+1)
      ind = mb
      DO  ii = 2, mb
        sum = sum * d / xmb + u(ind)
        ind = ind - 1
        xmb = xmb - one
      END DO
      zz = sum * d + u(ind)
    ELSE
!------------------------------------------------------------------
!   Accuracy test is based on Taylor's Series Expansion
!------------------------------------------------------------------
      u(2) = u(1) * u2(1)
      mb = 8
      j1 = mb
      xj1 = j1+1
      iexp = 0
!------------------------------------------------------------------
!   Accurate Summation
!------------------------------------------------------------------
      DO  ii = 1, mb
        j2 = 1
        50 j2 = j2 + 1
        IF (arr(j1,j2) /= x99) GO TO 50
        j2 = j2 - 1
        t1 = arr(j1,j2)
        j2 = j2 - 1
!------------------------------------------------------------------
!   Group I0 terms in the derivative
!------------------------------------------------------------------
        IF (j2 /= 0) THEN
          60 t1 = t1 / (y*y) + arr(j1,j2)
          j2 = j2 - 1
          IF (j2 >= 1) GO TO 60
        END IF
        IF (iexp == 1) t1 = t1 / y
        j2 = 6
        70 j2 = j2 - 1
        IF (arr(ii,j2) /= x99) GO TO 70
        j2 = j2 + 1
        t2 = arr(ii,j2)
        j2 = j2 + 1
        IF (iexp == 0) THEN
          iexp = 1
        ELSE
          iexp = 0
        END IF
!------------------------------------------------------------------
!   Group I1 terms in the derivative
!------------------------------------------------------------------
        IF (j2 /= 7) THEN
          80 t2 = t2 / (y*y) + arr(ii,j2)
          j2 = j2 + 1
          IF (j2 <= 6) GO TO 80
        END IF
        IF (iexp == 1) t2 = t2 / y
        IF (j1 == 8) THEN
          sum = u(1) * t1 + u(2) * t2
        ELSE
          sum = sum * (delta/xj1) + (u(1)*t1+u(2)*t2)
        END IF
        j1 = j1 - 1
        xj1 = j1+1
      END DO
      zz = sum * delta + u(1)
    END IF
    z = besi0(x)
!------------------------------------------------------------------
!   Accumulate Results
!------------------------------------------------------------------
    w = (z-zz) / z
    IF (w > zero) THEN
      k1 = k1 + 1
    ELSE IF (w < zero) THEN
      k3 = k3 + 1
    ELSE
      k2 = k2 + 1
    END IF
    w = ABS(w)
    IF (w > r6) THEN
      r6 = w
      x1 = x
    END IF
    r7 = r7 + w * w
    xl = xl + del
  END DO
!------------------------------------------------------------------
!   Gather and print statistics for test
!------------------------------------------------------------------
  n = k1 + k2 + k3
  r7 = SQRT(r7/xn)
  IF (j == 1) THEN
    WRITE (iout,5000)
  ELSE
    WRITE (iout,5100)
  END IF
  WRITE (iout,5200) n, a, b
  WRITE (iout,5300) k1, k2, k3
  WRITE (iout,5400) it, ibeta
  IF (r6 /= zero) THEN
    w = LOG(r6) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5500) r6, ibeta, w, x1
  w = MAX(ait+w,zero)
  WRITE (iout,5600) ibeta, w
  IF (r7 /= zero) THEN
    w = LOG(r7) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5700) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,5600) ibeta, w
!------------------------------------------------------------------
!   Initialize for next test
!------------------------------------------------------------------
  a = b
  b = b + b
  IF (j == 1) b = b + b - half
END DO
!-----------------------------------------------------------------
!   Test of error returns

!   Special tests
!-----------------------------------------------------------------
WRITE (iout,5800)
WRITE (iout,5900)
y = besi0(xmin)
WRITE (iout,6000) y
y = besi0(zero)
WRITE (iout,6100) 0, y
x = -one * ren()
y = besi0(x)
WRITE (iout,6200) x, y
x = -x
y = besi0(x)
WRITE (iout,6200) x, y
y = besei0(xmax)
WRITE (iout,6300) y
!-----------------------------------------------------------------
!   Determine largest safe argument for unscaled functions
!-----------------------------------------------------------------
WRITE (iout,6400)
xa = LOG(xmax)
120 xb = xa - (top(xa)-const) / bot(xa)
IF (ABS(xb-xa)/xb <= eps) THEN
  GO TO 130
ELSE
  xa = xb
  GO TO 120
END IF
130 xlarge = xb / xlam
y = besi0(xlarge)
WRITE (iout,6200) xlarge, y
xlarge = xb * xlam
y = besi0(xlarge)
WRITE (iout,6200) xlarge, y
WRITE (iout,6500)
STOP
!-----------------------------------------------------------------
5000 FORMAT ('1Test of I0(X) vs Multiplication Theorem'//)
5100 FORMAT ('1Test of I0(X) vs Taylor series'//)
5200 FORMAT (i7, ' Random arguments were tested from the interval (',  &
    f5.2, ',', f5.2, ')'//)
5300 FORMAT (' I0(X) was larger', i6, ' times,'/ t11, ' agreed', i6,  &
    ' times, and'/ t7, 'was smaller', i6, ' times.'//)
5400 FORMAT (' There are',i4,' base', i4,  &
    ' significant digits in a floating-point number'//)
5500 FORMAT (' The maximum relative error of', e15.4, ' = ', i4, ' **', f7.2  &
    /t5, 'occurred for X =', e13.6)
5600 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5700 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5800 FORMAT ('1Special Tests'//)
5900 FORMAT (' Test with extreme arguments'/)
6000 FORMAT (' I0(XMIN) = ', e24.17/)
6100 FORMAT (' I0(',i1,') = ', e24.17/)
6200 FORMAT (' I0(',e24.17,' ) = ', e24.17/)
6300 FORMAT (' E**-X * I0(XMAX) = ', e24.17/)
6400 FORMAT (' Tests near the largest argument for unscaled functions'/ )
6500 FORMAT (' This concludes the tests.')
!---------- Last line of BESI0 test program ----------
END PROGRAM i0test



FUNCTION top(x) RESULT(fn_val)
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

REAL (dp), PARAMETER   :: one = 1.0_dp, eight = 8.0_dp, one28 = 128.0_dp, &
                          half = 0.5_dp, xnine = 9.0_dp

fn_val = x - half * LOG(x) + LOG(one+(one/eight-xnine/one28/x)/x)
RETURN
END FUNCTION top



FUNCTION bot(x) RESULT(fn_val)
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

REAL (dp), PARAMETER   :: one = 1.0_dp, one28 = 128.0_dp, sixten = 16.0_dp, &
                          xnine = 9.0_dp, half = 0.5_dp, ateten = 18.0_dp

fn_val = -(sixten*x+ateten) / (((one28*x+sixten)*x+xnine)*x) + one - half / x
RETURN
END FUNCTION bot

