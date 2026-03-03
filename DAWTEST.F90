PROGRAM dawtst
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------
! FORTRAN 90 program to test DAW

!  Method:

!     Accuracy test compare function values against a local
!     Taylor's series expansion.  Derivatives are generated
!     from the recurrence relation.

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - An environmental inquiry program providing
!         information on the floating-point arithmetic
!         system.  Note that the call to MACHAR can
!         be deleted provided the following five
!         parameters are assigned the values indicated

!              IBETA  - the radix of the floating-point system
!              IT     - the number of base-IBETA digits in the
!                       significant of a floating-point number
!              XMIN   - the smallest positive floating-point number
!              XMAX   - the largest floating-point number

!     REN(K) - a function subprogram returning random real
!              numbers uniformly distributed over (0,1)


!  Intrinsic functions required are:

!         ABS, DBLE, LOG, MAX, REAL, SQRT

!  Reference: "The use of Taylor series to test accuracy of
!              function programs", W. J. Cody and L. Stoltz,
!              submitted for publication.

!  Latest modification: March 9, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!------------------------------------------------------------------

USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, ii, irnd, it, j, k1, k2, k3,  &
              machep, maxexp, minexp, n, ndum, negep, ngrd
REAL (dp)  :: a, ait, albeta, b, beta, del, eps, epsneg, r6, r7, t1, w,  &
              x, xbig, xkay, xl, xmax, xmin, xn, x1, y, z, zz
REAL (dp)  :: p(0:14)
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, forten = 14.0_dp, sixten = 16.0_dp,  &
                         x99 = -999.0_dp, delta = 0.0625_dp
INTEGER, PARAMETER    :: iout = 6

INTERFACE
  FUNCTION DAW(XX) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: xx
    REAL (dp)              :: fn_val
  END FUNCTION DAW
END INTERFACE

!------------------------------------------------------------------
!  Determine machine parameters and set constants
!------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
ait = it
albeta = LOG(beta)
a = delta
b = one
!-----------------------------------------------------------------
!  Random argument accuracy tests based on local Taylor expansion.
!-----------------------------------------------------------------
DO  j = 1, 4
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
!  Purify arguments
!------------------------------------------------------------------
    y = x - delta
    w = sixten * y
    t1 = w + y
    y = t1 - w
    x = y + delta
!------------------------------------------------------------------
!  Use Taylor's Series Expansion
!------------------------------------------------------------------
    p(0) = daw(y)
    z = y + y
    p(1) = one - z * p(0)
    xkay = two
    DO  ii = 2, 14
      p(ii) = -(z*p(ii-1)+xkay*p(ii-2))
      xkay = xkay + two
    END DO
    zz = p(14)
    xkay = forten
    DO  ii = 1, 14
      zz = zz * delta / xkay + p(14-ii)
      xkay = xkay - one
    END DO
    z = daw(x)
!------------------------------------------------------------------
!  Accumulate Results
!------------------------------------------------------------------
    w = (z-zz) / z
    IF (w > zero) THEN
      k1 = k1 + 1
    ELSE IF (w < zero) THEN
      k3 = k3 + 1
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
!  Gather and print statistics for test
!------------------------------------------------------------------
  k2 = n - k1 - k3
  r7 = SQRT(r7/xn)
  WRITE (iout,5000)
  WRITE (iout,5100) n, a, b
  WRITE (iout,5200) k1, k2, k3
  WRITE (iout,5300) it, ibeta
  w = x99
  IF (r6 /= zero) w = LOG(r6) / albeta
  WRITE (iout,5400) r6, ibeta, w, x1
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
  w = x99
  IF (r7 /= zero) w = LOG(r7) / albeta
  WRITE (iout,5600) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
!------------------------------------------------------------------
!  Initialize for next test
!------------------------------------------------------------------
  a = b
  b = b + b
  IF (j == 1) b = b + half
END DO
!-----------------------------------------------------------------
!  Special tests.  First check values for negative arguments.
!-----------------------------------------------------------------
WRITE (iout,5700)
WRITE (iout,5800) ibeta
DO  i = 1, 10
  x = ren() * (two+two)
  b = daw(x)
  a = b + daw(-x)
  IF (a*b /= zero) a = ait + LOG(ABS(a/b)) / albeta
  WRITE (iout,5900) x, a
  x = x + del
END DO
!-----------------------------------------------------------------
!  Next, test with special arguments
!-----------------------------------------------------------------
WRITE (iout,6000)
z = xmin
zz = daw(z)
WRITE (iout,6100) zz
!-----------------------------------------------------------------
!  Test of error return for arguments > xmax.  First, determine
!    xmax
!-----------------------------------------------------------------
IF (half < xmin*xmax) THEN
  xbig = half / xmin
ELSE
  xbig = xmax
END IF
WRITE (iout,6200)
z = xbig * (one-delta*delta)
WRITE (iout,6300) z
zz = daw(z)
WRITE (iout,6500) zz
z = xbig
WRITE (iout,6400) z
zz = daw(z)
WRITE (iout,6500) zz
w = one + delta * delta
IF (w < xmax/xbig) THEN
  z = xbig * w
  WRITE (iout,6400) z
  zz = daw(z)
  WRITE (iout,6500) zz
END IF
WRITE (iout,6600)
STOP
!-----------------------------------------------------------------
5000 FORMAT ('1Test of Dawson''S INTEGRAL VS TAYLOR EXPANSION'//)
5100 FORMAT (i7, ' Random arguments were tested from the interval (',  &
    f5.2, ',', f5.2, ')'//)
5200 FORMAT ('  F(X) was larger', i6, ' times, '/ t11, ' agreed', i6,  &
    ' times, and'/ t7, 'was smaller', i6, ' times.'//)
5300 FORMAT (' There are', i4, ' base', i4,  &
    ' significant digits in a floating-point number'//)
5400 FORMAT (' The maximum relative error of', e15.4, ' = ', i4, ' **', f7.2  &
    /t5, 'occurred for X =', e13.6)
5500 FORMAT (' The estimated loss of base', i4, ' significant digits is',  &
    f7.2//)
5600 FORMAT (' The root mean square relative error was', e15.4, ' = ', i4,  &
    ' **', f7.2)
5700 FORMAT ('1Special Tests'//)
5800 FORMAT (t8, 'Estimated loss of base', i3, ' significant digits in'//  &
    t9, 'X', t20, 'F(x)+F(-x)'/)
5900 FORMAT (t4, F7.3, f16.2)
6000 FORMAT (//' Test of special arguments'//)
6100 FORMAT ('  F(XMIN) = ', e24.17/)
6200 FORMAT (' Test of Error Returns'///)
6300 FORMAT (' DAW will be called with the argument', e13.6/  &
    ' This should not underflow'//)
6400 FORMAT (' DAW will be called with the argument', e13.6/  &
    ' This may underflow'//)
6500 FORMAT (' DAW returned the value', e13.6///)
6600 FORMAT (' This concludes the tests')
!---------- Last line of DAW test program ----------
END PROGRAM dawtst
