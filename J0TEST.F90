PROGRAM j0test
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!--------------------------------------------------------------------
!  Fortran 90 program to test BESJ0

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - an environmental inquiry program providing information on the
!              floating-point arithmetic system.  Note that the call to
!              MACHAR can be deleted provided the following six parameters
!              are assigned the values indicated

!                 IBETA  - the radix of the floating-point system
!                 IT     - the number of base-IBETA digits in the
!                          significand of a floating-point number
!                 MAXEXP - the smallest integer such that
!                          FLOAT(IBETA)**MAXEXP causes overflow
!                 EPS    - the smallest positive floating-point
!                          number such that 1.0+EPS .NE. 1.0
!                 XMIN   - the smallest positive normalized
!                          floating-point power of the radix
!                 XMAX   - the largest finite floating-point number

!     REN(K) - a function subprogram returning random real
!              numbers uniformly distributed over (0,1)


!  Intrinsic functions required are:

!         ABS, DBLE, FLOAT, LOG, MAX, SIGN, SQRT

!  Reference: "The use of Taylor series to test accuracy of function
!              programs", W. J. Cody and L. Stoltz, submitted for publication.

!  Latest modification: March 13, 1992

!  Authors: W. J. Cody
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439, USA

!--------------------------------------------------------------------

USE Bessel_JY0
USE Bessel_JY1
USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, ii, iii, irnd, it, j, jj, k1,  &
              k2, k3, machep, maxexp, minexp, n, negep, ngrd
REAL (dp)  :: a, ait, albeta, b, beta, bj0, bj1,  &
              cj0, cj1, d, del, eps, epsneg, r6, r7, sum, t, term, w, x,  &
              xl, xmax, xmin, xm, xn, x1, y, yinv, ysq, z, zz
!--------------------------------------------------------------------
!  Mathematical constants
!--------------------------------------------------------------------
INTEGER, PARAMETER    :: iout = 6
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, four = 4.0_dp,  &
                         delta = 0.0625_dp, eight = 8.0_dp, twenty = 20.0_dp, &
                         all9 = -999.0_dp, two56 = 256.0_dp,  &
                         sixten = 16.0_dp, elev = 11.0_dp
!--------------------------------------------------------------------
!  Coefficients for Taylor expansion
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: BJ0P(6,10) = RESHAPE(  &
           (/ 0.0_dp, 1.8144D6, -2.3688D5, 1.0845D4, -2.7D2, 5.0_dp,  &
             -1.8144D5, 2.394D4, -1.125D3, 30.0_dp, -1.0_dp, 0.0_dp,  &
              0.0_dp, 2.016D4, -2.7D3, 1.32D2, -4.0_dp, 0.0_dp,  &
             -2.52D3, 3.45D2, -18.0_dp, 1.0_dp, 0.0_dp, 0.0_dp,  &
              0.0_dp, 3.6D2, -51.0_dp, 3.0_dp, 0.0_dp, 0.0_dp,  &
            -60.0_dp, 9.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
              0.0_dp, 12.0_dp, -2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
             -3.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
              0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
             -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /), (/ 6, 10 /) )
REAL (dp), PARAMETER  :: BJ1P(6,10) = RESHAPE(  &
           (/ -3.6288D6, 9.2736D5, -6.201D4, 1.965D3, -40.0_dp, 1.0_dp,  &
               3.6288D5, -9.324D4, 6.345D3, -2.1D2, 5.0_dp, 0.0_dp,  &
              -4.032D4, 1.044D4, -7.29D2, 26.0_dp, -1.0_dp, 0.0_dp,  &
               5.04D3, -1.32D3, 96.0_dp, -4.0_dp, 0.0_dp, 0.0_dp,   &
              -7.2D2, 1.92D2, -15.0_dp, 1.0_dp, 0.0_dp, 0.0_dp,   &
               1.2D2, -33.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
             -24.0_dp, 7.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
               6.0_dp, -2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
              -2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,   &
               1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /), (/6, 10 /) )
!--------------------------------------------------------------------
!  Zeroes of J0
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: XI(2) = (/ 616.0_dp, 1413.0_dp /),  &
                         YX(2) = (/ -7.3927648221700192757D-4,  &
                                    -1.8608651797573879013D-4 /)
!--------------------------------------------------------------------

CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
b = zero
!--------------------------------------------------------------------
!  Random argument accuracy tests (based on a Taylor expansion)
!--------------------------------------------------------------------
DO  j = 1, 3
  k1 = 0
  k2 = 0
  k3 = 0
  x1 = zero
  r6 = zero
  r7 = zero
  n = 2000
  a = b
  IF (j == 1) THEN
    b = four
  ELSE IF (j == 2) THEN
    b = eight
  ELSE
    b = twenty
    n = 2000
  END IF
  xn = n
  del = (b-a) / xn
  xl = a
  DO  i = 1, n
    x = del * ren() + xl
!--------------------------------------------------------------------
!  Carefully purify arguments and evaluate identities
!--------------------------------------------------------------------
    y = x - delta
    w = sixten * y
    y = (w+y) - w
    x = y + delta
    sum = zero
    term = zero
    bj0 = besj0(y)
    z = besj0(x)
    d = delta
    IF (ABS(z) < ABS(bj0)) THEN
      cj0 = x
      x = y
      y = cj0
      cj0 = bj0
      bj0 = z
      z = cj0
      d = -d
    END IF
    bj1 = besj1(y)
    yinv = one / y
    ysq = one / (y*y)
    xm = elev
!--------------------------------------------------------------------
!  Evaluate (12-II)th derivative at Y.
!--------------------------------------------------------------------
    DO  ii = 1, 10
      cj0 = bj0p(1,ii)
      cj1 = bj1p(1,ii)
      jj = (11-ii) / 2 + 1
      DO  iii = 2, jj
        cj0 = cj0 * ysq + bj0p(iii,ii)
        cj1 = cj1 * ysq + bj1p(iii,ii)
      END DO
      IF ((ii/2)*2 /= ii) THEN
        cj0 = cj0 * yinv
      ELSE
        cj1 = cj1 * yinv
      END IF
      term = cj0 * bj0 + cj1 * bj1
      sum = (sum+term) * d / xm
      xm = xm - one
    END DO
    sum = (sum-bj1) * d + bj0
    zz = sum
!--------------------------------------------------------------------
!  Accumulate results
!--------------------------------------------------------------------
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
!--------------------------------------------------------------------
!  Process and output statistics
!--------------------------------------------------------------------
  n = k1 + k2 + k3
  r7 = SQRT(r7/xn)
  WRITE (iout,5000)
  WRITE (iout,5100) n, a, b
  WRITE (iout,5200) k1, k2, k3
  WRITE (iout,5300) it, ibeta
  IF (r6 /= zero) THEN
    w = LOG(ABS(r6)) / albeta
  ELSE
    w = all9
  END IF
  WRITE (iout,5400) r6, ibeta, w, x1
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
  IF (r7 /= zero) THEN
    w = LOG(ABS(r7)) / albeta
  ELSE
    w = all9
  END IF
  WRITE (iout,5600) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
END DO
!--------------------------------------------------------------------
!  Special tests
!--------------------------------------------------------------------
WRITE (iout,5700)
WRITE (iout,5800) ibeta
DO  i = 1, 2
  x = xi(i) / two56
  y = besj0(x)
  t = (y-yx(i)) / yx(i)
  IF (t /= zero) THEN
    w = LOG(ABS(t)) / albeta
  ELSE
    w = all9
  END IF
  w = MAX(ait+w,zero)
  WRITE (iout,5900) x, y, w
END DO
!--------------------------------------------------------------------
!  Test of error returns
!--------------------------------------------------------------------
WRITE (iout,6000)
x = xmax
WRITE (iout,6100) x
y = besj0(x)
WRITE (iout,6200) y
WRITE (iout,6300)
STOP
!--------------------------------------------------------------------
5000 FORMAT ('1Test of J0(X) VS Taylor expansion'//)
5100 FORMAT (i7,' random arguments were tested from the interval ','(',  &
    f5.1,',',f5.1,')'//)
5200 FORMAT (' ABS(J0(X)) was larger', i6, ' times',/ t16, ' agreed', i6,  &
    ' times, and'/ t12, 'was smaller', i6, ' times.'//)
5300 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number.'//)
5400 FORMAT (' The maximum relative error of', e15.4, ' = ', i4, ' **',f7.2  &
    / t5, 'occurred for X =', e13.6)
5500 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5600 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5700 FORMAT ('1Special Tests'//)
5800 FORMAT (' Accuracy near zeros'// t11, 'X', t27, 'BESJ0(X)', t48,  &
    'Loss of base', i3, ' digits'/)
5900 FORMAT (e20.10, e25.15, t54, f7.2/)
6000 FORMAT (//' Test with extreme arguments'///)
6100 FORMAT (' J0 will be called with the argument ',e17.10/  &
    ' This may stop execution.'//)
6200 FORMAT (' J0 returned the value', e25.17/)
6300 FORMAT (' This concludes the tests.')
!---------- Last card of BESJ0 test program ----------
END PROGRAM j0test
