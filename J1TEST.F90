PROGRAM j1test
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!--------------------------------------------------------------------
!  Fortran 90 program to test BESJ1

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - an environmental inquiry program providing
!              information on the floating-point arithmetic
!              system.  Note that the call to MACHAR can
!              be deleted provided the following six
!              parameters are assigned the values indicated

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
              cj0, cj1, d, del, eps, epsneg, r6, r7, sum, t, term,  &
              w, x, xl, xmax, xmin, xm, xn, x1, y, yinv, ysq, z, zz
!--------------------------------------------------------------------
!  Mathematical constants
!--------------------------------------------------------------------
INTEGER, PARAMETER    :: iout = 6
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, four = 4.0_dp,  &
                         delta = 0.0625_dp, eight = 8.0_dp, twenty = 20.0_dp, &
                         all9 = -999.0_dp, two56 = 256.0_dp, two = 2.0_dp,  &
                         thrten = 13.0_dp, sixten = 16.0_dp, elev = 11.0_dp
!--------------------------------------------------------------------
!  Coefficients for Taylor expansion
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: BJ0P(6, 10) = RESHAPE(  &
           (/ 1.99584D7, -2.58552D6, 1.16235D5, -2.775D3, 4.5D1, -1.0_dp,  &
              0.0_dp, -1.8144D6, 2.3688D5, -1.0845D4, 2.7D2, -5.0_dp,  &
              1.8144D5, -2.394D4, 1.125D3, -30.0_dp, 1.0_dp, 0.0_dp,  &
              0.0_dp, -2.016D4, 2.7D3, -1.32D2, 4.0_dp, 0.0_dp,  &
              2.52D3, -3.45D2, 18.0_dp, -1.0_dp, 0.0_dp, 0.0_dp,  &
              0.0_dp, -3.6D2, 51.0_dp, -3.0_dp, 0.0_dp, 0.0_dp,  &
              60.0_dp, -9.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
              0.0_dp, -12.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
              3.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
              0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /), (/ 6, 10 /) )
REAL (dp), PARAMETER  :: BJ1P(6, 10) = RESHAPE(  &
           (/ -3.99168D7, 1.016064D7, -6.7095D5, 2.067D4, -3.9D2, 6.0_dp,  &
               3.6288D6, -9.2736D5, 6.201D4, -1.965D3, 40.0_dp, -1.0_dp,  &
              -3.6288D5, 9.324D4, -6.345D3, 2.1D2, -5.0_dp, 0.0_dp,  &
               4.032D4, -1.044D4, 7.29D2, -26.0_dp, 1.0_dp, 0.0_dp,  &
              -5.04D3, 1.32D3, -96.0_dp, 4.0_dp, 0.0_dp, 0.0_dp,  &
               7.2D2, -1.92D2, 15.0_dp, -1.0_dp, 0.0_dp, 0.0_dp,  &
              -1.2D2, 33.0_dp, -3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
              24.0_dp, -7.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
              -6.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
               2.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp /), (/ 6, 10 /) )
!--------------------------------------------------------------------
!  Zeroes of J1
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: xi(2) = (/ 981.0_dp, 1796.0_dp /)
REAL (dp), PARAMETER  :: yx(2) = (/ -1.3100393001327972376D-4,  &
                                     1.1503460702301698285D-5 /)
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
DO  j = 1, 4
  k1 = 0
  k2 = 0
  k3 = 0
  x1 = zero
  r6 = zero
  r7 = zero
  n = 2000
  a = b
  IF (j == 1) THEN
    b = one
  ELSE IF (j == 2) THEN
    b = four
  ELSE IF (j == 3) THEN
    b = eight
  ELSE
    b = twenty
  END IF
  xn = n
  del = (b-a) / xn
  xl = a
  DO  i = 1, n
    x = del * ren() + xl
    IF (j == 1) THEN
!--------------------------------------------------------------------
!  Use traditional Maclaurin series for small arguments.
!--------------------------------------------------------------------
      y = x / two
      sum = y
      xm = thrten
      DO  ii = 1, 12
        sum = sum * y / xm
        xm = xm - one
        sum = (one-sum/xm) * y
      END DO
      zz = sum
      z = besj1(x)
    ELSE
!--------------------------------------------------------------------
!  Use local Taylor series elsewhere.  First, purify arguments.
!--------------------------------------------------------------------
      y = x - delta
      w = sixten * y
      y = (w+y) - w
      x = y + delta
      sum = zero
      term = zero
      bj1 = besj1(y)
      z = besj1(x)
      d = delta
      IF (ABS(z) < ABS(bj1)) THEN
        cj1 = x
        x = y
        y = cj1
        cj1 = bj1
        bj1 = z
        z = cj1
        d = -d
      END IF
      bj0 = besj0(y)
      yinv = one / y
      ysq = one / (y*y)
      xm = elev
!--------------------------------------------------------------------
!  Evaluate (12-II)th derivative at Y.
!--------------------------------------------------------------------
      DO  ii = 1, 10
        cj0 = bj0p(1,ii)
        cj1 = bj1p(1,ii)
        jj = (12-ii) / 2 + 1
        DO  iii = 2, jj
          cj0 = cj0 * ysq + bj0p(iii,ii)
          cj1 = cj1 * ysq + bj1p(iii,ii)
        END DO
        IF ((ii/2)*2 == ii) THEN
          cj0 = cj0 * yinv
        ELSE
          cj1 = cj1 * yinv
        END IF
        term = cj0 * bj0 + cj1 * bj1
        sum = (sum+term) * d / xm
        xm = xm - one
      END DO
      sum = (sum+bj0-bj1*yinv) * d + bj1
      zz = sum
    END IF
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
  IF (j == 1) THEN
    WRITE (iout,5000)
  ELSE
    WRITE (iout,5100)
  END IF
  WRITE (iout,5200) n, a, b
  WRITE (iout,5300) k1, k2, k3
  WRITE (iout,5400) it, ibeta
  IF (r6 /= zero) THEN
    w = LOG(ABS(r6)) / albeta
  ELSE
    w = all9
  END IF
  WRITE (iout,5500) r6, ibeta, w, x1
  w = MAX(ait+w,zero)
  WRITE (iout,5600) ibeta, w
  IF (r7 /= zero) THEN
    w = LOG(ABS(r7)) / albeta
  ELSE
    w = all9
  END IF
  WRITE (iout,5700) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,5600) ibeta, w
END DO
!--------------------------------------------------------------------
!  Special tests
!--------------------------------------------------------------------
WRITE (iout,5800)
WRITE (iout,5900) ibeta
DO  i = 1, 2
  x = xi(i) / two56
  y = besj1(x)
  t = (y-yx(i)) / yx(i)
  IF (t /= zero) THEN
    w = LOG(ABS(t)) / albeta
  ELSE
    w = all9
  END IF
  w = MAX(ait+w,zero)
  WRITE (iout,6000) x, y, w
END DO
!--------------------------------------------------------------------
!  Test of error returns
!--------------------------------------------------------------------
WRITE (iout,6100)
x = xmax
WRITE (iout,6200) x
y = besj1(x)
WRITE (iout,6300) y
WRITE (iout,6400)
STOP
!--------------------------------------------------------------------
5000 FORMAT ('1Test of J1(X) VS Maclaurin expansion'//)
5100 FORMAT ('1Test of J1(X) VS local Taylor expansion'//)
5200 FORMAT (i7,' random arguments were tested from the interval (',  &
    f5.1,',',f5.1,')'//)
5300 FORMAT (' ABS(J1(X)) was larger', i6, ' times',/ t16, ' agreed', i6,  &
    ' times, and'/ t12, 'was smaller', i6, ' times.'//)
5400 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number.'//)
5500 FORMAT (' The maximum relative error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =', e13.6)
5600 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5700 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5800 FORMAT ('1Special Tests'//)
5900 FORMAT (' Accuracy near zeros'// t11, 'X', t26, 'BESJ1(X)', t47,  &
    'Loss of base',i3,' digits'/)
6000 FORMAT (e20.10, e25.15, t54, f7.2/)
6100 FORMAT (//' Test with extreme arguments'///)
6200 FORMAT (' J1 will be called with the argument ', e17.10/  &
    ' This may stop execution.'//)
6300 FORMAT (' J1 returned the value', e25.17/)
6400 FORMAT (' This concludes the tests.')
!---------- Last card of BESJ1 test program ----------
END PROGRAM j1test
