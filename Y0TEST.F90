PROGRAM y0test
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!--------------------------------------------------------------------
!  Fortran 90 program to test BESY0

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - an environmental inquiry program providing information on the
!              floating-point arithmetic system.
!              Note that the call to MACHAR can be deleted provided the
!              following five parameters are assigned the values indicated

!                 IBETA  - the radix of the floating-point system
!                 IT     - the number of base-IBETA digits in the
!                          significand of a floating-point number
!                 MINEXP - the largest in magnitude negative
!                          integer such that  FLOAT(IBETA)**MINEXP
!                          is a positive floating-point number
!                 EPS    - the smallest positive floating-point
!                          number such that 1.0+EPS .NE. 1.0
!                 EPSNEG - the smallest positive floating-point
!                          number such that 1.0-EPSNEG .NE. 1.0

!     REN(K) - a function subprogram returning random real
!              numbers uniformly distributed over (0,1)


!  Intrinsic functions required are:

!         ABS, DBLE, LOG, MAX, REAL, SQRT

!  Reference: "Performance evaluation of programs for certain
!              Bessel functions", W. J. Cody and L. Stoltz,
!              ACM Trans. on Math. Software, Vol. 15, 1989, pp 41-48.

!  Latest modification: March 13, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439, USA

!--------------------------------------------------------------------
USE Bessel_JY0
USE Bessel_JY1
USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

LOGICAL    :: sflag, tflag
INTEGER    :: i, ibeta, iexp, ii, ind, irnd, it, j, k1, k2,  &
              k3, machep, maxexp, mb, minexp, n, negep, ngrd
REAL (dp)  :: a, ait, albeta, b, beta, c, del, eps, epsneg, r6, r7,  &
              sum, t, t1, w, x, xl, xmax, xmb, xmin, x1, y, z, zz
REAL (dp)  :: u(560)
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, three = 3.0_dp,  &
                         five5 = 5.5_dp, eight = 8.0_dp, twenty = 20.0_dp,  &
                         all9 = -999.0_dp, two56 = 256.0_dp, half = 0.5_dp, &
                         sixten = 16.0_dp, xlam = 0.9375_dp, hund = 100.0_dp
REAL (dp), PARAMETER  :: xi(3) = (/ 228.0D0, 1013.0D0, 1814.0D0 /)
REAL (dp), PARAMETER  :: yx(3) = (/ -2.6003142722933487915D-3,  &
                                     2.6053454911456774983D-4,  &
                                    -3.4079448714795552077D-5 /)
INTEGER, PARAMETER  :: iout = 6
!--------------------------------------------------------------------
!  Determine machine parameters and set constants
!--------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
a = eps
b = three
n = 2000
!--------------------------------------------------------------------
!  Random argument accuracy tests based on the multiplication theorem
!--------------------------------------------------------------------
DO  j = 1, 4
  sflag = (j == 1) .AND. (maxexp/it <= 5)
  k1 = 0
  k2 = 0
  k3 = 0
  x1 = zero
  r6 = zero
  r7 = zero
  del = (b-a) / n
  xl = a
  loop40:  DO  i = 1, n
    x = del * ren() + xl
!--------------------------------------------------------------------
!  Carefully purify arguments
!--------------------------------------------------------------------
    y = x / xlam
    w = sixten * y
    y = (w+y) - w
    x = y * xlam
!--------------------------------------------------------------------
!  Generate Bessel functions with forward recurrence
!--------------------------------------------------------------------
    u(1) = besy0(y)
    u(2) = besy1(y)
    tflag = sflag .AND. (y < half)
    IF (tflag) THEN
      u(1) = u(1) * eps
      u(2) = u(2) * eps
    END IF
    mb = 1
    xmb = one
    y = y * half
    w = (one-xlam) * (one+xlam)
    c = w * y
    t = ABS(u(1)+c*u(2))
    t1 = eps / hund
    DO  ii = 3, 60
      z = ABS(u(ii-1))
      IF (z/t1 < t) EXIT
      IF (y < xmb) THEN
        IF (z > xmax*(y/xmb)) THEN
          a = x
          xl = xl + del
          CYCLE loop40
        END IF
      END IF
      u(ii) = xmb / y * u(ii-1) - u(ii-2)
      IF (t1 > one/eps) THEN
        t = t * t1
        t1 = one
      END IF
      t1 = xmb * t1 / c
      xmb = xmb + one
      mb = mb + 1
    END DO
!--------------------------------------------------------------------
!  Evaluate Bessel series expansion
!--------------------------------------------------------------------
    sum = u(mb)
    ind = mb
    DO  ii = 2, mb
      ind = ind - 1
      xmb = xmb - one
      sum = sum * w * y / xmb + u(ind)
    END DO
    zz = sum
    IF (tflag) zz = zz / eps
    z = besy0(x)
    y = z
    IF (ABS(u(1)) > ABS(y)) y = u(1)
!--------------------------------------------------------------------
!  Accumulate results
!--------------------------------------------------------------------
    w = (z-zz) / y
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
  END DO loop40
!--------------------------------------------------------------------
!  Gather and print statistics
!--------------------------------------------------------------------
  n = k1 + k2 + k3
  r7 = SQRT(r7/n)
  WRITE (iout,5000)
  WRITE (iout,5100) n, a, b
  WRITE (iout,5200) k1, k2, k3
  WRITE (iout,5300) it, ibeta
  w = all9
  IF (r6 /= zero) w = LOG(ABS(r6)) / albeta
  IF (j == 4) THEN
    WRITE (iout,5700) r6, ibeta, w, x1
  ELSE
    WRITE (iout,5400) r6, ibeta, w, x1
  END IF
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
  w = all9
  IF (r7 /= zero) w = LOG(ABS(r7)) / albeta
  IF (j == 4) THEN
    WRITE (iout,5800) r7, ibeta, w
  ELSE
    WRITE (iout,5600) r7, ibeta, w
  END IF
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
!--------------------------------------------------------------------
!  Initialize for next test
!--------------------------------------------------------------------
  a = b
  IF (j == 1) THEN
    b = five5
    n = 2000
  ELSE IF (j == 2) THEN
    b = eight
    n = 2000
  ELSE
    b = twenty
    n = 500
  END IF
END DO
!--------------------------------------------------------------------
!  Special tests
!--------------------------------------------------------------------
WRITE (iout,5900)
WRITE (iout,6000) ibeta
DO  i = 1, 3
  x = xi(i) / two56
  y = besy0(x)
  w = all9
  t = (y-yx(i)) / yx(i)
  IF (t /= zero) w = LOG(ABS(t)) / albeta
  w = MAX(ait+w,zero)
  WRITE (iout,6100) x, y, w
END DO
!--------------------------------------------------------------------
!  Test of error returns
!--------------------------------------------------------------------
WRITE (iout,6200)
x = xmin
WRITE (iout,6400) x
y = besy0(x)
WRITE (iout,6500) y
x = zero
WRITE (iout,6300) x
y = besy0(x)
WRITE (iout,6500) y
x = xmax
WRITE (iout,6300) x
y = besy0(x)
WRITE (iout,6500) y
WRITE (iout,6600)
STOP
5000 FORMAT ('1Test of Y0(X) VS Multiplication Theorem'//)
5100 FORMAT (i7,' random arguments were tested from the interval (',  &
    f5.1,',',f5.1,')'//)
5200 FORMAT (' ABS(Y0(X)) was larger', i6,' times'/ t16,' agreed', i6,  &
    ' times, and'/ t12, 'was smaller', i6,' times.'//)
5300 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number.'//)
5400 FORMAT (' The maximum relative error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =', e13.6)
5500 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5600 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5700 FORMAT (' The maximum absolute error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =', e13.6)
5800 FORMAT (' The root mean square absolute error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5900 FORMAT ('1Special Tests'//)
6000 FORMAT (' Accuracy near zeros'// t11, 'X', t27, 'BESY0(X)', t48,  &
    'Loss of base',i3,' digits'/)
6100 FORMAT (e20.10, e25.15, t54, f7.2/)
6200 FORMAT (//' Test with extreme arguments'/)
6300 FORMAT (' Y0 will be called with the argument ',e17.10/  &
    ' This may stop execution.'//)
6400 FORMAT (' Y0 will be called with the argument ',e17.10/  &
    ' This should not stop execution.'//)
6500 FORMAT (' Y0 returned the value',e25.17/)
6600 FORMAT (' This concludes the tests.')
!---------- Last card of BESY0 test program ----------
END PROGRAM y0test
