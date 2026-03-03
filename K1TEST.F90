PROGRAM k1test
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!--------------------------------------------------------------------
!  FORTRAN 90 program to test BESK1

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - an environmental inquiry program providing information on the
!              floating-point arithmetic system.  Note that the call to
!              MACHAR can be deleted provided the following three
!              parameters are assigned the values indicated

!                 IBETA  - the radix of the floating-point system
!                 IT     - the number of base-IBETA digits in the
!                          significand of a floating-point number
!                 MAXEXP - the smallest positive power of BETA that overflows
!                 EPS    - the smallest positive floating-point
!                          number such that 1.0+EPS .NE. 1.0
!                 XMIN   - the smallest non-vanishing normalized
!                          floating-point power of the radix, i.e.,
!                          XMIN = FLOAT(IBETA) ** MINEXP
!                 XMAX   - the largest finite floating-point number.
!                          In particular XMAX = (1.0-EPSNEG) *
!                          FLOAT(IBETA) ** MAXEXP

!     REN(K) - a function subprogram returning random real
!              numbers uniformly distributed over (0,1)


!  Intrinsic functions required are:

!         ABS, DBLE, LOG, MAX, REAL, SQRT

!  User defined functions

!         BOT, TOP

!  Reference: "Performance evaluation of programs for certain
!              Bessel functions", W. J. Cody and L. Stoltz,
!              ACM Trans. on Math. Software, Vol. 15, 1989,
!              pp 41-48.

!  Latest modification: March 14, 1992

!  Author - Laura Stoltz
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439, USA
!--------------------------------------------------------------------
USE Bessel_K0
USE Bessel_K1
USE Toms715_Utilities

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

LOGICAL    :: sflag, tflag
INTEGER    :: i, ibeta, iexp, ii, ind, irnd, it, j, k1, k2,  &
              k3, machep, maxexp, mb, minexp, n, ndum, negep, ngrd
REAL (dp)  :: a, ait, albeta, amaxexp, b, beta, c, const, del, eps, epsneg,  &
              r6, r7, sum, t, t1, w, x, xa, xb,  &
              xl, xlam, xlarge, xmax, xmb, xmin, xn, x1, y, z, zz
REAL (dp)  :: u(0:559)
!--------------------------------------------------------------------
!  Mathematical constants
!--------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         eight = 8.0_dp, xnine = 9.0_dp, twenty = 20.0_dp, &
                         five = 5.0_dp, xden = 16.0_dp, all9 = -999.0_dp, &
                         pi = 3.141592653589793_dp, hund = 100.0_dp
!---------------------------------------------------------------------
!  Machine-dependent constant
!---------------------------------------------------------------------
REAL (dp), PARAMETER  :: XLEAST = 2.23D-308
INTEGER, PARAMETER    :: iout = 6
!---------------------------------------------------------------------

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

!--------------------------------------------------------------------
!  Determine machine parameters and set constants

CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
amaxexp = maxexp
b = eps
xlam = (xden-one) / xden
const = half * LOG(pi) - LOG(xmin)
!--------------------------------------------------------------------
!     Random argument accuracy tests
!--------------------------------------------------------------------
DO  j = 1, 3
  sflag = ((j == 1).AND.(amaxexp/ait <= five))
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
    b = eight
  ELSE
    b = twenty
  END IF
  xn = n
  del = (b-a) / xn
  xl = a
!---------------------------------------------------------------------
!   Accuracy test is based on the multiplication theorem
!---------------------------------------------------------------------
  loop40:  DO  i = 1, n
    x = del * ren() + xl
    y = x / xlam
    w = xden * y
    y = (w+y) - w
    x = y * xlam
    u(0) = besk0(y)
    u(1) = besk1(y)
    tflag = sflag .AND. (y < half)
    IF (tflag) THEN
      u(0) = u(0) * eps
      u(1) = u(1) * eps
    END IF
    mb = 1
    xmb = one
    y = y * half
    w = (one-xlam) * (one+xlam)
    c = w * y
    t = u(0) + c * u(1)
    t1 = eps / hund
    DO  ii = 2, 60
      z = u(ii-1)
      IF (z/t1 < t) THEN
        EXIT
      ELSE IF (u(ii-1) > one) THEN
        IF ((xmb/y) > (xmax/u(ii-1))) THEN
          xl = xl + del
          a = xl
          CYCLE loop40
        END IF
      END IF
      u(ii) = xmb / y * u(ii-1) + u(ii-2)
      IF (t1 > one/eps) THEN
        t = t * t1
        t1 = one
      END IF
      t1 = xmb * t1 / c
      xmb = xmb + one
      mb = mb + 1
    END DO
    sum = u(mb)
    ind = mb
    mb = mb - 1
    DO  ii = 1, mb
      xmb = xmb - one
      ind = ind - 1
      sum = sum * w * y / xmb + u(ind)
    END DO
    zz = sum
    IF (tflag) zz = zz / eps
    zz = zz * xlam
    z = besk1(x)
    y = z
    IF (u(0) > y) y = u(0)
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
  n = k1 + k2 + k3
  xn = n
  r7 = SQRT(r7/xn)
  WRITE (iout,5000)
  WRITE (iout,5100) n, a, b
  WRITE (iout,5200) k1, k2, k3
  WRITE (iout,5300) it, ibeta
  w = all9
  IF (r6 /= zero) w = LOG(r6) / albeta
  IF (j == 3) THEN
    WRITE (iout,5700) r6, ibeta, w, x1
  ELSE
    WRITE (iout,5400) r6, ibeta, w, x1
  END IF
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
  w = all9
  IF (r7 /= zero) w = LOG(r7) / albeta
  IF (j == 3) THEN
    WRITE (iout,5800) r7, ibeta, w
  ELSE
    WRITE (iout,5600) r7, ibeta, w
  END IF
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
END DO
!--------------------------------------------------------------------
!  Special tests
!--------------------------------------------------------------------
WRITE (iout,5900)
WRITE (iout,6000)
y = besk1(xleast)
WRITE (iout,6100) y
y = besk1(xmin)
WRITE (iout,6500) y
y = besk1(zero)
WRITE (iout,6200) 0, y
x = - ren()
y = besk1(x)
WRITE (iout,6300) x, y
y = besek1(xmax)
WRITE (iout,6400) y
xa = LOG(xmax)
60 xb = xa - (top(xa) + const) / bot(xa)
IF (ABS(xb-xa)/xb <= eps) THEN
  GO TO 70
ELSE
  xa = xb
  GO TO 60
END IF
70 xlarge = xb * xlam
y = besk1(xlarge)
WRITE (iout,6300) xlarge, y
xlarge = xb * (xnine/eight)
y = besk1(xlarge)
WRITE (iout,6300) xlarge, y
!--------------------------------------------------------------------
!  Test of error returns
!--------------------------------------------------------------------
STOP
5000 FORMAT ('1Test of K1(X) vs Multiplication Theorem'//)
5100 FORMAT (i7,' random arguments were tested from the interval (', f5.1,  &
    ',', f5.1,')'//)
5200 FORMAT (' ABS(K1(X)) was larger',i6, ' times,'/ t21,' agreed', i6,  &
    ' times, and'/ t17, 'was smaller', i6,' times.'//)
5300 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number.'//)
5400 FORMAT (' The maximum relative error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =',e13.6)
5500 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5600 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5700 FORMAT (' The maximum absolute error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =',e13.6)
5800 FORMAT (' The root mean square absolute error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5900 FORMAT ('1Special Tests'//)
6000 FORMAT (//' Test with extreme arguments'/)
6100 FORMAT (' K1(XLEAST) = ',e24.17/)
6200 FORMAT (' K1(',i1,') = ',e24.17/)
6300 FORMAT (' K1(',e24.17,' ) = ',e24.17/)
6400 FORMAT (' E**X * K1(XMAX) = ',e24.17/)
6500 FORMAT (' K1(XMIN) = ',e24.17/)
!---------- Last line of BESK0 test program ----------
END PROGRAM k1test



FUNCTION top(x) RESULT(fn_val)
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

REAL (dp), PARAMETER   :: one = 1.0_dp, eight = 8.0_dp, one28 = 128.0_dp, &
                          half = 0.5_dp, fiften = 15.0_dp, two = 2.0_dp,  &
                          three = 3.0_dp

fn_val = -x - half * LOG(two*x) + LOG(one + (three/eight-fiften/one28/x)/x)
RETURN
END FUNCTION top



FUNCTION bot(x) RESULT(fn_val)
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

REAL (dp), PARAMETER   :: one = 1.0_dp, one28 = 128.0_dp, fiften = 15.0_dp,  &
                          half = 0.5_dp, four8 = 48.0_dp, thirty = 30.0_dp

fn_val = -one - half / x +   &
         ((-four8*x + thirty)/(((one28*x + four8)*x - fiften)*x))
RETURN
END FUNCTION bot

