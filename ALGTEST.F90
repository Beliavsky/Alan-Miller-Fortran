PROGRAM algtst
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------
! FORTRAN 90 program to test ALGAMA

!   Method:

!      Accuracy tests compare function values against values
!      generated with the duplication formula.

!   Data required

!      None

!   Subprograms required from this package

!               IBETA  - The radix of the floating-point system
!               IT     - The number of base-ibeta digits in the
!                        significant of a floating-point number
!               EPS    - The smallest positive floating-point
!                        number such that 1.0+EPS .NE. 1.0
!               XMIN   - The smallest non-vanishing floating-point
!                        integral power of the radix
!               XMAX   - The largest finite floating-point number

!      REN(K) - A function subprogram returning random real
!               numbers uniformly distributed over (0,1)


!    Intrinsic functions required are:

!         ABS, ANINT, DBLE, LOG, MAX, REAL, SQRT

!  Reference: "Performance evaluation of programs related to the real gamma
!              function", W. J. Cody, submitted for publication.

!  Latest modification: March 9, 1992

!  Authors: W. J. Cody and L. Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439

!--------------------------------------------------------------------
USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, irnd, it, j, k1, k2, k3, machep,  &
              maxexp, minexp, n, negep, ngrd

INTERFACE
  FUNCTION ALGAMA(X) RESULT(res)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: res
  END FUNCTION ALGAMA
END INTERFACE

REAL (dp)  :: a, ait, albeta, b, beta, cl, del, eps, epsneg,  &
              r6, r7, u, v, w, x, xc, xl, xmax, xmin, xn, y, z, zz
!------------------------------------------------------------------
!   C1 = 0.5 - LN(SQRT(PI))
!   C2 = LN(2)
!   C3 = LN(2) - 11/16
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: C1 = -7.2364942924700087072D-2,   &
                         C2 = 6.9314718055994530942D-1,   &
                         C3 = 5.6471805599453094172D-3, zero = 0.0D0,   &
                         half = 0.5D0, one = 1.0_dp, two = 2.0_dp,  &
                         ten = 10.0_dp, sixten = 16.0_dp, p6875 = 0.6875_dp, &
                         p875 = 0.875_dp, p3125 = 1.3125_dp, p625 = 1.625_dp, &
                         all9 = -999.0_dp, xp99 = 0.99_dp
INTEGER, PARAMETER    :: iout = 6
!------------------------------------------------------------------
!  Determine machine parameters and set constants
!------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
a = zero
b = p875
n = 2000
xn = n
!-----------------------------------------------------------------
!  Determine largest argument for DLGAMA by iteration
!-----------------------------------------------------------------
cl = xp99 * xmax
z = -cl / all9
10 zz = cl / (LOG(z)-one)
IF (ABS(zz/z-one) > (two*beta*eps)) THEN
  z = zz
  GO TO 10
END IF
cl = zz
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
!-----------------------------------------------------------------
!  Use duplication formula
!-----------------------------------------------------------------
    IF (j /= 3) THEN
      IF (j == 1) THEN
        z = x + half
        x = z - half
        y = x + x
      ELSE
        x = x + x
        x = x * half
        y = (x+x) - one
        z = x - half
      END IF
      u = algama(x)
      w = (y-half) - half
      zz = ANINT(w*sixten) / sixten
      w = w - zz
      v = (((half - zz*p6875) - c1) - w*p6875) - c3 * (w+zz)
      v = ((v+algama(y)) - algama(z))
    ELSE
      z = x * half + half
      y = z - half
      x = y + y
      u = algama(x)
      v = (c1 + ((x-half) - half)*c2) + algama(y) + algama(z) - half
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
  ELSE IF (j == 2) THEN
    WRITE (iout,5100)
  ELSE
    WRITE (iout,5200)
  END IF
  WRITE (iout,5300) n, a, b
  WRITE (iout,5400) k1, k2, k3
  WRITE (iout,5500) it, ibeta
  IF (r6 /= zero) THEN
    w = LOG(ABS(r6)) / albeta
  ELSE
    w = all9
  END IF
  WRITE (iout,5600) r6, ibeta, w, xc
  w = MAX(ait+w,zero)
  WRITE (iout,5700) ibeta, w
  IF (r7 /= zero) THEN
    w = LOG(ABS(r7)) / albeta
  ELSE
    w = all9
  END IF
  WRITE (iout,5800) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,5700) ibeta, w
!------------------------------------------------------------------
!  Initialize for next test
!------------------------------------------------------------------
  a = p3125
  b = p625
  IF (j == 2) THEN
    a = two + two
    b = ten + ten
  END IF
END DO
!-----------------------------------------------------------------
!  Special tests
!  First test with special arguments
!-----------------------------------------------------------------
WRITE (iout,5900)
WRITE (iout,6000)
z = eps
zz = algama(z)
WRITE (iout,6100) z, zz
z = half
zz = algama(z)
WRITE (iout,6100) z, zz
z = one
zz = algama(z)
WRITE (iout,6100) z, zz
z = two
zz = algama(z)
WRITE (iout,6100) z, zz
!-----------------------------------------------------------------
!  Test of error returns
!-----------------------------------------------------------------
WRITE (iout,6200)
z = xmin
WRITE (iout,6400) z
zz = algama(z)
WRITE (iout,6500) zz
z = cl
WRITE (iout,6400) z
zz = algama(z)
WRITE (iout,6500) zz
z = -one
WRITE (iout,6300) z
zz = algama(z)
WRITE (iout,6500) zz
z = zero
WRITE (iout,6300) z
zz = algama(z)
WRITE (iout,6500) zz
z = xp99 * xmax
WRITE (iout,6300) z
zz = algama(z)
WRITE (iout,6500) zz
WRITE (iout,6600)
STOP

!-----------------------------------------------------------------
5000 FORMAT ('1Test of LGAMA(X) vs LN(2*SQRT(PI))-2X*LN(2)+',  &
    'LGAMA(2X)-LGAMA(X+1/2)'//)
5100 FORMAT ('1Test of LGAMA(X) vs LN(2*SQRT(PI))-(2X-1)*LN(2)+',  &
    'LGAMA(X-1/2)-LGAMA(2X-1)'//)
5200 FORMAT ('1Test of LGAMA(X) vs -LN(2*SQRT(PI))+X*LN(2)+',  &
    'LGAMA(X/2)+LGAMA(X/2+1/2)'//)
5300 FORMAT (i7,' Random arguments were tested from the interval (',f5  &
    .1,',',f5.1,')'//)
5400 FORMAT ('  LGAMA(X) was larger',i6,' times,'/ t15,' agreed',i6,  &
    ' times, and'/ t11,'was smaller',i6,' times.'//)
5500 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number'//)
5600 FORMAT (' The maximum relative error of',e15.4,' = ',i4,' **',f7.2  &
    /t5, 'occurred for X =', e13.6)
5700 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5800 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5900 FORMAT ('1Special Tests'//)
6000 FORMAT (//' Test of special arguments'//)
6100 FORMAT (' LGAMA (',e13.6,') = ',e13.6//)
6200 FORMAT ('1Test of Error Returns'///)
6300 FORMAT (' LGAMA will be called with the argument',e13.6,/  &
    ' This should trigger an error message'//)
6400 FORMAT (' LGAMA will be called with the argument',e13.6,/  &
    ' This should not trigger an error message'//)
6500 FORMAT (' LGAMA returned the value',e13.6///)
6600 FORMAT (' This concludes the tests')
!---------- Last line of ALGAMA test program ----------
END PROGRAM algtst
