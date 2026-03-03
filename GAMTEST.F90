PROGRAM gamtst
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------
! FORTRAN 90 program to test GAMMA or DGAMMA

!   Method:

!      Accuracy tests compare function values against values
!      generated with the duplication formula.

!   Data required

!      None

!   Subprograms required from this package

!     MACHAR - An environmental inquiry program providing
!              information on the floating-point arithmetic
!              system.  Note that the call to MACHAR can
!              be deleted provided the following five
!              parameters are assigned the values indicated

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

!      ABS, DBLE, INT, LOG, MAX, REAL, SQRT

!  Reference: "Performance evaluation of programs related
!              to the real gamma function", W. J. Cody,
!              submitted for publication.

!  Latest modification: March 12, 1992

!  Authors: W. J. Cody and L. Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439, USA

!--------------------------------------------------------------------

USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, irnd, it, j, k1, k2, k3, machep,  &
              maxexp, minexp, n, negep, ngrd, nx
REAL (dp)  :: a, ait, albeta, alnx, b, beta, c, cl, del, eps, epsneg, r6, r7, &
              w, x, xc, xl, xmax, xmin, xminv, xn, xnum, xxn, xp, xph, y, z, zz
REAL (dp), PARAMETER  :: c1 = 2.8209479177387814347D-1,  &
                         c2 = 9.1893853320467274178D-1,  &
                         zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,   &
                         two = 2.0_dp, ten = 10.0_dp, x99 = -999.0_dp,  &
                         xp99 = 0.99_dp
INTEGER, PARAMETER    :: iout = 6

INTERFACE
  FUNCTION GAMMA(X) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val
  END FUNCTION GAMMA
END INTERFACE
!------------------------------------------------------------------
!  Determine machine parameters and set constants
!------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
a = zero
b = two
n = 2000
xn = n
!-----------------------------------------------------------------
!  Determine smallest argument for GAMMA
!-----------------------------------------------------------------
IF (xmin*xmax < one) THEN
  xminv = one / xmax
ELSE
  xminv = xmin
END IF
!-----------------------------------------------------------------
!  Determine largest argument for GAMMA by Newton iteration
!-----------------------------------------------------------------
cl = LOG(xmax)
xp = half * cl
cl = c2 - cl
10 x = xp
alnx = LOG(x)
xnum = (x-half) * alnx - x + cl
xp = x - xnum / (alnx-half/x)
IF (ABS(xp-x)/x >= ten*eps) GO TO 10
cl = xp
!-----------------------------------------------------------------
!  Random argument accuracy tests
!-----------------------------------------------------------------
DO  j = 1, 4
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
!  Use duplication formula for X not close to the zero
!-----------------------------------------------------------------
    xph = x * half + half
    xp = xph - half
    x = xp + xp
    nx = INT(x)
    xxn = nx
    c = (two**nx) * (two**(x-xxn))
    z = gamma(x)
    zz = ((c*c1)*gamma(xp)) * gamma(xph)
!--------------------------------------------------------------------
!  Accumulate results
!--------------------------------------------------------------------
    w = (z-zz) / z
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
  WRITE (iout,5000)
  WRITE (iout,5100) n, a, b
  WRITE (iout,5200) k1, k2, k3
  WRITE (iout,5300) it, ibeta
  IF (r6 /= zero) THEN
    w = LOG(ABS(r6)) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5400) r6, ibeta, w, xc
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
  IF (r7 /= zero) THEN
    w = LOG(ABS(r7)) / albeta
  ELSE
    w = x99
  END IF
  WRITE (iout,5600) r7, ibeta, w
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
!------------------------------------------------------------------
!  Initialize for next test
!------------------------------------------------------------------
  a = b
  IF (j == 1) THEN
    b = ten
  ELSE IF (j == 2) THEN
    b = cl - half
  ELSE
    a = -(ten-half) * half
    b = a + half
  END IF
END DO
!-----------------------------------------------------------------
!  Special tests
!  First test with special arguments
!-----------------------------------------------------------------
WRITE (iout,5700)
WRITE (iout,5800)
x = -half
y = gamma(x)
WRITE (iout,5900) x, y
x = xminv / xp99
y = gamma(x)
WRITE (iout,5900) x, y
x = one
y = gamma(x)
WRITE (iout,5900) x, y
x = two
y = gamma(x)
WRITE (iout,5900) x, y
x = cl * xp99
y = gamma(x)
WRITE (iout,5900) x, y
!-----------------------------------------------------------------
!  Test of error returns
!-----------------------------------------------------------------
WRITE (iout,6000)
x = -one
WRITE (iout,6100) x
y = gamma(x)
WRITE (iout,6200) y
x = zero
WRITE (iout,6100) x
y = gamma(x)
WRITE (iout,6200) y
x = xminv * (one-eps)
WRITE (iout,6100) x
y = gamma(x)
WRITE (iout,6200) y
x = cl * (one+eps)
WRITE (iout,6100) x
y = gamma(x)
WRITE (iout,6200) y
WRITE (iout,6300)
STOP
!-----------------------------------------------------------------
5000 FORMAT ('1Test of GAMMA(X) vs Duplication Formula'//)
5100 FORMAT (i7,' Random arguments were tested from the interval (',f7.3,  &
    ',',f7.3,')'//)
5200 FORMAT (' GAMMA(X) was larger',i6,' times,'/ t14, ' agreed', i6,  &
    ' times, and'/ t10, 'was smaller',i6,' times.'//)
5300 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number'//)
5400 FORMAT (' The maximum relative error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =',e13.6)
5500 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5600 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5700 FORMAT ('1Special Tests'//)
5800 FORMAT (//' Test of special arguments'//)
5900 FORMAT (' GAMMA (',e13.6,') = ',e13.6//)
6000 FORMAT ('1Test of Error Returns'///)
6100 FORMAT (' GAMMA will be called with the argument',e13.6,/  &
    ' This should trigger an error message'//)
6200 FORMAT (' GAMMA returned the value',e13.6///)
6300 FORMAT (' This concludes the tests')
!---------- Last line of GAMMA test program ----------
END PROGRAM gamtst
