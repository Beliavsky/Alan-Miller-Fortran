PROGRAM psitst
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!--------------------------------------------------------------------
!  Fortran 77 program to test PSI

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - an environmental inquiry program providing information on
!              the floating-point arithmetic system.  Note that the call
!              system.  Note that the call to MACHAR can be deleted provided
!              the following five parameters are assigned the values indicated

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


!  Intrinsic Fortran functions required are:

!         ABS, DBLE, LOG, MAX, REAL, SQRT

!  Reference: "Performance evaluation of programs related to the real gamma
!              function", W. J. Cody, submitted for publication.

!  Latest modification: March 14, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!--------------------------------------------------------------------

USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, irnd, it, j, k1, k2, k3, machep,  &
              maxexp, minexp, n, negep, ngrd
REAL (dp)  :: a, ait, albeta, b, beta, del, eps, epsneg, r6, r7, y,  &
              w, x, xh, xl, xmax, xmin, xn, xx, x1, z, zh, zz
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, three = 3.0_dp,  &
                         half = 0.5_dp, eight = 8.0_dp, twenty = 20.0_dp,  &
                         all9 = -999.0_dp, xl2 = 6.9314718055994530942D-1,  &
                         one7 = -17.625_dp, one6 = -16.875_dp, x0 = 374.0_dp, &
                         x01 = 256.0_dp, v0 = -6.7240239024288040437D-04
INTEGER, PARAMETER  :: iout = 6

INTERFACE
  FUNCTION PSI(XX) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: xx
    REAL (dp)              :: fn_val
  END FUNCTION PSI
END INTERFACE

!--------------------------------------------------------------------
!  Determine machine parameters and set constants
!--------------------------------------------------------------------

CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
albeta = LOG(beta)
ait = it
!--------------------------------------------------------------------
!     Random argument accuracy tests
!--------------------------------------------------------------------
DO  j = 1, 4
  k1 = 0
  k3 = 0
  x1 = zero
  r6 = zero
  r7 = zero
  n = 2000
  IF (j == 1) THEN
    a = zero
    b = one
  ELSE IF (j == 2) THEN
    a = b + b
    b = eight
  ELSE IF (j == 3) THEN
    a = b
    b = twenty
  ELSE
    a = one7
    b = one6
    n = 500
  END IF
  xn = n
  del = (b-a) / xn
  xl = a
  DO  i = 1, n
    x = del * ren() + xl
!--------------------------------------------------------------------
!  Carefully purify arguments and evaluate identity
!--------------------------------------------------------------------
    xx = x * half
    xh = xx + half
    xx = xh - half
    x = xx + xx
    z = psi(x)
    zh = psi(xh)
    zz = psi(xx)
    zz = (zz+zh) * half + xl2
!--------------------------------------------------------------------
!  Accumulate results
!--------------------------------------------------------------------
    w = (zz-z) / zz
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
!--------------------------------------------------------------------
!  Process and output statistics
!--------------------------------------------------------------------
  k2 = n - k3 - k1
  r7 = SQRT(r7/xn)
  IF (2*(j/2) /= j) WRITE (iout,5000)
  WRITE (iout,5100)
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
x = x0 / x01
y = psi(x)
z = (y-v0) / v0
IF (z /= zero) THEN
  w = LOG(ABS(z)) / albeta
ELSE
  w = all9
END IF
w = MAX(ait+w,zero)
WRITE (iout,5900) x, y, ibeta, w
WRITE (iout,6000)
IF (xmax*xmin >= one) THEN
  x = xmin
ELSE
  x = one / xmax
END IF
WRITE (iout,6200) x
y = psi(x)
WRITE (iout,6300) y
x = xmax
WRITE (iout,6200) x
y = psi(x)
WRITE (iout,6300) y
!--------------------------------------------------------------------
!  Test of error returns
!--------------------------------------------------------------------
WRITE (iout,6400)
x = zero
WRITE (iout,6100) x
y = psi(x)
WRITE (iout,6300) y
x = -three / eps
WRITE (iout,6100) x
y = psi(x)
WRITE (iout,6300) y
WRITE (iout,6500)
STOP
!--------------------------------------------------------------------
5000 FORMAT ('1')
5100 FORMAT (' Test of PSI(X) vs (PSI(X/2)+PSI(X/2+1/2))/2 + ln(2)'//)
5200 FORMAT (i7, ' random arguments were tested from the interval (',  &
    f5.1, ',', f5.1, ')'//)
5300 FORMAT (' ABS(PSI(X)) was larger', i6, ' times'/ t22, ' agreed', i6,  &
    ' times, and'/ t18, 'was smaller', i6, ' times.'//)
5400 FORMAT (' There are', i4, ' base', i4,  &
    ' significant digits in a floating-point number.'//)
5500 FORMAT (' The maximum relative error of', e15.4, ' = ', i4, ' **', f7.2  &
    /t5, 'occurred for X =', e13.6)
5600 FORMAT (' The estimated loss of base', i4, ' significant digits is',  &
    f7.2//)
5700 FORMAT (' The root mean square relative error was', e15.4, ' = ', i4,  &
    ' **', f7.2)
5800 FORMAT ('1Special Tests'//)
5900 FORMAT (' Accuracy near positive zero'//' PSI(',e14.7,') = ',e24.17  &
    / t14, 'Loss of base', i3, ' digits = ', f7.2/)
6000 FORMAT (//' Test with extreme arguments'/)
6100 FORMAT (' PSI will be called with the argument ',e17.10/  &
    ' This may stop execution.'/)
6200 FORMAT (' PSI will be called with the argument ',e17.10/  &
    ' This should not stop execution.'/)
6300 FORMAT (' PSI returned the value',e25.17//)
6400 FORMAT (//' Test of error returns'//)
6500 FORMAT (' This concludes the tests.')
!---------- Last card of PSI test program ----------
END PROGRAM psitst
