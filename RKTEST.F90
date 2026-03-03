PROGRAM rktest
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!--------------------------------------------------------------------
! FORTRAN 90 program to test RKBESL

!  Data required

!     None

!  Subprograms required from this package

!     MACHAR - An environmental inquiry program providing information on the
!         floating-point arithmetic system.  Note that the call to MACHAR
!         can be deleted provided the following five parameters are
!         assigned the values indicated

!              IBETA  - the radix of the floating-point system
!              IT     - the number of base-IBETA digits in the
!                       significand of a floating-point number
!              MINEXP - the largest in magnitude negative
!                       integer such that FLOAT(IBETA)**MINEXP
!                       is a positive floating-point number
!              EPS    - the smallest positive floating-point
!                       number such that 1.0+EPS .NE. 1.0
!              EPSNEG - the smallest positive floating-point
!                       number such that 1.0-EPSNEG .NE. 1.0

!     REN(K) - a function subprogram returning random real
!              numbers uniformly distributed over (0,1)


!  Intrinsic functions required are:

!      ABS, DBLE, EXP, LOG, MAX, REAL, SQRT

!  Reference: "Performance evaluation of programs for certain
!              Bessel functions", W. J. Cody and L. Stoltz,
!              ACM Trans. on Math. Software, Vol. 15, 1989, pp 41-48.

!  Latest modification: March 14, 1992

!  Authors: W. J. Cody and L. Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439, USA

!--------------------------------------------------------------------

USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

LOGICAL    :: sflag, tflag
INTEGER    :: i, ibeta, iexp, ii, ind, irnd, it, ize, iz1, j,   &
              k1, k2, k3, machep, maxexp, mb, minexp, n, ncalc, negep, ngrd
REAL (dp)  :: a, ait, albeta, alpha, amaxexp, a1, b, beta, d, del,  &
              eps, epsneg, ovrchk, r6, r7, sum, t, t1, w, x, xl,  &
              xmax, xmb, xmin, xn, x1, y, z, zz
REAL (dp)  :: u(560), u2(560)
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
           five = 5.0_dp, eight = 8.0_dp, xnine = 9.0_dp, ten = 10.0_dp,  &
           hund = 100.0_dp, x99 = -999.0_dp, sixten = 16.0_dp,   &
           xlam = 0.9375_dp, c = 0.22579_dp, tinyx = 1.0D-10
INTEGER, PARAMETER  :: iout = 6

INTERFACE
  SUBROUTINE rkbesl(x, alpha, nb, ize, bk, ncalc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x
    REAL (dp), INTENT(IN)   :: alpha
    INTEGER, INTENT(IN)     :: nb
    INTEGER, INTENT(IN)     :: ize
    REAL (dp), INTENT(OUT)  :: bk(:)
    INTEGER, INTENT(OUT)    :: ncalc
  END SUBROUTINE rkbesl
END INTERFACE

!--------------------------------------------------------------------
!  Determine machine parameters and set constants
!--------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
ait = it
albeta = LOG(beta)
amaxexp = maxexp
a = eps
b = one
!-------------------------------------------------------------------
!  Random argument accuracy tests
!-------------------------------------------------------------------
DO  j = 1, 3
  sflag = (j == 1 .AND. amaxexp/ait <= five)
  n = 2000
  xn = n
  k1 = 0
  k2 = 0
  k3 = 0
  x1 = zero
  a1 = zero
  r6 = zero
  r7 = zero
  del = (b-a) / xn
  xl = a
  loop40:  DO  i = 1, n
    x = del * ren() + xl
    alpha = ren()
!--------------------------------------------------------------------
!   Accuracy test is based on the multiplication theorem
!--------------------------------------------------------------------
    ize = 1
    mb = 3
!--------------------------------------------------------------------
!   Carefully purify arguments
!--------------------------------------------------------------------
    y = x / xlam
    w = sixten * y
    t1 = w + y
    y = t1 - w
    x = y * xlam
    CALL rkbesl(y,alpha,mb,ize,u2,ncalc)
    tflag = sflag  .AND.  (y < half)
    IF (tflag) THEN
      u2(1) = u2(1) * eps
      u2(2) = u2(2) * eps
    END IF
    mb = 1
    xmb = zero
    w = (one-xlam) * (one+xlam) * half
    d = w * y
    t = u2(1) + d * u2(2)
    t1 = eps / hund
!--------------------------------------------------------------------
!   Generate terms using recurrence
!--------------------------------------------------------------------
    DO  ii = 3, 35
      xmb = xmb + one
      t1 = xmb * t1 / d
      z = u2(ii-1)
      ovrchk = (xmb+alpha) / (y*half)
      IF (z/t1 < t) THEN
        GO TO 20
      ELSE IF (u2(ii-1) > one) THEN
        IF (ovrchk > (xmax/u2(ii-1))) THEN
          xl = xl + del
          a = xl
          CYCLE loop40
        END IF
      END IF
      u2(ii) = ovrchk * u2(ii-1) + u2(ii-2)
      IF (t1 > one/eps) THEN
        t = t * t1
        t1 = one
      END IF
      mb = mb + 1
    END DO
!--------------------------------------------------------------------
!   Accurate Summation
!--------------------------------------------------------------------
    xmb = xmb + one
    20     sum = u2(mb+1)
    ind = mb
    DO  ii = 1, mb
      sum = sum * d / xmb + u2(ind)
      ind = ind - 1
      xmb = xmb - one
    END DO
    zz = sum
    IF (tflag) zz = zz / eps
    zz = zz * xlam ** alpha
    mb = 2
    CALL rkbesl(x,alpha,mb,ize,u,ncalc)
    z = u(1)
    y = z
    IF (u2(1) > y) y = u2(1)
!--------------------------------------------------------------------
!   Accumulate Results
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
      a1 = alpha
      iz1 = ize
    END IF
    r7 = r7 + w * w
    xl = xl + del
  END DO loop40
!--------------------------------------------------------------------
!   Gather and print statistics for test
!--------------------------------------------------------------------
  n = k1 + k2 + k3
  r7 = SQRT(r7/xn)
  WRITE (iout,5000)
  WRITE (iout,5100) n, a, b
  WRITE (iout,5200) k1, k2, k3
  WRITE (iout,5300) it, ibeta
  w = x99
  IF (r6 /= zero) w = LOG(r6) / albeta
  IF (j == 3) THEN
    WRITE (iout,5700) r6, ibeta, w, x1, a1, iz1
  ELSE
    WRITE (iout,5400) r6, ibeta, w, x1, a1, iz1
  END IF
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
  w = x99
  IF (r7 /= zero) w = LOG(r7) / albeta
  IF (j == 3) THEN
    WRITE (iout,5800) r7, ibeta, w
  ELSE
    WRITE (iout,5600) r7, ibeta, w
  END IF
  w = MAX(ait+w,zero)
  WRITE (iout,5500) ibeta, w
!--------------------------------------------------------------------
!   Initialize for next test
!--------------------------------------------------------------------
  a = b
  b = b + b
  IF (j == 1) b = ten
END DO
!-------------------------------------------------------------------
!   Test of error returns
!   First, test with bad parameters
!-------------------------------------------------------------------
WRITE (iout,5900)
x = -one
alpha = half
mb = 5
ize = 2
u(1) = zero
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
x = -x
alpha = one + half
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
alpha = half
mb = -mb
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
mb = -mb
ize = 5
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
!-------------------------------------------------------------------
!   Last tests are with extreme parameters
!-------------------------------------------------------------------
x = xmin
alpha = zero
mb = 2
ize = 2
u(1) = zero
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
x = tinyx * (one-SQRT(eps))
mb = 20
u(1) = zero
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
x = tinyx * (one+SQRT(eps))
mb = 20
u(1) = zero
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
!-------------------------------------------------------------------
!   Determine largest safe argument for unscaled functions
!-------------------------------------------------------------------
z = LOG(xmin)
w = z - c
zz = -z - ten
60 z = zz
zz = one / (eight*z)
a = z + LOG(z) * half + zz * (one-xnine*half*zz) + w
b = one + (half-zz*(one-xnine*zz)) / z
zz = z - a / b
IF (ABS(z-zz) > hund*eps*z) GO TO 60
x = zz * xlam
ize = 1
mb = 2
u(1) = zero
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
x = zz
mb = 2
u(1) = zero
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
x = ten / eps
ize = 2
u(1) = zero
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
x = xmax
u(1) = zero
CALL rkbesl(x,alpha,mb,ize,u,ncalc)
WRITE (iout,6000) x, alpha, mb, ize, u(1), ncalc
STOP
!-------------------------------------------------------------------
5000 FORMAT ('1Test of K(X,ALPHA) vs Multiplication Theorem'//)
5100 FORMAT (i7,' Random arguments were tested from the interval ','(',  &
    f5.2,',',f5.2,')'//)
5200 FORMAT (' K(X,ALPHA) was larger', i6, ' times,'/ t16, ' agreed', i6,  &
    ' times, and'/ t12, 'was smaller', i6, ' times.'//)
5300 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number'//)
5400 FORMAT (' The maximum relative error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =', e13.6, ', NU =',e13.6,' and IZE =',i2)
5500 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5600 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5700 FORMAT (' The maximum absolute error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =',e13.6,', NU =',e13.6,' and IZE =',i2)
5800 FORMAT (' The root mean square absolute error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5900 FORMAT ('1Check of Error Returns'///  &
    ' The following summarizes calls with indicated parameters'//  &
    ' NCALC different from MB indicates some form of error'//  &
    ' See documentation for RKBESL for details'//   &
    t8, 'ARG', t23, 'ALPHA      MB   IZ       RES      NCALC'//)
6000 FORMAT (2E15.7, 2I5, e15.7, i5//)
!---------- Last line of RKBESL test program ----------
END PROGRAM rktest
