PROGRAM ritest
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!----------------------------------------------------------------------
! FORTRAN 90 program to test RIBESL

!  Method:

!     Two different accuracy tests are used.  In the first interval,
!     function values are compared against values generated with the
!     multiplication formula, where the Bessel values used in the
!     multiplication formula are obtained from the function program.
!     In the remaining intervals, function values are compared against
!     values generated with a local Taylor series expansion.
!     Derivatives in the expansion are expressed in terms of the first two
!     Bessel functions, which are in turn obtained from the function program.

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
!              XMAX   - the largest finite floating-point number

!     REN(K) - a function subprogram returning random real
!              numbers uniformly distributed over (0,1)

!  Intrinsic functions required are:

!      ABS, DBLE, INT, LOG, MAX, REAL, SQRT

!  Reference: "Performance evaluation of programs for certain
!              Bessel functions", W. J. Cody and L. Stoltz,
!              ACM Trans. on Math. Software, Vol. 15, 1989,
!              pp 41-48.

!             "Use of Taylor series to test accuracy of function
!              programs," W. J. Cody and L. Stoltz, submitted for
!              publication.

!  Latest modification: March 14, 1992

!  Authors: W. J. Cody and L. Stoltz
!           Mathematics and Computer Science Division
!           Argonne National Laboratory
!           Argonne, IL 60439

!----------------------------------------------------------------------
USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, ii, iii, ind, irnd, it, ize, j, j1, j2, k, kk, &
              k1, k2, k3, last, m, machep, maxexp, mb, mborg, minexp, mvr,   &
              n, ncalc, negep, ngrd, nk, no1, num
REAL (dp)  :: a, ait, ak, akk, albeta, alpha, alphsq, a1, b, beta, d, del,  &
              delta, deriv, e, eps, epsneg, f, r6, r7, sum, t1, t2,   &
              w, x, xbad, xl, xmax, xmb, xmin, xj1, xn, x1, y, ysq, z, zz
REAL (dp)  :: g(5), u(560), u2(560)
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, ten = 10.0_dp, sixten = 16.0_dp,  &
                         hund = 100.0_dp, x99 = -999.0_dp, xlam = 1.03125_dp, &
                         xlarge = 1.0D4, c = 0.9189385332_dp
!----------------------------------------------------------------------
!  Arrays related to expansion of the derivatives in terms
!   of the first two Bessel functions.
!----------------------------------------------------------------------
INTEGER, PARAMETER    :: ndx(24) = (/ 9, 7, 5, 3, 1, 8, 6, 4, 2, 7, 5, 3, 1,  &
                                    6, 4, 2, 5, 3, 1, 4, 2, 3, 1, 2 /)
INTEGER, PARAMETER    :: ndx2(8) = (/ 5, 9, 13, 16, 19, 21, 23, 24 /)
REAL (dp), PARAMETER  :: ar1(11, 6) = RESHAPE(  &
           (/ 0.0D0, 1.0D0, 0.0D0, -1.0D0, 0.0D0, 1.0D0, 3.0D0, 0.0D0, -2.0D0, &
             -1.2D1, 0.0D0, 1.0D0, 0.0D0, -1.0D0, 1.0D0, 2.0D0, 0.0D0,         &
             -2.0D0, -6.0D0, 1.0D0, 7.0D0, 2.4D1, 0.0D0, 0.0D0, 1.0D0, 0.0D0,  &
             -3.0D0, 0.0D0, 2.0D0, 1.1D1, 0.0D0, -1.2D1, -5.0D1, 0.0D0,        &
              0.0D0, 0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, -6.0D0, 0.0D0, 2.0D0,  &
              3.5D1, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0,   &
              0.0D0, 0.0D0, -1.0D1, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,  &
              0.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0 /), (/ 11, 6 /) )
REAL (dp), PARAMETER  :: ar2(13, 9) = RESHAPE(  &
           (/ 1.0D0, 9.0D0, 6.0D1, 0.0D0, -3.0D0, -5.1D1, -3.6D2, 0.0D0,   &
              1.0D0, 1.8D1, 3.45D2, 2.52D3, 0.0D0, 0.0D0, -3.0D0, -3.3D1,  &
             -1.2D2, 1.0D0, 1.5D1, 1.92D2, 7.2D2, 0.0D0, -4.0D0, -9.6D1,   &
             -1.32D3, -5.04D3, 0.0D0, 3.0D0, 7.8D1, 2.74D2, 0.0D0, -2.7D1, &
             -5.7D2, -1.764D3, 0.0D0, 4.0D0, 2.46D2, 4.666D3, 1.3068D4,    &
              0.0D0, 0.0D0, -1.8D1, -2.25D2, 0.0D0, 3.0D0, 1.5D2, 1.624D3, &
              0.0D0, 0.0D0, -3.6D1, -1.32D3, -1.3132D4, 0.0D0, 0.0D0, 3.0D0,  &
              8.5D1, 0.0D0, 0.0D0, -4.5D1, -7.35D2, 0.0D0, 0.0D0, 6.0D0,   &
              5.5D2, 6.769D3, 0.0D0, 0.0D0, 0.0D0, -1.5D1, 0.0D0, 0.0D0,   &
              3.0D0, 1.75D2, 0.0D0, 0.0D0, 0.0D0, -6.0D1, -1.96D3, 0.0D0,  &
              0.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, -2.1D1, 0.0D0, 0.0D0, &
              0.0D0, 4.0D0, 3.22D2, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
              0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, -2.8D1, 0.0D0, 0.0D0, &
              0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,  &
              0.0D0, 1.0D0 /), (/ 13, 9 /) )
INTEGER, PARAMETER  :: iout = 6

INTERFACE
  SUBROUTINE ribesl(x, alpha, nb, ize, b, ncalc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x
    REAL (dp), INTENT(IN)   :: alpha
    INTEGER, INTENT(IN)     :: nb
    INTEGER, INTENT(IN)     :: ize
    REAL (dp), INTENT(OUT)  :: b(:)
    INTEGER, INTENT(OUT)    :: ncalc
  END SUBROUTINE ribesl
END INTERFACE

!----------------------------------------------------------------------
!  Determine machine parameters and set constants
!----------------------------------------------------------------------
CALL machar(ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp,  &
            eps, epsneg, xmin, xmax)
beta = ibeta
ait = it
albeta = LOG(beta)
a = zero
b = two
delta = xlam - one
f = (delta) * (xlam+one) * half
!----------------------------------------------------------------------
!  Random argument accuracy tests
!----------------------------------------------------------------------
DO  j = 1, 4
!----------------------------------------------------------------------
!  Determine the number of terms needed for convergence of the series
!  used in the multiplication theorem.  Use Newton iteration on the
!  asymptotic form of the convergence check for I0(X).
!----------------------------------------------------------------------
  xbad = b
  d = ait * albeta - c + one
  e = LOG(xbad*f) + one
  akk = one
  10   ak = akk
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
  a1 = zero
  r6 = zero
  r7 = zero
  del = (b-a) / xn
  xl = a
  DO  i = 1, n
    mb = mborg
    x = del * ren() + xl
    alpha = ren()
    ize = 1
!----------------------------------------------------------------------
!   Carefully purify arguments
!----------------------------------------------------------------------
    IF (j == 1) THEN
      y = x / xlam
    ELSE
      y = x - delta
    END IF
    w = sixten * y
    t1 = w + y
    t1 = w + t1
    y = t1 - w
    y = y - w
    IF (j == 1) THEN
      x = y * xlam
    ELSE
      x = y + delta
    END IF
    CALL ribesl(y, alpha, mb, ize, u2, ncalc)
    IF (j == 1) THEN
!----------------------------------------------------------------------
!   Accuracy test is based on the multiplication theorem
!----------------------------------------------------------------------
      d = f * y
      mb = ncalc - 2
      xmb = mb
      sum = u2(mb+1)
      ind = mb
      DO  ii = 2, mb
        sum = sum * d / xmb + u2(ind)
        ind = ind - 1
        xmb = xmb - one
      END DO
      zz = sum * d + u2(ind)
      zz = zz * xlam ** alpha
    ELSE
!----------------------------------------------------------------------
!   Accuracy test is based on local Taylor's series expansion
!----------------------------------------------------------------------
      ysq = y * y
      alphsq = alpha * alpha
      mb = 8
      j1 = mb
      xj1 = j1+1
      iexp = 0
      nk = 13
      num = 2
      DO  ii = 1, mb
        IF (nk == 0) THEN
          nk = 11
          num = 1
        END IF
        k = 9 - j1
        IF (k > 1) THEN
          no1 = ndx2(k-1) + 1
        ELSE
          no1 = 1
        END IF
        mvr = no1
        last = ndx2(k)
        k = last - no1 + 1
!----------------------------------------------------------------------
!         Group I(ALPHA) terms in the derivative
!----------------------------------------------------------------------
        DO  iii = 1, k
          j2 = ndx(mvr)
          IF (num == 1) THEN
            g(iii) = ar1(nk,j2)
          ELSE
            g(iii) = ar2(nk,j2)
          END IF
          IF (j2 > 1) THEN
            30             j2 = j2 - 1
            IF (num == 1) THEN
              g(iii) = g(iii) * alpha + ar1(nk,j2)
            ELSE
              g(iii) = g(iii) * alpha + ar2(nk,j2)
            END IF
            IF (j2 > 1) GO TO 30
          END IF
          mvr = mvr + 1
          nk = nk - 1
        END DO
        t1 = g(1)
        DO  iii = 2, k
          t1 = t1 / ysq + g(iii)
        END DO
        IF (iexp == 1) t1 = t1 / y
!----------------------------------------------------------------------
!         Group I(ALPHA+1) terms in the derivative
!----------------------------------------------------------------------
        iexp = 1 - iexp
        nk = nk + k
        mvr = no1
        kk = k
        DO  iii = 1, k
          j2 = ndx(mvr)
          m = MOD(j2,2)
          IF (m == 1) j2 = j2 - 1
          IF (j2 >= 2) THEN
            IF (num == 1) THEN
              g(iii) = ar1(nk,j2)
            ELSE
              g(iii) = ar2(nk,j2)
            END IF
            60             j2 = j2 - 2
            IF (j2 >= 2) THEN
              IF (num == 1) THEN
                g(iii) = g(iii) * alphsq + ar1(nk,j2)
              ELSE
                g(iii) = g(iii) * alphsq + ar2(nk,j2)
              END IF
              GO TO 60
            END IF
          ELSE
            kk = iii - 1
          END IF
          mvr = mvr + 1
          nk = nk - 1
        END DO
        t2 = g(1)
        DO  iii = 2, kk
          t2 = t2 / ysq + g(iii)
        END DO
        IF (iexp == 1) t2 = t2 / y
        deriv = u2(1) * t1 + u2(2) * t2
        IF (j1 == 8) THEN
          sum = deriv
        ELSE
          sum = sum * delta / xj1 + deriv
        END IF
        j1 = j1 - 1
        xj1 = xj1 - one
      END DO
      zz = sum * delta + u2(1)
    END IF
    mb = 2
    CALL ribesl(x, alpha, mb, ize, u, ncalc)
    z = u(1)
!----------------------------------------------------------------------
!   Accumulate Results
!----------------------------------------------------------------------
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
      a1 = alpha
    END IF
    r7 = r7 + w * w
    xl = xl + del
  END DO
!----------------------------------------------------------------------
!   Gather and print statistics for test
!----------------------------------------------------------------------
  k2 = n - k1 - k3
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
  WRITE (iout,5500) r6, ibeta, w, x1, a1
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
!----------------------------------------------------------------------
!   Initialize for next test
!----------------------------------------------------------------------
  a = b
  b = b + b
  IF (j == 2) b = ten
END DO
!----------------------------------------------------------------------
!   Test of error returns

!   First, test with bad parameters
!----------------------------------------------------------------------
WRITE (iout,5800)
x = one
alpha = one + half
mb = 5
ize = 2
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,5900) x, alpha, mb, ize, u(1), ncalc
alpha = half
mb = -mb
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,5900) x, alpha, mb, ize, u(1), ncalc
mb = -mb
ize = 5
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,5900) x, alpha, mb, ize, u(1), ncalc
!----------------------------------------------------------------------
!   Last tests are with extreme parameters
!----------------------------------------------------------------------
x = zero
alpha = ren()
mb = 2
ize = 1
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,5900) x, alpha, mb, ize, u(1), ncalc
alpha = zero
mb = 2
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,5900) x, alpha, mb, ize, u(1), ncalc
alpha = one
mb = 2
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,5900) x, alpha, mb, ize, u(1), ncalc
x = -one
alpha = half
mb = 5
ize = 2
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,6000) x
WRITE (iout,6100)
WRITE (iout,6200) ncalc, u(1)
!----------------------------------------------------------------------
!   Determine largest safe argument for scaled functions
!----------------------------------------------------------------------
WRITE (iout,6300)
x = xlarge * (one-SQRT(SQRT(eps)))
ize = 2
mb = 2
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,6000) x
WRITE (iout,6200) ncalc, u(1)
x = xlarge * (one+SQRT(SQRT(eps)))
mb = 2
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,6000) x
WRITE (iout,6100)
WRITE (iout,6200) ncalc, u(1)
!----------------------------------------------------------------------
!   Determine largest safe argument for unscaled functions
!----------------------------------------------------------------------
WRITE (iout,6400)
n = LOG(xmax)
z = n
x = z * (one - SQRT(SQRT(eps)))
ize = 1
mb = 2
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,6000) x
WRITE (iout,6200) ncalc, u(1)
x = z * (one + SQRT(SQRT(eps)))
mb = 2
u(1) = zero
CALL ribesl(x, alpha, mb, ize, u, ncalc)
WRITE (iout,6000) x
WRITE (iout,6100)
WRITE (iout,6200) ncalc, u(1)
WRITE (iout,6500)
STOP
!----------------------------------------------------------------------
5000 FORMAT ('1Test of I(X,ALPHA) vs Multiplication Theorem'//)
5100 FORMAT ('1Test of I(X,ALPHA) vs Taylor series'//)
5200 FORMAT (i7,' Random arguments were tested from the interval (',  &
    f5.2, ',', f5.2,')'//)
5300 FORMAT (' I(X,ALPHA) was larger',i6,' times,'/ t16, ' agreed', i6,  &
    ' times, and'/ t12, 'was smaller', i6,' times.'//)
5400 FORMAT (' There are', i4, ' base', i4,  &
    ' significant digits in a floating-point number'//)
5500 FORMAT (' The maximum relative error of', e15.4, ' = ', i4, ' **',f7.2  &
    / t5, 'occurred for X =', e13.6,' and NU =', e13.6)
5600 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5700 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5800 FORMAT ('1Check of Error Returns'///  &
    ' The following summarizes calls with indicated parameters'//  &
    ' NCALC different from MB indicates some form of error'//  &
    ' See documentation for RIBESL for details'// t8, 'ARG', t23,  &
    'ALPHA      MB   IZ       RES      NCALC'//)
5900 FORMAT (2E15.7, 2I5, e15.7, i5//)
6000 FORMAT (' RIBESL will be called with the argument',e13.6)
6100 FORMAT (' This should trigger an error message.')
6200 FORMAT (' NCALC returned the value',i5/  &
    ' and RIBESL returned the value',e13.6/)
6300 FORMAT (' Tests near the largest argument for scaled functions'/)
6400 FORMAT (' Tests near the largest argument for unscaled functions'/ )
6500 FORMAT (' This concludes the tests.')
!     ---------- Last line of RIBESL test program ----------
END PROGRAM ritest
