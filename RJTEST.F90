PROGRAM rjtest
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!----------------------------------------------------------------------
! FORTRAN 90 program to test RJBESL

!  Method:

!     Two different accuracy tests are used.  In the first interval,
!     function values are compared against values generated with the
!     multiplication formula, where the Bessel values used in the
!     multiplication formula are obtained from the function program.
!     In the remaining intervals, function values are compared
!     against values generated with a local Taylor series expansion.
!     Derivatives in the expansion are expressed in terms of the first
!     two Bessel functions, which are in turn obtained from the function
!     program.

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
!              ACM Trans. on Math. Software, Vol. 15, 1989, pp 41-48.

!             "Use of Taylor series to test accuracy of function programs,"
!              W. J. Cody and L. Stoltz, submitted for publication.

!  Latest modification: March 14, 1992

!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439, USA

! Acknowledgement: this program is a minor modification of the test
!          driver for RIBESL whose primary author was Laura Stoltz.

!----------------------------------------------------------------------

USE Toms715_Utilities
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER    :: i, ibeta, iexp, ii, iii, ind, irnd, it, j, j1,  &
              j2, k, kk, k1, k2, k3, last, m, machep, maxexp, mb, minexp,  &
              mvr, n, ncalc, negep, ngrd, nk, no1, num
REAL (dp)  :: a, ait, albeta, alpha, alphsq, a1, b, beta, d, del, delta,   &
              deriv, eps, epsneg, f, fxmx, r6, r7, sum, t1, t2, w, x, xl,  &
              xmax, xmb, xmin, xj1, xn, x1, y, ysq, z, zz
REAL (dp)  :: g(5), u(560), u2(560)
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
           two = 2.0_dp, ten = 10.0_dp, sixten = 16.0_dp, x99 = -999.0_dp, &
           xlam = 1.03125_dp, xlarge = 1.0D4, rtpi = 0.6366_dp
!----------------------------------------------------------------------
!  Arrays related to expansion of the derivatives in terms
!   of the first two Bessel functions.
!----------------------------------------------------------------------
INTEGER, PARAMETER  :: ndx(24) = (/ 9, 7, 5, 3, 1, 8, 6, 4, 2, 7, 5, 3,  &
                                    1, 6, 4, 2, 5, 3, 1, 4, 2, 3, 1, 2 /)
INTEGER, PARAMETER  :: ndx2(8) = (/ 5, 9, 13, 16, 19, 21, 23, 24 /)
REAL (dp), PARAMETER  :: AR1(11, 6) = RESHAPE(  &
  (/ 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, -3.0_dp, 0.0_dp, -2.0_dp, &
     1.2D1, 0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, -1.0_dp, 2.0_dp, 0.0_dp,           &
     2.0_dp, -6.0_dp, 1.0_dp, -7.0_dp, 2.4D1, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp,   &
     -3.0_dp, 0.0_dp, -2.0_dp, 1.1D1, 0.0_dp, 1.2D1, -5.0D1, 0.0_dp,            &
     0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, -6.0_dp, 0.0_dp, -2.0_dp,  &
     3.5D1, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp,     &
     0.0_dp, 0.0_dp, -1.0D1, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,    &
     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0D0 /), (/ 11, 6 /) )
REAL (dp), PARAMETER  :: AR2(13, 9) = RESHAPE(  &
  (/ -1.0_dp, 9.0_dp, -6.0D1, 0.0_dp, 3.0_dp, -5.1D1, 3.6D2, 0.0_dp,           &
      1.0_dp, -1.8D1, 3.45D2, -2.52D3, 0.0_dp, 0.0_dp, -3.0_dp, 3.3D1,         &
      -1.2D2, -1.0_dp, 1.5D1, -1.92D2, 7.2D2, 0.0_dp, 4.0_dp, -9.6D1,          &
      1.32D3, -5.04D3, 0.0_dp, 3.0_dp, -7.8D1, 2.74D2, 0.0_dp, -2.7D1,         &
      5.7D2, -1.764D3, 0.0_dp, -4.0_dp, 2.46D2, -4.666D3, 1.3068D4,            &
      0.0_dp, 0.0_dp, 1.8D1, -2.25D2, 0.0_dp, 3.0_dp, -1.5D2, 1.624D3,         &
      0.0_dp, 0.0_dp, -3.6D1, 1.32D3, -1.3132D4, 0.0_dp, 0.0_dp, -3.0_dp,      &
      8.5D1, 0.0_dp, 0.0_dp, 4.5D1, -7.35D2, 0.0_dp, 0.0_dp, 6.0_dp,           &
      -5.5D2, 6.769D3, 0.0_dp, 0.0_dp, 0.0_dp, -1.5D1, 0.0_dp, 0.0_dp,         &
      -3.0_dp, 1.75D2, 0.0_dp, 0.0_dp, 0.0_dp, 6.0D1, -1.96D3, 0.0_dp,         &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -2.1D1, 0.0_dp, 0.0_dp,  &
      0.0_dp, -4.0_dp, 3.22D2, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -2.8D1, 0.0_dp, 0.0_dp,  &
      0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
      0.0_dp, 1.0D0 /), (/ 13, 9 /) )
INTEGER, PARAMETER  :: iout = 6

INTERFACE
  SUBROUTINE rjbesl(x, alpha, nb, b, ncalc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)   :: x
    REAL (dp), INTENT(IN)   :: alpha
    INTEGER, INTENT(IN)     :: nb
    REAL (dp), INTENT(OUT)  :: b(:)
    INTEGER, INTENT(OUT)    :: ncalc
  END SUBROUTINE rjbesl
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
    x = del * ren() + xl
    10 alpha = ren()
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
      mb = 11
    ELSE
      x = y + delta
      mb = 2
    END IF
    CALL rjbesl(y, alpha, mb, u2, ncalc)
    CALL rjbesl(x, alpha, mb, u, ncalc)
    IF (j == 1) THEN
!----------------------------------------------------------------------
!   Accuracy test is based on the multiplication theorem
!----------------------------------------------------------------------
      d = -f * y
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
      IF (ABS(u(1)) < ABS(u2(1))) THEN
        z = x
        x = y
        y = z
        delta = x - y
        DO  ii = 1, mb
          z = u(ii)
          u(ii) = u2(ii)
          u2(ii) = z
        END DO
      END IF
!----------------------------------------------------------------------
!   Filter out cases where function values or derivatives are small
!----------------------------------------------------------------------
      w = SQRT(rtpi/x) / sixten
      z = MIN(ABS(u2(1)),ABS(u2(2)))
      IF (z < w) GO TO 10
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
!         Group J(ALPHA) terms in the derivative
!----------------------------------------------------------------------
        DO  iii = 1, k
          j2 = ndx(mvr)
          IF (num == 1) THEN
            g(iii) = ar1(nk,j2)
          ELSE
            g(iii) = ar2(nk,j2)
          END IF
          IF (j2 > 1) THEN
            40             j2 = j2 - 1
            IF (num == 1) THEN
              g(iii) = g(iii) * alpha + ar1(nk,j2)
            ELSE
              g(iii) = g(iii) * alpha + ar2(nk,j2)
            END IF
            IF (j2 > 1) GO TO 40
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
!         Group J(ALPHA+1) terms in the derivative
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
            70             j2 = j2 - 2
            IF (j2 >= 2) THEN
              IF (num == 1) THEN
                g(iii) = g(iii) * alphsq + ar1(nk,j2)
              ELSE
                g(iii) = g(iii) * alphsq + ar2(nk,j2)
              END IF
              GO TO 70
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
        deriv = u2(1) * t1 - u2(2) * t2
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
      fxmx = z
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
  WRITE (iout,5800) fxmx
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
  IF (j == 1) THEN
    b = ten
  ELSE IF (j == 2) THEN
    b = b + b
  ELSE IF (j == 3) THEN
    a = a + ten
    b = a + ten
  END IF
END DO
!----------------------------------------------------------------------
!   Test of error returns

!   First, test with bad parameters
!----------------------------------------------------------------------
WRITE (iout,5900)
x = one
alpha = one + half
mb = 5
CALL rjbesl(x, alpha, mb, u, ncalc)
WRITE (iout,6000) x, alpha, mb, u(1), ncalc
alpha = half
mb = -mb
CALL rjbesl(x, alpha, mb, u, ncalc)
WRITE (iout,6000) x, alpha, mb, u(1), ncalc
!----------------------------------------------------------------------
!   Last tests are with extreme parameters
!----------------------------------------------------------------------
x = zero
alpha = one
mb = 2
CALL rjbesl(x, alpha, mb, u, ncalc)
WRITE (iout,6000) x, alpha, mb, u(1), ncalc
x = -one
alpha = half
mb = 5
CALL rjbesl(x, alpha, mb, u, ncalc)
WRITE (iout,6000) x, alpha, mb, u(1), ncalc
!----------------------------------------------------------------------
!   Determine largest safe argument
!----------------------------------------------------------------------
WRITE (iout,6400)
x = xlarge * (one-SQRT(SQRT(eps)))
mb = 2
CALL rjbesl(x,alpha,mb,u,ncalc)
WRITE (iout,6100) x
WRITE (iout,6300) ncalc, u(1)
x = xlarge * (one+SQRT(SQRT(eps)))
mb = 2
CALL rjbesl(x,alpha,mb,u,ncalc)
WRITE (iout,6100) x
WRITE (iout,6200)
WRITE (iout,6300) ncalc, u(1)
WRITE (iout,6500)
STOP
!----------------------------------------------------------------------
5000 FORMAT ('1Test of J(X,ALPHA) vs Multiplication Theorem'//)
5100 FORMAT ('1Test of J(X,ALPHA) vs Taylor series'//)
5200 FORMAT (i7,' Random arguments were tested from the interval ','(',  &
    f5.2,',',f5.2,')'//)
5300 FORMAT (' J(X,ALPHA) was larger',i6,' times,'/ t16, ' agreed',i6,  &
    ' times, and'/ t12, 'was smaller',i6,' times.'//)
5400 FORMAT (' There are',i4,' base',i4,  &
    ' significant digits in a floating-point number'//)
5500 FORMAT (' The maximum relative error of',e15.4,' = ',i4,' **',f7.2  &
    / t5, 'occurred for X =',e13.6,' and NU =',e13.6)
5600 FORMAT (' The estimated loss of base',i4,' significant digits is',  &
    f7.2//)
5700 FORMAT (' The root mean square relative error was',e15.4,' = ',i4,  &
    ' **',f7.2)
5800 FORMAT ( t5, 'with J(X,ALPHA) = ',e13.6)
5900 FORMAT ('1Check of Error Returns'///  &
    ' The following summarizes calls with indicated parameters'//  &
    ' NCALC different from MB indicates some form of error'//  &
    ' See documentation for RJBESL for details'// t8, 'ARG', t23,  &
    'ALPHA      MB      B(1)      NCALC'//)
6000 FORMAT (2E15.7,i5,e15.7,i5//)
6100 FORMAT (' RJBESL will be called with the argument',e13.6)
6200 FORMAT (' This should trigger an error message.')
6300 FORMAT (' NCALC returned the value',i5/  &
    ' and RJBESL returned U(1) = ',e13.6/)
6400 FORMAT (' Tests near the largest acceptable argument for RJBESL'/)
6500 FORMAT (' This concludes the tests.')
!     ---------- Last line of RJBESL test program ----------
END PROGRAM rjtest
