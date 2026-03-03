SUBROUTINE rybesl(x, alpha, nb, by, ncalc)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!----------------------------------------------------------------------

!  This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
!  for non-negative argument X, and non-negative order N+ALPHA.

! Explanation of variables in the calling sequence

! X     - Working precision positive real argument for which
!         Y's are to be calculated.
! ALPHA - Working precision fractional part of order for which
!         Y's are to be calculated.  0 <= ALPHA < 1.0.
! NB    - Integer number of functions to be calculated, NB > 0.
!         The first function calculated is of order ALPHA, and the
!         last is of order (NB - 1 + ALPHA).
! BY    - Working precision output vector of length NB.  If the
!         routine terminates normally (NCALC=NB), the vector BY
!         contains the functions Y(ALPHA,X), ... , Y(NB-1+ALPHA,X),
!         If (0 < NCALC < NB), BY(I) contains correct function
!         values for I <= NCALC, and contains the ratios
!         Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
! NCALC - Integer output variable indicating possible errors.
!         Before using the vector BY, the user should check that
!         NCALC=NB, i.e., all orders have been calculated to
!         the desired accuracy.  See error returns below.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   beta   = Radix for the floating-point system
!   p      = Number of significant base-beta digits in the
!            significand of a floating-point number
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   EPS    = beta ** (-p)
!   DEL    = Machine number below which sin(x)/x = 1; approximately SQRT(EPS).
!   XMIN   = Smallest acceptable argument for RBESY; approximately
!            max(2*beta**minexp,2/XINF), rounded up
!   XINF   = Largest positive machine number; approximately beta**maxexp
!   THRESH = Lower bound for use of the asymptotic form; approximately
!            AINT(-LOG10(EPS/2.0))+1.0
!   XLARGE = Upper bound on X; approximately 1/DEL, because the sine
!            and cosine functions have lost about half of their
!            precision at that point.


!     Approximate values for some important machines are:

!                        beta    p     minexp      maxexp      EPS

!  CRAY-1        (S.P.)    2    48     -8193        8191    3.55E-15
!  Cyber 180/185
!    under NOS   (S.P.)    2    48      -975        1070    3.55E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2    24      -126         128    5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2    53     -1022        1024    1.11D-16
!  IBM 3033      (D.P.)   16    14       -65          63    1.39D-17
!  VAX           (S.P.)    2    24      -128         127    5.96E-8
!  VAX D-Format  (D.P.)    2    56      -128         127    1.39D-17
!  VAX G-Format  (D.P.)    2    53     -1024        1023    1.11D-16


!                         DEL      XMIN      XINF     THRESH  XLARGE

! CRAY-1        (S.P.)  5.0E-8  3.67E-2466 5.45E+2465  15.0E0  2.0E7
! Cyber 180/855
!   under NOS   (S.P.)  5.0E-8  6.28E-294  1.26E+322   15.0E0  2.0E7
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)  1.0E-4  2.36E-38   3.40E+38     8.0E0  1.0E4
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)  1.0D-8  4.46D-308  1.79D+308   16.0D0  1.0D8
! IBM 3033      (D.P.)  1.0D-8  2.77D-76   7.23D+75    17.0D0  1.0D8
! VAX           (S.P.)  1.0E-4  1.18E-38   1.70E+38     8.0E0  1.0E4
! VAX D-Format  (D.P.)  1.0D-9  1.18D-38   1.70D+38    17.0D0  1.0D9
! VAX G-Format  (D.P.)  1.0D-8  2.23D-308  8.98D+307   16.0D0  1.0D8

!*******************************************************************

! Error returns

!  In case of an error, NCALC .NE. NB, and not all Y's are
!  calculated to the desired accuracy.

!  NCALC <= -1:  An argument is out of range. For example,
!       NB <= 0, or ABS(X) >= XLARGE.  In this case, BY(1) = 0.0,
!       the remainder of the BY-vector is not calculated, and NCALC is set
!       to MIN(NB,0)-1  so that NCALC .NE. NB.
!  1 < NCALC < NB: Not all requested function values could be calculated
!       accurately.  BY(I) contains correct function values for I <= NCALC,
!       and and the remaining NB-NCALC array elements contain 0.0.

! Intrinsic functions required are:

!     ABS, AINT, COS, INT, LOG, SIN, SQRT


! Acknowledgement

!  This program draws heavily on Temme's Algol program for Y(a,x)
!  and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
!  scheme is used for  x < THRESH, and Campbell's scheme is used
!  in the asymptotic region.  Segments of code from both sources
!  have been translated into Fortran 77, merged, and heavily modified.
!  Modifications include parameterization of machine dependencies,
!  use of a new approximation for ln(gamma(x)), and built-in
!  protection against overflow and destructive underflow.

! References: "Bessel functions J_nu(x) and Y_nu(x) of real
!              order and real argument," Campbell, J. B.,
!              Comp. Phy. Comm. 18, 1979, pp. 133-142.

!             "On the numerical evaluation of the ordinary
!              Bessel function of the second kind," Temme,
!              N. M., J. Comput. Phys. 21, 1976, pp. 343-350.

!  Latest modification: March 15, 1992

!  Modified by: W. J. Cody
!               Applied Mathematics Division
!               Argonne National Laboratory
!               Argonne, IL  60439, USA

!----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: alpha
INTEGER, INTENT(IN)     :: nb
REAL (dp), INTENT(OUT)  :: by(:)
INTEGER, INTENT(OUT)    :: ncalc

INTEGER    :: i, k, na
REAL (dp)  :: alfa, aye, b, c, cosmu, d, den, ddiv, div, dmu, d1, d2,  &
              e, en, enu, en1, even, ex, f, g, gamma, h, odd, p, pa, pa1,   &
              q, qa, qa1, q0, r, s, sinmu, term, twobyx, xna, x2, ya, ya1
!----------------------------------------------------------------------
!  Mathematical constants
!    FIVPI = 5*PI
!    PIM5 = 5*PI - 15
!    ONBPI = 1/PI
!    PIBY2 = PI/2
!    SQ2BPI = SQUARE ROOT OF 2/PI
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
           two = 2.0_dp, three = 3.0_dp, eight = 8.0_dp, one5 = 15.0_dp,  &
           ten9 = 19.0_dp, fivpi = 1.5707963267948966192D1,  &
           piby2 = 1.5707963267948966192_dp, pi = 3.1415926535897932385_dp,  &
           sq2bpi = 7.9788456080286535588D-1, pim5 = 7.0796326794896619231D-1, &
           onbpi = 3.1830988618379067154D-1
!----------------------------------------------------------------------
!  Machine-dependent constants
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: del = 1.0D-8, xmin = 4.46D-308, xinf = 1.79D308,  &
           eps = EPSILON(0.0_dp), thresh = 16.0_dp, xlarge = 1.0D8
!----------------------------------------------------------------------
!  Coefficients for Chebyshev polynomial expansion of
!         1/gamma(1-x), abs(x) <= .5
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: ch(21) =  &
           (/ -0.67735241822398840964D-23, -0.61455180116049879894D-22,  &
               0.29017595056104745456D-20,  0.13639417919073099464D-18,  &
               0.23826220476859635824D-17, -0.90642907957550702534D-17,  &
              -0.14943667065169001769D-14, -0.33919078305362211264D-13,  &
              -0.17023776642512729175D-12,  0.91609750938768647911D-11,  &
               0.24230957900482704055D-09,  0.17451364971382984243D-08,  &
              -0.33126119768180852711D-07, -0.86592079961391259661D-06,  &
              -0.49717367041957398581D-05,  0.76309597585908126618D-04,  &
               0.12719271366545622927D-02,  0.17063050710955562222D-02,  &
              -0.76852840844786673690D-01, -0.28387654227602353814D+00,  &
               0.92187029365045265648D+00 /)
!----------------------------------------------------------------------
ex = x
enu = alpha
IF (nb > 0 .AND. x >= xmin .AND. ex < xlarge .AND. enu >= zero  &
       .AND. enu < one) THEN
  xna = INT(enu+half)
  na = xna
  IF (na == 1) enu = enu - xna
  IF (enu == -half) THEN
    p = sq2bpi / SQRT(ex)
    ya = p * SIN(ex)
    ya1 = -p * COS(ex)
  ELSE IF (ex < three) THEN
!----------------------------------------------------------------------
!  Use Temme's scheme for small X
!----------------------------------------------------------------------
    b = ex * half
    d = -LOG(b)
    f = enu * d
    e = b ** (-enu)
    IF (ABS(enu) < del) THEN
      c = onbpi
    ELSE
      c = enu / SIN(enu*pi)
    END IF
!----------------------------------------------------------------------
!  Computation of sinh(f)/f
!----------------------------------------------------------------------
    IF (ABS(f) < one) THEN
      x2 = f * f
      en = ten9
      s = one
      DO  i = 1, 9
        s = s * x2 / en / (en-one) + one
        en = en - two
      END DO
    ELSE
      s = (e-one/e) * half / f
    END IF
!----------------------------------------------------------------------
!  Computation of 1/gamma(1-a) using Chebyshev polynomials
!----------------------------------------------------------------------
    x2 = enu * enu * eight
    aye = ch(1)
    even = zero
    alfa = ch(2)
    odd = zero
    DO  i = 3, 19, 2
      even = -(aye+aye+even)
      aye = -even * x2 - aye + ch(i)
      odd = -(alfa+alfa+odd)
      alfa = -odd * x2 - alfa + ch(i+1)
    END DO
    even = (even*half+aye) * x2 - aye + ch(21)
    odd = (odd+alfa) * two
    gamma = odd * enu + even
!----------------------------------------------------------------------
!  End of computation of 1/gamma(1-a)
!----------------------------------------------------------------------
    g = e * gamma
    e = (e+one/e) * half
    f = two * c * (odd*e+even*s*d)
    e = enu * enu
    p = g * c
    q = onbpi / g
    c = enu * piby2
    IF (ABS(c) < del) THEN
      r = one
    ELSE
      r = SIN(c) / c
    END IF
    r = pi * c * r * r
    c = one
    d = -b * b
    h = zero
    ya = f + r * q
    ya1 = p
    en = zero
    30     en = en + one
    IF (ABS(g/(one+ABS(ya)))+ABS(h/(one+ABS(ya1))) > eps) THEN
      f = (f*en+p+q) / (en*en-e)
      c = c * d / en
      p = p / (en-enu)
      q = q / (en+enu)
      g = c * (f+r*q)
      h = c * p - en * g
      ya = ya + g
      ya1 = ya1 + h
      GO TO 30
    END IF
    ya = -ya
    ya1 = -ya1 / b
  ELSE IF (ex < thresh) THEN
!----------------------------------------------------------------------
!  Use Temme's scheme for moderate X
!----------------------------------------------------------------------
    c = (half-enu) * (half+enu)
    b = ex + ex
    e = (ex*onbpi*COS(enu*pi)/eps)
    e = e * e
    p = one
    q = -ex
    r = one + ex * ex
    s = r
    en = two
    40     IF (r*en*en < e) THEN
      en1 = en + one
      d = (en-one+c/en) / s
      p = (en+en-p*d) / en1
      q = (-b+q*d) / en1
      s = p * p + q * q
      r = r * s
      en = en1
      GO TO 40
    END IF
    f = p / s
    p = f
    g = -q / s
    q = g
    50     en = en - one
    IF (en > zero) THEN
      r = en1 * (two-p) - two
      s = b + en1 * q
      d = (en-one+c/en) / (r*r+s*s)
      p = d * r
      q = d * s
      e = f + one
      f = p * e - g * q
      g = q * e + p * g
      en1 = en
      GO TO 50
    END IF
    f = one + f
    d = f * f + g * g
    pa = f / d
    qa = -g / d
    d = enu + half - p
    q = q + ex
    pa1 = (pa*q-qa*d) / ex
    qa1 = (qa*q+pa*d) / ex
    b = ex - piby2 * (enu+half)
    c = COS(b)
    s = SIN(b)
    d = sq2bpi / SQRT(ex)
    ya = d * (pa*s+qa*c)
    ya1 = d * (qa1*s-pa1*c)
  ELSE
!----------------------------------------------------------------------
!  Use Campbell's asymptotic scheme.
!----------------------------------------------------------------------
    na = 0
    d1 = AINT(ex/fivpi)
    i = INT(d1)
    dmu = ((ex-one5*d1)-d1*pim5) - (alpha+half) * piby2
    IF (i-2*(i/2) == 0) THEN
      cosmu = COS(dmu)
      sinmu = SIN(dmu)
    ELSE
      cosmu = -COS(dmu)
      sinmu = -SIN(dmu)
    END IF
    ddiv = eight * ex
    dmu = alpha
    den = SQRT(ex)
    DO  k = 1, 2
      p = cosmu
      cosmu = sinmu
      sinmu = -p
      d1 = (two*dmu-one) * (two*dmu+one)
      d2 = zero
      div = ddiv
      p = zero
      q = zero
      q0 = d1 / div
      term = q0
      DO  i = 2, 20
        d2 = d2 + eight
        d1 = d1 - d2
        div = div + ddiv
        term = -term * d1 / div
        p = p + term
        d2 = d2 + eight
        d1 = d1 - d2
        div = div + ddiv
        term = term * d1 / div
        q = q + term
        IF (ABS(term) <= eps) EXIT
      END DO
      p = p + one
      q = q + q0
      IF (k == 1) THEN
        ya = sq2bpi * (p*cosmu-q*sinmu) / den
      ELSE
        ya1 = sq2bpi * (p*cosmu-q*sinmu) / den
      END IF
      dmu = dmu + one
    END DO
  END IF
  IF (na == 1) THEN
    h = two * (enu+one) / ex
    IF (h > one) THEN
      IF (ABS(ya1) > xinf/h) THEN
        h = zero
        ya = zero
      END IF
    END IF
    h = h * ya1 - ya
    ya = ya1
    ya1 = h
  END IF
!----------------------------------------------------------------------
!  Now have first one or two Y's
!----------------------------------------------------------------------
  by(1) = ya
  by(2) = ya1
  IF (ya1 == zero) THEN
    ncalc = 1
  ELSE
    aye = one + alpha
    twobyx = two / ex
    ncalc = 2
    DO  i = 3, nb
      IF (twobyx < one) THEN
        IF (ABS(by(i-1))*twobyx >= xinf/aye) GO TO 100
      ELSE
        IF (ABS(by(i-1)) >= xinf/aye/twobyx) GO TO 100
      END IF
      by(i) = twobyx * aye * by(i-1) - by(i-2)
      aye = aye + one
      ncalc = ncalc + 1
    END DO
  END IF
  100 DO  i = ncalc + 1, nb
    by(i) = zero
  END DO
ELSE
  by(1) = zero
  ncalc = MIN(nb,0) - 1
END IF
RETURN
!---------- Last line of RYBESL ----------
END SUBROUTINE rybesl
