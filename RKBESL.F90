SUBROUTINE rkbesl(x, alpha, nb, ize, bk, ncalc)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!-------------------------------------------------------------------

!  This FORTRAN 90 routine calculates modified Bessel functions of the
!  second kind, K SUB(N+ALPHA) (X), for non-negative argument X, and
!  non-negative order N+ALPHA, with or without exponential scaling.

!  Explanation of variables in the calling sequence

!  Description of output values ..

! X     - Working precision non-negative real argument for which
!         K's or exponentially scaled K's (K*EXP(X))
!         are to be calculated.  If K's are to be calculated,
!         X must not be greater than XMAX (see below).
! ALPHA - Working precision fractional part of order for which
!         K's or exponentially scaled K's (K*EXP(X)) are
!         to be calculated.  0 <= ALPHA < 1.0.
! NB    - Integer number of functions to be calculated, NB > 0.
!         The first function calculated is of order ALPHA, and the
!         last is of order (NB - 1 + ALPHA).
! IZE   - Integer type.  IZE = 1 if unscaled K's are to be calculated,
!         and 2 if exponentially scaled K's are to be calculated.
! BK    - Working precision output vector of length NB.  If the
!         routine terminates normally (NCALC=NB), the vector BK
!         contains the functions K(ALPHA,X), ... , K(NB-1+ALPHA,X),
!         or the corresponding exponentially scaled functions.
!         If (0 < NCALC < NB), BK(I) contains correct function
!         values for I <= NCALC, and contains the ratios
!         K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.
! NCALC - Integer output variable indicating possible errors.
!         Before using the vector BK, the user should check that
!         NCALC=NB, i.e., all orders have been calculated to
!         the desired accuracy.  See error returns below.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   EPS    = The smallest positive floating-point number such that
!            1.0+EPS > 1.0
!   XMAX   = Upper limit on the magnitude of X when IZE=1;  Solution
!            to equation:
!               W(X) * (1-1/8X+9/128X**2) = beta**minexp
!            where  W(X) = EXP(-X)*SQRT(PI/2X)
!   SQXMIN = Square root of beta**minexp
!   XINF   = Largest positive machine number; approximately
!            beta**maxexp
!   XMIN   = Smallest positive machine number; approximately beta**minexp

!     Approximate values for some important machines are:

!                          beta       minexp      maxexp      EPS

!  CRAY-1        (S.P.)      2        -8193        8191    7.11E-15
!  Cyber 180/185
!    under NOS   (S.P.)      2         -975        1070    3.55E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)      2         -126         128    1.19E-7
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)      2        -1022        1024    2.22D-16
!  IBM 3033      (D.P.)     16          -65          63    2.22D-16
!  VAX           (S.P.)      2         -128         127    5.96E-8
!  VAX D-Format  (D.P.)      2         -128         127    1.39D-17
!  VAX G-Format  (D.P.)      2        -1024        1023    1.11D-16


!                         SQXMIN       XINF        XMIN      XMAX

! CRAY-1        (S.P.)  6.77E-1234  5.45E+2465  4.59E-2467 5674.858
! Cyber 180/855
!   under NOS   (S.P.)  1.77E-147   1.26E+322   3.14E-294   672.788
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)  1.08E-19    3.40E+38    1.18E-38     85.337
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)  1.49D-154   1.79D+308   2.23D-308   705.342
! IBM 3033      (D.P.)  7.35D-40    7.23D+75    5.40D-79    177.852
! VAX           (S.P.)  5.42E-20    1.70E+38    2.94E-39     86.715
! VAX D-Format  (D.P.)  5.42D-20    1.70D+38    2.94D-39     86.715
! VAX G-Format  (D.P.)  7.46D-155   8.98D+307   5.57D-309   706.728

!*******************************************************************

! Error returns

!  In case of an error, NCALC .NE. NB, and not all K's are
!  calculated to the desired accuracy.

!  NCALC < -1:  An argument is out of range. For example,
!       NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) >= XMAX.
!       In this case, the B-vector is not calculated,
!       and NCALC is set to MIN0(NB,0)-2  so that NCALC .NE. NB.
!  NCALC = -1:  Either  K(ALPHA,X) >= XINF  or
!       K(ALPHA+NB-1,X)/K(ALPHA+NB-2,X) >= XINF.  In this case,
!       the B-vector is not calculated.  Note that again NCALC .NE. NB.

!  0 < NCALC < NB: Not all requested function values could
!       be calculated accurately.  BK(I) contains correct function
!       values for I <= NCALC, and contains the ratios
!       K(ALPHA+I-1,X)/K(ALPHA+I-2,X) for the rest of the array.


! Intrinsic functions required are:

!     ABS, AINT, EXP, INT, LOG, MAX, MIN, SINH, SQRT


! Acknowledgement

!  This program is based on a program written by J. B. Campbell (2) that
!  computes values of the Bessel functions K of real argument and real order.
!  Modifications include the addition of non-scaled functions,
!  parameterization of machine dependencies, and the use of more accurate
!  approximations for SINH and SIN.

! References: "On Temme's Algorithm for the Modified Bessel Functions of the
!              Third Kind," Campbell, J. B., TOMS 6(4), Dec. 1980, pp. 581-586.

!             "A FORTRAN IV Subroutine for the Modified Bessel Functions of
!              the Third Kind of Real Order and Real Argument"
!              Campbell, J. B., Report NRC/ERB-925, National Research Council,
!              Canada.

!  Latest modification: March 15, 1992

!  Modified by: W. J. Cody and L. Stoltz
!               Applied Mathematics Division
!               Argonne National Laboratory
!               Argonne, IL  60439, USA

!-------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: alpha
INTEGER, INTENT(IN)     :: nb
INTEGER, INTENT(IN)     :: ize
REAL (dp), INTENT(OUT)  :: bk(:)
INTEGER, INTENT(OUT)    :: ncalc

INTEGER    :: i, iend, itemp, j, k, m, mplus1
REAL (dp)  :: blpha, bk1, bk2, c, dm, d1, d2, d3, enu, ex, f0,  &
              f1, f2, p0, q0, ratio, twonu, twox, t1, t2, wminf, x2by4
!---------------------------------------------------------------------
!  Mathematical constants
!    A = LOG(2.D0) - Euler's constant
!    D = SQRT(2.D0/PI)
!---------------------------------------------------------------------
REAL (dp), PARAMETER  :: half = 0.5_dp, one = 1.0_dp, two = 2.0_dp,   &
           zero = 0.0_dp, four = 4.0_dp, tinyx = 1.0D-10,  &
           a = 0.11593151565841244881_dp, d = 0.797884560802865364_dp
!---------------------------------------------------------------------
!  Machine dependent parameters
!---------------------------------------------------------------------
REAL (dp), PARAMETER  :: eps = 2.22D-16, sqxmin = 1.49D-154, xinf = 1.79D+308, &
                         xmin = 2.23D-308, xmax = 705.342_dp
!---------------------------------------------------------------------
!  P, Q - Approximation for LOG(GAMMA(1+ALPHA))/ALPHA
!                                         + Euler's constant
!         Coefficients converted from hex to decimal and modified
!         by W. J. Cody, 2/26/82
!  R, S - Approximation for (1-ALPHA*PI/SIN(ALPHA*PI))/(2.D0*ALPHA)
!  T    - Approximation for SINH(Y)/Y
!---------------------------------------------------------------------
REAL (dp), PARAMETER  :: P(8) = (/   &
           0.805629875690432845D00, 0.204045500205365151D02,  &
           0.157705605106676174D03, 0.536671116469207504D03,  &
           0.900382759291288778D03, 0.730923886650660393D03,  &
           0.229299301509425145D03, 0.822467033424113231D00 /)
REAL (dp), PARAMETER  :: Q(7) = (/   &
           0.294601986247850434D02, 0.277577868510221208D03,  &
           0.120670325591027438D04, 0.276291444159791519D04,  &
           0.344374050506564618D04, 0.221063190113378647D04,  &
           0.572267338359892221D03 /)
REAL (dp), PARAMETER  :: R(5) = (/  &
          -0.48672575865218401848D+0, 0.13079485869097804016D+2,  &
          -0.10196490580880537526D+3, 0.34765409106507813131D+3,  &
           0.34958981245219347820D-3 /)
REAL (dp), PARAMETER  :: S(4) = (/  &
          -0.25579105509976461286D+2, 0.21257260432226544008D+3,  &
          -0.61069018684944109624D+3, 0.42269668805777760407D+3 /)
REAL (dp), PARAMETER  :: T(6) = (/   &
           0.16125990452916363814D-9, 0.25051878502858255354D-7,  &
           0.27557319615147964774D-5, 0.19841269840928373686D-3,  &
           0.83333333333334751799D-2, 0.16666666666666666446D+0 /)
REAL (dp), PARAMETER  :: ESTM(6) = (/  &
           5.20583D1, 5.7607D0, 2.7782D0, 1.44303D1, 1.853004D2, 9.3715D0 /)
REAL (dp), PARAMETER  :: ESTF(7) = (/  &
           4.18341D1, 7.1075D0, 6.4306D0, 4.25110D1, 1.35633D0,  &
           8.45096D1, 2.0D1 /)
!---------------------------------------------------------------------
ex = x
enu = alpha
ncalc = MIN(nb,0) - 2
IF (nb > 0 .AND. (enu >= zero .AND. enu < one) .AND.   &
    (ize >= 1 .AND. ize <= 2) .AND. (ize /= 1 .OR. ex <= xmax) .AND.  &
    ex > zero) THEN
  k = 0
  IF (enu < sqxmin) enu = zero
  IF (enu > half) THEN
    k = 1
    enu = enu - one
  END IF
  twonu = enu + enu
  iend = nb + k - 1
  c = enu * enu
  d3 = -c
  IF (ex <= one) THEN
!---------------------------------------------------------------------
!  Calculation of P0 = GAMMA(1+ALPHA) * (2/X)**ALPHA
!                 Q0 = GAMMA(1-ALPHA) * (X/2)**ALPHA
!---------------------------------------------------------------------
    d1 = zero
    d2 = p(1)
    t1 = one
    t2 = q(1)
    DO  i = 2, 7, 2
      d1 = c * d1 + p(i)
      d2 = c * d2 + p(i+1)
      t1 = c * t1 + q(i)
      t2 = c * t2 + q(i+1)
    END DO
    d1 = enu * d1
    t1 = enu * t1
    f1 = LOG(ex)
    f0 = a + enu * (p(8)-enu*(d1+d2)/(t1+t2)) - f1
    q0 = EXP(-enu*(a-enu*(p(8)+enu*(d1-d2)/(t1-t2))-f1))
    f1 = enu * f0
    p0 = EXP(f1)
!---------------------------------------------------------------------
!  Calculation of F0 =
!---------------------------------------------------------------------
    d1 = r(5)
    t1 = one
    DO  i = 1, 4
      d1 = c * d1 + r(i)
      t1 = c * t1 + s(i)
    END DO
    IF (ABS(f1) <= half) THEN
      f1 = f1 * f1
      d2 = zero
      DO  i = 1, 6
        d2 = f1 * d2 + t(i)
      END DO
      d2 = f0 + f0 * f1 * d2
    ELSE
      d2 = SINH(f1) / enu
    END IF
    f0 = d2 - enu * d1 / (t1*p0)
    IF (ex <= tinyx) THEN
!--------------------------------------------------------------------
!  X<=1.0E-10
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X)
!--------------------------------------------------------------------
      bk(1) = f0 + ex * f0
      IF (ize == 1) bk(1) = bk(1) - ex * bk(1)
      ratio = p0 / f0
      c = ex * xinf
      IF (k /= 0) THEN
!--------------------------------------------------------------------
!  Calculation of K(ALPHA,X) and X*K(ALPHA+1,X)/K(ALPHA,X),
!  ALPHA >= 1/2
!--------------------------------------------------------------------
        ncalc = -1
        IF (bk(1) >= c/ratio) GO TO 170
        bk(1) = ratio * bk(1) / ex
        twonu = twonu + two
        ratio = twonu
      END IF
      ncalc = 1
      IF (nb == 1) GO TO 170
!--------------------------------------------------------------------
!  Calculate  K(ALPHA+L,X)/K(ALPHA+L-1,X),  L  =  1, 2, ... , NB-1
!--------------------------------------------------------------------
      ncalc = -1
      DO  i = 2, nb
        IF (ratio >= c) GO TO 170
        bk(i) = ratio / ex
        twonu = twonu + two
        ratio = twonu
      END DO
      ncalc = 1
      GO TO 150
    ELSE
!--------------------------------------------------------------------
!  1.0E-10 < X <= 1.0
!--------------------------------------------------------------------
      c = one
      x2by4 = ex * ex / four
      p0 = half * p0
      q0 = half * q0
      d1 = -one
      d2 = zero
      bk1 = zero
      bk2 = zero
      f1 = f0
      f2 = p0
      50       d1 = d1 + two
      d2 = d2 + one
      d3 = d1 + d3
      c = x2by4 * c / d2
      f0 = (d2*f0+p0+q0) / d3
      p0 = p0 / (d2-enu)
      q0 = q0 / (d2+enu)
      t1 = c * f0
      t2 = c * (p0-d2*f0)
      bk1 = bk1 + t1
      bk2 = bk2 + t2
      IF ((ABS(t1/(f1+bk1)) > eps) .OR. (ABS(t2/(f2+bk2)) > eps)) GO TO 50
      bk1 = f1 + bk1
      bk2 = two * (f2+bk2) / ex
      IF (ize == 2) THEN
        d1 = EXP(ex)
        bk1 = bk1 * d1
        bk2 = bk2 * d1
      END IF
      wminf = estf(1) * ex + estf(2)
    END IF
  ELSE IF (eps*ex > one) THEN
!--------------------------------------------------------------------
!  X > ONE/EPS
!--------------------------------------------------------------------
    ncalc = nb
    bk1 = one / (d*SQRT(ex))
    bk(1:nb) = bk1
    GO TO 170
  ELSE
!--------------------------------------------------------------------
!  X > 1.0
!--------------------------------------------------------------------
    twox = ex + ex
    blpha = zero
    ratio = zero
    IF (ex <= four) THEN
!--------------------------------------------------------------------
!  Calculation of K(ALPHA+1,X)/K(ALPHA,X),  1.0 <= X <= 4.0
!--------------------------------------------------------------------
      d2 = INT(estm(1)/ex+estm(2))
      m = INT(d2)
      d1 = d2 + d2
      d2 = d2 - half
      d2 = d2 * d2
      DO  i = 2, m
        d1 = d1 - two
        d2 = d2 - d1
        ratio = (d3+d2) / (twox+d1-ratio)
      END DO
!--------------------------------------------------------------------
!  Calculation of I(|ALPHA|,X) and I(|ALPHA|+1,X) by backward
!    recurrence and K(ALPHA,X) from the wronskian
!--------------------------------------------------------------------
      d2 = INT(estm(3)*ex+estm(4))
      m = INT(d2)
      c = ABS(enu)
      d3 = c + c
      d1 = d3 - one
      f1 = xmin
      f0 = (two*(c+d2)/ex+half*ex/(c+d2+one)) * xmin
      DO  i = 3, m
        d2 = d2 - one
        f2 = (d3+d2+d2) * f0
        blpha = (one+d1/d2) * (f2+blpha)
        f2 = f2 / ex + f1
        f1 = f0
        f0 = f2
      END DO
      f1 = (d3+two) * f0 / ex + f1
      d1 = zero
      t1 = one
      DO  i = 1, 7
        d1 = c * d1 + p(i)
        t1 = c * t1 + q(i)
      END DO
      p0 = EXP(c*(a+c*(p(8)-c*d1/t1)-LOG(ex))) / ex
      f2 = (c+half-ratio) * f1 / ex
      bk1 = p0 + (d3*f0-f2+f0+blpha) / (f2+f1+f0) * p0
      IF (ize == 1) bk1 = bk1 * EXP(-ex)
      wminf = estf(3) * ex + estf(4)
    ELSE
!--------------------------------------------------------------------
!  Calculation of K(ALPHA,X) and K(ALPHA+1,X)/K(ALPHA,X), by backward
!  recurrence, for  X > 4.0
!--------------------------------------------------------------------
      dm = AINT(estm(5)/ex+estm(6))
      m = INT(dm)
      d2 = dm - half
      d2 = d2 * d2
      d1 = dm + dm
      DO  i = 2, m
        dm = dm - one
        d1 = d1 - two
        d2 = d2 - d1
        ratio = (d3+d2) / (twox+d1-ratio)
        blpha = (ratio+ratio*blpha) / dm
      END DO
      bk1 = one / ((d+d*blpha)*SQRT(ex))
      IF (ize == 1) bk1 = bk1 * EXP(-ex)
      wminf = estf(5) * (ex-ABS(ex-estf(7))) + estf(6)
    END IF
!--------------------------------------------------------------------
!  Calculation of K(ALPHA+1,X) from K(ALPHA,X) and
!    K(ALPHA+1,X)/K(ALPHA,X)
!--------------------------------------------------------------------
    bk2 = bk1 + bk1 * (enu+half-ratio) / ex
  END IF
!--------------------------------------------------------------------
!  Calculation of 'NCALC', K(ALPHA+I,X), I  =  0, 1, ... , NCALC-1,
!  K(ALPHA+I,X)/K(ALPHA+I-1,X), I  =  NCALC, NCALC+1, ... , NB-1
!--------------------------------------------------------------------
  ncalc = nb
  bk(1) = bk1
  IF (iend == 0) GO TO 170
  j = 2 - k
  IF (j > 0) bk(j) = bk2
  IF (iend == 1) GO TO 170
  m = MIN(INT(wminf-enu),iend)
  DO  i = 2, m
    t1 = bk1
    bk1 = bk2
    twonu = twonu + two
    IF (ex < one) THEN
      IF (bk1 >= (xinf/twonu)*ex) EXIT
      GO TO 110
    ELSE
      IF (bk1/ex >= xinf/twonu) EXIT
    END IF

    110 bk2 = twonu / ex * bk1 + t1
    itemp = i
    j = j + 1
    IF (j > 0) bk(j) = bk2
  END DO
  m = itemp
  IF (m == iend) GO TO 170
  ratio = bk2 / bk1
  mplus1 = m + 1
  ncalc = -1
  DO  i = mplus1, iend
    twonu = twonu + two
    ratio = twonu / ex + one / ratio
    j = j + 1
    IF (j > 1) THEN
      bk(j) = ratio
    ELSE
      IF (bk2 >= xinf/ratio) GO TO 170
      bk2 = ratio * bk2
    END IF
  END DO
  ncalc = MAX(mplus1-k,1)
  IF (ncalc == 1) bk(1) = bk2
  IF (nb == 1) GO TO 170
  150   j = ncalc + 1
  DO  i = j, nb
    IF (bk(ncalc) >= xinf/bk(i)) GO TO 170
    bk(i) = bk(ncalc) * bk(i)
    ncalc = i
  END DO
END IF
170 RETURN
!---------- Last line of RKBESL ----------
END SUBROUTINE rkbesl
