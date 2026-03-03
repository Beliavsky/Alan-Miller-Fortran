FUNCTION PSI(XX) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:01
 
!----------------------------------------------------------------------

! This function program evaluates the logarithmic derivative of the
!   gamma function,

!      psi(x) = d/dx (gamma(x)) / gamma(x) = d/dx (ln gamma(x))

!   for real x, where either

!          -xmax1 < x < -xmin (x not a negative integer), or
!            xmin < x.

!   The calling sequence for this function is

!                  Y = PSI(X)

!   The main computation uses rational Chebyshev approximations
!   published in Math. Comp. 27, 123-127 (1973) by Cody, Strecok and
!   Thacher.  This transportable program is patterned after the
!   machine-dependent FUNPACK program PSI(X), but cannot match that
!   version for efficiency or accuracy.  This version uses rational
!   approximations that are theoretically accurate to 20 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper selection
!   of the machine-dependent constants.

!*******************************************************************

! The following machine-dependent constants must be declared in
!   DATA statements.  IEEE values are provided as a default.

!   XINF   = largest positive machine number
!   XMAX1  = beta ** (p-1), where beta is the radix for the floating-point
!            system, and p is the number of base-beta digits in the
!            floating-point significand.  This is an upper bound on
!            non-integral floating-point numbers, and the negative of the
!            lower bound on acceptable negative arguments for PSI.
!            If rounding is necessary, round this value down.
!   XMIN1  = the smallest in magnitude acceptable argument.  We
!            recommend XMIN1 = MAX(1/XINF,xmin) rounded up, where
!            xmin is the smallest positive floating-point number.
!   XSMALL = absolute argument below which  PI*COTAN(PI*X)  may be
!            represented by 1/X.  We recommend XSMALL < sqrt(3 eps)/pi,
!            where eps is the smallest positive number such that
!            1 + eps > 1.
!   XLARGE = argument beyond which PSI(X) may be represented by
!            LOG(X).  The solution to the equation
!               x*ln(x) = beta ** p
!            is a safe value.

!     Approximate values for some important machines are

!                        beta  p     eps     xmin       XINF

!  CDC 7600      (S.P.)    2  48  7.11E-15  3.13E-294  1.26E+322
!  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
!  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
!  SUN 3/160     (D.P.)    2  53  1.11D-16  2.23D-308  1.79D+308
!  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
!                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
!   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307

!                         XMIN1      XMAX1     XSMALL    XLARGE

!  CDC 7600      (S.P.)  3.13E-294  1.40E+14  4.64E-08  9.42E+12
!  CRAY-1        (S.P.)  1.84E-2466 1.40E+14  4.64E-08  9.42E+12
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38   8.38E+06  1.90E-04  1.20E+06
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308  4.50D+15  5.80D-09  2.71D+14
!  IBM 3033      (D.P.)  1.39D-76   4.50D+15  5.80D-09  2.05D+15
!  SUN 3/160     (D.P.)  2.23D-308  4.50D+15  5.80D-09  2.71D+14
!  VAX 11/780    (S.P.)  5.89E-39   8.38E+06  1.35E-04  1.20E+06
!                (D.P.)  5.89D-39   3.60D+16  2.05D-09  2.05D+15
!   (G Format)   (D.P.)  1.12D-308  4.50D+15  5.80D-09  2.71D+14

!*******************************************************************

! Error Returns

!  The program returns XINF for  X < -XMAX1, for X zero or a negative
!    integer, or when X lies in (-XMIN1, 0), and returns -XINF
!    when X lies in (0, XMIN1).

! Intrinsic functions required are:

!     ABS, AINT, DBLE, INT, LOG, REAL, TAN


!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!  Latest modification: March 14, 1992

!----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: xx
REAL (dp)              :: fn_val

! Local variables

INTEGER    :: i, n, nq
REAL (dp)  :: aug, den, sgn, upper, w, x, z
!----------------------------------------------------------------------
!  Mathematical constants.  PIOV4 = pi / 4
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, fourth = 0.25_dp, half = 0.5_dp,  &
                         one = 1.0_dp, three = 3.0_dp, four = 4.0_dp,   &
                         piov4 = 7.8539816339744830962D-01
!----------------------------------------------------------------------
!  Machine-dependent constants
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: XINF = 1.79D+308, XMIN1 = 2.23D-308, XMAX1 = 4.50D+15, &
                         XSMALL = 5.80D-09, XLARGE = 2.71D+14
!----------------------------------------------------------------------
!  Zero of psi(x)
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: X01 = 187.0_dp, X01D = 128.0_dp,   &
                         X02 = 6.9464496836234126266D-04
!----------------------------------------------------------------------
!  Coefficients for approximation to  psi(x)/(x-x0)  over [0.5, 3.0]
!----------------------------------------------------------------------
REAL (dp), PARAMETER  ::  P1(9) =   &
        (/ 4.5104681245762934160D-03, 5.4932855833000385356D+00,  &
           3.7646693175929276856D+02, 7.9525490849151998065D+03,  &
           7.1451595818951933210D+04, 3.0655976301987365674D+05,  &
           6.3606997788964458797D+05, 5.8041312783537569993D+05,  &
           1.6585695029761022321D+05 /)
REAL (dp), PARAMETER  ::  Q1(8) =   &
        (/ 9.6141654774222358525D+01, 2.6287715790581193330D+03,  &
           2.9862497022250277920D+04, 1.6206566091533671639D+05,  &
           4.3487880712768329037D+05, 5.4256384537269993733D+05,  &
           2.4242185002017985252D+05, 6.4155223783576225996D-08 /)
!----------------------------------------------------------------------
!  Coefficients for approximation to  psi(x) - ln(x) + 1/(2x)
!     for  x > 3.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  ::  P2(7) =   &
        (/ -2.7103228277757834192D+00, -1.5166271776896121383D+01,  &
           -1.9784554148719218667D+01, -8.8100958828312219821D+00,  &
           -1.4479614616899842986D+00, -7.3689600332394549911D-02,  &
           -6.5135387732718171306D-21 /)
REAL (dp), PARAMETER  ::  Q2(6) =  &
        (/ 4.4992760373789365846D+01, 2.0240955312679931159D+02,  &
           2.4736979003315290057D+02, 1.0742543875702278326D+02,  &
           1.7463965060678569906D+01, 8.8427520398873480342D-01 /)
!----------------------------------------------------------------------

x = xx
w = ABS(x)
aug = zero
!----------------------------------------------------------------------
!  Check for valid arguments, then branch to appropriate algorithm
!----------------------------------------------------------------------
IF (-x >= xmax1 .OR. w < xmin1) THEN
  GO TO 50
ELSE IF (x >= half) THEN
  GO TO 20
!----------------------------------------------------------------------
!  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
!     Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
!----------------------------------------------------------------------
ELSE IF (w <= xsmall) THEN
  aug = -one / x
  GO TO 10
END IF
!----------------------------------------------------------------------
!  Argument reduction for cot
!----------------------------------------------------------------------
IF (x < zero) THEN
  sgn = piov4
ELSE
  sgn = -piov4
END IF
w = w - AINT(w)
nq = INT(w*four)
w = four * (w - nq*fourth)
!----------------------------------------------------------------------
!  W is now related to the fractional part of  4.0 * X.
!     Adjust argument to correspond to values in the first
!     quadrant and determine the sign.
!----------------------------------------------------------------------
n = nq / 2
IF ((n+n) /= nq) w = one - w
z = piov4 * w
IF (MOD(n,2) /= 0) sgn = -sgn
!----------------------------------------------------------------------
!  determine the final value for  -pi * cotan(pi*x)
!----------------------------------------------------------------------
n = (nq+1) / 2
IF (MOD(n,2) == 0) THEN
!----------------------------------------------------------------------
!  Check for singularity
!----------------------------------------------------------------------
  IF (z == zero) GO TO 50
  aug = sgn * (four/TAN(z))
ELSE
  aug = sgn * (four*TAN(z))
END IF
10 x = one - x
20 IF (x <= three) THEN
!----------------------------------------------------------------------
!  0.5 <= X <= 3.0
!----------------------------------------------------------------------
  den = x
  upper = p1(1) * x
  DO  i = 1, 7
    den = (den+q1(i)) * x
    upper = (upper+p1(i+1)) * x
  END DO
  den = (upper+p1(9)) / (den+q1(8))
  x = (x-x01/x01d) - x02
  fn_val = den * x + aug
  GO TO 60
END IF
!----------------------------------------------------------------------
!  3.0 < X
!----------------------------------------------------------------------
IF (x < xlarge) THEN
  w = one / (x*x)
  den = w
  upper = p2(1) * w
  DO  i = 1, 5
    den = (den+q2(i)) * w
    upper = (upper+p2(i+1)) * w
  END DO
  aug = (upper+p2(7)) / (den+q2(6)) - half / x + aug
END IF
fn_val = aug + LOG(x)
GO TO 60
!----------------------------------------------------------------------
!  Error return
!----------------------------------------------------------------------
50 fn_val = xinf
IF (x > zero) fn_val = -xinf
60 RETURN
!---------- Last card of PSI ----------
END FUNCTION PSI
