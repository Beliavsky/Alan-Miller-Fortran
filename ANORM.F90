FUNCTION ANORM(ARG) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------

! This function evaluates the normal distribution function:

!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo

!   The main computation evaluates near-minimax approximations
!   derived from those in "Rational Chebyshev approximations for
!   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
!   This transportable program uses rational functions that
!   theoretically approximate the normal distribution function to
!   at least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   XMIN  = the smallest positive floating-point number.

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   EPS   = argument below which anorm(x) may be represented by
!           0.5  and above which  x*x  will not underflow.
!           A conservative value is the largest machine number X
!           such that   1.0 + X = 1.0   to machine precision.
!   XLOW  = the most negative argument for which ANORM does not
!           vanish.  This is the negative of the solution to
!                    W(x) * (1-1/x**2) = XMIN,
!           where W(x) = exp(-x*x/2)/[x*sqrt(2*pi)].
!   XUPPR = positive argument beyond which anorm = 1.0.  A
!           conservative value is the solution to the equation
!                    exp(-x*x/2) = EPS,
!           i.e., XUPPR = sqrt[-2 ln(eps)].

!   Approximate values for some important machines are:

!                          XMIN        EPS        XLOW    XUPPR

!  CDC 7600      (S.P.)  3.13E-294   7.11E-15   -36.641   8.072
!  CRAY-1        (S.P.)  4.58E-246   7.11E-157 -106.521  26.816
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    5.96E-8    -12.949   5.768
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.11D-16   -37.519   8.572
!  IBM 195       (D.P.)  5.40D-79    1.39D-17   -18.781   8.811
!  VAX D-Format  (D.P.)  2.94D-39    1.39D-17   -13.055   8.811
!  VAX G-Format  (D.P.)  5.56D-309   1.11D-16   -37.556   8.572

!*******************************************************************

! Error returns

!  The program returns  ANORM = 0     for  ARG .LE. XLOW.


! Intrinsic functions required are:

!     ABS, AINT, EXP


!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!  Latest modification: March 15, 1992

!------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: arg
REAL (dp)              :: fn_val

! Local variables

INTEGER    :: i
REAL (dp)  :: del, x, xden, xnum, y, xsq
!------------------------------------------------------------------
!  Mathematical constants

!  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
!  THRSH is the argument for which anorm = 0.75.
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: one = 1.0_dp, half = 0.5_dp, zero = 0.0_dp,  &
                         sixten = 16.0_dp, sqrpi = 3.9894228040143267794D-1, &
                         thrsh = 0.66291_dp, root32 = 5.656854248_dp
!------------------------------------------------------------------
!  Machine-dependent constants
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: EPS = 1.11D-16, XLOW = -37.519_dp, XUPPR = 8.572_dp
!------------------------------------------------------------------
!  Coefficients for approximation in first interval
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: a(5) =  &
           (/ 2.2352520354606839287D00, 1.6102823106855587881D02,  &
              1.0676894854603709582D03, 1.8154981253343561249D04,  &
              6.5682337918207449113D-2 /)
REAL (dp), PARAMETER  :: b(4) =  &
           (/ 4.7202581904688241870D01, 9.7609855173777669322D02,  &
              1.0260932208618978205D04, 4.5507789335026729956D04 /)
!------------------------------------------------------------------
!  Coefficients for approximation in second interval
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: c(9) =  &
           (/ 3.9894151208813466764D-1, 8.8831497943883759412D00,  &
              9.3506656132177855979D01, 5.9727027639480026226D02,  &
              2.4945375852903726711D03, 6.8481904505362823326D03,  &
              1.1602651437647350124D04, 9.8427148383839780218D03,  &
              1.0765576773720192317D-8 /)
REAL (dp), PARAMETER  :: d(8) =  &
           (/ 2.2266688044328115691D01, 2.3538790178262499861D02,  &
              1.5193775994075548050D03, 6.4855582982667607550D03,  &
              1.8615571640885098091D04, 3.4900952721145977266D04,  &
              3.8912003286093271411D04, 1.9685429676859990727D04 /)
!------------------------------------------------------------------
!  Coefficients for approximation in third interval
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: p(6) =  &
           (/ 2.1589853405795699D-1, 1.274011611602473639D-1,  &
              2.2235277870649807D-2, 1.421619193227893466D-3,  &
              2.9112874951168792D-5, 2.307344176494017303D-2 /)
REAL (dp), PARAMETER  :: q(5) =  &
           (/ 1.28426009614491121D00, 4.68238212480865118D-1,  &
              6.59881378689285515D-2, 3.78239633202758244D-3,  &
              7.29751555083966205D-5 /)
!------------------------------------------------------------------

x = arg
y = ABS(x)
IF (y <= thrsh) THEN
!------------------------------------------------------------------
!  Evaluate  anorm  for  |X| <= 0.66291
!------------------------------------------------------------------
  xsq = zero
  IF (y > eps) xsq = x * x
  xnum = a(5) * xsq
  xden = xsq
  DO  i = 1, 3
    xnum = (xnum+a(i)) * xsq
    xden = (xden+b(i)) * xsq
  END DO
  fn_val = x * (xnum+a(4)) / (xden+b(4))
  fn_val = half + fn_val
!------------------------------------------------------------------
!  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
!------------------------------------------------------------------
ELSE IF (y <= root32) THEN
  xnum = c(9) * y
  xden = y
  DO  i = 1, 7
    xnum = (xnum+c(i)) * y
    xden = (xden+d(i)) * y
  END DO
  fn_val = (xnum+c(8)) / (xden+d(8))
  xsq = AINT(y*sixten) / sixten
  del = (y-xsq) * (y+xsq)
  fn_val = EXP(-xsq*xsq*half) * EXP(-del*half) * fn_val
  IF (x > zero) fn_val = one - fn_val
!------------------------------------------------------------------
!  Evaluate  anorm  for |X| > sqrt(32)
!------------------------------------------------------------------
ELSE
  fn_val = zero
  IF ((x >= xlow).AND.(x < xuppr)) THEN
    xsq = one / (x*x)
    xnum = p(6) * xsq
    xden = xsq
    DO  i = 1, 4
      xnum = (xnum+p(i)) * xsq
      xden = (xden+q(i)) * xsq
    END DO
    fn_val = xsq * (xnum+p(5)) / (xden+q(5))
    fn_val = (sqrpi-fn_val) / y
    xsq = AINT(x*sixten) / sixten
    del = (x-xsq) * (x+xsq)
    fn_val = EXP(-xsq*xsq*half) * EXP(-del*half) * fn_val
  END IF
  IF (x > zero) fn_val = one - fn_val
END IF
!------------------------------------------------------------------
!  Fix up for negative argument, erf, etc.
!------------------------------------------------------------------
RETURN
!---------- Last card of ANORM ----------
END FUNCTION ANORM
