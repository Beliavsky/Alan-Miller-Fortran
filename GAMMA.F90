FUNCTION GAMMA(X) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!----------------------------------------------------------------------

! This routine calculates the GAMMA function for a real argument X.
!   Computation is based on an algorithm outlined in reference 1.
!   The program uses rational functions that approximate the GAMMA
!   function to at least 20 significant decimal digits.  Coefficients
!   for the approximation over the interval (1,2) are unpublished.
!   Those for the approximation for X .GE. 12 are from reference 2.
!   The accuracy achieved depends on the arithmetic system, the
!   compiler, the intrinsic functions, and proper selection of the
!   machine-dependent constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

! XBIG   - the largest argument for which GAMMA(X) is representable
!          in the machine, i.e., the solution to the equation
!                  GAMMA(XBIG) = beta**maxexp
! XINF   - the largest machine representable floating-point number;
!          approximately beta**maxexp
! EPS    - the smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! XMININ - the smallest positive floating-point number such that
!          1/XMININ is machine representable

!     Approximate values for some important machines are:

!                            beta       maxexp        XBIG

! CRAY-1         (S.P.)        2         8191        966.961
! Cyber 180/855
!   under NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-Format   (D.P.)        2          127        34.844
! VAX G-Format   (D.P.)        2         1023        171.489

!                            XINF         EPS        XMININ

! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Cyber 180/855
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308

!*******************************************************************

! Error returns

!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.


!  Intrinsic functions required are:

!     INT, DBLE, EXP, LOG, REAL, SIN


! References: "An Overview of Software Development for Special
!              Functions", W. J. Cody, Lecture Notes in Mathematics,
!              506, Numerical Analysis Dundee, 1975, G. A. Watson
!              (ed.), Springer Verlag, Berlin, 1976.

!              Computer Approximations, Hart, Et. Al., Wiley and
!              sons, New York, 1968.

!  Latest modification: March 12, 1992

!  Authors: W. J. Cody and L. Stoltz
!           Applied Mathematics Division
!           Argonne National Laboratory
!           Argonne, IL 60439

!----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variables
INTEGER    :: i, n
LOGICAL    :: parity
REAL (dp)  :: fact, sum, xden, xnum, y, y1, ysq, z
!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: one = 1.0_dp, half = 0.5_dp, twelve = 12.0_dp,  &
                         two = 2.0_dp, zero = 0.0_dp,  &
                         sqrtpi = 0.9189385332046727417803297_dp,  &
                         pi = 3.1415926535897932384626434_dp
!----------------------------------------------------------------------
!  Machine dependent parameters
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: xbig = 171.624_dp, xminin = 2.23D-308,   &
                         eps = 2.22D-16, xinf = 1.79D308
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     approximation over (1,2).
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P(8) =  &
           (/ -1.71618513886549492533811D+0,  2.47656508055759199108314D+1,  &
              -3.79804256470945635097577D+2,  6.29331155312818442661052D+2,  &
               8.66966202790413211295064D+2, -3.14512729688483675254357D+4,  &
              -3.61444134186911729807069D+4,  6.64561438202405440627855D+4 /)
REAL (dp), PARAMETER  :: Q(8) =  &
           (/ -3.08402300119738975254353D+1,  3.15350626979604161529144D+2,  &
              -1.01515636749021914166146D+3, -3.10777167157231109440444D+3,  &
               2.25381184209801510330112D+4,  4.75584627752788110767815D+3,  &
              -1.34659959864969306392456D+5, -1.15132259675553483497211D+5 /)
!----------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: c(7) =  &
           (/ -1.910444077728D-03, 8.4171387781295D-04,  &
              -5.952379913043012D-04, 7.93650793500350248D-04,  &
              -2.777777777777681622553D-03, 8.333333333333333331554247D-02,  &
               5.7083835261D-03 /)
!----------------------------------------------------------------------

parity = .false.
fact = one
n = 0
y = x
IF (y <= zero) THEN
!----------------------------------------------------------------------
!  Argument is negative
!----------------------------------------------------------------------
  y = -x
  y1 = AINT(y)
  fn_val = y - y1
  IF (fn_val /= zero) THEN
    IF (y1 /= AINT(y1*half)*two) parity = .true.
    fact = -pi / SIN(pi*fn_val)
    y = y + one
  ELSE
    fn_val = xinf
    GO TO 900
  END IF
END IF
!----------------------------------------------------------------------
!  Argument is positive
!----------------------------------------------------------------------
IF (y < eps) THEN
!----------------------------------------------------------------------
!  Argument < EPS
!----------------------------------------------------------------------
  IF (y >= xminin) THEN
    fn_val = one / y
  ELSE
    fn_val = xinf
    GO TO 900
  END IF
ELSE IF (y < twelve) THEN
  y1 = y
  IF (y < one) THEN
!----------------------------------------------------------------------
!  0.0 < argument < 1.0
!----------------------------------------------------------------------
    z = y
    y = y + one
  ELSE
!----------------------------------------------------------------------
!  1.0 < argument < 12.0, reduce argument if necessary
!----------------------------------------------------------------------
    n = INT(y) - 1
    y = y - n
    z = y - one
  END IF
!----------------------------------------------------------------------
!  Evaluate approximation for 1.0 < argument < 2.0
!----------------------------------------------------------------------
  xnum = zero
  xden = one
  DO  i = 1, 8
    xnum = (xnum+p(i)) * z
    xden = xden * z + q(i)
  END DO
  fn_val = xnum / xden + one
  IF (y1 < y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  0.0 < argument < 1.0
!----------------------------------------------------------------------
    fn_val = fn_val / y1
  ELSE IF (y1 > y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  2.0 < argument < 12.0
!----------------------------------------------------------------------
    DO  i = 1, n
      fn_val = fn_val * y
      y = y + one
    END DO
  END IF
ELSE
!----------------------------------------------------------------------
!  Evaluate for argument .GE. 12.0,
!----------------------------------------------------------------------
  IF (y <= xbig) THEN
    ysq = y * y
    sum = c(7)
    DO  i = 1, 6
      sum = sum / ysq + c(i)
    END DO
    sum = sum / y - y + sqrtpi
    sum = sum + (y-half) * LOG(y)
    fn_val = EXP(sum)
  ELSE
    fn_val = xinf
    GO TO 900
  END IF
END IF
!----------------------------------------------------------------------
!  Final adjustments and return
!----------------------------------------------------------------------
IF (parity) fn_val = -fn_val
IF (fact /= one) fn_val = fact / fn_val
900 RETURN
! ---------- Last line of GAMMA ----------
END FUNCTION GAMMA
