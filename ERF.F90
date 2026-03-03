MODULE ErrorFunction

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

PRIVATE
PUBLIC  :: calerf, erf, erfc, erfcx

CONTAINS


SUBROUTINE calerf(arg, result, jint)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!------------------------------------------------------------------

! This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!   for a real argument  x.  It contains three FUNCTION type subprograms:
!   ERF, ERFC, and ERFCX and one SUBROUTINE type subprogram, CALERF.
!   The calling statements for the primary entries are:

!                   Y=ERF(X)     (or   Y=DERF(X)),

!                   Y=ERFC(X)    (or   Y=DERFC(X)),
!   and
!                   Y=ERFCX(X)   (or   Y=DERFCX(X)).

!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this routine.
!   The function subprograms invoke  CALERF  with the statement

!          CALL CALERF(ARG, RESULT, JINT)

!   where the parameter usage is as follows

!      Function                     Parameters for CALERF
!       call              ARG              Result          JINT

!     ERF(ARG)      ANY REAL ARGUMENT      ERF(ARG)          0
!     ERFC(ARG)     ABS(ARG) < XBIG        ERFC(ARG)         1
!     ERFCX(ARG)    XNEG < ARG < XMAX      ERFCX(ARG)        2

!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   XMIN   = the smallest positive floating-point number.

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].

!   Approximate values for some important machines are:

!                          XMIN       XINF        XNEG     XSMALL

!  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16


!                          XBIG       XHUGE       XMAX

!  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307

!*******************************************************************

! Error returns

!  The program returns  ERFC = 0      for  ARG >= XBIG;

!                       ERFCX = XINF  for  ARG < XNEG;
!      and
!                       ERFCX = 0     for  ARG >= XMAX.


! Intrinsic functions required are:

!     ABS, AINT, EXP


!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439, USA

!  Latest modification: March 12, 1992

!------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: arg
REAL (dp), INTENT(OUT)  :: result
INTEGER, INTENT(IN)     :: jint

INTEGER    :: i
REAL (dp)  :: del, x, xden, xnum, y, ysq
!------------------------------------------------------------------
!  Mathematical constants
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: four = 4.0_dp, one = 1.0_dp, half = 0.5_dp,   &
                         two = 2.0_dp, zero = 0.0_dp, sixten = 16.0_dp,  &
                         sqrpi = 5.6418958354775628695D-1, thresh = 0.46875_dp
!------------------------------------------------------------------
!  Machine-dependent constants
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: xinf = 1.79D308, xneg = -26.628_dp, xsmall = 1.11D-16, &
                         xbig = 26.543_dp, xhuge = 6.71D7, xmax = 2.53D307
!------------------------------------------------------------------
!  Coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: a(5) =  &
           (/ 3.16112374387056560D00, 1.13864154151050156D02,  &
              3.77485237685302021D02, 3.20937758913846947D03,  &
              1.85777706184603153D-1 /)
REAL (dp), PARAMETER  :: b(4) =  &
           (/ 2.36012909523441209D01, 2.44024637934444173D02,  &
              1.28261652607737228D03, 2.84423683343917062D03 /)
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: c(9) =  &
           (/ 5.64188496988670089D-1, 8.88314979438837594D0,   &
              6.61191906371416295D01, 2.98635138197400131D02,  &
              8.81952221241769090D02, 1.71204761263407058D03,  &
              2.05107837782607147D03, 1.23033935479799725D03,  &
              2.15311535474403846D-8 /)
REAL (dp), PARAMETER  :: d(8) =  &
           (/ 1.57449261107098347D01, 1.17693950891312499D02,  &
              5.37181101862009858D02, 1.62138957456669019D03,  &
              3.29079923573345963D03, 4.36261909014324716D03,  &
              3.43936767414372164D03, 1.23033935480374942D03 /)
!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
REAL (dp), PARAMETER  :: p(6) =  &
           (/ 3.05326634961232344D-1, 3.60344899949804439D-1,  &
              1.25781726111229246D-1, 1.60837851487422766D-2,  &
              6.58749161529837803D-4, 1.63153871373020978D-2 /)
REAL (dp), PARAMETER  :: q(5) =  &
           (/ 2.56852019228982242D00, 1.87295284992346047D00,  &
              5.27905102951428412D-1, 6.05183413124413191D-2,  &
              2.33520497626869185D-3 /)
!------------------------------------------------------------------

x = arg
y = ABS(x)
IF (y <= thresh) THEN
!------------------------------------------------------------------
!  Evaluate  erf  for  |X| <= 0.46875
!------------------------------------------------------------------
  ysq = zero
  IF (y > xsmall) ysq = y * y
  xnum = a(5) * ysq
  xden = ysq
  DO  i = 1, 3
    xnum = (xnum+a(i)) * ysq
    xden = (xden+b(i)) * ysq
  END DO
  result = x * (xnum+a(4)) / (xden+b(4))
  IF (jint /= 0) result = one - result
  IF (jint == 2) result = EXP(ysq) * result
  GO TO 50
!------------------------------------------------------------------
!  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!------------------------------------------------------------------
ELSE IF (y <= four) THEN
  xnum = c(9) * y
  xden = y
  DO  i = 1, 7
    xnum = (xnum+c(i)) * y
    xden = (xden+d(i)) * y
  END DO
  result = (xnum+c(8)) / (xden+d(8))
  IF (jint /= 2) THEN
    ysq = AINT(y*sixten) / sixten
    del = (y-ysq) * (y+ysq)
    result = EXP(-ysq*ysq) * EXP(-del) * result
  END IF
!------------------------------------------------------------------
!  Evaluate  erfc  for |X| > 4.0
!------------------------------------------------------------------
ELSE
  result = zero
  IF (y >= xbig) THEN
    IF ((jint /= 2).OR.(y >= xmax)) GO TO 40
    IF (y >= xhuge) THEN
      result = sqrpi / y
      GO TO 40
    END IF
  END IF
  ysq = one / (y*y)
  xnum = p(6) * ysq
  xden = ysq
  DO  i = 1, 4
    xnum = (xnum+p(i)) * ysq
    xden = (xden+q(i)) * ysq
  END DO
  result = ysq * (xnum+p(5)) / (xden+q(5))
  result = (sqrpi-result) / y
  IF (jint /= 2) THEN
    ysq = AINT(y*sixten) / sixten
    del = (y-ysq) * (y+ysq)
    result = EXP(-ysq*ysq) * EXP(-del) * result
  END IF
END IF
!------------------------------------------------------------------
!  Fix up for negative argument, erf, etc.
!------------------------------------------------------------------
40 IF (jint == 0) THEN
  result = (half-result) + half
  IF (x < zero) result = -result
ELSE IF (jint == 1) THEN
  IF (x < zero) result = two - result
ELSE
  IF (x < zero) THEN
    IF (x < xneg) THEN
      result = xinf
    ELSE
      ysq = AINT(x*sixten) / sixten
      del = (x-ysq) * (x+ysq)
      y = EXP(ysq*ysq) * EXP(del)
      result = (y+y) - result
    END IF
  END IF
END IF
50 RETURN
!---------- Last card of CALERF ----------
END SUBROUTINE calerf



FUNCTION ERF(X) RESULT(fn_val)
!--------------------------------------------------------------------

! This subprogram computes approximate values for erf(x).
!   (see comments heading CALERF).

!   Author/date: W. J. Cody, January 8, 1985

!--------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variable
INTEGER  :: jint
!------------------------------------------------------------------
jint = 0
CALL calerf(x, fn_val, jint)
RETURN
!---------- Last card of DERF ----------
END FUNCTION ERF



FUNCTION ERFC(X) RESULT(fn_val)
!--------------------------------------------------------------------

! This subprogram computes approximate values for erfc(x).
!   (see comments heading CALERF).

!   Author/date: W. J. Cody, January 8, 1985

!--------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variable
INTEGER  :: jint
!------------------------------------------------------------------
jint = 1
CALL calerf(x, fn_val, jint)
RETURN
!---------- Last card of DERFC ----------
END FUNCTION ERFC



FUNCTION ERFCX(X) RESULT(fn_val)
!------------------------------------------------------------------

! This subprogram computes approximate values for exp(x*x) * erfc(x).
!   (see comments heading CALERF).

!   Author/date: W. J. Cody, March 30, 1987

!------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

! Local variable
INTEGER  :: jint
!------------------------------------------------------------------
jint = 2
CALL calerf(x, fn_val, jint)
RETURN
!---------- Last card of DERFCX ----------
END FUNCTION ERFCX

END MODULE ErrorFunction
