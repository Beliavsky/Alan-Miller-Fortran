FUNCTION ALGAMA(X) RESULT(res)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!----------------------------------------------------------------------

! This routine calculates the LOG(GAMMA) function for a positive real
!   argument X.  Computation is based on an algorithm outlined in
!   references 1 and 2.  The program uses rational functions that
!   theoretically approximate LOG(GAMMA) to at least 18 significant
!   decimal digits.  The approximation for X > 12 is from reference
!   3, while approximations for X < 12.0 are similar to those in
!   reference 1, but are unpublished.  The accuracy achieved depends
!   on the arithmetic system, the compiler, the intrinsic functions,
!   and proper selection of the machine-dependent constants.

!*********************************************************************

! Explanation of machine-dependent constants.  Let

! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - largest argument for which LN(GAMMA(X)) is representable
!          in the machine, i.e., the solution to the equation
!                  LN(GAMMA(XBIG)) = beta**maxexp

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

! XINF   - largest machine representable floating-point number;
!          approximately beta**maxexp.
! EPS    - The smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! FRTBIG - Rough estimate of the fourth root of XBIG


!     Approximate values for some important machines are:

!                           beta      maxexp         XBIG

! CRAY-1        (S.P.)        2        8191       9.62E+2461
! Cyber 180/855
!   under NOS   (S.P.)        2        1070       1.72E+319
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)        2         128       4.08E+36
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)        2        1024       2.55D+305
! IBM 3033      (D.P.)       16          63       4.29D+73
! VAX D-Format  (D.P.)        2         127       2.05D+36
! VAX G-Format  (D.P.)        2        1023       1.28D+305


!                           XINF        EPS        FRTBIG

! CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
! Cyber 180/855
!   under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
! IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
! VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
! VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76

!**************************************************************

! Error returns

!  The program returns the value XINF for X .LE. 0.0 or when overflow would
!  occur.  The computation is believed to be free of underflow and overflow.


! Intrinsic functions required are:

!      LOG


! References:

!  1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for the Natural
!     Logarithm of the Gamma Function,' Math. Comp. 21, 1967, pp. 198-203.

!  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May, 1969.

!  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New York, 1968.


!  Authors: W. J. Cody and L. Stoltz
!           Argonne National Laboratory

!  Latest modification: March 9, 1992

!----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: res

! Local variables

INTEGER    :: i
REAL (dp)  :: corr, xden, xm1, xm2, xm4, xnum, y, ysq

!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: one = 1.0_dp, half = 0.5_dp, twelve = 12.0_dp,  &
                         zero = 0.0_dp, four = 4.0_dp, thrhal = 1.5_dp,  &
                         two = 2.0_dp, pnt68 = 0.6796875_dp,  &
                         sqrtpi = 0.9189385332046727417803297_dp
!----------------------------------------------------------------------
!  Machine dependent parameters
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: xbig = HUGE(one), xinf = HUGE(one),   &
                         eps = EPSILON(one), frtbig = 2.25e+76_dp
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     approximation over (0.5,1.5).
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: D1 = -0.5772156649015328605195174_dp
REAL (dp), PARAMETER  :: P1(8) =  &
         (/ 4.945235359296727046734888D0, 2.018112620856775083915565D2,  &
            2.290838373831346393026739D3, 1.131967205903380828685045D4,  &
            2.855724635671635335736389D4, 3.848496228443793359990269D4,  &
            2.637748787624195437963534D4, 7.225813979700288197698961D3 /)
REAL (dp), PARAMETER  :: Q1(8) =  &
         (/ 6.748212550303777196073036D1, 1.113332393857199323513008D3,  &
            7.738757056935398733233834D3, 2.763987074403340708898585D4,  &
            5.499310206226157329794414D4, 6.161122180066002127833352D4,  &
            3.635127591501940507276287D4, 8.785536302431013170870835D3 /)
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     Approximation over (1.5,4.0).
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: D2 = 0.4227843350984671393993777_dp
REAL (dp), PARAMETER  :: P2(8) =  &
         (/ 4.974607845568932035012064D0, 5.424138599891070494101986D2,  &
            1.550693864978364947665077D4, 1.847932904445632425417223D5,  &
            1.088204769468828767498470D6, 3.338152967987029735917223D6,  &
            5.106661678927352456275255D6, 3.074109054850539556250927D6 /)
REAL (dp), PARAMETER  :: Q2(8) =  &
         (/ 1.830328399370592604055942D2, 7.765049321445005871323047D3,  &
            1.331903827966074194402448D5, 1.136705821321969608938755D6,  &
            5.267964117437946917577538D6, 1.346701454311101692290052D7,  &
            1.782736530353274213975932D7, 9.533095591844353613395747D6 /)
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     Approximation over (4.0,12.0).
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: D4 = 1.791759469228055000094023_dp
REAL (dp), PARAMETER  :: P4(8) =  &
         (/ 1.474502166059939948905062D4, 2.426813369486704502836312D6,  &
            1.214755574045093227939592D8, 2.663432449630976949898078D9,  &
            2.940378956634553899906876D10, 1.702665737765398868392998D11,  &
            4.926125793377430887588120D11, 5.606251856223951465078242D11 /)
REAL (dp), PARAMETER  :: Q4(8) =  &
         (/ 2.690530175870899333379843D3, 6.393885654300092398984238D5,  &
            4.135599930241388052042842D7, 1.120872109616147941376570D9,  &
            1.488613728678813811542398D10, 1.016803586272438228077304D11, &
            3.417476345507377132798597D11, 4.463158187419713286462081D11 /)
!----------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: C(7) =  &
         (/ -1.910444077728D-03, 8.4171387781295D-04, -5.952379913043012D-04, &
             7.93650793500350248D-04, -2.777777777777681622553D-03,   &
             8.333333333333333331554247D-02, 5.7083835261D-03 /)
!----------------------------------------------------------------------
y = x
IF (y > zero .AND. y <= xbig) THEN
  IF (y <= eps) THEN
    res = -LOG(y)
  ELSE IF (y <= thrhal) THEN
!----------------------------------------------------------------------
!  EPS .LT. X .LE. 1.5
!----------------------------------------------------------------------
    IF (y < pnt68) THEN
      corr = -LOG(y)
      xm1 = y
    ELSE
      corr = zero
      xm1 = (y-half) - half
    END IF
    IF ((y <= half) .OR. (y >= pnt68)) THEN
      xden = one
      xnum = zero
      DO  i = 1, 8
        xnum = xnum * xm1 + p1(i)
        xden = xden * xm1 + q1(i)
      END DO
      res = corr + (xm1*(d1+xm1*(xnum/xden)))
    ELSE
      xm2 = (y-half) - half
      xden = one
      xnum = zero
      DO  i = 1, 8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
      END DO
      res = corr + xm2 * (d2+xm2*(xnum/xden))
    END IF
  ELSE IF (y <= four) THEN
!----------------------------------------------------------------------
!  1.5 .LT. X .LE. 4.0
!----------------------------------------------------------------------
    xm2 = y - two
    xden = one
    xnum = zero
    DO  i = 1, 8
      xnum = xnum * xm2 + p2(i)
      xden = xden * xm2 + q2(i)
    END DO
    res = xm2 * (d2+xm2*(xnum/xden))
  ELSE IF (y <= twelve) THEN
!----------------------------------------------------------------------
!  4.0 .LT. X .LE. 12.0
!----------------------------------------------------------------------
    xm4 = y - four
    xden = -one
    xnum = zero
    DO  i = 1, 8
      xnum = xnum * xm4 + p4(i)
      xden = xden * xm4 + q4(i)
    END DO
    res = d4 + xm4 * (xnum/xden)
  ELSE
!----------------------------------------------------------------------
!  Evaluate for argument .GE. 12.0,
!----------------------------------------------------------------------
    res = zero
    IF (y <= frtbig) THEN
      res = c(7)
      ysq = y * y
      DO  i = 1, 6
        res = res / ysq + c(i)
      END DO
    END IF
    res = res / y
    corr = LOG(y)
    res = res + sqrtpi - half * corr
    res = res + y * (corr-one)
  END IF
ELSE
!----------------------------------------------------------------------
!  Return for bad arguments
!----------------------------------------------------------------------
  res = xinf
END IF
!----------------------------------------------------------------------
!  Final adjustments and return
!----------------------------------------------------------------------
RETURN
! ---------- Last line of DLGAMA ----------
END FUNCTION ALGAMA
