FUNCTION DAW(XX) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!----------------------------------------------------------------------

! This function program evaluates Dawson's integral,

!                       2  / x   2
!                     -x   |    t
!             F(x) = e     |   e    dt
!                          |
!                          / 0

!   for a real argument x.

!   The calling sequence for this function is

!                   Y=DAW(X)

!   The main computation uses rational Chebyshev approximations published
!   in Math. Comp. 24, 171-178 (1970) by Cody, Paciorek and Thacher.
!   This transportable program is patterned after the machine-dependent
!   FUNPACK program DDAW(X), but cannot match that version for efficiency or
!   accuracy.  This version uses rational approximations that are
!   theoretically accurate to about 19 significant decimal digits.
!   The accuracy achieved depends on the arithmetic system, the compiler,
!   the intrinsic functions, and proper selection of the machine-dependent
!   constants.

!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   XINF   = largest positive machine number
!   XMIN   = the smallest positive machine number.
!   EPS    = smallest positive number such that 1+eps > 1.
!            Approximately  beta**(-p), where beta is the machine radix
!            and p is the number of significant base-beta digits in a
!            floating-point number.

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   XMAX   = absolute argument beyond which DAW(X) underflows.
!            XMAX = min(0.5/xmin, xinf).
!   XSMALL = absolute argument below DAW(X)  may be represented
!            by X.  We recommend XSMALL = sqrt(eps).
!   XLARGE = argument beyond which DAW(X) may be represented by
!            1/(2x).  We recommend XLARGE = 1/sqrt(eps).

!     Approximate values for some important machines are

!                        beta  p     eps     xmin       xinf

!  CDC 7600      (S.P.)    2  48  7.11E-15  3.14E-294  1.26E+322
!  CRAY-1        (S.P.)    2  48  7.11E-15  4.58E-2467 5.45E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2  24  1.19E-07  1.18E-38   3.40E+38
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2  53  1.11D-16  2.23E-308  1.79D+308
!  IBM 3033      (D.P.)   16  14  1.11D-16  5.40D-79   7.23D+75
!  VAX 11/780    (S.P.)    2  24  5.96E-08  2.94E-39   1.70E+38
!                (D.P.)    2  56  1.39D-17  2.94D-39   1.70D+38
!   (G Format)   (D.P.)    2  53  1.11D-16  5.57D-309  8.98D+307

!                         XSMALL     XLARGE     XMAX

!  CDC 7600      (S.P.)  5.96E-08   1.68E+07  1.59E+293
!  CRAY-1        (S.P.)  5.96E-08   1.68E+07  5.65E+2465
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)  2.44E-04   4.10E+03  4.25E+37
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)  1.05E-08   9.49E+07  2.24E+307
!  IBM 3033      (D.P.)  3.73D-09   2.68E+08  7.23E+75
!  VAX 11/780    (S.P.)  2.44E-04   4.10E+03  1.70E+38
!                (D.P.)  3.73E-09   2.68E+08  1.70E+38
!   (G Format)   (D.P.)  1.05E-08   9.49E+07  8.98E+307

!*******************************************************************

! Error Returns

!  The program returns 0.0 for |X| > XMAX.

! Intrinsic functions required are:

!     ABS


!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!  Latest modification: March 9, 1992

!----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: xx
REAL (dp)              :: fn_val

! Local variables

INTEGER    :: i
REAL (dp)  :: frac, sump, sumq, w2, x, y
!----------------------------------------------------------------------
!  Mathematical constants.
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         six25 = 6.25_dp, one225 = 12.25_dp, two5 = 25.0_dp
!----------------------------------------------------------------------
!  Machine-dependent constants
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: XSMALL = 1.05D-08, XLARGE = 9.49D+07,   &
                         XMAX = 2.24D+307
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation for  |x| < 2.5
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P1(10) = (/  &
        -2.69020398788704782410D-12, 4.18572065374337710778D-10,  &
        -1.34848304455939419963D-08, 9.28264872583444852976D-07,  &
        -1.23877783329049120592D-05, 4.07205792429155826266D-04,  &
        -2.84388121441008500446D-03, 4.70139022887204722217D-02,  &
        -1.38868086253931995101D-01, 1.00000000000000000004D+00 /)
REAL (dp), PARAMETER  :: Q1(10) = (/  &
         1.71257170854690554214D-10, 1.19266846372297253797D-08,  &
         4.32287827678631772231D-07, 1.03867633767414421898D-05,  &
         1.78910965284246249340D-04, 2.26061077235076703171D-03,  &
         2.07422774641447644725D-02, 1.32212955897210128811D-01,  &
         5.27798580412734677256D-01, 1.00000000000000000000D+00 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [2.5, 3.5)
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P2(10) = (/  &
        -1.70953804700855494930D+00,-3.79258977271042880786D+01,  &
         2.61935631268825992835D+01, 1.25808703738951251885D+01,  &
        -2.27571829525075891337D+01, 4.56604250725163310122D+00,  &
        -7.33080089896402870750D+00, 4.65842087940015295573D+01,  &
        -1.73717177843672791149D+01, 5.00260183622027967838D-01 /)
REAL (dp), PARAMETER  :: Q2(9) = (/  &
         1.82180093313514478378D+00, 1.10067081034515532891D+03,  &
        -7.08465686676573000364D+00, 4.53642111102577727153D+02,  &
         4.06209742218935689922D+01, 3.02890110610122663923D+02,  &
         1.70641269745236227356D+02, 9.51190923960381458747D+02,  &
         2.06522691539642105009D-01 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  x in [3.5, 5.0]
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P3(10) = (/  &
        -4.55169503255094815112D+00,-1.86647123338493852582D+01,  &
        -7.36315669126830526754D+00,-6.68407240337696756838D+01,  &
         4.84507265081491452130D+01, 2.69790586735467649969D+01,  &
        -3.35044149820592449072D+01, 7.50964459838919612289D+00,  &
        -1.48432341823343965307D+00, 4.99999810924858824981D-01 /)
REAL (dp), PARAMETER  :: Q3(9) = (/  &
         4.47820908025971749852D+01, 9.98607198039452081913D+01,  &
         1.40238373126149385228D+01, 3.48817758822286353588D+03,  &
        -9.18871385293215873406D+00, 1.24018500009917163023D+03,  &
        -6.88024952504512254535D+01,-2.31251575385145143070D+00,  &
         2.50041492369922381761D-01 /)
!----------------------------------------------------------------------
!  Coefficients for R(9,9) approximation in J-fraction form
!     for  |x| > 5.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P4(10) = (/  &
        -8.11753647558432685797D+00,-3.84043882477454453430D+01,  &
        -2.23787669028751886675D+01,-2.88301992467056105854D+01,  &
        -5.99085540418222002197D+00,-1.13867365736066102577D+01,  &
        -6.52828727526980741590D+00,-4.50002293000355585708D+00,  &
        -2.50000000088955834952D+00, 5.00000000000000488400D-01 /)
REAL (dp), PARAMETER  :: Q4(9) = (/  &
         2.69382300417238816428D+02, 5.04198958742465752861D+01,  &
         6.11539671480115846173D+01, 2.08210246935564547889D+02,  &
         1.97325365692316183531D+01,-1.22097010558934838708D+01,  &
        -6.99732735041547247161D+00,-2.49999970104184464568D+00,  &
         7.49999999999027092188D-01 /)
!----------------------------------------------------------------------
x = xx
IF (ABS(x) > xlarge) THEN
  IF (ABS(x) <= xmax) THEN
    fn_val = half / x
  ELSE
    fn_val = zero
  END IF
ELSE IF (ABS(x) < xsmall) THEN
  fn_val = x
ELSE
  y = x * x
  IF (y < six25) THEN
!----------------------------------------------------------------------
!  ABS(X) .LT. 2.5
!----------------------------------------------------------------------
    sump = p1(1)
    sumq = q1(1)
    DO  i = 2, 10
      sump = sump * y + p1(i)
      sumq = sumq * y + q1(i)
    END DO
    fn_val = x * sump / sumq
  ELSE IF (y < one225) THEN
!----------------------------------------------------------------------
!  2.5 .LE. ABS(X) .LT. 3.5
!----------------------------------------------------------------------
    frac = zero
    DO  i = 1, 9
      frac = q2(i) / (p2(i)+y+frac)
    END DO
    fn_val = (p2(10)+frac) / x
  ELSE IF (y < two5) THEN
!----------------------------------------------------------------------
!  3.5 .LE. ABS(X) .LT. 5.0
!---------------------------------------------------------------------
    frac = zero
    DO  i = 1, 9
      frac = q3(i) / (p3(i)+y+frac)
    END DO
    fn_val = (p3(10)+frac) / x
  ELSE
!----------------------------------------------------------------------
!  5.0 .LE. ABS(X) .LE. XLARGE
!------------------------------------------------------------------
    w2 = one / x / x
    frac = zero
    DO  i = 1, 9
      frac = q4(i) / (p4(i)+y+frac)
    END DO
    frac = p4(10) + frac
    fn_val = (half + half*w2*frac) / x
  END IF
END IF
RETURN
!---------- Last line of DAW ----------
END FUNCTION daw
