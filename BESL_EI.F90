MODULE Bessel_EI

! The Ei Bessel functions

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

PRIVATE
PUBLIC  :: calcei, ei, eone, expei, dsubn

CONTAINS


SUBROUTINE calcei(arg, result, INT)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!----------------------------------------------------------------------

! This Fortran 90 packet computes the exponential integrals Ei(x),
!  E1(x), and  exp(-x)*Ei(x)  for real arguments  x  where

!           integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
!  Ei(x) =
!          -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,

!  and where the first integral is a principal value integral.
!  The packet contains three function type subprograms: EI, EONE, and EXPEI;
!  and one subroutine type subprogram: CALCEI.  The calling statements for
!  the primary entries are

!                 Y = EI(X),            where  X .NE. 0,

!                 Y = EONE(X),          where  X .GT. 0,
!  and
!                 Y = EXPEI(X),         where  X .NE. 0,

!  and where the entry points correspond to the functions Ei(x),
!  E1(x), and exp(-x)*Ei(x), respectively.  The routine CALCEI
!  is intended for internal packet use only, all computations within
!  the packet being concentrated in this routine.  The function
!  subprograms invoke CALCEI with the Fortran statement
!         CALL CALCEI(ARG, RESULT, INT)
!  where the parameter usage is as follows

!     Function                  Parameters for CALCEI
!       Call                 ARG             RESULT         INT

!      EI(X)              X .NE. 0          Ei(X)            1
!      EONE(X)            X .GT. 0         -Ei(-X)           2
!      EXPEI(X)           X .NE. 0          exp(-X)*Ei(X)    3

!  The main computation involves evaluation of rational Chebyshev
!  approximations published in Math. Comp. 22, 641-649 (1968), and
!  Math. Comp. 23, 289-303 (1969) by Cody and Thacher.  This
!  transportable program is patterned after the machine-dependent
!  FUNPACK packet  NATSEI,  but cannot match that version for
!  efficiency or accuracy.  This version uses rational functions
!  that theoretically approximate the exponential integrals to
!  at least 18 significant decimal digits.  The accuracy achieved
!  depends on the arithmetic system, the compiler, the intrinsic
!  functions, and proper selection of the machine-dependent constants.


!*******************************************************************

! Explanation of machine-dependent constants.  Let

!   beta = radix for the floating-point system.
!   minexp = smallest representable power of beta.
!   maxexp = smallest power of beta that overflows.

! Then the following machine-dependent constants must be declared
!   in DATA statements.  IEEE values are provided as a default.

!   XBIG = largest argument acceptable to EONE; solution to
!          equation:
!                     exp(-x)/x * (1 + 1/x) = beta ** minexp.
!   XINF = largest positive machine number; approximately
!                     beta ** maxexp
!   XMAX = largest argument acceptable to EI; solution to
!          equation:  exp(x)/x * (1 + 1/x) = beta ** maxexp.

!     Approximate values for some important machines are:

!                           beta      minexp      maxexp

!  CRAY-1        (S.P.)       2       -8193        8191
!  Cyber 180/185
!    under NOS   (S.P.)       2        -975        1070
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)       2        -126         128
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)       2       -1022        1024
!  IBM 3033      (D.P.)      16         -65          63
!  VAX D-Format  (D.P.)       2        -128         127
!  VAX G-Format  (D.P.)       2       -1024        1023

!                           XBIG       XINF       XMAX

!  CRAY-1        (S.P.)    5670.31  5.45E+2465   5686.21
!  Cyber 180/185
!    under NOS   (S.P.)     669.31  1.26E+322     748.28
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)      82.93  3.40E+38       93.24
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)     701.84  1.79D+308     716.35
!  IBM 3033      (D.P.)     175.05  7.23D+75      179.85
!  VAX D-Format  (D.P.)      84.30  1.70D+38       92.54
!  VAX G-Format  (D.P.)     703.22  8.98D+307     715.66

!*******************************************************************

! Error returns

!  The following table shows the types of error that may be encountered
!  in this routine and the function value supplied in each case.

!       Error       Argument         Function values for
!                    Range         EI      EXPEI     EONE

!     UNDERFLOW  (-)X .GT. XBIG     0        -         0
!     OVERFLOW      X .GE. XMAX    XINF      -         -
!     ILLEGAL X       X = 0       -XINF    -XINF     XINF
!     ILLEGAL X      X .LT. 0       -        -     USE ABS(X)

! Intrinsic functions required are:

!     ABS, SQRT, EXP


!  Author: W. J. Cody
!          Mathematics abd Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439

!  Latest modification: March 9, 1992

!----------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: arg
REAL (dp), INTENT(OUT)  :: result
INTEGER, INTENT(IN)     :: INT

INTEGER    :: i
REAL (dp)  :: ei, frac, px(10), qx(10), sump, sumq, t, w, x, xmx0, y, ysq
!----------------------------------------------------------------------
!  Mathematical constants
!   EXP40 = exp(40)
!   X0 = zero of Ei
!   X01/X11 + X02 = zero of Ei to extra precision
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: zero = 0.0_dp, p037 = 0.037_dp, half = 0.5_dp,  &
                         one = 1.0_dp, two = 2.0_dp, three = 3.0_dp,   &
                         four = 4.0_dp, six = 6.0_dp, twelve = 12.0_dp,  &
                         two4 = 24.0_dp, fourty = 40.0_dp,   &
                         exp40 = 2.3538526683701998541D17,  &
                         x01 = 381.5_dp, x11 = 1024.0_dp,   &
                         x02 = -5.1182968633365538008D-5,   &
                         x0 = 3.7250741078136663466D-1
!----------------------------------------------------------------------
! Machine-dependent constants
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: XINF = 1.79D+308, XMAX = 716.351_dp, XBIG = 701.84_dp
!----------------------------------------------------------------------
! Coefficients  for -1.0 <= X < 0.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: A(7) =   &
           (/ 1.1669552669734461083368D2, 2.1500672908092918123209D3,  &
              1.5924175980637303639884D4, 8.9904972007457256553251D4,  &
              1.5026059476436982420737D5, -1.4815102102575750838086D5,  &
              5.0196785185439843791020_dp /)
REAL (dp), PARAMETER  :: B(6) =   &
           (/ 4.0205465640027706061433D1, 7.5043163907103936624165D2,  &
              8.1258035174768735759855D3, 5.2440529172056355429883D4,  &
              1.8434070063353677359298D5, 2.5666493484897117319268D5 /)
!----------------------------------------------------------------------
! Coefficients for -4.0 <= X < -1.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: C(9) =   &
           (/ 3.828573121022477169108D-1, 1.107326627786831743809D+1,  &
              7.246689782858597021199D+1, 1.700632978311516129328D+2,  &
              1.698106763764238382705D+2, 7.633628843705946890896D+1,  &
              1.487967702840464066613D+1, 9.999989642347613068437D-1,  &
              1.737331760720576030932D-8 /)
REAL (dp), PARAMETER  :: D(9) =   &
           (/ 8.258160008564488034698D-2, 4.344836335509282083360D+0,  &
              4.662179610356861756812D+1, 1.775728186717289799677D+2,  &
              2.953136335677908517423D+2, 2.342573504717625153053D+2,  &
              9.021658450529372642314D+1, 1.587964570758947927903D+1,  &
              1.000000000000000000000_dp /)
!----------------------------------------------------------------------
! Coefficients for X < -4.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: E(10) =   &
           (/ 1.3276881505637444622987D+2, 3.5846198743996904308695D+4,  &
              1.7283375773777593926828D+5, 2.6181454937205639647381D+5,  &
              1.7503273087497081314708D+5, 5.9346841538837119172356D+4,  &
              1.0816852399095915622498D+4, 1.0611777263550331766871D03,  &
              5.2199632588522572481039D+1, 9.9999999999999999087819D-1 /)
REAL (dp), PARAMETER  :: F(10) =   &
           (/ 3.9147856245556345627078D+4, 2.5989762083608489777411D+5,  &
              5.5903756210022864003380D+5, 5.4616842050691155735758D+5,  &
              2.7858134710520842139357D+5, 7.9231787945279043698718D+4,  &
              1.2842808586627297365998D+4, 1.1635769915320848035459D+3,  &
              5.4199632588522559414924D+1, 1.0_dp /)
!----------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: PLG(4) =   &
           (/ -2.4562334077563243311D+01, 2.3642701335621505212D+02,   &
              -5.4989956895857911039D+02, 3.5687548468071500413D+02 /)
REAL (dp), PARAMETER  :: QLG(4) =   &
           (/ -3.5553900764052419184D+01, 1.9400230218539473193D+02,   &
              -3.3442903192607538956D+02, 1.7843774234035750207D+02 /)
!----------------------------------------------------------------------
! Coefficients for  0.0 < X < 6.0,
!  ratio of Chebyshev polynomials
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P(10) =   &
           (/ -1.2963702602474830028590D01, -1.2831220659262000678155D03,  &
              -1.4287072500197005777376D04, -1.4299841572091610380064D06,  &
              -3.1398660864247265862050D05, -3.5377809694431133484800D08,  &
               3.1984354235237738511048D08, -2.5301823984599019348858D10,  &
               1.2177698136199594677580D10, -2.0829040666802497120940D11 /)
REAL (dp), PARAMETER  :: Q(10) =   &
           (/ 7.6886718750000000000000D01, -5.5648470543369082846819D03,  &
              1.9418469440759880361415D05, -4.2648434812177161405483D06,  &
              6.4698830956576428587653D07, -7.0108568774215954065376D08,  &
              5.4229617984472955011862D09, -2.8986272696554495342658D10,  &
              9.8900934262481749439886D10, -8.9673749185755048616855D10 /)
!----------------------------------------------------------------------
! J-fraction coefficients for 6.0 <= X < 12.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: R(10) =   &
           (/ -2.645677793077147237806_dp, -2.378372882815725244124_dp,  &
              -2.421106956980653511550D01,  1.052976392459015155422D01,  &
               1.945603779539281810439D01, -3.015761863840593359165D01,  &
               1.120011024227297451523D01, -3.988850730390541057912_dp,  &
               9.565134591978630774217_dp,  9.981193787537396413219D-1 /)
REAL (dp), PARAMETER  :: S(9) =   &
           (/ 1.598517957704779356479D-4,  4.644185932583286942650_dp,  &
              3.697412299772985940785D02, -8.791401054875438925029_dp,  &
              7.608194509086645763123D02,  2.852397548119248700147D01,  &
              4.731097187816050252967D02, -2.369210235636181001661D02,  &
              1.249884822712447891440_dp /)
!----------------------------------------------------------------------
! J-fraction coefficients for 12.0 <= X < 24.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P1(10) =   &
           (/ -1.647721172463463140042_dp, -1.860092121726437582253D01,  &
              -1.000641913989284829961D01, -2.105740799548040450394D01,  &
              -9.134835699998742552432D-1, -3.323612579343962284333D01,  &
               2.495487730402059440626D01,  2.652575818452799819855D01,  &
              -1.845086232391278674524_dp,  9.999933106160568739091D-1 /)
REAL (dp), PARAMETER  :: Q1(9) =   &
           (/ 9.792403599217290296840D01,  6.403800405352415551324D01,  &
              5.994932325667407355255D01,  2.538819315630708031713D02,  &
              4.429413178337928401161D01,  1.192832423968601006985D03,  &
              1.991004470817742470726D02, -1.093556195391091143924D01,  &
              1.001533852045342697818_dp /)
!----------------------------------------------------------------------
! J-fraction coefficients for  X .GE. 24.0
!----------------------------------------------------------------------
REAL (dp), PARAMETER  :: P2(10) =   &
           (/  1.75338801265465972390D02, -2.23127670777632409550D02,  &
              -1.81949664929868906455D01, -2.79798528624305389340D01,  &
              -7.63147701620253630855_dp, -1.52856623636929636839D01,  &
              -7.06810977895029358836_dp, -5.00006640413131002475_dp,  &
              -3.00000000320981265753_dp,  1.00000000000000485503_dp /)
REAL (dp), PARAMETER  :: Q2(9) =   &
           (/  3.97845977167414720840D04,  3.97277109100414518365_dp,  &
               1.37790390235747998793D02,  1.17179220502086455287D02,  &
               7.04831847180424675988D01, -1.20187763547154743238D01,  &
              -7.99243595776339741065_dp, -2.99999894040324959612_dp,  &
               1.99999999999048104167_dp /)
!----------------------------------------------------------------------

x = arg
IF (x == zero) THEN
  ei = -xinf
  IF (INT == 2) ei = -ei
ELSE IF (x < zero .OR. INT == 2) THEN
!----------------------------------------------------------------------
! Calculate EI for negative argument or for E1.
!----------------------------------------------------------------------
  y = ABS(x)
  IF (y <= one) THEN
    sump = a(7) * y + a(1)
    sumq = y + b(1)
    DO  i = 2, 6
      sump = sump * y + a(i)
      sumq = sumq * y + b(i)
    END DO
    ei = LOG(y) - sump / sumq
    IF (INT == 3) ei = ei * EXP(y)
  ELSE IF (y <= four) THEN
    w = one / y
    sump = c(1)
    sumq = d(1)
    DO  i = 2, 9
      sump = sump * w + c(i)
      sumq = sumq * w + d(i)
    END DO
    ei = -sump / sumq
    IF (INT /= 3) ei = ei * EXP(-y)
  ELSE
    IF (y > xbig .AND. INT < 3) THEN
      ei = zero
    ELSE
      w = one / y
      sump = e(1)
      sumq = f(1)
      DO  i = 2, 10
        sump = sump * w + e(i)
        sumq = sumq * w + f(i)
      END DO
      ei = -w * (one-w*sump/sumq)
      IF (INT /= 3) ei = ei * EXP(-y)
    END IF
  END IF
  IF (INT == 2) ei = -ei
ELSE IF (x < six) THEN
!----------------------------------------------------------------------
!  To improve conditioning, rational approximations are expressed
!    in terms of Chebyshev polynomials for 0 <= X < 6, and in
!    continued fraction form for larger X.
!----------------------------------------------------------------------
  t = x + x
  t = t / three - two
  px(1) = zero
  qx(1) = zero
  px(2) = p(1)
  qx(2) = q(1)
  DO  i = 2, 9
    px(i+1) = t * px(i) - px(i-1) + p(i)
    qx(i+1) = t * qx(i) - qx(i-1) + q(i)
  END DO
  sump = half * t * px(10) - px(9) + p(10)
  sumq = half * t * qx(10) - qx(9) + q(10)
  frac = sump / sumq
  xmx0 = (x-x01/x11) - x02
  IF (ABS(xmx0) >= p037) THEN
    ei = LOG(x/x0) + xmx0 * frac
    IF (INT == 3) ei = EXP(-x) * ei
  ELSE
!----------------------------------------------------------------------
! Special approximation to  ln(X/X0)  for X close to X0
!----------------------------------------------------------------------
    y = xmx0 / (x+x0)
    ysq = y * y
    sump = plg(1)
    sumq = ysq + qlg(1)
    DO  i = 2, 4
      sump = sump * ysq + plg(i)
      sumq = sumq * ysq + qlg(i)
    END DO
    ei = (sump/(sumq*(x+x0))+frac) * xmx0
    IF (INT == 3) ei = EXP(-x) * ei
  END IF
ELSE IF (x < twelve) THEN
  frac = zero
  DO  i = 1, 9
    frac = s(i) / (r(i)+x+frac)
  END DO
  ei = (r(10)+frac) / x
  IF (INT /= 3) ei = ei * EXP(x)
ELSE IF (x <= two4) THEN
  frac = zero
  DO  i = 1, 9
    frac = q1(i) / (p1(i)+x+frac)
  END DO
  ei = (p1(10)+frac) / x
  IF (INT /= 3) ei = ei * EXP(x)
ELSE
  IF (x >= xmax .AND. INT < 3) THEN
    ei = xinf
  ELSE
    y = one / x
    frac = zero
    DO  i = 1, 9
      frac = q2(i) / (p2(i)+x+frac)
    END DO
    frac = p2(10) + frac
    ei = y + y * y * frac
    IF (INT /= 3) THEN
      IF (x <= xmax-two4) THEN
        ei = ei * EXP(x)
      ELSE
!----------------------------------------------------------------------
! Calculation reformulated to avoid premature overflow
!----------------------------------------------------------------------
        ei = (ei*EXP(x-fourty)) * exp40
      END IF
    END IF
  END IF
END IF
result = ei
RETURN
!---------- Last line of CALCEI ----------
END SUBROUTINE calcei



FUNCTION EI(X) RESULT(fn_val)
!--------------------------------------------------------------------

! This function program computes approximate values for the
!   exponential integral  Ei(x), where  x  is real.

!  Author: W. J. Cody

!  Latest modification: March 9, 1992

!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

INTEGER  :: INT

!--------------------------------------------------------------------
INT = 1
CALL calcei(x, fn_val, INT)
RETURN
!---------- Last line of EI ----------
END FUNCTION EI



FUNCTION EXPEI(X) RESULT(fn_val)
!--------------------------------------------------------------------

! This function program computes approximate values for the
!   function  exp(-x) * Ei(x), where  Ei(x)  is the exponential
!   integral, and  x  is real.

!  Author: W. J. Cody

!  Latest modification: March 9, 1992

!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

INTEGER  :: INT
!--------------------------------------------------------------------
INT = 3
CALL calcei(x, fn_val, INT)
RETURN
!---------- Last line of EXPEI ----------
END FUNCTION EXPEI



FUNCTION EONE(X) RESULT(fn_val)
!--------------------------------------------------------------------

! This function program computes approximate values for the
!   exponential integral E1(x), where  x  is real.

!  Author: W. J. Cody

!  Latest modification: March 9, 1992

!--------------------------------------------------------------------
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: x
REAL (dp)              :: fn_val

INTEGER  :: INT
!--------------------------------------------------------------------
INT = 2
CALL calcei(x, fn_val, INT)
RETURN
!---------- Last line of EONE ----------
END FUNCTION EONE



SUBROUTINE dsubn(x, nmax, xmax, d)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-14  Time: 15:25:00
 
!-------------------------------------------------------------------
! Translation of Gautschi'f CACM Algorithm 282 for
!   derivatives of Ei(x).

!  Intrinsic functions required are:

!      ABS, EXP, INT, LOG, MIN

!-------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: x
INTEGER, INTENT(IN)     :: nmax
REAL (dp), INTENT(IN)   :: xmax
REAL (dp), INTENT(OUT)  :: d(0:nmax)

! Local variables

LOGICAL    :: bool1, bool2
INTEGER    :: j, lim, n0, mini, n, n1
REAL (dp)  :: e, en, p, q, t, x1, z

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp,       &
                         ten = 10.0_dp, c0 = 2.7183_dp, c1 = 46.7452_dp,  &
                         b0 = 5.7941D-5, B1 = -1.76148D-3, B2 = 2.08645D-2, &
                         B3 = -1.29013D-1, B4 = 8.5777D-1, B5 = 1.0125_dp,  &
                         B6 = 7.75D-1
!-------------------------------------------------------------------
x1 = ABS(x)
n0 = INT(x1)
e = EXP(x)
d(0) = e / x
bool1 = (x < zero) .OR. (x1 <= two)
bool2 = n0 < nmax
mini = MIN(n0,nmax)
IF (bool1) THEN
  lim = nmax
ELSE
  lim = mini
END IF
n = 1
en = one
10 d(n) = (e-en*d(n-1)) / x
n = n + 1
en = en + one
IF (x1 < one) THEN
  IF ((ABS(d(n-1)) < ABS(xmax*x/en)) .AND. (n <= lim)) GO TO 10
ELSE
  IF ((ABS(d(n-1)/x) < xmax/en) .AND. (n <= lim)) GO TO 10
END IF
DO  j = n, lim
  d(n) = zero
END DO
IF ((.NOT.bool1) .AND. bool2) THEN
  t = (x1+c1) / (c0*x1)
  IF (t < ten) THEN
    t = ((((b0*t+b1)*t+b2)*t+b3)*t+b4) * t + b5
  ELSE
    z = LOG(t) - b6
    p = (b6-LOG(z)) / (one+z)
    p = one / (one+p)
    t = t * p / z
  END IF
  n1 = c0 * x1 * t - one
  IF (n1 < nmax) n1 = nmax
  q = one / x
  en = one
  DO  n = 1, n1 + 1
    q = -en * q / x
    en = en + one
  END DO
  DO  n = n1, n0 + 1, -1
    en = en - one
    q = (e-x*q) / en
    IF (n <= nmax) d(n) = q
  END DO
END IF
RETURN
!---------- Last line of DSUBN ----------
END SUBROUTINE dsubn

END MODULE Bessel_EI
