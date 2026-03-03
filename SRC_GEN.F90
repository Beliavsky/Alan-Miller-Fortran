MODULE toms757

! ALGORITHM 757, COLLECTED ALGORITHMS FROM ACM.
! THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
! VOL. 22, NO. 3, September, 1996, P.  288--301.

! Code converted using TO_F90 by Alan Miller
! Date: 2001-10-10  Time: 12:42:51

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


FUNCTION abram0(xvalue) RESULT(fn_val)
 
! DESCRIPTION:
!    This function calculates the Abramowitz function of order 0, defined as

!     ABRAM0(x) = integral{ 0 to infinity } exp( -t*t - x/t ) dt

!     The code uses Chebyshev expansions with the coefficients
!     given to an accuracy of 20 decimal places.


! ERROR RETURNS:
!    If XVALUE < 0.0, the function prints a message and returns the value 0.0.


! MACHINE-DEPENDENT CONSTANTS:

!    NTERMF - INTEGER - No. of terms needed for the array AB0F.
!             Recommended value such that
!                   ABS( AB0F(NTERMF) ) < EPS/100

!    NTERMG - INTEGER - No. of terms needed for array AB0G.
!             Recommended value such that
!                   ABS( AB0G(NTERMG) ) < EPS/100

!    NTERMH - INTEGER - No. of terms needed for array AB0H.
!             Recommended value such that
!                   ABS( AB0H(NTERMH) ) < EPS/100

!    NTERMA - INTEGER - No. of terms needed for array AB0AS.
!             Recommended value such that
!                   ABS( AB0AS(NTERMA) ) < EPS/100

!   XLOW1 - REAL (dp) - The value below which
!            ABRAM0 = root(pi)/2 + X ( ln X - GVAL0 )
!           Recommended value is SQRT(2*EPSNEG)

!   LNXMIN - REAL (dp) - The value of ln XMIN. Used to prevent
!            exponential underflow for large X.

!   For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT.

!   The machine-dependent constants are computed internally by
!   using the D1MACH subroutine.


! INTRINSIC FUNCTIONS USED:
!   LOG, EXP, SQRT


! OTHER MISCFUN SUBROUTINES USED:
!        CHEVAL , ERRPRN, D1MACH

! AUTHOR:

!    DR. ALLAN J. MACLEOD,
!    DEPT. OF MATHEMATICS AND STATISTICS,
!    UNIVERSITY OF PAISLEY ,
!    HIGH ST.,
!    PAISLEY,
!    SCOTLAND.
!    PA1 2BE

!    ( e-mail: macl_ms0@paisley.ac.uk )

! LATEST REVISION:   23 January, 1996

REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterma, ntermf, ntermg, ntermh
REAL (dp)  :: asln, asval, fval, gval, hval, lnxmin, t, v, x, xlow1
CHARACTER (LEN=33)  :: errmsg = 'FUNCTION CALLED WITH ARGUMENT < 0'
CHARACTER (LEN=6)   :: fnname = 'ABRAM0'

REAL (dp), PARAMETER  :: ab0f(0:8) = (/  &
   -0.68121927093549469816D0, -0.78867919816149252495D0,  &
    0.5121581776818819543D-1, -0.71092352894541296D-3,  &
    0.368681808504287D-5, -0.917832337237D-8, 0.1270202563D-10,  &
   -0.1076888D-13, 0.599D-17 /)
REAL (dp), PARAMETER  :: ab0g(0:8) = (/  &
   -0.60506039430868273190D0, -0.41950398163201779803D0,  &
    0.1703265125190370333D-1, -0.16938917842491397D-3,  &
    0.67638089519710D-6, -0.135723636255D-8, 0.156297065D-11,  &
   -0.112887D-14, 0.55D-18 /)
REAL (dp), PARAMETER  :: ab0h(0:8) = (/  &
    1.38202655230574989705D0, -0.30097929073974904355D0,  &
    0.794288809364887241D-2, -0.6431910276847563D-4,  &
    0.22549830684374D-6, -0.41220966195D-9, 0.44185282D-12,  &
   -0.30123D-15, 0.14D-18 /)
REAL (dp), PARAMETER  :: ab0as(0:27) = (/  &
    1.97755499723693067407_dp, -0.1046024792004819485D-1,  &
    0.69680790253625366D-3, -0.5898298299996599D-4, 0.577164455305320D-5, &
   -0.61523013365756D-6, 0.6785396884767D-7, -0.723062537907D-8, &
    0.63306627365D-9, -0.989453793D-11, -0.1681980530D-10, 0.673799551D-11, &
   -0.200997939D-11, 0.54055903D-12, -0.13816679D-12, 0.3422205D-13, &
   -0.826686D-14, 0.194566D-14, -0.44268D-15, 0.9562D-16, -0.1883D-16, &
    0.301D-17, -0.19D-18, -0.14D-18, 0.11D-18, -0.4D-19, 0.2D-19, -0.1D-19 /)

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, two = 2.0_dp,  &
                         three = 3.0_dp, six = 6.0_dp, onehun = 100.0_dp,  &
                         rt3bpi = 0.97720502380583984317_dp, &
                         rtpib2 = 0.88622692545275801365_dp, &
                         gval0 = 0.13417650264770070909_dp,  &
                         onerpi = 0.56418958354775628695_dp

!   Start computation

x = xvalue

!   Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

t = 2*EPSILON(zero) / onehun
IF (x <= two) THEN
  DO  ntermf = 8, 0, -1
    IF (ABS(ab0f(ntermf)) > t) EXIT
  END DO
  DO  ntermg = 8, 0, -1
    IF (ABS(ab0g(ntermg)) > t) EXIT
  END DO
  DO  ntermh = 8, 0, -1
    IF (ABS(ab0h(ntermh)) > t) EXIT
  END DO
  xlow1 = SQRT(two*EPSILON(zero))
ELSE
  DO  nterma = 27, 0, -1
    IF (ABS(ab0as(nterma)) > t) EXIT
  END DO
  lnxmin = LOG(TINY(zero))
END IF

!   Code for 0 <= XVALUE <= 2

IF (x <= two) THEN
  IF (x == zero) THEN
    fn_val = rtpib2
    RETURN
  END IF
  IF (x < xlow1) THEN
    fn_val = rtpib2 + x * (LOG(x)-gval0)
    RETURN
  ELSE
    t = (x*x/two-half) - half
    fval = cheval(ntermf,ab0f,t)
    gval = cheval(ntermg,ab0g,t)
    hval = cheval(ntermh,ab0h,t)
    fn_val = fval / onerpi + x * (LOG(x)*hval-gval)
    RETURN
  END IF
ELSE
  
!   Code for XVALUE > 2
  
  v = three * ((x/two)**(two/three))
  t = (six/v-half) - half
  asval = cheval(nterma,ab0as,t)
  asln = LOG(asval/rt3bpi) - v
  IF (asln < lnxmin) THEN
    fn_val = zero
  ELSE
    fn_val = EXP(asln)
  END IF
END IF
RETURN
END FUNCTION abram0



FUNCTION abram1(xvalue) RESULT(fn_val)

!   DESCRIPTION:
!      This function calculates the Abramowitz function of order 1, defined as

!       ABRAM1(x) = integral{ 0 to infinity } t * exp( -t*t - x/t ) dt

!       The code uses Chebyshev expansions with the coefficients
!       given to an accuracy of 20 decimal places.


!   ERROR RETURNS:
!      If XVALUE < 0.0, the function prints a message and returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERMF - INTEGER - No. of terms needed for the array AB1F.
!               Recommended value such that
!                     ABS( AB1F(NTERMF) ) < EPS/100

!      NTERMG - INTEGER - No. of terms needed for array AB1G.
!               Recommended value such that
!                     ABS( AB1G(NTERMG) ) < EPS/100

!      NTERMH - INTEGER - No. of terms needed for array AB1H.
!               Recommended value such that
!                     ABS( AB1H(NTERMH) ) < EPS/100

!      NTERMA - INTEGER - No. of terms needed for array AB1AS.
!               Recommended value such that
!                     ABS( AB1AS(NTERMA) ) < EPS/100

!      XLOW - REAL (dp) - The value below which
!                ABRAM1(x) = 0.5 to machine precision.
!             The recommended value is EPSNEG/2

!      XLOW1 - REAL (dp) - The value below which
!                ABRAM1(x) = (1 - x ( sqrt(pi) + xln(x) ) / 2
!              Recommended value is SQRT(2*EPSNEG)

!      LNXMIN - REAL (dp) - The value of ln XMIN. Used to prevent
!              exponential underflow for large X.

!      For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT

!      The machine-dependent constants are computed internally by using
!      the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!     LOG, EXP, SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:

!      DR. ALLAN J. MACLEOD,
!      DEPT. OF MATHEMATICS AND STATISTICS,
!      UNIVERSITY OF PAISLEY,
!      HIGH ST.,
!      PAISLEY,
!      SCOTLAND.
!      PA1 2BE

!      ( e-mail: macl_ms0@paisley.ac.uk )

!   LATEST REVISION:   23 January, 1996

REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterma, ntermf, ntermg, ntermh
REAL (dp)  :: asln, asval, fval, gval, hval, lnxmin, t, v, x, xlow, xlow1
CHARACTER (LEN=33)  :: errmsg = 'FUNCTION CALLED WITH ARGUMENT < 0'
CHARACTER (LEN=6)   :: fnname = 'ABRAM1'

REAL (dp), PARAMETER  :: ab1f(0:9) = (/  &
    1.47285192577978807369_dp, 0.10903497570168956257_dp,  &
   -0.12430675360056569753_dp, 0.306197946853493315D-2, -0.2218410323076511D-4, 0.6989978834451D-7, -0.11597076444D-9,  &
    0.11389776D-12, -0.7173D-16, 0.3D-19 /)
REAL (dp), PARAMETER  :: ab1g(0:8) = (/  &
    0.39791277949054503528_dp, -0.29045285226454720849_dp,  &
    0.1048784695465363504D-1, -0.10249869522691336D-3, 0.41150279399110D-6,  &
   -0.83652638940D-9, 0.97862595D-12, -0.71868D-15, 0.35D-18 /)
REAL (dp), PARAMETER  :: ab1h(0:8) = (/  &
    0.84150292152274947030_dp, -0.7790050698774143395D-1,  &
    0.133992455878390993D-2, -0.808503907152788D-5, 0.2261858281728D-7,  &
   -0.3441395838D-10, 0.3159858D-13, -0.1884D-16, 0.1D-19 /)
REAL (dp), PARAMETER  :: ab1as(0:27) = (/  &
    2.13013643429065549448_dp, 0.6371526795218539933D-1,  &
   -0.129334917477510647D-2, 0.5678328753228265D-4, -0.279434939177646D-5,  &
    0.5600214736787D-7, 0.2392009242798D-7, -0.750984865009D-8,  &
    0.173015330776D-8, -0.36648877955D-9, 0.7520758307D-10, -0.1517990208D-10, &
    0.301713710D-11, -0.58596718D-12, 0.10914455D-12, -0.1870536D-13,  &
    0.262542D-14, -0.14627D-15, -0.9500D-16, 0.5873D-16, -0.2420D-16,  &
    0.868D-17, -0.290D-17, 0.93D-18, -0.29D-18, 0.9D-19, -0.3D-19, 0.1D-19 /)
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, three = 3.0_dp, six = 6.0_dp,  &
                         onehun = 100.0_dp,  &
                         rt3bpi = 0.97720502380583984317_dp, &
                         onerpi = 0.56418958354775628695_dp

!   Start calculation

x = xvalue

!   Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

t = 2*EPSILON(zero) / onehun
IF (x <= two) THEN
  DO  ntermf = 9, 0, -1
    IF (ABS(ab1f(ntermf)) > t) EXIT
  END DO
  DO  ntermg = 8, 0, -1
    IF (ABS(ab1g(ntermg)) > t) EXIT
  END DO
  DO  ntermh = 8, 0, -1
    IF (ABS(ab1h(ntermh)) > t) EXIT
  END DO
  t = EPSILON(zero)
  xlow1 = SQRT(two*t)
  xlow = t / two
ELSE
  DO  nterma = 27, 0, -1
    IF (ABS(ab1as(nterma)) > t) EXIT
  END DO
  lnxmin = LOG(TINY(zero))
END IF

!   Code for 0 <= XVALUE <= 2

IF (x <= two) THEN
  IF (x == zero) THEN
    fn_val = half
    RETURN
  END IF
  IF (x < xlow1) THEN
    IF (x < xlow) THEN
      fn_val = half
    ELSE
      fn_val = (one-x/onerpi-x*x*LOG(x)) * half
    END IF
    RETURN
  ELSE
    t = (x*x/two-half) - half
    fval = cheval(ntermf,ab1f,t)
    gval = cheval(ntermg,ab1g,t)
    hval = cheval(ntermh,ab1h,t)
    fn_val = fval - x * (gval/onerpi+x*LOG(x)*hval)
    RETURN
  END IF
ELSE
  
!   Code for XVALUE > 2
  
  v = three * ((x/two)**(two/three))
  t = (six/v-half) - half
  asval = cheval(nterma,ab1as,t)
  asln = LOG(asval*SQRT(v/three)/rt3bpi) - v
  IF (asln < lnxmin) THEN
    fn_val = zero
  ELSE
    fn_val = EXP(asln)
  END IF
END IF
RETURN
END FUNCTION abram1



FUNCTION abram2(xvalue) RESULT(fn_val)

!   DESCRIPTION:
!      This function calculates the Abramowitz function of order 2, defined as

!       ABRAM2(x) = integral{ 0 to infinity } (t**2) * exp( -t*t - x/t ) dt

!      The code uses Chebyshev expansions with the coefficients
!      given to an accuracy of 20 decimal places.


!   ERROR RETURNS:
!      If XVALUE < 0.0, the function prints a message and returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERMF - INTEGER - No. of terms needed for the array AB2F.
!               Recommended value such that
!                     ABS( AB2F(NTERMF) ) < EPS/100

!      NTERMG - INTEGER - No. of terms needed for array AB2G.
!               Recommended value such that
!                     ABS( AB2G(NTERMG) ) < EPS/100

!      NTERMH - INTEGER - No. of terms needed for array AB2H.
!               Recommended value such that
!                     ABS( AB2H(NTERMH) ) < EPS/100

!      NTERMA - INTEGER - No. of terms needed for array AB2AS.
!               Recommended value such that
!                     ABS( AB2AS(NTERMA) ) < EPS/100

!      XLOW - REAL (dp) - The value below which
!               ABRAM2 = root(pi)/4 to machine precision.
!             The recommended value is EPSNEG

!      XLOW1 - REAL (dp) - The value below which
!                ABRAM2 = root(pi)/4 - x/2 + x**3ln(x)/6
!              Recommended value is SQRT(2*EPSNEG)

!      LNXMIN - REAL (dp) - The value of ln XMIN. Used to prevent
!               exponential underflow for large X.

!     For values of EPS, EPSNEG, XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!     LOG, EXP

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:

!      DR. ALLAN J. MACLEOD,
!      DEPT. OF MATHEMATICS AND STATISTICS,
!      UNIVERSITY OF PAISLEY,
!      HIGH ST.,
!      PAISLEY,
!      SCOTLAND.
!      PA1 2BE

!      ( e-mail: macl_ms0@paisley.ac.uk )

!   LATEST REVISION:   23 January, 1996

REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterma, ntermf, ntermg, ntermh
REAL (dp)  :: asln, asval, fval, gval, hval, lnxmin, t, v, x, xlow, xlow1
CHARACTER (LEN=33)  :: errmsg = 'FUNCTION CALLED WITH ARGUMENT < 0'
CHARACTER (LEN= 6)  :: fnname = 'ABRAM2'

REAL (dp), PARAMETER  :: ab2f(0:9) = (/  &
    1.03612162804243713846_dp, 0.19371246626794570012_dp,  &
   -0.7258758839233007378D-1, 0.174790590864327399D-2, -0.1281223233756549D-4, &
    0.4115018153651D-7, -0.6971047256D-10, 0.6990183D-13, -0.4492D-16, 0.2D-19 /)
REAL (dp), PARAMETER  :: ab2g(0:8) = (/  &
    1.46290157198630741150D0, 0.20189466883154014317D0,  &
   -0.2908292087997129022D-1, 0.47061049035270050D-3, -0.257922080359333D-5, &
    0.656133712946D-8, -0.914110203D-11, 0.774276D-14, -0.429D-17 /)
REAL (dp), PARAMETER  :: ab2h(0:7) = (/  &
    0.30117225010910488881_dp, -0.1588667818317623783D-1,  &
    0.19295936935584526D-3, -0.90199587849300D-6, 0.206105041837D-8,  &
   -0.265111806D-11, 0.210864D-14, -0.111D-17 /)
REAL (dp), PARAMETER  :: ab2as(0:26) = (/  &
    2.46492325304334856893_dp, 0.23142797422248905432_dp,  &
   -0.94068173010085773D-3, 0.8290270038089733D-4, -0.883894704245866D-5,  &
    0.106638543567985D-5, -0.13991128538529D-6, 0.1939793208445D-7,  &
   -0.277049938375D-8, 0.39590687186D-9, -0.5408354342D-10, 0.635546076D-11, &
   -0.38461613D-12, -0.11696067D-12, 0.6896671D-13, -0.2503113D-13,  &
    0.785586D-14, -0.230334D-14, 0.64914D-15, -0.17797D-15, 0.4766D-16,  &
   -0.1246D-16, 0.316D-17, -0.77D-18, 0.18D-18, -0.4D-19, 0.1D-19 /)
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, two = 2.0_dp,  &
                         three = 3.0_dp, six = 6.0_dp, onehun = 100.0_dp,  &
                         rt3bpi = 0.97720502380583984317_dp,  &
                         rtpib4 = 0.44311346272637900682_dp,  &
                         onerpi = 0.56418958354775628695_dp

!   Start calculation

x = xvalue

!   Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

t = 2*EPSILON(zero) / onehun
IF (x <= two) THEN
  DO  ntermf = 9, 0, -1
    IF (ABS(ab2f(ntermf)) > t) EXIT
  END DO
  DO  ntermg = 8, 0, -1
    IF (ABS(ab2g(ntermg)) > t) EXIT
  END DO
  DO  ntermh = 7, 0, -1
    IF (ABS(ab2h(ntermh)) > t) EXIT
  END DO
  xlow = EPSILON(zero)
  xlow1 = SQRT(two*xlow)
ELSE
  DO  nterma = 26, 0, -1
    IF (ABS(ab2as(nterma)) > t) EXIT
  END DO
  lnxmin = LOG(TINY(zero))
END IF

!   Code for 0 <= XVALUE <= 2

IF (x <= two) THEN
  IF (x == zero) THEN
    fn_val = rtpib4
    RETURN
  END IF
  IF (x < xlow1) THEN
    IF (x < xlow) THEN
      fn_val = rtpib4
    ELSE
      fn_val = rtpib4 - half * x + x * x * x * LOG(x) / six
    END IF
    RETURN
  ELSE
    t = (x*x/two-half) - half
    fval = cheval(ntermf,ab2f,t)
    gval = cheval(ntermg,ab2g,t)
    hval = cheval(ntermh,ab2h,t)
    fn_val = fval / onerpi + x * (x*x*LOG(x)*hval-gval)
    RETURN
  END IF
ELSE
  
!   Code for XVALUE > 2
  
  v = three * ((x/two)**(two/three))
  t = (six/v-half) - half
  asval = cheval(nterma,ab2as,t)
  asln = LOG(asval/rt3bpi) + LOG(v/three) - v
  IF (asln < lnxmin) THEN
    fn_val = zero
  ELSE
    fn_val = EXP(asln)
  END IF
END IF
RETURN
END FUNCTION abram2



FUNCTION airint(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the integral of the Airy function Ai, defined as

!         AIRINT(x) = {integral 0 to x} Ai(t) dt

!      The program uses Chebyshev expansions, the coefficients of which
!      are given to 20 decimal places.


!   ERROR RETURNS:

!      If the argument is too large and negative, it is impossible
!      to accurately compute the necessary SIN and COS functions.
!      An error message is printed, and the program returns the
!      value -2/3 (the value at -infinity).


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - INTEGER - The no. of terms to be used from the array
!                          AAINT1. The recommended value is such that
!                             ABS(AAINT1(NTERM1)) < EPS/100,
!                          subject to 1 <= NTERM1 <= 25.

!      NTERM2 - INTEGER - The no. of terms to be used from the array
!                          AAINT2. The recommended value is such that
!                             ABS(AAINT2(NTERM2)) < EPS/100,
!                          subject to 1 <= NTERM2 <= 21.

!      NTERM3 - INTEGER - The no. of terms to be used from the array
!                          AAINT3. The recommended value is such that
!                             ABS(AAINT3(NTERM3)) < EPS/100,
!                          subject to 1 <= NTERM3 <= 40.

!      NTERM4 - INTEGER - The no. of terms to be used from the array
!                          AAINT4. The recommended value is such that
!                             ABS(AAINT4(NTERM4)) < EPS/100,
!                          subject to 1 <= NTERM4 <= 17.

!      NTERM5 - INTEGER - The no. of terms to be used from the array
!                          AAINT5. The recommended value is such that
!                             ABS(AAINT5(NTERM5)) < EPS/100,
!                          subject to 1 <= NTERM5 <= 17.

!      XLOW1 - REAL (dp) - The value such that, if |x| < XLOW1,
!                          AIRINT(x) = x * Ai(0)
!                     to machine precision. The recommended value is
!                          2 * EPSNEG.

!      XHIGH1 - REAL (dp) - The value such that, if x > XHIGH1,
!                          AIRINT(x) = 1/3,
!                      to machine precision. The recommended value is
!                          (-1.5*LOG(EPSNEG)) ** (2/3).

!      XNEG1 - REAL (dp) - The value such that, if x < XNEG1,
!                     the trigonometric functions in the asymptotic
!                     expansion cannot be calculated accurately.
!                     The recommended value is
!                          -(1/((EPS)**2/3))

!      For values of EPS and EPSNEG, refer to the file MACHCON.TXT.

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!                            COS, EXP, SIN, SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR: Dr. Allan J. MacLeod,
!           Dept. of Mathematics and Statistics,
!           Univ. of Paisley,
!           High St.,
!           Paisley,
!           SCOTLAND.
!           PA1 2BE

!           (e-mail:macl_ms0@paisley.ac.uk)

!   LATEST REVISION:  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3, nterm4, nterm5
REAL (dp)  :: arg, gval, hval, t, temp, x, xhigh1, xlow1, xneg1, z
CHARACTER (LEN=46)  :: errmsg = 'AIRINT'
CHARACTER (LEN= 6)  :: fnname = 'FUNCTION TOO NEGATIVE FOR ACCURATE COMPUTATION'

REAL (dp), PARAMETER  :: aaint1(0:25) = (/  &
    0.37713517694683695526_dp, -0.13318868432407947431_dp,  &
    0.3152497374782884809D-1, -0.318543076436574077D-2,  &
   -0.87398764698621915D-3, 0.46699497655396971D-3, -0.9544936738983692D-4,  &
    0.542705687156716D-5, 0.239496406252188D-5, -0.75690270205649D-6,  &
    0.9050138584518D-7, 0.320529456043D-8, -0.303825536444D-8,  &
    0.48900118596D-9, -0.1839820572D-10, -0.711247519D-11, 0.151774419D-11,  &
   -0.10801922D-12, -0.963542D-14, 0.313425D-14, -0.29446D-15, -0.477D-17,  &
    0.461D-17, -0.53D-18, 0.1D-19, 0.1D-19 /)

REAL (dp), PARAMETER  :: aaint2(0:21) = (/  &
    1.92002524081984009769_dp, -0.4220049417256287021D-1,  &
   -0.239457722965939223D-2, -0.19564070483352971D-3,  &
   -0.1547252891056112D-4, -0.140490186137889D-5, -0.12128014271367D-6,  &
   -0.1179186050192D-7, -0.104315578788D-8, -0.10908209293D-9,  &
   -0.929633045D-11, -0.110946520D-11, -0.7816483D-13, -0.1319661D-13,  &
   -0.36823D-15, -0.21505D-15, 0.1238D-16, -0.557D-17, 0.84D-18, -0.21D-18,  &
    0.4D-19, -0.1D-19 /)

REAL (dp), PARAMETER  :: aaint3(0:40) = (/  &
    0.47985893264791052053_dp, -0.19272375126169608863_dp,  &
    0.2051154129525428189D-1, 0.6332000070732488786D-1,  &
   -0.5093322261845754082D-1, 0.1284424078661663016D-1,  &
    0.2760137088989479413D-1, -0.1547066673866649507D-1,  &
   -0.1496864655389316026D-1, 0.336617614173574541D-2,  &
    0.530851163518892985D-2, 0.41371226458555081D-3, &
   -0.102490579926726266D-2, -0.32508221672025853D-3,  &
    0.8608660957169213D-4, 0.6671367298120775D-4, 0.449205999318095D-5,  &
   -0.670427230958249D-5, -0.196636570085009D-5, 0.22229677407226D-6,  &
    0.22332222949137D-6, 0.2803313766457D-7, -0.1155651663619D-7,  &
   -0.433069821736D-8, -0.6227777938D-10, 0.26432664903D-9, 0.5333881114D-10, &
   -0.522957269D-11, -0.382229283D-11, -0.40958233D-12, 0.11515622D-12,  &
    0.3875766D-13, 0.140283D-14, -0.141526D-14, -0.28746D-15, 0.923D-17,  &
    0.1224D-16, 0.157D-17, -0.19D-18, -0.8D-19, -0.1D-19 /)

REAL (dp), PARAMETER  :: aaint4(0:17) = (/  &
    1.99653305828522730048D0, -0.187541177605417759D-2,  &
   -0.15377536280305750D-3, -0.1283112967682349D-4, -0.108128481964162D-5, &
   -0.9182131174057D-7, -0.784160590960D-8, -0.67292453878D-9,  &
   -0.5796325198D-10, -0.501040991D-11, -0.43420222D-12, -0.3774305D-13,  &
   -0.328473D-14, -0.28700D-15, -0.2502D-16, -0.220D-17, -0.19D-18, -0.2D-19 /)
REAL (dp), PARAMETER  :: aaint5(0:17) = (/  &
    1.13024602034465716133D0, -0.464718064639872334D-2,  &
   -0.35137413382693203D-3, -0.2768117872545185D-4, -0.222057452558107D-5,  &
   -0.18089142365974D-6, -0.1487613383373D-7, -0.123515388168D-8,  &
   -0.10310104257D-9, -0.867493013D-11, -0.73080054D-12, -0.6223561D-13,   &
   -0.525128D-14, -0.45677D-15, -0.3748D-16, -0.356D-17, -0.23D-18, -0.4D-19 /)

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp,  &
                         three = 3.0_dp, four = 4.0_dp, eight = 8.0_dp,  &
                         nine = 9.0_dp, forty1 = 41.0_dp, onehun = 100.0_dp, &
                         ninhun = 900.0_dp, fr996 = 4996.0_dp,  &
                         piby4 = 0.78539816339744830962_dp,  &
                         pitim6 = 18.84955592153875943078_dp,  &
                         rt2b3p = 0.46065886596178063902_dp,  &
                         airzer = 0.35502805388781723926_dp

!   Start computation

x = xvalue

!   Compute the machine-dependent constants.

z = EPSILON(zero)
xlow1 = two * z
arg = 2*EPSILON(zero)
xneg1 = -one / (arg**(two/three))

!   Error test

IF (x < xneg1) THEN
  CALL errprn(fnname,errmsg)
  fn_val = -two / three
  RETURN
END IF

!  continue with machine-dependent constants

t = arg / onehun
IF (x >= zero) THEN
  DO  nterm1 = 25, 0, -1
    IF (ABS(aaint1(nterm1)) > t) EXIT
  END DO
  DO  nterm2 = 21, 0, -1
    IF (ABS(aaint2(nterm2)) > t) EXIT
  END DO
  xhigh1 = (-three*LOG(z)/two) ** (two/three)
ELSE
  DO  nterm3 = 40, 0, -1
    IF (ABS(aaint3(nterm3)) > t) EXIT
  END DO
  DO  nterm4 = 17, 0, -1
    IF (ABS(aaint4(nterm4)) > t) EXIT
  END DO
  DO  nterm5 = 17, 0, -1
    IF (ABS(aaint5(nterm5)) > t) EXIT
  END DO
END IF

!   Code for x >= 0

IF (x >= zero) THEN
  IF (x <= four) THEN
    IF (x < xlow1) THEN
      fn_val = airzer * x
    ELSE
      t = x / two - one
      fn_val = cheval(nterm1,aaint1,t) * x
    END IF
  ELSE
    IF (x > xhigh1) THEN
      temp = zero
    ELSE
      z = (x+x) * SQRT(x) / three
      temp = three * z
      t = (forty1-temp) / (nine+temp)
      temp = EXP(-z) * cheval(nterm2,aaint2,t) / SQRT(pitim6*z)
    END IF
    fn_val = one / three - temp
  END IF
ELSE
  
!   Code for x < 0
  
  IF (x >= -eight) THEN
    IF (x > -xlow1) THEN
      fn_val = airzer * x
    ELSE
      t = -x / four - one
      fn_val = x * cheval(nterm3,aaint3,t)
    END IF
  ELSE
    z = -(x+x) * SQRT(-x) / three
    arg = z + piby4
    temp = nine * z * z
    t = (fr996-temp) / (ninhun+temp)
    gval = cheval(nterm4,aaint4,t)
    hval = cheval(nterm5,aaint5,t)
    temp = gval * COS(arg) + hval * SIN(arg) / z
    fn_val = rt2b3p * temp / SQRT(z) - two / three
  END IF
END IF
RETURN
END FUNCTION airint



FUNCTION airygi(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This subroutine computes the modified Airy function Gi(x), defined as

!        AIRYGI(x) = [ Integral{0 to infinity} sin(x*t+t^3/3) dt ] / pi

!      The approximation uses Chebyshev expansions with the coefficients
!      given to 20 decimal places.


!   ERROR RETURNS:

!      If x < -XHIGH1*XHIGH1 (see below for definition of XHIGH1), then
!      the trig. functions needed for the asymptotic expansion of Bi(x)
!      cannot be computed to any accuracy. An error message is printed
!      and the code returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - INTEGER - The no. of terms to be used from the array
!                         ARGIP1. The recommended value is such that
!                                ABS(ARGIP1(NTERM1)) < EPS/100
!                         subject to 1 <= NTERM1 <= 30.

!      NTERM2 - INTEGER - The no. of terms to be used from the array
!                         ARGIP2. The recommended value is such that
!                                ABS(ARGIP2(NTERM2)) < EPS/100
!                         subject to 1 <= NTERM2 <= 29.

!      NTERM3 - INTEGER - The no. of terms to be used from the array
!                         ARGIN1. The recommended value is such that
!                                ABS(ARGIN1(NTERM3)) < EPS/100
!                         subject to 1 <= NTERM3 <= 42.

!      NTERM4 - INTEGER - The no. of terms to be used from the array
!                         ARBIN1. The recommended value is such that
!                                ABS(ARBIN1(NTERM4)) < EPS/100
!                         subject to 1 <= NTERM4 <= 10.

!      NTERM5 - INTEGER - The no. of terms to be used from the array
!                         ARBIN2. The recommended value is such that
!                                ABS(ARBIN2(NTERM5)) < EPS/100
!                         subject to 1 <= NTERM5 <= 11.

!      NTERM6 - INTEGER - The no. of terms to be used from the array
!                         ARGH2. The recommended value is such that
!                                ABS(ARHIN1(NTERM6)) < EPS/100
!                         subject to 1 <= NTERM6 <= 15.

!      XLOW1 - REAL (dp) - The value such that, if -XLOW1 < x < XLOW1,
!                     then AIRYGI = Gi(0) to machine precision.
!                     The recommended value is   EPS.

!      XHIGH1 - REAL (dp) - The value such that, if x > XHIGH1, then
!                      AIRYGI = 1/(Pi*x) to machine precision.
!                      Also used for error test - see above.
!                      The recommended value is
!                          cube root( 2/EPS ).

!      XHIGH2 - REAL (dp) - The value above which AIRYGI = 0.0.
!                      The recommended value is
!                          1/(Pi*XMIN).

!      XHIGH3 - REAL (dp) - The value such that, if x < XHIGH3,
!                      then the Chebyshev expansions for the
!                      asymptotic form of Bi(x) are not needed.
!                      The recommended value is
!                          -8 * cube root( 2/EPSNEG ).

!      For values of EPS, EPSNEG, and XMIN refer to the file
!      MACHCON.TXT.

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!                             COS , SIN , SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. Macleod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley,
!          High St.,
!          Paisley,
!          SCOTLAND.

!          (e-mail: macl_ms0@paisley.ac.uk)


!   LATEST UPDATE:
!                 23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3, nterm4, nterm5, nterm6
REAL (dp)  :: arg, bi, cheb1, cheb2, cosz, sinz, t, temp, x, xcube,  &
              xhigh1, xhigh2, xhigh3, xlow1, xminus, z, zeta
CHARACTER (LEN=46) :: errmsg = 'ARGUMENT TOO NEGATIVE FOR ACCURATE COMPUTATION'
CHARACTER (LEN= 6) :: fnname = 'AIRYGI'

REAL (dp), PARAMETER  :: argip1(0:30) = (/  &
    0.26585770795022745082_dp, -0.10500333097501922907_dp,  &
    0.841347475328454492D-2, 0.2021067387813439541D-1,  &
   -0.1559576113863552234D-1, 0.564342939043256481D-2,  &
   -0.59776844826655809D-3, -0.42833850264867728D-3, 0.22605662380909027D-3,  &
   -0.3608332945592260D-4, -0.785518988788901D-5, 0.473252480746370D-5,  &
   -0.59743513977694D-6, -0.15917609165602D-6, 0.6336129065570D-7,  &
   -0.276090232648D-8, -0.256064154085D-8, 0.47798676856D-9, 0.4488131863D-10, &
   -0.2346508882D-10, 0.76839085D-12, 0.73227985D-12, -0.8513687D-13,  &
   -0.1630201D-13, 0.356769D-14, 0.25001D-15, -0.10859D-15, -0.158D-17,  &
    0.275D-17, -0.5D-19, -0.6D-19 /)

REAL (dp), PARAMETER  :: argip2(0:29) = (/  &
    2.00473712275801486391_dp, 0.294184139364406724D-2,  &
    0.71369249006340167D-3, 0.17526563430502267D-3, 0.4359182094029882D-4,  &
    0.1092626947604307D-4, 0.272382418399029D-5, 0.66230900947687D-6,  &
    0.15425323370315D-6, 0.3418465242306D-7, 0.728157724894D-8,  &
    0.151588525452D-8, 0.30940048039D-9, 0.6149672614D-10, 0.1202877045D-10,  &
    0.233690586D-11, 0.43778068D-12, 0.7996447D-13, 0.1494075D-13,  &
    0.246790D-14, 0.37672D-15, 0.7701D-16, 0.354D-17, -0.49D-18, 0.62D-18,  &
   -0.40D-18, -0.1D-19, 0.2D-19, -0.3D-19, 0.1D-19 /)

REAL (dp), PARAMETER  :: argin1(0:42) = (/  &
   -0.20118965056732089130_dp, -0.7244175303324530499D-1,  &
    0.4505018923894780120D-1, -0.24221371122078791099_dp,  &
    0.2717884964361678294D-1, -0.5729321004818179697D-1,  &
   -0.18382107860337763587_dp, 0.7751546082149475511D-1,  &
    0.18386564733927560387_dp, 0.2921504250185567173D-1,  &
   -0.6142294846788018811D-1, -0.2999312505794616238D-1,  &
    0.585937118327706636D-2, 0.822221658497402529D-2, 0.132579817166846893D-2, &
   -0.96248310766565126D-3, -0.45065515998211807D-3, 0.772423474325474D-5,  &
    0.5481874134758052D-4, 0.1245898039742876D-4, -0.246196891092083D-5,  &
   -0.169154183545285D-5, -0.16769153169442D-6, 0.9636509337672D-7,  &
    0.3253314928030D-7, 0.5091804231D-10, -0.209180453553D-8,  &
   -0.41237387870D-9, 0.4163338253D-10, 0.3032532117D-10, 0.340580529D-11,  &
   -0.88444592D-12, -0.31639612D-12, -0.1505076D-13, 0.1104148D-13,  &
    0.246508D-14, -0.3107D-16, -0.9851D-16, -0.1453D-16, 0.118D-17,  &
    0.67D-18, 0.6D-19, -0.1D-19 /)

REAL (dp), PARAMETER  :: arbin1(0:10) = (/  &
    1.99983763583586155980D0, -0.8104660923669418D-4,  &
    0.13475665984689D-6, -0.70855847143D-9, 0.748184187D-11,  &
   -0.12902774D-12, 0.322504D-14, -0.10809D-15, 0.460D-17, -0.24D-18, 0.1D-19 /)
REAL (dp), PARAMETER  :: arbin2(0:11) = (/  &
    0.13872356453879120276D0, -0.8239286225558228D-4,  &
    0.26720919509866D-6, -0.207423685368D-8, 0.2873392593D-10,  &
   -0.60873521D-12, 0.1792489D-13, -0.68760D-15, 0.3280D-16,  &
   -0.188D-17, 0.13D-18, -0.1D-19 /)
REAL (dp), PARAMETER  :: arhin1(0:15) = (/  &
    1.99647720399779650525D0, -0.187563779407173213D-2,  &
   -0.12186470897787339D-3, -0.814021609659287D-5, -0.55050925953537D-6,  &
   -0.3763008043303D-7, -0.258858362365D-8, -0.17931829265D-9,  &
   -0.1245916873D-10, -0.87171247D-12, -0.6084943D-13, -0.431178D-14,  &
   -0.29787D-15, -0.2210D-16, -0.136D-17, -0.14D-18 /)

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, three = 3.0_dp,  &
                         four = 4.0_dp, five = 5.0_dp, seven = 7.0_dp, &
                         minate = -8.0_dp, nine = 9.0_dp, twent8 = 28.0_dp, &
                         seven2 = 72.0_dp, onehun = 100.0_dp, one76 = 176.0_dp, &
                         five14 = 514.0_dp, one024 = 1024.0_dp,  &
                         twelhu = 1200.0_dp,  &
                         gizero = 0.20497554248200024505_dp,  &
                         onebpi = 0.31830988618379067154_dp,  &
                         piby4 = 0.78539816339744830962_dp,   &
                         rtpiin = 0.56418958354775628695_dp

!   Start computation

x = xvalue

!   Compute the machine-dependent constants.

z = EPSILON(zero)
xlow1 = z
arg = 2*EPSILON(zero)
xhigh1 = one / arg
xhigh1 = (xhigh1+xhigh1) ** (one/three)

!   Error test

IF (x < -xhigh1*xhigh1) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!  continue with machine-dependent constants

t = arg / onehun
IF (x >= zero) THEN
  DO  nterm1 = 30, 0, -1
    IF (ABS(argip1(nterm1)) > t) EXIT
  END DO
  DO  nterm2 = 29, 0, -1
    IF (ABS(argip2(nterm2)) > t) EXIT
  END DO
  temp = four * piby4
  xhigh2 = one / (temp*TINY(zero))
ELSE
  DO  nterm3 = 42, 0, -1
    IF (ABS(argin1(nterm3)) > t) EXIT
  END DO
  DO  nterm4 = 10, 0, -1
    IF (ABS(arbin1(nterm4)) > t) EXIT
  END DO
  DO  nterm5 = 11, 0, -1
    IF (ABS(arbin2(nterm5)) > t) EXIT
  END DO
  DO  nterm6 = 15, 0, -1
    IF (ABS(arhin1(nterm6)) > t) EXIT
  END DO
  temp = one / z
  xhigh3 = minate * (temp+temp) ** (one/three)
END IF

!   Code for x >= 0.0

IF (x >= zero) THEN
  IF (x <= seven) THEN
    IF (x < xlow1) THEN
      fn_val = gizero
    ELSE
      t = (nine*x-twent8) / (x+twent8)
      fn_val = cheval(nterm1,argip1,t)
    END IF
  ELSE
    IF (x > xhigh1) THEN
      IF (x > xhigh2) THEN
        fn_val = zero
      ELSE
        fn_val = onebpi / x
      END IF
    ELSE
      xcube = x * x * x
      t = (twelhu-xcube) / (five14+xcube)
      fn_val = onebpi * cheval(nterm2,argip2,t) / x
    END IF
  END IF
ELSE
  
!   Code for x < 0.0
  
  IF (x >= minate) THEN
    IF (x > -xlow1) THEN
      fn_val = gizero
    ELSE
      t = -(x+four) / four
      fn_val = cheval(nterm3,argin1,t)
    END IF
  ELSE
    xminus = -x
    t = xminus * SQRT(xminus)
    zeta = (t+t) / three
    temp = rtpiin / SQRT(SQRT(xminus))
    cosz = COS(zeta+piby4)
    sinz = SIN(zeta+piby4) / zeta
    xcube = x * x * x
    IF (x > xhigh3) THEN
      t = -(one024/(xcube)+one)
      cheb1 = cheval(nterm4,arbin1,t)
      cheb2 = cheval(nterm5,arbin2,t)
      bi = (cosz*cheb1+sinz*cheb2) * temp
    ELSE
      bi = (cosz+sinz*five/seven2) * temp
    END IF
    t = (xcube+twelhu) / (one76-xcube)
    fn_val = bi + cheval(nterm6,arhin1,t) * onebpi / x
  END IF
END IF
RETURN
END FUNCTION airygi



FUNCTION airyhi(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This subroutine computes the modified Airy function Hi(x), defined as

!         AIRYHI(x) = [ Integral{0 to infinity} exp(x*t-t^3/3) dt ] / pi

!      The approximation uses Chebyshev expansions with the coefficients
!      given to 20 decimal places.


!   ERROR RETURNS:

!      If x > XHIGH1 (see below for definition of XHIGH1), then
!      the asymptotic expansion of Hi(x) will cause an overflow.
!      An error message is printed and the code returns the largest
!      floating-pt number as the result.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - INTEGER - The no. of terms to be used from the array
!                         ARHIP. The recommended value is such that
!                                ABS(ARHIP(NTERM1)) < EPS/100
!                         subject to 1 <= NTERM1 <= 31.

!      NTERM2 - INTEGER - The no. of terms to be used from the array
!                         ARBIP. The recommended value is such that
!                                ABS(ARBIP(NTERM2)) < EPS/100
!                         subject to 1 <= NTERM2 <= 23.

!      NTERM3 - INTEGER - The no. of terms to be used from the array
!                         ARGIP. The recommended value is such that
!                                ABS(ARGIP1(NTERM3)) < EPS/100
!                         subject to 1 <= NTERM3 <= 29.

!      NTERM4 - INTEGER - The no. of terms to be used from the array
!                         ARHIN1. The recommended value is such that
!                                ABS(ARHIN1(NTERM4)) < EPS/100
!                         subject to 1 <= NTERM4 <= 21.

!      NTERM5 - INTEGER - The no. of terms to be used from the array
!                         ARHIN2. The recommended value is such that
!                                ABS(ARHIN2(NTERM5)) < EPS/100
!                         subject to 1 <= NTERM5 <= 15.

!      XLOW1 - REAL (dp) - The value such that, if -XLOW1 < x < XLOW1,
!                     then AIRYGI = Hi(0) to machine precision.
!                     The recommended value is   EPS.

!      XHIGH1 - REAL (dp) - The value such that, if x > XHIGH1, then
!                      overflow might occur. The recommended value is
!                      computed as follows:
!                           compute Z = 1.5*LOG(XMAX)
!                        XHIGH1 = ( Z + LOG(Z)/4 + LOG(PI)/2 )**(2/3)

!      XNEG1 - REAL (dp) - The value below which AIRYHI = 0.0.
!                     The recommended value is
!                          -1/(Pi*XMIN).

!      XNEG2 - REAL (dp) - The value such that, if x < XNEG2, then
!                      AIRYHI = -1/(Pi*x) to machine precision.
!                      The recommended value is
!                          -cube root( 2/EPS ).

!      XMAX - REAL (dp) - The largest possible floating-pt. number.
!                    This is the value given to the function
!                    if x > XHIGH1.

!      For values of EPS, EPSNEG, XMIN  and XMAX refer to the file
!      MACHCON.TXT.

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!                            EXP , LOG , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. Macleod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley,
!          High St.,
!          Paisley,
!          SCOTLAND.

!          (e-mail: macl_ms0@paisley.ac.uk)


!   LATEST UPDATE:
!                  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3, nterm4, nterm5
REAL (dp)  :: bi, gi, t, temp, x, xcube, xhigh1,  &
              xlow1, xmax, xneg1, xneg2, z, zeta
CHARACTER (LEN=30)  :: errmsg = 'ARGUMENT TO FUNCTION TOO LARGE'
CHARACTER (LEN= 6)  :: fnname = 'AIRYHI'

REAL (dp), PARAMETER  :: arhip(0:31) = (/  &
    1.24013562561762831114_dp, 0.64856341973926535804_dp,  &
    0.55236252592114903246_dp, 0.20975122073857566794_dp,  &
    0.12025669118052373568_dp, 0.3768224931095393785D-1,  &
    0.1651088671548071651D-1, 0.455922755211570993D-2,  &
    0.161828480477635013D-2, 0.40841282508126663D-3, 0.12196479721394051D-3, &
    0.2865064098657610D-4, 0.742221556424344D-5, 0.163536231932831D-5,  &
    0.37713908188749D-6, 0.7815800336008D-7, 0.1638447121370D-7,  &
    0.319857665992D-8, 0.61933905307D-9, 0.11411161191D-9, 0.2064923454D-10, &
    0.360018664D-11, 0.61401849D-12, 0.10162125D-12, 0.1643701D-13,  &
    0.259084D-14, 0.39931D-15, 0.6014D-16, 0.886D-17, 0.128D-17, 0.18D-18, &
    0.3D-19 /)
REAL (dp), PARAMETER  :: arbip(0:23) = (/  &
    2.00582138209759064905_dp, 0.294478449170441549D-2,  &
    0.3489754514775355D-4, 0.83389733374343D-6, 0.3136215471813D-7,  &
    0.167865306015D-8, 0.12217934059D-9, 0.1191584139D-10, 0.154142553D-11,   &
    0.24844455D-12, 0.4213012D-13, 0.505293D-14, -0.60032D-15, -0.65474D-15,  &
   -0.22364D-15, -0.3015D-16, 0.959D-17, 0.616D-17, 0.97D-18, -0.37D-18,  &
   -0.21D-18, -0.1D-19, 0.2D-19, 0.1D-19 /)
REAL (dp), PARAMETER  :: argip1(0:29) = (/  &
    2.00473712275801486391_dp, 0.294184139364406724D-2, 0.71369249006340167D-3, &
    0.17526563430502267D-3, 0.4359182094029882D-4, 0.1092626947604307D-4,  &
    0.272382418399029D-5, 0.66230900947687D-6, 0.15425323370315D-6,  &
    0.3418465242306D-7, 0.728157724894D-8, 0.151588525452D-8, 0.30940048039D-9, &
    0.6149672614D-10, 0.1202877045D-10, 0.233690586D-11, 0.43778068D-12,  &
    0.7996447D-13, 0.1494075D-13, 0.246790D-14, 0.37672D-15, 0.7701D-16,  &
    0.354D-17, -0.49D-18, 0.62D-18, -0.40D-18, -0.1D-19, 0.2D-19, -0.3D-19,  &
    0.1D-19 /)
REAL (dp), PARAMETER  :: arhin1(0:21) = (/  &
    0.31481017206423404116_dp, -0.16414499216588964341_dp,  &
    0.6176651597730913071D-1, -0.1971881185935933028D-1,   &
    0.536902830023331343D-2, -0.124977068439663038D-2, 0.24835515596994933D-3, &
   -0.4187024096746630D-4, 0.590945437979124D-5, -0.68063541184345D-6,  &
    0.6072897629164D-7, -0.367130349242D-8, 0.7078017552D-10,  &
    0.1187894334D-10, -0.120898723D-11, 0.1189656D-13, 0.594128D-14,  &
   -0.32257D-15, -0.2290D-16, 0.253D-17, 0.9D-19, -0.2D-19 /)
REAL (dp), PARAMETER  :: arhin2(0:15) = (/  &
    1.99647720399779650525_dp, -0.187563779407173213D-2,  &
   -0.12186470897787339D-3, -0.814021609659287D-5, -0.55050925953537D-6,  &
   -0.3763008043303D-7, -0.258858362365D-8, -0.17931829265D-9,  &
   -0.1245916873D-10, -0.87171247D-12, -0.6084943D-13, -0.431178D-14,  &
   -0.29787D-15, -0.2210D-16, -0.136D-17, -0.14D-18 /)
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp,  &
                         three = 3.0_dp, four = 4.0_dp, seven = 7.0_dp,  &
                         minate = -8.0_dp, twelve = 12.0_dp, one76 = 176.0_dp, &
                         thre43 = 343.0_dp, five14 = 514.0_dp,  &
                         twelhu = 1200.0_dp, onehun = 100.0_dp,  &
                         hizero = 0.40995108496400049010_dp,  &
                         lnrtpi = 0.57236494292470008707_dp,  &
                         onebpi = 0.31830988618379067154_dp

!   Start computation

x = xvalue

!   Compute the machine-dependent constants.

xmax = HUGE(zero)
temp = three * LOG(xmax) / two
zeta = (temp + LOG(temp)/four - LOG(onebpi)/two)
xhigh1 = zeta ** (two/three)

!   Error test

IF (x > xhigh1) THEN
  CALL errprn(fnname,errmsg)
  fn_val = xmax
  RETURN
END IF

!  continue with machine-dependent constants

z = EPSILON(zero)
xlow1 = z
t = z / onehun
IF (x >= zero) THEN
  DO  nterm1 = 31, 0, -1
    IF (ABS(arhip(nterm1)) > t) EXIT
  END DO
  DO  nterm2 = 23, 0, -1
    IF (ABS(arbip(nterm2)) > t) EXIT
  END DO
  DO  nterm3 = 29, 0, -1
    IF (ABS(argip1(nterm3)) > t) EXIT
  END DO
ELSE
  DO  nterm4 = 21, 0, -1
    IF (ABS(arhin1(nterm4)) > t) EXIT
  END DO
  DO  nterm5 = 15, 0, -1
    IF (ABS(arhin2(nterm5)) > t) EXIT
  END DO
  temp = one / onebpi
  xneg1 = -one / (temp*TINY(zero))
  xneg2 = -((two/z)**(one/three))
END IF

!   Code for x >= 0.0

IF (x >= zero) THEN
  IF (x <= seven) THEN
    IF (x < xlow1) THEN
      fn_val = hizero
    ELSE
      t = (x+x) / seven - one
      temp = (x+x+x) / two
      fn_val = EXP(temp) * cheval(nterm1,arhip,t)
    END IF
  ELSE
    xcube = x * x * x
    temp = SQRT(xcube)
    zeta = (temp+temp) / three
    t = two * (SQRT(thre43/xcube)) - one
    temp = cheval(nterm2,arbip,t)
    temp = zeta + LOG(temp) - LOG(x) / four - lnrtpi
    bi = EXP(temp)
    t = (twelhu-xcube) / (xcube+five14)
    gi = cheval(nterm3,argip1,t) * onebpi / x
    fn_val = bi - gi
  END IF
ELSE
  
!   Code for x < 0.0
  
  IF (x >= minate) THEN
    IF (x > -xlow1) THEN
      fn_val = hizero
    ELSE
      t = (four*x+twelve) / (x-twelve)
      fn_val = cheval(nterm4,arhin1,t)
    END IF
  ELSE
    IF (x < xneg1) THEN
      fn_val = zero
    ELSE
      IF (x < xneg2) THEN
        temp = one
      ELSE
        xcube = x * x * x
        t = (xcube+twelhu) / (one76-xcube)
        temp = cheval(nterm5,arhin2,t)
      END IF
      fn_val = -temp * onebpi / x
    END IF
  END IF
END IF
RETURN
END FUNCTION airyhi



FUNCTION atnint(xvalue) RESULT(fn_val)

! DESCRIPTION:

!   The function ATNINT calculates the value of the
!   inverse-tangent integral defined by

!       ATNINT(x) = integral 0 to x ( (arctan t)/t ) dt

!   The approximation uses Chebyshev series with the coefficients
!   given to an accuracy of 20D.


! ERROR RETURNS:

!   There are no error returns from this program.


! MACHINE-DEPENDENT CONSTANTS:

!   NTERMS - INTEGER - The no. of terms of the array ATNINTT.
!                      The recommended value is such that
!                          ATNINA(NTERMS) < EPS/100

!   XLOW - REAL (dp) - A bound below which ATNINT(x) = x to machine
!                 precision. The recommended value is
!                     sqrt(EPSNEG/2).

!   XUPPER - REAL (dp) - A bound on x, above which, to machine precision
!                   ATNINT(x) = (pi/2)ln x
!                   The recommended value is 1/EPS.

!     For values of EPSNEG and EPS for various machine/compiler
!     combinations refer to the text file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.

! INTRINSIC FUNCTIONS USED:
!    ABS , LOG


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL ,  D1MACH


! AUTHOR: Dr. Allan J. MacLeod,
!         Dept. of Mathematics and Statistics,
!         University of Paisley,
!         High St.,
!         PAISLEY
!         SCOTLAND

!         (e-mail  macl_ms0@paisley.ac.uk )


! LATEST MODIFICATION:  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: ind, nterms
REAL (dp)  :: t, x, xlow, xupper

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         onehun = 100.0_dp, twobpi = 0.63661977236758134308_dp

REAL (dp), PARAMETER  :: atnina(0:22) = (/  &
    1.91040361296235937512_dp, -0.4176351437656746940D-1,  &
    0.275392550786367434D-2, -0.25051809526248881D-3, 0.2666981285121171D-4, &
   -0.311890514107001D-5, 0.38833853132249D-6, -0.5057274584964D-7,  &
    0.681225282949D-8, -0.94212561654D-9, 0.13307878816D-9,  &
   -0.1912678075D-10, 0.278912620D-11, -0.41174820D-12, 0.6142987D-13,  &
   -0.924929D-14, 0.140387D-14, -0.21460D-15, 0.3301D-16, -0.511D-17,  &
    0.79D-18, -0.12D-18, 0.2D-19 /)

!   Compute the machine-dependent constants.

t = 2*EPSILON(zero) / onehun
DO  nterms = 22, 0, -1
  IF (ABS(atnina(nterms)) > t) EXIT
END DO
t = EPSILON(zero)
xlow = SQRT(t/(one+one))
xupper = one / t

!   Start calculation

ind = 1
x = xvalue
IF (x < zero) THEN
  x = -x
  ind = -1
END IF

!   Code for X < =  1.0

IF (x <= one) THEN
  IF (x < xlow) THEN
    fn_val = x
  ELSE
    t = x * x
    t = (t-half) + (t-half)
    fn_val = x * cheval(nterms,atnina,t)
  END IF
ELSE
  
!   Code for X > 1.0
  
  IF (x > xupper) THEN
    fn_val = LOG(x) / twobpi
  ELSE
    t = one / (x*x)
    t = (t-half) + (t-half)
    fn_val = LOG(x) / twobpi + cheval(nterms,atnina,t) / x
  END IF
END IF
IF (ind < 0) fn_val = -fn_val
RETURN
END FUNCTION atnint



FUNCTION birint(xvalue) RESULT(fn_val)

!   DESCRIPTION:
!      This function calculates the integral of the Airy function Bi, defined

!          BIRINT(x) = integral{0 to x} Bi(t) dt

!      The program uses Chebyshev expansions, the coefficients of which
!      are given to 20 decimal places.


!   ERROR RETURNS:

!      If the function is too large and positive the correct
!      value would overflow. An error message is printed and the
!      program returns the value XMAX.

!      If the argument is too large and negative, it is impossible
!      to accurately compute the necessary SIN and COS functions,
!      for the asymptotic expansion.
!      An error message is printed, and the program returns the
!      value 0 (the value at -infinity).


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - INTEGER - The no. of terms to be used from the array
!                          ABINT1. The recommended value is such that
!                             ABS(ABINT1(NTERM1)) < EPS/100,
!                          subject to 1 <= NTERM1 <= 36.

!      NTERM2 - INTEGER - The no. of terms to be used from the array
!                          ABINT2. The recommended value is such that
!                             ABS(ABINT2(NTERM2)) < EPS/100,
!                          subject to 1 <= NTERM2 <= 37.

!      NTERM3 - INTEGER - The no. of terms to be used from the array
!                          ABINT3. The recommended value is such that
!                             ABS(ABINT3(NTERM3)) < EPS/100,
!                          subject to 1 <= NTERM3 <= 37.

!      NTERM4 - INTEGER - The no. of terms to be used from the array
!                          ABINT4. The recommended value is such that
!                             ABS(ABINT4(NTERM4)) < EPS/100,
!                          subject to 1 <= NTERM4 <= 20.

!      NTERM5 - INTEGER - The no. of terms to be used from the array
!                          ABINT5. The recommended value is such that
!                             ABS(ABINT5(NTERM5)) < EPS/100,
!                          subject to 1 <= NTERM5 <= 20.

!      XLOW1 - REAL (dp) - The value such that, if |x| < XLOW1,
!                          BIRINT(x) = x * Bi(0)
!                     to machine precision. The recommended value is
!                          2 * EPSNEG.

!      XHIGH1 - REAL (dp) - The value such that, if x > XHIGH1,
!                      the function value would overflow.
!                      The recommended value is computed as
!                          z = ln(XMAX) + 0.5ln(ln(XMAX)),
!                          XHIGH1 = (3z/2)^(2/3)

!      XNEG1 - REAL (dp) - The value such that, if x < XNEG1,
!                     the trigonometric functions in the asymptotic
!                     expansion cannot be calculated accurately.
!                     The recommended value is
!                          -(1/((EPS)**2/3))

!      XMAX - REAL (dp) - The value of the largest positive floating-pt
!                    number. Used in giving a value to the function
!                    if x > XHIGH1.

!      For values of EPS, EPSNEG, and XMAX see the file MACHCON.TXT.


!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.

!   INTRINSIC FUNCTIONS USED:
!                            COS, EXP, LOG, SIN, SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH

!   AUTHOR: Dr. Allan J. MacLeod,
!           Dept. of Mathematics and Statistics,
!           Univ. of Paisley,
!           High St.,
!           Paisley,
!           SCOTLAND.
!           PA1 2BE

!           (e-mail: macl_ms0@paisley.ac.uk )


!   LATEST REVISION:  23 January, 1996

REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3, nterm4, nterm5
REAL (dp)  :: arg, f1, f2, t, temp, x, xlow1, xhigh1, xmax, xneg1, z
CHARACTER (LEN=31) :: ermsg2 = 'ARGUMENT TOO LARGE AND NEGATIVE'
CHARACTER (LEN=31) :: ermsg1 = 'ARGUMENT TOO LARGE AND POSITIVE'
CHARACTER (LEN= 6) :: fnname = 'BIRINT'

REAL (dp), PARAMETER  :: abint1(0:36) = (/  &
    0.38683352445038543350_dp, -0.8823213550888908821D-1,  &
    0.21463937440355429239_dp, -0.4205347375891315126D-1,  &
    0.5932422547496086771D-1, -0.840787081124270210D-2,  &
    0.871824772778487955D-2, -0.12191600199613455D-3, 0.44024821786023234D-3, &
    0.27894686666386678D-3, -0.7052804689785537D-4, 0.5901080066770100D-4,  &
   -0.1370862587982142D-4, 0.505962573749073D-5, -0.51598837766735D-6,  &
    0.397511312349D-8, 0.9524985978055D-7, -0.3681435887321D-7,  &
    0.1248391688136D-7, -0.249097619137D-8, 0.31775245551D-9,  &
    0.5434365270D-10, -0.4024566915D-10, 0.1393855527D-10, -0.303817509D-11,  &
    0.40809511D-12, 0.1634116D-13, -0.2683809D-13, 0.896641D-14,  &
   -0.183089D-14, 0.21333D-15, 0.1108D-16, -0.1276D-16, 0.363D-17, -0.62D-18, &
    0.5D-19, 0.1D-19 /)

REAL (dp), PARAMETER  :: abint2(0:37) = (/  &
    2.04122078602516135181_dp, 0.2124133918621221230D-1,  &
    0.66617599766706276D-3, 0.3842047982808254D-4, 0.362310366020439D-5,  &
    0.50351990115074D-6, 0.7961648702253D-7, 0.717808442336D-8,  &
   -0.267770159104D-8, -0.168489514699D-8, -0.36811757255D-9, &
    0.4757128727D-10, 0.5263621945D-10, 0.778973500D-11, -0.460546143D-11,  &
   -0.183433736D-11, 0.32191249D-12, 0.29352060D-12, -0.1657935D-13,  &
   -0.4483808D-13, 0.27907D-15, 0.711921D-14, -0.1042D-16, -0.119591D-14,  &
    0.4606D-16, 0.20884D-15, -0.2416D-16, -0.3638D-16, 0.863D-17, 0.591D-17,  &
   -0.256D-17, -0.77D-18, 0.66D-18, 0.3D-19, -0.15D-18, 0.2D-19, 0.3D-19,  &
   -0.1D-19 /)

REAL (dp), PARAMETER  :: abint3(0:37) = (/  &
    0.31076961598640349251_dp, -0.27528845887452542718_dp,  &
    0.17355965706136543928_dp, -0.5544017909492843130D-1,  &
   -0.2251265478295950941D-1, 0.4107347447812521894D-1,  &
    0.984761275464262480D-2, -0.1555618141666041932D-1,  &
   -0.560871870730279234D-2, 0.246017783322230475D-2, 0.165740392292336978D-2, &
   -0.3277587501435402D-4, -0.24434680860514925D-3, -0.5035305196152321D-4,  &
    0.1630264722247854D-4, 0.851914057780934D-5, 0.29790363004664D-6,  &
   -0.64389707896401D-6, -0.15046988145803D-6, 0.1587013535823D-7,  &
    0.1276766299622D-7, 0.140578534199D-8, -0.46564739741D-9,  &
   -0.15682748791D-9, -0.403893560D-11, 0.666708192D-11, 0.128869380D-11,  &
   -0.6968663D-13, -0.6254319D-13, -0.718392D-14, 0.115296D-14, 0.42276D-15,  &
    0.2493D-16, -0.971D-17, -0.216D-17, -0.2D-19, 0.6D-19, 0.1D-19 /)

REAL (dp), PARAMETER  :: abint4(0:20) = (/  &
    1.99507959313352047614_dp, -0.273736375970692738D-2,  &
   -0.30897113081285850D-3, -0.3550101982798577D-4, -0.412179271520133D-5,  &
   -0.48235892316833D-6, -0.5678730727927D-7, -0.671874810365D-8,  &
   -0.79811649857D-9, -0.9514271478D-10, -0.1137468966D-10, -0.136359969D-11, &
   -0.16381418D-12, -0.1972575D-13, -0.237844D-14, -0.28752D-15, -0.3475D-16, &
   -0.422D-17, -0.51D-18, -0.6D-19, -0.1D-19 /)

REAL (dp), PARAMETER  :: abint5(0:20) = (/  &
    1.12672081961782566017_dp, -0.671405567525561198D-2,  &
   -0.69812918017832969D-3, -0.7561689886425276D-4, -0.834985574510207D-5,  &
   -0.93630298232480D-6, -0.10608556296250D-6, -0.1213128916741D-7,  &
   -0.139631129765D-8, -0.16178918054D-9, -0.1882307907D-10, -0.220272985D-11, &
   -0.25816189D-12, -0.3047964D-13, -0.358370D-14, -0.42831D-15, -0.4993D-16,  &
   -0.617D-17, -0.68D-18, -0.10D-18, -0.1D-19 /)

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, onept5 = 1.5_dp,  &
                         three = 3.0_dp, four = 4.0_dp, seven = 7.0_dp,  &
                         eight = 8.0_dp, nine = 9.0_dp, sixten = 16.0_dp,  &
                         onehun = 100.0_dp, ninhun = 900.0_dp,  &
                         thr644 = 3644.0_dp,  &
                         piby4 = 0.78539816339744830962_dp,  &
                         rt2b3p = 0.46065886596178063902_dp, &
                         birzer = 0.61492662744600073515_dp

!   Start computation

x = xvalue

!   Compute the machine-dependent constants.

t = EPSILON(zero)
f2 = one + one
xneg1 = -one / (t**(f2/three))
xmax = HUGE(zero)
f1 = LOG(xmax)
temp = f1 + LOG(f1) / f2
xhigh1 = (three*temp/f2) ** (f2/three)

!   Error test

IF (x > xhigh1) THEN
  CALL errprn(fnname,ermsg1)
  fn_val = xmax
  RETURN
END IF
IF (x < xneg1) THEN
  CALL errprn(fnname,ermsg2)
  fn_val = zero
  RETURN
END IF

!  continue with machine-dependent constants

xlow1 = f2 * t
t = t / onehun
IF (x >= zero) THEN
  DO  nterm1 = 36, 0, -1
    IF (ABS(abint1(nterm1)) > t) EXIT
  END DO
  DO  nterm2 = 37, 0, -1
    IF (ABS(abint2(nterm2)) > t) EXIT
  END DO
ELSE
  DO  nterm3 = 37, 0, -1
    IF (ABS(abint3(nterm3)) > t) EXIT
  END DO
  DO  nterm4 = 20, 0, -1
    IF (ABS(abint4(nterm4)) > t) EXIT
  END DO
  DO  nterm5 = 20, 0, -1
    IF (ABS(abint5(nterm5)) > t) EXIT
  END DO
END IF

!   Code for x >= 0.0

IF (x >= zero) THEN
  IF (x < xlow1) THEN
    fn_val = birzer * x
  ELSE
    IF (x <= eight) THEN
      t = x / four - one
      fn_val = x * EXP(onept5*x) * cheval(nterm1,abint1,t)
    ELSE
      t = sixten * SQRT(eight/x) / x - one
      z = (x+x) * SQRT(x) / three
      temp = rt2b3p * cheval(nterm2,abint2,t) / SQRT(z)
      temp = z + LOG(temp)
      fn_val = EXP(temp)
    END IF
  END IF
ELSE
  
!   Code for x < 0.0
  
  IF (x >= -seven) THEN
    IF (x > -xlow1) THEN
      fn_val = birzer * x
    ELSE
      t = -(x+x) / seven - one
      fn_val = x * cheval(nterm3,abint3,t)
    END IF
  ELSE
    z = -(x+x) * SQRT(-x) / three
    arg = z + piby4
    temp = nine * z * z
    t = (thr644-temp) / (ninhun+temp)
    f1 = cheval(nterm4,abint4,t) * SIN(arg)
    f2 = cheval(nterm5,abint5,t) * COS(arg) / z
    fn_val = (f2-f1) * rt2b3p / SQRT(z)
  END IF
END IF
RETURN
END FUNCTION birint



FUNCTION clausn(xvalue) RESULT(fn_val)

! DESCRIPTION:

!   This program calculates Clausen's integral defined by

!          CLAUSN(x) = integral 0 to x of (-ln(2*sin(t/2))) dt

!   The code uses Chebyshev expansions with the coefficients
!   given to 20 decimal places.


! ERROR RETURNS:

!   If |x| is too large it is impossible to reduce the argument
!   to the range [0,2*pi] with any precision. An error message
!   is printed and the program returns the value 0.0


! MACHINE-DEPENDENT CONSTANTS:

!   NTERMS - INTEGER - the no. of terms of the array ACLAUS
!                      to be used. The recommended value is
!                      such that ABS(ACLAUS(NTERMS)) < EPS/100
!                      subject to 1 <= NTERMS <= 15

!   XSMALL - REAL (dp) - the value below which Cl(x) can be
!                   approximated by x (1-ln x). The recommended
!                   value is pi*sqrt(EPSNEG/2).

!   XHIGH - REAL (dp) - The value of |x| above which we cannot
!                  reliably reduce the argument to [0,2*pi].
!                  The recommended value is   1/EPS.

!     For values of EPS and EPSNEG refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


! INTRINSIC FUNCTIONS USED:
!   AINT , LOG , SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


! AUTHOR:  Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley,
!          High St.
!          PAISLEY
!          SCOTLAND

!          ( e-mail: macl_ms0@paisley.ac.uk )


! LATEST MODIFICATION: 23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: indx, nterms
REAL (dp)  :: t, x, xhigh, xsmall
CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
CHARACTER (LEN= 6)  :: fnname = 'CLAUSN'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         onehun = 100.0_dp, pi = 3.1415926535897932385_dp,  &
                         pisq = 9.8696044010893586188_dp,  &
                         twopi = 6.2831853071795864769_dp,  &
                         twopia = 6.28125_dp, twopib = 0.19353071795864769253D-2
REAL (dp), PARAMETER  :: aclaus(0:15) = (/  &
    2.14269436376668844709_dp, 0.7233242812212579245D-1,  &
    0.101642475021151164D-2, 0.3245250328531645D-4,  &
    0.133315187571472D-5, 0.6213240591653D-7, 0.313004135337D-8,  &
    0.16635723056D-9, 0.919659293D-11, 0.52400462D-12,  &
    0.3058040D-13, 0.181969D-14, 0.11004D-15, 0.675D-17, 0.42D-18, 0.3D-19 /)

!  Start execution

x = xvalue

!   Compute the machine-dependent constants.

t = EPSILON(zero)
xhigh = one / t

!   Error test

IF (ABS(x) > xhigh) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Continue with machine-dependent constants

xsmall = pi * SQRT(half*t)
t = t / onehun
DO  nterms = 15, 0, -1
  IF (ABS(aclaus(nterms)) > t) GO TO 20
END DO

!  Continue with computation

20 indx = 1
IF (x < zero) THEN
  x = -x
  indx = -1
END IF

!  Argument reduced using simulated extra precision

IF (x > twopi) THEN
  t = AINT(x/twopi)
  x = (x-t*twopia) - t * twopib
END IF
IF (x > pi) THEN
  x = (twopia-x) + twopib
  indx = -indx
END IF

!  Set result to zero if X multiple of PI

IF (x == zero) THEN
  fn_val = zero
  RETURN
END IF

!  Code for X < XSMALL

IF (x < xsmall) THEN
  fn_val = x * (one-LOG(x))
ELSE
  
!  Code for XSMALL < =  X < =  PI
  
  t = (x*x) / pisq - half
  t = t + t
  IF (t > one) t = one
  fn_val = x * cheval(nterms,aclaus,t) - x * LOG(x)
END IF
IF (indx < 0) fn_val = -fn_val
RETURN
END FUNCTION clausn



FUNCTION debye1(xvalue) RESULT(fn_val)


!   DEFINITION:

!      This program calculates the Debye function of order 1, defined as

!            DEBYE1(x) = [Integral {0 to x} t/(exp(t)-1) dt] / x

!      The code uses Chebyshev series whose coefficients
!      are given to 20 decimal places.


!   ERROR RETURNS:

!      If XVALUE < 0.0 an error message is printed and the
!      function returns the value 0.0


!   MACHINE-DEPENDENT PARAMETERS:

!      NTERMS - INTEGER - The no. of elements of the array ADEB1.
!                         The recommended value is such that
!                             ABS(ADEB1(NTERMS)) < EPS/100 , with
!                                   1 <= NTERMS <= 18

!      XLOW - REAL (dp) - The value below which
!                    DEBYE1 = 1 - x/4 + x*x/36 to machine precision.
!                    The recommended value is
!                        SQRT(8*EPSNEG)

!      XUPPER - REAL (dp) - The value above which
!                      DEBYE1 = (pi*pi/(6*x)) - exp(-x)(x+1)/x.
!                      The recommended value is
!                          -LOG(2*EPS)

!      XLIM - REAL (dp) - The value above which DEBYE1 = pi*pi/(6*x)
!                    The recommended value is
!                          -LOG(XMIN)

!      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.

!   INTRINSIC FUNCTIONS USED:
!      AINT , EXP , INT , LOG , SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley
!          High St.
!          PAISLEY
!          SCOTLAND
!          PA1 2BE

!          (e-mail:  macl_ms0@paisley.ac.uk )


!   LATEST UPDATE:  23 january, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: i, nexp, nterms
REAL (dp)  :: expmx, rk, sum, t, x, xk, xlim, xlow, xupper
CHARACTER (LEN=17)  :: errmsg = 'ARGUMENT NEGATIVE'
CHARACTER (LEN= 6)  :: fnname = 'DEBYE1'

REAL (dp), PARAMETER  :: zero = 0.0_dp, quart = 0.25_dp, half = 0.5_dp,  &
                         one = 1.0_dp, four = 4.0_dp, eight = 8.0_dp,  &
                         nine = 9.0_dp, thirt6 = 36.0_dp, onehun = 100.0_dp, &
                         debinf = 0.60792710185402662866_dp

REAL (dp), PARAMETER  :: adeb1(0:18) = (/  &
    2.40065971903814101941_dp, 0.19372130421893600885_dp,  &
   -0.623291245548957703D-2, 0.35111747702064800D-3, -0.2282224667012310D-4, &
    0.158054678750300D-5, -0.11353781970719D-6, 0.835833611875D-8,  &
   -0.62644247872D-9, 0.4760334890D-10, -0.365741540D-11, 0.28354310D-12,  &
   -0.2214729D-13, 0.174092D-14, -0.13759D-15, 0.1093D-16, -0.87D-18,  &
    0.7D-19, -0.1D-19 /)

!   Start computation

x = xvalue

!   Check XVALUE >= 0.0

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

t = EPSILON(zero)
xlow = SQRT(t*eight)
xupper = -LOG(t+t)
xlim = -LOG(TINY(zero))
t = t / onehun
DO  nterms = 18, 0, -1
  IF (ABS(adeb1(nterms)) > t) GO TO 20
END DO

!   Code for x <= 4.0

20 IF (x <= four) THEN
  IF (x < xlow) THEN
    fn_val = ((x-nine)*x+thirt6) / thirt6
  ELSE
    t = ((x*x/eight)-half) - half
    fn_val = cheval(nterms,adeb1,t) - quart * x
  END IF
ELSE
  
!   Code for x > 4.0
  
  fn_val = one / (x*debinf)
  IF (x < xlim) THEN
    expmx = EXP(-x)
    IF (x > xupper) THEN
      fn_val = fn_val - expmx * (one+one/x)
    ELSE
      sum = zero
      rk = AINT(xlim/x)
      nexp = INT(rk)
      xk = rk * x
      DO  i = nexp, 1, -1
        t = (one+one/xk) / rk
        sum = sum * expmx + t
        rk = rk - one
        xk = xk - x
      END DO
      fn_val = fn_val - sum * expmx
    END IF
  END IF
END IF
RETURN
END FUNCTION debye1



FUNCTION debye2(xvalue) RESULT(fn_val)

!   DEFINITION:

!      This program calculates the Debye function of order 1, defined as

!            DEBYE2(x) = 2*[Integral {0 to x} t*t/(exp(t)-1) dt] / (x*x)

!      The code uses Chebyshev series whose coefficients
!      are given to 20 decimal places.


!   ERROR RETURNS:

!      If XVALUE < 0.0 an error message is printed and the
!      function returns the value 0.0


!   MACHINE-DEPENDENT PARAMETERS:

!      NTERMS - INTEGER - The no. of elements of the array ADEB2.
!                         The recommended value is such that
!                             ABS(ADEB2(NTERMS)) < EPS/100,
!                         subject to 1 <= NTERMS <= 18.

!      XLOW - REAL (dp) - The value below which
!                    DEBYE2 = 1 - x/3 + x*x/24 to machine precision.
!                    The recommended value is
!                        SQRT(8*EPSNEG)

!      XUPPER - REAL (dp) - The value above which
!                      DEBYE2 = (4*zeta(3)/x^2) - 2*exp(-x)(x^2+2x+1)/x^2.
!                      The recommended value is
!                          -LOG(2*EPS)

!      XLIM1 - REAL (dp) - The value above which DEBYE2 = 4*zeta(3)/x^2
!                     The recommended value is
!                          -LOG(XMIN)

!      XLIM2 - REAL (dp) - The value above which DEBYE2 = 0.0 to machine
!                     precision. The recommended value is
!                           SQRT(4.8/XMIN)

!      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT


!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      AINT , EXP , INT , LOG , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley
!          High St.
!          PAISLEY
!          SCOTLAND
!          PA1 2BE

!          (e-mail:  macl_ms0@paisley.ac.uk )


!   LATEST UPDATE:  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: i, nexp, nterms
REAL (dp)  :: expmx, rk, sum, t, x, xk, xlim1, xlim2, xlow, xupper
CHARACTER (LEN=17)  :: errmsg = 'ARGUMENT NEGATIVE'
CHARACTER (LEN= 6)  :: fnname = 'DEBYE2'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, three = 3.0_dp, four = 4.0_dp,  &
                         eight = 8.0_dp, twent4 = 24.0_dp, onehun = 100.0_dp, &
                         debinf = 4.80822761263837714160_dp
REAL (dp), PARAMETER  :: adeb2(0:18) = (/  &
    2.59438102325707702826_dp, 0.28633572045307198337_dp,  &
   -0.1020626561580467129D-1, 0.60491097753468435D-3, -0.4052576589502104D-4, &
    0.286338263288107D-5, -0.20863943030651D-6, 0.1552378758264D-7,  &
   -0.117312800866D-8, 0.8973585888D-10, -0.693176137D-11, 0.53980568D-12,  &
   -0.4232405D-13, 0.333778D-14, -0.26455D-15, 0.2106D-16, -0.168D-17,  &
    0.13D-18, -0.1D-19 /)

!   Start computation

x = xvalue

!   Check XVALUE >= 0.0

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

t = TINY(zero)
xlim1 = -LOG(t)
xlim2 = SQRT(debinf) / SQRT(t)
t = EPSILON(zero)
xlow = SQRT(t*eight)
xupper = -LOG(t+t)
t = t / onehun
DO  nterms = 18, 0, -1
  IF (ABS(adeb2(nterms)) > t) GO TO 20
END DO

!   Code for x <= 4.0

20 IF (x <= four) THEN
  IF (x < xlow) THEN
    fn_val = ((x-eight)*x+twent4) / twent4
  ELSE
    t = ((x*x/eight)-half) - half
    fn_val = cheval(nterms,adeb2,t) - x / three
  END IF
ELSE
  
!   Code for x > 4.0
  
  IF (x > xlim2) THEN
    fn_val = zero
  ELSE
    fn_val = debinf / (x*x)
    IF (x < xlim1) THEN
      expmx = EXP(-x)
      IF (x > xupper) THEN
        sum = ((x+two)*x+two) / (x*x)
      ELSE
        sum = zero
        rk = AINT(xlim1/x)
        nexp = INT(rk)
        xk = rk * x
        DO  i = nexp, 1, -1
          t = (one+two/xk+two/(xk*xk)) / rk
          sum = sum * expmx + t
          rk = rk - one
          xk = xk - x
        END DO
      END IF
      fn_val = fn_val - two * sum * expmx
    END IF
  END IF
END IF
RETURN
END FUNCTION debye2



FUNCTION debye3(xvalue) RESULT(fn_val)

!   DEFINITION:

!      This program calculates the Debye function of order 3, defined as

!            DEBYE3(x) = 3*[Integral {0 to x} t^3/(exp(t)-1) dt] / (x^3)

!      The code uses Chebyshev series whose coefficients
!      are given to 20 decimal places.


!   ERROR RETURNS:

!      If XVALUE < 0.0 an error message is printed and the
!      function returns the value 0.0


!   MACHINE-DEPENDENT PARAMETERS:

!      NTERMS - INTEGER - The no. of elements of the array ADEB3.
!                         The recommended value is such that
!                             ABS(ADEB3(NTERMS)) < EPS/100,
!                         subject to 1 <= NTERMS <= 18

!      XLOW - REAL (dp) - The value below which
!                    DEBYE3 = 1 - 3x/8 + x*x/20 to machine precision.
!                    The recommended value is
!                        SQRT(8*EPSNEG)

!      XUPPER - REAL (dp) - The value above which
!               DEBYE3 = (18*zeta(4)/x^3) - 3*exp(-x)(x^3+3x^2+6x+6)/x^3.
!                      The recommended value is
!                          -LOG(2*EPS)

!      XLIM1 - REAL (dp) - The value above which DEBYE3 = 18*zeta(4)/x^3
!                     The recommended value is
!                          -LOG(XMIN)

!      XLIM2 - REAL (dp) - The value above which DEBYE3 = 0.0 to machine
!                     precision. The recommended value is
!                          CUBE ROOT(19/XMIN)

!      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH

!   INTRINSIC FUNCTIONS USED:
!      AINT , EXP , INT , LOG , SQRT


!   AUTHOR:
!          Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley
!          High St.
!          PAISLEY
!          SCOTLAND
!          PA1 2BE

!          (e-mail:  macl_ms0@paisley.ac.uk )


!   LATEST UPDATE:  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: i, nexp, nterms
REAL (dp)  :: expmx, rk, sum, t, x, xk, xki, xlim1, xlim2, xlow, xupper
CHARACTER (LEN=17)  :: errmsg = 'ARGUMENT NEGATIVE'
CHARACTER (LEN= 6)  :: fnname = 'DEBYE3'

REAL (dp), PARAMETER  :: zero = 0.0_dp, pt375 = 0.375_dp, half = 0.5_dp,  &
                         one = 1.0_dp, three = 3.0_dp, four = 4.0_dp,  &
                         six = 6.0_dp, sevp5 = 7.5_dp, eight = 8.0_dp,  &
                         twenty = 20.0_dp, onehun = 100.0_dp,  &
                         debinf = 0.51329911273421675946D-1
REAL (dp), PARAMETER  :: adeb3(0:18) = (/  &
    2.70773706832744094526_dp, 0.34006813521109175100_dp,  &
   -0.1294515018444086863D-1, 0.79637553801738164D-3, -0.5463600095908238D-4, &
    0.392430195988049D-5, -0.28940328235386D-6, 0.2173176139625D-7,  &
   -0.165420999498D-8, 0.12727961892D-9, -0.987963459D-11, 0.77250740D-12,  &
   -0.6077972D-13, 0.480759D-14, -0.38204D-15, 0.3048D-16, -0.244D-17,  &
    0.20D-18, -0.2D-19 /)

!   Start computation

x = xvalue

!   Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

t = TINY(zero)
xlim1 = -LOG(t)
xk = one / three
xki = (one/debinf) ** xk
rk = t ** xk
xlim2 = xki / rk
t = EPSILON(zero)
xlow = SQRT(t*eight)
xupper = -LOG(t+t)
t = t / onehun
DO  nterms = 18, 0, -1
  IF (ABS(adeb3(nterms)) > t) GO TO 20
END DO

!   Code for x <= 4.0

20 IF (x <= four) THEN
  IF (x < xlow) THEN
    fn_val = ((x-sevp5)*x+twenty) / twenty
  ELSE
    t = ((x*x/eight)-half) - half
    fn_val = cheval(nterms,adeb3,t) - pt375 * x
  END IF
ELSE
  
!   Code for x > 4.0
  
  IF (x > xlim2) THEN
    fn_val = zero
  ELSE
    fn_val = one / (debinf*x*x*x)
    IF (x < xlim1) THEN
      expmx = EXP(-x)
      IF (x > xupper) THEN
        sum = (((x+three)*x+six)*x+six) / (x*x*x)
      ELSE
        sum = zero
        rk = AINT(xlim1/x)
        nexp = INT(rk)
        xk = rk * x
        DO  i = nexp, 1, -1
          xki = one / xk
          t = (((six*xki+six)*xki+three)*xki+one) / rk
          sum = sum * expmx + t
          rk = rk - one
          xk = xk - x
        END DO
      END IF
      fn_val = fn_val - three * sum * expmx
    END IF
  END IF
END IF
RETURN
END FUNCTION debye3



FUNCTION debye4(xvalue) RESULT(fn_val)

!   DEFINITION:

!      This program calculates the Debye function of order 4, defined as

!            DEBYE4(x) = 4*[Integral {0 to x} t^4/(exp(t)-1) dt] / (x^4)

!      The code uses Chebyshev series whose coefficients
!      are given to 20 decimal places.


!   ERROR RETURNS:

!      If XVALUE < 0.0 an error message is printed and the
!      function returns the value 0.0


!   MACHINE-DEPENDENT PARAMETERS:

!      NTERMS - INTEGER - The no. of elements of the array ADEB4.
!                         The recommended value is such that
!                             ABS(ADEB4(NTERMS)) < EPS/100,
!                         subject to 1 <= NTERMS <= 18

!      XLOW - REAL (dp) - The value below which
!                    DEBYE4 = 1 - 4x/10 + x*x/18 to machine precision.
!                    The recommended value is
!                        SQRT(8*EPSNEG)

!      XUPPER - REAL (dp) - The value above which
!               DEBYE4=(96*zeta(5)/x^4)-4*exp(-x)(x^4+4x^2+12x^2+24x+24)/x^4.
!                      The recommended value is
!                          -LOG(2*EPS)

!      XLIM1 - REAL (dp) - The value above which DEBYE4 = 96*zeta(5)/x^4
!                     The recommended value is
!                          -LOG(XMIN)

!      XLIM2 - REAL (dp) - The value above which DEBYE4 = 0.0 to machine
!                     precision. The recommended value is
!                          FOURTH ROOT(99/XMIN)

!      For values of EPS, EPSNEG, and XMIN see the file MACHCON.TXT


!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      AINT , EXP , INT , LOG , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley
!          High St.
!          PAISLEY
!          SCOTLAND
!          PA1 2BE

!          (e-mail:  macl_ms0@paisley.ac.uk )


!   LATEST UPDATE:  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: i, nexp, nterms
REAL (dp)  :: expmx, rk, sum, t, x, xk, xki, xlim1, xlim2, xlow, xupper
CHARACTER (LEN=17)  :: errmsg = 'ARGUMENT NEGATIVE'
CHARACTER (LEN= 6)  :: fnname = 'DEBYE4'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         twopt5 = 2.5_dp, four = 4.0_dp, five = 5.0_dp,  &
                         eight = 8.0_dp, twelve = 12.0_dp, eightn = 18.0_dp,  &
                         twent4 = 24.0_dp, forty5 = 45.0_dp, onehun = 100.0_dp, &
                         debinf = 99.54506449376351292781_dp
REAL (dp), PARAMETER  :: adeb4(0:18) = (/  &
    2.78186941502052346008_dp, 0.37497678352689286364_dp,  &
   -0.1494090739903158326D-1, 0.94567981143704274D-3,  &
   -0.6613291613893255D-4, 0.481563298214449D-5, -0.35880839587593D-6,  &
    0.2716011874160D-7, -0.208070991223D-8, 0.16093838692D-9,  &
   -0.1254709791D-10, 0.98472647D-12, -0.7772369D-13, 0.616483D-14,  &
   -0.49107D-15, 0.3927D-16, -0.315D-17, 0.25D-18, -0.2D-19 /)

!   Start computation

x = xvalue

!   Check XVALUE >= 0.0

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

t = TINY(zero)
xlim1 = -LOG(t)
rk = one / four
xk = debinf ** rk
xki = t ** rk
xlim2 = xk / xki
t = EPSILON(zero)
xlow = SQRT(t*eight)
xupper = -LOG(t+t)
t = t / onehun
DO  nterms = 18, 0, -1
  IF (ABS(adeb4(nterms)) > t) GO TO 20
END DO

!   Code for x <= 4.0

20 IF (x <= four) THEN
  IF (x < xlow) THEN
    fn_val = ((twopt5*x-eightn)*x+forty5) / forty5
  ELSE
    t = ((x*x/eight)-half) - half
    fn_val = cheval(nterms,adeb4,t) - (x+x) / five
  END IF
ELSE
  
!   Code for x > 4.0
  
  IF (x > xlim2) THEN
    fn_val = zero
  ELSE
    t = x * x
    fn_val = (debinf/t) / t
    IF (x < xlim1) THEN
      expmx = EXP(-x)
      IF (x > xupper) THEN
        sum = ((((x+four)*x+twelve)*x+twent4)*x+twent4) / (x*x*x*x )
      ELSE
        sum = zero
        rk = INT(xlim1/x)
        nexp = INT(rk)
        xk = rk * x
        DO  i = nexp, 1, -1
          xki = one / xk
          t = ((((twent4*xki+twent4)*xki+twelve)*xki+four)*xki+one ) / rk
          sum = sum * expmx + t
          rk = rk - one
          xk = xk - x
        END DO
      END IF
      fn_val = fn_val - four * sum * expmx
    END IF
  END IF
END IF
RETURN
END FUNCTION debye4



FUNCTION exp3(xvalue) RESULT(fn_val)

!   DESCRIPTION

!      This function calculates

!           EXP3(X) = integral 0 to X  (exp(-t*t*t)) dt

!      The code uses Chebyshev expansions, whose coefficients are
!      given to 20 decimal places.


!   ERROR RETURNS

!      If XVALUE < 0, an error message is printed and the function
!      returns the value 0.


!   MACHINE-DEPENDENT CONSTANTS

!      NTERM1 - INTEGER - The no. of terms of the array AEXP3,
!                         The recommended value is such that
!                               AEXP3(NTERM1) < EPS/100.

!      NTERM2 - INTEGER - The no. of terms of the array AEXP3A.
!                         The recommended value is such that
!                               AEXP3A(NTERM2) < EPS/100.

!      XLOW - REAL (dp) - The value below which EXP3(X) = X to machine
!                    precision. The recommended value is
!                          cube root(4*EPSNEG)

!      XUPPER - REAL (dp) - The value above which EXP3(X) = 0.89297...
!                      to machine precision. The recommended value is
!                           cube root(-ln(EPSNEG))

!      For values of EPS and EPSNEG for various machine/compiler
!      combinations refer to the file MACHCON.TXT.

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED
!      EXP, LOG

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR

!      DR. ALLAN J. MACLEOD,
!      DEPARTMENT OF MATHEMATICS AND STATISTICS,
!      UNIVERSITY OF PAISLEY,
!      HIGH ST.,
!      PAISLEY
!      SCOTLAND.

!      (e-mail  macl_ms0@paisley.ac.uk )


!   LATEST MODIFICATION:  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2
REAL (dp)  :: t, x, xlow, xupper
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'EXP3  '

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, three = 3.0_dp, four = 4.0_dp, &
                         sixten = 16.0_dp, onehun = 100.0_dp,  &
                         funinf = 0.89297951156924921122_dp
REAL (dp), PARAMETER  :: aexp3(0:24) = (/  &
    1.26919841422112601434_dp, -0.24884644638414098226_dp,  &
    0.8052622071723104125D-1, -0.2577273325196832934D-1,  &
    0.759987887307377429D-2, -0.203069558194040510D-2,  &
    0.49083458669932917D-3, -0.10768223914202077D-3, 0.2155172626428984D-4, &
   -0.395670513738429D-5, 0.66992409338956D-6, -0.10513218080703D-6, &
    0.1536258019825D-7, -0.209909603636D-8, 0.26921095381D-9,  &
   -0.3251952422D-10, 0.371148157D-11, -0.40136518D-12, 0.4123346D-13, &
   -0.403375D-14, 0.37658D-15, -0.3362D-16, 0.288D-17, -0.24D-18, 0.2D-19 /)

REAL (dp), PARAMETER  :: aexp3a(0:24) = (/  &
    1.92704649550682737293_dp, -0.3492935652048138054D-1, &
    0.145033837189830093D-2, -0.8925336718327903D-4, 0.705423921911838D-5, &
   -0.66717274547611D-6, 0.7242675899824D-7, -0.878258256056D-8, &
    0.116722344278D-8, -0.16766312812D-9, 0.2575501577D-10, -0.419578881D-11, &
    0.72010412D-12, -0.12949055D-12, 0.2428703D-13, -0.473311D-14,  &
    0.95531D-15, -0.19914D-15, 0.4277D-16, -0.944D-17, 0.214D-17, -0.50D-18, &
    0.12D-18, -0.3D-19, 0.1D-19 /)

!   Start calculation

x = xvalue

!   Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

t = EPSILON(zero)
xlow = (four*t) ** (one/three)
xupper = (-LOG(t)) ** (one/three)
t = t / onehun
IF (x <= two) THEN
  DO  nterm1 = 24, 0, -1
    IF (ABS(aexp3(nterm1)) > t) EXIT
  END DO
ELSE
  DO  nterm2 = 24, 0, -1
    IF (ABS(aexp3a(nterm2)) > t) EXIT
  END DO
END IF

!   Code for XVALUE < =  2

IF (x <= two) THEN
  IF (x < xlow) THEN
    fn_val = x
  ELSE
    t = ((x*x*x/four)-half) - half
    fn_val = x * cheval(nterm1,aexp3,t)
  END IF
ELSE
  
!   Code for XVALUE > 2
  
  IF (x > xupper) THEN
    fn_val = funinf
  ELSE
    t = ((sixten/(x*x*x))-half) - half
    t = cheval(nterm2,aexp3a,t)
    t = t * EXP(-x*x*x) / (three*x*x)
    fn_val = funinf - t
  END IF
END IF
RETURN
END FUNCTION exp3



FUNCTION goodst(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the function defined as

!        GOODST(x) = {integral 0 to inf} ( exp(-u*u)/(u+x) ) du

!      The code uses Chebyshev expansions whose coefficients are
!      given to 20 decimal places.


!   ERROR RETURNS:

!      If XVALUE <= 0.0, an error message is printed, and the
!      code returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - The no. of terms to be used in the array AGOST.
!                The recommended value is such that
!                    AGOST(NTERM1) < EPS/100,

!      NTERM2 - The no. of terms to be used in the array AGOSTA.
!                The recommended value is such that
!                    AGOSTA(NTERM2) < EPS/100,

!      XLOW - The value below which f(x) = -(gamma/2) - ln(x)
!             to machine precision. The recommended value is
!                EPSNEG

!      XHIGH - The value above which f(x) = sqrt(pi)/(2x) to
!              machine precision. The recommended value is
!                 2 / EPSNEG

!      For values of EPS and EPSNEG refer to the file MACHCON.TXT

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!       EXP , LOG

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:

!      Dr. Allan J. MacLeod,
!      Dept. of Mathematics and Statistics,
!      University of Paisley,
!      High St.,
!      Paisley.
!      SCOTLAND.

!      (e-mail: macl_ms0@paisley.ac.uk )


!   LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2
REAL (dp)  :: fval, t, x, xhigh, xlow
CHARACTER (LEN=15) :: errmsg = 'ARGUMENT <= 0.0'
CHARACTER (LEN= 6) :: fnname = 'GOODST'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp,  &
                         two = 2.0_dp, six = 6.0_dp, onehun = 100.0_dp,  &
                         gamby2 = 0.28860783245076643030_dp,  &
                         rtpib2 = 0.88622692545275801365_dp

REAL (dp), PARAMETER  :: agost(0:28) = (/  &
    0.63106560560398446247_dp, 0.25051737793216708827_dp,  &
   -0.28466205979018940757_dp, 0.8761587523948623552D-1,  &
    0.682602267221252724D-2, -0.1081129544192254677D-1,  &
    0.169101244117152176D-2, 0.50272984622615186D-3, -0.18576687204100084D-3, &
   -0.428703674168474D-5, 0.1009598903202905D-4, -0.86529913517382D-6,  &
   -0.34983874320734D-6, 0.6483278683494D-7, 0.757592498583D-8,  &
   -0.277935424362D-8, -0.4830235135D-10, 0.8663221283D-10, -0.394339687D-11, &
   -0.209529625D-11, 0.21501759D-12, 0.3959015D-13, -0.692279D-14,  &
   -0.54829D-15, 0.17108D-15, 0.376D-17, -0.349D-17, 0.7D-19, 0.6D-19 /)

REAL (dp), PARAMETER  :: agosta(0:23) = (/  &
    1.81775467984718758767_dp, -0.9921146570744097467D-1,  &
   -0.894058645254819243D-2, -0.94955331277726785D-3, -0.10971379966759665D-3, &
   -0.1346694539578590D-4, -0.172749274308265D-5, -0.22931380199498D-6,  &
   -0.3127844178918D-7, -0.436197973671D-8, -0.61958464743D-9,  &
   -0.8937991276D-10, -0.1306511094D-10, -0.193166876D-11, -0.28844270D-12,  &
   -0.4344796D-13, -0.659518D-14, -0.100801D-14, -0.15502D-15, -0.2397D-16,  &
   -0.373D-17, -0.58D-18, -0.9D-19, -0.1D-19 /)

!   Start computation

x = xvalue

!   Error test

IF (x <= zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

fval = EPSILON(zero)
t = fval / onehun
IF (x <= two) THEN
  DO  nterm1 = 28, 0, -1
    IF (ABS(agost(nterm1)) > t) EXIT
  END DO
  xlow = fval
ELSE
  DO  nterm2 = 23, 0, -1
    IF (ABS(agosta(nterm2)) > t) EXIT
  END DO
  xhigh = two / fval
END IF

!   Computation for 0 < x <= 2

IF (x <= two) THEN
  IF (x < xlow) THEN
    fn_val = -gamby2 - LOG(x)
  ELSE
    t = (x-half) - half
    fn_val = cheval(nterm1,agost,t) - EXP(-x*x) * LOG(x)
  END IF
ELSE
  
!   Computation for x > 2
  
  fval = rtpib2 / x
  IF (x > xhigh) THEN
    fn_val = fval
  ELSE
    t = (six-x) / (two+x)
    fn_val = fval * cheval(nterm2,agosta,t)
  END IF
END IF
RETURN
END FUNCTION goodst



FUNCTION i0int(xvalue) RESULT(fn_val)

!   DESCRIPTION:
!      This program computes the integral of the modified Bessel
!      function I0(x) using the definition

!         I0INT(x) = {integral 0 to x} I0(t) dt

!      The program uses Chebyshev expansions, the coefficients of
!      which are given to 20 decimal places.


!   ERROR RETURNS:
!      If |XVALUE| larger than a certain limit, the value of
!      I0INT would cause an overflow. If such a situation occurs
!      the programs prints an error message, and returns the
!      value sign(XVALUE)*XMAX, where XMAX is the largest
!      acceptable floating-pt. value.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - The no. of terms to be used from the array ARI01.
!                The recommended value is such that
!                    ABS(ARI01(NTERM1)) < EPS/100

!      NTERM2 - The no. of terms to be used from the array ARI0A.
!                The recommended value is such that
!                    ABS(ARI0A(NTERM2)) < EPS/100

!      XLOW - The value below which I0INT(x) = x, to machine precision.
!             The recommended value is
!                  sqrt(12*EPS).

!      XHIGH - The value above which overflow will occur. The
!              recommended value is
!                  ln(XMAX) + 0.5*ln(ln(XMAX)) + ln(2).

!      For values of EPS and XMAX refer to the file MACHCON.TXT.

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      EXP , LOG , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:

!      Dr. Allan J. MacLeod,
!      Dept. of Mathematics and Statistics,
!      University of Paisley,
!      High St.,
!      Paisley,
!      SCOTLAND
!      PA1 2BE

!      (e-mail :   macl_ms0@paisley.ac.uk )


!   LATEST REVISION:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: ind, nterm1, nterm2
REAL (dp)  :: t, temp, x, xhigh, xlow
CHARACTER (LEN=26) :: errmsg = 'SIZE OF ARGUMENT TOO LARGE'
CHARACTER (LEN= 6) :: fnname = 'I0INT '

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, three = 3.0_dp,  &
                         ateen = 18.0_dp, thirt6 = 36.0_dp, onehun = 100.0_dp, &
                         lnr2pi = 0.91893853320467274178_dp

REAL (dp), PARAMETER  :: ari01(0:28) = (/  &
    0.41227906926781516801_dp, -0.34336345150081519562_dp,  &
    0.22667588715751242585_dp, -0.12608164718742260032_dp,  &
    0.6012484628777990271D-1, -0.2480120462913358248D-1,  &
    0.892773389565563897D-2, -0.283253729936696605D-2, 0.79891339041712994D-3, &
   -0.20053933660964890D-3, 0.4416816783014313D-4, -0.822377042246068D-5,  &
    0.120059794219015D-5, -0.11350865004889D-6, 0.69606014466D-9,  &
    0.180622772836D-8, -0.26039481370D-9, -0.166188103D-11, 0.510500232D-11,  &
   -0.41515879D-12, -0.7368138D-13, 0.1279323D-13, 0.103247D-14, -0.30379D-15, &
   -0.1789D-16, 0.673D-17, 0.44D-18, -0.14D-18, -0.1D-19 /)

REAL (dp), PARAMETER  :: ari0a(0:33) = (/  &
    2.03739654571143287070_dp, 0.1917631647503310248D-1,  &
    0.49923334519288147D-3, 0.2263187103659815D-4, 0.158682108285561D-5,  &
    0.16507855636318D-6, 0.2385058373640D-7, 0.392985182304D-8,  &
    0.46042714199D-9, -0.7072558172D-10, -0.6747183961D-10, -0.2026962001D-10, &
   -0.87320338D-12, 0.175520014D-11, 0.60383944D-12, -0.3977983D-13,  &
   -0.8049048D-13, -0.1158955D-13, 0.827318D-14, 0.282290D-14, -0.77667D-15,  &
   -0.48731D-15, 0.7279D-16, 0.7873D-16, -0.785D-17, -0.1281D-16, 0.121D-17,  &
    0.214D-17, -0.27D-18, -0.36D-18, 0.7D-19, 0.6D-19, -0.2D-19, -0.1D-19 /)

!   Start computation

ind = 1
x = xvalue
IF (xvalue < zero) THEN
  ind = -1
  x = -x
END IF

!   Compute the machine-dependent constants.

t = LOG(HUGE(zero))
xhigh = t + LOG(t) * half - LOG(half)

!   Error test

IF (x > xhigh) THEN
  CALL errprn(fnname,errmsg)
  fn_val = EXP(xhigh-lnr2pi-half*LOG(xhigh))
  IF (ind == -1) fn_val = -fn_val
  RETURN
END IF

!   Continue with machine-constants

temp = EPSILON(zero)
t = temp / onehun
IF (x <= ateen) THEN
  DO  nterm1 = 28, 0, -1
    IF (ABS(ari01(nterm1)) > t) EXIT
  END DO
  xlow = SQRT(thirt6*temp/three)
ELSE
  DO  nterm2 = 33, 0, -1
    IF (ABS(ari0a(nterm2)) > t) EXIT
  END DO
END IF

!   Code for 0 <= |x| <= 18

IF (x <= ateen) THEN
  IF (x < xlow) THEN
    fn_val = x
  ELSE
    t = (three*x-ateen) / (x+ateen)
    fn_val = x * EXP(x) * cheval(nterm1,ari01,t)
  END IF
ELSE
  
!   Code for |x| > 18
  
  t = (thirt6/x-half) - half
  temp = x - half * LOG(x) - lnr2pi + LOG(cheval(nterm2,ari0a,t))
  fn_val = EXP(temp)
END IF
IF (ind == -1) fn_val = -fn_val
RETURN
END FUNCTION i0int



FUNCTION i0ml0(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This program calculates the function I0ML0 defined as

!                I0ML0(x) = I0(x) - L0(x)

!      where I0(x) is the modified Bessel function of the first kind of
!      order 0, and L0(x) is the modified Struve function of order 0.

!      The code uses Chebyshev expansions with the coefficients
!      given to an accuracy of 20D.


!   ERROR RETURNS:

!      The coefficients are only suitable for XVALUE >= 0.0. If
!      XVALUE < 0.0, an error message is printed and the function
!      returns the value 0.0


!   MACHINE-DEPENDENT PARAMETERS:

!      NTERM1 - INTEGER - The number of terms required for the array
!                         AI0L0. The recommended value is such that
!                              ABS(AI0L0(NTERM1)) < EPS/100

!      NTERM2 - INTEGER - The number of terms required for the array
!                         AI0L0A. The recommended value is such that
!                              ABS(AI0L0A(NTERM2)) < EPS/100

!      XLOW - REAL (dp) - The value below which I0ML0(x) = 1 to machine
!                    precision. The recommended value is
!                               EPSNEG

!      XHIGH - REAL (dp) - The value above which I0ML0(x) = 2/(pi*x) to
!                     machine precision. The recommended value is
!                               SQRT(800/EPS)

!      For values of EPS, and EPSNEG see the file MACHCON.TXT

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. MacLeod
!          Dept. of Mathematics and Statistics
!          University of Paisley
!          High St.
!          Paisley
!          SCOTLAND
!          PA1 2BE

!          ( e-mail: macl_ms0@paisley.ac.uk )


!   LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2
REAL (dp)  :: t, x, xhigh, xlow, xsq
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'I0ML0 '

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, six = 6.0_dp,  &
                         sixten = 16.0_dp, forty = 40.0_dp, two88 = 288.0_dp, &
                         onehun = 100.0_dp, atehun = 800.0_dp,  &
                         twobpi = 0.63661977236758134308_dp

REAL (dp), PARAMETER  :: ai0l0(0:23) = (/  &
    0.52468736791485599138_dp, -0.35612460699650586196_dp,  &
    0.20487202864009927687_dp, -0.10418640520402693629_dp,  &
    0.4634211095548429228D-1, -0.1790587192403498630D-1,  &
    0.597968695481143177D-2, -0.171777547693565429D-2,  &
    0.42204654469171422D-3, -0.8796178522094125D-4, 0.1535434234869223D-4,  &
   -0.219780769584743D-5, 0.24820683936666D-6, -0.2032706035607D-7,  &
    0.90984198421D-9, 0.2561793929D-10, -0.710609790D-11, 0.32716960D-12,  &
    0.2300215D-13, -0.292109D-14, -0.3566D-16, 0.1832D-16, -0.10D-18,  &
   -0.11D-18 /)

REAL (dp), PARAMETER  :: ai0l0a(0:23) = (/  &
    2.00326510241160643125_dp, 0.195206851576492081D-2,  &
    0.38239523569908328D-3, 0.7534280817054436D-4,  &
    0.1495957655897078D-4, 0.299940531210557D-5, 0.60769604822459D-6,  &
    0.12399495544506D-6, 0.2523262552649D-7, 0.504634857332D-8,  &
    0.97913236230D-9, 0.18389115241D-9, 0.3376309278D-10, 0.611179703D-11,  &
    0.108472972D-11, 0.18861271D-12, 0.3280345D-13, 0.565647D-14,  &
    0.93300D-15, 0.15881D-15, 0.2791D-16, 0.389D-17, 0.70D-18, 0.16D-18 /)

!   Start computation

x = xvalue

!   Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xsq = EPSILON(zero)
t = xsq / onehun
IF (x <= sixten) THEN
  DO  nterm1 = 23, 0, -1
    IF (ABS(ai0l0(nterm1)) > t) EXIT
  END DO
  xlow = xsq
ELSE
  DO  nterm2 = 23, 0, -1
    IF (ABS(ai0l0a(nterm2)) > t) EXIT
  END DO
  xhigh = SQRT(atehun/xsq)
END IF

!   Code for x <= 16

IF (x <= sixten) THEN
  IF (x < xlow) THEN
    fn_val = one
    RETURN
  ELSE
    t = (six*x-forty) / (x+forty)
    fn_val = cheval(nterm1,ai0l0,t)
    RETURN
  END IF
ELSE
  
!   Code for x > 16
  
  IF (x > xhigh) THEN
    fn_val = twobpi / x
  ELSE
    xsq = x * x
    t = (atehun-xsq) / (two88+xsq)
    fn_val = cheval(nterm2,ai0l0a,t) * twobpi / x
  END IF
END IF
RETURN
END FUNCTION i0ml0



FUNCTION i1ml1(xvalue) RESULT(fn_val)

! DESCRIPTION:

!    This program calculates the function I1ML1 defined as

!              I1ML1(x) = I1(x) - L1(x)

!    where I1(x) is the modified Bessel function of the first kind of
!    order 1, and L1(x) is the modified Struve function of order 1.

!    The code uses Chebyshev expansions with the coefficients
!    given to an accuracy of 20D.


! ERROR RETURNS:

!    The coefficients are only suitable for XVALUE >= 0.0. If
!    XVALUE < 0.0, an error message is printed and the function
!    returns the value 0.0


! MACHINE-DEPENDENT PARAMETERS:

!    NTERM1 - INTEGER - The number of terms required for the array
!                       AI1L1. The recommended value is such that
!                            ABS(AI1L1(NTERM1)) < EPS/100

!    NTERM2 - INTEGER - The number of terms required for the array
!                       AI1L1A. The recommended value is such that
!                            ABS(AI1L1A(NTERM2)) < EPS/100

!    XLOW - REAL (dp) - The value below which I1ML1(x) = x/2 to machine
!                  precision. The recommended value is
!                             2*EPSNEG

!    XHIGH - REAL (dp) - The value above which I1ML1(x) = 2/pi to
!                   machine precision. The recommended value is
!                             SQRT(800/EPS)

!    For values of EPS, and EPSNEG see the file MACHCON.TXT

!    The machine-dependent constants are computed internally by
!    using the D1MACH subroutine.


! INTRINSIC FUNCTIONS USED:
!    ABS , SQRT

! OTHER MISCFUN SUBROUTINES USED:
!        CHEVAL , ERRPRN, D1MACH


! AUTHOR:
!        Dr. Allan J. MacLeod
!        Dept. of Mathematics and Statistics
!        University of Paisley
!        High St.
!        Paisley
!        SCOTLAND
!        PA1 2BE

!        (e-mail: macl_ms0@paisley.ac.uk )


! LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2
REAL (dp)  :: t, x, xhigh, xlow, xsq
CHARACTER (LEN=14) :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6) :: fnname = 'I1ML1 '

REAL (dp), PARAMETER  :: zero = 0.0_dp, two = 2.0_dp, six = 6.0_dp,  &
                         sixten = 16.0_dp, forty = 40.0_dp, onehun = 100.0_dp, &
                         two88 = 288.0_dp, atehun = 800.0_dp, &
                         twobpi = 0.63661977236758134308_dp

REAL (dp), PARAMETER  :: ai1l1(0:23) = (/  &
    0.67536369062350576137_dp, -0.38134971097266559040_dp,  &
    0.17452170775133943559_dp, -0.7062105887235025061D-1,  &
    0.2517341413558803702D-1, -0.787098561606423321D-2,  &
    0.214814368651922006D-2, -0.50862199717906236D-3, 0.10362608280442330D-3, &
   -0.1795447212057247D-4, 0.259788274515414D-5, -0.30442406324667D-6,  &
    0.2720239894766D-7, -0.158126144190D-8, 0.1816209172D-10,  &
    0.647967659D-11, -0.54113290D-12, -0.308311D-14, 0.305638D-14,  &
   -0.9717D-16, -0.1422D-16, 0.84D-18, 0.7D-19, -0.1D-19 /)

REAL (dp), PARAMETER  :: ai1l1a(0:25) = (/  &
    1.99679361896789136501_dp, -0.190663261409686132D-2,  &
   -0.36094622410174481D-3, -0.6841847304599820D-4, -0.1299008228509426D-4,  &
   -0.247152188705765D-5, -0.47147839691972D-6, -0.9020819982592D-7,  &
   -0.1730458637504D-7, -0.332323670159D-8, -0.63736421735D-9,  &
   -0.12180239756D-9, -0.2317346832D-10, -0.439068833D-11, -0.82847110D-12,  &
   -0.15562249D-12, -0.2913112D-13, -0.543965D-14, -0.101177D-14,  &
   -0.18767D-15, -0.3484D-16, -0.643D-17, -0.118D-17, -0.22D-18, -0.4D-19,  &
   -0.1D-19 /)

!   Start computation

x = xvalue

!   Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xsq = EPSILON(zero)
t = xsq / onehun
IF (x <= sixten) THEN
  DO  nterm1 = 23, 0, -1
    IF (ABS(ai1l1(nterm1)) > t) EXIT
  END DO
  xlow = xsq + xsq
ELSE
  DO  nterm2 = 25, 0, -1
    IF (ABS(ai1l1a(nterm2)) > t) EXIT
  END DO
  xhigh = SQRT(atehun/xsq)
END IF

!   Code for x <= 16

IF (x <= sixten) THEN
  IF (x < xlow) THEN
    fn_val = x / two
    RETURN
  ELSE
    t = (six*x-forty) / (x+forty)
    fn_val = cheval(nterm1,ai1l1,t) * x / two
    RETURN
  END IF
ELSE
  
!   Code for x > 16
  
  IF (x > xhigh) THEN
    fn_val = twobpi
  ELSE
    xsq = x * x
    t = (atehun-xsq) / (two88+xsq)
    fn_val = cheval(nterm2,ai1l1a,t) * twobpi
  END IF
END IF
RETURN
END FUNCTION i1ml1



FUNCTION j0int(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the integral of the Bessel
!      function J0, defined as

!        J0INT(x) = {integral 0 to x} J0(t) dt

!      The code uses Chebyshev expansions whose coefficients are
!      given to 20 decimal places.


!   ERROR RETURNS:

!      If the value of |x| is too large, it is impossible to
!      accurately compute the trigonometric functions used. An
!      error message is printed, and the function returns the value 1.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - The no. of terms to be used from the array
!                ARJ01. The recommended value is such that
!                   ABS(ARJ01(NTERM1)) < EPS/100, provided that

!      NTERM2 - The no. of terms to be used from the array
!                ARJ0A1. The recommended value is such that
!                   ABS(ARJ0A1(NTERM2)) < EPS/100, provided that

!      NTERM3 - The no. of terms to be used from the array
!                ARJ0A2. The recommended value is such that
!                   ABS(ARJ0A2(NTERM3)) < EPS/100, provided that

!      XLOW - The value of |x| below which J0INT(x) = x to
!             machine-precision. The recommended value is
!                 sqrt(12*EPSNEG)

!      XHIGH - The value of |x| above which it is impossible
!              to calculate (x-pi/4) accurately. The recommended
!              value is      1/EPSNEG

!      For values of EPS and EPSNEG for various machine/compiler
!      combinations refer to the file MACHCON.TXT.

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      COS , SIN , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley,
!          Paisley,
!          SCOTLAND
!          PA1 2BE

!          (e-mail:   macl_ms0@paisley.ac.uk )


!   LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: ind, nterm1, nterm2, nterm3
REAL (dp)  :: pib41, t, temp, x, xhigh, xlow, xmpi4
CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
CHARACTER (LEN= 6)  :: fnname = 'J0INT '

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, twelve = 12.0_dp,  &
                         sixten = 16.0_dp, onehun = 100.0_dp, one28 = 128.0_dp, &
                         five12 = 512.0_dp,  &
                         rt2bpi = 0.79788456080286535588_dp,  &
                         pib411 = 201.0_dp, pib412 = 256.0_dp,  &
                         pib42 = 0.24191339744830961566D-3

REAL (dp), PARAMETER  :: arj01(0:23) = (/  &
    0.38179279321690173518_dp, -0.21275636350505321870_dp,  &
    0.16754213407215794187_dp, -0.12853209772196398954_dp,  &
    0.10114405455778847013_dp, -0.9100795343201568859D-1,  &
    0.6401345264656873103D-1, -0.3066963029926754312D-1,  &
    0.1030836525325064201D-1, -0.255670650399956918D-2,  &
    0.48832755805798304D-3, -0.7424935126036077D-4, 0.922260563730861D-5,  &
   -0.95522828307083D-6, 0.8388355845986D-7, -0.633184488858D-8,  &
    0.41560504221D-9, -0.2395529307D-10, 0.122286885D-11, -0.5569711D-13,  &
    0.227820D-14, -0.8417D-16, 0.282D-17, -0.9D-19 /)

REAL (dp), PARAMETER  :: arj0a1(0:21) = (/  &
    1.24030133037518970827_dp, -0.478125353632280693D-2,  &
    0.6613148891706678D-4, -0.186042740486349D-5, 0.8362735565080D-7,  &
   -0.525857036731D-8, 0.42606363251D-9, -0.4211761024D-10, 0.488946426D-11,  &
   -0.64834929D-12, 0.9617234D-13, -0.1570367D-13, 0.278712D-14, -0.53222D-15, &
    0.10844D-15, -0.2342D-16, 0.533D-17, -0.127D-17, 0.32D-18, -0.8D-19,  &
    0.2D-19, -0.1D-19 /)

REAL (dp), PARAMETER  :: arj0a2(0:18) = (/  &
    1.99616096301341675339_dp, -0.190379819246668161D-2,  &
    0.1539710927044226D-4, -0.31145088328103D-6, 0.1110850971321D-7,  &
   -0.58666787123D-9, 0.4139926949D-10, -0.365398763D-11, 0.38557568D-12,  &
   -0.4709800D-13, 0.650220D-14, -0.99624D-15, 0.16700D-15, -0.3028D-16,  &
    0.589D-17, -0.122D-17, 0.27D-18, -0.6D-19, 0.1D-19 /)

!   Start computation

x = xvalue
ind = 1
IF (x < zero) THEN
  x = -x
  ind = -1
END IF

!   Compute the machine-dependent constants.

temp = EPSILON(zero)
xhigh = one / temp

!   Error test

IF (x > xhigh) THEN
  CALL errprn(fnname,errmsg)
  fn_val = one
  IF (ind == -1) fn_val = -fn_val
  RETURN
END IF

!   continue with constants

t = temp / onehun
IF (x <= sixten) THEN
  DO  nterm1 = 23, 0, -1
    IF (ABS(arj01(nterm1)) > t) EXIT
  END DO
  xlow = SQRT(twelve*temp)
ELSE
  DO  nterm2 = 21, 0, -1
    IF (ABS(arj0a1(nterm2)) > t) EXIT
  END DO
  DO  nterm3 = 18, 0, -1
    IF (ABS(arj0a2(nterm3)) > t) EXIT
  END DO
END IF

!   Code for 0 <= |x| <= 16

IF (x <= sixten) THEN
  IF (x < xlow) THEN
    fn_val = x
  ELSE
    t = x * x / one28 - one
    fn_val = x * cheval(nterm1,arj01,t)
  END IF
ELSE
  
!   Code for |x| > 16
  
  t = five12 / (x*x) - one
  pib41 = pib411 / pib412
  xmpi4 = (x-pib41) - pib42
  temp = COS(xmpi4) * cheval(nterm2,arj0a1,t) / x
  temp = temp - SIN(xmpi4) * cheval(nterm3,arj0a2,t)
  fn_val = one - rt2bpi * temp / SQRT(x)
END IF
IF (ind == -1) fn_val = -fn_val
RETURN
END FUNCTION j0int



FUNCTION k0int(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the integral of the modified Bessel function
!      defined by

!         K0INT(x) = {integral 0 to x} K0(t) dt

!      The code uses Chebyshev expansions, whose coefficients are
!      given to 20 decimal places.


!   ERROR RETURNS:

!      If XVALUE < 0.0, the function is undefined. An error message is
!      printed and the function returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - The no. of terms to be used in the array AK0IN1. The
!                recommended value is such that
!                   ABS(AK0IN1(NTERM1)) < EPS/100,

!      NTERM2 - The no. of terms to be used in the array AK0IN2. The
!                recommended value is such that
!                   ABS(AK0IN2(NTERM2)) < EPS/100,

!      NTERM3 - The no. of terms to be used in the array AK0INA. The
!                recommended value is such that
!                   ABS(AK0INA(NTERM3)) < EPS/100,

!      XLOW - The value below which K0INT = x * ( const - ln(x) ) to
!             machine precision. The recommended value is
!                   sqrt (18*EPSNEG).

!      XHIGH - The value above which K0INT = pi/2 to machine precision.
!              The recommended value is
!                   - log (2*EPSNEG)

!      For values of EPS and EPSNEG refer to the file MACHCON.TXT.

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      EXP , LOG , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!         Dr. Allan J. MacLeod,
!         Dept. of Mathematics and Statistics,
!         University of Paisley,
!         High St.,
!         Paisley,
!         SCOTLAND

!         (e-mail: macl_ms0@paisley.ac.uk )


!   LATEST REVISION:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3
REAL (dp)  :: fval, t, temp, x, xhigh, xlow
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 8)  :: fnname = 'K0INT '

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, six = 6.0_dp,  &
                         twelve = 12.0_dp, eightn = 18.0_dp, onehun = 100.0_dp, &
                         const1 = 1.11593151565841244881_dp,   &
                         const2 =  -0.11593151565841244881_dp, &
                         piby2  = 1.57079632679489661923_dp,   &
                         rt2bpi = 0.79788456080286535588_dp

REAL (dp), PARAMETER  :: ak0in1(0:15) = (/  &
    16.79702714464710959477D0, 9.79134687676889407070D0,  &
     2.80501316044337939300D0, 0.45615620531888502068D0,  &
     0.4716224457074760784D-1, 0.335265148269698289D-2,  &
     0.17335181193874727D-3, 0.679951889364702D-5, 0.20900268359924D-6,  &
     0.516603846976D-8, 0.10485708331D-9, 0.177829320D-11, 0.2556844D-13,  &
     0.31557D-15, 0.338D-17, 0.3D-19 /)
REAL (dp), PARAMETER  :: ak0in2(0:15) = (/  &
    10.76266558227809174077D0, 5.62333479849997511550D0,  &
     1.43543664879290867158D0, 0.21250410143743896043D0,  &
     0.2036537393100009554D-1, 0.136023584095623632D-2,  &
     0.6675388699209093D-4, 0.250430035707337D-5, 0.7406423741728D-7,  &
     0.176974704314D-8, 0.3485775254D-10, 0.57544785D-12, 0.807481D-14,  &
     0.9747D-16, 0.102D-17, 0.1D-19 /)
REAL (dp), PARAMETER  :: ak0ina(0:27) = (/  &
     1.91172065445060453895D0, -0.4183064565769581085D-1,  &
     0.213352508068147486D-2, -0.15859497284504181D-3, 0.1497624699858351D-4, &
    -0.167955955322241D-5, 0.21495472478804D-6, -0.3058356654790D-7, &
     0.474946413343D-8, -0.79424660432D-9, 0.14156555325D-9,  &
    -0.2667825359D-10, 0.528149717D-11, -0.109263199D-11, 0.23518838D-12,  &
    -0.5247991D-13, 0.1210191D-13, -0.287632D-14, 0.70297D-15, -0.17631D-15, &
     0.4530D-16, -0.1190D-16, 0.319D-17, -0.87D-18, 0.24D-18, -0.7D-19,  &
     0.2D-19, -0.1D-19 /)

!   Start computation

x = xvalue

!   Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

temp = EPSILON(zero)
t = temp / onehun
IF (x <= six) THEN
  DO  nterm1 = 15, 0, -1
    IF (ABS(ak0in1(nterm1)) > t) EXIT
  END DO
  DO  nterm2 = 15, 0, -1
    IF (ABS(ak0in2(nterm2)) > t) EXIT
  END DO
  xlow = SQRT(eightn*temp)
ELSE
  DO  nterm3 = 27, 0, -1
    IF (ABS(ak0ina(nterm3)) > t) EXIT
  END DO
  xhigh = -LOG(temp+temp)
END IF

!   Code for 0 <= XVALUE <= 6

IF (x <= six) THEN
  IF (x < xlow) THEN
    fval = x
    IF (x > zero) THEN
      fval = fval * (const1-LOG(x))
    END IF
    fn_val = fval
  ELSE
    t = ((x*x)/eightn-half) - half
    fval = (const2+LOG(x)) * cheval(nterm2,ak0in2,t)
    fn_val = x * (cheval(nterm1,ak0in1,t)-fval)
  END IF
  
!   Code for x > 6
  
ELSE
  fval = piby2
  IF (x < xhigh) THEN
    t = (twelve/x-half) - half
    temp = EXP(-x) * cheval(nterm3,ak0ina,t)
    fval = fval - temp / (SQRT(x)*rt2bpi)
  END IF
  fn_val = fval
END IF
RETURN
END FUNCTION k0int



FUNCTION lobach(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the Lobachewsky function L(x), defined as

!         LOBACH(x) = {integral 0 to x} ( -ln ( | cos t | ) dt

!      The code uses Chebyshev expansions whose coefficients are given
!      to 20 decimal places.


!   ERROR RETURNS:

!      If |x| too large, it is impossible to accurately reduce the
!      argument to the range [0,pi]. An error message is printed
!      and the program returns the value 0.0


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - INTEGER - The no. of terms to be used of the array ARLOB1.
!                          The recommended value is such that
!                          ABS(ARLOB1(NTERM1)) < EPS/100

!      NTERM2 - INTEGER - The no. of terms to be used of the array ARLOB2.
!                          The recommended value is such that
!                          ABS(ARLOB2(NTERM2)) < EPS/100

!      XLOW1 - DOUBLE PRECISION - The value below which L(x) = 0.0 to
!                          machine-precision.
!                     The recommended value is
!                              cube-root ( 6*XMIN )

!      XLOW2 - REAL (dp) - The value below which L(x) = x**3/6 to
!                     machine-precision. The recommended value is
!                              sqrt ( 10*EPS )

!      XLOW3 - REAL (dp) - The value below which
!                         L(pi/2) - L(pi/2-x) = x ( 1 - log(x) )
!                     to machine-precision. The recommended value is
!                               sqrt ( 18*EPS )

!      XHIGH - REAL (dp) - The value of |x| above which it is impossible
!                     to accurately reduce the argument. The
!                     recommended value is   1 / EPS.

!      For values of EPS, and XMIN, refer to the file MACHCON.TXT

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      INT , LOG , SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:

!      Dr. Allan J. MacLeod,
!      Dept. of Mathematics and Statistics,
!      University of Paisley,
!      High St.,
!      Paisley,
!      SCOTLAND

!      ( e-mail: macl_ms0@paisley.ac.uk )


!   LATEST UPDATE:    23 January, 1996

REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: indpi2, indsgn, npi, nterm1, nterm2
REAL (dp)  :: fval, fval1, lbpb21, lobpi1, pi, piby2, piby21, piby4,  &
              pi1, t, x, xcub, xhigh, xlow1, xlow2, xlow3, xr
CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
CHARACTER (LEN= 6)  :: fnname = 'LOBACH'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, six = 6.0_dp, ten = 10.0_dp,  &
                         onehun = 100.0_dp, lobpia = 1115.0_dp,  &
                         lobpib = 512.0_dp,  &
                         lobpi2 = -1.48284696397869499311D-4, &
                         lbpb22 = -7.41423481989347496556D-5, &
                         pi11 = 201.0_dp, pi12 = 64.0_dp,  &
                         pi2 = 9.67653589793238462643D-4,  &
                         piby22 = 4.83826794896619231322D-4,  &
                         tcon = 3.24227787655480868620_dp

REAL (dp), PARAMETER  :: arlob1(0:15) = (/  &
    0.34464884953481300507_dp, 0.584198357190277669D-2,  &
    0.19175029694600330D-3, 0.787251606456769D-5,  &
    0.36507477415804D-6, 0.1830287272680D-7, 0.96890333005D-9,  &
    0.5339055444D-10, 0.303408025D-11, 0.17667875D-12,  &
    0.1049393D-13, 0.63359D-15, 0.3878D-16, 0.240D-17, 0.15D-18, 0.1D-19 /)
REAL (dp), PARAMETER  :: arlob2(0:10) = (/  &
    2.03459418036132851087_dp, 0.1735185882027407681D-1,  &
    0.5516280426090521D-4, 0.39781646276598D-6, 0.369018028918D-8,  &
    0.3880409214D-10, 0.44069698D-12, 0.527674D-14, 0.6568D-16, 0.84D-18,  &
    0.1D-19 /)

!   Start computation

x = ABS(xvalue)
indsgn = 1
IF (xvalue < zero) THEN
  indsgn = -1
END IF

!   Compute the machine-dependent constants.

xr = EPSILON(zero)
xhigh = one / xr

!   Error test

IF (x > xhigh) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   continue with constants

t = xr / onehun
DO  nterm1 = 15, 0, -1
  IF (ABS(arlob1(nterm1)) > t) EXIT
END DO
DO  nterm2 = 10, 0, -1
  IF (ABS(arlob2(nterm2)) > t) EXIT
END DO
xlow1 = (six*TINY(zero)) ** (two/six)
xlow2 = SQRT(ten*xr)
t = two * ten - two
xlow3 = SQRT(t*xr)

!   Reduce argument to [0,pi]

pi1 = pi11 / pi12
pi = pi1 + pi2
piby2 = pi / two
piby21 = pi1 / two
piby4 = piby2 / two
npi = INT(x/pi)
xr = (x-npi*pi1) - npi * pi2

!   Reduce argument to [0,pi/2]

indpi2 = 0
IF (xr > piby2) THEN
  indpi2 = 1
  xr = (pi1-xr) + pi2
END IF

!   Code for argument in [0,pi/4]

IF (xr <= piby4) THEN
  IF (xr < xlow1) THEN
    fval = zero
  ELSE
    xcub = xr * xr * xr
    IF (xr < xlow2) THEN
      fval = xcub / six
    ELSE
      t = (tcon*xr*xr-half) - half
      fval = xcub * cheval(nterm1,arlob1,t)
    END IF
  END IF
ELSE
  
!   Code for argument in [pi/4,pi/2]
  
  xr = (piby21-xr) + piby22
  IF (xr == zero) THEN
    fval1 = zero
  ELSE
    IF (xr < xlow3) THEN
      fval1 = xr * (one-LOG(xr))
    ELSE
      t = (tcon*xr*xr-half) - half
      fval1 = xr * (cheval(nterm2,arlob2,t)-LOG(xr))
    END IF
  END IF
  lbpb21 = lobpia / (lobpib+lobpib)
  fval = (lbpb21-fval1) + lbpb22
END IF
lobpi1 = lobpia / lobpib

!   Compute value for argument in [pi/2,pi]

IF (indpi2 == 1) THEN
  fval = (lobpi1-fval) + lobpi2
END IF
fn_val = fval

!   Scale up for arguments > pi

IF (npi > 0) THEN
  fn_val = (fval+npi*lobpi2) + npi * lobpi1
END IF
IF (indsgn == -1) THEN
  fn_val = -fn_val
END IF
RETURN
END FUNCTION lobach



FUNCTION strom(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates Stromgren's integral, defined as

!      STROM(X) = integral 0 to X { t**7 exp(2t)/[exp(t)-1]**3 } dt

!    The code uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ASTROM to be used.
!                       The recommended value is such that
!                             ASTROM(NTERMS) < EPS/100

!    XLOW0 - REAL (dp) - The value below which STROM = 0.0 to machine
!                    precision. The recommended value is
!                          5th root of (130*XMIN)

!    XLOW1 - REAL (dp) - The value below which STROM = 3*(X**5)/(4*(pi**4))
!                   to machine precision. The recommended value is
!                             2*EPSNEG

!    EPSLN - REAL (dp) - The value of ln(EPS). Used to determine the no.
!                   of exponential terms for large X.

!    EPNGLN - REAL (dp) - The value of ln(EPSNEG). Used to prevent
!                    overflow for large X.

!    XHIGH - REAL (dp) - The value above which
!                           STROM = 196.52 - 15*(x**7)*exp(-x)/(4pi**4)
!                   to machine precision. The recommended value is
!                             7 / EPS

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!    EXP, INT, LOG

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp
REAL (dp)  :: epngln, epsln, rk, sumexp, sum2, t, x, xhigh, xk, xk1, xlow0, &
              xlow1
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'STROM '

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, four = 4.0_dp, seven = 7.0_dp,  &
                         onehun = 100.0_dp, one30 = 130.0_dp, one5ln = 0.4055_dp, &
                         f15bp4 = 0.38497433455066256959D-1,  &
                         pi4b3 = 1.29878788045336582982D2,    &
                         valinf = 196.51956920868988261257_dp

REAL (dp), PARAMETER  :: astrom(0:26) = (/  &
    0.56556120872539155290_dp, 0.4555731969101785525D-1,  &
   -0.4039535875936869170D-1, -0.133390572021486815D-2,   &
    0.185862506250538030D-2, -0.4685555868053659D-4, -0.6343475643422949D-4,  &
    0.572548708143200D-5, 0.159352812216822D-5, -0.28884328431036D-6,  &
   -0.2446633604801D-7, 0.1007250382374D-7, -0.12482986104D-9,  &
   -0.26300625283D-9, 0.2490407578D-10, 0.485454902D-11, -0.105378913D-11,  &
   -0.3604417D-13, 0.2992078D-13, -0.163971D-14, -0.61061D-15, 0.9335D-16,  &
    0.709D-17, -0.291D-17, 0.8D-19, 0.6D-19, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 26, 0, -1
    IF (ABS(astrom(nterms)) > t) EXIT
  END DO
  xlow0 = (one30*TINY(zero)) ** (one/(seven-two))
  xlow1 = two * xk
ELSE
  epsln = LOG(2*EPSILON(zero))
  epngln = LOG(xk)
  xhigh = seven / xk
END IF

!   Code for x < =  4.0

IF (x <= four) THEN
  IF (x < xlow0) THEN
    fn_val = zero
  ELSE
    IF (x < xlow1) THEN
      fn_val = (x**5) / pi4b3
    ELSE
      t = ((x/two)-half) - half
      fn_val = (x**5) * cheval(nterms,astrom,t) * f15bp4
    END IF
  END IF
ELSE
  
!  Code for x > 4.0
  
  IF (x > xhigh) THEN
    sumexp = one
  ELSE
    numexp = INT(epsln/(one5ln-x)) + 1
    IF (numexp > 1) THEN
      t = EXP(-x)
    ELSE
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, 7
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sum2 = sum2 * (rk+one) / two
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = seven * LOG(x) - x + LOG(sumexp)
  IF (t < epngln) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t) * f15bp4
  END IF
END IF
RETURN
END FUNCTION strom



FUNCTION strvh0(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the value of the Struve function
!      of order 0, denoted H0(x), for the argument XVALUE, defined

!         STRVHO(x) = (2/pi) integral{0 to pi/2} sin(x cos(t)) dt

!      H0 also satisfies the second-order equation

!                 x*D(Df) + Df + x*f = 2x/pi

!      The code uses Chebyshev expansions whose coefficients are given to 20D.


!   ERROR RETURNS:

!      As the asymptotic expansion of H0 involves the Bessel function
!      of the second kind Y0, there is a problem for large x, since
!      we cannot accurately calculate the value of Y0. An error message
!      is printed and STRVH0 returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - The no. of terms to be used in the array ARRH0. The
!               recommended value is such that
!                      ABS(ARRH0(NTERM1)) < EPS/100.

!      NTERM2 - The no. of terms to be used in the array ARRH0A. The
!               recommended value is such that
!                      ABS(ARRH0A(NTERM2)) < EPS/100.

!      NTERM3 - The no. of terms to be used in the array AY0ASP. The
!               recommended value is such that
!                      ABS(AY0ASP(NTERM3)) < EPS/100.

!      NTERM4 - The no. of terms to be used in the array AY0ASQ. The
!               recommended value is such that
!                      ABS(AY0ASQ(NTERM4)) < EPS/100.

!      XLOW - The value for which H0(x) = 2*x/pi to machine precision, if
!             abs(x) < XLOW. The recommended value is
!                      XLOW = 3 * SQRT(EPSNEG)

!      XHIGH - The value above which we are unable to calculate Y0 with
!              any reasonable accuracy. An error message is printed and
!              STRVH0 returns the value 0.0. The recommended value is
!                      XHIGH = 1/EPS.

!      For values of EPS and EPSNEG refer to the file MACHCON.TXT.

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      ABS, COS, SIN, SQRT.

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          ALLAN J. MACLEOD
!          DEPT. OF MATHEMATICS AND STATISTICS
!          UNIVERSITY OF PAISLEY
!          HIGH ST.
!          PAISLEY
!          SCOTLAND
!          PA1 2BE

!          (e-mail: macl_ms0@paisley.ac.uk )


!   LATEST REVISION:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: indsgn, nterm1, nterm2, nterm3, nterm4
REAL (dp)  :: h0as, t, x, xhigh, xlow, xmp4, xsq, y0p, y0q, y0val
CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
CHARACTER (LEN= 6)  :: fnname = 'STRVH0'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         eight = 8.0_dp, eleven = 11.0_dp, twenty = 20.0_dp,  &
                         onehun = 100.0_dp, sixtp5 = 60.5_dp,  &
                         two62 = 262.0_dp, thr2p5 = 302.5_dp,  &
                         piby4 = 0.78539816339744830962_dp,   &
                         rt2bpi = 0.79788456080286535588_dp,  &
                         twobpi = 0.63661977236758134308_dp

REAL (dp), PARAMETER  :: arrh0(0:19) = (/  &
    0.28696487399013225740_dp, -0.25405332681618352305_dp,  &
    0.20774026739323894439_dp, -0.20364029560386585140_dp,  &
    0.12888469086866186016_dp, -0.4825632815622261202D-1,  &
    0.1168629347569001242D-1, -0.198118135642418416D-2,  &
    0.24899138512421286D-3, -0.2418827913785950D-4, 0.187437547993431D-5,  &
   -0.11873346074362D-6, 0.626984943346D-8, -0.28045546793D-9,  &
    0.1076941205D-10, -0.35904793D-12, 0.1049447D-13, -0.27119D-15,  &
    0.624D-17, -0.13D-18 /)

REAL (dp), PARAMETER  :: arrh0a(0:20) = (/  &
    1.99291885751992305515_dp, -0.384232668701456887D-2,  &
   -0.32871993712353050D-3, -0.2941181203703409D-4, -0.267315351987066D-5,  &
   -0.24681031075013D-6, -0.2295014861143D-7, -0.215682231833D-8,  &
   -0.20303506483D-9, -0.1934575509D-10, -0.182773144D-11, -0.17768424D-12,  &
   -0.1643296D-13, -0.171569D-14, -0.13368D-15, -0.2077D-16, 0.2D-19,  &
   -0.55D-18, 0.10D-18, -0.4D-19, 0.1D-19 /)

REAL (dp), PARAMETER  :: ay0asp(0:12) = (/  &
    1.99944639402398271568D0, -0.28650778647031958D-3,  &
   -0.1005072797437620D-4, -0.35835941002463D-6, -0.1287965120531D-7,  &
   -0.46609486636D-9, -0.1693769454D-10, -0.61852269D-12, -0.2261841D-13,  &
   -0.83268D-15, -0.3042D-16, -0.115D-17, -0.4D-19 /)
REAL (dp), PARAMETER  :: ay0asq(0:13) = (/  &
    1.99542681386828604092D0, -0.236013192867514472D-2,  &
   -0.7601538908502966D-4, -0.256108871456343D-5, -0.8750292185106D-7,  &
   -0.304304212159D-8, -0.10621428314D-9, -0.377371479D-11, -0.13213687D-12,  &
   -0.488621D-14, -0.15809D-15, -0.762D-17, -0.3D-19, -0.3D-19 /)

!   Start computation

x = xvalue
indsgn = 1
IF (x < zero) THEN
  x = -x
  indsgn = -1
END IF

!   Compute the machine-dependent constants.

h0as = EPSILON(zero)
xhigh = one / (2*EPSILON(zero))

!   Error test

IF (ABS(xvalue) > xhigh) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   continue with machine constants

t = h0as / onehun
IF (x <= eleven) THEN
  DO  nterm1 = 19, 0, -1
    IF (ABS(arrh0(nterm1)) > t) EXIT
  END DO
  y0p = SQRT(h0as)
  xlow = y0p + y0p + y0p
ELSE
  DO  nterm2 = 20, 0, -1
    IF (ABS(arrh0a(nterm2)) > t) EXIT
  END DO
  DO  nterm3 = 12, 0, -1
    IF (ABS(ay0asp(nterm3)) > t) EXIT
  END DO
  DO  nterm4 = 13, 0, -1
    IF (ABS(ay0asq(nterm4)) > t) EXIT
  END DO
END IF

!   Code for abs(x) <= 11

IF (x <= eleven) THEN
  IF (x < xlow) THEN
    fn_val = twobpi * x
  ELSE
    t = ((x*x)/sixtp5-half) - half
    fn_val = twobpi * x * cheval(nterm1,arrh0,t)
  END IF
ELSE
  
!   Code for abs(x) > 11
  
  xsq = x * x
  t = (two62-xsq) / (twenty+xsq)
  y0p = cheval(nterm3,ay0asp,t)
  y0q = cheval(nterm4,ay0asq,t) / (eight*x)
  xmp4 = x - piby4
  y0val = y0p * SIN(xmp4) - y0q * COS(xmp4)
  y0val = y0val * rt2bpi / SQRT(x)
  t = (thr2p5-xsq) / (sixtp5+xsq)
  h0as = twobpi * cheval(nterm2,arrh0a,t) / x
  fn_val = y0val + h0as
END IF
IF (indsgn == -1) fn_val = -fn_val
RETURN
END FUNCTION strvh0



FUNCTION strvh1(xvalue) RESULT(fn_val)

!   DESCRIPTION:
!      This function calculates the value of the Struve function
!      of order 1, denoted H1(x), for the argument XVALUE, defined as

!                                                                  2
!        STRVH1(x) = (2x/pi) integral{0 to pi/2} sin( x cos(t))*sin t dt

!      H1 also satisfies the second-order differential equation

!                    2   2                   2            2
!                   x * D f  +  x * Df  +  (x - 1)f  =  2x / pi

!      The code uses Chebyshev expansions with the coefficients given to 20D.


!   ERROR RETURNS:
!      As the asymptotic expansion of H1 involves the Bessel function
!      of the second kind Y1, there is a problem for large x, since
!      we cannot accurately calculate the value of Y1. An error message
!      is printed and STRVH1 returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - The no. of terms to be used in the array ARRH1. The
!               recommended value is such that
!                      ABS(ARRH1(NTERM1)) < EPS/100.

!      NTERM2 - The no. of terms to be used in the array ARRH1A. The
!               recommended value is such that
!                      ABS(ARRH1A(NTERM2)) < EPS/100.

!      NTERM3 - The no. of terms to be used in the array AY1ASP. The
!               recommended value is such that
!                      ABS(AY1ASP(NTERM3)) < EPS/100.

!      NTERM4 - The no. of terms to be used in the array AY1ASQ. The
!               recommended value is such that
!                      ABS(AY1ASQ(NTERM4)) < EPS/100.

!      XLOW1 - The value of x, below which H1(x) set to zero, if
!              abs(x)<XLOW1. The recommended value is
!                      XLOW1 = 1.5 * SQRT(XMIN)

!      XLOW2 - The value for which H1(x) = 2*x*x/pi to machine precision, if
!              abs(x) < XLOW2. The recommended value is
!                      XLOW2 = SQRT(15*EPSNEG)

!      XHIGH - The value above which we are unable to calculate Y1 with
!              any reasonable accuracy. An error message is printed and
!              STRVH1 returns the value 0.0.  The recommended value is
!                      XHIGH = 1/EPS.

!      For values of EPS, EPSNEG and XMIN refer to the file MACHCON.TXT.

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      ABS, COS, SIN, SQRT.

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          ALLAN J. MACLEOD
!          DEPT. OF MATHEMATICS AND STATISTICS
!          UNIVERSITY OF PAISLEY
!          HIGH ST.
!          PAISLEY
!          SCOTLAND
!          PA1 2BE

!          (e-mail: macl_ms0@paisley.ac.uk)


!   LATEST REVISION:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3, nterm4
REAL (dp)  :: h1as, t, x, xhigh, xlow1, xlow2, xm3p4, xsq, y1p, y1q, y1val
CHARACTER (LEN=26)  :: errmsg = 'ARGUMENT TOO LARGE IN SIZE'
CHARACTER (LEN= 6)  :: fnname = 'STRVH1'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, eight = 8.0_dp,  &
                         nine = 9.0_dp, fiften = 15.0_dp, twenty = 20.0_dp, &
                         onehun = 100.0_dp, fortp5 = 40.5_dp,  &
                         one82 = 182.0_dp, tw02p5 = 202.5_dp,  &
                         rt2bpi = 0.79788456080286535588_dp,  &
                         thpby4 = 2.35619449019234492885_dp,  &
                         twobpi = 0.63661977236758134308_dp

REAL (dp), PARAMETER  :: arrh1(0:17) = (/  &
    0.17319061083675439319_dp, -0.12606917591352672005_dp,  &
    0.7908576160495357500D-1, -0.3196493222321870820D-1,  &
    0.808040581404918834D-2, -0.136000820693074148D-2,  &
    0.16227148619889471D-3, -0.1442352451485929D-4,  &
    0.99219525734072D-6, -0.5441628049180D-7, 0.243631662563D-8,  &
    -0.9077071338D-10, 0.285926585D-11, -0.7716975D-13,  &
    0.180489D-14, -0.3694D-16, 0.67D-18, -0.1D-19 /)

REAL (dp), PARAMETER  :: arrh1a(0:21) = (/  &
    2.01083504951473379407_dp, 0.592218610036099903D-2,  &
    0.55274322698414130D-3, 0.5269873856311036D-4, 0.506374522140969D-5,  &
    0.49028736420678D-6, 0.4763540023525D-7, 0.465258652283D-8,  &
    0.45465166081D-9, 0.4472462193D-10, 0.437308292D-11, 0.43568368D-12,  &
    0.4182190D-13, 0.441044D-14, 0.36391D-15, 0.5558D-16, -0.4D-19,  &
    0.163D-17, -0.34D-18, 0.13D-18, -0.4D-19, 0.1D-19 /)

REAL (dp), PARAMETER  :: ay1asp(0:14) = (/  &
    2.00135240045889396402_dp, 0.71104241596461938D-3,  &
    0.3665977028232449D-4, 0.191301568657728D-5,  &
    0.10046911389777D-6, 0.530401742538D-8, 0.28100886176D-9,  &
    0.1493886051D-10, 0.79578420D-12, 0.4252363D-13, 0.227195D-14,  &
    0.12216D-15, 0.650D-17, 0.36D-18, 0.2D-19 /)
REAL (dp), PARAMETER  :: ay1asq(0:15) = (/  &
    5.99065109477888189116_dp, -0.489593262336579635D-2,  &
   -0.23238321307070626D-3, -0.1144734723857679D-4, -0.57169926189106D-6,  &
   -0.2895516716917D-7, -0.147513345636D-8, -0.7596537378D-10,  &
   -0.390658184D-11, -0.20464654D-12, -0.1042636D-13, -0.57702D-15,  &
   -0.2550D-16, -0.210D-17, 0.2D-19, -0.2D-19 /)

!   Start computation

x = ABS(xvalue)

!   Compute the machine-dependent constants.

xhigh = (half+half) / (2*EPSILON(zero))

!   Error test

IF (x > xhigh) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   continue with machine constants

h1as = EPSILON(zero)
t = h1as / onehun
IF (x <= nine) THEN
  DO  nterm1 = 17, 0, -1
    IF (ABS(arrh1(nterm1)) > t) EXIT
  END DO
  xlow1 = half * SQRT(TINY(zero))
  xlow1 = xlow1 + xlow1 + xlow1
  xlow2 = SQRT(fiften*h1as)
ELSE
  DO  nterm2 = 21, 0, -1
    IF (ABS(arrh1a(nterm2)) > t) EXIT
  END DO
  DO  nterm3 = 14, 0, -1
    IF (ABS(ay1asp(nterm3)) > t) EXIT
  END DO
  DO  nterm4 = 15, 0, -1
    IF (ABS(ay1asq(nterm4)) > t) EXIT
  END DO
END IF

!   Code for abs(x) <= 9

IF (x <= nine) THEN
  IF (x < xlow1) THEN
    fn_val = zero
  ELSE
    xsq = x * x
    IF (x < xlow2) THEN
      fn_val = twobpi * xsq
    ELSE
      t = (xsq/fortp5-half) - half
      fn_val = twobpi * xsq * cheval(nterm1,arrh1,t)
    END IF
  END IF
ELSE
  
!   Code for abs(x) > 9
  
  xsq = x * x
  t = (one82-xsq) / (twenty+xsq)
  y1p = cheval(nterm3,ay1asp,t)
  y1q = cheval(nterm4,ay1asq,t) / (eight*x)
  xm3p4 = x - thpby4
  y1val = y1p * SIN(xm3p4) + y1q * COS(xm3p4)
  y1val = y1val * rt2bpi / SQRT(x)
  t = (tw02p5-xsq) / (fortp5+xsq)
  h1as = twobpi * cheval(nterm2,arrh1a,t)
  fn_val = y1val + h1as
END IF
RETURN
END FUNCTION strvh1



FUNCTION strvl0(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the modified Struve function of
!      order 0, denoted L0(x), defined as the solution of the
!      second-order equation

!                  x*D(Df) + Df - x*f  =  2x/pi


!   ERROR RETURNS:

!      If the value of |XVALUE| is too large, the result
!      would cause an floating-pt overflow. An error message
!      is printed and the function returns the value of
!      sign(XVALUE)*XMAX where XMAX is the largest possible
!      floating-pt argument.


!   MACHINE-DEPENDENT PARAMETERS:

!      NTERM1 - INTEGER - The no. of terms for the array ARL0.
!                         The recommended value is such that
!                             ABS(ARL0(NTERM1)) < EPS/100

!      NTERM2 - INTEGER - The no. of terms for the array ARL0AS.
!                         The recommended value is such that
!                             ABS(ARL0AS(NTERM2)) < EPS/100

!      NTERM3 - INTEGER - The no. of terms for the array AI0ML0.
!                         The recommended value is such that
!                             ABS(AI0ML0(NTERM3)) < EPS/100

!      XLOW - REAL (dp) - The value of x below which L0(x) = 2*x/pi
!                    to machine precision. The recommended value is
!                             3*SQRT(EPS)

!      XHIGH1 - REAL (dp) - The value beyond which the Chebyshev series
!                      in the asymptotic expansion of I0 - L0 gives
!                      1.0 to machine precision. The recommended value
!                      is   SQRT( 30/EPSNEG )

!      XHIGH2 - REAL (dp) - The value beyond which the Chebyshev series
!                      in the asymptotic expansion of I0 gives 1.0
!                      to machine precision. The recommended value
!                      is   28 / EPSNEG

!      XMAX - REAL (dp) - The value of XMAX, where XMAX is the
!                    largest possible floating-pt argument.
!                    This is used to prevent overflow.

!      For values of EPS, EPSNEG and XMAX the user should refer
!      to the file MACHCON.TXT

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      EXP , LOG , SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          DR. ALLAN J. MACLEOD
!          DEPT. OF MATHEMATICS AND STATISTICS
!          UNIVERSITY OF PAISLEY
!          HIGH ST.
!          PAISLEY
!          SCOTLAND
!          PA1 2BE

!      (e-mail: macl_ms0@paisley.ac.uk )


!   LATEST REVISION:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: indsgn, nterm1, nterm2, nterm3
REAL (dp)  :: ch1, ch2, t, test, x, xhigh1, xhigh2, xlow, xmax, xsq
CHARACTER (LEN=24)  :: errmsg = 'ARGUMENT CAUSES OVERFLOW'
CHARACTER (LEN= 6)  :: fnname = 'STRVL0'

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp,  &
                         four = 4.0_dp, sixten = 16.0_dp, twent4 = 24.0_dp,  &
                         twent8 = 28.0_dp, onehun = 100.0_dp,  &
                         two88 = 288.0_dp, atehun = 800.0_dp,  &
                         lnr2pi = 0.91893853320467274178_dp,  &
                         twobpi = 0.63661977236758134308_dp

REAL (dp), PARAMETER  :: arl0(0:27) = (/  &
    0.42127458349979924863_dp, -0.33859536391220612188_dp,  &
    0.21898994812710716064_dp, -0.12349482820713185712_dp,  &
    0.6214209793866958440D-1, -0.2817806028109547545D-1,  &
    0.1157419676638091209D-1, -0.431658574306921179D-2,  &
    0.146142349907298329D-2, -0.44794211805461478D-3, 0.12364746105943761D-3, &
   -0.3049028334797044D-4, 0.663941401521146D-5, -0.125538357703889D-5,  &
    0.20073446451228D-6, -0.2588260170637D-7, 0.241143742758D-8,  &
   -0.10159674352D-9, -0.1202430736D-10, 0.262906137D-11, -0.15313190D-12,  &
   -0.1574760D-13, 0.315635D-14, -0.4096D-16, -0.3620D-16, 0.239D-17,  &
    0.36D-18, -0.4D-19 /)

REAL (dp), PARAMETER  :: arl0as(0:15) = (/  &
    2.00861308235605888600_dp, 0.403737966500438470D-2,  &
   -0.25199480286580267D-3, 0.1605736682811176D-4, -0.103692182473444D-5,  &
    0.6765578876305D-7, -0.444999906756D-8, 0.29468889228D-9,  &
   -0.1962180522D-10, 0.131330306D-11, -0.8819190D-13, 0.595376D-14,  &
   -0.40389D-15, 0.2651D-16, -0.208D-17, 0.11D-18 /)

REAL (dp), PARAMETER  :: ai0ml0(0:23) = (/  &
    2.00326510241160643125_dp, 0.195206851576492081D-2,  &
    0.38239523569908328D-3, 0.7534280817054436D-4, 0.1495957655897078D-4,  &
    0.299940531210557D-5, 0.60769604822459D-6, 0.12399495544506D-6,  &
    0.2523262552649D-7, 0.504634857332D-8, 0.97913236230D-9, 0.18389115241D-9, &
    0.3376309278D-10, 0.611179703D-11, 0.108472972D-11, 0.18861271D-12,  &
    0.3280345D-13, 0.565647D-14, 0.93300D-15, 0.15881D-15, 0.2791D-16,  &
    0.389D-17, 0.70D-18, 0.16D-18 /)

!   Start computation

x = xvalue
indsgn = 1
IF (x < zero) THEN
  x = -x
  indsgn = -1
END IF

!   Compute the machine-dependent constants.

test = EPSILON(zero)
t = test / onehun
IF (x <= sixten) THEN
  DO  nterm1 = 27, 0, -1
    IF (ABS(arl0(nterm1)) > t) EXIT
  END DO
  xlow = (one+two) * SQRT(test)
ELSE
  DO  nterm2 = 15, 0, -1
    IF (ABS(arl0as(nterm2)) > t) EXIT
  END DO
  DO  nterm3 = 23, 0, -1
    IF (ABS(ai0ml0(nterm3)) > t) EXIT
  END DO
  xmax = HUGE(zero)
  xhigh1 = SQRT((twent8+two)/test)
  xhigh2 = twent8 / test
END IF

!   Code for |xvalue| <= 16

IF (x <= sixten) THEN
  IF (x < xlow) THEN
    fn_val = twobpi * x
  ELSE
    t = (four*x-twent4) / (x+twent4)
    fn_val = twobpi * x * cheval(nterm1,arl0,t) * EXP(x)
  END IF
ELSE
  
!   Code for |xvalue| > 16
  
  IF (x > xhigh2) THEN
    ch1 = one
  ELSE
    t = (x-twent8) / (four-x)
    ch1 = cheval(nterm2,arl0as,t)
  END IF
  IF (x > xhigh1) THEN
    ch2 = one
  ELSE
    xsq = x * x
    t = (atehun-xsq) / (two88+xsq)
    ch2 = cheval(nterm3,ai0ml0,t)
  END IF
  test = LOG(ch1) - lnr2pi - LOG(x) / two + x
  IF (test > LOG(xmax)) THEN
    CALL errprn(fnname,errmsg)
    fn_val = xmax
  ELSE
    fn_val = EXP(test) - twobpi * ch2 / x
  END IF
END IF
IF (indsgn == -1) fn_val = -fn_val
RETURN
END FUNCTION strvl0



FUNCTION strvl1(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the modified Struve function of
!      order 1, denoted L1(x), defined as the solution of

!               x*x*D(Df) + x*Df - (x*x+1)f = 2*x*x/pi


!   ERROR RETURNS:

!      If the value of |XVALUE| is too large, the result
!      would cause an floating-pt overflow. An error message
!      is printed and the function returns the value of
!      sign(XVALUE)*XMAX where XMAX is the largest possible
!      floating-pt argument.


!   MACHINE-DEPENDENT PARAMETERS:

!      NTERM1 - INTEGER - The no. of terms for the array ARL1.
!                         The recommended value is such that
!                             ABS(ARL1(NTERM1)) < EPS/100

!      NTERM2 - INTEGER - The no. of terms for the array ARL1AS.
!                         The recommended value is such that
!                             ABS(ARL1AS(NTERM2)) < EPS/100

!      NTERM3 - INTEGER - The no. of terms for the array AI1ML1.
!                         The recommended value is such that
!                             ABS(AI1ML1(NTERM3)) < EPS/100

!      XLOW1 - REAL (dp) - The value of x below which L1(x) = 2*x*x/(3*pi)
!                     to machine precision. The recommended value is
!                              SQRT(15*EPS)

!      XLOW2 - REAL (dp) - The value of x below which L1(x) set to 0.0.
!                     This is used to prevent underflow. The
!                     recommended value is
!                              SQRT(5*XMIN)

!      XHIGH1 - REAL (dp) - The value of |x| above which the Chebyshev
!                      series in the asymptotic expansion of I1
!                      equals 1.0 to machine precision. The
!                      recommended value is  SQRT( 30 / EPSNEG ).

!      XHIGH2 - REAL (dp) - The value of |x| above which the Chebyshev
!                      series in the asymptotic expansion of I1 - L1
!                      equals 1.0 to machine precision. The recommended
!                      value is   30 / EPSNEG.

!      XMAX - REAL (dp) - The value of XMAX, where XMAX is the
!                    largest possible floating-pt argument.
!                    This is used to prevent overflow.

!      For values of EPS, EPSNEG, XMIN, and XMAX the user should refer
!      to the file MACHCON.TXT

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      EXP , LOG , SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          DR. ALLAN J. MACLEOD
!          DEPT. OF MATHEMATICS AND STATISTICS
!          UNIVERSITY OF PAISLEY
!          HIGH ST.
!          PAISLEY
!          SCOTLAND
!          PA1 2BE

!          (e-mail: macl_ms0@paisley.ac.uk )


!   LATEST UPDATE:
!                 23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3
REAL (dp)  :: ch1, ch2, t, test, x, xhigh1, xhigh2, xlow1, xlow2, xmax, xsq
CHARACTER (LEN=24)  :: errmsg = 'ARGUMENT CAUSES OVERFLOW'
CHARACTER (LEN= 6)  :: fnname = 'STRVL1'

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp,  &
                         four = 4.0_dp, sixten = 16.0_dp, twent4 = 24.0_dp,  &
                         thirty = 30.0_dp, onehun = 100.0_dp, two88 = 288.0_dp, &
                         atehun = 800.0_dp,  &
                         lnr2pi = 0.91893853320467274178_dp, &
                         pi3by2 = 4.71238898038468985769_dp, &
                         twobpi = 0.63661977236758134308_dp

REAL (dp), PARAMETER  :: arl1(0:26) = (/  &
    0.38996027351229538208_dp, -0.33658096101975749366_dp,  &
    0.23012467912501645616_dp, -0.13121594007960832327_dp,  &
    0.6425922289912846518D-1, -0.2750032950616635833D-1,  &
    0.1040234148637208871D-1, -0.350532294936388080D-2,  &
    0.105748498421439717D-2, -0.28609426403666558D-3, 0.6925708785942208D-4,  &
   -0.1489693951122717D-4, 0.281035582597128D-5, -0.45503879297776D-6,  &
    0.6090171561770D-7, -0.623543724808D-8, 0.38430012067D-9, 0.790543916D-11, &
   -0.489824083D-11, 0.46356884D-12, 0.684205D-14, -0.569748D-14,  &
    0.35324D-15, 0.4244D-16, -0.644D-17, -0.21D-18, 0.9D-19 /)

REAL (dp), PARAMETER  :: arl1as(0:16) = (/  &
    1.97540378441652356868_dp, -0.1195130555088294181D-1,  &
    0.33639485269196046D-3, -0.1009115655481549D-4, 0.30638951321998D-6,  &
   -0.953704370396D-8, 0.29524735558D-9, -0.951078318D-11, 0.28203667D-12,  &
   -0.1134175D-13, 0.147D-17, -0.6232D-16, -0.751D-17, -0.17D-18, 0.51D-18,  &
    0.23D-18, 0.5D-19 /)

REAL (dp), PARAMETER  :: ai1ml1(0:25) = (/  &
    1.99679361896789136501_dp, -0.190663261409686132D-2,  &
   -0.36094622410174481D-3, -0.6841847304599820D-4, -0.1299008228509426D-4,  &
   -0.247152188705765D-5, -0.47147839691972D-6, -0.9020819982592D-7,  &
   -0.1730458637504D-7, -0.332323670159D-8, -0.63736421735D-9,  &
   -0.12180239756D-9, -0.2317346832D-10, -0.439068833D-11, -0.82847110D-12,  &
   -0.15562249D-12, -0.2913112D-13, -0.543965D-14, -0.101177D-14,  &
   -0.18767D-15, -0.3484D-16, -0.643D-17, -0.118D-17, -0.22D-18, -0.4D-19,  &
   -0.1D-19 /)

!   START CALCULATION

x = ABS(xvalue)

!   Compute the machine-dependent constants.

test = EPSILON(zero)
t = test / onehun
IF (x <= sixten) THEN
  DO  nterm1 = 26, 0, -1
    IF (ABS(arl1(nterm1)) > t) EXIT
  END DO
  xlow1 = SQRT(thirty*test/two)
  xlow2 = SQRT((four+one)*TINY(zero))
ELSE
  DO  nterm2 = 16, 0, -1
    IF (ABS(arl1as(nterm2)) > t) EXIT
  END DO
  DO  nterm3 = 25, 0, -1
    IF (ABS(ai1ml1(nterm3)) > t) EXIT
  END DO
  xmax = HUGE(zero)
  xhigh2 = thirty / test
  xhigh1 = SQRT(xhigh2)
END IF

!   CODE FOR |XVALUE| <= 16

IF (x <= sixten) THEN
  IF (x <= xlow2) THEN
    fn_val = zero
  ELSE
    xsq = x * x
    IF (x < xlow1) THEN
      fn_val = xsq / pi3by2
    ELSE
      t = (four*x-twent4) / (x+twent4)
      fn_val = xsq * cheval(nterm1,arl1,t) * EXP(x) / pi3by2
    END IF
  END IF
ELSE
  
!   CODE FOR |XVALUE| > 16
  
  IF (x > xhigh2) THEN
    ch1 = one
  ELSE
    t = (x-thirty) / (two-x)
    ch1 = cheval(nterm2,arl1as,t)
  END IF
  IF (x > xhigh1) THEN
    ch2 = one
  ELSE
    xsq = x * x
    t = (atehun-xsq) / (two88+xsq)
    ch2 = cheval(nterm3,ai1ml1,t)
  END IF
  test = LOG(ch1) - lnr2pi - LOG(x) / two + x
  IF (test > LOG(xmax)) THEN
    CALL errprn(fnname,errmsg)
    fn_val = xmax
  ELSE
    fn_val = EXP(test) - twobpi * ch2
  END IF
END IF
RETURN
END FUNCTION strvl1



FUNCTION synch1(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the synchrotron radiation function defined as

!         SYNCH1(x) = x * Integral{x to inf} K(5/3)(t) dt,

!      where K(5/3) is a modified Bessel function of order 5/3.

!      The code uses Chebyshev expansions, the coefficients of which
!      are given to 20 decimal places.


!   ERROR RETURNS:

!      The function is undefined if x < 0.0. If XVALUE < 0.0,
!      an error message is printed and the function returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - INTEGER - The no. of terms needed from the array
!                         ASYNC1. The recommended value is such that
!                            ABS(ASYNC1(NTERM1)) < EPS/100.

!      NTERM2 - INTEGER - The no. of terms needed from the array
!                         ASYNC2. The recommended value is such that
!                            ABS(ASYNC2(NTERM2)) < EPS/100.

!      NTERM3 - INTEGER - The no. of terms needed from the array
!                         ASYNCA. The recommended value is such that
!                            ABS(ASYNCA(NTERM3)) < EPS/100.

!      XLOW - REAL (dp) - The value below which
!                        SYNCH1(x) = 2.14952.. * (x**(1/3))
!                    to machine precision. The recommended value
!                    is     sqrt (8*EPSNEG)

!      XHIGH1 - REAL (dp) - The value above which
!                          SYNCH1(x) = 0.0
!                      to machine precision. The recommended value
!                      is     -8*LN(XMIN)/7

!      XHIGH2 - REAL (dp) - The value of LN(XMIN). This is used
!                      to prevent underflow in calculations
!                      for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      EXP , LOG , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley,
!          Paisley,
!          SCOTLAND
!          PA1 2BE

!          ( e-mail:  macl_ms0@paisley.ac.uk )


!   LATEST UPDATE:
!                  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3
REAL (dp)  :: cheb1, cheb2, t, x, xhigh1, xhigh2, xlow, xpowth
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'SYNCH1'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         three = 3.0_dp, four = 4.0_dp, eight = 8.0_dp,  &
                         twelve = 12.0_dp, onehun = 100.0_dp,  &
                         conlow = 2.14952824153447863671_dp,  &
                         pibrt3 = 1.81379936423421785059_dp,  &
                         lnrtp2 = 0.22579135264472743236_dp

REAL (dp), PARAMETER  :: async1(0:13) = (/  &
    30.36468298250107627340_dp, 17.07939527740839457449_dp,  &
    4.56013213354507288887_dp, 0.54928124673041997963_dp,  &
    0.3729760750693011724D-1, 0.161362430201041242D-2,  &
    0.4819167721203707D-4, 0.105124252889384D-5,  &
    0.1746385046697D-7, 0.22815486544D-9, 0.240443082D-11,  &
    0.2086588D-13, 0.15167D-15, 0.94D-18 /)

REAL (dp), PARAMETER  :: async2(0:11) = (/  &
    0.44907216235326608443_dp, 0.8983536779941872179D-1,  &
    0.810445737721512894D-2, 0.42617169910891619D-3,  &
    0.1476096312707460D-4, 0.36286336153998D-6, 0.666348074984D-8,  &
    0.9490771655D-10, 0.107912491D-11, 0.1002201D-13, 0.7745D-16, 0.51D-18 /)

REAL (dp), PARAMETER  :: asynca(0:24) = (/  &
    2.13293051613550009848_dp, 0.7413528649542002401D-1,  &
    0.869680999099641978D-2, 0.117038262487756921D-2, 0.16451057986191915D-3, &
    0.2402010214206403D-4, 0.358277563893885D-5, 0.54477476269837D-6,  &
    0.8388028561957D-7, 0.1306988268416D-7, 0.205309907144D-8,  &
    0.32518753688D-9, 0.5179140412D-10, 0.830029881D-11, 0.133527277D-11,  &
    0.21591498D-12, 0.3499673D-13, 0.569942D-14, 0.92906D-15, 0.15222D-15,  &
    0.2491D-16, 0.411D-17, 0.67D-18, 0.11D-18, 0.2D-19 /)

!   Start calculation

x = xvalue
IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

cheb1 = EPSILON(zero)
t = cheb1 / onehun
IF (x <= four) THEN
  DO  nterm1 = 13, 0, -1
    IF (ABS(async1(nterm1)) > t) EXIT
  END DO
  DO  nterm2 = 11, 0, -1
    IF (ABS(async2(nterm2)) > t) EXIT
  END DO
  xlow = SQRT(eight*cheb1)
ELSE
  DO  nterm3 = 24, 0, -1
    IF (ABS(asynca(nterm3)) > t) EXIT
  END DO
  xhigh2 = LOG(TINY(zero))
  xhigh1 = -eight * xhigh2 / (eight-one)
END IF

!   Code for 0 <= x <= 4

IF (x <= four) THEN
  xpowth = x ** (one/three)
  IF (x < xlow) THEN
    fn_val = conlow * xpowth
  ELSE
    t = (x*x/eight-half) - half
    cheb1 = cheval(nterm1,async1,t)
    cheb2 = cheval(nterm2,async2,t)
    t = xpowth * cheb1 - (xpowth**11) * cheb2
    fn_val = t - pibrt3 * x
  END IF
ELSE
  IF (x > xhigh1) THEN
    fn_val = zero
  ELSE
    t = (twelve-x) / (x+four)
    cheb1 = cheval(nterm3,asynca,t)
    t = lnrtp2 - x + LOG(SQRT(x)*cheb1)
    IF (t < xhigh2) THEN
      fn_val = zero
    ELSE
      fn_val = EXP(t)
    END IF
  END IF
END IF
RETURN
END FUNCTION synch1



FUNCTION synch2(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the synchrotron radiation function
!      defined as

!         SYNCH2(x) = x * K(2/3)(x)

!      where K(2/3) is a modified Bessel function of order 2/3.

!      The code uses Chebyshev expansions, the coefficients of which
!      are given to 20 decimal places.


!   ERROR RETURNS:

!      The function is undefined if x < 0.0. If XVALUE < 0.0,
!      an error message is printed and the function returns the value 0.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - INTEGER - The no. of terms needed from the array
!                         ASYNC1. The recommended value is such that
!                            ABS(ASYN21(NTERM1)) < EPS/100.

!      NTERM2 - INTEGER - The no. of terms needed from the array
!                         ASYNC2. The recommended value is such that
!                            ABS(ASYN22(NTERM2)) < EPS/100.

!      NTERM3 - INTEGER - The no. of terms needed from the array
!                         ASYNCA. The recommended value is such that
!                            ABS(ASYN2A(NTERM3)) < EPS/100.

!      XLOW - REAL (dp) - The value below which
!                        SYNCH2(x) = 1.074764... * (x**(1/3))
!                    to machine precision. The recommended value
!                    is     sqrt (8*EPSNEG)

!      XHIGH1 - REAL (dp) - The value above which
!                          SYNCH2(x) = 0.0
!                      to machine precision. The recommended value
!                      is     -8*LN(XMIN)/7

!      XHIGH2 - REAL (dp) - The value of LN(XMIN). This is used
!                      to prevent underflow in calculations
!                      for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      EXP , LOG , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley,
!          Paisley,
!          SCOTLAND
!          PA1 2BE

!          (e-mail:  macl_ms0@paisley.ac.uk)


!   LATEST UPDATE:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3
REAL (dp)  :: cheb1, cheb2, t, x, xhigh1, xhigh2, xlow, xpowth
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'SYNCH2'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         two = 2.0_dp, three = 3.0_dp, four = 4.0_dp,  &
                         eight = 8.0_dp, ten = 10.0_dp, onehun = 100.0_dp,  &
                         conlow = 1.07476412076723931836_dp,  &
                         lnrtp2 = 0.22579135264472743236_dp

REAL (dp), PARAMETER  :: asyn21(0:14) = (/  &
    38.61783992384308548014_dp, 23.03771559496373459697_dp,  &
    5.38024998683357059676_dp, 0.61567938069957107760_dp,  &
    0.4066880046688955843D-1, 0.172962745526484141D-2,  &
    0.5106125883657699D-4, 0.110459595022012D-5,  &
    0.1823553020649D-7, 0.23707698034D-9, 0.248872963D-11,  &
    0.2152868D-13, 0.15607D-15, 0.96D-18, 0.1D-19 /)

REAL (dp), PARAMETER  :: asyn22(0:13) = (/  &
    7.90631482706608042875_dp, 3.13534636128534256841_dp,  &
    0.48548794774537145380_dp, 0.3948166758272372337D-1,  &
    0.196616223348088022D-2, 0.6590789322930420D-4,  &
    0.158575613498559D-5, 0.2868653011233D-7, 0.40412023595D-9,  &
    0.455684443D-11, 0.4204590D-13, 0.32326D-15, 0.210D-17, 0.1D-19 /)

REAL (dp), PARAMETER  :: asyn2a(0:18) = (/  &
    2.02033709417071360032_dp, 0.1095623712180740443D-1,  &
    0.85423847301146755D-3, 0.7234302421328222D-4,  &
    0.631244279626992D-5, 0.56481931411744D-6, 0.5128324801375D-7,  &
    0.471965329145D-8, 0.43807442143D-9, 0.4102681493D-10,  &
    0.386230721D-11, 0.36613228D-12, 0.3480232D-13, 0.333010D-14,  &
    0.31856D-15, 0.3074D-16, 0.295D-17, 0.29D-18, 0.3D-19 /)

!   Start calculation

x = xvalue
IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

cheb1 = EPSILON(zero)
t = cheb1 / onehun
IF (x <= four) THEN
  DO  nterm1 = 14, 0, -1
    IF (ABS(asyn21(nterm1)) > t) EXIT
  END DO
  DO  nterm2 = 13, 0, -1
    IF (ABS(asyn22(nterm2)) > t) EXIT
  END DO
  xlow = SQRT(eight*cheb1)
ELSE
  DO  nterm3 = 18, 0, -1
    IF (ABS(asyn2a(nterm3)) > t) EXIT
  END DO
  xhigh2 = LOG(TINY(zero))
  xhigh1 = -eight * xhigh2 / (eight-one)
END IF

!   Code for 0 <= x <= 4

IF (x <= four) THEN
  xpowth = x ** (one/three)
  IF (x < xlow) THEN
    fn_val = conlow * xpowth
  ELSE
    t = (x*x/eight-half) - half
    cheb1 = cheval(nterm1,asyn21,t)
    cheb2 = cheval(nterm2,asyn22,t)
    fn_val = xpowth * cheb1 - (xpowth**5) * cheb2
  END IF
ELSE
  IF (x > xhigh1) THEN
    fn_val = zero
  ELSE
    t = (ten-x) / (x+two)
    cheb1 = cheval(nterm3,asyn2a,t)
    t = lnrtp2 - x + LOG(SQRT(x)*cheb1)
    IF (t < xhigh2) THEN
      fn_val = zero
    ELSE
      fn_val = EXP(t)
    END IF
  END IF
END IF
RETURN
END FUNCTION synch2



FUNCTION tran02(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates the transport integral of order 2, defined as

!      TRAN02(X) = integral 0 to X { t**2 exp(t)/[exp(t)-1]**2 } dt

!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
!                       The recommended value is such that
!                             ATRAN(NTERMS) < EPS/100

!    XLOW1 - REAL (dp) - The value below which TRAN02 = x to
!                   machine precision. The recommended value is
!                             sqrt(8*EPSNEG)

!    XHIGH1 - REAL (dp) - The value above which the exponential series for
!                    large x contains only one term. The recommended value
!                    is        - ln(EPS).

!    XHIGH2 - REAL (dp) - The value above which
!                       TRAN02 = VALINF  -  x**2 exp(-x)
!                    The recommended value is 2/EPS

!    XHIGH3 - REAL (dp) - The value of ln(EPSNEG). Used to prevent overflow
!                    for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!     EXP, INT, LOG, SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY ,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp, numjn = 2
REAL (dp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1, xlow1
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'TRAN02'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, onehun = 100.0_dp,  &
                         rnumjn = 2.0_dp, valinf = 3.2898681336964528729_dp

REAL (dp), PARAMETER  :: atran(0:19) = (/  &
    1.67176044643453850301_dp, -0.14773535994679448986_dp,  &
    0.1482138199469363384D-1, -0.141953303263056126D-2,  &
    0.13065413244157083D-3, -0.1171557958675790D-4,  &
    0.103334984457557D-5, -0.9019113042227D-7, 0.781771698331D-8,  &
   -0.67445656840D-9, 0.5799463945D-10, -0.497476185D-11,  &
    0.42596097D-12, -0.3642189D-13, 0.311086D-14, -0.26547D-15,  &
    0.2264D-16, -0.193D-17, 0.16D-18, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 19, 0, -1
    IF (ABS(atran(nterms)) > t) EXIT
  END DO
  xlow1 = SQRT(eight*xk)
ELSE
  xhigh1 = -LOG(2*EPSILON(zero))
  xhigh2 = one / (half*xk)
  xhigh3 = LOG(xk)
END IF

!   Code for x < =  4.0

IF (x <= four) THEN
  IF (x < xlow1) THEN
    fn_val = (x**(numjn-1)) / (rnumjn-one)
  ELSE
    t = (((x*x)/eight)-half) - half
    fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
  END IF
ELSE
  
!  Code for x > 4.0
  
  IF (x > xhigh2) THEN
    sumexp = one
  ELSE
    IF (x <= xhigh1) THEN
      numexp = INT(xhigh1/x) + 1
      t = EXP(-x)
    ELSE
      numexp = 1
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, numjn
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = rnumjn * LOG(x) - x + LOG(sumexp)
  IF (t < xhigh3) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t)
  END IF
END IF
RETURN
END FUNCTION tran02



FUNCTION tran03(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates the transport integral of order 3, defined as

!      TRAN03(X) = integral 0 to X { t**3 exp(t)/[exp(t)-1]**2 } dt

!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
!                       The recommended value is such that
!                             ATRAN(NTERMS) < EPS/100

!    XLOW2 - REAL (dp) - The value below which TRAN03 = 0.0 to machine
!                    precision. The recommended value is
!                          square root of (2*XMIN)

!    XLOW1 - REAL (dp) - The value below which TRAN03 = X**2/2 to
!                   machine precision. The recommended value is
!                             sqrt(8*EPSNEG)

!    XHIGH1 - REAL (dp) - The value above which the exponential series for
!                    large X contains only one term. The recommended value
!                    is        - ln(EPS).

!    XHIGH2 - REAL (dp) - The value above which
!                       TRAN03 = VALINF  -  X**3 exp(-X)
!                    The recommended value is 3/EPS

!    XHIGH3 - REAL (dp) - The value of ln(EPSNEG). Used to prevent overflow
!                    for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!     EXP, INT, LOG, SQRT


!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY ,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp, numjn = 3
REAL (dp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
              xlow1, xlow2
CHARACTER (LEN=14) :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6) :: fnname = 'TRAN03'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, onehun = 100.0_dp,  &
                         rnumjn = 3.0_dp, valinf = 7.2123414189575657124_dp

REAL (dp), PARAMETER  :: atran(0:19) = (/  &
    0.76201254324387200657_dp, -0.10567438770505853250_dp,  &
    0.1197780848196578097D-1, -0.121440152036983073D-2,  &
    0.11550997693928547D-3, -0.1058159921244229D-4, 0.94746633853018D-6,  &
   -0.8362212128581D-7, 0.731090992775D-8, -0.63505947788D-9,  &
    0.5491182819D-10, -0.473213954D-11, 0.40676948D-12, -0.3489706D-13,  &
    0.298923D-14, -0.25574D-15, 0.2186D-16, -0.187D-17, 0.16D-18, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 19, 0, -1
    IF (ABS(atran(nterms)) > t) EXIT
  END DO
  xlow1 = SQRT(eight*xk)
  xlow2 = SQRT(TINY(zero)/half)
ELSE
  xhigh1 = -LOG(2*EPSILON(zero))
  xhigh2 = rnumjn / xk
  xhigh3 = LOG(xk)
END IF

!   Code for x < =  4.0

IF (x <= four) THEN
  IF (x < xlow2) THEN
    fn_val = zero
  ELSE
    IF (x < xlow1) THEN
      fn_val = (x**(numjn-1)) / (rnumjn-one)
    ELSE
      t = (((x*x)/eight)-half) - half
      fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
    END IF
  END IF
ELSE
  
!  Code for x > 4.0
  
  IF (x > xhigh2) THEN
    sumexp = one
  ELSE
    IF (x <= xhigh1) THEN
      numexp = INT(xhigh1/x) + 1
      t = EXP(-x)
    ELSE
      numexp = 1
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, numjn
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = rnumjn * LOG(x) - x + LOG(sumexp)
  IF (t < xhigh3) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t)
  END IF
END IF
RETURN
END FUNCTION tran03



FUNCTION tran04(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates the transport integral of order 4, defined as

!      TRAN04(X) = integral 0 to X { t**4 exp(t)/[exp(t)-1]**2 } dt

!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
!                       The recommended value is such that
!                             ATRAN(NTERMS) < EPS/100

!    XLOW2 - REAL (dp) - The value below which TRAN04 = 0.0 to machine
!                   precision. The recommended value is
!                          cube root of (3*XMIN)

!    XLOW1 - REAL (dp) - The value below which TRAN04 = X**3/3 to
!                   machine precision. The recommended value is
!                             sqrt(8*EPSNEG)

!    XHIGH1 - REAL (dp) - The value above which the exponential series for
!                    large X contains only one term. The recommended value
!                    is        - ln(EPS).

!    XHIGH2 - REAL (dp) - The value above which
!                       TRAN04 = VALINF  -  X**4 exp(-X)
!                    The recommended value is 4/EPS

!    XHIGH3 - REAL (dp) - The value of ln(EPSNEG). Used to prevent overflow
!                    for large x.


!    For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!     EXP, INT, LOG, SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY ,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp, numjn = 4
REAL (dp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
              xlow1, xlow2
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'TRAN04'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, onehun = 100.0_dp,  &
                         rnumjn = 4.0_dp, valinf = 25.975757609067316596_dp

REAL (dp), PARAMETER  :: atran(0:19) = (/  &
    0.48075709946151105786_dp, -0.8175378810321083956D-1,  &
    0.1002700665975162973D-1, -0.105993393598201507D-2,  &
    0.10345062450304053D-3, -0.964427054858991D-5, 0.87455444085147D-6,  &
   -0.7793212079811D-7, 0.686498861410D-8, -0.59995710764D-9,  &
    0.5213662413D-10, -0.451183819D-11, 0.38921592D-12, -0.3349360D-13,  &
    0.287667D-14, -0.24668D-15, 0.2113D-16, -0.181D-17, 0.15D-18, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 19, 0, -1
    IF (ABS(atran(nterms)) > t) EXIT
  END DO
  xlow1 = SQRT(eight*xk)
  xk1 = rnumjn - one
  xlow2 = (xk1*TINY(zero)) ** (one/xk1)
ELSE
  xhigh1 = -LOG(2*EPSILON(zero))
  xhigh2 = rnumjn / xk
  xhigh3 = LOG(xk)
END IF

!   Code for x < =  4.0

IF (x <= four) THEN
  IF (x < xlow2) THEN
    fn_val = zero
  ELSE
    IF (x < xlow1) THEN
      fn_val = (x**(numjn-1)) / (rnumjn-one)
    ELSE
      t = (((x*x)/eight)-half) - half
      fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
    END IF
  END IF
ELSE
  
!  Code for x > 4.0
  
  IF (x > xhigh2) THEN
    sumexp = one
  ELSE
    IF (x <= xhigh1) THEN
      numexp = INT(xhigh1/x) + 1
      t = EXP(-x)
    ELSE
      numexp = 1
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, numjn
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = rnumjn * LOG(x) - x + LOG(sumexp)
  IF (t < xhigh3) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t)
  END IF
END IF
RETURN
END FUNCTION tran04



FUNCTION tran05(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates the transport integral of order n, defined as

!      TRAN05(X) = integral 0 to X { t**5 exp(t)/[exp(t)-1]**2 } dt

!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
!                       The recommended value is such that
!                             ATRAN(NTERMS) < EPS/100

!    XLOW2 - REAL (dp) - The value below which TRAN05 = 0.0 to machine
!                   precision. The recommended value is
!                          4th root of (4*XMIN)

!    XLOW1 - REAL (dp) - The value below which TRAN05 = X**4/4 to
!                   machine precision. The recommended value is
!                             sqrt(8*EPSNEG)

!    XHIGH1 - REAL (dp) - The value above which the exponential series for
!                    large X contains only one term. The recommended value
!                    is        - ln(EPS).

!    XHIGH2 - REAL (dp) - The value above which
!                       TRAN05 = VALINF  -  X**5 exp(-X)
!                    The recommended value is 5/EPS

!    XHIGH3 - REAL (dp) - The value of ln(EPSNEG). Used to prevent overflow
!                    for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!     EXP, INT, LOG, SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY ,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp, numjn = 5
REAL (dp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
              xlow1, xlow2
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'TRAN05'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, onehun = 100.0_dp,  &
                         rnumjn = 5.0_dp,  &
                         valinf = 0.12443133061720439116D3
REAL (dp), PARAMETER  :: atran(0:19) = (/  &
    0.34777777713391078928D0, -0.6645698897605042801D-1,  &
    0.861107265688330882D-2, -0.93966822237555384D-3,  &
    0.9363248060815134D-4, -0.885713193408328D-5, 0.81191498914503D-6,  &
   -0.7295765423277D-7, 0.646971455045D-8, -0.56849028255D-9,  &
    0.4962559787D-10, -0.431093996D-11, 0.37310094D-12, -0.3219769D-13,  &
    0.277220D-14, -0.23824D-15, 0.2044D-16, -0.175D-17, 0.15D-18, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 19, 0, -1
    IF (ABS(atran(nterms)) > t) EXIT
  END DO
  xlow1 = SQRT(eight*xk)
  xk1 = rnumjn - one
  xlow2 = (xk1*TINY(zero)) ** (one/xk1)
ELSE
  xhigh1 = -LOG(2*EPSILON(zero))
  xhigh2 = rnumjn / xk
  xhigh3 = LOG(xk)
END IF

!   Code for x < =  4.0

IF (x <= four) THEN
  IF (x < xlow2) THEN
    fn_val = zero
  ELSE
    IF (x < xlow1) THEN
      fn_val = (x**(numjn-1)) / (rnumjn-one)
    ELSE
      t = (((x*x)/eight)-half) - half
      fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
    END IF
  END IF
ELSE
  
!  Code for x > 4.0
  
  IF (x > xhigh2) THEN
    sumexp = one
  ELSE
    IF (x <= xhigh1) THEN
      numexp = INT(xhigh1/x) + 1
      t = EXP(-x)
    ELSE
      numexp = 1
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, numjn
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = rnumjn * LOG(x) - x + LOG(sumexp)
  IF (t < xhigh3) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t)
  END IF
END IF
RETURN
END FUNCTION tran05



FUNCTION tran06(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates the transport integral of order 6, defined as

!      TRAN06(X) = integral 0 to X { t**6 exp(t)/[exp(t)-1]**2 } dt

!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
!                       The recommended value is such that
!                             ATRAN(NTERMS) < EPS/100

!    XLOW2 - REAL (dp) - The value below which TRAN06 = 0.0 to machine
!                   precision. The recommended value is
!                          5th root of (5*XMIN)

!    XLOW1 - REAL (dp) - The value below which TRAN06 = X**5/5 to
!                   machine precision. The recommended value is
!                             sqrt(8*EPSNEG)

!    XHIGH1 - REAL (dp) - The value above which the exponential series for
!                    large X contains only one term. The recommended value
!                    is        - ln(EPS).

!    XHIGH2 - REAL (dp) - The value above which
!                       TRAN06 = VALINF  -  X**6 exp(-X)
!                    The recommended value is 6/EPS

!    XHIGH3 - REAL (dp) - The value of ln(EPSNEG). Used to prevent overflow
!                    for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!     EXP, INT, LOG, SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY ,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp, numjn = 6
REAL (dp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
              xlow1, xlow2
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'TRAN06'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, onehun = 100.0_dp,  &
                         rnumjn = 6.0_dp, valinf = 732.48700462880338059_dp

REAL (dp), PARAMETER  :: atran(0:19) = (/  &
    0.27127335397840008227D0, -0.5588610553191453393D-1,  &
    0.753919513290083056D-2, -0.84351138579211219D-3, 0.8549098079676702D-4,  &
   -0.818715493293098D-5, 0.75754240427986D-6, -0.6857306541831D-7,  &
    0.611700376031D-8, -0.54012707024D-9, 0.4734306435D-10, -0.412701055D-11, &
    0.35825603D-12, -0.3099752D-13, 0.267501D-14, -0.23036D-15,  &
    0.1980D-16, -0.170D-17, 0.15D-18, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 19, 0, -1
    IF (ABS(atran(nterms)) > t) EXIT
  END DO
  xlow1 = SQRT(eight*xk)
  xk1 = rnumjn - one
  xlow2 = (xk1*TINY(zero)) ** (one/xk1)
ELSE
  xhigh1 = -LOG(2*EPSILON(zero))
  xhigh2 = rnumjn / xk
  xhigh3 = LOG(xk)
END IF

!   Code for x < =  4 .0

IF (x <= four) THEN
  IF (x < xlow2) THEN
    fn_val = zero
  ELSE
    IF (x < xlow1) THEN
      fn_val = (x**(numjn-1)) / (rnumjn-one)
    ELSE
      t = (((x*x)/eight)-half) - half
      fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
    END IF
  END IF
ELSE
  
!  Code for x > 4 .0
  
  IF (x > xhigh2) THEN
    sumexp = one
  ELSE
    IF (x <= xhigh1) THEN
      numexp = INT(xhigh1/x) + 1
      t = EXP(-x)
    ELSE
      numexp = 1
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, numjn
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = rnumjn * LOG(x) - x + LOG(sumexp)
  IF (t < xhigh3) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t)
  END IF
END IF
RETURN
END FUNCTION tran06



FUNCTION tran07(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates the transport integral of order 7, defined as

!      TRAN07(X) = integral 0 to X { t**7 exp(t)/[exp(t)-1]**2 } dt

!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
!                       The recommended value is such that
!                             ATRAN(NTERMS) < EPS/100

!    XLOW2 - REAL (dp) - The value below which TRAN07 = 0.0 to machine
!                   precision. The recommended value is
!                          6th root of (6*XMIN)

!    XLOW1 - REAL (dp) - The value below which TRAN07 = X**6/6 to
!                   machine precision. The recommended value is
!                             sqrt(8*EPSNEG)

!    XHIGH1 - REAL (dp) - The value above which the exponential series for
!                    large X contains only one term. The recommended value
!                    is        - ln(EPS).

!    XHIGH2 - REAL (dp) - The value above which
!                       TRAN07 = VALINF  -  X**7 exp(-X)
!                    The recommended value is 7/EPS

!    XHIGH3 - REAL (dp) - The value of ln(EPSNEG). Used to prevent overflow
!                    for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!     EXP, INT, LOG, SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY ,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp, numjn = 7
REAL (dp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
              xlow1, xlow2
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'TRAN07'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, onehun = 100.0_dp,  &
                         rnumjn = 7.0_dp,  &
                         valinf = 0.50820803580048910473D4
REAL (dp), PARAMETER  :: atran(0:19) = (/  &
    0.22189250734010404423_dp, -0.4816751061177993694D-1,  &
    0.670092448103153629D-2, -0.76495183443082557D-3,  &
    0.7863485592348690D-4, -0.761025180887504D-5,  &
    0.70991696299917D-6, -0.6468025624903D-7, 0.580039233960D-8,  &
    -0.51443370149D-9, 0.4525944183D-10, -0.395800363D-11,  &
    0.34453785D-12, -0.2988292D-13, 0.258434D-14, -0.22297D-15,  &
    0.1920D-16, -0.165D-17, 0.14D-18, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 19, 0, -1
    IF (ABS(atran(nterms)) > t) EXIT
  END DO
  xlow1 = SQRT(eight*xk)
  xk1 = rnumjn - one
  xlow2 = (xk1*TINY(zero)) ** (one/xk1)
ELSE
  xhigh1 = -LOG(2*EPSILON(zero))
  xhigh2 = rnumjn / xk
  xhigh3 = LOG(xk)
END IF

!   Code for x <= 4.0

IF (x <= four) THEN
  IF (x < xlow2) THEN
    fn_val = zero
  ELSE
    IF (x < xlow1) THEN
      fn_val = (x**(numjn-1)) / (rnumjn-one)
    ELSE
      t = (((x*x)/eight)-half) - half
      fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
    END IF
  END IF
ELSE
  
!  Code for x > 4.0
  
  IF (x > xhigh2) THEN
    sumexp = one
  ELSE
    IF (x <= xhigh1) THEN
      numexp = INT(xhigh1/x) + 1
      t = EXP(-x)
    ELSE
      numexp = 1
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, numjn
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = rnumjn * LOG(x) - x + LOG(sumexp)
  IF (t < xhigh3) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t)
  END IF
END IF
RETURN
END FUNCTION tran07



FUNCTION tran08(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates the transport integral of order 8, defined as

!      TRAN08(X) = integral 0 to X { t**8 exp(t)/[exp(t)-1]**2 } dt

!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
!                       The recommended value is such that
!                             ATRAN(NTERMS) < EPS/100

!    XLOW2 - REAL (dp) - The value below which TRAN08 = 0.0 to machine
!                   precision. The recommended value is
!                          7th root of (7*XMIN)

!    XLOW1 - REAL (dp) - The value below which TRAN08 = X**7/7 to
!                   machine precision. The recommended value is
!                             sqrt(8*EPSNEG)

!    XHIGH1 - REAL (dp) - The value above which the exponential series for
!                    large X contains only one term. The recommended value
!                    is        - ln(EPS).

!    XHIGH2 - REAL (dp) - The value above which
!                       TRAN08 = VALINF  -  X**8 exp(-X)
!                    The recommended value is 8/EPS

!    XHIGH3 - REAL (dp) - The value of ln(EPSNEG). Used to prevent overflow
!                    for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT


!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!     EXP, INT, LOG, SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY ,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:  23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp, numjn = 8
REAL (dp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1,  &
              xlow1, xlow2
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'TRAN08'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, onehun = 100.0_dp,  &
                         rnumjn = 8.0_dp, valinf = 0.40484399001901115764D5
REAL (dp), PARAMETER  :: atran(0:19) = (/  &
    0.18750695774043719233_dp, -0.4229527646093673337D-1,  &
    0.602814856929065592D-2, -0.69961054811814776D-3,  &
    0.7278482421298789D-4, -0.710846250050067D-5, 0.66786706890115D-6,  &
   -0.6120157501844D-7, 0.551465264474D-8, -0.49105307052D-9,  &
    0.4335000869D-10, -0.380218700D-11, 0.33182369D-12, -0.2884512D-13,  &
    0.249958D-14, -0.21605D-15, 0.1863D-16, -0.160D-17, 0.14D-18, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 19, 0, -1
    IF (ABS(atran(nterms)) > t) EXIT
  END DO
  xlow1 = SQRT(eight*xk)
  xk1 = rnumjn - one
  xlow2 = (xk1*TINY(zero)) ** (one/xk1)
ELSE
  xhigh1 = -LOG(2*EPSILON(zero))
  xhigh2 = rnumjn / xk
  xhigh3 = LOG(xk)
END IF

!   Code for x < =  4.0

IF (x <= four) THEN
  IF (x < xlow2) THEN
    fn_val = zero
  ELSE
    IF (x < xlow1) THEN
      fn_val = (x**(numjn-1)) / (rnumjn-one)
    ELSE
      t = (((x*x)/eight)-half) - half
      fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
    END IF
  END IF
ELSE
  
!  Code for x > 4.0
  
  IF (x > xhigh2) THEN
    sumexp = one
  ELSE
    IF (x <= xhigh1) THEN
      numexp = INT(xhigh1/x) + 1
      t = EXP(-x)
    ELSE
      numexp = 1
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, numjn
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = rnumjn * LOG(x) - x + LOG(sumexp)
  IF (t < xhigh3) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t)
  END IF
END IF
RETURN
END FUNCTION tran08



FUNCTION tran09(xvalue) RESULT(fn_val)

!  DESCRIPTION:

!    This program calculates the transport integral of order 9, defined as

!      TRAN09(X) = integral 0 to X { t**9 exp(t)/[exp(t)-1]**2 } dt

!    The program uses a Chebyshev series, the coefficients of which are
!    given to an accuracy of 20 decimal places.


!  ERROR RETURNS:

!    If XVALUE < 0.0, an error message is printed, and the program
!    returns the value 0.0.


!  MACHINE-DEPENDENT CONSTANTS:

!    NTERMS - INTEGER - The number of terms of the array ATRAN to be used.
!                       The recommended value is such that
!                             ATRAN(NTERMS) < EPS/100

!    XLOW2 - REAL (dp) - The value below which TRAN09 = 0.0 to machine
!                   precision. The recommended value is
!                          8th root of (8*XMIN)

!    XLOW1 - REAL (dp) - The value below which TRAN09 = X**8/8 to
!                   machine precision. The recommended value is
!                             sqrt(8*EPSNEG)

!    XHIGH1 - REAL (dp) - The value above which the exponential series for
!                    large X contains only one term. The recommended value
!                    is        - ln(EPS).

!    XHIGH2 - REAL (dp) - The value above which
!                       TRAN09 = VALINF  -  X**9 exp(-X)
!                    The recommended value is 9/EPS

!    XHIGH3 - REAL (dp) - The value of ln(EPSNEG). Used to prevent overflow
!                    for large x.

!     For values of EPS, EPSNEG, and XMIN refer to the file MACHCON.TXT

!     The machine-dependent constants are computed internally by
!     using the D1MACH subroutine.


!  INTRINSIC FUNCTIONS USED:
!     EXP, INT, LOG, SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!  AUTHOR:

!     DR. ALLAN J. MACLEOD,
!     DEPT. OF MATHEMATICS AND STATISTICS,
!     UNIVERSITY OF PAISLEY ,
!     HIGH ST.,
!     PAISLEY,
!     SCOTLAND.
!     PA1 2BE.

!     (e-mail: macl_ms0@paisley.ac.uk )


!  LATEST REVISION:   23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: k1, k2, nterms, numexp, numjn = 9
REAL (dp)  :: rk, sumexp, sum2, t, x, xhigh1, xhigh2, xhigh3, xk, xk1, xlow1, &
              xlow2
CHARACTER (LEN=14)  :: errmsg = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'TRAN09'

REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                         four = 4.0_dp, eight = 8.0_dp, onehun = 100.0_dp, &
                         rnumjn = 9.0_dp, valinf = 0.36360880558872871397D6

REAL (dp), PARAMETER  :: atran(0:19) = (/  &
    0.16224049991949846835D0, -0.3768351452195937773D-1,  &
    0.547669715917719770D-2, -0.64443945009449521D-3, 0.6773645285280983D-4,  &
   -0.666813497582042D-5, 0.63047560019047D-6, -0.5807478663611D-7,  &
    0.525551305123D-8, -0.46968861761D-9, 0.4159395065D-10, -0.365808491D-11, &
    0.32000794D-12, -0.2787651D-13, 0.242017D-14, -0.20953D-15, 0.1810D-16,  &
   -0.156D-17, 0.13D-18, -0.1D-19 /)

!  Start execution

x = xvalue

!  Error test

IF (x < zero) THEN
  CALL errprn(fnname,errmsg)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

xk = EPSILON(zero)
t = xk / onehun
IF (x <= four) THEN
  DO  nterms = 19, 0, -1
    IF (ABS(atran(nterms)) > t) EXIT
  END DO
  xlow1 = SQRT(eight*xk)
  xk1 = rnumjn - one
  xlow2 = (xk1*TINY(zero)) ** (one/xk1)
ELSE
  xhigh1 = -LOG(2*EPSILON(zero))
  xhigh2 = rnumjn / xk
  xhigh3 = LOG(xk)
END IF

!   Code for x < =  4.0

IF (x <= four) THEN
  IF (x < xlow2) THEN
    fn_val = zero
  ELSE
    IF (x < xlow1) THEN
      fn_val = (x**(numjn-1)) / (rnumjn-one)
    ELSE
      t = (((x*x)/eight)-half) - half
      fn_val = (x**(numjn-1)) * cheval(nterms,atran,t)
    END IF
  END IF
ELSE
  
!  Code for x > 4.0
  
  IF (x > xhigh2) THEN
    sumexp = one
  ELSE
    IF (x <= xhigh1) THEN
      numexp = INT(xhigh1/x) + 1
      t = EXP(-x)
    ELSE
      numexp = 1
      t = one
    END IF
    rk = zero
    DO  k1 = 1, numexp
      rk = rk + one
    END DO
    sumexp = zero
    DO  k1 = 1, numexp
      sum2 = one
      xk = one / (rk*x)
      xk1 = one
      DO  k2 = 1, numjn
        sum2 = sum2 * xk1 * xk + one
        xk1 = xk1 + one
      END DO
      sumexp = sumexp * t + sum2
      rk = rk - one
    END DO
  END IF
  t = rnumjn * LOG(x) - x + LOG(sumexp)
  IF (t < xhigh3) THEN
    fn_val = valinf
  ELSE
    fn_val = valinf - EXP(t)
  END IF
END IF
RETURN
END FUNCTION tran09



FUNCTION y0int(xvalue) RESULT(fn_val)

!   DESCRIPTION:

!      This function calculates the integral of the Bessel
!      function Y0, defined as

!        Y0INT(x) = {integral 0 to x} Y0(t) dt

!      The code uses Chebyshev expansions whose coefficients are
!      given to 20 decimal places.


!   ERROR RETURNS:

!      If x < 0.0, the function is undefined. An error message
!      is printed and the function returns the value 0.0.

!      If the value of x is too large, it is impossible to
!      accurately compute the trigonometric functions used. An
!      error message is printed, and the function returns the
!      value 1.0.


!   MACHINE-DEPENDENT CONSTANTS:

!      NTERM1 - The no. of terms to be used from the array
!                ARJ01. The recommended value is such that
!                   ABS(ARJ01(NTERM1)) < EPS/100

!      NTERM2 - The no. of terms to be used from the array
!                ARY01. The recommended value is such that
!                   ABS(ARY01(NTERM2)) < EPS/100

!      NTERM3 - The no. of terms to be used from the array
!                ARY0A1. The recommended value is such that
!                   ABS(ARY0A1(NTERM3)) < EPS/100

!      NTERM4 - The no. of terms to be used from the array
!                ARY0A2. The recommended value is such that
!                   ABS(ARY0A2(NTERM4)) < EPS/100

!      XLOW - The value of x below which
!                  Y0INT(x) = x*(ln(x) - 0.11593)*2/pi
!             to machine-precision. The recommended value is
!                 sqrt(9*EPSNEG)

!      XHIGH - The value of x above which it is impossible
!              to calculate (x-pi/4) accurately. The recommended
!              value is      1/EPSNEG

!      For values of EPS and EPSNEG, refer to the file MACHCON.TXT

!      The machine-dependent constants are computed internally by
!      using the D1MACH subroutine.


!   INTRINSIC FUNCTIONS USED:
!      COS , LOG , SIN , SQRT

!   OTHER MISCFUN SUBROUTINES USED:
!          CHEVAL , ERRPRN, D1MACH


!   AUTHOR:
!          Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley,
!          Paisley,
!          SCOTLAND
!          PA1 2BE

!          (e-mail: macl_ms0@paisley.ac.uk)


!   LATEST REVISION:    23 January, 1996


REAL (dp), INTENT(IN)  :: xvalue
REAL (dp)              :: fn_val

INTEGER    :: nterm1, nterm2, nterm3, nterm4
REAL (dp)  :: pib41, t, temp, x, xhigh, xlow, xmpi4
CHARACTER (LEN=18)  :: ermsg2 = 'ARGUMENT TOO LARGE'
CHARACTER (LEN=14)  :: ermsg1 = 'ARGUMENT < 0.0'
CHARACTER (LEN= 6)  :: fnname = 'Y0INT '

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, nine = 9.0_dp,  &
                         sixten = 16.0_dp, onehun = 100.0_dp,   &
                         one28 = 128.0_dp, five12 = 512.0_dp,   &
                         rt2bpi = 0.79788456080286535588_dp,    &
                         pib411 = 201.0_dp, pib412 = 256.0_dp,  &
                         pib42 = 0.24191339744830961566D-3,     &
                         twobpi = 0.63661977236758134308_dp,    &
                         gal2m1 =  -1.11593151565841244881_dp,  &
                         gamln2 =  -0.11593151565841244881_dp

REAL (dp), PARAMETER  :: arj01(0:23) = (/  &
    0.38179279321690173518_dp, -0.21275636350505321870_dp,  &
    0.16754213407215794187_dp, -0.12853209772196398954_dp,  &
    0.10114405455778847013_dp, -0.9100795343201568859D-1,   &
    0.6401345264656873103D-1, -0.3066963029926754312D-1,    &
    0.1030836525325064201D-1, -0.255670650399956918D-2,     &
    0.48832755805798304D-3, -0.7424935126036077D-4,         &
    0.922260563730861D-5, -0.95522828307083D-6, 0.8388355845986D-7, &
   -0.633184488858D-8, 0.41560504221D-9, -0.2395529307D-10, 0.122286885D-11,  &
   -0.5569711D-13, 0.227820D-14, -0.8417D-16, 0.282D-17, -0.9D-19 /)

REAL (dp), PARAMETER  :: ary01(0:24) = (/  &
    0.54492696302724365490_dp, -0.14957323588684782157_dp,  &
    0.11085634486254842337_dp, -0.9495330018683777109D-1,   &
    0.6820817786991456963D-1, -0.10324653383368200408_dp,   &
    0.10625703287534425491_dp, -0.6258367679961681990D-1,   &
    0.2385645760338293285D-1, -0.644864913015404481D-2,     &
    0.131287082891002331D-2, -0.20988088174989640D-3,       &
    0.2716042484138347D-4, -0.291199114014694D-5, 0.26344333093795D-6, &
   -0.2041172069780D-7, 0.137124781317D-8, -0.8070680792D-10, 0.419883057D-11, &
   -0.19459104D-12, 0.808782D-14, -0.30329D-15, 0.1032D-16, -0.32D-18, &
    0.1D-19 /)

REAL (dp), PARAMETER  :: ary0a1(0:21) = (/  &
    1.24030133037518970827_dp, -0.478125353632280693D-2,  &
    0.6613148891706678D-4, -0.186042740486349D-5, &
    0.8362735565080D-7, -0.525857036731D-8, 0.42606363251D-9, &
   -0.4211761024D-10, 0.488946426D-11, -0.64834929D-12, 0.9617234D-13, &
   -0.1570367D-13, 0.278712D-14, -0.53222D-15, 0.10844D-15, -0.2342D-16,  &
    0.533D-17, -0.127D-17, 0.32D-18, -0.8D-19, 0.2D-19, -0.1D-19 /)

REAL (dp), PARAMETER  :: ary0a2(0:18) = (/  &
    1.99616096301341675339_dp, -0.190379819246668161D-2, &
    0.1539710927044226D-4, -0.31145088328103D-6, 0.1110850971321D-7,  &
   -0.58666787123D-9, 0.4139926949D-10, -0.365398763D-11, 0.38557568D-12,  &
   -0.4709800D-13, 0.650220D-14, -0.99624D-15, 0.16700D-15, -0.3028D-16,  &
    0.589D-17, -0.122D-17, 0.27D-18, -0.6D-19, 0.1D-19 /)

!   Start computation

x = xvalue

!   First error test

IF (x < zero) THEN
  CALL errprn(fnname,ermsg1)
  fn_val = zero
  RETURN
END IF

!   Compute the machine-dependent constants.

temp = EPSILON(zero)
xhigh = one / temp

!   Second error test

IF (x > xhigh) THEN
  CALL errprn(fnname,ermsg2)
  fn_val = zero
  RETURN
END IF

!   continue with machine constants

t = temp / onehun
IF (x <= sixten) THEN
  DO  nterm1 = 23, 0, -1
    IF (ABS(arj01(nterm1)) > t) EXIT
  END DO
  DO  nterm2 = 24, 0, -1
    IF (ABS(ary01(nterm2)) > t) EXIT
  END DO
  xlow = SQRT(nine*temp)
ELSE
  DO  nterm3 = 21, 0, -1
    IF (ABS(ary0a1(nterm3)) > t) EXIT
  END DO
  DO  nterm4 = 18, 0, -1
    IF (ABS(ary0a2(nterm4)) > t) EXIT
  END DO
END IF

!   Code for 0 <= x <= 16

IF (x <= sixten) THEN
  IF (x < xlow) THEN
    IF (x == zero) THEN
      fn_val = zero
    ELSE
      fn_val = (LOG(x)+gal2m1) * twobpi * x
    END IF
  ELSE
    t = x * x / one28 - one
    temp = (LOG(x)+gamln2) * cheval(nterm1,arj01,t)
    temp = temp - cheval(nterm2,ary01,t)
    fn_val = twobpi * x * temp
  END IF
ELSE
  
!   Code for x > 16
  
  t = five12 / (x*x) - one
  pib41 = pib411 / pib412
  xmpi4 = (x-pib41) - pib42
  temp = SIN(xmpi4) * cheval(nterm3,ary0a1,t) / x
  temp = temp + COS(xmpi4) * cheval(nterm4,ary0a2,t)
  fn_val = -rt2bpi * temp / SQRT(x)
END IF
RETURN
END FUNCTION y0int



FUNCTION cheval(n, a, t) RESULT(fn_val)

!   This function evaluates a Chebyshev series, using the Clenshaw method
!   with Reinsch modification, as analysed in the paper by Oliver.

!   INPUT PARAMETERS

!       N - INTEGER - The no. of terms in the sequence

!       A - REAL (dp) ARRAY, dimension 0 to N - The coefficients of
!           the Chebyshev series

!       T - REAL (dp) - The value at which the series is to be evaluated


!   REFERENCES

!        "An error analysis of the modified Clenshaw method for
!         evaluating Chebyshev and Fourier series" J. Oliver,
!         J.I.M.A., vol. 20, 1977, pp379-391


! MACHINE-DEPENDENT CONSTANTS: NONE


! INTRINSIC FUNCTIONS USED;
!    ABS


! AUTHOR:  Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley ,
!          High St.,
!          PAISLEY,
!          SCOTLAND


! LATEST MODIFICATION:   21 December , 1992


INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: a(0:n)
REAL (dp), INTENT(IN)  :: t
REAL (dp)              :: fn_val

INTEGER    :: i
REAL (dp)  :: d1, d2, tt, u0, u1, u2
REAL (dp), PARAMETER  :: zero = 0.0_dp, half = 0.5_dp, test = 0.6_dp,  &
                         two = 2.0_dp

u1 = zero

!   If ABS ( T )  < 0.6 use the standard Clenshaw method

IF (ABS(t) < test) THEN
  u0 = zero
  tt = t + t
  DO  i = n, 0, -1
    u2 = u1
    u1 = u0
    u0 = tt * u1 + a(i) - u2
  END DO
  fn_val = (u0-u2) / two
ELSE
  
!   If ABS ( T )  > =  0.6 use the Reinsch modification
  
  d1 = zero
  
!   T > =  0.6 code
  
  IF (t > zero) THEN
    tt = (t-half) - half
    tt = tt + tt
    DO  i = n, 0, -1
      d2 = d1
      u2 = u1
      d1 = tt * u2 + a(i) + d2
      u1 = d1 + u2
    END DO
    fn_val = (d1+d2) / two
  ELSE
    
!   T < =  -0.6 code
    
    tt = (t+half) + half
    tt = tt + tt
    DO  i = n, 0, -1
      d2 = d1
      u2 = u1
      d1 = tt * u2 + a(i) - d2
      u1 = d1 - u2
    END DO
    fn_val = (d1-d2) / two
  END IF
END IF
RETURN
END FUNCTION cheval



SUBROUTINE errprn(fnname,errmsg)

! DESCRIPTION:
!    This subroutine prints out an error message if
!    an error has occurred in one of the MISCFUN functions.


! INPUT PARAMETERS:

!    FNNAME - CHARACTER - The name of the function with the error.

!    ERRMSG - CHARACTER - The message to be printed out.


! MACHINE-DEPENDENT PARAMETER:

!    OUTSTR - INTEGER - The numerical value of the output stream to be used
!                       for printing the error message.
!                       The subroutine has the default value   OUTSTR = 6.


! AUTHOR:
!       DR. ALLAN J. MACLEOD,
!       DEPT. OF MATHEMATICS AND STATISTICS,
!       UNIVERSITY OF PAISLEY,
!       HIGH ST.,
!       PAISLEY
!       SCOTLAND
!       PA1 2BE

!       (e-mail: macl_ms0@paisley.ac.uk )


! LATEST REVISION:    2 JUNE, 1995


CHARACTER (LEN= 6), INTENT(IN)  :: fnname
CHARACTER (LEN= *), INTENT(IN)  :: errmsg

INTEGER  :: outstr = 6

WRITE (outstr,5000) fnname
WRITE (outstr,5100) errmsg
RETURN

5000 FORMAT (/t6, 'ERROR IN MISCFUN FUNCTION  ', a6)
5100 FORMAT (/t6, a50)

END SUBROUTINE errprn

END MODULE toms757
