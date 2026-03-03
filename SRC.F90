!      ALGORITHM 819, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 28,NO. 3, September, 2002, P.  325--336.

! Code converted using TO_F90 by Alan Miller
! Date: 2002-11-04  Time: 15:13:16


MODULE AiryFunction
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

PRIVATE
PUBLIC  :: aiz, biz

! COMMON /param1/ pi, pihal

REAL (dp), PARAMETER  :: pi = 3.1415926535897932385_dp,  &
                         pihal = 1.5707963267948966192_dp

! COMMON /param2/ pih3, pisr, a, alf

REAL (dp), PARAMETER  :: pih3 = 4.71238898038469_dp,   &
                         pisr = 1.77245385090552_dp
REAL (dp), SAVE       :: a, alf

! COMMON /param3/ thet, r, th15, s1, c1, r32

REAL (dp), SAVE  :: thet, r, th15, s1, c1, r32

! COMMON /param4/ facto, th025, s3, c3

REAL (dp), SAVE  :: facto, th025, s3, c3


CONTAINS


SUBROUTINE aiz(ifun, ifac, x0, y0, gair, gaii, ierro)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COMPUTATION OF THE AIRY FUNCTION AI(Z) OR ITS DERIVATIVE AI'(Z)
! THE CODE USES:
!      1. MACLAURIN SERIES FOR |Y| < 3 AND -2.6 < X <1.3 (Z=X + I*Y)
!      2. GAUSS-LAGUERRE QUADRATURE  FOR |Z| < 15 AND  WHEN
!         MACLAURIN SERIES ARE NOT USED.
!      3. ASYMPTOTIC EXPANSION FOR |Z| > 15.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  INPUTS:
!    IFUN:
!         * IFUN=1, THE CODE COMPUTES AI(Z)
!         * IFUN=2, THE CODE COMPUTES AI'(Z)
!    IFAC:
!         * IFAC=1, THE CODE COMPUTES  AI(Z) OR AI'(Z)
!         * IFAC=2, THE CODE COMPUTES NORMALIZED AI(Z) OR AI'(Z)
!    X0:   REAL PART OF THE ARGUMENT Z
!    Y0:   IMAGINARY PART OF THE ARGUMENT  Z

!  OUTPUTS:
!    GAIR: REAL PART OF AI(Z) OR AI'(Z)
!    GAII: IMAGINARY PART OF AI(Z) OR AI'(Z)

!    IERRO: ERROR FLAG
!          * IERRO=0, SUCCESSFUL COMPUTATION
!          * IERRO=1, COMPUTATION OUT OF RANGE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   ACCURACY:

!     1) SCALED AIRY FUNCTIONS:
!        RELATIVE ACCURACY BETTER THAN 10**(-13) EXCEPT CLOSE TO THE ZEROS,
!        WHERE 10**(-13) IS THE ABSOLUTE PRECISION.
!        GRADUAL LOSS OF PRECISION TAKES PLACE FOR |Z|>1000
!        (REACHING 10**(-8) ABSOLUTE ACCURACY FOR |Z| CLOSE TO 10**(6))
!        IN THE CASE OF PHASE(Z) CLOSE TO PI OR -PI.
!     2) UNSCALED AIRY FUNCTIONS:
!        THE FUNCTION OVERFLOWS/UNDERFLOWS FOR
!        3/2*|Z|**(3/2) > LOG(OVER).
!        FOR |Z| < 30:
!        A) RELATIVE ACCURACY FOR THE MODULUS (EXCEPT AT THE ZEROS)
!           BETTER THAN 10**(-13).
!        B) ABSOLUTE ACCURACY FOR MIN(R(Z),1/R(Z)) BETTER THAN 10**(-13),
!           WHERE R(Z)=REAL(AI)/IMAG(AI) OR R(Z)=REAL(AI')/IMAG(AI').
!        FOR |Z| > 30, GRADUAL LOSS OF PRECISION TAKES PLACE AS |Z| INCREASES.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     AUTHORS:
!        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN).
!                      E-MAIL: AMPARO.GIL@UAM.ES
!        JAVIER SEGURA (U. CARLOS III DE MADRID, MADRID, SPAIN).
!                      E-MAIL: JSEGURA@MATH.UC3M.ES
!        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
!                      E-MAIL: NICO.TEMME@CWI.NL
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    REFERENCES:
!         COMPUTING AIRY FUNCTIONS BY NUMERICAL QUADRATURE.
!         NUMERICAL ALGORITHMS (2002).
!         A. GIL, J. SEGURA, N.M. TEMME
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)     :: ifun
INTEGER, INTENT(IN)     :: ifac
REAL (dp), INTENT(IN)   :: x0
REAL (dp), INTENT(IN)   :: y0
REAL (dp), INTENT(OUT)  :: gair
REAL (dp), INTENT(OUT)  :: gaii
INTEGER, INTENT(OUT)    :: ierro

REAL (dp) :: over, under, dl1, dl2, cover
REAL (dp) :: f23, pi23, sqrt3, xa, ya, f23r, df1, df2, s11, c11, dex, dre,  &
             dima, gar, gai, c, s, u, v, v0, ar, ai, ar1, ai1, ro, coe1,   &
             coe2, rex, dfr, dfi, ar11, ai11
INTEGER   :: iexpf, iexpf2, n
! COMMON /param1/ pi, pihal
! COMMON /param2/ pih3, pisr, a, alf
! COMMON /param3/ thet, r, th15, s1, c1, r32
! COMMON /param4/ facto, th025, s3, c3
REAL (dp), PARAMETER  :: x(1:25) = (/   &
    .283891417994567679D-1, .170985378860034935D0,  &
    .435871678341770460D0, .823518257913030858D0,  &
    1.33452543254227372D0, 1.96968293206435071D0,  &
    2.72998134002859938D0, 3.61662161916100897D0,  &
    4.63102611052654146D0, 5.77485171830547694D0,  &
    7.05000568630218682D0, 8.45866437513237792D0,  &
    10.0032955242749393D0, 11.6866845947722423D0,  &
    13.5119659344693551D0, 15.4826596959377140D0,  &
    17.6027156808069112D0, 19.8765656022785451D0,  &
    22.3091856773962780D0, 24.9061720212974207D0,  &
    27.6738320739497190D0, 30.6192963295084111D0,  &
    33.7506560850239946D0, 37.0771349708391198D0,  &
    40.6093049694341322D0 /)
REAL (dp), PARAMETER  :: w(1:25) = (/   &
    .143720408803313866D0,  &
    .230407559241880881D0, .242253045521327626D0,  &
    .203636639103440807D0, .143760630622921410D0,  &
    .869128834706078120D-1, .454175001832915883D-1,  &
    .206118031206069497D-1, .814278821268606972D-2,  &
    .280266075663377634D-2, .840337441621719716D-3,  &
    .219303732907765020D-3, .497401659009257760D-4,  &
    .978508095920717661D-5, .166542824603725563D-5,  &
    .244502736801316287D-6, .308537034236207072D-7,  &
    .333296072940112245D-8, .306781892316295828D-9,  &
    .239331309885375719D-10, .157294707710054952D-11,  &
    .864936011664392267D-13, .394819815638647111D-14,  &
    .148271173082850884D-15, .453390377327054458D-17 /)
REAL (dp), PARAMETER  :: xd(1:25) = (/   &
    .435079659953445D-1, .205779160144678D0,  &
    .489916161318751D0, .896390483211727D0, 1.42582496737580D0,  &
    2.07903190767599D0, 2.85702335104978D0, 3.76102058198275D0,  &
    4.79246521225895D0, 5.95303247470003D0, 7.24464710774066D0,  &
    8.66950223642504D0, 10.2300817341775D0, 11.9291866622602D0,  &
    13.7699665302828D0, 15.7559563095946D0, 17.8911203751898D0,  &
    20.1799048700978D0, 22.6273004064466D0, 25.2389175786164D0,  &
    28.0210785229929D0, 30.9809287996116D0, 34.1265753192057D0,  &
    37.4672580871163D0, 41.0135664833476D0 /)
REAL (dp), PARAMETER  :: wd(1:25) = (/   &
    .576354557898966D-1,  &
    .139560003272262D0, .187792315011311D0, .187446935256946D0,  &
    .150716717316301D0, .101069904453380D0, .575274105486025D-1,  &
    .280625783448681D-1, .117972164134041D-1, .428701743297432D-2,  &
    .134857915232883D-2, .367337337105948D-3, .865882267841931D-4,  &
    .176391622890609D-4, .309929190938078D-5, .468479653648208D-6,  &
    .607273267228907D-7, .672514812555074D-8, .633469931761606D-9,  &
    .504938861248542D-10, .338602527895834D-11,  &
    .189738532450555D-12, .881618802142698D-14,  &
    .336676636121976D-15, .104594827170761D-16 /)

! CONSTANTS
! pi = 3.1415926535897932385D0
! pihal = 1.5707963267948966192D0
! pih3 = 4.71238898038469D0
f23 = .6666666666666667D0
pi23 = 2.09439510239320D0
! pisr = 1.77245385090552D0
sqrt3 = 1.7320508075688772935D0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ya = y0
xa = x0
ierro = 0
iexpf = 0
iexpf2 = 0
IF (ya < 0.d0) ya = -ya
r = SQRT(xa*xa + ya*ya)
r32 = r * SQRT(r)
thet = phase(xa,ya)
cover = 2.d0 / 3.d0 * r32 * ABS(COS(1.5D0*thet))

! MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
over = HUGE(pi) * 1.d-4

! MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under = TINY(pi) * 1.d+4
dl1 = LOG(over)
dl2 = -LOG(under)
IF (dl1 > dl2) over = 1 / under
IF (ifac == 1) THEN
  IF (cover >= LOG(over)) THEN
! OVERFLOW/UNDERFLOW PROBLEMS.
!   CALCULATION ABORTED
    ierro = 1
    gair = 0
    gaii = 0
  END IF
  IF (cover >= (LOG(over)*0.2)) iexpf2 = 1
ELSE
  IF (cover >= (LOG(over)*0.2)) iexpf = 1
END IF
IF (ierro == 0) THEN
  IF (ifun == 1) THEN
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CC CALCULATION OF AI(Z) CCCCCCCCCCCCCCC
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CCCCCCCCCCC  SERIES, INTEGRALS OR EXPANSIONS
    IF ((ya < 3.d0) .AND. (xa < 1.3D0) .AND. (xa > -2.6D0)) THEN
! CC  SERIES
      CALL serai(xa, ya, gar, gai)
      IF (ifac == 2) THEN
        thet = phase(xa,ya)
        th15 = 1.5D0 * thet
        s1 = SIN(th15)
        c1 = COS(th15)
        f23r = f23 * r32
        df1 = f23r * c1
        df2 = f23r * s1
        s11 = SIN(df2)
        c11 = COS(df2)
        dex = EXP(df1)
        dre = dex * c11
        dima = dex * s11
        gair = dre * gar - dima * gai
        gaii = dre * gai + dima * gar
      ELSE
        gair = gar
        gaii = gai
        IF (y0 == 0.) gaii = 0.d0
      END IF
    ELSE
      IF (r > 15.d0) THEN
!CC ASYMPTOTIC EXPANSIONS CCC
        thet = phase(xa, ya)
        facto = 0.5D0 / pisr * r ** (-0.25D0)
        IF (thet > pi23) THEN
!CCCCCCCCCC CONNECTION FORMULAE CCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     N= 1: TRANSFORM Z TO W= U+IV=Z EXP( 2 PI I/3)
          n = 1
          c = -0.5D0
          s = n * 0.5 * sqrt3
          u = xa * c - ya * s
          v = xa * s + ya * c
          v0 = v
          IF (v < 0.d0) v = -v
          thet = phase(u, v)
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          th025 = thet * 0.25D0
          s3 = SIN(th025)
          c3 = COS(th025)
          CALL expai(ar1, ai1)
          IF (v0 < 0.d0) ai1 = -ai1
          ar = -(c*ar1 - s*ai1)
          ai = -(s*ar1 + c*ai1)
          IF (iexpf == 0) THEN
            IF (iexpf2 == 0) THEN
!     N=-1: TRANSFORM Z TO W= U+IV=Z EXP(-2 PI I/3)
              n = -1
              c = -0.5D0
              s = n * 0.5 * sqrt3
              u = xa * c - ya * s
              v = xa * s + ya * c
              v0 = v
              IF (v < 0.d0) v = -v
              thet = phase(u,v)
              th15 = 1.5D0 * thet
              s1 = SIN(th15)
              c1 = COS(th15)
              th025 = thet * 0.25D0
              s3 = SIN(th025)
              c3 = COS(th025)
              CALL expai(ar1, ai1)
              IF (v0 < 0.d0) ai1 = -ai1
              thet = phase(xa,ya)
              th15 = 1.5D0 * thet
              s1 = SIN(th15)
              c1 = COS(th15)
              ro = 1.333333333333333D0 * r32
              coe1 = ro * c1
              coe2 = ro * s1
              rex = EXP(coe1)
              dfr = rex * COS(coe2)
              dfi = rex * SIN(coe2)
              ar11 = dfr * ar1 - dfi * ai1
              ai11 = dfr * ai1 + dfi * ar1
              gair = ar - (c*ar11 - s*ai11)
              gaii = ai - (s*ar11 + c*ai11)
            ELSE
              thet = phase(xa,ya)
              th15 = 1.5D0 * thet
              s1 = SIN(th15)
              c1 = COS(th15)
              gair = ar
              gaii = ai
            END IF
          ELSE
            gair = ar
            gaii = ai
          END IF
        ELSE
!CCCCCC  ASYMPTOTIC EXPANSION CCCCCCCCCCCCCCC
          thet = phase(xa, ya)
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          th025 = thet * 0.25D0
          s3 = SIN(th025)
          c3 = COS(th025)
          CALL expai(gair, gaii)
        END IF
      ELSE
!CC INTEGRALS
        a = 0.1666666666666667D0
        alf = -a
        facto = 0.280514117723058D0 * r ** (-0.25D0)
        thet = phase(xa, ya)
        IF (thet <= pihal) THEN
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          th025 = thet * 0.25D0
          s3 = SIN(th025)
          c3 = COS(th025)
          CALL airy1(x, w, gair, gaii)
        END IF
        IF ((thet > pihal) .AND. (thet <= pi23)) THEN
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          th025 = thet * 0.25D0
          s3 = SIN(th025)
          c3 = COS(th025)
          CALL airy2(x, w, gair, gaii)
        END IF
        IF (thet > pi23) THEN
          n = 1
          c = -0.5D0
          s = n * 0.5 * sqrt3
          u = xa * c - ya * s
          v = xa * s + ya * c
          v0 = v
          IF (v < 0.d0) v = -v
          thet = phase(u,v)
          IF (thet <= pihal) THEN
            th15 = 1.5D0 * thet
            s1 = SIN(th15)
            c1 = COS(th15)
            th025 = thet * 0.25D0
            s3 = SIN(th025)
            c3 = COS(th025)
            CALL airy1(x,w,ar1,ai1)
          END IF
          IF ((thet > pihal) .AND. (thet <= pi23)) THEN
            th15 = 1.5D0 * thet
            s1 = SIN(th15)
            c1 = COS(th15)
            th025 = thet * 0.25D0
            s3 = SIN(th025)
            c3 = COS(th025)
            CALL airy2(x,w,ar1,ai1)
          END IF
          IF (v0 < 0.d0) ai1 = -ai1
          ar = -(c*ar1-s*ai1)
          ai = -(s*ar1+c*ai1)
          n = -1
          c = -0.5D0
          s = n * 0.5 * sqrt3
          u = xa * c - ya * s
          v = xa * s + ya * c
          v0 = v
          IF (v < 0.d0) v = -v
          thet = phase(u,v)
          IF (thet <= pihal) THEN
            th15 = 1.5D0 * thet
            s1 = SIN(th15)
            c1 = COS(th15)
            th025 = thet * 0.25D0
            s3 = SIN(th025)
            c3 = COS(th025)
            CALL airy1(x,w,ar1,ai1)
          END IF
          IF ((thet > pihal) .AND. (thet <= pi23)) THEN
            th15 = 1.5D0 * thet
            s1 = SIN(th15)
            c1 = COS(th15)
            th025 = thet * 0.25D0
            s3 = SIN(th025)
            c3 = COS(th025)
            CALL airy2(x,w,ar1,ai1)
          END IF
          IF (v0 < 0.d0) ai1 = -ai1
          thet = phase(xa,ya)
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          ro = 1.333333333333333D0 * r32
          coe1 = ro * c1
          coe2 = ro * s1
          rex = EXP(coe1)
          dfr = rex * COS(coe2)
          dfi = rex * SIN(coe2)
          ar11 = dfr * ar1 - dfi * ai1
          ai11 = dfr * ai1 + dfi * ar1
          gair = ar - (c*ar11 - s*ai11)
          gaii = ai - (s*ar11 + c*ai11)
        END IF
      END IF

      IF (ifac == 1) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF THE UNSCALED AI(Z) CCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        f23r = f23 * r32
        df1 = f23r * c1
        df2 = f23r * s1
        s11 = SIN(df2)
        c11 = COS(df2)
        dex = EXP(-df1)
        dre = dex * c11
        dima = -dex * s11
        gar = dre * gair - dima * gaii
        gai = dre * gaii + dima * gair
        gair = gar
        gaii = gai
        IF (y0 == 0.) gaii = 0.d0
      END IF
    END IF
  ELSE
!CCC CALCULATION OF AI(Z) CCCCCCCCCCC
    alf = 0.1666666666666667D0
    facto = -0.270898621247918D0 * r ** 0.25D0
!CCCCCCCCCCCCCC SERIES OR INTEGRALS CCCCCCCCCCCCCCCCCCCCCCCCCC
    IF (ya < 3.d0 .AND. xa < 1.3D0 .AND. xa > -2.6D0) THEN
!CC SERIES
      CALL seraid(xa, ya, gar, gai)
      IF (ifac == 2) THEN
        thet = phase(xa,ya)
        th15 = 1.5D0 * thet
        s1 = SIN(th15)
        c1 = COS(th15)
        f23r = f23 * r32
        df1 = f23r * c1
        df2 = f23r * s1
        s11 = SIN(df2)
        c11 = COS(df2)
        dex = EXP(df1)
        dre = dex * c11
        dima = dex * s11
        gair = dre * gar - dima * gai
        gaii = dre * gai + dima * gar
      ELSE
        gair = gar
        gaii = gai
        IF (y0 == 0.) gaii = 0.d0
      END IF
    ELSE
      IF (r > 15.d0) THEN
!CC  ASYMPTOTIC EXPANSIONS CCCCCCCCCCCCC
        thet = phase(xa, ya)
        facto = 0.5D0 / pisr * r ** 0.25D0
        IF (thet > pi23) THEN
!CCCCCC CONNECTION FORMULAE CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC     N= 1: TRANSFORM Z TO W= U+IV=Z EXP( 2 PI I/3)
          n = 1
          c = -0.5D0
          s = n * 0.5 * sqrt3
          u = xa * c - ya * s
          v = xa * s + ya * c
          v0 = v
          IF (v < 0.d0) v = -v
          thet = phase(u,v)
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          th025 = thet * 0.25D0
          s3 = SIN(th025)
          c3 = COS(th025)
          CALL expaid(ar1,ai1)
          IF (v0 < 0.d0) ai1 = -ai1
          ar = -(c*ar1+s*ai1)
          ai = -(-s*ar1+c*ai1)
          IF (iexpf == 0) THEN
            IF (iexpf2 == 0) THEN
!CC     N=-1: TRANSFORM Z TO W= U+IV=Z EXP(-2 PI I/3)
              n = -1
              c = -0.5D0
              s = n * 0.5 * sqrt3
              u = xa * c - ya * s
              v = xa * s + ya * c
              v0 = v
              IF (v < 0.d0) v = -v
              thet = phase(u,v)
              th15 = 1.5D0 * thet
              s1 = SIN(th15)
              c1 = COS(th15)
              th025 = thet * 0.25D0
              s3 = SIN(th025)
              c3 = COS(th025)
              CALL expaid(ar1,ai1)
              IF (v0 < 0.d0) ai1 = -ai1
              thet = phase(xa,ya)
              th15 = 1.5D0 * thet
              s1 = SIN(th15)
              c1 = COS(th15)
              ro = 1.333333333333333D0 * r32
              coe1 = ro * c1
              coe2 = ro * s1
              rex = EXP(coe1)
              dfr = rex * COS(coe2)
              dfi = rex * SIN(coe2)
              ar11 = dfr * ar1 - dfi * ai1
              ai11 = dfr * ai1 + dfi * ar1
              gair = ar - (c*ar11+s*ai11)
              gaii = ai - (-s*ar11+c*ai11)
            ELSE
              thet = phase(xa,ya)
              th15 = 1.5D0 * thet
              s1 = SIN(th15)
              c1 = COS(th15)
              gair = ar
              gaii = ai
            END IF
          ELSE
            gair = ar
            gaii = ai
          END IF
        ELSE
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          th025 = thet * 0.25D0
          s3 = SIN(th025)
          c3 = COS(th025)
          CALL expaid(gair, gaii)
        END IF
      ELSE
!CC INTEGRALS CCCCCCCCCCCCCCCC
        thet = phase(xa,ya)
        IF (thet <= pihal) THEN
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          th025 = thet * 0.25D0
          s3 = SIN(th025)
          c3 = COS(th025)
          CALL airy1d(xd, wd, gair, gaii)
        END IF
        IF (thet > pihal .AND. thet <= pi23) THEN
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          th025 = thet * 0.25D0
          s3 = SIN(th025)
          c3 = COS(th025)
          CALL airy2d(xd, wd, gair, gaii)
        END IF
        IF (thet > pi23) THEN
          n = 1
          c = -0.5D0
          s = n * 0.5 * sqrt3
          u = xa * c - ya * s
          v = xa * s + ya * c
          v0 = v
          IF (v < 0.d0) v = -v
          thet = phase(u,v)
          IF (thet <= pihal) THEN
            th15 = 1.5D0 * thet
            s1 = SIN(th15)
            c1 = COS(th15)
            th025 = thet * 0.25D0
            s3 = SIN(th025)
            c3 = COS(th025)
            CALL airy1d(xd, wd, ar1, ai1)
          END IF
          IF (thet > pihal .AND. thet <= pi23) THEN
            th15 = 1.5D0 * thet
            s1 = SIN(th15)
            c1 = COS(th15)
            th025 = thet * 0.25D0
            s3 = SIN(th025)
            c3 = COS(th025)
            CALL airy2d(xd, wd, ar1, ai1)
          END IF
          IF (v0 < 0.d0) ai1 = -ai1
          ar = -(c*ar1 + s*ai1)
          ai = -(-s*ar1 + c*ai1)
          n = -1
          c = -0.5D0
          s = n * 0.5 * sqrt3
          u = xa * c - ya * s
          v = xa * s + ya * c
          v0 = v
          IF (v < 0.d0) v = -v
          thet = phase(u,v)
          IF (thet <= pihal) THEN
            th15 = 1.5D0 * thet
            s1 = SIN(th15)
            c1 = COS(th15)
            th025 = thet * 0.25D0
            s3 = SIN(th025)
            c3 = COS(th025)
            CALL airy1d(xd, wd, ar1, ai1)
          END IF
          IF (thet > pihal .AND. thet <= pi23) THEN
            th15 = 1.5D0 * thet
            s1 = SIN(th15)
            c1 = COS(th15)
            th025 = thet * 0.25D0
            s3 = SIN(th025)
            c3 = COS(th025)
            CALL airy2d(xd, wd, ar1, ai1)
          END IF
          IF (v0 < 0.d0) ai1 = -ai1
          thet = phase(xa,ya)
          th15 = 1.5D0 * thet
          s1 = SIN(th15)
          c1 = COS(th15)
          ro = 1.333333333333333D0 * r32
          coe1 = ro * c1
          coe2 = ro * s1
          rex = EXP(coe1)
          dfr = rex * COS(coe2)
          dfi = rex * SIN(coe2)
          ar11 = dfr * ar1 - dfi * ai1
          ai11 = dfr * ai1 + dfi * ar1
          gair = ar - (c*ar11 + s*ai11)
          gaii = ai - (-s*ar11 + c*ai11)
        END IF
      END IF
      IF (ifac == 1) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC CALCULATION OF THE UNSCALED AI'(z) CCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        f23r = f23 * r32
        df1 = f23r * c1
        df2 = f23r * s1
        s11 = SIN(df2)
        c11 = COS(df2)
        dex = EXP(-df1)
        dre = dex * c11
        dima = -dex * s11
        gar = dre * gair - dima * gaii
        gai = dre * gaii + dima * gair
        gair = gar
        gaii = gai
        IF (y0 == 0) gaii = 0.d0
      END IF
    END IF
  END IF
END IF
IF (y0 < 0.d0) gaii = -gaii
RETURN
END SUBROUTINE aiz


SUBROUTINE airy1(x, w, gair, gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR    C
!CC              0 <= PHASE(Z) <= PI/2                           C
!CC                                                              C
!CC INPUTS:                                                      C
!CC      X,W,      NODES AND WEIGHTS FOR THE GAUSSIAN QUADRATURE C
!CC OUTPUTS:                                                     C
!CC      GAIR, GAII,  REAL AND IMAGINARY PARTS OF AI(Z)          C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(IN)   :: x(25)
REAL (dp), INTENT(IN)   :: w(25)
REAL (dp), INTENT(OUT)  :: gair
REAL (dp), INTENT(OUT)  :: gaii

REAL (dp) :: sumar, sumai, df1, df1c1, phi, phi6, s2,  &
             c2, dmodu, dmodu2, funr, funi, fac1, fac2
INTEGER   :: i

! COMMON /param2/ pih3, pisr, a, alf
! COMMON /param3/ thet, r, th15, s1, c1, r32
! COMMON /param4/ facto, th025, s3, c3

sumar = 0.d0
sumai = 0.d0
DO  i = 1, 25
  df1 = 1.5D0 * x(i) / r32
  df1c1 = df1 * c1
  phi = phase(2.d0+df1c1,df1*s1)
  phi6 = phi / 6.d0
  s2 = SIN(phi6)
  c2 = COS(phi6)
  dmodu = SQRT(4.d0+df1*df1+4.d0*df1c1)
  dmodu2 = dmodu ** alf
  funr = dmodu2 * c2
  funi = dmodu2 * s2
  sumar = sumar + w(i) * funr
  sumai = sumai + w(i) * funi
END DO
fac1 = facto * c3
fac2 = facto * s3
gair = fac1 * sumar + fac2 * sumai
gaii = fac1 * sumai - fac2 * sumar
RETURN
END SUBROUTINE airy1


SUBROUTINE airy2(x, w, gair, gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR      C
!CC              PI/2 < PHASE(Z) <= 2PI/3                          C
!CC                                                                C
!CC INPUTS:                                                        C
!CC      X,W,        NODES AND WEIGHTS FOR THE GAUSSIAN QUADRATURE C
!CC OUTPUTS:                                                       C
!CC      GAIR, GAII, REAL AND IMAGINARY PARTS OF AI(Z)             C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(IN)   :: x(25)
REAL (dp), INTENT(IN)   :: w(25)
REAL (dp), INTENT(OUT)  :: gair
REAL (dp), INTENT(OUT)  :: gaii

REAL (dp) :: sumar, sumai, df1, df1c1, phi, phi6, s2,  &
             c2, dmodu, dmodu2, funr, funi, fac1, fac2
REAL (dp) :: sqr2, sqr2i, tau, tgtau, b, ang, ctau, cfac, ct,  &
             st, sumr, sumi, ttau, beta
INTEGER   :: i

! COMMON /param2/ pih3, pisr, a, alf
! COMMON /param3/ thet, r, th15, s1, c1, r32
! COMMON /param4/ facto, th025, s3, c3

sqr2 = 1.41421356237310D0
sqr2i = 0.707106781186548D0
tau = th15 - pih3 * 0.5D0
tgtau = TAN(tau)
b = 5.d0 * a
ang = tau * b
ctau = COS(tau)
cfac = ctau ** (-b)
ct = COS(ang)
st = SIN(ang)
sumr = 0.d0
sumi = 0.d0
DO  i = 1, 25
  df1 = 3.d0 * x(i) / (ctau*r32)
  df1c1 = df1 * sqr2i * 0.5D0
  phi = phase(2.d0-df1c1, df1c1)
  phi6 = phi / 6.d0
  ttau = x(i) * tgtau
  beta = phi6 - ttau
  s2 = SIN(beta)
  c2 = COS(beta)
  dmodu = SQRT(4.d0 + df1*df1*0.25D0 - sqr2*df1)
  dmodu2 = dmodu ** alf
  funr = dmodu2 * c2
  funi = dmodu2 * s2
  sumr = sumr + w(i) * funr
  sumi = sumi + w(i) * funi
END DO
sumar = cfac * (ct*sumr - st*sumi)
sumai = cfac * (ct*sumi + st*sumr)
fac1 = facto * c3
fac2 = facto * s3
gair = fac1 * sumar + fac2 * sumai
gaii = fac1 * sumai - fac2 * sumar
RETURN
END SUBROUTINE airy2


SUBROUTINE airy1d(x, w, gair, gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI'(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR    C
!CC              0 <= PHASE(Z) <= PI/2                            C
!CC                                                               C
!CC INPUTS:                                                       C
!CC       X,W,      NODES AND WEIGHTS FOR THE GAUSSIAN QUADRATURE C
!CC OUTPUTS:                                                      C
!CC       GAIR,GAII, REAL AND IMAGINARY PARTS OF AI'(Z)           C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(IN)   :: x(25)
REAL (dp), INTENT(IN)   :: w(25)
REAL (dp), INTENT(OUT)  :: gair
REAL (dp), INTENT(OUT)  :: gaii

REAL (dp) :: sumar, sumai, df1, df1c1, phi, phi6, s2,  &
             c2, dmodu, dmodu2, funr, funi, fac1, fac2
INTEGER :: i

! COMMON /param2/ pih3, pisr, a, alf
! COMMON /param3/ thet, r, th15, s1, c1, r32
! COMMON /param4/ facto, th025, s3, c3

sumar = 0.d0
sumai = 0.d0
DO  i = 1, 25
  df1 = 1.5D0 * x(i) / r32
  df1c1 = df1 * c1
  phi = phase(2.d0+df1c1, df1*s1)
  phi6 = -phi * alf
  s2 = SIN(phi6)
  c2 = COS(phi6)
  dmodu = SQRT(4.d0 + df1*df1 + 4.d0*df1c1)
  dmodu2 = dmodu ** alf
  funr = dmodu2 * c2
  funi = dmodu2 * s2
  sumar = sumar + w(i) * funr
  sumai = sumai + w(i) * funi
END DO
fac1 = facto * c3
fac2 = facto * s3
gair = fac1 * sumar - fac2 * sumai
gaii = fac1 * sumai + fac2 * sumar
RETURN
END SUBROUTINE airy1d


SUBROUTINE airy2d(x, w, gair, gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC COMPUTES AI'(Z) BY GAUSS-LAGUERRE QUADRATURE IN THE SECTOR C
!CC              PI/2 < PHASE(Z) <= 3PI/2                      C
!CC                                                            C
!CC INPUTS:                                                    C
!CC      X,W,   NODES AND WEIGHTS FOR THE GAUSSIAN QUADRATURE  C
!CC OUTPUTS:                                                   C
!CC      GAIR,GAII, REAL AND IMAGINARY PARTS OF AI'(Z)         C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(IN)   :: x(25)
REAL (dp), INTENT(IN)   :: w(25)
REAL (dp), INTENT(OUT)  :: gair
REAL (dp), INTENT(OUT)  :: gaii

REAL (dp) :: sumar, sumai, df1, df1c1, phi, phi6, s2,  &
             c2, dmodu, dmodu2, funr, funi, fac1, fac2
REAL (dp) :: sqr2, sqr2i, tau, tgtau, b, ang, ctau, cfac, ct,  &
             st, sumr, sumi, ttau, beta
INTEGER   :: i

! COMMON /param2/ pih3, pisr, a, alf
! COMMON /param3/ thet, r, th15, s1, c1, r32
! COMMON /param4/ facto, th025, s3, c3

sqr2 = 1.41421356237310D0
sqr2i = 0.707106781186548D0
tau = th15 - pih3 * 0.5D0
tgtau = TAN(tau)
b = 7.d0 * alf
ang = tau * b
ctau = COS(tau)
cfac = ctau ** (-b)
ct = COS(ang)
st = SIN(ang)
sumr = 0.d0
sumi = 0.d0
DO  i = 1, 25
  df1 = 3.d0 * x(i) / (ctau*r32)
  df1c1 = df1 * sqr2i * 0.5D0
  phi = phase(2.d0-df1c1,df1c1)
  phi6 = -phi / 6.d0
  ttau = x(i) * tgtau
  beta = phi6 - ttau
  s2 = SIN(beta)
  c2 = COS(beta)
  dmodu = SQRT(4.d0 + df1*df1*0.25D0 - sqr2*df1)
  dmodu2 = dmodu ** alf
  funr = dmodu2 * c2
  funi = dmodu2 * s2
  sumr = sumr + w(i) * funr
  sumi = sumi + w(i) * funi
END DO
sumar = cfac * (ct*sumr - st*sumi)
sumai = cfac * (ct*sumi + st*sumr)
fac1 = facto * c3
fac2 = facto * s3
gair = fac1 * sumar - fac2 * sumai
gaii = fac1 * sumai + fac2 * sumar
RETURN
END SUBROUTINE airy2d


FUNCTION phase(x, y) RESULT(fn_val)

REAL (dp), INTENT(IN)  :: x
REAL (dp), INTENT(IN)  :: y
REAL (dp)              :: fn_val

REAL (dp) :: ay, p

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  COMPUTES THE PHASE OF Z = X + IY, IN (-PI,PI]
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COMMON /param1/ pi, pihal

IF (x == 0.0_dp .AND. y == 0.0_dp) THEN
  p = 0.d0
ELSE
  ay = ABS(y)
  IF (x >= ay) THEN
    p = ATAN(ay/x)
  ELSE IF (x+ay >= 0.d0) THEN
    p = pihal - ATAN(x/ay)
  ELSE
    p = pi + ATAN(ay/x)
  END IF
  IF (y < 0.d0) p = -p
END IF

fn_val = p
RETURN
END FUNCTION phase


SUBROUTINE fgp(x, y, eps, fr, fi, gr, gi)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    COMPUTES THE FUNCTIONS F AND G FOR THE SERIES OF AI'(Z).
!    THIS ROUTINE IS CALLED BY SERAID.
!
!    INPUTS:
!         X,Y,  REAL AND IMAGINARY PARTS OF Z
!         EPS,  PRECISION FOR THE COMPUTATION OF THE SERIES
!    OUTPUTS:
!         FR,FI, REAL AND IMAGINARY PARTS OF F
!         GR,GI, REAL AND IMAGINARY PARTS OF G
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: y
REAL (dp), INTENT(IN)   :: eps
REAL (dp), INTENT(OUT)  :: fr
REAL (dp), INTENT(OUT)  :: fi
REAL (dp), INTENT(OUT)  :: gr
REAL (dp), INTENT(OUT)  :: gi

INTEGER   :: a, b, k3
REAL (dp) :: x2, y2, u, v, p, q, cr, ci, dr, di

x2 = x * x
y2 = y * y
k3 = 0
u = x * (x2-3*y2)
v = y * (3*x2-y2)
cr = 0.5D0
ci = 0.d0
dr = 1.d0
di = 0.d0
fr = 0.5D0
fi = 0.d0
gr = 1.d0
gi = 0.d0
10 a = (k3+5) * (k3+3)
b = (k3+1) * (k3+3)
p = (u*cr-v*ci) / a
q = (v*cr+u*ci) / a
cr = p
ci = q
p = (u*dr-v*di) / b
q = (v*dr+u*di) / b
dr = p
di = q
fr = fr + cr
fi = fi + ci
gr = gr + dr
gi = gi + di
k3 = k3 + 3
IF (ABS(cr) + ABS(dr) + ABS(ci) + ABS(di) >= eps) GO TO 10
u = x2 - y2
v = 2.d0 * x * y
p = u * fr - v * fi
q = u * fi + v * fr
fr = p
fi = q
RETURN
END SUBROUTINE fgp


SUBROUTINE fg(x, y, eps, fr, fi, gr, gi)
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    COMPUTES THE FUNCTIONS F AND G IN EXPRESSION
!    10.4.2 OF ABRAMOWITZ & STEGUN FOR THE SERIES OF AI(Z).
!    THIS ROUTINE IS CALLED BY SERAI.
!
!    INPUTS:
!          X,Y,  REAL AND IMAGINARY PARTS OF Z
!          EPS,  PRECISION FOR THE COMPUTATION OF THE SERIES.
!    OUTPUTS:
!          FR,FI, REAL AND IMAGINARY PARTS OF F
!          GR,GI, REAL AND IMAGINARY PARTS OF G
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: y
REAL (dp), INTENT(IN)   :: eps
REAL (dp), INTENT(OUT)  :: fr
REAL (dp), INTENT(OUT)  :: fi
REAL (dp), INTENT(OUT)  :: gr
REAL (dp), INTENT(OUT)  :: gi

INTEGER   :: a, b, k3
REAL (dp) :: x2, y2, u, v, p, q, cr, ci, dr, di

x2 = x * x
y2 = y * y
k3 = 0
u = x * (x2 - 3.d0*y2)
v = y * (3.d0*x2 - y2)
cr = 1.d0
ci = 0.d0
dr = 1.d0
di = 0.d0
fr = 1.d0
fi = 0.d0
gr = 1.d0
gi = 0.d0
10 a = (k3+2) * (k3+3)
b = (k3+4) * (k3+3)
p = (u*cr-v*ci) / a
q = (v*cr+u*ci) / a
cr = p
ci = q
p = (u*dr-v*di) / b
q = (v*dr+u*di) / b
dr = p
di = q
fr = fr + cr
fi = fi + ci
gr = gr + dr
gi = gi + di
k3 = k3 + 3
IF (ABS(cr) + ABS(dr) + ABS(ci) + ABS(di) >= eps) GO TO 10
p = x * gr - y * gi
q = x * gi + y * gr
gr = p
gi = q
RETURN
END SUBROUTINE fg


SUBROUTINE serai(x, y, air, aii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI(Z), TAYLOR, COMPLEX Z      CCC
!CC                                      CCC
!CC   INPUTS:                            CCC
!CC        X,Y,    REAL AND IMAGINARY    CCC
!CC                PARTS OF Z            CCC
!CC   OUTPUTS:                           CCC
!CC        AIR,AII, REAL AND IMAGINARY   CCC
!CC                 PARTS OF AI(Z)       CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: y
REAL (dp), INTENT(OUT)  :: air
REAL (dp), INTENT(OUT)  :: aii

REAL (dp) :: eps
REAL (dp) :: fzr, fzi, gzr, gzi, cons1, cons2

eps = MAX(EPSILON(0.0D0), 1.D-15)
cons1 = 0.355028053887817239260D0
cons2 = 0.258819403792806798405D0
CALL fg(x, y, eps, fzr, fzi, gzr, gzi)
air = cons1 * fzr - cons2 * gzr
aii = cons1 * fzi - cons2 * gzi
RETURN
END SUBROUTINE serai


SUBROUTINE seraid(x, y, air, aii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI'(Z), TAYLOR, COMPLEX Z     CCC
!CC                                      CCC
!CC   INPUTS:                            CCC
!CC        X,Y,   REAL AND IMAGINARY     CCC
!CC               PARTS OF Z             CCC
!CC   OUTPUTS:                           CCC
!CC        AIR,AII, REAL AND IMAGINARY   CCC
!CC                 PARTS OF AI'(Z)      CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: y
REAL (dp), INTENT(OUT)  :: air
REAL (dp), INTENT(OUT)  :: aii

REAL (dp) :: eps
REAL (dp) :: fzr, fzi, gzr, gzi, cons1, cons2

eps = MAX(EPSILON(0.0D0), 1.D-15)
cons1 = 0.355028053887817239260D0
cons2 = 0.258819403792806798405D0
CALL fgp(x, y, eps, fzr, fzi, gzr, gzi)
air = cons1 * fzr - cons2 * gzr
aii = cons1 * fzi - cons2 * gzi
RETURN
END SUBROUTINE seraid


SUBROUTINE expai(gair, gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI(Z), ASYMPTOTIC EXPANSION, COMPLEX Z CCC
!CC                                               CCC
!CC   OUTPUTS:                                    CCC
!CC        GAIR, GAII,  REAL AND IMAGINARY        CCC
!CC                     PARTS OF AI(Z)            CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(OUT)  :: gair
REAL (dp), INTENT(OUT)  :: gaii

REAL (dp) :: eps
REAL (dp) :: df1, psiir, psiii, ck, dfrr, dfii, sumar, sumai,  &
             dfr, dfi, deltar, deltai, fac1, fac2
REAL (dp) :: df
INTEGER   :: k

! COMMON /param3/ thet, r, th15, s1, c1, r32
! COMMON /param4/ facto, th025, s3, c3

REAL (dp), PARAMETER  :: co(20) = (/  &
   -6.944444444444445D-2, 3.713348765432099D-2, -3.799305912780064D-2,  &
    5.764919041266972D-2, -0.116099064025515D0, 0.291591399230751D0,   &
   -0.877666969510017D0, 3.07945303017317D0, -12.3415733323452D0,  &
    55.6227853659171D0, -278.465080777603D0, 1533.16943201280D0,  &
    -9207.20659972641D0, 59892.5135658791D0, -419524.875116551D0,  &
    3148257.41786683D0, -25198919.8716024D0, 214288036.963680D0,  &
    -1929375549.18249D0, 18335766937.8906D0 /)

eps = MAX(EPSILON(0.0D0), 1.D-15)
df1 = 1.5D0 / r32
psiir = c1
psiii = -s1
k = 0
ck = 1.d0
df = 1.d0
dfrr = 1.d0
dfii = 0.d0
sumar = 1.d0
sumai = 0.d0

10 df = df * df1
ck = co(k+1) * df
dfr = dfrr
dfi = dfii
dfrr = dfr * psiir - dfi * psiii
dfii = dfr * psiii + dfi * psiir
deltar = dfrr * ck
deltai = dfii * ck
sumar = sumar + deltar
sumai = sumai + deltai
k = k + 1
IF (sumar /= 0.0_dp) THEN
  IF (ABS(deltar/sumar) > eps) GO TO 10
END IF
IF (sumai /= 0.0_dp) THEN
  IF (ABS(deltai/sumai) > eps) GO TO 10
END IF
fac1 = facto * c3
fac2 = facto * s3
gair = fac1 * sumar + fac2 * sumai
gaii = fac1 * sumai - fac2 * sumar
RETURN
END SUBROUTINE expai


SUBROUTINE expaid(gair, gaii)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   AIRY AI'(Z), ASYMPTOTIC EXPANSION, COMPLEX Z CCC
!CC                                                CCC
!CC   OUTPUTS:                                     CCC
!CC        GAIR, GAII,  REAL AND IMAGINARY         CCC
!CC                     PARTS OF AI'(Z)            CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

REAL (dp), INTENT(OUT)  :: gair
REAL (dp), INTENT(OUT)  :: gaii

REAL (dp) :: eps
REAL (dp) :: df1, psiir, psiii, vk, dfrr, dfii, sumar, sumai,  &
             dfr, dfi, deltar, deltai, fac1, fac2
REAL (dp) :: df
INTEGER   :: k

! COMMON /param3/ thet, r, th15, s1, c1, r32
! COMMON /param4/ facto, th025, s3, c3

REAL (dp), PARAMETER  :: co(20) = (/  &
     9.722222222222222D-2, -4.388503086419753D-2, 4.246283078989484D-2,  &
    -6.266216349203230D-2, 0.124105896027275D0, -0.308253764901079D0,   &
     0.920479992412945D0, -3.21049358464862D0, 12.8072930807356D0,  &
    -57.5083035139143D0, 287.033237109221D0, -1576.35730333710D0,   &
     9446.35482309593D0, -61335.7066638521D0, 428952.400400069D0,   &
    -3214536.52140086D0, 25697908.3839113D0, -218293420.832160D0,   &
     1963523788.99103D0, -18643931088.1072D0 /)

eps = MAX(EPSILON(0.0D0), 1.D-15)
df1 = 1.5D0 / r32
psiir = c1
psiii = -s1
k = 0
df = 1.d0
dfrr = 1.d0
dfii = 0.d0
sumar = 1.d0
sumai = 0.d0
10 df = df * df1
vk = co(k+1) * df
dfr = dfrr
dfi = dfii
dfrr = dfr * psiir - dfi * psiii
dfii = dfr * psiii + dfi * psiir
deltar = dfrr * vk
deltai = dfii * vk
sumar = sumar + deltar
sumai = sumai + deltai
k = k + 1
IF (sumar /= 0) THEN
  IF (ABS(deltar/sumar) > eps) GO TO 10
END IF
IF (sumai /= 0) THEN
  IF (ABS(deltai/sumai) > eps) GO TO 10
END IF
fac1 = facto * c3
fac2 = facto * s3
gair = -(fac1*sumar - fac2*sumai)
gaii = -(fac1*sumai + fac2*sumar)
RETURN
END SUBROUTINE expaid


SUBROUTINE biz(ifun, ifac, x0, y0, gbir, gbii, ierro)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! COMPUTATION OF THE AIRY FUNCTION BI(Z) OR ITS DERIVATIVE BI'(Z)
! THE CODE USES THE CONNECTION OF BI(Z) WITH AI(Z).
!                BIZ CALLS THE ROUTINE AIZ
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  INPUTS:
!    IFUN:
!         * IFUN=1, THE CODE COMPUTES BI(Z)
!         * IFUN=2, THE CODE COMPUTES BI'(Z)
!    IFAC:
!         * IFAC=1, THE CODE COMPUTES  BI(Z) OR BI'(Z)
!         * IFAC=2, THE CODE COMPUTES NORMALIZED BI(Z) OR BI'(Z)
!    X0:   REAL PART OF THE ARGUMENT Z
!    Y0:   IMAGINARY PART OF THE ARGUMENT  Z

!  OUTPUTS:
!    GBIR: REAL PART OF BI(Z) OR BI'(Z)
!    GBII: IMAGINARY PART OF BI(Z) OR BI'(Z)

!    IERRO: ERROR FLAG
!          * IERRO=0, SUCCESSFUL COMPUTATION
!          * IERRO=1, COMPUTATION OUT OF RANGE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!          MACHINE DEPENDENT CONSTANTS: FUNCTION D1MACH
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!   ACCURACY:

!     1) SCALED AIRY FUNCTIONS:
!        RELATIVE ACCURACY BETTER THAN 10**(-13) EXCEPT CLOSE TO THE ZEROS,
!        WHERE 10**(-13) IS THE ABSOLUTE PRECISION.
!        GRADUAL LOSS OF PRECISION TAKES PLACE FOR |Z| > 1000
!        IN THE CASE OF PHASE(Z) CLOSE TO +3*PI/2 OR -3*PI/2.
!     2) UNSCALED AIRY FUNCTIONS:
!        THE FUNCTION OVERFLOWS/UNDERFLOWS FOR
!        3/2*|Z|**(3/2) > LOG(OVER).
!        FOR |Z| < 30:
!        A) RELATIVE ACCURACY FOR THE MODULUS (EXCEPT AT THE ZEROS)
!           BETTER THAN 10**(-13).
!        B) ABSOLUTE ACCURACY FOR MIN(R(Z),1/R(Z)) BETTER
!           THAN 10**(-13), WHERE R(Z)=REAL(BI)/IMAG(BI)
!           OR R(Z)=REAL(BI')/IMAG(BI').
!        FOR |Z| > 30, GRADUAL LOSS OF PRECISION TAKES PLACE AS |Z| INCREASES.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     AUTHORS:
!        AMPARO GIL    (U. AUTONOMA DE MADRID, MADRID, SPAIN).
!                      E-MAIL: AMPARO.GIL@UAM.ES
!        JAVIER SEGURA (U. CARLOS III DE MADRID, MADRID, SPAIN).
!                      E-MAIL: JSEGURA@MATH.UC3M.ES
!        NICO M. TEMME (CWI, AMSTERDAM, THE NETHERLANDS).
!                      E-MAIL: NICO.TEMME@CWI.NL
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!    REFERENCES:
!         COMPUTING AIRY FUNCTIONS BY NUMERICAL QUADRATURE.
!         NUMERICAL ALGORITHMS (2002).
!         A. GIL, J. SEGURA, N.M. TEMME
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

INTEGER, INTENT(IN)     :: ifun
INTEGER, INTENT(IN)     :: ifac
REAL (dp), INTENT(IN)   :: x0
REAL (dp), INTENT(IN)   :: y0
REAL (dp), INTENT(OUT)  :: gbir
REAL (dp), INTENT(OUT)  :: gbii
INTEGER, INTENT(OUT)    :: ierro

REAL (dp) :: over, under, dl1, dl2, cover
REAL (dp) :: pi3, pi23, sqrt3, c, s, c1, s1, u, v, x, y,  &
             ar, ai, apr, api, br, bi, bpr, bpi, bbr, bbi, bbpr, bbpi
REAL (dp) :: thet, r, r32, thet32, a1, b1, df1, expo, expoi,  &
             etar, etai, etagr, etagi
INTEGER :: iexpf, ierr
! COMMON /param1/ pi, pihal

sqrt3 = 1.7320508075688772935D0
pi3 = pi / 3.d0
pi23 = 2.d0 * pi3
x = x0
c = 0.5D0 * sqrt3
s = 0.5D0
ierro = 0
iexpf = 0
IF (y0 < 0.d0) THEN
  y = -y0
ELSE
  y = y0
END IF
r = SQRT(x*x + y*y)
r32 = r * SQRT(r)
thet = phase(x,y)
cover = 2.d0 / 3.d0 * r32 * ABS(COS(1.5D0*thet))

!CC MACHINE DEPENDENT CONSTANT (OVERFLOW NUMBER)
over = HUGE(0.0D0) * 1.d-4
!CC MACHINE DEPENDENT CONSTANT (UNDERFLOW NUMBER)
under = TINY(0.0D0) * 1.d+4
dl1 = LOG(over)
dl2 = -LOG(under)
IF (dl1 > dl2) over = 1.0_dp / under
IF (ifac == 1) THEN
  IF (cover >= LOG(over)) THEN
!CC OVERFLOW/UNDERFLOW PROBLEMS.
!CC   CALCULATION ABORTED
    ierro = 1
    gbir = 0
    gbii = 0
  END IF
ELSE
  IF (cover >= LOG(over)*0.2) iexpf = 1
END IF

IF (ierro == 0) THEN
  IF (ifac == 1) THEN
    IF (y == 0.d0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC  TAKE TWICE THE REAL PART OF EXP(-PI I/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      c1 = -0.5D0
      s1 = -0.5D0 * sqrt3
      u = x * c1 - y * s1
      v = x * s1 + y * c1
      CALL aiz(ifun, ifac, u, v, ar, ai, ierr)
      IF (ifun == 1) THEN
        br = sqrt3 * ar + ai
        bi = 0.d0
      ELSE
        u = ar * c1 - ai * s1
        v = ar * s1 + ai * c1
        apr = u
        api = v
        bpr = sqrt3 * apr + api
        bpi = 0.d0
      END IF
    ELSE
      IF ((x < 0.d0) .AND. (y < -x*sqrt3)) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC      2 PI/3  < PHASE(Z) < PI
!CC      BI(Z)=EXP(I PI/6) AI_(-1)(Z) + EXP(-I PI/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        c1 = -0.5D0
        s1 = 0.5D0 * sqrt3
        u = x * c1 - y * s1
        v = x * s1 + y * c1
        CALL aiz(ifun, ifac, u, v, ar, ai, ierr)
        IF (ifun == 1) THEN
          br = c * ar - s * ai
          bi = c * ai + s * ar
        ELSE
          u = ar * c1 - ai * s1
          v = ar * s1 + ai * c1
          apr = u
          api = v
          bpr = c * apr - s * api
          bpi = c * api + s * apr
        END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   WE NEED ALSO AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        c1 = -0.5D0
        s1 = -0.5D0 * sqrt3
        u = x * c1 - y * s1
        v = x * s1 + y * c1
        CALL aiz(ifun, ifac, u, v, ar, ai, ierr)
        IF (ifun == 1) THEN
          s = -s
          br = br + c * ar - s * ai
          bi = bi + c * ai + s * ar
        ELSE
          u = ar * c1 - ai * s1
          v = ar * s1 + ai * c1
          apr = u
          api = v
          s = -s
          bpr = bpr + c * apr - s * api
          bpi = bpi + c * api + s * apr
        END IF
      ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   BI(Z) = I AI(Z) + 2 EXP(-I PI/6) AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        c1 = -0.5D0
        s1 = -0.5D0 * sqrt3
        u = x * c1 - y * s1
        v = x * s1 + y * c1
        CALL aiz(ifun, ifac, u, v, ar, ai, ierr)
        IF (ifun == 1) THEN
          br = sqrt3 * ar + ai
          bi = -ar + sqrt3 * ai
        ELSE
          u = ar * c1 - ai * s1
          v = ar * s1 + ai * c1
          apr = u
          api = v
          bpr = sqrt3 * apr + api
          bpi = -apr + sqrt3 * api
        END IF
        CALL aiz(ifun, ifac, x, y, ar, ai, ierr)
        IF (ifun == 1) THEN
          br = br - ai
          bi = bi + ar
        ELSE
          bpr = bpr - ai
          bpi = bpi + ar
        END IF
      END IF
    END IF
  ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC         SCALED BI AIRY FUNCTIONS         C
!CC   WE USE THE FOLLOWING NORMALIZATION:    C
!CC   LET ARGZ=ARG(Z), THEN:                 C
!CC   A) IF  0 <= ARGZ <= PI/3               C
!CC      BI=EXP(-2/3Z^3/2)BI                 C
!CC   B) IF  PI/3 <= ARGZ <= PI              C
!CC      BI=EXP(2/3Z^3/2)BI                  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    thet = phase(x, y)
    IF (thet <= pi3) THEN
      c1 = -0.5D0
      s1 = -0.5D0 * sqrt3
      u = x * c1 - y * s1
      v = x * s1 + y * c1
      CALL aiz(ifun, ifac, u, v, ar, ai, ierr)
      IF (ifun == 1) THEN
        br = sqrt3 * ar + ai
        bi = -ar + sqrt3 * ai
      ELSE
        u = ar * c1 - ai * s1
        v = ar * s1 + ai * c1
        apr = u
        api = v
        bpr = sqrt3 * apr + api
        bpi = -apr + sqrt3 * api
      END IF
      IF (iexpf == 0) THEN
        r = SQRT(x*x + y*y)
        r32 = r * SQRT(r)
        thet32 = thet * 1.5D0
        a1 = COS(thet32)
        b1 = SIN(thet32)
        df1 = 4.d0 / 3.d0 * r32
        expo = EXP(df1*a1)
        expoi = 1.d0 / expo
        etar = expo * COS(df1*b1)
        etai = expo * SIN(df1*b1)
        etagr = expoi * COS(-df1*b1)
        etagi = expoi * SIN(-df1*b1)
        CALL aiz(ifun, ifac, x, y, ar, ai, ierr)
        IF (ifun == 1) THEN
          br = br - ar * etagi - etagr * ai
          bi = bi + ar * etagr - etagi * ai
        ELSE
          bpr = bpr - ar * etagi - etagr * ai
          bpi = bpi + ar * etagr - etagi * ai
        END IF
      END IF
    END IF
    IF (thet > pi3 .AND. thet <= pi23) THEN
      IF (iexpf == 0) THEN
        r = SQRT(x*x+y*y)
        r32 = r * SQRT(r)
        thet32 = thet * 1.5D0
        a1 = COS(thet32)
        b1 = SIN(thet32)
        df1 = 4.d0 / 3.d0 * r32
        expo = EXP(df1*a1)
        expoi = 1.d0 / expo
        etar = expo * COS(df1*b1)
        etai = expo * SIN(df1*b1)
        etagr = expoi * COS(-df1*b1)
        etagi = expoi * SIN(-df1*b1)
        c1 = -0.5D0
        s1 = -0.5D0 * sqrt3
        u = x * c1 - y * s1
        v = x * s1 + y * c1
        CALL aiz(ifun, ifac, u, v, ar, ai, ierr)
        IF (ifun == 1) THEN
          bbr = sqrt3 * ar + ai
          bbi = -ar + sqrt3 * ai
          br = bbr * etar - bbi * etai
          bi = bbi * etar + bbr * etai
        ELSE
          u = ar * c1 - ai * s1
          v = ar * s1 + ai * c1
          apr = u
          api = v
          bbpr = sqrt3 * apr + api
          bbpi = -apr + sqrt3 * api
          bpr = bbpr * etar - bbpi * etai
          bpi = bbpi * etar + bbpr * etai
        END IF
      ELSE
        IF (ifun == 1) THEN
          br = 0.d0
          bi = 0.d0
        ELSE
          bpr = 0.d0
          bpi = 0.d0
        END IF
      END IF
      CALL aiz(ifun, ifac, x, y, ar, ai, ierr)
      IF (ifun == 1) THEN
        br = br - ai
        bi = bi + ar
      ELSE
        bpr = bpr - ai
        bpi = bpi + ar
      END IF
    END IF
    IF (thet > pi23) THEN
      c1 = -0.5D0
      s1 = 0.5D0 * sqrt3
      u = x * c1 - y * s1
      v = x * s1 + y * c1
      CALL aiz(ifun, ifac, u, v, ar, ai, ierr)
      IF (ifun == 1) THEN
        br = c * ar - s * ai
        bi = c * ai + s * ar
      ELSE
        u = ar * c1 - ai * s1
        v = ar * s1 + ai * c1
        apr = u
        api = v
        bpr = c * apr - s * api
        bpi = c * api + s * apr
      END IF
      IF (iexpf == 0) THEN
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   WE NEED ALSO AI_(1)(Z)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        r = SQRT(x*x+y*y)
        r32 = r * SQRT(r)
        thet32 = thet * 1.5D0
        a1 = COS(thet32)
        b1 = SIN(thet32)
        df1 = 4.d0 / 3.d0 * r32
        expo = EXP(df1*a1)
        expoi = 1.d0 / expo
        etar = expo * COS(df1*b1)
        etai = expo * SIN(df1*b1)
        etagr = expoi * COS(-df1*b1)
        etagi = expoi * SIN(-df1*b1)
        c1 = -0.5D0
        s1 = -0.5D0 * sqrt3
        u = x * c1 - y * s1
        v = x * s1 + y * c1
        CALL aiz(ifun, ifac, u, v, ar, ai, ierr)
        IF (ifun == 1) THEN
          s = -s
          bbr = c * ar - s * ai
          bbi = c * ai + s * ar
          br = br + etar * bbr - etai * bbi
          bi = bi + bbi * etar + etai * bbr
        ELSE
          u = ar * c1 - ai * s1
          v = ar * s1 + ai * c1
          apr = u
          api = v
          s = -s
          bbpr = c * apr - s * api
          bbpi = c * api + s * apr
          bpr = bpr + etar * bbpr - etai * bbpi
          bpi = bpi + bbpi * etar + etai * bbpr
        END IF
      END IF
    END IF
  END IF
  IF (y0 < 0) THEN
    bi = -bi
    bpi = -bpi
  END IF
  IF (ifun == 1) THEN
    gbir = br
    gbii = bi
  ELSE
    gbir = bpr
    gbii = bpi
  END IF
END IF
RETURN
END SUBROUTINE biz

END MODULE AiryFunction
