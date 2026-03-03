PROGRAM test_biv_norm
! Compare 3 algorithms for the bivariate normal

USE Owens_T
IMPLICIT NONE

INTERFACE
  FUNCTION bivnor(ah, ak, r) RESULT(b)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: ah
    REAL (dp), INTENT(IN)  :: ak
    REAL (dp), INTENT(IN)  :: r
    REAL (dp)              :: b
  END FUNCTION bivnor

  FUNCTION bvn ( lower, upper, infin, correl ) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: lower(:), upper(:), correl
    INTEGER, INTENT(IN)    :: infin(:)
    REAL (dp)              :: fn_val
  END FUNCTION bvn
END INTERFACE

INTEGER    :: infin(2), i
REAL (dp)  :: p1, p2, p3, lower(2), upper(2), correl, rand(3)
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

upper = 100.0_dp
infin = 1
WRITE(*, *) '   X       Y    Correl.      Donnelly        Genz        Patefield'
DO i = 1, 100
  CALL RANDOM_NUMBER( rand )
  lower(1) = 8.*(rand(1) - 0.5)
  lower(2) = 8.*(rand(2) - 0.5)
  correl   = 2.*(rand(3) - 0.5)
  p1 = bivnor(lower(1), lower(2), correl)
  p2 = bvn(lower, upper, infin, correl)
  p3 = biv_norm(lower(1), zero, one, lower(2), zero, one, correl)
  WRITE(*, '(3f8.4, "    ", 3g14.6)') lower, correl, p1, p2, p3
END DO

STOP
END PROGRAM test_biv_norm



FUNCTION bivnor(ah, ak, r) RESULT(b)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-04-02  Time: 15:44:21

!     BIVNOR is a controlled precision Fortran function to calculate the
!     bivariate normal upper right area, viz. the probability for two normal
!     variates X and Y whose correlation is R, that X > AH and Y > AK.
!     The accuracy is specified as the number of decimal digits, IDIG.

!     Reference:
!        Donnelly, T.G. (1973) Algorithm 462, Bivariate normal
!        distribution, Comm. A.C.M., vol.16, p. 638.

!     Corrected for the case AH = 0, AH non-zero, by reversing AH & AK
!     in such cases.

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: ah
REAL (dp), INTENT(IN)  :: ak
REAL (dp), INTENT(IN)  :: r
REAL (dp)              :: b

!     Local variables

REAL (dp) :: xah, xak, gh, gk, rr, h2, a2, h4, ex, w2, &
             ap, s2, sp, s1, sn, sgn, sqr, con, wh, wk, gw, t, g2, conex, cn
INTEGER   :: is

REAL (dp), PARAMETER :: two = 2._dp, zero = 0._dp, one = 1._dp, four = 4._dp,  &
                        quart = 0.25_dp, half = 0.5_dp

REAL (dp), PARAMETER :: twopi = 6.283185307179587_dp, explim = 80._dp
INTEGER, PARAMETER   :: idig = 15

b = zero
IF (ah == zero) THEN
  xah = ak
  xak = ah
ELSE
  xah = ah
  xak = ak
END IF

gh = gauss(-xah) / two
gk = gauss(-xak) / two
IF (r == zero) THEN
  b = four * gh * gk
  GO TO 350
END IF

rr = one - r*r
IF (rr < zero) THEN
  WRITE(*, *)'Error in BIVNOR, R = ', r
  GO TO 390
END IF
IF (rr > zero) GO TO 100

!     R^2 = 1.0

IF (r >= zero) GO TO 70
IF (xah + xak >= zero) GO TO 350
b = two * (gh + gk) - one
GO TO 350
70 IF (xah - xak < zero) THEN
  b = two * gk
ELSE
  b = two * gh
END IF
GO TO 350

!     Regular case, R^2 < 1

100 sqr = SQRT(rr)
IF (idig == 15) THEN
  con = twopi * 1.d-15 / two
ELSE
  con = twopi / two / 10**idig
END IF

IF (xah /= zero) GO TO 170
IF (xak /= zero) GO TO 190
b = ATAN(r/sqr) / twopi + quart
GO TO 350
170 b = gh

IF (xah*xak < 0.0) THEN
  GO TO   180
ELSE IF (xah*xak == 0.0) THEN
  GO TO   200
ELSE
  GO TO   190
END IF
180 b = b - half
190 b = b + gk
200 wh = -xah
wk = (xak/xah - r)/sqr
gw = two * gh
is = -1
210 sgn = -one
t = zero
IF (wk == zero) GO TO 320
IF (ABS(wk) - one < 0.0) THEN
  GO TO   270
ELSE IF (ABS(wk) - one == 0.0) THEN
  GO TO   230
ELSE
  GO TO   240
END IF
230 t = wk * gw * (one - gw) / two
GO TO 310
240 sgn = -sgn
wh = wh * wk
g2 = gauss(wh)
wk = one / wk
IF (wk < zero) b = b + half
b = b - (gw + g2)/two + gw*g2
270 h2 = wh * wh
a2 = wk * wk
h4 = h2 / two
IF (h4 < explim) THEN
  ex = EXP(-h4)
ELSE
  ex = zero
END IF
w2 = h4 * ex
ap = one
s2 = ap - ex
sp = ap
s1 = zero
sn = s1
conex = ABS(con / wk)
GO TO 290

280 sn = sp
sp = sp + one
s2 = s2 - w2
w2 = w2 * h4 / sp
ap = -ap*a2
290 cn = ap * s2 / (sn + sp)
s1 = s1 + cn
IF (ABS(cn) - conex > zero) GO TO 280

t = (ATAN(wk) - wk*s1) / twopi
310 b = b + sgn*t
320 IF (is >= 0) GO TO 350
IF (xak /= zero) THEN
  wh = -xak
  wk = (xah/xak - r) / sqr
  gw = two * gk
  is = 1
  GO TO 210
END IF

350 IF (b < zero) b = zero
IF (b > one) b = one

390 RETURN

CONTAINS


FUNCTION gauss(t) RESULT(fn_val)

!     GAUSS is a univariate lower normal tail area calculated here from
!     the central error function, DERF.
!     It may be replaced by the algorithm in Hill, I.D. and Joyce, S.A.
!     Algorithm 304, Normal curve integral (S15), Comm. A.C.M. (10),
!     (June 1967), p.374 or with Applied Statistics algorithm AS66.
!     Source: Owen, D.B., Ann. Math. Statist., vol.27 (1956), p.1075.

REAL (dp), INTENT(IN) :: t
REAL (dp)             :: fn_val

fn_val = (one + derf(t/SQRT(two))) / two

RETURN
END FUNCTION gauss



FUNCTION derf(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!        REAL (dp) EVALUATION OF THE ERROR FUNCTION
!
!        WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: ax, t, w
INTEGER   :: k
REAL (dp), PARAMETER :: a(21) = (/ .1283791670955125738961589031215D+00,  &
        -.3761263890318375246320529677070D+00,  &
         .1128379167095512573896158902931D+00,  &
        -.2686617064513125175943235372542D-01,  &
         .5223977625442187842111812447877D-02,  &
        -.8548327023450852832540164081187D-03,  &
         .1205533298178966425020717182498D-03,  &
        -.1492565035840625090430728526820D-04,  &
         .1646211436588924261080723578109D-05,  &
        -.1636584469123468757408968429674D-06,  &
         .1480719281587021715400818627811D-07,  &
        -.1229055530145120140800510155331D-08,  &
         .9422759058437197017313055084212D-10,  &
        -.6711366740969385085896257227159D-11,  &
         .4463222608295664017461758843550D-12,  &
        -.2783497395542995487275065856998D-13,  &
         .1634095572365337143933023780777D-14,  &
        -.9052845786901123985710019387938D-16,  &
         .4708274559689744439341671426731D-17,  &
        -.2187159356685015949749948252160D-18,  &
         .7043407712019701609635599701333D-20 /)
!-------------------------------

!                     ABS(X) <= 1

ax = ABS(x)
IF (ax <= one) THEN
  t = x * x
  w = a(21)
  DO k = 20, 1, -1
    w = t * w + a(k)
  END DO
  fn_val = x * (one + w)
  RETURN
END IF

!                     ABS(X) > 1

IF (ax < 8.5_dp) THEN
  fn_val = half + (half - EXP(-x*x)*derfc0(ax))
  IF (x < zero) fn_val = -fn_val
  RETURN
END IF

!                 LIMIT VALUE FOR LARGE X

fn_val = SIGN(one, x)
RETURN
END FUNCTION derf



FUNCTION derfc0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

!           EVALUATION OF EXP(X**2)*ERFC(X) FOR X >= 1

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!        APRIL 1992
!-------------------------------
REAL (dp)            :: t, u, v, z
REAL (dp), PARAMETER :: rpinv = .56418958354775628694807945156077259D0
REAL (dp), PARAMETER :: p0 = .16506148041280876191828601D-03,  &
                        p1 =  .15471455377139313353998665D-03,  &
                        p2 =  .44852548090298868465196794D-04,  &
                        p3 = -.49177280017226285450486205D-05,  &
                        p4 = -.69353602078656412367801676D-05,  &
                        p5 = -.20508667787746282746857743D-05,  &
                        p6 = -.28982842617824971177267380D-06,  &
                        p7 = -.17272433544836633301127174D-07,  &
                        q1 =  .16272656776533322859856317D+01,  &
                        q2 =  .12040996037066026106794322D+01,  &
                        q3 =  .52400246352158386907601472D+00,  &
                        q4 =  .14497345252798672362384241D+00,  &
                        q5 =  .25592517111042546492590736D-01,  &
                        q6 =  .26869088293991371028123158D-02,  &
                        q7 =  .13133767840925681614496481D-03
REAL (dp), PARAMETER :: r0 =  .145589721275038539045668824025D+00,  &
                        r1 = -.273421931495426482902320421863D+00,  &
                        r2 =  .226008066916621506788789064272D+00,  &
                        r3 = -.163571895523923805648814425592D+00,  &
                        r4 =  .102604312032193978662297299832D+00,  &
                        r5 = -.548023266949835519254211506880D-01,  &
                        r6 =  .241432239725390106956523668160D-01,  &
                        r7 = -.822062115403915116036874169600D-02,  &
                        r8 =  .180296241564687154310619200000D-02
REAL (dp), PARAMETER :: a0 = -.45894433406309678202825375D-03,   &
                        a1 = -.12281298722544724287816236D-01,  &
                        a2 = -.91144359512342900801764781D-01,  &
                        a3 = -.28412489223839285652511367D-01,  &
                        a4 =  .14083827189977123530129812D+01,  &
                        a5 =  .11532175281537044570477189D+01,  &
                        a6 = -.72170903389442152112483632D+01,  &
                        a7 = -.19685597805218214001309225D+01,  &
                        a8 =  .93846891504541841150916038D+01,  &
                        b1 =  .25136329960926527692263725D+02,  &
                        b2 =  .15349442087145759184067981D+03,  &
                        b3 = -.29971215958498680905476402D+03,  &
                        b4 = -.33876477506888115226730368D+04,  &
                        b5 =  .28301829314924804988873701D+04,  &
                        b6 =  .22979620942196507068034887D+05,  &
                        b7 = -.24280681522998071562462041D+05,  &
                        b8 = -.36680620673264731899504580D+05,  &
                        b9 =  .42278731622295627627042436D+05,  &
                        b10=  .28834257644413614344549790D+03,  &
                        b11=  .70226293775648358646587341D+03
REAL (dp), PARAMETER :: c0 = -.7040906288250128001000086D-04,   &
                        c1 = -.3858822461760510359506941D-02,  &
                        c2 = -.7708202127512212359395078D-01,  &
                        c3 = -.6713655014557429480440263D+00,  &
                        c4 = -.2081992124162995545731882D+01,  &
                        c5 =  .2898831421475282558867888D+01,  &
                        c6 =  .2199509380600429331650192D+02,  &
                        c7 =  .2907064664404115316722996D+01,  &
                        c8 = -.4766208741588182425380950D+02,  &
                        d1 =  .5238852785508439144747174D+02,  &
                        d2 =  .9646843357714742409535148D+03,  &
                        d3 =  .7007152775135939601804416D+04,  &
                        d4 =  .8515386792259821780601162D+04,  &
                        d5 = -.1002360095177164564992134D+06,  &
                        d6 = -.2065250031331232815791912D+06,  &
                        d7 =  .5695324805290370358175984D+06,  &
                        d8 =  .6589752493461331195697873D+06,  &
                        d9 = -.1192930193156561957631462D+07
REAL (dp), PARAMETER :: e0 = .540464821348814822409610122136D+00,  &
                        e1 = -.261515522487415653487049835220D-01, &
                        e2 = -.288573438386338758794591212600D-02, &
                        e3 = -.529353396945788057720258856000D-03
REAL (dp), PARAMETER :: s1 = .75000000000000000000D+00,   &
        s2  = -.18750000000000000000D+01, s3  = .65625000000000000000D+01,  &
        s4  = -.29531250000000000000D+02, s5  = .16242187500000000000D+03,  &
        s6  = -.10557421875000000000D+04, s7  = .79180664062500000000D+04,  &
        s8  = -.67303564453125000000D+05, s9  = .63938386230468750000D+06,  &
        s10 = -.67135305541992187500D+07, s11 = .77205601373291015625D+08
!-------------------------------
!     RPINV = 1/SQRT(PI)
!-------------------------------

!                     1 <= X <= 2

IF (x <= 2._dp) THEN
  u = ((((((p7*x + p6)*x + p5)*x + p4)*x + p3)*x + p2)*x + p1) * x + p0
  v = ((((((q7*x + q6)*x + q5)*x + q4)*x + q3)*x + q2)*x + q1) * x + 1._dp

  t = (x-3.75D0) / (x+3.75D0)
  fn_val = (((((((((u/v)*t + r8)*t + r7)*t + r6)*t + r5)*t + r4)*t + r3)*t + &
           r2)*t + r1) * t + r0
  RETURN
END IF

!                     2 < X <= 4

IF (x <= 4._dp) THEN
  z = 1._dp / (2.5_dp + x*x)
  u = (((((((a8*z + a7)*z + a6)*z + a5)*z + a4)*z + a3)*z + a2)*z + a1) * z + a0
  v = ((((((((((b11*z + b10)*z + b9)*z + b8)*z + b7)*z + b6)*z + b5)*z +  &
      b4)*z + b3)*z + b2)*z + b1) * z + 1._dp

  t = 13._dp * z - 1._dp
  fn_val = ((((u/v)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                     4 < X < 50

IF (x < 50._dp) THEN
  z = 1._dp / (2.5_dp + x*x)
  u = (((((((c8*z + c7)*z + c6)*z + c5)*z + c4)*z + c3)*z + c2)*z + c1) * z + &
      c0
  v = ((((((((d9*z + d8)*z + d7)*z + d6)*z + d5)*z + d4)*z + d3)*z + d2)*z +  &
      d1)*z + 1._dp
  t = 13._dp * z - 1._dp
  fn_val = (((((u/v)*t + e3)*t + e2)*t + e1)*t + e0) / x
  RETURN
END IF

!                        X >= 50

t = (1._dp/x) ** 2
z = (((((((((((s11*t + s10)*t + s9)*t + s8)*t + s7)*t + s6)*t + s5)*t +  &
    s4)*t + s3)*t + s2)*t + s1)*t - 0.5_dp) * t + 1._dp
fn_val = rpinv * (z/x)
RETURN
END FUNCTION derfc0

END FUNCTION bivnor



FUNCTION bvn ( lower, upper, infin, correl ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.
!     Extracted from Alan Genz's package for multivariate normal integration.

!  Parameters

!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, correlation coefficient.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN) :: lower(:), upper(:), correl
INTEGER, INTENT(IN)   :: infin(:)
REAL (dp)             :: fn_val

REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

IF ( infin(1) == 2  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )    &
            - bvnu ( upper(1), lower(2), correl )  &
            - bvnu ( lower(1), upper(2), correl )  &
            + bvnu ( upper(1), upper(2), correl )
ELSE IF ( infin(1) == 2  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )  &
            - bvnu ( upper(1), lower(2), correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )  &
            - bvnu ( lower(1), upper(2), correl )
ELSE IF ( infin(1) == 2  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )  &
            - bvnu ( -lower(1), -upper(2), correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 2 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )  &
            - bvnu ( -upper(1), -lower(2), correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( lower(1), -upper(2), -correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( -upper(1), lower(2), -correl )
ELSE IF ( infin(1) == 1  .AND. infin(2) == 1 ) THEN
  fn_val =  bvnu ( lower(1), lower(2), correl )
ELSE IF ( infin(1) == 0  .AND. infin(2) == 0 ) THEN
  fn_val =  bvnu ( -upper(1), -upper(2), correl )
END IF

RETURN


CONTAINS


FUNCTION bvnu( sh, sk, r ) RESULT(fn_val)

!     A function for computing bivariate normal probabilities.

!       Yihong Ge
!       Department of Computer Science and Electrical Engineering
!       Washington State University
!       Pullman, WA 99164-2752
!       Email : yge@eecs.wsu.edu
!     and
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu

! BVN - calculate the probability that X is larger than SH and Y is
!       larger than SK.

! Parameters

!   SH  REAL, integration limit
!   SK  REAL, integration limit
!   R   REAL, correlation coefficient
!   LG  INTEGER, number of Gauss Rule Points and Weights

REAL (dp), INTENT(IN) :: sh, sk, r
REAL (dp)             :: fn_val

! Local variables
INTEGER              :: i, lg, ng
REAL (dp), PARAMETER :: twopi = 6.283185307179586
REAL (dp)            :: as, a, b, c, d, rs, xs
REAL (dp)            :: bvn, sn, asr, h, k, bs, hs, hk
!     Gauss Legendre Points and Weights, N =  6
! DATA ( w(i,1), x(i,1), i = 1,3) /  &
! 0.1713244923791705D+00,-0.9324695142031522D+00,  &
! 0.3607615730481384D+00,-0.6612093864662647D+00,  &
! 0.4679139345726904D+00,-0.2386191860831970D+00/
!     Gauss Legendre Points and Weights, N = 12
! DATA ( w(i,2), x(i,2), i = 1,6) /  &
! 0.4717533638651177D-01,-0.9815606342467191D+00,  &
! 0.1069393259953183D+00,-0.9041172563704750D+00,  &
! 0.1600783285433464D+00,-0.7699026741943050D+00,  &
! 0.2031674267230659D+00,-0.5873179542866171D+00,  &
! 0.2334925365383547D+00,-0.3678314989981802D+00,  &
! 0.2491470458134029D+00,-0.1252334085114692D+00/
!     Gauss Legendre Points and Weights, N = 20
! DATA ( w(i,3), x(i,3), i = 1,10) /  &
! 0.1761400713915212D-01,-0.9931285991850949D+00,  &
! 0.4060142980038694D-01,-0.9639719272779138D+00,  &
! 0.6267204833410906D-01,-0.9122344282513259D+00,  &
! 0.8327674157670475D-01,-0.8391169718222188D+00,  &
! 0.1019301198172404D+00,-0.7463319064601508D+00,  &
! 0.1181945319615184D+00,-0.6360536807265150D+00,  &
! 0.1316886384491766D+00,-0.5108670019508271D+00,  &
! 0.1420961093183821D+00,-0.3737060887154196D+00,  &
! 0.1491729864726037D+00,-0.2277858511416451D+00,  &
! 0.1527533871307259D+00,-0.7652652113349733D-01/
REAL (dp), PARAMETER :: w(10,3) = RESHAPE( (/      &
      0.1713244923791705D+00, 0.3607615730481384D+00, 0.4679139345726904D+00, &
        0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,                  &
      0.4717533638651177D-01, 0.1069393259953183D+00, 0.1600783285433464D+00, &
      0.2031674267230659D+00, 0.2334925365383547D+00, 0.2491470458134029D+00, &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      0.1761400713915212D-01, 0.4060142980038694D-01, 0.6267204833410906D-01, &
      0.8327674157670475D-01, 0.1019301198172404D+00, 0.1181945319615184D+00, &
      0.1316886384491766D+00, 0.1420961093183821D+00, 0.1491729864726037D+00, &
      0.1527533871307259D+00 /), (/ 10, 3 /) )
REAL (dp), PARAMETER :: x(10,3) = RESHAPE( (/      &
      -0.9324695142031522D+00, -0.6612093864662647D+00,   &
      -0.2386191860831970D+00,  0.0D0, 0.0D0, 0.0D0,      &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      -0.9815606342467191D+00, -0.9041172563704750D+00,   &
      -0.7699026741943050D+00, -0.5873179542866171D+00,   &
      -0.3678314989981802D+00, -0.1252334085114692D+00,   &
       0.0D0, 0.0D0, 0.0D0, 0.0D0,                        &
      -0.9931285991850949D+00, -0.9639719272779138D+00,   &
      -0.9122344282513259D+00, -0.8391169718222188D+00,   &
      -0.7463319064601508D+00, -0.6360536807265150D+00,   &
      -0.5108670019508271D+00, -0.3737060887154196D+00,   &
      -0.2277858511416451D+00, -0.7652652113349733D-01 /), (/ 10, 3 /) )

IF ( ABS(r) < 0.3 ) THEN
  ng = 1
  lg = 3
ELSE IF ( ABS(r) < 0.75 ) THEN
  ng = 2
  lg = 6
ELSE
  ng = 3
  lg = 10
END IF
h = sh
k = sk
hk = h*k
bvn = zero
IF ( ABS(r) < 0.925 ) THEN
  hs = ( h*h + k*k )/2
  asr = ASIN(r)
  DO  i = 1, lg
    sn = SIN(asr*( x(i,ng)+1 )/2)
    bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/(one - sn*sn ) )
    sn = SIN(asr*(-x(i,ng)+1 )/2)
    bvn = bvn + w(i,ng)*EXP( ( sn*hk - hs )/(one - sn*sn ) )
  END DO
  bvn = bvn*asr/(2*twopi) + phi(-h)*phi(-k)
ELSE
  IF ( r < zero ) THEN
    k = -k
    hk = -hk
  END IF
  IF ( ABS(r) < one ) THEN
    as = ( one - r )*( one + r )
    a = SQRT(as)
    bs = ( h - k )**2
    c = ( 4. - hk )/8
    d = ( 12. - hk )/16.
    bvn = a*EXP( -(bs/as + hk)/2. )  &
    *( one - c*(bs - as)*(one - d*bs/5.)/3. + c*d*as*as/5. )
    IF ( hk > -160. ) THEN
      b = SQRT(bs)
      bvn = bvn - EXP(-hk/2)*SQRT(twopi)*phi(-b/a)*b  &
                      *( one - c*bs*( one - d*bs/5. )/3. )
    END IF
    a = a/2
    DO i = 1, lg
      xs = ( a*(x(i,ng) + one) )**2
      rs = SQRT( one - xs )
      bvn = bvn + a*w(i,ng)*( EXP( -bs/(2*xs) - hk/(1+rs) )/rs  &
                - EXP( -(bs/xs+hk)/2. )*( one + c*xs*( one + d*xs ) ) )
      xs = as*(-x(i,ng) + one)**2/4.
      rs = SQRT( 1 - xs )
      bvn = bvn + a*w(i,ng)*EXP( -(bs/xs + hk)/2 )  &
                * ( EXP( -hk*(one - rs)/(2*(one + rs)) )/rs - &
                       ( one + c*xs*( one + d*xs ) ) )
    END DO
    bvn = -bvn/twopi
  END IF
  IF ( r > 0 ) bvn =  bvn + phi( -MAX( h, k ) )
  IF ( r < 0 ) bvn = -bvn + MAX( zero, phi(-h) - phi(-k) )
END IF
fn_val = bvn

RETURN
END FUNCTION bvnu



SUBROUTINE normp(z, p, q, pdf)

! Normal distribution probabilities accurate to 1.e-15.
! Z = no. of standard deviations from the mean.
! P, Q = probabilities to the left & right of Z.   P + Q = 1.
! PDF = the probability density.

! Based upon algorithm 5666 for the error function, from:
! Hart, J.F. et al, 'Computer Approximations', Wiley 1968

! Programmer: Alan Miller

! Latest revision of Fortran 77 version - 30 March 1986
! Latest revision of Fortran 90 version - 12 August 1997

IMPLICIT NONE
REAL (dp), INTENT(IN)            :: z
REAL (dp), INTENT(OUT), OPTIONAL :: p, q, pdf

! Local variables
REAL (dp), PARAMETER :: p0 = 220.2068679123761D0, p1 = 221.2135961699311D0,  &
                        p2 = 112.0792914978709D0, p3 = 33.91286607838300D0,  &
                        p4 = 6.373962203531650D0, p5 = .7003830644436881D0,  &
                        p6 = .3526249659989109D-01,  &
                        q0 = 440.4137358247522D0, q1 = 793.8265125199484D0,  &
                        q2 = 637.3336333788311D0, q3 = 296.5642487796737D0,  &
                        q4 = 86.78073220294608D0, q5 = 16.06417757920695D0,  &
                        q6 = 1.755667163182642D0, q7 = .8838834764831844D-1, &
                        cutoff = 7.071D0, root2pi = 2.506628274631001D0
REAL (dp)            :: zabs, expntl, pp, qq, ppdf

zabs = ABS(z)

! |Z| > 37.

IF (zabs > 37.d0) THEN
  IF (PRESENT(pdf)) pdf = zero
  IF (z > zero) THEN
    IF (PRESENT(p)) p = one
    IF (PRESENT(q)) q = zero
  ELSE
    IF (PRESENT(p)) p = zero
    IF (PRESENT(q)) q = one
  END IF
  RETURN
END IF

! |Z| <= 37.

expntl = EXP(-0.5D0*zabs**2)
ppdf = expntl/root2pi
IF (PRESENT(pdf)) pdf = ppdf

! |Z| < CUTOFF = 10/sqrt(2).

IF (zabs < cutoff) THEN
  pp = expntl*((((((p6*zabs + p5)*zabs + p4)*zabs + p3)*zabs + p2)*zabs     &
                   + p1)*zabs + p0) / (((((((q7*zabs + q6)*zabs + q5)*zabs &
                   + q4)*zabs + q3)*zabs + q2)*zabs + q1)*zabs +q0)

! |Z| >= CUTOFF.

ELSE
  pp = ppdf/(zabs + one/(zabs + 2.d0/(zabs + 3.d0/(zabs + 4.d0/(zabs + 0.65D0)))))
END IF

IF (z < zero) THEN
  qq = one - pp
ELSE
  qq = pp
  pp = one - qq
END IF

IF (PRESENT(p)) p = pp
IF (PRESENT(q)) q = qq

RETURN
END SUBROUTINE normp



FUNCTION phi(z) RESULT(p)
REAL (dp), INTENT(IN) :: z
REAL (dp)             :: p

CALL normp(z, p)

RETURN
END FUNCTION phi

END FUNCTION bvn
