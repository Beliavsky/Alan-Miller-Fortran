MODULE NSWC_spec_func
USE constants_NSWC
IMPLICIT NONE

INTERFACE anorm
  MODULE PROCEDURE anorm
  MODULE PROCEDURE dnorm
  MODULE PROCEDURE znorm
END INTERFACE

CONTAINS


FUNCTION anorm(x, y) RESULT(fn_val)
! Replaces the statement function anorm in the F77 code.

REAL, INTENT(IN) :: x, y
REAL             :: fn_val

fn_val = MAX(ABS(x), ABS(y))
RETURN
END FUNCTION anorm


FUNCTION dnorm(x, y) RESULT(fn_val)
! Replaces the statement function anorm in the F77 code.

REAL (dp), INTENT(IN) :: x, y
REAL (dp)             :: fn_val

fn_val = MAX(ABS(x), ABS(y))
RETURN
END FUNCTION dnorm


FUNCTION znorm(z) RESULT(fn_val)
! Replaces the statement function anorm in the F77 code.

COMPLEX, INTENT(IN) :: z
REAL                :: fn_val

fn_val = MAX( ABS( REAL(z)), ABS(AIMAG(z) ) )
RETURN
END FUNCTION znorm


FUNCTION cpabs(x, y) RESULT(fn_val)
!     --------------------------------------
!     EVALUATION OF SQRT(X*X + Y*Y)
!     --------------------------------------
REAL, INTENT(IN) :: x, y
REAL             :: fn_val

! Local variable
REAL :: a

IF (ABS(x) > ABS(y)) THEN
  a = y / x
  fn_val = ABS(x) * SQRT(1.0 + a*a)
  RETURN
END IF
IF (y /= 0.0) THEN
  a = x / y
  fn_val = ABS(y) * SQRT(1.0 + a*a)
  RETURN
END IF
fn_val = 0.0
RETURN
END FUNCTION cpabs


FUNCTION dcpabs(x, y) RESULT(fn_val)
!     --------------------------------------
!     EVALUATION OF SQRT(X*X + Y*Y)
!     --------------------------------------
REAL (dp), INTENT(IN) :: x, y
REAL (dp)             :: fn_val

! Local variable
REAL (dp) :: a

IF (ABS(x) > ABS(y)) THEN
  a = y / x
  fn_val = ABS(x) * SQRT(1._dp+a*a)
  RETURN
END IF
IF (y /= 0._dp) THEN
  a = x / y
  fn_val = ABS(y) * SQRT(1._dp+a*a)
  RETURN
END IF
fn_val = 0._dp
RETURN
END FUNCTION dcpabs


SUBROUTINE crec(x, y, u, v)
!-----------------------------------------------------------------------
!             COMPLEX RECIPROCAL U + I*V = 1/(X + I*Y)
!-----------------------------------------------------------------------
REAL, INTENT(IN)  :: x, y
REAL, INTENT(OUT) :: u, v

! Local variables
REAL :: t, d

IF (ABS(x) <= ABS(y)) THEN
  t = x / y
  d = y + t * x
  u = t / d
  v = -1.0 / d
  RETURN
END IF
t = y / x
d = x + t * y
u = 1.0 / d
v = -t / d
RETURN
END SUBROUTINE crec


SUBROUTINE dcrec(x, y, u, v)
!-----------------------------------------------------------------------
!             COMPLEX RECIPROCAL U + I*V = 1/(X + I*Y)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x, y
REAL (dp), INTENT(OUT) :: u, v

! Local variables
REAL (dp) :: t, d

IF (ABS(x) <= ABS(y)) THEN
  t = x / y
  d = y + t * x
  u = t / d
  v = -1._dp / d
  RETURN
END IF
t = y / x
d = x + t * y
u = 1._dp / d
v = -t / d
RETURN
END SUBROUTINE dcrec


FUNCTION cdiv(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------
!              COMPLEX DIVISION A/B WHERE B IS NONZERO
!-----------------------------------------------------------------------
COMPLEX             :: fn_val
COMPLEX, INTENT(IN) :: a, b

! Local variables
REAL :: ar, ai, br, bi, t, d, u, v

ar = REAL(a)
ai = AIMAG(a)
br = REAL(b)
bi = AIMAG(b)

IF (ABS(br) >= ABS(bi)) THEN
  t = bi / br
  d = br + t * bi
  u = (ar+ai*t) / d
  v = (ai-ar*t) / d
  fn_val = CMPLX(u,v)
  RETURN
END IF
t = br / bi
d = bi + t * br
u = (ar*t+ai) / d
v = (ai*t-ar) / d
fn_val = CMPLX(u,v)
RETURN
END FUNCTION cdiv


SUBROUTINE cdivid(ar, ai, br, bi, cr, ci)
!-----------------------------------------------------------------------
!     REAL (dp) COMPLEX DIVISION C = A/B AVOIDING OVERFLOW
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: ar, ai, br, bi
REAL (dp), INTENT(OUT) :: cr, ci

! Local variables
REAL (dp) :: d, t, u, v

IF (ABS(br) > ABS(bi)) THEN
  t = bi / br
  d = br + t * bi
  u = (ar+ai*t) / d
  v = (ai-ar*t) / d
  cr = u
  ci = v
  RETURN
END IF

IF (bi /= 0._dp) THEN
  t = br / bi
  d = bi + t * br
  u = (ar*t+ai) / d
  v = (ai*t-ar) / d
  cr = u
  ci = v
  RETURN
END IF

!     DIVISION BY ZERO. C = INFINITY

cr = HUGE(1.0_dp)
ci = cr
RETURN
END SUBROUTINE cdivid


SUBROUTINE capo(x, y, r, theta)
REAL, INTENT(IN)  :: x, y
REAL, INTENT(OUT) :: r, theta

! Local variable
REAL :: a

IF (ABS(x) > ABS(y)) THEN
  a = y / x
  r = ABS(x) * SQRT(1.0+a*a)
  theta = ATAN2(y,x)
  RETURN
END IF
IF (y /= 0.) THEN
  a = x / y
  r = ABS(y) * SQRT(1.0+a*a)
  theta = ATAN2(y,x)
  RETURN
END IF
r = 0.0
theta = 0.0
RETURN
END SUBROUTINE capo


FUNCTION sin0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!          COMPUTATION OF SIN(X*PI/2) FOR ABS(X) <= 0.5
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL            :: t
REAL, PARAMETER :: a0 = .157079632679490E+01, a1 = -.645964097506244E+00,   &
                   a2 = .796926262460396E-01, a3 = -.468175413228242E-02,   &
                   a4 = .160441150291651E-03, a5 = -.359864175444606E-05,   &
                   a6 = .563372101191893E-07
!------------------------
t = x * x
fn_val = ((((((a6*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + a0) * x
RETURN
END FUNCTION sin0


FUNCTION cos0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!          COMPUTATION OF COS(X*PI/2) FOR ABS(X) <= 0.5
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL            :: t
REAL, PARAMETER :: a1 = -.123370055013615E+01, a2 = .253669507899753E+00,   &
                   a3 = -.208634807330586E-01, a4 = .919259935580283E-03,   &
                   a5 = -.252000841382533E-04, a6 = .465461768260405E-06
!------------------------
t = x * x
fn_val = (((((a6*t + a5)*t + a4)*t + a3)*t + a2)*t + a1) * t + 1.0
RETURN
END FUNCTION cos0


FUNCTION sin1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!                   EVALUATION  OF SIN(X*PI)
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: a0 = .314159265358979E+01, a1 = -.516771278004995E+01,   &
                   a2 = .255016403987327E+01, a3 = -.599264528932149E+00,   &
                   a4 = .821458689493251E-01, a5 = -.737001831310553E-02,   &
                   a6 = .461514425296398E-03, b1 = -.493480220054460E+01,   &
                   b2 = .405871212639605E+01, b3 = -.133526276691575E+01,   &
                   b4 = .235330543508553E+00, b5 = -.258048861575714E-01,   &
                   b6 = .190653140279462E-02
INTEGER         :: max, n
REAL            :: a, t
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                       MAX = IPMPAR(3)
max = HUGE(3)

!------------------------
a = ABS(x)
IF (a >= REAL(MAX)) THEN
  fn_val = 0.0
  RETURN
END IF

n = a
a = a - REAL(n)
IF (a <= 0.75) THEN
  IF (a < 0.25) GO TO 10

!                    0.25 <= A <= 0.75

  a = 0.25 + (0.25-a)
  t = a * a
  fn_val = ((((((b6*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t + 0.5) + 0.5
  GO TO 20
END IF

!                 A < 0.25  OR  A > 0.75

a = 0.25 + (0.75-a)
10 t = a * a
fn_val = ((((((a6*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + a0) * a

!                        TERMINATION

20 IF (x < 0.0) fn_val = -fn_val
IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION sin1


FUNCTION cos1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!                   EVALUATION  OF COS(X*PI)
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: a0 = .314159265358979E+01, a1 = -.516771278004995E+01,   &
                    a2 = .255016403987327E+01, a3 = -.599264528932149E+00,   &
                    a4 = .821458689493251E-01, a5 = -.737001831310553E-02,   &
                    a6 = .461514425296398E-03, b1 = -.493480220054460E+01,   &
                    b2 = .405871212639605E+01, b3 = -.133526276691575E+01,   &
                    b4 = .235330543508553E+00, b5 = -.258048861575714E-01,   &
                    b6 = .190653140279462E-02
INTEGER          :: max, n
REAL             :: a, t
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                       MAX = IPMPAR(3)
max = HUGE(3)

!------------------------
a = ABS(x)
IF (a >= REAL(MAX)) THEN
  fn_val = 1.0
  RETURN
END IF

n = a
a = a - REAL(n)
IF (a <= 0.75) THEN
  IF (a < 0.25) GO TO 10

!                    0.25 <= A <= 0.75

  a = 0.25 + (0.25-a)
  t = a * a
  fn_val = ((((((a6*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + a0) * a
  GO TO 20
END IF

!                 A < 0.25  OR  A > 0.75

a = 0.25 + (0.75-a)
n = n - 1
10 t = a * a
fn_val = ((((((b6*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t + 0.5) + 0.5

!                        TERMINATION

20 IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION cos1


FUNCTION dsin1(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!                REAL (dp) EVALUATION OF SIN(PI*X)

!                             --------------

!     THE EXPANSION FOR SIN(PI*A) (ABS(A) <= PI/4) USING A1,...,A13
!     IS ACCURATE TO WITHIN 2 UNITS OF THE 40-TH SIGNIFICANT DIGIT, AND
!     THE EXPANSION FOR COS(PI*A) (ABS(A) <= PI/4) USING B1,...,B13
!     IS ACCURATE TO WITHIN 4 UNITS OF THE 40-TH SIGNIFICANT DIGIT.

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp)            :: a, t, w
REAL (dp), PARAMETER :: pi = 3.141592653589793238462643383279502884197D+00
REAL (dp), PARAMETER :: a1 = -.1028083791780141522795259479153765743002D+00,   &
      a2  = .3170868848763100170457042079710451905600D-02,   &
      a3  = -.4657026956105571623449026167864697920000D-04,  &
      a4  = .3989844942879455643410226655783424000000D-06,   &
      a5  = -.2237397227721999776371894030796800000000D-08,  &
      a6  = .8847045483056962709715066675200000000000D-11,   &
      a7  = -.2598715447506450292885585920000000000000D-13,  &
      a8  = .5893449774331011070033920000000000000000D-16 ,  &
      a9  = -.1062975472045522550784000000000000000000D-18,   &
      a10 = .1561182648301780992000000000000000000000D-21,    &
      a11 = -.1903193516670976000000000000000000000000D-24,   &
      a12 = .1956617650176000000000000000000000000000D-27,    &
      a13 = -.1711276032000000000000000000000000000000D-30
REAL (dp), PARAMETER :: b1 = -.3084251375340424568385778437461297229882D+00, &
      b2  = .1585434424381550085228521039855226435920D-01,   &
      b3  = -.3259918869273900136414318317506279360000D-03,  &
      b4  = .3590860448591510079069203991239232000000D-05,   &
      b5  = -.2461136950494199754009084061808640000000D-07,  &
      b6  = .1150115912797405152263195572224000000000D-09,   &
      b7  = -.3898073171259675439899172864000000000000D-12,  &
      b8  = .1001886461636271969091584000000000000000D-14,   &
      b9  = -.2019653396886572027084800000000000000000D-17,  &
      b10 = .3278483561466560512000000000000000000000D-20,   &
      b11 = -.4377345082051788800000000000000000000000D-23,  &
      b12 = .4891532381388800000000000000000000000000D-26,   &
      b13 = -.4617089843200000000000000000000000000000D-29
INTEGER :: max, n
!------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                       MAX = IPMPAR(3)
max = HUGE(3)

!------------------------
a = ABS(x)
t = MAX
IF (a >= t) THEN
  fn_val = 0._dp
  RETURN
END IF

n = a
t = n
a = a - t
IF (a <= 0.75_dp) THEN
  IF (a < 0.25_dp) GO TO 10

!                    0.25 <= A <= 0.75

  a = 0.25_dp + (0.25_dp-a)
  t = 16._dp * a * a
  fn_val = (((((((((((((b13*t + b12)*t + b11)*t + b10)*t + b9)*t + b8)*t  &
           + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t +  &
           0.5_dp) + 0.5_dp
  GO TO 20
END IF

!                 A < 0.25  OR  A > 0.75

a = 0.25_dp + (0.75_dp-a)
10 t = 16._dp * a * a
w = (((((((((((((a13*t + a12)*t + a11)*t + a10)*t + a9)*t + a8)*t + a7)*t  &
    + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t + 0.5_dp) + 0.5_dp
fn_val = pi * a * w

!                        TERMINATION

20 IF (x < 0.0) fn_val = -fn_val
IF (MOD(n,2) /= 0) fn_val = -fn_val
RETURN
END FUNCTION dsin1


SUBROUTINE snhcsh(x, isw, sinhm, coshm)

REAL, INTENT(IN)            :: x
INTEGER, INTENT(IN)         :: isw
REAL, INTENT(OUT), OPTIONAL :: sinhm, coshm

! Local variables
REAL            :: ax, expx, xs, xx
REAL, PARAMETER :: cut(5) = (/ 1.65, 1.2, 1.2, 2.7, 1.65 /),  &
        sp5 =  .255251817302048E-09, sp4 =  .723809046696880E-07,  &
        sp3 =  .109233297700241E-04, sp2 =  .954811583154274E-03,  &
        sp1 =  .452867078563929E-01, sq1 = -.471329214363072E-02,  &
        cp5 =  .116744361560051E-08, cp4 =  .280407224259429E-06,  &
        cp3 =  .344417983443219E-04, cp2 =  .232293648552398E-02,  &
        cp1 =  .778752378267155E-01, cq1 = -.545809550662099E-02,  &
        zp3 =  5.59297116264720E-07, zp2 =  1.77943488030894E-04,  &
        zp1 =  1.69800461894792E-02, zq4 =  1.33412535492375E-09,  &
        zq3 = -5.80858944138663E-07, zq2 =  1.27814964403863E-04,  &
        zq1 = -1.63532871439181E-02

!                      FROM THE SPLINE UNDER TENSION PACKAGE
!                       CODED BY A. K. CLINE AND R. J. RENKA
!                            DEPARTMENT OF COMPUTER SCIENCES
!                              UNIVERSITY OF TEXAS AT AUSTIN
!                          MODIFIED BY A.H. MORRIS (NSWC/DL)

! THIS SUBROUTINE RETURNS APPROXIMATIONS TO
!       SINHM(X) = SINH(X) - X
!       COSHM(X) = COSH(X) - 1
! AND
!       COSHMM(X) = COSH(X) - 1 - X*X/2

! ON INPUT--

!   X CONTAINS THE VALUE OF THE INDEPENDENT VARIABLE.

!   ISW INDICATES THE FUNCTION DESIRED
!           = -1 IF ONLY SINHM IS DESIRED,
!           =  0 IF BOTH SINHM AND COSHM ARE DESIRED,
!           =  1 IF ONLY COSHM IS DESIRED,
!           =  2 IF ONLY COSHMM IS DESIRED,
!           =  3 IF BOTH SINHM AND COSHMM ARE DESIRED.

! ON OUTPUT--

!   SINHM CONTAINS THE VALUE OF SINHM(X) IF ISW <= 0 OR ISW = 3
!   (SINHM IS UNALTERED IF ISW =1 OR ISW = 2).

!   COSHM CONTAINS THE VALUE OF COSHM(X) IF ISW = 0 OR ISW = 1 AND CONTAINS
!   THE VALUE OF COSHMM(X) IF ISW >= 2 (COSHM IS UNALTERED IF ISW = -1).

! AND

!   X AND ISW ARE UNALTERED.

!-----------------------------------------------------------


xx = x
ax = ABS(xx)
xs = xx * xx
IF (ax >= cut(isw+2)) expx = EXP(ax)

! SINHM APPROXIMATION

IF (isw /= 1 .AND. isw /= 2) THEN
  IF (ax < 1.65) THEN
    sinhm = ((((((sp5*xs + sp4)*xs + sp3)*xs + sp2)*xs + sp1)*xs + 1.)*xs*xx)/ &
            ((sq1*xs + 1.)*6.)
  ELSE
    sinhm = -(((ax+ax) + 1./expx) - expx) / 2.
    IF (xx < 0.) sinhm = -sinhm
  END IF
END IF

! COSHM APPROXIMATION

IF (isw == 0.OR.isw == 1) THEN
  IF (ax < 1.2) THEN
    coshm = ((((((cp5*xs + cp4)*xs + cp3)*xs + cp2)*xs + cp1)*xs + 1.)*xs) /  &
            ((cq1*xs + 1.)*2.)
  ELSE
    coshm = ((1./expx - 2.) + expx) / 2.
  END IF
END IF

! COSHMM APPROXIMATION

IF (isw <= 1) RETURN
IF (ax < 2.70) THEN
  coshm = ((((zp3*xs + zp2)*xs + zp1)*xs + 1.)*xs*xs) /  &
          (((((zq4*xs + zq3)*xs + zq2)*xs + zq1)*xs + 1.)*24.)
  RETURN
END IF
coshm = (((1./expx - 2.) - xs) + expx) / 2.
RETURN
END SUBROUTINE snhcsh


FUNCTION esum(mu, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!                    EVALUATION OF EXP(MU + X)
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: mu, x
REAL             :: fn_val

! Local variable
REAL :: w

IF (x <= 0.0) THEN

  IF (mu < 0) GO TO 10
  w = mu + x
  IF (w > 0.0) GO TO 10
  fn_val = EXP(w)
  RETURN
END IF

IF (mu <= 0) THEN
  w = mu + x
  IF (w >= 0.0) THEN
    fn_val = EXP(w)
    RETURN
  END IF
END IF

10 w = mu
fn_val = EXP(w) * EXP(x)
RETURN
END FUNCTION esum


FUNCTION rexp(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: p1 =  .914041914819518E-09, p2 = .238082361044469E-01,  &
                    q1 = -.499999999085958E+00, q2 = .107141568980644E+00,  &
                    q3 = -.119041179760821E-01, q4 = .595130811860248E-03
REAL :: e
!-----------------------
IF (ABS(x) <= 0.15) THEN
  fn_val = x * (((p2*x + p1)*x + 1.0) / ((((q4*x+q3)*x+q2)*x+q1)*x + 1.0))
  RETURN
END IF

IF (x >= 0.0) THEN
  e = EXP(x)
  fn_val = e * (0.5 + (0.5 - 1.0/e))
  RETURN
END IF
IF (x >= -37.0) THEN
  fn_val = (EXP(x)-0.5) - 0.5
  RETURN
END IF
fn_val = -1.0
RETURN
END FUNCTION rexp


FUNCTION drexp(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: e, w, z
REAL (dp) :: a0 = .248015873015873015873016D-04,   &
    a1 = -.344452080605731005808147D-05, a2 = .206664230430046597475413D-06,  &
    a3 = -.447300111094328162971036D-08, a4 = .114734027080634968083920D-11,  &
    b1 = -.249994190011341852652396D+00, b2 = .249987228833107957725728D-01,  &
    b3 = -.119037506846942249362528D-02, b4 = .228908693387350391768682D-04
REAL (dp) :: c1 = .1666666666666666666666666666666667D+00,   &
             c2 = .4166666666666666666666666666666667D-01,   &
             c3 = .8333333333333333333333333333333333D-02,   &
             c4 = .1388888888888888888888888888888889D-02,   &
             c5 = .1984126984126984126984126984126984D-03
!---------------------------
IF (ABS(x) <= 0.15_dp) THEN

!        Z IS A MINIMAX APPROXIMATION OF THE SERIES

!                C6 + C7*X + C8*X**2 + ....

!        THIS APPROXIMATION IS ACCURATE TO WITHIN
!        1 UNIT OF THE 23-RD SIGNIFICANT DIGIT.
!        THE RESULTING VALUE FOR W IS ACCURATE TO
!        WITHIN 1 UNIT OF THE 33-RD SIGNIFICANT
!        DIGIT.

  z = ((((a4*x + a3)*x + a2)*x + a1)*x + a0) /  &
      ((((b4*x + b3)*x + b2)*x + b1)*x + 1._dp)
  w = ((((((z*x + c5)*x + c4)*x + c3)*x + c2)*x + c1)*x + 0.5_dp)*x + 1._dp
  fn_val = x * w
  RETURN
END IF

IF (x >= 0._dp) THEN
  e = EXP(x)
  fn_val = e * (0.5_dp+(0.5_dp-1._dp/e))
  RETURN
END IF
IF (x >= -77._dp) THEN
  fn_val = (EXP(x)-0.5_dp) - 0.5_dp
  RETURN
END IF
fn_val = -1._dp
RETURN
END FUNCTION drexp


FUNCTION alnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: p1 = -.129418923021993E+01, p2 = .405303492862024E+00,  &
        p3 = -.178874546012214E-01, q1 = -.162752256355323E+01,   &
        q2 =  .747811014037616E+00, q3 = -.845104217945565E-01,  &
        zero = 0._dp, half = 0.5_dp, one = 1._dp, two = 2._dp
REAL :: t, t2, w, x
!--------------------------
IF (ABS(a) <= 0.375) THEN
  t = a / (a + two)
  t2 = t * t
  w = (((p3*t2 + p2)*t2 + p1)*t2 + one) / (((q3*t2 + q2)*t2 + q1)*t2 + one)
  fn_val = two * t * w
  RETURN
END IF

x = one + a
IF (a < zero) x = (a + half) + half
fn_val = LOG(x)
RETURN
END FUNCTION alnrel


FUNCTION dlnrel(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF THE FUNCTION LN(1 + A)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: t, t2, w, z
REAL (dp) :: p0 = .7692307692307692307680D-01,   &
       p1 = -.1505958055914600184836D+00, p2 =  .9302355725278521726,   &
       p3 = -.1787900022182327735804D-01, q1 = -.2824412139355646910683D+01,  &
       q2 =  .2892424216041495392509D+01, q3 = -.1263560605948009364422D+01,  &
       q4 =  .1966769435894561313526D+00
REAL (dp) :: c1 = .3333333333333333333333333333333D+00,   &
             c2 = .2000000000000000000000000000000D+00,   &
             c3 = .1428571428571428571428571428571D+00,   &
             c4 = .1111111111111111111111111111111D+00,   &
             c5 = .9090909090909090909090909090909D-01
!-------------------------
IF (ABS(a) >= 0.375_dp) THEN
  t = 1._dp + a
  IF (a < 0._dp) t = 0.5_dp + (0.5_dp + a)
  fn_val = LOG(t)
  RETURN
END IF

!        W IS A MINIMAX APPROXIMATION OF THE SERIES

!               C6 + C7*T**2 + C8*T**4 + ...

!        THIS APPROXIMATION IS ACCURATE TO WITHIN 1.6 UNITS OF THE 21-ST
!        SIGNIFICANT DIGIT.
!        THE RESULTING VALUE FOR 1._dp + T2*Z IS ACCURATE TO WITHIN 1 UNIT OF
!        THE 30-TH SIGNIFICANT DIGIT.

t = a / (a + 2._dp)
t2 = t * t
w = (((p3*t2 + p2)*t2 + p1)*t2 + p0) /  &
    ((((q4*t2 + q3)*t2 + q2)*t2 + q1)*t2 + 1._dp)
z = ((((w*t2 + c5)*t2 + c4)*t2 + c3)*t2 + c2)*t2 + c1
fn_val = 2._dp * t * (1._dp + t2*z)
RETURN
END FUNCTION dlnrel


FUNCTION rlog(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE FUNCTION X - 1 - LN(X)
!-----------------------------------------------------------------------
!     A = RLOG (0.7)
!     B = RLOG (4/3)
!------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: a = .566749439387324E-01, b = .456512608815524E-01,   &
                   p0 = .333333333333333E+00, p1 = -.224696413112536E+00, &
                   p2 = .620886815375787E-02, q1 = -.127408923933623E+01, &
                   q2 = .354508718369557E+00
REAL :: r, t, u, up2, w, w1
!------------------------
IF (x >= 0.61 .AND. x <= 1.57) THEN
  IF (x >= 0.82) THEN
    IF (x > 1.18) GO TO 10

!                 ARGUMENT REDUCTION

    u = (x-0.5) - 0.5
    up2 = u + 2.0
    w1 = 0.0
    GO TO 20
  END IF

  u = (x-0.7) / 0.7
  up2 = u + 2.0
  w1 = a - u * 0.3
  GO TO 20

  10 t = 0.75 * (x-1.0)
  u = t - 0.25
  up2 = t + 1.75
  w1 = b + u / 3.0

!                  SERIES EXPANSION

  20 r = u / up2
  t = r * r
  w = ((p2*t + p1)*t + p0) / ((q2*t + q1)*t + 1.0)
  fn_val = r * (u - 2.0*t*w) + w1
  RETURN
END IF


r = (x-0.5) - 0.5
fn_val = r - LOG(x)
RETURN
END FUNCTION rlog


FUNCTION rlog1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE FUNCTION X - LN(1 + X)
!-----------------------------------------------------------------------
!     A = RLOG (0.7)
!     B = RLOG (4/3)
!------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: a = .566749439387324E-01, b = .456512608815524E-01,   &
                   p0 = .333333333333333E+00, p1 = -.224696413112536E+00, &
                   p2 = .620886815375787E-02, q1 = -.127408923933623E+01, &
                   q2 = .354508718369557E+00
REAL :: r, t, u, up2, w, w1
!------------------------
IF (x >= -0.39 .AND. x <= 0.57) THEN
  IF (x >= -0.18) THEN
    IF (x > 0.18) GO TO 10

!                 ARGUMENT REDUCTION

    u = x
    up2 = u + 2.0
    w1 = 0.0
    GO TO 20
  END IF

  u = (x+0.3) / 0.7
  up2 = u + 2.0
  w1 = a - u * 0.3
  GO TO 20

  10 t = 0.75 * x
  u = t - 0.25
  up2 = t + 1.75
  w1 = b + u / 3.0

!                  SERIES EXPANSION

  20 r = u / up2
  t = r * r
  w = ((p2*t + p1)*t + p0) / ((q2*t + q1)*t + 1.0)
  fn_val = r * (u - 2.0*t*w) + w1
  RETURN
END IF


w = (x+0.5) + 0.5
fn_val = x - LOG(w)
RETURN
END FUNCTION rlog1


FUNCTION drlog(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE FUNCTION X - 1 - LN(X)
!-----------------------------------------------------------------------

REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: a = .566749439387323789126387112411845D-01,   &
                        b = .456512608815524058941143273395059D-01
REAL (dp), PARAMETER :: p0 = .7692307692307692307680D-01,   &
       p1 = -.1505958055914600184836D+00, p2 = .9302355725278521726,   &
       p3 = -.1787900022182327735804D-01, q1 = -.2824412139355646910683D+01, &
       q2 = .2892424216041495392509D+01, q3 = -.1263560605948009364422D+01,  &
       q4 = .1966769435894561313526D+00
REAL (dp), PARAMETER :: c1 = .333333333333333333333333333333333D+00,   &
                        c2 = .200000000000000000000000000000000D+00,   &
                        c3 = .142857142857142857142857142857143D+00,   &
                        c4 = .111111111111111111111111111111111D+00,   &
                        c5 = .909090909090909090909090909090909D-01
REAL (dp) :: r, t, u, up2, w, w1, z
!-------------------------
!     A = DRLOG (0.7)
!     B = DRLOG (4/3)
!-------------------------
IF (x >= 0.61_dp .AND. x <= 1.57_dp) THEN
  IF (x >= 0.82_dp) THEN
    IF (x > 1.18_dp) GO TO 10

!                 ARGUMENT REDUCTION

    u = (x-0.5_dp) - 0.5_dp
    up2 = u + 2._dp
    w1 = 0._dp
    GO TO 20
  END IF

  u = (x-0.7_dp) / 0.7_dp
  up2 = u + 2._dp
  w1 = a - u * 0.3_dp
  GO TO 20

  10 t = 0.75_dp * (x-1._dp)
  u = t - 0.25_dp
  up2 = t + 1.75_dp
  w1 = b + u / 3._dp

!                  SERIES EXPANSION

  20 r = u / up2
  t = r * r

!        Z IS A MINIMAX APPROXIMATION OF THE SERIES

!               C6 + C7*R**2 + C8*R**4 + ...

!        FOR THE INTERVAL (0.0, 0.375). THE APPROX-
!        IMATION IS ACCURATE TO WITHIN 1.6 UNITS OF
!        THE 21-ST SIGNIFICANT DIGIT.

  z = (((p3*t + p2)*t + p1)*t + p0) / ((((q4*t+q3)*t+q2)*t+q1)*t+1._dp)

  w = ((((z*t + c5)*t + c4)*t + c3)*t + c2) * t + c1
  fn_val = r * (u-2._dp*t*w) + w1
  RETURN
END IF


r = (x-0.5_dp) - 0.5_dp
fn_val = r - LOG(x)
RETURN
END FUNCTION drlog


SUBROUTINE cerf(mo, z, w)
!-----------------------------------------------------------------------

!               COMPUTATION OF THE COMPLEX ERROR FUNCTION

!                           ----------------

!                        W = ERF(Z)   IF MO = 0
!                        W = ERFC(Z)  OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)  :: mo
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: cd(18) = (/ 0.000000000000000, 2.08605856013476E-2,  &
        8.29806940495687E-2, 1.85421653326079E-1, 3.27963479382361E-1,  &
        5.12675279912828E-1, 7.45412958045105E-1, 1.03695067418297,  &
        1.40378061255437, 1.86891662214001, 2.46314830523929,  &
        3.22719383737352, 4.21534348280013, 5.50178873151549,  &
        7.19258966683102, 9.45170208076408, 1.25710718314784E+1,  &
        1.72483537216334E+1 /), ce(18) = (/ 8.15723083324096E-2,  &
        1.59285285253437E-1, 1.48581625614499E-1, 1.33219670836245E-1,  &
        1.15690392878957E-1, 9.78580959447535E-2, 8.05908834297624E-2,  &
        6.40204538609872E-2, 4.81445242767885E-2, 3.33540658473295E-2,  &
        2.05548099470193E-2, 1.07847403887506E-2, 4.55634892214219E-3,  &
        1.43984458138925E-3, 3.07056139834171E-4, 3.78156541168541E-5,  &
        2.05173509616121E-6, 2.63564823682747E-8 /), c = .564189583547756
REAL    :: ef(2), qf(2), sm(2), sz(2), tm(2), ts(2), c2, dm, pm, qm, r,  &
           sn, ss, x, y
INTEGER :: i
!------------------------
!     C = 1/SQRT(PI)
!------------------------
x = REAL(z)
y = AIMAG(z)
sn = 1.0
IF (x < 0.0) THEN
  x = -x
  y = -y
  sn = -1.0
END IF

r = x * x + y * y
sz(1) = x * x - y * y
sz(2) = 2.0 * x * y

IF (r > 1.0) THEN
  IF (r >= 38.0) GO TO 50
  IF (sz(1) + 0.064*sz(2)*sz(2) > 0.0) GO TO 30
END IF

!                       TAYLOR SERIES

c2 = c + c
tm(1) = c2 * x
tm(2) = c2 * y
sm(1) = tm(1)
sm(2) = tm(2)
pm = 0.0
10 pm = pm + 1.0
dm = 2.0 * pm + 1.0
ts(1) = tm(1) * sz(1) - tm(2) * sz(2)
ts(2) = tm(1) * sz(2) + tm(2) * sz(1)
tm(1) = -ts(1) / pm
tm(2) = -ts(2) / pm
ts(1) = tm(1) / dm
ts(2) = tm(2) / dm
IF (ABS(sm(1))+ABS(ts(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(ts(2)) == ABS(sm(2))) GO TO 20
END IF
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
GO TO 10

!                       TERMINATION

20 IF (mo == 0) THEN
  w = CMPLX(sn*sm(1),sn*sm(2))
  RETURN
END IF
IF (sn >= 0.0) THEN
  sm(1) = 0.5 + (0.5-sm(1))
  sm(2) = -sm(2)
  w = CMPLX(sm(1),sm(2))
  RETURN
END IF
w = CMPLX(1.0+sm(1),sm(2))
RETURN

!              RATIONAL FUNCTION APPROXIMATION

30 sm(1) = 0.0
sm(2) = 0.0
qm = c * EXP(-sz(1))
ts(1) = qm * COS(-sz(2))
ts(2) = qm * SIN(-sz(2))
qf(1) = ts(1) * x - ts(2) * y
qf(2) = ts(1) * y + ts(2) * x
DO i = 1, 18
  ts(1) = sz(1) + cd(i)
  ts(2) = sz(2)
  ss = ts(1) * ts(1) + ts(2) * ts(2)
  tm(1) = ce(i) * ts(1) / ss
  tm(2) = -ce(i) * ts(2) / ss
  sm(1) = sm(1) + tm(1)
  sm(2) = sm(2) + tm(2)
END DO
ef(1) = qf(1) * sm(1) - qf(2) * sm(2)
ef(2) = qf(1) * sm(2) + qf(2) * sm(1)
GO TO 80

!                   ASYMPTOTIC EXPANSION

50 qf(1) = sz(1) / (r*r)
qf(2) = -sz(2) / (r*r)
qm = c * EXP(-sz(1))
ts(1) = qm * COS(-sz(2))
ts(2) = qm * SIN(-sz(2))
tm(1) = (ts(1)*x+ts(2)*y) / r
tm(2) = -(ts(1)*y-ts(2)*x) / r
sm(1) = tm(1)
sm(2) = tm(2)
pm = -0.5
60 pm = pm + 1.0
ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
tm(1) = -pm * ts(1)
tm(2) = -pm * ts(2)
IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 70
END IF
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (pm < 25.5) GO TO 60

70 IF (x < 0.01) THEN
  sn = -sn
  GO TO 20
END IF
ef(1) = sm(1)
ef(2) = sm(2)

!                       TERMINATION

80 IF (mo /= 0) THEN
  w = CMPLX(ef(1),ef(2))
  IF (sn == 1.0) RETURN
  w = CMPLX(2.0-ef(1),-ef(2))
  RETURN
END IF
ef(1) = sn * (1.0-ef(1))
ef(2) = -sn * ef(2)
w = CMPLX(ef(1),ef(2))
RETURN
END SUBROUTINE cerf


SUBROUTINE cerfc(mo, z, w)
!-----------------------------------------------------------------------

!              COMPUTATION OF THE COMPLEX COERROR FUNCTION

!                           ----------------

!           W = ERFC(Z)           IF MO = 0 OR REAL(Z) < 0
!           W = EXP(X*X)*ERFC(Z)  OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)  :: mo
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: cd(18) = (/ 0.00000000000000, 2.08605856013476E-2,  &
        8.29806940495687E-2, 1.85421653326079E-1, 3.27963479382361E-1,  &
        5.12675279912828E-1, 7.45412958045105E-1, 1.03695067418297,  &
        1.40378061255437, 1.86891662214001, 2.46314830523929,  &
        3.22719383737352, 4.21534348280013, 5.50178873151549,  &
        7.19258966683102, 9.45170208076408, 1.25710718314784E+1,  &
        1.72483537216334E+1 /), ce(18) = (/ 8.15723083324096E-2,  &
        1.59285285253437E-1, 1.48581625614499E-1, 1.33219670836245E-1,  &
        1.15690392878957E-1, 9.78580959447535E-2, 8.05908834297624E-2,  &
        6.40204538609872E-2, 4.81445242767885E-2, 3.33540658473295E-2,  &
        2.05548099470193E-2, 1.07847403887506E-2, 4.55634892214219E-3,  &
        1.43984458138925E-3, 3.07056139834171E-4, 3.78156541168541E-5,  &
        2.05173509616121E-6, 2.63564823682747E-8 /), c = .564189583547756
REAL    :: qf(2), sm(2), sz(2), tm(2), ts(2), c2, dm, pm, qm, r, sn, ss, x, y
INTEGER :: i
!------------------------
!     C = 1/SQRT(PI)
!------------------------
x = REAL(z)
y = AIMAG(z)
sn = 1.0
IF (x < 0.0) THEN
  x = -x
  y = -y
  sn = -1.0
END IF

IF (mo == 0.OR.sn /= 1.0.OR.MAX(x,ABS(y)) < 100.0) THEN
  r = x * x + y * y
  sz(1) = x * x - y * y
  sz(2) = 2.0 * x * y

  IF (r > 1.0) THEN
    IF (r >= 38.0) GO TO 50
    IF (sz(1) + 0.064*sz(2)*sz(2) > 0.0) GO TO 30
  END IF

!                       TAYLOR SERIES

  c2 = c + c
  tm(1) = c2 * x
  tm(2) = c2 * y
  sm(1) = tm(1)
  sm(2) = tm(2)
  pm = 0.0
  10 pm = pm + 1.0
  dm = 2.0 * pm + 1.0
  ts(1) = tm(1) * sz(1) - tm(2) * sz(2)
  ts(2) = tm(1) * sz(2) + tm(2) * sz(1)
  tm(1) = -ts(1) / pm
  tm(2) = -ts(2) / pm
  ts(1) = tm(1) / dm
  ts(2) = tm(2) / dm
  IF (ABS(sm(1))+ABS(ts(1)) == ABS(sm(1))) THEN
    IF (ABS(sm(2))+ABS(ts(2)) == ABS(sm(2))) GO TO 20
  END IF
  sm(1) = sm(1) + ts(1)
  sm(2) = sm(2) + ts(2)
  GO TO 10

!                       TERMINATION

  20 IF (sn /= 1.0) THEN
    w = CMPLX(1.0+sm(1),sm(2))
    RETURN
  END IF
  sm(1) = 0.5 + (0.5-sm(1))
  sm(2) = -sm(2)
  IF (mo == 0) GO TO 90

  qm = EXP(sz(1))
  qf(1) = qm * COS(sz(2))
  qf(2) = qm * SIN(sz(2))
  ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
  ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
  w = CMPLX(ts(1),ts(2))
  RETURN

!              RATIONAL FUNCTION APPROXIMATION

  30 sm(1) = 0.0
  sm(2) = 0.0
  DO i = 1, 18
    ts(1) = sz(1) + cd(i)
    ts(2) = sz(2)
    ss = ts(1) * ts(1) + ts(2) * ts(2)
    tm(1) = ce(i) * ts(1) / ss
    tm(2) = -ce(i) * ts(2) / ss
    sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
  END DO
  ts(1) = x * sm(1) - y * sm(2)
  ts(2) = x * sm(2) + y * sm(1)
  sm(1) = c * ts(1)
  sm(2) = c * ts(2)
  GO TO 80
END IF

!                   ASYMPTOTIC EXPANSION

50 CALL crec(x,y,tm(1),tm(2))
sm(1) = tm(1)
sm(2) = tm(2)
qf(1) = tm(1) * tm(1) - tm(2) * tm(2)
qf(2) = 2.0 * tm(1) * tm(2)
pm = -0.5
60 pm = pm + 1.0
ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
tm(1) = -pm * ts(1)
tm(2) = -pm * ts(2)
IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 70
END IF
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (pm < 25.5) GO TO 60
70 sm(1) = c * sm(1)
sm(2) = c * sm(2)
IF (x < 0.01) GO TO 100

!                       TERMINATION

80 IF (mo == 0 .OR. sn /= 1.0) THEN
  qm = EXP(-sz(1))
  qf(1) = qm * COS(-sz(2))
  qf(2) = qm * SIN(-sz(2))
  ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
  ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
  sm(1) = ts(1)
  sm(2) = ts(2)

  IF (sn /= 1.0) THEN
    w = CMPLX(2.0-sm(1),-sm(2))
    RETURN
  END IF
END IF
90 w = CMPLX(sm(1),sm(2))
RETURN

!               MODIFIED ASYMPTOTIC EXPANSION

100 IF (mo == 0.OR.sn /= 1.0) THEN
  qm = EXP(-sz(1))
  qf(1) = qm * COS(-sz(2))
  qf(2) = qm * SIN(-sz(2))
  ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
  ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
  sm(1) = 1.0 + sn * ts(1)
  sm(2) = sn * ts(2)
  w = CMPLX(sm(1),sm(2))
  RETURN
END IF

IF (ABS(y) >= 100.0) GO TO 90
IF (sz(1) <= exparg(1)) GO TO 90
qm = EXP(sz(1))
sm(1) = qm * COS(sz(2)) + sm(1)
sm(2) = qm * SIN(sz(2)) + sm(2)
w = CMPLX(sm(1),sm(2))
RETURN
END SUBROUTINE cerfc


FUNCTION erf (x) RESULT(fn_val)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL ERROR FUNCTION
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

REAL :: a(5) = (/ .771058495001320E-04, -.133733772997339E-02,   &
                  .323076579225834E-01,  .479137145607681E-01,   &
                  .128379167095513E+00 /),   &
        b(3) = (/ .301048631703895E-02,  .538971687740286E-01,   &
                  .375795757275549E+00 /),   &
        p(8) = (/-1.36864857382717E-07,  5.64195517478974E-01, &
                  7.21175825088309E+00,  4.31622272220567E+01, &
                  1.52989285046940E+02,  3.39320816734344E+02, &
                  4.51918953711873E+02,  3.00459261020162E+02 /),  &
        q(8) = (/ 1.00000000000000E+00,  1.27827273196294E+01, &
                  7.70001529352295E+01,  2.77585444743988E+02, &
                  6.38980264465631E+02,  9.31354094850610E+02, &
                  7.90950925327898E+02,  3.00459260956983E+02 /),  &
        r(5) = (/ 2.10144126479064E+00,  2.62370141675169E+01, &
                  2.13688200555087E+01,  4.65807828718470E+00, &
                  2.82094791773523E-01 /),   &
        s(4) = (/ 9.41537750555460E+01,  1.87114811799590E+02, &
                  9.90191814623914E+01,  1.80124575948747E+01 /)
!-------------------------
REAL :: ax, bot, c = .564189583547756, t, top, x2
!-------------------------
ax = ABS(x)
IF (ax < 0.5) THEN
  t = x*x
  top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0
  bot = ((b(1)*t + b(2))*t + b(3))*t + 1.0
  fn_val = x*(top/bot)
  RETURN

ELSE IF (ax < 4.0_dp) THEN
  top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax  &
        + p(6))*ax + p(7))*ax + p(8)
  bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax  &
        + q(6))*ax + q(7))*ax + q(8)
  fn_val = 0.5 + (0.5 - EXP(-x*x)*top/bot)
  IF (x < 0.0_dp) fn_val = -fn_val
  RETURN

ELSE IF (ax < 5.8) THEN
  x2 = x*x
  t = 1.0/x2
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4))*t + r(5)
  bot = (((s(1)*t + s(2))*t + s(3))*t + s(4))*t + 1.0
  fn_val = (c - top/(x2*bot)) / ax
  fn_val = 0.5 + (0.5 - EXP(-x2)*fn_val)
  IF (x < 0.0) fn_val = -fn_val
  RETURN

ELSE
  fn_val = SIGN(1.0,x)
END IF

RETURN
END FUNCTION erf



FUNCTION erfc(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL :: a(5) = (/ .771058495001320E-04, -.133733772997339E-02,  &
        .323076579225834E-01, .479137145607681E-01, .128379167095513E+00 /),  &
        b(3) = (/ .301048631703895E-02, .538971687740286E-01,  &
        .375795757275549E+00 /), p(8) = (/ -1.36864857382717E-07,  &
        5.64195517478974E-01, 7.21175825088309E+00, 4.31622272220567E+01,  &
        1.52989285046940E+02, 3.39320816734344E+02, 4.51918953711873E+02,  &
        3.00459261020162E+02 /), q(8) = (/ 1.00000000000000E+00,  &
        1.27827273196294E+01, 7.70001529352295E+01, 2.77585444743988E+02,  &
        6.38980264465631E+02, 9.31354094850610E+02, 7.90950925327898E+02,  &
        3.00459260956983E+02 /), r(5) = (/ 2.10144126479064E+00,  &
        2.62370141675169E+01, 2.13688200555087E+01, 4.65807828718470E+00,  &
        2.82094791773523E-01 /), s(4) = (/ 9.41537750555460E+01,  &
        1.87114811799590E+02, 9.90191814623914E+01, 1.80124575948747E+01 /),  &
        c = .564189583547756
REAL (dp) :: w
REAL      :: ax, t, top, bot, e
!-------------------------

!                     ABS(X) <= 0.5

ax = ABS(x)
IF (ax <= 0.5) THEN
  t = x * x
  top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0
  bot = ((b(1)*t + b(2))*t + b(3)) * t + 1.0
  fn_val = 0.5 + (0.5-x*(top/bot))
  RETURN
END IF

!                  0.5 < ABS(X) <= 4

IF (ax <= 4.0) THEN
  top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax + p(6))*ax  &
        + p(7)) * ax + p(8)
  bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax + q(6))*ax  &
        + q(7)) * ax + q(8)
  fn_val = top / bot
ELSE

!                       ABS(X) > 4

  IF (x <= -5.6) GO TO 10
  IF (x > 100.0) GO TO 20
  t = x * x
  IF (t > -exparg(1)) GO TO 20

  t = 1.0 / t
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4)) * t + r(5)
  bot = (((s(1)*t+s(2))*t+s(3))*t+s(4)) * t + 1.0
  fn_val = (c-t*top/bot) / ax
END IF

!                      FINAL ASSEMBLY

w = DBLE(x) * DBLE(x)
t = w
e = w - DBLE(t)
fn_val = ((0.5+(0.5-e))*EXP(-t)) * fn_val
IF (x < 0.0) fn_val = 2.0 - fn_val
RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

10 fn_val = 2.0
RETURN

!             LIMIT VALUE FOR LARGE POSITIVE X

20 fn_val = 0.0
RETURN
END FUNCTION erfc



FUNCTION erfc1(ind, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
!          ERFC1(IND,X) = ERFC(X)            IF IND = 0
!          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
!-----------------------------------------------------------------------
INTEGER, INTENT(IN) :: ind
REAL, INTENT(IN)    :: x
REAL                :: fn_val

! Local variables
REAL :: a(5) = (/ .771058495001320E-04, -.133733772997339E-02,  &
        .323076579225834E-01, .479137145607681E-01, .128379167095513E+00 /),  &
        b(3) = (/ .301048631703895E-02, .538971687740286E-01,  &
        .375795757275549E+00 /), p(8) = (/ -1.36864857382717E-07,  &
        5.64195517478974E-01, 7.21175825088309E+00, 4.31622272220567E+01,  &
        1.52989285046940E+02, 3.39320816734344E+02, 4.51918953711873E+02,  &
        3.00459261020162E+02 /), q(8) = (/ 1.00000000000000E+00,  &
        1.27827273196294E+01, 7.70001529352295E+01, 2.77585444743988E+02,  &
        6.38980264465631E+02, 9.31354094850610E+02, 7.90950925327898E+02,  &
        3.00459260956983E+02 /), r(5) = (/ 2.10144126479064E+00,  &
        2.62370141675169E+01, 2.13688200555087E+01, 4.65807828718470E+00,  &
        2.82094791773523E-01 /), s(4) = (/ 9.41537750555460E+01,  &
        1.87114811799590E+02, 9.90191814623914E+01, 1.80124575948747E+01 /),  &
        c = .564189583547756
REAL (dp) :: w
REAL      :: ax, t, top, bot, e
!-------------------------

!                     ABS(X) <= 0.5

ax = ABS(x)
IF (ax <= 0.5) THEN
  t = x * x
  top = ((((a(1)*t + a(2))*t + a(3))*t + a(4))*t + a(5)) + 1.0
  bot = ((b(1)*t + b(2))*t + b(3)) * t + 1.0
  fn_val = 0.5 + (0.5-x*(top/bot))
  IF (ind /= 0) fn_val = EXP(t) * fn_val
  RETURN
END IF

!                  0.5 < ABS(X) <= 4

IF (ax <= 4.0) THEN
  top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax + p(6))*ax  &
        + p(7)) * ax + p(8)
  bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax + q(6))*ax  &
        + q(7)) * ax + q(8)
  fn_val = top / bot
ELSE

!                      ABS(X) > 4

  IF (x <= -5.6) GO TO 10
  IF (ind == 0) THEN
    IF (x > 100.0) GO TO 20
    IF (x*x > -exparg(1)) GO TO 20
  END IF

  t = (1.0/x) ** 2
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4)) * t + r(5)
  bot = (((s(1)*t+s(2))*t+s(3))*t+s(4)) * t + 1.0
  fn_val = (c-t*top/bot) / ax
END IF

!                      FINAL ASSEMBLY

IF (ind /= 0) THEN
  IF (x < 0.0) fn_val = 2.0 * EXP(x*x) - fn_val
  RETURN
END IF
w = DBLE(x) * DBLE(x)
t = w
e = w - DBLE(t)
fn_val = ((0.5+(0.5-e))*EXP(-t)) * fn_val
IF (x < 0.0) fn_val = 2.0 - fn_val
RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

10 fn_val = 2.0
IF (ind /= 0) fn_val = 2.0 * EXP(x*x)
RETURN

!             LIMIT VALUE FOR LARGE POSITIVE X
!                       WHEN IND = 0

20 fn_val = 0.0
RETURN
END FUNCTION erfc1


SUBROUTINE dcerf(mo, z, w)
!-----------------------------------------------------------------------

!               COMPUTATION OF THE COMPLEX ERROR FUNCTION

!                          -----------------

!                       W = ERF(Z)    IF MO = 0
!                       W = ERFC(Z)   OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)    :: mo
REAL (dp), INTENT(IN)  :: z(2)
REAL (dp), INTENT(OUT) :: w(2)

! Local variables
REAL (dp) :: m, n, n2, n4, np1
REAL (dp) :: c = .56418958354775628694807945156077_dp, c2, d, d2, e,  &
             eps, r, sn, tol, x, y
REAL (dp) :: a0(2), an(2), b0(2), bn(2)
REAL (dp) :: g0(2), gn(2), h0(2), hn(2)
REAL (dp) :: qf(2), sm(2), sz(2), tm(2), ts(2), w0(2), wn(2)
!------------------------
!     C = 1/SQRT(PI)
!------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1._dp + EPS > 1._dp.

eps = EPSILON(1.0_dp)

!------------------------
x = z(1)
y = z(2)
sn = 1._dp
IF (x < 0._dp) THEN
  x = -x
  y = -y
  sn = -1._dp
END IF

r = x * x + y * y
sz(1) = x * x - y * y
sz(2) = 2._dp * x * y

IF (r > 1._dp) THEN
  IF (r >= 144._dp) GO TO 70
  IF (ABS(y) > 31.8_dp*x) GO TO 30
  IF (ABS(y) > 7.0_dp*x .AND. r < 64._dp) GO TO 30
  IF (ABS(y) > 3.2_dp*x .AND. r < 49._dp) GO TO 30
  IF (ABS(y) > 2.0_dp*x .AND. r < 36._dp) GO TO 30
  IF (ABS(y) > 1.2_dp*x .AND. r < 25._dp) GO TO 30
  IF (ABS(y) > 0.9_dp*x .AND. r < 16._dp) GO TO 30
  IF (r >= 6.25_dp) GO TO 50
  IF (ABS(y) > 0.6_dp*x) GO TO 30
  IF (r >= 4.0_dp) GO TO 20

  d = x - 2._dp
  IF (d*d + y*y < 1._dp) GO TO 20
  GO TO 30
END IF

!                          TAYLOR SERIES

c2 = c + c
tm(1) = c2 * x
tm(2) = c2 * y
sm(1) = tm(1)
sm(2) = tm(2)
tol = 2._dp * eps
m = 0._dp
10 m = m + 1._dp
d = m + m + 1._dp
ts(1) = tm(1) * sz(1) - tm(2) * sz(2)
ts(2) = tm(1) * sz(2) + tm(2) * sz(1)
tm(1) = -ts(1) / m
tm(2) = -ts(2) / m
ts(1) = tm(1) / d
ts(2) = tm(2) / d
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
IF (anorm(ts(1),ts(2)) > tol*anorm(sm(1),sm(2))) GO TO 10

IF (mo == 0) THEN
  w(1) = sn * sm(1)
  w(2) = sn * sm(2)
  RETURN
END IF
IF (sn /= 1._dp) THEN
  w(1) = 1._dp + sm(1)
  w(2) = sm(2)
  RETURN
END IF
w(1) = 0.5_dp + (0.5_dp-sm(1))
w(2) = -sm(2)
RETURN

!                  TAYLOR SERIES AROUND Z0 = 2

20 tm(1) = x
tm(2) = y
CALL erfcm2(0,tm,w)
IF (mo == 0) THEN
  w(1) = sn * (0.5_dp+(0.5_dp-w(1)))
  w(2) = -sn * w(2)
  RETURN
END IF
IF (sn > 0._dp) RETURN
w(1) = 2._dp - w(1)
w(2) = -w(2)
RETURN

!            PADE APPROXIMATION FOR THE TAYLOR SERIES
!                    FOR  (EXP(Z*Z)/Z)*ERF(Z)

30 d = 4._dp
IF (r > 16._dp) d = 16._dp
IF (r > 64._dp) d = 64._dp
d2 = d * d
CALL dcrec(sz(1),sz(2),w(1),w(2))
a0(1) = 1._dp
a0(2) = 0._dp
an(1) = (w(1)+4._dp/15._dp) * d
an(2) = w(2) * d
b0(1) = 1._dp
b0(2) = 0._dp
bn(1) = (w(1)-0.4_dp) * d
bn(2) = w(2) * d
CALL cdivid(an(1),an(2),bn(1),bn(2),wn(1),wn(2))
tol = 10._dp * eps
n4 = 0._dp

40 n4 = n4 + 4._dp
e = (n4+1._dp) * (n4+5._dp)
tm(1) = d * (w(1)-2._dp/e)
tm(2) = d * w(2)
e = d2 * (n4*(n4+2.0)) / ((n4-1.0)*(n4+3.0)*(n4+1.0)**2)

qf(1) = (tm(1)*an(1)-tm(2)*an(2)) + e * a0(1)
qf(2) = (tm(1)*an(2)+tm(2)*an(1)) + e * a0(2)
a0(1) = an(1)
a0(2) = an(2)
an(1) = qf(1)
an(2) = qf(2)
qf(1) = (tm(1)*bn(1)-tm(2)*bn(2)) + e * b0(1)
qf(2) = (tm(1)*bn(2)+tm(2)*bn(1)) + e * b0(2)
b0(1) = bn(1)
b0(2) = bn(2)
bn(1) = qf(1)
bn(2) = qf(2)

w0(1) = wn(1)
w0(2) = wn(2)
CALL cdivid(an(1),an(2),bn(1),bn(2),wn(1),wn(2))
IF (anorm(wn(1)-w0(1),wn(2)-w0(2)) > tol*anorm(wn(1),wn(2))) GO TO 40

c2 = c + c
sm(1) = c2 * (x*wn(1)-y*wn(2))
sm(2) = c2 * (x*wn(2)+y*wn(1))
e = EXP(-sz(1))
qf(1) = e * COS(-sz(2))
qf(2) = e * SIN(-sz(2))
tm(1) = qf(1) * sm(1) - qf(2) * sm(2)
tm(2) = qf(1) * sm(2) + qf(2) * sm(1)

w(1) = sn * tm(1)
w(2) = sn * tm(2)
IF (mo == 0) RETURN
w(1) = 1._dp - w(1)
w(2) = -w(2)
RETURN

!         PADE APPROXIMATION FOR THE ASYMPTOTIC EXPANSION
!                    FOR  Z*EXP(Z*Z)*ERFC(Z)

50 d = 4._dp * r
IF (r < 16._dp) d = 16._dp * r
d2 = d * d
tm(1) = sz(1) + sz(1)
tm(2) = sz(2) + sz(2)
g0(1) = 1._dp
g0(2) = 0._dp
gn(1) = (2._dp+tm(1)) / d
gn(2) = tm(2) / d
h0(1) = 1._dp
h0(2) = 0._dp
tm(1) = 3._dp + tm(1)
hn(1) = tm(1) / d
hn(2) = tm(2) / d
CALL cdivid(gn(1),gn(2),hn(1),hn(2),wn(1),wn(2))
np1 = 1._dp
tol = 10._dp * eps

60 n = np1
np1 = n + 1._dp
n2 = n + n
e = (n2*(n2+1._dp)) / d2
tm(1) = tm(1) + 4._dp
qf(1) = (tm(1)*gn(1)-tm(2)*gn(2)) / d - e * g0(1)
qf(2) = (tm(1)*gn(2)+tm(2)*gn(1)) / d - e * g0(2)
g0(1) = gn(1)
g0(2) = gn(2)
gn(1) = qf(1)
gn(2) = qf(2)
qf(1) = (tm(1)*hn(1)-tm(2)*hn(2)) / d - e * h0(1)
qf(2) = (tm(1)*hn(2)+tm(2)*hn(1)) / d - e * h0(2)
h0(1) = hn(1)
h0(2) = hn(2)
hn(1) = qf(1)
hn(2) = qf(2)

w0(1) = wn(1)
w0(2) = wn(2)
CALL cdivid(gn(1),gn(2),hn(1),hn(2),wn(1),wn(2))
IF (anorm(wn(1)-w0(1),wn(2)-w0(2)) > tol*anorm(wn(1),wn(2))) GO TO 60

tm(1) = x * hn(1) - y * hn(2)
tm(2) = x * hn(2) + y * hn(1)
CALL cdivid(c*gn(1),c*gn(2),tm(1),tm(2),sm(1),sm(2))
GO TO 90

!                      ASYMPTOTIC EXPANSION

70 CALL dcrec(x,y,tm(1),tm(2))
sm(1) = tm(1)
sm(2) = tm(2)
qf(1) = tm(1) * tm(1) - tm(2) * tm(2)
qf(2) = 2._dp * tm(1) * tm(2)
tol = 2._dp * eps
d = -0.5_dp
80 d = d + 1._dp
ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
tm(1) = -d * ts(1)
tm(2) = -d * ts(2)
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (anorm(tm(1),tm(2)) > tol*anorm(sm(1),sm(2))) GO TO 80
sm(1) = c * sm(1)
sm(2) = c * sm(2)
IF (x < 1.d-2) GO TO 100

!                       TERMINATION

90 e = EXP(-sz(1))
qf(1) = e * COS(-sz(2))
qf(2) = e * SIN(-sz(2))
ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)

IF (mo == 0) THEN
  w(1) = sn * (0.5_dp+(0.5_dp-sm(1)))
  w(2) = -sn * sm(2)
  RETURN
END IF
IF (sn /= 1._dp) THEN
  w(1) = 2._dp - sm(1)
  w(2) = -sm(2)
  RETURN
END IF
w(1) = sm(1)
w(2) = sm(2)
RETURN

!               MODIFIED ASYMPTOTIC EXPANSION

100 e = EXP(-sz(1))
qf(1) = e * COS(-sz(2))
qf(2) = e * SIN(-sz(2))
w(1) = qf(1) * sm(1) - qf(2) * sm(2)
w(2) = qf(1) * sm(2) + qf(2) * sm(1)
IF (mo /= 0) THEN
  w(1) = 1._dp + sn * w(1)
  w(2) = sn * w(2)
  RETURN
END IF
IF (sn < 0.0) RETURN
w(1) = -w(1)
w(2) = -w(2)
RETURN
END SUBROUTINE dcerf


SUBROUTINE dcerfc(mo, z, w)
!-----------------------------------------------------------------------

!              COMPUTATION OF THE COMPLEX COERROR FUNCTION

!                           ----------------

!           W = ERFC(Z)           IF MO = 0 OR REAL(Z) < 0
!           W = EXP(X*X)*ERFC(Z)  OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)    :: mo
REAL (dp), INTENT(IN)  :: z(2)
REAL (dp), INTENT(OUT) :: w(2)

! Local variables
REAL (dp) :: m, n, n2, n4, np1
REAL (dp) :: c = .56418958354775628694807945156077_dp, c2, d, d2, e,  &
             eps, r, sn, tol, x, y
REAL (dp) :: a0(2), an(2), b0(2), bn(2)
REAL (dp) :: g0(2), gn(2), h0(2), hn(2)
REAL (dp) :: qf(2), sm(2), sz(2), tm(2), ts(2), w0(2), wn(2)
!------------------------
!     C = 1/SQRT(PI)
!------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1._dp + EPS > 1._dp.

eps = EPSILON(1.0_dp)

!------------------------
x = z(1)
y = z(2)
sn = 1._dp
IF (x < 0._dp) THEN
  x = -x
  y = -y
  sn = -1._dp
END IF

IF (mo == 0.OR.sn /= 1._dp.OR.MAX(x,ABS(y)) < 144._dp) THEN
  r = x * x + y * y
  sz(1) = x * x - y * y
  sz(2) = 2._dp * x * y

  IF (r > 1._dp) THEN
    IF (r >= 144._dp) GO TO 70
    IF (ABS(y) > 31.8_dp*x) GO TO 30
    IF (ABS(y) > 7.0_dp*x .AND. r < 64._dp) GO TO 30
    IF (ABS(y) > 3.2_dp*x .AND. r < 49._dp) GO TO 30
    IF (ABS(y) > 2.0_dp*x .AND. r < 36._dp) GO TO 30
    IF (ABS(y) > 1.2_dp*x .AND. r < 25._dp) GO TO 30
    IF (ABS(y) > 0.9_dp*x .AND. r < 16._dp) GO TO 30
    IF (r >= 6.25_dp) GO TO 50
    IF (ABS(y) > 0.6_dp*x) GO TO 30
    IF (r >= 4.0_dp) GO TO 20

    d = x - 2._dp
    IF (d*d+y*y < 1._dp) GO TO 20
    GO TO 30
  END IF

!                          TAYLOR SERIES

  c2 = c + c
  tm(1) = c2 * x
  tm(2) = c2 * y
  sm(1) = tm(1)
  sm(2) = tm(2)
  tol = 2._dp * eps
  m = 0._dp
  10 m = m + 1._dp
  d = m + m + 1._dp
  ts(1) = tm(1) * sz(1) - tm(2) * sz(2)
  ts(2) = tm(1) * sz(2) + tm(2) * sz(1)
  tm(1) = -ts(1) / m
  tm(2) = -ts(2) / m
  ts(1) = tm(1) / d
  ts(2) = tm(2) / d
  sm(1) = sm(1) + ts(1)
  sm(2) = sm(2) + ts(2)
  IF (anorm(ts(1),ts(2)) > tol*anorm(sm(1),sm(2))) GO TO 10

  IF (sn /= 1._dp) THEN
    w(1) = 1._dp + sm(1)
    w(2) = sm(2)
    RETURN
  END IF
  sm(1) = 0.5_dp + (0.5_dp-sm(1))
  sm(2) = -sm(2)
  IF (mo == 0) GO TO 100

  e = EXP(sz(1))
  qf(1) = e * COS(sz(2))
  qf(2) = e * SIN(sz(2))
  w(1) = qf(1) * sm(1) - qf(2) * sm(2)
  w(2) = qf(1) * sm(2) + qf(2) * sm(1)
  RETURN

!                  TAYLOR SERIES AROUND Z0 = 2

  20 IF (sn >= 0._dp) THEN
    CALL erfcm2(mo,z,w)
    RETURN
  END IF
  tm(1) = x
  tm(2) = y
  CALL erfcm2(0,tm,w)
  w(1) = 2._dp - w(1)
  w(2) = -w(2)
  RETURN

!            PADE APPROXIMATION FOR THE TAYLOR SERIES
!                    FOR  (EXP(Z*Z)/Z)*ERF(Z)

  30 d = 4._dp
  IF (r > 16._dp) d = 16._dp
  IF (r > 64._dp) d = 64._dp
  d2 = d * d
  CALL dcrec(sz(1),sz(2),w(1),w(2))
  a0(1) = 1._dp
  a0(2) = 0._dp
  an(1) = (w(1)+4._dp/15._dp) * d
  an(2) = w(2) * d
  b0(1) = 1._dp
  b0(2) = 0._dp
  bn(1) = (w(1)-0.4_dp) * d
  bn(2) = w(2) * d
  CALL cdivid(an(1),an(2),bn(1),bn(2),wn(1),wn(2))
  tol = 10._dp * eps
  n4 = 0._dp

  40 n4 = n4 + 4._dp
  e = (n4+1._dp) * (n4+5._dp)
  tm(1) = d * (w(1)-2._dp/e)
  tm(2) = d * w(2)
  e = d2 * (n4*(n4+2.0)) / ((n4-1.0)*(n4+3.0)*(n4+1.0)**2)

  qf(1) = (tm(1)*an(1)-tm(2)*an(2)) + e * a0(1)
  qf(2) = (tm(1)*an(2)+tm(2)*an(1)) + e * a0(2)
  a0(1) = an(1)
  a0(2) = an(2)
  an(1) = qf(1)
  an(2) = qf(2)
  qf(1) = (tm(1)*bn(1)-tm(2)*bn(2)) + e * b0(1)
  qf(2) = (tm(1)*bn(2)+tm(2)*bn(1)) + e * b0(2)
  b0(1) = bn(1)
  b0(2) = bn(2)
  bn(1) = qf(1)
  bn(2) = qf(2)

  w0(1) = wn(1)
  w0(2) = wn(2)
  CALL cdivid(an(1),an(2),bn(1),bn(2),wn(1),wn(2))
  IF (anorm(wn(1)-w0(1),wn(2)-w0(2)) > tol*anorm(wn(1),wn(2)))  &
  GO TO 40

  c2 = c + c
  sm(1) = c2 * (x*wn(1)-y*wn(2))
  sm(2) = c2 * (x*wn(2)+y*wn(1))

  IF (mo /= 0 .AND. sn == 1._dp) THEN
    e = EXP(sz(1))
    w(1) = e * COS(sz(2)) - sm(1)
    w(2) = e * SIN(sz(2)) - sm(2)
    RETURN
  END IF
  e = EXP(-sz(1))
  qf(1) = e * COS(-sz(2))
  qf(2) = e * SIN(-sz(2))
  tm(1) = qf(1) * sm(1) - qf(2) * sm(2)
  tm(2) = qf(1) * sm(2) + qf(2) * sm(1)
  w(1) = 1._dp - sn * tm(1)
  w(2) = -sn * tm(2)
  RETURN

!         PADE APPROXIMATION FOR THE ASYMPTOTIC EXPANSION
!                    FOR  Z*EXP(Z*Z)*ERFC(Z)

  50 d = 4._dp * r
  IF (r < 16._dp) d = 16._dp * r
  d2 = d * d
  tm(1) = sz(1) + sz(1)
  tm(2) = sz(2) + sz(2)
  g0(1) = 1._dp
  g0(2) = 0._dp
  gn(1) = (2._dp+tm(1)) / d
  gn(2) = tm(2) / d
  h0(1) = 1._dp
  h0(2) = 0._dp
  tm(1) = 3._dp + tm(1)
  hn(1) = tm(1) / d
  hn(2) = tm(2) / d
  CALL cdivid(gn(1),gn(2),hn(1),hn(2),wn(1),wn(2))
  np1 = 1._dp
  tol = 10._dp * eps

  60 n = np1
  np1 = n + 1._dp
  n2 = n + n
  e = (n2*(n2+1._dp)) / d2
  tm(1) = tm(1) + 4._dp
  qf(1) = (tm(1)*gn(1)-tm(2)*gn(2)) / d - e * g0(1)
  qf(2) = (tm(1)*gn(2)+tm(2)*gn(1)) / d - e * g0(2)
  g0(1) = gn(1)
  g0(2) = gn(2)
  gn(1) = qf(1)
  gn(2) = qf(2)
  qf(1) = (tm(1)*hn(1)-tm(2)*hn(2)) / d - e * h0(1)
  qf(2) = (tm(1)*hn(2)+tm(2)*hn(1)) / d - e * h0(2)
  h0(1) = hn(1)
  h0(2) = hn(2)
  hn(1) = qf(1)
  hn(2) = qf(2)

  w0(1) = wn(1)
  w0(2) = wn(2)
  CALL cdivid(gn(1),gn(2),hn(1),hn(2),wn(1),wn(2))
  IF (anorm(wn(1)-w0(1),wn(2)-w0(2)) > tol*anorm(wn(1),wn(2)))  &
  GO TO 60

  tm(1) = x * hn(1) - y * hn(2)
  tm(2) = x * hn(2) + y * hn(1)
  CALL cdivid(c*gn(1),c*gn(2),tm(1),tm(2),sm(1),sm(2))
  GO TO 90
END IF

!                      ASYMPTOTIC EXPANSION

70 CALL dcrec(x,y,tm(1),tm(2))
sm(1) = tm(1)
sm(2) = tm(2)
qf(1) = tm(1) * tm(1) - tm(2) * tm(2)
qf(2) = 2._dp * tm(1) * tm(2)
tol = 2._dp * eps
d = -0.5_dp
80 d = d + 1._dp
ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
tm(1) = -d * ts(1)
tm(2) = -d * ts(2)
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (anorm(tm(1),tm(2)) > tol*anorm(sm(1),sm(2))) GO TO 80
sm(1) = c * sm(1)
sm(2) = c * sm(2)
IF (x < 1.d-2) GO TO 110

!                       TERMINATION

90 IF (mo == 0.OR.sn /= 1._dp) THEN
  e = EXP(-sz(1))
  qf(1) = e * COS(-sz(2))
  qf(2) = e * SIN(-sz(2))
  ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
  ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
  sm(1) = ts(1)
  sm(2) = ts(2)

  IF (sn /= 1._dp) THEN
    w(1) = 2._dp - sm(1)
    w(2) = -sm(2)
    RETURN
  END IF
END IF
100 w(1) = sm(1)
w(2) = sm(2)
RETURN

!               MODIFIED ASYMPTOTIC EXPANSION

110 IF (mo == 0.OR.sn /= 1._dp) THEN
  e = EXP(-sz(1))
  qf(1) = e * COS(-sz(2))
  qf(2) = e * SIN(-sz(2))
  ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
  ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
  w(1) = 1._dp + sn * ts(1)
  w(2) = sn * ts(2)
  RETURN
END IF

IF (ABS(y) >= 100._dp) GO TO 100
IF (sz(1) <= dxparg(1)) GO TO 100
e = EXP(sz(1))
w(1) = e * COS(sz(2)) + sm(1)
w(2) = e * SIN(sz(2)) + sm(2)
RETURN
END SUBROUTINE dcerfc


SUBROUTINE erfcm2(mo, z, w)
!-----------------------------------------------------------------------
!           CALCULATION OF ERFC(Z) USING THE TAYLOR SERIES
!                          AROUND Z0 = 2
!-----------------------------------------------------------------------
INTEGER, INTENT(IN)    :: mo
REAL (dp), INTENT(IN)  :: z(2)
REAL (dp), INTENT(OUT) :: w(2)

REAL (dp) :: a(63) = (/ .20000000000000000000000000000000000D+01,  &
    .23333333333333333333333333333333333D+01, .16666666666666666666666666666666667D+01,  &
    .63333333333333333333333333333333333D+00, -.22222222222222222222222222222222222D-01,  &
    -.16349206349206349206349206349206349D+00, -.76984126984126984126984126984126984D-01,  &
    -.24250440917107583774250440917107584D-02, .12716049382716049382716049382716049D-01,  &
    .50208433541766875100208433541766875D-02, -.25305969750414194858639303083747528D-03,  &
    -.78593217482106370995259884148773038D-03, -.19118154038788959423880058800693721D-03,  &
    .46324144207742091339974937858535742D-04, .33885549097189308829520469732109944D-04,  &
    .28637897646612243562134629672756034D-05, -.29071891082127275370004560446169188D-05,  &
    -.89674405786490646425523560263096103D-06, .96069103941908684338469767911200105D-07,  &
    .99432863129093191401848891268744113D-07, .97610310501460621303387795457283579D-08,  &
    -.65557500375673133822289344530892436D-08, -.18706782059105426900361744016236561D-08,  &
    .20329898993447386223176373714372370D-09, .16941915827254374668448114614201210D-09,  &
    .10619149520827430973786114446699534D-10, -.10136148256511788733365237088810952D-10,  &
    -.21042890133669970575386166675721692D-11, .37186985840699828780916522245407087D-12,  &
    .17921843632701679986488128324051002D-12, -.89823991804248069863542565948598397D-16,  &
    -.10533182313660970970232171410372199D-13, -.12340742690978398320850088252659714D-14,  &
    .44315624546581333350568244777175883D-15, .11584041639989442481950487524296214D-15,  &
    -.10765703619385988116658460442219647D-16, -.70653158723054941879586082239984222D-17,  &
    -.18708903154917138727191793341667090D-18, .32549879318817103966053527398133297D-18,  &
    .40654116689599228385911733319215613D-19, -.11250074516817311101947327325293424D-19,  &
    -.28923865378584966737386008432031980D-20, .23653053641701517160704870522922706D-21,  &
    .14665384680061888088099002254334292D-21, .26971039707314316218154193225264469D-23,  &
    -.58753834789274356433279671015522650D-23, -.59960357240498652932299485494869633D-24,  &
    .18586826578121663981412155416486531D-24, .38364131854692721721867481914852428D-25,  &
    -.41342210492630142578080062451711039D-26, -.17646283105274988992381528904600860D-26,  &
    .19828685934364181151988692232131608D-28, .65592252170840353572672782446212733D-28,  &
    .40626551379996340638338449938639730D-29, -.20097984104191034713653294173834095D-29,  &
    -.28104226475997460044096389060743131D-30, .48705319298749358709127987806547949D-31,  &
    .12664655832830787747161769929972617D-31, -.75168312488894341862391776330113688D-33,  &
    -.45760473722605993842481669806804415D-33, -.56725491529575395930156379514718000D-35,  &
    .13932664042920082608489441616061541D-34, .10452448992516358449586503951463322D-35 /)
REAL (dp), PARAMETER :: c = .20666985354092053857068941306585476D-01,  &
                        e = .46777349810472658379307436327470714D-02
REAL (dp)            :: eps, h(2), t(2), tol, x, y
INTEGER              :: n, j
!------------------------------
!     C = (2/SQRT(PI))*EXP(-4)
!     E = ERFC(2)
!------------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1._dp + EPS > 1._dp

eps = EPSILON(1.0_dp)

!------------------------------
tol = eps * 1.d+12
h(1) = 1._dp + (1._dp-z(1))
h(2) = -z(2)

x = 1._dp
y = 0._dp
w(1) = a(30)
w(2) = 0._dp
DO n = 31, 63
  t(1) = x * h(1) - y * h(2)
  t(2) = x * h(2) + y * h(1)
  x = t(1)
  y = t(2)
  t(1) = a(n) * x
  t(2) = a(n) * y
  w(1) = w(1) + t(1)
  w(2) = w(2) + t(2)
  IF (anorm(t(1),t(2)) <= tol*anorm(w(1),w(2))) GO TO 20
END DO

20 DO j = 1, 29
  n = 30 - j
  x = h(1) * w(1) - h(2) * w(2)
  w(2) = h(1) * w(2) + h(2) * w(1)
  w(1) = a(n) + x
END DO
x = h(1) * w(1) - h(2) * w(2)
w(2) = h(1) * w(2) + h(2) * w(1)
w(1) = 1._dp + x

x = c * (h(1)*w(1)-h(2)*w(2))
w(2) = c * (h(1)*w(2)+h(2)*w(1))
w(1) = e + x
IF (mo == 0) RETURN

!                     COMPUTE EXP(Z*Z)*ERFC(Z)

x = z(1) * z(1) - z(2) * z(2)
y = 2._dp * z(1) * z(2)
x = EXP(x)
t(1) = x * COS(y)
t(2) = x * SIN(y)
x = t(1) * w(1) - t(2) * w(2)
y = t(1) * w(2) + t(2) * w(1)
w(1) = x
w(2) = y
RETURN
END SUBROUTINE erfcm2


FUNCTION derf(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!        REAL (dp) EVALUATION OF THE ERROR FUNCTION
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: ax, t, w
INTEGER   :: i, k
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
IF (ax <= 1._dp) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = x * (1._dp+w)
  RETURN
END IF

!                     ABS(X) > 1

IF (ax < 8.5_dp) THEN
  fn_val = 0.5_dp + (0.5_dp-EXP(-x*x)*derfc0(ax))
  IF (x < 0._dp) fn_val = -fn_val
  RETURN
END IF

!                 LIMIT VALUE FOR LARGE X

fn_val = SIGN(1._dp,x)
RETURN
END FUNCTION derf


FUNCTION derfc(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: ax, t, w
INTEGER   :: i, k
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
IF (ax <= 1._dp) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = 0.5_dp + (0.5_dp-x*(1._dp+w))
  RETURN
END IF

!                       X < -1

IF (x <= 0._dp) THEN
  fn_val = 2._dp
  IF (x < -8.3_dp) RETURN
  t = x * x
  fn_val = 2._dp - EXP(-t) * derfc0(ax)
  RETURN
END IF

!                       X > 1

fn_val = 0._dp
IF (x > 100._dp) RETURN
t = x * x
IF (t > -dxparg(1)) RETURN
fn_val = EXP(-t) * derfc0(x)
RETURN
END FUNCTION derfc


FUNCTION derfc1(ind, x) RESULT(fn_val)
!-----------------------------------------------------------------------

!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

!          DERFC1(IND,X) = ERFC(X)           IF IND = 0
!          DERFC1(IND,X) = EXP(X*X)*ERFC(X)  OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)   :: ind
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: ax, t, w
INTEGER   :: i, k
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
IF (ax <= 1._dp) THEN
  t = x * x
  w = a(21)
  DO i = 1, 20
    k = 21 - i
    w = t * w + a(k)
  END DO
  fn_val = 0.5_dp + (0.5_dp-x*(1._dp+w))
  IF (ind /= 0) fn_val = EXP(t) * fn_val
  RETURN
END IF

!                       X < -1

IF (x <= 0._dp) THEN
  IF (x < -8.3_dp) GO TO 20
  IF (ind /= 0) THEN
    fn_val = 2._dp * EXP(x*x) - derfc0(ax)
    RETURN
  END IF
  fn_val = 2._dp - EXP(-x*x) * derfc0(ax)
  RETURN
END IF

!                       X > 1

IF (ind /= 0) THEN
  fn_val = derfc0(x)
  RETURN
END IF
fn_val = 0._dp
IF (x > 100._dp) RETURN
t = x * x
IF (t > -dxparg(1)) RETURN
fn_val = EXP(-t) * derfc0(x)
RETURN

!             LIMIT VALUE FOR LARGE NEGATIVE X

20 fn_val = 2._dp
IF (ind /= 0) fn_val = 2._dp * EXP(x*x)
RETURN
END FUNCTION derfc1


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
REAL (dp), PARAMETER :: rpinv = .56418958354775628694807945156077259_dp
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

  t = (x-3.75_dp) / (x+3.75_dp)
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


FUNCTION erfi(p, q) RESULT(fn_val)
!-----------------------------------------------------------------------

!              EVALUATION OF THE INVERSE ERROR FUNCTION

!                        ---------------

!     FOR 0 <= P < 1,  W = ERFI(P,Q) WHERE ERF(W) = P. IT IS
!     ASSUMED THAT Q = 1 - P.  IF P < 0, Q <= 0, OR P + Q IS
!     NOT 1, THEN ERFI(P,Q) IS SET TO A NEGATIVE VALUE.

!-----------------------------------------------------------------------
!     REFERENCE. MATHEMATICS OF COMPUTATION,OCT.1976,PP.827-830.
!                  J.M.BLAIR,C.A.EDWARDS,J.H.JOHNSON
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: p, q
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: a(6) = (/ .1400216916161353E+03, -.7204275515686407E+03,  &
                             .1296708621660511E+04, -.9697932901514031E+03,  &
                             .2762427049269425E+03, -.2012940180552054E+02 /),  &
                   b(6) = (/ .1291046303114685E+03, -.7312308064260973E+03,  &
                             .1494970492915789E+04, -.1337793793683419E+04,  &
                             .5033747142783567E+03, -.6220205554529216E+02 /),  &
                 a1(7)  = (/ -.1690478046781745E+00, .3524374318100228E+0,  &
                             -.2698143370550352E+02, .9340783041018743E+02,  &
                             -.1455364428646732E+03, .8805852004723659E+02,  &
                             -.1349018591231947E+02 /),  &
                  b1(7) = (/ -.1203221171313429E+00, .2684812231556632E+00,  &
                             -.2242485268704865E+02, .8723495028643494E+02,  &
                             -.1604352408444319E+03, .1259117982101525E+03,  &
                             -.3184861786248824E+02 /),   &
                  a2(9) = (/ .3100808562552958E-04, .4097487603011940E-02,  &
                             .1214902662897276E+00, .1109167694639028E+01,  &
                             .3228379855663924E+01, .2881691815651599E+01,  &
                             .2047972087262996E+01, .8545922081972148E+00,  &
                             .3551095884622383E-02 /),   &
                  b2(8) = (/ .3100809298564522E-04, .4097528678663915E-02,  &
                             .1215907800748757E+00, .1118627167631696E+01,  &
                             .3432363984305290E+01, .4140284677116202E+01,  &
                             .4119797271272204E+01, .2162961962641435E+01 /), &
                  a3(9) = (/ .3205405422062050E-08, .1899479322632128E-05,  &
                             .2814223189858532E-03, .1370504879067817E-01,  &
                             .2268143542005976E+00, .1098421959892340E+01,  &
                             .6791143397056208E+00, -.834334189167721E+00,  &
                             .3421951267240343E+00 /),   &
                  b3(6) = (/ .3205405053282398E-08, .1899480592260143E-05,  &
                             .2814349691098940E-03, .1371092249602266E-01,  &
                             .2275172815174473E+00, .1125348514036959E+01 /), &
                      c = .5625, c1 = .87890625, c2 = -.2302585092994046E+03
REAL :: eps, s, t, v, v1
!-----------------------------------------------------------------------
!     C2 = LN(1.E-100)
!-----------------------------------------------------------------------
IF (p >= 0.0 .AND. q > 0.0) THEN
  eps = MAX(EPSILON(1.0), 1.e-15)
  t = 0.5 + (0.5-(p+q))
  IF (ABS(t) > 3.0*eps) GO TO 10

!                      0 <= P <= 0.75

  IF (p <= 0.75) THEN
    v = p * p - c
    t = p * (((((a(6)*v + a(5))*v + a(4))*v + a(3))*v + a(2))*v + a(1))
    s = (((((v + b(6))*v + b(5))*v + b(4))*v + b(3))*v + b(2)) * v + b(1)
  ELSE

!                    0.75 < P <= 0.9375

    IF (p <= 0.9375) THEN
      v = p * p - c1
      t = p * ((((((a1(7)*v + a1(6))*v + a1(5))*v + a1(4))*v + a1(3))*v + a1(2))*v + a1(1))
      s = ((((((v + b1(7))*v + b1(6))*v + b1(5))*v + b1(4))*v + b1(3))*v +  &
          b1(2)) * v + b1(1)
    ELSE

!                  1.E-100 <= Q < 0.0625

      v1 = LOG(q)
      v = 1.0 / SQRT(-v1)
      IF (v1 >= c2) THEN
        t = (((((((a2(9)*v + a2(8))*v + a2(7))*v + a2(6))*v + a2(5))*v +  &
            a2(4))*v + a2(3))*v + a2(2)) * v + a2(1)
        s = v * ((((((((v + b2(8))*v + b2(7))*v + b2(6))*v + b2(5))*v +  &
            b2(4))*v + b2(3))*v + b2(2))*v + b2(1))
      ELSE

!                 1.E-10000 <= Q < 1.E-100

        t = (((((((a3(9)*v + a3(8))*v + a3(7))*v + a3(6))*v + a3(5))*v +  &
            a3(4))*v + a3(3))*v + a3(2)) * v + a3(1)
        s = v * ((((((v + b3(6))*v + b3(5))*v + b3(4))*v + b3(3))*v +  &
                b3(2))*v + b3(1))
      END IF
    END IF
  END IF
  fn_val = t / s
  RETURN
END IF

!                         ERROR RETURN

fn_val = -1.0
RETURN
10 fn_val = -2.0
RETURN
END FUNCTION erfi



FUNCTION derfi(p, q) RESULT(fn_val)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: p, q
REAL (dp)             :: fn_val

!                  REAL (dp) COMPUTATION OF
!                    THE INVERSE ERROR FUNCTION

!                         ----------------

!     FOR 0 <= P <= 1,  W = DERFI(P,Q) WHERE ERF(W) = P. IT
!     IS ASSUMED THAT Q = 1 - P. IF P < 0, Q <= 0, OR P + Q
!     IS NOT 1, THEN DERFI(P,Q) IS SET TO A NEGATIVE VALUE.

!-----------------------------------------------------------------------
!     REFERENCE. MATHEMATICS OF COMPUTATION,OCT.1976,PP.827-830.
!                  J.M.BLAIR,C.A.EDWARDS,J.H.JOHNSON
!-----------------------------------------------------------------------
REAL (dp) :: c = .5625_dp, c1 = .87890625_dp, c2  &
        = -.2302585092994045684017991454684364D+03, r  &
        = .8862269254527580136490837416705726D+00, eps, f, lnq, s, t, x
REAL (dp) :: a(7) = (/ .841467547194693616D-01,  &
        .160499904248262200D+01, .809451641478547505D+01,  &
        .164273396973002581D+02, .154297507839223692D+02,  &
        .669584134660994039D+01, .108455979679682472D+01 /), a1(7)  &
        = (/ .552755110179178015D+2, .657347545992519152D+3,  &
        .124276851197202733D+4, .818859792456464820D+3,  &
        .234425632359410093D+3, .299942187305427917D+2,  &
        .140496035731853946D+1 /), a2(7) = (/ .500926197430588206D+1,  &
        .111349802614499199D+3, .353872732756132161D+3,  &
        .356000407341490731D+3, .143264457509959760D+3,  &
        .240823237485307567D+2, .140496035273226366D+1 /), a3(11)  &
        = (/ .237121026548776092D4, .732899958728969905D6,  &
        .182063754893444775D7, .269191299062422172D7, .304817224671614253D7  &
       , .130643103351072345D7, .296799076241952125D6,  &
        .457006532030955554D5, .373449801680687213D4, .118062255483596543D3  &
       , .100000329157954960D1 /), a4(9) = (/ .154269429680540807D12,  &
        .430207405012067454D12, .182623446525965017D12,  &
        .248740194409838713D11, .133506080294978121D10,  &
        .302446226073105850D08, .285909602878724425D06,  &
        .101789226017835707D04, .100000004821118676D01 /),  &
        b(7) = (/ .352281538790042405D-02, .293409069065309557D+00,  &
        .326709873508963100D+01, .123611641257633210D+02,  &
        .207984023857547070D+02, .170791197367677668D+02,  &
        .669253523595376683D+01 /), b1(6) = (/ .179209835890172156D+3,  &
        .991315839349539886D+3, .138271033653003487D+4, .764020340925985926D+3,  &
        .194354053300991923D+3, .228139510050586581D+2 /),  &
        b2(6) = (/ .209004294324106981D+2, .198607335199741185D+3,  &
        .439311287748524270D+3, .355415991280861051D+3, .123303672628828521D+3,  &
        .186060775181898848D+2 /), b3(10) = (/ .851911109952055378D6,  &
        .194746720192729966D7, .373640079258593694D7, .397271370110424145D7,  &
        .339457682064283712D7 , .136888294898155938D7, .303357770911491406D6,  &
        .459721480357533823D5, .373762573565814355D4, .118064334590001264D3 /),  &
        b4(9) = (/ .220533001293836387D12, .347822938010402687D12,  &
        .468373326975152250D12, .185251723580351631D12,  &
        .249464490520921771D11, .133587491840784926D10,  &
        .302480682561295591D08, .285913799407861384D06,  &
        .101789250893050230D04 /)
!-----------------------------------------------------------------------
!     C2 = LN(1.E-100)
!     R  = SQRT(PI)/2
!-----------------------------------------------------------------------
IF (p >= 0._dp .AND. q > 0._dp) THEN
  eps = EPSILON(1.0_dp)
  t = 0.5_dp + (0.5_dp-(p+q))
  IF (ABS(t) > 3._dp*eps) GO TO 10

!                      0 <= P <= 0.75

  IF (p <= 0.75_dp) THEN
    x = c - p * p
    s = (((((a(1)*x+a(2))*x+a(3))*x+a(4))*x+a(5))*x+a(6)) * x +a(7)
    t = ((((((b(1)*x+b(2))*x+b(3))*x+b(4))*x+b(5))*x+b(6))*x+b(7)) * x + 1._dp
    fn_val = p * (s/t)
    IF (eps > 1.d-19) RETURN

    x = fn_val
    f = derf(x) - p
    fn_val = x - r * EXP(x*x) * f
    RETURN
  END IF

!                    0.75 < P <= 0.9375

  IF (p <= 0.9375_dp) THEN
    x = c1 - p * p
    IF (x <= 0.1_dp) THEN
      s = ((((((a1(1)*x+a1(2))*x+a1(3))*x+a1(4))*x+a1(5))*x+a1(6))*x+a1(7))
      t = ((((((b1(1)*x+b1(2))*x+b1(3))*x+b1(4))*x+b1(5))*x+b1(6))*x+1._dp)
    ELSE

      s = ((((((a2(1)*x+a2(2))*x+a2(3))*x+a2(4))*x+a2(5))*x+a2(6))*x+a2(7))
      t = ((((((b2(1)*x+b2(2))*x+b2(3))*x+b2(4))*x+b2(5))*x+b2(6))*x+1._dp)
    END IF

    fn_val = p * (s/t)
    IF (eps > 1.d-19) RETURN

    x = fn_val
    t = derfc1(1,x) - EXP(x*x) * q
    fn_val = x + r * t
    RETURN
  END IF

!                  1.E-100 <= Q < 0.0625

  lnq = LOG(q)
  x = 1._dp / SQRT(-lnq)
  IF (lnq >= c2) THEN
    s = (((((((((a3(1)*x+a3(2))*x+a3(3))*x+a3(4))*x+a3(5))*x+  &
    a3(6))*x+a3(7))*x+a3(8))*x+a3(9))*x+a3(10)) * x + a3(11)
    t = (((((((((b3(1)*x+b3(2))*x+b3(3))*x+b3(4))*x+b3(5))*x+  &
    b3(6))*x+b3(7))*x+b3(8))*x+b3(9))*x+b3(10)) * x + 1._dp
  ELSE

!                 1.E-10000 <= Q < 1.E-100

    s = (((((((a4(1)*x+a4(2))*x+a4(3))*x+a4(4))*x+a4(5))*x+a4(6))*  &
    x+a4(7))*x+a4(8)) * x + a4(9)
    t = ((((((((b4(1)*x+b4(2))*x+b4(3))*x+b4(4))*x+b4(5))*x+  &
    b4(6))*x+b4(7))*x+b4(8))*x+b4(9)) * x + 1._dp
  END IF

  fn_val = s / (x*t)
  IF (eps > 5.d-20) RETURN

  x = fn_val
  t = derfc1(1,x)
  f = (LOG(t)-lnq) - x * x
  fn_val = x + r * t * f
  RETURN
END IF

!                         ERROR RETURN

fn_val = -1._dp
RETURN
10 fn_val = -2._dp
RETURN
END FUNCTION derfi



FUNCTION aerf(x, h) RESULT(fn_val)
!-----------------------------------------------------------------------
!             COMPUTATION OF ERF(X + H) - ERF(X - H)
!-----------------------------------------------------------------------
!     C = 2/SQRT(PI)
!     P = LN(9*SQRT(PI))
!---------------------
REAL, INTENT(IN) :: x, h
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: c = 1.12837916709551257, p = 2.76959
REAL             :: ah, ax, a2, dn2, e, eps, hf, hg, h2, h3, s, st, t, u,  &
                    v, xmh, xph, x2, z
INTEGER          :: j, n, n1
!---------------------

!     **** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!          SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0 .

eps = EPSILON(1.0)

!---------------------
fn_val = 0.0
IF (h == 0.0) RETURN

ah = ABS(h)
ax = ABS(x)
xph = ax + ah
xmh = ax - ah

t = MAX(ah,ax)
t = t * t
IF (1.6*t*t >= eps) THEN
  IF ((ax*ah)**2 < 0.5*eps) GO TO 60
  IF (ax > ah) THEN

    IF (xmh >= 9.0) THEN
      IF (xmh*xmh+p > -exparg(1)) RETURN
    END IF
    IF (4.0*ah*ax > -epsln()) GO TO 40

    IF (ax <= 3.0*ah) THEN
      IF (xph < 1.0) GO TO 30
      GO TO 50
    END IF
!-------------------------------------------------
!        FOR (AX LESS THAN OR EQUAL TO .40)
!-------------------------------------------------
    e = MAX(1.e-15,eps)
    IF (ax <= 0.40) THEN
      h2 = xph * xph
      a2 = xmh * xmh
      x2 = ax + ax
      st = 1.
      hf = xmh
      n = 0
      n1 = 1
      dn2 = 1.
      s = 0.
      10   n = n + 1
      n1 = n1 + 2
      dn2 = -dn2 / n
      st = h2 * st + x2 * hf
      hf = a2 * hf
      t = st * dn2 / n1
      s = s + t
      IF (ABS(t) > e*ABS(s)) GO TO 10
      s = 0.5 + (s+0.5)
      fn_val = 2.0 * c * ah * s
    ELSE
!-------------------------------------------------
!        FOR (AX GREATER THAN .40)
!-------------------------------------------------
      n = 1
      j = 0
      h2 = 0.
      z = EXP(-0.5*ax*ax)
      u = 2.0 * ah * c * z
      h3 = z
      v = 2.0 * h * h
      hf = 2.0 * ax * ah
      s = 0.
      20 h2 = (hf*h3-v*h2) / n
      n = n + 1
      h3 = (hf*h2-v*h3) / n
      n = n + 1
      hg = h3 / n
      s = s + hg
      IF (ABS(hg) > e*ABS(s)) GO TO 20
      IF (j == 0) THEN
        j = 1
        GO TO 20
      END IF
      fn_val = u * (s+z)
    END IF
    IF (h < 0.0) fn_val = -fn_val
    RETURN
  END IF
!-------------------------------------------------
!        SPECIAL CASES
!-------------------------------------------------
  IF (xph >= 5.8) THEN
    IF (xmh > -5.6) GO TO 40
    fn_val = SIGN(2.0,h)
    RETURN
  END IF

  30 fn_val = erf(xph) - erf(xmh)
  IF (h < 0.0) fn_val = -fn_val
  RETURN

  40 fn_val = erfc(xmh)
  IF (h < 0.0) fn_val = -fn_val
  RETURN

  50 fn_val = erfc(xmh) - erfc(xph)
  IF (h < 0.0) fn_val = -fn_val
  RETURN
END IF

fn_val = 2.0 * c * h * (0.5+(0.5-(x*x+h*h/3.0)))
RETURN

!     THE VALUE IS  2.0*EXP(-X*X)*ERF(H)

60 t = 2.0
x2 = x * x
IF (x2 >= eps) t = 2.0 * EXP(-x2)
IF (h*h < 3.0*eps) THEN
  fn_val = c * h * t
  RETURN
END IF
fn_val = t * erf(h)
RETURN
END FUNCTION aerf


FUNCTION daerf(x, h) RESULT(fn_val)
!-----------------------------------------------------------------------
!             COMPUTATION OF ERF(X + H) - ERF(X - H)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x, h
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: c = 1.12837916709551257389615890312155_dp,  &
                        p = 2.76959_dp
REAL (dp) :: ah, ax, a2, dn2, e, eps, hf, hg, h2, h3, n, n1, s, st,  &
             t, u, v, xmh, xph, x2, z
INTEGER   :: j
!---------------------
!     C = 2/SQRT(PI)
!     P = LN(9*SQRT(PI))
!---------------------

!     **** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!          SMALLEST NUMBER SUCH THAT 1._dp + EPS > 1._dp .

eps = EPSILON(1.0_dp)

!---------------------
fn_val = 0._dp
IF (h == 0._dp) RETURN

ah = ABS(h)
ax = ABS(x)
xph = ax + ah
xmh = ax - ah

t = MAX(ah,ax)
t = t * t
IF (1.6_dp*t*t >= eps) THEN
  IF ((ax*ah)**2 < 0.5_dp*eps) GO TO 60
  IF (ax > ah) THEN

    IF (xmh >= 9._dp) THEN
      IF (xmh*xmh+p > -dxparg(1)) RETURN
    END IF
    IF (4._dp*ah*ax > -depsln()) GO TO 40

    IF (ax <= 3._dp*ah) THEN
      IF (xph < 1._dp) GO TO 30
      GO TO 50
    END IF
!-------------------------------------------------
!        FOR (AX LESS THAN OR EQUAL TO .40)
!-------------------------------------------------
    e = MAX(1.d-30,eps)
    IF (ax <= 0.4_dp) THEN
      h2 = xph * xph
      a2 = xmh * xmh
      x2 = ax + ax
      st = 1._dp
      hf = xmh
      n = 0._dp
      n1 = 1._dp
      dn2 = 1._dp
      s = 0._dp
      10   n = n + 1._dp
      n1 = n1 + 2._dp
      dn2 = -dn2 / n
      st = h2 * st + x2 * hf
      hf = a2 * hf
      t = st * dn2 / n1
      s = s + t
      IF (ABS(t) > e*ABS(s)) GO TO 10
      s = 0.5_dp + (0.5_dp+s)
      fn_val = 2._dp * c * ah * s
    ELSE
!-------------------------------------------------
!        FOR (AX GREATER THAN .40)
!-------------------------------------------------
      n = 1._dp
      j = 0
      h2 = 0._dp
      z = EXP(-0.5_dp*ax*ax)
      u = 2._dp * ah * c * z
      h3 = z
      v = 2._dp * h * h
      hf = 2._dp * ax * ah
      s = 0._dp
      20   h2 = (hf*h3-v*h2) / n
      n = n + 1._dp
      h3 = (hf*h2-v*h3) / n
      n = n + 1._dp
      hg = h3 / n
      s = s + hg
      IF (ABS(hg) > e*ABS(s)) GO TO 20
      IF (j == 0) THEN
        j = 1
        GO TO 20
      END IF
      fn_val = u * (s+z)
    END IF
    IF (h < 0._dp) fn_val = -fn_val
    RETURN
  END IF
!-------------------------------------------------
!        SPECIAL CASES
!-------------------------------------------------
  IF (xph >= 8.5_dp) THEN
    IF (xmh > -8.3_dp) GO TO 40
    fn_val = SIGN(2._dp,h)
    RETURN
  END IF

  30 fn_val = derf(xph) - derf(xmh)
  IF (h < 0._dp) fn_val = -fn_val
  RETURN

  40 fn_val = derfc(xmh)
  IF (h < 0._dp) fn_val = -fn_val
  RETURN

  50 fn_val = derfc(xmh) - derfc(xph)
  IF (h < 0._dp) fn_val = -fn_val
  RETURN
END IF

fn_val = 2._dp * c * h * (0.5_dp+(0.5_dp-(x*x+h*h/3._dp)))
RETURN

!     THE VALUE IS  2.0*EXP(-X*X)*ERF(H)

60 t = 2._dp
x2 = x * x
IF (x2 >= eps) t = 2._dp * EXP(-x2)
IF (h*h < 3._dp*eps) THEN
  fn_val = c * h * t
  RETURN
END IF
fn_val = t * derf(h)
RETURN
END FUNCTION daerf


FUNCTION pndf(x, ind) RESULT(fn_val)
!     ---------------
!     A = 1/SQRT(2)
!     C = SQRT(2/PI)
!     ---------------
REAL, INTENT(IN)    :: x
INTEGER, INTENT(IN) :: ind
REAL                :: fn_val

! Local variables
REAL, PARAMETER :: a = .707106781186548, c = .797884560802865
REAL            :: t
!     ---------------
t = a * x
IF (ind == 0) THEN
  IF (x >= -8.0) THEN
    fn_val = 0.5 * erfc1(0,-t)
    RETURN
  END IF
  fn_val = c / erfc1(1,-t)
  RETURN
END IF
IF (x <= 8.0) THEN
  fn_val = 0.5 * erfc1(0,t)
  RETURN
END IF
fn_val = c / erfc1(1,t)
RETURN
END FUNCTION pndf


SUBROUTINE pni(p, q, d, w, ierr)
!-----------------------------------------------------------------------

!         EVALUATION OF THE INVERSE NORMAL DISTRIBUTION FUNCTION

!                           ------------

!     LET F(T) = 1/(SQRT(2*PI)*EXP(-T*T/2)). THEN THE FUNCTION

!        PROB(X) = INTEGRAL FROM MINUS INFINITY TO X OF F(T)

!     IS THE NORMAL DISTRIBUTION FUNCTION OF ZERO MEAN AND UNIT
!     VARIANCE. IT IS ASSUMED THAT P > 0, Q > 0, P + Q = 1,
!     AND D = P - 0.5. THE VALUE W IS COMPUTED WHERE PROB(W) = P.

!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.

!       IERR = 0  NO INPUT ERRORS WERE DETECTED. W WAS COMPUTED.
!       IERR = 1  EITHER P OR Q IS INCORRECT.
!       IERR = 2  D IS INCORRECT.

!-----------------------------------------------------------------------
!     RT2 = SQRT(2)
!------------------------
REAL, INTENT(IN)     :: p, q, d
REAL, INTENT(OUT)    :: w
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL, PARAMETER :: rt2 = 1.414213562373095
REAL            :: eps, t, u, v
!------------------------
t = MIN(p,q)
IF (t > 0.0) THEN
  eps = MAX(EPSILON(1.0),1.e-15)
  w = 0.5 + (0.5-(p+q))
  IF (ABS(w) <= 2.0*eps) THEN

    u = ABS(d+d)
    v = t + t
    w = erfi(u,v)
    IF (w < 0.0) GO TO 10

    ierr = 0
    w = rt2 * w
    IF (d < 0.0) w = -w
    RETURN
  END IF
END IF

!                         ERROR RETURN

ierr = 1
RETURN
10 ierr = 2
RETURN
END SUBROUTINE pni


SUBROUTINE dpni(p, q, d, w, ierr)
!-----------------------------------------------------------------------

!         EVALUATION OF THE INVERSE NORMAL DISTRIBUTION FUNCTION

!                           ------------

!     LET F(T) = 1/(SQRT(2*PI)*EXP(-T*T/2)). THEN THE FUNCTION

!        PROB(X) = INTEGRAL FROM MINUS INFINITY TO X OF F(T)

!     IS THE NORMAL DISTRIBUTION FUNCTION OF ZERO MEAN AND UNIT
!     VARIANCE. IT IS ASSUMED THAT P > 0, Q > 0, P + Q = 1,
!     AND D = P - 0.5. THE VALUE W IS COMPUTED WHERE PROB(W) = P.

!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.

!       IERR = 0  NO INPUT ERRORS WERE DETECTED. W WAS COMPUTED.
!       IERR = 1  EITHER P OR Q IS INCORRECT.
!       IERR = 2  D IS INCORRECT.

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: p, q, d
REAL (dp), INTENT(OUT) :: w
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL (dp), PARAMETER :: rt2 = 1.4142135623730950488016887242097_dp
REAL (dp)            :: eps, t, u, v
!------------------------
!     RT2 = SQRT(2)
!------------------------
t = MIN(p,q)
IF (t > 0._dp) THEN
  eps = EPSILON(1.0_dp)
  w = 0.5_dp + (0.5_dp-(p+q))
  IF (ABS(w) <= 2._dp*eps) THEN

    u = ABS(d+d)
    v = t + t
    w = derfi(u,v)
    IF (w < 0._dp) GO TO 10

    ierr = 0
    w = rt2 * w
    IF (d < 0._dp) w = -w
    RETURN
  END IF
END IF

!                         ERROR RETURN

ierr = 1
RETURN
10 ierr = 2
RETURN
END SUBROUTINE dpni


FUNCTION daw(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!     THIS FUNCTION COMPUTES SINGLE PRECISION VALUES OF DAWSONS INTEGRAL,

!        EXP(-X*X) * INTEGRAL (FROM 0 TO X) EXP(T*T) DT,

!     DEFINED FOR ALL REAL ARGUMENTS.

!     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
!     APPROXIMATIONS PUBLISHED IN MATH. COMP. 24, 171-178(1970) BY
!     CODY, PACIOREK AND THACHER.

!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: p1(9) = (/ .100000000000000E+01, -.135599049815353E+00,  &
        .456738974064825E-01, -.258323495918050E-02, .360079463580992E-03,  &
        -.944375029163387E-05, .634674256878843E-06, -.711645839183817E-08,  &
        .977985913592343E-10 /),  &
                   q1(9) = (/ .100000000000000E+01, .531067616851310E+00,  &
        .133052308640737E+00, .206907491644210E-01, .220437428972266E-02,  &
        .166706801664365E-03, .887964712053131E-05, .311750854173480E-06,  &
        .574807177698046E-08 /),  &
                   p2(8) = (/ -.150695651187161E+01, .293365747395449E+02, &
       -.400000893643550E+02, -.757931918089369E-0, -.889106479747812E+01,  &
        .152644099623699E+02, -.597678086823489E+01, .500236896088668E+00 /), &
                   q2(7) = (/ -.673106069744813E+00, .124486788262252E+04,  &
        .721193217600229E+01, .112461662024575E+03, .729177556415532E+02,  &
        .115840292551888E+03, .226064666074309E+00 /),  &
                   p3(8) = (/ .476405645273229E+01, -.266167674896399E+02,  &
        -.916804879813552E+01, -.150507703496692E+0, .506460153742231E+01,  &
        -.498544802986608E+0, -.149838042036691E+01, .499999902705054E+00 /),  &
                   q3(7) = (/ .287776122973187E+03, .256105722342226E+02,  &
        .751701277744067E+02, .146515167783109E+03, .330707724676114E+02,  &
        -.148715811787195E+01, .250011459611839E+00 /),  &
                   p4(7) = (/ -.315576735766984E+02, -.100791496592972E+02,  &
        -.710713709224200E+01, -.596879853243925E+01, -.449773645376092E+01,  &
        -.249999965398199E+01, .499999999999330E+00 /),  &
                   q4(6) = (/ .168874162155616E+03, .698280748271071E+01,  &
        -.213029621139181E+02, -.712157348463305E+0, -.250005973192356E+01,  &
         .750000000715687E+00 /),  &
        xlarge = 16777216.0, xsmall = .59604644775391E-07
REAL    :: frac, sump, sumq, w2, y
INTEGER :: i
!-----------------------------------------------------------------------

IF (ABS(x) <= xlarge) THEN
  IF (ABS(x) < xsmall) GO TO 40
  y = x * x
  IF (y < 6.25) THEN

!     ---------- ABS(X) < 2.5 ----------

    sump = (((((((p1(9)*y+p1(8))*y+p1(7))*y+p1(6))*y+p1(5))*y+  &
           p1(4))*y+p1(3))*y+p1(2)) * y + p1(1)
    sumq = (((((((q1(9)*y+q1(8))*y+q1(7))*y+q1(6))*y+q1(5))*y+  &
           q1(4))*y+q1(3))*y+q1(2)) * y + q1(1)
    fn_val = x * sump / sumq
    GO TO 50
  END IF

!     ---------- 2.5 <= ABS(X) < 3.5 ----------

  IF (y < 12.25) THEN
    frac = 0.0

    DO i = 1, 7
      frac = q2(i) / (p2(i)+y+frac)
    END DO

    fn_val = (p2(8)+frac) / x
    GO TO 50
  END IF

!     ---------- 3.5 <= ABS(X) < 5.0 ----------

  IF (y < 25.0) THEN
    frac = 0.0

    DO i = 1, 7
      frac = q3(i) / (p3(i)+y+frac)
    END DO

    fn_val = (p3(8)+frac) / x
    GO TO 50
  END IF

!     ---------- 5.0 <= ABS(X) <= XLARGE ----------

  w2 = 1.0 / x / x
  frac = 0.0

  DO i = 1, 6
    frac = q4(i) / (p4(i)+y+frac)
  END DO

  frac = p4(7) + frac
  fn_val = (0.5+0.5*w2*frac) / x
  GO TO 50
END IF

!     ---------- XLARGE < ABS(X) ----------

fn_val = 0.5 / x
GO TO 50

!     ---------- RETURN FOR SMALL X ----------

40 fn_val = x

50 RETURN
END FUNCTION daw


FUNCTION dpdaw(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            REAL (dp) COMPUTATION OF DAWSONS INTEGRAL
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: a(20) = (/ -.6666666666666666666666666666657D+00,  &
        .2666666666666666666666666665302D+00,  &
        -.7619047619047619047619046823290D-01,  &
        .1693121693121693121693097101950D-01,  &
        -.3078403078403078403073750370528D-02,  &
        .4736004736004736004148385001356D-03,  &
        -.6314672981339647953899064401849D-04,  &
        .7429027036870170692270376716931D-05,  &
        -.7820028459863171536925638117632D-06,  &
        .7447646152244351276445219666744D-07,  &
        -.6476214045244314289022051868963D-08,  &
        .5180971231894821888421654461823D-09,  &
        -.3837756389504541092817011259727D-10,  &
        .2646727414301012080897585412600D-11,  &
        -.1707553348198261876085879075486D-12,  &
        .1034770025122653524451023758330D-13,  &
        -.5905667147861158816695814259561D-15,  &
        .3157018166820009192834256496230D-16,  &
        -.1501742103181747984387915732309D-17,  &
        .4921379778280206677674574916266D-19 /)
REAL (dp), PARAMETER :: b(45) = (/ -.56886544105215527114160533733674D-01,  &
  -.31811346996168131279322878048822D+00, .20873845413642236789741580198858D+00,  &
  -.12475409913779131214073498314784D+00, .67869305186676777092847516423676D-01,  &
  -.33659144895270939503068230966587D-01, .15260781271987971743682460381640D-01,  &
  -.63483709625962148230586094788535D-02, .24326740920748520596865966109343D-02,  &
  -.86219541491065032038526983549637D-03, .28376573336321625302857636538295D-03,  &
  -.87057549874170423699396581464335D-04, .24986849985481658331800044137276D-04,  &
  -.67319286764160294344603050339520D-05, .17078578785573543710504524047844D-05,  &
  -.40917551226475381271896592490038D-06, .92828292216755773260751785312273D-07,  &
  -.19991403610147617829845096332198D-07, .40963490644082195241210487868917D-08,  &
  -.80032409540993168075706781753561D-09, .14938503128761465059143225550110D-09,  &
  -.26687999885622329284924651063339D-10, .45712216985159458151405617724103D-11,  &
  -.75187305222043565872243727326771D-12, .11893100052629681879029828987302D-12,  &
  -.18116907933852346973490318263084D-13, .26611733684358969193001612199626D-14,  &
  -.37738863052129419795444109905930D-15, .51727953789087172679680082229329D-16,  &
  -.68603684084077500979419564670102D-17, .88123751354161071806469337321745D-18,  &
  -.10974248249996606292106299624652D-18, .13261199326367178513595545891635D-19,  &
  -.15562732768137380785488776571562D-20, .17751425583655720607833415570773D-21,  &
  -.19695006967006578384953608765439D-22, .21270074896998699661924010120533D-23,  &
  -.22375398124627973794182113962666D-24, .22942768578582348946971383125333D-25,  &
  -.22943788846552928693329592319999D-26, .22391702100592453618342297600000D-27,  &
  -.21338230616608897703678225066666D-28, .19866196585123531518028458666666D-29,  &
  -.18079295866694391771955199999999D-30, .16090686015283030305450666666666D-31 /)
REAL (dp) :: ax, eps = 1.d-31, t, w
INTEGER   :: i, k
!----------------------------
ax = ABS(x)
IF (ax < 4._dp) THEN
  t = x * x
  IF (ax <= 1._dp) THEN

!                       ABS(X) <= 1

    fn_val = x
    IF (t < eps) RETURN

    w = a(20)
    DO i = 1, 19
      k = 20 - i
      w = t * w + a(k)
    END DO
    fn_val = x * (0.75_dp+(0.25_dp+t*w))
    RETURN
  END IF

!                    1 < ABS(X) < 4

  fn_val = x * (.25_dp+dcsevl(.125_dp*t-1._dp,b,45))
  RETURN
END IF

!                       ABS(X) >= 4

fn_val = dpdaw0(ax) / x
RETURN
END FUNCTION dpdaw


FUNCTION dpdaw0(x) RESULT(fn_val)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

!                EVALUATION OF X*DAW(X) FOR X >= 4
!                WHERE DAW(X) IS THE DAWSON INTEGRAL

!-----------------------------------------------------------------------
REAL (dp) :: t, u, v, w
REAL (dp) :: a0 = .59682223611279114961181337D-06, a1  &
        = .17685355947137064277328544D-05, a2  &
        = .46539151619719879425847199D-05, a3  &
        = .60206549750426063518629015D-05, a4  &
        = .71968323029042065431569341D-05, a5  &
        = .30632314265730271259576310D-05, a6  &
        = .58023977792358623878717970D-06, a7  &
        = -.79009104459686847104040749D-06, a8  &
        = -.48136647436848802449955585D-06, a9  &
        = .15516701696125663593787151D-05, a10  &
        = .41907470564012404920368069D-05, a11  &
        = .58027971313506864533128271D-05, a12  &
        = .55370618938769991926278955D-05, a13  &
        = .38921890500515083447099834D-05, a14  &
        = .20373391058442140223632125D-05, a15  &
        = .77894504402267231707108862D-06, a16  &
        = .20671717275442647450228015D-06, a17  &
        = .34163346266402336952767687D-07, a18  &
        = .26550833614486808486653527D-08, b1  &
        = .24655153970246619491575782D+01, b2  &
        = .62943144539382480033771054D+01, b3  &
        = .62122304537339562619049238D+01, b4 = .72256247680648550388609993D+01
REAL (dp) :: c0 = .59682220964671964619429990D-06, c1  &
        = .15714896176335781761895194D-05, c2  &
        = .38969553030269966874001463D-05, c3  &
        = .52383725295173020317044495D-05, c4  &
        = .68105053752194312932660379D-05, c5  &
        = .54424292166005622934118953D-05, c6  &
        = .42426204827598626398351706D-05, c7  &
        = .16034175582037599078896885D-05, c8  &
        = .40869743822567069272296881D-06, c9  &
        = -.28877279517645392833652391D-06, c10  &
        = -.23330912423546272053845073D-06, c11  &
        = -.52391219684126122047201575D-07, c12  &
        = .21002079856304960793765863D-07, c13  &
        = .15999459808634652894884878D-07, c14  &
        = .30942877932880022278704223D-08, d1  &
        = .21353586303464124681525146D+01, d2  &
        = .51903635881846225790635588D+01, d3  &
        = .55426108503039000945912263D+01, d4  &
        = .72404445472735688881512590D+01, d5  &
        = .42043425514299588216548098D+01, d6  &
        = .39555198251457455938284071D+01, d7  &
        = .96558431899415447517392039D+00, d8 = .85802563214754700974536646D+00
REAL (dp) :: e0 = .59682220964671964619429989D-06, e1  &
        = .84093250039847132132545694D-06, e2  &
        = .21215362314612857942015465D-05, e3  &
        = .17544801458346344458699925D-05, e4  &
        = .23300467572238537581015267D-05, e5  &
        = .10789302963153438789721000D-05, e6  &
        = .84593690165753140097573765D-06, e7  &
        = .58921927492401852305814815D-07, e8  &
        = .88699252882391671775568852D-08, e9  &
        = -.45367333700483839577191791D-07, e10  &
        = -.57476692572774473717756092D-08, f1  &
        = .91128032791938035312580060D+00, f2  &
        = .28248447610354665352098820D+01, f3  &
        = .12208226691968579661669930D+01, f4  &
        = .26127044898109340807106578D+01, f5  &
        = .34850941500747816686718301D+00, f6  &
        = .11058424028893031129155960D+01, f7  &
        = -.58583432304677947215235382D-01, f8  &
        = .22235753204903993496120223D+00, f9  &
        = -.22116625272550591073419658D-01, f10  &
        = .14938281115881437851795141D-01
REAL (dp) :: g0 = .59682220962853689178993039D-06, g1  &
        = .10399202637216492921747166D-05, g2  &
        = .20011159669496258291170369D-05, g3  &
        = .19689301951834731308289277D-05, g4  &
        = .18861506874082130108404839D-05, g5  &
        = .10403672611266773427898830D-05, g6  &
        = .42700287233074049241986613D-06, g7  &
        = -.16199303908856237550380048D-07, g8  &
        = -.73900110466591484816182863D-07, g9  &
        = -.21369058482679700785965255D-07, g10  &
        = .42853587024761116737793072D-08, g11  &
        = -.27226352610679391576406666D-09, h1  &
        = .12446924554654251972427023D+01, h2  &
        = .24571243397858142676189669D+01, h3  &
        = .16710488988253927773151334D+01, h4  &
        = .17260859791081983044753913D+01, h5  &
        = .67406536368694046314885196D+00, h6  &
        = .47403526893401885333000527D+00, h7  &
        = .96251479080923959509108658D-01, h8  &
        = .53334517885765587426678624D-01, h9  &
        = .39450307003297031975491216D-02, h10  &
        = .20209221166462656808887976D-02
REAL (dp) :: s0 = .8210986449041747719684610504710D-02, s1  &
        = .8646073144815053170065230898334D-02, s2  &
        = .4768322737615973285410030924398D-03, s3  &
        = .4792593707378225992657970685396D-04, s4  &
        = .7507677744363576693833551190204D-05, s5  &
        = .1737929446861228512373840727547D-05
REAL (dp) :: p0 = .29531250000000000000002D+02, p1  &
        = -.14571781607273299440392D+04, p2 = .24285318385898860175073D+05,  &
        p3 = -.15843555052114168113822D+06, p4  &
        = .32969397422638395586636D+06, p5 = -.55331506994311967089636D+04,  &
        q1 = -.54843599093412230906873D+02, q2  &
        = .10882497826844164906477D+04, q3 = -.96578534881552358457185D+04,  &
        q4 = .37803382357862589384458D+05, q5 = -.51283783372259864777146D+05
!----------------------------
IF (x < 12._dp) THEN
  t = (32._dp/(x*x)-0.5_dp) - 0.5_dp
  IF (t < 0._dp) THEN

!                    -7/9 <= T <= -0.4

    IF (t <= -0.4_dp) THEN
      u = ((((((((a18*t + a17)*t + a16)*t + a15)*t + a14)*t + a13)*t + a12)*t+  &
      a11)*t + a10) * t + a9
      u = ((((((((u*t + a8)*t + a7)*t + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1) * t + a0
      v = (((b4*t + b3)*t + b2)*t + b1) * t + 1._dp
      GO TO 10
    END IF

!                     -0.4 < T < 0

    u = ((((((((c14*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*  &
    t + c6) * t + c5
    u = ((((u*t + c4)*t + c3)*t + c2)*t + c1) * t + c0
    v = (((((((d8*t+d7)*t+d6)*t+d5)*t+d4)*t+d3)*t+d2)*t+d1) * t +1._dp
  ELSE

!                      0 <= T <= 0.4

    IF (t <= 0.4_dp) THEN
      u = (((((((((e10*t + e9)*t + e8)*t + e7)*t + e6)*t + e5)*t + e4)*t + e3)*t  &
      +e2)*t + e1) * t + e0
      v = (((((((((f10*t+f9)*t+f8)*t+f7)*t+f6)*t+f5)*t+f4)*t+f3)*t  &
      +f2)*t+f1) * t + 1._dp
    ELSE

!                      0.4 < T <= 1

      u = ((((((((g11*t+g10)*t+g9)*t+g8)*t+g7)*t+g6)*t+g5)*t+g4)*t+g3) * t + g2
      u = (u*t+g1) * t + g0
      v = (((((((((h10*t+h9)*t+h8)*t+h7)*t+h6)*t+h5)*t+h4)*t+h3)*t  &
      +h2)*t+h1) * t + 1._dp
    END IF
  END IF

!          THE ABOVE FOUR MINIMAX APPROXIMATIONS U/V
!          ARE ACCURATE TO WITHIN 1 UNIT OF THE 25-TH
!          SIGNIFICANT DIGIT. THUS, THE APPROXIMATION
!          FOR W IS ACCURATE TO WITHIN 1 UNIT OF THE
!          29-TH SIGNIFICANT DIGIT.

  10 w = ((((((u/v)*t + s5)*t + s4)*t + s3)*t + s2)*t + s1) * t + s0
  fn_val = 0.5_dp + w
  RETURN
END IF

!                          X >= 12

t = (1._dp/x) ** 2
w = (((((p5*t + p4)*t + p3)*t + p2)*t + p1)*t + p0) / (((((q5*t+q4)*t+q3)*t+  &
q2)*t+q1)*t+1._dp)
w = ((((w*t+6.5625_dp)*t+1.875_dp)*t+0.75_dp)*t+0.5_dp) * t + 1._dp
fn_val = 0.5_dp * w
RETURN
END FUNCTION dpdaw0


SUBROUTINE cfrnli(mo, z, w)
!-----------------------------------------------------------------------

!            COMPUTATION OF THE COMPLEX FRESNEL INTEGRAL E(Z)

!                           ----------------

!                      W = E(Z)          IF MO = 0
!                      W = EXP(-Z)*E(Z)  OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)  :: mo
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: cd(18) = (/ 0.00000000000000, 2.08605856013476E-2,  &
        8.29806940495687E-2, 1.85421653326079E-1, 3.27963479382361E-1,  &
        5.12675279912828E-1, 7.45412958045105E-1, 1.03695067418297,  &
        1.40378061255437, 1.86891662214001, 2.46314830523929,  &
        3.22719383737352, 4.21534348280013, 5.50178873151549,  &
        7.19258966683102, 9.45170208076408, 1.25710718314784E+1,  &
        1.72483537216334E+1 /), ce(18) = (/ 8.15723083324096E-2,  &
        1.59285285253437E-1, 1.48581625614499E-1, 1.33219670836245E-1,  &
        1.15690392878957E-1, 9.78580959447535E-2, 8.05908834297624E-2,  &
        6.40204538609872E-2, 4.81445242767885E-2, 3.33540658473295E-2,  &
        2.05548099470193E-2, 1.07847403887506E-2, 4.55634892214219E-3,  &
        1.43984458138925E-3, 3.07056139834171E-4, 3.78156541168541E-5,  &
        2.05173509616121E-6, 2.63564823682747E-8 /),  &
        c = .564189583547756, c0 = -.707106781186548
REAL    :: qf(2), sm(2), tm(2), ts(2), zr(2), dm, pm, qm, r, ss, x, y
INTEGER :: i
!------------------------
!     C = 1/SQRT(PI)
!     C0 = -1/SQRT(2)
!------------------------
x = REAL(z)
y = AIMAG(z)
r = cpabs(x,y)
IF (r == 0.0) GO TO 110

!              EVALUATION OF ZR = SQRT(2*Z/PI)

IF (x < 0.0) THEN
  zr(2) = SQRT(r-x)
  zr(1) = y / zr(2)
ELSE
  zr(1) = SQRT(r+x)
  IF (y < 0.0) zr(1) = -zr(1)
  zr(2) = y / zr(1)
END IF
zr(1) = c * zr(1)
zr(2) = c * zr(2)

IF (r > 1.0) THEN
  IF (r >= 38.0) GO TO 50
  IF (x < 0.016*y*y) GO TO 30
END IF

!                       TAYLOR SERIES

sm(1) = 0.0
sm(2) = 0.0
tm(1) = zr(1)
tm(2) = zr(2)
pm = 0.0
10 pm = pm + 1.0
dm = 2.0 * pm + 1.0
ts(1) = tm(1) * x - tm(2) * y
ts(2) = tm(1) * y + tm(2) * x
tm(1) = ts(1) / pm
tm(2) = ts(2) / pm
ts(1) = tm(1) / dm
ts(2) = tm(2) / dm
IF (ABS(sm(1))+ABS(ts(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(ts(2)) == ABS(sm(2))) GO TO 20
END IF
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
GO TO 10
20 sm(1) = zr(1) + sm(1)
sm(2) = (c0+zr(2)) + sm(2)

IF (mo == 0) GO TO 100
qm = EXP(-x)
qf(1) = qm * COS(-y)
qf(2) = qm * SIN(-y)
GO TO 90

!              RATIONAL FUNCTION APPROXIMATION

30 sm(1) = 0.0
sm(2) = 0.0
DO i = 1, 18
  ts(1) = x - cd(i)
  ts(2) = y
  ss = ts(1) * ts(1) + ts(2) * ts(2)
  tm(1) = ce(i) * ts(1) / ss
  tm(2) = -ce(i) * ts(2) / ss
  sm(1) = sm(1) + tm(1)
  sm(2) = sm(2) + tm(2)
END DO
ts(1) = zr(1) * sm(1) - zr(2) * sm(2)
ts(2) = zr(1) * sm(2) + zr(2) * sm(1)
sm(1) = 0.5 * ts(1)
sm(2) = 0.5 * ts(2)
GO TO 80

!                   ASYMPTOTIC EXPANSION

50 qf(1) = (x/r) / r
qf(2) = -(y/r) / r
tm(1) = qf(1)
tm(2) = qf(2)
sm(1) = tm(1)
sm(2) = tm(2)
pm = -0.5
60 pm = pm + 1.0
ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
tm(1) = pm * ts(1)
tm(2) = pm * ts(2)
IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 70
END IF
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (pm < 25.5) GO TO 60
70 ts(1) = zr(1) * sm(1) - zr(2) * sm(2)
ts(2) = zr(1) * sm(2) + zr(2) * sm(1)
sm(1) = 0.5 * ts(1)
sm(2) = 0.5 * ts(2)
IF (zr(2) < 8.e-3) GO TO 120

!                       TERMINATION

80 IF (mo /= 0) GO TO 100
qm = EXP(x)
qf(1) = qm * COS(y)
qf(2) = qm * SIN(y)

90 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)

100 w = CMPLX(sm(1),sm(2))
RETURN

!                      CASE WHEN Z = 0

110 w = CMPLX(0.0,c0)
RETURN

!               MODIFIED ASYMPTOTIC EXPANSION

120 IF (mo == 0) THEN
  qm = EXP(x)
  qf(1) = qm * COS(y)
  qf(2) = qm * SIN(y)
  ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
  ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
  w = CMPLX(ts(1),c0+ts(2))
  RETURN
END IF

IF (-x <= exparg(1)) GO TO 100
qm = c0 * EXP(-x)
sm(1) = sm(1) + qm * SIN(y)
sm(2) = sm(2) + qm * COS(y)
w = CMPLX(sm(1),sm(2))
RETURN
END SUBROUTINE cfrnli


SUBROUTINE frnl(t, c, s)
!-----------------------------------------------------------------------
!             EVALUATION OF THE REAL FRESNEL INTEGRALS
!-----------------------------------------------------------------------
REAL, INTENT(IN)  :: t
REAL, INTENT(OUT) :: c, s

! Local variables
REAL, PARAMETER :: pi = 3.1415926535898
REAL, PARAMETER :: a(6) = (/ -.119278241233760E-05, .540730666359417E-04,  &
        -.160488306381990E-02, .281855008757077E-01, -.246740110027210E+00,  &
         .100000000000000E+01 /),  &
                   b(6) = (/ -.155653074871090E-06, .844415353045065E-05,  &
        -.312116934326082E-03, .724478420395276E-02, -.922805853580325E-01, &
         .523598775598300E+00 /),  &
                 cp(13) = (/ .114739945188034E-20, -.384444827287950E-18,  &
         .832125729394275E-16, -.142979507360076E-13, .198954961821465E-11,  &
        -.220226545457144E-09, .188434924092257E-07, -.120009722914157E-05,  &
         .540741337442140E-04, -.160488313553028E-02, .281855008777956E-01,  &
        -.246740110027196E+0, .999999999999996E+00 /),  &
                 sp(13) = (/ .705700784853927E-22, -.252757991492418E-19,  &
         .594117488940008E-17, -.112161631555448E-14, .173332189994074E-12,  &
        -.215742302078015E-10, .210821173208116E-08, -.156471443116560E-06,  &
         .844427287845253E-05, -.312116942346186E-03, .724478420418951E-02,  &
        -.922805853580323E-0, .523598775598300E+00 /)
REAL, PARAMETER :: pn(6) = (/ .318309816100920E+00, .134919391391516E+02,  &
        .158258097490377E+03, .598796451682535E+03, .632369782194966E+03,  &
        .967985390141920E+02 /), pd(6) = (/ .100000000000000E+01,  &
        .426900960480796E+02, .509085485682426E+03, .200034664144742E+04,  &
        .231910140792937E+04, .486678558201084E+03 /),  &
                   qn(6) = (/ .101320876178478E+00, .490534697099052E+01,  &
        .652095157811808E+02, .274183825747887E+03, .305040725009211E+03,  &
        .364566615872326E+02 /),  &
                   qd(6) = (/ .100000000000000E+01,  &
        .499330024470621E+02, .709854097670206E+03, .343470762861172E+04,  &
        .522213879312684E+04, .168801831831851E+04 /)
REAL, PARAMETER :: an(6) = (/ .318309885869756E+00, .254179177393500E+02,  &
        .575003792540838E+03, .426673405867140E+04, .891831887923938E+04,  &
        .267955736537967E+04 /),  &
                   ad(6) = (/ .100000000000000E+01,  &
        .801567066285184E+02, .182971463354850E+04, .138848884373420E+05,  &
        .309228411873207E+05, .120421274105856E+05 /),  &
                   bn(6) = (/ .101321181932417E+00, .925021984290547E+01,  &
        .240932023056602E+03, .206079616836437E+04, .484901973010149E+04,  &
        .130680669688315E+04 /),  &
                   bd(6) = (/ .100000000000000E+01,  &
        .928158182389149E+02, .250926840439955E+04, .233924458152954E+05,  &
        .685638896406835E+05, .418593101455019E+05 /)
REAL, PARAMETER :: cn(5) = (/ .318309886182000E+00, .299191968327887E+02,  &
        .691428839605668E+03, .394539800974744E+04, .290314254767015E+04 /),  &
        cd(5) = (/ .100000000000000E+01, .942978925136851E+02,  &
        .219977296283666E+04, .129726479671006E+05, .114991427758165E+05 /),  &
        dn(5) = (/ .101321183630876E+00, .110988033615242E+02,  &
        .306282306497228E+03, .213130259794164E+04, .171270676541694E+04 /),  &
        dd(5) = (/ .100000000000000E+01, .111060616085627E+03,  &
        .318197586347414E+04, .249342095714049E+05, .359241903823488E+05 /)
REAL, PARAMETER :: fp(7) = (/ .449763389301234E+05, -.188763642051836E+04,  &
        .669261097103246E+02, -.343966606879114E+01, .343112896133346E+00,  &
        -.967546019461500E-01, .318309886183465E+00 /),  &
                   gp(7) = (/ .316642183365360E+06, -.120618995106638E+05,  &
        .359164749179351E+03, -.142252603258172E+02, .982934118445454E+00,  &
        -.153989722912325E+00, .101321183639714E+00 /),  &
        p(6) = (/ -654729075.0, 2027025.0, -10395.0, 105.0, -3.0, 1.0 /),  &
        q(6) = (/ -13749310575.0, 34459425.0, -135135.0, 945.0, -15.0, 1.0 /)
INTEGER :: max, i, l, m
REAL    :: cy, f, fn, fd, g, gn, gd, n, pix, pixx, r, sy, x, xx, y
!--------------------------

!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                     MAX = IPMPAR(3)
max = HUGE(3)

!-----------------------------------------------------------------------
x = ABS(t)
IF (x <= 4.0) THEN
  xx = x * x
  y = xx * xx
!-----------------------------------------------------------------------
!             EVALUATION OF C(X) AND S(X) FOR X < 1.65
!                          WHERE X = ABS(T)
!-----------------------------------------------------------------------
  IF (x <= 0.6) THEN
    c = ((((a(1)*y+a(2))*y+a(3))*y+a(4))*y+a(5)) * y + a(6)
    s = ((((b(1)*y+b(2))*y+b(3))*y+b(4))*y+b(5)) * y + b(6)
    c = t * c
    s = t * xx * s
    RETURN
  END IF

  IF (x < 1.65) THEN
    c = cp(1)
    s = sp(1)
    DO i = 2, 13
      c = cp(i) + c * y
      s = sp(i) + s * y
    END DO
    c = t * c
    s = t * xx * s
    RETURN
  END IF
!-----------------------------------------------------------------------
!          EVALUATION OF THE AUXILIARY FUNCTIONS F(X) AND G(X)
!                        FOR X >= 1.65
!-----------------------------------------------------------------------
  IF (x < 2.0) THEN
    fn = ((((pn(1)*y+pn(2))*y+pn(3))*y+pn(4))*y+pn(5)) * y +pn(6)
    fd = ((((pd(1)*y+pd(2))*y+pd(3))*y+pd(4))*y+pd(5)) * y +pd(6)
    gn = ((((qn(1)*y+qn(2))*y+qn(3))*y+qn(4))*y+qn(5)) * y +qn(6)
    gd = ((((qd(1)*y+qd(2))*y+qd(3))*y+qd(4))*y+qd(5)) * y +qd(6)
    f = fn / (x*fd)
    g = gn / (x*xx*gd)
    y = 0.5 * xx
    GO TO 30
  END IF

  IF (x < 3.0) THEN
    fn = ((((an(1)*y+an(2))*y+an(3))*y+an(4))*y+an(5)) * y +an(6)
    fd = ((((ad(1)*y+ad(2))*y+ad(3))*y+ad(4))*y+ad(5)) * y +ad(6)
    gn = ((((bn(1)*y+bn(2))*y+bn(3))*y+bn(4))*y+bn(5)) * y +bn(6)
    gd = ((((bd(1)*y+bd(2))*y+bd(3))*y+bd(4))*y+bd(5)) * y +bd(6)
    f = fn / (x*fd)
    g = gn / (x*xx*gd)
    GO TO 20
  END IF

  fn = (((cn(1)*y+cn(2))*y+cn(3))*y+cn(4)) * y + cn(5)
  fd = (((cd(1)*y+cd(2))*y+cd(3))*y+cd(4)) * y + cd(5)
  gn = (((dn(1)*y+dn(2))*y+dn(3))*y+dn(4)) * y + dn(5)
  gd = (((dd(1)*y+dd(2))*y+dd(3))*y+dd(4)) * y + dd(5)
  f = fn / (x*fd)
  g = gn / (x*xx*gd)
ELSE

  IF (x < 6.0) THEN
    xx = x * x
    y = 1.0 / (xx*xx)
    f = (((((fp(1)*y+fp(2))*y+fp(3))*y+fp(4))*y+fp(5))*y+fp(6)) *y + fp(7)
    g = (((((gp(1)*y+gp(2))*y+gp(3))*y+gp(4))*y+gp(5))*y+gp(6)) *y + gp(7)
    f = f / x
    g = g / (x*xx)
  ELSE

    IF (x >= REAL(MAX)) GO TO 40
    pix = pi * x
    pixx = pix * x
    y = 1.0 / pixx
    y = y * y
    f = ((((p(1)*y+p(2))*y+p(3))*y+p(4))*y+p(5)) * y + p(6)
    g = ((((q(1)*y+q(2))*y+q(3))*y+q(4))*y+q(5)) * y + q(6)
    f = f / pix
    g = g / (pix*pixx)
  END IF
END IF
!-----------------------------------------------------------------------
!           EVALUATION OF SIN(0.5*PI*X*X) AND COS(0.5*PI*X*X)
!                 THE RESULTS ARE STORED IN SY AND CY
!-----------------------------------------------------------------------
20 m = x
l = MOD(m,2)
n = m - l
y = x - m
r = x - n

y = y * n
m = y
y = y - m
IF (MOD(m,2) /= 0) y = (y-0.5) - 0.5
y = y + 0.5 * r * r

30 sy = sin1(y)
cy = cos1(y)
!-----------------------------------------------------------------------
!                             TERMINATION
!-----------------------------------------------------------------------
c = 0.5 + (f*sy-g*cy)
s = 0.5 - (f*cy+g*sy)
IF (t >= 0.0) RETURN
c = -c
s = -s
RETURN

40 IF (t >= 0.0) THEN
  c = 0.5
  s = 0.5
  RETURN
END IF
c = -0.5
s = -0.5
RETURN
END SUBROUTINE frnl


SUBROUTINE cexpli(mo, z, w)
!-----------------------------------------------------------------------
!           EVALUATION OF THE COMPLEX EXPONENTIAL INTEGRAL
!-----------------------------------------------------------------------
INTEGER, INTENT(IN)  :: mo
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: pi = 3.14159265358979, euler = .577215664901533
REAL, PARAMETER :: cd(18) = (/ 0.00000000000000E+00, .311105957086528E-01,  &
        .103661260539112E+00, .216532335244554E+00, .369931427960192E+00,  &
        .566766259990589E+00, .814042066324748E+00, .112384247540813E+01,  &
        .151400478148512E+01, .200886795032284E+01, .264052411823592E+01,  &
        .345098449933392E+01, .449583360763202E+01, .585058263409822E+01,  &
        .762273501463380E+01, .997814501584578E+01, .132122064896408E+02,  &
        .180322948376021E+02 /), ce(18) = (/ .850156516121093E-02,  &
        .505037465849058E-01, .836817368956407E-01, .107047582417607E+00,  &
        .120424719029462E+00, .125096631582229E+00, .122314435224685E+00,  &
        .112621417553907E+00, .963419407392582E-01, .747398422757511E-01,  &
        .508596135953441E-01, .290822706773628E-01, .132201640530101E-01,  &
        .443802939829067E-02, .992612478987576E-03, .126579795112011E-03,  &
        .702150908253350E-05, .910281532564632E-07 /)
REAL    :: qf(2), sm(2), tm(2), ts(2), g0(2), gn(2), h0(2), hn(2), wn(2),  &
           c, cy, d, e, eps, n, np1, qm, r, ss, sy, tol, u, x, y
INTEGER :: i
LOGICAL :: ind

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1._dp + EPS > 1._dp.

eps = EPSILON(1.0)

!-------------------------

x = REAL(z)
y = AIMAG(z)
r = cpabs(x,y)
eps = MAX(eps,1.e-15)

IF (r > 1.0) THEN
  IF (r >= 40.0) GO TO 40
  IF (r >= 4.0) THEN
    IF (x <= 0.0.OR.ABS(y) > 8.0) GO TO 40
    IF (r < 10.0 .AND. ABS(y) > 1.8*x) GO TO 40
  ELSE
    IF (x < 0.09*y*y) GO TO 20
    IF (r > 3.6 .AND. ABS(y) > 1.8*x) GO TO 40
  END IF
END IF

!                        TAYLOR SERIES

sm(1) = 0.0
sm(2) = 0.0
tm(1) = x
tm(2) = y
n = 1.0
10 n = n + 1.0
ts(1) = tm(1) * x - tm(2) * y
ts(2) = tm(1) * y + tm(2) * x
tm(1) = ts(1) / n
tm(2) = ts(2) / n
ts(1) = tm(1) / n
ts(2) = tm(2) / n
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
IF (anorm(ts(1),ts(2)) > eps*anorm(sm(1),sm(2))) GO TO 10
sm(1) = x + sm(1)
sm(2) = y + sm(2)

sm(1) = (euler+LOG(r)) + sm(1)
sm(2) = ATAN2(-y,-x) + sm(2)
GO TO 70

!                      RATIONAL EXPANSION

20 sm(1) = 0.0
sm(2) = 0.0
DO i = 1, 18
  ts(1) = x - cd(i)
  ts(2) = y
  ss = ts(1) * ts(1) + ts(2) * ts(2)
  sm(1) = sm(1) + ce(i) * ts(1) / ss
  sm(2) = sm(2) - ce(i) * ts(2) / ss
END DO
GO TO 60

!         PADE APPROXIMATION FOR THE ASYMPTOTIC EXPANSION
!                       FOR EXP(-Z)*EI(Z)

40 x = -x
y = -y
d = 4.0 * r
IF (r < 10.0) d = 32.0
g0(1) = 1.0
g0(2) = 0.0
gn(1) = (1.0+x) / d
gn(2) = y / d
h0(1) = 1.0
h0(2) = 0.0
u = x + 2.0
hn(1) = u / d
hn(2) = gn(2)
w = CMPLX(1.0+x,y) / CMPLX(u,y)
wn(1) = REAL(w)
wn(2) = AIMAG(w)
np1 = 1.0
tol = 4.0 * eps

50 n = np1
np1 = n + 1.0
e = (n*np1) / d
u = u + 2.0
tm(1) = ((u*gn(1)-y*gn(2))-e*g0(1)) / d
tm(2) = ((u*gn(2)+y*gn(1))-e*g0(2)) / d
g0(1) = gn(1)
g0(2) = gn(2)
gn(1) = tm(1)
gn(2) = tm(2)
tm(1) = ((u*hn(1)-y*hn(2))-e*h0(1)) / d
tm(2) = ((u*hn(2)+y*hn(1))-e*h0(2)) / d
h0(1) = hn(1)
h0(2) = hn(2)
hn(1) = tm(1)
hn(2) = tm(2)

tm(1) = wn(1)
tm(2) = wn(2)
w = CMPLX(gn(1),gn(2)) / CMPLX(hn(1),hn(2))
wn(1) = REAL(w)
wn(2) = AIMAG(w)
IF (anorm(tm(1)-wn(1),tm(2)-wn(2)) > tol*anorm(wn(1),wn(2)))  &
GO TO 50

x = REAL(z)
y = AIMAG(z)
w = w / z
sm(1) = REAL(w)
sm(2) = AIMAG(w)

!                         TERMINATION

60 ind = x <= 0.0 .OR. ABS(y) > 1.e-2
IF (ind .AND. mo /= 0) GO TO 90
c = pi
IF (y > 0.0) c = -pi
qm = EXP(x)
cy = COS(y)
sy = SIN(y)
qf(1) = qm * cy
qf(2) = qm * sy
IF (mo == 0) GO TO 80

r = c / qm
sm(1) = sm(1) + r * sy
sm(2) = sm(2) + r * cy
GO TO 90

70 IF (mo == 0) GO TO 90
ind = .true.
qm = EXP(-x)
qf(1) = qm * COS(-y)
qf(2) = qm * SIN(-y)

80 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)
IF (.NOT.ind) sm(2) = sm(2) + c

90 w = CMPLX(sm(1),sm(2))
RETURN
END SUBROUTINE cexpli


SUBROUTINE expli(int, arg, result, ierr)

!-----------------------------------------------------------------------

!     THIS SUBROUTINE COMPUTES THE EXPONENTIAL INTEGRALS

!            EI(X),  E-SUB-1(X) = -EI(-X),  AND  EXP(-X)*EI(X)

!     WHERE

!             INTEGRAL (FROM T=-INFINITY TO T=X) (EXP(T)/T),  X > 0,
!     EI(X) =
!             -INTEGRAL (FROM T=-X TO T=INFINITY) (EXP(-T)/T), X < 0,

!     AND WHERE THE FIRST INTEGRAL IS A PRINCIPAL VALUE INTEGRAL. THE
!     ARGUMENTS INT, ARG, AND RESULT HAVE THE FOLLOWING USAGE ...

!            INT              ARG             RESULT
!             1            X .NE. 0          EI(X)
!             2            X > 0          E-SUB-1(X)
!             3            X .NE. 0          EXP(-X)*EI(X)

!     THE EXPANSION FOR 4 <= X <= 8 IS DUE TO WAYNE FULLERTON (LOS
!     ALAMOS). THE REMAINING EXPANSIONS ARE FROM MATH. COMP. 22, 641-649
!     (1968), AND MATH. COMP. 23, 289-303(1969) BY CODY AND THACHER.

!        ------------

!     ERROR MONITORING

!        THE PARAMETER IERR IS A VARIABLE THAT IS SET BY THE ROUTINE.
!     IF NO ERRORS ARE DETECTED THEN IERR IS SET TO 0. THE FOLLOWING
!     TABLE INDICATES THE TYPES OF ERRORS THAT MAY BE ENCOUNTERED IN
!     THE ROUTINE AND THE FUNCTION VALUES SUPPLIED IN EACH CASE.

!     IERR    ERROR     ARGUMENT          FUNCTION VALUES FOR
!                        RANGE        EI(X)  EXP(-X)*EI(X)  E-SUB-1(X)
!      1    UNDERFLOW   X < XMIN     0          -          0
!      2    OVERFLOW    X > XMAX     T          -          -
!      3    ILLEGAL X     X = 0         T          T          T
!      4    ILLEGAL X    X < 0       -          -          T

!     T INDICATES THAT THE ROUTINE TERMINATES WITHOUT ASSIGNING A VALUE
!     TO THE FUNCTION.

!        ----------

!     THIS SUBROUTINE WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR
!     THE FUNPACK PACKAGE OF SPECIAL FUNCTION SUBROUTINES. THE ROUTINE
!     WAS MODIFIED BY BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------

!     XMAX AND XMIN ARE MACHINE DEPENDENT CONSTANTS FOR DETECTING
!     UNDERFLOW AND OVERFLOW. XMAX AND XMIN ARE GIVEN APPROXIMATE
!     VALUES IN STATEMENTS 240 AND 340.

!-----------------------------------------------------------------------

!     VALUE OF EXP(40.0)

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)  :: int
REAL, INTENT(IN)     :: arg
REAL, INTENT(OUT)    :: result
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL, PARAMETER :: a(6) = (/ -.577215664901531E+00, .758833087029943E+00,  &
        .125660818982053E+00, .204158408934305E-01, .825035122466538E-03,  &
        .962949813453924E-05 /), b(5) = (/ .100000000000000E+01,  &
        .417810755380398E+00, .730228560396799E-01, .642720224671078E-02,  &
        .245134203588369E-03 /), c(7) = (/ .465627107975096E-06,  &
        .999979577051595E+00, .904161556946328E+01, .243784088791317E+02,  &
        .230192559391334E+02, .690522522784443E+01, .430967839469389E+00 /)
REAL, PARAMETER :: d(7) = (/ .100000000000000E+01, .100411643829054E+02,  &
        .324264210695138E+02, .412807841891424E+02, .204494785013794E+02,  &
        .331909213593302E+01, .103400130404874E+00 /), e(7)  &
        = (/ -.999999999998447E+00, -.266271060431811E+02,  &
        -.241055827097015E+03, -.895927957772937E+03, -.129885688746484E+04,  &
        -.545374158883133E+03, -.566575206533869E+01 /), f(7)  &
        = (/ .100000000000000E+01, .286271060422192E+02,  &
        .292310039388533E+03, .133278537748257E+04, .277761949509163E+04,  &
        .240401713225909E+04, .631657483280800E+03 /), p1(8)  &
        = (/ -.866937339951070E+01, -.549142265521085E+0,  &
        -.421001615357070E+04, -.249301393458648E+06, -.119623669349247E+06  &
       , -.221744627758845E+08, .389280421311201E+07, -.195773036904548E+09  &
        /), q1(8) = (/ .341718750000000E+02, -.160708926587221E+04,  &
        .357300298058508E+05, -.483547436162164E+06, .428559624611749E+07,  &
        -.249033375740540E+08, .891925767575612E+08, -.826271498626055E+08  &
        /), p2(8) = (/ -.218086381520723E+01, -.219010233854881E+0,  &
        .930816385662165E+01, .250762811293560E+02, -.331842531997221E+02,  &
        .601217990830080E+02, -.432531132878135E+02, .100443109228078E+01 /)
REAL, PARAMETER :: q2(7) = (/ .393707701852715E+01, .300892648372915E+03,  &
        -.625041161671876E+01, .100367439516726E+04, .143256738121938E+02,  &
        .273624119889328E+04, .527468851962908E+00 /), p3(8)  &
        = (/ -.348334653602852E+01, -.186545454883399E+0,  &
        -.828561994140641E+01, -.323467330305403E+02, .179601688769252E+02,  &
        .175656315469614E+01, -.195022321289660E+01, .999994296074708E+00  &
        /), q3(7) = (/ .695000655887434E+02, .572837193837324E+02,  &
        .257776384238440E+02, .760761148007735E+03, .289516727925135E+02,  &
        -.343942266899870E+01, .100083867402639E+01 /), p4(8)  &
        = (/ -.531686623494482E+02, .891263822573708E+01,  &
        -.139381360364405E+01, -.308336269051763E+0, -.749289167792884E+01,  &
        -.500140345515924E+01, -.300000016782086E+01, .100000000000058E+01  &
        /), q4(7) = (/ .104745362652468E+04, -.674704580465832E+01,  &
        .295999399486831E+03, -.431325836146628E+01, -.790404992298926E+01,  &
        -.299996432944446E+0, .199999999924131E+01 /)
REAL, PARAMETER :: r(20) = (/ .636295897967470E+00, -.130811686750676E+00,  &
        -.843674102130539E-02, .265684915310067E-02, .328227217816581E-03,  &
        -.237834477714302E-04, -.114398043081001E-04, -.144059434332383E-05, &
        .524159566511488E-08, .384073064078443E-07, .858802448602672E-08,  &
        .102192266258550E-08, .217491323232897E-10, -.220902381426231E-10,  &
        -.634575335449288E-11, -.108377465668577E-1, -.119098228722226E-12 ,  &
        -.284386823892656E-14, .250803270266868E-14, .787296415285598E-15 /), &
        x0 = .372507410781366, EXP40 = .235385266837020E+18
INTEGER   :: i, m
REAL      :: px(9), qx(9), ei, frac, sump, sumq, t, w, x, xx0, xmx0,  &
             xmax, xmin, y
REAL (dp), PARAMETER :: dx0 = .37250741078136663446199186658_dp

x = arg
ierr = 0
IF (INT == 2) GO TO 100
IF (x < 0.0) THEN
  GO TO 60
ELSE IF (x == 0.0) THEN
  GO TO 130
END IF
IF (x < 12.) THEN
  IF (x <= 8.) THEN
    IF (x < 4.) THEN
!     ---------- 0.0 < X < 4.0.
!                RATIONAL APPROXIMATION USED IS EXPRESSED
!                IN TERMS OF CHEBYSHEV POLYNOMIALS TO
!                IMPROVE CONDITIONING  ----------
      t = x + x
      t = t / 3.0 - 2.0
      px(1) = 0.0
      qx(1) = 0.0
      px(2) = p1(1)
      qx(2) = q1(1)

      DO i = 2, 7
        px(i+1) = t * px(i) - px(i-1) + p1(i)
        qx(i+1) = t * qx(i) - qx(i-1) + q1(i)
      END DO

      sump = .5 * t * px(8) - px(7) + p1(8)
      sumq = .5 * t * qx(8) - qx(7) + q1(8)
      frac = sump / sumq
      xmx0 = DBLE(x) - dx0
      IF (ABS(xmx0) >= 0.07) THEN
        xx0 = x / x0
        ei = LOG(xx0) + xmx0 * frac
        IF (INT == 3) ei = EXP(-x) * ei
        GO TO 90
      END IF
!     ---------- EVALUATE APPROXIMATION FOR LN(X/X0)
!                FOR X CLOSE TO X0 ----------
      y = xmx0 / x0
      ei = alnrel(y) + xmx0 * frac
      IF (INT == 3) ei = EXP(-x) * ei
      GO TO 90
    END IF
!     ---------- 4.0 <= X <= 8.0 ----------
    m = 20
    ei = (1.0+csevl(3.0-16.0/x,r,m)) / x
    IF (INT == 3) GO TO 90
    ei = ei * EXP(x)
    GO TO 90
  END IF
!     ---------- 8.0 < X < 12.0 ----------
  frac = 0.0

  DO i = 1, 7
    frac = q2(i) / (p2(i)+x+frac)
  END DO

  ei = (p2(8)+frac) / x
  IF (INT == 3) GO TO 90
  ei = ei * EXP(x)
  GO TO 90
END IF
!     ---------- 12.0 <= X < 24.0 ----------
IF (x < 24.) THEN
  frac = 0.0

  DO i = 1, 7
    frac = q3(i) / (p3(i)+x+frac)
  END DO

  ei = (p3(8)+frac) / x
  IF (INT == 3) GO TO 90
  ei = ei * EXP(x)
  GO TO 90
END IF
!     ---------- 24.0 <= X ----------
xmax = exparg(0)
IF ((x > xmax) .AND. (INT < 3)) GO TO 120
y = 1.0 / x
frac = 0.0

DO i = 1, 7
  frac = q4(i) / (p4(i)+x+frac)
END DO

frac = p4(8) + frac
ei = y + y * y * frac
IF (INT == 3) GO TO 90
IF (x <= 150.0) THEN
  ei = ei * EXP(x)
  GO TO 90
END IF
!     ---------- CALCULATION REFORMULATED TO AVOID
!                PREMATURE OVERFLOW ----------
ei = (ei*EXP(x-40.0)) * EXP40
GO TO 90
!     ---------- ORIGINAL X WAS NEGATIVE.  CALCULATION OF
!                E-SUB-1 JOINS AT LABEL 300 ----------
60 y = -x
70 w = 1.0 / y
IF (y <= 4.0) THEN
  IF (y <= 1.0) THEN
!     ---------- 0.0 < -X <= 1.0 ----------
    ei = LOG(y) - (((((a(6)*y+a(5))*y+a(4))*y+a(3))*y+a(2))*y+  &
    a(1)) / ((((b(5)*y+b(4))*y+b(3))*y+b(2))*y+b(1))
    IF (INT == 3) ei = ei * EXP(y)
    GO TO 80
  END IF
!     ---------- 1.0 < -X <= 4.0 ----------
  ei = -((((((c(7)*w+c(6))*w+c(5))*w+c(4))*w+c(3))*w+c(2))*w+  &
  c(1)) / ((((((d(7)*w+d(6))*w+d(5))*w+d(4))*w+d(3))*w+d(2))*w+d(1))
  IF (INT == 3) GO TO 90
  ei = ei * EXP(-y)
ELSE
!     ---------- 4.0 < -X ----------
  xmin = exparg(1)
  IF ((-ABS(x) < xmin) .AND. (INT < 3)) GO TO 110
  ei = -w * (1.0+w*((((((e(7)*w+e(6))*w+e(5))*w+e(4))*w+e(3))*w+  &
  e(2))*w+e(1))/((((((f(7)*w+f(6))*w+f(5))*w+f(4))*w+f(3))*w+f(2))*w+f(1)))
  IF (INT == 3) GO TO 90
  ei = ei * EXP(-y)
  t = 0.5 * ei
  IF (t == 0.0) GO TO 110
END IF
80 IF (INT == 2) ei = -ei
90 result = ei
RETURN

100 y = x
IF (y < 0.0) THEN
  GO TO 140
ELSE IF (y > 0.0) THEN
  GO TO 70
ELSE
  GO TO 130
END IF
!     ---------- ERROR RETURN FOR X < XMIN,
!                CAUSING UNDERFLOW ----------
110 ei = 0.0
ierr = 1
GO TO 90
!     ---------- ERROR RETURN FOR X > XMAX,
!                CAUSING OVERFLOW ----------
120 ierr = 2
RETURN
!     ---------- ERROR RETURN FOR ILLEGAL
!                ARGUMENT, X = 0 ----------
130 ierr = 3
RETURN
!     ---------- ERROR RETURN FOR NEGATIVE
!                ARGUMENT IN E-SUB-1 ----------
140 ierr = 4
RETURN
END SUBROUTINE expli


FUNCTION dei(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!        REAL (dp) EVALUATION OF THE EXPONENTIAL INTEGRAL
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

fn_val = -de1e(-x)
IF (x > 4._dp.OR.x < -1._dp) fn_val = EXP(x) * fn_val
RETURN
END FUNCTION dei


FUNCTION dei1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!          REAL (dp) EVALUATION OF EXP(-X)*EI(X)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

fn_val = -de1e(-x)
IF (x > 4._dp.OR.x < -1._dp) RETURN
fn_val = EXP(-x) * fn_val
RETURN
END FUNCTION dei1


FUNCTION de1e(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!     LET E1(X) DENOTE THE EXPONENTIAL INTEGRAL FOR POSITIVE X AND
!     THE CAUCHY PRINCIPAL VALUE FOR NEGATIVE X. IF X IS NONZERO
!     THEN DE1E HAS THE VALUE ...

!          DE1E(X) = E1(X)          IF -4 <= X <= 1
!          DE1E(X) = EXP(X)*E1(X)   OTHERWISE

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------

!     THE FOLLOWING SERIES FOR E1 WERE DEVELOPED BY WAYNE FULLERTON
!     (LOS ALAMOS NATIONAL LABORATORY).

!     SERIES A             ON THE INTERVAL -3.12500E-02 TO  0.
!                                        WITH WEIGHTED ERROR   4.62E-32
!                                         LOG WEIGHTED ERROR  31.34
!                               SIGNIFICANT FIGURES REQUIRED  29.70
!                                    DECIMAL PLACES REQUIRED  32.18


!     SERIES B             ON THE INTERVAL -1.25000E-01 TO -3.12500E-02
!                                        WITH WEIGHTED ERROR   2.22E-32
!                                         LOG WEIGHTED ERROR  31.65
!                               SIGNIFICANT FIGURES REQUIRED  30.75
!                                    DECIMAL PLACES REQUIRED  32.54


!     SERIES D             ON THE INTERVAL -2.50000E-01 TO -1.25000E-01
!                                        WITH WEIGHTED ERROR   5.19E-32
!                                         LOG WEIGHTED ERROR  31.28
!                               SIGNIFICANT FIGURES REQUIRED  30.82
!                                    DECIMAL PLACES REQUIRED  32.09


!     SERIES E             ON THE INTERVAL -4.00000E+00 TO -1.00000E+00
!                                        WITH WEIGHTED ERROR   8.49E-34
!                                         LOG WEIGHTED ERROR  33.07
!                               SIGNIFICANT FIGURES REQUIRED  34.13
!                                    DECIMAL PLACES REQUIRED  33.80


!     SERIES R             ON THE INTERVAL -1.00000E+00 TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   8.08E-33
!                                         LOG WEIGHTED ERROR  32.09
!                        APPROX SIGNIFICANT FIGURES REQUIRED  30.4
!                                    DECIMAL PLACES REQUIRED  32.79


!     SERIES P             ON THE INTERVAL  2.50000E-01 TO  1.00000E+00
!                                        WITH WEIGHTED ERROR   6.65E-32
!                                         LOG WEIGHTED ERROR  31.18
!                               SIGNIFICANT FIGURES REQUIRED  30.69
!                                    DECIMAL PLACES REQUIRED  32.03


!     SERIES Q             ON THE INTERVAL  0.          TO  2.50000E-01
!                                        WITH WEIGHTED ERROR   5.07E-32
!                                         LOG WEIGHTED ERROR  31.30
!                               SIGNIFICANT FIGURES REQUIRED  30.40
!                                    DECIMAL PLACES REQUIRED  32.20

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: a(50) = (/ .3284394579616699087873844201881D-01,  &
    -.1669920452031362851476184343387D-01, .2845284724361346807424899853252D-03,  &
    -.7563944358516206489487866938533D-05, .2798971289450859157504843180879D-06,  &
    -.1357901828534531069525563926255D-07, .8343596202040469255856102904906D-09,  &
    -.6370971727640248438275242988532D-10, .6007247608811861235760831561584D-11,  &
    -.7022876174679773590750626150088D-12, .1018302673703687693096652346883D-12,  &
    -.1761812903430880040406309966422D-13, .3250828614235360694244030353877D-14,  &
    -.5071770025505818678824872259044D-15, .1665177387043294298172486084156D-16,  &
    .3166753890797514400677003536555D-16, -.1588403763664141515133118343538D-16,  &
    .4175513256138018833003034618484D-17, -.2892347749707141906710714478852D-18,  &
    -.2800625903396608103506340589669D-18, .1322938639539270903707580023781D-18,  &
    -.1804447444177301627283887833557D-19, -.7905384086522616076291644817604D-20,  &
    .4435711366369570103946235838027D-20, -.4264103994978120868865309206555D-21,  &
    -.3920101766937117541553713162048D-21, .1527378051343994266343752326971D-21,  &
    .1024849527049372339310308783117D-22, -.2134907874771433576262711405882D-22,  &
    .3239139475160028267061694700366D-23,  &
    .2142183762299889954762643168296D-23, -.8234609419601018414700348082312D-24,  &
    -.1524652829645809479613694401140D-24, .1378208282460639134668480364325D-24,  &
    .2131311202833947879523224999253D-26, -.2012649651526484121817466763127D-25,  &
    .1995535662263358016106311782673D-26, .2798995808984003464948686520319D-26,  &
    -.5534511845389626637640819277823D-27,  &
    -.3884995396159968861682544026146D-27, .1121304434507359382850680354679D-27,  &
    .5566568152423740948256563833514D-28, -.2045482929810499700448533938176D-28,  &
    -.8453813992712336233411457493674D-29, .3565758433431291562816111116287D-29,  &
    .1383653872125634705539949098871D-29, -.6062167864451372436584533764778D-30,  &
    -.2447198043989313267437655119189D-30, .1006850640933998348011548180480D-30,  &
    .4623685555014869015664341461674D-31 /)
REAL (dp), PARAMETER :: b(60) = (/ .20263150647078889499401236517381D+00,  &
    -.73655140991203130439536898728034D-01, .63909349118361915862753283840020D-02,  &
    -.60797252705247911780653153363999D-03, -.73706498620176629330681411493484D-04,  &
    .48732857449450183453464992488076D-04, -.23837064840448290766588489460235D-05,  &
    -.30518612628561521027027332246121D-05, .17050331572564559009688032992907D-06,  &
    .23834204527487747258601598136403D-06, .10781772556163166562596872364020D-07,  &
    -.17955692847399102653642691446599D-07, -.41284072341950457727912394640436D-08,  &
    .68622148588631968618346844526664D-09, .53130183120506356147602009675961D-09,  &
    .78796880261490694831305022893515D-10, -.26261762329356522290341675271232D-10,  &
    -.15483687636308261963125756294100D-10, -.25818962377261390492802405122591D-11,  &
    .59542879191591072658903529959352D-12, .46451400387681525833784919321405D-12,  &
    .11557855023255861496288006203731D-12, -.10475236870835799012317547189670D-14,  &
    -.11896653502709004368104489260929D-13, -.47749077490261778752643019349950D-14,  &
    -.81077649615772777976249734754135D-15, .13435569250031554199376987998178D-15,  &
    .14134530022913106260248873881287D-15, .49451592573953173115520663232883D-16,  &
    .79884048480080665648858587399367D-17, -.14008632188089809829248711935393D-17,  &
    -.14814246958417372107722804001680D-17, -.55826173646025601904010693937113D-18,  &
    -.11442074542191647264783072544598D-18, .25371823879566853500524018479923D-20,  &
    .13205328154805359813278863389097D-19, .62930261081586809166287426789485D-20,  &
    .17688270424882713734999261332548D-20, .23266187985146045209674296887432D-21,  &
    -.67803060811125233043773831844113D-22, -.59440876959676373802874150531891D-22,  &
    -.23618214531184415968532592503466D-22, -.60214499724601478214168478744576D-23,  &
    -.65517906474348299071370444144639D-24, .29388755297497724587042038699349D-24,  &
    .22601606200642115173215728758510D-24, .89534369245958628745091206873087D-25,  &
    .24015923471098457555772067457706D-25, .34118376888907172955666423043413D-26,  &
    -.71617071694630342052355013345279D-27, -.75620390659281725157928651980799D-27,  &
    -.33774612157467324637952920780800D-27, -.10479325703300941711526430332245D-27,  &
    -.21654550252170342240854880201386D-28, -.75297125745288269994689298432000D-30,  &
    .19103179392798935768638084000426D-29, .11492104966530338547790728833706D-29,  &
    .43896970582661751514410359193600D-30, .12320883239205686471647157725866D-30,  &
    .22220174457553175317538581162666D-31 /)
REAL (dp), PARAMETER :: d(41) = (/ .63629589796747038767129887806803D+00,  &
    -.13081168675067634385812671121135D+00, -.84367410213053930014487662129752D-02,  &
    .26568491531006685413029428068906D-02, .32822721781658133778792170142517D-03,  &
    -.23783447771430248269579807851050D-04, -.11439804308100055514447076797047D-04,  &
    -.14405943433238338455239717699323D-05, .52415956651148829963772818061664D-08,  &
    .38407306407844323480979203059716D-07, .85880244860267195879660515759344D-08,  &
    .10219226625855003286339969553911D-08, .21749132323289724542821339805992D-10,  &
    -.22090238142623144809523503811741D-10, -.63457533544928753294383622208801D-11,  &
    -.10837746566857661115340539732919D-11, -.11909822872222586730262200440277D-12,  &
    -.28438682389265590299508766008661D-14, .25080327026686769668587195487546D-14,  &
    .78729641528559842431597726421265D-15, .15475066347785217148484334637329D-15,  &
    .22575322831665075055272608197290D-16, .22233352867266608760281380836693D-17,  &
    .16967819563544153513464194662399D-19, -.57608316255947682105310087304533D-19,  &
    -.17591235774646878055625369408853D-19, -.36286056375103174394755328682666D-20,  &
    -.59235569797328991652558143488000D-21, -.76030380926310191114429136895999D-22,  &
    -.62547843521711763842641428479999D-23, .25483360759307648606037606400000D-24,  &
    .25598615731739857020168874666666D-24, .71376239357899318800207052800000D-25,  &
    .14703759939567568181578956800000D-25, .25105524765386733555198634666666D-26,  &
    .35886666387790890886583637333333D-27, .39886035156771301763317759999999D-28,  &
    .21763676947356220478805333333333D-29, -.46146998487618942367607466666666D-30,  &
    -.20713517877481987707153066666666D-30, -.51890378563534371596970666666666D-31 /)
REAL (dp), PARAMETER :: e(29) = (/ -.16113461655571494025720663927566180D+02,  &
    .77940727787426802769272245891741497D+01,  &
    -.19554058188631419507127283812814491D+01,  &
    .37337293866277945611517190865690209D+00,  &
    -.56925031910929019385263892220051166D-01,  &
    .72110777696600918537847724812635813D-02,  &
    -.78104901449841593997715184089064148D-03,  &
    .73880933562621681878974881366177858D-04,  &
    -.62028618758082045134358133607909712D-05,  &
    .46816002303176735524405823868362657D-06,  &
    -.32092888533298649524072553027228719D-07,  &
    .20151997487404533394826262213019548D-08,  &
    -.11673686816697793105356271695015419D-09,  &
    .62762706672039943397788748379615573D-11,  &
    -.31481541672275441045246781802393600D-12,  &
    .14799041744493474210894472251733333D-13,  &
    -.65457091583979673774263401588053333D-15,  &
    .27336872223137291142508012748799999D-16,  &
    -.10813524349754406876721727624533333D-17,  &
    .40628328040434303295300348586666666D-19,  &
    -.14535539358960455858914372266666666D-20,  &
    .49632746181648636830198442666666666D-22,  &
    -.16208612696636044604866560000000000D-23,  &
    .50721448038607422226431999999999999D-25,  &
    -.15235811133372207813973333333333333D-26,  &
    .44001511256103618696533333333333333D-28,  &
    -.12236141945416231594666666666666666D-29,  &
    .32809216661066001066666666666666666D-31,  &
    -.84933452268306432000000000000000000D-33 /)
REAL (dp), PARAMETER :: r(25) = (/ -.3739021479220279511668698204827D-01,  &
    .4272398606220957726049179176528D-01,  &
    -.130318207984970054415392055219726D+00,  &
    .144191240246988907341095893982137D-01,  &
    -.134617078051068022116121527983553D-02,  &
    .107310292530637799976115850970073D-03,  &
    -.742999951611943649610283062223163D-05,  &
    .453773256907537139386383211511827D-06,  &
    -.247641721139060131846547423802912D-07,  &
    .122076581374590953700228167846102D-08,  &
    -.548514148064092393821357398028261D-10,  &
    .226362142130078799293688162377002D-11,  &
    -.863589727169800979404172916282240D-13,  &
    .306291553669332997581032894881279D-14,  &
    -.101485718855944147557128906734933D-15,  &
    .315482174034069877546855328426666D-17,  &
    -.923604240769240954484015923200000D-19,  &
    .255504267970814002440435029333333D-20,  &
    -.669912805684566847217882453333333D-22,  &
    .166925405435387319431987199999999D-23,  &
    -.396254925184379641856000000000000D-25,  &
    .898135896598511332010666666666666D-27,  &
    -.194763366993016433322666666666666D-28,  &
    .404836019024630033066666666666666D-30,  &
    -.807981567699845120000000000000000D-32 /)
REAL (dp), PARAMETER :: p(50) = (/ -.60577324664060345999319382737747D+00,  &
    -.11253524348366090030649768852718D+00, .13432266247902779492487859329414D-01,  &
    -.19268451873811457249246838991303D-02, .30911833772060318335586737475368D-03,  &
    -.53564132129618418776393559795147D-04, .98278128802474923952491882717237D-05,  &
    -.18853689849165182826902891938910D-05, .37494319356894735406964042190531D-06,  &
    -.76823455870552639273733465680556D-07, .16143270567198777552956300060868D-07,  &
    -.34668022114907354566309060226027D-08, .75875420919036277572889747054114D-09,  &
    -.16886433329881412573514526636703D-09, .38145706749552265682804250927272D-10,  &
    -.87330266324446292706851718272334D-11, .20236728645867960961794311064330D-11,  &
    -.47413283039555834655210340820160D-12, .11221172048389864324731799928920D-12,  &
    -.26804225434840309912826809093395D-13, .64578514417716530343580369067212D-14,  &
    -.15682760501666478830305702849194D-14, .38367865399315404861821516441408D-15,  &
    -.94517173027579130478871048932556D-16, .23434812288949573293896666439133D-16,  &
    -.58458661580214714576123194419882D-17, .14666229867947778605873617419195D-17,  &
    -.36993923476444472706592538274474D-18, .93790159936721242136014291817813D-19,  &
    -.23893673221937873136308224087381D-19, .61150624629497608051934223837866D-20,  &
    -.15718585327554025507719853288106D-20, .40572387285585397769519294491306D-21,  &
    -.10514026554738034990566367122773D-21, .27349664930638667785806003131733D-22,  &
    -.71401604080205796099355574271999D-23, .18705552432235079986756924211199D-23,  &
    -.49167468166870480520478020949333D-24, .12964988119684031730916087125333D-24,  &
    -.34292515688362864461623940437333D-25, .90972241643887034329104820906666D-26,  &
    -.24202112314316856489934847999999D-26, .64563612934639510757670475093333D-27,  &
    -.17269132735340541122315987626666D-27, .46308611659151500715194231466666D-28,  &
    -.12448703637214131241755170133333D-28, .33544574090520678532907007999999D-29,  &
    -.90598868521070774437543935999999D-30, .24524147051474238587273216000000D-30,  &
    -.66528178733552062817107967999999D-31 /)
REAL (dp), PARAMETER :: q(64) = (/ -.1892918000753016825495679942820D+00,  &
    -.8648117855259871489968817056824D-01, .7224101543746594747021514839184D-02,  &
    -.8097559457557386197159655610181D-03, .1099913443266138867179251157002D-03,  &
    -.1717332998937767371495358814487D-04, .2985627514479283322825342495003D-05,  &
    -.5659649145771930056560167267155D-06, .1152680839714140019226583501663D-06,  &
    -.2495030440269338228842128765065D-07, .5692324201833754367039370368140D-08,  &
    -.1359957664805600338490030939176D-08, .3384662888760884590184512925859D-09,  &
    -.8737853904474681952350849316580D-10, .2331588663222659718612613400470D-10,  &
    -.6411481049213785969753165196326D-11, .1812246980204816433384359484682D-11,  &
    -.5253831761558460688819403840466D-12, .1559218272591925698855028609825D-12,  &
    -.4729168297080398718476429369466D-13, .1463761864393243502076199493808D-13,  &
    -.4617388988712924102232173623604D-14, .1482710348289369323789239660371D-14,  &
    -.4841672496239229146973165734417D-15, .1606215575700290408116571966188D-15,  &
    -.5408917538957170947895023784252D-16, .1847470159346897881370231402310D-16,  &
    -.6395830792759094470500610425050D-17, .2242780721699759457250233276170D-17,  &
    -.7961369173983947552744555308646D-18, .2859308111540197459808619929272D-18,  &
    -.1038450244701137145900697137446D-18, .3812040607097975780866841008319D-19,  &
    -.1413795417717200768717562723696D-19, .5295367865182740958305442594815D-20,  &
    -.2002264245026825902137211131439D-20, .7640262751275196014736848610918D-21,  &
    -.2941119006868787883311263523362D-21, .1141823539078927193037691483586D-21,  &
    -.4469308475955298425247020718489D-22, .1763262410571750770630491408520D-22,  &
    -.7009968187925902356351518262340D-23, .2807573556558378922287757507515D-23,  &
    -.1132560944981086432141888891562D-23, .4600574684375017946156764233727D-24,  &
    -.1881448598976133459864609148108D-24, .7744916111507730845444328478037D-25,  &
    -.3208512760585368926702703826261D-25, .1337445542910839760619930421384D-25,  &
    -.5608671881802217048894771735210D-26, .2365839716528537483710069473279D-26,  &
    -.1003656195025305334065834526856D-26, .4281490878094161131286642556927D-27,  &
    -.1836345261815318199691326958250D-27, .7917798231349540000097468678144D-28,  &
    -.3431542358742220361025015775231D-28, .1494705493897103237475066008917D-28,  &
    -.6542620279865705439739042420053D-29, .2877581395199171114340487353685D-29,  &
    -.1271557211796024711027981200042D-29, .5644615555648722522388044622506D-30,  &
    -.2516994994284095106080616830293D-30, .1127259818927510206370368804181D-30,  &
    -.5069814875800460855562584719360D-31 /)
REAL (dp) :: c, eps, t, w
INTEGER   :: m
!------------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1._dp + EPS > 1._dp.

eps = EPSILON(1.0_dp)

!------------------------------
IF (ABS(x) < 90._dp) THEN
  IF (x <= -1._dp) THEN

!                    -90 < X < -4

    IF (x <= -32._dp) THEN
      m = 50
      IF (eps >= 1.d-20) m = 25
      fn_val = (1._dp+dcsevl(64._dp/x+1._dp,a,m)) / x
      RETURN
    END IF

    IF (x <= -8._dp) THEN
      m = 60
      IF (eps >= 1.d-20) m = 37
      fn_val = (1._dp+dcsevl((64._dp/x+5._dp)/3._dp,b,m)) / x
      RETURN
    END IF

    IF (x < -4._dp) THEN
      m = 41
      IF (eps >= 1.d-20) m = 27
      fn_val = (1._dp+dcsevl(16._dp/x+3._dp,d,m)) / x
      RETURN
    END IF

!                     -4 <= X <= 1

    m = 29
    IF (eps >= 1.d-20) m = 20
    fn_val = -LOG(-x) + dcsevl((2._dp*x+5._dp)/3._dp,e,m)
    RETURN
  END IF

  IF (x <= 1.0_dp) THEN
    IF (x >= -0.4_dp .AND. x <= -0.35_dp) THEN
      fn_val = -dei0(-x,eps)
      RETURN
    END IF
    m = 25
    IF (eps >= 1.d-20) m = 18
    fn_val = (-LOG(ABS(x))-0.6875_dp+x) + dcsevl(x,r,m)
    RETURN
  END IF

!                     1 < X < 90

  IF (x <= 4.0_dp) THEN
    m = 50
    IF (eps >= 1.d-20) m = 31
    fn_val = (1._dp+dcsevl((8._dp/x-5._dp)/3._dp,p,m)) / x
    RETURN
  END IF

  m = 64
  IF (eps >= 1.d-20) m = 35
  fn_val = (1._dp+dcsevl(8._dp/x-1._dp,q,m)) / x
  RETURN
END IF

!                   ASYMPTOTIC EXPANSION

t = -1._dp / x
c = t
w = c
m = 1
10 m = m + 1
c = (m*t) * c
w = c + w
IF (ABS(c) > eps) GO TO 10
fn_val = (1._dp+w) / x
RETURN
END FUNCTION de1e


FUNCTION dei0(x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------

!            TAYLOR SERIES EXPANSION OF EI(X) AROUND X0,
!                  WHERE X0 IS THE ZERO OF EI(X).
!                    EPS IS THE TOLERANCE USED.

!-------------------------
!     WRITTEN BY A.H. MORRIS
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x, eps
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: a(40) = (/ .3896215733907167310156502703593482682018D+01,  &
        -.3281607866398561670879044070702055058438D+01,  &
        .6522376145438925697728352767902339522245D+01,  &
        -.1296969738353651703636356975116693457132D+02,  &
        .2788629796294204997855360701398702087604D+02,  &
        -.6237880152891541873078526672920295283143D+02,  &
        .1435349488096750987841265647073135861344D+03,  &
        -.3371558271787468916821364466977375583658D+03,  &
        .8045318399821382506595322457265602778098D+03,  &
        -.1943796645723498840655451915157946462648D+04,  &
        .4743765650402430835228269085129777320454D+04,  &
        -.1167346399116716364394668734600584330571D+05,  &
        .2892695530543545087445160311373446386859D+05,  &
        -.7210794586837158996878001987822898188198D+05,  &
        .1806695585893919626172098163733836311447D+06,  &
        -.4546962188544665746524572520110515778526D+06,  &
        .1148834546817744310374556891236080193473D+07,  &
        -.2912721663850837498392234670693435881386D+07,  &
        .7407692958000587759747953639495510375408D+07,  &
        -.1889172700038153127288849726417780854730D+08,  &
        .4830003493086024720868271496253148288055D+08,  &
        -.1237682190024917092137405370407916520821D+09,  &
        .3178111056663621852260468265468336116367D+09,  &
        -.8176185693184928170769596413793786736279D+09,  &
        .2107109291864363741291032089276432438927D+10,  &
        -.5438996831077284596300440196418865401363D+10,  &
        .1406026390995585037838894210474681627693D+11,  &
        -.3639689100149205333626392754168250384373D+11,  &
        .9433859509219164512733865811047107199412D+11,  &
        -.2448111705066430130314746602027041462835D+12,  &
        .6359981818273706257655285041587660739835D+12,  &
        -.1653989211524391716301960841541179924503D+13,  &
        .4305601123377464671923711939926758523701D+13,  &
        -.1121847693567642152208443795868288937687D+14,  &
        .2925565695557339262045727352754930716608D+14,  &
        -.7635552741959392076619218035480359307499D+14,  &
        .1994372792759425025753893705017248674884D+15,  &
        -.5213021921201092276891722450906568692592D+15,  &
        .1363558024737805584657706536660107687818D+16,  &
        -.3568973490569445692988895507297245137908D+16 /)
REAL (dp) :: c, db2, h, t, w
INTEGER   :: n
REAL (dp) :: dk1 = 25598514349._dp, dk2 = 12212826724._dp,  &
             dk3 = 52346020729._dp, db = 68719476736._dp,  &
              dx = .64725688445954142292644880487403537155379408215561D-33
!-------------------------

!     SET  H = X - X0  WHERE X0 IS THE ZERO OF EI(X). X0 HAS THE
!     APPROXIMATE 60 DIGIT VALUE ...

!      .37250741078136663446 19918665801191335356 89497771654051555657

!     A MORE ACCURATE VALUE IS GIVEN BY ...

!            X0 = DK1/8**12 + DK2/8**24 + DK3/8**36 + DX

!     THE FOLLOWING CODE SHOULD YIELD THE CORRECT VALUE FOR H IF A
!     BINARY, OCTAL, OR HEXADECIMAL REAL (dp) ARITHMETIC IS BEING USED.

db2 = db * db
h = (((x-dk1/db)-dk2/db2)-dk3/(db*db2)) - dx

!-------------------------
t = h
w = 0._dp
DO n = 2, 40
  c = a(n) * t
  w = w + c
  IF (ABS(c) < eps) GO TO 20
  t = h * t
END DO

20 fn_val = h * (a(1)+w)
RETURN
END FUNCTION dei0


FUNCTION si(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!              EVALUATION OF THE SINE INTEGRAL FUNCTION

!                         ------------------

!     ALGORITHM.  A MINIMAX APPROXIMATION OBTAINED BY A. H. MORRIS
!     IS USED WHEN ABS(X) <= 5, AND THE CHEBYSHEV EXPANSION GIVEN
!     ON PAGE 326 OF THE REFERENCE IS USED WHEN ABS(X) > 5.

!     REFERENCE.  LUKE, YUDELL L., THE SPECIAL FUNCTIONS AND THEIR
!     APPROXIMATIONS, VOL. 2, ACADEMIC PRESS, NEW YORK, 1969.

!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!        FEB 1993

!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: c(46) = (/ 9.76155271128712E-01, 8.96845854916423E-02,  &
        -3.04656658030696E-02, 8.50892472922945E-02, -5.78073683148386E-03  &
       , -5.07182677775691E-03, 8.38643256650893E-04, -3.34223415981738E-04  &
       , -2.15746207281216E-05, 1.28560650086065E-04, -1.56456413510232E-05  &
       , -1.52025513597262E-0, 4.04001013843204E-06, -5.95896122752160E-0,  &
        -4.34985305974340E-07, 7.13472533530840E-07, -5.34302186061100E-08  &
       , -1.76003581156610E-07, 3.85028855125900E-08, 1.92576544441700E-08  &
       , -1.00735358217200E-08, 3.36359194377000E-09, 1.28049619406000E-09  &
       , -2.42546870827000E-09, 1.86917288950000E-10, 7.13431298340000E-10  &
       , -1.70673483710000E-10, -1.14604070350000E-1, 5.88004411500000E-11  &
       , -6.78417843000000E-1, -1.21572380900000E-11, 1.26561248700000E-11  &
       , 4.74814180000000E-13, -5.32309477000000E-1, 9.05903810000000E-13,  &
        1.40046450000000E-12, -5.00968320000000E-13, -1.80458040000000E-1,  &
        1.66162910000000E-13, -5.02616400000000E-14, -3.48453600000000E-14  &
       , 4.60056600000000E-14, 5.74000000000000E-16, -1.95310700000000E-14  &
       , 3.68837000000000E-15, 5.62862000000000E-15 /),  &
        pihalf = 1.5707963267949, a1 = -.480279472444504E-01,  &
        a2 = .127177528378855E-02, a3 = -.170630463362755E-04,  &
        a4 = .129975549721579E-06, a5 = -.582322888431340E-09,  &
        a6 = .148011790132481E-11, a7 = -.172103429855786E-14,  &
        b1 = .752760831110726E-02, b2 = .233090788469112E-04,  &
        b3 = .305598403979701E-07
INTEGER, PARAMETER :: n = 46, m = 21
REAL    :: ax, eps, p, q, t, t1, t2, w, z
INTEGER :: i, j
!-------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0.
!            IT IS ASSUMED THAT SIN(X) AND COS(X) ARE DEFINED
!            FOR ABS(X) < 1.0/EPS.

!-------------------------
ax = ABS(x)
IF (ax <= 5.0) THEN

!                          ABS(X) <= 5

  t = x * x
  w = ((((((((a7*t + a6)*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t+0.5)+0.5) /  &
      (((b3*t + b2)*t + b1)*t+1.0)
  fn_val = x * w
  RETURN
END IF

!                          ABS(X) > 5

IF (ax >= 1.e+5) THEN
  eps = EPSILON(1.0)
  IF (ax*eps >= 1.0) GO TO 30
END IF

z = 10.0 / ax - 1.0
w = z + z
j = n - 1
t1 = c(j)
t2 = 0.0
DO i = 1, m
  j = j - 2
  t = t1
  t1 = w * t1 - t2 + c(j)
  t2 = t
END DO
p = z * t1 - t2 + c(1)

j = n
t1 = c(j)
t2 = 0.0
DO i = 1, m
  j = j - 2
  t = t1
  t1 = w * t1 - t2 + c(j)
  t2 = t
END DO
q = z * t1 - t2 + c(2)

fn_val = (p*COS(ax)+q*SIN(ax)) / ax
fn_val = pihalf - fn_val
IF (x < 0.0) fn_val = -fn_val
RETURN

!                        ABS(X) >= 1/EPS

30 fn_val = SIGN(pihalf,x)
RETURN
END FUNCTION si


FUNCTION cin(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!             EVALUATION OF THE COSINE INTEGRAL FUNCTION

!                        --------------------

!     ALGORITHM.  A MINIMAX APPROXIMATION OBTAINED BY A. H. MORRIS
!     IS USED WHEN ABS(X) <= 5, AND THE CHEBYSHEV EXPANSION GIVEN
!     ON PAGE 326 OF THE REFERENCE IS USED WHEN ABS(X) > 5.

!     REFERENCE.  LUKE, YUDELL L., THE SPECIAL FUNCTIONS AND THEIR
!     APPROXIMATIONS, VOL. 2, ACADEMIC PRESS, NEW YORK, 1969.

!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!        FEB 1993

!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: c(46) = (/ 9.76155271128712E-01, 8.96845854916423E-02,  &
        -3.04656658030696E-02, 8.50892472922945E-02, -5.78073683148386E-03  &
       , -5.07182677775691E-03, 8.38643256650893E-04, -3.34223415981738E-04  &
       , -2.15746207281216E-05, 1.28560650086065E-04, -1.56456413510232E-05  &
       , -1.52025513597262E-05, 4.04001013843204E-06, -5.95896122752160E-07,  &
        -4.34985305974340E-07, 7.13472533530840E-07, -5.34302186061100E-08  &
       , -1.76003581156610E-07, 3.85028855125900E-08, 1.92576544441700E-08  &
       , -1.00735358217200E-08, 3.36359194377000E-09, 1.28049619406000E-09  &
       , -2.42546870827000E-09, 1.86917288950000E-10, 7.13431298340000E-10  &
       , -1.70673483710000E-10, -1.14604070350000E-1, 5.88004411500000E-11  &
       , -6.78417843000000E-12, -1.21572380900000E-11, 1.26561248700000E-11  &
       , 4.74814180000000E-13, -5.32309477000000E-12, 9.05903810000000E-13,  &
        1.40046450000000E-12, -5.00968320000000E-13, -1.80458040000000E-13,  &
        1.66162910000000E-13, -5.02616400000000E-14, -3.48453600000000E-14  &
       , 4.60056600000000E-14, 5.74000000000000E-16, -1.95310700000000E-14  &
       , 3.68837000000000E-15, 5.62862000000000E-15 /),  &
        euler = .577215664901533, a1 = -.3204778460606743E-01,  &
        a2 = .5675344472695174E-03, a3 = -.5156904589230793E-05,  &
        a4 = .2599346476124840E-07, a5 = -.7022103360527070E-10,  &
        a6 = .8201072014709134E-13, b1 = .9618882060598870E-02,  &
        b2 = .4239527386910357E-04, b3 = .1039865274359493E-06,  &
        b4 = .1238483219145955E-09
INTEGER, PARAMETER :: n = 46, m = 21
REAL    :: ax, eps, p, q, t, t1, t2, w, z
INTEGER :: i, j
!-------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0.
!            IT IS ASSUMED THAT SIN(X) AND COS(X) ARE DEFINED
!            FOR ABS(X) < 1.0/EPS.

!-------------------------
ax = ABS(x)
IF (ax <= 5.0) THEN

!                          ABS(X) <= 5

  t = x * x
  w = (((((((a6*t + a5)*t + a4)*t + a3)*t + a2)*t + a1)*t+0.5)+0.5) /  &
      ((((b4 *t + b3)*t + b2)*t + b1)*t+1.0)
  fn_val = 0.25 * t * w
  RETURN
END IF

!                          ABS(X) > 5

IF (ax >= 1.e+5) THEN
  eps = EPSILON(1.0)
  IF (ax*eps >= 1.0) GO TO 30
END IF

z = 10.0 / ax - 1.0
w = z + z
j = n - 1
t1 = c(j)
t2 = 0.0
DO i = 1, m
  j = j - 2
  t = t1
  t1 = w * t1 - t2 + c(j)
  t2 = t
END DO
p = z * t1 - t2 + c(1)

j = n
t1 = c(j)
t2 = 0.0
DO i = 1, m
  j = j - 2
  t = t1
  t1 = w * t1 - t2 + c(j)
  t2 = t
END DO
q = z * t1 - t2 + c(2)

fn_val = (p*SIN(ax)-q*COS(ax)) / ax
fn_val = (euler-fn_val) + LOG(ax)
RETURN

!                        ABS(X) >= 1/EPS

30 fn_val = euler + LOG(ax)
RETURN
END FUNCTION cin


SUBROUTINE cexexi(z, w)
!-----------------------------------------------------------------------
!        COMPUTATION OF THE EXPONENTIAL EXPONENTIAL INTEGRAL
!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: const = 2.46740110027234, euler = .577215664901533,  &
                   zeta2 = 1.64493406684823
REAL, PARAMETER :: ed(18) = (/ 0.00000000000000E+00, .311105957086528E-01,  &
        .103661260539112E+00, .216532335244554E+00, .369931427960192E+00,  &
        .566766259990589E+00, .814042066324748E+00, .112384247540813E+01,  &
        .151400478148512E+01, .200886795032284E+01, .264052411823592E+01,  &
        .345098449933392E+01, .449583360763202E+01, .585058263409822E+01,  &
        .762273501463380E+01, .997814501584578E+01, .132122064896408E+02,  &
        .180322948376021E+02 /)
REAL, PARAMETER :: ee(18) = (/ .850156516121093E-02,  &
        .505037465849058E-01, .836817368956407E-01, .107047582417607E+00,  &
        .120424719029462E+00, .125096631582229E+00, .122314435224685E+00,  &
        .112621417553907E+00, .963419407392582E-01, .747398422757511E-01,  &
        .508596135953441E-01, .290822706773628E-01, .132201640530101E-01,  &
        .443802939829067E-02, .992612478987576E-03, .126579795112011E-03,  &
        .702150908253350E-05, .910281532564632E-07 /)
REAL, PARAMETER :: dd(19) = (/ 0.00000000000000E+00, .419556678374293E-01,  &
        .117533661648665E+00, .228560237455987E+00, .375667350161240E+00,  &
        .791594846276672E+00, .107546889623058E+01, .142659208030841E+01,  &
        .186290554952377E+01, .240730009509856E+01, .308854035607524E+01,  &
        .394277605155259E+01, .501593196543981E+01, .636759180748651E+01,  &
        .807776193096055E+01, .102598961138887E+02, .130896768422610E+02,  &
        .168832169085916E+02, .224083240941713E+02 /)
REAL, PARAMETER :: de(19) = (/ -.346911733535892E-03, -.603787732461745E-02,  &
        -.152461305949249E-01, -.210582169827291E-01, -.171894208720754E-01  &
       , .314323467033032E-01, .750898531566972E-01, .124689787807260E+00,  &
        .168579075090035E+00, .191715080699511E+00, .182600794550836E+00,  &
        .142345674307147E+00, .874862222419327E-01, .402175083288425E-01,  &
        .128575005680180E-01, .257673782598441E-02, .275955003784349E-03,  &
        .119139315517122E-04, .107292980199386E-06 /)
REAL, PARAMETER :: cd(18) = (/ 0.00000000000000E+00, .237286128313683E-01,  &
        .854113210668760E-01, .185276627282059E+00, .323741526616688E+00,  &
        .503045460381267E+00, .728806607587188E+00, .101122770102872E+01,  &
        .136598448171249E+01, .181510139038929E+01, .238824701955419E+01,  &
        .312490532008812E+01, .407802489894445E+01, .532033545554865E+01,  &
        .695624307290579E+01, .914759902547031E+01, .121829074388544E+02,  &
        .167511311969873E+02 /)
REAL, PARAMETER :: ce(18) = (/ .349517258926827E-01,  &
        .135849105925897E+00, .158850581552296E+00, .153001434535435E+00,  &
        .134520752856461E+00, .111913051619671E+00, .892008386656190E-01,  &
        .679227205472067E-01, .486723197887211E-01, .320170976532266E-01,  &
        .187008965021111E-01, .929708414427865E-02, .372604763161087E-02,  &
        .111989537559823E-02, .228057496872353E-03, .269596227781453E-04,  &
        .141255430224301E-05, .176352326808806E-07 /)
REAL    :: ts(2), sr(2), sm(2), tm(2), qf(2), zl(2), pm, qm, r, ss, x, y
INTEGER :: i
!-------------------------
!     EULER = EULER CONSTANT
!     CONST = (PI*PI)/4
!     ZETA2 = THE RIEMANN ZETA FUNCTION EVALUATED AT 2
!-------------------------
x = REAL(z)
y = AIMAG(z)
r = x * x + y * y

IF (r > 1.0) THEN
  IF (r < 1296.0) THEN
    IF (x < 0.07*y*y) GO TO 20
    GO TO 40
  END IF
!-----------------------------------------------------------------------
!                     ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
  qf(1) = x / r
  qf(2) = -y / r
  sm(1) = -qf(1)
  sm(2) = -qf(2)
  tm(1) = sm(1)
  tm(2) = sm(2)
  pm = 1.0
  10 ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
  ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
  tm(1) = pm * ts(1)
  tm(2) = pm * ts(2)
  pm = pm + 1.0
  ts(1) = tm(1) / pm
  ts(2) = tm(2) / pm
  IF (ABS(sm(1))+ABS(ts(1)) == ABS(sm(1))) THEN
    IF (ABS(sm(2))+ABS(ts(2)) == ABS(sm(2))) GO TO 110
  END IF
  sm(1) = sm(1) + ts(1)
  sm(2) = sm(2) + ts(2)
  IF (pm < 36.0) GO TO 10
  GO TO 110
!-----------------------------------------------------------------------
!                      RATIONAL EXPANSION
!-----------------------------------------------------------------------
  20 sm(1) = 0.0
  sm(2) = 0.0
  DO i = 1, 18
    ts(1) = x - cd(i)
    ts(2) = y
    ss = ts(1) * ts(1) + ts(2) * ts(2)
    tm(1) = -ce(i) * ts(1) / ss
    tm(2) = ce(i) * ts(2) / ss
    sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
  END DO
  GO TO 110
!-----------------------------------------------------------------------
!                 EXPANSION INVOLVING EI AND DI
!-----------------------------------------------------------------------
  40 zl(1) = euler + 0.5 * LOG(r)
  zl(2) = ATAN2(-y,-x)

!                    SET SM = EXP(Z)*EI(-Z)

  sm(1) = 0.0
  sm(2) = 0.0
  DO i = 1, 18
    ts(1) = -x - ed(i)
    ts(2) = -y
    ss = ts(1) * ts(1) + ts(2) * ts(2)
    tm(1) = ee(i) * ts(1) / ss
    tm(2) = -ee(i) * ts(2) / ss
    sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
  END DO
  sr(1) = zl(1) * sm(1) - zl(2) * sm(2)
  sr(2) = zl(1) * sm(2) + zl(2) * sm(1)

!                    SET SM = EXP(Z)*DI(-Z)

  sm(1) = 0.0
  sm(2) = 0.0
  DO i = 1, 19
    ts(1) = -x - dd(i)
    ts(2) = -y
    ss = ts(1) * ts(1) + ts(2) * ts(2)
    tm(1) = de(i) * ts(1) / ss
    tm(2) = -de(i) * ts(2) / ss
    sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
  END DO
  ts(1) = -(x*sm(1)+y*sm(2)) / r
  ts(2) = -(x*sm(2)-y*sm(1)) / r

!                    COMPUTE THE EXPANSION

  sr(1) = sr(1) - ts(1) - zeta2
  sr(2) = sr(2) - ts(2)
  sm(1) = 0.0
  sm(2) = 0.0
  tm(1) = 1.0
  tm(2) = 0.0
  pm = 0.0
  qm = zeta2
  70 pm = pm + 1.0
  qm = qm - 1.0 / (pm*pm)
  ts(1) = tm(1) * x - tm(2) * y
  ts(2) = tm(1) * y + tm(2) * x
  tm(1) = ts(1) / pm
  tm(2) = ts(2) / pm
  ts(1) = qm * tm(1)
  ts(2) = qm * tm(2)
  IF (ABS(sm(1))+ABS(ts(1)) == ABS(sm(1))) THEN
    IF (ABS(sm(2))+ABS(ts(2)) == ABS(sm(2))) GO TO 80
  END IF
  sm(1) = sm(1) - ts(1)
  sm(2) = sm(2) - ts(2)
  GO TO 70
  80 sm(1) = sr(1) + sm(1)
  sm(2) = sr(2) + sm(2)

!                      COMPUTE EXP(-Z)*SM

  qm = EXP(-x)
  qf(1) = qm * COS(-y)
  qf(2) = qm * SIN(-y)
  ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
  ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
  sm(1) = ts(1)
  sm(2) = ts(2)
ELSE
!-----------------------------------------------------------------------
!                  SERIES FOR X*X + Y*Y <= 1
!-----------------------------------------------------------------------
  sm(1) = 0.0
  sm(2) = 0.0
  zl(1) = euler + 0.5 * LOG(r)
  zl(2) = ATAN2(-y,-x)
  sr(1) = const + 0.5 * (zl(1)*zl(1)-zl(2)*zl(2))
  sr(2) = zl(1) * zl(2)
  tm(1) = 1.0
  tm(2) = 0.0
  pm = 0.0
  qm = 0.0
  90 pm = pm + 1.0
  qm = qm + 1.0 / pm
  ts(1) = -tm(1) * x + tm(2) * y
  ts(2) = -tm(1) * y - tm(2) * x
  tm(1) = ts(1) / pm
  tm(2) = ts(2) / pm
  r = (zl(1)-1.0/pm) - qm
  ts(1) = (r*tm(1)-zl(2)*tm(2)) / pm
  ts(2) = (r*tm(2)+zl(2)*tm(1)) / pm
  IF (ABS(sm(1))+ABS(ts(1)) == ABS(sm(1))) THEN
    IF (ABS(sm(2))+ABS(ts(2)) == ABS(sm(2))) GO TO 100
  END IF
  sm(1) = sm(1) + ts(1)
  sm(2) = sm(2) + ts(2)
  GO TO 90
  100 sm(1) = sr(1) + sm(1)
  sm(2) = sr(2) + sm(2)
END IF
!-----------------------------------------------------------------------
!                         TERMINATION
!-----------------------------------------------------------------------
110 w = CMPLX(sm(1),sm(2))
RETURN
END SUBROUTINE cexexi



FUNCTION cli(z) RESULT(fn_val)
!     ******************************************************************
!     COMPUTATION OF THE COMPLEX LOGARITHMIC INTEGRAL
!     ******************************************************************
COMPLEX, INTENT(IN) :: z
COMPLEX             :: fn_val

! Local variables
REAL, PARAMETER :: qb(25) = (/ 2.77777777777778E-2, -1.00000000000000E-2,  &
        -1.70068027210884E-2, -1.94444444444444E-2, -2.06611570247934E-2,  &
        -2.14173006480699E-2, -2.19488663772311E-2, -2.23492338111715E-2,  &
        -2.26636891351914E-2, -2.29178211549926E-2, -2.31276449354844E-2,  &
        -2.33038680700203E-2, -2.34539766464373E-2, -2.35833786876607E-2,  &
        -2.36960832049849E-2, -2.37951264448373E-2, -2.38828504258091E-2,  &
        -2.39610907251825E-2, -2.40313063764460E-2, -2.40946717197585E-2,  &
        -2.41521426124012E-2, -2.42045049812210E-2, -2.42524109782181E-2,  &
        -2.42964062815807E-2, -2.43369509729144E-2 /),  &
        c = 1.64493406684823
REAL    :: az(2), pm, r, sm(2), tm(2), ts(2), sr(2), qf(2), dl(2), ds,   &
           zd(2), zl(2)
INTEGER :: n
!     ---------------------
!     C = PI**2/6
!     ---------------------
az(1) = REAL(z)
az(2) = AIMAG(z)
r = cpabs(az(1),az(2))
IF (r <= 0.5) THEN
  sr(1) = 0.0
  sr(2) = 0.0
  qf(1) = -az(1)
  qf(2) = -az(2)
  tm(1) = az(1)
  tm(2) = az(2)
ELSE

  IF (r >= 3.0) THEN
    zl(1) = LOG(r)
    zl(2) = ATAN2(az(2),az(1))
    sr(1) = c + 0.5 * (zl(1)*zl(1)-zl(2)*zl(2))
    sr(2) = zl(1) * zl(2)
    qf(1) = (-az(1)/r) / r
    qf(2) = (az(2)/r) / r
    tm(1) = qf(1)
    tm(2) = qf(2)
  ELSE

    zd(1) = 1.0 + az(1)
    zd(2) = az(2)
    ds = zd(1) * zd(1) + zd(2) * zd(2)
    IF (ds == 0.0) GO TO 50
    dl(1) = 0.5 * LOG(ds)
    dl(2) = ATAN2(zd(2),zd(1))
    IF (ds > 0.25) GO TO 20
    zl(1) = LOG(r)
    zl(2) = ATAN2(-az(2),-az(1))
    sr(1) = -c + (dl(1)*zl(1)-dl(2)*zl(2))
    sr(2) = dl(1) * zl(2) + dl(2) * zl(1)
    qf(1) = zd(1)
    qf(2) = zd(2)
    tm(1) = qf(1)
    tm(2) = qf(2)
  END IF
END IF

!             EVALUATION OF THE TAYLOR SERIES

sr(1) = sr(1) + tm(1)
sr(2) = sr(2) + tm(2)
sm(1) = 0.0
sm(2) = 0.0
pm = 1.0
10 pm = pm + 1.0
ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
tm(1) = ts(1)
tm(2) = ts(2)
ts(1) = tm(1) / (pm*pm)
ts(2) = tm(2) / (pm*pm)
IF (ABS(sm(1))+ABS(ts(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(ts(2)) == ABS(sm(2))) GO TO 40
END IF
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
GO TO 10

!        EVALUATION OF THE SERIES IN  U = -LN(1 + Z)

20 qf(1) = dl(1) * dl(1) - dl(2) * dl(2)
qf(2) = 2.0 * dl(1) * dl(2)
sr(1) = dl(1) + 0.25 * qf(1)
sr(2) = dl(2) + 0.25 * qf(2)
sm(1) = 0.0
sm(2) = 0.0
tm(1) = dl(1)
tm(2) = dl(2)
DO n = 1, 25
  ts(1) = qb(n) * (tm(1)*qf(1)-tm(2)*qf(2))
  ts(2) = qb(n) * (tm(1)*qf(2)+tm(2)*qf(1))
  tm(1) = ts(1)
  tm(2) = ts(2)
  IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
    IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 40
  END IF
  sm(1) = sm(1) + tm(1)
  sm(2) = sm(2) + tm(2)
END DO

40 fn_val = CMPLX(sr(1)+sm(1),sr(2)+sm(2))
RETURN

!                  EVALUATION AT Z = -1

50 fn_val = CMPLX(-c,0.0)
RETURN
END FUNCTION cli


FUNCTION ali(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!            COMPUTATION OF THE REAL DILOGARITHM FUNCTION
!-----------------------------------------------------------------------
!     WRITTEN BY
!        ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WEAPONS CENTER
!        DAHLGREN VIRGINIA
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: a(5) = (/ .217590467528373E+01, .165569610692639E+01,  &
        .522944061702389E+00, .626073688152965E-01, .187280204672313E-02 /),  &
        b(6) = (/ .242590467528371E+01, .215106116463796E+01,  &
        .853664388896516E+00, .148635712775060E+00, .936304016023909E-02,  &
        .115362459229893E-03 /), c(10) = (/ -.139792925233661E+01,  &
        .368504569727477E+00, .467406917183686E-01, .113795257294490E-01,  &
        .369638462505741E-02, .140888669464352E-02, .580505641503297E-03,  &
        .279065584075104E-03, .727678355839120E-04, .941452067850052E-04 /),  &
        d(2) = (/ -.164792925233634E+01, .669375771675355E+00 /), e(18)  &
        = (/ -.194565741631859E+00, -.430017756528812E-02,  &
        -.129188263110634E-03, -.344864872694838E-05, .566899694553089E-09,  &
        .126641834906132E-07, .163966793864421E-08, .164221074630109E-09,  &
        .149644905021032E-10, .130214292886747E-11, .110415518123737E-12,  &
        .921674760163207E-14, .761646464974859E-15, .625216733700975E-16,  &
        .510910937990370E-17, .416215390793180E-18, .338357379188308E-19,  &
        .274674744366340E-20 /), p(4) = (/ -.124827318209942E+01,  &
        -.593706951284264E-01, .368603360394688E-01, .243497524184253E-02 /),  &
        q(8) = (/ .100000000000000E+01, .252618047164349E+00,  &
        .171618729068655E-01, .234444792844727E-03, -.174928841869743E-05,  &
        .347369010951250E-07, -.713275908929482E-09, .958397514026421E-11 /),  &
        r(7) = (/ .265189940015693E+00, .230201018075415E+00,  &
        .315999623504943E-01, .154066621939470E-02, .286697611038892E-04,  &
        .163031291368652E-06, .838957807732251E-10 /), s(6)  &
        = (/ .100000000000000E+01, .177195068872258E+00,  &
        .110559275223905E-01, .291916852717175E-03, .304793254397420E-05,  &
        .882114921507386E-08 /), const = 1.64493406684823
REAL (dp), PARAMETER :: x0 = -12.5951703698450161286398965_dp
REAL (dp) :: t, w
INTEGER   :: i, l
!--------------------------
!     CONST = PI**2/6
!     X0 = ZERO OF THE REAL DILOGARITHM FUNCTION
!--------------------------

IF (x <= 1.0) THEN
  IF (x < 0.0) THEN
    IF (x < -0.5) THEN
      IF (x < -1.0) THEN
        IF (x < -2.0) THEN

!                            X < -2

          IF (x >= -26.63 .AND. x <= -6.97) GO TO 10
          t = -1.0 / x
          w = (((((((((((c(10)*t + c(9))*t + c(8))*t + c(7))*t + c(6))*t +  &
              c(5))*t + c(4))*t + c(3))*t + c(2))*t + c(1))*t+0.5)+0.5) /  &
              ((d(2)*t+d(1))*t+1.0)
          fn_val = 0.5 * LOG(-x) ** 2 - (2.0*const+w/x)
          RETURN
        END IF

!                        -2 <= X < -1

        t = -(1.0+x)
        w = (((((a(5)*t + a(4))*t + a(3))*t + a(2))*t + a(1))*t+1.0) / (((  &
            (((b(6)*t + b(5))*t + b(4))*t + b(3))*t + b(2))*t + b(1))*t+1.0)
        fn_val = -(const+t*w) + LOG(-x) * LOG(t)
        RETURN
      END IF

!                       -1 <= X < -1/2

      t = 0.5 + (0.5+x)
      fn_val = -const
      IF (t == 0.0) RETURN
      w = (((((((((((c(10)*t + c(9))*t + c(8))*t + c(7))*t + c(6))*t +  &
          c(5))*t + c(4))*t + c(3))*t + c(2))*t + c(1))*t+0.5)+0.5) /  &
          ((d(2)*t+d(1))*t+1.0)
      fn_val = (-const+t*w) + LOG(-x) * LOG(t)
      RETURN
    END IF

!                       -1/2 <= X < 0

    t = -x
    w = (((((((((((c(10)*t + c(9))*t + c(8))*t + c(7))*t + c(6))*t + c(5))*t  &
        + c(4))*t + c(3))*t + c(2))*t + c(1))*t+0.5)+0.5) / &
        ((d(2)*t+d(1))*t+1.0)
    fn_val = x * w
    RETURN
  END IF

!                         0 <= X <= 1

  w = (((((a(5)*x + a(4))*x + a(3))*x + a(2))*x + a(1))*x+1.0) /   &
      ((((((b(6)*x + b(5))*x + b(4))*x + b(3))*x + b(2))*x + b(1))*x + 1.0)
  fn_val = x * w
  RETURN
END IF

!                            X > 1

t = 1.0 / x
w = (((((a(5)*t + a(4))*t + a(3))*t + a(2))*t + a(1))*t+1.0) /   &
    ((((((b(6)*t + b(5))*t + b(4))*t + b(3))*t + b(2))*t + b(1))*t + 1.0)
fn_val = (const-w/x) + 0.5 * LOG(x) ** 2
RETURN

!-----------------------------------------------------------------------
!              EVALUATION FOR  -26.63 <= X <= -6.97
!-----------------------------------------------------------------------

10 IF (x > -14.0) THEN
  IF (x > -11.1) THEN

!                    -11.1 < X <= -6.97

    t = -(x+7.0)
    fn_val = (((p(4)*t + p(3))*t + p(2))*t + p(1)) / (((((((q(8)*t+q(7))*t + &
             q(6))*t+q(5))*t+q(4))*t+q(3))*t+q(2))*t+q(1))
    RETURN
  END IF


!                     -14 < X <= -11.1

  t = DBLE(x) - x0
  w = e(14)
  DO l = 1, 13
    i = 14 - l
    w = w * t + e(i)
  END DO
  fn_val = t * w
  RETURN
END IF

!                     -26.63 <= X <= -14

t = -(x+14.0)
fn_val = ((((((r(7)*t + r(6))*t + r(5))*t + r(4))*t + r(3))*t + r(2))*t + r(1)) /  &
         (((((s(6)*t + s(5))*t + s(4))*t + s(3))*t + s(2))*t + s(1))
RETURN
END FUNCTION ali


SUBROUTINE cgamma(mo, z, w)
!-----------------------------------------------------------------------

!        EVALUATION OF THE COMPLEX GAMMA AND LOGGAMMA FUNCTIONS

!                        ---------------

!     MO IS AN INTEGER, Z A COMPLEX ARGUMENT, AND W A COMPLEX VARIABLE.

!                 W = GAMMA(Z)       IF MO = 0
!                 W = LN(GAMMA(Z))   OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)  :: mo
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w

! Local variables
COMPLEX :: eta, eta2, sum
REAL, PARAMETER :: c0(12) = (/ .833333333333333E-01, -.277777777777778E-02,  &
        .793650793650794E-03, -.595238095238095E-03, .841750841750842E-03,  &
        -.191752691752692E-02, .641025641025641E-02, -.295506535947712E-01,  &
        .179644372368831E+00, -.139243221690590E+01, .134028640441684E+02,  &
        -.156848284626002E+03 /), pi = 3.14159265358979,  &
        pi2 = 6.28318530717959, alpi = 1.14472988584940,  &
        hl2p = .918938533204673
REAL    :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, s, sn, &
           s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
INTEGER :: j, k, l, m, max, n, nm1
!---------------------------
!     ALPI = LOG(PI)
!     HL2P = 0.5 * LOG(2*PI)
!---------------------------

!     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
!            BE USED, AND EPS IS THE SMALLEST REAL NUMBER
!            SUCH THAT 1.0 + EPS > 1.0.

!                      MAX = IPMPAR(3)
max = HUGE(3)
eps = EPSILON(1.0)

!---------------------------
x = REAL(z)
y = AIMAG(z)
IF (x < 0.0) THEN
!-----------------------------------------------------------------------
!            CASE WHEN THE REAL PART OF Z IS NEGATIVE
!-----------------------------------------------------------------------
  y = ABS(y)
  t = -pi * y
  et = EXP(t)
  e2t = et * et

!     SET  A1 = (1 + E2T)/2  AND  A2 = (1 - E2T)/2

  a1 = 0.5 * (1.0+e2t)
  t2 = t + t
  IF (t2 >= -0.15) THEN
    a2 = -0.5 * rexp(t2)
  ELSE
    a2 = 0.5 * (0.5+(0.5-e2t))
  END IF

!     COMPUTE SIN(PI*X) AND COS(PI*X)

  IF (ABS(x) >= MIN(REAL(MAX),1.0/eps)) GO TO 70
  k = ABS(x)
  u = x + k
  k = MOD(k,2)
  IF (u <= -0.5) THEN
    u = 0.5 + (0.5+u)
    k = k + 1
  END IF
  u = pi * u
  sn = SIN(u)
  cn = COS(u)
  IF (k == 1) THEN
    sn = -sn
    cn = -cn
  END IF

!     SET  H1 + H2*I  TO  PI/SIN(PI*Z)  OR  LOG(PI/SIN(PI*Z))

  a1 = sn * a1
  a2 = cn * a2
  a = a1 * a1 + a2 * a2
  IF (a == 0.0) GO TO 70
  IF (mo == 0) THEN

    h1 = a1 / a
    h2 = -a2 / a
    c = pi * et
    h1 = c * h1
    h2 = c * h2
  ELSE

    h1 = (alpi+t) - 0.5 * LOG(a)
    h2 = -ATAN2(a2,a1)
  END IF
  IF (AIMAG(z) >= 0.0) THEN
    x = 1.0 - x
    y = -y
  ELSE
    h2 = -h2
    x = 1.0 - x
  END IF
END IF
!-----------------------------------------------------------------------
!           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
!-----------------------------------------------------------------------
w1 = 0.0
w2 = 0.0
n = 0
t = x
y2 = y * y
a = t * t + y2
cut = 36.0
IF (eps > 1.e-8) cut = 16.0
IF (a < cut) THEN
  IF (a == 0.0) GO TO 70
  10 n = n + 1
  t = t + 1.0
  a = t * t + y2
  IF (a < cut) GO TO 10

!     LET S1 + S2*I BE THE PRODUCT OF THE TERMS (Z+J)/(Z+N)

  u1 = (x*t+y2) / a
  u2 = y / a
  s1 = u1
  s2 = n * u2
  IF (n >= 2) THEN
    u = t / a
    nm1 = n - 1
    DO j = 1, nm1
      v1 = u1 + j * u
      v2 = (n-j) * u2
      c = s1 * v1 - s2 * v2
      d = s1 * v2 + s2 * v1
      s1 = c
      s2 = d
    END DO
  END IF

!     SET  W1 + W2*I = LOG(S1 + S2*I)  WHEN MO IS NONZERO

  s = s1 * s1 + s2 * s2
  IF (mo /= 0) THEN
    w1 = 0.5 * LOG(s)
    w2 = ATAN2(s2,s1)
  END IF
END IF

!     SET  V1 + V2*I = (Z - 0.5) * LOG(Z + N) - Z

t1 = 0.5 * LOG(a) - 1.0
t2 = ATAN2(y,t)
u = x - 0.5
v1 = (u*t1-0.5) - y * t2
v2 = u * t2 + y * t1

!     LET A1 + A2*I BE THE ASYMPTOTIC SUM

eta = CMPLX(t/a,-y/a)
eta2 = eta * eta
m = 12
IF (a >= 289.0) m = 6
IF (eps > 1.e-8) m = m / 2
sum = CMPLX(c0(m),0.0)
l = m
DO j = 2, m
  l = l - 1
  sum = CMPLX(c0(l),0.0) + sum * eta2
END DO
sum = sum * eta
a1 = REAL(sum)
a2 = AIMAG(sum)
!-----------------------------------------------------------------------
!                 GATHERING TOGETHER THE RESULTS
!-----------------------------------------------------------------------
w1 = (((a1+hl2p)-w1)+v1) - n
w2 = (a2-w2) + v2
IF (REAL(z) < 0.0) GO TO 50
IF (mo == 0) THEN

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO = 0

  a = EXP(w1)
  w1 = a * COS(w2)
  w2 = a * SIN(w2)
  IF (n == 0) GO TO 60
  c = (s1*w1+s2*w2) / s
  d = (s1*w2-s2*w1) / s
  w1 = c
  w2 = d
  GO TO 60
END IF

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO IS NONZERO.
!     THE ANGLE W2 IS REDUCED TO THE INTERVAL -PI < W2 <= PI.

40 IF (w2 <= pi) THEN
  k = 0.5 - w2 / pi2
  w2 = w2 + pi2 * k
  GO TO 60
END IF
k = w2 / pi2 - 0.5
w2 = w2 - pi2 * REAL(k+1)
IF (w2 <= -pi) w2 = pi
GO TO 60

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO IS NONZERO

50 IF (mo /= 0) THEN
  w1 = h1 - w1
  w2 = h2 - w2
  GO TO 40
END IF

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO = 0

a = EXP(-w1)
t1 = a * COS(-w2)
t2 = a * SIN(-w2)
w1 = h1 * t1 - h2 * t2
w2 = h1 * t2 + h2 * t1
IF (n /= 0) THEN
  c = w1 * s1 - w2 * s2
  d = w1 * s2 + w2 * s1
  w1 = c
  w2 = d
END IF

!     TERMINATION

60 w = CMPLX(w1,w2)
RETURN
!-----------------------------------------------------------------------
!             THE REQUESTED VALUE CANNOT BE COMPUTED
!-----------------------------------------------------------------------
70 w = (0.0,0.0)
RETURN
END SUBROUTINE cgamma


FUNCTION cgam0(z) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF 1/GAMMA(1 + Z)  FOR ABS(Z) < 1.0
!-----------------------------------------------------------------------
COMPLEX, INTENT(IN) :: z
COMPLEX             :: fn_val

! Local variables
COMPLEX :: w
REAL, PARAMETER :: a(25) = (/ .577215664901533E+00, -.655878071520254E+00,  &
        -.420026350340952E-01, .166538611382291E+00, -.421977345555443E-01,  &
        -.962197152787697E-02, .721894324666310E-02, -.116516759185907E-02,  &
        -.215241674114951E-03, .128050282388116E-03, -.201348547807882E-04,  &
        -.125049348214267E-0, .113302723198170E-05, -.205633841697761E-0,  &
        .611609510448142E-08, .500200764446922E-08, -.118127457048702E-08,  &
        .104342671169110E-09, .778226343990507E-11, -.369680561864221E-11,  &
        .510037028745448E-12, -.205832605356651E-13, -.534812253942302E-14,  &
        .122677862823826E-14, -.118125930169746E-15 /)
INTEGER :: i, k, n
REAL    :: x, y
!-----------------------
n = 25
x = REAL(z)
y = AIMAG(z)
IF (x*x+y*y <= 0.04) n = 14

k = n
w = a(n)
DO i = 2, n
  k = k - 1
  w = a(k) + z * w
END DO
fn_val = 1.0 + z * w
RETURN
END FUNCTION cgam0


FUNCTION gamma(a) RESULT(fn_val)
!-----------------------------------------------------------------------

!         EVALUATION OF THE GAMMA FUNCTION FOR REAL ARGUMENTS

!                           -----------

!     GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
!     BE COMPUTED.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!          NAVAL SURFACE WEAPONS CENTER
!          DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: p(7) = (/ .539637273585445E-03, .261939260042690E-02,  &
        .204493667594920E-01, .730981088720487E-01, .279648642639792E+00,  &
        .553413866010467E+00, 1.0 /), q(7) = (/ -.832979206704073E-03,  &
        .470059485860584E-02, .225211131035340E-01, -.170458969313360E+00,  &
        -.567902761974940E-01, .113062953091122E+01, 1.0 /),  &
        pi = 3.1415926535898, r1 = .820756370353826E-03, r2  &
        = -.595156336428591E-03, r3 = .793650663183693E-03, r4  &
        = -.277777777770481E-02, r5 = .833333333333333E-01
REAL (dp), PARAMETER :: d = .41893853320467274178_dp
REAL    :: bot, g, lnx, s, t, top, w, x, z
INTEGER :: i, j, m, n
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
fn_val = 0.0
x = a
IF (ABS(a) < 15.0) THEN
!-----------------------------------------------------------------------
!            EVALUATION OF GAMMA(A) FOR ABS(A) < 15
!-----------------------------------------------------------------------
  t = 1.0
  m = INT(a) - 1

!     LET T BE THE PRODUCT OF A-J WHEN A >= 2

  IF (m < 0) GO TO 40

  DO j = 1, m
    x = x - 1.0
    t = x * t
  END DO
  x = x - 1.0
  GO TO 60

!     LET T BE THE PRODUCT OF A+J WHEN A < 1

  40 t = a
  IF (a <= 0.0) THEN
    m = -m - 1
    IF (m /= 0) THEN
      DO j = 1, m
        x = x + 1.0
        t = x * t
      END DO
    END IF
    x = (x+0.5) + 0.5
    t = x * t
    IF (t == 0.0) RETURN
  END IF


!     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!     CODE MAY BE OMITTED IF DESIRED.

  IF (ABS(t) < 1.e-30) THEN
    IF (ABS(t)*HUGE(1.0) <= 1.0001) RETURN
    fn_val = 1.0 / t
    RETURN
  END IF

!     COMPUTE GAMMA(1 + X) FOR  0 <= X < 1

  60 top = p(1)
  bot = q(1)
  DO i = 2, 7
    top = p(i) + x * top
    bot = q(i) + x * bot
  END DO
  fn_val = top / bot

!     TERMINATION

  IF (a >= 1.0) THEN
    fn_val = fn_val * t
    RETURN
  END IF
  fn_val = fn_val / t
  RETURN
END IF
!-----------------------------------------------------------------------
!            EVALUATION OF GAMMA(A) FOR ABS(A) >= 15
!-----------------------------------------------------------------------
IF (ABS(a) >= 1.e3) RETURN
IF (a <= 0.0) THEN
  x = -a
  n = x
  t = x - n
  IF (t > 0.9) t = 1.0 - t
  s = SIN(pi*t) / pi
  IF (MOD(n,2) == 0) s = -s
  IF (s == 0.0) RETURN
END IF

!     COMPUTE THE MODIFIED ASYMPTOTIC SUM

t = 1.0 / (x*x)
g = ((((r1*t + r2)*t + r3)*t + r4)*t + r5) / x

!     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = LOG(X)
!     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.

lnx = log(x)

!     FINAL ASSEMBLY

z = x
g = (d+g) + (z-0.5_dp) * (lnx-1._dp)
w = g
t = g - DBLE(w)
IF (w > 0.99999*exparg(0)) RETURN
fn_val = EXP(w) * (1.0+t)
IF (a < 0.0) fn_val = (1.0/(fn_val*s)) / x
RETURN
END FUNCTION gamma


FUNCTION glog(x) RESULT(fn_val)
!     -------------------
!     EVALUATION OF LN(X) FOR X >= 15
!     -------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: c1 = .286228750476730_dp, c2 = .399999628131494+dp,  &
                        c3 = .666666666752663_dp
REAL (dp), PARAMETER :: w(163) = (/ .270805020110221007D+01, .277258872223978124D+01,  &
  .283321334405621608D+01, .289037175789616469D+01, .294443897916644046D+01, .299573227355399099D+01, .304452243772342300D+01, &
  .309104245335831585D+01, .313549421592914969D+01, .317805383034794562D+01, .321887582486820075D+01, .325809653802148205D+01, &
  .329583686600432907D+01, .333220451017520392D+01, .336729582998647403D+01, .340119738166215538D+01, .343398720448514625D+01, &
  .346573590279972655D+01, .349650756146648024D+01, .352636052461616139D+01, .355534806148941368D+01, .358351893845611000D+01, &
  .361091791264422444D+01, .363758615972638577D+01, .366356164612964643D+01, .368887945411393630D+01, .371357206670430780D+01, &
  .373766961828336831D+01, .376120011569356242D+01, .378418963391826116D+01, .380666248977031976D+01, .382864139648909500D+01, &
  .385014760171005859D+01, .387120101090789093D+01, .389182029811062661D+01, .391202300542814606D+01, .393182563272432577D+01, &
  .395124371858142735D+01, .397029191355212183D+01, .398898404656427438D+01, .400733318523247092D+01, .402535169073514923D+01, &
  .404305126783455015D+01, .406044301054641934D+01, .407753744390571945D+01, .409434456222210068D+01, .411087386417331125D+01, &
  .412713438504509156D+01, .414313472639153269D+01, .415888308335967186D+01, .417438726989563711D+01, .418965474202642554D+01, &
  .420469261939096606D+01, .421950770517610670D+01, .423410650459725938D+01, .424849524204935899D+01, .426267987704131542D+01, &
  .427666611901605531D+01, .429045944114839113D+01, .430406509320416975D+01, .431748811353631044D+01, .433073334028633108D+01, &
  .434380542185368385D+01, .435670882668959174D+01, .436944785246702149D+01, .438202663467388161D+01, .439444915467243877D+01, &
  .440671924726425311D+01, .441884060779659792D+01, .443081679884331362D+01, .444265125649031645D+01, .445434729625350773D+01, &
  .446590811865458372D+01, .447733681447820647D+01, .448863636973213984D+01, .449980967033026507D+01, .451085950651685004D+01, &
  .452178857704904031D+01, .453259949315325594D+01, .454329478227000390D+01, .455387689160054083D+01, .456434819146783624D+01, &
  .457471097850338282D+01, .458496747867057192D+01, .459511985013458993D+01, .460517018598809137D+01, .461512051684125945D+01, &
  .462497281328427108D+01, .463472898822963577D+01, .464439089914137266D+01, .465396035015752337D+01, .466343909411206714D+01, &
  .467282883446190617D+01, .468213122712421969D+01, .469134788222914370D+01, .470048036579241623D+01, .470953020131233414D+01, &
  .471849887129509454D+01, .472738781871234057D+01, .473619844839449546D+01, .474493212836325007D+01, .475359019110636465D+01, &
  .476217393479775612D+01, .477068462446566476D+01, .477912349311152939D+01, .478749174278204599D+01, .479579054559674109D+01, &
  .480402104473325656D+01, .481218435537241750D+01, .482028156560503686D+01, .482831373730230112D+01, .483628190695147800D+01, &
  .484418708645859127D+01, .485203026391961717D+01, .485981240436167211D+01, .486753445045558242D+01, .487519732320115154D+01, &
  .488280192258637085D+01, .489034912822175377D+01, .489783979995091137D+01, .490527477843842945D+01, .491265488573605201D+01, &
  .491998092582812492D+01, .492725368515720469D+01, .493447393313069176D+01, .494164242260930430D+01, .494875989037816828D+01, &
  .495582705760126073D+01, .496284463025990728D+01, .496981329957600062D+01, .497673374242057440D+01, .498360662170833644D+01, &
  .499043258677873630D+01, .499721227376411506D+01, .500394630594545914D+01, .501063529409625575D+01, .501727983681492433D+01, &
  .502388052084627639D+01, .503043792139243546D+01, .503695260241362916D+01, .504342511691924662D+01, .504985600724953705D+01, &
  .505624580534830806D+01, .506259503302696680D+01, .506890420222023153D+01, .507517381523382692D+01, .508140436498446300D+01, &
  .508759633523238407D+01, .509375020080676233D+01, .509986642782419842D+01, .510594547390058061D+01, .511198778835654323D+01, &
  .511799381241675511D+01, .512396397940325892D+01, .512989871492307347D+01, .513579843705026176D+01, .514166355650265984D+01, &
  .514749447681345304D+01, .515329159449777895D+01, .515905529921452903D+01, .516478597392351405D+01, .517048399503815178D+01, &
  .517614973257382914D+01 /)
REAL (dp) :: t, t2, z
INTEGER   :: n
!     -------------------

IF (x < 178.0_dp) THEN
  n = x
  t = (x-n) / (x+n)
  t2 = t * t
  z = (((c1*t2 + c2)*t2 + c3)*t2 + 2.0) * t
  fn_val = w(n-14) + z
  RETURN
END IF

fn_val = LOG(x)
RETURN
END FUNCTION glog


FUNCTION gam1(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 <= A <= 1.5
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: p(7) = (/ .577215664901533E+00, -.409078193005776E+00,  &
        -.230975380857675E+00, .597275330452234E-01, .766968181649490E-02,  &
        -.514889771323592E-02, .589597428611429E-03 /),  &
        q(5) = (/ .100000000000000E+01, .427569613095214E+00,  &
        .158451672430138E+00, .261132021441447E-01, .423244297896961E-02 /),  &
        r(9) = (/ -.422784335098468E+00, -.771330383816272E+00,  &
        -.244757765222226E+00, .118378989872749E+00, .930357293360349E-03,  &
        -.118290993445146E-01, .223047661158249E-02, .266505979058923E-03,  &
        -.132674909766242E-03 /), s1 = .273076135303957E+00,   &
        s2 = .559398236957378E-01
REAL :: bot, d, t, top, w
!------------------------
t = a
d = a - 0.5
IF (d > 0.0) t = d - 0.5
IF (t < 0.0) THEN
  GO TO 30
ELSE IF (t > 0.0) THEN
  GO TO 20
END IF

fn_val = 0.0
RETURN

20 top = (((((p(7)*t + p(6))*t + p(5))*t + p(4))*t + p(3))*t + p(2)) * t + p(1)
bot = (((q(5)*t+q(4))*t+q(3))*t+q(2)) * t + 1.0
w = top / bot
IF (d <= 0.0) THEN
  fn_val = a * w
  RETURN
END IF
fn_val = (t/a) * ((w-0.5)-0.5)
RETURN

30 top = (((((((r(9)*t + r(8))*t + r(7))*t + r(6))*t + r(5))*t + r(4))*t + r(3))*  &
         t + r(2)) * t + r(1)
bot = (s2*t + s1) * t + 1.0
w = top / bot
IF (d <= 0.0) THEN
  fn_val = a * ((w+0.5)+0.5)
  RETURN
END IF
fn_val = t * w / a
RETURN
END FUNCTION gam1


FUNCTION gamln(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!          NAVAL SURFACE WARFARE CENTER
!          DAHLGREN, VIRGINIA
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------
REAL, INTENT(IN) :: a
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: d = .418938533204673, c0 = .833333333333333E-01,   &
              c1 = -.277777777760991E-02, c2 = .793650666825390E-03,   &
              c3 = -.595202931351870E-03, c4 = .837308034031215E-03,   &
              c5 = -.165322962780713E-02
REAL    :: t, w
INTEGER :: i, n
!-----------------------------------------------------------------------
IF (a <= 0.8) THEN
  fn_val = gamln1(a) - LOG(a)
  RETURN
END IF
IF (a <= 2.25) THEN
  t = (a-0.5) - 0.5
  fn_val = gamln1(t)
  RETURN
END IF

IF (a < 10.0) THEN
  n = a - 1.25
  t = a
  w = 1.0
  DO i = 1, n
    t = t - 1.0
    w = t * w
  END DO
  fn_val = gamln1(t-1.0) + LOG(w)
  RETURN
END IF

t = (1.0/a) ** 2
w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0) / a
fn_val = (d+w) + (a-0.5) * (LOG(a)-1.0)
RETURN
END FUNCTION gamln


FUNCTION gamln1(a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: p0 = .577215664901533E+00, p1 = .844203922187225E+00,  &
        p2 = -.168860593646662E+00, p3 = -.780427615533591E+00,   &
        p4 = -.402055799310489E+00, p5 = -.673562214325671E-01,   &
        p6 = -.271935708322958E-02, q1 = .288743195473681E+01,   &
        q2  = .312755088914843E+01, q3 = .156875193295039E+01,   &
        q4  = .361951990101499E+00, q5 = .325038868253937E-01,   &
        q6  = .667465618796164E-03, r0 = .422784335098467E+00,   &
        r1  = .848044614534529E+00, r2 = .565221050691933E+00,   &
        r3  = .156513060486551E+00, r4 = .170502484022650E-01,   &
        r5  = .497958207639485E-03, s1 = .124313399877507E+01,   &
        s2  = .548042109832463E+00, s3 = .101552187439830E+00,   &
        s4  = .713309612391000E-02, s5 = .116165475989616E-03
REAL :: w, x
!----------------------
IF (a < 0.6) THEN
  w = ((((((p6*a+p5)*a+p4)*a+p3)*a+p2)*a+p1)*a+p0) /   &
      ((((((q6*a+q5)*a+q4)*a+q3)*a+q2)*a+q1)*a+1.0)
  fn_val = -a * w
  RETURN
END IF

x = (a-0.5) - 0.5
w = (((((r5*x+r4)*x+r3)*x+r2)*x+r1)*x+r0) / (((((s5*x+s4)*x+s3)*x +  &
    s2)*x+s1)*x+1.0)
fn_val = x * w
RETURN
END FUNCTION gamln1


SUBROUTINE dcgama(mo, z, w)
!-----------------------------------------------------------------------

!        EVALUATION OF THE COMPLEX GAMMA AND LOGGAMMA FUNCTIONS

!                        ---------------

!     MO IS AN INTEGER. Z AND W ARE INTERPRETED AS REAL (dp)
!     COMPLEX NUMBERS. IT IS ASSUMED THAT Z(1) AND Z(2) ARE THE REAL
!     AND IMAGINARY PARTS OF THE COMPLEX NUMBER Z, AND THAT W(1) AND
!     W(2) ARE THE REAL AND IMAGINARY PARTS OF W.

!                 W = GAMMA(Z)       IF MO = 0
!                 W = LN(GAMMA(Z))   OTHERWISE

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)    :: mo
REAL (dp), INTENT(IN)  :: z(:)
REAL (dp), INTENT(OUT) :: w(:)

! Local variables
REAL (dp), PARAMETER :: c0(30)  &
        = (/ .8333333333333333333333333333333333333333D-01,  &
        -.2777777777777777777777777777777777777778D-02,  &
        .7936507936507936507936507936507936507937D-03,  &
        -.5952380952380952380952380952380952380952D-03,  &
        .8417508417508417508417508417508417508418D-03,  &
        -.1917526917526917526917526917526917526918D-02,  &
        .6410256410256410256410256410256410256410D-02,  &
        -.2955065359477124183006535947712418300654D-01,  &
        .1796443723688305731649384900158893966944D+00,  &
        -.1392432216905901116427432216905901116427D+01,  &
        .1340286404416839199447895100069013112491D+02,  &
        -.1568482846260020173063651324520889738281D+03,  &
        .2193103333333333333333333333333333333333D+04,  &
        -.3610877125372498935717326521924223073648D+05,  &
        .6914722688513130671083952507756734675533D+06,  &
        -.1523822153940741619228336495888678051866D+08,  &
        .3829007513914141414141414141414141414141D+09,  &
        -.1088226603578439108901514916552510537473D+11,  &
        .3473202837650022522522522522522522522523D+12,  &
        -.1236960214226927445425171034927132488108D+14,  &
        .4887880647930793350758151625180229021085D+15,  &
        -.2132033396091937389697505898213683855747D+17,  &
        .1021775296525700077565287628053585500394D+19,  &
        -.5357547217330020361082770919196920448485D+20,  &
        .3061578263704883415043151051329622758194D+22,  &
        -.1899991742639920405029371429306942902947D+24,  &
        .1276337403382883414923495137769782597654D+26,  &
        -.9252847176120416307230242348347622779519D+27,  &
        .7218822595185610297836050187301637922490D+29,  &
        -.6045183405995856967743148238754547286066D+31 /),  &
        dlpi = 1.144729885849400174143427351353058711647_dp,  &
        hl2p = .9189385332046727417803297364056176398614_dp,  &
        pi = 3.141592653589793238462643383279502884197_dp,  &
        pi2 = 6.283185307179586476925286766559005768394_dp
REAL (dp) :: a, a1, a2, c, cn, cut, d, eps, et, e2t, h1, h2, q1, q2, s, sn,  &
             s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, w2, x, y, y2
INTEGER   :: j, k, max, n, nm1
!---------------------------
!     DLPI = LOG(PI)
!     HL2P = 0.5 * LOG(2*PI)
!---------------------------

!     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
!            BE USED, AND EPS IS THE SMALLEST NUMBER SUCH
!            THAT  1._dp + EPS > 1._dp.

!                      MAX = IPMPAR(3)
max = HUGE(3)
eps = EPSILON(1.0_dp)

!---------------------------
x = z(1)
y = z(2)
IF (x < 0._dp) THEN
!-----------------------------------------------------------------------
!            CASE WHEN THE REAL PART OF Z IS NEGATIVE
!-----------------------------------------------------------------------
  y = ABS(y)
  t = -pi * y
  et = EXP(t)
  e2t = et * et

!     SET  A1 = (1 + E2T)/2  AND  A2 = (1 - E2T)/2

  a1 = 0.5_dp * (1._dp+e2t)
  t2 = t + t
  IF (t2 >= -0.15_dp) THEN
    a2 = -0.5_dp * drexp(t2)
  ELSE
    a2 = 0.5_dp * (0.5_dp+(0.5_dp-e2t))
  END IF

!     COMPUTE SIN(PI*X) AND COS(PI*X)

  u = MAX
  IF (ABS(x) >= MIN(u,1._dp/eps)) GO TO 80
  k = ABS(x)
  u = x + k
  k = MOD(k,2)
  IF (u <= -0.5_dp) THEN
    u = 0.5_dp + (0.5_dp+u)
    k = k + 1
  END IF
  u = pi * u
  sn = SIN(u)
  cn = COS(u)
  IF (k == 1) THEN
    sn = -sn
    cn = -cn
  END IF

!     SET  H1 + H2*I  TO  PI/SIN(PI*Z)  OR  LOG(PI/SIN(PI*Z))

  a1 = sn * a1
  a2 = cn * a2
  a = a1 * a1 + a2 * a2
  IF (a == 0._dp) GO TO 80
  IF (mo == 0) THEN

    h1 = a1 / a
    h2 = -a2 / a
    c = pi * et
    h1 = c * h1
    h2 = c * h2
  ELSE

    h1 = (dlpi+t) - 0.5_dp * LOG(a)
    h2 = -ATAN2(a2,a1)
  END IF
  IF (z(2) >= 0._dp) THEN
    x = 1.0 - x
    y = -y
  ELSE
    h2 = -h2
    x = 1.0 - x
  END IF
END IF
!-----------------------------------------------------------------------
!           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
!-----------------------------------------------------------------------
w1 = 0._dp
w2 = 0._dp
n = 0
t = x
y2 = y * y
a = t * t + y2
cut = 225._dp
IF (eps > 1.d-30) cut = 144._dp
IF (eps > 1.d-20) cut = 64._dp
IF (a < cut) THEN
  IF (a == 0._dp) GO TO 80
  10 n = n + 1
  t = t + 1._dp
  a = t * t + y2
  IF (a < cut) GO TO 10

!     LET S1 + S2*I BE THE PRODUCT OF THE TERMS (Z+J)/(Z+N)

  u1 = (x*t+y2) / a
  u2 = y / a
  s1 = u1
  s2 = n * u2
  IF (n >= 2) THEN
    u = t / a
    nm1 = n - 1
    DO j = 1, nm1
      v1 = u1 + j * u
      v2 = (n-j) * u2
      c = s1 * v1 - s2 * v2
      d = s1 * v2 + s2 * v1
      s1 = c
      s2 = d
    END DO
  END IF

!     SET  W1 + W2*I = LOG(S1 + S2*I)  WHEN MO IS NONZERO

  s = s1 * s1 + s2 * s2
  IF (mo /= 0) THEN
    w1 = 0.5_dp * LOG(s)
    w2 = ATAN2(s2,s1)
  END IF
END IF

!     SET  V1 + V2*I = (Z - 0.5) * LOG(Z + N) - Z

t1 = 0.5_dp * LOG(a) - 1._dp
t2 = ATAN2(y,t)
u = x - 0.5_dp
v1 = (u*t1-0.5_dp) - y * t2
v2 = u * t2 + y * t1

!     LET A1 + A2*I BE THE ASYMPTOTIC SUM

u1 = t / a
u2 = -y / a
q1 = u1 * u1 - u2 * u2
q2 = 2._dp * u1 * u2
a1 = 0._dp
a2 = 0._dp
DO j = 1, 30
  t1 = a1
  t2 = a2
  a1 = a1 + c0(j) * u1
  a2 = a2 + c0(j) * u2
  IF (a1 == t1) THEN
    IF (a2 == t2) GO TO 40
  END IF
  t1 = u1 * q1 - u2 * q2
  t2 = u1 * q2 + u2 * q1
  u1 = t1
  u2 = t2
END DO
!-----------------------------------------------------------------------
!                 GATHERING TOGETHER THE RESULTS
!-----------------------------------------------------------------------
40 w1 = (((a1+hl2p)-w1)+v1) - n
w2 = (a2-w2) + v2
IF (z(1) < 0._dp) GO TO 60
IF (mo == 0) THEN

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO = 0

  a = EXP(w1)
  w1 = a * COS(w2)
  w2 = a * SIN(w2)
  IF (n == 0) GO TO 70
  c = (s1*w1+s2*w2) / s
  d = (s1*w2-s2*w1) / s
  w1 = c
  w2 = d
  GO TO 70
END IF

!     CASE WHEN THE REAL PART OF Z IS NONNEGATIVE AND MO IS NONZERO.
!     THE ANGLE W2 IS REDUCED TO THE INTERVAL -PI < W2 <= PI.

50 IF (w2 <= pi) THEN
  k = 0.5_dp - w2 / pi2
  w2 = w2 + pi2 * k
  GO TO 70
END IF
k = w2 / pi2 - 0.5_dp
u = k + 1
w2 = w2 - pi2 * u
IF (w2 <= -pi) w2 = pi
GO TO 70

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO IS NONZERO

60 IF (mo /= 0) THEN
  w1 = h1 - w1
  w2 = h2 - w2
  GO TO 50
END IF

!     CASE WHEN THE REAL PART OF Z IS NEGATIVE AND MO = 0

a = EXP(-w1)
t1 = a * COS(-w2)
t2 = a * SIN(-w2)
w1 = h1 * t1 - h2 * t2
w2 = h1 * t2 + h2 * t1
IF (n /= 0) THEN
  c = w1 * s1 - w2 * s2
  d = w1 * s2 + w2 * s1
  w1 = c
  w2 = d
END IF

!     TERMINATION

70 w(1) = w1
w(2) = w2
RETURN
!-----------------------------------------------------------------------
!             THE REQUESTED VALUE CANNOT BE COMPUTED
!-----------------------------------------------------------------------
80 w(1) = 0._dp
w(2) = 0._dp
RETURN
END SUBROUTINE dcgama


FUNCTION dgamma(a) RESULT(fn_val)
!-----------------------------------------------------------------------

!                EVALUATION OF THE GAMMA FUNCTION FOR
!                     REAL (dp) ARGUMENTS

!                           -----------

!     DGAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT
!     BE COMPUTED.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!          NAVAL SURFACE WEAPONS CENTER
!          DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: d = 0.41893853320467274178032973640562_dp,  &
                        pi = 3.14159265358979323846264338327950_dp
REAL (dp) :: s, t, x, w
INTEGER   :: j, n
!-----------------------------------------------------------------------
!     D = 0.5*(LN(2*PI) - 1)
!-----------------------------------------------------------------------
fn_val = 0._dp
x = a
IF (ABS(a) <= 20._dp) THEN
!-----------------------------------------------------------------------
!             EVALUATION OF DGAMMA(A) FOR ABS(A) <= 20
!-----------------------------------------------------------------------
  t = 1._dp
  n = x
  n = n - 1

!     LET T BE THE PRODUCT OF A-J WHEN A >= 2

  IF (n < 0) THEN
    GO TO 40
  ELSE IF (n == 0) THEN
    GO TO 30
  END IF

  DO j = 1, n
    x = x - 1._dp
    t = x * t
  END DO
  30 x = x - 1._dp
  GO TO 60

!     LET T BE THE PRODUCT OF A+J WHEN A < 1

  40 t = a
  IF (a <= 0._dp) THEN
    n = -n - 1
    IF (n /= 0) THEN
      DO j = 1, n
        x = x + 1._dp
        t = x * t
      END DO
    END IF
    x = (x+0.5_dp) + 0.5_dp
    t = x * t
    IF (t == 0._dp) RETURN
  END IF


!     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!     CODE MAY BE OMITTED IF DESIRED.

  IF (ABS(t) < 1.d-33) THEN
    IF (ABS(t)*HUGE(1.0_dp) <= 1.000000001_dp) RETURN
    fn_val = 1._dp / t
    RETURN
  END IF

!     COMPUTE DGAMMA(1 + X) FOR 0 <= X < 1

  60 fn_val = 1._dp / (1._dp + dgam1(x))

!     TERMINATION

  IF (a >= 1._dp) THEN
    fn_val = fn_val * t
    RETURN
  END IF
  fn_val = fn_val / t
  RETURN
END IF
!-----------------------------------------------------------------------
!           EVALUATION OF DGAMMA(A) FOR ABS(A) > 20
!-----------------------------------------------------------------------
IF (ABS(a) >= 1.d3) RETURN
IF (a <= 0._dp) THEN
  s = dsin1(a) / pi
  IF (s == 0._dp) RETURN
  x = -a
END IF

!     COMPUTE THE MODIFIED ASYMPTOTIC SUM

w = dpdel(x)

!     FINAL ASSEMBLY

w = (d+w) + (x-0.5_dp) * (LOG(x)-1._dp)
IF (w > dxparg(0)) RETURN
fn_val = EXP(w)
IF (a < 0._dp) fn_val = (1._dp/(fn_val*s)) / x

RETURN
END FUNCTION dgamma


FUNCTION dpdel(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF THE FUNCTION DEL(X) FOR  X >= 10  WHERE
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X)

!                         --------

!     THE SERIES FOR DPDEL ON THE INTERVAL 0.0 TO 1.0 DERIVED BY
!     A.H. MORRIS FROM THE CHEBYSHEV SERIES IN THE SLATEC LIBRARY
!     OBTAINED BY WAYNE FULLERTON (LOS ALAMOS).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: a(15) = (/ .833333333333333333333333333333D-01,  &
        -.277777777777777777777777752282D-04,  &
        .793650793650793650791732130419D-07,  &
        -.595238095238095232389839236182D-09,  &
        .841750841750832853294451671990D-11,  &
        -.191752691751854612334149171243D-12,  &
        .641025640510325475730918472625D-14,  &
        -.295506514125338232839867823991D-15,  &
        .179643716359402238723287696452D-16,  &
        -.139228964661627791231203060395D-17,  &
        .133802855014020915603275339093D-18,  &
        -.154246009867966094273710216533D-19,  &
        .197701992980957427278370133333D-20,  &
        -.234065664793997056856992426667D-21,  &
        .171348014966398575409015466667D-22 /)
REAL (dp) :: t, w
INTEGER   :: i, k
!-----------------------------------------------------------------------
t = (10._dp/x) ** 2
w = a(15)
DO i = 1, 14
  k = 15 - i
  w = t * w + a(k)
END DO
fn_val = w / x
RETURN
END FUNCTION dpdel


FUNCTION dgam1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF 1/GAMMA(1 + X) - 1  FOR -0.5 <= X <= 1.5
!-----------------------------------------------------------------------

!     THE FOLLOWING ARE THE FIRST 49 COEFFICIENTS OF THE MACLAURIN
!     EXPANSION FOR 1/GAMMA(1 + X) - 1. THE COEFFICIENTS ARE
!     CORRECT TO 40 DIGITS. THE COEFFICIENTS WERE OBTAINED BY
!     ALFRED H. MORRIS JR. (NAVAL SURFACE WARFARE CENTER) AND ARE
!     GIVEN HERE FOR REFERENCE. ONLY THE FIRST 14 COEFFICIENTS ARE
!     USED IN THIS CODE.

!                           -----------

!     DATA A(1)  / .5772156649015328606065120900824024310422D+00/,
!    *     A(2)  /-.6558780715202538810770195151453904812798D+00/,
!    *     A(3)  /-.4200263503409523552900393487542981871139D-01/,
!    *     A(4)  / .1665386113822914895017007951021052357178D+00/,
!    *     A(5)  /-.4219773455554433674820830128918739130165D-01/,
!    *     A(6)  /-.9621971527876973562114921672348198975363D-02/,
!    *     A(7)  / .7218943246663099542395010340446572709905D-02/,
!    *     A(8)  /-.1165167591859065112113971084018388666809D-02/,
!    *     A(9)  /-.2152416741149509728157299630536478064782D-03/,
!    *     A(10) / .1280502823881161861531986263281643233949D-03/
!     DATA A(11) /-.2013485478078823865568939142102181838229D-04/,
!    *     A(12) /-.1250493482142670657345359473833092242323D-05/,
!    *     A(13) / .1133027231981695882374129620330744943324D-05/,
!    *     A(14) /-.2056338416977607103450154130020572836513D-06/,
!    *     A(15) / .6116095104481415817862498682855342867276D-08/,
!    *     A(16) / .5002007644469222930055665048059991303045D-08/,
!    *     A(17) /-.1181274570487020144588126565436505577739D-08/,
!    *     A(18) / .1043426711691100510491540332312250191401D-09/,
!    *     A(19) / .7782263439905071254049937311360777226068D-11/,
!    *     A(20) /-.3696805618642205708187815878085766236571D-11/
!     DATA A(21) / .5100370287454475979015481322863231802727D-12/,
!    *     A(22) /-.2058326053566506783222429544855237419746D-13/,
!    *     A(23) /-.5348122539423017982370017318727939948990D-14/,
!    *     A(24) / .1226778628238260790158893846622422428165D-14/,
!    *     A(25) /-.1181259301697458769513764586842297831212D-15/,
!    *     A(26) / .1186692254751600332579777242928674071088D-17/,
!    *     A(27) / .1412380655318031781555803947566709037086D-17/,
!    *     A(28) /-.2298745684435370206592478580633699260285D-18/,
!    *     A(29) / .1714406321927337433383963370267257066813D-19/,
!    *     A(30) / .1337351730493693114864781395122268022875D-21/
!     DATA A(31) /-.2054233551766672789325025351355733796682D-21/,
!    *     A(32) / .2736030048607999844831509904330982014865D-22/,
!    *     A(33) /-.1732356445910516639057428451564779799070D-23/,
!    *     A(34) /-.2360619024499287287343450735427531007926D-25/,
!    *     A(35) / .1864982941717294430718413161878666898946D-25/,
!    *     A(36) /-.2218095624207197204399716913626860379732D-26/,
!    *     A(37) / .1297781974947993668824414486330594165619D-27/,
!    *     A(38) / .1180697474966528406222745415509971518560D-29/,
!    *     A(39) /-.1124584349277088090293654674261439512119D-29/,
!    *     A(40) / .1277085175140866203990206677751124647749D-30/
!     DATA A(41) /-.7391451169615140823461289330108552823711D-32/,
!    *     A(42) / .1134750257554215760954165259469306393009D-34/,
!    *     A(43) / .4639134641058722029944804907952228463058D-34/,
!    *     A(44) /-.5347336818439198875077418196709893320905D-35/,
!    *     A(45) / .3207995923613352622861237279082794391090D-36/,
!    *     A(46) /-.4445829736550756882101590352124643637401D-38/,
!    *     A(47) /-.1311174518881988712901058494389922190237D-38/,
!    *     A(48) / .1647033352543813886818259327906394145400D-39/,
!    *     A(49) /-.1056233178503581218600561071538285049997D-40/

!                           -----------

!     C = A(1) - 1 IS ALSO FREQUENTLY NEEDED. C HAS THE VALUE ...

!     DATA C /-.4227843350984671393934879099175975689578D+00/

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: d, t, w, z
REAL (dp), PARAMETER :: a0 = .611609510448141581788D-08, a1  &
        = .624730830116465516210D-08, b1 = .203610414066806987300D+00, b2  &
        = .266205348428949217746D-01, b3 = .493944979382446875238D-03, b4  &
        = -.851419432440314906588D-05, b5 = -.643045481779353022248D-05, b6  &
        = .992641840672773722196D-06, b7 = -.607761895722825260739D-07, b8  &
        = .195755836614639731882D-09
REAL (dp), PARAMETER :: p0 = .6116095104481415817861D-08, p1  &
        = .6871674113067198736152D-08, p2 = .6820161668496170657, p3  &
        = .4686843322948848031080D-10, p4 = .1572833027710446286995D-11, p5  &
        = -.1249441572276366213222D-12, p6 = .4343529937408594255178D-14, q1  &
        = .3056961078365221025009D+00, q2 = .5464213086042296536016D-01, q3  &
        = .4956830093825887312, q4 = .2692369466186361192876D-03
REAL (dp), PARAMETER :: c = -.422784335098467139393487909917598D+00, c0  &
        = .577215664901532860606512090082402D+00, c1  &
        = -.655878071520253881077019515145390D+00, c2  &
        = -.420026350340952355290039348754298D-01, c3  &
        = .166538611382291489501700795102105D+00, c4  &
        = -.421977345555443367482083012891874D-01, c5  &
        = -.962197152787697356211492167234820D-02, c6  &
        = .721894324666309954239501034044657D-02, c7  &
        = -.116516759185906511211397108401839D-02, c8  &
        = -.215241674114950972815729963053648D-03, c9  &
        = .128050282388116186153198626328164D-03, c10  &
        = -.201348547807882386556893914210218D-04, c11  &
        = -.125049348214267065734535947383309D-05, c12  &
        = .113302723198169588237412962033074D-05, c13  &
        = -.205633841697760710345015413002057D-06
!----------------------------
t = x
d = x - 0.5_dp
IF (d > 0._dp) t = d - 0.5_dp
IF (t < 0.0_dp) THEN
  GO TO 30
ELSE IF (t > 0.0_dp) THEN
  GO TO 20
END IF

fn_val = 0._dp
RETURN
!------------

!                CASE WHEN 0 < T <= 0.5

!              W IS A MINIMAX APPROXIMATION FOR
!              THE SERIES A(15) + A(16)*T + ...

!------------
20 w = ((((((p6*t + p5)*t + p4)*t + p3)*t + p2)*t + p1)*t + p0) /   &
       ((((q4*t+q3)*t + q2)*t + q1)*t + 1._dp)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c0

IF (d <= 0._dp) THEN
  fn_val = x * z
  RETURN
END IF
fn_val = (t/x) * ((z-0.5_dp)-0.5_dp)
RETURN
!------------

!                CASE WHEN -0.5 <= T < 0

!              W IS A MINIMAX APPROXIMATION FOR
!              THE SERIES A(15) + A(16)*T + ...

!------------
30 w = (a1*t + a0) / ((((((((b8*t + b7)*t + b6)*t + b5)*t + b4)*t + b3)*t + b2)*t + b1)*t+1._dp)
z = (((((((((((((w*t + c13)*t + c12)*t + c11)*t + c10)*t + c9)*t + c8)*t + c7)*t  &
    + c6)*t + c5)*t + c4)*t + c3)*t + c2)*t + c1) * t + c

IF (d <= 0._dp) THEN
  fn_val = x * ((z+0.5_dp)+0.5_dp)
  RETURN
END IF
fn_val = t * z / x
RETURN
END FUNCTION dgam1


FUNCTION dgamln(a) RESULT(fn_val)
!-----------------------------------------------------------------------

!           EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS
!          NAVAL SURFACE WEAPONS CENTER
!          DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
!     D = 0.5*(LN(2*PI) - 1)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: d = 0.41893853320467274178032973640562_dp
REAL (dp) :: w, x
INTEGER   :: i, n
!--------------------------
IF (a < 0.5_dp) THEN
  fn_val = dgmln1(a) - LOG(a)
  RETURN
END IF
IF (a <= 2.5_dp) THEN
  x = a - 1._dp
  IF (a < 1._dp) x = (a-0.5_dp) - 0.5_dp
  fn_val = dgmln1(x)
  RETURN
END IF

IF (a < 10._dp) THEN
  n = a - 1.5_dp
  x = a
  w = 1._dp
  DO i = 1, n
    x = x - 1._dp
    w = x * w
  END DO
  fn_val = dgmln1(x-1._dp) + LOG(w)
  RETURN
END IF

w = dpdel(a)
fn_val = (d+w) + (a-0.5_dp) * (LOG(a)-1._dp)
RETURN
END FUNCTION dgamln


FUNCTION dgmln1(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + X)) FOR -0.5 <= X <= 1.5
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: w
!-----------------------
w = dgam1(x)
fn_val = -dlnrel(w)
RETURN
END FUNCTION dgmln1


SUBROUTINE cpsi(z, w)
!-----------------------------------------------------------------------
!           EVALUATION OF THE COMPLEX DIGAMMA FUNCTION
!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: c0(12) = (/ .833333333333333E-01, -.833333333333333E-02,  &
        .396825396825397E-02, -.416666666666667E-02, .757575757575758E-02,  &
        -.210927960927961E-01, .833333333333333E-01, -.443259803921569E+00,  &
        .305395433027012E+01, -.264562121212121E+02, .281460144927536E+03,  &
        -.360751054639805E+04 /), pi = 3.14159265358979324,  &
        pi2 = 6.28318530717958648
COMPLEX   :: eta, eta2, sum
REAL (dp) :: ds1, ds2
REAL      :: a, a1, a2, c1, c2, cn, eps, et, h1, h2, s, s1, s2, sn, t,  &
             u, w1, w2, x, y, y2
INTEGER   :: j, k, l, m, max
!----------------------------
!     PI2 = 2*PI
!----------------------------

!     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
!            BE USED, AND EPS IS THE SMALLEST REAL NUMBER
!            SUCH THAT 1.0 + EPS > 1.0.

!                      MAX = IPMPAR(3)
max = HUGE(3)
eps = EPSILON(1.0)

!----------------------------
x = REAL(z)
y = AIMAG(z)
IF (x < 0.0) THEN
!-----------------------------------------------------------------------
!            CASE WHEN THE REAL PART OF Z IS NEGATIVE
!-----------------------------------------------------------------------
  y = ABS(y)
  t = -pi2 * y
  et = EXP(t)

!     SET  A1 = (1 + ET)/2  AND  A2 = (1 - ET)/2

  a1 = 0.5 * (1.0+et)
  IF (t >= -0.15) THEN
    a2 = -0.5 * rexp(t)
  ELSE
    a2 = 0.5 * (0.5+(0.5-et))
  END IF

!     COMPUTE SIN(PI*X) AND COS(PI*X), OR -SIN(PI*X) AND -COS(PI*X)

  IF (ABS(x) >= MIN(REAL(MAX),1.0/eps)) GO TO 30
  k = ABS(x)
  u = x + k
  IF (u <= -0.5) u = 0.5 + (0.5+u)
  u = pi * u
  sn = SIN(u)
  cn = COS(u)

!     SET H1 + H2*I = PI*COT(PI*Z)

  s1 = a1 * sn
  s2 = a2 * cn
  c1 = a1 * cn
  c2 = -a2 * sn
  s = s1 * s1 + s2 * s2
  h1 = pi * (s1*c1+s2*c2) / s
  h2 = pi * (s1*c2-s2*c1) / s

  IF (AIMAG(z) >= 0.0) THEN
    x = 1.0 - x
    y = -y
  ELSE
    h2 = -h2
    x = 1.0 - x
  END IF
END IF
!-----------------------------------------------------------------------
!           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
!-----------------------------------------------------------------------
t = x
y2 = y * y
a = x * x + y2
IF (a /= 0.0) THEN

!     LET S1 + S2*I BE THE SUM OF THE TERMS 1/(Z+J) FOR J = 0,1,...,N-1

  ds1 = 0._dp
  ds2 = 0._dp
  10 IF (a < 36.0) THEN
    ds1 = ds1 + DBLE(t/a)
    ds2 = ds2 - DBLE(y/a)
    t = t + 1.0
    a = t * t + y2
    GO TO 10
  END IF
  s1 = ds1
  s2 = ds2

!     SET W1 + W2*I = LOG(Z+N)

  w1 = 0.5 * LOG(a)
  w2 = ATAN2(y,t)

!     LET A1 + A2*I BE THE ASYMPTOTIC SUM

  eta = CMPLX(t/a,-y/a)
  eta2 = eta * eta
  m = 12
  l = m
  sum = CMPLX(c0(m),0.0)
  DO j = 2, m
    l = l - 1
    sum = CMPLX(c0(l),0.0) + sum * eta2
  END DO
  sum = CMPLX(0.5,0.0) * eta + eta2 * sum
  a1 = REAL(sum)
  a2 = AIMAG(sum)
!-----------------------------------------------------------------------
!                 GATHERING TOGETHER THE RESULTS
!-----------------------------------------------------------------------
  w1 = (w1-s1) - a1
  w2 = (w2-a2) - s2
  w = CMPLX(w1,w2)
  IF (REAL(z) >= 0.0) RETURN
  w = CMPLX(w1-h1,w2-h2)
  RETURN
END IF
!-----------------------------------------------------------------------
!             THE REQUESTED VALUE CANNOT BE COMPUTED
!-----------------------------------------------------------------------
30 w = (0.0,0.0)
RETURN
END SUBROUTINE cpsi


FUNCTION psi(xx) RESULT(fn_val)
!---------------------------------------------------------------------

!                 EVALUATION OF THE DIGAMMA FUNCTION

!                           -----------

!     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
!     BE COMPUTED.

!     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
!     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
!     CODY, STRECOK AND THACHER.

!---------------------------------------------------------------------
!     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
!     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY
!     A.H. MORRIS (NSWC).
!---------------------------------------------------------------------
REAL, INTENT(IN) :: xx
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: p1(7) = (/ .895385022981970E-02, .477762828042627E+01,  &
        .142441585084029E+03, .118645200713425E+04, .363351846806499E+04,  &
        .413810161269013E+04, .130560269827897E+04 /), p2(4)  &
        = (/ -.212940445131011E+01, -.701677227766759E+0,  &
        -.448616543918019E+01, -.648157123766197E+00 /),  &
        q1(6) = (/ .448452573429826E+02, .520752771467162E+03,  &
        .221000799247830E+04, .364127349079381E+04, .190831076596300E+04,  &
        .691091682714533E-05 /), q2(4) = (/ .322703493791143E+02,  &
        .892920700481861E+02, .546117738103215E+02, .777788548522962E+01 /),  &
        piov4 = .785398163397448
REAL (dp), PARAMETER :: dx0 = 1.461632144968362341262659542325721325_dp
REAL    :: aug, den, sgn, upper, w, x, xmax1, xmx0, xsmall, z
INTEGER :: i, m, n, nq
!---------------------------------------------------------------------

!     PIOV4 = PI/4
!     DX0 = ZERO OF PSI TO EXTENDED PRECISION

!---------------------------------------------------------------------

!     MACHINE DEPENDENT CONSTANTS ...

!        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
!                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED
!                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
!                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
!                 PSI MAY BE REPRESENTED AS LOG(X).

!        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
!                 MAY BE REPRESENTED BY 1/X.

!---------------------------------------------------------------------
!      XMAX1 = IPMPAR(3)
xmax1 = HUGE(3)
xmax1 = MIN(xmax1, 1.0/EPSILON(1.0))
xsmall = 1.e-9
!---------------------------------------------------------------------
x = xx
aug = 0.0
IF (x < 0.5) THEN
!---------------------------------------------------------------------
!     X < 0.5,  USE REFLECTION FORMULA
!     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!---------------------------------------------------------------------
  IF (ABS(x) <= xsmall) THEN
    IF (x == 0.0) GO TO 30
!---------------------------------------------------------------------
!     0 < ABS(X) <= XSMALL.  USE 1/X AS A SUBSTITUTE
!     FOR  PI*COTAN(PI*X)
!---------------------------------------------------------------------
    aug = -1.0 / x
  ELSE
!---------------------------------------------------------------------
!     REDUCTION OF ARGUMENT FOR COTAN
!---------------------------------------------------------------------
    w = -x
    sgn = piov4
    IF (w <= 0.0) THEN
      w = -w
      sgn = -sgn
    END IF
!---------------------------------------------------------------------
!     MAKE AN ERROR EXIT IF X <= -XMAX1
!---------------------------------------------------------------------
    IF (w >= xmax1) GO TO 30
    nq = INT(w)
    w = w - REAL(nq)
    nq = INT(w*4.0)
    w = 4.0 * (w-REAL(nq)*.25)
!---------------------------------------------------------------------
!     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
!     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
!     QUADRANT AND DETERMINE SIGN
!---------------------------------------------------------------------
    n = nq / 2
    IF ((n+n) /= nq) w = 1.0 - w
    z = piov4 * w
    m = n / 2
    IF ((m+m) /= n) sgn = -sgn
!---------------------------------------------------------------------
!     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
!---------------------------------------------------------------------
    n = (nq+1) / 2
    m = n / 2
    m = m + m
    IF (m == n) THEN
!---------------------------------------------------------------------
!     CHECK FOR SINGULARITY
!---------------------------------------------------------------------
      IF (z == 0.0) GO TO 30
!---------------------------------------------------------------------
!     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
!     SIN/COS AS A SUBSTITUTE FOR TAN
!---------------------------------------------------------------------
      aug = sgn * ((COS(z)/SIN(z))*4.0)
    ELSE
      aug = sgn * ((SIN(z)/COS(z))*4.0)
    END IF
  END IF
  x = 1.0 - x
END IF
IF (x <= 3.0) THEN
!---------------------------------------------------------------------
!     0.5 <= X <= 3.0
!---------------------------------------------------------------------
  den = x
  upper = p1(1) * x

  DO i = 1, 5
    den = (den+q1(i)) * x
    upper = (upper+p1(i+1)) * x
  END DO

  den = (upper+p1(7)) / (den+q1(6))
  xmx0 = DBLE(x) - dx0
  fn_val = den * xmx0 + aug
  RETURN
END IF
!---------------------------------------------------------------------
!     IF X >= XMAX1, PSI = LN(X)
!---------------------------------------------------------------------
IF (x < xmax1) THEN
!---------------------------------------------------------------------
!     3.0 < X < XMAX1
!---------------------------------------------------------------------
  w = 1.0 / (x*x)
  den = w
  upper = p2(1) * w

  DO i = 1, 3
    den = (den+q2(i)) * w
    upper = (upper+p2(i+1)) * w
  END DO

  aug = upper / (den+q2(4)) - 0.5 / x + aug
END IF
fn_val = aug + LOG(x)
RETURN
!---------------------------------------------------------------------
!     ERROR RETURN
!---------------------------------------------------------------------
30 fn_val = 0.0
RETURN
END FUNCTION psi



SUBROUTINE dcpsi(z, w)
!-----------------------------------------------------------------------
!           EVALUATION OF THE COMPLEX DIGAMMA FUNCTION
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: z(:)
REAL (dp), INTENT(OUT) :: w(:)

! Local variables
REAL (dp), PARAMETER :: c0(30)  &
        = (/ .8333333333333333333333333333333333333333D-01,  &
        -.8333333333333333333333333333333333333333D-02,  &
        .3968253968253968253968253968253968253968D-02,  &
        -.4166666666666666666666666666666666666667D-02,  &
        .7575757575757575757575757575757575757576D-02,  &
        -.2109279609279609279609279609279609279609D-01,  &
        .8333333333333333333333333333333333333333D-01,  &
        -.4432598039215686274509803921568627450980D+00,  &
        .3053954330270119743803954330270119743804D+01,  &
        -.2645621212121212121212121212121212121212D+02,  &
        .2814601449275362318840579710144927536232D+03,  &
        -.3607510546398046398046398046398046398046D+04,  &
        .5482758333333333333333333333333333333333D+05,  &
        -.9749368238505747126436781609195402298851D+06,  &
        .2005269579668807894614346227249453055905D+08,  &
        -.4723848677216299019607843137254901960784D+09,  &
        .1263572479591666666666666666666666666667D+11,  &
        -.3808793112524536881155302207933786881155D+12,  &
        .1285085049930508333333333333333333333333D+14,  &
        -.4824144835485017037158167036215816703622D+15,  &
        .2004031065651625273810842166323893898645D+17,  &
        -.9167743603195330775699275362318840579710D+18,  &
        .4597988834365650349043794326241134751773D+20,  &
        -.2518047192145109569708902332022552610788D+22,  &
        .1500173349215392873371144015151515151515D+24,  &
        -.9689957887463594065649794289465408805031D+25,  &
        .6764588237929282099094524230179847767567D+27,  &
        -.5089065946866228968976633291591192528736D+29,  &
        .4114728879255797869766548606761933615819D+31,  &
        -.3566658209537555610968457460865182898779D+33 /),  &
        pi = 3.141592653589793238462643383279502884197_dp, pi2 = 6.283185307179586476925286766
REAL (dp) :: a, a1, a2, cn, cut, c1, c2, eps, et, h1, h2, q1,  &
             q2, s, sn, s1, s2, t, t1, t2, u, u1, u2, v1, v2, w1, &
             w2, x, y,y2
INTEGER   :: j, k, m, max
!----------------------------
!     PI2 = 2*PI
!----------------------------

!     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
!            BE USED, AND EPS IS THE SMALLEST REAL NUMBER
!            SUCH THAT 1._dp + EPS > 1._dp.

!                      MAX = IPMPAR(3)
max = HUGE(3)
eps = EPSILON(1.0_dp)

!----------------------------
x = z(1)
y = z(2)
IF (x < 0._dp) THEN
!-----------------------------------------------------------------------
!            CASE WHEN THE REAL PART OF Z IS NEGATIVE
!-----------------------------------------------------------------------
  y = ABS(y)
  t = -pi2 * y
  et = EXP(t)

!     SET  A1 = (1 + ET)/2  AND  A2 = (1 - ET)/2

  a1 = 0.5_dp * (1._dp+et)
  IF (t >= -0.15_dp) THEN
    a2 = -0.5_dp * drexp(t)
  ELSE
    a2 = 0.5_dp * (0.5_dp+(0.5_dp-et))
  END IF

!     COMPUTE SIN(PI*X) AND COS(PI*X), OR -SIN(PI*X) AND -COS(PI*X)

  u = MAX
  IF (ABS(x) >= MIN(u,1._dp/eps)) GO TO 40
  k = ABS(x)
  u = x + k
  IF (u <= -0.5_dp) u = 0.5_dp + (0.5_dp+u)
  u = pi * u
  sn = SIN(u)
  cn = COS(u)

!     SET H1 + H2*I = PI*COT(PI*Z)

  s1 = a1 * sn
  s2 = a2 * cn
  c1 = a1 * cn
  c2 = -a2 * sn
  s = s1 * s1 + s2 * s2
  h1 = pi * (s1*c1+s2*c2) / s
  h2 = pi * (s1*c2-s2*c1) / s

  IF (z(2) >= 0._dp) THEN
    x = 1._dp - x
    y = -y
  ELSE
    h2 = -h2
    x = 1._dp - x
  END IF
END IF
!-----------------------------------------------------------------------
!           CASE WHEN THE REAL PART OF Z IS NONNEGATIVE
!-----------------------------------------------------------------------
t = x
y2 = y * y
a = x * x + y2
IF (a /= 0._dp) THEN
  cut = 225._dp
  IF (eps > 1.d-30) cut = 144._dp

!     LET S1 + S2*I BE THE SUM OF THE TERMS 1/(Z+J) FOR J = 0,1,...,N-1

  s1 = 0._dp
  s2 = 0._dp
  10 IF (a < cut) THEN
    s1 = s1 + t / a
    s2 = s2 - y / a
    t = t + 1._dp
    a = t * t + y2
    GO TO 10
  END IF

!     SET W1 + W2*I = LOG(Z+N)

  w1 = 0.5_dp * LOG(a)
  w2 = ATAN2(y,t)

!     LET A1 + A2*I BE THE ASYMPTOTIC SUM

  u1 = t / a
  u2 = -y / a
  q1 = u1 * u1 - u2 * u2
  q2 = 2._dp * u1 * u2
  v1 = q1
  v2 = q2
  a1 = 0._dp
  a2 = 0._dp
  m = 30
  IF (eps > 1.d-30) m = 25
  DO j = 1, m
    t1 = a1
    t2 = a2
    a1 = a1 + c0(j) * v1
    a2 = a2 + c0(j) * v2
    IF (a1 == t1) THEN
      IF (a2 == t2) GO TO 30
    END IF
    t1 = v1 * q1 - v2 * q2
    t2 = v1 * q2 + v2 * q1
    v1 = t1
    v2 = t2
  END DO
!-----------------------------------------------------------------------
!                 GATHERING TOGETHER THE RESULTS
!-----------------------------------------------------------------------
  30 a1 = a1 + 0.5_dp * u1
  a2 = a2 + 0.5_dp * u2
  w(1) = (w1-a1) - s1
  w(2) = (w2-a2) - s2
  IF (z(1) >= 0._dp) RETURN
  w(1) = w(1) - h1
  w(2) = w(2) - h2
  RETURN
END IF
!-----------------------------------------------------------------------
!             THE REQUESTED VALUE CANNOT BE COMPUTED
!-----------------------------------------------------------------------
40 w(1) = 0._dp
w(2) = 0._dp
RETURN
END SUBROUTINE dcpsi



FUNCTION dpsi(a) RESULT(fn_val)
!-----------------------------------------------------------------------

!               EVALUATION OF THE DIGAMMA FUNCTION FOR
!                     REAL (dp) ARGUMENTS

!                           -----------

!     DPSI(A) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
!     BE COMPUTED.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!          NAVAL SURFACE WARFARE CENTER
!          DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------

!     THE SERIES FOR DPSI ON THE INTERVAL 0.0 TO 1.0 WAS DERIVED
!     BY WAYNE FULLERTON (LOS ALAMOS NATIONAL LABORATORY).

!                                  WITH WEIGHTED ERROR   5.79E-32
!                                   LOG WEIGHTED ERROR  31.24
!                         SIGNIFICANT FIGURES REQUIRED  30.93
!                              DECIMAL PLACES REQUIRED  32.05

!     THE SERIES FOR  A >= 10  WAS DERIVED BY A.H. MORRIS FROM
!     THE CHEBYSHEV SERIES IN THE SLATEC LIBRARY OBTAINED BY WAYNE
!     FULLERTON (LOS ALAMOS).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: c(42) = (/ -.38057080835217921520437677667039D-01,  &
    .49141539302938712748204699654277D+00,  &
    -.56815747821244730242892064734081D-01, .83578212259143131362775650747862D-02,  &
    -.13332328579943425998079274172393D-02, .22031328706930824892872397979521D-03,  &
    -.37040238178456883592889086949229D-04, .62837936548549898933651418717690D-05,  &
    -.10712639085061849855283541747074D-05, .18312839465484165805731589810378D-06,  &
    -.31353509361808509869005779796885D-07, .53728087762007766260471919143615D-08,  &
    -.92116814159784275717880632624730D-09, .15798126521481822782252884032823D-09,  &
    -.27098646132380443065440589409707D-10, .46487228599096834872947319529549D-11,  &
    -.79752725638303689726504797772737D-12, .13682723857476992249251053892838D-12,  &
    -.23475156060658972717320677980719D-13, .40276307155603541107907925006281D-14,  &
    -.69102518531179037846547422974771D-15, .11856047138863349552929139525768D-15,  &
    -.20341689616261559308154210484223D-16, .34900749686463043850374232932351D-17,  &
    -.59880146934976711003011081393493D-18, .10273801628080588258398005712213D-18,  &
    -.17627049424561071368359260105386D-19, .30243228018156920457454035490133D-20,  &
    -.51889168302092313774286088874666D-21, .89027730345845713905005887487999D-22,  &
    -.15274742899426728392894971904000D-22, .26207314798962083136358318079999D-23,  &
    -.44964642738220696772598388053333D-24, .77147129596345107028919364266666D-25,  &
    -.13236354761887702968102638933333D-25, .22709994362408300091277311999999D-26,  &
    -.38964190215374115954491391999999D-27, .66851981388855302310679893333333D-28,  &
    -.11469986654920864872529919999999D-28, .19679385886541405920515413333333D-29,  &
    -.33764488189750979801907200000000D-30, .57930703193214159246677333333333D-31 /)
REAL (dp), PARAMETER :: p(15) = (/ .833333333333333333333333333147D-03,  &
    -.833333333333333333333317475057D-06,  &
    .396825396825396825343072884056D-08,  &
    -.416666666666666570859890514548D-10,  &
    .757575757575654146210665696401D-12,  &
    -.210927960920616064592099772274D-13,  &
    .833333329719356554828382131321D-15,  &
    -.443259676504784387819140445894D-16,  &
    .305392145578967948828783519552D-17,  &
    -.264499326810660590871410866039D-18,  &
    .280568932535744579536244004181D-19,  &
    -.351388195869099967789469969066D-20,  &
    .476233402067211507540059750399D-21,  &
    -.575024569953144855161645738666D-22,  &
    .416180125797657207803740160000D-23 /),  &
    pi = 3.1415926535897932384626433832795_dp,  &
    x0 = .46163214496836234126265954232572_dp
REAL (dp) :: eps, s, t, x, xmax, w
INTEGER   :: j, k, l, max, n
!----------------------------

!     ****** XMAX, MAX, AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            XMAX IS THE LARGEST POSITIVE REAL NUMBER THAT MAY
!            BE USED, MAX IS THE LARGEST POSITIVE INTEGER THAT
!            MAY BE USED, AND EPS IS THE SMALLEST REAL NUMBER
!            SUCH THAT 1._dp + EPS > 1._dp.

!                      MAX  = IPMPAR(3)
max = HUGE(3)
eps = EPSILON(1.0_dp)
xmax = HUGE(1.0_dp)

!----------------------------
fn_val = 0._dp
x = a
IF (ABS(a) < 10._dp) THEN
!-----------------------------------------------------------------------
!             EVALUATION OF DPSI(A) FOR ABS(A) < 10
!-----------------------------------------------------------------------
  t = 0._dp
  n = x
  n = n - 1

!     LET T BE THE SUM OF 1/(A-J) WHEN A >= 2

  IF (n < 0) GO TO 40
  DO j = 1, n
    x = x - 1._dp
    t = 1._dp / x + t
  END DO
  x = x - 1._dp
  GO TO 60

!     CHECK IF 1/A CAN OVERFLOW

  40 IF (ABS(a) < 1.d-35) THEN
    IF (ABS(a)*xmax <= 1.000000001_dp) RETURN
  END IF

!     LET T BE THE SUM OF -1/(A+J) WHEN A < 1

  t = -1._dp / a
  IF (a <= 0._dp) THEN
    n = -n - 1
    IF (n /= 0) THEN
      DO j = 1, n
        x = x + 1._dp
        IF (x == 0._dp) RETURN
        t = t - 1._dp / x
      END DO
    END IF
    x = (x+0.5_dp) + 0.5_dp
    IF (x == 0._dp) RETURN
    t = t - 1._dp / x
  END IF

!     COMPUTE  T + DPSI(1 + X)  FOR 0 <= X < 1

  60 IF (ABS(x-x0) <= 2.d-2) THEN
    fn_val = t + dpsi0(1._dp+x)
    RETURN
  END IF
  k = 42
  IF (eps > 1.d-20) k = 28
  fn_val = t + dcsevl(2._dp*x-1._dp,c,k)
  RETURN
END IF
!-----------------------------------------------------------------------
!           EVALUATION OF DPSI(A) FOR ABS(A) >= 10
!-----------------------------------------------------------------------
IF (a <= 0._dp) THEN
  t = MAX
  IF (ABS(a) >= MIN(t,1._dp/eps)) RETURN

!     SET W = PI*COT(PI*A) WHEN A IS NEGATIVE

  k = ABS(a)
  t = a + k
  IF (t == 0._dp) RETURN
  IF (t <= -0.5_dp) t = 1._dp + t
  t = pi * t
  w = pi * (COS(t)/SIN(t))
  x = 1._dp - x
END IF

!     COMPUTE THE MODIFIED ASYMPTOTIC SUM

t = (10._dp/x) ** 2
s = p(15)
DO j = 1, 14
  l = 15 - j
  s = p(l) + t * s
END DO
s = 0.5_dp / x + t * s

!     FINAL ASSEMBLY

fn_val = LOG(x) - s
IF (a < 0._dp) fn_val = fn_val - w
RETURN
END FUNCTION dpsi


FUNCTION dpsi0(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!            TAYLOR SERIES EXPANSION OF PSI(X) AROUND X0,
!                  WHERE X0 IS THE ZERO OF PSI(X).

!-------------------------
!     WRITTEN BY A.H. MORRIS
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: a(20) = (/ .967672245447621170427444761710D+00,  &
        -.442763168983592106092865281853D+00,  &
        .258499760955651010624401385701D+00,  &
        -.163942705442406527504251292747D+00,  &
        .107824050691262365757182948867D+00,  &
        -.721995612564547109261217836051D-01,  &
        .488042881641431072250925255079D-01,  &
        -.331611264748473592922583984045D-01,  &
        .225976482322181046596248251178D-01,  &
        -.154247659049489591388003168412D-01,  &
        .105387916166121753881240498824D-01,  &
        -.720453438635686824097047437040D-02,  &
        .492678139572985344635426640268D-02,  &
        -.336980165543932808279285672353D-02,  &
        .230512632673492783693838028298D-02,  &
        -.157693677143019725927093497173D-02,  &
        .107882520191629658069191777474D-02,  &
        -.738070938996005129566047389379D-03,  &
        .504953265834602035177398177463D-03,  &
        -.345468025106307699555567970882D-03 /)
REAL (dp), PARAMETER :: dk1 = 100442596182._dp, dk2 = 51069247913._dp,  &
        dk3 = 53827985572._dp, db = 68719476736._dp,   &
        dx = .28939299282041499433886199389507989269636D-32
REAL (dp) :: db2, h, w
INTEGER   :: i, l, n, nm1
!-------------------------

!     SET  H = X - X0  WHERE X0 IS THE ZERO OF PSI(X). X0 HAS THE
!     APPROXIMATE 60 DIGIT VALUE ...

!      1.4616321449683623412 62659542325721328468 19620400644635129598

!     A MORE ACCURATE VALUE IS GIVEN BY ...

!            X0 = DK1/8**12 + DK2/8**24 + DK3/8**36 + DX

!     THE FOLLOWING CODE SHOULD YIELD THE CORRECT VALUE FOR H IF A
!     BINARY, OCTAL, OR HEXADECIMAL REAL (dp) ARITHMETIC IS
!     BEING USED.

db2 = db * db
h = (((x-dk1/db)-dk2/db2)-dk3/(db*db2)) - dx

!-------------------------

n = 20
nm1 = n - 1
w = a(n)
DO i = 1, nm1
  l = n - i
  w = a(l) + h * w
END DO
fn_val = h * w
RETURN
END FUNCTION dpsi0


SUBROUTINE psidf(x, n, m, ans, iflag)
!-----------------------------------------------------------------------

!         PSIDF COMPUTES M MEMBER SEQUENCES OF SCALED DERIVATIVES OF
!         THE PSI FUNCTION

!                W(K,X)=(-1)**(K+1)*PSI(K,X)/GAMMA(K+1)

!         K=N,...,N+M-1 WHERE PSI(K,X) IS THE K-TH DERIVATIVE OF THE
!         PSI FUNCTION.

!         THE BASIC METHOD OF EVALUATION IS THE ASYMPTOTIC EXPANSION
!         FOR LARGE X>=XMIN FOLLOWED BY BACKWARD RECURSION ON A TWO
!         TERM RECURSION RELATION

!                  W(X+1) + X**(-N-1) = W(X).

!         THIS IS SUPPLEMENTED BY A SERIES

!                  SUM( (X+K)**(-N-1) , K=0,1,2,... )

!         WHICH CONVERGES RAPIDLY FOR LARGE N. BOTH XMIN AND THE
!         NUMBER OF TERMS OF THE SERIES ARE CALCULATED FROM THE UNIT
!         ROUND OFF OF THE MACHINE ENVIRONMENT.

!         THE NOMINAL COMPUTATIONAL ACCURACY IS THE MAXIMUM OF UNIT
!         ROUNDOFF (=EPSILON(1.0)) AND 1.0E-18 SINCE CRITICAL CONSTANTS
!         ARE GIVEN TO ONLY 18 DIGITS.

!     DESCRIPTION OF ARGUMENTS

!         INPUT

!           X      - ARGUMENT, X > 0.0

!           N      - FIRST MEMBER OF THE SEQUENCE, N >= 0

!           M      - NUMBER OF MEMBERS OF THE SEQUENCE, M >= 1

!         OUTPUT

!           ANS    - A VECTOR OF LENGTH AT LEAST M WHOSE FIRST M
!                    COMPONENTS ARE THE SCALED DERIVATIVES.

!           IFLAG  - A VARIABLE WHICH REPORTS THE STATUS OF THE
!                    RESULTS.
!                    IFLAG = 0 THE DESIRED VALUES WERE OBTAINED.
!                    IFLAG = 1 AN INPUT ERROR WAS DETECTED.
!                    IFLAG = 2 OVERFLOW. X TOO SMALL OR N+M-1
!                              TOO LARGE.
!                    IFLAG = 3 UNDERFLOW. X TOO LARGE OR N+M-1
!                              TOO LARGE.
!                    IFLAG = 4 N+M-1 IS TOO LARGE FOR THE CURRENT
!                              VALUE OF X. THIS SETTING WILL NOT
!                              OCCUR WHEN N+M-1 <= 100.

!-----------------------------------------------------------------------
!     WRITTEN BY DONALD E. AMOS
!         SANDIA LABORATORIES
!         JUNE 1982
!     MODIFIED BY A. H. MORRIS (NSWC), 1990.

!     REFERENCES ...

!     (1) ACM TRANS. MATH SOFTWARE, 1983.
!     (2) HANDBOOK OF MATHEMATICAL FUNCTIONS, AMS 55, NATIONAL BUREAU
!         OF STANDARDS BY M. ABRAMOWITZ AND I.A. STEGUN, 1964, PP.
!         258-260, EQUATIONS 6.3.5, 6.3.18, 6.4.6, 6.4.9, 6.4.10.
!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: x
INTEGER, INTENT(IN)  :: n, m
REAL, INTENT(OUT)    :: ans(:)
INTEGER, INTENT(OUT) :: iflag

! Local variables
REAL :: alpha, arg, beta, c = .43429, den, elim, eps, fln, fn, fnp, fns,  &
        fx, nd, rxsq, s, t, ta, tk, tol, tols, tss, tst, tt, t1, t2,  &
        wdtol, xdmln, xdmy, xinc, xln, xm, xmin
REAL, PARAMETER :: b(22) = (/ 1.0, -5.00000000000000000E-01,  &
        1.66666666666666667E-01, -3.33333333333333333E-02,  &
        2.38095238095238095E-02, -3.33333333333333333E-02,  &
        7.57575757575757576E-02, -2.53113553113553114E-01,  &
        1.16666666666666667E+00, -7.09215686274509804E+00,  &
        5.49711779448621554E+01, -5.29124242424242424E+02,  &
        6.19212318840579710E+03, -8.65802531135531136E+04,  &
        1.42551716666666667E+06, -2.72982310678160920E+07,  &
        6.01580873900642368E+08, -1.51163157670921569E+10,  &
        4.29614643061166667E+11, -1.37116552050883328E+13,  &
        4.88332318973593167E+14, -1.92965793419400681E+16 /)
REAL    :: trm(22), trmr(100)
INTEGER :: i, j, k, mx, nn, np, nx
INTEGER, PARAMETER :: nmax = 100
!------------------------
!     C = 1/LN(10)
!------------------------
iflag = 0
IF (x > 0.0 .AND. n >= 0 .AND. m >= 1) THEN

  nn = n + m - 1
  fn = REAL(nn)
  fnp = fn + 1.0
  eps = EPSILON(1.0)
  wdtol = MAX(0.5*eps,0.5E-18)
!-----------------------------------------------------------------------
!     ELIM = APPROXIMATE EXPONENTIAL OVER AND UNDERFLOW LIMIT
!-----------------------------------------------------------------------
  elim = MIN(exparg(0),ABS(exparg(1))) - 6.906
  xln = LOG(x)
  t = fnp * xln
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST FOR SMALL AND LARGE X
!-----------------------------------------------------------------------
  IF (ABS(t) > elim) GO TO 170
  IF (x >= wdtol) THEN
!-----------------------------------------------------------------------
!     COMPUTE XMIN AND THE NUMBER OF TERMS OF THE SERIES, FLN+1
!-----------------------------------------------------------------------
    nd = -c * epsln()
    nd = MIN(nd,18.0)
    fln = nd - 3.0
    alpha = 3.5 + 0.4 * fln
    beta = 0.21 + fln * (0.0006038*fln+0.008677)
    xm = alpha + beta * fn
    mx = INT(xm) + 1
    xmin = REAL(mx)

    IF (n /= 0) THEN
      xm = -2.302 * nd - MIN(0.0,xln)
      fns = REAL(n)
      arg = xm / fns
      arg = MIN(0.0,arg)
      eps = EXP(arg)
      xm = 1.0 - eps
      IF (ABS(arg) < 1.0E-3) xm = -arg
      fln = x * xm / eps
      xm = xmin - x
      IF (xm > 7.0 .AND. fln < 15.0) GO TO 110
    END IF

    xdmy = x
    xdmln = xln
    xinc = 0.0
    IF (x < xmin) THEN
      nx = INT(x)
      xinc = xmin - REAL(nx)
      xdmy = x + xinc
      xdmln = LOG(xdmy)
    END IF
!-----------------------------------------------------------------------
!     GENERATE W(N+M-1,X) BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
    t = fn * xdmln
    t1 = xdmln + xdmln
    t2 = t + xdmln
    tk = MAX(ABS(t),ABS(t1),ABS(t2))
    IF (tk > elim) GO TO 180

    tss = EXP(-t)
    tt = 0.5 / xdmy
    t1 = tt
    tst = wdtol * tt
    IF (nn /= 0) t1 = tt + 1.0 / fn
    rxsq = 1.0 / (xdmy*xdmy)
    ta = 0.5 * rxsq
    t = fnp * ta
    s = t * b(3)
    IF (ABS(s) >= tst) THEN

      tk = 2.0
      DO k = 4, 22
        t = t * ((tk+fn+1.0)/(tk+1.0)) * ((tk+fn)/(tk+2.0)) * rxsq
        trm(k) = t * b(k)
        IF (ABS(trm(k)) < tst) GO TO 20
        s = s + trm(k)
        tk = tk + 2.0
      END DO
    END IF

    20   s = (s+t1) * tss
    IF (xinc /= 0.0) THEN
!-----------------------------------------------------------------------
!     BACKWARD RECUR FROM XDMY TO X
!-----------------------------------------------------------------------
      nx = INT(xinc)
      np = nn + 1
      IF (nx > nmax) GO TO 190
      IF (nn == 0) GO TO 80
      xm = xinc - 1.0
      fx = x + xm
!-----------------------------------------------------------------------
!     THIS LOOP SHOULD NOT BE CHANGED. FX IS ACCURATE WHEN X IS SMALL
!-----------------------------------------------------------------------
      DO i = 1, nx
        trmr(i) = fx ** (-np)
        s = s + trmr(i)
        xm = xm - 1.0
        fx = x + xm
      END DO
    END IF

    ans(m) = s
    IF (fn == 0.0) GO TO 100
!-----------------------------------------------------------------------
!     GENERATE LOWER DERIVATIVES, J<N+M-1
!-----------------------------------------------------------------------
    IF (m == 1) RETURN
    DO j = 2, m
      fnp = fn
      fn = fn - 1.0
      tss = tss * xdmy
      t1 = tt
      IF (fn /= 0.0) t1 = tt + 1.0 / fn
      t = fnp * ta
      s = t * b(3)
      IF (ABS(s) >= tst) THEN

        tk = 3.0 + fnp
        DO k = 4, 22
          trm(k) = trm(k) * fnp / tk
          IF (ABS(trm(k)) < tst) GO TO 50
          s = s + trm(k)
          tk = tk + 2.0
        END DO
      END IF

      50   s = (s+t1) * tss
      IF (xinc /= 0.0) THEN
        IF (fn == 0.0) GO TO 80
        xm = xinc - 1.0
        fx = x + xm
        DO i = 1, nx
          trmr(i) = trmr(i) * fx
          s = s + trmr(i)
          xm = xm - 1.0
          fx = x + xm
        END DO
      END IF

      mx = m - j + 1
      ans(mx) = s
      IF (fn == 0.0) GO TO 100
    END DO
    RETURN
!-----------------------------------------------------------------------
!     RECURSION FOR N = 0
!-----------------------------------------------------------------------
    80 DO i = 1, nx
      s = s + 1.0 / (x+REAL(nx-i))
    END DO

    100 ans(1) = s - xdmln
    RETURN
!-----------------------------------------------------------------------
!     COMPUTE BY SERIES (X+K)**(-(N+1)) , K=0,1,2,...
!-----------------------------------------------------------------------
    110 nn = INT(fln) + 1
    np = n + 1
    t1 = (fns+1.0) * xln
    t = EXP(-t1)
    s = t
    den = x
    DO i = 1, nn
      den = den + 1.0
      trm(i) = den ** (-np)
      s = s + trm(i)
    END DO
    ans(1) = s
    IF (m == 1) RETURN
!-----------------------------------------------------------------------
!     GENERATE HIGHER DERIVATIVES, J > N
!-----------------------------------------------------------------------
    tol = wdtol / 5.0
    DO j = 2, m
      t = t / x
      s = t
      tols = t * tol
      den = x
      DO i = 1, nn
        den = den + 1.0
        trm(i) = trm(i) / den
        s = s + trm(i)
        IF (trm(i) < tols) GO TO 140
      END DO
      140   ans(j) = s
    END DO
    RETURN
  END IF
!-----------------------------------------------------------------------
!     SMALL X < UNIT ROUND OFF
!-----------------------------------------------------------------------
  ans(1) = x ** (-n-1)
  IF (m == 1) RETURN
  k = 1
  DO i = 2, m
    ans(k+1) = ans(k) / x
    k = k + 1
  END DO
  RETURN
END IF
!-----------------------------------------------------------------------
!     ERROR RETURN
!-----------------------------------------------------------------------
iflag = 1
RETURN

170 IF (t <= 0.0) THEN
  iflag = 2
  RETURN
END IF

180 iflag = 3
RETURN
!                            INCREASE THE DIMENSION OF TRMR(NMAX)
190 iflag = 4
RETURN
END SUBROUTINE psidf


FUNCTION betaln(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
REAL, INTENT(IN) :: a0, b0
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: e = .918938533204673
REAL             :: a, b, c, h, u, v, w, z
INTEGER          :: i, n
!--------------------------
a = MIN(a0,b0)
b = MAX(a0,b0)
IF (a < 8.0) THEN
  IF (a < 1.0) THEN
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A < 1
!-----------------------------------------------------------------------
    IF (b < 8.0) THEN
      fn_val = gamln(a) + (gamln(b)-gamln(a+b))
      RETURN
    END IF
    fn_val = gamln(a) + algdiv(a,b)
    RETURN
  END IF
!-----------------------------------------------------------------------
!                PROCEDURE WHEN 1 <= A < 8
!-----------------------------------------------------------------------
  IF (a <= 2.0) THEN
    IF (b <= 2.0) THEN
      fn_val = gamln(a) + gamln(b) - gsumln(a,b)
      RETURN
    END IF
    w = 0.0
    IF (b < 8.0) GO TO 20
    fn_val = gamln(a) + algdiv(a,b)
    RETURN
  END IF

!                REDUCTION OF A WHEN B <= 1000

  IF (b > 1000.0) GO TO 40
  n = a - 1.0
  w = 1.0
  DO i = 1, n
    a = a - 1.0
    h = a / b
    w = w * (h/(1.0+h))
  END DO
  w = LOG(w)
  IF (b >= 8.0) THEN
    fn_val = w + gamln(a) + algdiv(a,b)
    RETURN
  END IF

!                 REDUCTION OF B WHEN B < 8

  20 n = b - 1.0
  z = 1.0
  DO i = 1, n
    b = b - 1.0
    z = z * (b/(a+b))
  END DO
  fn_val = w + LOG(z) + (gamln(a)+(gamln(b)-gsumln(a,b)))
  RETURN

!                REDUCTION OF A WHEN B > 1000

  40 n = a - 1.0
  w = 1.0
  DO i = 1, n
    a = a - 1.0
    w = w * (a/(1.0+a/b))
  END DO
  fn_val = (LOG(w)-n*LOG(b)) + (gamln(a)+algdiv(a,b))
  RETURN
END IF
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A >= 8
!-----------------------------------------------------------------------
w = bcorr(a,b)
h = a / b
c = h / (1.0+h)
u = -(a-0.5) * LOG(c)
v = b * alnrel(h)
IF (u > v) THEN
  fn_val = (((-0.5*LOG(b)+e)+w)-v) - u
  RETURN
END IF
fn_val = (((-0.5*LOG(b)+e)+w)-u) - v
RETURN
END FUNCTION betaln


FUNCTION gsumln(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
!          FOR 1 <= A <= 2  AND  1 <= B <= 2
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a, b
REAL             :: fn_val

! Local variables
REAL :: x

x = DBLE(a) + DBLE(b) - 2._dp
IF (x <= 0.25) THEN
  fn_val = gamln1(1.0+x)
  RETURN
END IF
IF (x <= 1.25) THEN
  fn_val = gamln1(x) + alnrel(x)
  RETURN
END IF
fn_val = gamln1(x-1.0) + LOG(x*(1.0+x))
RETURN
END FUNCTION gsumln


FUNCTION bcorr(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------

!     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
!     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
!     IT IS ASSUMED THAT A0 >= 8 AND B0 >= 8.

!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a0, b0
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: c0 = .833333333333333E-01, c1 = -.277777777760991E-02,  &
                   c2 = .793650666825390E-03, c3 = -.595202931351870E-03,  &
                   c4 = .837308034031215E-03, c5 = -.165322962780713E-02
REAL :: a, b, c, h, s3, s5, s7, s9, s11, t, w, x, x2
!------------------------
a = MIN(a0,b0)
b = MAX(a0,b0)

h = a / b
c = h / (1.0+h)
x = 1.0 / (1.0+h)
x2 = x * x

!                SET SN = (1 - X**N)/(1 - X)

s3 = 1.0 + (x+x2)
s5 = 1.0 + (x+x2*s3)
s7 = 1.0 + (x+x2*s5)
s9 = 1.0 + (x+x2*s7)
s11 = 1.0 + (x+x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (1.0/b) ** 2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3) * t + c0
w = w * (c/b)

!                   COMPUTE  DEL(A) + W

t = (1.0/a) ** 2
fn_val = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0) / a + w
RETURN
END FUNCTION bcorr


FUNCTION algdiv(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8

!                         --------

!     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).

!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a, b
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: c0 = .833333333333333E-01, c1 = -.277777777760991E-02,  &
                    c2 = .793650666825390E-03, c3 = -.595202931351870E-03,  &
                    c4 = .837308034031215E-03, c5 = -.165322962780713E-02
REAL :: c, d, h, s3, s5, s7, s9, s11, t, u, v, w, x, x2
!------------------------
IF (a > b) THEN
  h = b / a
  c = 1.0 / (1.0+h)
  x = h / (1.0+h)
  d = a + (b-0.5)
ELSE
  h = a / b
  c = h / (1.0+h)
  x = 1.0 / (1.0+h)
  d = b + (a-0.5)
END IF

!                SET SN = (1 - X**N)/(1 - X)

x2 = x * x
s3 = 1.0 + (x+x2)
s5 = 1.0 + (x+x2*s3)
s7 = 1.0 + (x+x2*s5)
s9 = 1.0 + (x+x2*s7)
s11 = 1.0 + (x+x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (1.0/b) ** 2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3) * t + c0
w = w * (c/b)

!                    COMBINE THE RESULTS

u = d * alnrel(a/b)
v = a * (LOG(b)-1.0)
IF (u > v) THEN
  fn_val = (w-v) - u
  RETURN
END IF
fn_val = (w-u) - v
RETURN
END FUNCTION algdiv


FUNCTION dbetln(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a0, b0
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: e = .9189385332046727417803297364056_dp
REAL (dp)            :: a, b, c, h, sn, u, v, w, z
INTEGER              :: i, n
!--------------------------
!     E = 0.5*LN(2*PI)
!--------------------------
a = MIN(a0,b0)
b = MAX(a0,b0)
IF (a < 10._dp) THEN
  IF (a < 1._dp) THEN
!-----------------------------------------------------------------------
!                   PROCEDURE WHEN A < 1
!-----------------------------------------------------------------------
    IF (b < 10._dp) THEN
      fn_val = dgamln(a) + (dgamln(b)-dgamln(a+b))
      RETURN
    END IF
    fn_val = dgamln(a) + dlgdiv(a,b)
    RETURN
  END IF
!-----------------------------------------------------------------------
!               PROCEDURE WHEN 1 <= A < 10
!-----------------------------------------------------------------------
  IF (a <= 2._dp) THEN
    IF (b <= 2._dp) THEN
      fn_val = dgamln(a) + dgamln(b) - dgsmln(a,b)
      RETURN
    END IF
    w = 0._dp
    IF (b < 10._dp) GO TO 20
    fn_val = dgamln(a) + dlgdiv(a,b)
    RETURN
  END IF

!               REDUCTION OF A WHEN B <= 1000

  IF (b > 1.d3) GO TO 40
  n = a - 1._dp
  w = 1._dp
  DO i = 1, n
    a = a - 1._dp
    h = a / b
    w = w * (h/(1._dp+h))
  END DO
  w = LOG(w)
  IF (b >= 10._dp) THEN
    fn_val = w + dgamln(a) + dlgdiv(a,b)
    RETURN
  END IF

!                REDUCTION OF B WHEN B < 10

  20 n = b - 1._dp
  z = 1._dp
  DO i = 1, n
    b = b - 1._dp
    z = z * (b/(a+b))
  END DO
  fn_val = w + LOG(z) + (dgamln(a)+(dgamln(b)-dgsmln(a,b)))
  RETURN

!               REDUCTION OF A WHEN B > 1000

  40 n = a - 1._dp
  w = 1._dp
  DO i = 1, n
    a = a - 1._dp
    w = w * (a/(1._dp+a/b))
  END DO
  sn = n
  fn_val = (LOG(w)-sn*LOG(b)) + (dgamln(a)+dlgdiv(a,b))
  RETURN
END IF
!-----------------------------------------------------------------------
!                  PROCEDURE WHEN A >= 10
!-----------------------------------------------------------------------
w = dbcorr(a,b)
h = a / b
c = h / (1._dp+h)
u = -(a-0.5_dp) * LOG(c)
v = b * dlnrel(h)
IF (u > v) THEN
  fn_val = (((-0.5_dp*LOG(b)+e)+w)-v) - u
  RETURN
END IF
fn_val = (((-0.5_dp*LOG(b)+e)+w)-u) - v
RETURN
END FUNCTION dbetln


FUNCTION dgsmln(a, b) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
!          FOR 1 <= A <= 2  AND  1 <= B <= 2
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a, b
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: x

x = (a-1._dp) + (b-1._dp)
IF (x <= 0.5_dp) THEN
  fn_val = dgmln1(1._dp+x)
  RETURN
END IF
IF (x <= 1.5_dp) THEN
  fn_val = dgmln1(x) + dlnrel(x)
  RETURN
END IF
fn_val = dgmln1(x-1._dp) + LOG(x*(1._dp+x))
RETURN
END FUNCTION dgsmln


FUNCTION dbcorr(a0, b0) RESULT(fn_val)
!-----------------------------------------------------------------------

!     EVALUATION OF DEL(A) + DEL(B0) - DEL(A) + B0) WHERE
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
!     IT IS ASSUMED THAT A0 >= 10 AND B0 >= 10.

!                         --------

!     THE SERIES FOR DEL(X), WHICH APPLIES FOR X >= 10, WAS
!     DERIVED BY A.H. MORRIS FROM THE CHEBYSHEV SERIES IN THE
!     SLATEC LIBRARY OBTAINED BY WAYNE FULLERTON (LOS ALAMOS).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a0, b0
REAL (dp)             :: fn_val

! Local variables
REAL (dp)            :: a, b, c, h, s(15), t, w, x, x2, z
INTEGER              :: j, k
REAL (dp), PARAMETER :: e(15) = (/ .833333333333333333333333333333D-01,  &
   -.277777777777777777777777752282D-04, .793650793650793650791732130419D-07,  &
   -.595238095238095232389839236182D-09, .841750841750832853294451671990D-11,  &
   -.191752691751854612334149171243D-12, .641025640510325475730918472625D-14,  &
   -.295506514125338232839867823991D-15, .179643716359402238723287696452D-16,  &
   -.139228964661627791231203060395D-17, .133802855014020915603275339093D-18,  &
   -.154246009867966094273710216533D-19, .197701992980957427278370133333D-20,  &
   -.234065664793997056856992426667D-21, .171348014966398575409015466667D-22 /)
!--------------------------
a = MIN(a0,b0)
b = MAX(a0,b0)

h = a / b
c = h / (1._dp+h)
x = 1._dp / (1._dp+h)
x2 = x * x

!        COMPUTE (1 - X**N)/(1 - X) FOR N = 1,3,5,...
!            STORE THESE VALUES IN S(1),S(2),...

s(1) = 1._dp
DO j = 1, 14
  s(j+1) = 1._dp + (x+x2*s(j))
END DO

!                SET W = DEL(B) - DEL(A + B)

t = (10._dp/b) ** 2
w = e(15) * s(15)
DO j = 1, 14
  k = 15 - j
  w = t * w + e(k) * s(k)
END DO
w = w * (c/b)

!                    COMPUTE  DEL(A) + W

t = (10._dp/a) ** 2
z = e(15)
DO j = 1, 14
  k = 15 - j
  z = t * z + e(k)
END DO
fn_val = z / a + w
RETURN
END FUNCTION dbcorr


FUNCTION dlgdiv (a, b) RESULT(fn_val)
!-----------------------------------------------------------------------

!     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B >= 8

!                         --------

!     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
!     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a, b
REAL (dp)             :: fn_val

! EXTERNAL   alnrel
REAL (dp) :: c, d, h, s11, s3, s5, s7, s9, t, u, v, w, x, x2
REAL (dp), PARAMETER :: c0 = .833333333333333D-01,  &
             one = 1._dp, half = 0.5_dp, c1 = -.277777777760991D-02, &
             c2 = .793650666825390D-03, c3 = -.595202931351870D-03, &
             c4 = .837308034031215D-03, c5 = -.165322962780713D-02
!------------------------
IF (a > b) THEN
  h = b/a
  c = one/(one + h)
  x = h/(one + h)
  d = a + (b - half)
ELSE
  h = a/b
  c = h/(one + h)
  x = one/(one + h)
  d = b + (a - half)
END IF

!                SET SN = (1 - X**N)/(1 - X)

x2 = x*x
s3 = one + (x + x2)
s5 = one + (x + x2*s3)
s7 = one + (x + x2*s5)
s9 = one + (x + x2*s7)
s11 = one + (x + x2*s9)

!                SET W = DEL(B) - DEL(A + B)

t = (one/b)**2
w = ((((c5*s11*t + c4*s9)*t + c3*s7)*t + c2*s5)*t + c1*s3)*t + c0
w = w*(c/b)

!                    COMBINE THE RESULTS

u = d*dlnrel(a/b)
v = a*(LOG(b) - one)
IF (u > v) THEN
  fn_val = (w - v) - u
ELSE
  fn_val = (w - u) - v
END IF

RETURN
END FUNCTION dlgdiv



SUBROUTINE gratio(a, x, ans, qans, ind)
!-----------------------------------------------------------------------

!        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
!                      P(A,X) AND Q(A,X)

!                        ----------

!     IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X
!     ARE NOT BOTH 0.

!     ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE
!     P(A,X) AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER.
!     IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS
!     POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE, IF
!     IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE
!     6-TH SIGNIFICANT DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY
!     IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT.

!     ERROR RETURN ...

!        ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE,
!     WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT.
!     P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN
!     X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!     REVISED ... DEC 1991
!-------------------------
REAL, INTENT(IN)    :: a, x
INTEGER, INTENT(IN) :: ind
REAL, INTENT(OUT)   :: ans, qans

! Local variables
REAL, PARAMETER :: acc0(3) = (/ 5.e-15, 5.e-7, 5.e-4 /),  &
        big(3) = (/ 25.0, 14.0 , 10.0 /), e0(3) = (/ .25E-3, .25E-1, .14 /),  &
        x0(3) = (/ 31.0, 17.0 , 9.7 /), ALOG10 = 2.30258509299405,  &
        rt2pin = .398942280401433, rtpi = 1.77245385090552,  &
        d00 = -.333333333333333E+00, d10 = -.185185185185185E-02,   &
        d20 =  .413359788359788E-02, d30 = .649434156378601E-03,   &
        d40 = -.861888290916712E-03, d50 = -.336798553366358E-03,   &
        d60 =  .531307936463992E-03, d70 = .344367606892378E-03,   &
        d80 = -.652623918595309E-03
REAL, PARAMETER :: a0(4) = (/ -.231272501940775E-02, -.335378520024220E-01,  &
                           -.159840143443990E+00, -.333333333333333E+00 /),  &
        a1(4) = (/ -.398783924370770E-05, -.587926036018402E-03,  &
                   -.491687131726920E-02, -.185185185184291E-02 /),  &
        a2(2) = (/ .669564126155663E-03, .413359788442192E-02 /),  &
        a3(2) = (/ .810586158563431E-03, .649434157619770E-03 /),  &
        a4(2) = (/ -.105014537920131E-03, -.861888301199388E-03 /),  &
        a5(2) = (/ -.435211415445014E-03, -.336806989710598E-03 /),  &
        a6(2) = (/ -.182503596367782E-03, .531279816209452E-03 /),  &
        a7(2) = (/ .443219646726422E-03, .344430064306926E-03 /),  &
        a8(2) = (/ .878371203603888E-03, -.686013280418038E-03 /)
REAL, PARAMETER :: b0(6) = (/ .633763414209504E-06, -.939001940478355E-05,  &
        .239521354917408E-02, .376245718289389E-01, .238549219145773E+00,  &
        .729520430331981E+00 /),  &
        b1(4) = (/ .386325038602125E-02,  &
        .506042559238939E-01, .283344278023803E+00, .780110511677243E+00 /),  &
        b2(5) = (/ -.421924263980656E-03, .650837693041777E-02,  &
        .682034997401259E-01, .339173452092224E+00, .810647620703045E+00 /),  &
        b3(5) = (/ -.632276587352120E-03, .905375887385478E-02,  &
        .906610359762969E-01, .406288930253881E+00, .894800593794972E+00 /),  &
        b4(4) = (/ .322609381345173E-01, .178295773562970E+00,  &
        .591353097931237E+00, .103151890792185E+01 /),  &
        b5(3) = (/ .178716720452422E+00, .600380376956324E+00,  &
        .108515217314415E+01 /),  &
        b6(2) = (/ .345608222411837E+00, .770341682526774E+00 /),  &
        b7(2) = (/ .821824741357866E+00, .115029088777769E+01 /)
REAL :: d0(6) = (/ .833333333333333E-01, -.148148148148148E-01,  &
        .115740740740741E-02, .352733686067019E-03, -.178755144032922E-03,  &
        .391926317852244E-04 /), d1(4) = (/ -.347222222222222E-02,  &
        .264550264550265E-02, -.990226337448560E-03, .205761316872428E-03 /),  &
        d2(2) = (/ -.268132716049383E-02, .771604938271605E-03 /),  &
        d3(2) = (/ .229472093621399E-03, -.469189494395256E-03 /),  &
        d4(1) = (/ .784039221720067E-03 /), d5(1) = (/ -.697281375836586E-04 /),  &
        d6(1) = (/ -.592166437353694E-03 /)
REAL    :: acc, amn, apn, a2n, a2nm1, b2n, b2nm1, c, c0, c1, c2, c3, c4, c5, &
           c6, c7, c8, e, g, h, j, l, r, rta, rtx, s, sum, t, tol, twoa, u,  &
           w, wk(20), y, z
INTEGER :: i, iop, m, n, nl1
!-------------------------

!     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + E > 1.0 .

e = EPSILON(1.0)

!-------------------------
IF (a < 0.0.OR.x < 0.0) GO TO 320
IF (a == 0.0 .AND. x == 0.0) GO TO 320
IF (a*x == 0.0) GO TO 310

iop = ind + 1
IF (iop /= 1 .AND. iop /= 2) iop = 3
acc = MAX(acc0(iop),e)

!            SELECT THE APPROPRIATE ALGORITHM

IF (a < 1.0) THEN
  IF (a == 0.5) GO TO 290
  IF (x < 1.1) GO TO 100
  r = rcomp(a,x)
  IF (r == 0.0) GO TO 280
  GO TO 160
END IF

IF (a < big(iop)) THEN
  IF (a > x.OR.x >= x0(iop)) GO TO 10
  twoa = a + a
  m = INT(twoa)
  IF (twoa /= REAL(m)) GO TO 10
  i = m / 2
  IF (a == REAL(i)) GO TO 130
  GO TO 140
END IF

l = x / a
IF (l == 0.0) GO TO 270
s = 0.5 + (0.5-l)
z = rlog(l)
IF (z >= 700.0/a) GO TO 300
y = a * z
rta = SQRT(a)
IF (ABS(s) <= e0(iop)/rta) GO TO 230
IF (ABS(s) <= 0.4) GO TO 180

10 r = rcomp(a,x)
IF (r == 0.0) GO TO 310
IF (x > MAX(a,ALOG10)) THEN
  IF (x < x0(iop)) GO TO 160
ELSE

!                 TAYLOR SERIES FOR P/R

  apn = a + 1.0
  t = x / apn
  wk(1) = t
  DO n = 2, 20
    apn = apn + 1.0
    t = t * (x/apn)
    IF (t <= 1.e-3) GO TO 30
    wk(n) = t
  END DO
  n = 20

  30 sum = t
  tol = 0.5 * acc
  40 apn = apn + 1.0
  t = t * (x/apn)
  sum = sum + t
  IF (t > tol) GO TO 40

  nl1 = n - 1
  DO m = 1, nl1
    n = n - 1
    sum = sum + wk(n)
  END DO
  ans = (r/a) * (1.0+sum)
  qans = 0.5 + (0.5-ans)
  RETURN
END IF

!                 ASYMPTOTIC EXPANSION

amn = a - 1.0
t = amn / x
wk(1) = t
DO n = 2, 20
  amn = amn - 1.0
  t = t * (amn/x)
  IF (ABS(t) <= 1.e-3) GO TO 70
  wk(n) = t
END DO
n = 20

70 sum = t
80 IF (ABS(t) >= acc) THEN
  amn = amn - 1.0
  t = t * (amn/x)
  sum = sum + t
  GO TO 80
END IF

nl1 = n - 1
DO m = 1, nl1
  n = n - 1
  sum = sum + wk(n)
END DO
qans = (r/x) * (1.0+sum)
ans = 0.5 + (0.5-qans)
RETURN

!             TAYLOR SERIES FOR P(A,X)/X**A

100 l = 3.0
c = x
sum = x / (a+3.0)
tol = 3.0 * acc / (a+1.0)
110 l = l + 1.0
c = -c * (x/l)
t = c / (a+l)
sum = sum + t
IF (ABS(t) > tol) GO TO 110
j = a * x * ((sum/6.0-0.5/(a+2.0))*x+1.0/(a+1.0))

z = a * LOG(x)
h = gam1(a)
g = 1.0 + h
IF (x >= 0.25) THEN
  IF (a < x/2.59) GO TO 120
ELSE
  IF (z > -.13394) GO TO 120
END IF

w = EXP(z)
ans = w * g * (0.5+(0.5-j))
qans = 0.5 + (0.5-ans)
RETURN

120 l = rexp(z)
w = 0.5 + (0.5+l)
qans = (w*j-l) * g - h
IF (qans < 0.0) GO TO 280
ans = 0.5 + (0.5-qans)
RETURN

!             FINITE SUMS FOR Q WHEN A >= 1
!                 AND 2*A IS AN INTEGER

130 sum = EXP(-x)
t = sum
n = 1
c = 0.0
GO TO 150

140 rtx = SQRT(x)
sum = erfc1(0,rtx)
t = EXP(-x) / (rtpi*rtx)
n = 0
c = -0.5

150 IF (n /= i) THEN
  n = n + 1
  c = c + 1.0
  t = (x*t) / c
  sum = sum + t
  GO TO 150
END IF
qans = sum
ans = 0.5 + (0.5-qans)
RETURN

!              CONTINUED FRACTION EXPANSION

160 tol = MAX(8.0*e,4.0*acc)
a2nm1 = 1.0
a2n = 1.0
b2nm1 = x
b2n = x + (1.0-a)
c = 1.0
170 a2nm1 = x * a2n + c * a2nm1
b2nm1 = x * b2n + c * b2nm1
c = c + 1.0
t = c - a
a2n = a2nm1 + t * a2n
b2n = b2nm1 + t * b2n

a2nm1 = a2nm1 / b2n
b2nm1 = b2nm1 / b2n
a2n = a2n / b2n
b2n = 1.0
IF (ABS(a2n-a2nm1/b2nm1) >= tol*a2n) GO TO 170

qans = r * a2n
ans = 0.5 + (0.5-qans)
RETURN

180 IF (ABS(s) <= 2.0*e .AND. a*e*e > 3.28E-3) GO TO 320
c = EXP(-y)
w = 0.5 * erfc1(1,SQRT(y))
u = 1.0 / a
z = SQRT(z+z)
IF (l < 1.0) z = -z
IF (iop == 2) THEN
  GO TO 200
ELSE IF (iop > 2) THEN
  GO TO 210
END IF

IF (ABS(s) <= 1.e-3) GO TO 240

!            USING THE MINIMAX APPROXIMATIONS

c0 = (((a0(1)*z + a0(2))*z + a0(3))*z + a0(4)) / ((((((b0(1)*z + b0(2))*z +  &
     b0(3))*z + b0(4))*z + b0(5))*z + b0(6))*z + 1.0)
c1 = (((a1(1)*z + a1(2))*z + a1(3))*z + a1(4)) / ((((b1(1)*z + b1(2))*z +   &
     b1(3))*z + b1(4))*z + 1.0)
c2 = (a2(1)*z + a2(2)) / (((((b2(1)*z + b2(2))*z + b2(3))*z + b2(4))*z +  &
     b2(5))*z + 1.0)
c3 = (a3(1)*z + a3(2)) / (((((b3(1)*z + b3(2))*z + b3(3))*z + b3(4))*z +  &
     b3(5))*z + 1.0)
c4 = (a4(1)*z + a4(2)) / ((((b4(1)*z + b4(2))*z + b4(3))*z + b4(4))*z + 1.0)
c5 = (a5(1)*z + a5(2)) / (((b5(1)*z + b5(2))*z + b5(3))*z + 1.0)
c6 = (a6(1)*z + a6(2)) / ((b6(1)*z + b6(2))*z + 1.0)
c7 = (a7(1)*z + a7(2)) / ((b7(1)*z + b7(2))*z + 1.0)
c8 = a8(1) * z + a8(2)
t = (((((((c8*u + c7)*u + c6)*u + c5)*u + c4)*u + c3)*u + c2)*u + c1)*u + c0
GO TO 220

!                    TEMME EXPANSION

200 c0 = (((((d0(6)*z+d0(5))*z+d0(4))*z+d0(3))*z+d0(2))*z+d0(1)) * z +d00
c1 = (((d1(4)*z+d1(3))*z+d1(2))*z+d1(1)) * z + d10
c2 = d2(1) * z + d20
t = (c2*u+c1) * u + c0
GO TO 220

210 t = ((d0(3)*z+d0(2))*z+d0(1)) * z + d00

220 IF (l >= 1.0) THEN
  qans = c * (w+rt2pin*t/rta)
  ans = 0.5 + (0.5-qans)
  RETURN
END IF
ans = c * (w-rt2pin*t/rta)
qans = 0.5 + (0.5-ans)
RETURN

!               TEMME EXPANSION FOR L = 1

230 IF (a*e*e > 3.28E-3) GO TO 320
c = 0.5 + (0.5-y)
w = (0.5 - SQRT(y)*(0.5 + (0.5 - y/3.0))/rtpi) / c
u = 1.0 / a
z = SQRT(z+z)
IF (l < 1.0) z = -z
IF (iop == 2) THEN
  GO TO 250
ELSE IF (iop > 2) THEN
  GO TO 260
END IF

240 c0 = ((d0(3)*z+d0(2))*z+d0(1)) * z + d00
c1 = ((d1(3)*z+d1(2))*z+d1(1)) * z + d10
c2 = (d2(2)*z+d2(1)) * z + d20
c3 = (d3(2)*z+d3(1)) * z + d30
c4 = d4(1) * z + d40
c5 = d5(1) * z + d50
c6 = d6(1) * z + d60
t = (((((((d80*u+d70)*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1) * u + c0
GO TO 220

250 c0 = (d0(2)*z+d0(1)) * z + d00
c1 = d1(1) * z + d10
t = (d20*u+c1) * u + c0
GO TO 220

260 t = d0(1) * z + d00
GO TO 220

!                     SPECIAL CASES

270 ans = 0.0
qans = 1.0
RETURN

280 ans = 1.0
qans = 0.0
RETURN

290 IF (x < 0.25) THEN
  ans = erf(SQRT(x))
  qans = 0.5 + (0.5-ans)
  RETURN
END IF
qans = erfc1(0,SQRT(x))
ans = 0.5 + (0.5-qans)
RETURN

300 IF (ABS(s) <= 2.0*e) GO TO 320
310 IF (x <= a) GO TO 270
GO TO 280

!                     ERROR RETURN

320 ans = 2.0
RETURN
END SUBROUTINE gratio



FUNCTION rcomp(a, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!                EVALUATION OF EXP(-X)*X**A/GAMMA(A)
!-----------------------------------------------------------------------
!     RT2PIN = 1/SQRT(2*PI)
!------------------------
REAL, INTENT(IN) :: a, x
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: rt2pin = .398942280401433
REAL             :: t, t1, u
!------------------------
fn_val = 0.0
IF (x == 0.0) RETURN
IF (a < 20.0) THEN

  t = a * LOG(x) - x
  IF (t < exparg(1)) RETURN
  IF (a < 1.0) THEN
    fn_val = (a*EXP(t)) * (1.0+gam1(a))
    RETURN
  END IF
  fn_val = EXP(t) / gamma(a)
  RETURN
END IF

u = x / a
IF (u == 0.0) RETURN
t = (1.0/a) ** 2
t1 = (((0.75*t-1.0)*t+3.5)*t-105.0) / (a*1260.0)
t1 = t1 - a * rlog(u)
IF (t1 >= exparg(1)) fn_val = rt2pin * SQRT(a) * EXP(t1)
RETURN
END FUNCTION rcomp


FUNCTION drcomp(a, x) RESULT(fn_val)
!-----------------------------------------------------------------------
!              EVALUATION OF EXP(-X)*X**A/GAMMA(A)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a, x
REAL (dp)             :: fn_val

! Local variables
REAL (dp), PARAMETER :: c = .398942280401432677939946059934_dp
REAL (dp)            :: t, w
!--------------------------
!     C = 1/SQRT(2*PI)
!--------------------------
fn_val = 0._dp
IF (x == 0._dp) RETURN
IF (a <= 20._dp) THEN

  t = a * LOG(x) - x
  IF (t < dxparg(1)) RETURN
  IF (a < 1._dp) THEN
    fn_val = (a*EXP(t)) * (1._dp+dgam1(a))
    RETURN
  END IF
  fn_val = EXP(t) / dgamma(a)
  RETURN
END IF

t = x / a
IF (t == 0._dp) RETURN
w = -(dpdel(a)+a*drlog(t))
IF (w >= dxparg(1)) fn_val = c * SQRT(a) * EXP(w)
RETURN
END FUNCTION drcomp


SUBROUTINE gaminv(a, x, x0, p, q, ierr)
!-----------------------------------------------------------------------

!             INVERSE INCOMPLETE GAMMA RATIO FUNCTION

!     GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1.
!     THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
!     ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X
!     TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE
!     PARTICULAR COMPUTER ARITHMETIC BEING USED.

!                        ------------

!     X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
!     AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
!     NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN
!     A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
!     IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.

!     X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER DOES NOT
!     WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET X0 <= 0.

!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!     WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
!     VALUES ...

!       IERR =  0    THE SOLUTION WAS OBTAINED. ITERATION WAS
!                    NOT USED.
!       IERR>0    THE SOLUTION WAS OBTAINED. IERR ITERATIONS
!                    WERE PERFORMED.
!       IERR = -2    (INPUT ERROR) A <= 0
!       IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A
!                    IS TOO LARGE.
!       IERR = -4    (INPUT ERROR) P OR Q IS NEGATIVE, OR
!                    P + Q .NE. 1.
!       IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST
!                    RECENT VALUE OBTAINED FOR X IS GIVEN.
!                    THIS CANNOT OCCUR IF X0 <= 0.
!       IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
!                    THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
!       IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE
!                    ROUTINE IS NOT CERTAIN OF ITS ACCURACY.
!                    ITERATION CANNOT BE PERFORMED IN THIS
!                    CASE. IF X0 <= 0, THIS CAN OCCUR ONLY
!                    WHEN P OR Q IS APPROXIMATELY 0. IF X0 IS
!                    POSITIVE THEN THIS CAN OCCUR WHEN A IS
!                    EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY
!                    LARGE (SAY A >= 1.E20).

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!     REVISED ... JANUARY 1992
!------------------------
REAL, INTENT(IN)     :: a, x0, p, q
REAL, INTENT(OUT)    :: x
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL, PARAMETER :: ln10 = 2.302585, bmin(2) = (/ 1.e-28, 1.e-13 /),  &
        emin(2) = (/ 2.e-03, 6.e-03 /), c = .577215664901533, tol = 1.e-5
REAL    :: amax, amin, am1, ap1, ap2, ap3, apn, b, c1, c2, c3, c4, c5, d,  &
           e, eps, e2, g, h, pn, qg, qn, r, rta, s, sum, s2, t, u, w, xmin, &
           xn, y, z
INTEGER :: ier, iop
!------------------------
!     LN10 = LN(10)
!     C = EULER CONSTANT
!------------------------

!     ****** E AND XMIN ARE MACHINE DEPENDENT CONSTANTS. E IS THE
!            SMALLEST NUMBER FOR WHICH 1.0 + E > 1.0, AND XMIN
!            IS THE SMALLEST POSITIVE NUMBER.

e = EPSILON(1.0)
xmin = TINY(1.0)

!------------------------
x = 0.0
IF (a > 0.0) THEN
  IF (p < 0.0.OR.q < 0.0) GO TO 120
  t = ((p+q)-0.5) - 0.5
  IF (ABS(t) > 5.0*MAX(e,1.e-15)) GO TO 120

  ierr = 0
  xmin = xmin / e
  IF ((p/e) > xmin) THEN
    IF ((q/e) <= xmin) GO TO 160
    IF (a == 1.0) GO TO 100

    e2 = e + e
    amax = 0.4E-10 / (e*e)
    eps = MAX(100.0*e,1.e-10)
    iop = 1
    IF (e > 1.e-10) iop = 2
    xn = x0
    IF (x0 <= 0.0) THEN

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A < 1

      IF (a > 1.0) GO TO 30
      g = gamma(a+1.0)
      qg = q * g
      IF (qg == 0.0) GO TO 160
      b = qg / a
      IF (qg > 0.6*a) GO TO 20
      IF (a < 0.30 .AND. b >= 0.35) THEN
        t = EXP(-(b+c))
        u = t * EXP(t)
        xn = t * EXP(u)
        GO TO 50
      END IF

      IF (b >= 0.45) GO TO 20
      IF (b == 0.0) GO TO 160
      y = -LOG(b)
      s = 0.5 + (0.5-a)
      z = LOG(y)
      t = y - s * z
      IF (b >= 0.15) THEN
        xn = y - s * LOG(t) - LOG(1.0+s/(t+1.0))
        GO TO 80
      END IF
      IF (b > 1.e-2) THEN
        u = ((t+2.0*(3.0-a))*t+(2.0-a)*(3.0-a)) / ((t+(5.0-a))*t+2.0)
        xn = y - s * LOG(t) - LOG(u)
        GO TO 80
      END IF
      10 c1 = -s * z
      c2 = -s * (1.0+c1)
      c3 = s * ((0.5*c1+(2.0-a))*c1+(2.5-1.5*a))
      c4 = -s * (((c1/3.0+(2.5-1.5*a))*c1+((a-6.0)*a+7.0))*c1 +  &
                ((11.0*a-46.0)*a+47.0)/6.0)
      c5 = -s * ((((-c1/4.0+(11.0*a-17.0)/6.0)*c1+((-3.0*a+13.0)*a  &
           -13.0))*c1+0.5*(((2.0*a-25.0)*a+72.0)*a-61.0))*c1 +  &
           (((25.0*a-195.0)*a+477.0)*a-379.0)/12.0)
      xn = ((((c5/y+c4)/y+c3)/y+c2)/y+c1) + y
      IF (a > 1.0) GO TO 80
      IF (b > bmin(iop)) GO TO 80
      x = xn
      RETURN

      20 IF (b*q <= 1.e-8) THEN
        xn = EXP(-(q/a+c))
      ELSE
        IF (p > 0.9) THEN
          xn = EXP((alnrel(-q)+gamln1(a))/a)
        ELSE
          xn = EXP(LOG(p*g)/a)
        END IF
      END IF

      IF (xn == 0.0) GO TO 110
      t = 0.5 + (0.5-xn/(a+1.0))
      xn = xn / t
      GO TO 50

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A > 1

      30 t = p - 0.5
      IF (q < 0.5) t = 0.5 - q
      CALL pni(p, q, t, s, ier)
      IF (ier /= 0) WRITE(*, *) '** Error in call to PNI from GAMINV **'

      rta = SQRT(a)
      s2 = s * s
      xn = (((12.0*s2 - 243.0)*s2 - 923.0)*s2 + 1472.0) / 204120.0
      xn = (xn/a + s*((9.0*s2 + 256.0)*s2 - 433.0)/(38880.0*rta)) -   &
           ((3.0*s2 + 7.0)*s2 - 16.0) / 810.0
      xn = a + s * rta + (s2-1.0) / 3.0 + s * (s2-7.0) / (36.0*rta) + xn / a
      xn = MAX(xn, 0.0)

      amin = 20.0
      IF (e < 1.e-8) amin = 250.0
      IF (a >= amin) THEN
        x = xn
        d = 0.5 + (0.5-x/a)
        IF (ABS(d) <= 1.e-1) RETURN
      END IF

      IF (p > 0.5) THEN
        IF (xn < 3.0*a) GO TO 80
        w = LOG(q)
        y = -(w+gamln(a))
        d = MAX(2.0,a*(a-1.0))
        IF (y >= ln10*d) THEN
          s = 1.0 - a
          z = LOG(y)
          GO TO 10
        END IF
        t = a - 1.0
        xn = y + t * LOG(xn) - alnrel(-t/(xn+1.0))
        xn = y + t * LOG(xn) - alnrel(-t/(xn+1.0))
        GO TO 80
      END IF

      ap1 = a + 1.0
      IF (xn > 0.70*ap1) GO TO 60
      w = LOG(p) + gamln(ap1)
      IF (xn <= 0.15*ap1) THEN
        ap2 = a + 2.0
        ap3 = a + 3.0
        x = EXP((w+x)/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+x/ap2)))/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+x/ap2)))/a)
        x = EXP((w+x-LOG(1.0+(x/ap1)*(1.0+(x/ap2)*(1.0+x/ap3))))/a)
        xn = x
        IF (xn <= 1.e-2*ap1) THEN
          IF (xn <= emin(iop)*ap1) RETURN
          GO TO 60
        END IF
      END IF

      apn = ap1
      t = xn / apn
      sum = 1.0 + t
      40 apn = apn + 1.0
      t = t * (xn/apn)
      sum = sum + t
      IF (t > 1.e-4) GO TO 40
      t = w - LOG(sum)
      xn = EXP((xn+t)/a)
      xn = xn * (1.0-(a*LOG(xn)-xn-t)/(a-xn))
      GO TO 60
    END IF

!                 SCHRODER ITERATION USING P

    50 IF (p > 0.5) GO TO 80
    60 IF (p <= xmin) GO TO 150
    am1 = (a-0.5) - 0.5
    70 IF (a > amax) THEN
      d = 0.5 + (0.5-xn/a)
      IF (ABS(d) <= e2) GO TO 150
    END IF

    IF (ierr >= 20) GO TO 130
    ierr = ierr + 1
    CALL gratio(a,xn,pn,qn,0)
    IF (pn == 0.0.OR.qn == 0.0) GO TO 150
    r = rcomp(a,xn)
    IF (r < xmin) GO TO 150
    t = (pn-p) / r
    w = 0.5 * (am1-xn)
    IF (ABS(t) > 0.1.OR.ABS(w*t) > 0.1) THEN
      x = xn * (1.0-t)
      IF (x <= 0.0) GO TO 140
      d = ABS(t)
    ELSE

      h = t * (1.0+w*t)
      x = xn * (1.0-h)
      IF (x <= 0.0) GO TO 140
      IF (ABS(w) >= 1.0 .AND. ABS(w)*t*t <= eps) RETURN
      d = ABS(h)
    END IF
    xn = x
    IF (d > tol) GO TO 70
    IF (d <= eps) RETURN
    IF (ABS(p-pn) <= tol*p) RETURN
    GO TO 70

!                 SCHRODER ITERATION USING Q

    80 IF (q <= xmin) GO TO 150
    am1 = (a-0.5) - 0.5
    90 IF (a > amax) THEN
      d = 0.5 + (0.5-xn/a)
      IF (ABS(d) <= e2) GO TO 150
    END IF

    IF (ierr >= 20) GO TO 130
    ierr = ierr + 1
    CALL gratio(a,xn,pn,qn,0)
    IF (pn == 0.0.OR.qn == 0.0) GO TO 150
    r = rcomp(a,xn)
    IF (r < xmin) GO TO 150
    t = (q-qn) / r
    w = 0.5 * (am1-xn)
    IF (ABS(t) > 0.1.OR.ABS(w*t) > 0.1) THEN
      x = xn * (1.0-t)
      IF (x <= 0.0) GO TO 140
      d = ABS(t)
    ELSE

      h = t * (1.0+w*t)
      x = xn * (1.0-h)
      IF (x <= 0.0) GO TO 140
      IF (ABS(w) >= 1.0 .AND. ABS(w)*t*t <= eps) RETURN
      d = ABS(h)
    END IF
    xn = x
    IF (d > tol) GO TO 90
    IF (d <= eps) RETURN
    IF (ABS(q-qn) <= tol*q) RETURN
    GO TO 90
  END IF

!                       SPECIAL CASES

  ierr = -8
  RETURN

  100 IF (q >= 0.9) THEN
    x = -alnrel(-p)
    RETURN
  END IF
  x = -LOG(q)
  RETURN
END IF

!                       ERROR RETURN

ierr = -2
RETURN

110 ierr = -3
RETURN

120 ierr = -4
RETURN

130 ierr = -6
RETURN

140 ierr = -7
RETURN

150 x = xn
ierr = -8
RETURN

160 x = HUGE(1.0)
ierr = -8
RETURN
END SUBROUTINE gaminv


SUBROUTINE dgrat(a, x, ans, qans, ierr)
!-----------------------------------------------------------------------

!        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
!                      P(A,X) AND Q(A,X)

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!     REVISED ... JAN 1992
!-------------------------
REAL (dp), INTENT(IN)  :: a, x
REAL (dp), INTENT(OUT) :: ans, qans
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL (dp) :: amn, apn, a2n, a2nm1, big, b2n, b2nm1, c, e, g, h, j, l, r, rta, &
             rtx, s, sum, t, tol, twoa, u, x0, y, z, wk(20)
INTEGER   :: i, m, n, nl1
REAL (dp), PARAMETER :: ALOG10 = 2.30258509299404568401799145468_dp,  &
                        rtpi = 1.77245385090551602729816748334_dp
!-------------------------
!     ALOG10 = LN(10)
!     RTPI   = SQRT(PI)
!-------------------------

!     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + E > 1.0 .

e = EPSILON(1.0_dp)

!-------------------------
IF (a >= 0._dp .AND. x >= 0._dp) THEN
  IF (a == 0._dp .AND. x == 0._dp) GO TO 230
  ierr = 0
  IF (a*x == 0._dp) GO TO 220

  e = MAX(e,1.d-30)
  big = 30._dp
  IF (e < 1.d-17) big = 50._dp
  x0 = 45._dp
  IF (e < 1.d-17) x0 = 68._dp

!            SELECT THE APPROPRIATE ALGORITHM

  IF (a < 1._dp) THEN
    IF (a == 0.5_dp) GO TO 200
    IF (x <= 2._dp) GO TO 100
    r = drcomp(a,x)
    IF (r == 0._dp) GO TO 190
    GO TO 150
  END IF

  IF (a < big) THEN
    IF (a > x.OR.x >= x0) GO TO 10
    twoa = a + a
    m = twoa
    l = m
    IF (twoa /= l) GO TO 10
    i = m / 2
    l = i
    IF (a == l) GO TO 120
    GO TO 130
  END IF

  l = x / a
  IF (l == 0._dp) GO TO 180
  s = 0.5_dp + (0.5_dp-l)
  z = drlog(l)
  IF (z >= 700._dp/a) GO TO 210
  y = a * z
  rta = SQRT(a)
  IF (ABS(s) <= 0.4_dp) GO TO 170

  10 r = drcomp(a,x)
  IF (r == 0._dp) GO TO 220
  IF (x > MAX(a,ALOG10)) THEN
    IF (x < x0) GO TO 150
  ELSE

!                 TAYLOR SERIES FOR P/R

    apn = a + 1._dp
    t = x / apn
    wk(1) = t
    DO n = 2, 20
      apn = apn + 1._dp
      t = t * (x/apn)
      IF (t < 1.d-3) GO TO 30
      wk(n) = t
    END DO
    n = 20

    30 sum = t
    tol = 0.5_dp * e
    40 apn = apn + 1._dp
    t = t * (x/apn)
    sum = sum + t
    IF (t > tol) GO TO 40

    nl1 = n - 1
    DO m = 1, nl1
      n = n - 1
      sum = sum + wk(n)
    END DO
    ans = (r/a) * (1._dp+sum)
    qans = 0.5_dp + (0.5_dp-ans)
    RETURN
  END IF

!                 ASYMPTOTIC EXPANSION

  amn = a - 1._dp
  t = amn / x
  wk(1) = t
  DO n = 2, 20
    amn = amn - 1._dp
    t = t * (amn/x)
    IF (ABS(t) <= 1.d-3) GO TO 70
    wk(n) = t
  END DO
  n = 20

  70 sum = t
  80 IF (ABS(t) >= e) THEN
    amn = amn - 1._dp
    t = t * (amn/x)
    sum = sum + t
    GO TO 80
  END IF

  nl1 = n - 1
  DO m = 1, nl1
    n = n - 1
    sum = sum + wk(n)
  END DO
  qans = (r/x) * (1._dp+sum)
  ans = 0.5_dp + (0.5_dp-qans)
  RETURN

!             TAYLOR SERIES FOR P(A,X)/X**A

  100 l = 3._dp
  c = x
  sum = x / (a+3._dp)
  tol = 3._dp * e / (a+1._dp)
  110 l = l + 1._dp
  c = -c * (x/l)
  t = c / (a+l)
  sum = sum + t
  IF (ABS(t) > tol) GO TO 110
  j = a * x * ((sum/6._dp-0.5_dp/(a+2._dp))*x+1._dp/(a+1._dp))

  z = a * LOG(x)
  u = EXP(z)
  h = dgam1(a)
  g = 1._dp + h
  ans = u * g * (0.5_dp+(0.5_dp-j))
  qans = 0.5_dp + (0.5_dp-ans)
  IF (ans <= 0.9_dp) RETURN

  l = drexp(z)
  qans = (u*j-l) * g - h
  IF (qans <= 0._dp) GO TO 190
  ans = 0.5_dp + (0.5_dp-qans)
  RETURN

!             FINITE SUMS FOR Q WHEN A >= 1
!                 AND 2*A IS AN INTEGER

  120 sum = EXP(-x)
  t = sum
  n = 1
  c = 0._dp
  GO TO 140

  130 rtx = SQRT(x)
  sum = derfc1(0,rtx)
  t = EXP(-x) / (rtpi*rtx)
  n = 0
  c = -0.5_dp

  140 IF (n /= i) THEN
    n = n + 1
    c = c + 1._dp
    t = (x*t) / c
    sum = sum + t
    GO TO 140
  END IF
  qans = sum
  ans = 0.5_dp + (0.5_dp-qans)
  RETURN

!              CONTINUED FRACTION EXPANSION

  150 tol = 8._dp * e
  a2nm1 = 1._dp
  a2n = 1._dp
  b2nm1 = x
  b2n = x + (1._dp-a)
  c = 1._dp
  160 a2nm1 = x * a2n + c * a2nm1
  b2nm1 = x * b2n + c * b2nm1
  c = c + 1._dp
  t = c - a
  a2n = a2nm1 + t * a2n
  b2n = b2nm1 + t * b2n

  a2nm1 = a2nm1 / b2n
  b2nm1 = b2nm1 / b2n
  a2n = a2n / b2n
  b2n = 1._dp
  IF (ABS(a2n-a2nm1/b2nm1) >= tol*a2n) GO TO 160

  qans = r * a2n
  ans = 0.5_dp + (0.5_dp-qans)
  RETURN

!                 MINIMAX APPROXIMATIONS

  170 IF (ABS(s) <= 2._dp*e .AND. a*e*e > 3.28D-3) GO TO 240
  IF (e >= 1.d-17) THEN
    CALL dgr17(a,y,l,z,rta,ans,qans)
    RETURN
  END IF
  CALL dgr29(a,y,l,z,rta,ans,qans)
  RETURN

!                     SPECIAL CASES

  180 ans = 0._dp
  qans = 1._dp
  RETURN

  190 ans = 1._dp
  qans = 0._dp
  RETURN

  200 IF (x < 0.25_dp) THEN
    ans = derf(SQRT(x))
    qans = 0.5_dp + (0.5_dp-ans)
    RETURN
  END IF
  qans = derfc1(0,SQRT(x))
  ans = 0.5_dp + (0.5_dp-qans)
  RETURN

  210 IF (ABS(s) <= 2._dp*e) GO TO 240
  220 IF (x <= a) GO TO 180
  GO TO 190
END IF

!                     ERROR RETURN

ierr = 1
ans = 2._dp
RETURN

230 ierr = 2
ans = 2._dp
RETURN

240 ierr = 3
ans = 2._dp
RETURN
END SUBROUTINE dgrat



SUBROUTINE dgr29(a,  y,  l,  z,  rta,  ans,  qans)
!-----------------------------------------------------------------------

!            ALGORITHM USING MINIMAX APPROXIMATIONS

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)     :: a, y, l, rta
REAL (dp), INTENT(IN OUT) :: z
REAL (dp), INTENT(OUT)    :: ans, qans

! Local variables
REAL (dp) :: t, u, w
REAL (dp) :: a0(7) = (/ -.234443848930188413698825870D-08,  &
        -.408902435641223939887180303D-07,  &
        -.327874000161065050049103731D-06,  &
        -.145717031728609218851588740D-05,  &
        -.372722892959910688597417881D-05,  &
        -.490033281596113358850307112D-05,  &
        -.218544851067999216147364227D-05 /), a1(7)  &
        = (/ -.162671127226300802902860047D-05,  &
        -.359791514993122440319624428D-04,  &
        -.334816794629374699945489443D-03,  &
        -.167787748352827199882047653D-02,  &
        -.462960105006279850867332060D-02,  &
        -.627269388216833251971110268D-02,  &
        -.185185185185185185185185200D-02 /), a2(7)  &
        = (/ .100841467329617467204527243D-06,  &
        .261809837060522545971782889D-05, .351658023234640143803014403D-04  &
       , .287368655528567495658887760D-03, .138385867950361368914038461D-02  &
       , .365985331203490698463644329D-02, .413359788359788359788359644D-02  &
        /), a3(7) = (/ .352304123782956092061364635D-06,  &
        .695396758348887902366951353D-05, .620467118988901865955998784D-04  &
       , .331552280167649130371474456D-03, .987931909328964685388525477D-03  &
       , .141844584435355290321010006D-02, .649434156378600823045102236D-03  &
        /), a4(7) = (/ -.260879135093022176005540138D-07,  &
        -.470448694272734954500324169D-06,  &
        -.487392507564453824976295590D-05,  &
        -.337525643163070607393381432D-04,  &
        -.173138093150706317400323103D-03,  &
        -.619343030286408407629007048D-03,  &
        -.861888290916711698604710684D-03 /)
REAL (dp) :: a5(7) = (/ -.116166342948098688243985652D-07,  &
        .506465072067030007394288471D-08, -.556701576804390213081214801D-05,  &
        -.332229941748769925615918550D-04,  &
        -.171902547619915856635305717D-03,  &
        -.548868487607991087508092013D-03,  &
        -.336798553366358151161633777D-03 /), a6(4)  &
        = (/ .118384620224413424936260301D-04,  &
        .694345283181981060040314140D-05, .209213745619758030399432459D-03,  &
        .531307936463992224884286210D-03 /), a7(5)  &
        = (/ .972342656522493967167788395D-05,  &
        .462793722775687016808279009D-04, .208913588225005764102252127D-03,  &
        .605983804794748515383615779D-03, .344367606892381545765962366D-03  &
        /), a8(5) = (/ -.231069438570167401077137510D-05,  &
        -.192877995065652524742879002D-04,  &
        -.282551884312564905942488077D-04,  &
        -.353272052089782073130912603D-03,  &
        -.652623918595320914510590273D-03 /), a9(5)  &
        = (/ -.203007139532451428594124139D-04,  &
        -.120148495117517992204095691D-03,  &
        -.377126645910917006921076652D-03,  &
        -.109151697941931403194363814D-02,  &
        -.596761290192642722092337263D-03 /), a10(4)  &
        = (/ .475862254251166503473724173D-04,  &
        -.352503880413640910997936559D-04, .580375987713106460207815603D-03  &
       , .133244544950730832649306319D-02 /)
REAL (dp) :: a11(4) = (/ .121185049262809526794966703D-03,  &
        .717725173388339108430635016D-05, .246371734409638623215800502D-02,  &
        .157972766214718575927904484D-02 /), a12(4)  &
        = (/ -.246294151509758620837749269D-03,  &
        .650624975008642297405944869D-03, -.214376520139497301154749750D-03,  &
        -.407251199495291398243480255D-02 /), a13(3)  &
        = (/ -.159520095187034545391135461D-02,  &
        -.109727312966041723997078734D-01,  &
        -.594758070915055362667114240D-02 /), a14(3)  &
        = (/ .245543970647383469794050102D-02,  &
        -.119636668153843644820445054D-01, .175722793448246103440764372D-01  &
        /), a15(2) = (/ .588261033368548917447688791D-01,  &
        .400765463491067514929787780D-01 /), a16(2)  &
        = (/ .119522261141925960204472459D+00,  &
        -.100326700196947262548667584D+00 /), a17(1)  &
        = (/ -.259949826752497731336860753D+00 /), a18(1)  &
        = (/ .724036968309299822373280436D+00 /)
REAL (dp) :: b0(9) = (/ -.129786815987713980865910767D-09,  &
        .319268409139858531586963150D-08, .597739416777031660496708557D-04  &
       , .131659965062389880196860991D-02, .138263099503103838517015533D-01  &
       , .866750030433403450681521877D-01, .349373447613102956696810725D+00  &
       , .902581259032419042347458484D+00, .139388806936391316154237713D+01  &
        /), b1(9) = (/ .361538770500640888027927000D-09,  &
        .974094440943696092434381137D-05, .275463718595762102271929980D-03  &
       , .356903970692700621824901511D-02, .276755209895072417713430394D-01  &
       , .140741499324744724262767201D+00, .482173396010404307346794795D+00  &
       , .109307843990990308990473663D+01, .151225469637089956064399494D+01  &
        /), b2(8) = (/ .144996224602847932479320241D-04,  &
        .378705615967233119938297206D-03, .457258679387716305283282667D-02  &
       , .333036784835643463383606186D-01, .160392471625881407829191009D+00  &
       , .524238095721639512312120765D+00, .114320896084982707537755002D+01  &
       , .153405837991415136438992306D+01 /), b3(8)  &
        = (/ .656342109234806261144233394D-04,  &
        .130398975231883219976260776D-02, .126418031281256648240652355D-01  &
       , .760733201461716525855765749D-01, .308149284260387354956024487D+00  &
       , .856743428738899911100227393D+00, .159678625605457556492814589D+01  &
       , .183078413578083710405050462D+01 /), b4(8)  &
        = (/ .561738585657138771286755470D-04,  &
        .104553622856827932853059322D-02, .990129468337836044520381371D-02  &
       , .590964360473404599955095091D-01, .241580582651643837306299024D+00  &
       , .686949677014349678482109368D+00, .133507902144433100426436242D+01  &
       , .162826466816694512158165085D+01 /)
REAL (dp) :: b5(7) = (/ .106576106868815233442641444D-03,  &
        .280714123386276098548285440D-02, .254669201041872409738119341D-01  &
       , .136071713023783507468096673D+00, .462890328922621047510807887D+00  &
       , .103913867517817784825064299D+01, .142263185288429590449288300D+01  &
        /), b6(9) = (/ -.633002360430352916354621750D-05,  &
        -.248639208901374031411609873D-04, .151734058829700925162000373D-03  &
       , .477475914272399601740818883D-02, .384410125775084107229541456D-01  &
       , .184699876959596092801262547D+00, .571784440733980642101712125D+00  &
       , .118432122801495778365352945D+01, .150831585220968267709550582D+01  &
        /), b7(7) = (/ .215964480325937088444595990D-03,  &
        .621296161441756044580440529D-02, .497403555098433701440032746D-01  &
       , .230812334251394761909158355D+00, .682159830165959997577293001D+00  &
       , .133753662990343866552766613D+01, .160951809815647533045690195D+01  &
        /), b8(7) = (/ .156052480203446255774109882D-02,  &
        .189231675289329563916597032D-01, .110127834209242088316741250D+00  &
       , .407929996207245634766606879D+00, .101702505946784412105505734D+01  &
       , .172269407630659768618234623D+01, .182765408802230546887514255D+01 /)
REAL (dp) :: b9(6) = (/ .108808775028021530146610124D-01,  &
        .803149717787956717154553908D-01, .335555306170768573903990019D+00  &
       , .881575022436158946373557744D+00, .156222230858412078350692234D+01  &
       , .170833470935668756293234818D+01 /), b10(6)  &
        = (/ .161103572271541189817119144D-01,  &
        .114651544043625219459951640D+00, .448280675300097555552484502D+00  &
       , .110810715319704031415255670D+01, .183146436130501918547134176D+01  &
       , .187235769169449339141968881D+01 /), b11(5)  &
        = (/ .794610889405176143379963912D-02,  &
        .131627017265860324219513170D+00, .505939635317477779328000706D+00  &
       , .116082103318559904744144217D+01, .145670749780693850410866175D+01  &
        /), b12(4) = (/ .168390445944818504703640731D+00,  &
        .653453590771198550320727688D+00, .140298208333879535577602171D+01  &
       , .162497775209192630951344224D+01 /)
REAL (dp) :: b13(4) = (/ .207815761771742289849225339D+00,  &
        .790935125477975506817064616D+00, .158706682625067673596619095D+01  &
       , .175409273929961597148916309D+01 /), b14(2)  &
        = (/ .676925518749829493412063599D+00,  &
        .100158659226079685399214158D+01 /), b15(2)  &
        = (/ .124266359850901469771032599D+01,  &
        .149189509890654955611528542D+01 /), b16(1)  &
        = (/ .536462039767059451769400255D+00 /)
REAL (dp) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,  &
                    c12, c13, c14, c15, c16
REAL (dp) :: d0(7) = (/ -.333333333333333333333333333333D+00,  &
        .833333333333333333333333333333D-01,  &
        -.148148148148148148148148148148D-01,  &
        .115740740740740740740740740741D-02,  &
        .352733686067019400352733686067D-03,  &
        -.178755144032921810699588477366D-03,  &
        .391926317852243778169704095630D-04 /), e,  &
        rt2pin = .398942280401432677939946059934_dp
!---------------------------
!     RT2PIN = 1/SQRT(2*PI)
!---------------------------
e = EXP(-y)
w = 0.5_dp * derfc1(1,SQRT(y))
u = 1._dp / a
z = SQRT(z+z)
IF (l < 1._dp) z = -z

t = ((((((a0(1)*z + a0(2))*z + a0(3))*z + a0(4))*z + a0(5))*z + a0(6))*z+  &
a0(7)) / (((((((((b0(1)*z + b0(2))*z + b0(3))*z + b0(4))*z + b0(5))*z+  &
b0(6))*z + b0(7))*z + b0(8))*z + b0(9))*z+1._dp)
c0 = ((((((t*z + d0(7))*z + d0(6))*z + d0(5))*z + d0(4))*z + d0(3))*z + d0(2)) * z + d0(1)
c1 = ((((((a1(1)*z + a1(2))*z + a1(3))*z + a1(4))*z + a1(5))*z + a1(6))*z+  &
a1(7)) / (((((((((b1(1)*z + b1(2))*z + b1(3))*z + b1(4))*z + b1(5))*z+  &
b1(6))*z + b1(7))*z + b1(8))*z + b1(9))*z+1._dp)
c2 = ((((((a2(1)*z + a2(2))*z + a2(3))*z + a2(4))*z + a2(5))*z + a2(6))*z+  &
a2(7)) / ((((((((b2(1)*z + b2(2))*z + b2(3))*z + b2(4))*z + b2(5))*z+  &
b2(6))*z + b2(7))*z + b2(8))*z+1._dp)
c3 = ((((((a3(1)*z + a3(2))*z + a3(3))*z + a3(4))*z + a3(5))*z + a3(6))*z+  &
a3(7)) / ((((((((b3(1)*z + b3(2))*z + b3(3))*z + b3(4))*z + b3(5))*z+  &
b3(6))*z + b3(7))*z + b3(8))*z+1._dp)
c4 = ((((((a4(1)*z + a4(2))*z + a4(3))*z + a4(4))*z + a4(5))*z + a4(6))*z+  &
a4(7)) / ((((((((b4(1)*z + b4(2))*z + b4(3))*z + b4(4))*z + b4(5))*z+  &
b4(6))*z + b4(7))*z + b4(8))*z+1._dp)
c5 = ((((((a5(1)*z + a5(2))*z + a5(3))*z + a5(4))*z + a5(5))*z + a5(6))*z+  &
a5(7)) / (((((((b5(1)*z + b5(2))*z + b5(3))*z + b5(4))*z + b5(5))*z+  &
b5(6))*z + b5(7))*z+1._dp)
c6 = (((a6(1)*z + a6(2))*z + a6(3))*z + a6(4)) / (((((((((b6(1)*z+  &
b6(2))*z + b6(3))*z + b6(4))*z + b6(5))*z + b6(6))*z + b6(7))*z + b6(8))*z + b6(9))*z+1._dp)
c7 = ((((a7(1)*z + a7(2))*z + a7(3))*z + a7(4))*z + a7(5)) / (((((((  &
b7(1)*z + b7(2))*z + b7(3))*z + b7(4))*z + b7(5))*z + b7(6))*z + b7(7))*z+1._dp)
c8 = ((((a8(1)*z + a8(2))*z + a8(3))*z + a8(4))*z + a8(5)) / (((((((  &
b8(1)*z + b8(2))*z + b8(3))*z + b8(4))*z + b8(5))*z + b8(6))*z + b8(7))*z+1._dp)
c9 = ((((a9(1)*z + a9(2))*z + a9(3))*z + a9(4))*z + a9(5)) / ((((((b9(1)*z  &
+b9(2))*z + b9(3))*z + b9(4))*z + b9(5))*z + b9(6))*z+1._dp)
c10 = (((a10(1)*z + a10(2))*z + a10(3))*z + a10(4)) / ((((((b10(1)*z+  &
b10(2))*z + b10(3))*z + b10(4))*z + b10(5))*z + b10(6))*z+1._dp)
c11 = (((a11(1)*z + a11(2))*z + a11(3))*z + a11(4)) / (((((b11(1)*z+  &
b11(2))*z + b11(3))*z + b11(4))*z + b11(5))*z+1._dp)
c12 = (((a12(1)*z + a12(2))*z + a12(3))*z + a12(4)) / ((((b12(1)*z+  &
b12(2))*z + b12(3))*z + b12(4))*z+1._dp)
c13 = ((a13(1)*z + a13(2))*z + a13(3)) / ((((b13(1)*z + b13(2))*z+  &
b13(3))*z + b13(4))*z+1._dp)
c14 = ((a14(1)*z + a14(2))*z + a14(3)) / ((b14(1)*z + b14(2))*z+1._dp)
c15 = (a15(1)*z + a15(2)) / ((b15(1)*z + b15(2))*z+1._dp)
c16 = (a16(1)*z + a16(2)) / (b16(1)*z+1._dp)

t = (a18(1)*u+a17(1)) * u + c16
t = (((((((((((((((t*u+c15)*u+c14)*u+c13)*u+c12)*u+c11)*u+c10)*u+  &
c9)*u+c8)*u+c7)*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1) * u + c0

IF (l >= 1._dp) THEN
  qans = e * (w + rt2pin*t/rta)
  ans = 0.5_dp + (0.5_dp-qans)
  RETURN
END IF
ans = e * (w - rt2pin*t/rta)
qans = 0.5_dp + (0.5_dp-ans)
RETURN
END SUBROUTINE dgr29



SUBROUTINE dgr17(a, y, l, z, rta, ans, qans)
!-----------------------------------------------------------------------

!            ALGORITHM USING MINIMAX APPROXIMATIONS
!                        FOR C0,...,C10

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)     :: a, y, l, rta
REAL (dp), INTENT(IN OUT) :: z
REAL (dp), INTENT(OUT)    :: ans, qans

! Local variables
REAL (dp) :: e, rt2pin = .398942280401432678_dp, t, u, w
REAL (dp) :: c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10
REAL (dp) :: a0(5) = (/ -.73324404807556026D-03,  &
        -.11758531313175796D-01, -.76816029947195974D-01,  &
        -.24232172943558393D+00, -.33333333333333333D+00 /), a1(4)  &
        = (/ -.16746784557475121D-03, -.16090334014223031D-02,  &
        -.52949366601406939D-02, -.18518518518518417D-02 /), a2(4)  &
        = (/ .12049855113125238D-04, .13743853858711134D-03,  &
        .15067356806896441D-02, .41335978835983393D-02 /), a3(4)  &
        = (/ .46318872971699924D-05, .13012396979747783D-04,  &
        .81804333975935872D-03, .64943415637082551D-03 /), a4(3)  &
        = (/ -.37567394580525597D-05, -.82794205648271314D-04,  &
        -.86188829773520181D-03 /), a5(2) = (/ -.43263341886764011D-03,  &
        -.33679854644784478D-03 /), a6(2) = (/ -.12962670089753501D-03,  &
        .53130115408837152D-03 /), a7(2) = (/ .47861364421780889D-03,  &
        .34438428473168988D-03 /), a8(2) = (/ .27086391808339115D-03,  &
        -.65256615574219131D-03 /), a9(3) = (/ .84725086921921823D-03,  &
        -.14838721516118744D-03, -.60335050249571475D-03 /), a10(2)  &
        = (/ -.19144384985654775D-02, .13324454494800656D-02 /)
REAL (dp) :: b0(6) = (/ .10555647473018528D-06,  &
        .73121701584237188D-03, .13250270182342259D-01,  &
        .10288837674434487D+00, .43024494247383254D+00,  &
        .97696518830675185D+00 /), b1(6) = (/ .12328086517283227D-05,  &
        .98671953445602142D-03, .15954049115266936D-01,  &
        .11439610256504704D+00, .45195109694529839D+00,  &
        .98426579647613593D+00 /), b2(5) = (/ .15927093345670077D-02,  &
        .22316881460606523D-01, .14009848931638062D+00,  &
        .50379606871703058D+00, .10131761625405203D+01 /), b3(4)  &
        = (/ .12414068921653593D-01, .10044290377295469D+00,  &
        .42226789458984594D+00, .90628317147366376D+00 /), b4(4)  &
        = (/ .31290397554562032D-01, .16988291247058802D+00,  &
        .57225859400072754D+00, .10057375981227881D+01 /), b5(4)  &
        = (/ .22714615451529335D-01, .17081504060220639D+00,  &
        .60019022026983067D+00, .10775200414676195D+01 /), b6(3)  &
        = (/ .65929776650152292D-01, .45957439582639129D+00,  &
        .87058903334443855D+00 /), b7(3) = (/ .27176241899664174D+00,  &
        .78991370162247144D+00, .12396875725833093D+01 /), b8(2)  &
        = (/ .44207055629598579D+00, .87002402612484571D+00 /)
!------------------------
!     RT2PIN = 1/SQRT(2*PI)
!------------------------
e = EXP(-y)
w = 0.5_dp * derfc1(1,SQRT(y))
u = 1._dp / a
z = SQRT(z+z)
IF (l < 1._dp) z = -z

c0 = ((((a0(1)*z + a0(2))*z + a0(3))*z + a0(4))*z + a0(5)) / ((((((b0(1)*z  &
      + b0(2))*z + b0(3))*z + b0(4))*z + b0(5))*z + b0(6))*z + 1._dp)
c1 = (((a1(1)*z + a1(2))*z + a1(3))*z + a1(4)) / ((((((b1(1)*z + b1(2))*z + &
     b1(3))*z + b1(4))*z + b1(5))*z + b1(6))*z + 1._dp)
c2 = (((a2(1)*z + a2(2))*z + a2(3))*z + a2(4)) / (((((b2(1)*z + b2(2))*z +  &
     b2(3))*z + b2(4))*z + b2(5))*z + 1._dp)
c3 = (((a3(1)*z + a3(2))*z + a3(3))*z + a3(4)) / ((((b3(1)*z + b3(2))*z +   &
     b3(3))*z + b3(4))*z + 1._dp)
c4 = ((a4(1)*z + a4(2))*z + a4(3)) / ((((b4(1)*z + b4(2))*z + b4(3))*z +    &
     b4(4))*z + 1._dp)
c5 = (a5(1)*z + a5(2)) / ((((b5(1)*z + b5(2))*z + b5(3))*z + b5(4))*z + 1._dp)
c6 = (a6(1)*z + a6(2)) / (((b6(1)*z + b6(2))*z + b6(3))*z + 1._dp)
c7 = (a7(1)*z + a7(2)) / (((b7(1)*z + b7(2))*z + b7(3))*z + 1._dp)
c8 = (a8(1)*z + a8(2)) / ((b8(1)*z + b8(2))*z + 1._dp)
c9 = (a9(1)*z + a9(2)) * z + a9(3)
c10 = a10(1) * z + a10(2)

t = (((((((((c10*u + c9)*u + c8)*u + c7)*u + c6)*u + c5)*u + c4)*u + c3)*u + &
    c2)*u + c1)*u + c0

IF (l >= 1._dp) THEN
  qans = e * (w + rt2pin*t/rta)
  ans = 0.5_dp + (0.5_dp - qans)
  RETURN
END IF
ans = e * (w - rt2pin*t/rta)
qans = 0.5_dp + (0.5_dp - ans)

RETURN
END SUBROUTINE dgr17



SUBROUTINE dginv(a, x, p, q, ierr)
!-----------------------------------------------------------------------

!                        REAL (dp)
!             INVERSE INCOMPLETE GAMMA RATIO FUNCTION

!     GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1.
!     THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER
!     ITERATION IS EMPLOYED.

!                        ------------

!     X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
!     AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT
!     NUMBER AVAILABLE. OTHERWISE, DGINV ATTEMPTS TO OBTAIN
!     A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE
!     IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.

!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!     WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING
!     VALUES ...

!       IERR =  0    THE SOLUTION WAS OBTAINED. ITERATION WAS
!                    NOT USED.
!       IERR>0    THE SOLUTION WAS OBTAINED. IERR ITERATIONS
!                    WERE PERFORMED.
!       IERR = -2    (INPUT ERROR) A <= 0
!       IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A
!                    IS TOO LARGE.
!       IERR = -4    (INPUT ERROR) P OR Q IS NEGATIVE, OR
!                    P + Q .NE. 1.
!       IERR = -6    10 ITERATIONS WERE PERFORMED. THE MOST
!                    RECENT VALUE OBTAINED FOR X IS GIVEN.
!                    (THIS SETTING SHOULD NEVER OCCUR.)
!       IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X.
!                    THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
!       IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE
!                    ROUTINE IS NOT CERTAIN OF ITS ACCURACY.
!                    ITERATION CANNOT BE PERFORMED IN THIS
!                    CASE. THIS SETTING CAN OCCUR ONLY WHEN
!                    P OR Q IS APPROXIMATELY 0.

!-----------------------------------------------------------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!     WRITTEN ... JANUARY 1992
!------------------------
REAL (dp), INTENT(IN)  :: a, p, q
REAL (dp), INTENT(OUT) :: x
INTEGER, INTENT(OUT)   :: ierr

REAL      :: p0, q0, x0
REAL (dp) :: am1, apn, ap1, ap2, ap3, b, c1, c2, c3, c4, c5, d, e, eps, &
             g, h, pn, qg, qn, r, rta, s, sum, s2, t, u, w, xmin, xn, y, &
             z, amin
INTEGER   :: ier, ind
REAL (dp), PARAMETER :: c = .577215664901533_dp, ln10 = 2.302585_dp, tol = 1.d-10
!------------------------
!     LN10 = LN(10)
!     C = EULER CONSTANT
!------------------------

!     ****** E AND XMIN ARE MACHINE DEPENDENT CONSTANTS. E IS THE
!            SMALLEST NUMBER FOR WHICH 1.0 + E > 1.0, AND XMIN
!            IS THE SMALLEST POSITIVE NUMBER.

e = EPSILON(1.0_dp)
xmin = TINY(1.0_dp)

!------------------------
x = 0._dp
IF (a > 0._dp) THEN
  IF (p < 0._dp.OR.q < 0._dp) GO TO 120
  t = ((p+q)-0.5_dp) - 0.5_dp
  IF (ABS(t) > 5._dp*MAX(e,1.d-30)) GO TO 120

  ierr = 0
  xmin = xmin / e
  IF ((p/e) > xmin) THEN
    IF ((q/e) <= xmin) GO TO 160
    IF (a == 1._dp) GO TO 100

    e = MAX(e,1.d-30)
    eps = 1.d3 * e
    amin = 5.d3
    IF (e < 1.d-17) amin = 2.d6
    IF (a >= amin) GO TO 30

!        GET AN INITIAL APPROXIMATION USING THE SINGLE
!         PRECISION ARITHMETIC (IF THIS IS POSSIBLE)

    p0 = p
    q0 = q
    IF (p0 /= 0.0 .AND. q0 /= 0.0) THEN
      CALL gaminv(REAL(a),x0,0.0,p0,q0,ier)
      IF (ier >= 0.0 .OR. ier == -8) THEN
        ierr = MAX(ier,0)
        IF (x0 <= 1.e34) THEN
          xn = x0
          GO TO 50
        END IF
      END IF
    END IF

    IF (a > 1._dp) GO TO 30
    xn = 0._dp

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A < 1

    g = dgamma(a+1._dp)
    qg = q * g
    IF (qg == 0._dp) GO TO 160
    b = qg / a
    IF (qg > 0.6_dp*a) GO TO 20
    IF (a < 0.30_dp .AND. b >= 0.35_dp) THEN
      t = EXP(-(b+c))
      u = t * EXP(t)
      xn = t * EXP(u)
      GO TO 50
    END IF

    IF (b >= 0.45_dp) GO TO 20
    IF (b == 0._dp) GO TO 160
    y = -LOG(b)
    s = 0.5_dp + (0.5_dp-a)
    z = LOG(y)
    t = y - s * z
    IF (b >= 0.15_dp) THEN
      xn = y - s * LOG(t) - LOG(1._dp+s/(t+1._dp))
      GO TO 80
    END IF
    IF (b > 1.d-2) THEN
      u = ((t+2._dp*(3._dp-a))*t+(2._dp-a)*(3._dp-a)) / ((t+(5._dp-a))*t+2._dp)
      xn = y - s * LOG(t) - LOG(u)
      GO TO 80
    END IF
    10 c1 = -s * z
    c2 = -s * (1._dp+c1)
    c3 = s * ((0.5_dp*c1 + (2._dp-a))*c1 + (2.5_dp-1.5_dp*a))
    c4 = -s * (((c1/3._dp+(2.5_dp-1.5_dp*a))*c1 + ((a-6._dp)*a + 7._dp))*c1  &
         +((11._dp*a - 46._dp)*a + 47._dp)/6._dp)
    c5 = -s * ((((-c1/4._dp + (11._dp*a-17._dp)/6._dp)*c1+((-3._dp*a+  &
         13._dp)*a - 13._dp))*c1 + 0.5_dp*(((2._dp*a - 25._dp)*a + 72._dp)*a-  &
         61._dp))*c1 + (((25._dp*a - 195._dp)*a + 477._dp)*a - 379._dp)/12._dp)
    xn = ((((c5/y+c4)/y+c3)/y+c2)/y+c1) + y
    GO TO 80

    20 IF (b*q <= 1.d-8) THEN
      xn = EXP(-(q/a+c))
    ELSE
      IF (p > 0.9_dp) THEN
        xn = EXP((dlnrel(-q)+dgmln1(a))/a)
      ELSE
        xn = EXP(LOG(p*g)/a)
      END IF
    END IF

    IF (xn == 0._dp) GO TO 110
    t = 0.5_dp + (0.5_dp-xn/(a+1._dp))
    xn = xn / t
    GO TO 50

!        SELECTION OF THE INITIAL APPROXIMATION XN OF X
!                       WHEN A > 1

    30 t = p - 0.5_dp
    IF (q < 0.5_dp) t = 0.5_dp - q
    CALL dpni(p,q,t,s,ier)

    rta = SQRT(a)
    s2 = s * s
    xn = (((12._dp*s2-243._dp)*s2-923._dp)*s2+1472._dp) / 204120._dp -  &
    s * (((3753._dp*s2+4353._dp)*s2-289517._dp)*s2-289717._dp) / (146966400._dp*rta)
    xn = (xn/a+s*((9._dp*s2+256._dp)*s2-433._dp)/(38880._dp*rta)) - ((  &
    3._dp*s2+7._dp)*s2-16._dp) / 810._dp
    xn = a + s * rta + (s2-1._dp) / 3._dp + s * (s2-7._dp) / (36._dp*rta) + xn / a
    xn = MAX(xn,0._dp)

    IF (a >= amin) THEN
      x = xn
      d = 0.5_dp + (0.5_dp-x/a)
      IF (ABS(d) <= 1.d-1) THEN
        IF (ABS(d) > 1.d-3) GO TO 50
        RETURN
      END IF
    END IF

    IF (p > 0.5_dp) THEN
      IF (xn < 3._dp*a) GO TO 80
      w = LOG(q)
      y = -(w+dgamln(a))
      d = MAX(2._dp,a*(a-1._dp))
      IF (y >= ln10*d) THEN
        s = 1._dp - a
        z = LOG(y)
        GO TO 10
      END IF
      t = a - 1._dp
      xn = y + t * LOG(xn) - dlnrel(-t/(xn+1._dp))
      xn = y + t * LOG(xn) - dlnrel(-t/(xn+1._dp))
      GO TO 80
    END IF

    ap1 = a + 1._dp
    IF (xn > 0.7_dp*ap1) GO TO 60
    w = LOG(p) + dgamln(ap1)
    IF (xn <= 0.15_dp*ap1) THEN
      ap2 = a + 2._dp
      ap3 = a + 3._dp
      x = EXP((w+x)/a)
      x = EXP((w+x-LOG(1.0+(x/ap1)*(1._dp+x/ap2)))/a)
      x = EXP((w+x-LOG(1.0+(x/ap1)*(1._dp+x/ap2)))/a)
      x = EXP((w+x-LOG(1.0+(x/ap1)*(1._dp+(x/ap2)*(1._dp+x/ap3))))/a)
      xn = x
      IF (xn <= 1.d-2*ap1) GO TO 60
    END IF

    apn = ap1
    t = xn / apn
    sum = 1._dp + t
    40 apn = apn + 1._dp
    t = t * (xn/apn)
    sum = sum + t
    IF (t > 1.d-4) GO TO 40
    t = w - LOG(sum)
    xn = EXP((xn+t)/a)
    xn = xn * (1._dp - (a*LOG(xn) - xn - t) / (a - xn))
    GO TO 60

!                 SCHRODER ITERATION USING P

    50 IF (p > 0.5_dp) GO TO 80
    60 IF (p <= xmin) GO TO 150
    am1 = (a-0.5_dp) - 0.5_dp

    70 IF (ierr >= 10) GO TO 130
    ierr = ierr + 1
    CALL dgrat(a, xn, pn, qn, ind)
    IF (ind /= 0) ierr = ierr + 1
    IF (pn == 0._dp .OR. qn == 0._dp) GO TO 150
    r = drcomp(a,xn)
    IF (r < xmin) GO TO 150
    t = (pn-p) / r
    w = 0.5_dp * (am1 - xn)
    IF (ABS(t) > 0.1_dp .OR. ABS(w*t) > 0.1_dp) THEN
      x = xn * (1._dp-t)
      IF (x <= 0._dp) GO TO 140
      d = ABS(t)
    ELSE

      h = t * (1._dp + w*t)
      x = xn * (1._dp - h)
      IF (x <= 0._dp) GO TO 140
      IF (ABS(w) >= 1._dp .AND. ABS(w)*t*t <= eps) RETURN
      d = ABS(h)
    END IF
    xn = x
    IF (d > tol) GO TO 70
    IF (d <= eps) RETURN
    IF (ABS(p - pn) <= tol*p) RETURN
    GO TO 70

!                 SCHRODER ITERATION USING Q

    80 IF (q <= xmin) GO TO 150
    am1 = (a - 0.5_dp) - 0.5_dp

    90 IF (ierr >= 10) GO TO 130
    ierr = ierr + 1
    CALL dgrat(a, xn, pn, qn, ind)
    IF (pn == 0._dp .OR. qn == 0._dp) GO TO 150
    r = drcomp(a,xn)
    IF (r < xmin) GO TO 150
    t = (q-qn) / r
    w = 0.5_dp * (am1-xn)
    IF (ABS(t) > 0.1_dp .OR. ABS(w*t) > 0.1_dp) THEN
      x = xn * (1._dp - t)
      IF (x <= 0._dp) GO TO 140
      d = ABS(t)
    ELSE

      h = t * (1._dp + w*t)
      x = xn * (1._dp - h)
      IF (x <= 0._dp) GO TO 140
      IF (ABS(w) >= 1._dp .AND. ABS(w)*t*t <= eps) RETURN
      d = ABS(h)
    END IF
    xn = x
    IF (d > tol) GO TO 90
    IF (d <= eps) RETURN
    IF (ABS(q-qn) <= tol*q) RETURN
    GO TO 90
  END IF

!                       SPECIAL CASES

  ierr = -8
  RETURN

  100 IF (q >= 0.9_dp) THEN
    x = -dlnrel(-p)
    RETURN
  END IF
  x = -LOG(q)
  RETURN
END IF

!                       ERROR RETURN

ierr = -2
RETURN

110 ierr = -3
RETURN

120 ierr = -4
RETURN

130 ierr = -6
RETURN

140 ierr = -7
RETURN

150 x = xn
ierr = -8
RETURN

160 x = HUGE(1.0_dp)
ierr = -8
RETURN
END SUBROUTINE dginv


SUBROUTINE bratio(a, b, x, y, w, w1, ierr)
!-----------------------------------------------------------------------

!            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)

!                     --------------------

!     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X <= 1
!     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES

!                      W  = IX(A,B)
!                      W1 = 1 - IX(A,B)

!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
!     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
!     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
!     ONE OF THE FOLLOWING VALUES ...

!        IERR = 1  IF A OR B IS NEGATIVE
!        IERR = 2  IF A = B = 0
!        IERR = 3  IF X < 0 OR X > 1
!        IERR = 4  IF Y < 0 OR Y > 1
!        IERR = 5  IF X + Y .NE. 1
!        IERR = 6  IF X = A = 0
!        IERR = 7  IF Y = B = 0

!--------------------
!     WRITTEN BY ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN, VIRGINIA
!     REVISED ... APRIL 1993
!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: a, b, x, y
REAL, INTENT(OUT)    :: w, w1
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL    :: a0, b0, eps, lambda, t, x0, y0, z
INTEGER :: ierr1, ind, n
!-----------------------------------------------------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
!            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS > 1.0

eps = EPSILON(1.0)

!-----------------------------------------------------------------------
w = 0.0
w1 = 0.0
IF (a >= 0.0 .AND. b >= 0.0) THEN
  IF (a == 0.0 .AND. b == 0.0) GO TO 150
  IF (x < 0.0.OR.x > 1.0) GO TO 160
  IF (y < 0.0.OR.y > 1.0) GO TO 170
  z = ((x+y)-0.5) - 0.5
  IF (ABS(z) > 3.0*eps) GO TO 180

  ierr = 0
  IF (x /= 0.0) THEN
    IF (y == 0.0) GO TO 110
    IF (a == 0.0) GO TO 120
    IF (b == 0.0) GO TO 100

    eps = MAX(eps,1.e-15)
    IF (MAX(a,b) < 1.e-3*eps) GO TO 140

    ind = 0
    a0 = a
    b0 = b
    x0 = x
    y0 = y
    IF (MIN(a0,b0) <= 1.0) THEN

!             PROCEDURE FOR A0 <= 1 OR B0 <= 1

      IF (x > 0.5) THEN
        ind = 1
        a0 = b
        b0 = a
        x0 = y
        y0 = x
      END IF

      IF (b0 < MIN(eps,eps*a0)) GO TO 10
      IF (a0 < MIN(eps,eps*b0) .AND. b0*x0 <= 1.0) GO TO 20
      IF (MAX(a0,b0) <= 1.0) THEN
        IF (a0 >= MIN(0.2,b0)) GO TO 30
        IF (x0**a0 <= 0.9) GO TO 30
        IF (x0 >= 0.3) GO TO 40
        n = 20
        GO TO 60
      END IF

      IF (b0 <= 1.0) GO TO 30
      IF (x0 >= 0.3) GO TO 40
      IF (x0 < 0.1) THEN
        IF ((x0*b0)**a0 <= 0.7) GO TO 30
      END IF
      IF (b0 > 15.0) GO TO 70
      n = 20
      GO TO 60
    END IF

!             PROCEDURE FOR A0 > 1 AND B0 > 1

    IF (a <= b) THEN
      lambda = a - (a+b) * x
    ELSE
      lambda = (a+b) * y - b
    END IF
    IF (lambda < 0.0) THEN
      ind = 1
      a0 = b
      b0 = a
      x0 = y
      y0 = x
      lambda = ABS(lambda)
    END IF

    IF (b0 < 40.0 .AND. b0*x0 <= 0.7) GO TO 30
    IF (b0 < 40.0) GO TO 80
    IF (a0 <= b0) THEN
      IF (a0 <= 100.0) GO TO 50
      IF (lambda > 0.03*a0) GO TO 50
      GO TO 90
    END IF
    IF (b0 <= 100.0) GO TO 50
    IF (lambda > 0.03*b0) GO TO 50
    GO TO 90

!            EVALUATION OF THE APPROPRIATE ALGORITHM

    10 w = fpser(a0,b0,x0,eps)
    w1 = 0.5 + (0.5-w)
    GO TO 130

    20 w1 = apser(a0,b0,x0,eps)
    w = 0.5 + (0.5-w1)
    GO TO 130

    30 w = bpser(a0,b0,x0,eps)
    w1 = 0.5 + (0.5-w)
    GO TO 130

    40 w1 = bpser(b0, a0, y0, eps)
    w = 0.5 + (0.5-w1)
    GO TO 130

    50 w = bfrac(a0, b0, x0, y0, lambda, 15.0*eps)
    w1 = 0.5 + (0.5-w)
    GO TO 130

    60 w1 = bup(b0, a0, y0, x0, n, eps)
    b0 = b0 + n
    70 CALL bgrat(b0, a0, y0, x0, w1, eps, ierr1)
    IF (ierr1 /= 0) WRITE(*, *) '** Error in call to BGRAT from BRATIO **'
    w = 0.5 + (0.5-w1)
    GO TO 130

    80 n = b0
    b0 = b0 - n
    IF (b0 == 0.0) THEN
      n = n - 1
      b0 = 1.0
    END IF
    w = bup(b0, a0, y0, x0, n, eps)
    IF (x0 <= 0.7) THEN
      w = w + bpser(a0, b0, x0, eps)
      w1 = 0.5 + (0.5-w)
      GO TO 130
    END IF

    IF (a0 <= 15.0) THEN
      n = 20
      w = w + bup(a0, b0, x0, y0, n, eps)
      a0 = a0 + n
    END IF
    CALL bgrat(a0, b0, x0, y0, w, eps, ierr1)
    w1 = 0.5 + (0.5-w)
    GO TO 130

    90 w = basym(a0, b0, lambda, 100.0*eps)
    w1 = 0.5 + (0.5-w)
    GO TO 130
  END IF

!               TERMINATION OF THE PROCEDURE

  IF (a == 0.0) GO TO 190
  100 w = 0.0
  w1 = 1.0
  RETURN

  110 IF (b == 0.0) GO TO 200
  120 w = 1.0
  w1 = 0.0
  RETURN

  130 IF (ind == 0) RETURN
  t = w
  w = w1
  w1 = t
  RETURN

!           PROCEDURE FOR A AND B < 1.E-3*EPS

  140 w = b / (a+b)
  w1 = a / (a+b)
  RETURN
END IF

!                       ERROR RETURN

ierr = 1
RETURN
150 ierr = 2
RETURN
160 ierr = 3
RETURN
170 ierr = 4
RETURN
180 ierr = 5
RETURN
190 ierr = 6
RETURN
200 ierr = 7
RETURN
END SUBROUTINE bratio


FUNCTION fpser(a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------

!                 EVALUATION OF I (A,B)
!                                X

!          FOR B < MIN(EPS,EPS*A) AND X <= 0.5.

!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a, b, x, eps
REAL             :: fn_val

! Local variables
REAL :: an, c, s, t, tol

!                  SET  FPSER = X**A

fn_val = 1.0
IF (a > 1.e-3*eps) THEN
  fn_val = 0.0
  t = a * LOG(x)
  IF (t < exparg(1)) RETURN
  fn_val = EXP(t)
END IF

!                NOTE THAT 1/B(A,B) = B

fn_val = (b/a) * fn_val
tol = eps / a
an = a + 1.0
t = x
s = t / an
10 an = an + 1.0
t = x * t
c = t / an
s = s + c
IF (ABS(c) > tol) GO TO 10

fn_val = fn_val * (1.0+a*s)
RETURN
END FUNCTION fpser


FUNCTION apser(a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
!     A <= MIN(EPS,EPS*B), B*X <= 1, AND X <= 0.5. USED WHEN
!     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a, b, x, eps
REAL             :: fn_val

! Local variables
REAL            :: aj, bx, c, j, s, t, tol
REAL, PARAMETER :: g = .577215664901533
!--------------------
bx = b * x
t = x - bx
IF (b*eps <= 2.e-2) THEN
  c = LOG(x) + psi(b) + g + t
ELSE
  c = LOG(bx) + g + t
END IF

tol = 5.0 * eps * ABS(c)
j = 1.0
s = 0.0
10 j = j + 1.0
t = t * (x-bx/j)
aj = t / j
s = s + aj
IF (ABS(aj) > tol) GO TO 10

fn_val = -a * (c+s)
RETURN
END FUNCTION apser


FUNCTION bpser(a, b, x, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B <= 1
!     OR B*X <= 0.7.  EPS IS THE TOLERANCE USED.
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a, b, x, eps
REAL             :: fn_val

! Local variables
REAL    :: apb, a0, b0, c, n, sum, t, tol, u, w, z
INTEGER :: i, m

fn_val = 0.0
IF (x == 0.0) RETURN
!-----------------------------------------------------------------------
!            COMPUTE THE FACTOR X**A/(A*BETA(A,B))
!-----------------------------------------------------------------------
a0 = MIN(a,b)
IF (a0 >= 1.0) THEN
  z = a * LOG(x) - betaln(a,b)
  fn_val = EXP(z) / a
ELSE
  b0 = MAX(a,b)
  IF (b0 < 8.0) THEN
    IF (b0 <= 1.0) THEN

!            PROCEDURE FOR A0 < 1 AND B0 <= 1

      fn_val = x ** a
      IF (fn_val == 0.0) RETURN

      apb = a + b
      IF (apb <= 1.0) THEN
        z = 1.0 + gam1(apb)
      ELSE
        u = DBLE(a) + DBLE(b) - 1._dp
        z = (1.0+gam1(u)) / apb
      END IF

      c = (1.0+gam1(a)) * (1.0+gam1(b)) / z
      fn_val = fn_val * c * (b/apb)
      GO TO 20
    END IF

!         PROCEDURE FOR A0 < 1 AND 1 < B0 < 8

    u = gamln1(a0)
    m = b0 - 1.0
    IF (m >= 1) THEN
      c = 1.0
      DO i = 1, m
        b0 = b0 - 1.0
        c = c * (b0/(a0+b0))
      END DO
      u = LOG(c) + u
    END IF

    z = a * LOG(x) - u
    b0 = b0 - 1.0
    apb = a0 + b0
    IF (apb <= 1.0) THEN
      t = 1.0 + gam1(apb)
    ELSE
      u = DBLE(a0) + DBLE(b0) - 1._dp
      t = (1.0+gam1(u)) / apb
    END IF
    fn_val = EXP(z) * (a0/a) * (1.0+gam1(b0)) / t
  ELSE

!            PROCEDURE FOR A0 < 1 AND B0 >= 8

    u = gamln1(a0) + algdiv(a0,b0)
    z = a * LOG(x) - u
    fn_val = (a0/a) * EXP(z)
  END IF
END IF
20 IF (fn_val == 0.0.OR.a <= 0.1*eps) RETURN
!-----------------------------------------------------------------------
!                     COMPUTE THE SERIES
!-----------------------------------------------------------------------
sum = 0.0
n = 0.0
c = 1.0
tol = eps / a
30 n = n + 1.0
c = c * (0.5+(0.5-b/n)) * x
w = c / (a+n)
sum = sum + w
IF (ABS(w) > tol) GO TO 30
fn_val = fn_val * (1.0+a*sum)
RETURN
END FUNCTION bpser


FUNCTION bup(a, b, x, y, n, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
!     EPS IS THE TOLERANCE USED.
!-----------------------------------------------------------------------
REAL, INTENT(IN)    :: a, b, x, y, eps
INTEGER, INTENT(IN) :: n
REAL                :: fn_val

! Local variables
REAL    :: apb, ap1, d, l, r, t, w
INTEGER :: i, k, kp1, mu, nm1

!          OBTAIN THE SCALING FACTOR EXP(-MU) AND
!             EXP(MU)*(X**A*Y**B/BETA(A,B))/A

apb = a + b
ap1 = a + 1.0
mu = 0
d = 1.0
IF (n /= 1 .AND. a >= 1.0) THEN
  IF (apb >= 1.1*ap1) THEN
    mu = ABS(exparg(1))
    k = exparg(0)
    IF (k < mu) mu = k
    t = mu
    d = EXP(-t)
  END IF
END IF

fn_val = brcmp1(REAL(mu), a, b, x, y) / a
IF (n == 1 .OR. fn_val == 0.0) RETURN
nm1 = n - 1
w = d

!          LET K BE THE INDEX OF THE MAXIMUM TERM

k = 0
IF (b > 1.0) THEN
  IF (y <= 1.e-4) THEN
    k = nm1
  ELSE
    r = (b-1.0) * x / y - a
    IF (r < 1.0) GO TO 20
    k = nm1
    t = nm1
    IF (r < t) k = r
  END IF

!          ADD THE INCREASING TERMS OF THE SERIES

  DO i = 1, k
    l = i - 1
    d = ((apb+l)/(ap1+l)) * x * d
    w = w + d
  END DO
  IF (k == nm1) GO TO 40
END IF

!          ADD THE REMAINING TERMS OF THE SERIES

20 kp1 = k + 1
DO i = kp1, nm1
  l = i - 1
  d = ((apb+l)/(ap1+l)) * x * d
  w = w + d
  IF (d <= eps*w) GO TO 40
END DO

!               TERMINATE THE PROCEDURE

40 fn_val = fn_val * w
RETURN
END FUNCTION bup


FUNCTION bfrac(a, b, x, y, lambda, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B > 1.
!     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a, b, x, y, lambda, eps
REAL             :: fn_val

! Local variables
REAL :: alpha, an, anp1, beta, bn, bnp1, c, c0, c1, e, n, yp1, p, r, r0, &
        s, t, w

fn_val = brcomp(a,b,x,y)
IF (fn_val == 0.0) RETURN

c = 1.0 + lambda
c0 = b / a
c1 = 1.0 + 1.0 / a
yp1 = y + 1.0

n = 0.0
p = 1.0
s = a + 1.0
an = 0.0
bn = 1.0
anp1 = 1.0
bnp1 = c / c1
r = c1 / c

!        CONTINUED FRACTION CALCULATION

10 n = n + 1.0
t = n / a
w = n * (b-n) * x
e = a / s
alpha = (p*(p+c0)*e*e) * (w*x)
IF (alpha > 0.0) THEN
  e = (1.0+t) / (c1+t+t)
  beta = n + w / s + e * (c+n*yp1)
  p = 1.0 + t
  s = s + 2.0

!        UPDATE AN, BN, ANP1, AND BNP1

  t = alpha * an + beta * anp1
  an = anp1
  anp1 = t
  t = alpha * bn + beta * bnp1
  bn = bnp1
  bnp1 = t
  r0 = r
  r = anp1 / bnp1
  IF (ABS(r-r0) > eps*r) THEN

!        RESCALE AN, BN, ANP1, AND BNP1

    an = an / bnp1
    bn = bn / bnp1
    anp1 = r
    bnp1 = 1.0
    GO TO 10
  END IF
END IF

!                 TERMINATION

fn_val = fn_val * r
RETURN
END FUNCTION bfrac


FUNCTION brcomp(a, b, x, y) RESULT(fn_val)
!-----------------------------------------------------------------------
!               EVALUATION OF X**A*Y**B/BETA(A,B)
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: a, b, x, y
REAL             :: fn_val

! Local variables
REAL            :: apb, a0, b0, c, e, h, lambda, lnx, lny, t, u, v, x0, y0, z
REAL, PARAMETER :: const = .398942280401433
INTEGER         :: i, n
!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------

fn_val = 0.0
IF (x == 0.0.OR.y == 0.0) RETURN
a0 = MIN(a,b)
IF (a0 < 8.0) THEN

  IF (x <= 0.375) THEN
    lnx = LOG(x)
    lny = alnrel(-x)
  ELSE
    IF (y <= 0.375) THEN
      lnx = alnrel(-y)
      lny = LOG(y)
    ELSE
      lnx = LOG(x)
      lny = LOG(y)
    END IF
  END IF

  z = a * lnx + b * lny
  IF (a0 >= 1.0) THEN
    z = z - betaln(a,b)
    fn_val = EXP(z)
    RETURN
  END IF
!-----------------------------------------------------------------------
!              PROCEDURE FOR A < 1 OR B < 1
!-----------------------------------------------------------------------
  b0 = MAX(a,b)
  IF (b0 < 8.0) THEN
    IF (b0 <= 1.0) THEN

!                   ALGORITHM FOR B0 <= 1

      fn_val = EXP(z)
      IF (fn_val == 0.0) RETURN

      apb = a + b
      IF (apb <= 1.0) THEN
        z = 1.0 + gam1(apb)
      ELSE
        u = DBLE(a) + DBLE(b) - 1._dp
        z = (1.0 + gam1(u)) / apb
      END IF

      c = (1.0+gam1(a)) * (1.0+gam1(b)) / z
      fn_val = fn_val * (a0*c) / (1.0 + a0/b0)
      RETURN
    END IF

!                ALGORITHM FOR 1 < B0 < 8

    u = gamln1(a0)
    n = b0 - 1.0
    IF (n >= 1) THEN
      c = 1.0
      DO i = 1, n
        b0 = b0 - 1.0
        c = c * (b0/(a0+b0))
      END DO
      u = LOG(c) + u
    END IF

    z = z - u
    b0 = b0 - 1.0
    apb = a0 + b0
    IF (apb <= 1.0) THEN
      t = 1.0 + gam1(apb)
    ELSE
      u = DBLE(a0) + DBLE(b0) - 1._dp
      t = (1.0 + gam1(u)) / apb
    END IF
    fn_val = a0 * EXP(z) * (1.0 + gam1(b0)) / t
    RETURN
  END IF

!                   ALGORITHM FOR B0 >= 8

  u = gamln1(a0) + algdiv(a0,b0)
  fn_val = a0 * EXP(z-u)
  RETURN
END IF
!-----------------------------------------------------------------------
!              PROCEDURE FOR A >= 8 AND B >= 8
!-----------------------------------------------------------------------
IF (a <= b) THEN
  h = a / b
  x0 = h / (1.0+h)
  y0 = 1.0 / (1.0+h)
  lambda = a - (a+b) * x
ELSE
  h = b / a
  x0 = 1.0 / (1.0+h)
  y0 = h / (1.0+h)
  lambda = (a+b) * y - b
END IF

e = -lambda / a
IF (ABS(e) <= 0.6) THEN
  u = rlog1(e)
ELSE
  u = e - LOG(x/x0)
END IF

e = lambda / b
IF (ABS(e) <= 0.6) THEN
  v = rlog1(e)
ELSE
  v = e - LOG(y/y0)
END IF

z = EXP(-(a*u+b*v))
fn_val = const * SQRT(b*x0) * z * EXP(-bcorr(a,b))
RETURN
END FUNCTION brcomp


FUNCTION brcmp1(mu, a, b, x, y) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
!-----------------------------------------------------------------------
REAL, INTENT(IN)    :: mu, a, b, x, y
REAL                :: fn_val

! Local variables
REAL            :: apb, a0, b0, c, e, h, lambda, lnx, lny, t, u, v, x0, y0, z
REAL, PARAMETER :: const = .398942280401433
INTEGER         :: i, n
!-----------------
!     CONST = 1/SQRT(2*PI)
!-----------------

a0 = MIN(a,b)
IF (a0 < 8.0) THEN

  IF (x <= 0.375) THEN
    lnx = LOG(x)
    lny = alnrel(-x)
  ELSE
    IF (y <= 0.375) THEN
      lnx = alnrel(-y)
      lny = LOG(y)
    ELSE
      lnx = LOG(x)
      lny = LOG(y)
    END IF
  END IF

  z = a * lnx + b * lny
  IF (a0 >= 1.0) THEN
    z = z - betaln(a,b)
    fn_val = esum(mu,z)
    RETURN
  END IF
!-----------------------------------------------------------------------
!              PROCEDURE FOR A < 1 OR B < 1
!-----------------------------------------------------------------------
  b0 = MAX(a,b)
  IF (b0 < 8.0) THEN
    IF (b0 <= 1.0) THEN

!                   ALGORITHM FOR B0 <= 1

      fn_val = esum(mu,z)
      IF (fn_val == 0.0) RETURN

      apb = a + b
      IF (apb <= 1.0) THEN
        z = 1.0 + gam1(apb)
      ELSE
        u = DBLE(a) + DBLE(b) - 1._dp
        z = (1.0 + gam1(u)) / apb
      END IF

      c = (1.0+gam1(a)) * (1.0 + gam1(b)) / z
      fn_val = fn_val * (a0*c) / (1.0 + a0/b0)
      RETURN
    END IF

!                ALGORITHM FOR 1 < B0 < 8

    u = gamln1(a0)
    n = b0 - 1.0
    IF (n >= 1) THEN
      c = 1.0
      DO i = 1, n
        b0 = b0 - 1.0
        c = c * (b0/(a0+b0))
      END DO
      u = LOG(c) + u
    END IF

    z = z - u
    b0 = b0 - 1.0
    apb = a0 + b0
    IF (apb <= 1.0) THEN
      t = 1.0 + gam1(apb)
    ELSE
      u = DBLE(a0) + DBLE(b0) - 1._dp
      t = (1.0+gam1(u)) / apb
    END IF
    fn_val = a0 * esum(mu,z) * (1.0 + gam1(b0)) / t
    RETURN
  END IF

!                   ALGORITHM FOR B0 >= 8

  u = gamln1(a0) + algdiv(a0,b0)
  fn_val = a0 * esum(mu,z-u)
  RETURN
END IF
!-----------------------------------------------------------------------
!              PROCEDURE FOR A >= 8 AND B >= 8
!-----------------------------------------------------------------------
IF (a <= b) THEN
  h = a / b
  x0 = h / (1.0+h)
  y0 = 1.0 / (1.0+h)
  lambda = a - (a+b) * x
ELSE
  h = b / a
  x0 = 1.0 / (1.0+h)
  y0 = h / (1.0+h)
  lambda = (a+b) * y - b
END IF

e = -lambda / a
IF (ABS(e) <= 0.6) THEN
  u = rlog1(e)
ELSE
  u = e - LOG(x/x0)
END IF

e = lambda / b
IF (ABS(e) <= 0.6) THEN
  v = rlog1(e)
ELSE
  v = e - LOG(y/y0)
END IF

z = esum(mu,-(a*u+b*v))
fn_val = const * SQRT(b*x0) * z * EXP(-bcorr(a,b))
RETURN
END FUNCTION brcmp1


SUBROUTINE bgrat(a, b, x, y, w, eps, ierr)
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
!     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
!     THAT A >= 15 AND B <= 1.  EPS IS THE TOLERANCE USED.
!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: a, b, x, y, eps
REAL, INTENT(OUT)    :: w
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL    :: bm1, bp2n, c(30), cn, coef, d(30), dj, j, l, lnx, nu, n2, q, &
           r, s, sum, t, tol, t2, u, v, z
INTEGER :: i, n, nm1

bm1 = (b-0.5) - 0.5
nu = a + 0.5 * bm1
IF (y <= 0.375) THEN
  lnx = alnrel(-y)
ELSE
  lnx = LOG(x)
END IF
z = -nu * lnx
IF (b*z /= 0.0) THEN

!                 COMPUTATION OF THE EXPANSION
!                 SET R = EXP(-Z)*Z**B/GAMMA(B)

  r = b * (1.0 + gam1(b)) * EXP(b*LOG(z))
  r = r * EXP(a*lnx) * EXP(0.5*bm1*lnx)
  u = algdiv(b,a) + b * LOG(nu)
  u = r * EXP(-u)
  IF (u /= 0.0) THEN
    CALL grat1(b, z, r, t, q, eps)
    tol = 15.0 * eps
    v = 0.25 * (1.0/nu) ** 2
    t2 = 0.25 * lnx * lnx
    l = w / u
    j = q / r
    sum = j
    t = 1.0
    cn = 1.0
    n2 = 0.0
    DO n = 1, 30
      bp2n = b + n2
      j = (bp2n*(bp2n + 1.0)*j + (z + bp2n + 1.0)*t) * v
      n2 = n2 + 2.0
      t = t * t2
      cn = cn / (n2*(n2 + 1.0))
      c(n) = cn
      s = 0.0
      IF (n /= 1) THEN
        nm1 = n - 1
        coef = b - n
        DO i = 1, nm1
          s = s + coef * c(i) * d(n-i)
          coef = coef + b
        END DO
      END IF
      d(n) = bm1 * cn + s / n
      dj = d(n) * j
      sum = sum + dj
      IF (sum <= 0.0) GO TO 40
      IF (ABS(dj) <= tol*(sum+l)) GO TO 30
    END DO

!                    ADD THE RESULTS TO W

    30 ierr = 0
    w = w + u * sum
    RETURN
  END IF
END IF

!               THE EXPANSION CANNOT BE COMPUTED

40 ierr = 1
RETURN
END SUBROUTINE bgrat


SUBROUTINE grat1(a, x, r, p, q, eps)
!-----------------------------------------------------------------------
!           EVALUATION OF P(A,X) AND Q(A,X) WHERE A <= 1 AND
!        THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A)
!-----------------------------------------------------------------------
REAL, INTENT(IN)  :: a, x, r, eps
REAL, INTENT(OUT) :: p, q

! Local variables
REAL :: a2n, a2nm1, an, b2n, b2nm1, c, g, h, j, l, sum, t, tol, w, z

!--------------------
IF (a*x == 0.0) GO TO 70
IF (a == 0.5) GO TO 60
IF (x >= 1.1) THEN
ELSE

!             TAYLOR SERIES FOR P(A,X)/X**A

  an = 3.0
  c = x
  sum = x / (a+3.0)
  tol = 3.0 * eps / (a+1.0)
  10 an = an + 1.0
  c = -c * (x/an)
  t = c / (a+an)
  sum = sum + t
  IF (ABS(t) > tol) GO TO 10
  j = a * x * ((sum/6.0 - 0.5/(a + 2.0))*x + 1.0/(a + 1.0))

  z = a * LOG(x)
  h = gam1(a)
  g = 1.0 + h
  IF (x >= 0.25) THEN
    IF (a < x/2.59) GO TO 20
  ELSE
    IF (z > -.13394) GO TO 20
  END IF

  w = EXP(z)
  p = w * g * (0.5+(0.5-j))
  q = 0.5 + (0.5-p)
  RETURN

  20 l = rexp(z)
  q = ((0.5 + (0.5+l))*j-l) * g - h
  IF (q <= 0.0) GO TO 50
  p = 0.5 + (0.5-q)
  RETURN
END IF

!              CONTINUED FRACTION EXPANSION

tol = 8.0 * eps
a2nm1 = 1.0
a2n = 1.0
b2nm1 = x
b2n = x + (1.0-a)
c = 1.0
30 a2nm1 = x * a2n + c * a2nm1
b2nm1 = x * b2n + c * b2nm1
c = c + 1.0
a2n = a2nm1 + (c-a) * a2n
b2n = b2nm1 + (c-a) * b2n
a2nm1 = a2nm1 / b2n
b2nm1 = b2nm1 / b2n
a2n = a2n / b2n
b2n = 1.0
IF (ABS(a2n-a2nm1/b2nm1) >= tol*a2n) GO TO 30

q = r * a2n
p = 0.5 + (0.5-q)
RETURN

!                SPECIAL CASES

40 p = 0.0
q = 1.0
RETURN

50 p = 1.0
q = 0.0
RETURN

60 IF (x < 0.25) THEN
  p = erf(SQRT(x))
  q = 0.5 + (0.5-p)
  RETURN
END IF
q = erfc1(0,SQRT(x))
p = 0.5 + (0.5-q)
RETURN

70 IF (x <= a) GO TO 40
GO TO 50
END SUBROUTINE grat1


FUNCTION basym(a, b, lambda, eps) RESULT(fn_val)
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
!     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
!     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
!     A AND B ARE GREATER THAN OR EQUAL TO 15.
!-----------------------------------------------------------------------
REAL, INTENT(IN)  :: a, b, lambda, eps
REAL              :: fn_val

! Local variables
REAL, PARAMETER    :: e0 = 1.12837916709551, e1 = .353553390593274
REAL               :: a0(21), b0(21), bsum, c(21), d(21), dsum, f, h, h2, hn, &
                      j0, j1, r, r0, r1, s, sum, t, t0, t1, u, w, w0,  &
                      z, z0, z2, zn, znm1
INTEGER, PARAMETER :: num = 20
INTEGER            :: i, im1, imj, j, m, mmj, mm1, n, np1
!------------------------
!     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
!            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
!            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.

!------------------------
fn_val = 0.0
IF (a < b) THEN
  h = a / b
  r0 = 1.0 / (1.0+h)
  r1 = (b-a) / b
  w0 = 1.0 / SQRT(a*(1.0+h))
ELSE
  h = b / a
  r0 = 1.0 / (1.0+h)
  r1 = (b-a) / a
  w0 = 1.0 / SQRT(b*(1.0+h))
END IF

f = a * rlog1(-lambda/a) + b * rlog1(lambda/b)
t = EXP(-f)
IF (t == 0.0) RETURN
z0 = SQRT(f)
z = 0.5 * (z0/e1)
z2 = f + f

a0(1) = (2.0/3.0) * r1
c(1) = -0.5 * a0(1)
d(1) = -c(1)
j0 = (0.5/e0) * erfc1(1,z0)
j1 = e1
sum = j0 + d(1) * w0 * j1

s = 1.0
h2 = h * h
hn = 1.0
w = w0
znm1 = z
zn = z2
DO n = 2, num, 2
  hn = h2 * hn
  a0(n) = 2.0 * r0 * (1.0+h*hn) / (n+2.0)
  np1 = n + 1
  s = s + hn
  a0(np1) = 2.0 * r1 * s / (n+3.0)

  DO i = n, np1
    r = -0.5 * (i+1.0)
    b0(1) = r * a0(1)
    DO m = 2, i
      bsum = 0.0
      mm1 = m - 1
      DO j = 1, mm1
        mmj = m - j
        bsum = bsum + (j*r-mmj) * a0(j) * b0(mmj)
      END DO
      b0(m) = r * a0(m) + bsum / m
    END DO
    c(i) = b0(i) / (i+1.0)

    dsum = 0.0
    im1 = i - 1
    DO j = 1, im1
      imj = i - j
      dsum = dsum + d(imj) * c(j)
    END DO
    d(i) = -(dsum+c(i))
  END DO

  j0 = e1 * znm1 + (n-1.0) * j0
  j1 = e1 * zn + n * j1
  znm1 = z2 * znm1
  zn = z2 * zn
  w = w0 * w
  t0 = d(n) * w * j0
  w = w0 * w
  t1 = d(np1) * w * j1
  sum = sum + (t0+t1)
  IF ((ABS(t0)+ABS(t1)) <= eps*sum) GO TO 60
END DO

60 u = EXP(-bcorr(a,b))
fn_val = e0 * t * u * sum
RETURN
END FUNCTION basym


SUBROUTINE isubx(a0, b0, x0, p, ierr, eps)
REAL, INTENT(IN)     :: a0, b0, x0, eps
REAL, INTENT(OUT)    :: p
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL    :: a, afrac, am1, an, b, bfrac, c, gamrat, i, j, k, lambda, m, n, &
           p0, rty, sum, t, temp, tol, x, y
INTEGER :: imax, ind, l
REAL, PARAMETER :: w(10) = (/ 6.6671344308688E-2,  &
          1.4945134915058E-1, 2.1908636251598E-1, 2.6926671931000E-1,  &
          2.9552422471475E-1, 2.9552422471475E-1, 2.6926671931000E-1,  &
          2.1908636251598E-1, 1.4945134915058E-1, 6.6671344308688E-2 /),  &
    z(10) = (/ 1.3046735791414E-2, 6.7468316655507E-2, 1.6029521585049E-1,  &
               2.8330230293538E-1, 4.2556283050918E-1, 5.7443716949081E-1,  &
               7.1669769706462E-1, 8.3970478414951E-1, 9.3253168334449E-1,  &
               9.8695326420859E-1 /),  &
    pihalf = 1.5707963267949, rpinv = .56418958354776
!     -------------------
!     ****** MAX IS A MACHINE DEPENDENT CONSTANT. MAX IS THE
!            LARGEST POSITIVE INTEGER THAT MAY BE USED.

!                     MAX = IPMPAR(3)
imax = HUGE(3)

!     -------------------
a = a0
b = b0
x = x0
y = 0.5 + (0.5-x)

!                        CHECK THE ARGUMENTS

p = 0.0
ierr = 1
IF (a < 0.5.OR.b < 0.5) GO TO 180
IF (x /= 0.0 .AND. x /= 1.0) THEN
  IF (x < 0.0) GO TO 170
  m = imax
  IF (a >= m.OR.b > 70.0.OR.y < 0.0) GO TO 190
  k = INT(a)
  j = INT(b)
  afrac = a - k
  bfrac = b - j
  IF ((afrac /= 0.0 .AND. afrac /= 0.5).OR.(bfrac /= 0.0 .AND. bfrac /= 0.5)) &
                  GO TO 200
  IF (a >= 5000.0 .AND. x < 0.96) RETURN

!                      CHECK IF B IS AN INTEGER

  ind = 0
  tol = 0.5 * MAX(eps,1.e-11)
  IF (bfrac /= 0.0) GO TO 50
  IF (afrac /= 0.0) GO TO 20
  IF (a >= b) GO TO 20

!                        INTERCHANGE A AND B

  10 ind = 1
  t = b
  b = a
  a = t
  t = y
  y = x
  x = t
  t = j
  j = k
  k = t

!                        COMPUTE EXPANSION 14

  20 am1 = a - 1.0
  n = 1.0
  IF (am1 >= 0.5) THEN
    n = j
    IF (y < 2.0*j*x) THEN
      t = am1 * y / x + 1.0
      IF (t < j) n = INT(t)
    END IF
  END IF

  i = n - 1.0
  c = (a*LOG(x)+i*alnrel(-x)) - blnd(a,i)
  IF (c > -30) THEN
    tol = tol / j
    an = EXP(c)
    IF (an > tol) THEN
      IF (an >= 1.0-tol) GO TO 160

      c = an
      sum = 0.0
      30 i = i + 1.0
      IF (i < j) THEN
        c = ((am1+i)/i) * y * c
        sum = sum + c
        IF (c > tol) GO TO 30
      END IF

      i = n
      c = an
      40 i = i - 1.0
      IF (i /= 0.0) THEN
        c = i * c / ((i+am1)*y)
        sum = sum + c
        IF (c > tol) GO TO 40
      END IF
      p = an + sum
    END IF
  END IF

  IF (p >= 1.0) p = 1.0
  IF (ind == 0) RETURN
  p = 0.5 + (0.5-p)
  IF (p < 0.0) p = 0.0
  RETURN

!                SELECTION OF THE APPROPRIATE ALGORITHM

  50 am1 = a - 1.0
  IF (a <= 70.0) THEN
    IF (afrac == 0.0) GO TO 10

!                COMPUTE P0 = IX(A,1/2) OR P0 = IX(1/2,B)
!                           USING FORMULA 22

    temp = SQRT(x)
    rty = SQRT(y)
    c = ATAN(temp/rty) / pihalf
    IF (k == 0.0) GO TO 80
    ind = j
    m = k + k
    temp = -temp

    60 i = 0.0
    t = 1.0
    sum = 0.0
    70 i = i + 2.0
    IF (i /= m) THEN
      t = x * (i/(i+1.0)) * t
      sum = t + sum
      GO TO 70
    END IF

    p0 = (sum+1.0) * temp * rty / pihalf + c
    IF (ind /= 0) GO TO 110
    p = p0
    RETURN

    80 IF (j == 0.0) GO TO 140
    m = j + j
    x = y
    GO TO 60
  END IF

!                  COMPUTE P0 = IX(A,1/2) FOR A > 70
!                       USING EXPANSION 52 OR 53

  p0 = 0.0
  IF (x >= 0.7) THEN
    t = tol ** (1.0/am1)
    IF (x > t) THEN

      t = 0.5 + (0.5-t)
      lambda = SQRT(t)
      rty = SQRT(y)
      gamrat = rpinv * EXP(-algdiv(0.5,a))
      IF (t < 4.0*y) THEN

        c = lambda - rty
        temp = 2.0 * rty
        sum = 0.0
        DO l = 1, 10
          t = c * z(l)
          sum = sum + w(l) * (x-t*(t+temp)) ** am1
        END DO
        p0 = c * gamrat * sum + 0.5 * tol
      ELSE

        sum = 0.0
        DO l = 1, 10
          t = 1.0 - y * z(l) * z(l)
          sum = sum + w(l) * t ** am1
        END DO
        p0 = 1.0 - rty * gamrat * sum
      END IF
    END IF
  END IF

!                     COMPUTE P USING EXPANSION 21

  110 IF (j /= 0.0) THEN
    n = j
    IF (y < 2.0*j*x) THEN
      t = am1 * y / x + 0.5
      IF (t < 2.0) THEN
        n = 1.0
      ELSE
        IF (t < j) n = INT(t)
      END IF
    END IF

    t = n - 0.5
    c = (a*LOG(x)+t*alnrel(-x)) - blnd(a,t)
    IF (c > -30.0) THEN
      c = EXP(c)
      IF (c > tol/j) THEN
        IF (p0+c >= 1.0-tol) GO TO 150

        tol = tol / j
        lambda = c
        sum = 0.0
        120 t = t + 1.0
        IF (t <= j) THEN
          lambda = (am1+t) * y * lambda / t
          sum = sum + lambda
          IF (lambda > tol) GO TO 120
        END IF

        lambda = c
        t = a - 0.5
        i = n
        130 i = i - 1.0
        IF (i > 0.0) THEN
          lambda = ((i+0.5)/(i+t)) * lambda / y
          sum = lambda + sum
          IF (lambda > tol) GO TO 130
        END IF

        p = c + sum
      END IF
    END IF
  END IF
  p = p + p0
  IF (p >= 1.0) p = 1.0
  RETURN
END IF

!                           SPECIAL CASES

p = x
RETURN

140 p = c
RETURN

150 p = 1.0
RETURN

160 p = 1 - ind
RETURN

!                           ERROR RETURN

170 ierr = 2
RETURN

180 IF (a <= 0.0.OR.b <= 0.0) GO TO 170
190 ierr = 3
RETURN

200 ierr = 4
RETURN
END SUBROUTINE isubx


FUNCTION blnd(a, b) RESULT(fn_val)
REAL, INTENT(IN) :: a, b
REAL             :: fn_val

IF (a <= 20.0) THEN
  fn_val = (logam(a) - logam(a+b)) + logam(b+1.0)
  RETURN
END IF
fn_val = algdiv(b,a) + logam(b+1.0)
RETURN
END FUNCTION blnd


FUNCTION logam(x) RESULT(fn_val)
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: w(200) = (/ .57236494292470E+00, 0.0, -.12078223763525E+00, 0.0,  &
  .28468287047292E+00, .69314718055995E+00, .12009736023471E+01, .17917594692281E+01, .24537365708424E+01, .31780538303479E+01, &
  .39578139676187E+01, .47874917427820E+01, .56625620598571E+01, .65792512120101E+01, .75343642367587E+01, .85251613610654E+01, &
  .95492672573010E+01, .10604602902745E+02, .11689333420797E+02, .12801827480081E+02, .13940625219404E+02, .15104412573076E+02, &
  .16292000476567E+02, .17502307845874E+02, .18734347511936E+02, .19987214495662E+02, .21260076156245E+02, .22552163853123E+02, &
  .23862765841689E+02, .25191221182739E+02, .26536914491116E+02, .27899271383841E+02, .29277754515041E+02, .30671860106081E+02, &
  .32081114895947E+02, .33505073450137E+02, .34943315776877E+02, .36395445208033E+02, .37861086508961E+02, .39339884187199E+02, &
  .40831500974531E+02, .42335616460753E+02, .43851925860675E+02, .45380138898477E+02, .46919978795809E+02, .48471181351835E+02, &
  .50033494105019E+02, .51606675567764E+02, .53190494526169E+02, .54784729398112E+02, .56389167643720E+02, .58003605222981E+02, &
  .59627846095884E+02, .61261701761002E+02, .62904990828877E+02, .64557538627006E+02, .66219176833549E+02, .67889743137182E+02, &
  .69569080920824E+02, .71257038967168E+02, .72953471184169E+02, .74658236348830E+02, .76371197867783E+02, .78092223553315E+02, &
  .79821185413614E+02, .81557959456115E+02, .83302425502950E+02, .85054467017582E+02, .86813970941781E+02, .88580827542198E+02, &
  .90354930265818E+02, .92136175603687E+02, .93924462962300E+02, .95719694542143E+02, .97521775222888E+02, .99330612454787E+02, &
  .10114611615586E+03, .10296819861451E+03, .10479677439716E+03, .10663176026064E+03, .10847307506907E+03, .11032063971476E+03, &
  .11217437704318E+03, .11403421178146E+03, .11590007047041E+03, .11777188139975E+03, .11964957454634E+03, .12153308151544E+03, &
  .12342233548444E+03, .12531727114936E+03, .12721782467361E+03, .12912393363913E+03, .13103553699957E+03, .13295257503562E+03, &
  .13487498931216E+03, .13680272263733E+03, .13873571902320E+03, .14067392364823E+03, .14261728282115E+03, .14456574394634E+03, &
  .14651925549072E+03, .14847776695177E+03, .15044122882700E+03, .15240959258450E+03, .15438281063467E+03, .15636083630308E+03, &
  .15834362380427E+03, .16033112821663E+03, .16232330545817E+03, .16432011226320E+03, .16632150615984E+03, .16832744544843E+03, &
  .17033788918059E+03, .17235279713916E+03, .17437212981875E+03, .17639584840700E+03, .17842391476655E+03, .18045629141754E+03, &
  .18249294152079E+03, .18453382886145E+03, .18657891783334E+03, .18862817342367E+03, .19068156119837E+03, .19273904728784E+03, &
  .19480059837319E+03, .19686618167289E+03, .19893576492993E+03, .20100931639928E+03, .20308680483583E+03, .20516819948264E+03, &
  .20725347005963E+03, .20934258675254E+03, .21143552020227E+03, .21353224149456E+03, .21563272214993E+03, .21773693411395E+03, &
  .21984484974781E+03, .22195644181913E+03, .22407168349308E+03, .22619054832373E+03, .22831301024565E+03, .23043904356578E+03, &
  .23256862295547E+03, .23470172344282E+03, .23683832040517E+03, .23897838956183E+03, .24112190696703E+03, .24326884900298E+03, &
  .24541919237325E+03, .24757291409619E+03, .24972999149863E+03, .25189040220972E+03, .25405412415489E+03, .25622113555001E+03, &
  .25839141489572E+03, .26056494097186E+03, .26274169283208E+03, .26492164979855E+03, .26710479145687E+03, .26929109765102E+03, &
  .27148054847853E+03, .27367312428569E+03, .27586880566295E+03, .27806757344037E+03, .28026940868320E+03, .28247429268763E+03, &
  .28468220697654E+03, .28689313329543E+03, .28910705360840E+03, .29132395009427E+03, .29354380514276E+03, .29576660135076E+03, &
  .29799232151870E+03, .30022094864701E+03, .30245246593264E+03, .30468685676567E+03, .30692410472600E+03, .30916419358015E+03, &
  .31140710727802E+03, .31365282994988E+03, .31590134590330E+03, .31815263962021E+03, .32040669575401E+03, .32266349912673E+03, &
  .32492303472629E+03, .32718528770378E+03, .32945024337081E+03, .33171788719693E+03, .33398820480710E+03, .33626118197920E+03, &
  .33853680464160E+03, .34081505887080E+03, .34309593088909E+03, .34537940706227E+03, .34766547389743E+03, .34995411804077E+03, &
  .35224532627544E+03, .35453908551944E+03, .35683538282361E+03, .35913420536958E+03 /)
!     ------------------------------------------------------------------
!     D = 0.5*(LN(2*PI) - 1)
!     ------------------------------------------------------------------
REAL, PARAMETER :: d = .41893853320467
REAL            :: t, z
INTEGER         :: n
!     ------------------------------------------------------------------
!     COMPUTATION OF LN(GAMMA(X)) FOR X = N/2  WHERE N IS AN INTEGER
!     ------------------------------------------------------------------
IF (x <= 100.0) THEN
  n = 2.0 * x + 0.1
  fn_val = w(n)
  RETURN
END IF

t = (1.0/x) ** 2
z = (((-0.75*t + 1.0)*t - 3.5)*t + 105.0) / (x*1260.0)
fn_val = (d + z) + (x - 0.5) * (LOG(x)-1.0)
RETURN
END FUNCTION logam


SUBROUTINE cbsslj(z, cnu, w)
!-----------------------------------------------------------------------

!         EVALUATION OF THE COMPLEX BESSEL FUNCTION J   (Z)
!                                                    CNU
!-----------------------------------------------------------------------

!     WRITTEN BY
!         ANDREW H. VAN TUYL AND ALFRED H. MORRIS, JR.
!         NAVAL SURFACE WARFARE CENTER
!         OCTOBER, 1991

!     A MODIFICATION OF THE PROCEDURE DEVELOPED BY ALLEN V. HERSHEY
!     (NAVAL SURFACE WARFARE CENTER) IN 1978 FOR HANDLING THE DEBYE
!     APPROXIMATION IS EMPLOYED.

!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, cnu
COMPLEX, INTENT(OUT) :: w

! Local variables
COMPLEX         :: c, nu, s, sm1, sm2, t, tsc, w0, w1, zn, zz
REAL, PARAMETER :: pi = 3.14159265358979
REAL            :: a, cn1, cn2, e, fn, pn, qm, qn, qnp1, r, rn2, r2, sn,  &
                   t1, t2, u, v, x, y
INTEGER         :: i, k, m, n
!-----------------------
x = REAL(z)
y = AIMAG(z)
r = cpabs(x,y)
cn1 = REAL(cnu)
cn2 = AIMAG(cnu)
rn2 = cn1 * cn1 + cn2 * cn2
pn = AINT(cn1)
fn = cn1 - pn
sn = 1.0

!          CALCULATION WHEN ORDER IS AN INTEGER

IF (fn == 0.0 .AND. cn2 == 0.0) THEN
  n = pn
  pn = ABS(pn)
  cn1 = pn
  IF (n < 0 .AND. n /= (n/2)*2) sn = -1.0
END IF

!          SELECTION OF METHOD

IF (r > 17.5) THEN
  IF (r > 17.5+0.5*rn2) GO TO 10
  GO TO 20
END IF

!          USE MACLAURIN EXPANSION AND RECURSION

IF (cn1 < 0.0) THEN
  qn = -1.25 * (r+0.5*ABS(cn2)-ABS(y-0.5*cn2))
  IF (cn1 < qn) THEN
    qn = 1.25 * (r-MAX(1.2*r,ABS(y-cn2)))
    IF (cn1 < qn) THEN
      qn = MIN(pn,-AINT(1.25*(r-ABS(cn2))))
      GO TO 60
    END IF
  END IF
END IF

r2 = r * r
qm = 0.0625 * r2 * r2 - cn2 * cn2
qn = MAX(pn,AINT(SQRT(MAX(0.0,qm))))
GO TO 60

!          USE ASYMPTOTIC EXPANSION

10 CALL cbja(z,cnu,w)
RETURN

!          CALCULATION FOR 17.5 < ABS(Z) <= 17.5 + 0.5*ABS(CNU)**2

20 n = 0
IF (ABS(cn2) < 0.8*ABS(y)) THEN
  qm = -1.25 * (r+0.5*ABS(cn2)-ABS(y-0.5*cn2))
  IF (cn1 < qm) THEN
    qm = 1.25 * (r-MAX(1.2*r,ABS(y-cn2)))
    IF (cn1 < qm) n = 1
  END IF
END IF

qn = pn
a = 4.e-3 * r * r
zz = z
IF (x < 0.0) zz = -z

!          CALCULATION OF ZONE OF EXCLUSION OF DEBYE APPROXIMATION

30 nu = CMPLX(qn+fn,cn2)
zn = nu / z
t2 = AIMAG(zn) * AIMAG(zn)
u = 1.0 - REAL(zn)
t1 = u * u + t2
u = 1.0 + REAL(zn)
t2 = u * u + t2
u = t1 * t2
v = a * u / (t1*t1+t2*t2)
IF (u*v*v <= 1.0) THEN

!          THE ARGUMENT LIES INSIDE THE ZONE OF EXCLUSION

  qn = qn + 1.0
  IF (n == 0) GO TO 30

!          USE MACLAURIN EXPANSION WITH FORWARD RECURRENCE

  qn = MIN(pn,-AINT(1.25*(r-ABS(cn2))))
ELSE

!          USE BACKWARD RECURRENCE STARTING FROM THE
!          ASYMPTOTIC EXPANSION

  qnp1 = qn + 1.0
  IF (ABS(qn) < ABS(pn)) THEN
    IF (r >= 17.5+0.5*(qnp1*qnp1+cn2*cn2)) THEN

      nu = CMPLX(qn+fn,cn2)
      CALL cbja(zz,nu,sm1)
      nu = CMPLX(qnp1+fn,cn2)
      CALL cbja(zz,nu,sm2)
      GO TO 40
    END IF
  END IF

!          USE BACKWARD RECURRENCE STARTING FROM THE DEBYE APPROXIMATION

  nu = CMPLX(qn+fn,cn2)
  CALL cbdb(zz,nu,fn,sm1)
  IF (qn == pn) GO TO 50
  nu = CMPLX(qnp1+fn,cn2)
  CALL cbdb(zz,nu,fn,sm2)

  40 nu = CMPLX(qn+fn,cn2)
  tsc = 2.0 * nu * sm1 / zz - sm2
  sm2 = sm1
  sm1 = tsc
  qn = qn - 1.0
  IF (qn /= pn) GO TO 40

  50 w = sm1
  IF (sn < 0.0) w = -w
  IF (x >= 0.0) RETURN

  nu = pi * CMPLX(-cn2,cn1)
  IF (y < 0.0) nu = -nu
  w = EXP(nu) * w
  RETURN
END IF

!          USE MACLAURIN EXPANSION WITH FORWARD OR BACKWARD RECURRENCE.

60 m = qn - pn
IF (ABS(m) <= 1) THEN
  nu = CMPLX(cn1,cn2)
  CALL cbjm(z,nu,w)
ELSE
  nu = CMPLX(qn+fn,cn2)
  CALL cbjm(z,nu,w1)
  w0 = 0.25 * z * z
  IF (m <= 0) THEN

!          FORWARD RECURRENCE

    m = ABS(m)
    nu = nu + 1.0
    CALL cbjm(z,nu,w)
    DO i = 2, m
      c = nu * (nu+1.0)
      t = (c/w0) * (w-w1)
      w1 = w
      w = t
      nu = nu + 1.0
    END DO
  ELSE

!          BACKWARD RECURRENCE

    nu = nu - 1.0
    CALL cbjm(z,nu,w)
    DO i = 2, m
      c = nu * (nu+1.0)
      t = (w0/c) * w1
      w1 = w
      w = w - t
      nu = nu - 1.0
    END DO
  END IF
END IF

!          FINAL ASSEMBLY

IF (fn == 0.0 .AND. cn2 == 0.0) THEN
  k = pn
  IF (k == 0) RETURN
  e = sn / gamma(pn+1.0)
  w = e * w * (0.5*z) ** k
  RETURN
END IF

s = cnu * LOG(0.5*z)
w = EXP(s) * w
IF (rn2 <= 0.81) THEN
  w = w * cgam0(cnu)
  RETURN
END IF
CALL cgamma(0, cnu, t)
w = cdiv(w, cnu*t)
RETURN
END SUBROUTINE cbsslj


SUBROUTINE cbjm(z, cnu, w)
!-----------------------------------------------------------------------

!       COMPUTATION OF  (Z/2)**(-CNU) * GAMMA(CNU + 1) * J(CNU,Z)

!                           -----------------

!     THE MACLAURIN EXPANSION IS USED. IT IS ASSUMED THAT CNU IS NOT
!     A NEGATIVE INTEGER.

!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, cnu
COMPLEX, INTENT(OUT) :: w

! Local variables
COMPLEX :: nu, nup1, p, s, sn, t, ti
REAL    :: a, a0, eps, inu, m, rnu
INTEGER :: i, imin, k, km1, km2
!--------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0 .

eps = EPSILON(1.0)

!--------------------------
s = -0.25 * (z*z)
nu = cnu
rnu = REAL(nu)
inu = AIMAG(nu)
a = 0.5 + (0.5+rnu)
nup1 = CMPLX(a,inu)

IF (a > 0.0) THEN
  m = 1.0
  t = s / nup1
  w = 1.0 + t
ELSE

!     ADD 1.0 AND THE FIRST K-1 TERMS

  k = INT(-a) + 2
  km1 = k - 1
  w = (1.0,0.0)
  t = w
  DO i = 1, km1
    m = i
    t = t * (s/(m*(nu+m)))
    w = w + t
    IF (anorm(t) <= eps*anorm(w)) GO TO 20
  END DO
  GO TO 50

!     CHECK IF THE (K-1)-ST AND K-TH TERMS CAN BE IGNORED.
!     IF SO THEN THE SUMMATION IS COMPLETE.

  20 IF (i /= km1) THEN
    imin = i + 1
    IF (imin < k-5) THEN
      ti = t

      m = km1
      t = s / (nu+m)
      a0 = anorm(t) / m
      t = t * (s/(nu+(m+1.0)))
      a = anorm(t) / (m*(m+1.0))
      a = MAX(a,a0)

      t = (1.0,0.0)
      km2 = k - 2
      DO i = imin, km2
        m = i
        t = t * (s/(m*(nu+m)))
        IF (a*anorm(t) < 0.5) RETURN
      END DO
      t = t * ti
      imin = km2
    END IF

!     ADD THE (K-1)-ST TERM

    a = 1.0
    p = (1.0,0.0)
    sn = p
    DO i = imin, km1
      m = i
      a = a * m
      p = p * (nu+m)
      sn = s * sn
    END DO
    t = t * (cdiv(sn,p)/a)
    w = w + t
  END IF
END IF

!     ADD THE REMAINING TERMS

50 m = m + 1.0
t = t * (s/(m*(nu+m)))
w = w + t
IF (anorm(t) > eps*anorm(w)) GO TO 50

RETURN
END SUBROUTINE cbjm


SUBROUTINE cbdb(cz, cnu, fn, w)
!-----------------------------------------------------------------------

!         CALCULATION OF J   (CZ) BY THE DEBYE APPROXIMATION
!                         CNU
!                         ------------------

!     IT IS ASSUMED THAT REAL(CZ) >= 0 AND THAT REAL(CNU) = FN + K
!     WHERE K IS AN INTEGER.

!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: cz, cnu
REAL, INTENT(IN)     :: fn
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL :: a(136) = (/ 1.0,-.208333333333333E+00, .125000000000000E+00, .334201388888889E+00, &
  -.401041666666667E+00, .703125000000000E-01,-.102581259645062E+01, .184646267361111E+01, &
  -.891210937500000E+00, .732421875000000E-01, .466958442342625E+01,-.112070026162230E+02, &
   .878912353515625E+01,-.236408691406250E+01, .112152099609375E+00,-.282120725582002E+02, &
   .846362176746007E+02,-.918182415432400E+02, .425349987453885E+02,-.736879435947963E+01, &
   .227108001708984E+00, .212570130039217E+03,-.765252468141182E+03, .105999045252800E+04, &
  -.699579627376133E+03, .218190511744212E+03,-.264914304869516E+02, .572501420974731E+00, &
  -.191945766231841E+04, .806172218173731E+04,-.135865500064341E+05, .116553933368645E+05, &
  -.530564697861340E+04, .120090291321635E+04,-.108090919788395E+03, .172772750258446E+01, &
   .202042913309661E+05,-.969805983886375E+05, .192547001232532E+06,-.203400177280416E+06, &
   .122200464983017E+06,-.411926549688976E+05, .710951430248936E+04,-.493915304773088E+03, &
   .607404200127348E+01,-.242919187900551E+06, .131176361466298E+07,-.299801591853811E+07, &
   .376327129765640E+07,-.281356322658653E+07, .126836527332162E+07,-.331645172484564E+06, &
   .452187689813627E+05,-.249983048181121E+04, .243805296995561E+02, .328446985307204E+07, &
  -.197068191184322E+08, .509526024926646E+08,-.741051482115327E+08, .663445122747290E+08, &
  -.375671766607634E+08, .132887671664218E+08,-.278561812808645E+07, .308186404612662E+06, &
  -.138860897537170E+05, .110017140269247E+03,-.493292536645100E+08, .325573074185766E+09, &
  -.939462359681578E+09, .155359689957058E+10,-.162108055210834E+10, .110684281682301E+10, &
  -.495889784275030E+09, .142062907797533E+09,-.244740627257387E+08, .224376817792245E+07, &
  -.840054336030241E+05, .551335896122021E+03, .814789096118312E+09,-.586648149205185E+10, &
   .186882075092958E+11,-.346320433881588E+11, .412801855797540E+11,-.330265997498007E+11, &
   .179542137311556E+11,-.656329379261928E+10, .155927986487926E+10,-.225105661889415E+09, &
   .173951075539782E+08,-.549842327572289E+06, .303809051092238E+04,-.146792612476956E+11, &
   .114498237732026E+12,-.399096175224466E+12, .819218669548577E+12,-.109837515608122E+13, &
   .100815810686538E+13,-.645364869245377E+12, .287900649906151E+12,-.878670721780233E+11, &
   .176347306068350E+11,-.216716498322380E+10, .143157876718889E+09,-.387183344257261E+07, &
   .182577554742932E+05, .286464035717679E+12,-.240629790002850E+13, .910934118523990E+13, &
  -.205168994109344E+14, .305651255199353E+14,-.316670885847852E+14, .233483640445818E+14, &
  -.123204913055983E+14, .461272578084913E+13,-.119655288019618E+13, .205914503232410E+12, &
  -.218229277575292E+11, .124700929351271E+10,-.291883881222208E+08, .118838426256783E+06, &
  -.601972341723401E+13, .541775107551060E+14,-.221349638702525E+15, .542739664987660E+15, &
  -.889496939881026E+15, .102695519608276E+16,-.857461032982895E+15, .523054882578445E+15, &
  -.232604831188940E+15, .743731229086791E+14,-.166348247248925E+14, .248500092803409E+13, &
  -.229619372968246E+12, .114657548994482E+11,-.234557963522252E+09, .832859304016289E+06 /)
REAL               :: alpha, am, aq, ar, is, inu, izn, phi, sgn, theta, u,  &
                      v, x, y
REAL, PARAMETER    :: c = .398942280401433, pi = 3.14159265358979,  &
                      pi2 = 6.28318530717959, bnd = 1.04719755119660
COMPLEX            :: c1, c2, eta, nu, p, p1, q, r, s, s1, s2, sm, t, z, zn
COMPLEX, PARAMETER :: j = (0.0, 1.0)
INTEGER            :: ind, k, l, m
!----------------------
!     C = 1/SQRT(2*PI)
!     BND = PI/3
!----------------------
z = cz
nu = cnu
inu = AIMAG(cnu)
IF (inu < 0.0) THEN
  z = CONJG(z)
  nu = CONJG(nu)
END IF
x = REAL(z)
y = AIMAG(z)

!          TANH(GAMMA) = SQRT(1 - (Z/NU)**2) = W/NU
!          T = EXP(NU*(TANH(GAMMA) - GAMMA))

zn = z / nu
izn = AIMAG(zn)
IF (ABS(izn) <= 0.1*ABS(REAL(zn))) THEN

  s = (1.0-zn) * (1.0+zn)
  eta = 1.0 / s
  q = SQRT(s)
  s = 1.0 / (nu*q)
  t = zn / (1.0 + q)
  t = EXP(nu*(q + LOG(t)))
ELSE

  s = (nu-z) * (nu+z)
  eta = (nu*nu) / s
  w = SQRT(s)
  q = w / nu
  IF (REAL(q) < 0.0) w = -w
  s = 1.0 / w
  t = z / (nu+w)
  t = EXP(w + nu*LOG(t))
END IF

is = AIMAG(s)
r = SQRT(s)
c1 = r * t
ar = REAL(r) * REAL(r) + AIMAG(r) * AIMAG(r)
aq = -1.0 / (REAL(q)*REAL(q) + AIMAG(q)*AIMAG(q))

phi = ATAN2(y,x) / 3.0
q = nu - z
theta = ATAN2(AIMAG(q),REAL(q)) - phi
ind = 0
IF (ABS(theta) >= 2.0*bnd) THEN

  ind = 1
  CALL crec(REAL(t), AIMAG(t), u, v)
  c2 = -j * r * CMPLX(u,v)
  IF (is >= 0.0) THEN
    IF (is > 0.0) GO TO 10
    IF (REAL(s) <= 0.0) GO TO 10
  END IF
  c2 = -c2
END IF

!          SUMMATION OF THE SERIES S1 AND S2

10 sm = s * s
p = (a(2)*eta + a(3)) * s
p1 = ((a(4)*eta + a(5))*eta + a(6)) * sm
s1 = (1.0 + p) + p1
IF (ind /= 0) s2 = (1.0-p) + p1
sgn = 1.0
am = ar * ar
m = 4
l = 6

!          P = VALUE OF THE M-TH POLYNOMIAL

20 l = l + 1
alpha = a(l)
p = CMPLX(a(l), 0.0)
DO k = 2, m
  l = l + 1
  alpha = a(l) + aq * alpha
  p = a(l) + eta * p
END DO

!          ONLY THE S1 SUM IS FORMED WHEN IND = 0

sm = s * sm
p = p * sm
s1 = s1 + p
IF (ind /= 0) THEN
  sgn = -sgn
  s2 = s2 + sgn * p
END IF
am = ar * am
IF (1.0 + alpha*am /= 1.0) THEN
  m = m + 1
  IF (m <= 16) GO TO 20
END IF

!          FINAL ASSEMBLY

s1 = c * c1 * s1
IF (ind == 0) THEN
  w = s1
ELSE

  s2 = c * c2 * s2
  q = nu + z
  theta = ATAN2(AIMAG(q),REAL(q)) - phi
  IF (ABS(theta) <= bnd) THEN
    w = s1 + s2
  ELSE

    alpha = pi2
    IF (izn < 0.0) alpha = -alpha
    t = alpha * CMPLX(ABS(inu), -fn)
    alpha = EXP(REAL(t))
    u = AIMAG(t)
    r = CMPLX(COS(u),SIN(u))
    t = s1 - (alpha*r) * s1
    IF (x == 0.0 .AND. inu == 0.0) t = -t

    IF (y < 0.0) THEN
      IF (izn >= 0.0 .AND. theta <= SIGN(pi,theta)) s2 = s2 * (  &
      CONJG(r)/alpha)
      IF (x == 0.0) GO TO 40
      IF (izn >= 0.0) THEN
        IF (is < 0.0) GO TO 40
      END IF
    END IF

    w = s2 + t
    GO TO 50
    40 w = s2 - t
  END IF
END IF

50 IF (inu < 0.0) w = CONJG(w)
RETURN
END SUBROUTINE cbdb


SUBROUTINE cbja(cz, cnu, w)
!-----------------------------------------------------------------------
!        COMPUTATION OF J(NU,Z) BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: cz, cnu
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL               :: eps, inu, m, r, rnu, tol, u, v, x, y
REAL, PARAMETER    :: pihalf = 1.5707963267949, c = 1.12837916709551
COMPLEX            :: a, a1, arg, e, eta, nu, p, q, t, z, zr, zz
COMPLEX, PARAMETER :: j = (0.0,1.0)
INTEGER            :: i, ind
!--------------------------
!     PIHALF = PI/2
!     C = 2*PI**(-1/2)
!--------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0 .

eps = EPSILON(1.0)

!--------------------------
z = cz
x = REAL(z)
y = AIMAG(z)
nu = cnu
ind = 0
IF (ABS(x) <= 1.e-2*ABS(y)) THEN
  IF (AIMAG(nu) < 0.0 .AND. ABS(REAL(nu)) < 1.e-2*ABS(AIMAG(nu))) THEN
    ind = 1
    nu = CONJG(nu)
    z = CONJG(z)
    y = -y
  END IF
END IF

IF (x < -1.e-2*y) z = -z
zz = z + z
CALL crec(REAL(zz),AIMAG(zz),u,v)
zr = CMPLX(u,v)
eta = -zr * zr

p = (0.0,0.0)
q = (0.0,0.0)
a1 = nu * nu - 0.25
a = a1
t = a1
m = 1.0
tol = eps * anorm(a1)
DO i = 1, 16
  a = a - 2.0 * m
  m = m + 1.0
  t = t * a * eta / m
  p = p + t
  a = a - 2.0 * m
  m = m + 1.0
  t = t * a / m
  q = q + t
  IF (anorm(t) <= tol) GO TO 20
END DO

20 p = p + 1.0
q = (q+a1) * zr
w = z - pihalf * nu
IF (ABS(AIMAG(w)) <= 1.0) THEN
  arg = w - 0.5 * pihalf
  w = c * SQRT(zr) * (p*COS(arg)-q*SIN(arg))
ELSE
  e = EXP(-j*w)
  t = q - j * p
  IF (AIMAG(z) > 0.0 .AND. REAL(z) <= 1.e-2*AIMAG(z) .AND.   &
  ABS(REAL(nu)) < 1.e-2*AIMAG(nu)) t = 0.5 * t
  CALL crec(REAL(e),AIMAG(e),u,v)
  w = 0.5 * c * SQRT(j*zr) * ((p-j*q)*e+t*CMPLX(u,v))
END IF

IF (x < -1.e-2*y) THEN
  IF (y < 0.0) nu = -nu

!     COMPUTATION OF EXP(I*PI*NU)

  rnu = REAL(nu)
  inu = AIMAG(nu)
  r = EXP(-2.0*pihalf*inu)
  u = r * cos1(rnu)
  v = r * sin1(rnu)
  w = w * CMPLX(u,v)
END IF

IF (ind /= 0) w = CONJG(w)
RETURN
END SUBROUTINE cbja


SUBROUTINE bsslj(a, in, w)
!     ******************************************************************
!     FORTRAN SUBROUTINE FOR ORDINARY BESSEL FUNCTION OF INTEGRAL ORDER
!     ******************************************************************
!     A  = ARGUMENT (COMPLEX NUMBER)
!     IN = ORDER (INTEGER)
!     W  = FUNCTION OF FIRST KIND (COMPLEX NUMBER)
!     -------------------
COMPLEX, INTENT(IN)  :: a
INTEGER, INTENT(IN)  :: in
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: cd(30) = (/ 0.00000000000000,  &
        -1.64899505142212E-2, -7.18621880068536E-2, -1.67086878124866E-1,  &
        -3.02582250219469E-1, -4.80613945245927E-1, -7.07075239357898E-1,  &
        -9.92995790539516E-1, -1.35583925612592, -1.82105907899132,  &
        -2.42482175310879, -3.21956655708750, -4.28658077248384,  &
        -5.77022816798128, -8.01371260952526, 0.00000000000000,  &
        -5.57742429879505E-3, -4.99112944172476E-2, -1.37440911652397E-1,  &
        -2.67233784710566E-1, -4.40380166808682E-1, -6.61813614872541E-1,  &
        -9.41861077665017E-1, -1.29754130468326, -1.75407696719816,  &
        -2.34755299882276, -3.13041332689196, -4.18397120563729,  &
        -5.65251799214994, -7.87863959810677 /), ce(30)  &
        = (/ 0.00000000000000, -4.80942336387447E-3,  &
        -1.31366200347759E-2, -1.94843834008458E-2, -2.19948900032003E-2,  &
        -2.09396625676519E-2, -1.74600268458650E-2, -1.27937813362085E-2,  &
        -8.05234421796592E-3, -4.15817375002760E-3, -1.64317738747922E-3,  &
        -4.49175585314709E-4, -7.28594765574007E-5, -5.38265230658285E-6,  &
        -9.93779048036289E-8, 0.00000000000000, 7.53805779200591E-2,  &
        7.12293537403464E-2, 6.33116224228200E-2, 5.28240264523301E-2,  &
        4.13305359441492E-2, 3.01350573947510E-2, 2.01043439592720E-2,  &
        1.18552223068074E-2, 5.86055510956010E-3, 2.25465148267325E-3,  &
        6.08173041536336E-4, 9.84215550625747E-5, 7.32139093038089E-6,  &
        1.37279667384666E-7 /)
REAL    :: an, aq(2), az(2), fj(2), pm, pn, qf(2), qm, qn, qz(2), rz(2),  &
           sz(2), rm(4), sm(4), sn, ss, ts(2), tm(2), zm, zr(2), zs
INTEGER :: i, m, n
!     -------------------
az(1) = REAL(a)
az(2) = AIMAG(a)
zs = az(1) * az(1) + az(2) * az(2)
zm = SQRT(zs)
pn = ABS(in)
sn = 1.0
IF (in < 0) THEN
  IF (in /= in/2*2) THEN
    sn = -1.0
  END IF
END IF
IF (az(1) < 0) THEN
  qz(1) = -az(1)
  qz(2) = -az(2)
  IF (in == in/2*2) GO TO 10
  sn = -sn
ELSE
  qz(1) = +az(1)
  qz(2) = +az(2)
END IF
10 IF (zm > 17.5+0.5*pn*pn) THEN
  qn = pn
ELSE
  qn = 0.5 * zm - 0.5 * ABS(qz(2)) + 0.5 *ABS(0.5*zm-ABS(qz(2)))
  IF (pn > qn) THEN
    qn = +AINT(0.0625*zs)
    IF (pn <= qn) GO TO 130
    qn = pn
    GO TO 130
  END IF
  IF (zm > 17.5) THEN
    qn = +AINT(SQRT(2.0*(zm-17.5)))
  ELSE
    IF (zs >= 1.0) THEN
      IF (-ABS(az(2))+0.096*az(1)*az(1) >= 0) GO TO 20
    END IF
    qn = +AINT(0.0625*zs)
    IF (pn <= qn) GO TO 130
    qn = pn
    GO TO 130
    20 qn = 0.0
  END IF
END IF
sz(1) = qz(1)
sz(2) = qz(2)
qm = sn * 0.797884560802865
zr(1) = SQRT(sz(1)+zm)
zr(2) = sz(2) / zr(1)
zr(1) = 0.707106781186548 * zr(1)
zr(2) = 0.707106781186548 * zr(2)
qf(1) = +qm * zr(1) / zm
qf(2) = -qm * zr(2) / zm
IF (zm > 17.5) THEN
  rz(1) = +0.5 * qz(1) / zs
  rz(2) = -0.5 * qz(2) / zs
  an = qn * qn - 0.25
  sm(1) = 0.0
  sm(2) = 0.0
  sm(3) = 0.0
  sm(4) = 0.0
  tm(1) = 1.0
  tm(2) = 0.0
  pm = 0.0
  GO TO 40
  30 an = an - 2.0 * pm
  pm = pm + 1.0
  ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
  ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
  tm(1) = -an * ts(1) / pm
  tm(2) = -an * ts(2) / pm
  40 sm(1) = sm(1) + tm(1)
  sm(2) = sm(2) + tm(2)
  an = an - 2.0 * pm
  pm = pm + 1.0
  ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
  ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
  tm(1) = +an * ts(1) / pm
  tm(2) = +an * ts(2) / pm
  IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
    IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 60
  END IF
  sm(3) = sm(3) + tm(1)
  sm(4) = sm(4) + tm(2)
  IF (pm < 35.0) GO TO 30
ELSE
  sm(1) = 1.0
  sm(2) = 0.0
  sm(3) = 1.0
  sm(4) = 0.0
  m = 15.0 * qn + 2.0
  n = 15.0 * qn + 15.0
  DO i = m, n
    ts(1) = +qz(2) - cd(i)
    ts(2) = -qz(1)
    ss = ts(1) * ts(1) + ts(2) * ts(2)
    tm(1) = +ce(i) * ts(1) / ss
    tm(2) = -ce(i) * ts(2) / ss
    sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
    ts(1) = -qz(2) - cd(i)
    ts(2) = +qz(1)
    ss = ts(1) * ts(1) + ts(2) * ts(2)
    tm(1) = +ce(i) * ts(1) / ss
    tm(2) = -ce(i) * ts(2) / ss
    sm(3) = sm(3) + tm(1)
    sm(4) = sm(4) + tm(2)
  END DO
  ts(1) = +0.5 * (sm(2)-sm(4))
  ts(2) = -0.5 * (sm(1)-sm(3))
  sm(1) = +0.5 * (sm(1)+sm(3))
  sm(2) = +0.5 * (sm(2)+sm(4))
  sm(3) = ts(1)
  sm(4) = ts(2)
END IF
60 aq(1) = qz(1) - 1.57079632679490 * (qn+0.5)
aq(2) = qz(2)
ts(1) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
ts(2) = -SIN(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
tm(1) = sm(1) * ts(1) - sm(2) * ts(2)
tm(2) = sm(1) * ts(2) + sm(2) * ts(1)
ts(1) = +SIN(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
ts(2) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
rm(1) = tm(1) - sm(3) * ts(1) + sm(4) * ts(2)
rm(2) = tm(2) - sm(3) * ts(2) - sm(4) * ts(1)
IF (qn /= pn) THEN
  rm(3) = rm(1)
  rm(4) = rm(2)
  qn = qn + 1.0
  IF (zm > 17.5) THEN
    an = qn * qn - 0.25
    sm(1) = 0.0
    sm(2) = 0.0
    sm(3) = 0.0
    sm(4) = 0.0
    tm(1) = 1.0
    tm(2) = 0.0
    pm = 0.0
    GO TO 80
    70 an = an - 2.0 * pm
    pm = pm + 1.0
    ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
    ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
    tm(1) = -an * ts(1) / pm
    tm(2) = -an * ts(2) / pm
    80 sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
    an = an - 2.0 * pm
    pm = pm + 1.0
    ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
    ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
    tm(1) = +an * ts(1) / pm
    tm(2) = +an * ts(2) / pm
    IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
      IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 100
    END IF
    sm(3) = sm(3) + tm(1)
    sm(4) = sm(4) + tm(2)
    IF (pm < 35.0) GO TO 70
  ELSE
    sm(1) = 1.0
    sm(2) = 0.0
    sm(3) = 1.0
    sm(4) = 0.0
    m = 15.0 * qn + 2.0
    n = 15.0 * qn + 15.0
    DO i = m, n
      ts(1) = +qz(2) - cd(i)
      ts(2) = -qz(1)
      ss = ts(1) * ts(1) + ts(2) * ts(2)
      tm(1) = +ce(i) * ts(1) / ss
      tm(2) = -ce(i) * ts(2) / ss
      sm(1) = sm(1) + tm(1)
      sm(2) = sm(2) + tm(2)
      ts(1) = -qz(2) - cd(i)
      ts(2) = +qz(1)
      ss = ts(1) * ts(1) + ts(2) * ts(2)
      tm(1) = +ce(i) * ts(1) / ss
      tm(2) = -ce(i) * ts(2) / ss
      sm(3) = sm(3) + tm(1)
      sm(4) = sm(4) + tm(2)
    END DO
    ts(1) = +0.5 * (sm(2)-sm(4))
    ts(2) = -0.5 * (sm(1)-sm(3))
    sm(1) = +0.5 * (sm(1)+sm(3))
    sm(2) = +0.5 * (sm(2)+sm(4))
    sm(3) = ts(1)
    sm(4) = ts(2)
  END IF
  100 aq(1) = qz(1) - 1.57079632679490 * (qn+0.5)
  aq(2) = qz(2)
  ts(1) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
  ts(2) = -SIN(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
  tm(1) = sm(1) * ts(1) - sm(2) * ts(2)
  tm(2) = sm(1) * ts(2) + sm(2) * ts(1)
  ts(1) = +SIN(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
  ts(2) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
  rm(1) = tm(1) - sm(3) * ts(1) + sm(4) * ts(2)
  rm(2) = tm(2) - sm(3) * ts(2) - sm(4) * ts(1)
  GO TO 120
  110 tm(1) = +2.0 * qn * qz(1) / zs
  tm(2) = -2.0 * qn * qz(2) / zs
  ts(1) = tm(1) * rm(1) - tm(2) * rm(2) - rm(3)
  ts(2) = tm(1) * rm(2) + tm(2) * rm(1) - rm(4)
  rm(3) = rm(1)
  rm(4) = rm(2)
  rm(1) = ts(1)
  rm(2) = ts(2)
  qn = qn + 1.0
  120 IF (qn < pn) GO TO 110
END IF
fj(1) = qf(1) * rm(1) - qf(2) * rm(2)
fj(2) = qf(1) * rm(2) + qf(2) * rm(1)
w = CMPLX(fj(1), fj(2))
RETURN

130 sz(1) = 0.25 * (qz(1)*qz(1) - qz(2)*qz(2))
sz(2) = 0.5 * qz(1) * qz(2)
an = qn
sm(1) = 0.0
sm(2) = 0.0
sm(3) = 0.0
sm(4) = 0.0
tm(1) = 1.0
tm(2) = 0.0
pm = 0.0
140 an = an + 1.0
ts(1) = +tm(1) / an
ts(2) = +tm(2) / an
sm(3) = sm(3) + ts(1)
sm(4) = sm(4) + ts(2)
tm(1) = -ts(1) * sz(1) + ts(2) * sz(2)
tm(2) = -ts(1) * sz(2) - ts(2) * sz(1)
pm = pm + 1.0
tm(1) = tm(1) / pm
tm(2) = tm(2) / pm
IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 150
END IF
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
GO TO 140

150 sm(1) = sm(1) + 1.0
an = qn + 1.0
sm(3) = an * sm(3)
sm(4) = an * sm(4)
GO TO 170
160 an = qn * (qn+1.0)
tm(1) = sz(1) / an
tm(2) = sz(2) / an
ts(1) = -tm(1) * sm(3) + tm(2) * sm(4)
ts(2) = -tm(1) * sm(4) - tm(2) * sm(3)
sm(3) = sm(1)
sm(4) = sm(2)
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
qn = qn - 1.0
170 IF (qn > pn) GO TO 160
qf(1) = sn
qf(2) = 0.0
qn = 0.0
GO TO 190
180 qn = qn + 1.0
tm(1) = qf(1) * qz(1) - qf(2) * qz(2)
tm(2) = qf(1) * qz(2) + qf(2) * qz(1)
qf(1) = 0.5 * tm(1) / qn
qf(2) = 0.5 * tm(2) / qn
190 IF (qn < pn) GO TO 180
fj(1) = qf(1) * sm(1) - qf(2) * sm(2)
fj(2) = qf(1) * sm(2) + qf(2) * sm(1)
w = CMPLX(fj(1),fj(2))
RETURN
END SUBROUTINE bsslj



SUBROUTINE besj(x, alpha, n, y, nz)

!     WRITTEN BY D.E. AMOS, S.L. DANIEL AND M.K. WESTON, JANUARY, 1975.

!     REFERENCES
!         SAND-75-0147

!         CDC 6600 SUBROUTINES IBESS AND JBESS FOR BESSEL FUNCTIONS
!         I(NU,X) AND J(NU,X), X>=0, NU>=0  BY D.E. AMOS, S.L.
!         DANIEL, M.K. WESTON. ACM TRANS MATH SOFTWARE,3,PP 76-92 (1977)

!         TABLES OF BESSEL FUNCTIONS OF MODERATE OR LARGE ORDERS,
!         NPL MATHEMATICAL TABLES, VOL. 6, BY F.W.J. OLVER, HER
!         MAJESTY-S STATIONERY OFFICE, LONDON, 1962.

!     ABSTRACT
!         BESJ COMPUTES AN N MEMBER SEQUENCE OF J BESSEL FUNCTIONS
!         J/SUB(ALPHA+K-1)/(X), K=1,...,N FOR NON-NEGATIVE ALPHA AND X.
!         A COMBINATION OF THE POWER SERIES, THE ASYMPTOTIC EXPANSION
!         FOR X TO INFINITY AND THE UNIFORM ASYMPTOTIC EXPANSION FOR
!         NU TO INFINITY ARE APPLIED OVER SUBDIVISIONS OF THE (NU,X)
!         PLANE. FOR VALUES OF (NU,X) NOT COVERED BY ONE OF THESE
!         FORMULAE, THE ORDER IS INCREMENTED OR DECREMENTED BY INTEGER
!         VALUES INTO A REGION WHERE ONE OF THE FORMULAE APPLY. BACKWARD
!         RECURSION IS APPLIED TO REDUCE ORDERS BY INTEGER VALUES EXCEPT
!         WHERE THE ENTIRE SEQUENCE LIES IN THE OSCILLATORY REGION. IN
!         THIS CASE FORWARD RECURSION IS STABLE AND VALUES FROM THE
!         ASYMPTOTIC EXPANSION FOR X TO INFINITY START THE RECURSION
!         WHEN IT IS EFFICIENT TO DO SO. LEADING TERMS OF THE SERIES AND
!         UNIFORM EXPANSION ARE TESTED FOR UNDERFLOW. IF A SEQUENCE IS
!         REQUESTED AND THE LAST MEMBER WOULD UNDERFLOW, THE RESULT IS
!         SET TO ZERO AND THE NEXT LOWER ORDER TRIED, ETC., UNTIL A
!         MEMBER COMES ON SCALE OR ALL MEMBERS ARE SET TO ZERO. OVERFLOW
!         CANNOT OCCUR.

!         BESJ CALLS ASJY, JAIRY, GAMLN, SPMPAR, AND IPMPAR

!     DESCRIPTION OF ARGUMENTS

!         INPUT
!           X      - X >= 0.0
!           ALPHA  - ORDER OF FIRST MEMBER OF THE SEQUENCE,
!                    ALPHA >= 0.0
!           N      - NUMBER OF MEMBERS IN THE SEQUENCE, N >= 1

!         OUTPUT
!           Y      - A VECTOR WHOSE FIRST N COMPONENTS CONTAIN
!                    VALUES FOR J/SUB(ALPHA+K-1)/(X), K=1,...,N
!           NZ     - ERROR INDICATOR
!                    NZ=0      NORMAL RETURN - COMPUTATION COMPLETED
!                    NZ=-1     X IS LESS THAN 0.0
!                    NZ=-2     ALPHA IS LESS THAN 0.0
!                    NZ=-3     N IS LESS THAN 1
!                    NZ>0   LAST NZ COMPONENTS OF Y SET TO 0.0
!                              BECAUSE OF UNDERFLOW

!     ERROR CONDITIONS
!         IMPROPER INPUT ARGUMENTS - A FATAL ERROR
!         UNDERFLOW  - A NON-FATAL ERROR (NZ>0)

REAL, INTENT(IN)     :: x, alpha
INTEGER, INTENT(IN)  :: n
REAL, INTENT(OUT)    :: y(:)
INTEGER, INTENT(OUT) :: nz

! Local variables
INTEGER :: i, ialp, idalp, iflw, in, inlim = 150, is, i1, i2, k, kk, km,  &
           kt, nn, ns
REAL    :: ak, akm, ans, ap, arg, coef, dalpha, dfn, dtm, earg, elim, etx,  &
           fidal, flgjy, fn, fnf, fni, fnp1, fnu, gln, rden, relb,   &
           rtx, rzden, s, sa, sb, sxo2, s1, s2, t, ta, tau, tb, temp(3), tfn, tm,  &
           tol, tolln, trx, tx, t1, t2, wk(7), xo2, xo2l
REAL, PARAMETER :: fnulim(2) = (/ 100.0, 60.0 /),  &
                   pdf = 7.85398163397448E-01, pidt = 1.57079632679490E+00, &
                   pp(4) = (/ 8.72909153935547E+00, 2.65693932265030E-01, &
                   1.24578576865586E-01, 7.70133747430388E-04 /),   &
                   rttp = 7.97884560802865E-01, rtwo = 1.34839972492648E+00
!     -------------------
!     IPMPAR(8) REPLACES IPMPAR(5) IN A REAL (dp) CODE
!     IPMPAR(9) REPLACES IPMPAR(6) IN A REAL (dp) CODE

!     DEFINITION OF THE TOLERANCES TOL AND ELIM

tb = RADIX(1.0)
ta = EPSILON(1.0) / tb
IF (tb /= 2.0) THEN
  IF (tb == 8.0) GO TO 10
  IF (tb == 16.0) GO TO 20
  tb = LOG(tb)
  GO TO 30
END IF
tb = .69315
GO TO 30
10 tb = 2.07944
GO TO 30
20 tb = 2.77259

30 tol = MAX(ta,1.e-15)
i1 = DIGITS(1.0)
i2 = MINEXPONENT(1.0)
!     LN(10**3) = 6.90776
elim = REAL(-i2) * tb - 6.90776
!     TOLLN = -LN(TOL)
tolln = REAL(i1) * tb
tolln = MIN(tolln, 34.5388)

nz = 0
kt = 1
IF (n < 1) THEN
  GO TO 480
ELSE IF (n > 1) THEN
  GO TO 50
END IF
kt = 2
50 nn = n
IF (x < 0.0) THEN
  GO TO 490
ELSE IF (x > 0.0) THEN
  GO TO 110
END IF
IF (alpha < 0.0) THEN
  GO TO 470
ELSE IF (alpha > 0.0) THEN
  GO TO 80
END IF
y(1) = 1.0
IF (n == 1) RETURN
i1 = 2
GO TO 90
80 i1 = 1
90 DO i = i1, n
  y(i) = 0.0
END DO
RETURN

110 IF (alpha < 0.0) GO TO 470

ialp = INT(alpha)
fni = REAL(ialp+n-1)
fnf = alpha - REAL(ialp)
dfn = fni + fnf
fnu = dfn
xo2 = x * 0.5
sxo2 = xo2 * xo2

!     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
!     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE APPLIED.

IF (sxo2 > (fnu+1.0)) THEN
  ta = MAX(20.0,fnu)
  IF (x > ta) GO TO 130
  IF (x > 12.0) GO TO 120
  xo2l = LOG(xo2)
  ns = INT(sxo2-fnu) + 1
ELSE
  fn = fnu
  fnp1 = fn + 1.0
  xo2l = LOG(xo2)
  is = kt
  IF (x <= 0.50) GO TO 170
  ns = 0
END IF
fni = fni + REAL(ns)
dfn = fni + fnf
fn = dfn
fnp1 = fn + 1.0
is = kt
IF (n-1+ns > 0) is = 3
GO TO 170

120 ans = MAX(36.0-fnu,0.0)
ns = INT(ans)
fni = fni + REAL(ns)
dfn = fni + fnf
fn = dfn
is = kt
IF (n-1+ns > 0) is = 3
GO TO 140

130 rtx = SQRT(x)
tau = rtwo * rtx
ta = tau + fnulim(kt)
IF (fnu <= ta) GO TO 310
fn = fnu
is = kt

!     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY

140 i1 = ABS(3-is)
i1 = MAX(i1,1)
flgjy = 1.0
CALL asjy(jairy, x, fn, flgjy, i1, tol, elim, temp(is:), wk, iflw)
IF (iflw /= 0) GO TO 220
SELECT CASE (is)
  CASE (1)
    GO TO 160
  CASE (2)
    GO TO 280
  CASE (3)
    GO TO 390
END SELECT

150 temp(1) = temp(3)
kt = 1
160 is = 2
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
IF (i1 == 2) GO TO 280
GO TO 140

!     SERIES FOR (X/2)**2<=NU+1

170 gln = gamln(fnp1)
arg = fn * xo2l - gln
IF (arg < (-elim)) GO TO 240
earg = EXP(arg)

180 s = 1.0
IF (x >= tol) THEN
  ak = 3.0
  t2 = 1.0
  t = 1.0
  s1 = fn
  DO k = 1, 17
    s2 = t2 + s1
    t = -t * sxo2 / s2
    s = s + t
    IF (ABS(t) < tol) EXIT
    t2 = t2 + ak
    ak = ak + 2.0
    s1 = s1 + fn
  END DO
END IF

temp(is) = s * earg
IF (is == 2) THEN
  GO TO 280
ELSE IF (is == 3) THEN
  GO TO 380
END IF
earg = earg * fn / xo2
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
is = 2
GO TO 180

!     SET UNDERFLOW VALUE AND UPDATE PARAMETERS

220 y(nn) = 0.0
nn = nn - 1
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
IF (nn < 1) THEN
  GO TO 270
ELSE IF (nn > 1) THEN
  GO TO 140
END IF
kt = 2
is = 2
GO TO 140

240 y(nn) = 0.0
nn = nn - 1
fnp1 = fn
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
IF (nn < 1) THEN
  GO TO 270
ELSE IF (nn > 1) THEN
  GO TO 260
END IF
kt = 2
is = 2
260 IF (sxo2 > fnp1) THEN
  GO TO 140
END IF
arg = arg - xo2l + LOG(fnp1)
IF (arg < -elim) GO TO 240
GO TO 170

270 nz = n - nn
RETURN

!     BACKWARD RECURSION SECTION

280 nz = n - nn
IF (kt /= 2) THEN
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
  y(nn) = temp(1)
  y(nn-1) = temp(2)
  IF (nn == 2) RETURN
  trx = 2.0 / x
  dtm = fni
  tm = (dtm+fnf) * trx
  k = nn + 1
  DO i = 3, nn
    k = k - 1
    y(k-2) = tm * y(k-1) - y(k)
    dtm = dtm - 1.0
    tm = (dtm+fnf) * trx
  END DO
  RETURN
END IF
300 y(1) = temp(2)
RETURN

!     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN
!     OSCILLATORY REGION X>MAX(20, NU), PROVIDED THE LAST MEMBER
!     OF THE SEQUENCE IS ALSO IN THE REGION.

310 in = INT(alpha-tau+2.0)
IF (in > 0) THEN
  idalp = ialp - in - 1
  kt = 1
ELSE
  idalp = ialp
  in = 0
END IF
is = kt
fidal = REAL(idalp)
dalpha = fidal + fnf
arg = x - pidt * dalpha - pdf
sa = SIN(arg)
sb = COS(arg)
coef = rttp / rtx
etx = 8.0 * x

320 dtm = fidal + fidal
dtm = dtm * dtm
tm = 0.0
IF (fidal /= 0.0.OR.ABS(fnf) >= tol) THEN
  tm = 4.0 * fnf * (fidal+fidal+fnf)
END IF
trx = dtm - 1.0
t2 = (trx+tm) / etx
s2 = t2
relb = tol * ABS(t2)
t1 = etx
s1 = 1.0
fn = 1.0
ak = 8.0
DO k = 1, 13
  t1 = t1 + etx
  fn = fn + ak
  trx = dtm - fn
  ap = trx + tm
  t2 = -t2 * ap / t1
  s1 = s1 + t2
  t1 = t1 + etx
  ak = ak + 8.0
  fn = fn + ak
  trx = dtm - fn
  ap = trx + tm
  t2 = t2 * ap / t1
  s2 = s2 + t2
  IF (ABS(t2) <= relb) GO TO 340
  ak = ak + 8.0
END DO
340 temp(is) = coef * (s1*sb-s2*sa)
IF (is /= 2) THEN
  fidal = fidal + 1.0
  dalpha = fidal + fnf
  is = 2
  tb = sa
  sa = -sb
  sb = tb
  GO TO 320
END IF

!     FORWARD RECURSION SECTION

IF (kt == 2) GO TO 300
s1 = temp(1)
s2 = temp(2)
tx = 2.0 / x
tm = dalpha * tx
IF (in /= 0) THEN

!     FORWARD RECUR TO INDEX ALPHA

  DO i = 1, in
    s = s2
    s2 = tm * s2 - s1
    tm = tm + tx
    s1 = s
  END DO
  IF (nn == 1) GO TO 370
  s = s2
  s2 = tm * s2 - s1
  tm = tm + tx
  s1 = s
END IF

!     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1

y(1) = s1
y(2) = s2
IF (nn == 2) RETURN
DO i = 3, nn
  y(i) = tm * y(i-1) - y(i-2)
  tm = tm + tx
END DO
RETURN
370 y(1) = s2
RETURN

!     BACKWARD RECURSION WITH NORMALIZATION BY
!     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.

!     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
380 akm = MAX(3.0-fn,0.0)
km = INT(akm)
tfn = fn + REAL(km)
ta = (gln+tfn-0.9189385332-0.0833333333/tfn) / (tfn+0.5)
ta = xo2l - ta
tb = -(1.0-1.5/tfn) / tfn
akm = tolln / (-ta+SQRT(ta*ta-tolln*tb)) + 1.5
in = km + INT(akm)
GO TO 430

!     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
390 gln = wk(3) + wk(2)
IF (wk(6) > 30.0) GO TO 410
rden = (pp(4)*wk(6)+pp(3)) * wk(6) + 1.0
rzden = pp(1) + pp(2) * wk(6)
ta = rzden / rden
IF (wk(1) >= 0.10) THEN
  tb = gln / wk(5)
  GO TO 420
END IF
400 tb = (1.259921049 + (0.1679894730 + 0.0887944358*wk(1))*wk(1)) /wk(7)
GO TO 420

410 ta = 0.5 * tolln / wk(4)
ta = ((0.0493827160*ta - 0.1111111111)*ta + 0.6666666667) * ta *wk(6)
IF (wk(1) < 0.10) GO TO 400
tb = gln / wk(5)
420 in = INT(ta/tb+1.5)
IF (in > inlim) GO TO 150

430 dtm = fni + REAL(in)
trx = 2.0 / x
tm = (dtm+fnf) * trx
ta = 0.0
tb = tol
kk = 1

!     BACKWARD RECUR UNINDEXED

440 DO i = 1, in
  s = tb
  tb = tm * tb - ta
  ta = s
  dtm = dtm - 1.0
  tm = (dtm+fnf) * trx
END DO
!     NORMALIZATION
IF (kk == 1) THEN
  ta = (ta/tb) * temp(3)
  tb = temp(3)
  kk = 2
  in = ns
  IF (ns /= 0) GO TO 440
END IF
y(nn) = tb
nz = n - nn
IF (nn == 1) RETURN
k = nn - 1
y(k) = tm * tb - ta
IF (nn == 2) RETURN
dtm = dtm - 1.0
tm = (dtm+fnf) * trx
km = k - 1

!     BACKWARD RECUR INDEXED

DO i = 1, km
  y(k-1) = tm * y(k) - y(k+1)
  dtm = dtm - 1.0
  tm = (dtm + fnf) * trx
  k = k - 1
END DO
RETURN


470 nz = -2
RETURN

480 nz = -3
RETURN

490 nz = -1
RETURN
END SUBROUTINE besj


SUBROUTINE asjy(funjy, x, fnu, flgjy, in, tol, elim, y, wk, iflw)

!            ASJY COMPUTES BESSEL FUNCTIONS J AND Y
!            FOR ARGUMENTS X > 0.0 AND ORDERS FNU >= 35.0
!            ON FLGJY = 1 AND FLGJY = -1 RESPECTIVELY

!                               INPUT

!   FUNJY - EXTERNAL FUNCTION JAIRY OR YAIRY
!       X - ARGUMENT, X > 0.0
!     FNU - ORDER OF THE FIRST BESSEL FUNCTION
!   FLGJY - SELECTION FLAG
!           FLGJY =  1.0 GIVES THE J FUNCTION
!           FLGJY = -1.0 GIVES THE Y FUNCTION
!      IN - NUMBER OF FUNCTIONS DESIRED, IN = 1 OR 2
!     TOL - TOLERANCE SPECIFIED BY BESJ OR BESY
!    ELIM - TOLERANCE SPECIFIED BY BESJ OR BESY

!                               OUTPUT

!      Y  - A VECTOR WHOSE FIRST IN COMPONENTS CONTAIN THE SEQUENCE
!    IFLW - A FLAG INDICATING UNDERFLOW OR OVERFLOW
!                 RETURN VARIABLES FOR BESJ ONLY
!   WK(1) = 1 - (X/FNU)**2 = W**2
!   WK(2) = SQRT(ABS(WK(1)))
!   WK(3) = ABS(WK(2) - ATAN(WK(2)))  OR
!           ABS(LN((1 + WK(2))/(X/FNU)) - WK(2))
!         = ABS((2/3)*ZETA**(3/2))
!   WK(4) = FNU*WK(3)
!   WK(5) = (1.5*WK(3)*FNU)**(1/3) = SQRT(ZETA)*FNU**(1/3)
!   WK(6) = SIGN(1.,W**2)*WK(5)**2 = SIGN(1.,W**2)*ZETA*FNU**(2/3)
!   WK(7) = FNU**(1/3)

!                               WRITTEN BY
!                               D. E. AMOS

!  ABSTRACT
!      ASJK IMPLEMENTS THE UNIFORM ASYMPTOTIC EXPANSION OF THE J AND Y
!      BESSEL FUNCTIONS FOR FNU >= 35 AND REAL X > 0.0.
!      THE FORMS ARE IDENTICAL EXCEPT FOR A CHANGE IN SIGN OF SOME OF
!      THE TERMS.  THIS CHANGE IN SIGN IS ACCOMPLISHED BY MEANS OF THE
!      FLAG FLGJY = 1 OR -1.  ON FLGJY = 1 THE AIRY FUNCTIONS AI(X) AND
!      DAI(X) ARE SUPPLIED BY THE EXTERNAL FUNCTION JAIRY, AND ON FLGJY = -1
!      THE AIRY FUNCTIONS BI(X) AND DBI(X) ARE SUPPLIED BY THE EXTERNAL
!      FUNCTION YAIRY.
!  N.B. JAIRY & YAIRY are actually subroutines, not functions.

REAL, INTENT(IN)     :: x, fnu, flgjy, tol, elim
INTEGER, INTENT(IN)  :: in
REAL, INTENT(OUT)    :: y(:), wk(:)
INTEGER, INTENT(OUT) :: iflw

INTERFACE
  SUBROUTINE funjy(x, rx, c, ai, dai)
    REAL, INTENT(IN)  :: x, rx, c
    REAL, INTENT(OUT) :: ai, dai
  END SUBROUTINE funjy
END INTERFACE

! Local variables
INTEGER :: i, j, jn, jr, ju, k, kb, klast, kmax(5), kp1, ks,  &
           ksp1, kstemp, l, lr, lrp1
REAL    :: abw2, akm, ap, asum, az, bsum, cr(10), crz32, dfi, dr(10), fi,  &
           fn, fn2, phi, rcz, rden, relb, rfn2, rtz, rzden, sa, sb, suma,  &
           sumb, s1, ta, tau, tb, tfn, t2, upol(10), xx, z, z32
REAL, PARAMETER :: alfa(26,4) = RESHAPE( (/ -4.44444444444444E-03,  &
      -9.22077922077922E-04, -8.84892884892885E-05,  1.65927687832450E-04,  &
       2.46691372741793E-04,  2.65995589346255E-04,  2.61824297061501E-04,  &
       2.48730437344656E-04,  2.32721040083232E-04,  2.16362485712365E-04,  &
       2.00738858762752E-04,  1.86267636637545E-04,  1.73060775917876E-04,  &
       1.61091705929016E-04,  1.50274774160908E-04,  1.40503497391270E-04,  &
       1.31668816545923E-04,  1.23667445598253E-04,  1.16405271474738E-04,  &
       1.09798298372713E-04,  1.03772410422993E-04,  9.82626078369363E-05,  &
       9.32120517249503E-05,  8.85710852478712E-05,  8.42963105715700E-05,  &
       8.03497548407791E-05,  6.93735541354589E-04,  2.32241745182922E-04,  &
      -1.41986273556691E-05, -1.16444931672049E-04, -1.50803558053049E-04,  &
      -1.55121924918096E-04, -1.46809756646466E-04, -1.33815503867491E-04,  &
      -1.19744975684254E-04, -1.06184319207974E-04, -9.37699549891194E-05,  &
      -8.26923045588193E-05, -7.29374348155221E-05, -6.44042357721016E-05,  &
      -5.69611566009369E-05, -5.04731044303562E-05, -4.48134868008883E-05,  &
      -3.98688727717599E-05, -3.55400532972042E-05, -3.17414256609022E-05,  &
      -2.83996793904175E-05, -2.54522720634871E-05, -2.28459297164725E-05,  &
      -2.05352753106481E-05, -1.84816217627666E-05, -1.66519330021394E-05,  &
      -3.54211971457744E-04, -1.56161263945159E-04,  3.04465503594936E-05,  &
       1.30198655773243E-04,  1.67471106699712E-04,  1.70222587683593E-04,  &
       1.56501427608595E-04,  1.36339170977445E-04,  1.14886692029825E-04,  &
       9.45869093034688E-05,  7.64498419250898E-05,  6.07570334965197E-05,  &
       4.74394299290509E-05,  3.62757512005344E-05,  2.69939714979225E-05,  &
       1.93210938247939E-05,  1.30056674793963E-05,  7.82620866744497E-06,  &
       3.59257485819352E-06,  1.44040049814252E-07, -2.65396769697939E-06,  &
      -4.91346867098486E-06, -6.72739296091248E-06, -8.17269379678658E-06,  &
      -9.31304715093561E-06, -1.02011418798016E-05,  3.78194199201773E-04,  &
       2.02471952761816E-04, -6.37938506318862E-05, -2.38598230603006E-04,  &
      -3.10916256027362E-04, -3.13680115247576E-04, -2.78950273791323E-04,  &
      -2.28564082619141E-04, -1.75245280340847E-04, -1.25544063060690E-04,  &
      -8.22982872820208E-05, -4.62860730588116E-05, -1.72334302366962E-05,  &
       5.60690482304602E-06,  2.31395443148287E-05,  3.62642745856794E-05,  &
       4.58006124490189E-05,  5.24595294959114E-05,  5.68396208545815E-05,  &
       5.94349820393104E-05,  6.06478527578422E-05,  6.08023907788436E-05,  &
       6.01577894539460E-05,  5.89199657344698E-05,  5.72515823777593E-05,  &
       5.52804375585853E-05 /), (/ 26, 4 /) )
REAL, PARAMETER :: ar(8) = (/ 8.35503472222222E-02, 1.28226574556327E-01,  &
        2.91849026464140E-01, 8.81627267443758E-01, 3.32140828186277E+00,  &
        1.49957629868626E+01, 7.89230130115865E+01, 4.74451538868264E+02 /)
REAL, PARAMETER :: beta(26,5) = RESHAPE( (/    1.79988721413553E-02, 5.59964911064388E-03, &
   2.88501402231133E-03, 1.80096606761054E-03, 1.24753110589199E-03, 9.22878876572938E-04, &
   7.14430421727287E-04, 5.71787281789705E-04, 4.69431007606482E-04, 3.93232835462917E-04, &
   3.34818889318298E-04, 2.88952148495752E-04, 2.52211615549573E-04, 2.22280580798883E-04, &
   1.97541838033063E-04, 1.76836855019718E-04, 1.59316899661821E-04, 1.44347930197334E-04, &
   1.31448068119965E-04, 1.20245444949303E-04, 1.10449144504599E-04, 1.01828770740567E-04, &
   9.41998224204238E-05, 8.74130545753834E-05, 8.13466262162801E-05, 7.59002269646219E-05, &
  -1.49282953213429E-03,-8.78204709546389E-04,-5.02916549572035E-04,-2.94822138512746E-04, &
  -1.75463996970783E-04,-1.04008550460816E-04,-5.96141953046458E-05,-3.12038929076098E-05, &
  -1.26089735980230E-05,-2.42892608575730E-07, 8.05996165414274E-06, 1.36507009262147E-05, &
   1.73964125472926E-05, 1.98672978842134E-05, 2.14463263790823E-05, 2.23954659232457E-05, &
   2.28967783814713E-05, 2.30785389811178E-05, 2.30321976080909E-05, 2.28236073720349E-05, &
   2.25005881105292E-05, 2.20981015361991E-05, 2.16418427448104E-05, 2.11507649256221E-05, &
   2.06388749782171E-05, 2.01165241997082E-05, 5.52213076721293E-04, 4.47932581552385E-04, &
   2.79520653992021E-04, 1.52468156198447E-04, 6.93271105657044E-05, 1.76258683069991E-05, &
  -1.35744996343269E-05,-3.17972413350427E-05,-4.18861861696693E-05,-4.69004889379141E-05, &
  -4.87665447413787E-05,-4.87010031186735E-05,-4.74755620890087E-05,-4.55813058138628E-05, &
  -4.33309644511266E-05,-4.09230193157750E-05,-3.84822638603221E-05,-3.60857167535411E-05, &
  -3.37793306123367E-05,-3.15888560772110E-05,-2.95269561750807E-05,-2.75978914828336E-05, &
  -2.58006174666884E-05,-2.41308356761280E-05,-2.25823509518346E-05,-2.11479656768913E-05, &
  -4.74617796559960E-04,-4.77864567147321E-04,-3.20390228067038E-04,-1.61105016119962E-04, &
  -4.25778101285435E-05, 3.44571294294968E-05, 7.97092684075675E-05, 1.03138236708272E-04, &
   1.12466775262204E-04, 1.13103642108481E-04, 1.08651634848774E-04, 1.01437951597662E-04, &
   9.29298396593364E-05, 8.40293133016090E-05, 7.52727991349134E-05, 6.69632521975731E-05, &
   5.92564547323195E-05, 5.22169308826976E-05, 4.58539485165361E-05, 4.01445513891487E-05, &
   3.50481730031328E-05, 3.05157995034347E-05, 2.64956119950516E-05, 2.29363633690998E-05, &
   1.97893056664022E-05, 1.70091984636413E-05, 7.36465810572578E-04, 8.72790805146194E-04, &
   6.22614862573135E-04, 2.85998154194304E-04, 3.84737672879366E-06,-1.87906003636972E-04, &
  -2.97603646594555E-04,-3.45998126832656E-04,-3.53382470916038E-04,-3.35715635775049E-04, &
  -3.04321124789040E-04,-2.66722723047613E-04,-2.27654214122820E-04,-1.89922611854562E-04, &
  -1.55058918599094E-04,-1.23778240761874E-04,-9.62926147717644E-05,-7.25178327714425E-05, &
  -5.22070028895634E-05,-3.50347750511901E-05,-2.06489761035552E-05,-8.70106096849767E-06, &
   1.13698686675100E-06, 9.16426474122779E-06, 1.56477785428873E-05, 2.08223629482467E-05 /), &
   (/ 26, 5 /) )
REAL, PARAMETER :: br(10) = (/ -1.45833333333333E-01,  &
        -9.87413194444444E-02, -1.43312053915895E-01, -3.17227202678414E-01, &
        -9.42429147957120E-01, -3.51120304082635E+00, -1.57272636203680E+01, &
        -8.22814390971859E+01, -4.92355370523671E+02, -3.31621856854797E+03 /)
REAL, PARAMETER :: c(65) = (/  -2.08333333333333E-01,  1.25000000000000E-01,  &
         3.34201388888889E-01, -4.01041666666667E-01,  7.03125000000000E-02,  &
        -1.02581259645062E+00,  1.84646267361111E+00, -8.91210937500000E-01,  &
         7.32421875000000E-02,  4.66958442342625E+00, -1.12070026162230E+01,  &
         8.78912353515625E+00, -2.36408691406250E+00,  1.12152099609375E-01,  &
        -2.82120725582002E+01,  8.46362176746007E+01, -9.18182415432400E+01,  &
         4.25349987453885E+01, -7.36879435947963E+00,  2.27108001708984E-01,  &
         2.12570130039217E+02, -7.65252468141182E+02,  1.05999045252800E+03,  &
        -6.99579627376133E+02,  2.18190511744212E+02, -2.64914304869516E+01,  &
         5.72501420974731E-01, -1.91945766231841E+03,  8.06172218173731E+03,  &
        -1.35865500064341E+04,  1.16553933368645E+04, -5.30564697861340E+03,  &
         1.20090291321635E+03, -1.08090919788395E+02,  1.72772750258446E+00,  &
         2.02042913309661E+04, -9.69805983886375E+04,  1.92547001232532E+05,  &
        -2.03400177280416E+05,  1.22200464983017E+05, -4.11926549688976E+04,  &
         7.10951430248936E+03, -4.93915304773088E+02,  6.07404200127348E+00,  &
        -2.42919187900551E+05,  1.31176361466298E+06, -2.99801591853811E+06,  &
         3.76327129765640E+06, -2.81356322658653E+06,  1.26836527332162E+06,  &
        -3.31645172484564E+05,  4.52187689813627E+04, -2.49983048181121E+03,  &
         2.43805296995561E+01,  3.28446985307204E+06, -1.97068191184322E+07,  &
         5.09526024926646E+07, -7.41051482115327E+07,  6.63445122747290E+07,  &
        -3.75671766607634E+07,  1.32887671664218E+07, -2.78561812808645E+06,  &
         3.08186404612662E+05, -1.38860897537170E+04,  1.10017140269247E+02 /)
REAL, PARAMETER :: con1 = 6.66666666666667E-01, con2 = 3.33333333333333E-01,  &
        con548 = 1.04166666666667E-01,  &
        gama(26) = (/ 6.29960524947437E-01, 2.51984209978975E-01,  &
        1.54790300415656E-01, 1.10713062416159E-01, 8.57309395527395E-02,  &
        6.97161316958684E-02, 5.86085671893714E-02, 5.04698873536311E-02,  &
        4.42600580689155E-02, 3.93720661543510E-02, 3.54283195924455E-02,  &
        3.21818857502098E-02, 2.94646240791158E-02, 2.71581677112934E-02,  &
        2.51768272973862E-02, 2.34570755306079E-02, 2.19508390134907E-02,  &
        2.06210828235646E-02, 1.94388240897881E-02, 1.83810633800683E-02,  &
        1.74293213231963E-02, 1.65685837786612E-02, 1.57865285987918E-02,  &
        1.50729501494096E-02, 1.44193250839955E-02, 1.38184805735342E-02 /)
REAL, PARAMETER :: tols = -6.90775527898214E+00

fn = fnu
iflw = 0
DO jn = 1, in
  xx = x / fn
  wk(1) = 1.0 - xx * xx
  abw2 = ABS(wk(1))
  wk(2) = SQRT(abw2)
  wk(7) = fn ** con2
  IF (abw2 > 0.27750) GO TO 80

!     ASYMPTOTIC EXPANSION
!     CASES NEAR X=FN, ABS(1.-(X/FN)**2) <= 0.2775
!     COEFFICIENTS OF ASYMPTOTIC EXPANSION BY SERIES

!     ZETA AND TRUNCATION FOR A(ZETA) AND B(ZETA) SERIES

!     KMAX IS TRUNCATION INDEX FOR A(ZETA) AND B(ZETA) SERIES=MAX(2,SA)

  sa = 0.0
  IF (abw2 /= 0.0) THEN
    sa = tols / LOG(abw2)
  END IF
  sb = sa
  DO i = 1, 5
    akm = MAX(sa,2.0)
    kmax(i) = INT(akm)
    sa = sa + sb
  END DO
  kb = kmax(5)
  klast = kb - 1
  sa = gama(kb)
  DO k = 1, klast
    kb = kb - 1
    sa = sa * wk(1) + gama(kb)
  END DO
  z = wk(1) * sa
  az = ABS(z)
  rtz = SQRT(az)
  wk(3) = con1 * az * rtz
  wk(4) = wk(3) * fn
  wk(5) = rtz * wk(7)
  wk(6) = -wk(5) * wk(5)
  IF (z > 0.0) THEN
    IF (wk(4) > elim) GO TO 70
    wk(6) = -wk(6)
  END IF
  phi = SQRT(SQRT(sa+sa+sa+sa))

!     B(ZETA) FOR S=0

  kb = kmax(5)
  klast = kb - 1
  sb = beta(kb,1)
  DO k = 1, klast
    kb = kb - 1
    sb = sb * wk(1) + beta(kb,1)
  END DO
  ksp1 = 1
  fn2 = fn * fn
  rfn2 = 1.0 / fn2
  rden = 1.0
  asum = 1.0
  relb = tol * ABS(sb)
  bsum = sb
  DO ks = 1, 4
    ksp1 = ksp1 + 1
    rden = rden * rfn2

!     A(ZETA) AND B(ZETA) FOR S=1,2,3,4

    kstemp = 5 - ks
    kb = kmax(kstemp)
    klast = kb - 1
    sa = alfa(kb,ks)
    sb = beta(kb,ksp1)
    DO k = 1, klast
      kb = kb - 1
      sa = sa * wk(1) + alfa(kb,ks)
      sb = sb * wk(1) + beta(kb,ksp1)
    END DO
    ta = sa * rden
    tb = sb * rden
    asum = asum + ta
    bsum = bsum + tb
    IF (ABS(ta) <= tol .AND. ABS(tb) <= relb) GO TO 60
  END DO

  60 bsum = bsum / (fn*wk(7))
  GO TO 140

  70 iflw = 1
  RETURN

  80 upol(1) = 1.0
  tau = 1.0 / wk(2)
  t2 = 1.0 / wk(1)
  IF (wk(1) < 0.0) THEN

!     CASES FOR (X/FN)>SQRT(1.2775)

    wk(3) = ABS(wk(2)-ATAN(wk(2)))
    wk(4) = wk(3) * fn
    rcz = -con1 / wk(4)
    z32 = 1.5 * wk(3)
    rtz = z32 ** con2
    wk(5) = rtz * wk(7)
    wk(6) = -wk(5) * wk(5)
  ELSE

!     CASES FOR (X/FN)<SQRT(0.7225)

    wk(3) = ABS(LOG((1.0+wk(2))/xx)-wk(2))
    wk(4) = wk(3) * fn
    rcz = con1 / wk(4)
    IF (wk(4) > elim) GO TO 70
    z32 = 1.5 * wk(3)
    rtz = z32 ** con2
    wk(7) = fn ** con2
    wk(5) = rtz * wk(7)
    wk(6) = wk(5) * wk(5)
  END IF
  phi = SQRT((rtz+rtz)*tau)
  tb = 1.0
  asum = 1.0
  tfn = tau / fn
  upol(2) = (c(1)*t2+c(2)) * tfn
  crz32 = con548 * rcz
  bsum = upol(2) + crz32
  relb = tol * ABS(bsum)
  ap = tfn
  ks = 0
  kp1 = 2
  rzden = rcz
  l = 2
  DO lr = 2, 8, 2

!     COMPUTE TWO U POLYNOMIALS FOR NEXT A(ZETA) AND B(ZETA)

    lrp1 = lr + 1
    DO k = lr, lrp1
      ks = ks + 1
      kp1 = kp1 + 1
      l = l + 1
      s1 = c(l)
      DO j = 2, kp1
        l = l + 1
        s1 = s1 * t2 + c(l)
      END DO
      ap = ap * tfn
      upol(kp1) = ap * s1
      cr(ks) = br(ks) * rzden
      rzden = rzden * rcz
      dr(ks) = ar(ks) * rzden
    END DO
    suma = upol(lrp1)
    sumb = upol(lr+2) + upol(lrp1) * crz32
    ju = lrp1
    DO jr = 1, lr
      ju = ju - 1
      suma = suma + cr(jr) * upol(ju)
      sumb = sumb + dr(jr) * upol(ju)
    END DO
    tb = -tb
    IF (wk(1) > 0.0) tb = ABS(tb)
    asum = asum + suma * tb
    bsum = bsum + sumb * tb
    IF (ABS(suma) <= tol .AND. ABS(sumb) <= relb) GO TO 130
  END DO
  130 tb = wk(5)
  IF (wk(1) > 0.0) tb = -tb
  bsum = bsum / tb

  140 CALL funjy(wk(6), wk(5), wk(4), fi, dfi)
  y(jn) = flgjy * phi * (fi*asum+dfi*bsum) / wk(7)
  fn = fn - flgjy
END DO
RETURN
END SUBROUTINE asjy


SUBROUTINE jairy(x, rx, c, ai, dai)

!                  JAIRY COMPUTES THE AIRY FUNCTION AI(X)
!                   AND ITS DERIVATIVE DAI(X) FOR JBESS

!                                   INPUT

!         X - ARGUMENT, COMPUTED BY JBESS, X UNRESTRICTED
!        RX - RX=SQRT(ABS(X)), COMPUTED BY JBESS
!         C - C=2.*(ABS(X)**1.5)/3., COMPUTED BY JBESS

!                                  OUTPUT

!        AI - VALUE OF FUNCTION AI(X)
!       DAI - VALUE OF THE DERIVATIVE DAI(X)

!                                WRITTEN BY

!                                D. E. AMOS
!                               S. L. DANIEL
!                               M. K. WESTON

REAL, INTENT(IN)  :: x, rx, c
REAL, INTENT(OUT) :: ai, dai

! Local variables
REAL    :: ak1(14) = (/ 2.20423090987793E-01, -1.25290242787700E-01,  &
        1.03881163359194E-02, 8.22844152006343E-04, -2.34614345891226E-04,  &
        1.63824280172116E-05, 3.06902589573189E-07, -1.29621999359332E-07,  &
        8.22908158823668E-09, 1.53963968623298E-11, -3.39165465615682E-11,  &
        2.03253257423626E-12, -1.10679546097884E-14, -5.16169497785080E-15  &
        /), ak2(23) = (/ 2.74366150869598E-01, 5.39790969736903E-03,  &
        -1.57339220621190E-03, 4.27427528248750E-04, -1.12124917399925E-04  &
       , 2.88763171318904E-05, -7.36804225370554E-06, 1.87290209741024E-06  &
       , -4.75892793962291E-07, 1.21130416955909E-07, -3.09245374270614E-08  &
       , 7.92454705282654E-09, -2.03902447167914E-09, 5.26863056595742E-10  &
       , -1.36704767639569E-10, 3.56141039013708E-11, -9.31388296548430E-12  &
       , 2.44464450473635E-12, -6.43840261990955E-13, 1.70106030559349E-13  &
       , -4.50760104503281E-14, 1.19774799164811E-14, -3.19077040865066E-15  &
        /), ak3(14) = (/ 2.80271447340791E-01, -1.78127042844379E-03,  &
        4.03422579628999E-05, -1.63249965269003E-06, 9.21181482476768E-08,  &
        -6.52294330229155E-09, 5.47138404576546E-10, -5.24408251800260E-11  &
       , 5.60477904117209E-12, -6.56375244639313E-13, 8.31285761966247E-14  &
       , -1.12705134691063E-14, 1.62267976598129E-15, -2.46480324312426E-16 /)
REAL    :: ajp(19) = (/ 7.78952966437581E-02, -1.84356363456801E-01,  &
        3.01412605216174E-02, 3.05342724277608E-02, -4.95424702513079E-03,  &
        -1.72749552563952E-03, 2.43137637839190E-04, 5.04564777517082E-05,  &
        -6.16316582695208E-06, -9.03986745510768E-07, 9.70243778355884E-08  &
       , 1.09639453305205E-08, -1.04716330588766E-09, -9.60359441344646E-11  &
       , 8.25358789454134E-12, 6.36123439018768E-13, -4.96629614116015E-14  &
       , -3.29810288929615E-15, 2.35798252031104E-16 /), ajn(19)  &
        = (/ 3.80497887617242E-02, -2.45319541845546E-01,  &
        1.65820623702696E-01, 7.49330045818789E-02, -2.63476288106641E-02,  &
        -5.92535597304981E-03, 1.44744409589804E-03, 2.18311831322215E-04,  &
        -4.10662077680304E-05, -4.66874994171766E-06, 7.15218807277160E-07  &
       , 6.52964770854633E-08, -8.44284027565946E-09, -6.44186158976978E-10  &
       , 7.20802286505285E-11, 4.72465431717846E-12, -4.66022632547045E-13  &
       , -2.67762710389189E-14, 2.36161316570019E-15 /), a(15)  &
        = (/ 4.90275424742791E-01, 1.57647277946204E-03,  &
        -9.66195963140306E-05, 1.35916080268815E-07, 2.98157342654859E-07,  &
        -1.86824767559979E-08, -1.03685737667141E-09, 3.28660818434328E-10  &
       , -2.57091410632780E-11, -2.32357655300677E-12, 9.57523279048255E-13  &
       , -1.20340828049719E-13, -2.90907716770715E-15, 4.55656454580149E-15  &
       , -9.99003874810259E-16 /), b(15) = (/ 2.78593552803079E-01,  &
        -3.52915691882584E-03, -2.31149677384994E-05, 4.71317842263560E-06  &
       , -1.12415907931333E-07, -2.00100301184339E-08, 2.60948075302193E-09  &
       , -3.55098136101216E-11, -3.50849978423875E-11, 5.83007187954202E-12  &
       , -2.04644828753326E-13, -1.10529179476742E-13, 2.87724778038775E-14  &
       , -2.88205111009939E-15, -3.32656311696166E-16 /)
REAL    :: dak1(14) = (/ 2.04567842307887E-01, -6.61322739905664E-02,  &
        -8.49845800989287E-03, 3.12183491556289E-03, -2.70016489829432E-04  &
       , -6.35636298679387E-06, 3.02397712409509E-06, -2.18311195330088E-07  &
       , -5.36194289332826E-10, 1.13098035622310E-09, -7.43023834629073E-11  &
       , 4.28804170826891E-13, 2.23810925754539E-13, -1.39140135641182E-14  &
        /), dak2(24) = (/ 2.93332343883230E-01, -8.06196784743112E-03,  &
        2.42540172333140E-03, -6.82297548850235E-04, 1.85786427751181E-04,  &
        -4.97457447684059E-05, 1.32090681239497E-05, -3.49528240444943E-06  &
       , 9.24362451078835E-07, -2.44732671521867E-07, 6.49307837648910E-08  &
       , -1.72717621501538E-08, 4.60725763604656E-09, -1.23249055291550E-09  &
       , 3.30620409488102E-10, -8.89252099772401E-11, 2.39773319878298E-11  &
       , -6.48013921153450E-12, 1.75510132023731E-12, -4.76303829833637E-13  &
       , 1.29498241100810E-13, -3.52679622210430E-14, 9.62005151585923E-15  &
       , -2.62786914342292E-15 /), dak3(14) = (/ 2.84675828811349E-01,  &
        2.53073072619080E-03, -4.83481130337976E-05, 1.84907283946343E-06,  &
        -1.01418491178576E-07, 7.05925634457153E-09, -5.85325291400382E-10  &
       , 5.56357688831339E-11, -5.90889094779500E-12, 6.88574353784436E-13  &
       , -8.68588256452194E-14, 1.17374762617213E-14, -1.68523146510923E-15  &
       , 2.55374773097056E-16 /)
REAL    :: dajp(19) = (/ 6.53219131311457E-02,  &
        -1.20262933688823E-01, 9.78010236263823E-03, 1.67948429230505E-02,  &
        -1.97146140182132E-03, -8.45560295098867E-04, 9.42889620701976E-05  &
       , 2.25827860945475E-05, -2.29067870915987E-06, -3.76343991136919E-07  &
       , 3.45663933559565E-08, 4.29611332003007E-09, -3.58673691214989E-10  &
       , -3.57245881361895E-11, 2.72696091066336E-12, 2.26120653095771E-13  &
       , -1.58763205238303E-14, -1.12604374485125E-15, 7.31327529515367E-17  &
        /), dajn(19) = (/ 1.08594539632967E-02, 8.53313194857091E-02,  &
        -3.15277068113058E-01, -8.78420725294257E-02, 5.53251906976048E-02  &
       , 9.41674060503241E-03, -3.32187026018996E-03, -4.11157343156826E-04  &
       , 1.01297326891346E-04, 9.87633682208396E-06, -1.87312969812393E-06  &
       , -1.50798500131468E-07, 2.32687669525394E-08, 1.59599917419225E-09  &
       , -2.07665922668385E-10, -1.24103350500302E-11, 1.39631765331043E-12  &
       , 7.39400971155740E-14, -7.32887475627500E-15 /),  &
       da(15) = (/ 4.91627321104601E-01, 3.11164930427489E-03,  &
        8.23140762854081E-05, -4.61769776172142E-06, -6.13158880534626E-08  &
       , 2.87295804656520E-08, -1.81959715372117E-09, -1.44752826642035E-10  &
       , 4.53724043420422E-11, -3.99655065847223E-12, -3.24089119830323E-13  &
       , 1.62098952568741E-13, -2.40765247974057E-14, 1.69384811284491E-16  &
       , 8.17900786477396E-16 /), db(15) = (/ -2.77571356944231E-01,  &
        4.44212833419920E-03, -8.42328522190089E-05, -2.58040318418710E-06  &
       , 3.42389720217621E-07, -6.24286894709776E-09, -2.36377836844577E-09  &
       , 3.16991042656673E-10, -4.40995691658191E-12, -5.18674221093575E-12  &
       , 9.64874015137022E-13, -4.90190576608710E-14, -1.77253430678112E-14  &
       , 5.55950610442662E-15, -7.11793337579530E-16 /),  &
        fpi12 = 1.30899693899575E+00,  &
        con2 = 5.03154716196777E+00, con3 = 3.80004589867293E-01,  &
        con4 = 8.33333333333333E-01, con5 = 8.66025403784439E-01
INTEGER, PARAMETER :: n1 = 14, n2 = 23, n3 = 19, n4 = 15, m1 = 12, m2 = 21, &
                      m3 = 17, m4 = 13, n1d = 14, n2d = 24, n3d = 19,  &
                      n4d = 15, m1d = 12, m2d = 22, m3d = 17, m4d = 13
REAL     :: ccv, cv, ec, e1, e2, f1, f2, rtrx, scv, t, temp1, temp2, tt
INTEGER  :: i, j
!     -------------------
IF (x >= 0.) THEN
  IF (c <= 5.) THEN
    IF (x <= 1.2) THEN
      t = (x+x-1.2) * con4
      tt = t + t
      j = n1
      f1 = ak1(j)
      f2 = 0.
      DO i = 1, m1
        j = j - 1
        temp1 = f1
        f1 = tt * f1 - f2 + ak1(j)
        f2 = temp1
      END DO
      ai = t * f1 - f2 + ak1(1)

      j = n1d
      f1 = dak1(j)
      f2 = 0.
      DO i = 1, m1d
        j = j - 1
        temp1 = f1
        f1 = tt * f1 - f2 + dak1(j)
        f2 = temp1
      END DO
      dai = -(t*f1 - f2 + dak1(1))
      RETURN
    END IF

    t = (x + x - con2) * con3
    tt = t + t
    j = n2
    f1 = ak2(j)
    f2 = 0.
    DO i = 1, m2
      j = j - 1
      temp1 = f1
      f1 = tt * f1 - f2 + ak2(j)
      f2 = temp1
    END DO
    rtrx = SQRT(rx)
    ec = EXP(-c)
    ai = ec * (t*f1-f2+ak2(1)) / rtrx
    j = n2d
    f1 = dak2(j)
    f2 = 0.
    DO i = 1, m2d
      j = j - 1
      temp1 = f1
      f1 = tt * f1 - f2 + dak2(j)
      f2 = temp1
    END DO
    dai = -ec * (t*f1-f2+dak2(1)) * rtrx
    RETURN
  END IF

  t = 10. / c - 1.
  tt = t + t
  j = n1
  f1 = ak3(j)
  f2 = 0.
  DO i = 1, m1
    j = j - 1
    temp1 = f1
    f1 = tt * f1 - f2 + ak3(j)
    f2 = temp1
  END DO
  rtrx = SQRT(rx)
  ec = EXP(-c)
  ai = ec * (t*f1-f2+ak3(1)) / rtrx
  j = n1d
  f1 = dak3(j)
  f2 = 0.
  DO i = 1, m1d
    j = j - 1
    temp1 = f1
    f1 = tt * f1 - f2 + dak3(j)
    f2 = temp1
  END DO
  dai = -rtrx * ec * (t*f1-f2+dak3(1))
  RETURN
END IF

IF (c <= 5.) THEN
  t = .4 * c - 1.
  tt = t + t
  j = n3
  f1 = ajp(j)
  e1 = ajn(j)
  f2 = 0.
  e2 = 0.
  DO i = 1, m3
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt * f1 - f2 + ajp(j)
    e1 = tt * e1 - e2 + ajn(j)
    f2 = temp1
    e2 = temp2
  END DO
  ai = (t*e1-e2+ajn(1)) - x * (t*f1-f2+ajp(1))
  j = n3d
  f1 = dajp(j)
  e1 = dajn(j)
  f2 = 0.
  e2 = 0.
  DO i = 1, m3d
    j = j - 1
    temp1 = f1
    temp2 = e1
    f1 = tt * f1 - f2 + dajp(j)
    e1 = tt * e1 - e2 + dajn(j)
    f2 = temp1
    e2 = temp2
  END DO
  dai = x * x * (t*f1-f2+dajp(1)) + (t*e1-e2+dajn(1))
  RETURN
END IF

t = 10. / c - 1.
tt = t + t
j = n4
f1 = a(j)
e1 = b(j)
f2 = 0.
e2 = 0.
DO i = 1, m4
  j = j - 1
  temp1 = f1
  temp2 = e1
  f1 = tt * f1 - f2 + a(j)
  e1 = tt * e1 - e2 + b(j)
  f2 = temp1
  e2 = temp2
END DO
temp1 = t * f1 - f2 + a(1)
temp2 = t * e1 - e2 + b(1)
rtrx = SQRT(rx)
cv = c - fpi12
ccv = COS(cv)
scv = SIN(cv)
ai = (temp1*ccv - temp2*scv) / rtrx
j = n4d
f1 = da(j)
e1 = db(j)
f2 = 0.
e2 = 0.
DO i = 1, m4d
  j = j - 1
  temp1 = f1
  temp2 = e1
  f1 = tt * f1 - f2 + da(j)
  e1 = tt * e1 - e2 + db(j)
  f2 = temp1
  e2 = temp2
END DO
temp1 = t * f1 - f2 + da(1)
temp2 = t * e1 - e2 + db(1)
e1 = ccv * con5 + .5 * scv
e2 = scv * con5 - .5 * ccv
dai = (temp1*e1-temp2*e2) * rtrx
RETURN
END SUBROUTINE jairy


SUBROUTINE bssly(a, in, w)
!     ******************************************************************
!     FORTRAN SUBROUTINE FOR ORDINARY BESSEL FUNCTION OF INTEGRAL ORDER
!     ******************************************************************
!     A  = ARGUMENT (COMPLEX NUMBER)
!     IN = ORDER (INTEGER)
!     W  = FUNCTION OF SECOND KIND (COMPLEX NUMBER)
!     -------------------
COMPLEX, INTENT(IN)  :: a
INTEGER, INTENT(IN)  :: in
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL    :: an, az(2), aq(2), pm, pn, qm, qn, rm, rn, qz(2), rz(2), sn, ss,  &
           sz(2), ts(2), tm(4), sm(4), sl(2), sq(2), sr(2), qf(2), zl(2), zs
INTEGER :: i, la, lr, ly, m, n
REAL, PARAMETER :: cd(30) = (/ 0.00000000000000, -1.64899505142212E-2,  &
        -7.18621880068536E-2, -1.67086878124866E-1, -3.02582250219469E-1,  &
        -4.80613945245927E-1, -7.07075239357898E-1, -9.92995790539516E-1,  &
        -1.35583925612592, -1.82105907899132, -2.42482175310879,  &
        -3.21956655708750, -4.28658077248384, -5.77022816798128,  &
        -8.01371260952526,  0.00000000000000, -5.57742429879505E-3,  &
        -4.99112944172476E-2, -1.37440911652397E-1, -2.67233784710566E-1,  &
        -4.40380166808682E-1, -6.61813614872541E-1, -9.41861077665017E-1,  &
        -1.29754130468326, -1.75407696719816, -2.34755299882276,  &
        -3.13041332689196, -4.18397120563729, -5.65251799214994,  &
        -7.87863959810677 /),   &
        ce(30) = (/ 0.00000000000000,  &
        -4.80942336387447E-3, -1.31366200347759E-2, -1.94843834008458E-2,  &
        -2.19948900032003E-2, -2.09396625676519E-2, -1.74600268458650E-2,  &
        -1.27937813362085E-2, -8.05234421796592E-3, -4.15817375002760E-3,  &
        -1.64317738747922E-3, -4.49175585314709E-4, -7.28594765574007E-5,  &
        -5.38265230658285E-6, -9.93779048036289E-8, 0.00000000000000,  &
         7.53805779200591E-2,  7.12293537403464E-2, 6.33116224228200E-2,  &
         5.28240264523301E-2,  4.13305359441492E-2, 3.01350573947510E-2,  &
         2.01043439592720E-2,  1.18552223068074E-2, 5.86055510956010E-3,  &
         2.25465148267325E-3,  6.08173041536336E-4, 9.84215550625747E-5,  &
         7.32139093038089E-6,  1.37279667384666E-7 /)

!     -------------------
az(1) = REAL(a)
az(2) = AIMAG(a)
zs = az(1) * az(1) + az(2) * az(2)
zl(1) = 0.5 * LOG(zs)
zl(2) = ATAN2(az(2),az(1))
an = ABS(in)
sn = +1.0
IF (in < 0) THEN
  IF (in /= in/2*2) THEN
    sn = -1.0
  END IF
END IF
IF (az(1) < 0) THEN
  qz(1) = -az(1)
  qz(2) = -az(2)
ELSE
  qz(1) = +az(1)
  qz(2) = +az(2)
END IF
IF (zs <= 1.0) GO TO 100
IF (zs < 289.0) THEN
  IF (-ABS(az(2))+0.096*az(1)*az(1) > 0.0) THEN
    GO TO 50
  ELSE
    GO TO 100
  END IF
END IF
qm = sn * 0.797884560802865 * EXP(-0.5*zl(1))
qf(1) = qm * COS(-0.5*zl(2))
qf(2) = qm * SIN(-0.5*zl(2))
IF (an > 1.0) GO TO 20
pn = an
la = 10
GO TO 190

10 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)
GO TO 180
20 pn = 1.0
la = 30
GO TO 190

30 sq(1) = -qf(1) * sm(1) + qf(2) * sm(2)
sq(2) = -qf(1) * sm(2) - qf(2) * sm(1)
pn = 0.0
la = 40
GO TO 190

40 sr(1) = +qf(1) * sm(1) - qf(2) * sm(2)
sr(2) = +qf(1) * sm(2) + qf(2) * sm(1)
GO TO 150
50 qm = sn * 0.3989422804014327 * EXP(-0.5*zl(1))
qf(1) = qm * COS(-0.5*zl(2))
qf(2) = qm * SIN(-0.5*zl(2))
IF (an > 1.0) GO TO 70
pn = an
lr = 60
GO TO 240

60 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)
GO TO 180
70 pn = 1.0
lr = 80
GO TO 240

80 sq(1) = -qf(1) * sm(1) + qf(2) * sm(2)
sq(2) = -qf(1) * sm(2) - qf(2) * sm(1)
pn = 0.0
lr = 90
GO TO 240

90 sr(1) = +qf(1) * sm(1) - qf(2) * sm(2)
sr(2) = +qf(1) * sm(2) + qf(2) * sm(1)
GO TO 150
100 qf(1) = sn * 0.6366197723675813
qf(2) = 0.0
IF (an > 1.0) GO TO 120
pn = an
ly = 110
GO TO 270

110 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)
GO TO 180
120 pn = 1.0
ly = 130
GO TO 270

130 sq(1) = -qf(1) * sm(1) + qf(2) * sm(2)
sq(2) = -qf(1) * sm(2) - qf(2) * sm(1)
pn = 0.0
ly = 140
GO TO 270

140 sr(1) = +qf(1) * sm(1) - qf(2) * sm(2)
sr(2) = +qf(1) * sm(2) + qf(2) * sm(1)
150 rz(1) = +az(1) / zs
rz(2) = -az(2) / zs
pn = 0.0
GO TO 170

160 sq(1) = sr(1)
sq(2) = sr(2)
sr(1) = sm(1)
sr(2) = sm(2)
170 sm(1) = 2.0 * pn * (rz(1)*sr(1) - rz(2)*sr(2)) - sq(1)
sm(2) = 2.0 * pn * (rz(1)*sr(2) + rz(2)*sr(1)) - sq(2)
pn = pn + 1.0
IF (pn < an) GO TO 160
180 w = CMPLX(sm(1), sm(2))
RETURN

190 sm(1) = 0.0
sm(2) = 0.0
sm(3) = 0.0
sm(4) = 0.0
rz(1) = +0.5 * qz(1) / zs
rz(2) = -0.5 * qz(2) / zs
qn = pn * pn - 0.25
tm(1) = 1.0
tm(2) = 0.0
pm = 0.0
GO TO 210

200 qn = qn - 2.0 * pm
pm = pm + 1.0
ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
tm(1) = -qn * ts(1) / pm
tm(2) = -qn * ts(2) / pm
210 sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
qn = qn - 2.0 * pm
pm = pm + 1.0
ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
tm(1) = +qn * ts(1) / pm
tm(2) = +qn * ts(2) / pm
IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
  IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 220
END IF
sm(3) = sm(3) + tm(1)
sm(4) = sm(4) + tm(2)
IF (pm < 35.0) GO TO 200
220 aq(1) = qz(1) - 1.57079632679490 * (pn+0.5)
aq(2) = qz(2)
ts(1) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
ts(2) = -SIN(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
tm(1) = sm(1) * ts(1) - sm(2) * ts(2)
tm(2) = sm(1) * ts(2) + sm(2) * ts(1)
tm(3) = sm(3) * ts(1) - sm(4) * ts(2)
tm(4) = sm(3) * ts(2) + sm(4) * ts(1)
ts(1) = +SIN(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
ts(2) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
tm(1) = tm(1) - sm(3) * ts(1) + sm(4) * ts(2)
tm(2) = tm(2) - sm(3) * ts(2) - sm(4) * ts(1)
tm(3) = tm(3) + sm(1) * ts(1) - sm(2) * ts(2)
tm(4) = tm(4) + sm(1) * ts(2) + sm(2) * ts(1)
IF (az(1) < 0) THEN
  IF (az(2) < 0) THEN
    sm(1) = -2.0 * tm(1) + tm(4)
    sm(2) = -2.0 * tm(2) - tm(3)
  ELSE
    sm(1) = -2.0 * tm(1) - tm(4)
    sm(2) = -2.0 * tm(2) + tm(3)
  END IF
  IF (pn == 0.0) GO TO 230
  sm(1) = -sm(1)
  sm(2) = -sm(2)
ELSE
  sm(1) = tm(3)
  sm(2) = tm(4)
END IF

230 SELECT CASE (la)
  CASE (10)
    GO TO 10
  CASE (30)
    GO TO 30
  CASE (40)
    GO TO 40
END SELECT

240 sm(1) = 1.0
sm(2) = 0.0
sm(3) = 1.0
sm(4) = 0.0
m = 15.0 * pn + 2.0
n = 15.0 * pn + 15.0
DO i = m, n
  ts(1) = +qz(2) - cd(i)
  ts(2) = -qz(1)
  ss = ts(1) * ts(1) + ts(2) * ts(2)
  tm(1) = +ce(i) * ts(1) / ss
  tm(2) = -ce(i) * ts(2) / ss
  sm(1) = sm(1) + tm(1)
  sm(2) = sm(2) + tm(2)
  ts(1) = -qz(2) - cd(i)
  ts(2) = +qz(1)
  ss = ts(1) * ts(1) + ts(2) * ts(2)
  tm(1) = +ce(i) * ts(1) / ss
  tm(2) = -ce(i) * ts(2) / ss
  sm(3) = sm(3) + tm(1)
  sm(4) = sm(4) + tm(2)
END DO
aq(1) = qz(1) - 1.57079632679490 * (pn+0.5)
aq(2) = qz(2)
ts(1) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
ts(2) = -SIN(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
tm(1) = +ts(1) * sm(1) - ts(2) * sm(2) + ts(1) * sm(3) - ts(2) *sm(4)
tm(2) = +ts(1) * sm(2) + ts(2) * sm(1) + ts(1) * sm(4) + ts(2) *sm(3)
tm(3) = +ts(1) * sm(2) + ts(2) * sm(1) - ts(1) * sm(4) - ts(2) *sm(3)
tm(4) = -ts(1) * sm(1) + ts(2) * sm(2) + ts(1) * sm(3) - ts(2) *sm(4)
ts(1) = +SIN(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
ts(2) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
tm(1) = tm(1) - ts(1) * sm(2) - ts(2) * sm(1) + ts(1) * sm(4) +ts(2) * sm(3)
tm(2) = tm(2) + ts(1) * sm(1) - ts(2) * sm(2) - ts(1) * sm(3) +ts(2) * sm(4)
tm(3) = tm(3) + ts(1) * sm(1) - ts(2) * sm(2) + ts(1) * sm(3) -ts(2) * sm(4)
tm(4) = tm(4) + ts(1) * sm(2) + ts(2) * sm(1) + ts(1) * sm(4) +ts(2) * sm(3)
IF (az(1) < 0) THEN
  IF (az(2) < 0) THEN
    sm(1) = -2.0 * tm(1) + tm(4)
    sm(2) = -2.0 * tm(2) - tm(3)
  ELSE
    sm(1) = -2.0 * tm(1) - tm(4)
    sm(2) = -2.0 * tm(2) + tm(3)
  END IF
  IF (pn == 0.0) GO TO 260
  sm(1) = -sm(1)
  sm(2) = -sm(2)
ELSE
  sm(1) = tm(3)
  sm(2) = tm(4)
END IF
260 SELECT CASE (lr)
  CASE (60)
    GO TO 60
  CASE (80)
    GO TO 80
  CASE (90)
    GO TO 90
END SELECT

270 aq(1) = 1.0
aq(2) = 0.0
rn = 0.0
pm = 0.0
GO TO 290
280 pm = pm + 1.0
rn = rn + 0.5 / pm
ts(1) = 0.5 * (az(1)*aq(1)-az(2)*aq(2))
ts(2) = 0.5 * (az(1)*aq(2)+az(2)*aq(1))
aq(1) = ts(1) / pm
aq(2) = ts(2) / pm
290 IF (pm < pn) GO TO 280
sz(1) = 0.25 * (az(1)-az(2)) * (az(1)+az(2))
sz(2) = 0.5 * az(1) * az(2)
sr(1) = 0.0
sr(2) = 0.0
ss = aq(1) * aq(1) + aq(2) * aq(2)
tm(1) = +aq(1) / ss
tm(2) = -aq(2) / ss
pm = 0.0
GO TO 310

300 tm(1) = tm(1) / (pn-pm)
tm(2) = tm(2) / (pn-pm)
sr(1) = sr(1) - 0.5 * tm(1)
sr(2) = sr(2) - 0.5 * tm(2)
pm = pm + 1.0
ts(1) = sz(1) * tm(1) - sz(2) * tm(2)
ts(2) = sz(1) * tm(2) + sz(2) * tm(1)
tm(1) = +ts(1) / pm
tm(2) = +ts(2) / pm
310 IF (pm < pn) GO TO 300
sm(1) = 0.0
sm(2) = 0.0
rm = 1.0
qm = 0.0
sl(1) = -0.115931515658412 + zl(1) - rn
sl(2) = +zl(2)
pm = 0.0
GO TO 330

320 qm = qm + rm
pm = pm + 1.0
rm = 0.25 * zs * rm / (pm*(pn+pm))
ts(1) = sz(1) * aq(1) - sz(2) * aq(2)
ts(2) = sz(1) * aq(2) + sz(2) * aq(1)
aq(1) = -ts(1) / (pm*(pn+pm))
aq(2) = -ts(2) / (pm*(pn+pm))
sl(1) = sl(1) - 0.5 / pm - 0.5 / (pn+pm)
330 tm(1) = aq(1) * sl(1) - aq(2) * sl(2)
tm(2) = aq(1) * sl(2) + aq(2) * sl(1)
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (qm+rm > qm) GO TO 320
sm(1) = sr(1) + sm(1)
sm(2) = sr(2) + sm(2)
SELECT CASE (ly)
  CASE (110)
    GO TO 110
  CASE (130)
    GO TO 130
  CASE (140)
    GO TO 140
END SELECT

RETURN
END SUBROUTINE bssly



SUBROUTINE cbssli(z, cnu, w)
!------------------------------------------------------------
!     CALCULATION OF THE MODIFIED BESSEL FUNCTION OF THE
!     FIRST KIND FOR COMPLEX ORDER CNU AND COMPLEX ARGUMENT Z.
!     IT IS ASSUMED THAT -PI < ARG Z <= PI.
!------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, cnu
COMPLEX, INTENT(OUT) :: w

! Local variables
COMPLEX         :: nu, z0
REAL            :: t, x, y
REAL, PARAMETER :: pihalf = 1.5707963267949

x = REAL(z)
y = AIMAG(z)
IF (y >= 0.0) THEN
  z0 = CMPLX(y,-x)
  nu = cnu
ELSE
  z0 = CMPLX(-y,x)
  nu = -cnu
END IF

t = -pihalf * AIMAG(nu)
IF (t <= exparg(1)) THEN
  w = (0.0,0.0)
  RETURN
END IF

CALL cbsslj(z0,cnu,w)
w = EXP(t) * w
t = 0.5 * REAL(nu)
w = w * CMPLX(cos1(t),sin1(t))
RETURN
END SUBROUTINE cbssli


SUBROUTINE bssli(mo, a, in, w)
!     ******************************************************************
!     FORTRAN SUBROUTINE FOR MODIFIED BESSEL FUNCTION OF INTEGRAL ORDER
!     ******************************************************************
!     MO = MODE OF OPERATION
!     A  = ARGUMENT (COMPLEX NUMBER)
!     IN = ORDER (INTEGER)
!     W  = FUNCTION OF FIRST KIND (COMPLEX NUMBER)
!     -------------------
COMPLEX, INTENT(IN)  :: a
INTEGER, INTENT(IN)  :: mo, in
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: cd(30) = (/ 0.00000000000000,  &
        -1.64899505142212E-2, -7.18621880068536E-2, -1.67086878124866E-1,  &
        -3.02582250219469E-1, -4.80613945245927E-1, -7.07075239357898E-1,  &
        -9.92995790539516E-1, -1.35583925612592, -1.82105907899132,  &
        -2.42482175310879, -3.21956655708750, -4.28658077248384,  &
        -5.77022816798128, -8.01371260952526, 0.00000000000000,  &
        -5.57742429879505E-3, -4.99112944172476E-2, -1.37440911652397E-1,  &
        -2.67233784710566E-1, -4.40380166808682E-1, -6.61813614872541E-1,  &
        -9.41861077665017E-1, -1.29754130468326, -1.75407696719816,  &
        -2.34755299882276, -3.13041332689196, -4.18397120563729,  &
        -5.65251799214994, -7.87863959810677 /), ce(30)  &
        = (/ 0.00000000000000, -4.80942336387447E-3,  &
        -1.31366200347759E-2, -1.94843834008458E-2, -2.19948900032003E-2,  &
        -2.09396625676519E-2, -1.74600268458650E-2, -1.27937813362085E-2,  &
        -8.05234421796592E-3, -4.15817375002760E-3, -1.64317738747922E-3,  &
        -4.49175585314709E-4, -7.28594765574007E-5, -5.38265230658285E-6,  &
        -9.93779048036289E-8, 0.00000000000000, 7.53805779200591E-2,  &
        7.12293537403464E-2, 6.33116224228200E-2, 5.28240264523301E-2,  &
        4.13305359441492E-2, 3.01350573947510E-2, 2.01043439592720E-2,  &
        1.18552223068074E-2, 5.86055510956010E-3, 2.25465148267325E-3,  &
        6.08173041536336E-4, 9.84215550625747E-5, 7.32139093038089E-6,  &
        1.37279667384666E-7 /)
REAL    :: an, aq(2), az(2), fi(2), pm, pn, qm, qn, qz(2), rz(2), sz(2), &
           zr(2), ts(2), tm(2), rm(4), sm(4), sn, ss, qf(2), zm, zs
INTEGER :: i, m, n
!     -------------------
az(1) = REAL(a)
az(2) = AIMAG(a)
zs = az(1) * az(1) + az(2) * az(2)
zm = SQRT(zs)
pn = ABS(in)
sn = +1.0
IF (az(1) < 0) THEN
  qz(1) = -az(1)
  qz(2) = -az(2)
  IF (in == in/2*2) GO TO 10
  sn = -1.0
ELSE
  qz(1) = az(1)
  qz(2) = az(2)
END IF
10 IF (zm > 17.5+0.5*pn*pn) THEN
  qn = pn
ELSE
  qn = 0.5 * zm - 0.5 * ABS(qz(1)) + 0.5 *ABS(0.5*zm-ABS(qz(1)))
  IF (pn > qn) THEN
    qn = +AINT(0.0625*zs)
    IF (pn <= qn) GO TO 130
    qn = pn
    GO TO 130
  END IF
  IF (zm > 17.5) THEN
    qn = +AINT(SQRT(2.0*(zm-17.5)))
  ELSE
    IF (zs >= 1.0) THEN
      IF (-ABS(az(1))+0.096*az(2)*az(2) >= 0) GO TO 20
    END IF
    qn = AINT(0.0625*zs)
    IF (pn <= qn) GO TO 130
    qn = pn
    GO TO 130
    20 qn = 0.0
  END IF
END IF
sz(1) = qz(1)
sz(2) = qz(2)
qm = sn * 0.398942280401433
zr(1) = SQRT(sz(1)+zm)
zr(2) = sz(2) / zr(1)
zr(1) = 0.707106781186548 * zr(1)
zr(2) = 0.707106781186548 * zr(2)
qf(1) = +qm * zr(1) / zm
qf(2) = -qm * zr(2) / zm
IF (zm > 17.5) THEN
  rz(1) = +0.5 * qz(1) / zs
  rz(2) = -0.5 * qz(2) / zs
  an = qn * qn - 0.25
  sm(1) = 0.0
  sm(2) = 0.0
  sm(3) = 0.0
  sm(4) = 0.0
  tm(1) = 1.0
  tm(2) = 0.0
  pm = 0.0
  GO TO 40
  30 an = an - 2.0 * pm
  pm = pm + 1.0
  ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
  ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
  tm(1) = an * ts(1) / pm
  tm(2) = an * ts(2) / pm
  40 sm(1) = sm(1) + tm(1)
  sm(2) = sm(2) + tm(2)
  an = an - 2.0 * pm
  pm = pm + 1.0
  ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
  ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
  tm(1) = an * ts(1) / pm
  tm(2) = an * ts(2) / pm
  IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
    IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 50
  END IF
  sm(3) = sm(3) + tm(1)
  sm(4) = sm(4) + tm(2)
  IF (pm < 35.0) GO TO 30
  50 ts(1) = sm(1) + sm(3)
  ts(2) = sm(2) + sm(4)
  sm(1) = sm(1) - sm(3)
  sm(2) = sm(2) - sm(4)
  sm(3) = ts(1)
  sm(4) = ts(2)
ELSE
  sm(1) = 1.0
  sm(2) = 0.0
  sm(3) = 1.0
  sm(4) = 0.0
  m = 15.0 * qn + 2.0
  n = 15.0 * qn + 15.0
  DO i = m, n
    ts(1) = -qz(1) - cd(i)
    ts(2) = -qz(2)
    ss = ts(1) * ts(1) + ts(2) * ts(2)
    tm(1) = +ce(i) * ts(1) / ss
    tm(2) = -ce(i) * ts(2) / ss
    sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
    ts(1) = qz(1) - cd(i)
    ts(2) = qz(2)
    ss = ts(1) * ts(1) + ts(2) * ts(2)
    tm(1) = +ce(i) * ts(1) / ss
    tm(2) = -ce(i) * ts(2) / ss
    sm(3) = sm(3) + tm(1)
    sm(4) = sm(4) + tm(2)
  END DO
END IF
rm(1) = sm(1)
rm(2) = sm(2)
IF (qz(1) < 17.5) THEN
  aq(1) = -2.0 * qz(1)
  IF (qz(2) < 0) THEN
    aq(2) = -2.0 * qz(2) - 3.14159265358979 * (qn+0.5)
  ELSE
    aq(2) = -2.0 * qz(2) + 3.14159265358979 * (qn+0.5)
  END IF
  qm = EXP(aq(1))
  ts(1) = qm * COS(aq(2))
  ts(2) = qm * SIN(aq(2))
  rm(1) = rm(1) + ts(1) * sm(3) - ts(2) * sm(4)
  rm(2) = rm(2) + ts(1) * sm(4) + ts(2) * sm(3)
END IF
IF (qn /= pn) THEN
  rm(3) = rm(1)
  rm(4) = rm(2)
  qn = qn + 1.0
  IF (zm > 17.5) THEN
    an = qn * qn - 0.25
    sm(1) = 0.0
    sm(2) = 0.0
    sm(3) = 0.0
    sm(4) = 0.0
    tm(1) = 1.0
    tm(2) = 0.0
    pm = 0.0
    GO TO 80
    70 an = an - 2.0 * pm
    pm = pm + 1.0
    ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
    ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
    tm(1) = an * ts(1) / pm
    tm(2) = an * ts(2) / pm
    80 sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
    an = an - 2.0 * pm
    pm = pm + 1.0
    ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
    ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
    tm(1) = an * ts(1) / pm
    tm(2) = an * ts(2) / pm
    IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
      IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 90
    END IF
    sm(3) = sm(3) + tm(1)
    sm(4) = sm(4) + tm(2)
    IF (pm < 35.0) GO TO 70
    90 ts(1) = sm(1) + sm(3)
    ts(2) = sm(2) + sm(4)
    sm(1) = sm(1) - sm(3)
    sm(2) = sm(2) - sm(4)
    sm(3) = ts(1)
    sm(4) = ts(2)
  ELSE
    sm(1) = 1.0
    sm(2) = 0.0
    sm(3) = 1.0
    sm(4) = 0.0
    m = 15.0 * qn + 2.0
    n = 15.0 * qn + 15.0
    DO i = m, n
      ts(1) = -qz(1) - cd(i)
      ts(2) = -qz(2)
      ss = ts(1) * ts(1) + ts(2) * ts(2)
      tm(1) = +ce(i) * ts(1) / ss
      tm(2) = -ce(i) * ts(2) / ss
      sm(1) = sm(1) + tm(1)
      sm(2) = sm(2) + tm(2)
      ts(1) = +qz(1) - cd(i)
      ts(2) = +qz(2)
      ss = ts(1) * ts(1) + ts(2) * ts(2)
      tm(1) = +ce(i) * ts(1) / ss
      tm(2) = -ce(i) * ts(2) / ss
      sm(3) = sm(3) + tm(1)
      sm(4) = sm(4) + tm(2)
    END DO
  END IF
  rm(1) = sm(1)
  rm(2) = sm(2)
  IF (qz(1) >= 17.5) GO TO 120
  aq(1) = -2.0 * qz(1)
  IF (qz(2) < 0) THEN
    aq(2) = -2.0 * qz(2) - 3.14159265358979 * (qn+0.5)
  ELSE
    aq(2) = -2.0 * qz(2) + 3.14159265358979 * (qn+0.5)
  END IF
  qm = EXP(aq(1))
  ts(1) = qm * COS(aq(2))
  ts(2) = qm * SIN(aq(2))
  rm(1) = rm(1) + ts(1) * sm(3) - ts(2) * sm(4)
  rm(2) = rm(2) + ts(1) * sm(4) + ts(2) * sm(3)
  GO TO 120
  110 tm(1) = -2.0 * qn * qz(1) / zs
  tm(2) = +2.0 * qn * qz(2) / zs
  ts(1) = tm(1) * rm(1) - tm(2) * rm(2) + rm(3)
  ts(2) = tm(1) * rm(2) + tm(2) * rm(1) + rm(4)
  rm(3) = rm(1)
  rm(4) = rm(2)
  rm(1) = ts(1)
  rm(2) = ts(2)
  qn = qn + 1.0
  120 IF (qn < pn) GO TO 110
END IF
IF (mo == 0) THEN
  qm = EXP(qz(1))
  tm(1) = qm * COS(qz(2))
  tm(2) = qm * SIN(qz(2))
  ts(1) = tm(1) * rm(1) - tm(2) * rm(2)
  ts(2) = tm(1) * rm(2) + tm(2) * rm(1)
  rm(1) = ts(1)
  rm(2) = ts(2)
END IF
fi(1) = qf(1) * rm(1) - qf(2) * rm(2)
fi(2) = qf(1) * rm(2) + qf(2) * rm(1)
w = CMPLX(fi(1),fi(2))
RETURN
130 sz(1) = 0.25 * (qz(1)*qz(1)-qz(2)*qz(2))
sz(2) = 0.5 * qz(1) * qz(2)
an = qn
sm(1) = 0.0
sm(2) = 0.0
sm(3) = 0.0
sm(4) = 0.0
tm(1) = 1.0
tm(2) = 0.0
pm = 0.0
140 an = an + 1.0
ts(1) = tm(1) / an
ts(2) = tm(2) / an
sm(3) = sm(3) + ts(1)
sm(4) = sm(4) + ts(2)
tm(1) = ts(1) * sz(1) - ts(2) * sz(2)
tm(2) = ts(1) * sz(2) + ts(2) * sz(1)
pm = pm + 1.0
tm(1) = tm(1) / pm
tm(2) = tm(2) / pm
IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 150
END IF
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
GO TO 140
150 sm(1) = sm(1) + 1.0
an = qn + 1.0
sm(3) = an * sm(3)
sm(4) = an * sm(4)
GO TO 170
160 an = qn * (qn+1.0)
tm(1) = sz(1) / an
tm(2) = sz(2) / an
ts(1) = +tm(1) * sm(3) - tm(2) * sm(4)
ts(2) = +tm(1) * sm(4) + tm(2) * sm(3)
sm(3) = sm(1)
sm(4) = sm(2)
sm(1) = sm(1) + ts(1)
sm(2) = sm(2) + ts(2)
qn = qn - 1.0
170 IF (qn > pn) GO TO 160
qf(1) = sn
qf(2) = 0.0
qn = 0.0
GO TO 190
180 qn = qn + 1.0
tm(1) = qf(1) * qz(1) - qf(2) * qz(2)
tm(2) = qf(1) * qz(2) + qf(2) * qz(1)
qf(1) = 0.5 * tm(1) / qn
qf(2) = 0.5 * tm(2) / qn
190 IF (qn < pn) GO TO 180
IF (mo /= 0) THEN
  qm = EXP(-qz(1))
  tm(1) = qm * COS(-qz(2))
  tm(2) = qm * SIN(-qz(2))
  ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
  ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
  qf(1) = ts(1)
  qf(2) = ts(2)
END IF
fi(1) = qf(1) * sm(1) - qf(2) * sm(2)
fi(2) = qf(1) * sm(2) + qf(2) * sm(1)
w = CMPLX(fi(1),fi(2))
RETURN
END SUBROUTINE bssli



SUBROUTINE besi(x, alpha, kode, n, y, nz)

!     WRITTEN BY D. E. AMOS AND S. L. DANIEL, JANUARY,1975.

!     REFERENCE
!         SAND-75-0152

!         CDC 6600 SUBROUTINES IBESS AND JBESS FOR BESSEL FUNCTIONS
!         I(NU,X) AND J(NU,X), X >= 0, NU >= 0  BY D.E. AMOS, S.L.
!         DANIEL, M.K. WESTON. ACM TRANS MATH SOFTWARE,3,PP 76-92 (1977)

!         TABLES OF BESSEL FUNCTIONS OF MODERATE OR LARGE ORDERS,
!         NPL MATHEMATICAL TABLES, VOL. 6, BY F.W.J. OLVER, HER
!         MAJESTY'S STATIONERY OFFICE, LONDON, 1962.

!     ABSTRACT
!         BESI COMPUTES AN N MEMBER SEQUENCE OF I BESSEL FUNCTIONS
!         I/SUB(ALPHA+K-1)/(X), K=1,...,N OR SCALED BESSEL FUNCTIONS
!         EXP(-X)*I/SUB(ALPHA+K-1)/(X), K=1,...,N FOR NON-NEGATIVE ALPHA
!         AND X. A COMBINATION OF THE POWER SERIES, THE ASYMPTOTIC
!         EXPANSION FOR X TO INFINITY, AND THE UNIFORM ASYMPTOTIC
!         EXPANSION FOR NU TO INFINITY ARE APPLIED OVER SUBDIVISIONS OF
!         THE (NU,X) PLANE. FOR VALUES NOT COVERED BY ONE OF THESE
!         FORMULAE, THE ORDER IS INCREMENTED BY AN INTEGER SO THAT ONE
!         OF THESE FORMULAE APPLY. BACKWARD RECURSION IS USED TO REDUCE
!         ORDERS BY INTEGER VALUES. THE ASYMPTOTIC EXPANSION FOR X TO
!         INFINITY IS USED ONLY WHEN THE ENTIRE SEQUENCE (SPECIFICALLY
!         THE LAST MEMBER) LIES WITHIN THE REGION COVERED BY THE
!         EXPANSION. LEADING TERMS OF THESE EXPANSIONS ARE USED TO TEST
!         FOR OVER OR UNDERFLOW WHERE APPROPRIATE. IF A SEQUENCE IS
!         REQUESTED AND THE LAST MEMBER WOULD UNDERFLOW, THE RESULT IS
!         SET TO ZERO AND THE NEXT LOWER ORDER TRIED, ETC., UNTIL A
!         MEMBER COMES ON SCALE OR ALL ARE SET TO ZERO. AN OVERFLOW
!         CANNOT OCCUR WITH SCALING.

!         BESI CALLS ASIK, GAMLN, SPMPAR, AND IPMPAR

!     DESCRIPTION OF ARGUMENTS

!         INPUT
!           X      - X >= 0.0
!           ALPHA  - ORDER OF FIRST MEMBER OF THE SEQUENCE,
!                    ALPHA >= 0.0
!           KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!                    KODE=1 RETURNS
!                           Y(K)=        I/SUB(ALPHA+K-1)/(X),
!                                K=1,...,N
!                    KODE=2 RETURNS
!                           Y(K)=EXP(-X)*I/SUB(ALPHA+K-1)/(X),
!                                K=1,...,N
!           N      - NUMBER OF MEMBERS IN THE SEQUENCE, N>=1

!         OUTPUT
!           Y      - A VECTOR WHOSE FIRST N COMPONENTS CONTAIN
!                    VALUES FOR I/SUB(ALPHA+K-1)/(X) OR SCALED
!                    VALUES FOR EXP(-X)*I/SUB(ALPHA+K-1)/(X),
!                    K=1,...,N DEPENDING ON KODE
!           NZ     - ERROR INDICATOR
!                    NZ= 0     NORMAL RETURN-COMPUTATION COMPLETED
!                    NZ=-1     X IS LESS THAN 0.0
!                    NZ=-2     ALPHA IS LESS THAN 0.0
!                    NZ=-3     N IS LESS THAN 1
!                    NZ=-4     KODE IS NOT 1 OR 2
!                    NZ=-5     X IS TOO LARGE FOR KODE=1
!                    NZ>0   LAST NZ COMPONENTS OF Y SET TO 0.0
!                              BECAUSE OF UNDERFLOW

!     ERROR CONDITIONS
!         IMPROPER INPUT ARGUMENTS - A FATAL ERROR
!         OVERFLOW WITH KODE=1 - A FATAL ERROR
!         UNDERFLOW - A NON-FATAL ERROR(NZ>0)

REAL, INTENT(IN)     :: x, alpha
INTEGER, INTENT(IN)  :: kode, n
REAL, INTENT(OUT)    :: y(:)
INTEGER, INTENT(OUT) :: nz

! Local variables
INTEGER :: i, ialp, in, inlim = 80, is, i1, i2, k, kk, km, kt, nn, ns
REAL    :: ain, ak, akm, ans, ap, arg, atol, tolln, dfn, dtm, dx,  &
           earg, elim, etx, flgik, fn, fnf, fni, fnp1, fnu, gln, ra,  &
           rttpi = 3.98942280401433E-01, s, sx, sxo2, s1, s2, t, ta, tb,  &
           temp(3), tfn, tm, tol, trx, t2, xo2, xo2l, z
!     -------------------
!     IPMPAR(8) REPLACES IPMPAR(5) IN A REAL (dp) CODE
!     IPMPAR(9) REPLACES IPMPAR(6) IN A REAL (dp) CODE

!     DEFINITION OF THE TOLERANCES TOL AND ELIM

tb = RADIX(1.0)
ta = EPSILON(1.0) / tb
IF (tb /= 2.0) THEN
  IF (tb == 8.0) GO TO 10
  IF (tb == 16.0) GO TO 20
  tb = LOG(tb)
  GO TO 30
END IF
tb = .69315
GO TO 30
10 tb = 2.07944
GO TO 30
20 tb = 2.77259

30 tol = MAX(ta,1.e-15)
i1 = DIGITS(1.0)
i2 = MINEXPONENT(1.0)
!     LN(10**3) = 6.90776
elim = REAL(-i2) * tb - 6.90776
!     TOLLN = -LN(TOL)
tolln = REAL(i1) * tb
tolln = MIN(tolln,34.5388)



nz = 0
kt = 1
IF (n < 1) THEN
  GO TO 480
ELSE IF (n > 1) THEN
  GO TO 50
END IF
kt = 2
50 nn = n
IF (kode < 1.OR.kode > 2) GO TO 460

IF (x < 0.0) THEN
  GO TO 490
ELSE IF (x > 0.0) THEN
  GO TO 110
END IF

IF (alpha < 0.0) THEN
  GO TO 470
ELSE IF (x > 0.0) THEN
  GO TO 80
END IF
y(1) = 1.0
IF (n == 1) RETURN
i1 = 2
GO TO 90
80 i1 = 1
90 DO i = i1, n
  y(i) = 0.0
END DO
RETURN

110 IF (alpha < 0.0) GO TO 470

ialp = INT(alpha)
fni = REAL(ialp+n-1)
fnf = alpha - REAL(ialp)
dfn = fni + fnf
fnu = dfn
in = 0
xo2 = x * 0.5
sxo2 = xo2 * xo2
etx = REAL(kode-1)
sx = etx * x

!     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
!     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE APPLIED.

IF (sxo2 > (fnu+1.0)) THEN
  IF (x <= 12.0) GO TO 130
  fn = 0.55 * fnu * fnu
  fn = MAX(17.0,fn)
  IF (x >= fn) GO TO 360
  ans = MAX(36.0-fnu,0.0)
  ns = INT(ans)
  fni = fni + REAL(ns)
  dfn = fni + fnf
  fn = dfn
  is = kt
  km = n - 1 + ns
  IF (km > 0) is = 3
  GO TO 140
END IF
fn = fnu
fnp1 = fn + 1.0
xo2l = LOG(xo2)
is = kt
IF (x <= 0.5) GO TO 210
ns = 0
120 fni = fni + REAL(ns)
dfn = fni + fnf
fn = dfn
fnp1 = fn + 1.0
is = kt
IF (n-1+ns > 0) is = 3
GO TO 210
130 xo2l = LOG(xo2)
ns = INT(sxo2-fnu)
GO TO 120

!     OVERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION

140 IF (kode /= 2) THEN
  IF (alpha < 1.0) GO TO 170
  z = x / alpha
  ra = SQRT(1.0+z*z)
  gln = LOG((1.0+ra)/z)
  t = ra * (1.0-etx) + etx / (z+ra)
  arg = alpha * (t-gln)
  IF (arg > elim) GO TO 500
  IF (km == 0) GO TO 160
END IF

!     UNDERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION

150 z = x / fn
ra = SQRT(1.0+z*z)
gln = LOG((1.0+ra)/z)
t = ra * (1.0-etx) + etx / (z+ra)
arg = fn * (t-gln)
160 IF (arg < (-elim)) GO TO 260
GO TO 200
170 IF (x > elim) GO TO 500
GO TO 150

!     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY

180 IF (km == 0) THEN
  y(1) = temp(3)
  RETURN
END IF
temp(1) = temp(3)
in = ns
kt = 1
i1 = 0

190 is = 2
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
IF (i1 == 2) GO TO 320
z = x / fn
ra = SQRT(1.0+z*z)
gln = LOG((1.0+ra)/z)
t = ra * (1.0-etx) + etx / (z+ra)
arg = fn * (t-gln)

200 i1 = ABS(3-is)
i1 = MAX(i1,1)
flgik = 1.0
CALL asik(x, fn, kode, flgik, ra, arg, i1, tol, temp(is:))
SELECT CASE (is)
  CASE (1)
    GO TO 190
  CASE (2)
    GO TO 320
  CASE (3)
    GO TO 410
END SELECT

!     SERIES FOR (X/2)**2<=NU+1

210 gln = gamln(fnp1)
arg = fn * xo2l - gln - sx
IF (arg < (-elim)) GO TO 280
earg = EXP(arg)

220 s = 1.0
IF (x >= tol) THEN
  ak = 3.0
  t2 = 1.0
  t = 1.0
  s1 = fn
  DO k = 1, 17
    s2 = t2 + s1
    t = t * sxo2 / s2
    s = s + t
    IF (ABS(t) < tol) GO TO 240
    t2 = t2 + ak
    ak = ak + 2.0
    s1 = s1 + fn
  END DO
END IF

240 temp(is) = s * earg
SELECT CASE (is)
  CASE (1)
    GO TO 250
  CASE (2)
    GO TO 320
  CASE (3)
    GO TO 400
END SELECT

250 earg = earg * fn / xo2
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
is = 2
GO TO 220

!     SET UNDERFLOW VALUE AND UPDATE PARAMETERS

260 y(nn) = 0.0
nn = nn - 1
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
IF (nn < 1) THEN
  GO TO 310
ELSE IF (nn > 1) THEN
  GO TO 150
END IF
kt = 2
is = 2
GO TO 150

280 y(nn) = 0.0
nn = nn - 1
fnp1 = fn
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
IF (nn < 1) THEN
  GO TO 310
ELSE IF (nn > 1) THEN
  GO TO 300
END IF
kt = 2
is = 2
300 IF (sxo2 > fnp1) THEN
  GO TO 150
END IF
arg = arg - xo2l + LOG(fnp1)
IF (arg < (-elim)) GO TO 280
GO TO 210
310 nz = n - nn
RETURN

!     BACKWARD RECURSION SECTION

320 nz = n - nn
330 IF (kt /= 2) THEN
  s1 = temp(1)
  s2 = temp(2)
  trx = 2.0 / x
  dtm = fni
  tm = (dtm+fnf) * trx
  IF (in /= 0) THEN
!     BACKWARD RECUR TO INDEX ALPHA+NN-1
    DO i = 1, in
      s = s2
      s2 = tm * s2 + s1
      s1 = s
      dtm = dtm - 1.0
      tm = (dtm+fnf) * trx
    END DO
    y(nn) = s1
    IF (nn == 1) RETURN
    y(nn-1) = s2
    IF (nn == 2) RETURN
  ELSE
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
    y(nn) = s1
    y(nn-1) = s2
    IF (nn == 2) RETURN
  END IF
  k = nn + 1
  DO i = 3, nn
    k = k - 1
    y(k-2) = tm * y(k-1) + y(k)
    dtm = dtm - 1.0
    tm = (dtm+fnf) * trx
  END DO
  RETURN
END IF
y(1) = temp(2)
RETURN

!     ASYMPTOTIC EXPANSION FOR X TO INFINITY

360 earg = rttpi / SQRT(x)
IF (kode /= 2) THEN
  IF (x > elim) GO TO 500
  earg = earg * EXP(x)
END IF
etx = 8.0 * x
is = kt
in = 0
fn = fnu
370 dx = fni + fni
tm = 0.0
IF (fni /= 0.0.OR.ABS(fnf) >= tol) THEN
  tm = 4.0 * fnf * (fni+fni+fnf)
END IF
dtm = dx * dx
s1 = etx
trx = dtm - 1.0
dx = -(trx+tm) / etx
t = dx
s = 1.0 + dx
atol = tol * ABS(s)
s2 = 1.0
ak = 8.0
DO k = 1, 25
  s1 = s1 + etx
  s2 = s2 + ak
  dx = dtm - s2
  ap = dx + tm
  t = -t * ap / s1
  s = s + t
  IF (ABS(t) <= atol) GO TO 390
  ak = ak + 8.0
END DO
390 temp(is) = s * earg
IF (is == 2) GO TO 330
is = 2
fni = fni - 1.0
dfn = fni + fnf
fn = dfn
GO TO 370

!     BACKWARD RECURSION WITH NORMALIZATION BY
!     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.

!     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
400 akm = MAX(3.0-fn,0.0)
km = INT(akm)
tfn = fn + REAL(km)
ta = (gln + tfn - 0.9189385332 - 0.0833333333/tfn) / (tfn + 0.5)
ta = xo2l - ta
tb = -(1.0 - 1.0/tfn) / tfn
ain = tolln / (-ta + SQRT(ta*ta-tolln*tb)) + 1.5
in = INT(ain)
in = in + km
GO TO 420

!     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
410 t = 1.0 / (fn*ra)
ain = tolln / (gln + SQRT(gln*gln + t*tolln)) + 1.5
in = INT(ain)
IF (in > inlim) GO TO 180

420 trx = 2.0 / x
dtm = fni + REAL(in)
tm = (dtm+fnf) * trx
ta = 0.0
tb = tol
kk = 1

!     BACKWARD RECUR UNINDEXED

430 DO i = 1, in
  s = tb
  tb = tm * tb + ta
  ta = s
  dtm = dtm - 1.0
  tm = (dtm+fnf) * trx
END DO
!     NORMALIZATION
IF (kk == 1) THEN
  ta = (ta/tb) * temp(3)
  tb = temp(3)
  kk = 2
  in = ns
  IF (ns /= 0) GO TO 430
END IF
y(nn) = tb
nz = n - nn
IF (nn == 1) RETURN
tb = tm * tb + ta
k = nn - 1
y(k) = tb
IF (nn == 2) RETURN
dtm = dtm - 1.0
tm = (dtm+fnf) * trx
km = k - 1

!     BACKWARD RECUR INDEXED

DO i = 1, km
  y(k-1) = tm * y(k) + y(k+1)
  dtm = dtm - 1.0
  tm = (dtm+fnf) * trx
  k = k - 1
END DO
RETURN


460 nz = -4
RETURN

470 nz = -2
RETURN

480 nz = -3
RETURN

490 nz = -1
RETURN

500 nz = -5
RETURN
END SUBROUTINE besi


SUBROUTINE asik(x, fnu, kode, flgik, ra, arg, in, tol, y)

!                  ASIK COMPUTES BESSEL FUNCTIONS I AND K
!                  FOR ARGUMENTS X > 0.0 AND ORDERS FNU >= 35
!                  ON FLGIK = 1 AND FLGIK = -1 RESPECTIVELY.

!                                    INPUT

!      X    - ARGUMENT, X > 0.0
!      FNU  - ORDER OF FIRST BESSEL FUNCTION
!      KODE - A PARAMETER TO INDICATE THE SCALING OPTION
!             KODE=1 RETURNS Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN
!                    OR      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN
!                    ON FLGIK = 1.0 OR FLGIK = -1.0
!             KODE=2 RETURNS Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN
!                    OR      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN
!                    ON FLGIK = 1.0 OR FLGIK = -1.0
!     FLGIK - SELECTION PARAMETER FOR I OR K FUNCTION
!             FLGIK =  1.0 GIVES THE I FUNCTION
!             FLGIK = -1.0 GIVES THE K FUNCTION
!        RA - SQRT(1.+Z*Z), Z=X/FNU
!       ARG - ARGUMENT OF THE LEADING EXPONENTIAL
!        IN - NUMBER OF FUNCTIONS DESIRED, IN=1 OR 2
!       TOL - TOLERANCE SPECIFIED BY BESI OR BESK

!                                    OUTPUT

!         Y - A VECTOR WHOSE FIRST IN COMPONENTS CONTAIN THE SEQUENCE

!                                 WRITTEN BY
!                                 D. E. AMOS

!     ABSTRACT
!         ASIK IMPLEMENTS THE UNIFORM ASYMPTOTIC EXPANSION OF
!         THE I AND K BESSEL FUNCTIONS FOR FNU>=35 AND REAL
!         X > 0.0.  THE FORMS ARE IDENTICAL EXCEPT FOR A CHANGE
!         IN SIGN OF SOME OF THE TERMS. THIS CHANGE IN SIGN IS
!         ACCOMPLISHED BY MEANS OF THE FLAG FLGIK = 1 OR -1.

REAL, INTENT(IN)     :: x, fnu, flgik, tol
REAL, INTENT(IN OUT) :: arg, ra
INTEGER, INTENT(IN)  :: kode, in
REAL, INTENT(OUT)    :: y(:)

! Local variables
INTEGER         :: j, jn, k, kk, l
REAL, PARAMETER :: c(65) = (/ -2.08333333333333E-01, 1.25000000000000E-01  &
       , 3.34201388888889E-01, -4.01041666666667E-01, 7.03125000000000E-02  &
       , -1.02581259645062E+00, 1.84646267361111E+00, -8.91210937500000E-01  &
       , 7.32421875000000E-02, 4.66958442342625E+00, -1.12070026162230E+01  &
       , 8.78912353515625E+00, -2.36408691406250E+00, 1.12152099609375E-01  &
       , -2.82120725582002E+01, 8.46362176746007E+01, -9.18182415432400E+01  &
       , 4.25349987453885E+01, -7.36879435947963E+00, 2.27108001708984E-01  &
       , 2.12570130039217E+02, -7.65252468141182E+02, 1.05999045252800E+03  &
       , -6.99579627376133E+02, 2.18190511744212E+02, -2.64914304869516E+01  &
       , 5.72501420974731E-01, -1.91945766231841E+03, 8.06172218173731E+03  &
       , -1.35865500064341E+04, 1.16553933368645E+04, -5.30564697861340E+03  &
       , 1.20090291321635E+03, -1.08090919788395E+02, 1.72772750258446E+00  &
       , 2.02042913309661E+04, -9.69805983886375E+04, 1.92547001232532E+05  &
       , -2.03400177280416E+05, 1.22200464983017E+05, -4.11926549688976E+04  &
       , 7.10951430248936E+03, -4.93915304773088E+02, 6.07404200127348E+00  &
       , -2.42919187900551E+05, 1.31176361466298E+06, -2.99801591853811E+06  &
       , 3.76327129765640E+06, -2.81356322658653E+06, 1.26836527332162E+06  &
       , -3.31645172484564E+05, 4.52187689813627E+04, -2.49983048181121E+03  &
       , 2.43805296995561E+01, 3.28446985307204E+06, -1.97068191184322E+07  &
       , 5.09526024926646E+07, -7.41051482115327E+07, 6.63445122747290E+07  &
       , -3.75671766607634E+07, 1.32887671664218E+07, -2.78561812808645E+06  &
       , 3.08186404612662E+05, -1.38860897537170E+04, 1.10017140269247E+02  &
        /), con(2) = (/ 3.98942280401432678E-01, 1.25331413731550025E+00 /)
REAL            :: ak, ap, coef, etx, fn, gln, s1, s2, t, t2, z
!     ---------------------
fn = fnu
z = (3.0-flgik) / 2.0
kk = INT(z)
DO jn = 1, in
  IF (jn /= 1) THEN
    fn = fn - flgik
    z = x / fn
    ra = SQRT(1.0+z*z)
    gln = LOG((1.0+ra)/z)
    etx = REAL(kode-1)
    t = ra * (1.0-etx) + etx / (z+ra)
    arg = fn * (t-gln) * flgik
  END IF
  coef = EXP(arg)
  t = 1.0 / ra
  t2 = t * t
  t = t / fn
  t = SIGN(t,flgik)
  s2 = 1.0
  ap = 1.0
  l = 0
  DO k = 2, 11
    l = l + 1
    s1 = c(l)
    DO j = 2, k
      l = l + 1
      s1 = s1 * t2 + c(l)
    END DO
    ap = ap * t
    ak = ap * s1
    s2 = s2 + ak
    IF (MAX(ABS(ak),ABS(ap)) < tol) GO TO 30
  END DO

  30 t = ABS(t)
  y(jn) = s2 * coef * SQRT(t) * con(kk)
END DO
RETURN
END SUBROUTINE asik



SUBROUTINE cbesk(z, cnu, w)
!-----------------------------------------------------------------------

!        CALCULATION OF THE MODIFIED BESSEL FUNCTION OF THE
!        SECOND KIND FOR COMPLEX ORDER CNU AND COMPLEX ARGUMENT Z.
!        IT IS ASSUMED THAT -PI < ARG Z <= PI.

!-----------------------------------------------------------------------
!     WRITTEN BY
!        ALFRED H. MORRIS, JR. AND ANDREW H. VAN TUYL
!        NAVAL SURFACE WARFARE CENTER
!        OCT 1992
!--------------------------
COMPLEX, INTENT(IN)  :: z, cnu
COMPLEX, INTENT(OUT) :: w

! Local variables
COMPLEX :: nu
REAL    :: a, b, c, e, rn, rz, s, tau, x, y

x = ABS(REAL(z))
y = ABS(AIMAG(z))
nu = cnu
IF (REAL(nu) < 0.0) nu = -nu
a = REAL(nu)
b = ABS(AIMAG(nu))
rn = cpabs(a,b)
rz = cpabs(x,y)

!     ASYMPTOTIC EXPANSION

IF (rz >= 17.5+0.5*rn*rn) THEN
  CALL cbka(z, nu, w)
  RETURN
END IF

IF (b > 1.5E-2*a) THEN
  tau = rn / rz
  IF (tau < 1.5) THEN
    IF (tau <= 0.05) GO TO 60
    IF (tau >= 0.70.OR.REAL(z) > 0.0) THEN
      IF (AIMAG(nu) > 0.0 .AND. b < 0.07*a) GO TO 60

      IF (tau >= 1.291) THEN
        IF (AIMAG(nu) < 0.0 .AND. b < 0.127*a) GO TO 10
        s = 0.5 * (1.5-tau)
        IF (b < s*a) GO TO 60
        GO TO 50
      END IF

      IF (tau < 0.639) GO TO 30
      IF (tau <= 0.691) THEN
        IF (b > 0.5*a) GO TO 40
        e = b / a
        GO TO 20
      END IF

      IF (b >= 0.191*a) THEN
        IF (tau > 0.99 .AND. b > 0.257*a .AND. b < 0.64*a) GO TO 40
        IF (tau <= 1.16 .AND. b < 0.727*a .AND. y < 0.727*x) GO TO 60
        IF (tau <= 0.91 .AND. a > 0.45*b .AND. y < 0.325*x) GO TO 60
        c = 0.471
        IF (tau < 0.75) c = 0.55
        IF (tau < 0.844 .AND. a < 0.55*b .AND. y < c*x) GO TO 60
      END IF

      10 s = 1.65 * (1.54-tau) ** 2
      IF (tau < 0.91) s = 0.82 - 1.5 * (tau-0.8)
      e = 2.25
      IF (tau < 0.78) e = 2.90
      IF (b > e*s*a) GO TO 50
      IF (b > 0.5*a) GO TO 40

      e = b / (s*a)
      IF (e >= 0.50) THEN
        c = 2.83 - 1.66 * e
        IF (y > c*x) GO TO 50
        GO TO 60
      END IF
      20 c = 7.0 - 10.0 * e
      IF (tau > 0.86) c = 8.0 - 12.0 * e
      IF (y > c*x) GO TO 50
      GO TO 60

      30 IF (b <= 0.191*a) GO TO 60
      40 IF (x >= 0.64*(tau-0.2)*y) THEN
        s = 1.5 * b / (a+1.e-7)
        e = 0.95
        IF (tau > 0.95 .AND. tau < 1.16 .AND. b > 0.471*a .AND. b <=  &
        a) e = 0.75
        IF (tau > 0.85 .AND. tau <= 0.95 .AND. b <= a) e = 0.80
        IF (tau > 0.71 .AND. tau <= 0.85 .AND. b < 1.21*a) e = 0.85
        IF (tau > 0.61 .AND. tau <= 0.71 .AND. b > 0.63*a .AND. b <  &
        1.15*a) e = 0.80
        IF (tau > 0.50 .AND. tau <= 0.61 .AND. b > 0.7*a .AND. b <= a)  &
        e = 0.70

        IF (tau > 0.68 .AND. tau <= 0.77 .AND. a < 0.75*b) e = 1.15
        IF (tau > 0.77 .AND. tau < 0.95 .AND. a < 0.83*b) e = 1.10
        c = (1.0+e*tau) * tau * TANH(s) ** 2
        IF (x >= c*y) GO TO 60
      END IF
    END IF
  END IF

!     CALCULATION IN TERMS OF THE MODIFIED
!     BESSEL FUNCTION I

  50 CALL cbki(z,nu,w)
  RETURN
END IF

!     POWER SERIES AND MILLER ALGORITHM

60 CALL cbkm(z, rz, nu, w)
RETURN
END SUBROUTINE cbesk


SUBROUTINE cbki(z, cnu, w)
!------------------------------------------------------------
!     CALCULATION OF THE MODIFIED BESSEL FUNCTION OF THE
!     SECOND KIND WITH COMPLEX ORDER AND ARGUMENT IN TERMS
!     OF THE MODIFIED BESSEL FUNCTION OF THE FIRST KIND.
!------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, cnu
COMPLEX, INTENT(OUT) :: w

! Local variables
COMPLEX         :: w1, w2
REAL, PARAMETER :: pi = 3.14159265358979E+00, pihalf = 1.5707963267949E+00
REAL            :: a, b, u1, u2

CALL cbssli(z, -cnu, w1)
CALL cbssli(z, cnu, w2)
a = REAL(cnu)
b = pi * AIMAG(cnu)
u1 = sin1(a) * COSH(b)
u2 = cos1(a) * SINH(b)
w = pihalf * cdiv(w1-w2, CMPLX(u1,u2))
RETURN
END SUBROUTINE cbki


SUBROUTINE cbka(z, cnu, w)
!-----------------------------------------------------------------------
!        COMPUTATION OF THE BESSEL FUNCTION K FOR COMPLEX ORDER
!        CNU AND COMPLEX ARGUMENT Z BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, cnu
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL            :: eps, m, u, v
REAL, PARAMETER :: c = 1.77245385090552E+00
COMPLEX         :: a, p, t, zr
INTEGER         :: i
!--------------------------
!     C = PI**(1/2)
!--------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0 .

eps = EPSILON(1.0)

!--------------------------
CALL crec(REAL(z),AIMAG(z),u,v)
zr = CMPLX(0.5*u,0.5*v)
a = cnu * cnu - 0.25
m = 1.0
t = a * zr
p = t

DO i = 1, 16
  a = a - 2.0 * m
  m = m + 1.0
  t = t * a * zr / m
  p = p + t
  IF (anorm(t) <= eps*anorm(p)) GO TO 20
END DO

20 p = p + 1.0
t = SQRT(zr)
IF (AIMAG(z) == 0.0) t = CONJG(t)
w = c * t * p * EXP(-z)
RETURN
END SUBROUTINE cbka


SUBROUTINE cbkm(z, rz, cnu, w)
!-------------------------------------------------------------------
!     CALCULATION OF THE MODIFIED BESSEL FUNCTION OF THE SECOND KIND
!     WITH COMPLEX ORDER AND ARGUMENT BY MEANS OF MACLAURIN EXPANSIONS
!     AND THE MILLER ALGORITHM.  IT IS ASSUMED THAT RZ = ABS(Z).
!-------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, cnu
REAL, INTENT(IN)     :: rz
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: pi = 3.14159265358979, cpi = 1.25331413731550
REAL            :: a, ai, b, c, nu, t, x, y, znorm, zr1, zr2
COMPLEX         :: cz, ex, r, r1, u1, u2, u3, w1, w2, zr
INTEGER         :: i, ind, m, n, n1
!---------------------------
!     CPI = SQRT(PI/2)
!---------------------------

!        REDUCTION OF CNU TO R = NU + B*I WHERE
!                -0.5 < NU <= 0.5

r = cnu
IF (REAL(r) < 0.0) r = -r

a = REAL(r)
b = AIMAG(r)
n = a
nu = a - REAL(n)
t = nu - 0.5
IF (t > 0.0) THEN
  nu = t - 0.5
  n = n + 1
END IF
m = n
r1 = CMPLX(nu, b)

ind = 1
cz = z
x = 0.5 * REAL(z)
y = 0.5 * AIMAG(z)
CALL crec(x, y, zr1, zr2)
zr = CMPLX(zr1, zr2)
IF (t == 0.0) THEN
  IF (b == 0.0) THEN

!               CALCULATION FOR NU = 0.5

    w = CMPLX(cpi,0.0)
    w1 = (1.0,0.0)
    IF (n == 0) GO TO 30
    IF (x < 0.0) THEN
      ind = -1
      cz = -z
      zr = -zr
    END IF
    u1 = w
    u2 = w
    n = n + 1
    r1 = (-0.5,0.0)
    GO TO 10
  END IF
END IF

!           CALCULATION FOR ABS(NU) < 0.5

znorm = 0.5 * rz
IF (znorm <= 1.0) THEN
  ind = 0
  CALL ckps(z, znorm, zr, r1, u1, u2)
ELSE
  u1 = (1.0,0.0)
  IF (x >= 0.0) THEN
    CALL cbkml(z,znorm,zr,r1,n,w1,u2)
  ELSE
    ind = -1
    cz = -z
    zr = -zr
    CALL cbkml(cz,znorm,zr,r1,n,w1,u2)
  END IF
END IF

IF (n <= 1) THEN
  w = u1
  IF (n /= 0) w = u2
  GO TO 30
END IF

!                 RECURSION

10 n1 = n - 1
DO i = 1, n1
  ai = i
  u3 = (r1+ai) * zr * u2 + u1
  u1 = u2
  u2 = u3
END DO
w = u3

30 IF (ind == 0) RETURN
w = w * w1 * EXP(-cz) / SQRT(cz)
IF (ind > 0) RETURN

!            ANALYTIC CONTINUATION

c = EXP(0.5*b*pi)
ex = cxp(m, nu)
IF (y >= 0.0) THEN
  ex = c * ex
  w1 = CMPLX(AIMAG(z), -REAL(z))
  CALL cbsslj(w1, r, w2)
  w2 = CMPLX(pi*AIMAG(w2), -pi*REAL(w2))
  w = ex * (ex*w+w2)
  RETURN
END IF
ex = CONJG(ex) / c
w1 = CMPLX(-AIMAG(z), REAL(z))
CALL cbsslj(w1, r, w2)
w2 = CMPLX(-pi*AIMAG(w2),pi*REAL(w2))
w = ex * (ex*w+w2)
RETURN
END SUBROUTINE cbkm


SUBROUTINE ckps(z, r, zr, nu, w1, w2)
!-----------------------------------------------------------------------

!        CALCULATION OF THE MODIFIED BESSEL FUNCTIONS

!                  W1 = K  (Z)  AND  W2 = K    (Z)
!                        NU                NU+1

!     FOR A COMPLEX ARGUMENT Z WHERE ABS(Z) <= 2 AND A COMPLEX
!     ORDER NU WHERE ABS(REAL(NU)) <= 0.5.  IT IS ASSUMED THAT
!     -PI < ARG Z <= PI, R = ABS(Z/2), AND ZR = 2/Z.  POWER
!     SERIES EXPANSIONS ARE USED.

!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, zr, nu
REAL, INTENT(IN)     :: r
COMPLEX, INTENT(OUT) :: w1, w2

! Local variables
REAL, PARAMETER :: d(7) = (/ .577215664901533E+00, -.420026350340952E-01,  &
        -.421977345555443E-01, .721894324666310E-02, -.215241674114951E-03,  &
        -.201348547807882E-04, .113302723198170E-05 /), tol = 1.e-10,  &
        pi = 3.14159265358979
COMPLEX :: a, c, ch, cl, f, g1, g2, gm1, gm2, mu, p, q, sh, t, t1, t2
INTEGER :: k
REAL    :: ak, eps, eps0, phi, s, u, x, y
!-----------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

!-----------------------------------------------------------------------
eps0 = MAX(eps,5.e-15)

!                     CL = - LN(Z/2)

phi = ATAN2(AIMAG(z),REAL(z))
cl = CMPLX(-LOG(r),-phi)
mu = nu * cl
c = EXP(mu)

!                   G1 = GAMMA(1 + NU)
!                   G2 = GAMMA(1 - NU)

t = pi * nu
IF (anorm(nu) <= tol) THEN
  a = 1.0 + (t*t) / 6.0
ELSE
  a = t / SIN(t)
END IF

s = REAL(nu) ** 2 + AIMAG(nu) ** 2
IF (s < 1.0) THEN
  t = cgam0(nu)
  g1 = 1.0 / t
  g2 = a * t
ELSE
  t = 0.5 + (0.5+nu)
  CALL cgamma(0,t,g1)
  g2 = cdiv(a,g1)
END IF

gm2 = 0.5 * (g1+g2)
IF (s <= 0.04) THEN

!         THE FOLLOWING IS THE TAYLOR SERIES FOR
!     W1 = (1/G2 - 1/G1)/(2*NU). NOTE THAT G1*G2 = A.

  t = nu * nu
  w1 = -((((((d(7)*t+d(6))*t+d(5))*t+d(4))*t+d(3))*t+d(2))*t+d(1))
  gm1 = a * w1
ELSE
  gm1 = 0.5 * (g1-g2) / nu
END IF

!            INITIALIZATION OF THE SUMMATION

p = 0.5 * c * g1
q = 0.5 * cdiv(g2,c)
x = REAL(mu)
y = AIMAG(mu)
IF (anorm(mu) <= tol) THEN
  u = x * y
  sh = CMPLX(1.0,u/3.0)
  ch = CMPLX(1.0,u)
ELSE
  t = CMPLX(-y,x)
  sh = SIN(t) / t
  ch = COS(t)
END IF

f = gm1 * ch + gm2 * cl * sh
c = (1.0,0.0)
w1 = f
w2 = p

!                 SUMMATION OF SERIES

t = 0.25 * (z*z)
DO k = 1, 50
  ak = k
  f = (ak*f+p+q) / ((ak-nu)*(ak+nu))
  p = p / (ak-nu)
  q = q / (ak+nu)
  c = c * t / ak
  t1 = c * f
  w1 = w1 + t1
  t2 = c * (p-ak*f)
  w2 = w2 + t2
  IF (anorm(t1) <= eps0*anorm(w1)) GO TO 20
END DO

20 w2 = w2 * zr
RETURN
END SUBROUTINE ckps


SUBROUTINE cbkml(z, r, zr, nu, n, w, w0)
!-----------------------------------------------------------------------

!             COMPUTATION OF THE SCALED BESSEL FUNCTION

!                   W = (EXP(Z)*SQRT(Z)) K  (Z)
!                                         NU
!     AND THE VALUE

!                   W0 = K    (Z) /K  (Z)
!                         NU+1      NU

!     FOR COMPLEX ORDERS NU AND NU + 1 AND FOR COMPLEX ARGUMENT Z
!     BY USE OF THE MILLER ALGORITHM.  FOR THE GREATEST ACCURACY,
!     Z SHOULD LIE IN A SECTOR SLIGHTLY LARGER THAN THE RIGHT HALF
!     PLANE. IT IS ASSUMED THAT ABS(REAL(NU)) < 0.5, AND THAT
!     R = ABS(Z/2) AND ZR = 2/Z.

!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, zr, nu
REAL, INTENT(IN)     :: r
COMPLEX, INTENT(OUT) :: w, w0
INTEGER, INTENT(IN)  :: n

! Local variables
COMPLEX         :: al, bl, s, u1, u2, u3
REAL            :: a, b, c, e, eps, eps0, f, inu, l, rnu, t, th, x, y
INTEGER         :: i, m, num
REAL, PARAMETER :: c1 = 1.25331413731559E+00, c2 = 1.77245385090552E+00
!---------------------
!     C1 = SQRT(PI/2)
!     C2 = SQRT(PI)
!---------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

!-----------------------------------------------------------------------
eps0 = MAX(eps,5.e-15)

x = REAL(z)
y = AIMAG(z)
rnu = REAL(nu)
inu = AIMAG(nu)

!     CALCULATION OF M FOR USE IN MILLER ALGORITHM.

th = ATAN2(y,x)
a = 3.0 / (1.0+r)
b = 14.7 / (28.0+r)
f = (2.0*r) ** (0.25)
IF (rnu /= 0.5 .OR. inu == 0.0) THEN
  c = 4.0 * cos1(rnu) / (c1*eps0*f)
ELSE
  c = 4.0 * c2 / f
END IF
m = (0.485/r) * (LOG(c)+r*COS(a*th)/(1.0+0.008*r)) ** 2 / (2.0*  &
    COS(b*th)) ** 2 + 1.5
t = 0.0
IF (r > 22.5) t = 0.01 * (r-22.5) ** 2
num = 35.0 + 0.8 * (r-22.5) + t
IF (inu /= 0.0) m = m + num
c = n
IF (c+rnu < 0.327*ABS(inu)) m = m + 10

!     BACKWARD RECURRENCE IN MILLER ALGORITHM.

a = (0.5-rnu) * (0.5+rnu) + inu * inu
e = 2.0 * rnu * inu
u2 = (0.0,0.0)
u1 = (1.0,0.0)
s = u1
l = m
DO i = 2, m
  u3 = u2
  u2 = u1
  c = a / l + (l-1.0)
  al = CMPLX(c,-e/l)
  bl = 2.0 * (l+z)
  u1 = (bl*u2-(l+1.0)*u3) / al
  s = s + u1
  c = ABS(REAL(u1)) + ABS(AIMAG(u1))
  IF (i /= m .AND. c >= 1.e+8) THEN

!           RESCALE TO AVOID OVERFLOW

    u2 = u2 / c
    u1 = u1 / c
    s = s / c
  END IF

  l = l - 1.0
END DO

!     LAST STEP IN THE MILLER ALGORITHM.

IF (c >= 2.0) THEN
  u2 = u2 / c
  u1 = u1 / c
  s = s / c
END IF
u3 = u2
u2 = u1
al = CMPLX(0.5*a,-0.5*e)
bl = 1.0 + z
u1 = (bl*u2-u3) / al
s = s + u1

!     FINAL ASSEMBLY

w = c1 * (u1/s)
w0 = 1.0 + 0.5 * (nu+0.5-u2/u1) * zr
RETURN
END SUBROUTINE cbkml



SUBROUTINE cka(z, nu, w)
!-----------------------------------------------------------------------
!        COMPUTATION OF THE BESSEL FUNCTION K FOR REAL ORDER NU
!        AND COMPLEX ARGUMENT Z BY THE ASYMPTOTIC EXPANSION
!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z
REAL, INTENT(IN)     :: nu
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL            :: a, eps, m, u, v
REAL, PARAMETER :: c = 1.77245385090552
COMPLEX         :: p, t, zr
INTEGER         :: i
!--------------------------
!     C = PI**(1/2)
!--------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0 .

eps = EPSILON(1.0)

!--------------------------
CALL crec(REAL(z), AIMAG(z), u, v)
zr = CMPLX(0.5*u,0.5*v)
a = nu * nu - 0.25
m = 1.0
t = a * zr
p = t

DO i = 1, 16
  a = a - 2.0 * m
  m = m + 1.0
  t = (a/m) * t * zr
  p = p + t
  IF (anorm(t) <= eps*anorm(p)) GO TO 20
END DO

20 p = p + 1.0
t = SQRT(zr)
IF (AIMAG(z) == 0.0) t = CONJG(t)
w = c * t * p * EXP(-z)
RETURN
END SUBROUTINE cka


SUBROUTINE ckm(z, r, zr, nu, w1, w2)
!-----------------------------------------------------------------------

!        CALCULATION OF THE MODIFIED BESSEL FUNCTIONS

!                  W1 = K  (Z)  AND  W2 = K    (Z)
!                        NU                NU+1

!     FOR A COMPLEX ARGUMENT Z WHERE ABS(Z) <= 2  AND A REAL
!     ORDER NU WHERE ABS(NU) <= 0.5.  IT IS ASSUMED THAT
!     -PI < ARG Z <= PI, R = ABS(Z/2), AND ZR = 2/Z. POWER
!     SERIES EXPANSIONS ARE USED.

!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, zr
REAL, INTENT(IN)     :: r, nu
COMPLEX, INTENT(OUT) :: w1, w2

! Local variables
REAL, PARAMETER  :: d(7) = (/ .577215664901533E+00, -.420026350340952E-01,  &
        -.421977345555443E-01, .721894324666310E-02, -.215241674114951E-03,  &
        -.201348547807882E-04, .113302723198170E-05 /), tol = 1.e-10,  &
        pi = 3.14159265358979
REAL    :: a, ak, eps, eps0, gm1, gm2, g1, g2, phi, t, x, y
INTEGER :: k
COMPLEX :: c, ch, cl, f, mu, p, q, sh, t1, t2, w
!-----------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

!-----------------------------------------------------------------------
eps0 = MAX(eps,5.e-15)

!                     CL = - LN(Z/2)

phi = ATAN2(AIMAG(z),REAL(z))
cl = CMPLX(-LOG(r),-phi)
mu = nu * cl
c = EXP(mu)

!                   G1 = GAMMA(1 + NU)
!                   G2 = GAMMA(1 - NU)

t = pi * nu
IF (ABS(nu) <= tol) THEN
  a = 1.0 + (t*t) / 6.0
ELSE
  a = t / SIN(t)
END IF

t = 0.5 + (0.5+gam1(nu))
g1 = 1.0 / t
g2 = a * t
gm2 = 0.5 * (g1+g2)
IF (ABS(nu) <= 0.2) THEN

!         THE FOLLOWING IS THE TAYLOR SERIES FOR
!     T = (1/G2 - 1/G1)/(2*NU). NOTE THAT G1*G2 = A.

  t = nu * nu
  t = -((((((d(7)*t+d(6))*t+d(5))*t+d(4))*t+d(3))*t+d(2))*t+d(1))
  gm1 = a * t
ELSE
  gm1 = 0.5 * (g1-g2) / nu
END IF

!            INITIALIZATION OF THE SUMMATION

p = (0.5*g1) * c
q = (0.5*g2) / c
x = REAL(mu)
y = AIMAG(mu)
IF (anorm(mu) <= tol) THEN
  t = x * y
  sh = CMPLX(1.0,t/3.0)
  ch = CMPLX(1.0,t)
ELSE
  w = CMPLX(-y,x)
  sh = SIN(w) / w
  ch = COS(w)
END IF

f = gm1 * ch + gm2 * cl * sh
c = (1.0,0.0)
w1 = f
w2 = p

!                 SUMMATION OF SERIES

w = 0.25 * (z*z)
DO k = 1, 50
  ak = k
  f = (ak*f+p+q) / ((ak-nu)*(ak+nu))
  p = p / (ak-nu)
  q = q / (ak+nu)
  c = c * w / ak
  t1 = c * f
  w1 = w1 + t1
  t2 = c * (p-ak*f)
  w2 = w2 + t2
  IF (anorm(t1) <= eps0*anorm(w1)) GO TO 20
END DO

20 w2 = w2 * zr
RETURN
END SUBROUTINE ckm


SUBROUTINE ckml(z, r, zr, nu, k1, k2)
!-----------------------------------------------------------------------

!             COMPUTATION OF THE SCALED BESSEL FUNCTIONS

!                   K1 = (EXP(Z)*SQRT(Z)) K  (Z)
!                                          NU
!                   K2 = (EXP(Z)*SQRT(Z)) K    (Z)
!                                          NU+1

!     FOR REAL ORDERS NU AND NU + 1 AND FOR COMPLEX ARGUMENT Z
!     BY USE OF THE MILLER ALGORITHM. FOR THE GREATEST ACCURACY,
!     Z SHOULD LIE IN A SECTOR SLIGHTLY LARGER THAN THE RIGHT
!     HALF PLANE. IT IS ASSUMED THAT ABS(NU) < 0.5, AND THAT
!     R = ABS(Z/2) AND ZR = 2/Z.

!-----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z, zr
COMPLEX, INTENT(OUT) :: k1, k2
REAL, INTENT(IN)     :: r, nu

! Local variables
COMPLEX         :: bl, s, u1, u2, u3
REAL            :: a, al, b, c, e, eps, eps0, l, nu2, th, x, y
INTEGER         :: i, m
REAL, PARAMETER :: c1 = 1.25331413731559
!---------------------
!     C1 = SQRT(PI/2)
!---------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

!-----------------------------------------------------------------------
eps0 = MAX(eps,5.e-15)
x = REAL(z)
y = AIMAG(z)
nu2 = nu * nu

!     CALCULATION OF M FOR USE IN MILLER ALGORITHM.

th = ATAN2(y,x)
a = 3.0 / (1.0+r)
b = 14.7 / (28.0+r)
c = 4.0 * cos1(nu) / (c1*eps0*(2.0*r)**(0.25))
m = (0.485/r) * (LOG(c)+r*COS(a*th)/(1.0+0.008*r)) ** 2 / (2.0*  &
COS(b*th)) ** 2 + 1.5

!     BACKWARD RECURRENCE IN MILLER ALGORITHM.

u2 = (0.0,0.0)
u1 = (1.0,0.0)
s = u1
l = m
DO i = 2, m
  u3 = u2
  u2 = u1
  e = l - 0.5
  al = (e*e-nu2) / (l*(l+1.0))
  bl = (2.0/(l+1.0)) * (l+z)
  u1 = (bl*u2-u3) / al
  s = s + u1
  c = ABS(REAL(u1)) + ABS(AIMAG(u1))
  IF (i /= m .AND. c >= 1.e+6) THEN

!           RESCALE TO AVOID OVERFLOW

    u2 = u2 / c
    u1 = u1 / c
    s = s / c
  END IF

  l = l - 1.0
END DO

!     LAST STEP IN THE MILLER ALGORITHM.

IF (c >= 2.0) THEN
  u2 = u2 / c
  u1 = u1 / c
  s = s / c
END IF
u3 = u2
u2 = u1
al = 0.5 * (0.5-nu) * (0.5+nu)
bl = 1.0 + z
u1 = (bl*u2-u3) / al
s = s + u1

!     FINAL ASSEMBLY

k1 = c1 * (u1/s)
k2 = k1 * (1.0+0.5*(nu+0.5-u2/u1)*zr)
RETURN
END SUBROUTINE ckml


FUNCTION cxp(n, nu) RESULT(fn_val)
!-----------------------------------------------------------------------
!                 COMPUTATION OF EXP(-R*(PI/2)*I)
!              WHERE R = N + NU FOR ABS(NU) <= 0.5
!-----------------------------------------------------------------------
REAL, INTENT(IN)    :: nu
INTEGER, INTENT(IN) :: n
COMPLEX             :: fn_val

! Local variables
REAL    :: c, s
INTEGER :: k

c = cos0(nu)
s = sin0(nu)
k = MOD(n,4)
IF (k /= 0) THEN
  IF (k == 1) GO TO 10
  IF (k == 2) GO TO 20
  GO TO 30
END IF

fn_val = CMPLX(c,-s)
RETURN

10 fn_val = CMPLX(-s,-c)
RETURN

20 fn_val = CMPLX(-c,s)
RETURN

30 fn_val = CMPLX(s,c)
RETURN
END FUNCTION cxp


SUBROUTINE bsslk(mo, a, in, w)
!     ******************************************************************
!     FORTRAN SUBROUTINE FOR MODIFIED BESSEL FUNCTION OF INTEGRAL ORDER
!     ******************************************************************
!     MO = MODE OF OPERATION
!     A  = ARGUMENT (COMPLEX NUMBER)
!     IN = ORDER (INTEGER)
!     W  = FUNCTION OF SECOND KIND (COMPLEX NUMBER)
!     -------------------
INTEGER, INTENT(IN)  :: mo, in
COMPLEX, INTENT(IN)  :: a
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL, PARAMETER :: cd(30) = (/ 0.00000000000000, -1.64899505142212E-2,  &
        -7.18621880068536E-2, -1.67086878124866E-1, -3.02582250219469E-1,  &
        -4.80613945245927E-1, -7.07075239357898E-1, -9.92995790539516E-1,  &
        -1.35583925612592, -1.82105907899132, -2.42482175310879,  &
        -3.21956655708750, -4.28658077248384, -5.77022816798128,  &
        -8.01371260952526, 0.00000000000000, -5.57742429879505E-3,  &
        -4.99112944172476E-2, -1.37440911652397E-1, -2.67233784710566E-1,  &
        -4.40380166808682E-1, -6.61813614872541E-1, -9.41861077665017E-1,  &
        -1.29754130468326, -1.75407696719816, -2.34755299882276,  &
        -3.13041332689196, -4.18397120563729, -5.65251799214994,  &
        -7.87863959810677 /),  &
        ce(30) = (/ 0.00000000000000,  &
        -4.80942336387447E-3, -1.31366200347759E-2, -1.94843834008458E-2,  &
        -2.19948900032003E-2, -2.09396625676519E-2, -1.74600268458650E-2,  &
        -1.27937813362085E-2, -8.05234421796592E-3, -4.15817375002760E-3,  &
        -1.64317738747922E-3, -4.49175585314709E-4, -7.28594765574007E-5,  &
        -5.38265230658285E-6, -9.93779048036289E-8, 0.00000000000000,  &
        7.53805779200591E-2, 7.12293537403464E-2, 6.33116224228200E-2,  &
        5.28240264523301E-2, 4.13305359441492E-2, 3.01350573947510E-2,  &
        2.01043439592720E-2, 1.18552223068074E-2, 5.86055510956010E-3,  &
        2.25465148267325E-3, 6.08173041536336E-4, 9.84215550625747E-5,  &
        7.32139093038089E-6, 1.37279667384666E-7 /)
REAL    :: az(2), sz(2), rz(2), zl(2), ts(2), tm(2), sm(2), sl(2), sq(2), &
           sr(2), aq(2), qf(2)
REAL    :: an, pm, pn, qm, qn, rm, rn, sn, ss, zs
INTEGER :: i, la, lr, lk, m, n
!     -------------------
az(1) = REAL(a)
az(2) = AIMAG(a)
zs = az(1) * az(1) + az(2) * az(2)
zl(1) = 0.5 * LOG(zs)
zl(2) = ATAN2(az(2),az(1))
an = ABS(in)
tm(1) = 0.0
tm(2) = 0.0
IF (mo == 0) THEN
  tm(1) = az(1)
  tm(2) = az(2)
END IF
IF (zs <= 1.0) GO TO 100
IF (zs < 289.0) THEN
  IF (az(1) + 0.096*az(2)*az(2) > 0.0) THEN
    GO TO 50
  ELSE
    GO TO 100
  END IF
END IF
qm = 1.25331413731550 * EXP(-0.5*zl(1)-tm(1))
qf(1) = qm * COS(-0.5*zl(2)-tm(2))
qf(2) = qm * SIN(-0.5*zl(2)-tm(2))
IF (an > 1.0) GO TO 20
pn = an
la = 10
GO TO 190
10 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)
GO TO 180
20 pn = 1.0
la = 30
GO TO 190
30 sq(1) = qf(1) * sm(1) - qf(2) * sm(2)
sq(2) = qf(1) * sm(2) + qf(2) * sm(1)
pn = 0.0
la = 40
GO TO 190
40 sr(1) = qf(1) * sm(1) - qf(2) * sm(2)
sr(2) = qf(1) * sm(2) + qf(2) * sm(1)
GO TO 150
50 qm = 1.25331413731550 * EXP(-0.5*zl(1)-tm(1))
qf(1) = qm * COS(-0.5*zl(2)-tm(2))
qf(2) = qm * SIN(-0.5*zl(2)-tm(2))
IF (an > 1.0) GO TO 70
pn = an
lr = 60
GO TO 230
60 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)
GO TO 180
70 pn = 1.0
lr = 80
GO TO 230
80 sq(1) = qf(1) * sm(1) - qf(2) * sm(2)
sq(2) = qf(1) * sm(2) + qf(2) * sm(1)
pn = 0.0
lr = 90
GO TO 230
90 sr(1) = qf(1) * sm(1) - qf(2) * sm(2)
sr(2) = qf(1) * sm(2) + qf(2) * sm(1)
GO TO 150
100 qf(1) = 1.0
qf(2) = 0.0
IF (mo /= 0) THEN
  qm = EXP(az(1))
  qf(1) = qm * COS(az(2))
  qf(2) = qm * SIN(az(2))
END IF
IF (an > 1.0) GO TO 120
pn = an
lk = 110
GO TO 250
110 ts(1) = qf(1) * sm(1) - qf(2) * sm(2)
ts(2) = qf(1) * sm(2) + qf(2) * sm(1)
sm(1) = ts(1)
sm(2) = ts(2)
GO TO 180
120 pn = 1.0
lk = 130
GO TO 250
130 sq(1) = qf(1) * sm(1) - qf(2) * sm(2)
sq(2) = qf(1) * sm(2) + qf(2) * sm(1)
pn = 0.0
lk = 140
GO TO 250
140 sr(1) = qf(1) * sm(1) - qf(2) * sm(2)
sr(2) = qf(1) * sm(2) + qf(2) * sm(1)
150 rz(1) = +az(1) / zs
rz(2) = -az(2) / zs
pn = 0.0
GO TO 170
160 sq(1) = sr(1)
sq(2) = sr(2)
sr(1) = sm(1)
sr(2) = sm(2)
170 sm(1) = 2.0 * pn * (rz(1)*sr(1)-rz(2)*sr(2)) + sq(1)
sm(2) = 2.0 * pn * (rz(1)*sr(2)+rz(2)*sr(1)) + sq(2)
pn = pn + 1.0
IF (pn < an) GO TO 160
180 w = CMPLX(sm(1),sm(2))
RETURN

190 sm(1) = 0.0
sm(2) = 0.0
rz(1) = +0.5 * az(1) / zs
rz(2) = -0.5 * az(2) / zs
qn = (pn-0.5) * (pn+0.5)
tm(1) = 1.0
tm(2) = 0.0
pm = 0.0
GO TO 210
200 qn = qn - 2.0 * pm
pm = pm + 1.0
ts(1) = rz(1) * tm(1) - rz(2) * tm(2)
ts(2) = rz(1) * tm(2) + rz(2) * tm(1)
tm(1) = qn * ts(1) / pm
tm(2) = qn * ts(2) / pm
IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
  IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 220
END IF
210 sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (pm < 36.0) GO TO 200
220 SELECT CASE (la)
  CASE (10)
    GO TO 10
  CASE (30)
    GO TO 30
  CASE (40)
    GO TO 40
END SELECT

230 sm(1) = 1.0
sm(2) = 0.0
m = 15.0 * pn + 2.0
n = 15.0 * pn + 15.0
DO i = m, n
  ts(1) = az(1) - cd(i)
  ts(2) = az(2)
  ss = ts(1) * ts(1) + ts(2) * ts(2)
  tm(1) = +ce(i) * ts(1) / ss
  tm(2) = -ce(i) * ts(2) / ss
  sm(1) = sm(1) + tm(1)
  sm(2) = sm(2) + tm(2)
END DO
SELECT CASE (lr)
  CASE (60)
    GO TO 60
  CASE (80)
    GO TO 80
  CASE (90)
    GO TO 90
END SELECT

250 aq(1) = 1.0
aq(2) = 0.0
rn = 0.0
sn = -1.0
pm = 0.0
GO TO 270
260 pm = pm + 1.0
rn = rn + 0.5 / pm
sn = -sn
ts(1) = 0.5 * (az(1)*aq(1)-az(2)*aq(2))
ts(2) = 0.5 * (az(1)*aq(2)+az(2)*aq(1))
aq(1) = ts(1) / pm
aq(2) = ts(2) / pm
270 IF (pm < pn) GO TO 260
sz(1) = 0.25 * (az(1)-az(2)) * (az(1)+az(2))
sz(2) = 0.5 * az(1) * az(2)
sr(1) = 0.0
sr(2) = 0.0
ss = aq(1) * aq(1) + aq(2) * aq(2)
tm(1) = +aq(1) / ss
tm(2) = -aq(2) / ss
pm = 0.0
GO TO 290
280 tm(1) = tm(1) / (pn-pm)
tm(2) = tm(2) / (pn-pm)
sr(1) = sr(1) + 0.5 * tm(1)
sr(2) = sr(2) + 0.5 * tm(2)
pm = pm + 1.0
ts(1) = sz(1) * tm(1) - sz(2) * tm(2)
ts(2) = sz(1) * tm(2) + sz(2) * tm(1)
tm(1) = -ts(1) / pm
tm(2) = -ts(2) / pm
290 IF (pm < pn) GO TO 280
sm(1) = 0.0
sm(2) = 0.0
rm = 1.0
qm = 0.0
aq(1) = sn * aq(1)
aq(2) = sn * aq(2)
sl(1) = -0.115931515658412 + zl(1) - rn
sl(2) = +zl(2)
pm = 0.0
GO TO 310

300 qm = qm + rm
pm = pm + 1.0
rm = 0.25 * zs * rm / (pm*(pn+pm))
ts(1) = sz(1) * aq(1) - sz(2) * aq(2)
ts(2) = sz(1) * aq(2) + sz(2) * aq(1)
aq(1) = ts(1) / (pm*(pn+pm))
aq(2) = ts(2) / (pm*(pn+pm))
sl(1) = sl(1) - 0.5 / pm - 0.5 / (pn+pm)
310 tm(1) = aq(1) * sl(1) - aq(2) * sl(2)
tm(2) = aq(1) * sl(2) + aq(2) * sl(1)
sm(1) = sm(1) + tm(1)
sm(2) = sm(2) + tm(2)
IF (qm+rm > qm) GO TO 300
sm(1) = sr(1) + sm(1)
sm(2) = sr(2) + sm(2)

SELECT CASE (lk)
  CASE (110)
    GO TO 110
  CASE (130)
    GO TO 130
  CASE (140)
    GO TO 140
END SELECT

RETURN
END SUBROUTINE bsslk


SUBROUTINE cai(ind, z, ai, aip, ierr)
!------------------------------------------------------------
!     CALCULATES THE AIRY FUNCTION AI AND ITS DERIVATIVE AIP
!     FOR COMPLEX ARGUMENT Z.
!------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: ai, aip
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL    :: a, b, r

ierr = 0
a = REAL(z)
b = AIMAG(z)
r = cpabs(a,b)
IF (r <= 1.0) THEN

!     MACLAURIN EXPANSION

  CALL airm(ind, z, ai=ai, aip=aip)
  RETURN
END IF
IF (r <= 10.0) THEN

!     INTERMEDIATE RANGE CALCULATION

  CALL aii(ind, z, ai, aip, ierr)
  RETURN
END IF

!     ASYMPTOTIC EXPANSION

CALL aia(ind, z, ai, aip, ierr)
RETURN
END SUBROUTINE cai


SUBROUTINE airm(ind, z, ai, aip, bi, bip)
!--------------------------------------------------------------
!     CALCULATES THE AIRY FUNCTIONS AI AND BI AND THEIR DERIVATIVES
!     AIP AND BIP BY USE OF THEIR MACLAURIN EXPANSIONS.
!--------------------------------------------------------------
INTEGER, INTENT(IN)            :: ind
COMPLEX, INTENT(IN)            :: z
COMPLEX, INTENT(OUT), OPTIONAL :: ai, aip, bi, bip

! Local variables
COMPLEX :: z1, z2, z3, zz, f, f1, g, g1, e, e1
REAL, PARAMETER :: a(8) = (/ .166666666666667E+00, .555555555555556E-02,  &
        .771604938271605E-04, .584549195660307E-06, .278356759838241E-08,  &
        .909662613850462E-11, .216586336631062E-13, .392366551867867E-16 /),  &
        b(8) = (/ .833333333333333E-01, .198412698412698E-02,  &
        .220458553791887E-04, .141319585764030E-06, .588831607350126E-09,  &
        .172172984605300E-11, .372668797846970E-14, .621114663078283E-17 /),  &
        c(8) = (/ .333333333333333E-01, .694444444444444E-03,  &
        .701459034792368E-05, .417535139757362E-07, .163739270493083E-09,  &
        .454831306925231E-12, .941679724482880E-15, .150910212256872E-17 /),  &
        d(8) = (/ .333333333333333E+00, .138888888888889E-01,  &
        .220458553791887E-03, .183715461493239E-05, .942130571760201E-08,  &
        .327128670750070E-10, .819871355263333E-13, .155278665769571E-15 /),  &
        c1 = 3.55028053887817E-01, c2 = 2.58819403792807E-01,  &
        sqt3 = 1.73205080756888E+00
REAL    :: x, y
INTEGER :: i, n
!-----------------------
!     C1 = 3**(-2/3)/GAMMA(2/3)
!     C2 = 3**(-1/3)/GAMMA(1/3)
!-----------------------
z2 = z * z
z3 = z * z2

!     SUMMATION OF F AND G

f = CMPLX(a(8),0.0)
g = CMPLX(b(8),0.0)
DO n = 1, 7
  i = 8 - n
  f = a(i) + z3 * f
  g = b(i) + z3 * g
END DO
f = 1.0 + z3 * f
g = z + z2 * z2 * g

!     SUMMATION OF F1 AND G1

f1 = CMPLX(c(8),0.0)
g1 = CMPLX(d(8),0.0)
DO n = 1, 7
  i = 8 - n
  f1 = c(i) + z3 * f1
  g1 = d(i) + z3 * g1
END DO
f1 = z2 * (0.5+z3*f1)
g1 = 1.0 + z3 * g1

!     FINAL ASSEMBLY

IF (PRESENT(ai)) ai = c1 * f - c2 * g
IF (PRESENT(bi)) bi = sqt3 * (c1*f + c2*g)
IF (PRESENT(aip)) aip = c1 * f1 - c2 * g1
IF (PRESENT(bip)) bip = sqt3 * (c1*f1 + c2*g1)
IF (ind == 0) RETURN

x = REAL(z)
y = AIMAG(z)
z1 = SQRT(z)
zz = z * z1 / 1.5
e = EXP(zz)
IF (PRESENT(ai)) ai = ai * e
IF (PRESENT(aip)) aip = aip * e
IF (PRESENT(bi) .OR. PRESENT(bip)) THEN
  IF (ABS(y) <= x*sqt3) THEN
    e1 = 1.0 / e
    IF (PRESENT(bi)) bi = bi * e1
    IF (PRESENT(bip)) bip = bip * e1
    RETURN
  END IF
END IF
bi = bi * e
bip = bip * e

RETURN
END SUBROUTINE airm


SUBROUTINE aii(ind, z, ai, aip, ierr)
!------------------------------------------------------------
!     CALCULATES THE AIRY FUNCTION AI AND ITS DERIVATIVE AIP FOR
!     COMPLEX ARGUMENT Z IN THE INTERMEDIATE RANGE 1 <= ABS(Z) <= 10.0.
!------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: ai, aip
INTEGER, INTENT(OUT) :: ierr

! Local variables
COMPLEX         :: z1, z2, z3, zm, w1, w2, w1m, w2m, e
REAL, PARAMETER :: c1 = 1.83776298473931E-01
REAL            :: a, b, r

!     C1 = 1/(PI*SQRT(3))

ierr = 0
a = REAL(z)
b = AIMAG(z)
r = cpabs(a,b)
z1 = SQRT(z)
z2 = z1 * z / 1.5
IF (ABS(b) >= -5.0*a) THEN

!           ----  ABS(B) >= -5.0*A  ----

  CALL ka(ind, z2, w1, w2)
  ai = c1 * z1 * w1
  aip = -c1 * z * w2
  RETURN
END IF

!           ----  ABS(B) < -5.0*A  ----

IF (ABS(b) < -1.74*a) GO TO 20
IF (r >= 8.2) GO TO 30
10 zm = -z
z1 = SQRT(zm)
z3 = z1 * zm / 1.5
CALL ja(z3, w1, w2, w1m, w2m)
ai = (z1/3.0) * (w1m+w1)
aip = (z/3.0) * (w2m-w2)
IF (ind == 0) RETURN
e = EXP(z2)
ai = ai * e
aip = aip * e
RETURN
20 IF (r < 7.4) GO TO 10
30 CALL aia(ind, z, ai, aip, ierr)
RETURN
END SUBROUTINE aii


SUBROUTINE aia(ind, z, ai, aip, ierr)
!-----------------------------------------------------------------------
!     CALCULATES THE AIRY FUNCTION AI AND ITS DERIVATIVE AIP FOR
!     COMPLEX ARGUMENT Z BY MEANS OF ASYMPTOTIC EXPANSIONS.
!-----------------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: ai, aip
INTEGER, INTENT(OUT) :: ierr

! Local variables
COMPLEX :: z1, z2, z2r, zz, w, w2, s1, s2, s3, s4, e,zeta, si, cn
COMPLEX :: alpha, beta, j
REAL, PARAMETER :: c(30) = (/ .100000000000000E+01, .694444444444444E-01,  &
        .371334876543210E-01, .379930591278006E-01, .576491904126697E-01,  &
        .116099064025515E+00, .291591399230751E+00, .877666969510017E+00,  &
        .307945303017317E+01, .123415733323452E+02, .556227853659171E+02,  &
        .278465080777603E+03, .153316943201280E+04, .920720659972641E+04,  &
        .598925135658791E+05, .419524875116551E+06, .314825741786683E+07,  &
        .251989198716024E+08, .214288036963680E+09, .192937554918249E+10,  &
        .183357669378906E+11, .183418303528833E+12, .192647115897045E+13,  &
        .211969993886476E+14, .243826826879716E+15, .292659921929793E+16,  &
        .365903070126431E+17, .475768102036307E+18, .642404935790194E+19,  &
        .899520742705838E+20 /), d(30) = (/ .100000000000000E+01,  &
        -.972222222222222E-01, -.438850308641975E-01, -.424628307898948E-01,  &
        -.626621634920323E-01, -.124105896027275E+00, -.308253764901079E+00,  &
        -.920479992412945E+00, -.321049358464862E+01, -.128072930807356E+02,  &
        -.575083035139143E+02, -.287033237109221E+03, -.157635730333710E+04,  &
        -.944635482309593E+04, -.613357066638521E+05, -.428952400400069E+06,  &
        -.321453652140086E+07, -.256979083839113E+08, -.218293420832160E+09,  &
        -.196352378899103E+10, -.186439310881072E+11, -.186352996385294E+12,  &
        -.195588293238984E+13, -.215064446351972E+14, -.247236992290621E+15,  &
        -.296588243029521E+16, -.370624400063547E+17, -.481678264794522E+18,  &
        -.650098408075106E+19, -.909919826436541E+20 /),  &
        c1 = .564189583547756, c2 = .398942280401433
REAL    :: eps, r, t, u, u1, v, v1, xm, xpos, xneg
INTEGER :: i, k, m, m2
!------------------------

!     EPS, XPOS, AND XNEG ARE MACHINE DEPENDENT CONSTANTS. EPS IS
!     THE SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0, XPOS IS THE
!     THE LARGEST POSTIVE NUMBER FOR WHICH EXP(X) CAN BE COMPUTED,
!     AND XNEG IS THE LARGEST NEGATIVE NUMBER FOR WHICH EXP(X) DOES
!     NOT UNDERFLOW.

eps = EPSILON(1.0)
xpos = exparg(0)
xneg = exparg(1)

!------------------------
ierr = 0
IF (REAL(z) >= 0.0) THEN

!             -----   REAL(Z) >= 0    -----

  z1 = SQRT(z)
  z2 = SQRT(z1)
  z2r = 1.0 / z2
  CALL crec(REAL(z),AIMAG(z),u,v)
  w = -1.5 * CMPLX(u,v) / z1
  u = ABS(REAL(w))
  v = ABS(AIMAG(w))
  t = MAX(u,v)
  IF (ind == 0) THEN

    IF (t == 0.0) GO TO 30
    u1 = u / t
    v1 = v / t
    r = u * u1 + v * v1
    xm = xpos
    IF (REAL(w) < 0.0) xm = -xneg
    IF (u1 >= r*xm.OR.v1 >= 0.1*r/eps) GO TO 30
    zeta = z1 * z / 1.5
    e = EXP(-zeta)
  END IF

  m = 20
  IF (t > 30.0) m = 8
  s1 = CMPLX(c(m),0.0)
  s2 = CMPLX(d(m),0.0)
  i = m
  DO k = 2, m
    i = i - 1
    s1 = c(i) + w * s1
    s2 = d(i) + w * s2
  END DO

  ai = 0.5 * c1 * z2r * s1
  aip = -0.5 * c1 * z2 * s2
  IF (ind /= 0) RETURN
  ai = e * ai
  aip = e * aip
  RETURN
END IF

!             -----   REAL(Z) < 0    -----

zz = -z
z1 = SQRT(zz)
z2 = SQRT(z1)
z2r = 1.0 / z2
CALL crec(REAL(zz),AIMAG(zz),u,v)
w = 1.5 * CMPLX(u,v) / z1
u = ABS(REAL(w))
v = ABS(AIMAG(w))
t = MAX(u,v)

IF (t /= 0.0) THEN
  u1 = u / t
  v1 = v / t
  r = u * u1 + v * v1
  IF (ind == 0) THEN
    IF (v1 >= r*xpos.OR.u1 >= 0.1*r/eps) GO TO 30
    zeta = z1 * zz / 1.5
  ELSE
    e = (0.0,0.0)
    j = (0.0,-1.0)
    IF (AIMAG(z) < 0.0) j = (0.0,1.0)
    IF (v1 <= 0.5*r*ABS(xneg)) THEN
      IF (u1 >= 0.05*r/eps) GO TO 30
      zeta = z1 * zz / 1.5
      e = EXP(2.0*j*zeta)
    END IF
  END IF

  w2 = w * w
  m = 15
  IF (t > 30.0) m = 5
  m2 = m + m
  i = m2 - 1
  s1 = CMPLX(c(i),0.0)
  s2 = CMPLX(c(m2),0.0)
  s3 = CMPLX(d(i),0.0)
  s4 = CMPLX(d(m2),0.0)
  DO k = 2, m
    i = i - 1
    s2 = c(i) - s2 * w2
    s4 = d(i) - s4 * w2
    i = i - 1
    s1 = c(i) - s1 * w2
    s3 = d(i) - s3 * w2
  END DO
  s2 = w * s2
  s4 = w * s4

  IF (ind == 0) THEN
    cn = COS(zeta)
    si = SIN(zeta)
  ELSE
    cn = 0.5 * (1.0+e)
    si = 0.5 * (1.0-e) * j
  END IF

  alpha = s1 - s2
  beta = s1 + s2
  ai = c2 * z2r * (alpha*cn+beta*si)
  alpha = s3 - s4
  beta = s3 + s4
  aip = c2 * z2 * (alpha*si-beta*cn)
  RETURN
END IF

!         RETURN WITH ZERO VALUES IF SCALING IS NEEDED

30 ai = (0.0,0.0)
aip = (0.0,0.0)
ierr = 1
RETURN
END SUBROUTINE aia



SUBROUTINE cbi(ind, z, bi, bip, ierr)
!------------------------------------------------------------
!     CALCULATES THE AIRY FUNCTION BI AND ITS DERIVATIVE BIP
!     FOR COMPLEX ARGUMENT Z.
!------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: bi, bip
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL    :: a, b, r

ierr = 0
a = REAL(z)
b = AIMAG(z)
r = cpabs(a,b)
IF (r <= 1.0) THEN

!     MACLAURIN EXPANSION

  CALL airm(ind, z, bi=bi, bip=bip)
  RETURN
END IF
IF (r <= 9.6) THEN

!     INTERMEDIATE RANGE CALCULATION

  CALL bii(ind, z, bi, bip, ierr)
  RETURN
END IF

!     ASYMPTOTIC EXPANSION

CALL bia(ind, z, bi, bip, ierr)
RETURN
END SUBROUTINE cbi


SUBROUTINE bii(ind, z, bi, bip, ierr)
!------------------------------------------------------------
!     CALCULATES THE AIRY FUNCTION BI AND ITS DERIVATIVE BIP FOR
!     COMPLEX ARGUMENT Z IN THE INTERMEDIATE RANGE 1 <= ABS(Z) <= 10.0.
!------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: bi, bip
INTEGER, INTENT(OUT) :: ierr

! Local variables
COMPLEX         :: z1, z2, zm, w1, w2, w1m, w2m, e, e1
REAL, PARAMETER :: c1 = 5.77350269189626E-01
REAL            :: a, r, x, y

!     C1 = 1/SQRT(3)

ierr = 0
x = REAL(z)
y = AIMAG(z)
r = cpabs(x,y)
z1 = SQRT(z)
z2 = z1 * z / 1.5
e = EXP(-z2)
e1 = 1.0 / e
IF (REAL(z) >= 0.0) THEN

!               ----  REAL(Z) >= 0  ----

  IF (r >= 8.9) THEN
    a = 0.156 * r - 0.913
    IF (ABS(y) < a*x.OR.ABS(y) > 0.58*x) GO TO 10
  END IF
  CALL ia(z2,w1,w2,w1m,w2m)
  bi = c1 * z1 * (w1+w1m)
  bip = c1 * z * (w2+w2m)
  IF (ind == 0) RETURN
ELSE

!               ----  REAL(Z) < 0  ----

  IF (r >= 8.1) THEN
    IF (ABS(y) < 3.89*ABS(x)) GO TO 10
  END IF
  zm = -z
  z1 = SQRT(zm)
  z2 = z1 * zm / 1.5
  CALL ja(z2,w1,w2,w1m,w2m)
  bi = c1 * z1 * (w1m-w1)
  bip = c1 * zm * (w2m+w2)
  IF (ind == 0) RETURN
END IF
IF (x < c1*ABS(y)) THEN
  bi = bi * e1
  bip = bip * e1
  RETURN
END IF
bi = bi * e
bip = bip * e
RETURN

10 CALL bia(ind, z, bi, bip, ierr)
RETURN
END SUBROUTINE bii


SUBROUTINE bia(ind, z, bi, bip, ierr)
!---------------------------------------------------------------
!     CALCULATES THE AIRY FUNCTION BI AND ITS DERIVATIVE BIP FOR
!     COMPLEX ARGUMENT Z BY MEANS OF ASYMPTOTIC EXPANSIONS.
!---------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: bi, bip
INTEGER, INTENT(OUT) :: ierr

! Local variables
COMPLEX :: z1, z2, z2r, zz, w, w2, s1, s2, s3, s4, e, zeta, si, cn, cf1, cf2, &
           alpha, beta, j, cz
COMPLEX, PARAMETER :: ex3c = (5.e-01,-8.66025403784439E-01),  &
                      ex6 = (8.66025403784439E-01,5.e-01),    &
                      ex6c = (8.66025403784439E-01, -5.e-01), &
                      cln2 = (0.0,3.46573590279973E-01)
REAL, PARAMETER :: c(30) = (/ .100000000000000E+01, .694444444444444E-01,  &
        .371334876543210E-01, .379930591278006E-01, .576491904126697E-01,  &
        .116099064025515E+00, .291591399230751E+00, .877666969510017E+00,  &
        .307945303017317E+01, .123415733323452E+02, .556227853659171E+02,  &
        .278465080777603E+03, .153316943201280E+04, .920720659972641E+04,  &
        .598925135658791E+05, .419524875116551E+06, .314825741786683E+07,  &
        .251989198716024E+08, .214288036963680E+09, .192937554918249E+10,  &
        .183357669378906E+11, .183418303528833E+12, .192647115897045E+13,  &
        .211969993886476E+14, .243826826879716E+15, .292659921929793E+16,  &
        .365903070126431E+17, .475768102036307E+18, .642404935790194E+19,  &
        .899520742705838E+20 /),  &
        d(30) = (/ .100000000000000E+01,  &
        -.972222222222222E-01, -.438850308641975E-01, -.424628307898948E-01,  &
        -.626621634920323E-01, -.124105896027275E+00, -.308253764901079E+00,  &
        -.920479992412945E+00, -.321049358464862E+01, -.128072930807356E+02,  &
        -.575083035139143E+02, -.287033237109221E+03, -.157635730333710E+04,  &
        -.944635482309593E+04, -.613357066638521E+05, -.428952400400069E+06,  &
        -.321453652140086E+07, -.256979083839113E+08, -.218293420832160E+09,  &
        -.196352378899103E+10, -.186439310881072E+11, -.186352996385294E+12,  &
        -.195588293238984E+13, -.215064446351972E+14, -.247236992290621E+15,  &
        -.296588243029521E+16, -.370624400063547E+17, -.481678264794522E+18,  &
        -.650098408075106E+19, -.909919826436541E+20 /),  &
        sqt3 = 1.73205080756888, c1 = 5.64189583547756E-01,   &
        c2 = 3.98942280401433E-01, c3 = 7.07106781186548E-01
REAL    :: ce, cf, eps, r, s, t, u, u1, v, v1, x, xpos, xneg, y
INTEGER :: i, k, m, m2
!------------------------

!     EPS AND XM ARE MACHINE DEPENDENT CONSTANTS. EPS IS THE
!     SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0, XPOS IS THE
!     LARGEST POSITIVE NUMBER FOR WHICH EXP(XM) CAN BE COMPUTED,
!     AND XNEG IS THE NEGATIVE NUMBER OF LARGEST MAGNITUDE FOR
!     WHICH EXP(X) DOES NOT UNDERFLOW.

eps = EPSILON(1.0)
xpos = exparg(0)
xneg = exparg(1)

!------------------------
ierr = 0
x = REAL(z)
y = AIMAG(z)
IF (x >= ABS(y)*sqt3) THEN

!          -----  ABS(ARG(Z)) <= PI/6  ----

  z1 = SQRT(z)
  z2 = SQRT(z1)
  z2r = 1.0 / z2
  CALL crec(x,y,u,v)
  w = 1.5 * CMPLX(u,v) / z1
  u = ABS(REAL(w))
  v = ABS(AIMAG(w))
  t = MAX(u,v)
  IF (ind == 0) THEN
    IF (t == 0.0) GO TO 40
    u1 = u / t
    v1 = v / t
    r = u * u1 + v * v1
    IF (u1 >= r*xpos.OR.v1 >= 0.1*r/eps) GO TO 40
    zeta = z1 * z / 1.5
    e = EXP(zeta)
  END IF

  m = 20
  t = MAX(x,ABS(y))
  IF (t > 30.0) m = 8
  s1 = CMPLX(c(m),0.0)
  s2 = CMPLX(d(m),0.0)
  i = m
  DO k = 2, m
    i = i - 1
    s1 = c(i) + w * s1
    s2 = d(i) + w * s2
  END DO

  bi = c1 * z2r * s1
  bip = c1 * z2 * s2
  IF (ind /= 0) RETURN
  bi = e * bi
  bip = e * bip
  RETURN
END IF
IF (x >= 0.0) THEN

!          ----  PI/6 < ABS(ARG(Z)) <= PI/2  ----

  cz = z
  IF (y < 0.0) cz = CONJG(cz)
  zz = cz * ex3c
  z1 = SQRT(zz)
  z2 = SQRT(z1)
  z2r = 1.0 / z2
  cf1 = c1 * z2r * ex6
  cf2 = c1 * z2 * ex6c
  CALL crec(REAL(zz),AIMAG(zz),u,v)
  w = 1.5 * CMPLX(u,v) / z1
  u = ABS(REAL(w))
  v = ABS(AIMAG(w))
  t = MAX(u,v)

  IF (t == 0.0) GO TO 40
  u1 = u / t
  v1 = v / t
  r = u * u1 + v * v1
  IF (ind == 0) THEN
    IF (v1 >= r*xpos.OR.u1 >= 0.1*r/eps) GO TO 40
    zeta = z1 * zz / 1.5
    cn = COS(zeta-cln2)
    si = SIN(zeta-cln2)
    GO TO 20
  END IF

!        E = EXP(-2*I*(ZETA - CLN2)) IF ABS(ARG(ZZ)) <= PI/3
!        E = EXP( 2*I*(ZETA - CLN2)) IF ABS(ARG(ZZ)) > PI/3

  e = (0.0,0.0)
  j = (0.0,-1.0)
  s = 1.0
  ce = 1.0
  cf = 0.5
  IF (AIMAG(zz) > 0.0) THEN
    s = -1.0
    ce = 0.5
    cf = 2.0
  END IF
  IF (v1 < 0.5*r*ABS(xneg)) THEN
    IF (u1 >= 0.05*r/eps) GO TO 40
    zeta = z1 * zz / 1.5
    e = cf * EXP(2.0*s*j*zeta)
  END IF
  cn = ce * c3 * (1+e)
  si = ce * s * c3 * (1-e) * j
ELSE

!                  ----  REAL(Z) < 0  ----

  zz = -z
  IF (y < 0.0) zz = CONJG(zz)
  z1 = SQRT(zz)
  z2 = SQRT(z1)
  z2r = 1.0 / z2
  cf1 = c2 * z2r
  cf2 = c2 * z2
  CALL crec(REAL(zz),AIMAG(zz),u,v)
  w = 1.5 * CMPLX(u,v) / z1
  u = ABS(REAL(w))
  v = ABS(AIMAG(w))
  t = MAX(u,v)

  IF (t == 0.0) GO TO 40
  u1 = u / t
  v1 = v / t
  r = u * u1 + v * v1
  IF (ind == 0) THEN
    IF (v1 >= r*xpos.OR.u1 >= 0.1*r/eps) GO TO 40
    zeta = z1 * zz / 1.5
    cn = COS(zeta)
    si = SIN(zeta)
  ELSE
    e = (0.0,0.0)
    j = (0.0,-1.0)
    IF (v1 < 0.5*r*ABS(xneg)) THEN
      IF (u1 >= 0.05*r/eps) GO TO 40
      zeta = z1 * zz / 1.5
      e = EXP(2.0*j*zeta)
    END IF
    cn = 0.5 * (1.0+e)
    si = 0.5 * (1.0-e) * j
  END IF
END IF

20 w2 = w * w
m = 15
t = MAX(ABS(x),ABS(y))
IF (t > 30.0) m = 5
m2 = m + m
i = m2 - 1
s1 = CMPLX(c(i),0.0)
s2 = CMPLX(c(m2),0.0)
s3 = CMPLX(d(i),0.0)
s4 = CMPLX(d(m2),0.0)
DO k = 2, m
  i = i - 1
  s2 = c(i) - s2 * w2
  s4 = d(i) - s4 * w2
  i = i - 1
  s1 = c(i) - s1 * w2
  s3 = d(i) - s3 * w2
END DO
s2 = w * s2
s4 = w * s4
IF (x < 0.0) THEN
  alpha = s1 + s2
  beta = s2 - s1
ELSE
  alpha = s1 - s2
  beta = s1 + s2
END IF
bi = cf1 * (alpha*cn+beta*si)
IF (x < 0.0) THEN
  alpha = s3 - s4
  beta = s3 + s4
ELSE
  alpha = s3 + s4
  beta = s4 - s3
END IF
bip = cf2 * (alpha*cn+beta*si)
IF (y >= 0.0) RETURN
bi = CONJG(bi)
bip = CONJG(bip)
RETURN

!            RETURN WITH ZERO VALUES IF SCALING IS NEEDED.

40 bi = (0.0,0.0)
bip = (0.0,0.0)
ierr = 1
RETURN
END SUBROUTINE bia



SUBROUTINE ia(z, i1, i2, i1m, i2m)
!-------------------------------------------------------------
!     CALCULATES THE MODIFIED BESSEL FUNCTION OF THE FIRST
!     KIND FOR ORDERS 1/3, 2/3, -1/3, AND -2/3 AND FOR COMPLEX
!     ARGUMENT Z, WHERE -PI < ARG(Z) <= PI.  I1 AND I2
!     ARE REPLACED BY THE FUNCTIONS OF ORDERS 1/3 AND 2/3,
!     RESPECTIVELY, AND I1M AND I2M BY THE FUNCTIONS OF ORDERS
!     -1/3 AND -2/3, RESPECTIVELY.
!-------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: i1, i2, i1m, i2m

! Local variables
COMPLEX, PARAMETER :: ex13 = (5.0E-01, 8.66025403784439E-01),   &
                      ex13c = (5.0E-01, -8.66025403784439E-01), &
                      ex23 = (-5.0E-01, 8.66025403784439E-01),  &
                      ex23c = (-5.0E-01, -8.66025403784439E-01)
COMPLEX            :: cz

!     EX13 = EXP(I*PI/3)
!     EX13C = EXP(-I*PI/3)
!     EX23 = EXP(2*I*PI/3)
!     EX23C = EXP(-2*I*PI/3)

IF (REAL(z) < 0.0) THEN
  cz = -z

!     CALCULATION OF I1, I2, I1M, AND I2M WHEN REAL(CZ) > 0.

  CALL imc(cz, i1, i2, i1m, i2m)

!     FINAL ASSEMBLY

  IF (AIMAG(z) >= 0.0) THEN
    i1 = ex13 * i1
    i2 = ex23 * i2
    i1m = ex13c * i1m
    i2m = ex23c * i2m
    RETURN
  END IF
  i1 = ex13c * i1
  i2 = ex23c * i2
  i1m = ex13 * i1m
  i2m = ex23 * i2m
  RETURN
END IF

CALL imc(z, i1, i2, i1m, i2m)
RETURN
END SUBROUTINE ia


SUBROUTINE imc(z, i1, i2, i1m, i2m)
!----------------------------------------------------------------
!     CALCULATES THE MODIFIED BESSEL FUNCTION OF THE FIRST KIND
!     FOR ORDERS 1/3, 2/3, -1/3, AND -2/3 AND FOR COMPLEX ARGUMENT Z.
!     THE MACLAURIN EXPANSION AND BACKWARD RECURRENCE ARE USED.
!     I1 AND I2 ARE REPLACED BY THE FUNCTIONS OF ORDERS 1/3 AND 2/3,
!     RESPECTIVELY, AND I1M AND I2M BY THE FUNCTIONS OF ORDERS -1/3
!     AND -2/3, RESPECTIVELY.  FOR GREATEST ACCURACY,
!     Z SHOULD LIE IN THE REGION REAL(Z) >= 0.
!----------------------------------------------------------------
COMPLEX, INTENT(IN)            :: z
COMPLEX, INTENT(OUT)           :: i1, i2
COMPLEX, INTENT(OUT), OPTIONAL :: i1m, i2m

! Local variables
COMPLEX         :: ia1, ia2, ia3, ib1, ib2, ib3, sz, zh, e, cf1, cf2, cf3, cf4
REAL, PARAMETER :: c1 = .333333333333333E+00, c2 = .666666666666667E+00,  &
                   gm1 = .892979511569248E+00, gm2 = .902745292950932E+00
REAL            :: a, an, b, cfa, cfb, cn1, cn2, m
INTEGER         :: i, n, n1

!     GM1 = GAMMA(4.0/3.0)
!     GM2 = GAMMA(5.0/3.0)

zh = 0.5 * z
sz = zh * zh
a = REAL(zh)
b = AIMAG(zh)
an = AINT(a*a+b*b)
cn1 = c1 + an
cn2 = c2 + an

!     CALCULATION OF INITIAL VALUES FOR BACKWARD RECURRENCE BY
!     USE OF THE MACLAURIN EXPANSION.

CALL bim(z, cn1, ia1)
CALL bim(z, cn1+1.0, ia2)
CALL bim(z, cn2, ib1)
CALL bim(z, cn2+1.0, ib2)

!     BACKWARD RECURRENCE

n = an
m = an
n1 = n + 1
DO i = 1, n1
  ia3 = ia2
  ia2 = ia1
  ib3 = ib2
  ib2 = ib1
  cfa = (m+c1) * (m+c1+1.0)
  cfb = (m+c2) * (m+c2+1.0)
  m = m - 1.0
  ia1 = ia2 + (sz/cfa) * ia3
  ib1 = ib2 + (sz/cfb) * ib3
END DO
e = EXP(c1*LOG(zh))
cf1 = e / gm1
cf2 = e * e / gm2
i1 = cf1 * ia2
i2 = cf2 * ib2

IF (PRESENT(i1m)) THEN
  cf3 = c2 * cf2 / zh
  i1m = cf3 * ib1
END IF

IF (PRESENT(i2m)) THEN
  cf4 = c1 * cf1 / zh
  i2m = cf4 * ia1
END IF

RETURN
END SUBROUTINE imc


SUBROUTINE bim(z, cn, w)
!-------------------------------------------------------------
!     CALCULATES THE MODIFIED BESSEL FUNCTION OF THE FIRST KIND
!     FOR REAL ORDER CN > -1 AND COMPLEX ARGUMENT Z BY MEANS OF THE
!     MACLAURIN EXPANSION.  W IS REPLACED BY THE CALCULATED VALUE.
!-------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z
REAL, INTENT(IN)     :: cn
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL    :: d, eps, m
COMPLEX :: sz, t
!------------------
!     EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE SMALLEST
!     NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

!-------------------------------------------------------------
sz = 0.25 * z * z

!     INITIALIZATION OF MACLAURIN EXPANSION

m = 1.0
t = sz / (cn+1.0)
w = t

!     SUMMATION OF MACLAURIN EXPANSION

10 m = m + 1.0
d = m * (cn+m)
t = t * (sz/d)
w = w + t
IF (anorm(t) > eps*anorm(w)) GO TO 10

w = w + 1.0
RETURN
END SUBROUTINE bim


SUBROUTINE ja(z, i1, i2, i1m, i2m)
!------------------------------------------------------------
!     CALCULATES THE BESSEL FUNCTION OF THE FIRST KIND FOR
!     ORDERS 1/3, 2/3, -1/3, AND -2/3 AND FOR COMPLEX ARGUMENT
!     Z, WHERE -PI < ARG(Z) <= PI.  I1 AND I2 ARE REPLACED
!     BY THE FUNCTIONS OF ORDERS 1/3 AND 2/3, RESPECTIVELY, AND
!     I1M AND I2M BY FUNCTIONS OF ORDERS -1/3 AND -2/3, RESPECTIVELY.
!--------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: i1, i2, i1m, i2m

! Local variables
COMPLEX, PARAMETER :: ex13  = ( 5.0E-01, 8.66025403784439E-01),  &
                      ex13c = ( 5.0E-01,-8.66025403784439E-01),  &
                      ex23  = (-5.0E-01, 8.66025403784439E-01),  &
                      ex23c = (-5.0E-01,-8.66025403784439E-01)
COMPLEX :: cz

!     EX13 = EXP(I*PI/3)
!     EX13C = EXP(-I*PI/3)
!     EX23 = EXP(2*I*PI/3)
!     EX23C = EXP(-2*I*PI/3)

IF (REAL(z) < 0.0) THEN
  cz = -z

!     CALCULATION OF I1, I2, I1M, AND I2M WHEN REAL(CZ) > 0.0

  CALL jmc(cz,i1,i2,i1m,i2m)

!     FINAL ASSEMBLY

  IF (AIMAG(z) >= 0.0) THEN
    i1 = ex13 * i1
    i2 = ex23 * i2
    i1m = ex13c * i1m
    i2m = ex23c * i2m
    RETURN
  END IF
  i1 = ex13c * i1
  i2 = ex23c * i2
  i1m = ex13 * i1m
  i2m = ex23 * i2m
  RETURN
END IF
CALL jmc(z,i1,i2,i1m,i2m)
RETURN
END SUBROUTINE ja


SUBROUTINE jmc(z, i1, i2, i1m, i2m)
!----------------------------------------------------------------
!     CALCULATES THE BESSEL FUNCTION OF THE FIRST
!     KIND FOR ORDERS 1/3, 2/3, -1/3, AND -2/3 AND FOR COMPLEX
!     ARGUMENT Z. THE MACLAURIN EXPANSION AND BACKWARD RECURRENCE
!     ARE USED. I1 AND I2 ARE REPLACED BY THE FUNCTIONS OF ORDERS
!     1/3 AND 2/3, RESPECTIVELY, AND I1M AND I2M BY THE FUNCTIONS
!     OF ORDERS -1/3 AND -2/3, RESPECTIVELY.  FOR GREATEST
!     ACCURACY, Z SHOULD LIE IN THE REGION REAL(Z) >= 0.
!----------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: i1, i2, i1m, i2m

! Local variables
COMPLEX :: ia1, ia2, ia3, ib1, ib2, ib3, sz, zh, e, cf1, cf2, cf3, cf4
REAL    :: a, an, b, cfa, cfb, cn1, cn2, m
INTEGER :: i, n, n1
REAL, PARAMETER :: c1 = .333333333333333E+00, c2 = .666666666666667E+00, &
                   gm1 = .892979511569248E+00, gm2 = .902745292950932E+00

!     GM1 = GAMMA(4.0/3.0)
!     GM2 = GAMMA(5.0/3.0)

zh = 0.5 * z
sz = zh * zh
a = REAL(zh)
b = AIMAG(zh)
an = AINT(a*a + b*b)
cn1 = c1 + an
cn2 = c2 + an

!     CALCULATION OF INITIAL VALUES FOR BACKWARD RECURRENCE BY
!     USE OF THE MACLAURIN EXPANSION.

CALL bjm(z,cn1,ia1)
CALL bjm(z,cn1+1.0,ia2)
CALL bjm(z,cn2,ib1)
CALL bjm(z,cn2+1.0,ib2)

!     BACKWARD RECURRENCE

n = an
n1 = n + 1
m = an
DO i = 1, n1
  ia3 = ia2
  ia2 = ia1
  ib3 = ib2
  ib2 = ib1
  cfa = (m+c1) * (m+c1+1.0)
  cfb = (m+c2) * (m+c2+1.0)
  m = m - 1.0
  ia1 = ia2 - (sz/cfa) * ia3
  ib1 = ib2 - (sz/cfb) * ib3
END DO
e = EXP(c1*LOG(zh))
cf1 = e / gm1
cf2 = e * e / gm2
cf3 = c2 * cf2 / zh
cf4 = c1 * cf1 / zh
i1 = cf1 * ia2
i2 = cf2 * ib2
i1m = cf3 * ib1
i2m = cf4 * ia1
RETURN
END SUBROUTINE jmc


SUBROUTINE bjm(z, cn, w)
!-------------------------------------------------------------
!     CALCULATES THE BESSEL FUNCTION OF THE FIRST KIND
!     FOR REAL ORDER CN > -1 AND COMPLEX ARGUMENT Z BY MEANS
!     OF THE MACLAURIN EXPANSION.  W IS REPLACED BY THE
!     CALCULATED VALUE.
!-------------------------------------------------------------
COMPLEX, INTENT(IN)  :: z
REAL, INTENT(IN)     :: cn
COMPLEX, INTENT(OUT) :: w

! Local variables
REAL    :: d, eps, m
COMPLEX :: sz, t
!------------------
!     EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE SMALLEST
!     NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

!-------------------------------------------------------------
sz = -0.25 * z * z

!     INITIALIZATION OF MACLAURIN EXPANSION

m = 1.0
t = sz / (cn+1.0)
w = t

!     SUMMATION OF MACLAURIN EXPANSION

10 m = m + 1.0
d = m * (cn+m)
t = t * (sz/d)
w = w + t
IF (anorm(t) > eps*anorm(w)) GO TO 10

w = w + 1.0
RETURN
END SUBROUTINE bjm


SUBROUTINE ka(ind, z, k1, k2)
!------------------------------------------------------------
!     CALCULATES THE MODIFIED BESSEL FUNCTION OF THE SECOND
!     KIND FOR ORDERS 1/3 AND 2/3 AND FOR COMPLEX ARGUMENT Z,
!     WHERE -PI < ARG(Z) <= PI.  K1 IS REPLACED BY THE
!     FUNCTION OF ORDER 1/3, AND K2 BY THE FUNCTION OF ORDER 2/3.
!------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: k1, k2

! Local variables
COMPLEX            :: i1, i2, cz, e
COMPLEX, PARAMETER :: ex13c = (5.0E-01, -8.66025403784439E-01),  &
                      ex23c = (-5.0E-01,-8.66025403784439E-01), j = (0.0,1.0)
REAL, PARAMETER    :: pi = 3.14159265358979E+00
REAL               :: a, b
INTEGER            :: ind1

!     EX13C = EXP(-PI*I/3)
!     EX23C = EXP(-2*PI*I/3)

a = REAL(z)
b = AIMAG(z)
IF (ABS(b) >= -0.5*a) THEN

!            ----  ABS(B) >= -0.5*A  ----

  CALL kml(ind,z,k1,k2)
  RETURN
END IF

!            ----  ABS(B) < -0.5*A  ----

cz = -z
IF (AIMAG(z) < 0.0) cz = CONJG(cz)
ind1 = 0
CALL kml(ind1, cz, k1, k2)
CALL imc(cz, i1, i2)
k1 = ex13c * k1 - j * pi * i1
k2 = ex23c * k2 - j * pi * i2
IF (ind /= 0) THEN
  e = EXP(z)
  k1 = k1 * e
  k2 = k2 * e
END IF
IF (AIMAG(z) >= 0.0) RETURN
k1 = CONJG(k1)
k2 = CONJG(k2)

RETURN
END SUBROUTINE ka


SUBROUTINE kml(ind, z, k1, k2)
!------------------------------------------------------------
!     CALCULATES THE MODIFIED BESSEL FUNCTION OF THE SECOND
!     KIND FOR ORDERS 1/3 AND 2/3 AND FOR COMPLEX ARGUMENT Z
!     BY USE OF THE MILLER ALGORITHM.  K1 IS REPLACED BY THE
!     FUNCTION OF ORDER 1/3, AND K2 BY THE FUNCTION OF ORDER
!     2/3.  FOR GREATEST ACCURACY, Z SHOULD LIE IN THE REGION
!     REAL(Z) >= 0.
!------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: k1, k2

! Local variables
COMPLEX         :: bi, u1, u2, u3, s, e
REAL, PARAMETER :: c1 = 1.25331413731550E+00
REAL            :: a, ai, al, b, c, eps, r, th, x1, x2
INTEGER         :: i, l, m

!     C1 = SQRT(PI/2)

eps = EPSILON(1.0)
x1 = REAL(z)
x2 = AIMAG(z)

!     CALCULATION OF M FOR USE IN MILLER ALGORITHM.

CALL capo(x1,x2,r,th)
a = 3.0 / (1.0+r)
b = 14.7 / (28.0+r)
c = 2.0 / (c1*eps*(2.0*r)**(0.25))
m = (0.485/r) * (LOG(c)+r*COS(a*th)/(1.0+0.008*r)) ** 2 / (2.0*  &
    COS(b*th)) ** 2 + 1.5

!     BACKWARD RECURRENCE IN MILLER ALGORITHM.

s = 0.0
u2 = 0.0
u1 = eps
l = m
DO i = 1, m
  al = l
  u3 = u2
  u2 = u1
  ai = ((al-0.5)**2-1.0/9.0) / (al*(al+1.0))
  bi = 2.0 * (al+z) / (al+1.0)
  u1 = (bi*u2-u3) / ai
  s = s + u1
  l = l - 1
END DO

!     FINAL ASSEMBLY

k1 = c1 * u1 / (s*SQRT(z))
k2 = k1 * (z+1.0/6.0-u2/u1) / z
IF (ind /= 0) RETURN
e = EXP(-z)
k1 = k1 * e
k2 = k2 * e
RETURN
END SUBROUTINE kml


FUNCTION ai(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF THE AIRY FUNCTION AI(X)
!-----------------------------------------------------------------------
!     X0 = 2**(2/3)
!     C = EXP(2/3)
!-----------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER  :: x0 = 1.58740105196820, c = 1.94773404105468, an0  &
        = .355028053887818E+00, an1 = -.187394912983414E+00, an2  &
        = -.383735973881972E-01, an3 = .491952571236878E-01, an4  &
        = -.967017625191329E-02, an5 = -.205648610308316E-02, an6  &
        = .114176040526844E-02, an7 = -.117114823456866E-03, an8  &
        = -.270165470074755E-04, an9 = .789002965889206E-05, ad0  &
        = .100000000000000E+01, ad1 = .201179850513612E+00, ad2  &
        = .385762517106249E-01, ad3 = .230887443780120E-04, bn0  &
        = .355028053887817E+00, bn1 = -.997169317338190E-01, bn2  &
        = -.602216060213075E-01, bn3 = .297705337630730E-01, bn4  &
        = -.152969932286570E-02, bn5 = -.147868368189372E-02, bn6  &
        = .350518617006107E-03, bn7 = -.257766924610873E-04, bd0  &
        = .100000000000000E+01, bd1 = .448140563306831E+00, bd2  &
        = .157074537566686E+00, bd3 = .316964519364865E-01, bd4  &
        = .485922740843953E-02, bd5 = .423326964456309E-03, pn0  &
        = .282094378896566E+00, pn1 = .807868561687271E-01, pn2  &
        = .630644564152247E-02, pn3 = .147116711467936E-03, pn4  &
        = .750490748341483E-06, pd0 = .100000000000000E+01, pd1  &
        = .292890323271551E+00, pd2 = .239376862143358E-01, pd3  &
        = .612353984250624E-03, pd4 = .384461189764830E-05, pd5  &
        = .123247804102182E-08, qn0 = .282094791017188E+00, qn1  &
        = .149585822742689E+00, qn2 = .241876418864958E-01, qn3  &
        = .138190913282142E-02, qn4 = .241862862465003E-04, qn5  &
        = .709733720554615E-07, qd0 = .100000000000000E+01, qd1  &
        = .536778341756648E+00, qd2 = .889112579703465E-01, qd3  &
        = .533368703697049E-02, qd4 = .103812739863315E-03, qd5  &
        = .408838544650398E-06, rn0 = .282094791773878E+00, rn1  &
        = .203731967781874E+00, rn2 = .436660479870037E-01, rn3  &
        = .306595563073142E-02, rn4 = .517398800281618E-04, rd0  &
        = .100000000000000E+01, rd1 = .728721438361672E+00, rd2  &
        = .159210021472267E+00, rd3 = .116985268534248E-01, rd4  &
        = .225973894323078E-03, rd5 = .232707159780478E-06
REAL    :: phi, r, rtx, t, w
INTEGER :: n, n2
!-----------------------------------------------------------------------
IF (x < -1.0) THEN
  CALL aimp(-x,r,phi)
  fn_val = r * SIN(phi)
  RETURN
END IF

IF (x < 0.0) THEN
  fn_val = (((((((((an9*x + an8)*x + an7)*x + an6)*x + an5)*x + an4)*x + &
           an3)*x + an2)*x + an1)*x + an0) / (((ad3*x + ad2)*x + ad1)*x + ad0)
  RETURN
END IF

IF (x < 1.0) THEN
  fn_val = (((((((bn7*x + bn6)*x + bn5)*x + bn4)*x + bn3)*x + bn2)*x + bn1)*x + bn0)  &
           / (((((bd5*x + bd4)*x + bd3)*x + bd2)*x + bd1)*x + bd0)
  RETURN
END IF

rtx = SQRT(x)
IF (x <= x0) THEN
  t = 16.0 / (x*rtx)
  w = ((((pn4*t + pn3)*t + pn2)*t + pn1)*t + pn0) /  &
      (((((pd5*t + pd4)*t + pd3)*t + pd2)*t + pd1)*t + pd0)
  fn_val = (w/SQRT(rtx)) * EXP(-2.0*x*rtx/3.0)
  RETURN
END IF

IF (x <= 4.0D0) THEN
  t = 16.0 / (x*rtx)
  w = (((((qn5*t+qn4)*t+qn3)*t+qn2)*t+qn1)*t+qn0) /   &
      (((((qd5*t+qd4)*t+qd3)*t+qd2)*t+qd1)*t+qd0)
  fn_val = (w/SQRT(rtx)) * EXP(-2.0*x*rtx/3.0)
  RETURN
END IF

IF (x*rtx <= 1.5*exparg(0)) THEN
  t = 16.0 / (x*rtx)
  w = ((((rn4*t + rn3)*t + rn2)*t + rn1)*t + rn0) /   &
      (((((rd5*t + rd4)*t + rd3)*t + rd2)*t + rd1)*t + rd0)
  n = rtx
  n2 = n * n
  t = (x-n2) / (rtx+n)
  fn_val = ((w/SQRT(rtx))/c**(n2*n)) * EXP(-2.0*t*(n*rtx+t*t/3.0))
  RETURN
END IF

fn_val = 0.0
RETURN
END FUNCTION ai


FUNCTION aie(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!                   SCALED AIRY FUNCTION AI(X)


!             AIE(X) = EXP(ZETA)*AI(X)  WHEN X >= 0
!             AIE(X) = AI(X)            WHEN X < 0

!             ZETA = (2/3) * X**(3/2)

!-----------------------------------------------------------------------
!     X0 = 2**(2/3)
!-----------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: x0 = .158740105196820E+01, an0 = .355028053887818E+00, an1  &
        = -.187394912983414E+00, an2 = -.383735973881972E-01, an3  &
        = .491952571236878E-01, an4 = -.967017625191329E-02, an5  &
        = -.205648610308316E-02, an6 = .114176040526844E-02, an7  &
        = -.117114823456866E-03, an8 = -.270165470074755E-04, an9  &
        = .789002965889206E-05, ad0 = .100000000000000E+01, ad1  &
        = .201179850513612E+00, ad2 = .385762517106249E-01, ad3  &
        = .230887443780120E-04, bn0 = .355028053887817E+00, bn1  &
        = -.997169317338190E-01, bn2 = -.602216060213075E-01, bn3  &
        = .297705337630730E-01, bn4 = -.152969932286570E-02, bn5  &
        = -.147868368189372E-02, bn6 = .350518617006107E-03, bn7  &
        = -.257766924610873E-04, bd0 = .100000000000000E+01, bd1  &
        = .448140563306831E+00, bd2 = .157074537566686E+00, bd3  &
        = .316964519364865E-01, bd4 = .485922740843953E-02, bd5  &
        = .423326964456309E-03, pn0 = .282094378896566E+00, pn1  &
        = .807868561687271E-01, pn2 = .630644564152247E-02, pn3  &
        = .147116711467936E-03, pn4 = .750490748341483E-06, pd0  &
        = .100000000000000E+01, pd1 = .292890323271551E+00, pd2  &
        = .239376862143358E-01, pd3 = .612353984250624E-03, pd4  &
        = .384461189764830E-05, pd5 = .123247804102182E-08, qn0  &
        = .282094791017188E+00, qn1 = .149585822742689E+00, qn2  &
        = .241876418864958E-01, qn3 = .138190913282142E-02, qn4  &
        = .241862862465003E-04, qn5 = .709733720554615E-07, qd0  &
        = .100000000000000E+01, qd1 = .536778341756648E+00, qd2  &
        = .889112579703465E-01, qd3 = .533368703697049E-02, qd4  &
        = .103812739863315E-03, qd5 = .408838544650398E-06, rn0  &
        = .282094791773878E+00, rn1 = .203731967781874E+00, rn2  &
        = .436660479870037E-01, rn3 = .306595563073142E-02, rn4  &
        = .517398800281618E-04, rd0 = .100000000000000E+01, rd1  &
        = .728721438361672E+00, rd2 = .159210021472267E+00, rd3  &
        = .116985268534248E-01, rd4 = .225973894323078E-03, rd5  &
        = .232707159780478E-06
REAL   :: phi, r, rtx, t, w
!-----------------------------------------------------------------------
IF (x < -1.0) THEN
  CALL aimp(-x,r,phi)
  fn_val = r * SIN(phi)
  RETURN
END IF

IF (x < 0.0) THEN
  fn_val = (((((((((an9*x + an8)*x + an7)*x + an6)*x + an5)*x + an4)*x +  &
           an3)*x + an2)*x + an1)*x + an0) / (((ad3*x + ad2)*x + ad1)*x + ad0)
  RETURN
END IF

IF (x < 1.0) THEN
  fn_val = (((((((bn7*x + bn6)*x + bn5)*x + bn4)*x + bn3)*x + bn2)*x +  &
           bn1)*x + bn0) /  &
           (((((bd5*x + bd4)*x + bd3)*x + bd2)*x + bd1)*x + bd0)
  IF (x > 1.e-20) fn_val = fn_val * EXP(2.0*x*SQRT(x)/3.0)
  RETURN
END IF

rtx = SQRT(x)
IF (x <= x0) THEN
  t = 16.0 / (x*rtx)
  w = ((((pn4*t + pn3)*t + pn2)*t + pn1)*t + pn0) /   &
      (((((pd5*t + pd4)*t + pd3)*t + pd2)*t + pd1)*t + pd0)
  fn_val = w / SQRT(rtx)
  RETURN
END IF

IF (x <= 4.0D0) THEN
  t = 16.0 / (x*rtx)
  w = (((((qn5*t+qn4)*t+qn3)*t+qn2)*t+qn1)*t+qn0) /   &
      (((((qd5*t+qd4)*t+qd3)*t+qd2)*t+qd1)*t+qd0)
  fn_val = w / SQRT(rtx)
  RETURN
END IF

IF (x <= 1.e20) THEN
  t = 16.0 / (x*rtx)
  w = ((((rn4*t + rn3)*t + rn2)*t + rn1)*t + rn0) / (((((rd5*t + rd4)*t + &
      rd3)*t + rd2)*t + rd1)*t + rd0)
  fn_val = w / SQRT(rtx)
  RETURN
END IF

fn_val = rn0 / SQRT(rtx)
RETURN
END FUNCTION aie


FUNCTION bi(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!             EVALUATION OF THE AIRY FUNCTION BI(X)


!     NOTE... IF X IS A POSITIVE NUMBER WHERE BI(X) IS TOO LARGE
!     TO BE COMPUTED, THEN BI(X) IS SET TO 0.

!-----------------------------------------------------------------------
!     X0 = 16**(2/3)
!     C = EXP(2/3)
!-----------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL    :: x0 = 6.3496042078728, c = 1.94773404105468, an0  &
        = .614926627446001E+00, an1 = .462726943978834E+00, an2  &
        = .867811386408974E-02, an3 = .974670609357959E-01, an4  &
        = .370856545413908E-01, an5 = .569193415071716E-03, an6  &
        = .269172131237236E-02, an7 = .746473849872868E-03, an8  &
        = .105638036899269E-04, an9 = .242726195973978E-04, an10  &
        = .557260250681542E-05, ad0 = .100000000000000E+01, ad1  &
        = .234801779278695E-01, ad2 = -.300487317759152E-02, ad3  &
        = -.597414466459612E-02, bn0 = .614926627446001E+00, bn1  &
        = .548653374523520E+00, bn2 = .582684047163842E-01, bn3  &
        = .871954925712688E-01, bn4 = .508547058449004E-01, bn5  &
        = .361412623711710E-02, bn6 = .177269722794511E-02, bn7  &
        = .117774184027185E-02, bn8 = .627004834186143E-04, bn9  &
        = .774782269814080E-06, bn10 = .118116474369315E-04, bd0  &
        = .100000000000000E+01, bd1 = .163214622184402E+00, bd2  &
        = -.242285981710408E-01, bd3 = -.720554280297616E-02
REAL    :: pn0 = .619911943572678E+00, pn1 = .100411558489626E+01, pn2  &
        = .563659963795768E+00, pn3 = .274925508033015E+00, pn4  &
        = .115641822943246E+00, pn5 = .120048517441127E-01, pn6  &
        = .501838091254330E-02, pd0 = .100000000000000E+01, pd1  &
        = .159751878026937E+01, pd2 = .104664867034140E+01, pd3  &
        = .512560333664022E+00, pd4 = .159144727666995E+00, pd5  &
        = .394456748956258E-01, pd6 = .529926873250079E-02, pd7  &
        = .288921845412576E-03, qn0 = .595123543430856E+00, qn1  &
        = .652692120245803E+00, qn2 = .436851872835894E+00, qn3  &
        = .201626141057807E+00, qn4 = .649535170626944E-01, qn5  &
        = .171798867787816E-01, qn6 = .287998748038892E-02, qn7  &
        = .359634362348937E-03, qd0 = .100000000000000E+01, qd1  &
        = .114259871204893E+01, qd2 = .766390439057101E+00, qd3  &
        = .348287281255683E+00, qd4 = .117049276946157E+00, qd5  &
        = .294545450289541E-01, qd6 = .523951773968125E-02, qd7  &
        = .622692248774973E-03, qd8 = .674811395957744E-06, rn0  &
        = .568067636505865E+00, rn1 = .462183136291541E-01, rn2  &
        = .268519638203645E+00, rn3 = .199427104235673E-02, rn4  &
        = .135599161332010E-03, rn5 = .229937707171804E-04, rn6  &
        = .697888081361175E-05, rn7 = .153277172934286E-05, rn8  &
        = -.149322381877245E-05, rn9 = -.113533571972859E-05, rn10  &
        = .740721412702102E-06, rn11 = -.120160431596119E-06, rd0  &
        = .100000000000000E+01, rd1 = .741293424676788E-01, rd2  &
        = .471695968238457E+00, sn0 = .564189583547757E+00, sn1  &
        = .112605519585866E+00, sn2 = .893329124921909E-03, sn3  &
        = .532139134120350E-04, sn4 = .592725458717738E-05, sn5  &
        = .921448923850546E-06, sn6 = .404558310611815E-06, sn7  &
        = -.660517686759109E-06, sn8 = .174667472383815E-05, sn9  &
        = -.287037710548882E-05, sn10 = .322304072982791E-05, sn11  &
        = -.231569499551950E-05, sn12 = .963478964685941E-06, sn13  &
        = -.173784488565533E-06, sd0 = .100000000000000E+01, sd1  &
        = .193077670156841E+00
REAL    :: phi, r, rtx, t, w
INTEGER :: n, n2
!-----------------------------------------------------------------------
IF (x < -1.0) THEN
  CALL aimp(-x,r,phi)
  fn_val = r * COS(phi)
  RETURN
END IF

IF (x < 0.0) THEN
  fn_val = ((((((((((an10*x + an9)*x + an8)*x + an7)*x + an6)*x + an5)*x +  &
           an4)*x + an3)*x + an2)*x + an1)*x + an0) /  &
           (((ad3*x + ad2)*x + ad1)*x + ad0)
  RETURN
END IF

IF (x <= 1.0) THEN
  fn_val = ((((((((((bn10*x + bn9)*x + bn8)*x + bn7)*x + bn6)*x + bn5)*x +  &
           bn4)*x + bn3)*x + bn2)*x + bn1)*x + bn0) /  &
           (((bd3*x + bd2)*x + bd1)*x + bd0)
  RETURN
END IF

rtx = SQRT(x)
IF (x <= 2.0) THEN
  t = x - 1.0
  w = ((((((pn6*t + pn5)*t + pn4)*t + pn3)*t + pn2)*t + pn1)*t + pn0) /  &
      (((((((pd7*t + pd6)*t + pd5)*t + pd4)*t + pd3)*t + pd2)*t + pd1)*t + pd0)
  fn_val = (w/SQRT(rtx)) * EXP(2.0*x*rtx/3.0)
  RETURN
END IF

IF (x <= 4.0) THEN
  t = x - 2.0
  w = (((((((qn7*t+qn6)*t+qn5)*t+qn4)*t+qn3)*t+qn2)*t+qn1)*t+qn0)  &
      / ((((((((qd8*t+qd7)*t+qd6)*t+qd5)*t+qd4)*t+qd3)*t+qd2)*t+qd1)*t+qd0)
  fn_val = (w/SQRT(rtx)) * EXP(2.0*x*rtx/3.0)
  RETURN
END IF

IF (x <= x0) THEN
  t = 16.0 / (x*rtx) - 1.0
  w = (((((((((((rn11*t + rn10)*t + rn9)*t + rn8)*t + rn7)*t + rn6)*t + rn5)*t  &
      + rn4)*t + rn3)*t + rn2)*t + rn1)*t + rn0) / ((rd2*t + rd1)*t + rd0)
  fn_val = (w/SQRT(rtx)) * EXP(2.0*x*rtx/3.0)
  RETURN
END IF

IF (x*rtx <= 1.5*exparg(0)) THEN
  t = 16.0 / (x*rtx)
  w = (((((((((((((sn13*t + sn12)*t + sn11)*t + sn10)*t + sn9)*t + sn8)*t +  &
      sn7)*t + sn6)*t + sn5)*t + sn4)*t + sn3)*t + sn2)*t + sn1)*t + sn0) / (sd1*t + sd0)
  n = rtx
  n2 = n * n
  t = (x-n2) / (rtx+n)
  fn_val = (w/SQRT(rtx)) * c ** (n2*n) * EXP(2.0*t*(n*rtx+t*t/3.0))
  RETURN
END IF

fn_val = 0.0
RETURN
END FUNCTION bi


FUNCTION bie(x) RESULT(fn_val)
!-----------------------------------------------------------------------

!                    SCALED AIRY FUNCTION BI(X)


!             BIE(X) = EXP(-ZETA)*BI(X)  WHEN X >= 0
!             BIE(X) = BI(X)             WHEN X < 0

!             ZETA = (2/3) * X**(3/2)

!-----------------------------------------------------------------------
!     X0 = 16**(2/3)
!-----------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: x0 = 6.3496042078728, an0 = .614926627446001E+00, an1  &
        = .462726943978834E+00, an2 = .867811386408974E-02, an3  &
        = .974670609357959E-01, an4 = .370856545413908E-01, an5  &
        = .569193415071716E-03, an6 = .269172131237236E-02, an7  &
        = .746473849872868E-03, an8 = .105638036899269E-04, an9  &
        = .242726195973978E-04, an10 = .557260250681542E-05, ad0  &
        = .100000000000000E+01, ad1 = .234801779278695E-01, ad2  &
        = -.300487317759152E-02, ad3 = -.597414466459612E-02, bn0  &
        = .614926627446001E+00, bn1 = .548653374523520E+00, bn2  &
        = .582684047163842E-01, bn3 = .871954925712688E-01, bn4  &
        = .508547058449004E-01, bn5 = .361412623711710E-02, bn6  &
        = .177269722794511E-02, bn7 = .117774184027185E-02, bn8  &
        = .627004834186143E-04, bn9 = .774782269814080E-06, bn10  &
        = .118116474369315E-04, bd0 = .100000000000000E+01, bd1  &
        = .163214622184402E+00, bd2 = -.242285981710408E-01, bd3  &
        = -.720554280297616E-02
REAL, PARAMETER :: pn0 = .619911943572678E+00, pn1  &
        = .100411558489626E+01, pn2 = .563659963795768E+00, pn3  &
        = .274925508033015E+00, pn4 = .115641822943246E+00, pn5  &
        = .120048517441127E-01, pn6 = .501838091254330E-02, pd0  &
        = .100000000000000E+01, pd1 = .159751878026937E+01, pd2  &
        = .104664867034140E+01, pd3 = .512560333664022E+00, pd4  &
        = .159144727666995E+00, pd5 = .394456748956258E-01, pd6  &
        = .529926873250079E-02, pd7 = .288921845412576E-03, qn0  &
        = .595123543430856E+00, qn1 = .652692120245803E+00, qn2  &
        = .436851872835894E+00, qn3 = .201626141057807E+00, qn4  &
        = .649535170626944E-01, qn5 = .171798867787816E-01, qn6  &
        = .287998748038892E-02, qn7 = .359634362348937E-03, qd0  &
        = .100000000000000E+01, qd1 = .114259871204893E+01, qd2  &
        = .766390439057101E+00, qd3 = .348287281255683E+00, qd4  &
        = .117049276946157E+00, qd5 = .294545450289541E-01, qd6  &
        = .523951773968125E-02, qd7 = .622692248774973E-03, qd8  &
        = .674811395957744E-06, rn0 = .568067636505865E+00, rn1  &
        = .462183136291541E-01, rn2 = .268519638203645E+00, rn3  &
        = .199427104235673E-02, rn4 = .135599161332010E-03, rn5  &
        = .229937707171804E-04, rn6 = .697888081361175E-05, rn7  &
        = .153277172934286E-05, rn8 = -.149322381877245E-05, rn9  &
        = -.113533571972859E-05, rn10 = .740721412702102E-06, rn11  &
        = -.120160431596119E-06, rd0 = .100000000000000E+01, rd1  &
        = .741293424676788E-01, rd2 = .471695968238457E+00, sn0  &
        = .564189583547757E+00, sn1 = .112605519585866E+00, sn2  &
        = .893329124921909E-03, sn3 = .532139134120350E-04, sn4  &
        = .592725458717738E-05, sn5 = .921448923850546E-06, sn6  &
        = .404558310611815E-06, sn7 = -.660517686759109E-06, sn8  &
        = .174667472383815E-05, sn9 = -.287037710548882E-05, sn10  &
        = .322304072982791E-05, sn11 = -.231569499551950E-05, sn12  &
        = .963478964685941E-06, sn13 = -.173784488565533E-06, sd0  &
        = .100000000000000E+01, sd1 = .193077670156841E+00
REAL :: phi, r, rtx, t, w
!-----------------------------------------------------------------------
IF (x < -1.0) THEN
  CALL aimp(-x,r,phi)
  fn_val = r * COS(phi)
  RETURN
END IF

IF (x < 0.0) THEN
  fn_val = ((((((((((an10*x + an9)*x + an8)*x + an7)*x + an6)*x + an5)*x + &
           an4)*x + an3)*x + an2)*x + an1)*x + an0) /   &
           (((ad3*x + ad2)*x + ad1)*x + ad0)
  RETURN
END IF

IF (x <= 1.0) THEN
  fn_val = ((((((((((bn10*x + bn9)*x + bn8)*x + bn7)*x + bn6)*x + bn5)*x + &
           bn4)*x + bn3)*x + bn2)*x + bn1)*x + bn0) /  &
           (((bd3*x + bd2)*x + bd1)*x + bd0)
  IF (x > 1.e-20) fn_val = fn_val * EXP(-2.0*x*SQRT(x)/3.0)
  RETURN
END IF

rtx = SQRT(x)
IF (x <= 2.0) THEN
  t = x - 1.0
  w = ((((((pn6*t + pn5)*t + pn4)*t + pn3)*t + pn2)*t + pn1)*t + pn0) /  &
      (((((((pd7*t + pd6)*t + pd5)*t + pd4)*t + pd3)*t + pd2)*t + pd1)*t + pd0)
  fn_val = w / SQRT(rtx)
  RETURN
END IF

IF (x <= 4.0) THEN
  t = x - 2.0
  w = (((((((qn7*t+qn6)*t+qn5)*t+qn4)*t+qn3)*t+qn2)*t+qn1)*t+qn0)  &
      / ((((((((qd8*t+qd7)*t+qd6)*t+qd5)*t+qd4)*t+qd3)*t+qd2)*t+qd1)*t+qd0)
  fn_val = w / SQRT(rtx)
  RETURN
END IF

IF (x <= x0) THEN
  t = 16.0 / (x*rtx) - 1.0
  w = (((((((((((rn11*t + rn10)*t + rn9)*t + rn8)*t + rn7)*t + rn6)*t + rn5)*t  &
      + rn4)*t + rn3)*t + rn2)*t + rn1)*t + rn0) / ((rd2*t + rd1)*t + rd0)
  fn_val = w / SQRT(rtx)
  RETURN
END IF

IF (x <= 1.e20) THEN
  t = 16.0 / (x*rtx)
  w = (((((((((((((sn13*t + sn12)*t + sn11)*t + sn10)*t + sn9)*t + sn8)*t +  &
      sn7)*t + sn6)*t + sn5)*t + sn4)*t + sn3)*t + sn2)*t + sn1)*t + sn0) / &
      (sd1*t + sd0)
  fn_val = w / SQRT(rtx)
  RETURN
END IF

fn_val = sn0 / SQRT(rtx)
RETURN
END FUNCTION bie


SUBROUTINE aimp(x, r, phi)
!-----------------------------------------------------------------------
!     COMPUTATION OF THE AIRY MODULUS AND PHASE FOR X >= 1
!-----------------------------------------------------------------------
REAL, INTENT(IN)  :: x
REAL, INTENT(OUT) :: r, phi

! Local variables
REAL            :: rtx, z
REAL, PARAMETER :: pi4 = .785398163397448, an0 = .297640916735064E+00,   &
               an1 = .772796814419809E+00, an2 = .764990563560236E+00,   &
               an3 = .375694096095838E+00, an4 = .978661044870204E-01,   &
               an5 = .110446639522696E-01, an6 = .145271249611697E-05,   &
               ad0 = .100000000000000E+01, ad1 = .247380029946443E+01,   &
               ad2 = .240125897828762E+01, ad3 = .118267264172257E+01,   &
               ad4 = .306942883081787E+00, ad5 = .347670057203535E-01,   &
               bn0 = .593601051670149E+00, bn1 = .223281495955754E+01,   &
               bn2 = .317718143418600E+01, bn3 = .229890914530923E+01,   &
               bn4 = .933580623665765E+00, bn5 = .209164380960390E+00,   &
               bn6 = .207910965366403E-01, bd0 = .100000000000000E+01,   &
               bd1 = .345985556561483E+01, bd2 = .479629661187354E+01,   &
               bd3 = .345429311552596E+01, bd4 = .140017214942186E+01,   &
               bd5 = .313770549939860E+00, bd6 = .311852186700025E-01,   &
               cn0 = .313541841678871E+00, cn1 = .470104287134296E+00,   &
               cn2 = .291795874641314E+00, cn3 = .962250689852768E-01,   &
               cn4 = .171024484244850E-01, cn5 = .134933201907052E-02,   &
               cd0 = .100000000000000E+01, cd1 = .148070947673639E+01,   &
               cd2 = .917484386216329E+00, cd3 = .302281922152536E+00,   &
               cd4 = .537309296828367E-01, cd5 = .423890576557513E-02,   &
               cd6 = .525954318463502E-08, dn0 = .654836896032068E+00,   &
               dn1 = .117099614856528E+01, dn2 = .831899010444840E+00,   &
               dn3 = .301060337976575E+00, dn4 = .564712748150658E-01,   &
               dn5 = .444134415666317E-02, dd0 = .100000000000000E+01,   &
               dd1 = .176306543768126E+01, dd2 = .124897609613487E+01,   &
               dd3 = .451576491257036E+00, dd4 = .847085955634988E-01,   &
               dd5 = .666188176245820E-02, dd6 = .537600060708764E-08,   &
               pn0 = .318309886183791E+00, pn1 = .100996327221962E+01,   &
               pn2 = .902315148591491E+00, pn3 = .259820640977615E+00,   &
               pn4 = .203717769716282E-01, pn5 = .216893438784765E-03,   &
               pd0 = .100000000000000E+01, pd1 = .317533460265059E+01,   &
               pd2 = .284232123705698E+01, pd3 = .822777439238360E+00,   &
               pd4 = .656865942543526E-01, pd5 = .775376048996392E-03,   &
               qn0 = .666666666666667E+00, qn1 = .141905542385598E+01,   &
               qn2 = .772778148352443E+00, qn3 = .115170415082442E+00,   &
               qn4 = .326457319318373E-02, qd0 = .100000000000000E+01,   &
               qd1 = .213102454203392E+01, qd2 = .116432601041188E+01,   &
               qd3 = .175509465791633E+00, qd4 = .528319849831061E-02,   &
               qd5 = .867802002275824E-05
!-----------------------------------------------------------------------

IF (x <= 2.0) THEN
  z = x - 1.0
  r = ((((((an6*z + an5)*z + an4)*z + an3)*z + an2)*z + an1)*z + an0) /   &
      (((((ad5*z + ad4)*z + ad3)*z + ad2)*z + ad1)*z + ad0)
  phi = ((((((bn6*z + bn5)*z + bn4)*z + bn3)*z + bn2)*z + bn1)*z + bn0) /   &
        ((((((bd6*z + bd5)*z + bd4)*z + bd3)*z + bd2)*z + bd1)*z + bd0)
ELSE

  IF (x < 4.0) THEN
    z = x - 2.0
    r = (((((cn5*z + cn4)*z + cn3)*z + cn2)*z + cn1)*z + cn0) /  &
        ((((((cd6*z+ cd5)*z + cd4)*z + cd3)*z + cd2)*z + cd1)*z + cd0)
    phi = (((((dn5*z + dn4)*z + dn3)*z + dn2)*z + dn1)*z + dn0) /   &
          ((((((dd6*z + dd5)*z + dd4)*z + dd3)*z + dd2)*z + dd1)*z + dd0)
  ELSE

    IF (x <= 1.e10) THEN
      z = 64.0 / x ** 3
      r = (((((pn5*z+pn4)*z+pn3)*z+pn2)*z+pn1)*z+pn0) /  &
          (((((pd5*z +pd4)*z+pd3)*z+pd2)*z+pd1)*z+pd0)
      phi = ((((qn4*z+qn3)*z+qn2)*z+qn1)*z+qn0) /   &
            (((((qd5*z+qd4)*z+qd3)*z+qd2)*z+qd1)*z+qd0)
    ELSE

      r = pn0
      phi = qn0
    END IF
  END IF
END IF

rtx = SQRT(x)
r = SQRT(r/rtx)
phi = pi4 + x * rtx * phi
RETURN
END SUBROUTINE aimp


FUNCTION ck(k, l) RESULT(fn_val)
!     ------------------------------------------------------------------
!     THIS FUNCTION CALCULATES THE COMPLETE ELLIPTIC INTEGRAL F(K)
!     FOR COMPLEX VALUES OF THE MODULUS K. IT IS ASSUMED THAT L.NE.0
!     AND THAT K**2 + L**2 = 1.
!     ------------------------------------------------------------------
COMPLEX, INTENT(IN) :: k, l
COMPLEX             :: fn_val

! Local variables
COMPLEX :: ak, al, ak1, al1, al2, ckk, ckp, f1, f2, f3, fxk,  &
           aktemp, ck1, j = (0.0,1.0), z
REAL    :: ln4 = 1.3862943611199, x1(12) = (/ 6.5487222790801E-03,  &
        3.8946809560450E-02, 9.8150263106007E-02, 1.8113858159063E-01,  &
        2.8322006766737E-01, 3.9843443516344E-01, 5.1995262679235E-01,  &
        6.4051091671611E-01, 7.5286501205183E-01, 8.5024002416230E-01,  &
        9.2674968322391E-01, 9.7775612969000E-01 /),  &
        x2(12) = (/ -9.8156063424672E-01, -9.0411725637048E-01,  &
        -7.6990267419431E-01, -5.8731795428662E-01, -3.6783149899818E-01,  &
        -1.2523340851147E-01, 1.2523340851147E-01, 3.6783149899818E-01,  &
        5.8731795428662E-01, 7.6990267419431E-01, 9.0411725637048E-01,  &
        9.8156063424672E-01 /)
REAL    :: w1(12) = (/ 9.3192691443932E-02,  &
        1.4975182757632E-01, 1.6655745436459E-01, 1.5963355943699E-01,  &
        1.3842483186484E-01, 1.1001657063572E-01, 7.9961821770829E-02,  &
        5.2406954824642E-02, 3.0071088873761E-02, 1.4249245587998E-02,  &
        4.8999245823217E-03, 8.3402903805690E-04 /),  &
        w2(12) = (/ 4.7175336386512E-02, 1.0693932599532E-01,  &
        1.6007832854335E-01, 2.0316742672307E-01, 2.3349253653836E-01,  &
        2.4914704581340E-01, 2.4914704581340E-01, 2.3349253653836E-01,  &
        2.0316742672307E-01, 1.6007832854335E-01, 1.0693932599532E-01,  &
        4.7175336386512E-02 /), fl(12) = (/ 1.5708005371203E+00,  &
        1.5709452753591E+00, 1.5717433742881E+00 , 1.5740325056162E+00,  &
        1.5787613653341E+00, 1.5867393901613E+00, 1.5983969635617E+00,  &
        1.6135762587884E+00, 1.6313677113831E+00, 1.6500349733510E+00,  &
        1.6671202200919E+00, 1.6798403417359E+00 /),  &
        fa(12) = (/ 2.0794472764428E+00, 2.0795966441739E+00,  &
        2.0803359313463E+00, 2.0823286205438E+00, 2.0862633195105E+00,  &
        2.0926508621232E+00, 2.1016440761258E+00, 2.1128974786197E+00,  &
        2.1254857173540E+00, 2.1379218133017E+00, 2.1483404506064E+00,  &
        2.1548934173960E+00 /), fb(12) = (/ 1.5744273529551E+00,  &
        1.5899097325063E+00, 1.6176685384410E+00, 1.6574605448620E+00,  &
        1.7087245795822E+00, 1.7703459462057E+00, 1.8403280188791E+00,  &
        1.9154060277115E+00, 1.9907093877047E+00, 2.0596975322636E+00,  &
        2.1146977530430E+00, 2.1482986855683E+00 /),   &
        c1 = .20264236728467, c2 = .15915494309189
LOGICAL :: branch
REAL    :: eps, phi, r, t, tol, u, v, x, xx, y
INTEGER :: i, ind
!     --------------------------------------------------------------

!     EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE
!     SMALLEST NUMBER SUCH THAT 1. + EPS > 1.

eps = EPSILON(1.0)

!     ---------------------------------------------------
IF (.NOT.(l == (0.0,0.0))) THEN
  ind = 0
  branch = .true.
  tol = 8.0 * MAX(eps,1.e-14)

  ak1 = cflect(k)
  al1 = cflect(l)
  ak = ak1
  al = al1

  x = REAL(ak)
  y = AIMAG(ak)
  u = REAL(al)
  v = AIMAG(al)
  IF (MAX(x,ABS(y)) < 1.0/eps) THEN
    IF (MAX(u,ABS(v)) >= 1.1/eps) GO TO 110

!     CHECK THAT K**2 + L**2 = 1

    IF (x >= u) THEN
      t = u / x
      IF (ABS(x*x/(v*v+1.0/(1.0+t*t))-1.0) > tol) GO TO 110
      IF (ABS(y+t*v) > tol*MAX(1.0,ABS(v))) GO TO 110
    ELSE
      t = x / u
      IF (ABS(u*u/(y*y+1.0/(1.0+t*t))-1.0) > tol) GO TO 110
      IF (ABS(v+t*y) > tol*MAX(1.0,ABS(y))) GO TO 110
    END IF

!     USES  LOGARITHMIC SERIES WHEN ABS(AL) <= 0.55

    IF (u > 1.42.OR.ABS(v) > 1.42) GO TO 80
    10 IF (ABS(al) <= 0.55) THEN
      CALL kl(al,ckk,ckp)
      IF (branch) GO TO 20
      ck1 = ckk
      fn_val = ckp
      al = ak
      GO TO 100
    END IF

!     USES MACLAURIN EXPANSION WHEN THE ABSOLUTE VALUE OF
!     THE MODULUS AK IS LESS THAN OR EQUAL TO 0.55

    r = ABS(ak)
    IF (r > 0.55) GO TO 30
    IF (.NOT.branch) THEN
      CALL kl(ak,ckp,ck1)
      fn_val = ckp
      al = ak
      GO TO 100
    END IF
    ckk = km(ak*ak)
    20 fn_val = ckk
    GO TO 90

!     NUMERICAL QUADRATURE APPROXIMATION

    30 IF (ind /= 0.OR.r <= 1.0) THEN
      40 al2 = al * al

      f1 = (0.0,0.0)
      DO i = 1, 12
        xx = x1(i) / 2.
        fxk = ak * xx
        f1 = f1 + w1(i) * fl(i) / (al2+fxk*fxk)
      END DO
      f2 = (0.0,0.0)
      DO i = 1, 12
        xx = .25 * (1.+x2(i))
        fxk = ak * xx
        f2 = f2 + w2(i) * fa(i) / (al2+fxk*fxk)
      END DO
      f3 = (0.0,0.0)
      DO i = 1, 12
        xx = .25 * (3.-x2(i))
        fxk = ak * xx
        f3 = f3 + w2(i) * fb(i) / (al2+fxk*fxk)
      END DO

      fn_val = al * (c1*f1+c2*(f2+f3))

!     END OF NUMERICAL QUADRATURE APPROXIMATION

      IF (branch) GO TO 90
      ck1 = fn_val
      branch = .true.

!     INTERCHANGE AK AND AL

      aktemp = ak
      ak = al
      al = aktemp
      GO TO 40
    END IF

!     USES INVERSE MODULUS TRANSFORMATION WHEN ABS(AK) IS GREATER
!     THAN 1 AND REAL(AK**2) IS GREATER THAN 0.5.

    80 IF (x*x > y*y+0.5) THEN
      ind = 1
      branch = .false.
      ak = 1.0 / ak1
      al = cflect(j*al1/ak1)
      GO TO 10
    END IF

!     USES COMPLEMENTARY INVERSE MODULUS TRANSFORMATION WHEN ABS(AK)
!     IS GREATER THAN 1 AND REAL(AK**2) IS LESS THAN OR EQUAL TO 0.5

    ind = 2
    ak = cflect(j*ak1/al1)
    al = 1.0 / al1
    GO TO 10

!     RETURN IF NO TRANSFORMATIONS HAVE BEEN PERFORMED

    90 IF (ind == 0) RETURN
    IF (ind /= 1) THEN

!     COMPLEMENTARY INVERSE MODULUS TRANSFORMATION

      fn_val = al * fn_val
      RETURN
    END IF

!     INVERSE MODULUS TRANSFORMATION

    100 IF (AIMAG(ak1) < 0.0) THEN
      fn_val = al * (ck1 - j*fn_val)
      RETURN
    END IF
    fn_val = al * (ck1 + j*fn_val)
    RETURN
  END IF

!     CALCULATION OF F(K) FOR LARGE K AND L

  IF (x > ABS(y)) THEN
    IF (ABS(ABS(v/x)-1.0) > tol) GO TO 110
    IF (ABS(u/x+y/v) > tol) GO TO 110
    t = y / x
    phi = ATAN2(x,ABS(y))
    r = (ln4 + 0.5*alnrel(t*t)) + LOG(x)
    IF (y < 0.0) r = -r
    fn_val = (CMPLX(phi,r)/CMPLX(1.0,t)) / x
    RETURN
  END IF

  IF (ABS(ABS(u/y)-1.0) <= tol) THEN
    IF (ABS(x/u+v/y) <= tol) THEN
      t = v / u
      z = CMPLX((ln4+0.5*alnrel(t*t))+LOG(u),ATAN2(v,u))
      fn_val = (z/CMPLX(1.0,t)) / u
      RETURN
    END IF
  END IF
END IF

!     ERROR RETURN

110 fn_val = (0.0,0.0)
RETURN
END FUNCTION ck



FUNCTION cflect(z) RESULT(fn_val)
!---------------------------------------------------------
!     REFLECTS Z WITH RESPECT TO THE ORIGIN IF REAL(Z)
!     < 0.0 OR IF Z IS ON THE NEGATIVE IMAGINARY AXIS.
!---------------------------------------------------------
COMPLEX, INTENT(IN) :: z
COMPLEX             :: fn_val
!     ----------
IF (REAL(z) < 0.0) THEN
  fn_val = -z
ELSE IF (REAL(z) > 0.0) THEN
  fn_val = z
ELSE
  fn_val = CMPLX(0.0, ABS(AIMAG(z)))
END IF

RETURN
END FUNCTION cflect


FUNCTION km(k2) RESULT(fn_val)
!---------------------------------------------------------------------
!     KM COMPUTES THE COMPLETE ELLIPTIC INTEGRAL F(K) FOR A GIVEN
!     VALUE OF K2 = K**2 BY USE OF THE MACLAURIN EXPANSION.
!---------------------------------------------------------------------
COMPLEX, INTENT(IN) :: k2
COMPLEX             :: fn_val

! Local variables
COMPLEX         :: an, s1
REAL, PARAMETER :: hpi = 1.5707963267949
REAL            :: c, eps, ri, tol
INTEGER         :: i
! ---------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

! ---------------

tol = MAX(eps,1.e-14)
s1 = (1.0,0.0)
an = (1.0,0.0)
DO i = 1, 50
  ri = i
  c = ((ri-0.5)/ri) ** 2
  an = c * (an*k2)
  s1 = s1 + an
  IF (ABS(REAL(an))+ABS(AIMAG(an)) < tol) GO TO 20
END DO

20 fn_val = hpi * s1
RETURN
END FUNCTION km


SUBROUTINE kl(l, fk, fl)
! ----------------------------------------------------------------------
!     KL COMPUTES THE COMPLETE ELLIPTIC INTEGRALS F(K) AND F(L) FOR
!     A GIVEN VALUE OF L, WHERE ABS(L) < 1 AND K**2 + L**2 = 1.
!     IT IS ASSUMED THAT -PI <= ARG(L**2) < PI FOR THE RESULTING
!     VALUE FOR F(K) TO BE MEANINGFUL.
! ----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: l
COMPLEX, INTENT(OUT) :: fk, fl

! Local variables
COMPLEX          :: an, l2, s1, s2, w
REAL, PARAMETER  :: ln4 = 1.3862943611199, hpi = 1.5707963267949
REAL             :: bn, c, eps, ri, tol, u, x, y
INTEGER          :: i
! ---------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

! ---------------

!            THE LOGARITHMIC EXPANSION IS USED FOR F(K)
!              AND THE MACLAURIN EXPANSION FOR F(L)

tol = MAX(eps,1.e-14)
l2 = l * l
s1 = (0.0,0.0)
s2 = (0.0,0.0)
an = (1.0,0.0)
bn = 0.0
DO i = 1, 50
  ri = i
  c = ((ri-0.5)/ri) ** 2
  an = c * (an*l2)
  bn = bn + 1.0 / (ri*(2.0*ri-1.0))
  s1 = s1 + an
  s2 = s2 + an * bn
  IF (ABS(REAL(an))+ABS(AIMAG(an)) < tol) GO TO 20
END DO
20 s1 = s1 + (1.0,0.0)

!                   SET W = 0.5*LOG(16.0/L2)

x = REAL(l)
y = AIMAG(l)
IF (x == 0.0) THEN
  w = CMPLX(ln4-LOG(ABS(y)),hpi)
ELSE

  IF (ABS(x) <= ABS(y)) THEN
    u = (ln4-0.5*alnrel((x/y)**2)) - LOG(ABS(y))
  ELSE
    u = (ln4-0.5*alnrel((y/x)**2)) - LOG(ABS(x))
  END IF

  IF (x <= 0.0) THEN
    w = CMPLX(u,-ATAN2(-y,-x))
  ELSE
    w = CMPLX(u,-ATAN2(y,x))
  END IF
END IF

!                        FINAL ASSEMBLY

fk = w * s1 - s2
fl = hpi * s1
RETURN
END SUBROUTINE kl


SUBROUTINE cke(k, l, ck, ce, ierr)
!     ------------------------------------------------------------------
!     THIS FUNCTION CALCULATES THE COMPLETE ELLIPTIC INTEGRALS F(K)
!     AND E(K) FOR COMPLEX VALUES OF THE MODULUS K.  IT IS ASSUMED
!     THAT L.NE.0 AND THAT K**2 + L**2 = 1.
!     ------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: k, l
COMPLEX, INTENT(OUT) :: ck, ce
INTEGER, INTENT(OUT) :: ierr

COMPLEX :: ak, al, ak1, al1, ckk, ckp, f1, f2, f3, aktemp, ck1
COMPLEX, PARAMETER :: j = (0.0, 1.0)
COMPLEX :: cee, cep, ce1, e1, e2, e3, at, fx, fxk
COMPLEX :: k1, l1, ak2, al2, z, g, g1, gg, gp
REAL, PARAMETER :: ln4 = 1.3862943611199, x1(12) = (/ 6.5487222790801E-03,  &
        3.8946809560450E-02, 9.8150263106007E-02, 1.8113858159063E-01,  &
        2.8322006766737E-01, 3.9843443516344E-01, 5.1995262679235E-01,  &
        6.4051091671611E-01, 7.5286501205183E-01, 8.5024002416230E-01,  &
        9.2674968322391E-01, 9.7775612969000E-01 /),  &
        x2(12) = (/ -9.8156063424672E-01, -9.0411725637048E-01,  &
        -7.6990267419431E-01, -5.8731795428662E-01, -3.6783149899818E-01,  &
        -1.2523340851147E-01, 1.2523340851147E-01, 3.6783149899818E-01,  &
        5.8731795428662E-01, 7.6990267419431E-01, 9.0411725637048E-01,  &
        9.8156063424672E-01 /), w1(12) = (/ 9.3192691443932E-02,  &
        1.4975182757632E-01, 1.6655745436459E-01, 1.5963355943699E-01,  &
        1.3842483186484E-01, 1.1001657063572E-01, 7.9961821770829E-02,  &
        5.2406954824642E-02, 3.0071088873761E-02, 1.4249245587998E-02,  &
        4.8999245823217E-03, 8.3402903805690E-04 /),  &
        w2(12) = (/ 4.7175336386512E-02, 1.0693932599532E-01,  &
        1.6007832854335E-01, 2.0316742672307E-01, 2.3349253653836E-01,  &
        2.4914704581340E-01, 2.4914704581340E-01, 2.3349253653836E-01,  &
        2.0316742672307E-01, 1.6007832854335E-01, 1.0693932599532E-01,  &
        4.7175336386512E-02 /), fl(12) = (/ 1.5708005371203E+00,  &
        1.5709452753591E+00, 1.5717433742881E+00 , 1.5740325056162E+00,  &
        1.5787613653341E+00, 1.5867393901613E+00, 1.5983969635617E+00,  &
        1.6135762587884E+00, 1.6313677113831E+00, 1.6500349733510E+00,  &
        1.6671202200919E+00, 1.6798403417359E+00 /),  &
        fa(12) = (/ 2.0794472764428E+00, 2.0795966441739E+00,  &
        2.0803359313463E+00, 2.0823286205438E+00, 2.0862633195105E+00,  &
        2.0926508621232E+00, 2.1016440761258E+00, 2.1128974786197E+00,  &
        2.1254857173540E+00, 2.1379218133017E+00, 2.1483404506064E+00,  &
        2.1548934173960E+00 /), fb(12) = (/ 1.5744273529551E+00,  &
        1.5899097325063E+00, 1.6176685384410E+00, 1.6574605448620E+00,  &
        1.7087245795822E+00, 1.7703459462057E+00, 1.8403280188791E+00,  &
        1.9154060277115E+00, 1.9907093877047E+00, 2.0596975322636E+00,  &
        2.1146977530430E+00, 2.1482986855683E+00 /),   &
        c1 = .20264236728467, c2 = .15915494309189
REAL    :: c, eps, phi, r, t, tol, u, v, x, xx, y
LOGICAL :: branch
INTEGER :: i, ind
!     --------------------------------------------------------------

!     EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE
!     SMALLEST NUMBER SUCH THAT 1. + EPS > 1.

eps = EPSILON(1.0)

!     ---------------------------------------------------
IF (.NOT.(l == (0.0,0.0))) THEN
  ind = 0
  branch = .true.
  tol = 8.0 * MAX(eps,1.e-14)

  ak1 = cflect(k)
  al1 = cflect(l)
  ak = ak1
  al = al1
  ierr = 0

  x = REAL(ak)
  y = AIMAG(ak)
  u = REAL(al)
  v = AIMAG(al)
  IF (MAX(x,ABS(y)) < 1.0/eps) THEN
    IF (MAX(u,ABS(v)) >= 1.1/eps) GO TO 110

!     CHECK THAT K**2 + L**2 = 1

    IF (x >= u) THEN
      t = u / x
      IF (ABS(x*x/(v*v+1.0/(1.0+t*t))-1.0) > tol) GO TO 110
      IF (ABS(y+t*v) > tol*MAX(1.0,ABS(v))) GO TO 110
    ELSE
      t = x / u
      IF (ABS(u*u/(y*y+1.0/(1.0+t*t))-1.0) > tol) GO TO 110
      IF (ABS(v+t*y) > tol*MAX(1.0,ABS(y))) GO TO 110
    END IF

!     USES  LOGARITHMIC SERIES WHEN ABS(AL)
!     IS LESS THAN OR EQUAL TO 0.55

    IF (u > 1.42.OR.ABS(v) > 1.42) GO TO 80
    10 IF (ABS(al) <= 0.55) THEN
      CALL ekl(al,ckk,ckp,cee,cep,gg,gp)
      IF (branch) GO TO 20
      ck1 = ckk
      ck = ckp
      ce1 = cee
      ce = cep
      g1 = gg
      g = gp
      ak2 = al * al
      al = ak
      al2 = al * al
      GO TO 100
    END IF

!     USES MACLAURIN EXPANSION WHEN THE ABSOLUTE VALUE OF
!     THE MODULUS AK IS LESS THAN OR EQUAL TO 0.55

    r = ABS(ak)
    IF (r > 0.55) GO TO 30
    IF (.NOT.branch) THEN
      CALL ekl(ak,ckp,ck1,cep,ce1,gp,g1)
      ck = ckp
      ce = cep
      g = gp
      ak2 = al * al
      al = ak
      al2 = al * al
      GO TO 100
    END IF
    CALL ekm(ak*ak,ckk,cee)
    20 ck = ckk
    ce = cee
    GO TO 90

!     NUMERICAL QUADRATURE APPROXIMATION

    30 IF (ind /= 0.OR.r <= 1.0) THEN
      40 al2 = al * al
      ak2 = ak * ak

      f1 = (0.0,0.0)
      e1 = (0.0,0.0)
      DO i = 1, 12
        xx = x1(i) / 2.
        fx = ak * xx / al
        fxk = ak * xx
        at = atn(fx)
        e1 = e1 + w1(i) * fl(i) * (1.0+at)
        f1 = f1 + w1(i) * fl(i) / (al2+fxk*fxk)
      END DO
      f2 = (0.0,0.0)
      e2 = (0.0,0.0)
      DO i = 1, 12
        xx = .25 * (1.+x2(i))
        fx = ak * xx / al
        fxk = ak * xx
        at = atn(fx)
        e2 = e2 + w2(i) * fa(i) * (1.0 + at)
        f2 = f2 + w2(i) * fa(i) / (al2 + fxk*fxk)
      END DO
      f3 = (0.0,0.0)
      e3 = (0.0,0.0)
      DO i = 1, 12
        xx = .25 * (3.-x2(i))
        fx = ak * xx / al
        fxk = ak * xx
        at = atn(fx)
        e3 = e3 + w2(i) * fb(i) * (1.0 + at)
        f3 = f3 + w2(i) * fb(i) / (al2 + fxk*fxk)
      END DO

      ck = al * (c1*f1 + c2*(f2+f3))
      ce = al * (c1*e1 + c2*(e2+e3))

!     END OF NUMERICAL QUADRATURE APPROXIMATION

      IF (branch) GO TO 90
      ck1 = ck
      ce1 = ce
      branch = .true.

!     INTERCHANGE AK AND AL

      aktemp = ak
      ak = al
      al = aktemp
      GO TO 40
    END IF

!     USES INVERSE MODULUS TRANSFORMATION WHEN ABS(AK) IS GREATER
!     THAN 1 AND REAL(AK**2) IS GREATER THAN 0.5.

    80 IF (x*x > y*y+0.5) THEN
      ind = 1
      branch = .false.
      ak = 1.0 / ak1
      al = cflect(j*al1/ak1)
      GO TO 10
    END IF

!     USES COMPLEMENTARY INVERSE MODULUS TRANSFORMATION WHEN ABS(AK)
!     IS GREATER THAN 1 AND REAL(AK**2) IS LESS THAN OR EQUAL TO 0.5

    ind = 2
    ak = cflect(j*ak1/al1)
    al = 1.0 / al1
    GO TO 10

!     RETURN IF NO TRANSFORMATIONS HAVE BEEN PERFORMED

    90 IF (ind == 0) RETURN
    IF (ind /= 1) THEN

!     COMPLEMENTARY INVERSE MODULUS TRANSFORMATION

      ck = al * ck
      ce = ce / al
      RETURN
    END IF

!     INVERSE MODULUS TRANSFORMATION

    g = ce - al2 * ck
    g1 = ce1 - ak2 * ck1
    100 IF (AIMAG(ak2) < 0.0) THEN
      ce = (g1+j*g) / al
      ck = al * (ck1-j*ck)
      RETURN
    END IF
    ce = (g1-j*g) / al
    ck = al * (ck1+j*ck)
    RETURN
  END IF

!     CALCULATION OF F(K) AND E(K) FOR LARGE K AND L

  IF (x > ABS(y)) THEN
    IF (ABS(ABS(v/x)-1.0) > tol) GO TO 110
    IF (ABS(u/x+y/v) > tol) GO TO 110
    t = y / x
    k1 = CMPLX(1.0,t)
    phi = ATAN2(x,ABS(y))
    r = (ln4 + 0.5*alnrel(t*t)) + LOG(x)
    c = 0.5 * r + 0.25
    z = CMPLX(y,-x)
    IF (y < 0.0) THEN
      r = -r
      c = -c
      z = -z
    END IF
    ck = (CMPLX(phi,r)/k1) / x
    ce = z + (CMPLX(0.5*phi,c)/k1) / x
    RETURN
  END IF

  IF (ABS(ABS(u/y)-1.0) > tol) GO TO 110
  IF (ABS(x/u+v/y) > tol) GO TO 110
  t = v / u
  l1 = CMPLX(1.0,t)
  r = (ln4 + 0.5*alnrel(t*t)) + LOG(u)
  phi = ATAN2(v,u)
  ck = (CMPLX(r,phi)/l1) / u
  ce = al + (CMPLX(0.5*r-0.25, 0.5*phi)/l1) / u
  RETURN
END IF

!     ERROR RETURN

ierr = 1
RETURN
110 ierr = 2
RETURN
END SUBROUTINE cke



FUNCTION atn(z) RESULT(fn_val)
!---------------------------------------------------
!     CALCULATES COMPLEX FUNCTION ATN(Z) = Z*ATAN(Z)
!     USING REAL (dp).
!---------------------------------------------------
COMPLEX, INTENT(IN) :: z
COMPLEX             :: fn_val

! Local variables
REAL (dp) :: dx, dy
REAL      :: atn1, atn2, d, da, db, t, x, y

x = REAL(z)
y = AIMAG(z)
dx = x
dy = y
t = 1.d0 - dx * dx - dy * dy
da = -0.5 * ATAN2(-2.0*x,t)
d = (1.0-dy) ** 2 + dx * dx
db = 0.25 * alnrel(4.0*y/d)
atn1 = da * x - db * y
atn2 = da * y + db * x
fn_val = CMPLX(atn1,atn2)
RETURN
END FUNCTION atn


SUBROUTINE ekm(k2, fk, ek)
! ----------------------------------------------------------------------
!     EKM COMPUTES THE COMPLETE ELLIPTIC INTEGRALS F(K) AND E(K) FOR
!     A GIVEN VALUE OF K2 = K**2 BY USE OF THE MACLAURIN EXPANSIONS.
! ----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: k2
COMPLEX, INTENT(OUT) :: fk, ek

! Local variables
COMPLEX         :: an, cn, s1, s2
REAL, PARAMETER :: hpi = 1.5707963267949
REAL            :: c, eps, ri, tol
INTEGER         :: i
! ---------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

! ---------------

tol = MAX(eps,1.e-14)
s1 = (1.0,0.0)
s2 = (1.0,0.0)
an = (1.0,0.0)
DO i = 1, 50
  ri = i
  c = ((ri-0.5)/ri) ** 2
  an = c * (an*k2)
  cn = an / (2.0*ri-1.0)
  s1 = s1 + an
  s2 = s2 - cn
  IF (ABS(REAL(an))+ABS(AIMAG(an)) < tol) GO TO 20
END DO

20 fk = hpi * s1
ek = hpi * s2
RETURN
END SUBROUTINE ekm


SUBROUTINE ekl(l, fk, fl, ek, el, gk, gl)
! ----------------------------------------------------------------------
!     EKL COMPUTES THE COMPLETE ELLIPTIC INTEGRALS F(K), F(L), E(K),
!     E(L) FOR A GIVEN VALUE OF L2, WHERE L2 = L**2 AND K**2 + L**2 = 1.
!     IT IS ASSUMED THAT -PI < ARG(L2) <= PI FOR THE RESULTING
!     VALUE FOR F(K) TO BE MEANINGFUL.  THE COMBINATIONS OF FUNCTIONS
!     G(K) = E(K) - L**2*F(K) AND G(L) = E(L) - K**2*F(L) ARE ALSO
!     CALCULATED.
! ----------------------------------------------------------------------
COMPLEX, INTENT(IN)  :: l
COMPLEX, INTENT(OUT) :: fk, fl, ek, el, gk, gl

! Local variables
COMPLEX         :: an, cn, en, l2, s1, s2, s3, s4, s5, s6, s7, w
REAL, PARAMETER :: ln4 = 1.3862943611199, hpi = 1.5707963267949
REAL            :: bn, c, dn, eps, fn, gn, ri, tol, u, x, y
INTEGER         :: i
! ---------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1 + EPS > 1.

eps = EPSILON(1.0)

! ---------------

!            THE LOGARITHMIC EXPANSIONS ARE USED FOR F(K) AND E(K)
!            AND THE MACLAURIN EXPANSIONS FOR F(L) AND E(L)

tol = MAX(eps,1.e-14)
l2 = l * l
s1 = (0.0,0.0)
s2 = (0.0,0.0)
s3 = (0.0,0.0)
s4 = (0.0,0.0)
s5 = (0.0,0.0)
s6 = (0.0,0.0)
s7 = (0.0,0.0)
an = (1.0,0.0)
bn = 0.0
DO i = 1, 300
  ri = i
  c = ((ri-0.5)/ri) ** 2
  an = c * (an*l2)
  bn = bn + 1.0 / (ri*(2.0*ri-1.0))
  cn = an / (2.0*ri-1.0)
  dn = bn * ri / (ri-0.5)
  en = cn / (2.0*ri-1.0)
  fn = ri / (ri-0.5)
  gn = 0.5 / (ri+1.0)
  s1 = s1 + an
  s2 = s2 + an * bn
  s3 = s3 - cn
  s4 = s4 + an * dn
  s5 = s5 + en
  s6 = s6 + an * fn
  s7 = s7 + an * gn
  IF (ABS(REAL(an))+ABS(AIMAG(an)) < tol) GO TO 20
END DO
20 s1 = s1 + (1.0,0.0)
s3 = s3 + (1.0,0.0)
s5 = s5 + (1.0,0.0)
s7 = s7 + (0.5,0.0)

!              SET W = 0.5*LOG(16.0/L2)

x = REAL(l)
y = AIMAG(l)
IF (x == 0.0) THEN
  w = CMPLX(ln4-LOG(ABS(y)),hpi)
ELSE

  IF (ABS(x) <= ABS(y)) THEN
    u = (ln4-0.5*alnrel((x/y)**2)) - LOG(ABS(y))
  ELSE
    u = (ln4-0.5*alnrel((y/x)**2)) - LOG(ABS(x))
  END IF

  IF (x <= 0.0) THEN
    w = CMPLX(u,-ATAN2(-y,-x))
  ELSE
    w = CMPLX(u,-ATAN2(y,x))
  END IF
END IF

!              FINAL ASSEMBLY

fk = w * s1 - s2
fl = hpi * s1
ek = w * s6 - s4 + s5
el = hpi * s3
gk = -w * s7 * l2 - s4 + s5 + s2 * l2
gl = hpi * s7 * l2
RETURN
END SUBROUTINE ekl


SUBROUTINE ellpi(phi, cphi, k, l, f, e, ierr)
!-----------------------------------------------------------------------

!          REAL ELLIPTIC INTEGRALS OF THE FIRST AND SECOND KINDS

!                        -----------------

!     PHI = ARGUMENT                    (0.0 <= PHI  <= PI/2)
!     CPHI = PI/2 - PHI                 (0.0 <= CPHI <= PI/2)
!     K = MODULUS                       (ABS(K) <= 1.0)
!     L = COMODULUS = SQRT (1 - K*K)    (ABS(L) <= 1.0)
!     F = ELLIPTIC INTEGRAL OF FIRST KIND = F(PHI, K)
!     E = ELLIPTIC INTEGRAL OF SECOND KIND = E(PHI, K)
!     IERR = ERROR INDICATOR (IERR = 0  IF NO ERRORS)
!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: phi, cphi, k, l
REAL, INTENT(OUT)    :: f, e
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL, PARAMETER :: ln4 = 1.3862943611199, th1 = .76159415595576
REAL            :: an, cn, dn, d1, d2, d3, d4, hn, k2, l2, p, pn, px, qn, qx,  &
                   r, ri, rj, rk, rm, rn, r0, r2, si, sj, sk, sn, ss,  &
                   s0, s1, s2, s3, s4, td, tr, ts, t1, t2, w
!------------------------
!     LN4 = LN(4)
!     TH1 = TANH(1)
!------------------------
IF (phi >= 0.0 .AND. cphi >= 0.0) THEN
  IF (ABS(k) > 1.0.OR.ABS(l) > 1.0) GO TO 40
  ierr = 0
  IF (phi == 0.0) THEN
    f = 0.0
    e = 0.0
    RETURN
  END IF

  IF (phi >= 0.79) THEN
    sn = COS(cphi)
    cn = SIN(cphi)
  ELSE
    sn = SIN(phi)
    cn = COS(phi)
  END IF

  k2 = k * k
  l2 = l * l
  ss = sn * sn
  px = ABS(k*sn)
  qx = ABS(k*cn)
  IF (px < th1) THEN

!     SERIES EXPANSION FOR ABS(K*SIN(PHI)) <= TANH(1)

    pn = 1.0
    qn = 2.0
    an = phi
    hn = 1.0
    s1 = 0.0
    s2 = 0.0
    tr = phi * ss
    ts = sn * cn

    10 an = (pn*an-ts) / qn
    r = k2 * hn / qn
    s2 = s2 + r * an
    hn = pn * r
    s0 = s1
    s1 = s1 + hn * an
    IF (ABS(tr) >= ABS(an)) THEN
      IF (ABS(s1) > ABS(s0)) THEN
        pn = qn + 1.0
        qn = pn + 1.0
        tr = ss * tr
        ts = ss * ts
        GO TO 10
      END IF
    END IF

    f = phi + s1
    e = phi - s2
    RETURN
  END IF

!     SERIES EXPANSION FOR ABS(K*SIN(PHI)) > TANH(1)

  r = cpabs(l,qx)
  IF (r == 0.0) GO TO 50
  r2 = r * r
  si = 1.0
  sj = 1.0
  sk = 0.0
  rm = 0.0
  rn = 0.0
  s1 = 0.0
  s2 = 0.0
  s3 = 0.0
  s4 = 0.0
  td = qx * r
  dn = 2.0
  GO TO 30

  20 si = ri
  sj = rj
  sk = rk
  dn = dn + 2.0
  td = r2 * td
  30 pn = (dn-1.0) / dn
  qn = (dn+1.0) / (dn+2.0)
  ri = pn * si
  rj = pn * pn * l2 * sj
  rk = sk + 2.0 / (dn*(dn-1.0))
  r0 = td / dn
  rm = qn * qn * l2 * (rm-r0*ri)
  rn = pn * qn * l2 * (rn-r0*si)
  d1 = rj
  d2 = qn * rj
  d3 = rm - rj * rk
  d4 = rn - pn * l2 * sj * rk + l2 * sj / (dn*dn)
  r0 = s3
  s1 = s1 + d1
  s2 = s2 + d2
  s3 = s3 + d3
  s4 = s4 + d4
  IF (s3 < r0) GO TO 20

  w = 1.0 + px
  p = ln4 - LOG(r+qx)
  t1 = (1.0+s1) * p + qx / r * alnrel(-0.5*r2/w)
  t2 = (0.5+s2) * l2 * p + (1.0-qx*r/w)
  f = t1 + s3
  e = t2 + s4
  RETURN
END IF

!     ERROR RETURN

ierr = 1
RETURN
40 ierr = 2
RETURN
50 ierr = 3
RETURN
END SUBROUTINE ellpi


SUBROUTINE dellpi(phi, cphi, k, l, f, e, ierr)
!-----------------------------------------------------------------------

!                  REAL (dp) COMPUTATION OF THE
!          REAL ELLIPTIC INTEGRALS OF THE FIRST AND SECOND KINDS

!                        -----------------

!     PHI = ARGUMENT                    (0.0 <= PHI  <= PI/2)
!     CPHI = PI/2 - PHI                 (0.0 <= CPHI <= PI/2)
!     K = MODULUS                       (ABS(K) <= 1.0)
!     L = COMODULUS = SQRT (1 - K*K)    (ABS(L) <= 1.0)
!     F = ELLIPTIC INTEGRAL OF FIRST KIND = F(PHI, K)
!     E = ELLIPTIC INTEGRAL OF SECOND KIND = E(PHI, K)
!     IERR = ERROR INDICATOR (IERR = 0  IF NO ERRORS)
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: phi, cphi, k, l
REAL (dp), INTENT(OUT) :: f, e
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL (dp), PARAMETER :: ln4 = .1386294361119890618834464242916353136151D+01, &
                        th1 = .7615941559557648881194582826047935904128D+00
REAL (dp) :: an, cn, dn, d1, d2, d3, d4, hn, k2, l2, p, pn, px, qn, qx, r,   &
             ri, rj, rk, rm, rn, r0, r2, si, sj, sk, sn, ss, s0, s1, s2, s3, &
             s4, td, tr, ts, t1, t2, w
!------------------------
!     LN4 = LN(4)
!     TH1 = TANH(1)
!------------------------
IF (phi >= 0.d0 .AND. cphi >= 0.d0) THEN
  IF (ABS(k) > 1.d0.OR.ABS(l) > 1.d0) GO TO 40
  ierr = 0
  IF (phi == 0.d0) THEN
    f = 0.d0
    e = 0.d0
    RETURN
  END IF

  IF (phi >= 0.79D0) THEN
    sn = COS(cphi)
    cn = SIN(cphi)
  ELSE
    sn = SIN(phi)
    cn = COS(phi)
  END IF

  k2 = k * k
  l2 = l * l
  ss = sn * sn
  px = ABS(k*sn)
  qx = ABS(k*cn)
  IF (px < th1) THEN

!     SERIES EXPANSION FOR ABS(K*SIN(PHI)) <= TANH(1)

    pn = 1.d0
    qn = 2.d0
    an = phi
    hn = 1.d0
    s1 = 0.d0
    s2 = 0.d0
    tr = phi * ss
    ts = sn * cn

    10 an = (pn*an-ts) / qn
    r = k2 * hn / qn
    s2 = s2 + r * an
    hn = pn * r
    s0 = s1
    s1 = s1 + hn * an
    IF (ABS(tr) >= ABS(an)) THEN
      IF (ABS(s1) > ABS(s0)) THEN
        pn = qn + 1.d0
        qn = pn + 1.d0
        tr = ss * tr
        ts = ss * ts
        GO TO 10
      END IF
    END IF

    f = phi + s1
    e = phi - s2
    RETURN
  END IF

!     SERIES EXPANSION FOR ABS(K*SIN(PHI)) > TANH(1)

  r = dcpabs(l,qx)
  IF (r == 0.d0) GO TO 50
  r2 = r * r
  si = 1.d0
  sj = 1.d0
  sk = 0.d0
  rm = 0.d0
  rn = 0.d0
  s1 = 0.d0
  s2 = 0.d0
  s3 = 0.d0
  s4 = 0.d0
  td = qx * r
  dn = 2.d0
  GO TO 30

  20 si = ri
  sj = rj
  sk = rk
  dn = dn + 2.d0
  td = r2 * td
  30 pn = (dn-1.d0) / dn
  qn = (dn+1.d0) / (dn+2.d0)
  ri = pn * si
  rj = pn * pn * l2 * sj
  rk = sk + 2.d0 / (dn*(dn-1.d0))
  r0 = td / dn
  rm = qn * qn * l2 * (rm-r0*ri)
  rn = pn * qn * l2 * (rn-r0*si)
  d1 = rj
  d2 = qn * rj
  d3 = rm - rj * rk
  d4 = rn - pn * l2 * sj * rk + l2 * sj / (dn*dn)
  r0 = s3
  s1 = s1 + d1
  s2 = s2 + d2
  s3 = s3 + d3
  s4 = s4 + d4
  IF (s3 < r0) GO TO 20

  w = 1.d0 + px
  p = ln4 - LOG(r+qx)
  t1 = (1.d0+s1) * p + qx / r * dlnrel(-0.5D0*r2/w)
  t2 = (0.5D0+s2) * l2 * p + (1.d0-qx*r/w)
  f = t1 + s3
  e = t2 + s4
  RETURN
END IF

!     ERROR RETURN

ierr = 1
RETURN
40 ierr = 2
RETURN
50 ierr = 3
RETURN
END SUBROUTINE dellpi


SUBROUTINE epi(phi, cphi, k2, l2, n, m, p, ierr)
!-----------------------------------------------------------------------
!             REAL ELLIPTIC INTEGRAL OF THE THIRD KIND
!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: phi, cphi, k2, l2, n, m
REAL, INTENT(OUT)    :: p
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL, PARAMETER :: pihalf = 1.5707963267948966192
REAL            :: a, b, c, eps, r, rf, s, s2, tol
!---------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.0 + EPS > 1.0.

eps = EPSILON(1.0)

!-----------------------------------------------------------------------
tol = 4.0 * eps
IF (MIN(phi,cphi) >= 0.0) THEN
  IF (ABS((phi+cphi)-pihalf) <= tol*pihalf) THEN
    IF (ABS(n) > 1.0) GO TO 10
    IF (k2 < 0.0 .OR. l2 < 0.0) GO TO 20
    IF (ABS((k2+l2)-1.0) > tol) GO TO 20

    IF (phi >= 0.79) THEN
      s = COS(cphi)
      c = SIN(cphi)
    ELSE
      s = SIN(phi)
      c = COS(phi)
    END IF
    a = c * c
    b = l2 + k2 * a
    s2 = s * s

    IF (n <= 0.0) THEN
      r = 1.0 - n * s2
    ELSE
      IF (m < 0.0 .OR. m > 1.0) GO TO 10
      IF (ABS((m+n)-1.0) > tol) GO TO 10
      r = m + n * a
    END IF

    CALL rjval(a,b,1.0,r,p,ierr)
    IF (ierr /= 0) GO TO 30
    p = p * (s*s2) * n / 3.0
    CALL rfval(a,b,1.0,rf,ierr)
    p = p + s * rf
    RETURN
  END IF
END IF

!     ERROR RETURN

ierr = 1
RETURN
10 ierr = 2
RETURN
20 ierr = 3
RETURN
30 ierr = 4
RETURN
END SUBROUTINE epi


SUBROUTINE depi(phi, cphi, k2, l2, n, m, p, ierr)
!-----------------------------------------------------------------------
!               REAL (dp) COMPUTATION OF THE
!             REAL ELLIPTIC INTEGRAL OF THE THIRD KIND
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: phi, cphi, k2, l2, n, m
REAL (dp), INTENT(OUT) :: p
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL (dp)            :: a, b, c, eps, r, rf, s, s2, tol
REAL (dp), PARAMETER :: pihalf = 1.570796326794896619231321691639751442099D0
!-----------------------------------------------------------------------

!     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE
!            SMALLEST NUMBER SUCH THAT 1.D0 + EPS > 1.D0.

eps = EPSILON(1.0D0)

!-----------------------------------------------------------------------
tol = 4.d0 * eps
IF (MIN(phi,cphi) >= 0.d0) THEN
  IF (ABS((phi+cphi)-pihalf) <= tol*pihalf) THEN
    IF (ABS(n) > 1.d0) GO TO 10
    IF (k2 < 0.d0.OR.l2 < 0.d0) GO TO 20
    IF (ABS((k2+l2)-1.d0) > tol) GO TO 20

    IF (phi >= 0.79D0) THEN
      s = COS(cphi)
      c = SIN(cphi)
    ELSE
      s = SIN(phi)
      c = COS(phi)
    END IF
    a = c * c
    b = l2 + k2 * a
    s2 = s * s

    IF (n <= 0.d0) THEN
      r = 1.d0 - n * s2
    ELSE
      IF (m < 0.d0.OR.m > 1.d0) GO TO 10
      IF (ABS((m+n)-1.d0) > tol) GO TO 10
      r = m + n * a
    END IF

    CALL drjval(a,b,1.d0,r,p,ierr)
    IF (ierr /= 0) GO TO 30
    p = p * (s*s2) * n / 3.d0
    CALL drfval(a,b,1.d0,rf,ierr)
    p = p + s * rf
    RETURN
  END IF
END IF

!     ERROR RETURN

ierr = 1
RETURN
10 ierr = 2
RETURN
20 ierr = 3
RETURN
30 ierr = 4
RETURN
END SUBROUTINE depi


SUBROUTINE rfval(x, y, z, rf, ierr)
!-----------------------------------------------------------------------

!          THIS SUBROUTINE COMPUTES THE INCOMPLETE ELLIPTIC INTEGRAL
!          OF THE FIRST KIND

!          RF(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF

!                                -1/2     -1/2     -1/2
!                      (1/2)(T+X)    (T+Y)    (T+Z)    DT,

!          WHERE X, Y, AND Z ARE NONNEGATIVE AND AT MOST ONE OF THEM
!          IS ZERO.  IF ONE OF THEM IS ZERO, THE INTEGRAL IS COMPLETE.
!          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
!          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
!          SERIES TO FIFTH ORDER.  REFERENCE. B. C. CARLSON, COMPUTING
!          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
!          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
!          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
!          MARCH 1, 1980. MODIFIED BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: x, y, z
REAL, INTENT(OUT)    :: rf
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL :: c1, c2, c3, e2, e3, epslon, errtol, lamda, lolim, mu, s, uplim, xn, &
        xndev, xnroot, yn, yndev, ynroot, zn, zndev, znroot
!-----------------------------------------------------------------------

!          INPUT ...

!          X, Y, AND Z ARE THE VARIABLES IN THE INTEGRAL RF(X,Y,Z).

!          OUTPUT ...

!          RF IS THE VALUE OF THE INCOMPLETE ELLIPTIC INTEGRAL.

!          IERR IS THE RETURN ERROR CODE.
!               IERR = 0  FOR NORMAL COMPLETION OF THE SUBROUTINE.
!               IERR = 1  X, Y, OR Z IS NEGATIVE.
!               IERR = 2  X+Y, X+Z, OR Y+Z IS TOO SMALL.
!               IERR = 3  X, Y, OR Z IS TOO LARGE.

!-----------------------------------------------------------------------

!          MACHINE DEPENDENT PARAMETERS ...

!          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!          LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
!          UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.

lolim = 5.0 * TINY(1.0)
uplim = 0.2 * HUGE(1.0)

!          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
!          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
!          ERRTOL ** 6 / (4 * (1 - ERRTOL)).

errtol = (3.6*EPSILON(1.0)) ** (1.0/6.0)

!-----------------------------------------------------------------------
!          WARNING. CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
!          EXPENSE OF ROBUSTNESS.
!-----------------------------------------------------------------------

IF (MIN(x,y,z) >= 0.0) THEN
  IF (MIN(x+y,x+z,y+z) < lolim) GO TO 20
  IF (MAX(x,y,z) > uplim) GO TO 30

  ierr = 0
  xn = x
  yn = y
  zn = z

  10 mu = (xn+yn+zn) / 3.0
  xndev = 2.0 - (mu+xn) / mu
  yndev = 2.0 - (mu+yn) / mu
  zndev = 2.0 - (mu+zn) / mu
  epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev))
  IF (epslon >= errtol) THEN
    xnroot = SQRT(xn)
    ynroot = SQRT(yn)
    znroot = SQRT(zn)
    lamda = xnroot * (ynroot+znroot) + ynroot * znroot
    xn = (xn+lamda) * 0.25
    yn = (yn+lamda) * 0.25
    zn = (zn+lamda) * 0.25
    GO TO 10
  END IF

  c1 = 1.0 / 24.0
  c2 = 3.0 / 44.0
  c3 = 1.0 / 14.0
  e2 = xndev * yndev - zndev * zndev
  e3 = xndev * yndev * zndev
  s = 1.0 + (c1*e2-0.1-c2*e3) * e2 + c3 * e3
  rf = s / SQRT(mu)
  RETURN
END IF

!                      ERROR RETURN

rf = 0.0
ierr = 1
RETURN
20 rf = 0.0
ierr = 2
RETURN
30 rf = 0.0
ierr = 3
RETURN
END SUBROUTINE rfval


SUBROUTINE drfval(x, y, z, rf, ierr)
!-----------------------------------------------------------------------

!          THIS SUBROUTINE COMPUTES THE INCOMPLETE ELLIPTIC INTEGRAL
!          OF THE FIRST KIND

!          RF(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF

!                                -1/2     -1/2     -1/2
!                      (1/2)(T+X)    (T+Y)    (T+Z)    DT,

!          WHERE X, Y, AND Z ARE NONNEGATIVE AND AT MOST ONE OF THEM
!          IS ZERO.  IF ONE OF THEM IS ZERO, THE INTEGRAL IS COMPLETE.
!          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
!          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
!          SERIES TO FIFTH ORDER.  REFERENCE. B. C. CARLSON, COMPUTING
!          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
!          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
!          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
!          MARCH 1, 1980. MODIFIED BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x, y, z
REAL (dp), INTENT(OUT) :: rf
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL (dp) :: c1, c2, c3, e2, e3, epslon, errtol, lamda, lolim, mu, s,  &
             uplim, xn, xndev, xnroot, yn, yndev, ynroot, zn, zndev, znroot
!-----------------------------------------------------------------------

!          INPUT ...

!          X, Y, AND Z ARE THE VARIABLES IN THE INTEGRAL RF(X,Y,Z).

!          OUTPUT ...

!          RF IS THE VALUE OF THE INCOMPLETE ELLIPTIC INTEGRAL.

!          IERR IS THE RETURN ERROR CODE.
!               IERR = 0  FOR NORMAL COMPLETION OF THE SUBROUTINE.
!               IERR = 1  X, Y, OR Z IS NEGATIVE.
!               IERR = 2  X+Y, X+Z, OR Y+Z IS TOO SMALL.
!               IERR = 3  X, Y, OR Z IS TOO LARGE.

!-----------------------------------------------------------------------

!          MACHINE DEPENDENT PARAMETERS ...

!          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!          LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
!          UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.

lolim = 5.0D0 * TINY(1.0D0)
uplim = 0.2D0 * HUGE(1.0D0)

!          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
!          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
!          ERRTOL ** 6 / (4 * (1 - ERRTOL)).

errtol = (3.6 * REAL(EPSILON(1.0D0))) ** (1.0/6.0)

!-----------------------------------------------------------------------
!          WARNING. CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
!          EXPENSE OF ROBUSTNESS.
!-----------------------------------------------------------------------

IF (MIN(x,y,z) >= 0.d0) THEN
  IF (MIN(x+y,x+z,y+z) < lolim) GO TO 20
  IF (MAX(x,y,z) > uplim) GO TO 30

  ierr = 0
  xn = x
  yn = y
  zn = z

  10 mu = (xn+yn+zn) / 3.d0
  xndev = 2.d0 - (mu+xn) / mu
  yndev = 2.d0 - (mu+yn) / mu
  zndev = 2.d0 - (mu+zn) / mu
  epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev))
  IF (epslon >= errtol) THEN
    xnroot = SQRT(xn)
    ynroot = SQRT(yn)
    znroot = SQRT(zn)
    lamda = xnroot * (ynroot+znroot) + ynroot * znroot
    xn = (xn+lamda) * 0.25D0
    yn = (yn+lamda) * 0.25D0
    zn = (zn+lamda) * 0.25D0
    GO TO 10
  END IF

  c1 = 1.d0 / 24.d0
  c2 = 3.d0 / 44.d0
  c3 = 1.d0 / 14.d0
  e2 = xndev * yndev - zndev * zndev
  e3 = xndev * yndev * zndev
  s = 1.d0 + (c1*e2-0.1D0-c2*e3) * e2 + c3 * e3
  rf = s / SQRT(mu)
  RETURN
END IF

!                      ERROR RETURN

rf = 0.d0
ierr = 1
RETURN
20 rf = 0.d0
ierr = 2
RETURN
30 rf = 0.d0
ierr = 3
RETURN
END SUBROUTINE drfval


SUBROUTINE rdval(x, y, z, rd, ierr)
!-----------------------------------------------------------------------

!          THIS SUBROUTINE COMPUTES THE INCOMPLETE ELLIPTIC INTEGRAL
!          OF THE SECOND KIND

!          RD(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF

!                                -1/2     -1/2     -3/2
!                      (3/2)(T+X)    (T+Y)    (T+Z)    DT,

!          WHERE X AND Y ARE NONNEGATIVE, X + Y IS POSITIVE, AND Z IS
!          POSITIVE.  IF X OR Y IS ZERO, THE INTEGRAL IS COMPLETE.
!          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
!          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
!          SERIES TO FIFTH ORDER.  REFERENCE. B. C. CARLSON, COMPUTING
!          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
!          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
!          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
!          MARCH 1, 1980. MODIFIED BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: x, y, z
REAL, INTENT(OUT)    :: rd
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL    :: c1, c2, c3, c4, ea, eb, ec, ed, ef, epslon, errtol, lamda
REAL    :: lolim, mu, power4, sigma, s1, s2, uplim, xn, xndev
REAL    :: xnroot, yn, yndev, ynroot, zn, zndev, znroot
!-----------------------------------------------------------------------

!          INPUT ...

!          X, Y, AND Z ARE THE VARIABLES IN THE INTEGRAL RD(X,Y,Z).

!          OUTPUT ...

!          RD IS THE VALUE OF THE INCOMPLETE ELLIPTIC INTEGRAL.

!          IERR IS THE RETURN ERROR CODE.
!               IERR = 0  FOR NORMAL COMPLETION OF THE SUBROUTINE.
!               IERR = 1  X, Y, OR Z IS NEGATIVE.
!               IERR = 2  X+Y OR Z IS TOO SMALL.
!               IERR = 3  X, Y, OR Z IS TOO LARGE.

!-----------------------------------------------------------------------

!          MACHINE DEPENDENT PARAMETERS ...

!          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
!          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
!          3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2.

errtol = (.28*EPSILON(1.0)) ** (1.0/6.0)

!          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!          LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ** (2/3).
!          UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
!          MINIMUM) ** (2/3).

mu = -2.0 / 3.0
lolim = 2.0001 * HUGE(1.0) ** mu
uplim = (10.0*TINY(1.0)/errtol) ** mu

!-----------------------------------------------------------------------
!          WARNING. CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
!          EXPENSE OF ROBUSTNESS.
!-----------------------------------------------------------------------

IF (MIN(x,y,z) >= 0.0) THEN
  IF (MIN(x+y,z) < lolim) GO TO 20
  IF (MAX(x,y,z) > uplim) GO TO 30

  ierr = 0
  xn = x
  yn = y
  zn = z
  sigma = 0.0
  power4 = 1.0

  10 mu = (xn+yn+3.0*zn) * 0.2
  xndev = (mu-xn) / mu
  yndev = (mu-yn) / mu
  zndev = (mu-zn) / mu
  epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev))
  IF (epslon >= errtol) THEN
    xnroot = SQRT(xn)
    ynroot = SQRT(yn)
    znroot = SQRT(zn)
    lamda = xnroot * (ynroot+znroot) + ynroot * znroot
    sigma = sigma + power4 / (znroot*(zn+lamda))
    power4 = power4 * 0.25
    xn = (xn+lamda) * 0.25
    yn = (yn+lamda) * 0.25
    zn = (zn+lamda) * 0.25
    GO TO 10
  END IF

  c1 = 3.0 / 14.0
  c2 = 1.0 / 6.0
  c3 = 9.0 / 22.0
  c4 = 3.0 / 26.0
  ea = xndev * yndev
  eb = zndev * zndev
  ec = ea - eb
  ed = ea - 6.0 * eb
  ef = ed + ec + ec
  s1 = ed * (-c1+0.25*c3*ed-1.5*c4*zndev*ef)
  s2 = zndev * (c2*ef+zndev*(-c3*ec+zndev*c4*ea))
  rd = 3.0 * sigma + power4 * (1.0+s1+s2) / (mu*SQRT(mu))
  RETURN
END IF

!                        ERROR RETURN

rd = 0.0
ierr = 1
RETURN
20 rd = 0.0
ierr = 2
RETURN
30 rd = 0.0
ierr = 3
RETURN
END SUBROUTINE rdval


SUBROUTINE drdval(x, y, z, rd, ierr)
!-----------------------------------------------------------------------

!          THIS SUBROUTINE COMPUTES THE INCOMPLETE ELLIPTIC INTEGRAL
!          OF THE SECOND KIND

!          RD(X,Y,Z) = INTEGRAL FROM ZERO TO INFINITY OF

!                                -1/2     -1/2     -3/2
!                      (3/2)(T+X)    (T+Y)    (T+Z)    DT,

!          WHERE X AND Y ARE NONNEGATIVE, X + Y IS POSITIVE, AND Z IS
!          POSITIVE.  IF X OR Y IS ZERO, THE INTEGRAL IS COMPLETE.
!          THE DUPLICATION THEOREM IS ITERATED UNTIL THE VARIABLES ARE
!          NEARLY EQUAL, AND THE FUNCTION IS THEN EXPANDED IN TAYLOR
!          SERIES TO FIFTH ORDER.  REFERENCE. B. C. CARLSON, COMPUTING
!          ELLIPTIC INTEGRALS BY DUPLICATION, NUMER. MATH. 33 (1979),
!          1-16.  CODED BY B. C. CARLSON AND ELAINE M. NOTIS, AMES
!          LABORATORY-DOE, IOWA STATE UNIVERSITY, AMES, IOWA 50011.
!          MARCH 1, 1980. MODIFIED BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x, y, z
REAL (dp), INTENT(OUT) :: rd
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL (dp) :: c1, c2, c3, c4, ea, eb, ec, ed, ef, epslon, errtol, lamda, &
             lolim, mu, power4, sigma, s1, s2, uplim, xn, xndev, xnroot, &
             yn, yndev, ynroot, zn, zndev, znroot
!-----------------------------------------------------------------------

!          INPUT ...

!          X, Y, AND Z ARE THE VARIABLES IN THE INTEGRAL RD(X,Y,Z).

!          OUTPUT ...

!          RD IS THE VALUE OF THE INCOMPLETE ELLIPTIC INTEGRAL.

!          IERR IS THE RETURN ERROR CODE.
!               IERR = 0  FOR NORMAL COMPLETION OF THE SUBROUTINE.
!               IERR = 1  X, Y, OR Z IS NEGATIVE.
!               IERR = 2  X+Y OR Z IS TOO SMALL.
!               IERR = 3  X, Y, OR Z IS TOO LARGE.

!-----------------------------------------------------------------------

!          MACHINE DEPENDENT PARAMETERS ...

!          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
!          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
!          3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2.

errtol = (.28 * REAL(EPSILON(1.0D0))) ** (1.0/6.0)

!          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!          LOLIM IS NOT LESS THAN 2 / (MACHINE MAXIMUM) ** (2/3).
!          UPLIM IS NOT GREATER THAN (0.1 * ERRTOL / MACHINE
!          MINIMUM) ** (2/3).

mu = -2.d0 / 3.d0
lolim = 2.00000000001D0 * HUGE(1.0D0) ** mu
uplim = (10.d0*TINY(1.0D0)/errtol) ** mu

!-----------------------------------------------------------------------
!          WARNING. CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
!          EXPENSE OF ROBUSTNESS.
!-----------------------------------------------------------------------

IF (MIN(x,y,z) >= 0.d0) THEN
  IF (MIN(x+y,z) < lolim) GO TO 20
  IF (MAX(x,y,z) > uplim) GO TO 30

  ierr = 0
  xn = x
  yn = y
  zn = z
  sigma = 0.d0
  power4 = 1.d0

  10 mu = (xn+yn+3.d0*zn) * 0.2D0
  xndev = (mu-xn) / mu
  yndev = (mu-yn) / mu
  zndev = (mu-zn) / mu
  epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev))
  IF (epslon >= errtol) THEN
    xnroot = SQRT(xn)
    ynroot = SQRT(yn)
    znroot = SQRT(zn)
    lamda = xnroot * (ynroot+znroot) + ynroot * znroot
    sigma = sigma + power4 / (znroot*(zn+lamda))
    power4 = power4 * 0.25D0
    xn = (xn+lamda) * 0.25D0
    yn = (yn+lamda) * 0.25D0
    zn = (zn+lamda) * 0.25D0
    GO TO 10
  END IF

  c1 = 3.d0 / 14.d0
  c2 = 1.d0 / 6.d0
  c3 = 9.d0 / 22.d0
  c4 = 3.d0 / 26.d0
  ea = xndev * yndev
  eb = zndev * zndev
  ec = ea - eb
  ed = ea - 6.d0 * eb
  ef = ed + ec + ec
  s1 = ed * (-c1+0.25D0*c3*ed-1.5D0*c4*zndev*ef)
  s2 = zndev * (c2*ef+zndev*(-c3*ec+zndev*c4*ea))
  rd = 3.d0 * sigma + power4 * (1.d0+s1+s2) / (mu*SQRT(mu))
  RETURN
END IF

!                        ERROR RETURN

rd = 0.d0
ierr = 1
RETURN
20 rd = 0.d0
ierr = 2
RETURN
30 rd = 0.d0
ierr = 3
RETURN
END SUBROUTINE drdval


SUBROUTINE rjval(x, y, z, p, rj, ierr)
!-----------------------------------------------------------------------

!          THIS SUBROUTINE COMPUTES THE INCOMPLETE ELLIPTIC INTEGRAL
!          OF THE THIRD KIND

!          RJ(X,Y,Z,P) = INTEGRAL FROM ZERO TO INFINITY OF

!                                  -1/2     -1/2     -1/2     -1
!                        (3/2)(T+X)    (T+Y)    (T+Z)    (t + p)  DT,

!          WHERE X, Y, AND Z ARE NONNEGATIVE, AT MOST ONE OF THEM IS
!          ZERO, AND P IS POSITIVE.  IF X OR Y OR Z IS ZERO, THE
!          INTEGRAL IS COMPLETE.  THE DUPLICATION THEOREM IS ITERATED
!          UNTIL THE VARIABLES ARE NEARLY EQUAL, AND THE FUNCTION IS
!          THEN EXPANDED IN TAYLOR SERIES TO FIFTH ORDER.  REFERENCE.
!          B. C. CARLSON, COMPUTING ELLIPTIC INTEGRALS BY DUPLICATION,
!          NUMER. MATH. 33 (1979), 1-16.  CODED BY B. C. CARLSON AND
!          ELAINE M. NOTIS, AMES LABORATORY-DOE, IOWA STATE UNIVERSITY,
!          AMES, IOWA 50011.  MARCH 1, 1980. MODIFIED BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: x, y, z, p
REAL, INTENT(OUT)    :: rj
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL :: rc, alfa, beta, c1, c2, c3, c4, ea, eb, ec, e2, e3, epslon,  &
        errtol, etolrc, lamda, lolim, mu, pn, pndev, power4, sigma, s1, &
        s2, s3, uplim, xn, xndev, xnroot, yn, yndev, ynroot, zn, zndev, &
        znroot
!-----------------------------------------------------------------------

!          INPUT ...

!          X, Y, Z, AND P ARE THE VARIABLES IN THE INTEGRAL RJ(X,Y,Z,P).

!          OUTPUT ...

!          RJ IS THE VALUE OF THE INCOMPLETE ELLIPTIC INTEGRAL.

!          IERR IS THE RETURN ERROR CODE.
!               IERR = 0  FOR NORMAL COMPLETION OF THE SUBROUTINE.
!               IERR = 1  X, Y, Z, OR P IS NEGATIVE.
!               IERR = 2  X+Y, X+Z, Y+Z, OR P IS TOO SMALL.
!               IERR = 3  X, Y, Z, OR P IS TOO LARGE.

!-----------------------------------------------------------------------

!          MACHINE DEPENDENT PARAMETERS ...

!          RC IS A FUNCTION COMPUTED BY THE SUBROUTINE RCVAL1.
!          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!          LOLIM IS NOT LESS THAN THE CUBE ROOT OF THE VALUE
!          OF LOLIM USED IN THE CODE FOR RC, AND
!          UPLIM IS NOT GREATER THAN 0.3 TIMES THE CUBE ROOT OF
!          THE VALUE OF UPLIM USED IN THE CODE FOR RC.

mu = 1.0 / 3.0
lolim = 1.0001 * (5.0*TINY(1.0)) ** mu
uplim = .29999 * (0.2*HUGE(1.0)) ** mu

!          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE. THE
!          RELATIVE ERROR DUE TO TRUNCATION OF THE SERIES FOR RJ
!          IS LESS THAN 3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2.
!          AN ERROR TOLERANCE (ETOLRC) WILL BE PASSED TO THE CODE FOR
!          RC TO MAKE THE TRUNCATION ERROR FOR RC LESS THAN FOR RJ.

errtol = (.28*EPSILON(1.0)) ** (1.0/6.0)

!-----------------------------------------------------------------------
!          WARNING. CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
!          EXPENSE OF ROBUSTNESS.
!-----------------------------------------------------------------------

IF (MIN(x,y,z,p) >= 0.0) THEN
  IF (MIN(x+y,x+z,y+z,p) < lolim) GO TO 20
  IF (MAX(x,y,z,p) > uplim) GO TO 30

  ierr = 0
  xn = x
  yn = y
  zn = z
  pn = p
  sigma = 0.0
  power4 = 1.0
  etolrc = 0.5 * errtol

  10 mu = (xn+yn+zn+pn+pn) * 0.2
  xndev = (mu-xn) / mu
  yndev = (mu-yn) / mu
  zndev = (mu-zn) / mu
  pndev = (mu-pn) / mu
  epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev),ABS(pndev))
  IF (epslon >= errtol) THEN
    xnroot = SQRT(xn)
    ynroot = SQRT(yn)
    znroot = SQRT(zn)
    lamda = xnroot * (ynroot+znroot) + ynroot * znroot
    alfa = pn * (xnroot+ynroot+znroot) + xnroot * ynroot * znroot
    alfa = alfa * alfa
    beta = pn * (pn+lamda) * (pn+lamda)
    CALL rcval1(alfa,beta,etolrc,rc,ierr)
    IF (ierr /= 0) RETURN
    sigma = sigma + power4 * rc
    power4 = power4 * 0.25
    xn = (xn+lamda) * 0.25
    yn = (yn+lamda) * 0.25
    zn = (zn+lamda) * 0.25
    pn = (pn+lamda) * 0.25
    GO TO 10
  END IF

  c1 = 3.0 / 14.0
  c2 = 1.0 / 3.0
  c3 = 3.0 / 22.0
  c4 = 3.0 / 26.0
  ea = xndev * (yndev+zndev) + yndev * zndev
  eb = xndev * yndev * zndev
  ec = pndev * pndev
  e2 = ea - 3.0 * ec
  e3 = eb + 2.0 * pndev * (ea-ec)
  s1 = 1.0 + e2 * (-c1+0.75*c3*e2-1.5*c4*e3)
  s2 = eb * (0.5*c2+pndev*(-c3-c3+pndev*c4))
  s3 = pndev * ea * (c2-pndev*c3) - c2 * pndev * ec
  rj = 3.0 * sigma + power4 * (s1+s2+s3) / (mu*SQRT(mu))
  RETURN
END IF

!                        ERROR RETURN

rj = 0.0
ierr = 1
RETURN
20 rj = 0.0
ierr = 2
RETURN
30 rj = 0.0
ierr = 3
RETURN
END SUBROUTINE rjval


SUBROUTINE drjval(x, y, z, p, rj, ierr)
!-----------------------------------------------------------------------

!          THIS SUBROUTINE COMPUTES THE INCOMPLETE ELLIPTIC INTEGRAL
!          OF THE THIRD KIND

!          RJ(X,Y,Z,P) = INTEGRAL FROM ZERO TO INFINITY OF

!                                  -1/2     -1/2     -1/2     -1
!                        (3/2)(T+X)    (T+Y)    (T+Z)    (t + p)  DT,

!          WHERE X, Y, AND Z ARE NONNEGATIVE, AT MOST ONE OF THEM IS
!          ZERO, AND P IS POSITIVE.  IF X OR Y OR Z IS ZERO, THE
!          INTEGRAL IS COMPLETE.  THE DUPLICATION THEOREM IS ITERATED
!          UNTIL THE VARIABLES ARE NEARLY EQUAL, AND THE FUNCTION IS
!          THEN EXPANDED IN TAYLOR SERIES TO FIFTH ORDER.  REFERENCE.
!          B. C. CARLSON, COMPUTING ELLIPTIC INTEGRALS BY DUPLICATION,
!          NUMER. MATH. 33 (1979), 1-16.  CODED BY B. C. CARLSON AND
!          ELAINE M. NOTIS, AMES LABORATORY-DOE, IOWA STATE UNIVERSITY,
!          AMES, IOWA 50011.  MARCH 1, 1980. MODIFIED BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x, y, z, p
REAL (dp), INTENT(OUT) :: rj
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL (dp) :: rc, alfa, beta, c1, c2, c3, c4, ea, eb, ec, e2, e3, epslon, &
             errtol, etolrc, lamda, lolim, mu, pn, pndev, power4, sigma, &
             s1, s2, s3, uplim, xn, xndev, xnroot, yn, yndev, ynroot,    &
             zn, zndev, znroot
!-----------------------------------------------------------------------

!          INPUT ...

!          X, Y, Z, AND P ARE THE VARIABLES IN THE INTEGRAL RJ(X,Y,Z,P).

!          OUTPUT ...

!          RJ IS THE VALUE OF THE INCOMPLETE ELLIPTIC INTEGRAL.

!          IERR IS THE RETURN ERROR CODE.
!               IERR = 0  FOR NORMAL COMPLETION OF THE SUBROUTINE.
!               IERR = 1  X, Y, Z, OR P IS NEGATIVE.
!               IERR = 2  X+Y, X+Z, Y+Z, OR P IS TOO SMALL.
!               IERR = 3  X, Y, Z, OR P IS TOO LARGE.

!-----------------------------------------------------------------------

!          MACHINE DEPENDENT PARAMETERS ...

!          RC IS A FUNCTION COMPUTED BY THE SUBROUTINE DRCVL1.
!          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!          LOLIM IS NOT LESS THAN THE CUBE ROOT OF THE VALUE
!          OF LOLIM USED IN THE CODE FOR RC, AND
!          UPLIM IS NOT GREATER THAN 0.3 TIMES THE CUBE ROOT OF
!          THE VALUE OF UPLIM USED IN THE CODE FOR RC.

mu = 1._dp / 3._dp
lolim = 1.00000000001_dp * (5.0_dp*TINY(1.0_dp)) ** mu
uplim = .299999999999_dp * (0.2_dp*HUGE(1.0_dp)) ** mu

!          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
!          RELATIVE ERROR DUE TO TRUNCATION OF THE SERIES FOR RJ
!          IS LESS THAN 3 * ERRTOL ** 6 / (1 - ERRTOL) ** 3/2.
!          AN ERROR TOLERANCE (ETOLRC) WILL BE PASSED TO THE CODE FOR
!          RC TO MAKE THE TRUNCATION ERROR FOR RC LESS THAN FOR RJ.

errtol = (.28 * REAL(EPSILON(1.0_dp))) ** (1.0/6.0)

!-----------------------------------------------------------------------
!          WARNING. CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
!          EXPENSE OF ROBUSTNESS.
!-----------------------------------------------------------------------

IF (MIN(x,y,z,p) >= 0._dp) THEN
  IF (MIN(x+y,x+z,y+z,p) < lolim) GO TO 20
  IF (MAX(x,y,z,p) > uplim) GO TO 30

  ierr = 0
  xn = x
  yn = y
  zn = z
  pn = p
  sigma = 0._dp
  power4 = 1._dp
  etolrc = 0.5_dp * errtol

  10 mu = (xn+yn+zn+pn+pn) * 0.2_dp
  xndev = (mu-xn) / mu
  yndev = (mu-yn) / mu
  zndev = (mu-zn) / mu
  pndev = (mu-pn) / mu
  epslon = MAX(ABS(xndev),ABS(yndev),ABS(zndev),ABS(pndev))
  IF (epslon >= errtol) THEN
    xnroot = SQRT(xn)
    ynroot = SQRT(yn)
    znroot = SQRT(zn)
    lamda = xnroot * (ynroot+znroot) + ynroot * znroot
    alfa = pn * (xnroot+ynroot+znroot) + xnroot * ynroot * znroot
    alfa = alfa * alfa
    beta = pn * (pn+lamda) * (pn+lamda)
    CALL drcvl1(alfa,beta,etolrc,rc,ierr)
    IF (ierr /= 0) RETURN
    sigma = sigma + power4 * rc
    power4 = power4 * 0.25_dp
    xn = (xn+lamda) * 0.25_dp
    yn = (yn+lamda) * 0.25_dp
    zn = (zn+lamda) * 0.25_dp
    pn = (pn+lamda) * 0.25_dp
    GO TO 10
  END IF

  c1 = 3._dp / 14._dp
  c2 = 1._dp / 3._dp
  c3 = 3._dp / 22._dp
  c4 = 3._dp / 26._dp
  ea = xndev * (yndev+zndev) + yndev * zndev
  eb = xndev * yndev * zndev
  ec = pndev * pndev
  e2 = ea - 3._dp * ec
  e3 = eb + 2._dp * pndev * (ea-ec)
  s1 = 1._dp + e2 * (-c1+0.75_dp*c3*e2-1.5_dp*c4*e3)
  s2 = eb * (0.5_dp*c2+pndev*(-c3-c3+pndev*c4))
  s3 = pndev * ea * (c2-pndev*c3) - c2 * pndev * ec
  rj = 3._dp * sigma + power4 * (s1+s2+s3) / (mu*SQRT(mu))
  RETURN
END IF

!                        ERROR RETURN

rj = 0._dp
ierr = 1
RETURN
20 rj = 0._dp
ierr = 2
RETURN
30 rj = 0._dp
ierr = 3
RETURN
END SUBROUTINE drjval


SUBROUTINE rcval1(x, y, errtol, rc, ierr)
!-----------------------------------------------------------------------

!          THIS SUBROUTINE COMPUTES THE INTEGRAL

!          RC(X,Y) = INTEGRAL FROM ZERO TO INFINITY OF

!                              -1/2     -1
!                    (1/2)(T+X)    (T+Y)  DT,

!          WHERE X IS NONNEGATIVE AND Y IS POSITIVE.  THE DUPLICATION
!          THEOREM IS ITERATED UNTIL THE VARIABLES ARE NEARLY EQUAL,
!          AND THE FUNCTION IS THEN EXPANDED IN TAYLOR SERIES TO FIFTH
!          ORDER.  LOGARITHMIC, INVERSE CIRCULAR, AND INVERSE HYPERBOLIC
!          FUNCTIONS CAN BE EXPRESSED IN TERMS OF RC.   REFERENCE.
!          B. C. CARLSON, COMPUTING ELLIPTIC INTEGRALS BY DUPLICATION,
!          NUMER. MATH. 33 (1979), 1-16.  CODED BY B. C. CARLSON AND
!          ELAINE M. NOTIS, AMES LABORATORY-DOE, IOWA STATE UNIVERSITY,
!          AMES, IOWA 50011.  MARCH 1, 1980. MODIFIED BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------
REAL, INTENT(IN)     :: x, y, errtol
REAL, INTENT(OUT)    :: rc
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL :: c1, c2, lamda, lolim, mu, s, sn, uplim, xn, yn
!-----------------------------------------------------------------------

!          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!          LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
!          UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.

lolim = 5.0 * TINY(1.0)
uplim = 0.2 * HUGE(1.0)

!          INPUT ...

!          X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).

!          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
!          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
!          16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).

!          SAMPLE CHOICES   ERRTOL   RELATIVE TRUNCATION
!                                    ERROR LESS THAN
!                           1.E-3    2.E-17
!                           3.E-3    2.E-14
!                           1.E-2    2.E-11
!                           3.E-2    2.E-8
!                           1.E-1    2.E-5

!          OUTPUT ...

!          RC IS THE VALUE OF THE INTEGRAL.

!          IERR IS THE RETURN ERROR CODE.
!               IERR = 0  FOR NORMAL COMPLETION OF THE SUBROUTINE.
!               IERR = 1  X OR Y IS NEGATIVE, OR Y = 0.
!               IERR = 2  X+Y IS TOO SMALL.
!               IERR = 3  X OR Y IS TOO LARGE.

!-----------------------------------------------------------------------
!          WARNING. CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
!          EXPENSE OF ROBUSTNESS.
!-----------------------------------------------------------------------

IF (x >= 0.0 .AND. y > 0.0) THEN
  IF ((x+y) < lolim) GO TO 20
  IF (MAX(x,y) > uplim) GO TO 30

  ierr = 0
  xn = x
  yn = y

  10 mu = (xn+yn+yn) / 3.0
  sn = (yn+mu) / mu - 2.0
  IF (ABS(sn) >= errtol) THEN
    lamda = 2.0 * SQRT(xn) * SQRT(yn) + yn
    xn = (xn+lamda) * 0.25
    yn = (yn+lamda) * 0.25
    GO TO 10
  END IF

  c1 = 1.0 / 7.0
  c2 = 9.0 / 22.0
  s = sn * sn * (0.3+sn*(c1+sn*(0.375+sn*c2)))
  rc = (1.0+s) / SQRT(mu)
  RETURN
END IF

!                      ERROR RETURN

rc = 0.0
ierr = 1
RETURN
20 rc = 0.0
ierr = 2
RETURN
30 rc = 0.0
ierr = 3
RETURN
END SUBROUTINE rcval1


SUBROUTINE drcvl1(x, y, errtol, rc, ierr)
!-----------------------------------------------------------------------

!          THIS SUBROUTINE COMPUTES THE INTEGRAL

!          RC(X,Y) = INTEGRAL FROM ZERO TO INFINITY OF

!                              -1/2     -1
!                    (1/2)(T+X)    (T+Y)  DT,

!          WHERE X IS NONNEGATIVE AND Y IS POSITIVE.  THE DUPLICATION THEOREM
!          IS ITERATED UNTIL THE VARIABLES ARE NEARLY EQUAL, AND THE FUNCTION
!          IS THEN EXPANDED IN TAYLOR SERIES TO FIFTH ORDER.
!          LOGARITHMIC, INVERSE CIRCULAR, AND INVERSE HYPERBOLIC FUNCTIONS
!          CAN BE EXPRESSED IN TERMS OF RC.   REFERENCE.
!          B. C. CARLSON, COMPUTING ELLIPTIC INTEGRALS BY DUPLICATION,
!          NUMER. MATH. 33 (1979), 1-16.  CODED BY B. C. CARLSON AND
!          ELAINE M. NOTIS, AMES LABORATORY-DOE, IOWA STATE UNIVERSITY,
!          AMES, IOWA 50011.  MARCH 1, 1980. MODIFIED BY A.H. MORRIS (NSWC).

!-----------------------------------------------------------------------
REAL (dp), INTENT(IN)  :: x, y, errtol
REAL (dp), INTENT(OUT) :: rc
INTEGER, INTENT(OUT)   :: ierr

! Local variables
REAL (dp) :: c1, c2, lamda, lolim, mu, s, sn, uplim, xn, yn
!-----------------------------------------------------------------------

!          LOLIM AND UPLIM DETERMINE THE RANGE OF VALID ARGUMENTS.
!          LOLIM IS NOT LESS THAN THE MACHINE MINIMUM MULTIPLIED BY 5.
!          UPLIM IS NOT GREATER THAN THE MACHINE MAXIMUM DIVIDED BY 5.

lolim = 5.0_dp * TINY(1.0_dp)
uplim = 0.2_dp * HUGE(1.0_dp)

!          INPUT ...

!          X AND Y ARE THE VARIABLES IN THE INTEGRAL RC(X,Y).

!          ERRTOL IS SET TO THE DESIRED ERROR TOLERANCE.
!          RELATIVE ERROR DUE TO TRUNCATION IS LESS THAN
!          16 * ERRTOL ** 6 / (1 - 2 * ERRTOL).

!          OUTPUT ...

!          RC IS THE VALUE OF THE INTEGRAL.

!          IERR IS THE RETURN ERROR CODE.
!               IERR = 0  FOR NORMAL COMPLETION OF THE SUBROUTINE.
!               IERR = 1  X OR Y IS NEGATIVE, OR Y = 0.
!               IERR = 2  X+Y IS TOO SMALL.
!               IERR = 3  X OR Y IS TOO LARGE.

!-----------------------------------------------------------------------
!          WARNING. CHANGES IN THE PROGRAM MAY IMPROVE SPEED AT THE
!          EXPENSE OF ROBUSTNESS.
!-----------------------------------------------------------------------

IF (x >= 0._dp .AND. y > 0._dp) THEN
  IF ((x+y) < lolim) GO TO 20
  IF (MAX(x,y) > uplim) GO TO 30

  ierr = 0
  xn = x
  yn = y

  10 mu = (xn+yn+yn) / 3._dp
  sn = (yn+mu) / mu - 2._dp
  IF (ABS(sn) >= errtol) THEN
    lamda = 2._dp * SQRT(xn) * SQRT(yn) + yn
    xn = (xn+lamda) * 0.25_dp
    yn = (yn+lamda) * 0.25_dp
    GO TO 10
  END IF

  c1 = 1._dp / 7._dp
  c2 = 9._dp / 22._dp
  s = sn * sn * (0.3_dp+sn*(c1+sn*(0.375_dp+sn*c2)))
  rc = (1._dp+s) / SQRT(mu)
  RETURN
END IF

!                      ERROR RETURN

rc = 0._dp
ierr = 1
RETURN
20 rc = 0._dp
ierr = 2
RETURN
30 rc = 0._dp
ierr = 3
RETURN
END SUBROUTINE drcvl1


SUBROUTINE ellpf(u, k, l, s, c, d, ierr)
!     -------------------------------------------------------------
!     ELLPF CALCULATES THE ELLIPTIC FUNCTIONS SN(U,K), CN(U,K), AND
!     DN(U,K) FOR REAL U AND REAL MODULUS K.  IT IS ASSUMED THAT
!     ABS(K) <= 1. AND K**2 + L**2 = 1.
!     -------------------------------------------------------------
REAL, INTENT(IN)     :: u, k, l
REAL, INTENT(OUT)    :: s, c, d
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL, PARAMETER :: pihalf = 1.5707963267949
INTEGER :: max, n
REAL    :: e, eps, f, f2, r, sg, tol, u1, u2, z
!     ----------------
!     ****** MAX AND EPS ARE MACHINE DEPENDENT CONSTANTS.
!            MAX IS THE LARGEST POSITIVE INTEGER THAT MAY
!            BE USED, AND EPS IS THE SMALLEST REAL NUMBER
!            FOR WHICH 1 + EPS > 1.

!                    MAX = IPMPAR(3)
MAX = huge(3)
eps = EPSILON(1.0)

!     ----------------
!     CALCULATION FOR L = 0.0

IF (l == 0.0) THEN
  s = TANH(u)
  e = EXP(-ABS(u))
  c = 2.0 * e / (1.0 + e*e)
  d = c
  ierr = 0
  RETURN
END IF

!     CHECK THAT K**2 + L**2 = 1

tol = 2.0 * eps
z = DBLE(k*k) + (DBLE(l*l)-1._dp)
IF (ABS(z) <= tol) THEN

  f = pihalf
  IF (k /= 0.0) CALL ellpi(pihalf, 0.0, k, l, f, e, ierr)
  f2 = 2.0 * f

!                   ARGUMENT REDUCTION

  u1 = ABS(u)
  r = u1 / f2
  IF (r >= MIN(REAL(MAX), 1.0/eps)) GO TO 20
  n = INT(r)
  u1 = u1 - REAL(n) * f2
  sg = 1.0
  IF (MOD(n,2) /= 0) sg = -1.0

  IF (u1 > 0.0) THEN
    IF (u1 > f) THEN
      u1 = u1 - f2
      sg = -sg
      IF (u1 >= 0.0) GO TO 10
    END IF

!     CALCULATION OF ELLIPTIC FUNCTIONS FOR 0.0 <= U2 <= F(K)

    u2 = ABS(u1)
    CALL scd(u2,ABS(k),ABS(l),f,s,c,d)
    ierr = 0
    IF (u1 < 0.0) s = -s

!     FINAL ASSEMBLY

    s = sg * s
    c = sg * c
    IF (u < 0.0) s = -s
    RETURN
  END IF

!             U IS AN INTEGER MULTIPLE OF F2

  10 s = 0.0
  c = sg
  d = 1.0
  ierr = 0
  RETURN
END IF

!                      ERROR RETURN

ierr = 1
RETURN
20 ierr = 2
RETURN
END SUBROUTINE ellpf


SUBROUTINE scd(u, k, l, f, s, c, d)
!     --------------------------------------------------------
!     SCD COMPUTES THE ELLIPTIC FUNCTIONS SN(U,K), CN(U,K),
!     AND DN(U,K) FOR REAL U AND REAL MODULUS K SUCH THAT
!     0.0 <= U <= F AND 0.0 <= K < 1.0, WHERE
!     F = F(K) IS THE COMPLETE ELLIPTIC INTEGRAL OF THE
!     FIRST KIND, AND F1 = F(L) IS THE COMPLEMENTARY INTEGRAL.
!     IT IS ASSUMED THAT K**2 + L**2 = 1.
!     --------------------------------------------------------
REAL, INTENT(IN)  :: u, k, l, f
REAL, INTENT(OUT) :: s, c, d

! Local variables
REAL, PARAMETER :: pihalf = 1.5707963267949
REAL            :: c1, d1, e1, f1, s1, v
INTEGER         :: ierr
!     ------------------------
IF (k /= 0.0) THEN
  v = f - u

!     USES MACLAURIN EXPANSION WHEN U OR V <= 0.01

  IF (u <= 0.01) THEN
    CALL scdm(u, k, s, c, d)
    RETURN
  END IF
  IF (v <= 0.01) THEN
    CALL scdm(v, k, s1, c1, d1)
    s = c1 / d1
    c = l * s1 / d1
    d = l / d1
    RETURN
  END IF

!     USES FOURIER EXPANSION WHEN K <= L

  CALL ellpi(pihalf, 0.0, l, k, f1, e1, ierr)
  IF (ierr /= 0) THEN
    WRITE(*, *) '** Error in call to ELLPI from SCD **'
    e1 = 1.0 + e1            ! Redundant; to ensure that e1 is used
  END IF

  IF (k <= l) THEN
    CALL scdf(u, k, l, f, f1, s, c, d)
    RETURN
  END IF

!     USES IMAGINARY TRANSFORMATION OF JACOBI & FOURIER EXPANSION WHEN K > L

  CALL scdj(u, k, l, f, f1, s, c, d)
  RETURN
END IF

!     COMPUTATION FOR K = 0.0

s = SIN(u)
c = COS(u)
d = 1.0
RETURN
END SUBROUTINE scd


SUBROUTINE scdm(u, k, s, c, d)
!     -------------------------------------------------
!     CALCULATES SN(U,K), CN(U,K), AND DN(U,K) FOR
!     0.0 <= U <= 0.01 AND FOR 0.0 <= K <= 1.0
!     BY USE OF THE MACLAURIN EXPANSION FOR SN(U,K)
!     -------------------------------------------------
REAL, INTENT(IN)  :: u, k
REAL, INTENT(OUT) :: s, c, d

! Local variables
REAL :: c1, c2, c3, c4, k2, u2

k2 = k * k
u2 = u * u
c1 = -(1.0+k2) / 6.0
c2 = (1.0+k2*(14.0+k2)) / 120.0
c3 = -(1.0+k2*(135.0+k2*(135.0+k2))) / 5040.0
c4 = (1.0+k2*(1228.0+k2*(5478.0+k2*(1228.0+k2)))) / 362880.0
s = u * (1.0+u2*(c1+u2*(c2+u2*(c3+c4*u2))))
c = SQRT(1.0-s*s)
d = SQRT(1.0-(k*s)**2)
RETURN
END SUBROUTINE scdm


SUBROUTINE scdf(u, k, l, f, f1, s, c, d)
!     -------------------------------------------------------------
!     SCDF COMPUTES SN(U,K), CN(U,K), AND DN(U,K) FOR REAL U AND
!     K BY USE OF THE FOURIER EXPANSION FOR SN(U,K).  IT IS
!     ASSUMED THAT 0.0 <= K < 1.0 AND 0.0 <= U <= F,
!     WHERE F = F(K) IS THE COMPLETE ELLIPTIC INTEGRAL OF THE
!     FIRST KIND AND F1 = F(L) IS THE COMPLEMENTARY INTEGRAL, WITH
!     L .NE. 0. AND K**2 + L**2 = 1.
!     -------------------------------------------------------------
REAL, INTENT(IN)  :: u, k, l, f, f1
REAL, INTENT(OUT) :: s, c, d

! Local variables
REAL, PARAMETER :: pihalf = 1.5707963267949
REAL :: a, ai, coef, eps, i, qh, q1, q2, qn, qd, sum, temp, tol, v, w, x
!     -------------------------------------------------
!     EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE
!     SMALLEST NUMBER SUCH THAT 1. + EPS > 1.

eps = EPSILON(1.0)

!     -------------------------------------------------
tol = eps / 10.0
v = f - u
qh = EXP(-pihalf*f1/f)
q1 = qh * qh
q2 = q1 * q1
coef = 4. * pihalf * qh / (k*f)
qn = 1.0
qd = q1
w = MIN(u,v)
x = pihalf * w / f

!     CALCULATION OF SERIES FOR W = MIN(U,V)

i = 1.0
sum = 0.0
10 ai = qn / (1.0-qd)
a = ai * SIN(i*x)
sum = sum + a
IF (ABS(ai) >= tol*ABS(sum)) THEN
  qn = qn * q1
  qd = qd * q2
  i = i + 2.0
  GO TO 10
END IF

!     ASSEMBLY FOR U <= V

s = coef * sum
c = SQRT(1.0-s*s)
d = SQRT(1.0-(k*s)**2)
IF (u == w) RETURN

!     ASSEMBLY FOR U > V

temp = s
s = c / d
c = l * temp / d
d = l / d
RETURN
END SUBROUTINE scdf


SUBROUTINE scdj(u, k, l, f, f1, s, c, d)
!     ----------------------------------------------------------------
!     SCDJ COMPUTES SN(U,K), CN(U,K), AND DN(U,K) FOR REAL U AND
!     K USING THE IMAGINARY TRANSFORMATION OF JACOBI AND A
!     FOURIER EXPANSION.  IT IS ASSUMED THAT 0.0 <= K < 1.0
!     AND 0.0 <= U <= F, WHERE F = F(K) IS THE COMPLETE ELLIPTIC
!     INTEGRAL OF THE FIRST KIND AND F1 = F(L) IS THE COMPLEMENTARY
!     INTEGRAL, AND THAT L .NE. 0. AND K**2 + L**2 = 1.
!     ----------------------------------------------------------------
REAL, INTENT(IN)  :: u, k, l, f, f1
REAL, INTENT(OUT) :: s, c, d

! Local variables
REAL            :: a, coef, e1, e2, e1n, e2n, eps, n, q1, q2, q1n, q2n,  &
                   sh, sum, temp, tol, v, w, x, x2, xn
REAL, PARAMETER :: pihalf = 1.5707963267949, pi = 3.1415926535898
!     ------------------------------------------------
!     EPS IS A MACHINE DEPENDENT CONSTANT.  EPS IS THE
!     SMALLEST NUMBER SUCH THAT 1. + EPS > 1.

eps = EPSILON(1.0)

!     ------------------------------------------------
tol = eps / 10.0
v = f - u
q1 = -EXP(-pi*f/f1)
q2 = q1 * q1

w = MIN(u,v)
e1 = pi * MAX(u,v) / f1
e2 = pi * (f+w) / f1
e1 = EXP(-e1)
e2 = EXP(-e2)

coef = pihalf / (k*f1)
x = pihalf * w / f1
x2 = 2.0 * x

!     CALCULATION OF SERIES FOR W = MIN(U,V)

n = 1.0
q1n = q1
q2n = q2
e1n = e1
e2n = e2
sum = 0.0

10 xn = n * x2
IF (xn <= 1.0) THEN
  CALL snhcsh(xn, -1, sinhm=sh)
  sh = sh + xn
  a = 2.0 * q1n * ABS(q1n) * sh / (1.0 + q2n)
ELSE
  a = q1n * (e1n-e2n) / (1.0+q2n)
END IF
sum = sum + a
IF (ABS(a) >= tol*ABS(sum)) THEN
  q1n = q1n * q1
  q2n = q2n * q2
  e1n = e1n * e1
  e2n = e2n * e2
  n = n + 1.0
  GO TO 10
END IF

!     ASSEMBLY FOR U <= V

s = coef * (TANH(x)+2.0*sum)
c = SQRT(1.0-s*s)
d = SQRT(1.0-(k*s)**2)
IF (u == w) RETURN

!     ASSEMBLY FOR U > V

temp = s
s = c / d
c = l * temp / d
d = l / d
RETURN
END SUBROUTINE scdj


SUBROUTINE elpfc1(u, k, l, s, c, d, ierr)
!     -------------------------------------------------------------
!     ELPFC1 CALCULATES THE ELLIPTIC FUNCTIONS SN(U,K), CN(U,K),
!     DN(U,K) FOR COMPLEX U AND REAL MODULUS K.  IT IS ASSUMED THAT
!     ABS(K) <= 1. AND K**2 + L**2 = 1.
!     -------------------------------------------------------------
COMPLEX, INTENT(IN)  :: u
COMPLEX, INTENT(OUT) :: s, c, d
REAL, INTENT(IN)     :: k, l
INTEGER, INTENT(OUT) :: ierr

! Local variables
REAL :: ck, cl, coef, c1, c2, dk, dl, d1, d2, k2, r, sk, sl, s1, s2,  &
        t, td1, td2, t1, t2, u1, u2

u1 = REAL(u)
u2 = AIMAG(u)
k2 = k * k
IF (u1 /= 0.0) THEN
  IF (u2 /= 0.0) GO TO 10

!     CALCULATION FOR U2 = 0.

  CALL ellpf(u1, k, l, s1, c1, d1, ierr)
  IF (ierr /= 0) RETURN
  s2 = 0.0
  c2 = 0.0
  d2 = 0.0
  GO TO 20
END IF

!     CALCULATION FOR U1 = 0.

CALL ellpf(u2,l,k,s2,c2,d2,ierr)
IF (ierr /= 0) RETURN
IF (c2 == 0.0) GO TO 30
s1 = 0.0
s2 = s2 / c2
d1 = d2 / c2
d2 = 0.0
c1 = 1.0 / c2
c2 = 0.0
GO TO 20

!     CALCULATION FOR U1 AND U2 .NE. 0.

10 CALL ellpf(u1,k,l,sk,ck,dk,ierr)
IF (ierr /= 0) RETURN
CALL ellpf(u2,l,k,sl,cl,dl,ierr)
IF (ierr /= 0) RETURN
coef = ABS(k) * sl
t1 = cl
t2 = coef * sk
td1 = coef * t1
td2 = coef * t2
IF (ABS(t2) > ABS(t1)) THEN
  IF (t2 == 0.0) GO TO 30
  IF (td2 == 0.0) GO TO 30
  t = t1 / t2
  r = 1.0 / (1.0+t*t)
  s1 = dl * r / td2
  s2 = ck * dk * sl * t * r / t2
  c1 = ck * t * r / t2
  c2 = -dk * sl * dl * r / td2
  d1 = dk * dl * t * r / t2
  d2 = -k2 * ck * sl * r / td2
ELSE
  IF (t1 == 0.0) GO TO 30
  IF (td1 == 0.0) GO TO 30
  t = t2 / t1
  r = 1.0 / (1.0+t*t)
  s1 = dl * t * r / td1
  s2 = ck * dk * sl * r / t1
  c1 = ck * r / t1
  c2 = -dk * sl * dl * t * r / td1
  d1 = dk * dl * r / t1
  d2 = -k2 * ck * sl * t * r / td1
END IF

!     FINAL ASSEMBLY

20 s = CMPLX(s1,s2)
c = CMPLX(c1,c2)
d = CMPLX(d1,d2)
RETURN

!     ERROR RETURN

30 ierr = 3
RETURN
END SUBROUTINE elpfc1


SUBROUTINE peq(z, w, ierr)

!     WEIERSTRASS P-FUNCTION IN THE EQUIANHARMONIC CASE
!     FOR COMPLEX ARGUMENT WITH UNIT PERIOD PARALLELOGRAM

COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w
INTEGER, INTENT(OUT) :: ierr

! Local variables
COMPLEX :: z1, z4, z6
REAL    :: zr, zi
INTEGER :: m, n

!     REDUCTION TO FUNDAMENTAL PARALLELOGRAM

zi = 1.1547005383792515E0 * AIMAG(z) + 0.5E0
m = INT(zi)
IF (zi < 0.E0) m = m - 1
zr = REAL(z) - 0.5E0 * REAL(m) + 0.5E0
n = INT(zr)
IF (zr < 0.E0) n = n - 1
z1 = z - REAL(n) - (0.5E0,0.86602540378443865E0) * REAL(m)

!     IF Z1=0 THEN Z COINCIDES WITH A LATTICE POINT.
!     THE LATTICE POINTS ARE POLES FOR P.

w = z1 * z1
zr = ABS(REAL(w)) + ABS(AIMAG(w))
IF (zr == 0E0) THEN
  ierr = 1
  RETURN
END IF

!     EVALUATION OF P(Z1)

ierr = 0
z4 = w * w
z6 = z4 * w
w = 1E0 / w + 6E0 * z4 * (5E0+z6) / (1E0-z6) ** 2 + z4 * (((((  &
    -2.6427662E-10*z6 + 1.610954818E-8)*z6 + 7.38610752879E-6)*z6 +  &
    4.3991444671178E-4)*z6 + 7.477288220490697E-2)*z6 -  &
    6.8484153287299201E-1) / (((((6.2252191E-10*z6 + 2.553314573E-7)  &
    *z6 - 2.619832920421E-5)*z6 - 5.6444801847646E-4)*z6 +  &
    4.565553484820106E-2)*z6 + 1.0)
RETURN
END SUBROUTINE peq


SUBROUTINE peq1(z, w, ierr)

!     FIRST DERIVATIVE OF WEIERSTRASS P-FUNCTION IN THE EQUIANHARMONIC CASE
!     FOR COMPLEX ARGUMENT WITH UNIT PERIOD PARALLELOGRAM

COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w
INTEGER, INTENT(OUT) :: ierr

! Local variables
COMPLEX :: z1, z3, z6
REAL    :: zr, zi
INTEGER :: m, n

!     REDUCTION TO FUNDAMENTAL PARALLELOGRAM

zi = 1.1547005383792515 * AIMAG(z) + 0.5
m = INT(zi)
IF (zi < 0.0) m = m - 1
zr = REAL(z) - 0.5 * REAL(m) + 0.5
n = INT(zr)
IF (zr < 0.0) n = n - 1
z1 = z - REAL(n) - (0.5, 0.86602540378443865) * REAL(m)

!     IF Z1=0 THEN Z COINCIDES WITH A LATTICE POINT.
!     THE LATTICE POINTS ARE POLES FOR DP.

z3 = z1 * z1 * z1
z6 = z3 * z3
w = z3 * (1.0 - z6) ** 3
zr = ABS(REAL(w)) + ABS(AIMAG(w))
IF (zr == 0.0) THEN
  ierr = 1
  RETURN
END IF

!     EVALUATION OF DP(Z1)

ierr = 0
w = (((14*z6+294)*z6 + 126)*z6 - 2.0) / w + z3 *   &
    ((((((-2.95539175E-9*z6 - 2.6764693031E-7)*z6 + 2.402192743346E-5)*z6 + &
    1.9656661451391E-4)*z6 + 1.760135529461036E-2)*z6 +  &
    8.1026243498822636E-1)*z6 - 2.73936613149196804) /   &
    ((((((4.6397763E-10*z6 + 5.413482233E-8)*z6 - 1.56293298374E-6)*z6 -  &
    1.0393701076352E-4)*z6 + 9.5553182532237E-4)*z6 + 9.131106969640212E-2)*z6 + &
    1.0)
RETURN
END SUBROUTINE peq1


SUBROUTINE plem(z, w, ierr)

!     WEIERSTRASS P-FUNCTION IN THE LEMNISCATIC CASE
!     FOR COMPLEX ARGUMENT WITH UNIT PERIOD PARALLELOGRAM

COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w
INTEGER, INTENT(OUT) :: ierr

! Local variables
COMPLEX :: z1, z4
REAL    :: zr, zi
INTEGER :: m, n

!     REDUCTION TO FUNDAMENTAL PARALLELOGRAM

zr = REAL(z) + 0.5
zi = AIMAG(z) + 0.5
m = INT(zr)
n = INT(zi)
IF (zr < 0.0) m = m - 1
IF (zi < 0.0) n = n - 1
z1 = z - REAL(m) - (0.0,1.0) * REAL(n)

!     IF Z1=0 THEN Z COINCIDES WITH A LATTICE POINT.
!     THE LATTICE POINTS ARE POLES FOR P.

w = z1 * z1
zr = ABS(REAL(w)) + ABS(AIMAG(w))
IF (zr == 0.0) THEN
  ierr = 1
  RETURN
END IF

!     EVALUATION OF P(Z1)

ierr = 0
z4 = w * w
w = 1.0 / w + 4.0 * w * (3.0+z4) / (1.0-z4) ** 2 + w * ((((((((  &
    -7.233108E-11*z4 + 1.714197273E-8)*z4 - 2.5369036492E-7)*z4 -  &
    7.98710206868E-6)*z4 + 6.4850606909737E-4)*z4 +  &
    7.39624629362938E-3)*z4 + 2.012382768497244E-2)*z4 +  &
    7.1177297543136598E-1)*z4 - 2.54636399353830738) / ((((((((  &
    5.1161516E-10*z4 + 6.61289408E-9)*z4 + 4.4618987048E-7)*z4 -  &
    8.42694918892E-6)*z4 + 4.42886829095E-6)*z4 - 4.22629935217101E-3)  &
    *z4 + 2.577496871700433E-2)*z4 + 4.2359940482277074E-1)*z4 + 1.0)
RETURN
END SUBROUTINE plem


SUBROUTINE plem1(z, w, ierr)

!     FIRST DERIVATIVE OF WEIERSTRASS P-FUNCTION IN THE LEMNISCATIC CASE
!     FOR COMPLEX ARGUMENT WITH UNIT PERIOD PARALLELOGRAM

COMPLEX, INTENT(IN)  :: z
COMPLEX, INTENT(OUT) :: w
INTEGER, INTENT(OUT) :: ierr

! Local variables
COMPLEX :: z1, z3, z4
REAL    :: zr, zi
INTEGER :: m, n

!     REDUCTION TO FUNDAMENTAL PARALLELOGRAM

zr = REAL(z) + 0.5
zi = AIMAG(z) + 0.5
m = INT(zr)
n = INT(zi)
IF (zr < 0.0) m = m - 1
IF (zi < 0.0) n = n - 1
z1 = z - REAL(m) - (0.0,1.0) * REAL(n)

!     IF Z1=0 THEN Z COINCIDES WITH A LATTICE POINT.
!     THE LATTICE POINTS ARE POLES FOR DP.

z3 = z1 * z1 * z1
z4 = z3 * z1
w = (z1*(1.0-z4)) ** 3
zr = ABS(REAL(w)) + ABS(AIMAG(w))
IF (zr == 0.0) THEN
  ierr = 1
  RETURN
END IF

!     EVALUATION OF DP(Z1)

ierr = 0
w = (((10.*z4 + 90.)*z4 + 30.)*z4 - 20.) /  &
    w + z1 * ((((((((((  &
    -3.9046302E-9*z4 - 1.001487137E-8)*z4 + 5.9573043092E-7)*z4 -  &
    2.482518130524E-5)*z4 + 1.4557266595395E-4)*z4 +   &
    4.56633655643206E-3)*z4 + 6.224782572111135E-2)*z4 +   &
    1.038527937794269E-2)*z4 + 1.19804620802637942)*z4 +   &
    6.42791439683811718)*z4 - 5.09272798707661477) / ((((((((((  &
    4.726888E-11*z4 - 3.0667983E-9)*z4 + 1.0087596089E-7)*z4-  &
    8.060683451E-8)*z4 + 1.184299251664E-5)*z4 - 2.3096723361547E-4)*  &
    z4 - 2.90730903142055E-3)*z4 + 1.338392411135511E-2)*z4 +   &
    2.3098639320021426E-1)*z4 + 8.4719880964554148E-1)*z4 + 1.0)
RETURN
END SUBROUTINE plem1


SUBROUTINE valr2(x, y, n0, p, iop, a, ind, ko)
!     -------------------
! Integral of the uncorrelated bivariate normal density function over arbitrary
! polygons and semi-infinite angular regions (pie-shaped wedges).  That is, it
! is the integral of:
!                      [1/(2.pi)].exp[-(x^2 + y^2)/2]
!
! Input
! If n0 = 1, the integral is over a wedge-shaped region.
!            (x(1), y(1)) is the point of the wedge, (x(2), y(2)) is any point
!            on one edge and (x(3), y(3)) is any point on the other edge.
! If n0 >= 3, the integral is over a polygon.
!            The polygon is defined by the points (x(i), y(i)), i=1, ..., n.
!            The routine adds (x(n+1), y(n+1)) = (x(1), y(1)).
! iop = accuracy required.   The integral (p) is calculated to an accuracy
!       of at least 3 times iop decimal digits for iop = 1, 2 or 3.
!
! Output
! p   = estimate of the integral.
! a   = area of the polygon if n0 >= 3, otherwise zero.
! ind = error indicator
!     = 0 if integral has been calculated to the required accuracy
!     = 1 if n0 = 1 and (x(1), y(1)) was equal to, or too close to one of
!         the other two points to determine the angle.  p = 5.0 is returned.
!     = 2 if n0 = 1 and the angle = pi (180 degrees), or n0 >= 3 and one of
!         the angles of the polygon is equal to pi.  The value of p returned
!         in this case is correct.
!     = 3 illegal value for n0 (either < 1 or = 2)
!
! Fortran 77 version by Armido DiDonato, Richard Hageman and Alfred Morris
! from the Naval Surface Warfare Mathematics Library.
!     -------------------

REAL, INTENT(IN OUT) :: x(:), y(:)
INTEGER, INTENT(IN)  :: n0, iop
REAL, INTENT(OUT)    :: p, a
INTEGER, INTENT(OUT) :: ind, ko

! Local variables
INTEGER :: n, k, m
REAL    :: g(2), h(2), kom, l, tau, tausq, w, z, v, xk,     &
           psi1, t1, u, yk, d1sq, d2sq, bgd1, bgd2, aj0, b, c, cee, p1,  &
           cape, caph, f, aj1, circm, sum, t, capv, ykm1
!     -------------------
REAL    :: pi = 3.1415926535898, twopi = 6.28318530717958,  &
           alnpi = 1.14472988584940, rtpi = 1.77245385090552,  &
           rtpii = .56418958354776, &
           e(5) = (/ .885777518572895E+00, -.981151952778050E+00,  &
                     .759305502082485E+00, -.353644980686977E+00,  &
                     .695232092435207E-01 /),  &
           e2(10) = (/ .886226470016632E+00, -.999950714561036E+00,  &
                       .885348820003892E+00, -.660611239043357E+00,  &
                       .421821197160099E+00, -.222898055667208E+00,  &
                       .905057384150449E-01, -.254906111884287E-01,  &
                       .430895168984138E-02, -.323377239693247E-03 /),  &
           e3(15) = (/ .886226924931465E+00, -.999999899776252E+00,  &
                       .886223733186722E+00, -.666626670510907E+00,  &
                       .442851899328569E+00, -.265638206366025E+00,  &
                       .145060043403014E+00, -.714909837799889E-01,  &
                       .309199295521210E-01, -.112323532148441E-01,  &
                       .324944543171185E-02, -.704260243309096E-03,  &
                       .105787574480633E-03, -.971864864160461E-05,  &
                       .408335517232165E-06 /),  &
           aph1(3) = (/ 2.02E-7, 2.08E-13, 2.71E-19 /),  &
           aph2(3) = (/ 1.22E-2, 1.23E-4,  1.34E-6 /),   &
           aph4(3) = (/ .6962E-1, .6990E-2, .7311E-3 /), &
           rsq(3)  = (/ 6.0516, 12.60605, 19.201924 /),  &
           a3d8(3) = (/ 0.28125E-4, 0.285E-7, 0.32625E-10 /)

!     -------------------
!     TAU IS A MACHINE DEPENDENT TOLERANCE. IT IS ASSUMED THAT A
!     7 OR MORE DIGIT FLOATING POINT ARITHMETIC IS BEING USED.
!     -------------------
tau=2.0*spmpar(1)
IF (tau < 3.e-11) tau=MAX(5.0*tau, 1.e-14)
!     -------------------
n=n0
IF (n == 2.OR.n < 1) GO TO 4021
tausq=tau*tau

p=0.0
ind=0
a=0.0
kom=0.0
k=1
IF (n /= 1) GO TO 10

w=x(2)-x(1)
z=y(2)-y(1)
u=x(3)-x(1)
v=y(3)-y(1)
xk=0.0
psi1=v*w - u*z
IF (psi1 >= 0.0) GO TO 21

p=-1.0
t1=w
w=u
u=t1
t1=v
v=z
z=t1
GO TO 21

10 ko=0
x(n+1)=x(1)
y(n+1)=y(1)
u=x(2)-x(1)
v=y(2)-y(1)
xk=x(1)
yk=y(1)

20 w=x(1)-x(n)
z=y(1)-y(n)
21 d1sq=w*w + z*z
IF (d1sq > tausq) GO TO 30
IF (n == 1) GO TO 4011
n=n-1
IF (n == 2) RETURN
GO TO 20

30 d2sq=u*u+v*v
IF (d2sq > tausq) GO TO 40
IF (n == 1) GO TO 4011
31 k=k+1
u=x(k+1)-xk
v=y(k+1)-yk
d2sq=u*u + v*v
IF (d2sq <= tausq) GO TO 31
IF (k == n-1) RETURN

40 a=xk*(y(k+1)-y(n))
bgd1=SQRT(d1sq+d1sq)
bgd2=SQRT(d2sq+d2sq)

!             PROCESSING VERTEX (XK,YK)

50 psi1=v*w-u*z
cee=u*w+v*z
aj0=ATAN2(psi1,cee)
kom=kom+aj0
l=0.0
b=.5*(x(k)*x(k)+y(k)*y(k))
IF (b > aph1(iop)) GO TO 60
p1=aj0/twopi
GO TO 3621

60 g(1)=(w*x(k)+z*y(k))/bgd1
g(2)=(u*x(k)+v*y(k))/bgd2
h(1)=(-y(k)*w+x(k)*z)/bgd1
h(2)=(-y(k)*u+x(k)*v)/bgd2
IF (ABS(psi1) > bgd1*bgd2*a3d8(iop)) GO TO 80
IF (cee < 0.0) GO TO 70
IF (ABS(aj0) > tau .AND. g(1) < 0.0) GO TO 80
p1=0.0
GO TO 3621

70 IF (ABS(psi1) <= (.5*tau*bgd1*bgd2)) ind=2
IF (psi1 < 0.0) GO TO 71
p1=.5*erfc1(0,h(2))
GO TO 3621
71 p1=-.5*erfc1(0,h(1))
GO TO 3621

80 IF (b > aph2(iop)) GO TO 90
c=rtpi*(h(2)-h(1))-(g(2)*h(2)-g(1)*h(1))
p1=(aj0-c)/twopi
GO TO 3621

!                 COMPUTATION OF L

90 IF (g(1) < 0.0) GO TO 100
IF (g(2) >= 0.0) GO TO 130
g(2)=-g(2)
h(2)=-h(2)
IF (ABS(h(2)) <= aph4(iop)) GO TO 91
l=.5*erfc1(0,-h(2))
GO TO 120
91 l=.5+rtpii*h(2)
GO TO 120

100 g(1)=-g(1)
h(1)=-h(1)
IF (g(2) < 0.0) GO TO 110
IF (ABS(h(1)) <= aph4(iop)) GO TO 101
l=.5*erfc1(0,h(1))
GO TO 120
101 l=.5-rtpii*h(1)
GO TO 120

110 g(2)=-g(2)
h(2)=-h(2)
IF (ABS(h(1)) <= aph4(iop)) GO TO 112
IF (ABS(h(2)) <= aph4(iop)) GO TO 111
l=.5*(erfc1(0,h(1))-erfc1(0,h(2)))
GO TO 130
111 l=rtpii*h(2)-.5*erf(h(1))
GO TO 130
112 IF (ABS(h(2)) <= aph4(iop)) GO TO 113
l=.5*erf(h(2))-rtpii*h(1)
GO TO 130
113 l=rtpii*(h(2)-h(1))
GO TO 130

120 psi1=-psi1
IF (psi1 <= 0.0) GO TO 121
l=l-1.0
aj0=aj0+pi
GO TO 130
121 aj0=aj0-pi

!                 SERIES EVALUATION

130 IF (b >= rsq(iop)) GO TO 171
cape=aj0
caph=.5*aj0
m=1
f=0.0
aj1=h(2)-h(1)
circm=aj1
IF (iop-2 < 0) THEN
  GO TO 140
ELSE IF (iop-2 == 0) THEN
  GO TO 150
ELSE
  GO TO 160
END IF

140 sum=e(m)*aj1
141 m=m+1
h(1)=h(1)*g(1)
h(2)=h(2)*g(2)
t=h(2)-h(1)
f=f+b
capv=(f*cape+t)/m
sum=sum+e(m)*capv
IF (m >= 5) GO TO 170
cape=circm
circm=capv
GO TO 141

150 sum=e2(m)*aj1
151 m=m+1
h(1)=h(1)*g(1)
h(2)=h(2)*g(2)
t=h(2)-h(1)
f=f+b
capv=(f*cape+t)/m
sum=sum+e2(m)*capv
IF (m >= 10) GO TO 170
cape=circm
circm=capv
GO TO 151

160 sum=e3(m)*aj1
161 m=m+1
h(1)=h(1)*g(1)
h(2)=h(2)*g(2)
t=h(2)-h(1)
f=f+b
capv=(f*cape+t)/m
sum=sum+e3(m)*capv
IF (m >= 15) GO TO 170
cape=circm
circm=capv
GO TO 161

170 p1=l+EXP(-(b+alnpi))*(caph-sum)
GO TO 3621
171 p1=l

!               STANDARD TERMINATION

3621 IF (k /= n) GO TO 3651
IF (n /= 1) GO TO 3631
p=ABS(p+ABS(p1))
RETURN

3631 p=p-p1
kom=kom/twopi
a=.5*a
IF (kom < 0.0) GO TO 3641
ko=INT(kom+.125)
GO TO 3645
3641 ko=INT(kom-.125)
3645 p=p+REAL(ko)
RETURN

!              SET UP THE NEXT VERTEX

3651 w=u
z=v
bgd1=bgd2
xk=x(k+1)
yk=y(k+1)
ykm1=y(k)
3661 k=k+1
u=x(k+1)-xk
v=y(k+1)-yk
d2sq=u*u+v*v
IF (d2sq <= tausq) GO TO 3661
bgd2=SQRT(d2sq+d2sq)
p=p-p1
a=a+xk*(y(k+1)-ykm1)
GO TO 50

!                  ERROR RETURN

4011 ind=1
p=5.0
RETURN
4021 ind=3
RETURN
END SUBROUTINE valr2



SUBROUTINE circv(r, d, j, p, ierr)
!-----------------------------------------------------------------------
!       IF J .NE. 0, OUTPUT IS P = CIRCULAR COVERAGE FUNCTION. P GIVES
!       THE PROBABILITY OF A SHOT FALLING, UNDER A NORMAL DISTRIBUTION
!       WITH MEAN (0,0) AND EQUAL STANDARD DEVIATIONS, S, IN A CIRCLE
!       OF RADIUS R0, OFFSET A DISTANCE D0 FROM (0,0).
!       INPUT IS R = R0/S, D = D0/S.

!       IF J = 0, OUTPUT IS P = GENERALIZED CIRCULAR ERROR FUNCTION.
!       P GIVES THE PROBABILITY OF A SHOT FALLING ,UNDER A NORMAL
!       BIVARIATE DISTRIBUTION WITH MEAN (0,0) AND STANDARD DEVIATIONS
!       SMIN AND S, IN A CIRCLE OF RADIUS R0 CENTERED AT (0,0).
!       INPUT FOR J = 0, R = R0/S, D = SMIN/S <= 1.
!       IF SMIN = 0, S .NE. 0, P = ERF(R/(SQR(2)).

!       IF IERR .NE. 0, SOME PORTION OF THE INPUT IS UNACCEPTABLE.
!       IF R < 0., THEN CIRCV SETS IERR = 1.
!       IF D < 0, OR J = 0 AND D > 1., THEN CIRCV SETS IERR = 2.

!       REFERENCES
!       MATH OF COMP APRIL 1961,PP169,173 AND OCT.1961, PP 375, 382.
!       NWL REPORT N0.1768, JAN. 1962. NSWC REPORT N0.83-13, NOV. 1982.
!       IEEE TRANS. INFO. TH. APRIL 1965, P. 312.
!-----------------------------------------------------------------------
!       NEGATIVE R AND D ARE NOT PERMITTED.
!--------------------------------------------
REAL, INTENT(IN)     :: r, d
INTEGER, INTENT(IN)  :: j
REAL, INTENT(OUT)    :: p
INTEGER, INTENT(OUT) :: ierr

! Local variable
REAL            :: a, a1, an, anm1, bk2, const, c0, delta, d2, eps, eps0, &
                   expt, f, m, m0, s, s0, s1, s2, sum, t, td, tn, tr,  &
                   t0, t2, w, x, x0, y1, z, zm, zp, zr, zrs
INTEGER         :: i, n
REAL, PARAMETER :: e = 2.71828182845905, c1 = .707106781186548,   &
                   c2 = .564189583547756, c3 = 6.28318530717959
!-----------------------
!     C1 = 1/SQRT(2)
!     C2 = 1/SQRT(PI)
!     C3 = 2*PI
!-----------------------
p = 0.0
IF (r >= 0.0) THEN
  IF (d < 0.0) GO TO 110
  ierr = 0
  IF (r == 0.0) RETURN

  eps0 = EPSILON(1.0)
  z = -LOG(eps0)
  IF (j == 0) THEN
!------------------------------------------------------------------
!   FOR J = 0, ERROR IN D IF D < 0 OR D > 1
!------------------------------------------------------------------
!                J = 0
!--------------------------------------------------
    IF (d > 1.0) GO TO 110
    IF (d == 1.0) GO TO 30
    IF (d == 0.0) THEN
!---------------------------------------
!          J = 0,     D = 0
!---------------------------------------
      p = erf(r*c1)
      RETURN
    END IF
!-------------------------------------------------
!   J = 0,  (R*R > -2*LOG(EPS0))  P = 1
!-------------------------------------------------
    IF (r*r >= 2.0*z) THEN
      p = 1.0
      RETURN
    END IF
!---------------------------------------------
    x = r
    d2 = d * d
    y1 = r / d
    t = 0.5 * y1
    zm = (0.5-d2) + 0.5
    t2 = t * t
    t = t2 * zm

    eps = 10.0 * eps0
    IF (t <= 14.0) THEN
!-------------------------------------------------
!    J = 0,   T <= 14
!-------------------------------------------------
      zp = 1.0 + d2
      zr = zm / zp
      zrs = zr * zr
      bk2 = t2 * zp
      c0 = 2.0 * d / zp
      s0 = EXP(-bk2)
      t0 = 0.5 + (0.5-s0)
      IF (bk2 <= 0.15) t0 = -rexp(-bk2)
      t0 = c0 * t0
      s0 = c0 * s0

      p = t0
      an = 0.0
      10 an = an + 2.0
      f = (an-1.0) / an
      w = t / an
      x = s0 * w
      t0 = f * zrs * t0 - (w+zr) * x
      s0 = w * x
      p = p + t0
      IF (t0 > eps*p) GO TO 10
      RETURN
    END IF
!------------------------------------------------
!   J = 0,    T > 14
!------------------------------------------------
    t = 0.25 / t
    const = d2 / zm
    delta = SQRT(zm)
    x = 2.0 * (c1*c2) / (r*delta)
    expt = EXP(-0.5*r*r)
    CALL erfc0(1,c1*r,expt,m)
    m = m / delta

    p = 1.0
    sum = m
    an = 0.0
    IF (expt*m >= 5.e-3) THEN

!          SET  P = 1 - EXP(-R*R/2)*M

      p = (erf(c1*r)-d2/(1.0+delta)) / delta
      sum = 0.0
    END IF

!          COMPUTE THE ASYMPTOTIC EXPANSION

    20 an = an + 2.0
    anm1 = an - 1.0
    f = anm1 / an
    m0 = m
    m = const * f * (x-m)
    IF (m < m0 .AND. m >= 0.0) THEN
      x = anm1 * f * t * x
      sum = sum + m
      IF (m > eps*sum) GO TO 20
    END IF
    p = p - expt * sum
    GO TO 100
  END IF
!------------------------------------
!              J .NE. 0
!------------------------------------
  a1 = r - d
  t = r * d
  IF (d /= 0.0) GO TO 40
!----------------------------------------------
!     J = 0 AND D = 1, OR J .NE. 0 AND D = 0
!----------------------------------------------
  30 p = -rexp(-0.5*r*r)
  RETURN
!----------------------------------------
  40 a = 0.5 * (a1*a1)
  IF (a1 > 5.0) THEN
    IF (a <= z) GO TO 50
    p = 1.0
    RETURN
  END IF
  IF (a1 < -5.0 .AND. a > -exparg(1)) RETURN
!----------------------------------------
  50 eps = 1.5E2 * eps0
  IF (r < 1.7.OR.t <= 16.0) THEN
!------------------------------------------------------------------
!      J .NE. 0,    R < 1.7  OR  R*D <= 16
!------------------------------------------------------------------
    tr = 0.5 * (r*r)
    td = 0.5 * (d*d)

!       FIND THE NUMBER N OF TERMS TO BE USED IN THE SERIES

    z = 0.5 * t
    n = z * e + 1.0
    tn = ((e*z/n)**(n+n)) / (c3*n)
    60 n = n + 1
    w = z / n
    tn = tn * (w*w)
    IF (tn > eps) GO TO 60

!       COMPUTE THE SERIES

    m = n
    CALL gratio(m+1.0,tr,s,w,0)
    w = rcomp(m,tr) / m
    p = s
    DO i = 1, n
      s = s + w
      w = (m/tr) * w
      p = s + (td/m) * p
      m = m - 1.0
    END DO
    p = EXP(-td) * p
  ELSE
!------------------------------------------------------------------
!      J .NE. 0,    R >= 1.7  AND  R*D > 16
!------------------------------------------------------------------
    z = c1 * ABS(a1)
    w = EXP(-a)
    CALL erfc0(1,z,w,s1)

    a = 0.5 * z / t
    t = 0.25 / t
    m = c2 - z * s1
    IF (z >= 4.0) m = erfcr(z)
    m = 0.5 * a * m
    x = 0.5 * c2 * t
    s1 = s1 + m
    s2 = c2 + x

    an = 2.0
    80 an = an + 2.0
    anm1 = an - 1.0
    f = anm1 / an
    m0 = m
    m = f * a * (x-z*m)
    x = f * (anm1*t) * x
    s2 = s2 + x
    IF (m > 0.0 .AND. m < m0) THEN
      s1 = s1 + m
      IF (m > eps*s1) GO TO 80
    END IF

    90 an = an + 2.0
    anm1 = an - 1.0
    x0 = x
    x = (anm1/an) * (anm1*t) * x
    IF (x < x0) THEN
      s2 = s2 + x
      IF (x > eps*s2) GO TO 90
    END IF

    s1 = 0.5 * (r+d) * s1
    s2 = c1 * s2
    w = w / SQRT(r*d)
    p = 0.5 * w * ABS(s1-s2)
    IF (a1 > 0.0) p = ABS(1.0-0.5*w*(s1+s2))
  END IF

!       TERMINATION

  100 IF (p > 1.0) p = 1.0
  RETURN
END IF
ierr = 1
p = -1.0
RETURN
110 ierr = 2
p = -1.0
RETURN
END SUBROUTINE circv


SUBROUTINE erfc0(ind, x, e, y)
!-----------------------------------------------------------------------
!         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

!               Y = ERFC(X)            IF IND = 0
!               Y = EXP(X*X)*ERFC(X)   OTHERWISE

!     E IS AN INPUT/OUTPUT VARIABLE. IF E >= 0 THEN IT IS ASSUMED
!     THAT E = EXP(-X*X). IN THIS CASE E IS NOT MODIFIED. IF E IS
!     NEGATIVE THEN E IS SET TO EXP(-X*X) WHEN THIS VALUE IS NEEDED.

!-----------------------------------------------------------------------
INTEGER, INTENT(IN)  :: ind
REAL, INTENT(IN)     :: x
REAL, INTENT(IN OUT) :: e
REAL, INTENT(OUT)    :: y

! Local variables

REAL, PARAMETER :: a(4) = (/ -1.65581836870402E-4, 3.25324098357738E-2,    &
                              1.02201136918406E-1, 1.12837916709552 /), &
                   b(4) = (/ 4.64988945913179E-3, 7.01333417158511E-2,     &
                             4.23906732683201E-1, 1.00000000000000 /),  &
                   p(8) = (/ -1.36864857382717E-7, 5.64195517478974E-1,    &
                              7.21175825088309, 4.31622272220567E01,    &
                              1.52989285046940E02, 3.39320816734344E02,    &
                              4.51918953711873E02, 3.00459261020162E02 /), &
                   q(8) = (/  1.00000000000000, 1.27827273196294E01,    &
                              7.70001529352295E01, 2.77585444743988E02,    &
                              6.38980264465631E02, 9.31354094850610E02,    &
                              7.90950925327898E02, 3.00459260956983E02 /), &
                   r(5) = (/  2.10144126479064, 2.62370141675169E01,    &
                              2.13688200555087E01, 4.65807828718470,    &
                              2.82094791773523E-1 /),  &
                   s(5) = (/  9.41537750555460E01, 1.87114811799590E02,    &
                              9.90191814623914E01, 1.80124575948747E01,    &
                              1.00000000000000 /), c = .564189583547756
REAL (dp) :: w
REAL      :: ax, t, top, bot, eps
!-------------------------

!                     ABS(X) < 0.47

ax = ABS(x)
IF (ax < 0.47) THEN
  t = x * x
  top = ((a(1)*t + a(2))*t + a(3)) * t + a(4)
  bot = ((b(1)*t + b(2))*t + b(3)) * t + b(4)
  y = 0.5 + (0.5 - x*top/bot)
  IF (ind == 0) RETURN

  IF (e < 0.0) e = EXP(-t)
  y = y / e
  RETURN
END IF

!                  0.47 <= ABS(X) <= 4

IF (ax > 4.0) GO TO 20
top = ((((((p(1)*ax + p(2))*ax + p(3))*ax + p(4))*ax + p(5))*ax + p(6))*ax +  &
      p(7)) * ax + p(8)
bot = ((((((q(1)*ax + q(2))*ax + q(3))*ax + q(4))*ax + q(5))*ax + q(6))*ax +  &
      q(7)) * ax + q(8)
y = top / bot

10 IF (ind /= 0) THEN
  IF (x >= 0.0) RETURN
  IF (e < 0.0) e = EXP(-x*x)
  y = 2.0 / e - y
  RETURN
END IF
w = DBLE(x) * DBLE(x)
t = w
eps = w - DBLE(t)
IF (e < 0.0) e = EXP(-t)
y = ((0.5 + (0.5-eps))*e) * y
IF (x < 0.0) y = 2.0 - y
RETURN

!                      ABS(X) > 4

20 IF (x > -5.5) THEN
  IF (ind == 0 .AND. x > 50.0) GO TO 30
  t = (1.0/x) ** 2
  top = (((r(1)*t + r(2))*t + r(3))*t + r(4)) * t + r(5)
  bot = (((s(1)*t + s(2))*t + s(3))*t + s(4)) * t + s(5)
  y = (c-t*top/bot) / ax
  GO TO 10
END IF

!             LIMIT VALUE FOR LARGE NEGATIVE X

y = 2.0
IF (ind == 0) RETURN
IF (e < 0.0) e = EXP(-x*x)
y = 2.0 / e
RETURN

!             LIMIT VALUE FOR LARGE POSITIVE X
!                       WHEN IND = 0

30 y = 0.0
RETURN
END SUBROUTINE erfc0


FUNCTION erfcr(x) RESULT(fn_val)
!-----------------------------------------------------------------------
!     COMPUTATION OF  1/SQRT(PI) - X*EXP(X*X)*ERFC(X)  FOR X >= 4
!-----------------------------------------------------------------------
REAL, INTENT(IN) :: x
REAL             :: fn_val

! Local variables
REAL, PARAMETER :: r(5) = (/ 2.10144126479064E+00, 2.62370141675169E+01,  &
                             2.13688200555087E+01, 4.65807828718470E+00,  &
                             2.82094791773523E-01 /),  &
                   s(4) = (/ 9.41537750555460E+01, 1.87114811799590E+02,  &
                             9.90191814623914E+01, 1.80124575948747E+01 /)
REAL :: t, top, bot
!------------------------
t = (1.0/x) ** 2
top = (((r(1)*t + r(2))*t + r(3))*t + r(4)) * t + r(5)
bot = (((s(1)*t + s(2))*t + s(3))*t + s(4)) * t + 1.0
fn_val = t * top / bot
RETURN
END FUNCTION erfcr


SUBROUTINE pkill(r0, sx, sy, h, k, p)
!-----------------------------------------------------------------------

!           COMPUTATION OF THE ELLIPTICAL COVERAGE FUNCTION

!                         ---------------

!     THE RESULT P IS ACCURATE TO AT LEAST 6 SIGNIFICANT DIGITS IF
!     P >= 1.E-20 AND MAX(H/S,K/S,SX/S,SY/S) <= 10/SQRT(EPSILON(1.0))
!     FOR S = MIN(SX,SY).

!-----------------------------------------------------------------------
REAL, INTENT(IN)  :: r0, sx, sy, h, k
REAL, INTENT(OUT) :: p

! Local variables
REAL :: v(16) = (/ .4830766568773832E-01, .1444719615827965E+00,  &
                   .2392873622521371E+00, .3318686022821276E+00,  &
                   .4213512761306353E+00, .5068999089322294E+00,  &
                   .5877157572407623E+00, .6630442669302152E+00,  &
                   .7321821187402897E+00, .7944837959679424E+00,  &
                   .8493676137325700E+00, .8963211557660521E+00,  &
                   .9349060759377397E+00, .9647622555875064E+00,  &
                   .9856115115452683E+00, .9972638618494816E+00 /), &
        w(16) = (/ .9654008851472780E-01, .9563872007927486E-01,  &
                   .9384439908080457E-01, .9117387869576388E-01,  &
                   .8765209300440381E-01, .8331192422694676E-01,  &
                   .7819389578707031E-01, .7234579410884851E-01,  &
                   .6582222277636185E-01, .5868409347853555E-01,  &
                   .5099805926237618E-01, .4283589802222668E-01,  &
                   .3427386291302143E-01, .2539206530926206E-01,  &
                   .1627439473090567E-01, .7018610009470097E-02 /),  &
        a3 = 20.0 , c = .4638648042895004, rpi = .5641895835477563,  &
        rt2 = 1.414213562373095
REAL :: a, a1, a2, a4, c8, d, dm, d1, d2, eps, e3, e4, &
        f, g2, g3, h1, h2, h8, h3, h5, h6, k8, q, q1, r, r8, sa, seps1,   &
        sqeps, s0, s1, s2, s9, t, t1, t2, t3, t4, t5, t9, u, x1, x5, yy, &
        z, z1, z2, z3, z8, z9
INTEGER :: j, iy, iz, iz3, i, ii, n1, nt, j8
!-----------------------------------------------------------------------
!   V(*), W(*)-- GAUSSIAN ABSCISSAS AND WEIGHTS OF ORDER 32, ON (-1,1).
!-----------------------------------------------------------------------

p = 0.0
IF (sx <= 0.0 .OR. sy <= 0.0) THEN
  p = -1.e10
  RETURN
END IF
j = 0
eps = EPSILON(1.0)
sqeps = SQRT(eps)
a = 6.5
!----------------------------------------
!     TEST NO.1 (REPORT 83-13).
!----------------------------------------
z3 = MIN(SQRT(SQRT(TINY(1.0))),1.e-30) * sx * sy
IF (r0*r0 <= z3) RETURN

h2 = h * h + k * k
h8 = ABS(h)
k8 = ABS(k)
dm = MAX(sx,sy)
!-----------------------------------------
!     TEST NO.2 (REPORT 83-13).
!-----------------------------------------
IF ((r0-h8 + a3*sx) <= 0.0) RETURN
IF ((r0-k8 + a3*sy) <= 0.0) RETURN
!-----------------------------------------
!     TEST NO.3 (REPORT 83-13).
!-----------------------------------------
c8 = -epsln()
!-----------EXP(-10.3026) = 3.35E(-5)---------
a4 = MAX(c8-1.74249E1,10.3026)
a4 = SQRT(a4+a4)
h3 = (r0-h) * (r0+h)
h5 = (r0-k) * (r0+k)
s0 = a4 * dm
t = r0 - s0
IF (t >= 0.0) THEN
  IF (t*t >= h2) THEN
    IF (r0 >= 1/sqeps) THEN
      t = h3 - 2.0 * r0 * s0 + s0 * s0 - k * k
      IF (ABS(h3) > ABS(h5)) THEN
        t = h5 - 2.0 * r0 * s0 + s0 * s0 - h * h
      END IF
      IF (t < 0.0) GO TO 10
    END IF
    p = 1.0
    RETURN
  END IF
END IF
!-----------------------------------------
!     TEST NO.4 (REPORT 83-13).
!-----------------------------------------
10 s0 = SQRT(h2)
g2 = 0.0
IF (s0 > r0) THEN
  d = ((s0-r0)/dm) ** 2
  IF (r0*r0*EXP(-0.5*d) <= z3) RETURN
END IF
!-----------------------------------------
!     SX - SY < EPS
!-----------------------------------------
IF (ABS(sx-sy) <= 20.0*dm*eps) THEN
  h8 = s0
  k8 = 0.0
  IF ((r0-h8+a3*dm) <= 0.0) RETURN
  IF (r0 >= dm/sqeps) THEN
    j = 1
    g2 = (k*k-h3) / ((r0+h8)*dm)
    IF (ABS(h3) > ABS(h5)) THEN
      g2 = (h*h-h5) / ((r0+h8)*dm)
    END IF
  END IF
END IF
!--------------------------
!     SMALL R
!--------------------------
t1 = r0 / sx
t2 = r0 / sy
t3 = t1 * t2
t1 = t1 * t1
t2 = t2 * t2
z1 = (h/sx) ** 2
z2 = (k/sy) ** 2
t = t1 * (z1-1.0) + t2 * (z2-1.0)
IF (ABS(t) <= 1.e-3) THEN
  t9 = (t*t-4.0*(t1*t1*(z1-0.5)+t2*t2*(z2-0.5))) / 192.0
  IF (ABS(t9) <= MAX(10.0*eps,1.e-10)) THEN
    p = 0.5 * t3 * (1.0+0.125*t) * EXP(-0.5*(z1+z2))
    RETURN
  END IF
END IF
!-----------------------------------------------
!     NORMALIZE MAX(SX,SY) = 1.
!-----------------------------------------------
r = r0 / dm
s1 = sx / dm
s2 = sy / dm
h8 = h8 / dm
k8 = k8 / dm
h2 = h2 / (dm*dm)
!---------------------------------
!     S1 = 1 >= S2
!---------------------------------
IF (s1 < s2) THEN
  h1 = s1
  s1 = s2
  s2 = h1
  h1 = h8
  h8 = k8
  k8 = h1
END IF
!-----------------------------------------------
!     LIMITING RESULTS FOR MIN(S1,S2) = 0
!-----------------------------------------------
seps1 = MIN(6.71*sqeps,1.e-5)
!-----------------------------------------------
!            R = K,  S2 SMALL
!-----------------------------------------------
IF (k8 == r) THEN
  yy = .166484 * (r*(h8*h8+1.0)+1.0/r) * s2
  IF (ABS(yy) > seps1) GO TO 20
  h1 = h8 / rt2
  p = c * EXP(-h1*h1) * SQRT(k8*s2)
  RETURN
END IF
!-----------------------------------------------
!            R > K,  S2 SMALL
!-----------------------------------------------
IF (r >= k8) THEN
  IF (k8 == 0.0) THEN
    h1 = s2 * s2 / (4.0*rt2*r)
    g2 = g2 / rt2
    IF (j == 0) g2 = ABS(h8-r) / rt2
    IF (ABS(g2) >= 4.0) THEN
      IF (h1*ABS(g2) > seps1) GO TO 20
      p = 0.5 * aerf(h8/rt2,r/rt2)
      RETURN
    END IF
    IF (j == 0) THEN
      p = 0.5 * aerf(h8/rt2,r/rt2)
    ELSE
      IF (h8+r >= rt2) THEN
        p = 0.5 * (erfc(g2)-erfc((h8+r)/rt2))
      ELSE
        p = 0.5 * (erf((h8+r)/rt2)-erf(g2))
      END IF
    END IF
    IF (h1*EXP(-g2*g2) > p*seps1) GO TO 20
    RETURN
  END IF

  z = (r-k8) * (r+k8)
  g2 = SQRT(z)
  h1 = ABS(h8-g2)
  j = 0
  IF (h1 <= 5.0) THEN
    j = 1
    z8 = aerf(h8/rt2,g2/rt2)
    IF (z8 == 0.0) GO TO 20
    h1 = EXP(-0.5*(h8-g2)**2) / z8
  END IF
  h1 = 1.0 / (8.0*rt2*z) * (k8*k8*ABS(h8-g2)+r*r/g2) * s2 * s2 *h1
  IF (ABS(h1) <= seps1) THEN
    u = 0.5
    IF (j == 0) z8 = aerf(h8/rt2,g2/rt2)
    s2 = rt2 * s2
    IF (k8-r > -13.0*s2) u = 0.25 * aerf(k8/s2,r/s2)
    p = u * z8
    RETURN
  END IF
END IF
!----------------------------
!     FIND THE VALUE OF A
!----------------------------
20 s0 = s1 * s1
s9 = s2 * s2
g3 = h8 * h8
g2 = k8 * k8
z8 = s0 + s9
z = s0 * s0 + s9 * s9
h3 = g3 * s0 + g2 * s9
t1 = 2.0 * (z+2*h3)
yy = r * r * (h2+z8) / t1
t1 = (h2+z8) * (h2+z8) / t1
CALL gratio(t1,yy,z,z8,0)
z2 = r / (rt2*s2)
r8 = r / (rt2*s1)
h2 = h8 / (rt2*s1)
h3 = k8 / (rt2*s2)
s0 = 0.0
s9 = 0.0
IF (z >= 1.e-13.OR.s2 <= 5.e-10) THEN
  IF (h2 > 50.0.OR.h3 > 50.0) s0 = 1.5
END IF
u = aerf(h3,z2)
yy = 0.25 * u * aerf(h2,r8)
IF (yy > 0.1) s9 = 0.5
IF (yy < 5.e-15) THEN
  z = yy
ELSE
  z = MIN(yy,z)
END IF

IF (z < 0.5) THEN
  a = a + .5
  IF (z < 5.e-4) THEN
    a = a + .5 + s0
    IF (z < 1.e-6) THEN
      a = a + .5 + s9
      IF (z < 5.e-9) THEN
        a = a + .5
        IF (z < 5.e-10) THEN
          a = a + .25
          IF (z < 5.e-11) THEN
            a = a + .25
            IF (z < 5.e-12) THEN
              a = a + .25
              IF (z < 5.e-13) THEN
                a = a + .5 + .5 * s9
                IF (z <= 5.e-15) THEN
                  a = a + .5
                  IF (z <= 1.e-18) THEN
                    a = a + .25
                    IF (z <= 1.e-20) THEN
                      a = a + .25
                      IF (z <= 1.e-25) THEN
                        a = a + 1.
                        IF (z <= 1.e-30) THEN
                          a = a + 2.
                        END IF
                      END IF
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF
END IF

IF (s2 < 5.e-2 .AND. h3 > z2) THEN
  t9 = r8 * u * EXP(-h2*h2)
  IF (t9 < 5.e-2*z) a = a + 0.5
END IF
!---------------------------------------
!     START INTEGRATION PROCEDURE
!---------------------------------------
s0 = s1
s9 = s2
z8 = s2
g2 = k8
g3 = h8
z9 = s1
j8 = 0
!-------------------------------------------------
!     DETERMINE INTERVAL OF INTEGRATION, (A,RT2).
!           E3 = (RT2-A)/2, D1 = (RT2+A)/2.
!-------------------------------------------------
30 z = g2 + a * z8
h3 = g2 - a * z8
h5 = 0.0
t3 = -1.0
a1 = (g3-a*z9) / r
IF (ABS(a1-0.5) < 0.5) THEN
  a2 = MAX(((1.0-g3/r)+a*z9/r)*(1.0+a1),0.0)
  sa = SQRT(a2)
  IF (sa >= h3/r) THEN
    t3 = a1 / SQRT(1.0+sa)
  END IF
END IF
!----------------------------------------
!            T9 = 1.0 - D1
!----------------------------------------
IF (h3 <= 0.0) THEN
!----------------------------------------
!            H3 <= 0.0
!----------------------------------------
  IF (t3 >= 0.0) THEN
    d2 = 0.5 * (1.0+t3)
    e4 = 0.25 * sa / d2
    t5 = e4
  END IF
  IF (z >= r) THEN
    e3 = 0.5
    d1 = e3
    t9 = d1
  ELSE
    d1 = 0.5 * (1.0+SQRT(1.0-z/r))
    e3 = 0.25 * z / (r*d1)
    t9 = e3
  END IF
  IF (t3 <= d1-e3) GO TO 40
ELSE
!----------------------------------------
!            H3 > 0.0
!----------------------------------------
  h5 = 1.0
  q = h3 / r
  IF (q > 1.0) THEN
    e3 = 0.5
    d1 = e3
    t9 = d1
    GO TO 40
  END IF
  IF (t3 >= 0.0) THEN
    e4 = MAX((1.0-g2/r)+a*z8/r,0.0)
    e4 = SQRT(e4)
    d2 = 0.5 * (e4+t3)
    t5 = 0.5 * (q/(1.0+e4)+sa/(1.0+t3))
    e4 = 0.25 * (sa-q) / d2
  END IF
  e3 = MAX((1.0-g2/r)+a*z8/r,0.0)
  IF (z >= r) THEN
    e3 = 0.5 * SQRT(e3)
    d1 = e3
    t9 = 0.25 * (3.0+h3/r) / (1.0+e3)
  ELSE
    t9 = SQRT(e3)
    t2 = SQRT(1.0-z/r)
    d1 = 0.5 * (t9+t2)
    e3 = a * z8 / (r*2.0*d1)
    t9 = 0.5 * ((h3/r)/(1.0+t9)+(z/r)/(1.0+t2))
  END IF
  IF (t3 <= d1-e3) GO TO 40
END IF
e3 = e4
d1 = d2
t9 = t5
40 IF (j8 == 0) THEN
  j8 = 1
  f = e3
  t = d1
  t1 = t9
  h6 = h5
  z8 = s1
  g2 = h8
  g3 = k8
  z9 = s2
  GO TO 30
END IF
!-----------------------------------------------------------------------
!     DETERMINE IN WHICH ORDER THE X AMD Y INTEGRATIONS ARE CARRIED OUT.
!-----------------------------------------------------------------------
IF (s2 <= 2.e-2.OR.h8+k8 >= 2.e2) THEN
  IF (ABS(e3-f) <= 0.4*f) THEN
    IF (d1 < t) THEN
      GO TO 90
    ELSE IF (d1 > t) THEN
      GO TO 80
    ELSE
      GO TO 50
    END IF
  END IF
END IF
IF (e3 < 2.e4*eps) GO TO 80
IF (f < 2.e4*eps) GO TO 90
IF (MAX(h8/s1,k8/s2) >= 2.0) THEN
  IF (s2 >= 1.e-5) THEN
    IF (s2 >= 5.e-4) GO TO 50
    IF (s1 /= s2) GO TO 80
  END IF
END IF

IF (MIN(e3,f) < 2.5E-2*sqeps) GO TO 60
50 IF (e3 < f) THEN
  GO TO 90
ELSE IF (e3 > f) THEN
  GO TO 80
ELSE
  GO TO 60
END IF

60 IF (e3 == f .AND. s0 >= s9) GO TO 90
IF (e3 > f) GO TO 90
80 e3 = f
d1 = t
t9 = t1
s9 = s1
s0 = s2
z8 = h8
h8 = k8
k8 = z8
h5 = h6

90 z2 = r / (rt2*s9)
r8 = r / (rt2*s0)
h2 = h8 / (rt2*s0)
h3 = k8 / (rt2*s9)
n1 = 16
iz = 0
iz3 = 0
iy = 0
p = 0.0
t = h2 - r8 + r8 * (d1-e3*v(16)) ** 2
IF (t > 0.0 .AND. t*t > -exparg(1)) RETURN

q1 = rpi * e3 * r8
g3 = 0.0
nt = 2 * n1 + 1
z = (.5+d1) + .5
DO ii = 1, nt
  i = ii - (n1+1)
  IF (i /= 0) THEN
    j = ABS(i)
    q = e3 * SIGN(1,i) * v(j)
    t = q + d1
    h6 = z + q
    q = t9 - q
    f = h6 * q
    t1 = r8 * f
    t2 = (h2-t1) * (h2-t1)
    IF (h2-r8 >= 0..OR.t < eps) t2 = ((h2-r8)+r8*t*t) ** 2
    t4 = EXP(-t2)
    IF (h2 == 0.0) THEN
      t4 = t4 + t4
    ELSE
      IF (h5 == 0.0) THEN
        t2 = 4.0 * h2 * t1
        IF (t2 <= c8) THEN
          t4 = t4 * (1.0+EXP(-t2))
        END IF
      END IF
    END IF
    IF (iz == 0) THEN
      g2 = SQRT(1.0+f)
      z1 = z2 * t * g2
      x1 = h3 - z1
      IF (x1 <= -a) THEN
        iz = 1
        x5 = 2.0
      ELSE
        IF (ABS(x1) <= 1.e-2*z1) THEN
          iy = 1
          x1 = ((k8-r)+r*f*f/(1.0+t*g2)) / (rt2*s9)
        END IF
        IF (iz3 == 0) THEN
          sa = h3 + z1
          IF (sa <= a3) THEN
            IF (iy /= 0) THEN
              IF (x1 <= rt2) THEN
                x5 = erf(sa) - erf(x1)
                GO TO 100
              END IF
              x5 = erfc(x1) - erfc(sa)
              GO TO 100
            END IF
            x5 = aerf(h3,z1)
            GO TO 100
          END IF
          iz3 = 1
        END IF
        x5 = erfc(x1)
      END IF
    END IF
    100 g3 = g3 + x5 * t4 * t * w(j)
  END IF
END DO
p = q1 * g3
IF (p > yy) p = yy
IF (p > (1.0-MIN(1.e6*eps,1.e-5))) p = 1.0
IF (p < 0.0) p = 0.0
RETURN
END SUBROUTINE pkill



FUNCTION csevl(x, a, n) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATE THE N TERM CHEBYSHEV SERIES A AT X.
!          ONLY HALF OF THE FIRST COEFFICIENT IS USED.
!-----------------------------------------------------------------------
REAL, INTENT(IN)    :: a(:), x
INTEGER, INTENT(IN) :: n
REAL                :: fn_val

! Local variables
REAL    :: x2, s0, s1, s2
INTEGER :: i, k

IF (n <= 1) THEN
  fn_val = 0.5 * a(1)
  RETURN
END IF

x2 = x + x
s0 = a(n)
s1 = 0.0
DO i = 2, n
  s2 = s1
  s1 = s0
  k = n - i + 1
  s0 = x2 * s1 - s2 + a(k)
END DO
fn_val = 0.5 * (s0-s2)
RETURN
END FUNCTION csevl


FUNCTION dcsevl(x, a, n) RESULT(fn_val)
!-----------------------------------------------------------------------
!          EVALUATE THE N TERM CHEBYSHEV SERIES A AT X.
!          ONLY HALF OF THE FIRST COEFFICIENT IS USED.
!-----------------------------------------------------------------------
REAL (dp), INTENT(IN) :: a(:), x
INTEGER, INTENT(IN)   :: n
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: x2, s0, s1, s2
INTEGER   :: i, k

IF (n <= 1) THEN
  fn_val = 0.5_dp * a(1)
  RETURN
END IF

x2 = x + x
s0 = a(n)
s1 = 0._dp
DO i = 2, n
  s2 = s1
  s1 = s0
  k = n - i + 1
  s0 = x2 * s1 - s2 + a(k)
END DO
fn_val = 0.5_dp * (s0-s2)
RETURN
END FUNCTION dcsevl

END MODULE NSWC_spec_func
