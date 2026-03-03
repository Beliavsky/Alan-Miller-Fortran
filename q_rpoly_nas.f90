MODULE Solve_Real_Poly
! CACM Algorithm 493 by Jenkins & Traub

! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-08  Time: 16:02:50

! Latest revision - 10 February 2000
! amiller @ bigpond.net.au

USE quadruple_precision
IMPLICIT NONE

! COMMON /global/ p, qp, k, qk, svk, sr, si, u, v, a, b, c, d, a1,  &
!    a2, a3, a6, a7, e, f, g, h, szr, szi, lzr, lzi, eta, are, mre, n, nn

TYPE (quad), ALLOCATABLE, SAVE  :: p(:), qp(:), k(:), qk(:), svk(:)
TYPE (quad), SAVE               :: sr, si, u, v, a, b, c, d, a1, a2, a3, a6,  &
                                   a7, ee, f, g, h, szr, szi, lzr, lzi
REAL (r10), SAVE                :: eta, are, mre
INTEGER, SAVE                   :: n, nn

PRIVATE
PUBLIC  :: rpoly


CONTAINS


SUBROUTINE rpoly(op, degree, zeror, zeroi, fail)

! FINDS THE ZEROS OF A REAL POLYNOMIAL
! OP  - TYPE (quad) VECTOR OF COEFFICIENTS IN ORDER OF DECREASING POWERS.
! DEGREE   - INTEGER DEGREE OF POLYNOMIAL.
! ZEROR, ZEROI - OUTPUT TYPE (quad) VECTORS OF
!                REAL AND IMAGINARY PARTS OF THE ZEROS.
! FAIL  - OUTPUT LOGICAL PARAMETER, TRUE ONLY IF THE LEADING COEFFICIENT
!         IS ZERO, OR IF RPOLY HAS FOUND FEWER THAN DEGREE ZEROS.
!         IN THE LATTER CASE DEGREE IS RESET TO THE NUMBER OF ZEROS FOUND.
! TO CHANGE THE SIZE OF POLYNOMIALS WHICH CAN BE SOLVED, RESET THE DIMENSIONS
! OF THE ARRAYS IN THE COMMON AREA AND IN THE FOLLOWING DECLARATIONS.
! THE SUBROUTINE USES SINGLE PRECISION CALCULATIONS FOR SCALING, BOUNDS AND
! ERROR CALCULATIONS.  ALL CALCULATIONS FOR THE ITERATIONS ARE DONE IN DOUBLE
! PRECISION.

TYPE (quad), INTENT(IN)   :: op(:)
INTEGER, INTENT(IN OUT)   :: degree
TYPE (quad), INTENT(OUT)  :: zeror(:)
TYPE (quad), INTENT(OUT)  :: zeroi(:)
LOGICAL, INTENT(OUT)      :: fail

TYPE (quad), ALLOCATABLE  :: temp(:)
REAL (r10), ALLOCATABLE   :: pt(:)

TYPE (quad)  :: t, aa, bb, cc, factor
REAL (r10)   :: lo, MAX, MIN, xx, yy, cosr, sinr, xxx, x, sc, bnd,  &
                xm, ff, df, dx, infin, smalno, base
INTEGER      :: cnt, nz, i, j, jj, l, nm1
LOGICAL      :: zerok

! THE FOLLOWING STATEMENTS SET MACHINE CONSTANTS USED IN VARIOUS PARTS OF THE
! PROGRAM.  THE MEANING OF THE FOUR CONSTANTS ARE...
! ETA     THE MAXIMUM RELATIVE REPRESENTATION ERROR WHICH CAN BE DESCRIBED AS
!         THE SMALLEST POSITIVE FLOATING POINT NUMBER SUCH THAT 1._r10+ETA > 1.
! INFINY  THE LARGEST FLOATING-POINT NUMBER.
! SMALNO  THE SMALLEST POSITIVE FLOATING-POINT NUMBER IF THE EXPONENT RANGE
!         DIFFERS IN SINGLE AND TYPE (quad) THEN SMALNO AND INFIN
!         SHOULD INDICATE THE SMALLER RANGE.
! BASE    THE BASE OF THE FLOATING-POINT NUMBER SYSTEM USED.

base = RADIX(0.0)
eta = EPSILON(1.0_r10) * 1.0E-18
infin = HUGE(0.0)
smalno = TINY(0.0)

! ARE AND MRE REFER TO THE UNIT ERROR IN + AND * RESPECTIVELY.
! THEY ARE ASSUMED TO BE THE SAME AS ETA.
are = eta
mre = eta
lo = smalno / eta

! INITIALIZATION OF CONSTANTS FOR SHIFT ROTATION
xx = .70710678
yy = -xx
cosr = -.069756474
sinr = .99756405
fail = .false.
n = degree
nn = n + 1

! ALGORITHM FAILS IF THE LEADING COEFFICIENT IS ZERO.
IF (op(1)%hi == 0._r10) THEN
  fail = .true.
  degree = 0
  RETURN
END IF

! REMOVE THE ZEROS AT THE ORIGIN IF ANY
10 IF (op(nn)%hi == 0.0_r10) THEN
  j = degree - n + 1
  zeror(j) = 0._r10
  zeroi(j) = 0._r10
  nn = nn - 1
  n = n - 1
  GO TO 10
END IF

! Allocate various arrays

IF (ALLOCATED(p)) DEALLOCATE(p, qp, k, qk, svk, temp, pt)
ALLOCATE( p(nn), qp(nn), k(nn), qk(nn), svk(nn), temp(nn), pt(nn) )

! MAKE A COPY OF THE COEFFICIENTS
p(1:nn) = op(1:nn)

! START THE ALGORITHM FOR ONE ZERO
30 IF (n <= 2) THEN
  IF (n < 1) RETURN

! CALCULATE THE FINAL ZERO OR PAIR OF ZEROS
  IF (n /= 2) THEN
    zeror(degree) = -p(2) / p(1)
    zeroi(degree) = 0.0_r10
    RETURN
  END IF
  CALL quadratic(p(1), p(2), p(3), zeror(degree-1), zeroi(degree-1),  &
                 zeror(degree), zeroi(degree))
  RETURN
END IF

! FIND LARGEST AND SMALLEST MODULI OF COEFFICIENTS.
MAX = 0.
MIN = infin
DO  i = 1, nn
  x = ABS(p(i))
  IF (x > MAX) MAX = x
  IF (x /= 0. .AND. x < MIN) MIN = x
END DO

! SCALE IF THERE ARE LARGE OR VERY SMALL COEFFICIENTS COMPUTES A SCALE FACTOR
! TO MULTIPLY THE COEFFICIENTS OF THE POLYNOMIAL.  THE SCALING IS DONE TO
! AVOID OVERFLOW AND TO AVOID UNDETECTED UNDERFLOW INTERFERING WITH THE
! CONVERGENCE CRITERION.  THE FACTOR IS A POWER OF THE BASE.

sc = lo / MIN
IF (sc <= 1.0) THEN
  IF (MAX < 10.) GO TO 60
  IF (sc == 0.) sc = smalno
ELSE
  IF (infin/sc < MAX) GO TO 60
END IF
l = LOG(sc) / LOG(base) + .5
factor = (base) ** l
IF (factor%hi /= 1._r10) THEN
  DO i = 1, nn
    p(i) = factor * p(i)
  END DO
END IF

! COMPUTE LOWER BOUND ON MODULI OF ZEROS.
60 DO i = 1, nn
  pt(i) = ABS(p(i))
END DO
pt(nn) = -pt(nn)

! COMPUTE UPPER ESTIMATE OF BOUND
x = EXP( (LOG(-pt(nn)) - LOG(pt(1)) ) / n )
IF (pt(n) /= 0.) THEN

! IF NEWTON STEP AT THE ORIGIN IS BETTER, USE IT.
  xm = -pt(nn) / pt(n)
  IF (xm < x) x = xm
END IF

! CHOP THE INTERVAL (0,X) UNTIL FF <= 0
80 xm = x * .1
ff = pt(1)
DO  i = 2, nn
  ff = ff * xm + pt(i)
END DO
IF (ff > 0.) THEN
  x = xm
  GO TO 80
END IF
dx = x

! DO NEWTON ITERATION UNTIL X CONVERGES TO TWO DECIMAL PLACES
100 IF (ABS(dx/x) > .005) THEN
  ff = pt(1)
  df = ff
  DO  i = 2, n
    ff = ff * x + pt(i)
    df = df * x + ff
  END DO
  ff = ff * x + pt(nn)
  dx = ff / df
  x = x - dx
  GO TO 100
END IF
bnd = x

! COMPUTE THE DERIVATIVE AS THE INTIAL K POLYNOMIAL
! AND DO 5 STEPS WITH NO SHIFT
nm1 = n - 1
DO  i = 2, n
  k(i) = (nn-i) * p(i) / REAL(n)
END DO
k(1) = p(1)
aa = p(nn)
bb = p(n)
zerok = k(n)%hi == 0._r10
DO  jj = 1, 5
  cc = k(n)
  IF (.NOT.zerok) THEN
! USE SCALED FORM OF RECURRENCE IF VALUE OF K AT 0 IS NONZERO
    t = -aa / cc
    DO  i = 1, nm1
      j = nn - i
      k(j) = t * k(j-1) + p(j)
    END DO
    k(1) = p(1)
    zerok = ABS(k(n)) <= ABS(bb) * eta * 10.
  ELSE
! USE UNSCALED FORM OF RECURRENCE
    DO  i = 1, nm1
      j = nn - i
      k(j) = k(j-1)
    END DO
    k(1) = 0._r10
    zerok = k(n)%hi == 0._r10
  END IF
END DO

! SAVE K FOR RESTARTS WITH NEW SHIFTS
temp(1:n) = k(1:n)

! LOOP TO SELECT THE QUADRATIC  CORRESPONDING TO EACH NEW SHIFT
DO  cnt = 1, 20
! QUADRATIC CORRESPONDS TO A DOUBLE SHIFT TO A NON-REAL POINT AND ITS COMPLEX
! CONJUGATE.  THE POINT HAS MODULUS BND AND AMPLITUDE ROTATED BY 94 DEGREES
! FROM THE PREVIOUS SHIFT.
  xxx = cosr * xx - sinr * yy
  yy = sinr * xx + cosr * yy
  xx = xxx
  sr = bnd * xx
  si = bnd * yy
  u = -2.0_r10 * sr
  v = bnd

! SECOND STAGE CALCULATION, FIXED QUADRATIC
  CALL fxshfr(20*cnt,nz)
  IF (nz /= 0) THEN

! THE SECOND STAGE JUMPS DIRECTLY TO ONE OF THE THIRD STAGE ITERATIONS AND
! RETURNS HERE IF SUCCESSFUL.
! DEFLATE THE POLYNOMIAL, STORE THE ZERO OR ZEROS AND RETURN TO THE MAIN
! ALGORITHM.
    j = degree - n + 1
    zeror(j) = szr
    zeroi(j) = szi
    nn = nn - nz
    n = nn - 1
    p(1:nn) = qp(1:nn)
    IF (nz == 1) GO TO 30
    zeror(j+1) = lzr
    zeroi(j+1) = lzi
    GO TO 30
  END IF

! IF THE ITERATION IS UNSUCCESSFUL ANOTHER QUADRATIC
! IS CHOSEN AFTER RESTORING K
  k(1:n) = temp(1:n)
END DO

! RETURN WITH FAILURE IF NO CONVERGENCE WITH 20 SHIFTS
fail = .true.
degree = degree - n
RETURN
END SUBROUTINE rpoly



SUBROUTINE fxshfr(l2, nz)
! COMPUTES UP TO  L2  FIXED SHIFT K-POLYNOMIALS, TESTING FOR CONVERGENCE IN
! THE LINEAR OR QUADRATIC CASE.  INITIATES ONE OF THE VARIABLE SHIFT
! ITERATIONS AND RETURNS WITH THE NUMBER OF ZEROS FOUND.
! L2 - LIMIT OF FIXED SHIFT STEPS
! NZ - NUMBER OF ZEROS FOUND

INTEGER, INTENT(IN)   :: l2
INTEGER, INTENT(OUT)  :: nz

TYPE (quad)  :: svu, svv, ui, vi, s
REAL (r10)   :: betas, betav, oss, ovv, ss, vv, ts, tv, ots, otv, tvv, tss
INTEGER      :: TYPE, j, iflag
LOGICAL      :: vpass, spass, vtry, stry

nz = 0
betav = .25
betas = .25
oss = sr
ovv = v

! EVALUATE POLYNOMIAL BY SYNTHETIC DIVISION
CALL quadsd(nn, u, v, p, qp, a, b)
CALL calcsc(TYPE)
DO  j = 1, l2
! CALCULATE NEXT K POLYNOMIAL AND ESTIMATE V
  CALL nextk(TYPE)
  CALL calcsc(TYPE)
  CALL newest(TYPE,ui,vi)
  vv = vi

! ESTIMATE S
  ss = 0.
  IF (k(n)%hi /= 0._r10) ss = -p(nn) / k(n)
  tv = 1.
  ts = 1.
  IF (j /= 1 .AND. TYPE /= 3) THEN

! COMPUTE RELATIVE MEASURES OF CONVERGENCE OF S AND V SEQUENCES
    IF (vv /= 0.) tv = ABS((vv-ovv)/vv)
    IF (ss /= 0.) ts = ABS((ss-oss)/ss)

! IF DECREASING, MULTIPLY TWO MOST RECENT CONVERGENCE MEASURES
    tvv = 1.
    IF (tv < otv) tvv = tv * otv
    tss = 1.
    IF (ts < ots) tss = ts * ots

! COMPARE WITH CONVERGENCE CRITERIA
    vpass = tvv < betav
    spass = tss < betas
    IF ((spass.OR.vpass)) THEN

! AT LEAST ONE SEQUENCE HAS PASSED THE CONVERGENCE TEST.
! STORE VARIABLES BEFORE ITERATING
      svu = u
      svv = v
      svk(1:n) = k(1:n)
      s = ss

! CHOOSE ITERATION ACCORDING TO THE FASTEST CONVERGING SEQUENCE
      vtry = .false.
      stry = .false.
      IF (spass .AND. ((.NOT.vpass) .OR. tss < tvv)) GO TO 40
      20 CALL quadit(ui, vi, nz)
      IF (nz > 0) RETURN

! QUADRATIC ITERATION HAS FAILED.  FLAG THAT IT HAS
! BEEN TRIED AND DECREASE THE CONVERGENCE CRITERION.
      vtry = .true.
      betav = betav * .25

! TRY LINEAR ITERATION IF IT HAS NOT BEEN TRIED AND
! THE S SEQUENCE IS CONVERGING
      IF (stry .OR. (.NOT.spass)) GO TO 50
      k(1:n) = svk(1:n)
      40 CALL realit(s, nz, iflag)
      IF (nz > 0) RETURN

! LINEAR ITERATION HAS FAILED.  FLAG THAT IT HAS BEEN
! TRIED AND DECREASE THE CONVERGENCE CRITERION
      stry = .true.
      betas = betas * .25
      IF (iflag /= 0) THEN

! IF LINEAR ITERATION SIGNALS AN ALMOST DOUBLE REAL
! ZERO ATTEMPT QUADRATIC INTERATION
        ui = -(s+s)
        vi = s * s
        GO TO 20
      END IF

! RESTORE VARIABLES
      50 u = svu
      v = svv
      k(1:n) = svk(1:n)

! TRY QUADRATIC ITERATION IF IT HAS NOT BEEN TRIED
! AND THE V SEQUENCE IS CONVERGING
      IF (vpass .AND. (.NOT.vtry)) GO TO 20

! RECOMPUTE QP AND SCALAR VALUES TO CONTINUE THE SECOND STAGE
      CALL quadsd(nn, u, v, p, qp, a, b)
      CALL calcsc(TYPE)
    END IF
  END IF
  ovv = vv
  oss = ss
  otv = tv
  ots = ts
END DO
RETURN
END SUBROUTINE fxshfr


SUBROUTINE quadit(uu, vv, nz)
! VARIABLE-SHIFT K-POLYNOMIAL ITERATION FOR A QUADRATIC FACTOR, CONVERGES
! ONLY IF THE ZEROS ARE EQUIMODULAR OR NEARLY SO.
! UU,VV - COEFFICIENTS OF STARTING QUADRATIC
! NZ - NUMBER OF ZERO FOUND

TYPE (quad), INTENT(IN)  :: uu
TYPE (quad), INTENT(IN)  :: vv
INTEGER, INTENT(OUT)     :: nz

TYPE (quad)  :: ui, vi
REAL (r10)   :: mp, omp, ee, relstp, t, zm
INTEGER      :: TYPE, i, j
LOGICAL      :: tried

nz = 0
tried = .false.
u = uu
v = vv
j = 0

! MAIN LOOP
10 CALL quadratic(quad(1._r10,0.0_r10), u, v, szr, szi, lzr, lzi)

! RETURN IF ROOTS OF THE QUADRATIC ARE REAL AND NOT
! CLOSE TO MULTIPLE OR NEARLY EQUAL AND  OF OPPOSITE SIGN
IF (ABS(ABS(szr)-ABS(lzr)) > .01_r10*ABS(lzr)) RETURN

! EVALUATE POLYNOMIAL BY QUADRATIC SYNTHETIC DIVISION
CALL quadsd(nn, u, v, p, qp, a, b)
mp = ABS(a-szr*b) + ABS(szi*b)

! COMPUTE A RIGOROUS  BOUND ON THE ROUNDING ERROR IN EVALUTING P
zm = SQRT(ABS(v))
ee = 2. * ABS(qp(1))
t = -szr * b
DO  i = 2, n
  ee = ee * zm + ABS(qp(i))
END DO
ee = ee * zm + ABS(a + t)
ee = (5.*mre + 4.*are) * ee - (5.*mre + 2.*are) * (ABS(a) + t) +  &
     ABS(b)*zm + 2. * are * ABS(t)

! ITERATION HAS CONVERGED SUFFICIENTLY IF THE
! POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND
IF (mp <= 20.*ee) THEN
  nz = 2
  RETURN
END IF
j = j + 1

! STOP ITERATION AFTER 20 STEPS
IF (j > 20) RETURN
IF (j >= 2) THEN
  IF (.NOT.(relstp > .01 .OR. mp < omp .OR. tried)) THEN

! A CLUSTER APPEARS TO BE STALLING THE CONVERGENCE.
! FIVE FIXED SHIFT STEPS ARE TAKEN WITH A U,V CLOSE TO THE CLUSTER
    IF (relstp < eta) relstp = eta
    relstp = SQRT(relstp)
    u = u - u * relstp
    v = v + v * relstp
    CALL quadsd(nn, u, v, p, qp, a, b)
    DO  i = 1, 5
      CALL calcsc(TYPE)
      CALL nextk(TYPE)
    END DO
    tried = .true.
    j = 0
  END IF
END IF
omp = mp

! CALCULATE NEXT K POLYNOMIAL AND NEW U AND V
CALL calcsc(TYPE)
CALL nextk(TYPE)
CALL calcsc(TYPE)
CALL newest(TYPE, ui, vi)

! IF VI IS ZERO THE ITERATION IS NOT CONVERGING
IF (vi%hi == 0._r10) RETURN
relstp = ABS((vi-v)/vi)
u = ui
v = vi
GO TO 10
END SUBROUTINE quadit


SUBROUTINE realit(sss, nz, iflag)
! VARIABLE-SHIFT H POLYNOMIAL ITERATION FOR A REAL ZERO.
! SSS   - STARTING ITERATE
! NZ    - NUMBER OF ZERO FOUND
! IFLAG - FLAG TO INDICATE A PAIR OF ZEROS NEAR REAL AXIS.

TYPE (quad), INTENT(IN OUT)  :: sss
INTEGER, INTENT(OUT)         :: nz
INTEGER, INTENT(OUT)         :: iflag

TYPE (quad)  :: pv, kv, t, s
REAL (r10)   :: ms, mp, omp, ee
INTEGER      :: i, j

nz = 0
s = sss
iflag = 0
j = 0

! MAIN LOOP
10 pv = p(1)

! EVALUATE P AT S
qp(1) = pv
DO  i = 2, nn
  pv = pv * s + p(i)
  qp(i) = pv
END DO
mp = ABS(pv)

! COMPUTE A RIGOROUS BOUND ON THE ERROR IN EVALUATING P
ms = ABS(s)
ee = (mre/(are+mre)) * ABS(qp(1))
DO  i = 2, nn
  ee = ee * ms + ABS(qp(i))
END DO

! ITERATION HAS CONVERGED SUFFICIENTLY IF THE
! POLYNOMIAL VALUE IS LESS THAN 20 TIMES THIS BOUND
IF (mp <= 20.*((are+mre)*ee - mre*mp)) THEN
  nz = 1
  szr = s
  szi = 0._r10
  RETURN
END IF
j = j + 1

! STOP ITERATION AFTER 10 STEPS
IF (j > 10) RETURN
IF (j >= 2) THEN
  IF (ABS(t) <= .001*ABS(s-t) .AND. mp > omp) THEN

! A CLUSTER OF ZEROS NEAR THE REAL AXIS HAS BEEN ENCOUNTERED RETURN WITH
! IFLAG SET TO INITIATE A QUADRATIC ITERATION
    iflag = 1
    sss = s
    RETURN
  END IF
END IF

! RETURN IF THE POLYNOMIAL VALUE HAS INCREASED SIGNIFICANTLY
omp = mp

! COMPUTE T, THE NEXT POLYNOMIAL, AND THE NEW ITERATE
kv = k(1)
qk(1) = kv
DO  i = 2, n
  kv = kv * s + k(i)
  qk(i) = kv
END DO
IF (ABS(kv) > ABS(k(n))*10.*eta) THEN

! USE THE SCALED FORM OF THE RECURRENCE IF THE VALUE OF K AT S IS NONZERO
  t = -pv / kv
  k(1) = qp(1)
  DO  i = 2, n
    k(i) = t * qk(i-1) + qp(i)
  END DO
ELSE
! USE UNSCALED FORM
  k(1) = 0.0_r10
  DO  i = 2, n
    k(i) = qk(i-1)
  END DO
END IF
kv = k(1)
DO  i = 2, n
  kv = kv * s + k(i)
END DO
t = 0._r10
IF (ABS(kv) > ABS(k(n))*10.*eta) t = -pv / kv
s = s + t
GO TO 10
END SUBROUTINE realit


SUBROUTINE calcsc(TYPE)
! THIS ROUTINE CALCULATES SCALAR QUANTITIES USED TO COMPUTE THE NEXT K
! POLYNOMIAL AND NEW ESTIMATES OF THE QUADRATIC COEFFICIENTS.
! TYPE - INTEGER VARIABLE SET HERE INDICATING HOW THE CALCULATIONS ARE
! NORMALIZED TO AVOID OVERFLOW

INTEGER, INTENT(OUT)  :: TYPE

! SYNTHETIC DIVISION OF K BY THE QUADRATIC 1,U,V
CALL quadsd(n, u, v, k, qk, c, d)
IF (ABS(c) <= ABS(k(n))*100.*eta) THEN
  IF (ABS(d) <= ABS(k(n-1))*100.*eta) THEN
    TYPE = 3

! TYPE=3 INDICATES THE QUADRATIC IS ALMOST A FACTOR OF K
    RETURN
  END IF
END IF
IF (ABS(d) >= ABS(c)) THEN
  TYPE = 2

! TYPE=2 INDICATES THAT ALL FORMULAS ARE DIVIDED BY D
  ee = a / d
  f = c / d
  g = u * b
  h = v * b
  a3 = (a+g) * ee + h * (b/d)
  a1 = b * f - a
  a7 = (f+u) * a + h
  RETURN
END IF
TYPE = 1

! TYPE=1 INDICATES THAT ALL FORMULAS ARE DIVIDED BY C
ee = a / c
f = d / c
g = u * ee
h = v * b
a3 = a * ee + (h/c+g) * b
a1 = b - a * (d/c)
a7 = a + g * d + h * f
RETURN
END SUBROUTINE calcsc


SUBROUTINE nextk(TYPE)
! COMPUTES THE NEXT K POLYNOMIALS USING SCALARS COMPUTED IN CALCSC

INTEGER, INTENT(IN)  :: TYPE

TYPE (quad) :: temp
INTEGER     :: i

IF (TYPE /= 3) THEN
  temp = a
  IF (TYPE == 1) temp = b
  IF (ABS(a1) <= ABS(temp)*eta*10.) THEN

! IF A1 IS NEARLY ZERO THEN USE A SPECIAL FORM OF THE RECURRENCE
    k(1) = 0._r10
    k(2) = -a7 * qp(1)
    DO  i = 3, n
      k(i) = a3 * qk(i-2) - a7 * qp(i-1)
    END DO
    RETURN
  END IF

! USE SCALED FORM OF THE RECURRENCE
  a7 = a7 / a1
  a3 = a3 / a1
  k(1) = qp(1)
  k(2) = qp(2) - a7 * qp(1)
  DO  i = 3, n
    k(i) = a3 * qk(i-2) - a7 * qp(i-1) + qp(i)
  END DO
  RETURN
END IF

! USE UNSCALED FORM OF THE RECURRENCE IF TYPE IS 3
k(1) = 0._r10
k(2) = 0._r10
DO  i = 3, n
  k(i) = qk(i-2)
END DO
RETURN
END SUBROUTINE nextk


SUBROUTINE newest(TYPE, uu, vv)
! COMPUTE NEW ESTIMATES OF THE QUADRATIC COEFFICIENTS
! USING THE SCALARS COMPUTED IN CALCSC.

INTEGER, INTENT(IN)       :: TYPE
TYPE (quad), INTENT(OUT)  :: uu
TYPE (quad), INTENT(OUT)  :: vv

TYPE (quad) :: a4, a5, b1, b2, c1, c2, c3, c4, temp

! USE FORMULAS APPROPRIATE TO SETTING OF TYPE.
IF (TYPE /= 3) THEN
  IF (TYPE /= 2) THEN
    a4 = a + u * b + h * f
    a5 = c + (u+v*f) * d
  ELSE
    a4 = (a+g) * f + h
    a5 = (f+u) * c + v * d
  END IF

! EVALUATE NEW QUADRATIC COEFFICIENTS.
  b1 = -k(n) / p(nn)
  b2 = -(k(n-1) + b1*p(n)) / p(nn)
  c1 = v * b2 * a1
  c2 = b1 * a7
  c3 = b1 * b1 * a3
  c4 = c1 - c2 - c3
  temp = a5 + b1 * a4 - c4
  IF (temp%hi /= 0._r10) THEN
    uu = u - (u*(c3 + c2) + v*(b1*a1 + b2*a7)) / temp
    vv = v * (1. + c4/temp)
    RETURN
  END IF
END IF

! IF TYPE=3 THE QUADRATIC IS ZEROED
uu = 0._r10
vv = 0._r10
RETURN
END SUBROUTINE newest


SUBROUTINE quadsd(nn, u, v, p, q, a, b)
! DIVIDES P BY THE QUADRATIC  1,U,V  PLACING THE
! QUOTIENT IN Q AND THE REMAINDER IN A,B

INTEGER, INTENT(IN)       :: nn
TYPE (quad), INTENT(IN)   :: u
TYPE (quad), INTENT(IN)   :: v
TYPE (quad), INTENT(IN)   :: p(:)
TYPE (quad), INTENT(OUT)  :: q(:)
TYPE (quad), INTENT(OUT)  :: a
TYPE (quad), INTENT(OUT)  :: b

TYPE (quad) :: c
INTEGER     :: i

b = p(1)
q(1) = b
a = p(2) - u * b
q(2) = a
DO  i = 3, nn
  c = p(i) - u * a - v * b
  q(i) = c
  b = a
  a = c
END DO
RETURN
END SUBROUTINE quadsd


SUBROUTINE quadratic(a, b1, c, sr, si, lr, li)
! CALCULATE THE ZEROS OF THE QUADRATIC A*Z**2 + B1*Z + C.
! THE QUADRATIC FORMULA, MODIFIED TO AVOID OVERFLOW, IS USED TO FIND THE
! LARGER ZERO IF THE ZEROS ARE REAL AND BOTH ZEROS ARE COMPLEX.
! THE SMALLER REAL ZERO IS FOUND DIRECTLY FROM THE PRODUCT OF THE ZEROS C/A.

TYPE (quad), INTENT(IN)   :: a
TYPE (quad), INTENT(IN)   :: b1
TYPE (quad), INTENT(IN)   :: c
TYPE (quad), INTENT(OUT)  :: sr
TYPE (quad), INTENT(OUT)  :: si
TYPE (quad), INTENT(OUT)  :: lr
TYPE (quad), INTENT(OUT)  :: li

TYPE (quad) :: b, d, ee

IF (a%hi /= 0._r10) GO TO 20
sr = 0._r10
IF (b1%hi /= 0._r10) sr = -c / b1
lr = 0._r10
10 si = 0._r10
li = 0._r10
RETURN

20 IF (c%hi == 0._r10) THEN
  sr = 0._r10
  lr = -b1 / a
  GO TO 10
END IF

! COMPUTE DISCRIMINANT AVOIDING OVERFLOW
b = b1 / 2._r10
IF (ABS(b) >= ABS(c)) THEN
  ee = 1._r10 - (a/b) * (c/b)
  d = SQRT(ABS(ee)) * ABS(b)
ELSE
  ee = a
  IF (c%hi < 0._r10) ee = -a
  ee = b * (b/ABS(c)) - ee
  d = SQRT(ABS(ee)) * SQRT(ABS(c))
END IF

IF (ee%hi >= 0._r10) THEN

! REAL ZEROS
  IF (b%hi >= 0._r10) d = -d
  lr = (-b+d) / a
  sr = 0._r10
  IF (lr%hi /= 0._r10) sr = (c/lr) / a
  GO TO 10
END IF

! COMPLEX CONJUGATE ZEROS
sr = -b / a
lr = sr
si = ABS(d/a)
li = -si
RETURN
END SUBROUTINE quadratic

END MODULE Solve_Real_Poly



PROGRAM test_rpoly
USE quadruple_precision
USE Solve_Real_Poly
IMPLICIT NONE

TYPE (quad)  :: p(50), zr(50), zi(50)
INTEGER      :: degree, i
LOGICAL      :: fail

WRITE(*, 5000)
degree = 10
p(1) = 1.
p(2) = -55.
p(3) = 1320.
p(4) = -18150.
p(5) = 157773.
p(6) = -902055.
p(7) = 3416930.
p(8) = -8409500.
p(9) = 12753576.
p(10) = -10628640.
p(11) = 3628800.

CALL rpoly(p, degree, zr, zi, fail)
IF (fail) THEN
  WRITE(*, *) ' ** Failure by RPOLY **'
ELSE
  WRITE(*, '(a/ (2g26.18))') ' Real part              Imaginary part',  &
                             (zr(i)%hi, zi(i)%hi, i=1,degree)
END IF

STOP

5000 FORMAT (' EXAMPLE 1. POLYNOMIAL WITH ZEROS 1,2,...,10.')

END PROGRAM test_rpoly
