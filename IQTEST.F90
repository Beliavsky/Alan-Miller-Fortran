PROGRAM iqtest
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-01-18  Time: 22:19:39
 
!     TEST PROGRAM FOR IQPACK: DOUBLE PRECISION VERSION

USE iqpack
IMPLICIT NONE

REAL (dp) :: a, alpha, b, beta, frq, qfsum, qfsx
INTEGER   :: i, ier, ik, ind, key, kind, mtk, niwf, nt, nwf, nwts
INTEGER   :: mlt(100), ndx(100), iwf(500)
REAL (dp) :: t(100), wts(100), wf(10000)

REAL (dp), PARAMETER  :: four = 4.0_dp
REAL (dp), PARAMETER  :: pi = 3.14159265358979323846264338327950_dp

!     FUNCTIONS AND SUBROUTINES REFERENCED - CEGQF,CEGQFS,CEIQF,
!         CEIQFS,CGQF,CIQF,CIQFS,CLIQF,CLIQFS,EIQFS,F,COS,SIN

INTERFACE
  FUNCTION f(x, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: x
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION f
END INTERFACE

!     LOGICAL UNIT NUMBERS: LO,LU=OUTPUT  LI=INPUT (NOT USED)
INTEGER, PARAMETER  :: lo = 6, lu = 6

!     SIZE OF KNOTS ARRAY AND WORKSPACE IN DIMENSION STATEMENTS
nwf = 10000
niwf = 500

!     GENERATE FEJER TYPE RULES - I.E. GAUSS CHEBYSHEV KNOTS AND
!     LEGENDRE WEIGHT FUNCTION
!     SET WEIGHT FUNCTION TYPE AND PARAMETERS
kind = 1
a = one
b = four
!     ALPHA, BETA NOT USED IN LEGENDRE WEIGHT FUNCTION BUT SET ANYWAY
alpha = half
beta = -half/two

!     SET MULTIPLICITY OF KNOTS
mtk = 2
key = 1

!     NUMBER OF KNOTS.
nt = 5

!     IF NT IS CHANGED IT MUST BE <= NTMAX
!     SET SIZE OF WEIGHTS ARRAY
nwts = mtk*nt
frq = pi/two/nt

!     TEST THE DRIVERS
DO  ind = 1,7
  
!        SET UP THE KNOTS FOR NT-POINT FEJER TYPE RULE
!        FOR IND = 2,4,6 SCALE THE KNOTS TO GIVEN A,B
!        OTHERWISE - DEFAULT A,B - DON'T SCALE
  DO  i = 1,nt
    t(i) = COS((two*i-1)*frq)
    IF(MOD(ind,2) == 0) THEN
      IF(kind == 5 .OR. kind == 6 .OR. kind == 8) THEN
        t(i) = t(i)+a
      ELSE
        t(i) = ((b-a)*t(i) + (a+b))/two
      END IF
    END IF
  END DO
!        SET KNOT MULTIPLICITIES
  mlt(1:nt) = mtk
  WRITE(lo,*)'----------------------------------------'
  SELECT CASE ( ind )
    CASE (    1)
      GO TO 30
    CASE (    2)
      GO TO 40
    CASE (    3)
      GO TO 50
    CASE (    4)
      GO TO 60
    CASE (    5)
      GO TO 70
    CASE (    6)
      GO TO 80
    CASE (    7)
      GO TO 90
  END SELECT
  
  30 WRITE(lo,*)' BEGINNING TEST OF CIQFS'
  CALL ciqfs(nt, t, mlt, nwts, wts, ndx, key, kind, alpha, beta, lu,  &
             nwf, wf, ier)
  WRITE(lo,*)' RETURN FROM CIQFS: ERROR INDICATOR IER= ', ier
  CYCLE
  
!        CIQF CALLS CIQFS, CGQF CALLS CGQFS
!        SO THIS TESTS BOTH
  40 WRITE(lo,*)' BEGINNING TEST OF CIQF,CIQFS,CGQF AND CGQFS'
  
  WRITE(lo,*)' WITH ALL CLASSICAL WEIGHT FUNCTIONS'
  a = -half
  b = two
  DO  ik = 1,8
    kind = ik
    beta = two
    IF(kind == 8) beta = -four*four
    WRITE(lo,*)' KNOTS AND WEIGHTS OF GAUSS QF'
    WRITE(lo,*)' COMPUTED BY CGQF(S)'
    CALL cgqf(nt, t, wts, kind, alpha, beta, a, b, lu, ier)
    WRITE(lo,*)
    WRITE(lo,*)' RETURN FROM CGQF(S): ERROR INDICATOR IER= ', ier

!           NOW COMPUTE THE WEIGHTS FOR THE SAME KNOTS BY CIQF
    WRITE(lo,*)' WEIGHTS OF GAUSS QF COMPUTED FROM THE'
    WRITE(lo,*)' KNOTS BY CIQF(S)'
    CALL ciqf(nt, t, mlt, nwts, wts, ndx, key, kind, alpha, beta, a, b, lu, &
              nwf, wf, ier)
    WRITE(lo,*)
    WRITE(lo,*)' RETURN FROM CIQF(S): ERROR INDICATOR IER= ', ier
  END DO
  CYCLE
  
  50 kind = 1
  WRITE(lo,*)' BEGINNING TEST OF CEIQFS'
  qfsx = COS(-one) - COS(one)

!        QFSX IS EXACT VALUE OF INTEGRAL
  CALL ceiqfs(nt, t, mlt, kind, alpha, beta, f, qfsum, nwf, wf, niwf, iwf, ier)
  WRITE(lo,*)
  WRITE(lo,*)' RETURN FROM CEIQFS: ERROR INDICATOR IER= ', ier
  WRITE(lo,*)' INTEGRAL OF SIN(X) ON Õ-1,1å BY FEJER TYPE RULE'
  WRITE(lo,*)' WITH ', nt, ' POINTS OF MULTIPLICITY ', mtk
  WRITE(lo,*)' QUADRATURE FORMULA:', qfsum
  WRITE(lo,*)' EXACT VALUE       :', qfsx
  CYCLE
  
  60 WRITE(lo,*)' BEGINNING TEST OF CEIQF'
  kind = 1
  qfsx = COS(a) - COS(b)

!        QFSX IS EXACT VALUE OF INTEGRAL
  CALL ceiqf(nt, t, mlt, kind, alpha, beta, a, b, f, qfsum, nwf, wf, niwf, &
             iwf, ier)
  WRITE(lo,*)
  WRITE(lo,*)' RETURN FROM CEIQF: ERROR INDICATOR IER= ', ier
  WRITE(lo,*)' INTEGRAL OF SIN(X) FROM ', a, ' TO ', b
  WRITE(lo,*)' BY FEJER TYPE RULE WITH ', nt, ' POINTS'
  WRITE(lo,*)' OF MULTIPLICITY ', mtk
  WRITE(lo,*)' QUADRATURE FORMULA:', qfsum
  WRITE(lo,*)' EXACT VALUE       :', qfsx
  CYCLE
  
  70 WRITE(lo,*)' BEGINNING TEST OF CLIQFS'
  CALL cliqfs(nt, t, wts, kind, alpha, beta, lu, nwf, wf, niwf, iwf, ier)
  WRITE(lo,*)
  WRITE(lo,*)'RETURN FROM CLIQFS. IER= ', ier
  CYCLE
  
  80 WRITE(lo,*)' BEGINNING TEST OF CLIQF AND EIQFS'
  CALL cliqf(nt, t, wts, kind, alpha, beta, a, b, lu, nwf, wf, niwf, iwf, ier)
  WRITE(lo,*)
  WRITE(lo,*)' RETURN FROM CLIQF: ERROR INDICATOR IER= ', ier
  CALL eiqfs(nt, t, wts, f, qfsum, ier)
  WRITE(lo,*)' RETURN FROM EIQFS: ERROR INDICATOR IER= ', ier
  qfsx = COS(a) - COS(b)
  WRITE(lo,*)' INTEGRAL OF SIN(X) FROM ', a, ' TO ', b
  WRITE(lo,*)' BY FEJER TYPE RULE WITH ', nt, ' POINTS'
  WRITE(lo,*)' OF MULTIPLICITY ONE'
  WRITE(lo,*)' QUADRATURE FORMULA:', qfsum
  WRITE(lo,*)' EXACT VALUE       :', qfsx
  CYCLE

  90 WRITE(lo,*)' BEGINNING TEST OF CEGQF'
  nt = 12
  
!        EXPONENTIAL WEIGHT FUNCTION ON (A,B) CHOSEN
  
  kind = 7
  alpha = one
  beta = zero
  a = one
  b = four
  
!        COMPUTE AND EVALUATE THE GAUSS QF
  CALL cegqf(nt, kind, alpha, beta, a, b, f, qfsum, ier)
  WRITE(lo,*)' RETURN FROM CEGQF: ERROR INDICATOR IER= ', ier
  IF(ier /= 0) CYCLE
  
  qfsx = (b-a)*half*(COS(a)-COS(b)) + SIN(b) + SIN(a) - two*SIN((a+b)/two)
  WRITE(lo,*)' INTEGRAL OF X*SIN(X) FROM ', a, ' TO ', b
  WRITE(lo,*)' BY GAUSS-EXPONENTIAL RULE WITH ', nt, ' POINTS'
  WRITE(lo,*)' QUADRATURE FORMULA:', qfsum
  WRITE(lo,*)' EXACT VALUE       :', qfsx
  
  WRITE(lo,*)' BEGINNING TEST OF CEGQFS'
  nt = 12
  
!        EXPONENTIAL WEIGHT FUNCTION, DEFAULT VALUES OF A,B
  
  kind = 7
  alpha = one
  beta = zero
  
!        COMPUTE AND EVALUATE THE GAUSS QF
  CALL cegqfs(nt, kind, alpha, beta, f, qfsum, ier)
  WRITE(lo,*)' RETURN FROM CEGQFS: ERROR INDICATOR IER= ', ier
  IF(ier /= 0) CYCLE
  
  a = -one
  b = one
  qfsx = (b-a)*half*(COS(a)-COS(b)) + SIN(b) + SIN(a) - two*SIN((a+b)/two)
  WRITE(lo,*)' INTEGRAL OF X*SIN(X) FROM ', a, ' TO ', b
  WRITE(lo,*)' BY GAUSS-EXPONENTIAL RULE WITH ', nt, ' POINTS'
  WRITE(lo,*)' QUADRATURE FORMULA:', qfsum
  WRITE(lo,*)' EXACT VALUE       :', qfsx
  
END DO
STOP
END PROGRAM iqtest



FUNCTION f(x, i) RESULT(fn_val)

!     THIS FUNCTION GENERATES VALUES OF THE INTEGRAND
!     AND POSSIBLY ITS DERIVATIVES. SINCE IT APPEARS AS AN ACTUAL
!     PARAMETER IN THE ROUTINE WHICH USES IT, IT SHOULD BE DECLARED
!     IN AN EXTERNAL STATEMENT IN ANY CALLING PROGRAM

USE iqpack, ONLY: dp
IMPLICIT NONE

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: i
REAL (dp)              :: fn_val

INTEGER :: k, l

!     FUNCTIONS AND SUBROUTINES REFERENCED - COS,SIN,MOD

!     RETURNS THE I-TH DERIVATIVE OF SIN(T) AT T = X
!     (ZERO-TH DERIVATIVE IS THE FUNCTION VALUE)
k = MOD(i,4)
l = MOD(k,2)
IF(l == 0) fn_val = SIN(x)
IF(l == 1) fn_val = COS(x)
IF(k > 1) fn_val = - fn_val

RETURN
END FUNCTION f
