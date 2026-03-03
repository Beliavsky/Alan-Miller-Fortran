MODULE common_refnum
IMPLICIT NONE

! COMMON /refnum/ nprob,nfev,njev
INTEGER, SAVE :: nprob, nfev, njev

END MODULE common_refnum



PROGRAM test_lmder
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 23:25:02
 
!  **********

!  THIS PROGRAM TESTS CODES FOR THE LEAST-SQUARES SOLUTION OF
!  M NONLINEAR EQUATIONS IN N VARIABLES. IT CONSISTS OF A DRIVER
!  AND AN INTERFACE SUBROUTINE FCN. THE DRIVER READS IN DATA,
!  CALLS THE NONLINEAR LEAST-SQUARES SOLVER, AND FINALLY PRINTS
!  OUT INFORMATION ON THE PERFORMANCE OF THE SOLVER. THIS IS
!  ONLY A SAMPLE DRIVER, MANY OTHER DRIVERS ARE POSSIBLE. THE
!  INTERFACE SUBROUTINE FCN IS NECESSARY TO TAKE INTO ACCOUNT THE
!  FORMS OF CALLING SEQUENCES USED BY THE FUNCTION AND JACOBIAN
!  SUBROUTINES IN THE VARIOUS NONLINEAR LEAST-SQUARES SOLVERS.

!  SUBPROGRAMS CALLED

!    USER-SUPPLIED ...... FCN

!    MINPACK-SUPPLIED ... DPMPAR,ENORM,INITPT,LMDER1,SSQFCN

!    FORTRAN-SUPPLIED ... SQRT

!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!  **********

USE Levenberg_Marquardt
USE common_refnum
IMPLICIT NONE
INTEGER   :: i, ic, info, k, m, n, ntries
INTEGER   :: iwa(40), ma(60), na(60), nf(60), nj(60), np(60), nx(60)
REAL (dp) :: factor, fnorm1, fnorm2, tol
REAL (dp) :: fjac(65,40), fnm(60), fvec(65), x(40)

!     LOGICAL INPUT UNIT IS ASSUMED TO BE NUMBER 5.
!     LOGICAL OUTPUT UNIT IS ASSUMED TO BE NUMBER 6.

INTEGER, PARAMETER :: nread = 5, nwrite = 6

! EXTERNAL fcn

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, fjac, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    REAL (dp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn

  SUBROUTINE ssqfcn (m, n, x, fvec, nprob)
    USE Levenberg_Marquardt
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: m, n
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: fvec(:)
    INTEGER, INTENT(IN)     :: nprob
  END SUBROUTINE ssqfcn

  SUBROUTINE ssqjac (m, n, x, fjac, nprob)
    USE Levenberg_Marquardt
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: m, n, nprob
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: fjac(:,:)
  END SUBROUTINE ssqjac
END INTERFACE

REAL (dp), PARAMETER :: one = 1.0_dp, ten = 10.0_dp

tol = SQRT( EPSILON(one) )
ic = 0
10 READ (nread,50) nprob, n, m, ntries
IF (nprob <= 0) GO TO 30
factor = one
DO  k=1,ntries
  ic = ic + 1
  CALL initpt (n, x, nprob, factor)
  CALL ssqfcn (m, n, x, fvec, nprob)
  fnorm1 = enorm(m, fvec)
  WRITE (nwrite,60) nprob, n, m
  nfev = 0
  njev = 0
  CALL lmder1 (fcn, m, n, x, fvec, fjac, tol, info, iwa)
  CALL ssqfcn (m, n, x, fvec, nprob)
  fnorm2 = enorm(m, fvec)
  np(ic) = nprob
  na(ic) = n
  ma(ic) = m
  nf(ic) = nfev
  nj(ic) = njev
  nx(ic) = info
  fnm(ic) = fnorm2
  WRITE (nwrite,70) fnorm1, fnorm2, nfev, njev, info, x(1:n)
  factor = ten*factor
END DO
GO TO 10

30 WRITE (nwrite,80) ic
WRITE (nwrite,90)
DO  i=1,ic
  WRITE (nwrite,100) np(i), na(i), ma(i), nf(i), nj(i), nx(i), fnm(i)
END DO
STOP

50 FORMAT (4I5)
60 FORMAT (////'      PROBLEM', i5, '      DIMENSIONS', 2I5//)
70 FORMAT ('      INITIAL L2 NORM OF THE RESIDUALS', g15.7//  &
           '      FINAL L2 NORM OF THE RESIDUALS  ', g15.7//  &
           '      NUMBER OF FUNCTION EVALUATIONS  ', i10//  &
           '      NUMBER OF JACOBIAN EVALUATIONS  ', i10//  &
           '      EXIT PARAMETER', t39, i10//  &
           '      FINAL APPROXIMATE SOLUTION'// (t6, 5g15.7))
80 FORMAT (' SUMMARY OF ', i3, ' CALLS TO LMDER1'/)
90 FORMAT (' NPROB   N    M   NFEV  NJEV  INFO  FINAL L2 NORM'/)
100 FORMAT (3I5, 3I6, ' ', g15.7)


CONTAINS


SUBROUTINE initpt (n, x, nprob, factor)
INTEGER, INTENT(IN)     :: n, nprob
REAL (dp), INTENT(IN)   :: factor
REAL (dp), INTENT(OUT)  :: x(:)
!     **********

!     SUBROUTINE INITPT

!     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR THE
!     FUNCTIONS DEFINED BY SUBROUTINE SSQFCN. THE SUBROUTINE RETURNS
!     IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
!     THE 11TH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
!     THIS CASE, IF FACTOR IS NOT UNITY, THEN THE SUBROUTINE RETURNS
!     THE VECTOR  X(J) = FACTOR, J=1,...,N.

!     THE SUBROUTINE STATEMENT IS

!       SUBROUTINE INITPT(N,X,NPROB,FACTOR)

!     WHERE

!       N IS A POSITIVE INTEGER INPUT VARIABLE.

!       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE STANDARD
!         STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY FACTOR.

!       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
!         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.

!       FACTOR IS AN INPUT VARIABLE WHICH SPECIFIES THE MULTIPLE OF
!         THE STANDARD STARTING POINT.  IF FACTOR IS UNITY, NO
!         MULTIPLICATION IS PERFORMED.

!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!     **********
INTEGER   :: j
REAL (dp) :: h
REAL (dp), PARAMETER :: zero = 0.0_dp, half = 0.5_dp, one = 1.0_dp,  &
                        two = 2.0_dp, three = 3.0_dp, five = 5.0_dp, &
                        seven = 7.0_dp, ten = 10.0_dp, twenty = 20.0_dp, &
                        twntf = 25.0_dp
REAL (dp), PARAMETER :: c1 = 1.2_dp, c2 = 0.25_dp, c3 = 0.39_dp,    &
                        c4 = 0.415_dp, c5 = 0.02_dp, c6 = 4000._dp, &
                        c7 = 250._dp, c8 = 0.3_dp, c9 = 0.4_dp,     &
                        c10 = 1.5_dp, c11 = 0.01_dp, c12 = 1.3_dp,  &
                        c13 = 0.65_dp, c14 = 0.7_dp, c15 = 0.6_dp,  &
                        c16 = 4.5_dp, c17 = 5.5_dp

!     SELECTION OF INITIAL POINT.

SELECT CASE ( nprob )
  CASE (    1:3)   !     LINEAR FUNCTION - FULL RANK OR RANK 1.

    x(1:n) = one

  CASE (    4)     !     ROSENBROCK FUNCTION.

    x(1) = -c1
    x(2) = one

  CASE (    5)     !     HELICAL VALLEY FUNCTION.

    x(1) = -one
    x(2) = zero
    x(3) = zero

  CASE (    6)     !     POWELL SINGULAR FUNCTION.

    x(1) = three
    x(2) = -one
    x(3) = zero
    x(4) = one

  CASE (    7)     !     FREUDENSTEIN AND ROTH FUNCTION.

x(1) = half
x(2) = -two

  CASE (    8)     !     BARD FUNCTION.

x(1) = one
x(2) = one
x(3) = one

  CASE (    9)     !     KOWALIK AND OSBORNE FUNCTION.

x(1) = c2
x(2) = c3
x(3) = c4
x(4) = c3

  CASE (   10)     !     MEYER FUNCTION.

    x(1) = c5
    x(2) = c6
    x(3) = c7

  CASE (   11)     !     WATSON FUNCTION.

    x(1:n) = zero

  CASE (   12)     !     BOX 3-DIMENSIONAL FUNCTION.

    x(1) = zero
    x(2) = ten
    x(3) = twenty

  CASE (   13)     !     JENNRICH AND SAMPSON FUNCTION.

    x(1) = c8
    x(2) = c9

  CASE (   14)     !     BROWN AND DENNIS FUNCTION.

    x(1) = twntf
    x(2) = five
    x(3) = -five
    x(4) = -one

  CASE (   15)     !     CHEBYQUAD FUNCTION.

    h = one / DBLE(n+1)
    DO  j=1,n
      x(j) = DBLE(j)*h
    END DO

  CASE (   16)     !     BROWN ALMOST-LINEAR FUNCTION.

    x(1:n) = half

  CASE (   17)     !     OSBORNE 1 FUNCTION.

    x(1) = half
    x(2) = c10
    x(3) = -one
    x(4) = c11
    x(5) = c5

  CASE (   18)     !     OSBORNE 2 FUNCTION.

    x(1) = c12
    x(2) = c13
    x(3) = c13
    x(4) = c14
    x(5) = c15
    x(6) = three
    x(7) = five
    x(8) = seven
    x(9) = two
    x(10) = c16
    x(11) = c17

END SELECT

!     COMPUTE MULTIPLE OF INITIAL POINT.

IF (factor == one) GO TO 250
IF (nprob == 11) GO TO 230
x(1:n) = factor*x(1:n)
GO TO 250

230 x(1:n) = factor

250 RETURN

!     LAST CARD OF SUBROUTINE INITPT.

END SUBROUTINE initpt

!     LAST CARD OF DRIVER.

END PROGRAM test_lmder



SUBROUTINE fcn (m, n, x, fvec, fjac, iflag)
USE common_refnum
IMPLICIT NONE
INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)        :: m, n
REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN OUT)  :: fvec(:)
REAL (dp), INTENT(OUT)     :: fjac(:,:)
INTEGER, INTENT(IN OUT)    :: iflag

INTERFACE
  SUBROUTINE ssqfcn (m, n, x, fvec, nprob)
    USE Levenberg_Marquardt
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: m, n
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: fvec(:)
    INTEGER, INTENT(IN)     :: nprob
  END SUBROUTINE ssqfcn

  SUBROUTINE ssqjac (m, n, x, fjac, nprob)
    USE Levenberg_Marquardt
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: m, n, nprob
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: fjac(:,:)
  END SUBROUTINE ssqjac
END INTERFACE

!  **********

!  THE CALLING SEQUENCE OF FCN SHOULD BE IDENTICAL TO THE
!  CALLING SEQUENCE OF THE FUNCTION SUBROUTINE IN THE NONLINEAR
!  LEAST-SQUARES SOLVER. FCN SHOULD ONLY CALL THE TESTING
!  FUNCTION AND JACOBIAN SUBROUTINES SSQFCN AND SSQJAC WITH
!  THE APPROPRIATE VALUE OF PROBLEM NUMBER (NPROB).

!  SUBPROGRAMS CALLED

!    MINPACK-SUPPLIED ... SSQFCN,SSQJAC

!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!  **********

IF (iflag == 1) CALL ssqfcn (m, n, x, fvec, nprob)
IF (iflag == 2) CALL ssqjac (m, n, x, fjac, nprob)
IF (iflag == 1) nfev = nfev + 1
IF (iflag == 2) njev = njev + 1
RETURN

!     LAST CARD OF INTERFACE SUBROUTINE FCN.

END SUBROUTINE fcn



SUBROUTINE ssqjac (m, n, x, fjac, nprob)
USE Levenberg_Marquardt
IMPLICIT NONE
INTEGER, INTENT(IN)     :: m, n, nprob
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: fjac(:,:)

!  **********

!  SUBROUTINE SSQJAC

!  THIS SUBROUTINE DEFINES THE JACOBIAN MATRICES OF EIGHTEEN
!  NONLINEAR LEAST SQUARES PROBLEMS. THE PROBLEM DIMENSIONS ARE
!  AS DESCRIBED IN THE PROLOGUE COMMENTS OF SSQFCN.

!  THE SUBROUTINE STATEMENT IS

!    SUBROUTINE SSQJAC(M,N,X,FJAC,LDFJAC,NPROB)

!  WHERE

!    M AND N ARE POSITIVE INTEGER INPUT VARIABLES. N MUST NOT EXCEED M.

!    X IS AN INPUT ARRAY OF LENGTH N.

!    FJAC IS AN M BY N OUTPUT ARRAY WHICH CONTAINS THE JACOBIAN
!      MATRIX OF THE NPROB FUNCTION EVALUATED AT X.

!    LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
!      WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.

!    NPROB IS A POSITIVE INTEGER VARIABLE WHICH DEFINES THE
!      NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.

!  SUBPROGRAMS CALLED

!    FORTRAN-SUPPLIED ... DATAN,DCOS,DEXP,DSIN,SQRT

!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!  **********

INTEGER   :: i, j, k, mm1, nm1
REAL (dp) :: div, dx, prod, s2, temp, ti, tmp1, tmp2, tmp3, tmp4, tpi
REAL (dp), PARAMETER :: v(11) = (/ 4.0D0, 2.0D0, 1.0D0, 5.0D-1, 2.5D-1, &
                                   1.67D-1, 1.25D-1, 1.0D-1, 8.33D-2,  &
                                   7.14D-2, 6.25D-2 /)

REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp,  &
                        three = 3.0_dp, four = 4.0_dp, five = 5.0_dp, &
                        eight = 8.0_dp, ten = 10.0_dp, c14 = 14.0_dp, &
                        c20 = 20.0_dp, c29 = 29.0_dp, c45 = 45.0_dp,  &
                        c100 = 100.0_dp

!     JACOBIAN ROUTINE SELECTOR.

SELECT CASE ( nprob )
  CASE (    1)     !     LINEAR FUNCTION - FULL RANK.

    temp = two / DBLE(m)
    fjac(1:m,1:n) = -temp
    DO  j=1,n
      fjac(j,j) = fjac(j,j) + one
    END DO

  CASE (    2)     !     LINEAR FUNCTION - RANK 1.

    DO  j=1,n
      DO  i=1,m
        fjac(i,j) = i * j
      END DO
    END DO

  CASE (    3)     !     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.

    fjac(1:m,1:n) = zero
    nm1 = n-1
    mm1 = m-1
    DO  j=2,nm1
      DO  i=2,mm1
        fjac(i,j) = (i-1) * j
      END DO
    END DO

  CASE (    4)     !     ROSENBROCK FUNCTION.

    fjac(1,1) = -c20*x(1)
    fjac(1,2) = ten
    fjac(2,1) = -one
    fjac(2,2) = zero

  CASE (    5)     !     HELICAL VALLEY FUNCTION.

    tpi = eight*ATAN(one)
    temp = x(1)**2 + x(2)**2
    tmp1 = tpi*temp
    tmp2 = SQRT(temp)
    fjac(1,1) = c100*x(2)/tmp1
    fjac(1,2) = -c100*x(1)/tmp1
    fjac(1,3) = ten
    fjac(2,1) = ten*x(1)/tmp2
    fjac(2,2) = ten*x(2)/tmp2
    fjac(2,3) = zero
    fjac(3,1) = zero
    fjac(3,2) = zero
    fjac(3,3) = one

  CASE (    6)     !     POWELL SINGULAR FUNCTION.

    fjac(1:4,1:4) = zero
    fjac(1,1) = one
    fjac(1,2) = ten
    fjac(2,3) = SQRT(five)
    fjac(2,4) = -fjac(2,3)
    fjac(3,2) = two*(x(2) - two*x(3))
    fjac(3,3) = -two*fjac(3,2)
    fjac(4,1) = two*SQRT(ten)*(x(1) - x(4))
    fjac(4,4) = -fjac(4,1)

  CASE (    7)     !     FREUDENSTEIN AND ROTH FUNCTION.

    fjac(1,1) = one
    fjac(1,2) = x(2)*(ten-three*x(2)) - two
    fjac(2,1) = one
    fjac(2,2) = x(2)*(two+three*x(2)) - c14

  CASE (    8)     !     BARD FUNCTION.

DO  i=1,15
  tmp1 = i
  tmp2 = 16 - i
  tmp3 = tmp1
  IF (i > 8) tmp3 = tmp2
  tmp4 = (x(2)*tmp2 + x(3)*tmp3)**2
  fjac(i,1) = -one
  fjac(i,2) = tmp1*tmp2 / tmp4
  fjac(i,3) = tmp1*tmp3 / tmp4
END DO

  CASE (    9)     !     KOWALIK AND OSBORNE FUNCTION.

    DO  i=1,11
      tmp1 = v(i)*(v(i) + x(2))
      tmp2 = v(i)*(v(i) + x(3)) + x(4)
      fjac(i,1) = -tmp1 / tmp2
      fjac(i,2) = -v(i)*x(1) / tmp2
      fjac(i,3) = fjac(i,1)*fjac(i,2)
      fjac(i,4) = fjac(i,3) / v(i)
    END DO

  CASE (   10)     !     MEYER FUNCTION.

DO  i=1,16
  temp = five*i + c45 + x(3)
  tmp1 = x(2) / temp
  tmp2 = EXP(tmp1)
  fjac(i,1) = tmp2
  fjac(i,2) = x(1)*tmp2 / temp
  fjac(i,3) = -tmp1*fjac(i,2)
END DO

  CASE (   11)     !     WATSON FUNCTION.

    DO  i=1,29
      div = i / c29
      s2 = zero
      dx = one
      DO  j=1,n
        s2 = s2 + dx*x(j)
        dx = div*dx
      END DO
      temp = two*div*s2
      dx = one / div
      DO  j=1,n
        fjac(i,j) = dx*(DBLE(j-1) - temp)
        dx = div*dx
      END DO
    END DO
    fjac(30:31,1:n) = zero
    fjac(30,1) = one
    fjac(31,1) = -two*x(1)
    fjac(31,2) = one

  CASE (   12)     !     BOX 3-DIMENSIONAL FUNCTION.

DO  i=1,m
  temp = DBLE(i)
  tmp1 = temp/ten
  fjac(i,1) = -tmp1*EXP(-tmp1*x(1))
  fjac(i,2) = tmp1*EXP(-tmp1*x(2))
  fjac(i,3) = EXP(-temp)-EXP(-tmp1)
END DO

  CASE (   13)     !     JENNRICH AND SAMPSON FUNCTION.

    DO  i=1,m
      temp = DBLE(i)
      fjac(i,1) = -temp*EXP(temp*x(1))
      fjac(i,2) = -temp*EXP(temp*x(2))
    END DO

  CASE (   14)     !     BROWN AND DENNIS FUNCTION.

DO  i=1,m
  temp = DBLE(i) / five
  ti = SIN(temp)
  tmp1 = x(1) + temp*x(2) - EXP(temp)
  tmp2 = x(3) + ti*x(4) - COS(temp)
  fjac(i,1) = two*tmp1
  fjac(i,2) = temp*fjac(i,1)
  fjac(i,3) = two*tmp2
  fjac(i,4) = ti*fjac(i,3)
END DO

  CASE (   15)     !     CHEBYQUAD FUNCTION.

    dx = one/DBLE(n)
    DO  j=1,n
      tmp1 = one
      tmp2 = two*x(j) - one
      temp = two*tmp2
      tmp3 = zero
      tmp4 = two
      DO  i=1,m
        fjac(i,j) = dx*tmp4
        ti = four*tmp2 + temp*tmp4 - tmp3
        tmp3 = tmp4
        tmp4 = ti
        ti = temp*tmp2 - tmp1
        tmp1 = tmp2
        tmp2 = ti
      END DO
    END DO

  CASE (   16)     !     BROWN ALMOST-LINEAR FUNCTION.

    prod = one
    DO  j=1,n
      prod = x(j)*prod
      fjac(1:n,j) = one
      fjac(j,j) = two
    END DO
    DO  j=1,n
      temp = x(j)
      IF (temp /= zero) GO TO 440
      temp = one
      prod = one
      DO  k=1,n
        IF (k /= j) prod = x(k)*prod
      END DO
      440 fjac(n,j) = prod / temp
    END DO

  CASE (   17)     !     OSBORNE 1 FUNCTION.

    DO  i=1,33
      temp = ten*DBLE(i-1)
      tmp1 = EXP(-x(4)*temp)
      tmp2 = EXP(-x(5)*temp)
      fjac(i,1) = -one
      fjac(i,2) = -tmp1
      fjac(i,3) = -tmp2
      fjac(i,4) = temp*x(2)*tmp1
      fjac(i,5) = temp*x(3)*tmp2
    END DO

  CASE (   18)     !     OSBORNE 2 FUNCTION.

    DO  i=1,65
      temp = DBLE(i-1) / ten
      tmp1 = EXP(-x(5)*temp)
      tmp2 = EXP(-x(6)*(temp - x(9))**2)
      tmp3 = EXP(-x(7)*(temp - x(10))**2)
      tmp4 = EXP(-x(8)*(temp - x(11))**2)
      fjac(i,1) = -tmp1
      fjac(i,2) = -tmp2
      fjac(i,3) = -tmp3
      fjac(i,4) = -tmp4
      fjac(i,5) = temp*x(1)*tmp1
      fjac(i,6) = x(2)*(temp - x(9))**2*tmp2
      fjac(i,7) = x(3)*(temp - x(10))**2*tmp3
      fjac(i,8) = x(4)*(temp - x(11))**2*tmp4
      fjac(i,9) = -two*x(2)*x(6)*(temp - x(9))*tmp2
      fjac(i,10) = -two*x(3)*x(7)*(temp - x(10))*tmp3
      fjac(i,11) = -two*x(4)*x(8)*(temp - x(11))*tmp4
    END DO

END SELECT

RETURN

!     LAST CARD OF SUBROUTINE SSQJAC.

END SUBROUTINE ssqjac



SUBROUTINE ssqfcn (m, n, x, fvec, nprob)
USE Levenberg_Marquardt
IMPLICIT NONE
INTEGER, INTENT(IN)     :: m, n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: fvec(:)
INTEGER, INTENT(IN)     :: nprob

!  **********

!  SUBROUTINE SSQFCN

!  THIS SUBROUTINE DEFINES THE FUNCTIONS OF EIGHTEEN NONLINEAR
!  LEAST SQUARES PROBLEMS. THE ALLOWABLE VALUES OF (M,N) FOR
!  FUNCTIONS 1,2 AND 3 ARE VARIABLE BUT WITH M >= N.
!  FOR FUNCTIONS 4,5,6,7,8,9 AND 10 THE VALUES OF (M,N) ARE
!  (2,2),(3,3),(4,4),(2,2),(15,3),(11,4) AND (16,3), RESPECTIVELY.
!  FUNCTION 11 (WATSON) HAS M = 31 WITH N USUALLY 6 OR 9.
!  HOWEVER, ANY N, N = 2,...,31, IS PERMITTED.
!  FUNCTIONS 12,13 AND 14 HAVE N = 3,2 AND 4, RESPECTIVELY, BUT
!  ALLOW ANY M >= N, WITH THE USUAL CHOICES BEING 10,10 AND 20.
!  FUNCTION 15 (CHEBYQUAD) ALLOWS M AND N VARIABLE WITH M >= N.
!  FUNCTION 16 (BROWN) ALLOWS N VARIABLE WITH M = N.
!  FOR FUNCTIONS 17 AND 18, THE VALUES OF (M,N) ARE
!  (33,5) AND (65,11), RESPECTIVELY.

!  THE SUBROUTINE STATEMENT IS

!    SUBROUTINE SSQFCN(M,N,X,FVEC,NPROB)

!  WHERE

!    M AND N ARE POSITIVE INTEGER INPUT VARIABLES. N MUST NOT
!      EXCEED M.

!    X IS AN INPUT ARRAY OF LENGTH N.

!    FVEC IS AN OUTPUT ARRAY OF LENGTH M WHICH CONTAINS THE NPROB
!      FUNCTION EVALUATED AT X.

!    NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
!      NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.

!  SUBPROGRAMS CALLED

!    FORTRAN-SUPPLIED ... DATAN,DCOS,EXP,DSIN,SQRT,DSIGN

!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
!  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE

!  **********

INTEGER   :: i, iev, j, nm1
REAL (dp) :: div, dx, prod, total, s1, s2, temp, ti, tmp1, tmp2, tmp3, tmp4, tpi

REAL (dp), PARAMETER :: zero = 0.0_dp, zp25 = 0.25_dp, zp5 = 0.5_dp,  &
    one = 1.0_dp, two = 2.0_dp, five = 5.0_dp, eight = 8.0_dp, ten = 10._dp, &
    c13 = 13._dp, c14 = 14._dp, c29 = 29._dp, c45 = 45._dp
REAL (dp), PARAMETER :: v(11) = (/  &
    4.0D0, 2.0D0, 1.0D0, 5.0E-1_dp, 2.5E-1_dp, 1.67E-1_dp, 1.25E-1_dp,  &
    1.0E-1_dp, 8.33E-2_dp, 7.14E-2_dp, 6.25E-2_dp /)
REAL (dp), PARAMETER :: y1(15) = (/  &
    1.4E-1_dp, 1.8E-1_dp, 2.2E-1_dp, 2.5E-1_dp, 2.9E-1_dp, 3.2E-1_dp,  &
    3.5E-1_dp, 3.9E-1_dp, 3.7E-1_dp, 5.8E-1_dp, 7.3E-1_dp, 9.6E-1_dp,  &
    1.34_dp, 2.1_dp, 4.39_dp /)
REAL (dp), PARAMETER :: y2(11) = (/   &
    1.957E-1_dp, 1.947E-1_dp, 1.735E-1_dp, 1.6E-1_dp, 8.44E-2_dp, 6.27E-2_dp, &
    4.56E-2_dp, 3.42E-2_dp, 3.23E-2_dp, 2.35E-2_dp, 2.46E-2_dp /)
REAL (dp), PARAMETER :: y3(16) = (/   &
    3.478D4, 2.861D4, 2.365D4, 1.963D4, 1.637D4, 1.372D4, 1.154D4, 9.744D3,  &
    8.261D3, 7.03D3, 6.005D3, 5.147D3, 4.427D3, 3.82D3, 3.307D3, 2.872D3 /)
REAL (dp), PARAMETER :: y4(33) = (/   &
    8.44E-1_dp, 9.08E-1_dp, 9.32E-1_dp, 9.36E-1_dp, 9.25E-1_dp, 9.08E-1_dp,  &
    8.81E-1_dp, 8.5E-1_dp, 8.18E-1_dp, 7.84E-1_dp, 7.51E-1_dp, 7.18E-1_dp,   &
    6.85E-1_dp, 6.58E-1_dp, 6.28E-1_dp, 6.03E-1_dp, 5.8E-1_dp, 5.58E-1_dp,   &
    5.38E-1_dp, 5.22E-1_dp, 5.06E-1_dp, 4.9E-1_dp, 4.78E-1_dp, 4.67E-1_dp,   &
    4.57E-1_dp, 4.48E-1_dp, 4.38E-1_dp, 4.31E-1_dp, 4.24E-1_dp, 4.2E-1_dp,   &
    4.14E-1_dp, 4.11E-1_dp, 4.06E-1_dp /)
REAL (dp), PARAMETER :: y5(65) = (/   &
    1.366_dp, 1.191_dp, 1.112_dp, 1.013_dp, 9.91E-1_dp, 8.85E-1_dp, 8.31E-1_dp,   &
    8.47E-1_dp, 7.86E-1_dp, 7.25E-1_dp, 7.46E-1_dp, 6.79E-1_dp, 6.08E-1_dp,  &
    6.55E-1_dp, 6.16E-1_dp, 6.06E-1_dp, 6.02E-1_dp, 6.26E-1_dp, 6.51E-1_dp,  &
    7.24E-1_dp, 6.49E-1_dp, 6.49E-1_dp, 6.94E-1_dp, 6.44E-1_dp, 6.24E-1_dp,  &
    6.61E-1_dp, 6.12E-1_dp, 5.58E-1_dp, 5.33E-1_dp, 4.95E-1_dp, 5.0E-1_dp,   &
    4.23E-1_dp, 3.95E-1_dp, 3.75E-1_dp, 3.72E-1_dp, 3.91E-1_dp, 3.96E-1_dp,  &
    4.05E-1_dp, 4.28E-1_dp, 4.29E-1_dp, 5.23E-1_dp, 5.62E-1_dp, 6.07E-1_dp,  &
    6.53E-1_dp, 6.72E-1_dp, 7.08E-1_dp, 6.33E-1_dp, 6.68E-1_dp, 6.45E-1_dp,  &
    6.32E-1_dp, 5.91E-1_dp, 5.59E-1_dp, 5.97E-1_dp, 6.25E-1_dp, 7.39E-1_dp,  &
    7.1E-1_dp, 7.29E-1_dp, 7.2E-1_dp, 6.36E-1_dp, 5.81E-1_dp, 4.28E-1_dp,    &
    2.92E-1_dp, 1.62E-1_dp,  9.8E-2_dp, 5.4E-2_dp /)

!     FUNCTION ROUTINE SELECTOR.

SELECT CASE ( nprob )
  CASE (    1)     !     LINEAR FUNCTION - FULL RANK.

    total = SUM( x(1:n) )
    temp = two*total/DBLE(m) + one
    DO  i=1,m
      fvec(i) = -temp
      IF (i <= n) fvec(i) = fvec(i) + x(i)
    END DO

  CASE (    2)     !     LINEAR FUNCTION - RANK 1.

    total = zero
    DO  j=1,n
      total = total + DBLE(j)*x(j)
    END DO
    DO  i=1,m
      fvec(i) = DBLE(i)*total - one
    END DO

  CASE (    3)     !     LINEAR FUNCTION - RANK 1 WITH ZERO COLUMNS AND ROWS.

    total = zero
    nm1 = n-1
    DO  j=2,nm1
      total = total + DBLE(j)*x(j)
    END DO
    DO  i=1,m
      fvec(i) = DBLE(i-1)*total - one
    END DO
    fvec(m) = -one

  CASE (    4)     !     ROSENBROCK FUNCTION.

    fvec(1) = ten*(x(2) - x(1)**2)
    fvec(2) = one - x(1)

  CASE (    5)     !     HELICAL VALLEY FUNCTION.

    tpi = eight*ATAN(one)
    tmp1 = SIGN(zp25,x(2))
    IF (x(1) > zero) tmp1 = ATAN(x(2)/x(1))/tpi
    IF (x(1) < zero) tmp1 = ATAN(x(2)/x(1))/tpi + zp5
    tmp2 = SQRT(x(1)**2 + x(2)**2)
    fvec(1) = ten*(x(3) - ten*tmp1)
    fvec(2) = ten*(tmp2 - one)
    fvec(3) = x(3)

  CASE (    6)     !     POWELL SINGULAR FUNCTION.

    fvec(1) = x(1) + ten*x(2)
    fvec(2) = SQRT(five)*(x(3) - x(4))
    fvec(3) = (x(2) - two*x(3))**2
    fvec(4) = SQRT(ten)*(x(1) - x(4))**2

  CASE (    7)     !     FREUDENSTEIN AND ROTH FUNCTION.

    fvec(1) = -c13 + x(1) + ((five - x(2))*x(2) - two)*x(2)
    fvec(2) = -c29 + x(1) + ((one + x(2))*x(2) - c14)*x(2)

  CASE (    8)     !     BARD FUNCTION.

    DO  i=1,15
      tmp1 = DBLE(i)
      tmp2 = DBLE(16-i)
      tmp3 = tmp1
      IF (i > 8) tmp3 = tmp2
      fvec(i) = y1(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
    END DO

  CASE (    9)     !     KOWALIK AND OSBORNE FUNCTION.

    DO  i=1,11
      tmp1 = v(i)*(v(i) + x(2))
      tmp2 = v(i)*(v(i) + x(3)) + x(4)
      fvec(i) = y2(i) - x(1)*tmp1/tmp2
    END DO

  CASE (   10)     !     MEYER FUNCTION.

    DO  i=1,16
      temp = five*DBLE(i) + c45 + x(3)
      tmp1 = x(2)/temp
      tmp2 = EXP(tmp1)
      fvec(i) = x(1)*tmp2 - y3(i)
    END DO

  CASE (   11)     !     WATSON FUNCTION.

    DO  i=1,29
      div = DBLE(i)/c29
      s1 = zero
      dx = one
      DO  j=2,n
        s1 = s1 + DBLE(j-1)*dx*x(j)
        dx = div*dx
      END DO
      s2 = zero
      dx = one
      DO  j=1,n
        s2 = s2 + dx*x(j)
        dx = div*dx
      END DO
      fvec(i) = s1 - s2**2 - one
    END DO
    fvec(30) = x(1)
    fvec(31) = x(2) - x(1)**2 - one

  CASE (   12)     !     BOX 3-DIMENSIONAL FUNCTION.

    DO  i=1,m
      temp = DBLE(i)
      tmp1 = temp/ten
      fvec(i) = EXP(-tmp1*x(1)) - EXP(-tmp1*x(2)) + (EXP(-temp) -  &
                EXP(- tmp1))*x(3)
    END DO

  CASE (   13)     !     JENNRICH AND SAMPSON FUNCTION.

    DO  i=1,m
      temp = DBLE(i)
      fvec(i) = two + two*temp - EXP(temp*x(1)) - EXP(temp*x(2))
    END DO

  CASE (   14)     !     BROWN AND DENNIS FUNCTION.

    DO  i=1,m
      temp = DBLE(i)/five
      tmp1 = x(1) + temp*x(2) - EXP(temp)
      tmp2 = x(3) + SIN(temp)*x(4) - COS(temp)
      fvec(i) = tmp1**2 + tmp2**2
    END DO

  CASE (   15)     !     CHEBYQUAD FUNCTION.

    fvec(1:m) = zero
    DO  j=1,n
      tmp1 = one
      tmp2 = two*x(j) - one
      temp = two*tmp2
      DO  i=1,m
        fvec(i) = fvec(i) + tmp2
        ti = temp*tmp2 - tmp1
        tmp1 = tmp2
        tmp2 = ti
      END DO
    END DO
    dx = one/DBLE(n)
    iev = -1
    DO  i=1,m
      fvec(i) = dx*fvec(i)
      IF (iev > 0) fvec(i) = fvec(i) + one/(DBLE(i)**2 - one)
      iev = -iev
    END DO

  CASE (   16)     !     BROWN ALMOST-LINEAR FUNCTION.

    total = -DBLE(n+1)
    prod = one
    DO  j=1,n
      total = total + x(j)
      prod = x(j)*prod
    END DO
    DO  i=1,n
      fvec(i) = x(i) + total
    END DO
    fvec(n) = prod - one

  CASE (   17)     !     OSBORNE 1 FUNCTION.

DO  i=1,33
  temp = ten*DBLE(i-1)
  tmp1 = EXP(-x(4)*temp)
  tmp2 = EXP(-x(5)*temp)
  fvec(i) = y4(i) - (x(1) + x(2)*tmp1 + x(3)*tmp2)
END DO

  CASE (   18)     !     OSBORNE 2 FUNCTION.

    DO  i=1,65
      temp = DBLE(i-1)/ten
      tmp1 = EXP(-x(5)*temp)
      tmp2 = EXP(-x(6)*(temp-x(9))**2)
      tmp3 = EXP(-x(7)*(temp-x(10))**2)
      tmp4 = EXP(-x(8)*(temp-x(11))**2)
      fvec(i) = y5(i) - (x(1)*tmp1 + x(2)*tmp2 + x(3)*tmp3 + x(4)*tmp4)
    END DO

END SELECT

RETURN

!     LAST CARD OF SUBROUTINE SSQFCN.

END SUBROUTINE ssqfcn
