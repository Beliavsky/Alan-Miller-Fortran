MODULE common_tun
IMPLICIT NONE
INTEGER, SAVE  :: nprob
END MODULE common_tun



PROGRAM Test677

! (ALGORITHM SIGMA)
! MAIN PROGRAM (TEST VERSION)
! (CALL  SIGMA  VIA THE DRIVER  SIGMA1 )

USE toms667
USE common_tun
IMPLICIT NONE

! COMMON AREA TO PASS TEST-PROBLEM NUMBER  NPROB
! TO THE FUNCTION  FUNCT  WHICH WILL COMPUTE
! THE FUNCTION VALUES OF TEST-PROBLEM  NPROB
! BY CALLING THE TEST-PROBLEM COLLECTION SUBROUTINE  GLOMTF

! COMMON /tun/ nprob

!  X0     INITIAL POINT
!  XMIN   FINAL ESTIMATE OF GLOBAL MINIMUM
!  XMINGL, XMAXGL  MUST BE DIMENSIONED HERE IN ORDER TO CALL
!         THE PRE-EXISTING SUBROUTINE  GLOMIP.
REAL (dp)  :: x0(100), xmin(100), xmingl(100), xmaxgl(100)

INTERFACE
  SUBROUTINE glomip(nprob, n, x0, xmin, xmax)
    USE toms667, ONLY: dp
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: nprob
    INTEGER, INTENT(OUT)    :: n
    REAL (dp), INTENT(OUT)  :: x0(:)
    REAL (dp), INTENT(OUT)  :: xmin(:)
    REAL (dp), INTENT(OUT)  :: xmax(:)
  END SUBROUTINE glomip
END INTERFACE

INTEGER    :: iout, iprint, n, nfev, nsuc
REAL (dp)  :: fmin

! INPUT PROBLEM NUMBER
10 WRITE (6,5000, ADVANCE='NO')
READ (5,5100) nprob
WRITE (6,5200) nprob

! TERMINATE OR CONTINUE
IF (nprob > 0 .AND. nprob < 38) THEN

! CALL  GLOMIP  TO GET PROBLEM DIMENSION  N  AND INITIAL POINT  X0
! NOTE THAT  GLOMIP  RETURNS ALSO THE BOUNDARIES  XMINGL ,  XMAXGL
! OF THE OBSERVATION REGION (NOT NEEDED HERE)
  CALL glomip(nprob,n,x0,xmingl,xmaxgl)

! SET  NSUC  SO AS TO HAVE GOOD CHANCES, WITHOUT PROHIBITIVE
! COMPUTATIONAL EFFORT
  nsuc = 5

! SET  IPRINT  SO AS TO HAVE A MODERATE OUTPUT
  iprint = 0

! CALL DRIVER SUBROUTINE  SIGMA1
  CALL sigma1(n,x0,nsuc,iprint,xmin,fmin,nfev,iout)

! GO TO THE NEXT PROBLEM
  GO TO 10
END IF

! END OF TEST PROBLEMS
WRITE (6,5300)

STOP
5000 FORMAT (///' INPUT PROBLEM NUMBER (1 TO 37, 0 = STOP): ')
5100 FORMAT (i2)
5200 FORMAT (//' PROBLEM NUMBER = ',i2//)
5300 FORMAT (/' END OF TEST PROBLEMS '/)
END PROGRAM Test677



FUNCTION funct(n,x) RESULT(fn_val)

! COMPUTES THE FUNCTION VALUES OF TEST-PROBLEM  NPROB
! BY CALLING THE SUBROUTINE  GLOMTF .

USE common_tun
USE toms667, ONLY: dp
IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(:)
REAL (dp)              :: fn_val

INTERFACE
  SUBROUTINE glomtf(nprob, n, x, funz)
    USE toms667, ONLY: dp
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: nprob
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: funz
  END SUBROUTINE glomtf
END INTERFACE

REAL (dp)  :: f

! COMMON /tun/ nprob

CALL glomtf(nprob,n,x,f)
fn_val = f

RETURN
END FUNCTION funct




SUBROUTINE glomtf(nprob, n, x, funz)

!  THE SUBROUTINE  GLOMTF  PROVIDES THE CODING OF 37 REAL-VALUED
!  FUNCTIONS OF  N  REAL VARIABLES, TO BE USED, TOGETHER WITH THE
!  SUBROUTINE  GLOMIP , TO DEFINE 37 TEST PROBLEMS FOR GLOBAL
!  MINIMIZATION SOFTWARE.

!  THE SUBROUTINE  GLOMTF  RETURNS IN  FUNZ  THE FUNCTION VALUE
!  AT THE POINT  X = (X1,X2,...,XN)  FOR THE FUNCTION DEFINED BY
!  PROBLEM NUMBER  NPROB .

!  CALLING STATEMENT

!     CALL GLOMTF (NPROB, N, X, FUNZ)

!  DESCRIPTION OF THE CALL PARAMETERS
!     (THE FORTRAN IMPLICIT TYPE DEFINITION FOR INTEGERS IS USED.
!     ALL NON-INTEGER PARAMETERS ARE DOUBLE-PRECISION).

!  NPROB   IS THE (INPUT) TEST-PROBLEM NUMBER

!  N   IS THE (INPUT) DIMENSION OF THE PROBLEM

!  X   IS AN (INPUT) N-VECTOR CONTAINING THE INDEPENDENT VARIABLES

!  FUNZ   IS THE (OUTPUT) VALUE AT  X  OF THE FUNCTION DEFINED BY
!     PROBLEM NUMBER  NPROB .

USE toms667, ONLY: dp
IMPLICIT NONE

INTEGER, INTENT(IN)     :: nprob
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: funz

REAL (dp)  :: xx, di, yy, beta, s1, s2, a, b, r2, r4, r8, u, v, uu, vv, s
REAL (dp)  :: rg, rp, h, pund, range

REAL (dp)  :: y(10)

REAL (dp)  :: pi = 3.14159265358979323846D0, p1a = 0.1D0, p4a = 0.1D0,  &
              p5a = 0.1D0, p6a = 0.1D0
REAL (dp)  :: p9a = -1.4251284283197609708D0, p9b = -0.80032110047197312466D0
REAL (dp)  :: p17a = 1.275D0

REAL (dp)  :: p20a(10,4) = RESHAPE( (/  &
    4.d0, 1.d0, 8.d0, 6.d0, 3.d0, 2.d0, 5.d0, 8.d0, 6.d0, 7.d0,  &
    4.d0, 1.d0, 8.d0, 6.d0, 7.d0, 9.d0, 5.d0, 1.d0, 2.d0, 3.6D0, &
    4.d0, 1.d0, 8.d0, 6.d0, 3.d0, 2.d0, 3.d0, 8.d0, 6.d0, 7.d0,  &
    4.d0, 1.d0, 8.d0, 6.d0, 7.d0, 9.d0, 3.d0, 1.d0, 2.d0, 3.6D0 /), (/ 10, 4 /) )
REAL (dp)  :: p20b(10) = (/  &
    0.1D0, 0.2D0, 0.2D0, 0.4D0, 0.4D0, 0.6D0, 0.3D0, 0.7D0, 0.5D0, 0.5D0 /)
REAL (dp)  :: p21a(4,3) = RESHAPE( (/  3.d0,  0.1D0, 3.d0, 0.1D0,  &
                                      10.d0, 10.d0, 10.d0, 10.d0,  &
                                      30.d0, 35.d0, 30.d0, 35.d0 /), (/ 4, 3 /) )
REAL (dp)  :: p21b(4,3) = RESHAPE( (/ 0.3689D0, 0.4699D0, 0.1091D0, 0.03815D0, &
                                      0.1170D0, 0.4387D0, 0.8732D0, 0.5743D0, &
                                      0.2673D0, 0.7470D0, 0.5547D0, 0.8828D0 /), &
                                      (/ 4, 3 /) )
REAL (dp)  :: p21c(4) = (/ 1.0d0, 1.2D0, 3.d0, 3.2D0 /)
REAL (dp)  :: p21d =  -69.d0
REAL (dp)  :: p22a(4,6) = RESHAPE( (/ 10.d0, 0.05D0, 3.d0, 17.d0,  &
                                       3.d0, 10.d0, 3.5D0, 8.d0,   &
                                      17.d0, 17.d0, 1.7D0, 0.05D0, &
                                      3.5D0, 0.1D0, 10.d0, 10.d0,  &
                                      1.7D0, 8.d0,  17.d0, 0.1D0,  &
                                      8.d0, 14.d0,  8.d0, 14.d0 /), (/ 4, 6 /) )
REAL (dp)  :: p22b(4,6) = RESHAPE( (/ 0.1312D0, 0.2329D0, 0.2348D0, 0.4047D0, &
                                      0.1696D0, 0.4135D0, 0.1451D0, 0.8828D0, &
                                      0.5569D0, 0.8307D0, 0.3522D0, 0.8732D0, &
                                      0.0124D0, 0.3736D0, 0.2883D0, 0.5743D0, &
                                      0.8283D0, 0.1004D0, 0.3047D0, 0.1091D0, &
                                      0.5886D0, 0.9991D0, 0.6650D0, 0.0381D0 /), &
                                      (/ 4, 6 /) )
REAL (dp)  :: p22c(4) = (/ 1.d0, 1.2D0, 3.d0, 3.2D0 /)
REAL (dp)  :: p22d = -69.d0
REAL (dp)  :: p36a = 100.0_dp, p36b = 1.0_dp, p36c = 10.0_dp, p36d = 0.98_dp
REAL (dp)  :: p37a = 10.0_dp, p37b = 1.0_dp, p37c = 10.0_dp, p37d = 0.98_dp

INTEGER    :: i, j, m

SELECT CASE (nprob)
  CASE (1)
    !    1  A FOURTH ORDER POLYNOMIAL  (N = 1)

    funz = ((0.25D0*x(1)*x(1) - 0.5D0)*x(1) + p1a) * x(1)

  CASE (2)
    !    2  GOLDSTEIN SIXTH ORDER POLYNOMIAL  (N = 1)

    xx = x(1) * x(1)
    funz = ((xx-15.d0)*xx+27.d0) * xx + 250.d0

  CASE (3)
    !    3  ONE-DIMENSIONAL PENALIZED SHUBERT FUNCTION  (N = 1)

    funz = 0.d0
    DO  i = 1, 5
      di = i
      funz = funz + di * COS((di+1.d0)*x(1)+di)
    END DO
    IF (ABS(x(1)) > 10.d0) funz = funz + penfun(x(1),10.d0,100.d0,2)

  CASE (4)
    !    4  A FOURTH ORDER POLYNOMIAL IN TWO VARIABLES  (N = 2)

    xx = x(1) * x(1)
    yy = x(2) * x(2)
    funz = 0.25D0 * xx * xx - 0.5D0 * xx + p4a * x(1) + 0.5D0 * yy

  CASE (5)
    !    5  A FUNCTION WITH A SINGLE ROW OF LOCAL MINIMA  (N = 2)

    funz = 0.5D0 * (p5a*x(1)*x(1) + 1.d0 - COS(2.d0*x(1))) + x(2)*x(2)

  CASE (6)
    !    6  SIX-HUMP CAMEL FUNCTION  (N = 2)

    xx = x(1) * x(1)
    yy = x(2) * x(2)
    funz = ((xx/3.d0 - (2.d0 + p6a))*xx + 4.d0)*xx + x(1)*x(2) + 4.d0*(yy-1.d0)*yy

  CASE (7)
    !    7  TWO-DIMENSIONAL PENALIZED SHUBERT FUNCTION, BETA = 0  (N = 2)

    beta = 0.d0
    GO TO 110

  CASE (8)
    !    8  TWO-DIMENSIONAL PENALIZED SHUBERT FUNCTION, BETA = 1/2  (N = 2)

    beta = 0.5D0
    GO TO 110

  CASE (9)
    !    9  TWO-DIMENSIONAL PENALIZED SHUBERT FUNCTION, BETA = 1  (N = 2)

    beta = 1.d0
    GO TO 110

  CASE (10)
    !   10  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10  (N = 2)

    a = 1.d1
    b = 1.d-1
    GO TO 190

  CASE (11)
    !   11  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**2  (N = 2)

    a = 1.d2
    b = 1.d-2
    GO TO 190

  CASE (12)
    !   12  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**3  (N = 2)

    a = 1.d3
    b = 1.d-3
    GO TO 190

  CASE (13)
    !   13  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**4  (N = 2)

    a = 1.d4
    b = 1.d-4
    GO TO 190

  CASE (14)
    !   14  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**5  (N = 2)

    a = 1.d5
    b = 1.d-5
    GO TO 190

  CASE (15)
    !   15  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**6  (N = 2)

    a = 1.d6
    b = 1.d-6
    GO TO 190

  CASE (16)
    !   16  GOLDSTEIN-PRICE FUNCTION  (N = 2)

    u = x(1) + x(2) + 1.d0
    v = 2.d0 * x(1) - 3.d0 * x(2)
    uu = u * u
    vv = v * v
    funz = (1.d0+uu*(36.d0-20.d0*u+3.d0*uu)) * (30.d0+vv*(18.d0-16.d0*v+3.d0*vv))

  CASE (17)
    !   17  PENALIZED BRANIN FUNCTION  (N = 2)

    funz = (x(2)-p17a*(x(1)/pi)**2+(5.d0/pi)*x(1)-6.d0) ** 2 + 10.d0 *  &
           (1.d0-1.d0/(8.d0*pi)) * COS(x(1)) + 10.d0
    IF (ABS(x(1)-2.5D0) > 7.5D0) funz = funz + penfun(x(1)-2.5D0,7.5D0,100.d0,2)
    IF (ABS(x(2)-7.5D0) > 7.5D0) funz = funz + penfun(x(2)-7.5D0,7.5D0,100.d0,2)

  CASE (18)
    !   18  PENALIZED SHEKEL FUNCTION, M = 5  (N = 4)

    m = 5
    GO TO 250

  CASE (19)
    !   19  PENALIZED SHEKEL FUNCTION, M = 7  (N = 4)

    m = 7
    GO TO 250

  CASE (20)
    !   20  PENALIZED SHEKEL FUNCTION, M = 10  (N = 4)

    m = 10
    GO TO 250

  CASE (21)
    !   21  PENALIZED THREE-DIMENSIONAL HARTMAN FUNCTION  (N = 3)

    m = 4
    funz = 0.d0
    DO  i = 1, m
      s = 0.d0
      DO  j = 1, n
        s = s - p21a(i,j) * (x(j)-p21b(i,j)) ** 2
      END DO
      IF (s >= p21d) funz = funz - p21c(i) * EXP(s)
    END DO
    GO TO 350

  CASE (22)
    !   22  PENALIZED SIX-DIMENSIONAL HARTMAN FUNCTION  (N = 6)

    m = 4
    funz = 0.d0
    DO  i = 1, m
      s = 0.d0
      DO  j = 1, n
        s = s - p22a(i,j) * (x(j)-p22b(i,j)) ** 2
      END DO
      IF (s >= p22d) funz = funz - p22c(i) * EXP(s)
    END DO
    GO TO 350

  CASE (23:25)
    !   23  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 1  (N = 2)

    !   24  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 1  (N = 3)

    !   25  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 1  (N = 4)

    y(1:n) = 1.d0 + 0.25D0 * (x(1:n) - 1.d0)
    GO TO 450

  CASE (26:28)
    !   26  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 2  (N = 5)

    !   27  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 2  (N = 8)

    !   28  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 2  (N = 10)

    y(1:n) = x(1:n)
    GO TO 450

  CASE (29:31)
    !   29  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 10   (N = 2)

    !   30  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 10   (N = 3)

    !   31  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 10   (N = 4)

    range = 10.d0
    GO TO 530

  CASE (32:34)
    !   32  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 5   (N = 5)

    !   33  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 5   (N = 6)

    !   34  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 5  (N = 7)

    range = 5.d0
    GO TO 530

  CASE (35)
    !   35  A FUNCTION WITH A SINGLE CUSP-SHAPED MINIMUM  (N = 5)

    funz = 0.d0
    DO  i = 1, n
      funz = funz + i * x(i) * x(i)
    END DO
    funz = SQRT(SQRT(funz))

  CASE (36)
    !   36  A FUNCTION WITH A SMALL-ATTRACTION-REGION GLOBAL MINIMUM (N = 2)

    rg = p36a
    rp = p36b
    h = p36c
    pund = p36d
    GO TO 610

  CASE (37)
    !   37  A FUNCTION WITH A SMALL-ATTRACTION-REGION GLOBAL MINIMUM (N = 5)

    rg = p37a
    rp = p37b
    h = p37c
    pund = p37d
    GO TO 610

END SELECT
RETURN

110 funz = ((x(1)-p9a)**2 + (x(2)-p9b)**2) * beta
s1 = 0.d0
s2 = 0.d0
DO  i = 1, 5
  di = i
  s1 = s1 + di * COS((di+1.d0)*x(1) + di)
  s2 = s2 + di * COS((di+1.d0)*x(2) + di)
END DO
funz = funz + s1 * s2
IF (ABS(x(1)) > 10.d0) funz = funz + penfun(x(1),10.d0,100.d0,2)
IF (ABS(x(2)) > 10.d0) funz = funz + penfun(x(2),10.d0,100.d0,2)
RETURN

190 xx = x(1) * x(1)
yy = x(2) * x(2)
r2 = xx + yy
r4 = r2 * r2
r8 = r4 * r4
funz = a * xx + yy - r4 + b * r8
RETURN

250 funz = 0.d0
DO  i = 1, m
  s = p20b(i)
  DO  j = 1, n
    s = s + (x(j)-p20a(i,j)) ** 2
  END DO
  funz = funz - 1.d0 / s
END DO
DO  i = 1, n
  IF (ABS(x(i)-5.d0) > 5.d0) funz = funz + penfun(x(i)-5.d0,5.d0,100.d0,2)
END DO
RETURN

350 DO  i = 1, n
  IF (ABS(x(i)-0.5D0) > 0.5D0) funz = funz + penfun(x(i)-0.5D0,0.5D0,100.d0,2)
END DO
RETURN

450 funz = 10.d0 * SIN(pi*y(1)) ** 2 + (y(n)-1.d0) ** 2
DO  i = 2, n
  funz = funz + (y(i-1)-1.d0) ** 2 * (1.d0 + 10.d0*SIN(pi*y(i))**2)
END DO
funz = funz * pi / n
range = 10.d0
GO TO 550

530 funz = SIN(3.d0*pi*x(1))**2 + (x(n) - 1.d0)**2 * (1.d0 + SIN(2.d0*pi*x(n))**2)
DO  i = 2, n
  funz = funz + (x(i-1)-1.d0)**2 * (1.d0 + SIN(3.d0*pi*x(i))**2)
END DO
funz = funz * 0.1D0

550 DO  i = 1, n
  IF (ABS(x(i)) > range) funz = funz + penfun(x(i),range,100.d0,4)
END DO
RETURN

610 s = SUM( x(2:n)**2 )
funz = s + x(1) * x(1)
s = s + (x(1)-rg) ** 2
IF (s < rp*rp*pund) funz = funz - (rg*rg+h) * EXP(-s/(rp*rp-s))
RETURN

CONTAINS


FUNCTION penfun(var, ran, fact, iexp) RESULT(fn_val)
!     PENALIZATION FUNCTION

REAL (dp), INTENT(IN)  :: var, ran, fact
INTEGER, INTENT(IN)    :: iexp
REAL (dp)              :: fn_val

fn_val = fact * (ABS(var) - ran) ** iexp
RETURN
END FUNCTION penfun

END SUBROUTINE glomtf



SUBROUTINE glomip(nprob, n, x0, xmin, xmax)

!  THE SUBROUTINE  GLOMIP  PROVIDES THE CODING FOR THE NUMBER OF VARIABLES,
!  THE INITIAL POINT, AND THE OBSERVATION REGION TO BE USED, TOGETHER WITH
!  THE 37 TEST FUNCTIONS GIVEN BY SUBROUTINE  GLOMTF , TO DEFINE 37 TEST
!  PROBLEMS FOR GLOBAL MINIMIZATION SOFTWARE.

!  THE SUBROUTINE  GLOMIP  RETURNS IN  N , X0 , AND  XMIN ,
!  XMAX  THE NUMBER OF VARIABLES, THE INITIAL POINT, AND THE
!  BOUNDARIES OF THE OBSERVATION REGION.

!  CALLING STATMENT
!     CALL GLOMIP (NPROB, N, X0, XMIN, XMAX)

!  DESCRIPTION OF THE CALL PARAMETERS
!     (THE FORTRAN IMPLICIT TYPE DEFINITION FOR INTEGERS IS USED.
!     ALL NON-INTEGER PARAMETERS ARE DOUBLE-PRECISION).

!  NPROB   IS THE (INPUT) NUMBER OF THE TEST PROBLEM TO BE CONSIDERED
!  N       IS THE (OUTPUT) NUMBER OF VARIABLES (DIMENSION) OF THE PROBLEM
!  XMIN , XMAX   ARE THE (OUTPUT) N-VECTORS CONTAINING THE LEFT AND RIGHT
!          BOUNDARIES OF THE OBSERVATION REGION DEFINED BY THE POINTS
!          X = (X1,...,XN) SUCH THAT XMIN(I) <= X(I) <= XMAX(I), I = 1,...,N

USE toms667, ONLY: dp
IMPLICIT NONE

INTEGER, INTENT(IN)     :: nprob
INTEGER, INTENT(OUT)    :: n
REAL (dp), INTENT(OUT)  :: x0(:)
REAL (dp), INTENT(OUT)  :: xmin(:)
REAL (dp), INTENT(OUT)  :: xmax(:)

REAL (dp) :: v0, vmin, vmax
INTEGER   :: i

SELECT CASE (nprob)

  CASE (1)
    !    1  A FOURTH-ORDER POLYNOMIAL  (N = 1)

    n = 1
    v0 = 1.0D0
    vmin = -10.d0
    vmax = 10.d0
    GO TO 430

  CASE (2)
    !    2  GOLDSTEIN SIXTH ORDER POLYNOMIAL  (N = 1)

    n = 1
    v0 = 0.d0
    vmin = -4.d0
    vmax = 4.d0
    GO TO 430

  CASE (3)
    !    3  ONE-DIMENSIONAL PENALIZED SHUBERT FUNCTION  (N = 1)

    n = 1
    v0 = 0.d0
    vmin = -10.d0
    vmax = 10.d0
    GO TO 430

  CASE (4)
    !    4  A FOURTH ORDER POLYNOMIAL IN TWO VARIABLES  (N = 2)

    n = 2
    x0(1) = 1.d0
    x0(2) = 0.d0
    vmin = -10.d0
    vmax = 10.d0
    GO TO 410

  CASE (5)
    !    5  A FUNCTION WITH A SINGLE ROW OF LOCAL MINIMA  (N = 2)

    n = 2
    x0(1) = -3.d0
    x0(2) = 0.d0
    xmin(1) = -15.d0
    xmax(1) = 25.d0
    xmin(2) = -5.d0
    xmax(2) = 15.d0
    RETURN

  CASE (6)
    !    6  SIX-HUMP CAMEL FUNCTION  (N = 2)

    n = 2
    v0 = 0.d0
    xmin(1) = -3.d0
    xmax(1) = 3.d0
    xmin(2) = -2.d0
    xmax(2) = 2.d0
    GO TO 450

  CASE (7:9)
    !    7  TWO-DIMENSIONAL PENALIZED SHUBERT FUNCTION, BETA = 0  (N = 2)

    !    8  TWO-DIMENSIONAL PENALIZED SHUBERT FUNCTION, BETA = 1/2  (N = 2)

    !    9  TWO-DIMENSIONAL PENALIZED SHUBERT FUNCTION, BETA = 1  (N = 2)

    n = 2
    v0 = 0.d0
    vmin = -10.d0
    vmax = 10.d0
    GO TO 430

  CASE (10:15)
    !   10  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10  (N = 2)

    !   11  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**2  (N = 2)

    !   12  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**3  (N = 2)

    !   13  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**4  (N = 2)

    !   14  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**5  (N = 2)

    !   15  A FUNCTION WITH THREE ILL-CONDITIONED MINIMA, A=10**6  (N = 2)

    n = 2
    x0(1) = 0.d0
    x0(2) = 0.d0
    xmin(1) = -10.d0
    xmax(1) = 10.d0
    xmin(2) = -100.d0
    xmax(2) = 100.d0
    RETURN

  CASE (16)
    !   16  GOLDSTEIN-PRICE FUNCTION  (N = 2)

    n = 2
    v0 = 1.d0
    vmin = -2.d0
    vmax = 2.d0
    GO TO 430

  CASE (17)
    !   17  PENALIZED BRANIN FUNCTION  (N = 2)

    n = 2
    x0(1) = 2.5D0
    x0(2) = 7.5D0
    xmin(1) = -5.d0
    xmax(1) = 10.d0
    xmin(2) = 0.d0
    xmax(2) = 15.d0
    RETURN

  CASE (18:20)
    !   18  PENALIZED SHEKEL FUNCTION, M = 5  (N = 4)

    !   19  PENALIZED SHEKEL FUNCTION, M = 7  (N = 4)

    !   20  PENALIZED SHEKEL FUNCTION, M = 10  (N = 4)

    n = 4
    v0 = 9.d0
    vmin = 0.d0
    vmax = 10.d0
    GO TO 430

  CASE (21)
    !   21  PENALIZED THREE-DIMENSIONAL HARTMAN FUNCTION  (N = 3)

    n = 3
    GO TO 230

  CASE (22)
    !   22  PENALIZED SIX-DIMENSIONAL HARTMAN FUNCTION  (N = 6)

    n = 6
    GO TO 230

  CASE (23)
    !   23  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 1  (N = 2)

    n = 2
    GO TO 330

  CASE (24)
    !   24  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 1  (N = 3)

    n = 3
    GO TO 330

  CASE (25)
    !   25  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 1  (N = 4)

    n = 4
    GO TO 330

  CASE (26)
    !   26  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 2  (N = 5)

    n = 5
    GO TO 330

  CASE (27)
    !   27  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 2  (N = 8)

    n = 8
    GO TO 330

  CASE (28)
    !   28  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 2  (N = 10)

    n = 10
    GO TO 330

  CASE (29)
    !   29  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 10   (N = 2)

    n = 2
    GO TO 330

  CASE (30)
    !   30  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 10   (N = 3)

    n = 3
    GO TO 330

  CASE (31)
    !   31  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 10   (N = 4)

    n = 4
    GO TO 330

  CASE (32)
    !   32  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 5   (N = 5)

    n = 5
    GO TO 370

  CASE (33)
    !   33  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 5   (N = 6)

    n = 6
    GO TO 370

  CASE (34)
    !   34  PENALIZED LEVY-MONTALVO FUNCTION, TYPE 3, RANGE = 5  (N = 7)

    n = 7
    GO TO 370

  CASE (35)
    !   35  A FUNCTION WITH A SINGLE CUSP-SHAPED MINIMUM  (N = 5)

    n = 5
    v0 = 1000.d0
    vmin = -20000.d0
    vmax = 10000.d0
    GO TO 430

  CASE (36)
    !   36  A FUNCTION WITH A SMALL-ATTRACTION-REGION GLOBAL MINIMUM (N = 2)

    n = 2
    x0(1) = 0.d0
    x0(2) = 100.d0
    vmin = -1000.d0
    vmax = 1000.d0
    GO TO 410

  CASE (37)
    !   37  A FUNCTION WITH A SMALL-ATTRACTION-REGION GLOBAL MINIMUM (N = 5)

    n = 5
    x0(1) = 0.d0
    x0(2) = 0.d0
    x0(3) = 0.d0
    x0(4) = 0.d0
    x0(5) = 10.d0
    vmin = -100.d0
    vmax = 100.d0
    GO TO 410

END SELECT

230 v0 = 0.5D0
vmin = 0.d0
vmax = 1.d0
GO TO 430

330 v0 = 0.d0
vmin = -10.d0
vmax = 10.d0
GO TO 430

370 v0 = 0.d0
vmin = -5.d0
vmax = 5.d0
GO TO 430

410 DO  i = 1, n
  xmin(i) = vmin
  xmax(i) = vmax
END DO
RETURN

430 DO  i = 1, n
  xmin(i) = vmin
  xmax(i) = vmax
END DO

450 x0(1:n) = v0

RETURN
END SUBROUTINE glomip



SUBROUTINE ptseg(n,xpfmin,fpfmin,fpfmax,kp,nfev)

!  SAMPLE VERSION OF THE OUTPUT SUBROUTINE  PTSEG
!  (PERFORMS END-OF-SEGMENT OUTPUT)
!  WHICH MUST BY SUPPLIED BY THE USER AND WHICH IS DESCRIBED IN
!  DETAIL WITHIN THE SUBROUTINE  SIGMA.

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: xpfmin(:)
REAL (dp), INTENT(IN)  :: fpfmin
REAL (dp), INTENT(IN)  :: fpfmax
INTEGER, INTENT(IN)    :: kp
INTEGER, INTENT(IN)    :: nfev

WRITE (6,5000) kp, nfev, fpfmin, fpfmax
WRITE (6,5100) xpfmin

RETURN

5000 FORMAT (' OBSERVATION PERIOD  KP= ', i4,  &
             ',   FUNCTION EVALUATIONS  NFEV= ', i7/  &
             '  FINAL BEST, WORST FUNCT. VALUES  FPFMIN= ', g13.5,  &
             ',  FPFMAX= ', g13.5)
5100 FORMAT ('  BEST FINAL POINT  XPFMIN'/(' ', 6G13.5))
END SUBROUTINE ptseg



SUBROUTINE ptrial(n,xopt,fopt,ftfmin,ftfmax,ftfopt,istop,istopt,nfev,kp,iprint)

!  SAMPLE VERSION OF THE OUTPUT SUBROUTINE  PTRIAL
!  (PERFORMS END-OF-TRIAL OUTPUT)
!  WHICH MUST BY SUPPLIED BY THE USER AND WHICH IS DESCRIBED IN
!  DETAIL WITHIN THE SUBROUTINE  SIGMA.

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: xopt(:)
REAL (dp), INTENT(IN)  :: fopt
REAL (dp), INTENT(IN)  :: ftfmin
REAL (dp), INTENT(IN)  :: ftfmax
REAL (dp), INTENT(IN)  :: ftfopt
INTEGER, INTENT(IN)    :: istop
INTEGER, INTENT(IN)    :: istopt
INTEGER, INTENT(IN)    :: nfev
INTEGER, INTENT(IN)    :: kp
INTEGER, INTENT(IN)    :: iprint

WRITE (6,5000)

IF (iprint == 0) WRITE (6,5100) kp, nfev, ftfmin, ftfmax
WRITE (6,5200) istop, istopt, ftfopt, fopt
WRITE (6,5300) xopt(1:n)

RETURN
5000 FORMAT (//' END OF A TRIAL ')
5100 FORMAT ('  OBSERVATION PERIOD  KP= ', i4,  &
             ',   FUNCTION EVALUATIONS  NFEV= ', i7/  &
             '  FINAL BEST, WORST FUNCT. VALUES  FTFMIN= ', g13.5,  &
             ',  FTFMAX= ', g13.5)
5200 FORMAT (/'  TRIAL STOP INDICATOR  ISTOP= ', i2,  &
             ',   PAST STOPS INDICATOR  ISTOPT= ', i2/  &
             '  END-OF-TRIAL BEST FUNCTION VALUE   FTFOPT= ', g14.6/  &
             '  BEST CURRENT MINIMUM FUNCTION VALUE  FOPT= ', g14.6)
5300 FORMAT ('  BEST CURRENT MINIMIZER   XOPT'/(' ', 6G13.5))
END SUBROUTINE ptrial



SUBROUTINE ptksuc(ksuc)

!  SAMPLE VERSION OF THE OUTPUT SUBROUTINE  PTKSUC
!  (PERFORMS END-OF-TRIAL OUTPUT RELATED TO THE COUNT OF
!    SUCCESSFUL TRIALS)
!  WHICH MUST BY SUPPLIED BY THE USER AND WHICH IS DESCRIBED IN
!  DETAIL WITHIN THE SUBROUTINE  SIGMA.

IMPLICIT NONE

INTEGER, INTENT(IN)    :: ksuc

WRITE (6,5000) ksuc, ksuc

RETURN
5000 FORMAT (// ' THE CURRENT COUNT  KSUC  OF SUCCESSFUL TRIALS HAS ',  &
             'REACHED FOR THE FIRST'/  &
             '  TIME THE VALUE ', i2,  &
             ' (IF THE REQUESTED COUNT  NSUC  OF SUCCESSFUL TRIALS '/  &
             '  HAD BEEN GIVEN THE VALUE ', i2, ' THE ALGORITHM WOULD ',  &
             'HAVE STOPPED HERE)'//)
END SUBROUTINE ptksuc
