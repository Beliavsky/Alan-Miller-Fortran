!   PROGRAM NAME    - FIT
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-27  Time: 13:16:33


! N.B. The files INPUT & OUTPUT follow at the end of the code.

!-----------------------------------------------------------------------

!   LATEST REVISION - DECEMBER 27, 1986

!   PURPOSE     - PARAMETER ESTIMATION OF RESPIRATORY MODELS
!                 BY A GLOBAL OPTIMIZATION ALGORITHM

!   REQUIRED ROUTINES   - GLOBAL, LOCAL, URDMN, FUNCT, FUN

!-----------------------------------------------------------------------

MODULE fit_common
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

! COMMON /fu/ om(50), rezz(50), irel, zre(50), zim(50)

REAL (dp), SAVE :: om(50), rezz(50), zre(50), zim(50)
INTEGER, SAVE   :: irel

END MODULE fit_common



PROGRAM fit
USE fit_common
USE global_minimum
IMPLICIT NONE

REAL (dp)             :: x0(15,20), f00(20), MIN(18), MAX(18),   &
                         df, do, f0, o0, o1, tpi
CHARACTER (LEN=80)    :: label
INTEGER               :: i, in, ipr, k, m, nc, ng, nparm, npt, nsampl, nsig
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp, onep5 = 1.5_dp,  &
                         eight = 8.0_dp, tausend = 1000._dp
INTEGER, ALLOCATABLE  :: seed(:)

INTERFACE
  SUBROUTINE ladder(x, f, f0, df, o0, DO, npt, ipr, np)
    USE fit_common
    IMPLICIT NONE
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN)      :: f
    REAL (dp), INTENT(IN)      :: f0
    REAL (dp), INTENT(IN)      :: df
    REAL (dp), INTENT(IN)      :: o0
    REAL (dp), INTENT(IN OUT)  :: DO
    INTEGER, INTENT(IN OUT)    :: npt
    INTEGER, INTENT(IN OUT)    :: ipr
    INTEGER, INTENT(IN OUT)    :: np
  END SUBROUTINE ladder
END INTERFACE

tpi = eight*ATAN(one)
in = 5
ipr = 6
OPEN(5, FILE='INPUT')
OPEN(6, FILE='OUTPUT')
READ(in, 901) label
901 FORMAT(A)
WRITE(ipr, 902) label
WRITE(*, 902) label
902 FORMAT(' ', A/)
nparm = 5
READ(in, *) MIN(1), MAX(1), MIN(2), MAX(2), MIN(3), MAX(3), MIN(4), MAX(4), &
            MIN(5), MAX(5)
920 FORMAT(' ', f9.4, ' ', f9.4)
WRITE(ipr, 920) MIN(1), MAX(1), MIN(2), MAX(2), MIN(3), MAX(3),  &
                MIN(4), MAX(4), MIN(5), MAX(5)
WRITE(*, 920) MIN(1), MAX(1), MIN(2), MAX(2), MIN(3), MAX(3),  &
              MIN(4), MAX(4), MIN(5), MAX(5)
WRITE(ipr, 905)
WRITE(*, 905)
905 FORMAT(// ' RUN PARAMETERS'/)
READ(in, *) irel, nsampl, ng, nsig
WRITE(ipr, 922) irel, nsampl, ng, nsig
WRITE(*, 922) irel, nsampl, ng, nsig
922 FORMAT(' ', i2, ' ', i4, ' ', i2, ' ', i1 /)
i = 1
WRITE(ipr, 906)
WRITE(*, 906)
906 FORMAT(' SAMPLE'/)

3 READ(in, *, END=6, ERR=9) om(i), zre(i), zim(i)
WRITE(ipr, 925) om(i), zre(i), zim(i)
WRITE(*, 925) om(i), zre(i), zim(i)
925 FORMAT(' ', f6.3, ' ', f7.2, ' ', f8.2)
i = i + 1
GO TO 3

6 m = i - 1
WRITE(ipr, 926)
WRITE(*, 926)
926 FORMAT(//)
o0 = tausend
o1 = zero
DO = (om(2) - om(1))*tpi
DO  i=1,m
  om(i) = om(i)*tpi
  IF (om(i) < o0) o0 = om(i)
  IF (om(i) > o1) o1 = om(i)
  IF (i <= 1) CYCLE
  IF (om(i) - om(i-1) < DO) DO = om(i) - om(i-1)
END DO

! Set the random number seed

CALL RANDOM_SEED(size=k)
ALLOCATE (seed(k))
WRITE(*, '(1x, a, i4, a)') 'Enter ', k, ' integers as random number seeds: '
READ(*, *) seed
CALL RANDOM_SEED(put=seed)
WRITE(ipr, '(a / (7i11) )') ' Random number seed(s): ', seed
WRITE(ipr, * )

f0 = o0/tpi
df = DO/tpi
npt = INT((o1-o0)/DO + onep5)
DO  i=1,m
  rezz(i) = SQRT(zre(i)*zre(i) + zim(i)*zim(i))
END DO
CALL global(MIN, MAX, nparm, m, nsampl, ng, ipr, nsig, x0, nc, f00)
DO  i=1,nc
  CALL ladder(x0(1:,i), f00(i), f0, df, o0, DO, npt, ipr, nparm)
END DO

9 STOP
END PROGRAM fit
    

SUBROUTINE funct(x, value, nparm, m)

USE fit_common
IMPLICIT NONE

REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: value
INTEGER, INTENT(IN)     :: nparm
INTEGER, INTENT(IN)     :: m

REAL (dp) :: f(100), zimi, zrei
INTEGER   :: i, j, kk, mm

DO  kk=1,m
  zrei = x(1) + (x(2)/(om(kk)**x(3)))
  zimi = om(kk)*x(4) - (x(5)/(om(kk)**x(3)))
  zrei = zre(kk) - zrei
  zimi = zim(kk) - zimi
  j = kk*2
  i = j-1
  IF (irel /= 0) GO TO 100
  f(i) = zrei/rezz(kk)
  f(j) = zimi/rezz(kk)
  CYCLE
  100 f(i) = zrei
  f(j) = zimi
END DO

mm = m+m
value = SUM( f(1:mm)**2 )
value = SQRT(value/m)
RETURN
END SUBROUTINE funct


SUBROUTINE ladder(x, f, f0, df, o0, DO, npt, ipr, np)

USE fit_common
IMPLICIT NONE

REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN)      :: f
REAL (dp), INTENT(IN)      :: f0
REAL (dp), INTENT(IN)      :: df
REAL (dp), INTENT(IN)      :: o0
REAL (dp), INTENT(IN OUT)  :: DO
INTEGER, INTENT(IN OUT)    :: npt
INTEGER, INTENT(IN OUT)    :: ipr
INTEGER, INTENT(IN OUT)    :: np

REAL (dp) :: zbr(50), zbi(50), fok, oj, tr1, tr2, tr3, tr4, tr5, tr6, tr7
INTEGER   :: i, kk
REAL (dp), PARAMETER :: zero = 0.0_dp, one = 1.0_dp, eight = 8.0_dp,  &
                        ts = 360._dp

oj = o0 - DO
fok = ts/(eight*ATAN(one))
WRITE(ipr, 901) f, x(1:np)
901 FORMAT(/////' ', g14.8, 3(/'    ', 5(g14.8, ' ')))
WRITE(ipr, 903)
903 FORMAT(//'  FREQ    REAL       IMAG       ABS       PHASE    DRE ',  &
             '   DIM   DABS    FREQ'/)
DO  kk=1,npt
  oj = oj + DO
  zbr(kk) = x(1) + (x(2)/(oj**x(3)))
  zbi(kk) = oj*x(4) - (x(5)/(oj**x(3)))
END DO
oj = f0 - df
DO  i=1,npt
  oj = oj + df
  tr1 = zbr(i)
  tr2 = zbi(i)
  tr3 = SQRT(tr1*tr1 + tr2*tr2)
  tr4 = ATAN2(tr2, tr1)*fok
  IF (zre(i) /= zero) tr5 = (zre(i)-tr1)/ABS(zre(i))
  IF (zre(i) == zero) tr5 = zre(i)-tr1
  IF (zim(i) /= zero) tr6 = (zim(i)-tr2)/ABS(zim(i))
  IF (zim(i) == zero) tr6 = zim(i)-tr2
  tr7 = SQRT(tr5*tr5 + tr6*tr6)
  WRITE(ipr, 904) oj, tr1, tr2, tr3, tr4, tr5, tr6, tr7, oj
  904 FORMAT(' ', f5.2, ' ', 3(' ', g10.4), 4(' ', f6.2), '   ', f5.2)
END DO
RETURN
END SUBROUTINE ladder


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The file: INPUT

 ez itt a cimke helye
   0.0000   1.0000
   0.0000   1.0000
   0.0000   2.0000
   0.0000   1.0000
   0.0000   1.0000
 0  100  2  3
  0.025  5.00   -5.00
  0.050  3.00   -2.00
  0.075  2.00   -1.00
  0.100  1.50   -0.50
  0.125  1.20   -0.20
  0.150  1.10   -0.10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! An output file follows:

  ez itt a cimke helye

    0.0000    1.0000
    0.0000    1.0000
    0.0000    2.0000
    0.0000    1.0000
    0.0000    1.0000


 RUN PARAMETERS

  0  100  2 3

 SAMPLE

  0.025    5.00    -5.00
  0.050    3.00    -2.00
  0.075    2.00    -1.00
  0.100    1.50    -0.50
  0.125    1.20    -0.20
  0.150    1.10    -0.10



 Random number seed(s): 
      87632     389707         32     643890      32987    4365908     243987
        689      43897     654987


   100 FUNCTION EVALUATIONS USED FOR SAMPLING
 *** THE LOCAL MINIMUM NO.  1: 0.67609502E-01, NFEV=  200
 0.67609502E-01
    0.51391532     0.55396685      1.1950652     0.56849758     0.53076047     
    
 SAMPLE POINT ADDED TO THE CLUSTER NO.  1
 0.21660555    
    0.45845671     0.51020943      1.3761022     0.59273621     0.57636474     
    

   100 FUNCTION EVALUATIONS USED FOR SAMPLING
 NEW SEED POINT ADDED TO THE CLUSTER NO.  1, NFEV=  163
 0.67609502E-01
    0.10385758     0.85375662      1.0035828     0.90094354     0.75231245     
    
 *** IMPROVEMENT ON THE LOCAL MINIMUM NO.  1:0.67609502E-01 FOR 0.67609502E-01
 0.67609502E-01
    0.51391543     0.55396676      1.1950652     0.56849753     0.53076042     
    
 SAMPLE POINT ADDED TO THE CLUSTER NO.  1
 0.16335736    
    0.14014181     0.93676996     0.78321739     0.85749589     0.73506345     
    
 NEW SEED POINT ADDED TO THE CLUSTER NO.  1, NFEV=  184
 0.67609502E-01
    0.71011390     0.62601696      1.0769925     0.95468350     0.69247156     
    
 *** IMPROVEMENT ON THE LOCAL MINIMUM NO.  1:0.67609502E-01 FOR 0.67609502E-01
 0.67609502E-01
    0.51391535     0.55396682      1.1950652     0.56849757     0.53076046     
    





 LOCAL MINIMA FOUND:


 0.67609502E-01
    0.51391535     0.55396682      1.1950652     0.56849757     0.53076046     
    



 NORMAL TERMINATION AFTER   747 FUNCTION EVALUATIONS








 0.67609502E-01
    0.51391535     0.55396682      1.1950652     0.56849757     0.53076046     
    


  FREQ    REAL       IMAG       ABS       PHASE    DRE    DIM   DABS    FREQ

  0.03   5.574     -4.759      7.329     -40.49  -0.11  -0.05   0.12    0.03
  0.05   2.724     -1.939      3.344     -35.44   0.09  -0.03   0.10    0.05
  0.08   1.875     -1.036      2.143     -28.93   0.06   0.04   0.07    0.08
  0.10   1.479     -.5677      1.584     -21.00   0.01   0.14   0.14    0.10
  0.13   1.253     -.2619      1.280     -11.80  -0.04   0.31   0.31    0.13
  0.15   1.109     -.3390E-01  1.109      -1.75  -0.01  -0.66   0.66    0.15
