MODULE Kernel_Regression

! General Remarks
! The subroutines glkern.f and lokern.f use an efficient and fast algorithm for
! automatically adaptive nonparametric regression estimation with a kernel method.
! Roughly speaking, the method performs a local averaging of the observations when
! estimating the regression function.  Analogously, one can estimate derivatives
! of small order of the regression function.

! Crucial for the kernel regression estimation used here is the choice of global
! or local bandwidths. Too small ones will lead to a wiggly curve, too large ones
! will smooth away important details.

! The subroutine glkern.f calculates an estimator of the regression function or
! derivatives of the regression function with an automatically chosen global
! plugin bandwidth.  The subroutine lokern.f calculates such an estimator with an
! automatically chosen bandwidth function.  It is also possible to use global and
! local bandwidths, respectively, which are specified by the user.
! The main idea of the plugin method is to estimate the optimal bandwidths by
! estimating the asymptotically optimal mean (integrated) squared error optimal
! bandwidths.  Therefore, one has to estimate the variance for homoscedastic error
! variables and a functional of a smooth variance function for heteroscedastic
! error variables, respectively.  Also, one has to estimate an integral functional
! of the squared k-th derivative of the regression function (k=KORD) for the
! global bandwidth and the squared k-th derivative itself for the local
! bandwidths.  Here, a further kernel estimator for this derivative is used with a
! bandwidth which is adapted iteratively to the regression function.
! A convolution form of the kernel estimator for the regression function and its
! derivatives is used.  Thereby, one can adapt the S-array (which is
! S(I)=(T(I)+T(I+1))/2 in the standard convolution form) to be a smoothed grid
! more suitable for random design, see Herrmann (1996). Using this estimator leads
! to an asymptotically minimax efficient estimator for fixed and random design.
! Polynomial kernels and boundary kernels are used with a fast and stable updating
! algorithm for kernel regression estimation.

! More details can be found in the papers referred to in the references.

! References
! On the global iterative plugin bandwidth estimator:
! T. Gasser, A. Kneip, and W. Köhler (1991). A flexible and fast method for
! automatic smoothing. Journal of the American Statistical Association, 86,
! 643-652.
! On the local plugin bandwidth estimator:
! M. Brockmann, T. Gasser, and E. Herrmann (1993). Locally adaptive bandwidth
! choice for kernel regression estimators. Journal of the American Statistical
! Association, 88, 1302-1309.
! On nonparametric variance estimation:
! T. Gasser, L. Sroka, and C. Jennen-Steinmetz (1986). Residual and residual
! pattern in nonlinear regression. Biometrika, 73, 625-633.
! On adapting heteroscedasticity:
! E. Herrmann (1997). Local bandwidth choice in kernel regression estimation.
! Journal of Graphical and Computational Statistics, 6, 35-54.
! On the fast algorithm for kernel regression estimator:
! T. Gasser and A. Kneip (1989) discussion of Buja, A., Hastie, T. and Tibshirani,
! R.: Linear smoothers and additive models, The Annals of Statistics, 17, 532-535.
!
! B. Seifert, M. Brockmann, J. Engel, and T. Gasser (1994). Fast algorithms for
! nonparametric curve estimation. J. Computational and Graphical Statistics 3,
! 192-213.
! On the special kernel estimator for random design:
! E. Herrmann (1996). On the convolution type kernel regression estimator.
! Preprint 1833, FB Mathematik, Technische Universität Darmstadt (available from
! the preprint server of the mathematical department of the technical university
! of Darmstadt )

! Comments to Eva Herrmann: eherrmann@mathematik.tu-darmstadt.de

! Last update (Fortran 77): 10-December-98 / mz

! Code converted using TO_F90 by Alan Miller
! Date: 2000-10-04  Time: 16:41:52

IMPLICIT NONE

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


SUBROUTINE glkern(t, x, n, tt, m, ihom, nue, kord, irnd,  &
                  ismo, m1, tl, tu, s, sig, b, y)

! N.B. Arguments WN & W1 have been removed.
 
!------------------------------------------------------------------*
!   SHORT-VERSION: OCT 1996

!   PURPOSE:

!   GENERAL SUBROUTINE FOR KERNEL SMOOTHING:
!   COMPUTATION OF ITERATIVE PLUG-IN ALGORITHM FOR GLOBAL BANDWIDTH
!   SELECTION FOR KERNELS WITH (NUE,KORD) = (0,2),(0,4),(1,3) OR (2,4).

!------------------------------------------------------------------*
!   THE RAW DATA SHOULD BE GIVEN BY THE POINTS
!   (T(1),X(1)),...,(T(N),X(N))

!   THE RESULTING ESTIMATOR OF THE NUE-TH DERIVATIVE OF THE
!   REGRESSION CURVE IS GIVEN THROUGH THE POINTS
!   (TT(1),Y(1)),...,(TT(M),Y(M))

!   THE PLUG-IN BANDWIDTH IS GIVEN BY B
!------------------------------------------------------------------*

!  PARAMETERS :

!  INPUT    T(N)         INPUT GRID (T(1)<T(2)<...<T(N))
!  INPUT    X(N)         DATA
!  INPUT    N            LENGTH OF X

!  INPUT    TT(M)        OUTPUT GRID, SHOULD BE ORDERED
!  INPUT    M            LENGTH OF TT

!  INPUT    IHOM         HOMOSKEDASTICY OF VARIANCE
!                        0: HOMOSKEDASTIC ERROR VARIABLES,
!                        <> 0: IF THE VARIANCE SHOULD ESTIMATED AS
!                              SMOOTH FUNCTION.
!                        ****** DEFAULT VALUE: IHOM=0

!  INPUT    NUE          ORDER OF DERIVATIVE (0-4) OF THE REGRESSION
!                        FUNCTION WHICH SHALL BE ESTIMATED
!                        ****** DEFAULT VALUE: NUE=0

!  INPUT    KORD         ORDER OF KERNEL (<=6), FOR ISMO=0  ONLY
!                        NUE=0, KORD=2 OR KORD=4
!                        OR NUE=1, KORD=3 OR NUE=2, KORD=4 ARE ALLOWED
!                        ****** DEFAULT VALUE: KORD=NUE+2

!  INPUT    IRND         0: IF RANDOM GRID POINTS T MAY OCCUR
!                        <>0 ELSE  (ONLY NECESSARY IF S SHOULD BE
!                            COMPUTED)
!                        ****** DEFAULT VALUE IRND=0

!  INPUT    ISMO         0:ESTIMATING THE OPTIMAL GLOBAL BANDWIDTH
!                        <>0 USING GLOBAL INPUT BANDWIDTH B
!                        ****** DEFAULT VALUE ISMO=0

!  INPUT    M1           >=10, LENGTH OF W1, LARGE VALUES WILL INCREASE
!                        THE ACCURACY OF THE INTEGRAL APPROXIMATION
!                        ****** DEFAULT VALUE: M1=400

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! IN/OUTPUT TL/TU        LOWER/UPPER BOUND FOR INTEGRAL APPROXIMATION
!                        AND VARIANCE ESTIMATION (IF SIG=0 AND IHOM=0),
!                        IF TU<=TL, [TL,TU] ARE COMPUTED AS ABOUT
!                        THE 87% MIDDLE PART OF [T(1),T(N)]
!                        ****** DEFAULT VALUES: TL=1.0, TU=0.0

! IN/OUTPUT S(0:N)       IF S(N)<=S(0) THIS ARRAY IS COMPUTED AS
!                        MIDPOINTS OF T, FOR NON-RANDOM DESIGN AND AS
!                        SMOOTHED QUANTILES FOR RANDOM DESIGN
!                        ****** DEFAULT VALUES: S(0)=1.0, S(N)=0.0
!                               AND THE OTHER S(I) CAN BE UNDEFINED

! IN/OUTPUT SIG          RESIDUAL VARIANCE, ESTIMATED FOR SIG=0 OR
!                        IHOM<>0, ELSE GIVEN BY INPUT
!                        ****** DEFAULT VALUE: SIG=-1.0

! IN/OUTPUT B            GLOBAL PLUG-IN BANDWIDTH
!                        ****** B CAN BE UNDIFINED IF ISMO=0


! WORK     WN(0:N,5)     WORK ARRAY FOR KERNEL SMOOTHING ROUTINE
!                        OR NUE=1, KORD=3 OR NUE=2, KORD=4 ARE ALLOWED
!                        ****** WILL BE SET IN SUBROUTINE
! WORK     W1(M1,3)      WORK ARRAY FOR INTEGRAL APPROXIMATION
!                        ****** WILL BE SET IN SUBROUTINE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! OUTPUT   Y(M)          KERNEL ESTIMATE WITH BOP (=B0 FOR ISMO<>0)
!                        ****** WILL BE SET IN SUBROUTINE
!-----------------------------------------------------------------------
!  USED SUBROUTINES: COFF, RESEST, KERNEL WITH FURTHER SUBROUTINES
!                    WHICH ARE CONTAINED IN THE FILE subs.f
!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: tt(:)
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: ihom
INTEGER, INTENT(IN)        :: nue
INTEGER, INTENT(IN OUT)    :: kord
INTEGER, INTENT(IN)        :: irnd
INTEGER, INTENT(IN OUT)    :: ismo
INTEGER, INTENT(IN)        :: m1
REAL (dp), INTENT(IN OUT)  :: tl
REAL (dp), INTENT(IN OUT)  :: tu
REAL (dp), INTENT(IN OUT)  :: s(0:)
REAL (dp), INTENT(IN OUT)  :: sig
REAL (dp), INTENT(IN OUT)  :: b
REAL (dp), INTENT(OUT)     :: y(:)

! Local variables

REAL (dp)  :: w1(m1,3), wn(0:n,5)
INTEGER    :: i, ii, iil, il, inputs, iprint, isort, it, itende, itt, iu,  &
              j, kk, kk2, nn, nyg
REAL (dp)  :: alpha, b2, bmax, bmin, bres, bs, const, ex, exs, exsvi, fac,  &
              q, r2, rvar, s0, sn, snr, ssi, tll, tuu, vi, xi, xmy2
!-
!-------- 1. INITIALISATIONS AND SOME ERROR-CHECKS
REAL (dp), SAVE  :: bias(2,0:2) = RESHAPE(  &
                    (/ 0.2, 0.04762, 0.4286, 0.1515, 1.33, 0.6293 /), (/ 2,3 /))
REAL (dp), SAVE  :: vark(2,0:2) = RESHAPE(  &
                    (/ 0.6, 1.250, 2.143, 11.93, 35.0, 381.6 /), (/ 2,3 /))
REAL (dp), SAVE  :: fak2(2:4) = (/ 4., 36., 576. /)

nyg=0
inputs=0
!-------- IF NO ERRORS SHOULD BE WRITTEN ON STANDARD OUTPUT, SET IPRINT=1
!-------- IF ERRORS AND VERY DETAILED WARNINGS SHOULD BE WRITTEN ON
!--------           STANDARD OUTPUT, SET IPRINT < 0
iprint=0

IF(nue > 4 .OR. nue < 0) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: Order of derivative not allowed'
  STOP
END IF
IF(nue > 2 .AND. ismo == 0) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: Order of derivative not allowed'
  STOP
END IF
IF(n <= 2) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: Number of data too small'
  STOP
END IF
IF(m < 1) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: No output points'
  STOP
END IF
IF(m1 < 10) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: Variable M1 is choosen too small'
  STOP
END IF

kk=(kord-nue)/2
IF(2*kk+nue /= kord) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: Kernel order not allowed, set to ',nue+2
  kord=nue+2
END IF
IF(kord > 4 .AND. ismo == 0) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: Kernel order not allowed, set to ',nue+2
  kord=nue+2
END IF
IF(kord > 6 .OR. kord <= nue) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: Kernel order not allowed, set to ',nue+2
  kord=nue+2
END IF
IF(ismo /= 0 .AND. b <= 0) THEN
  IF(iprint == 0) WRITE(*, *) 'glkern: Plug-in bandwidth is used'
  ismo=0
END IF
rvar=sig
!-
!-------- 2. COMPUTATION OF S-SEQUENCE
s0=1.5*t(1) - 0.5*t(2)
sn=1.5*t(n) - 0.5*t(n-1)
IF(s(n) <= s(0)) THEN
  inputs=1
  DO  i=1,n-1
    s(i)=.5*(t(i)+t(i+1))
  END DO
  s(0)=s0
  s(n)=sn
  IF(ismo /= 0 .AND. irnd /= 0) GO TO 160
ELSE
  IF(ismo /= 0) GO TO 160
END IF
!-
!-------- 3. COMPUTATION OF MINIMAL, MAXIMAL ALLOWED BANDWIDTH
bmax=(sn-s0)*.5
bmin=(sn-s0)/DBLE(n)*DBLE(kord-1)*.6
!-
!-------- 4. WARNINGS IF TT-GRID LARGER THAN T-GRID
IF(tt(1) < s0 .AND. tt(m) > sn .AND. iprint < 0) WRITE(*, *)   &
    'glkern: Extrapolation at both boundaries not optimized'
IF(tt(1) < s0 .AND. tt(m) <= sn .AND. iprint < 0) WRITE(*, *)   &
    'glkern: Extrapolation at left boundary not optimized'
IF(tt(1) >= s0 .AND. tt(m) > sn .AND. iprint < 0) WRITE(*, *)   &
    'glkern: Extrapolation at right boundary not optimized'

!-
!-------- 5. COMPUTE TL,TU AND THEIR T-GRID AS INNER PART FOR
!            INTEGRAL APPROXIMATION IN THE ITERATIONS
itt=0
51 IF (tu <= tl) THEN
  tl=.933*s0 +.067*sn
  tu=.067*s0 +.933*sn
  itt=itt+1
END IF
tl=MAX(tl,s0)
tu=MIN(tu,sn)
il=1
iu=n
wn(1,1)=0.0
wn(n,1)=0.0
DO  i=1,n
  IF(t(i) <= tl .OR. t(i) >= tu) wn(i,1)=0.0
  IF(t(i) > tl .AND. t(i) < tu) wn(i,1)=1.0
  IF(t(i) < tl) il=i+1
  IF(t(i) <= tu) iu=i
END DO
nn=iu-il+1
IF(nn == 0 .AND. itt == 0) THEN
  tu=tl-1.0
  GO TO 51
END IF
IF(nn == 0 .AND. itt == 1) THEN
  tu=sn
  tl=s0
  GO TO 51
END IF
!-
!-------- 6. COMPUTE T-GRID FOR INTEGRAL APPROXIMATION
DO  i=1,m1
  w1(i,2)=1.0
  w1(i,1)=tl+(tu-tl)*DBLE(i-1)/DBLE(m1-1)
END DO
!-
!-------- 7. CALCULATION OF WEIGHT FUNCTION
alpha=1.d0/DBLE(13)
DO  i=il,iu
  xi=(t(i) - tl)/alpha/(tu-tl)
  IF(xi > 1) GO TO 71
  wn(i,1)=(10.0 - 15*xi + 6*xi*xi)*xi*xi*xi
END DO
71 DO  i=iu,il,-1
  xi=(tu-t(i))/alpha/(tu-tl)
  IF(xi > 1) GO TO 73
  wn(i,1)=(10.0 - 15*xi + 6*xi*xi)*xi*xi*xi
END DO
73 DO  i=1,m1
  xi=(w1(i,1)-tl)/alpha/(tu-tl)
  IF(xi > 1) GO TO 75
  w1(i,2)=(10.0 - 15*xi + 6*xi*xi)*xi*xi*xi
END DO
75 DO  i=m1,1,-1
  xi=(tu-w1(i,1))/alpha/(tu-tl)
  IF(xi > 1) GO TO 77
  w1(i,2)=(10.0 - 15*xi + 6*xi*xi)*xi*xi*xi
END DO
!-
!-------- 8. COMPUTE CONSTANTS FOR ITERATION
77 ex=1./DBLE(kord+kord+1)
kk2=(kord-nue)
kk=kk2/2
!-
!-------- 9. ESTIMATING VARIANCE AND SMOOTHED PSEUDORESIDUALS
IF(sig <= .0 .AND. ihom == 0) CALL resest(t(il:),x(il:),nn,wn(il:,2),r2,sig)
IF(ihom /= 0) THEN
  CALL resest(t,x,n,wn(1:,2),snr,sig)
  bres=MAX(bmin, 0.2*nn**(-.2)*(s(iu)-s(il-1)))
  DO  i=1,n
    wn(i,3)=t(i)
    wn(i,2)=wn(i,2)*wn(i,2)
  END DO
  CALL kernel(t,wn(1:,2),n,bres,0,kk2,nyg,s,wn(il:,3),nn,wn(il:,4))
ELSE
  CALL coff(wn(1:,4),n,sig)
END IF
!-
!-------- 10. ESTIMATE/COMPUTE INTEGRAL CONSTANT
100 vi=0.
DO  i=il,iu
  vi=vi + wn(i,1)*n*(s(i)-s(i-1))**2*wn(i,4)
END DO
!-
!-------- 11. REFINEMENT OF S-SEQUENCE FOR RANDOM DESIGN
IF(inputs == 1 .AND. irnd == 0) THEN
  DO  i=0,n
    wn(i,5)=DBLE(i)/DBLE(n+1)
    wn(i,2)=(DBLE(i)+.5)/DBLE(n+1)
    wn(i,3)=wn(i,2)
  END DO
  exs=-DBLE(3*kord+1)/DBLE(6*kord+3)
  exsvi=DBLE(kord)/DBLE(6*kord+3)
  bs=0.1*(vi/(sn-s0)**2)**exsvi*n**exs
  CALL kernel(wn(1:,5),t,n,bs,0,2,nyg,wn(0:,3),wn(0:,2),n+1,s(0:))
  111 isort=0
  vi=0.0
  DO  i=1,n
    vi=vi + wn(i,1)*n*(s(i)-s(i-1))**2*wn(i,4)
    IF(s(i) < s(i-1)) THEN
      ssi=s(i-1)
      s(i-1)=s(i)
      s(i)=ssi
      isort=1
    END IF
  END DO
  IF(isort == 1) GO TO 111
  IF(ismo /= 0) GO TO 160
END IF
b=bmin*2.
!-
!-------- 12. COMPUTE INFLATION CONSTANT AND EXPONENT AND LOOP OF ITERATIONS
const=DBLE(2*nue+1)*fak2(kord)*vark(kk,nue)*vi  &
      /(DBLE(2*kord-2*nue)*bias(kk,nue)**2*DBLE(n))
fac=1.1*(1.+(nue/10.)+0.05*(kord-nue-2.)) *n**(2./DBLE((2*kord+1)*(2*kord+3)))
itende=1 + 2*kord + kord*(2*kord+1)

DO  it=1,itende
!-
!-------- 13. ESTIMATE DERIVATIVE OF ORDER KORD IN ITERATIONS
  b2=b*fac
  b2=MAX(b2,bmin/DBLE(kord-1)*DBLE(kord+1))
  b2=MIN(b2,bmax)
  CALL kernel(t,x,n,b2,kord,kord+2,nyg,s,w1(1:,1),m1,w1(1:,3))
!-
!-------- 14. ESTIMATE INTEGRALFUNCTIONAL IN ITERATIONS
  xmy2=.75*(w1(1,2)*w1(1,3)*w1(1,3) + w1(m1,2)*w1(m1,3)*w1(m1,3))
  DO  i=2,m1-1
    xmy2=xmy2 + w1(i,2)*w1(i,3)*w1(i,3)
  END DO
  xmy2=xmy2*(tu-tl)/m1
!-
!-------- 15. FINISH OF ITERATIONS
  b=(const/xmy2)**ex
  b=MAX(bmin,b)
  b=MIN(bmax,b)
END DO
!-
!-------- 16  COMPUTE SMOOTHED FUNCTION WITH PLUG-IN BANDWIDTH
160 CALL kernel(t, x, n, b, nue, kord, nyg, s, tt, m, y)
!-
!-------- 17. VARIANCE CHECK
IF(ihom /= 0) sig=rvar
IF(rvar == sig .OR. r2 < .88 .OR. ihom /= 0 .OR. nue > 0) RETURN
ii=0
iil=0
j=2
tll=MAX(tl, tt(1))
tuu=MIN(tu, tt(m))
DO  i=il,iu
  IF(t(i) < tll .OR. t(i) > tuu) CYCLE
  ii=ii+1
  IF(iil == 0) iil=i
  171 IF(tt(j) < t(i)) THEN
    j=j+1
    IF(j <= m) GO TO 171
  END IF
  wn(ii,3)=x(i)-y(j) + (y(j)-y(j-1))*(tt(j)-t(i))/(tt(j)-tt(j-1))
END DO
IF(iil == 0 .OR. ii-iil < 10) THEN
  CALL resest(t(il:),wn(1:,3),nn,wn(1:,4),snr,rvar)
ELSE
  CALL resest(t(iil:), wn(1:,3), ii, wn(1:,4), snr, rvar)
END IF
q=sig/rvar
IF(q <= 2.) RETURN
IF(q > 5. .AND. r2 > .95) rvar=rvar*.5
sig=rvar
CALL coff(wn(1:,4), n, sig)
GO TO 100
END SUBROUTINE glkern



SUBROUTINE lokern(t, x, n, tt, m, ihom, nue, kord, irnd,  &
                  ismo, m1, tl, tu, s, sig, ban, y)

! N.B. Arguments WN, W1 & WM have been removed.
 
!-------------------------------------------------------------------*
!--------------------------------------------------------------------
!    SHORT-VERSION: JANUARY 1997

!    PURPOSE:

!    GENERAL SUBROUTINE FOR KERNEL SMOOTHING:
!    COMPUTATION OF ITERATIVE PLUG-IN ALGORITHM FOR LOCAL BANDWIDTH
!    SELECTION FOR KERNELS WITH (NUE,KORD) = (0,2),(0,4),(1,3) OR (2,4).

!-------------------------------------------------------------------*
!    THE RAW DATA SHOULD BE GIVEN BY THE POINTS
!    (T(1),X(1)),...,(T(N),X(N))

!    THE RESULTING ESTIMATOR OF THE NUE-TH DERIVATIVE OF THE
!    REGRESSION CURVE IS GIVEN THROUGH THE POINTS
!    (TT(1),Y(1)),...,(TT(M),Y(M))

!    THE LOCAL PLUG-IN BANDWIDTH ARRAY IS GIVEN BY BAN(1),...,BAN(M)
!-------------------------------------------------------------------*

!  PARAMETERS :

!  INPUT    T(N)         INPUT GRID (T(1)<T(2)<...<T(N))
!  INPUT    X(N)         DATA
!  INPUT    N            LENGTH OF X

!  INPUT    TT(M)        OUTPUT GRID
!  INPUT    M            LENGTH OF TT

!  INPUT    IHOM         HOMOSZEDASTICY OF VARIANCE
!                        0: HOMOSZEDASTIC ERROR VARIABLES
!                        <> 0: IF THE VARIANCE SHOULD ESTIMATED AS
!                              SMOOTH FUNCTION.
!                        ****** DEFAULT VALUE: IHOM=0

!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!                        ****** DEFAULT VALUE: NUE=0

!  INPUT    KORD         ORDER OF KERNEL (<=6), FOR ISMO=0  ONLY
!                        NUE=0, KORD=2 OR KORD=4
!                        OR NUE=1, KORD=3 OR NUE=2, KORD=4 ARE ALLOWED
!                        ****** DEFAULT VALUE: KORD=NUE+2

!  INPUT    IRND         0: IF RANDOM GRID POINTS T MAY OCCUR
!                        <>0 ELSE  (ONLY NECESSARY IF S SHOULD BE
!                            COMPUTED)
!                        ****** DEFAULT VALUE IRND=0

!  INPUT    ISMO         0:ESTIMATING THE OPTIMAL LOCAL BANDWIDTH
!                        <>0 USING LOCAL INPUT BANDWIDTH-ARRAY IN BAN
!                        ****** DEFAULT VALUE ISMO=0

!  INPUT    M1           >=10, LENGTH OF W1, LARGE VALUES WILL INCREASE
!                        THE ACCURACY OF THE INTEGRAL APPROXIMATION
!                        ****** DEFAULT VALUE: M1=400

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! IN/OUTPUT TL/TU        LOWER/UPPER BOUND FOR INTEGRAL APPROXIMATION
!                        AND VARIANCE ESTIMATION (IF SIG=0 AND IHOM=0),
!                        IF TU<=TL, [TL,TU] ARE COMPUTED AS ABOUT
!                        THE 87% MIDDLE PART OF [T(1),T(N)]
!                        ****** DEFAULT VALUES: TL=1.0, TU=0.0

! IN/OUTPUT S(0:N)       IF S(N)<=S(0) THIS ARRAY IS COMPUTED AS
!                        MIDPOINTS OF T, FOR NON-RANDOM DESIGN AND AS
!                        SMOOTHED QUANTILES FOR RANDOM DESIGN
!                        ****** DEFAULT VALUES: S(0)=1.0, S(N)=0.0
!                               AND THE OTHER S(I) CAN BE UNDEFINED

! IN/OUTPUT SIG          RESIDUAL VARIANCE, ESTIMATED FOR SIG=0 OR
!                        IHOM<>0, ELSE GIVEN BY INPUT
!                        ****** DEFAULT VALUE: SIG=-1.0

! IN/OUTPUT BAN(M)       LOCAL PLUG-IN BANDWIDTH ARRAY
!                        ****** WILL BE SET IN SUBROUTINE FOR ISMO=0

! WORK     WN(0:N,5)     WORK ARRAY FOR KERNEL SMOOTHING ROUTINE
!                        ****** WILL BE SET IN SUBROUTINE
! WORK     W1(M1,3)      WORK ARRAY FOR INTEGRAL APPROXIMATION
!                        ****** WILL BE SET IN SUBROUTINE
! WORK     WM(M)         WORK ARRAY FOR INTEGRAL APPROXIMATION
!                        ****** WILL BE SET IN SUBROUTINE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! OUTPUT   Y(M)          KERNEL ESTIMATE WITH BOP(=B0 FOR ISMO<>0)
!                        ****** WILL BE SET IN SUBROUTINE
!-----------------------------------------------------------------------
!  USED SUBROUTINES: COFF, RESEST, KERNEL WITH FURTHER SUBROUTINES
!                     WHICH ARE CONTAINED IN THE FILE subs.f
!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: tt(:)
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: ihom
INTEGER, INTENT(IN)        :: nue
INTEGER, INTENT(IN OUT)    :: kord
INTEGER, INTENT(IN)        :: irnd
INTEGER, INTENT(IN)        :: ismo
INTEGER, INTENT(IN)        :: m1
REAL (dp), INTENT(IN OUT)  :: tl
REAL (dp), INTENT(IN OUT)  :: tu
REAL (dp), INTENT(IN OUT)  :: s(0:)
REAL (dp), INTENT(IN OUT)  :: sig
REAL (dp), INTENT(IN OUT)  :: ban(:)
REAL (dp), INTENT(OUT)     :: y(:)

! Local variables

REAL (dp)  :: wn(0:n,5), w1(m1,3), wm(m)
INTEGER    :: i, ii, iil, il, inputs, iprint, isort, it, itende, itt, iu,  &
              j, kk, kk2, kordv, nn, nuev, nyg, nyl
REAL (dp)  :: alpha, b, b2, bmax, bmin, bres, bs, bvar, const, dist,  &
              ex, exs, exsvi, fac, g1, g2, q, r2, rvar, s0, sn, snr, ssi, &
              tll, tuu, vi, wstep, xh, xi, xmy2, xxh
!-
!-------- 1. INITIALISATIONS AND SOME ERROR-CHECKS
REAL (dp), SAVE  :: bias(2,0:2) = RESHAPE(  &
                    (/ 0.2, 0.04762, 0.4286, 0.1515, 1.33, 0.6293 /), &
                    (/ 2,3 /) )
REAL (dp), SAVE  :: vark(2,0:2) = RESHAPE(  &
                    (/ 0.6, 1.250, 2.143, 11.93, 35.0, 381.6 /), (/ 2,3 /) )
REAL (dp), SAVE  :: fak2(2:4) = (/4., 36., 576. /)

nyg=0
inputs=0
!-------- IF NO ERRORS SHOULD BE WRITTEN ON STANDARD OUTPUT, SET IPRINT > 0
!-------- IF ERRORS AND VERY DETAILED WARNINGS SHOULD BE WRITTEN ON
!--------           STANDARD OUTPUT, SET IPRINT < 0
iprint=0

IF(nue > 4 .OR. nue < 0) THEN
  IF(iprint == 0) WRITE(*, *) 'lokern: Order of derivative not allowed'
  STOP
END IF
IF(nue > 2 .AND. ismo == 0) THEN
  IF(iprint == 0) WRITE(*, *) 'lokern: Order of derivative not allowed'
  STOP
END IF
IF(n <= 2) THEN
  IF(iprint == 0) WRITE(*, *) 'lokern: Number of data too small'
  STOP
END IF
IF(m < 1) THEN
  IF(iprint == 0) WRITE(*, *) 'lokern: No output points'
  STOP
END IF
IF(m1 < 10) THEN
  IF(iprint == 0) WRITE(*, *) 'lokern: Variable M1 is choosen too small'
  STOP
END IF
kk=(kord-nue)/2
IF(2*kk+nue /= kord) THEN
  IF(iprint == 0) WRITE(*, *) 'lokern: Kernel order not allowed, set to ',nue+2
  kord=nue+2
END IF
IF(kord > 4 .AND. ismo == 0) THEN
  IF(iprint == 0) WRITE(*, *) 'lokern: Kernel order not allowed, set to ',nue+2
  kord=nue+2
END IF
IF(kord > 6 .OR. kord <= nue) THEN
  IF(iprint == 0) WRITE(*, *) 'lokern: Kernel order not allowed, set to ',nue+2
  kord=nue+2
END IF
rvar=sig
!-
!-------- 2. COMPUTATION OF S-SEQUENCE
s0=1.5*t(1)-0.5*t(2)
sn=1.5*t(n)-0.5*t(n-1)
IF(s(n) <= s(0)) THEN
  inputs=1
  DO  i=1,n-1
    s(i)=.5*(t(i)+t(i+1))
  END DO
  s(0)=s0
  s(n)=sn
  IF(ismo /= 0 .AND. irnd /= 0) GO TO 230
ELSE
  IF(ismo /= 0) GO TO 230
END IF
!-
!-------- 3. COMPUTATION OF MINIMAL, MAXIMAL ALLOWED GLOBAL BANDWIDTH
bmax=(sn-s0)*.5
bmin=(sn-s0)/DBLE(n)*DBLE(kord-1)*.6
!-
!-------- 4. WARNINGS IF TT-GRID LARGER THAN T-GRID
IF(tt(1) < s0 .AND. tt(m) > sn .AND. iprint < 0) WRITE(*, *)   &
    'lokern: Extrapolation at both boundaries not optimized'
IF(tt(1) < s0 .AND. tt(m) <= sn .AND. iprint < 0) WRITE(*, *)   &
    'lokern: Extrapolation at left boundary not optimized'
IF(tt(1) >= s0 .AND. tt(m) > sn .AND. iprint < 0) WRITE(*, *)   &
    'lokern: Extrapolation at right boundary not optimized'
!-
!-------- 5. COMPUTE TL,TU AND THEIR T-GRID AS INNER PART FOR
!            INTEGRAL APPROXIMATION IN THE ITERATIONS
itt=0
51 IF (tu <= tl) THEN
  tl=.933*s0+.067*sn
  tu=.067*s0+.933*sn
  itt=itt+1
END IF
tl=MAX(s0,tl)
tu=MIN(sn,tu)
il=1
iu=n
wn(1,1)=0.0
wn(n,1)=0.0
DO  i=1,n
  IF(t(i) <= tl .OR. t(i) >= tu) wn(i,1)=0.0
  IF(t(i) > tl .AND. t(i) < tu) wn(i,1)=1.0
  IF(t(i) < tl) il=i+1
  IF(t(i) <= tu) iu=i
END DO
nn=iu-il+1
IF(nn == 0 .AND. itt == 0) THEN
  tu=tl-1.0
  GO TO 51
END IF
IF(nn == 0 .AND. itt == 1) THEN
  tu=sn
  tl=s0
  GO TO 51
END IF
!-
!-------- 6. COMPUTE T-GRID FOR INTEGRAL APPROXIMATION
DO  i=1,m1
  w1(i,2)=1.0
  w1(i,1)=tl+(tu-tl)*DBLE(i-1)/DBLE(m1-1)
END DO
!-
!-------- 7. CALCULATION OF WEIGHT FUNCTION
alpha=1.d0/DBLE(13)
DO  i=il,iu
  xi=(t(i) - tl)/alpha/(tu-tl)
  IF(xi > 1) GO TO 71
  wn(i,1)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
END DO
71 DO  i=iu,il,-1
  xi=(tu-t(i))/alpha/(tu-tl)
  IF(xi > 1) GO TO 73
  wn(i,1)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
END DO
73 DO  i=1,m1
  xi=(w1(i,1)-tl)/alpha/(tu-tl)
  IF(xi > 1) GO TO 75
  w1(i,2)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
END DO
75 DO  i=m1,1,-1
  xi=(tu-w1(i,1))/alpha/(tu-tl)
  IF(xi > 1) GO TO 77
  w1(i,2)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
END DO
!-
!-------- 8. COMPUTE CONSTANTS FOR ITERATION
77 ex=1./DBLE(kord+kord+1)
kk2=(kord-nue)
kk=kk2/2
!-
!-------- 9. ESTIMATING VARIANCE AND SMOOTHED PSEUDORESIDUALS
IF(sig <= .0 .AND. ihom == 0) CALL resest(t(il:),x(il:),nn,wn(il:,2),r2,sig)
IF(ihom /= 0) THEN
  CALL resest(t,x,n,wn(1:,2),snr,sig)
  bres=MAX(bmin,.2*nn**(-.2)*(s(iu)-s(il-1)))
  DO  i=1,n
    wn(i,3)=t(i)
    wn(i,2)=wn(i,2)*wn(i,2)
  END DO
  CALL kernel(t,wn(1:,2),n,bres,0,kk2,nyg,s, wn(1:,3),n,wn(1:,4))
ELSE
  CALL coff(wn(1:,4),n,sig)
END IF
!-
!-------- 10. ESTIMATE/COMPUTE INTEGRAL CONSTANT
100 vi=0.
DO  i=il,iu
  vi=vi + wn(i,1)*n*(s(i)-s(i-1))**2*wn(i,4)
END DO
!-
!-------- 11. REFINEMENT OF S-SEQUENCE FOR RANDOM DESIGN
IF(inputs == 1 .AND. irnd == 0) THEN
  DO  i=0,n
    wn(i,5)=DBLE(i)/DBLE(n+1)
    wn(i,2)=(DBLE(i)+.5)/DBLE(n+1)
    wn(i,3)=wn(i,2)
  END DO
  exs=-DBLE(3*kord+1)/DBLE(6*kord+3)
  exsvi=DBLE(kord)/DBLE(6*kord+3)
  bs=0.1*(vi/(sn-s0)**2)**exsvi*n**exs
  CALL kernel(wn(1:,5),t,n,bs,0,2,nyg,wn(0:,3),wn(0:,2),n+1,s(0:))
  111 isort=0
  vi=0.0
  DO  i=1,n
    vi=vi + wn(i,1)*n*(s(i)-s(i-1))**2*wn(i,4)
    IF(s(i) < s(i-1)) THEN
      ssi=s(i-1)
      s(i-1)=s(i)
      s(i)=ssi
      isort=1
    END IF
  END DO
  IF(isort == 1) GO TO 111
  IF(ismo /= 0) GO TO 230
END IF
b=bmin*2.
!-
!-------- 12. COMPUTE INFLATION CONSTANT AND EXPONENT AND LOOP OF ITERATIONS
const=DBLE(2*nue+1)*fak2(kord)*vark(kk,nue)*vi  &
      /(DBLE(2*kord-2*nue)*bias(kk,nue)**2*DBLE(n))
fac=1.1*(1. + (nue/10.) + 0.05*(kord-nue-2.))*n**(2./DBLE((2*kord+1)*(2*kord+3)))
itende=1 + 2*kord + kord*(2*kord+1)

DO  it=1,itende
!-
!-------- 13. ESTIMATE DERIVATIVE OF ORDER KORD IN ITERATIONS
  b2=b*fac
  b2=MAX(b2,bmin/(kord-1) * DBLE(kord+1))
  b2=MIN(b2,bmax)
  CALL kernel(t,x,n,b2,kord,kord+2,nyg,s,w1(1:,1),m1,w1(1:,3))
!-
!-------- 14. ESTIMATE INTEGRALFUNCTIONAL IN ITERATIONS
  xmy2=.75*(w1(1,2)*w1(1,3)*w1(1,3) + w1(m1,2)*w1(m1,3)*w1(m1,3))
  DO  i=2,m1-1
    xmy2=xmy2 + w1(i,2)*w1(i,3)*w1(i,3)
  END DO
  xmy2=xmy2*(tu-tl)/m1
!-
!-------- 15. FINISH OF ITERATIONS
  b=(const/xmy2)**ex
  b=MAX(bmin,b)
  b=MIN(bmax,b)
  
END DO
!-------- 16  COMPUTE SMOOTHED FUNCTION WITH GLOBAL PLUG-IN BANDWIDTH
CALL kernel(t, x, n, b, nue, kord, nyg, s, tt, m, y)
!-
!-------- 17. VARIANCE CHECK
IF(ihom /= 0) sig=rvar
IF(rvar == sig .OR. r2 < .88 .OR. ihom /= 0 .OR. nue > 0) GO TO 180
ii=0
iil=0
j=2
tll=MAX(tl, tt(1))
tuu=MIN(tu, tt(m))
DO  i=il,iu
  IF(t(i) < tll .OR. t(i) > tuu) CYCLE
  ii=ii+1
  IF(iil == 0) iil=i
  171 IF(tt(j) < t(i)) THEN
    j=j+1
    IF(j <= m) GO TO 171
  END IF
  wn(ii,3)=x(i)-y(j) + (y(j)-y(j-1))*(tt(j)-t(i))/(tt(j)-tt(j-1))
END DO
CALL resest(t(iil:), wn(1:,3), ii, wn(1:,4), snr, rvar)
q=sig/rvar
CALL coff(wn(1:,4), n, sig)
IF(q <= 2.) GO TO 180
IF(q > 5. .AND. r2 > .95) rvar=rvar*.5
sig=rvar
CALL coff(wn(1:,4), n, sig)
GO TO 100
!-
!-------- 18. LOCAL INITIALIZATIONS
180 bvar=b
nuev=0
kordv=2
nyl=1
!-
!-------- 19. COMPUTE INNER BANDWIDTHS
g1=0.86*(1.0 + DBLE(kord-nue-2)*.05)*b
g1=g1*DBLE(n)**(4./DBLE(2*kord+1)/(2*kord+5))
g1=MAX(g1, bmin/DBLE(kord-1)*DBLE(kord+1))
g1=MIN(g1, bmax)

g2=1.4*(1.0 + DBLE(kord-nue-2)*0.05)*b
g2=g2*DBLE(n)**(2./DBLE(2*kord+1)/DBLE(2*kord+3))
g2=MAX(g2, bmin)
g2=MIN(g2, bmax)
!-
!-------- 20. ESTIMATE/COMPUTE INTEGRAL CONSTANT VI LOCALLY
DO  i=1,n
  wn(i,4)=DBLE(n)*wn(i,4)*(s(i)-s(i-1))
END DO
DO  j=1,m
  ban(j)=bvar
  wm(j)=tt(j)
  IF(tt(j) < s(0)+g1) THEN
    dist=((tt(j)-g1-s(0))/g1)**2
    ban(j)=bvar*(1.0 + dist)
    ban(j)=MIN(ban(j), bmax)
    wm(j)=tt(j) + 0.5*dist*g1
  ELSE IF(tt(j) > s(n)-g1) THEN
    dist=((tt(j)-s(n)+g1)/g1)**2
    ban(j)=bvar*(1.0 + dist)
    ban(j)=MIN(ban(j), bmax)
    wm(j)=tt(j) - .5*dist*g1
  END IF
END DO
CALL kernel(t, wn(1:,4), n, bvar, nuev, kordv, nyl, s, wm, m, ban)
!-
!-------- 21. ESTIMATION OF KORDTH DERIVATIVE LOCALLY
wstep=(tt(m)-tt(1))/(m1-2)
DO  j=2,m1
  w1(j,2)=tt(1) + (j-2)*wstep
  w1(j,1)=tt(1) + (j-1.5)*wstep
END DO
w1(1,1)=tt(1) + 0.5*wstep

CALL kernel(t, x, n, g1, kord, kord+2, nyg, s, w1(2:,2), m1-1, w1(2:,3))

DO  j=2,m1
  w1(j,3)=w1(j,3)*w1(j,3)
END DO
DO  j=1,m
  y(j)=g2
  IF(tt(j) < s(0)+g1) THEN
    y(j)=g2*(1.0 + ((tt(j)-g1-s(0))/g1)**2)
    y(j)=MIN(y(j), bmax)
  ELSE IF(tt(j) > s(n)-g1) THEN
    y(j)=g2*(1.0 + ((tt(j)-s(n)+g1)/g1)**2)
    y(j)=MIN(y(j), bmax)
  END IF
END DO

CALL kernp(w1(2:,2), w1(2:,3), m1-1, g2, nuev, kordv, nyl, w1(1:,1), wm, m, y)
!-
!-------- 22. FINISH
DO  j=1,m
  xh=bmin**(2*kord+1)*ABS(y(j))*vi/const
  xxh=const*ABS(ban(j))/vi/bmax**(2*kord+1)
  IF(ban(j) < xh) THEN
    ban(j)=bmin
  ELSE
    IF(y(j) < xxh) THEN
      ban(j)=bmax
    ELSE
      ban(j)=(const*ban(j)/y(j)/vi)**ex
    END IF
  END IF
END DO
!-
!-------- 23. COMPUTE SMOOTHED FUNCTION WITH LOCAL PLUG-IN BANDWIDTH
230 y(1:m)=ban(1:m)
CALL kernel(t, x, n, b, nue, kord, nyl, s, tt, m, y)

RETURN
END SUBROUTINE lokern



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
!C     SEVERAL KERNEL SMOOTHING SUBROUTINES WHICH ARE USED BY GLKERN.F
!C     AND LOKERN.F, VERSION JANUARY 1997
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     THIS FILE CONTAINS:
!C
!C     SUBROUTINE RESEST(T,X,N,RES,SNR,SIGMA2)
!C                FOR VARIANCE ESTIMATION
!C
!C     SUBROUTINE KERNEL(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
!C                DRIVER SUBROUTINE FOR KERNEL REGRESSION ESTIMATION
!C                CALLS FAST OR CONVENTIAL KERNEL ROUTINE
!C
!C     SUBROUTINE KERNP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
!C                DRIVER SUBROUTINE FOR KERNEL REGRESSION ESTIMATION
!C                WITHOUT USE OF BOUNDARY KERNELS
!C---------------------------------------------------------------------
!C     SUBROUTINE KERNFA(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
!C                FAST ALGORITHM FOR KERNEL ESTIMATION
!C     SUBROUTINE DREG(SW,A1,A2,IORD,X,SL,SR,T,B,IFLOP)
!C                USED BY SUBROUTINE KERNFA,KERNFP
!C     SUBROUTINE LREG(SW,A3,IORD,D,DOLD,Q,C)
!C                USED BY SUBROUTINE KERNFA,KERNFP
!C     SUBROUTINE FREG(SW,NUE,KORD,IBOUN,Y,C,ICALL,A)
!C                USED BY SUBROUTINE KERNFA,KERNFP
!C     SUBROUTINE KERNFP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
!C                FAST ALGORITHM FOR KERNP ESTIMATION WITHOUT BOUNDARY
!C---------------------------------------------------------------------
!C     SUBROUTINE KERNCL(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
!C                CONVENTIONAL ALGORITHM FOR KERNEL ESTIMATION
!C     SUBROUTINE SMO(S,X,N,TAU,WID,NUE,IORD,IBOUN,IST,S1,C,Y)
!C                SINGLE ESTIMATION STEP, USED BY KERNCL
!C     SUBROUTINE KERNCP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
!C                CONVENTIONAL ALGORITHM WITHOUT BOUNDARY KERNELS
!C     SUBROUTINE SMOP(S,X,N,TAU,WID,NUE,IORD,IBOUN,IST,S1,C,Y)
!C                SINGLE ESTIMATION STEP, USED BY KERNCP
!C---------------------------------------------------------------------
!C     SUBROUTINE COFFI(NUE,KORD,C)
!C                KERNEL COEFFICIENT OF POLYNOMIAL KERNELS USED BY
!C                KERNCL,KERNCP AND KERNFP
!C     SUBROUTINE COFFB(NUE,KORD,Q,IBOUN,C)
!C                KERNEL COEFFICIENT OF POLYNOMIAL BOUNDARY KERNELS
!C                USED BY KERNFA AND KERNCL
!C---------------------------------------------------------------------
!C     SUBROUTINE COFF(X,N,FA)
!C                SIMPLE SUBROUTINE FOR ARRAY INITIALIZATION
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


SUBROUTINE resest(t, x, n, res, snr, sigma2)
!--------------------------------------------------------------------
!    VERSION: JUNE, 1996

!    PURPOSE:

!    COMPUTES ONE-LEAVE-OUT RESIDUALS FOR NONPARAMETRIC ESTIMATION
!    OF RESIDUAL VARIANCE (LOCAL LINEAR APPROXIMATION FOLLOWED BY
!    REWEIGHTING)

!  PARAMETERS:

!  INPUT   T(N)      ABSCISSAE (ORDERED: T(I)<=T(I+1))
!  INPUT   X(N)      DATA
!  INPUT   N         LENGTH OF DATA ( >2 )
!  OUTPUT  RES(N)    RESIDUALS AT T(1),...,T(N)
!  OUTPUT  SNR       EXPLAINED VARIANCE OF THE TRUE CURVE
!  OUTPUT  SIGMA2    ESTIMATION OF SIGMA**2 (RESIDUAL VARIANCE)

!--------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: t(:)
REAL (dp), INTENT(IN)   :: x(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: res(:)
REAL (dp), INTENT(OUT)  :: snr
REAL (dp), INTENT(OUT)  :: sigma2

! Local variables

INTEGER    :: i
REAL (dp)  :: dn, g1, g2, ex, ex2, sx, tt

!-
sigma2=0.
ex=x(1)*(t(2)-t(1))
ex2=x(1)*ex
DO  i=2,n-1
  tt=t(i+1)-t(i-1)
  IF(tt /= 0.) g1=(t(i+1)-t(i))/tt
  IF(tt == 0.) g1=.5
  g2=1.-g1
  res(i)=(x(i)-g1*x(i-1)-g2*x(i+1))/SQRT(1.+g1*g1+g2*g2)
  sigma2=sigma2+res(i)*res(i)
  sx=x(i)*tt
  ex=ex+sx
  ex2=ex2+x(i)*sx
END DO
tt=t(3)-t(2)
IF(tt /= 0.) g1=(t(1)-t(2))/tt
IF(tt == 0.) g1=.5
g2=1.-g1
res(1)=(x(1)-g1*x(3)-g2*x(2))/SQRT(1.+g1*g1+g2*g2)
tt=t(n-1)-t(n-2)
IF(tt /= 0.) g1=(t(n-1)-t(n))/tt
IF(tt == 0.) g1=.5
g2=1.-g1
res(n)=(x(n)-g1*x(n-2)-g2*x(n-1))/SQRT(1.+g1*g1+g2*g2)
sigma2=(sigma2+res(1)*res(1)+res(n)*res(n))/n
!-
sx=x(n)*(t(n)-t(n-1))
dn=2.*(t(n)-t(1))
ex=(ex+sx)/dn
ex2=(ex2+x(n)*sx)/dn
IF(ex2 == 0) snr=0.
IF(ex2 > 0) snr=1-sigma2/(ex2-ex*ex)
RETURN
END SUBROUTINE resest



SUBROUTINE kernel(t, x, n, b, nue, kord, ny, s, tt, m, y)
!-----------------------------------------------------------------------
!       SHORT-VERSION MAY, 1995

!       DRIVER SUBROUTINE FOR KERNEL SMOOTHING, CHOOSES BETWEEN
!       STANDARD AND O(N) ALGORITHM

!  PARAMETERS :

!  INPUT    T(N)         INPUT GRID (REGRESSION DESIGN)
!  INPUT    X(N)         DATA, GIVEN ON T(N)
!  INPUT    N            LENGTH OF X
!  INPUT    B            ONE SIDED BANDWIDTH (FOR NY=1 MEAN BANDWIDTH)
!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!  INPUT    KORD         ORDER OF KERNEL (<=6); DEFAULT IS KORD=NUE+2
!  INPUT    NY           0: GLOBAL BANDWIDTH (DEFAULT)
!                        1: VARIABLE BANDWIDTHS, GIVEN IN Y AS INPUT
!  INPUT    S(0:N)       INTERPOLATION SEQUENCE
!  INPUT    TT(M)        OUTPUT GRID.  MUST BE PART OF INPUT GRID FOR IEQ=0
!  INPUT    M            NUMBER OF POINTS WHERE FUNCTION IS ESTIMATED,
!                         OR  LENGTH OF TT. DEFAULT IS M=400
!  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
!  OUTPUT   Y(M)         ESTIMATED REGRESSION FUNCTION

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: b
INTEGER, INTENT(IN)        :: nue
INTEGER, INTENT(IN)        :: kord
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN)      :: s(0:)
REAL (dp), INTENT(IN)      :: tt(:)
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN OUT)  :: y(:)

! Local variable

REAL (dp) :: chan

!-
!------  COMPUTING CHANGE POINT
chan=(5. + kord)*MAX(1.0, SQRT(REAL(n)/REAL(m)))
!------
IF(b*(n-1)/(t(n)-t(1)) < chan) THEN
  CALL kerncl(t, x, n, b, nue, kord, ny, s, tt, m, y)
ELSE
  CALL kernfa(t, x, n, b, nue, kord, ny, s, tt, m, y)
END IF
!-
RETURN
END SUBROUTINE kernel



SUBROUTINE kernp(t, x, n, b, nue, kord, ny, s, tt, m, y)
!-----------------------------------------------------------------------
!       SHORT-VERSION JANUARY, 1997

!       DRIVER SUBROUTINE FOR KERNEL SMOOTHING, CHOOSES BETWEEN
!       STANDARD AND O(N) ALGORITHM WITHOUT USING BOUNDARY KERNELS

!  PARAMETERS :

!  INPUT    T(N)         INPUT GRID (REGRESSION DESIGN)
!  INPUT    X(N)         DATA, GIVEN ON T(N)
!  INPUT    N            LENGTH OF X
!  INPUT    B            ONE SIDED BANDWIDTH (FOR NY=1 MEAN BANDWIDTH)
!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!  INPUT    KORD         ORDER OF KERNEL (<=6); DEFAULT IS KORD=NUE+2
!  INPUT    NY           0: GLOBAL BANDWIDTH (DEFAULT)
!                        1: VARIABLE BANDWIDTHS, GIVEN IN Y AS INPUT
!  INPUT    S(0:N)       INTERPOLATION SEQUENCE
!  INPUT    TT(M)        OUTPUT GRID. MUST BE PART OF INPUT GRID FOR IEQ=0
!  INPUT    M            NUMBER OF POINTS WHERE FUNCTION IS ESTIMATED,
!                         OR  LENGTH OF TT.  DEFAULT IS M=400
!  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
!  OUTPUT   Y(M)         ESTIMATED REGRESSION FUNCTION

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: b
INTEGER, INTENT(IN)        :: nue
INTEGER, INTENT(IN)        :: kord
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN)      :: s(0:)
REAL (dp), INTENT(IN)      :: tt(:)
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN OUT)  :: y(:)

! Local variable

REAL (dp) :: chan

!-
!------  COMPUTING CHANGE POINT
chan=(5.+kord)*MAX(1., SQRT(REAL(n)/REAL(m)))
!------
IF(b*(n-1)/(t(n)-t(1)) < chan) THEN
  CALL kerncp(t, x, n, b, nue, kord, ny, s, tt, m, y)
ELSE
  CALL kernfp(t, x, n, b, nue, kord, ny, s, tt, m, y)
END IF
!-
RETURN
END SUBROUTINE kernp



SUBROUTINE kernfa(t, x, n, b, nue, kord, ny, s, tt, m, y)
!-----------------------------------------------------------------------
!       SHORT-VERSION: MAY, 1995

!       PURPOSE:

!       COMPUTATION OF KERNEL ESTIMATE USING O(N) ALGORITHM BASED ON
!       LEGENDRE POLYNOMIALS, GENERAL SPACED DESIGN AND LOCAL
!       BANDWIDTH ALLOWED. (NEW INITIALISATIONS OF THE LEGENDRE SUMS
!       FOR NUMERICAL REASONS)

!  PARAMETERS :

!  INPUT    T(N)         INPUT GRID
!  INPUT    X(N)         DATA, GIVEN ON T(N)
!  INPUT    N            LENGTH OF X
!  INPUT    B            ONE SIDED BANDWIDTH
!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!  INPUT    KORD         ORDER OF KERNEL (<=6)
!  INPUT    NY           0, GLOBAL BANDWIDTH; 1, LOCAL BANDWIDTH IN Y
!  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
!  INPUT    TT(M)        OUTPUT GRID
!  INPUT    M            NUMBER OF POINTS TO ESTIMATE
!  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
!  OUTPUT   Y(M)         ESTIMATED FUNCTION

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: b
INTEGER, INTENT(IN)        :: nue
INTEGER, INTENT(IN)        :: kord
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN)      :: s(0:n)
REAL (dp), INTENT(IN)      :: tt(:)
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN OUT)  :: y(:)


INTEGER   :: j, k, iord, init, icall, i, iboun
INTEGER   :: jl, jr, jnr, jnl

REAL (dp) :: c(7), sw(7), xf(7), dold
REAL (dp) :: a(7,7), a1(7), a2(7), a3(7,7), cm(7,6)
REAL (dp) :: bmin, bmax, bb, s0, sn, wwl, wwr, wid, wr, wido
!-
!------ COMPUTE CONSTANTS FOR LATER USE
s0=1.5*t(1)-0.5*t(2)
sn=1.5*t(n)-0.5*t(n-1)
bmin=(sn-s0)*.6D0/DBLE(n)*DBLE(kord-1)
bmax=(s(n)-s(0))*.5
IF(kord == 2) bmin=bmin*0.1D0
iord=kord+1
DO  k=3,iord
  a1(k)=DBLE(2*k-1)/DBLE(k)
  a2(k)=DBLE(1-k)/DBLE(k)
END DO
!-
init=0
icall=0
dold=0.d0
!-
!------ SMOOTHING LOOP
DO  i=1,m
  bb=b
  IF (ny == 1) bb=y(i)
  IF(bb < bmin) bb=bmin
  IF(bb > bmax) bb=bmax
  iboun=0
!-
!------ COMPUTE LEFT BOUNDARY KERNEL
  IF(tt(i) < s(0)+bb) THEN
    wwl=s(0)
    wwr=s(0)+bb+bb
    wid=wwr-tt(i)
    iboun=1
    CALL coffb(nue,kord,(tt(i)-s(0))/wid,iboun,c)
  END IF
!-
!------ COMPUTE RIGHT BOUNDARY KERNEL
  IF(tt(i)+bb > s(n)) THEN
    wwl=s(n)-(bb+bb)
    wwr=s(n)
    wid=tt(i)-wwl
    iboun=-1
    CALL coffb(nue,kord,(s(n)-tt(i))/wid,iboun,c)
  END IF
!-
!------ NO BOUNDARY
  IF(iboun == 0) THEN
    wid=bb
    wwl=tt(i)-bb
    wwr=tt(i)+bb
  END IF
!-
!------ INITIALISATION FOR INIT=0
  IF(init == 0) THEN
    sw(1:iord)=0.
    jl=1
    DO  j=1,n
      IF(s(j-1) < wwl) THEN
        jl=j+1
      ELSE
        IF(s(j) > wwr) EXIT
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
      END IF
    END DO
    jr=j-1
    wr=wwr
    init=1
    GO TO 6666
  ELSE
    init=init+1
  END IF
!-
!------ COMPARE OLD SUM WITH NEW SMOOTHING INTERVALL TT(I)-B,TT(I)+B
  IF(s(jr-1) >= wwl) THEN
    jnr=jr
    jnl=jl
    IF(s(jr) > wwr) THEN
      DO  j=jr,jl,-1
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i-1),wido,-1)
        jnr=j-1
        IF(s(jnr) <= wwr) EXIT
      END DO
    END IF
    IF(s(jl-1) < wwl) THEN
      DO  j=jl,jr
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i-1),wido,-1)
        jnl=j+1
        IF(s(j) >= wwl) EXIT
      END DO
    END IF
!-
!------ UPDATING OF SW
    CALL lreg(sw,a3,iord,(tt(i)-tt(i-1))/wid,dold,wido/wid,cm)
    IF(jnr == jr) THEN
      DO  j=jr+1,n
        IF(s(j) > wwr) EXIT
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
        jnr=j
      END DO
    END IF
    jr=jnr
    IF(jl == jnl) THEN
      DO   j=jl-1,1,-1
        IF(s(j-1) < wwl) EXIT
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
        jnl=j
      END DO
    END IF
    jl=jnl
  ELSE
!-
!------ NEW INITIALISATION OF SW
    sw(1:iord)=0.
    DO  j=jr,n
      IF(s(j-1) < wwl) THEN
        jl=j+1
      ELSE
        IF(s(j) > wwr) EXIT
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
      END IF
    END DO
    jr=j-1
    wr=wwr
  END IF
!-
!------ IF BANDWIDTH IS TOO SMALL NO SMOOTHING
  6666 IF(wwl >= s(jr-1) .AND. wwr <= s(jr)) THEN
    y(i)=x(jr)
    IF(nue > 0) y(i)=0.d0
  ELSE
!-
!------ ADD FIRST AND LAST POINT OF THE SMOOTHING INTERVAL
    DO  k=1,iord
      xf(k)=sw(k)
    END DO
    IF(jl /= 1) CALL dreg(xf,a1,a2,iord,x(jl-1),wwl,s(jl-1),tt(i),wid,1)
    IF(jr /= n) CALL dreg(xf,a1,a2,iord,x(jr+1),s(jr),wwr,tt(i),wid,1)
!-
!------ NOW THE SUMS ARE BUILT THAT ARE NEEDED TO COMPUTE THE ESTIMATE
    CALL freg(xf,nue,kord,iboun,y(i),c,icall,a)
    IF(nue > 0) y(i)=y(i)/(wid**nue)
  END IF
!-
!------ NEW INITIALISATION ?
  IF(jl > jr .OR. wwl > wr .OR. init > 100) init=0
  wido=wid
!-
END DO
!-
RETURN
END SUBROUTINE kernfa



SUBROUTINE kernfp(t, x, n, b, nue, kord, ny, s, tt, m, y)
!-----------------------------------------------------------------------
!       SHORT-VERSION: JANUARY, 1997

!       PURPOSE:

!       COMPUTATION OF KERNEL ESTIMATE USING O(N) ALGORITHM BASED ON
!       LEGENDRE POLYNOMIALS, GENERAL SPACED DESIGN AND LOCAL
!       BANDWIDTH ALLOWED. (NEW INITIALISATIONS OF THE LEGENDRE SUMS
!       FOR NUMERICAL REASONS) WITHOUT BOUNDARY KERNELS, JUST NORMALIZING

!  PARAMETERS :

!  INPUT    T(N)         INPUT GRID
!  INPUT    X(N)         DATA, GIVEN ON T(N)
!  INPUT    N            LENGTH OF X
!  INPUT    B            ONE SIDED BANDWIDTH
!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!  INPUT    KORD         ORDER OF KERNEL (<=6)
!  INPUT    NY           0, GLOBAL BANDWIDTH; 1, LOCAL BANDWIDTH IN Y
!  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
!  INPUT    TT(M)        OUTPUT GRID
!  INPUT    M            NUMBER OF POINTS TO ESTIMATE
!  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
!  OUTPUT   Y(M)         ESTIMATED FUNCTION

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: b
INTEGER, INTENT(IN)        :: nue
INTEGER, INTENT(IN)        :: kord
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN)      :: s(0:n)
REAL (dp), INTENT(IN)      :: tt(:)
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN OUT)  :: y(:)


INTEGER    :: j, k, iord, init, icall, i, iboun
INTEGER    :: jl, jr, jnr, jnl

REAL (dp)  :: c(7), sw(7), xf(7), dold, qq, q, xnor
REAL (dp)  :: a(7,7), a1(7), a2(7), a3(7,7), cm(7,6)
REAL (dp)  :: bmin, bmax, bb, s0, sn, wwl, wwr, wid, wr, wido
!-
!------ COMPUTE CONSTANTS FOR LATER USE
s0=1.5*t(1)-0.5*t(2)
sn=1.5*t(n)-0.5*t(n-1)
bmin=(sn-s0)*.6D0/DBLE(n)*DBLE(kord-1)
bmax=(s(n)-s(0))*.5
IF(kord == 2) bmin=bmin*0.1D0
iord=kord+1
CALL coffi(nue,kord,c)
DO  k=3,iord
  a1(k)=DBLE(2*k-1)/DBLE(k)
  a2(k)=DBLE(1-k)/DBLE(k)
END DO
!-
init=0
icall=0
dold=0.d0
!-
!------ SMOOTHING LOOP
DO  i=1,m
  bb=b
  IF (ny == 1) bb=y(i)
  IF(bb < bmin) bb=bmin
  IF(bb > bmax) bb=bmax
  iboun=0
!-
!------ COMPUTE LEFT BOUNDARY
  IF(tt(i) < s(0)+bb) THEN
    wwl=s(0)
    wwr=s(0)+bb+bb
    wid=wwr-tt(i)
    iboun=1
  END IF
!-
!------ COMPUTE RIGHT BOUNDARY
  IF(tt(i)+bb > s(n)) THEN
    wwl=s(n)-(bb+bb)
    wwr=s(n)
    wid=tt(i)-wwl
    iboun=-1
  END IF
!-
!------ NO BOUNDARY
  IF(iboun == 0) THEN
    wid=bb
    wwl=tt(i)-bb
    wwr=tt(i)+bb
    xnor=1.d0
  END IF
!-
!------ COMPUTE NORMALIZING CONSTANT
  IF(iboun /= 0) THEN
    IF(iboun == 1) q=(tt(i)-s(0))/wid
    IF(iboun == -1) q=(s(n)-tt(i))/wid
    qq=q*q
    xnor=c(1)*(1.d0+q)
    DO  k=3,iord,2
      q=q*qq
      xnor=xnor+c(k)*(1.d0+q)
    END DO
    iboun=0
  END IF
!-
!------ INITIALISATION FOR INIT=0
  IF(init == 0) THEN
    DO  k=1,iord
      sw(k)=0.
    END DO
    jl=1
    DO  j=1,n
      IF(s(j-1) < wwl) THEN
        jl=j+1
      ELSE
        IF(s(j) > wwr) EXIT
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
      END IF
    END DO
    jr=j-1
    wr=wwr
    init=1
    GO TO 6666
  ELSE
    init=init+1
  END IF
!-
!------ COMPARE OLD SUM WITH NEW SMOOTHING INTERVALL TT(I)-B,TT(I)+B
  IF(s(jr-1) >= wwl) THEN
    jnr=jr
    jnl=jl
    IF(s(jr) > wwr) THEN
      DO  j=jr,jl,-1
        CALL dreg(sw, a1, a2, iord, x(j), s(j-1), s(j), tt(i-1), wido, -1)
        jnr=j-1
        IF(s(jnr) <= wwr) EXIT
      END DO
    END IF
    IF(s(jl-1) < wwl) THEN
      DO  j=jl,jr
        CALL dreg(sw, a1, a2, iord, x(j), s(j-1), s(j), tt(i-1), wido, -1)
        jnl=j+1
        IF(s(j) >= wwl) EXIT
      END DO
    END IF
!-
!------ UPDATING OF SW
    CALL lreg(sw, a3, iord, (tt(i)-tt(i-1))/wid, dold, wido/wid, cm)
    IF(jnr == jr) THEN
      DO  j=jr+1,n
        IF(s(j) > wwr) EXIT
        CALL dreg(sw, a1, a2, iord, x(j), s(j-1), s(j), tt(i), wid, 1)
        jnr=j
      END DO
    END IF
    jr=jnr
    IF(jl == jnl) THEN
      DO   j=jl-1,1,-1
        IF(s(j-1) < wwl) EXIT
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
        jnl=j
      END DO
    END IF
    jl=jnl
  ELSE
!-
!------ NEW INITIALISATION OF SW
    sw(1:iord)=0.
    DO  j=jr,n
      IF(s(j-1) < wwl) THEN
        jl=j+1
      ELSE
        IF(s(j) > wwr) EXIT
        CALL dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
      END IF
    END DO
    jr=j-1
    wr=wwr
  END IF
!-
!------ IF BANDWIDTH IS TOO SMALL NO SMOOTHING
  6666 IF(wwl >= s(jr-1) .AND. wwr <= s(jr)) THEN
    y(i)=x(jr)
    IF(nue > 0) y(i)=0.d0
  ELSE
!-
!------ ADD FIRST AND LAST POINT OF THE SMOOTHING INTERVAL
    DO  k=1,iord
      xf(k)=sw(k)
    END DO
    IF(jl /= 1) CALL dreg(xf,a1,a2,iord,x(jl-1),wwl,s(jl-1),tt(i),wid,1)
    IF(jr /= n) CALL dreg(xf,a1,a2,iord,x(jr+1),s(jr),wwr,tt(i),wid,1)
!-
!------ NOW THE SUMS ARE BUILT THAT ARE NEEDED TO COMPUTE THE ESTIMATE
    CALL freg(xf,nue,kord,iboun,y(i),c,icall,a)
    IF(nue > 0) y(i)=y(i)/(wid**nue)
    IF(xnor /= 1.d0) y(i)=y(i)/xnor
  END IF
!-
!------ NEW INITIALISATION ?
  IF(jl > jr .OR. wwl > wr .OR. init > 100) init=0
  wido=wid
!-
END DO
!-
RETURN
END SUBROUTINE kernfp



SUBROUTINE dreg(sw, a1, a2, iord, x, sl, sr, t, b, iflop)
!-----------------------------------------------------------------------
!    VERSION: MAY, 1995

!    PURPOSE:

!    COMPUTES NEW LEGENDRE SUMS (FOR REGRESSION)

!    PARAMETERS:
!                 **************   INPUT   *******************

!     SW(IORD)  :   OLD SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.
!     A1(7)     :   CONSTANTS OF RECURSIVE FORMULA FOR LEGENDRE POL.
!     A2(7)     :                             "
!     IORD      :   ORDER OF KERNEL POLYNOMIAL
!     X         :   DATA POINT
!     SL        :   LEFT S-VALUE
!     SR        :   RIGHT S-VALUE
!     T         :   POINT WHERE THE SMOOTHED VALUE IS TO BE ESTIMATED
!     B         :   BANDWIDTH
!     IFLOP     :   1: ADDITION, ELSE SUBTRACTION


!                 **************   OUTPUT   *******************

!     SW(IORD)  :   NEW SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.

!-----------------------------------------------------------------------

REAL (dp), INTENT(OUT)  :: sw(7)
REAL (dp), INTENT(IN)   :: a1(7)
REAL (dp), INTENT(IN)   :: a2(7)
INTEGER, INTENT(IN)     :: iord
REAL (dp), INTENT(IN)   :: x
REAL (dp), INTENT(IN)   :: sl
REAL (dp), INTENT(IN)   :: sr
REAL (dp), INTENT(IN)   :: t
REAL (dp), INTENT(IN)   :: b
INTEGER, INTENT(IN)     :: iflop


INTEGER    :: k
REAL (dp)  :: p(7,2)
!-
!------  COMPUTE LEGENDRE POLYNOMIALS
p(1,1)=(t-sl)/b
p(1,2)=(t-sr)/b
p(2,1)=1.5D0*p(1,1)*p(1,1)-.5D0
p(2,2)=1.5D0*p(1,2)*p(1,2)-.5D0
DO  k=3,iord
  p(k,1)=a1(k)*p(k-1,1)*p(1,1)+a2(k)*p(k-2,1)
  p(k,2)=a1(k)*p(k-1,2)*p(1,2)+a2(k)*p(k-2,2)
END DO
!-
!------  COMPUTE NEW LEGENDRE SUMS
IF(iflop == 1) THEN
  DO  k=1,iord
    sw(k)=sw(k)+(p(k,1)-p(k,2))*x
  END DO
ELSE
  DO  k=1,iord
    sw(k)=sw(k)+(p(k,2)-p(k,1))*x
  END DO
END IF
RETURN
END SUBROUTINE dreg



SUBROUTINE lreg(sw, a3, iord, d, dold, q, c)
!------------------------------------------------------------------
!    VERSION: MAY, 1995

!    PURPOSE:

!    UPDATE OF SW-SEQUENCE ACCORDING TO NEW BANDWIDTH AND NEW DATA
!    (VERSION FOR REGRESSION)

!    PARAMETERS :
!                       **************   INPUT   *******************

!            SW(IORD)  :  SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.
!            IORD      :  ORDER OF KERNEL POLYNOMIAL
!            D         :  DIST. TO THE NEXT POINT DIVIDED BY BANDW.
!            DOLD      :  D       PREVIOUS STEP
!            Q         :  NEW BANDWIDTH DIVIDED BY OLD BANDWIDTH
!                       **************    WORK   *******************

!            A3(7,7)   :  MATRIX (P*Q*P)**(-1)
!            C(7,6)    :  MATRIX OF COEFFICIENTS
!                       **************   OUTPUT   ******************

!            SW(IORD)  :  UPDATED VERSION OF SW

!---------------------------------------------------------------------

REAL (dp), INTENT(IN OUT)  :: sw(7)
REAL (dp), INTENT(OUT)     :: a3(7,7)
INTEGER, INTENT(IN)        :: iord
REAL (dp), INTENT(IN)      :: d
REAL (dp), INTENT(IN OUT)  :: dold
REAL (dp), INTENT(IN)      :: q
REAL (dp), INTENT(OUT)     :: c(7,6)


INTEGER    :: k, i, l
REAL (dp)  :: dd, ww, qq, xx

!-
!- BUILD UP MATRIX
IF(dold /= d .OR. dold == 0) THEN
  dold=d
  dd=d*d
!-
  IF(iord == 7) THEN
    c(7,6)=13.d0*d
    c(7,5)=71.5D0*dd
    c(7,4)=(214.5D0*dd+9.d0)*d
    c(7,3)=(375.375D0*dd+77.d0)*dd
    c(7,2)=((375.375D0*dd+247.5D0)*dd+5.d0)*d
    c(7,1)=((187.6875D0*dd+346.5D0)*dd+40.5D0)*dd
  END IF
!-
  IF(iord >= 6) THEN
    c(6,5)=11.d0*d
    c(6,4)=49.5D0*dd
    c(6,3)=(115.5D0*dd+7.d0)*d
    c(6,2)=(144.375D0*dd+45.d0)*dd
    c(6,1)=((86.625D0*dd+94.5D0)*dd+3.d0)*d
  END IF
!-
  IF(iord >= 5) THEN
    c(5,4)=9.d0*d
    c(5,3)=31.5D0*dd
    c(5,2)=(52.5D0*dd+5.d0)*d
    c(5,1)=(39.375D0*dd+21.d0)*dd
  END IF
!-
  IF(iord >= 4) THEN
    c(4,3)=7.d0*d
    c(4,2)=17.5D0*dd
    c(4,1)=(17.5D0*dd+3.d0)*d
  END IF
!-
  IF(iord >= 3) THEN
    c(3,2)=5.d0*d
    c(3,1)=7.5D0*dd
  END IF
!-
  c(2,1)=3.d0*d
END IF
IF(q < .9999 .OR. q > 1.0001) THEN
!-
!------- BUILT UP MATRIX A3=P*Q*P**-1
  a3(1,1)=q
  DO  k=2,iord
    a3(k,k)=a3(k-1,k-1)*q
  END DO
  ww=q*q-1.d0
  DO  k=1,iord-2
    ww=ww*q
    a3(k+2,k)=(k+.5D0)*ww
  END DO
!-
  IF(iord >= 5) THEN
    qq=a3(2,2)
    a3(5,1)=q*(1.875D0+qq*(-5.25D0+qq*3.375D0))
  END IF
  IF(iord >= 6) a3(6,2)=qq*(4.375D0+qq*(-11.25D0+qq*6.875D0))
  IF(iord == 7) THEN
    a3(7,1)=q*(-2.1875D0+qq*(11.8125D0+qq*(-18.5625D0+qq* 8.9375D0)))
    a3(7,3)=q*qq*(7.875D0+qq*(-19.25D0+qq*11.375D0))
  END IF
!-
!------- COMPUTE A*C AND NEW LEGENDRE SUMS
  DO  i=iord,2,-1
    xx=0.
    DO  k=1,i
      ww=0.
      DO  l=k,i-1,2
        ww=ww+a3(l,k)*c(i,l)
      END DO
      IF(MOD(i-k,2) == 0) ww=ww+a3(i,k)
      xx=xx+ww*sw(k)
    END DO
    sw(i)=xx
  END DO
  sw(1)=a3(1,1)*sw(1)
ELSE
  DO  i=iord,2,-1
    DO  k=1,i-1
      sw(i)=sw(i)+c(i,k)*sw(k)
    END DO
  END DO
END IF
RETURN
END SUBROUTINE lreg



SUBROUTINE freg(sw, nue, kord, iboun, y, c, icall, a)
!------------------------------------------------------------------
!    SHORT-VERSION: MAY, 1995

!    PURPOSE:

!    FINAL COMPUTATION OF A SMOOTHED VALUE VIA LEGENDRE POLYNOMIALS

!    PARAMETERS :
!                       **************   INPUT   *******************

!            SW(KORD+1):  SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.
!            NUE       :  ORDER OF DERIVATIVE
!            KORD      :  ORDER OF KERNEL
!            IBOUN     :  0: INTERIOR KERNEL, ELSE BOUNDARY KERNEL
!            C(KORD+1) :  SEQUENCE OF POLYN. COEFF. FOR BOUND. KERNEL
!            ICALL     :  PARAMETER USED TO INITIALISE COMPUTATION
!                      :   OF A MATRIX
!                       **************    WORK    ******************

!            A(7,7)    :  MATRIX OF COEFFICIENTS
!                       **************   OUTPUT   ******************

!            Y          :  COMPUTED ESTIMATE

!--------------------------------------------------------------------

REAL (dp), INTENT(IN)    :: sw(7)
INTEGER, INTENT(IN)      :: nue
INTEGER, INTENT(IN)      :: kord
INTEGER, INTENT(IN)      :: iboun
REAL (dp), INTENT(OUT)   :: y
REAL (dp), INTENT(IN)    :: c(7)
INTEGER, INTENT(IN OUT)  :: icall
REAL (dp), INTENT(OUT)   :: a(7,7)


INTEGER    :: i, j
REAL (dp)  :: ww
!-
!------- DEFINITION OF LEGENDRE COEFFICIENTS FOR BOUNDARY
IF(icall == 0 .AND. iboun /= 0) THEN
  a(2,2)=2./3.
  a(1,3)=.6
  a(3,3)=.4
  a(2,4)=4./7.
  a(4,4)=8./35.
  a(1,5)=27./63.
  a(3,5)=28./63.
  a(5,5)=8./63.
  a(2,6)=110./231.
  a(4,6)=72./231.
  a(6,6)=16./231.
  a(1,7)=143./429.
  a(3,7)=182./429.
  a(5,7)=88./429.
  a(7,7)=16./429.
  icall=1
END IF
IF(iboun /= 0) THEN
!-
!------- COMPUTATION OF THE SMOOTHED VALUE AT BOUNDARY
  y=c(1)*sw(1)+c(2)*a(2,2)*sw(2)
  DO  j=3,kord+1
    ww=a(j,j)*sw(j)
    DO  i=j-2,1,-2
      ww=ww+a(i,j)*sw(i)
    END DO
    y=y+c(j)*ww
  END DO
ELSE
!-
!------- COMPUTATION OF THE SMOOTHED VALUE AT INTERIOR
  IF(nue == 0) THEN
    IF(kord == 2) y=-.1*sw(3)+.6*sw(1)
    IF(kord == 4) y=(sw(5)-4.*sw(3)+9.*sw(1))/12.
    IF(kord == 6) y=-7.2115379E-02*sw(7)+.25961537*sw(5)  &
        -.4375*sw(3)+.75*sw(1)
  END IF
  IF(nue == 1) THEN
    IF(kord == 3) y=(3.*sw(4)-10.*sw(2))/14.
    IF(kord == 5) y=(-15.*sw(6)+48.*sw(4)-55.*sw(2))/44.
  END IF
  IF(nue == 2) THEN
    IF(kord == 4) y=(-5.*sw(5)+14.*sw(3)-9.*sw(1))/6.
    IF(kord == 6) y=2.01923*sw(7)-5.76923*sw(5)+5.25*sw(3) -1.5*sw(1)
  END IF
  IF(nue == 3) y=4.772727*sw(6)-12.272727*sw(4)+7.5*sw(2)
  IF(nue == 4) y=-36.34615*sw(7)+88.84615*sw(5)-52.5*sw(3)
END IF
!-
RETURN
END SUBROUTINE freg



SUBROUTINE kerncl(t, x, n, b, nue, kord, ny, s, tt, m, y)
!-----------------------------------------------------------------------
!       SHORT-VERSION JANUARY 1995

!       KERNEL SMOOTHING, CONVENTIONAL ALGORITHM,GENERAL
!       DESIGN, LOCAL BANDWIDTH ALLOWED

!  PARAMETERS :

!  INPUT    T(N)         INPUT GRID
!  INPUT    X(N)         DATA, GIVEN ON T(N)
!  INPUT    N            LENGTH OF X
!  INPUT    B            ONE SIDED BANDWIDTH (FOR NY=1 MEAN BANDWIDTH)
!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!  INPUT    KORD         ORDER OF KERNEL (<=6)
!  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
!  INPUT    NY           0: GLOBAL, 1: LOCAL BANDWIDTH INPUTED IN Y
!  INPUT    TT(M)        OUTPUT GRID
!  INPUT    M            NUMBER OF POINTS TO ESTIMATE
!  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
!  OUTPUT   Y(M)         ESTIMATED REGRESSION FUNCTION

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: b
INTEGER, INTENT(IN)        :: nue
INTEGER, INTENT(IN)        :: kord
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN)      :: s(0:n)
REAL (dp), INTENT(IN)      :: tt(:)
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN OUT)  :: y(:)


REAL (dp)  :: c(7), c1(7)
INTEGER    :: ist, i, iboun, iord
REAL (dp)  :: bb, bmax, wid, s0, s1, sn, bmin
!-
!------  COMPUTE KERNEL COEFFICIENTS FOR INTERIOR AND SOME CONSTANTS
CALL coffi(nue, kord, c)
iord=kord+1
bb=b
s0=1.5*t(1) - 0.5*t(2)
sn=1.5*t(n)-0.5*t(n-1)
bmin=(sn-s0)*.6D0/DBLE(n)*DBLE(kord-1)
bmax=(s(n)-s(0))*.5
IF(kord == 2) bmin=0.1D0*bmin
ist=1
!-
!-------  LOOP OVER OUTPUT GRID
DO  i=1,m
  IF(ny /= 0) bb=y(i)
  IF(bb > bmax) bb=bmax
  IF(bb < bmin) bb=bmin
  wid=bb
  s1=tt(i)-bb
  iboun=0
!-
!-------  COMPUTE LEFT BOUNDARY KERNEL
  IF(s1 < s(0)) THEN
    s1=s(0)
    wid=bb+bb+s(0)-tt(i)
    CALL coffb(nue,kord,(tt(i)-s(0))/wid,1,c1)
    iboun=1
  END IF
!-
!-------  COMPUTE RIGHT BOUNDARY KERNEL
  IF(tt(i)+bb > s(n)) THEN
    s1=s(n)-(bb+bb)
    wid=tt(i)-s1
    CALL coffb(nue,kord,(s(n)-tt(i))/wid,-1,c1)
    iboun=-1
  END IF
!-
!------  SEARCH FIRST S-POINT OF SMOOTHING INTERVAL
  2 IF(s(ist) <= s1) THEN
    ist=ist+1
    GO TO 2
  END IF
  3 IF(s(ist-1) > s1) THEN
    ist=ist-1
    GO TO 3
  END IF
!-
!-------  IF BANDWIDTH IS TOO SMALL NO SMOOTHING
  IF(s(ist) >= tt(i)+wid .OR. ist == n) THEN
    y(i)=x(ist)
    IF(nue > 0) y(i)=0.
  ELSE
!-
!-----  COMPUTE SMOOTHED DATA AT TT(I)
    IF (iboun /= 0) THEN
      CALL smo(s, x, n, tt(i), wid, nue, iord, iboun, ist, s1, c1, y(i))
    ELSE
      CALL smo(s, x, n, tt(i), wid, nue, iord, iboun, ist, s1, c, y(i))
    END IF
  END IF
END DO
!-
RETURN
END SUBROUTINE kerncl



SUBROUTINE smo(s, x, n, tau, wid, nue, iord, iboun, ist, s1, c, y)
!-----------------------------------------------------------------------
!       SHORT-VERSION JANUARY 1995

!       PERFORMS ONE SMOOTHING STEP

!  PARAMETERS :

!  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
!  INPUT    X(N)         DATA
!  INPUT    N            LENGTH OF X
!  INPUT    TAU          POINT WHERE FUNCTION IS ESTIMATED
!  INPUT    WID          ONE SIDED BANDWIDTH
!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!  INPUT    IORD          ORDER OF KERNEL POLYNOMIAL
!  INPUT    IBOUN         TYPE OF BOUNDARY
!  INPUT    IST          INDEX OF FIRST POINT OF SMOOTHING INTERVAL
!  INPUT    S1           LEFT BOUNDARY OF SMOOTHING INTERVAL
!  INPUT    C(7)         KERNEL COEFFICIENTS
!  OUTPUT   Y            SMOOTHED VALUE AT TAU

!  WORK     WO(7)        WORK ARRAY

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: s(0:)
REAL (dp), INTENT(IN)   :: x(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: tau
REAL (dp), INTENT(IN)   :: wid
INTEGER, INTENT(IN)     :: nue
INTEGER, INTENT(IN)     :: iord
INTEGER, INTENT(IN)     :: iboun
INTEGER, INTENT(IN)     :: ist
REAL (dp), INTENT(IN)   :: s1
REAL (dp), INTENT(IN)   :: c(7)
REAL (dp), INTENT(OUT)  :: y


REAL (dp)  :: wo(7)
INTEGER    :: jend, ibeg, incr, i, j
REAL (dp)  :: yy, yyy, w, widnue
!-
y=0.
jend=0
ibeg=2
IF(iboun /= 0 .OR. (nue /= 1 .AND. nue /= 3)) ibeg=1
incr=2
IF(iboun /= 0) incr=1
!-
!------  COMPUTE INITIAL KERNEL VALUES
IF(iboun > 0) THEN
  yy=(tau-s1)/wid
  wo(ibeg)=yy
  DO  i=ibeg,iord-incr,incr
    wo(i+incr)=wo(i)*yy
  END DO
ELSE
  DO  i=ibeg,iord,incr
    wo(i)=1.
  END DO
END IF
!-
!------  LOOP OVER SMOOTHING INTERVAL
DO  j=ist,n
  yy=(tau-s(j))/wid
  IF(yy < -1.) THEN
    yy=-1.
    jend=1
  END IF
  yyy=yy
  IF(iboun == 0) THEN
    yy=yy*yy
    IF(nue == 1 .OR. nue == 3) yyy=yy
  END IF
!-
!------  LOOP FOR COMPUTING WEIGHTS
  w=0.
  DO  i=ibeg,iord,incr
    w=w+c(i)*(wo(i)-yyy)
    wo(i)=yyy
    yyy=yyy*yy
  END DO
  y=y+w*x(j)
  IF(jend == 1) EXIT
END DO
!-
!-------  NORMALIZING FOR NUE>0
IF(nue > 0) THEN
  widnue=wid**nue
  y=y/widnue
END IF
!-
RETURN
END SUBROUTINE smo



SUBROUTINE kerncp(t, x, n, b, nue, kord, ny, s, tt, m, y)
!-----------------------------------------------------------------------
!       SHORT-VERSION JANUARY 1997

!       KERNEL SMOOTHING, CONVENTIONAL ALGORITHM,GENERAL
!       DESIGN, LOCAL BANDWIDTH ALLOWED, WITHOUT BOUNDARY KERNELS,
!       JUST NORMALIZING

!  PARAMETERS :

!  INPUT    T(N)         INPUT GRID
!  INPUT    X(N)         DATA, GIVEN ON T(N)
!  INPUT    N            LENGTH OF X
!  INPUT    B            ONE SIDED BANDWIDTH (FOR NY=1 MEAN BANDWIDTH)
!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!  INPUT    KORD         ORDER OF KERNEL (<=6)
!  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
!  INPUT    NY           0: GLOBAL, 1: LOCAL BANDWIDTH INPUTED IN Y
!  INPUT    TT(M)        OUTPUT GRID
!  INPUT    M            NUMBER OF POINTS TO ESTIMATE
!  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
!  OUTPUT   Y(M)         ESTIMATED REGRESSION FUNCTION

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)      :: t(:)
REAL (dp), INTENT(IN)      :: x(:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: b
INTEGER, INTENT(IN)        :: nue
INTEGER, INTENT(IN)        :: kord
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN)      :: s(0:n)
REAL (dp), INTENT(IN)      :: tt(:)
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN OUT)  :: y(:)

REAL (dp)  :: c(7), c1(7)
INTEGER    :: ist, i, iboun, iord
REAL (dp)  :: bb, bmax, wid, s0, s1, sn, bmin
!-
!------  COMPUTE KERNEL COEFFICIENTS FOR INTERIOR AND SOME CONSTANTS
CALL coffi(nue, kord, c)
iord=kord+1
bb=b
s0=1.5*t(1)-0.5*t(2)
sn=1.5*t(n)-0.5*t(n-1)
bmin=(sn-s0)*.6D0/DBLE(n)*DBLE(kord-1)
bmax=(s(n)-s(0))*.5
IF(kord == 2) bmin=0.1D0*bmin
ist=1
!-
!-------  LOOP OVER OUTPUT GRID
DO  i=1,m
  IF(ny /= 0) bb=y(i)
  IF(bb > bmax) bb=bmax
  IF(bb < bmin) bb=bmin
  wid=bb
  s1=tt(i)-bb
  iboun=0
!-
!-------  COMPUTE LEFT BOUNDARY KERNEL
  IF(s1 < s(0)) THEN
    s1=s(0)
    wid=bb+bb+s(0)-tt(i)
    CALL coffb(nue,kord,(tt(i)-s(0))/wid,1,c1)
    iboun=1
  END IF
!-
!-------  COMPUTE RIGHT BOUNDARY KERNEL
  IF(tt(i)+bb > s(n)) THEN
    s1=s(n)-(bb+bb)
    wid=tt(i)-s1
    iboun=-1
  END IF
!-
!------  SEARCH FIRST S-POINT OF SMOOTHING INTERVAL
  2 IF(s(ist) <= s1) THEN
    ist=ist+1
    GO TO 2
  END IF
  3 IF(s(ist-1) > s1) THEN
    ist=ist-1
    GO TO 3
  END IF
!-
!-------  IF BANDWIDTH IS TOO SMALL NO SMOOTHING
  IF(s(ist) >= tt(i)+wid .OR. ist == n) THEN
    y(i)=x(ist)
    IF(nue > 0) y(i)=0.
  ELSE
!-
!-----  COMPUTE SMOOTHED DATA AT TT(I)
    CALL smop(s, x, n, tt(i), wid, nue, iord, iboun, ist, s1, c, y(i))
  END IF
END DO
!-
RETURN
END SUBROUTINE kerncp



SUBROUTINE smop(s, x, n, tau, wid, nue, iord, iboun, ist, s1, c, y)
!-----------------------------------------------------------------------
!       SHORT-VERSION JANUARY 1997

!       PERFORMS ONE SMOOTHING STEP

!  PARAMETERS :

!  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
!  INPUT    X(N)         DATA
!  INPUT    N            LENGTH OF X
!  INPUT    TAU          POINT WHERE FUNCTION IS ESTIMATED
!  INPUT    WID          ONE SIDED BANDWIDTH
!  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
!  INPUT    IORD          ORDER OF KERNEL POLYNOMIAL
!  INPUT    IBOUN         TYPE OF BOUNDARY
!  INPUT    IST          INDEX OF FIRST POINT OF SMOOTHING INTERVAL
!  INPUT    S1           LEFT BOUNDARY OF SMOOTHING INTERVAL
!  INPUT    C(7)         KERNEL COEFFICIENTS
!  OUTPUT   Y            SMOOTHED VALUE AT TAU

!  WORK     WO(7)        WORK ARRAY

!-----------------------------------------------------------------------

REAL (dp), INTENT(IN)   :: s(0:)
REAL (dp), INTENT(IN)   :: x(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: tau
REAL (dp), INTENT(IN)   :: wid
INTEGER, INTENT(IN)     :: nue
INTEGER, INTENT(IN)     :: iord
INTEGER, INTENT(IN)     :: iboun
INTEGER, INTENT(IN)     :: ist
REAL (dp), INTENT(IN)   :: s1
REAL (dp), INTENT(IN)   :: c(7)
REAL (dp), INTENT(OUT)  :: y

REAL (dp)  :: wo(7)
INTEGER    :: jend, ibeg, incr, i, j
REAL (dp)  :: yy, yyy, w, widnue, ww
!-
y=0.
ww=0.
jend=0
ibeg=2
IF(nue /= 1 .AND. nue /= 3) ibeg=1
incr=2
!-
!------  COMPUTE INITIAL KERNEL VALUES
IF(iboun > 0) THEN
  yy=(tau-s1)/wid
  wo(ibeg)=yy
  yy=yy*yy
  IF(nue == 1 .OR. nue == 3) wo(ibeg)=yy
  DO  i=ibeg,iord-incr,incr
    wo(i+incr)=wo(i)*yy
  END DO
ELSE
  DO  i=ibeg,iord,incr
    wo(i)=1.
  END DO
END IF
!-
!------  LOOP OVER SMOOTHING INTERVAL
DO  j=ist,n
  yy=(tau-s(j))/wid
  IF(yy < -1.) THEN
    yy=-1.
    jend=1
  END IF
  yyy=yy
  yy=yy*yy
  IF(nue == 1 .OR. nue == 3) yyy=yy
!-
!------  LOOP FOR COMPUTING WEIGHTS
  w=0.
  DO  i=ibeg,iord,incr
    w=w+c(i)*(wo(i)-yyy)
    wo(i)=yyy
    yyy=yyy*yy
  END DO
  y=y+w*x(j)
  ww=ww+w
  IF(jend == 1) EXIT
END DO
!-
!-------  NORMALIZING FOR NUE>0
IF(ww /= 0) y=y/ww
IF(nue > 0) THEN
  widnue=wid**nue
  y=y/widnue
END IF
!-
RETURN
END SUBROUTINE smop



SUBROUTINE coffi(nue, kord, c)
!-----------------------------------------------------------------------
!       SHORT-VERSION: JANUARY 1995

!       PURPOSE:

!       DEFINES POLYNOMIAL KERNEL COEFFICIENTS FOR INTERIOR.

!  PARAMETERS:

!  INPUT  NUE        ORDER OF DERIVATIVE (0-4)
!  INPUT  KORD       ORDER OF KERNEL (NUE+I, I=2,4,6;  KORD<=6)
!  OUTPUT C(7)       POLYNOMIAL KERNEL COEFFICIENTS

!-----------------------------------------------------------------------

INTEGER, INTENT(IN)     :: nue
INTEGER, INTENT(IN)     :: kord
REAL (dp), INTENT(OUT)  :: c(7)

!-
c(1:7)=0.
IF(nue == 0 .AND. kord == 2) THEN
  c(1)=0.75_dp
  c(3)=-0.25_dp
END IF

IF(nue == 0 .AND. kord == 4) THEN
  c(1)=1.40625D0
  c(3)=-1.5625D0
  c(5)=0.65625D0
END IF

IF(nue == 0 .AND. kord == 6) THEN
  c(1)=2.05078125D0
  c(3)=-4.78515625D0
  c(5)=5.16796875D0
  c(7)=-1.93359375D0
END IF

IF(nue == 1 .AND. kord == 3) THEN
  c(2)=-1.875D0
  c(4)=0.9375D0
END IF

IF(nue == 1 .AND. kord == 5) THEN
  c(2)=-8.203125D0
  c(4)=11.484375D0
  c(6)=-4.921875D0
END IF

IF(nue == 2 .AND. kord == 4) THEN
  c(1)=-6.5625D0
  c(3)=13.125D0
  c(5)=-6.5625D0
END IF

IF(nue == 2 .AND. kord == 6) THEN
  c(1)=-24.609375D0
  c(3)=103.359375D0
  c(5)=-132.890625D0
  c(7)=54.140625D0
END IF

IF(nue == 3 .AND. kord == 5) THEN
  c(2)=88.59375D0
  c(4)=-147.65625D0
  c(6)=68.90625D0
END IF

IF(nue == 4 .AND. kord == 6) THEN
  c(1)=324.84375D0
  c(3)=-1624.21875D0
  c(5)=2273.90625D0
  c(7)=-974.53125D0
END IF

RETURN
END SUBROUTINE coffi



SUBROUTINE coffb(nue, kord, q, iboun, c)
!-----------------------------------------------------------------------
!       SHORT-VERSION: JANUARY 1995

!       PURPOSE:

!       COMPUTES COEFFICIENTS OF POLYNOMIAL BOUNDARY KERNELS,
!       FOLLOWING GASSER + MUELLER PREPRINT 38 SFB 123 HEIDELBERG
!       AND UNPUBLISHED RESULTS

!  PARAMETERS:

!  INPUT  NUE        ORDER OF DERIVATIVE (0-4)
!  INPUT  KORD       ORDER OF KERNEL (NUE+I, I=2,4,6;  KORD<=6)
!  INPUT  Q          PERCENTAGE OF WID AT BOUNDARY
!  INPUT  IBOUN      < 0 RIGHT BOUNDARY OF DATA
!                    > 0 LEFT BOUNDARY OF DATA
!  OUTPUT C(7)       POLYNOMIAL KERNEL COEFFICIENTS

!-----------------------------------------------------------------------

INTEGER, INTENT(IN)     :: nue
INTEGER, INTENT(IN)     :: kord
REAL (dp), INTENT(IN)   :: q
INTEGER, INTENT(IN)     :: iboun
REAL (dp), INTENT(OUT)  :: c(7)

! Local variables

REAL (dp)  :: d, p, p1, p3, p6, p12
INTEGER    :: i, j

c(1:7)=0.
p=-q
p1=1.+q
p3=p1*p1*p1

IF(nue == 0 .AND. kord == 2) THEN
  d=1./(p3*p1)
  c(1)=(6.+p*(12.+p*18.))*d
  c(2)=9.*(1.+p)*(1.+p)*d
  c(3)=(4.+p*8.)*d
END IF

IF(nue == 0 .AND. kord == 4) THEN
  d=p1/(p3*p3*p3)
  p12=(1.+p)*(1.+p)*d
  c(1)=20.*(1.+p*(12.+p*(78.+p*(164.+p*(165.+p*(60.+p*10.))))))*d
  c(2)=100.*(1.+p*(5.+p))**2*p12
  c(3)=200.*(1.+p*(12.+p*(33.+p*(36.+p*(14.+p+p)))))*d
  c(4)=175.*(1.+p*(10.+p*3.))*p12
  c(5)=56.*(1.+p*(12.+p*(18.+p*4.)))*d
END IF

IF(nue == 0 .AND. kord == 6) THEN
  p6=p3*p3
  d=1./(p6*p6)
  p12=(1.+p)*(1.+p)*d
  c(1)=42.*(1.+p*(30.+p*(465.+p*(3000.+p*(10050.+p*(17772.+p  &
      *(17430.+p*(9240.+p*(2625.+p*(350.+p*21.))))))))))*d
  c(2)=441.*(1.+p*(14.+p*(36.+p*(14.+p))))**2*p12
  c(3)=1960.*(1.+p*(30.+p*(255.+p*(984.+p*(1902.+p*(1956.+p  &
      *(1065.+p*(300.+p*(39.+p+p)))))))))*d
  c(4)=4410.*(1.+p*(28.+p*(156.+p*(308.+p*(188.+p*(42.+p*3.)))))) *p12
  c(5)=5292.*(1.+p*(30.+p*(185.+p*(440.+p*(485.+p*(250.+p*(57. +p*4.)))))))*d
  c(6)=3234.*(1.+p*(28.+p*(108.+p*(56.+p*5.))))*p12
  c(7)=792.*(1.+p*(30.+p*(150.+p*(200.+p*(75.+p*6.)))))*d
END IF

IF(nue == 1 .AND. kord == 3) THEN
  d=-1./(p3*p3)
  p12=(1.+p)*(1.+p)*d
  c(1)=(60.+p*240.)*p12
  c(2)=120.*(2.+p*(6.+p*(6.+p)))*d
  c(3)=300.*p12
  c(4)=(120.+p*180.)*d
END IF

IF(nue == 1 .AND. kord == 5) THEN
  d=-1./(p3*p3*p3*p1)
  p12=(1.+p)*(1.+p)*d
  c(1)=420.*(1.+p*(18.+p*(98.+p*(176.+p*(75.+p*10.)))))*p12
  c(2)=2100.*(2.+p*(25.+p*(120.+p*(245.+p*(238.+p*(105.+p*(20. +p)))))))*d
  c(3)=14700.*(1.+p*(4.+p))**2*p12
  c(4)=5880.*(4.+p*(35.+p*(90.+p*(95.+p*(40.+p*6.)))))*d
  c(5)=17640.*(1.+p*(6.+p+p))*p12
  c(6)=2520.*(2.+p*(15.+p*(20.+p*5.)))*d
END IF

IF(nue == 2 .AND. kord == 4) THEN
  d=p1/(p3*p3*p3)
  p12=(1.+p)*(1.+p)*d
  c(1)=840.*(1.+p*(12.+p*(28.+p*(24.+p*5.))))*d
  c(2)=2100.*(3.+p*(10.+p))*p12
  c(3)=1680.*(9.+p*(28.+p*(27.+p*6.)))*d
  c(4)=14700.*p12
  c(5)=(5040.+p*6720.)*d
END IF

IF(nue == 2 .AND. kord == 6) THEN
  p6=p3*p3
  d=1./(p6*p6)
  p12=(1.+p)*(1.+p)*d
  c(1)=5040.*(2.+p*(60.+p*(489.+p*(1786.+p*(3195.+p*(2952.+p  &
      *(1365.+p*(294.+p*21.))))))))*d
  c(2)=52920.*(3.+p*(42.+p*(188.+p*(308.+p*(156.+p*(28. +p))))))*p12
  c(3)=141120.*(6.+p*(68.+p*(291.+p*(570.+p*(555.+p*(264.+p  &
      *(57.+p*4.)))))))*d
  c(4)=529200.*(2.+p*(7.+p+p))**2*p12
  c(5)=90720.*(30.+p*(228.+p*(559.+p*(582.+p*(255. +p*40.)))))*d
  c(6)=582120.*(3.+p*(14.+p*5.))*p12
  c(7)=221760.*(2.+p*(12.+p*(15.+p*4.)))*d
END IF

IF(nue == 3 .AND. kord == 5) THEN
  d=-1./(p3*p3*p3*p1)
  p12=(1.+p)*(1.+p)*d
  c(1)=15120.*(1.+p*(18.+p*(38.+p*6.)))*p12
  c(2)=45360.*(4.+p*(35.+p*(80.+p*(70.+p*(20.+p)))))*d
  c(3)=352800.*(2.+p*(6.+p))*p12
  c(4)=151200.*(8.+p*(25.+p*(24.+p*6.)))*d
  c(5)=952560.*p12
  c(6)=70560.*(4.+p*5.)*d
END IF

IF(nue == 4 .AND. kord == 6) THEN
  p6=p3*p3
  d=1./(p6*p6)
  p12=(1.+p)*(1.+p)*d
  c(1)=332640.*(1.+p*(30.+p*(171.+p*(340.+p*(285.+p*(90.+p*7. ))))))*d
  c(2)=1164240.*(5.+p*(56.+p*(108.+p*(28.+p))))*p12
  c(3)=6652800.*(5.+p*(38.+p*(85.+p*(76.+p*(25.+p+p)))))*d
  c(4)=17463600.*(5.+p*(14.+p*3.))*p12
  c(5)=4656960.*(25.+p*(78.+p*(75.+p*20.)))*d
  c(6)=76839840.*p12
  c(7)=3991680.*(5.+p*6.)*d
END IF

IF(iboun > 0) RETURN
j=2
IF(nue == 1 .OR. nue == 3) j=1
DO  i=j,kord,2
  c(i)=-c(i)
END DO

RETURN
END SUBROUTINE coffb



SUBROUTINE coff(x, n, fa)

REAL (dp), INTENT(OUT)  :: x(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: fa

x(1:n)=fa

RETURN
END SUBROUTINE coff

END MODULE Kernel_Regression
