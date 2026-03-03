MODULE MV_Normal

!      ALGORITHM 725, COLLECTED ALGORITHMS FROM ACM.
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-05  Time: 09:51:33
 
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 19, NO. 4, DECEMBER, 1993, P. 546.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)


CONTAINS


SUBROUTINE dmv(m, k, h, r, prob, eps, ier, ERR, cut)

INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN OUT)  :: k
REAL (dp), INTENT(IN)    :: h(20)
REAL (dp), INTENT(IN)    :: r(20,20)
REAL (dp), INTENT(OUT)   :: prob
REAL (dp), INTENT(IN)    :: eps
INTEGER, INTENT(OUT)     :: ier
REAL (dp), INTENT(OUT)   :: ERR
REAL (dp), INTENT(IN)    :: cut


!   DMV CALCULATES THE MULTIVARIATE NORMAL INTEGRAL.

!   THE INTEGRAL LOWER LIMITS ARE -INFINITY FOR ALL VARIABLES.
!   THE UPPER LIMITS ARE GIVEN BY THE VECTOR H (INPUT).
!   THE CORRELATION MATRIX IS R (INPUT).

!   PARAMETERS:

!   M IS THE NUMBER OF VARIABLES.  M <= 20 (INPUT).

!   K IS A FLAG CONTROLLING THE INTEGRATION METHOD (INPUT).
!     2<=K<=10 ONE GAUSIAN QUADRATURE WITH K POINTS FOR EACH
!              VARIABLE IS PERFORMED.  EPS IS DISREGARDED.
!    12<=K<=20 PROGRESSIVE GAUSIAN QUADRATURES ARE PERFORMED FOR
!              K=2,3,...,KMAX, WHERE KMAX IS 10 LESS THAN THE INPUT
!              VALUE OF K, UNTIL THE DIFFERENCE BETWEEN THE CALCULATIONS
!              FOR TWO SUCCESSIVE K'S DOES NOT EXCEED EPS.
!              THE LAST VALUE OF K IS RETURNED IN IER.
!              IF KMAX IS REACHED BEFORE A DIFFERENCE OF EPS IS ACHIEVED,
!              THEN IER=0 AND THIS DIFFERENCE IS RETURNED IN ERR.
!     EVERY OTHER K IS EQUIVALENT TO K=20.

!   H IS THE VECTOR OF UPPER LIMITS FOR THE INTEGRATION (INPUT).

!   R IS THE CORRELATION MATRIX DIMENSIONED R(20,20) (INPUT).
!     IT MUST BE A SYMMETRIC, POSITIVE SEMI-DEFINITE MATRIX.
!     THE CODE DOES NOT CHECK FOR THAT.
!     IF R IS SINGULAR, NUMERICAL PROBLEMS MAY OCCUR.

!   PROB IS THE CALCULATED PROBABILITY (OUTPUT).

!   EPS IS THE REQUIRED ERROR LIMIT (INPUT).

!   IER IS AN ERROR INDICATOR. ON NORMAL TERMINATION 0 < IER <= 10 IS THE
!       VALUE OF K CORRESPONDING TO THE RETURNED VALUE PROB, (OUTPUT).
!       IER= 0 IS RETURNED WHEN THE ERROR IS GREATER THAN EPS.
!       IER=98 MEANS AN ILLEGAL M VALUE (M SHOULD NOT EXCEED 20).
!       IER=99 MEANS THAT R IS SINGULAR OR ITS DETERMINANT IS NEGATIVE.

!   ERR IS THE DIFFERENCE IN CALCULATING BY K-1 AND K. CALCULATED
!       ONLY WHEN THE INPUT K > 10 (OUTPUT).

!   CUT - IF X<-CUT THEN EXP(X) IS IGNORED (ASSUMED ZERO).
!         SUGGESTED VALUE CUT=15 (INPUT).

!   TWO SUBROUTINES (DPOFA, DPODI) ARE TAKEN FROM THE NAAS:LINPACK
!   LIBRARY. THEY CAN BE REPLACED BY A SUBROUTINE THAT INVERTS A
!   SYMMETRIC MATRIX AND CALCULATES ITS DETERMINANT.

INTEGER    :: ker, kmax
REAL (dp)  :: pold

IF(k < 2 .OR. k > 10) GO TO 10

!     DMV1 CALCULATES THE PROBABILITY FOR A GIVEN K

CALL dmv1(m, k, h, r, cut, prob, ier)
IF(ier == 0) ier=k
RETURN

10   kmax=10
IF(k > 11 .AND. k < 21) kmax = k - 10
pold=-1.
DO  ker=2,kmax
  CALL dmv1(m, ker, h, r, cut, prob, ier)
  IF(ier > 0) RETURN
  ier=ker
  ERR=prob-pold
  IF(ABS(ERR) <= eps) RETURN
  pold=prob
END DO
ier=0
RETURN
END SUBROUTINE dmv


SUBROUTINE dmv1(m, k, h, r, cut, res, ier)

!   DMV1 CALCULATES THE PROBABILITY FOR A GIVEN K BY APPLYING THEOREM 1.
!   IT CALLS SUBROUTINE SOLVE TO CALCULATE INDIVIDUAL VALUES BY (1) and (2).
!   DMV1 CAN BE CALLED DIRECTLY.

!   K IS THE NUMBER OF POINTS IN THE GAUSS QUADRATURE (INPUT).
!       IF K < 2 OR K > 10 THEN K IS SET TO 10.

!   IER IS AN ERROR INDICATOR.
!       IER=0 IS NORMAL TERMINATION.
!       IER=98 INDICATES THAT M IS OUT OF RANGE.
!       IER=99 MEANS THAT R IS SINGULAR OR ITS DETERMINANT IS NEGATIVE.

!   THE OTHER VARIABLES ARE DESCRIBED ABOVE IN DMV.

INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN OUT)  :: k
REAL (dp), INTENT(IN)    :: h(20)
REAL (dp), INTENT(IN)    :: r(20,20)
REAL (dp), INTENT(IN)    :: cut
REAL (dp), INTENT(OUT)   :: res
INTEGER, INTENT(OUT)     :: ier

REAL (dp)  :: h1(20), r1(20,20), s
INTEGER    :: i, ipos, j, l, l1, l2, list(20), list1(20)

ier=0
IF(m < 0 .OR. m > 20) ier=98
IF(ier > 0) RETURN
res=0.

!   WE EXPAND UP TO 2 TO THE POWER OF M TERMS DESCRIBED IN THEOREM 1.
!   THE VECTOR LIST REPRESENTS THE "+" "-" AND INFINITY IN THE TERM.
!   LIST(I) HAS THE VALUE OF 0 IF H(I) <= 0, 1 IF H(I) > 0, AND
!   THE VALUE OF 2 IF INFINITY IS USED AS THE UPPER LIMIT AND THE
!   VARIABLE IS DROPPED FROM THE INTEGRATION.
!   LIST1 IS THE LIST OF VARIABLES THAT PARTICIPATE IN THE
!   INTEGRATION (THAT DO NOT HAVE INFINITY AS AN UPPER LIMIT).

DO  i=1,m
  list(i)=0
  IF(h(i) > 0.) list(i)=1
END DO
20   ipos=0
l=0
DO  i=1,m
  IF(list(i) == 2) CYCLE
  l=l+1
  list1(l)=i
END DO
IF(l == 0) GO TO 50
DO  i=1,l
  l1=list1(i)
  h1(i)=h(l1)
  IF(list(l1) == 1) h1(i)=-h1(i)
  IF(list(l1) == 1) ipos=1-ipos
  DO  j=1,l
    l2=list1(j)
    r1(i,j)=r(l1,l2)
    IF(list(l1) == 1) r1(i,j)=-r1(i,j)
    IF(list(l2) == 1) r1(i,j)=-r1(i,j)
  END DO
END DO

50 CALL solve(l, k, h1, r1, cut, s, ier)
IF(ier > 0) RETURN
IF(ipos == 1) s = -s
res = res + s
DO  i=1,m
  IF(list(i) == 0) CYCLE
  list(i) = list(i) + 1
  IF(list(i) == 2) GO TO 20
  list(i) = 1
END DO
RETURN
END SUBROUTINE dmv1


SUBROUTINE solve(m, k, h, r, cut, solution, ier)

!     SOLVE CALCULATES THE MULTIVARIATE PROBABILITY OF A GIVEN VECTOR H
!     AND GIVEN K.  SOLVE CAN BE CALLED DIRECTLY.

!     ALL THE ARGUMENTS ARE AS DESCRIBED ABOVE IN DMV1 OR DMV.

!     THE INTEGRATION PARAMETERS TAKEN FROM [8] HAVE 15 DIGITS ACCURACY.

INTEGER, INTENT(IN)      :: m
INTEGER, INTENT(IN OUT)  :: k
REAL (dp), INTENT(IN)    :: h(20)
REAL (dp), INTENT(IN)    :: r(20,20)
REAL (dp), INTENT(IN)    :: cut
REAL (dp), INTENT(OUT)   :: solution
INTEGER, INTENT(OUT)     :: ier

REAL (dp)  :: d(2), det, h1(20), h2(20), pi, prod, r1(20,20), sum, y(20)
INTEGER    :: list(20)
INTEGER    :: i, j, l

!     THE VECTORS COEF AND X CONTAIN THE GAUSS QUADRATURE COEFFICIENTS
!     AND POINTS RESPECTIVELY AND ARE OBTAINED FROM [6].

REAL(dp), PARAMETER  :: coef(10,10) = RESHAPE( (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, .6405291796843790D0,  &
    .245697745768379D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,  &
    0.0D0, .446029770466658D0, .396468266998335D0, 4.37288879877644D-2,  &
    0.d0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, .325302999756919D0,  &
    .421107101852062D0, .133442500357520D0, .637432348625728D-2, 0.d0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, .248406152028443D0, .392331066652399D0,   &
    .211418193076057D0, .332466603513439D-1, .824853344515628D-3, 0.0D0,  &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, .196849675488598D0, .349154201525395D0, &
    .257259520584421D0, .0760131375840057D0, .685191862513596D-2,  &
    9.84716452019267D-5, 0.0D0, 0.0D0, 0.0D0, 0.0D0, .160609965149261D0,  &
    .306319808158099D0, .275527141784905D0, .120630193130784D0,   &
    .0218922863438067D0, .00123644672831056D0, 1.10841575911059D-5, 0.0D0,  &
    0.0D0, 0.0D0, .13410918845336D0, .26833075447264D0, .275953397988422D0,  &
    .15744828261879D0, .0448141099174625D0, .00536793575602526D0,   &
    2.02063649132407D-4, 1.19259692659532D-6, 0.d0, 0.0D0, .114088970242118D0, &
    .235940791223685D0, .266425473630253D0, .183251679101663D0,   &
    .0713440493066916D0, .0139814184155604D0, .00116385272078519D0,  &
    .305670214897831D-4, 1.23790511337496D-7, 0.d0, 0.0985520975191087D0,  &
    .208678066608185D0, .252051688403761D0, .198684340038387D0,  &
    .097198422760062D0, .0270244164355446D0, .00380464962249537D0,   &
    2.28886243044656D-4, 4.34534479844469D-6, 1.24773714817825D-8 /),  &
    (/ 10, 10 /) )
REAL(dp), PARAMETER  :: x(10,10) = RESHAPE( (/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, .300193931060839D0,  &
    1.25242104533372D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,  &
    0.0D0, .190554149798192D0, .848251867544577D0, 1.79977657841573D0,  &
    0.d0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, .133776446996068D0,  &
    .624324690187190D0, 1.34253782564499D0, 2.26266447701036D0, 0.d0,  &
    0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, .100242151968216D0,  &
    .482813966046201D0, 1.06094982152572D0, 1.77972941852026D0,  &
    2.66976035608766D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, .0786006594130979D0,&
    .386739410270631D0, .866429471682044D0, 1.46569804966352D0,   &
    2.172707796939D0, 3.03682016932287D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,   &
    .0637164846067008D0, .318192018888619D0, .724198989258373D0,   &
    1.23803559921509D0, 1.83852822027095D0, 2.53148815132768D0,   &
    3.37345643012458D0, 0.0D0, 0.0D0, 0.0D0, .0529786439318514D0,  &
    .267398372167767D0, .616302884182402D0, 1.06424631211623D0,  &
    1.58885586227006D0, 2.18392115309586D0, 2.86313388370808D0,  &
    3.6860071627244D0, 0.d0, 0.0D0, .0449390308011934D0, .228605305560535D0, &
    .532195844331646D0, .927280745338081D0, 1.39292385519588D0,  &
    1.91884309919743D0, 2.50624783400574D0, 3.17269213348124D0,  &
    3.97889886978978D0, 0.d0, .0387385243257289D0, .198233304013083D0,  &
    .465201111814767D0, .816861885592273D0, 1.23454132402818D0,  &
    1.70679814968913D0, 2.22994008892494D0, 2.80910374689875D0,  &
    3.46387241949586D0, 4.25536180636608D0 /), (/ 10, 10 /) )

pi = 4.0D0*ATAN(1.d0)
ier=0
IF(m < 0 .OR. m > 20) ier=98
IF(ier > 0) RETURN
IF(k < 2 .OR. k > 10) k=10
IF(m > 0) GO TO 10
solution=1.
ier=0
RETURN

10 r1(1:m,1:m) = r(1:m,1:m)
CALL dpofa(r1, 20, m, ier)
IF(ier > 0) ier=99
IF(ier == 99) RETURN
CALL dpodi(r1, 20, m, d, 11)
IF(d(1) <= 0. .OR. d(2) <= -30) ier=99
IF(ier == 99) RETURN

!   NOTE THAT ANOTHER SUBROUTINE MAY GIVE THE DETERMINANT DIRECTLY
!   IN DET AND THEREFORE THE NEXT STATEMENT MAY BE UNNECESSARY.
!   R1 CONTAINS THE INVERSE OF R.  LOOP 30 CATERS TO THE POSSIBILITY
!   THAT THE INVERSE IS STORED ONLY IN THE LOWER TRIANGLE.

det=d(1)*10.d0**d(2)
DO  i=1,m
  DO  j=i,m
    r1(j,i)=r1(i,j)
  END DO
END DO
prod=1.d0
DO  i=1,m
  h1(i)=SQRT(2./r1(i,i))
  h2(i)=h(i)/h1(i)
  prod=prod*pi*r1(i,i)
  list(i)=1
END DO
prod=1.d0 / SQRT(det*prod)
DO  i=1,m
  DO  j=1,m
    r1(i,j)=r1(i,j)*h1(i)*h1(j)*0.5D0
  END DO
END DO

solution=0.
60   sum=0.
DO  i=1,m
  l=list(i)
  sum=sum + x(l,k)*x(l,k)
  y(i) = x(l,k) - h2(i)
END DO
DO  i=1,m
  DO  j=1,m
    sum = sum - y(i)*y(j)*r1(i,j)
  END DO
END DO
IF(sum < -cut) GO TO 100
sum=EXP(sum)
DO  i=1,m
  sum=sum*coef(list(i),k)
END DO
solution = solution + sum
100  DO  i=1,m
  list(i)=list(i)+1
  IF(list(i) <= k) GO TO 60
  list(i)=1
END DO
solution = solution*prod
RETURN
END SUBROUTINE solve


SUBROUTINE dpofa(a, lda, n, info)

INTEGER, INTENT(IN)        :: lda
REAL (dp), INTENT(IN OUT)  :: a(lda,:)
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(OUT)       :: info

!   dpofa factors a double precision symmetric positive definite matrix.

!   dpofa is usually called by dpoco, but it can be called
!   directly with a saving in time if  rcond  is not needed.
!   (time for dpoco) = (1 + 18/n)*(time for dpofa) .

!   on entry

!      a       double precision(lda, n)
!              the symmetric matrix to be factored.  Only the
!              diagonal and upper triangle are used.

!      lda     integer
!              the leading dimension of the array  a .

!      n       integer
!              the order of the matrix  a .

!   on return

!      a       an upper triangular matrix  r  so that  a = trans(r)*r
!              where  trans(r)  is the transpose.
!              the strict lower triangle is unaltered.
!              if  info .ne. 0 , the factorization is not complete.

!      info    integer
!              = 0  for normal return.
!              = k  signals an error condition.  The leading minor
!                   of order  k  is not positive definite.

!   linpack.  this version dated 08/14/78 .
!   cleve moler, university of new mexico, argonne national lab.

!   subroutines and functions

!   blas ddot
!   fortran dsqrt

!   internal variables

REAL (dp)  :: s, t
INTEGER    :: j, jm1, k

!     begin block with ...exits to 40

DO  j = 1, n
  info = j
  s = 0.0D0
  jm1 = j - 1
  IF (jm1 < 1) GO TO 20
  DO  k = 1, jm1
    t = a(k,j) - ddot(k-1, a(1:,k), 1, a(1:,j), 1)
    t = t/a(k,k)
    a(k,j) = t
    s = s + t*t
  END DO
  20 s = a(j,j) - s
!     ......exit
  IF (s <= 0.0D0) GO TO 40
  a(j,j) = SQRT(s)
END DO
info = 0
40 RETURN
END SUBROUTINE dpofa


SUBROUTINE dpodi(a, lda, n, det, job)

INTEGER, INTENT(IN)        :: lda
REAL (dp), INTENT(IN OUT)  :: a(lda,:)
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(OUT)     :: det(2)
INTEGER, INTENT(IN)        :: job

!   dpodi computes the determinant and inverse of a certain
!   double precision symmetric positive definite matrix (see below)
!   using the factors computed by dpoco, dpofa or dqrdc.

!   on entry

!      a       double precision(lda, n)
!              the output  a  from dpoco or dpofa
!              or the output  x  from dqrdc.

!      lda     integer
!              the leading dimension of the array  a .

!      n       integer
!              the order of the matrix  a .

!      job     integer
!              = 11   both determinant and inverse.
!              = 01   inverse only.
!              = 10   determinant only.

!   on return

!      a       if dpoco or dpofa was used to factor  a  then
!              dpodi produces the upper half of inverse(a) .
!              if dqrdc was used to decompose  x  then
!              dpodi produces the upper half of inverse(trans(x)*x)
!              where trans(x) is the transpose.
!              elements of  a  below the diagonal are unchanged.
!              if the units digit of job is zero,  a  is unchanged.

!      det     double precision(2)
!              determinant of  a  or of  trans(x)*x  if requested.
!              otherwise not referenced.
!              determinant = det(1) * 10.0**det(2)
!              with  1.0 .le. det(1) .lt. 10.0
!              or  det(1) .eq. 0.0 .

!   error condition

!      a division by zero will occur if the input factor contains
!      a zero on the diagonal and the inverse is requested.
!      it will not occur if the subroutines are called correctly
!      and if dpoco or dpofa has set info .eq. 0 .

!   linpack.  this version dated 08/14/78 .
!   cleve moler, university of new mexico, argonne national lab.

!   subroutines and functions

!   blas daxpy,dscal
!   fortran mod

!   internal variables

REAL (dp)  :: s, t
INTEGER    :: i, j, jm1, k, kp1

!     compute determinant

IF (job/10 == 0) GO TO 70
det(1) = 1.0D0
det(2) = 0.0D0
s = 10.0D0
DO  i = 1, n
  det(1) = a(i,i)**2*det(1)
!        ...exit
  IF (det(1) == 0.0D0) EXIT
  10  IF (det(1) >= 1.0D0) GO TO 30
  det(1) = s*det(1)
  det(2) = det(2) - 1.0D0
  GO TO 10
  30  IF (det(1) < s) CYCLE
  det(1) = det(1)/s
  det(2) = det(2) + 1.0D0
  GO TO 30
END DO

!     compute inverse(r)

70 IF (MOD(job,10) == 0) GO TO 140
DO  k = 1, n
  a(k,k) = 1.0D0/a(k,k)
  t = -a(k,k)
  CALL dscal(k-1, t, a(1:,k), 1)
  kp1 = k + 1
  IF (n >= kp1) THEN
    DO  j = kp1, n
      t = a(k,j)
      a(k,j) = 0.0D0
      CALL daxpy(k, t, a(1:,k), 1, a(1:,j), 1)
    END DO
  END IF
END DO

!        form  inverse(r) * trans(inverse(r))

DO  j = 1, n
  jm1 = j - 1
  IF (jm1 >= 1) THEN
    DO  k = 1, jm1
      t = a(k,j)
      CALL daxpy(k, t, a(1:,j), 1, a(1:,k), 1)
    END DO
  END IF
  t = a(j,j)
  CALL dscal(j, t, a(1:,j), 1)
END DO

140 RETURN
END SUBROUTINE dpodi


SUBROUTINE daxpy(n, da, dx, incx, dy, incy)

!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: da
REAL (dp), INTENT(IN)   :: dx(:)
INTEGER, INTENT(IN)     :: incx
REAL (dp), INTENT(OUT)  :: dy(:)
INTEGER, INTENT(IN)     :: incy

INTEGER  :: i, ix, iy, m, mp1

IF(n <= 0) RETURN
IF (da == 0.0D0) RETURN
IF(incx == 1 .AND. incy == 1) GO TO 20

!        code for unequal increments or equal increments not equal to 1

ix = 1
iy = 1
IF(incx < 0) ix = (-n+1)*incx + 1
IF(incy < 0) iy = (-n+1)*incy + 1
DO  i = 1,n
  dy(iy) = dy(iy) + da*dx(ix)
  ix = ix + incx
  iy = iy + incy
END DO
RETURN

!        code for both increments equal to 1

!        clean-up loop

20 m = MOD(n,4)
IF( m == 0 ) GO TO 40
DO  i = 1,m
  dy(i) = dy(i) + da*dx(i)
END DO
IF( n < 4 ) RETURN
40 mp1 = m + 1
DO  i = mp1,n,4
  dy(i) = dy(i) + da*dx(i)
  dy(i + 1) = dy(i + 1) + da*dx(i + 1)
  dy(i + 2) = dy(i + 2) + da*dx(i + 2)
  dy(i + 3) = dy(i + 3) + da*dx(i + 3)
END DO
RETURN
END SUBROUTINE daxpy


SUBROUTINE dscal(n, da, dx, incx)

!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: da
REAL (dp), INTENT(OUT)  :: dx(:)
INTEGER, INTENT(IN)     :: incx

INTEGER :: i, m, mp1, nincx

IF(n <= 0) RETURN
IF(incx == 1) GO TO 20

!        code for increment not equal to 1

nincx = n*incx
DO  i = 1,nincx,incx
  dx(i) = da*dx(i)
END DO
RETURN

!        code for increment equal to 1


!        clean-up loop

20 m = MOD(n,5)
IF( m == 0 ) GO TO 40
DO  i = 1,m
  dx(i) = da*dx(i)
END DO
IF( n < 5 ) RETURN
40 mp1 = m + 1
DO  i = mp1,n,5
  dx(i) = da*dx(i)
  dx(i + 1) = da*dx(i + 1)
  dx(i + 2) = da*dx(i + 2)
  dx(i + 3) = da*dx(i + 3)
  dx(i + 4) = da*dx(i + 4)
END DO
RETURN
END SUBROUTINE  dscal


FUNCTION ddot(n, dx, incx, dy, incy) RESULT(fn_val)

!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: dx(:)
INTEGER, INTENT(IN)    :: incx
REAL (dp), INTENT(IN)  :: dy(:)
INTEGER, INTENT(IN)    :: incy
REAL (dp)              :: fn_val

REAL (dp)  :: dtemp
INTEGER    :: i, ix, iy, m, mp1

fn_val = 0.0D0
dtemp = 0.0D0
IF(n <= 0) RETURN
IF(incx == 1 .AND. incy == 1) GO TO 20

!        code for unequal increments or equal increments
!          not equal to 1

ix = 1
iy = 1
IF(incx < 0)ix = (-n+1)*incx + 1
IF(incy < 0)iy = (-n+1)*incy + 1
DO  i = 1,n
  dtemp = dtemp + dx(ix)*dy(iy)
  ix = ix + incx
  iy = iy + incy
END DO
fn_val = dtemp
RETURN

!        code for both increments equal to 1

!        clean-up loop

20 m = MOD(n,5)
IF( m == 0 ) GO TO 40
DO  i = 1,m
  dtemp = dtemp + dx(i)*dy(i)
END DO
IF( n < 5 ) GO TO 60
40 mp1 = m + 1
DO  i = mp1,n,5
  dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +  &
      dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
END DO
60 fn_val = dtemp
RETURN
END FUNCTION ddot

END MODULE MV_Normal



PROGRAM test725

!   THIS TEST PROGRAM CHECKS DMV ON A SET OF PRE-SPECIFIED PROBLEMS.
!   FOR EACH PROBLEM, THE RESULT CALCULATED BY THE PROGRAM IS PRINTED
!   NEXT TO THE CORRECT RESULT SO THAT ONE CAN CHECK HIS MACHINE.
!   NO INPUT IS REQUIRED.

USE MV_Normal
IMPLICIT NONE

REAL (dp)  :: diff, err, h(20), prob, r(20,20)
INTEGER    :: ier, ih, ir, j, k, l, m, number
REAL(dp), PARAMETER  :: result(12) = (/ .72437, .74520, .95618, .95855,   &
                                        .63439, .67778, .93656, .94253,   &
                                        .56299, .62670, .91819, .92845 /)

!       All output is to the terminal.  To modify this, change all
!       following "WRITE(*" to "WRITE(LUN" and precede this comment
!       by code to initialize unit LUN for output.

WRITE(*,101)
101  FORMAT('   M    CALCULATED     TRUE       DIFF.    IER')
NUMBER=0
k=0
DO  m=2,4
  DO  ih=1,2
    DO  ir=1,2
      DO  l=1,m
        h(l)=ih
        DO  j=1,m
          r(l,j)=0.25*ir
        END DO
        r(l,l)=1.
      END DO
      NUMBER=NUMBER+1
      CALL dmv(m, k, h, r, prob, 1.d-5, ier, ERR, 15.d0)
      diff=prob - result(NUMBER)
      WRITE(*, 34) m, prob, result(NUMBER), diff, ier
      34  FORMAT(i4, 2F12.5, ' ', g12.2, i5)
    END DO
  END DO
END DO
STOP
END PROGRAM test725
