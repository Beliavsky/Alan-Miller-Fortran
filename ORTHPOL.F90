MODULE orthpol
!      ALGORITHM 726, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 20, NO. 1, MARCH, 1994, PP. 21-62.
!
! This is a package of routines, called ORTHPOL, for generating orthogonal
! polynomials and Gauss-type quadrature rules developed by Walter Gautschi.
! A description of the underlying methods can be found in a companion paper
! published in ``ACM Transactions on Mathematical Software''.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:51
! Latest revision - 7 June 2001

IMPLICIT NONE

INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)


CONTAINS


SUBROUTINE drecur(n, ipoly, dal, dbe, da, db, iderr)
 
! This subroutine generates the coefficients  a(k),b(k), k=0,1,...,n-1,
! in the recurrence relation

!       p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
!                            k=0,1,...,n-1,

!       p(-1)(x)=0,  p(0)(x)=1,

! for some classical (monic) orthogonal polynomials, and sets  b(0)
! equal to the total mass of the weight distribution. The results are
! stored in the arrays  a,b,  which hold, respectively, the coefficients
! a(k-1),b(k-1), k=1,2,...,n.

!       Input:  n - - the number of recursion coefficients desired
!               ipoly-integer identifying the polynomial as follows:
!                     1=Legendre polynomial on (-1,1)
!                     2=Legendre polynomial on (0,1)
!                     3=Chebyshev polynomial of the first kind
!                     4=Chebyshev polynomial of the second kind
!                     5=Jacobi polynomial with parameters  al=-.5,be=.5
!                     6=Jacobi polynomial with parameters  al,be
!                     7=generalized Laguerre polynomial with
!                       parameter  al
!                     8=Hermite polynomial
!               al,be-input parameters for Jacobi and generalized
!                     Laguerre polynomials

!       Output: a,b - arrays containing, respectively, the recursion
!                     coefficients  a(k-1),b(k-1), k=1,2,...,n.
!               ierr -an error flag, equal to  0  on normal return,
!                     equal to  1  if  al  or  be  are out of range
!                     when  ipoly=6  or  ipoly=7, equal to  2  if  b(0)
!                     overflows when  ipoly=6  or  ipoly=7, equal to  3
!                     if  n  is out of range, and equal to  4  if  ipoly
!                     is not an admissible integer. In the case  ierr=2,
!                     the coefficient  b(0)  is set equal to the largest
!                     machine-representable number.

! The subroutine calls for the function subroutines  r1mach,gamma  and
! alga. The routines  gamma  and  alga, which are included in this file,
! evaluate respectively the gamma function and its logarithm for positive
! arguments.  They are used only in the cases  ipoly=6  and ipoly=7.

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: ipoly
REAL (dp), INTENT(IN)   :: dal
REAL (dp), INTENT(IN)   :: dbe
REAL (dp), INTENT(OUT)  :: da(:)
REAL (dp), INTENT(OUT)  :: db(:)
INTEGER, INTENT(OUT)    :: iderr

REAL (dp) :: dlmach, dkm1, dalpbe, dt, dal2, dbe2
INTEGER   :: k

IF(n < 1) THEN
  iderr = 3
  RETURN
END IF

dlmach = LOG(HUGE(0.0_dp))
iderr = 0
da(1:n) = 0.0_dp
IF(ipoly == 1) THEN
  db(1) = 2._dp
  IF (n == 1) RETURN
  DO  k=2,n
    dkm1 = k-1
    db(k) = 1._dp/(4._dp - 1._dp/(dkm1*dkm1))
  END DO
  RETURN
ELSE IF(ipoly == 2) THEN
  da(1) = .5_dp
  db(1) = 1._dp
  IF(n == 1) RETURN
  DO  k=2,n
    da(k) = .5_dp
    dkm1 = k-1
    db(k) = .25_dp/(4._dp - 1._dp/(dkm1*dkm1))
  END DO
  RETURN
ELSE IF(ipoly == 3) THEN
  db(1) = 4._dp*ATAN(1._dp)
  IF(n == 1) RETURN
  db(2) = .5_dp
  IF(n == 2) RETURN
  db(3:n) = .25_dp
  RETURN
ELSE IF(ipoly == 4) THEN
  db(1) = 2._dp*ATAN(1._dp)
  IF(n == 1) RETURN
  db(2:n) = .25_dp
  RETURN
ELSE IF(ipoly == 5) THEN
  db(1) = 4._dp*ATAN(1._dp)
  da(1) = .5_dp
  IF(n == 1) RETURN
  db(2:n) = .25_dp
  RETURN
ELSE IF(ipoly == 6) THEN
  IF(dal <= -1._dp .OR. dbe <= -1._dp) THEN
    iderr = 1
    RETURN
  ELSE
    dalpbe = dal + dbe
    da(1) = (dbe - dal)/(dalpbe + 2._dp)
    dt = (dalpbe+1._dp)*LOG(2._dp) + lngamma(dal+1._dp) + lngamma(dbe+1._dp) -  &
          lngamma(dalpbe + 2._dp)
    IF(dt > dlmach) THEN
      iderr = 2
      db(1) = HUGE(0.0_dp)
    ELSE
      db(1) = EXP(dt)
    END IF
    IF(n == 1) RETURN
    dal2 = dal*dal
    dbe2 = dbe*dbe
    da(2) = (dbe2-dal2)/((dalpbe + 2._dp)*(dalpbe + 4._dp))
    db(2) = 4._dp*(dal + 1._dp)*(dbe + 1._dp) /  &
            ((dalpbe + 3._dp)*(dalpbe + 2._dp)**2)
    IF(n == 2) RETURN
    DO  k=3,n
      dkm1 = k-1
      da(k) = .25_dp*(dbe2 - dal2)/(dkm1*dkm1*(1._dp + .5_dp*dalpbe/dkm1)  &
              *(1._dp + .5_dp*(dalpbe + 2._dp)/dkm1))
      db(k) = .25_dp*(1._dp + dal/dkm1)*(1._dp + dbe/dkm1)*(1._dp + dalpbe/  &
              dkm1)/((1._dp + .5_dp*(dalpbe + 1._dp)/dkm1)*(1._dp + .5_dp*(dalpbe  &
              -1._dp)/dkm1)*(1._dp + .5_dp*dalpbe/dkm1)**2)
    END DO
    RETURN
  END IF
ELSE IF(ipoly == 7) THEN
  IF(dal <= -1._dp) THEN
    iderr = 1
    RETURN
  ELSE
    da(1) = dal + 1._dp
    db(1) = dgamma(dal + 1._dp)
    IF(iderr == 2) db(1) = HUGE(0.0_dp)
    IF(n == 1) RETURN
    DO  k=2,n
      dkm1 = k-1
      da(k) = 2._dp*dkm1 + dal + 1._dp
      db(k) = dkm1*(dkm1 + dal)
    END DO
    RETURN
  END IF
ELSE IF(ipoly == 8) THEN
  db(1) = SQRT(4._dp*ATAN(1._dp))
  IF(n == 1) RETURN
  DO  k=2,n
    db(k) = .5_dp*(k-1)
  END DO
ELSE
  iderr = 4
END IF

RETURN
END SUBROUTINE drecur



FUNCTION lngamma(z) RESULT(lanczos)

!  Uses Lanczos-type approximation to ln(gamma) for z > 0.
!  Reference:
!       Lanczos, C. 'A precision approximation of the gamma
!               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
!  Accuracy: About 14 significant digits except for small regions
!            in the vicinity of 1 and 2.

!  Programmer: Alan Miller
!              1 Creswick Street, Brighton, Vic. 3187, Australia
!  Latest revision - 14 October 1996

REAL(dp), INTENT(IN) :: z
REAL(dp)             :: lanczos

! Local variables

REAL(dp)  :: a(9) = (/ 0.9999999999995183D0, 676.5203681218835D0, &
                      -1259.139216722289D0, 771.3234287757674D0, &
                      -176.6150291498386D0, 12.50734324009056D0, &
                      -0.1385710331296526D0, 0.9934937113930748D-05, &
                       0.1659470187408462D-06 /), zero = 0.D0,   &
                       one = 1.d0, lnsqrt2pi = 0.9189385332046727D0, &
                       half = 0.5d0, sixpt5 = 6.5d0, seven = 7.d0, tmp
INTEGER  :: j

IF (z <= zero) THEN
  WRITE(*, *) 'Error: zero or -ve argument for lngamma'
  RETURN
END IF

lanczos = zero
tmp = z + seven
DO j = 9, 2, -1
  lanczos = lanczos + a(j)/tmp
  tmp = tmp - one
END DO
lanczos = lanczos + a(1)
lanczos = LOG(lanczos) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
RETURN

END FUNCTION lngamma



FUNCTION dgamma(dx) RESULT(fn_val)

! This evaluates the gamma function for real positive  dx, using the
! function  lngamma.

REAL (dp), INTENT(IN)  :: dx
REAL (dp)              :: fn_val

REAL (dp) :: dlmach, dt

dlmach = LOG(HUGE(0.0_dp))
dt = lngamma(dx)
IF(dt >= dlmach) THEN
  fn_val = HUGE(0.0_dp)
  WRITE(*, *) '*** Overflow in calculating gamma function ***'
ELSE
  fn_val = EXP(dt)
END IF

RETURN
END FUNCTION dgamma



SUBROUTINE dcheb(n, da, db, dnu, dalpha, dbeta, ds, iderr)

! N.B. Arguments DS0, DS1 & DS2 have been removed.
 
! Given a set of polynomials  p(0),p(1),...,p(2*n-1)  satisfying

!        p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
!                        k=0,1,...,2*n-2,

!        p(-1)(x)=0,  p(0)(x)=1,

! and associated modified moments

!           fnu(k)=integral of p(k)(x)*dlambda(x),
!                        k=0,1,...,2*n-1,

! this subroutine uses the modified Chebyshev algorithm (see, e.g.,
! Section 2.4 of W. Gautschi,``On generating orthogonal polynomials'',
! SIAM J. Sci. Statist. Comput. 3, 1982, 289-317) to generate the
! recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the
! polynomials  pi(k)  orthogonal with respect to the integration
! measure  dlambda(x), i.e.,

!        pi(k+1)(x)=(x-alpha(k))*pi(k)(x)-beta(k)*pi(k-1)(x),
!                          k=0,1,...,n-1,

!        pi(-1)(x)=0,  pi(0)(x)=1.

!     Input:    n - - the number of recursion coefficients desired
!               a,b-- arrays of dimension 2*n-1 to be filled with the
!                     values of  a(k-1),b(k-1), k=1,2,...,2*n-1
!               fnu-- array of dimension  2*n  to be filled with the values
!                     of the modified moments  fnu(k-1), k=1,2, ...,2*n
!     Output:   alpha,beta-- arrays containing, respectively, the
!                     recursion coefficients  alpha(k-1),beta(k-1),
!                     k=1,2,...,n, where  beta(0)  is the total mass.
!               s - - array containing the normalization factors
!                     s(k)=integral [pi(k)(x)]**2 dlambda(x), k=0,1,
!                     2,...,n-1.
!               ierr- an error flag, equal to  0  on normal return,
!                     equal to  1  if  abs(fnu(0))  is less than the machine
!                     zero, equal to  2  if  n  is out of range, equal to  -k
!                     if  s(k), k=0,1,2,...,n-1, is about to underflow,
!                     and equal to  +k  if it is about to overflow.

! The arrays  ds0,ds1,ds2  are needed for working space.

! On machines with limited exponent range, the occurrence of underflow
! [overflow] in the computation of the  alpha's  and  beta's  can often
! be avoided by multiplying all modified moments by a sufficiently large
! [small] scaling factor and dividing the new  beta(0)  by the same
! scaling factor.

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: da(:)
REAL (dp), INTENT(IN)   :: db(:)
REAL (dp), INTENT(IN)   :: dnu(:)
REAL (dp), INTENT(OUT)  :: dalpha(:)
REAL (dp), INTENT(OUT)  :: dbeta(:)
REAL (dp), INTENT(OUT)  :: ds(:)
INTEGER, INTENT(OUT)    :: iderr

REAL (dp) :: dtiny, dhuge, ds0(2*n), ds1(2*n), ds2(2*n)
INTEGER   :: k, l, lk, nd

! The arrays  da, db  are assumed to have dimension  2*n-1, the arrays
! dnu, ds0, ds1, ds2  dimension  2*n.

nd = 2*n
dtiny = 10._dp*TINY(0.0_dp)
dhuge = .1_dp*HUGE(0.0_dp)
iderr = 0
IF(ABS(dnu(1)) < dtiny) THEN
  iderr = 1
  RETURN
END IF
IF(n < 1) THEN
  iderr = 2
  RETURN
END IF

dalpha(1) = da(1) + dnu(2)/dnu(1)
dbeta(1) = dnu(1)
IF(n == 1) RETURN
ds(1) = dnu(1)
DO  l=1,nd
  ds0(l) = 0._dp
  ds1(l) = dnu(l)
END DO
DO  k=2,n
  lk = nd - k + 1
  DO  l=k,lk
    ds2(l) = ds1(l+1) - (dalpha(k-1)-da(l))*ds1(l) - dbeta(k-1)*ds0(l)  &
             + db(l)*ds1(l-1)
    IF(l == k) ds(k) = ds2(k)
  END DO
  IF(ABS(ds(k)) < dtiny) THEN
    iderr = -(k-1)
    RETURN
  ELSE IF(ABS(ds(k)) > dhuge) THEN
    iderr = k-1
    RETURN
  END IF
  dalpha(k) = da(k) + (ds2(k+1)/ds2(k)) - (ds1(k)/ds1(k-1))
  dbeta(k) = ds2(k)/ds1(k-1)
  DO  l=k,lk
    ds0(l) = ds1(l)
    ds1(l) = ds2(l)
  END DO
END DO

RETURN
END SUBROUTINE dcheb



SUBROUTINE dsti(n, ncap, dx, dw, dalpha, dbeta, ierr)

! N.B. Arguments DP0, DP1 & DP2 have been removed.
 
! This routine applies ``Stieltjes's procedure'' (cf. Section 2.1 of
! W. Gautschi,``On generating orthogonal polynomials'', SIAM J. Sci.
! Statist. Comput. 3, 1982, 289-317) to generate the recursion
! coefficients  alpha(k), beta(k) , k=0,1,...,n-1, for the discrete
! (monic) orthogonal polynomials associated with the inner product

!     (f,g)=sum over k from 1 to ncap of w(k)*f(x(k))*g(x(k)).

! The integer  n  must be between  1  and  ncap, inclusive; otherwise,
! there is an error exit with  ierr=1. The results are stored in the
! arrays  alpha, beta; the arrays  p0, p1, p2  are working arrays.

! If there is a threat of underflow or overflow in the calculation
! of the coefficients  alpha(k)  and  beta(k), the routine exits with
! the error flag  ierr  set equal to  -k  (in the case of underflow)
! or  +k  (in the case of overflow), where  k  is the recursion index
! for which the problem occurs. The former [latter] can often be avoided
! by multiplying all weights  w(k)  by a sufficiently large [small]
! scaling factor prior to entering the routine, and, upon exit, divide
! the coefficient  beta(0)  by the same factor.

! This routine should be used with caution if  n  is relatively close
! to  ncap, since there is a distinct possibility of numerical
! instability developing. (See W. Gautschi,``Is the recurrence relation
! for orthogonal polynomials always stable?'', BIT, 1993, to appear.)
! In that case, the routine  lancz  should be used.

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: ncap
REAL (dp), INTENT(IN)   :: dx(:)
REAL (dp), INTENT(IN)   :: dw(:)
REAL (dp), INTENT(OUT)  :: dalpha(:)
REAL (dp), INTENT(OUT)  :: dbeta(:)
INTEGER, INTENT(OUT)    :: ierr

REAL (dp) :: dp0(ncap), dp1(ncap), dp2(ncap)
REAL (dp) :: dtiny, dhuge, dsum0, dsum1, dsum2, dt
INTEGER   :: k, m, nm1

dtiny = 10._dp*TINY(0.0_dp)
dhuge = .1_dp*HUGE(0.0_dp)
ierr = 0
IF(n <= 0 .OR. n > ncap) THEN
  ierr = 1
  RETURN
END IF

nm1 = n-1
dsum0 = 0._dp
dsum1 = 0._dp
DO  m=1,ncap
  dsum0 = dsum0 + dw(m)
  dsum1 = dsum1 + dw(m)*dx(m)
END DO
dalpha(1) = dsum1/dsum0
dbeta(1) = dsum0
IF(n == 1) RETURN
DO  m=1,ncap
  dp1(m) = 0._dp
  dp2(m) = 1._dp
END DO
DO  k=1,nm1
  dsum1 = 0._dp
  dsum2 = 0._dp
  DO  m=1,ncap
    IF(dw(m) == 0._dp) CYCLE
    dp0(m) = dp1(m)
    dp1(m) = dp2(m)
    dp2(m) = (dx(m) - dalpha(k))*dp1(m) - dbeta(k)*dp0(m)
    IF(ABS(dp2(m)) > dhuge .OR. ABS(dsum2) > dhuge) THEN
      ierr = k
      RETURN
    END IF
    dt = dw(m)*dp2(m)*dp2(m)
    dsum1 = dsum1 + dt
    dsum2 = dsum2 + dt*dx(m)
  END DO
  IF(ABS(dsum1) < dtiny) THEN
    ierr = -k
    RETURN
  END IF
  dalpha(k+1) = dsum2/dsum1
  dbeta(k+1) = dsum1/dsum0
  dsum0 = dsum1
END DO

RETURN
END SUBROUTINE dsti



SUBROUTINE dlancz(n, ncap, dx, dw, dalpha, dbeta, ierr)

! N.B. Arguments DP0 & DP1 have been removed.
 
! This routine carries out the same task as the routine  sti, but
! uses the more stable Lanczos method. The meaning of the input
! and output parameters is the same as in the routine  sti. (This
! routine is adapted from the routine RKPW in W.B. Gragg and
! W.J. Harrod,``The numerically stable reconstruction of Jacobi
! matrices from spectral data'', Numer. Math. 44, 1984, 317-335.)

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: ncap
REAL (dp), INTENT(IN)   :: dx(:)
REAL (dp), INTENT(IN)   :: dw(:)
REAL (dp), INTENT(OUT)  :: dalpha(:)
REAL (dp), INTENT(OUT)  :: dbeta(:)
INTEGER, INTENT(OUT)    :: ierr

REAL (dp) :: dp0(ncap), dp1(ncap)
REAL (dp) :: dpi, dgam, dsig, dt, dxlam, drho, dtmp, dtsig, dtk
INTEGER   :: i, k

IF(n <= 0 .OR. n > ncap) THEN
  ierr = 1
  RETURN
ELSE
  ierr = 0
END IF
DO  i=1,ncap
  dp0(i) = dx(i)
  dp1(i) = 0._dp
END DO
dp1(1) = dw(1)
DO  i=1,ncap-1
  dpi = dw(i+1)
  dgam = 1._dp
  dsig = 0._dp
  dt = 0._dp
  dxlam = dx(i+1)
  DO  k=1,i+1
    drho = dp1(k) + dpi
    dtmp = dgam*drho
    dtsig = dsig
    IF(drho <= 0._dp) THEN
      dgam = 1._dp
      dsig = 0._dp
    ELSE
      dgam = dp1(k)/drho
      dsig = dpi/drho
    END IF
    dtk = dsig*(dp0(k)-dxlam) - dgam*dt
    dp0(k) = dp0(k) - (dtk-dt)
    dt = dtk
    IF(dsig <= 0._dp) THEN
      dpi = dtsig*dp1(k)
    ELSE
      dpi = (dt**2)/dsig
    END IF
    dtsig = dsig
    dp1(k) = dtmp
  END DO
END DO
DO  k=1,n
  dalpha(k) = dp0(k)
  dbeta(k) = dp1(k)
END DO

RETURN
END SUBROUTINE dlancz



SUBROUTINE dmcdis(n, ncapm, mc, mp, dxp, dyp, dquad, deps, iq, idelta,  &
                  irout, finld, finrd, dendl, dendr, dalpha, &
                  dbeta, ncap, kount, ierrd, ied, dbe, dx, dw, dwf)

! N.B. The name dwf of the user's function has been added as an extra argument.
! N.B. Arguments DXFER, DWFER, DXM, DWM, DP0, DP1 & DP2 have been removed.

! This is a multiple-component discretization procedure as described in
! Section 4.3 of the companion paper. It generates to a relative
! accuracy of  eps  the recursion coefficients  alpha(k), beta(k),
! k=0,1,...,n-1, for the polynomials orthogonal with respect to a
! weight distribution consisting of the sum of  mc  continuous
! components and a discrete component with  mp  points. The continuous
! part of the spectrum is made up of  mc  weight functions, each
! supported on its own interval. These intervals may or may not be
! disjoint. The discretization of the inner product on the i-th
! interval is furnished either by a user-supplied subroutine  quad,
! or by the general-purpose subroutine  qgp  provided in this package,
! depending on whether  iq  is equal, or not equal, to  1, respectively.
! The user-supplied routine must have the form  quad(n,x,w,i,ierr)  and
! is assumed to supply the abscissas  x(k)  and weights  w(k), k=1,2,
! ...,n, to be used in approximating the i-th inner product

!               integral of p(x)*q(x)*wf(x,i)dx

! by the

!       sum over k from 1 to n of w(k)*p(x(k))*q(x(k)),

!                                        i=1,2,...,mc.

! The desired recurrence coefficients are then approximated by the
! recursion coefficients of the discrete orthogonal polynomials
! belonging to the discretized inner product, which in turn are
! computed by either the Stieltjes procedure or the Lanczos algorithm
! according as  irout  is equal to, or not equal to  1, respectively.
! Two error flags  ierr,ie  are provided which signal the occurrence
! of an error condition in the quadrature process, or in the routine
! sti  or  lancz  (whichever is used), respectively. The point spectrum
! is given through its abscissas  xp  and jumps  yp.

! If the quadrature routine  quad  has polynomial degree of exactness
! at least  id(n)  for each i, and if  id(n)/n = idelta + O(1/n)  as
! n  goes to infinity, then the procedure is designed to converge after
! one iteration, provided  idelta  is set with the appropriate
! integer. Normally,  idelta=1 (for interpolatory rules) or  idelta=2
! (for Gaussian rules). The default value is  idelta=1.

!    Input:  n    - - the number of recursion coefficients desired;
!                     type integer
!            ncapm  - a discretization parameter indicating an upper
!                     limit of the fineness of the discretization;
!                     ncapm=500  will usually be satisfactory; type
!                     integer
!            mc  - -  the number of disjoint intervals in the
!                     continuous part of the spectrum; type integer
!            mp  - -  the number of points in the discrete part of
!                     the spectrum; type integer. If there is no
!                     point spectrum, set  mp=0.
!            xp  - -  an array of dimension  mp  containing the
!                     abscissas of the point spectrum
!            yp  - -  an array of dimension  mp  containing the jumps
!                     of the point spectrum
!            quad  -  a subroutine determining the discretization of
!                     the inner product on each component interval,
!                     or a dummy routine if  iq  is not equal to  1
!                     (see below)
!            eps  - - the desired relative accuracy of the nonzero
!                     recursion coefficients; type real
!            iq   - - an integer selecting a user-supplied quadrature
!                     routine  quad  if  iq=1  or the ORTHPOL routine
!                     qgp  otherwise
!            idelta - a nonzero integer, typically  1  or  2, inducing
!                     fast convergence in the case of special quadrature
!                     routines
!            irout  - an integer selecting the routine for generating
!                     the recursion coefficients from the discrete
!                     inner product. Specifically,  irout=1  selects the
!                     routine  sti, whereas any other value selects the
!                     routine  lancz

! The logical variables  finl,finr, the arrays  endl,endr  of
! dimension  mc, and the arrays  xfer,wfer  of dimension  ncapm  are
! input variables to the subroutine  qgp  and are used (and hence need
! to be properly dimensioned) only if  iq  is not equal to  1.

!    Output:  alpha,beta - arrays of dimension n, holding as k-th
!                     element  alpha(k-1), beta(k-1), k=1,2,...,n,
!                     respectively
!             ncap  - an integer indicating the fineness of the
!                     discretization that yields convergence within
!                     the eps-tolerance
!             kount - the number of iterations used
!             ierr  - an error flag, equal to  0  on normal return,
!                     equal to  -1  if  n  is not in the proper range,
!                     equal to  i  if there is an error condition in
!                     the discretization of the i-th interval,
!                     and equal to  ncapm  if the discretized
!                     Stieltjes procedure does not converge within the
!                     discretization resolution specified by  ncapm
!             ie - -  an error flag inherited from the routine  sti
!                     or  lancz  (whichever is used)

! The array  be  of dimension  n, the arrays  x,w  of dimension  ncapm,
! and the arrays  xm,wm,p0,p1,p2  of dimension mc*ncapm + mp  are used
! for working space.  The routine calls upon the subroutine  sti  or
! lancz, depending on the choice of  irout.

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: ncapm
INTEGER, INTENT(IN)      :: mc
INTEGER, INTENT(IN)      :: mp
REAL (dp), INTENT(IN)    :: dxp(:)
REAL (dp), INTENT(IN)    :: dyp(:)
REAL (dp), INTENT(IN)    :: deps
INTEGER, INTENT(IN)      :: iq
INTEGER, INTENT(IN OUT)  :: idelta
INTEGER, INTENT(IN)      :: irout
LOGICAL, INTENT(IN)      :: finld
LOGICAL, INTENT(IN)      :: finrd
REAL (dp), INTENT(IN)    :: dendl(:)
REAL (dp), INTENT(IN)    :: dendr(:)
REAL (dp), INTENT(OUT)   :: dalpha(:)
REAL (dp), INTENT(OUT)   :: dbeta(:)
INTEGER, INTENT(OUT)     :: ncap
INTEGER, INTENT(OUT)     :: kount
INTEGER, INTENT(OUT)     :: ierrd
INTEGER, INTENT(OUT)     :: ied
REAL (dp), INTENT(OUT)   :: dbe(:)
REAL (dp), INTENT(OUT)   :: dx(:)
REAL (dp), INTENT(OUT)   :: dw(:)

INTERFACE
  SUBROUTINE dquad(n, dx, dw, i, ierr)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN OUT)  :: dx(:), dw(:)
    INTEGER, INTENT(IN)        :: i
    INTEGER, INTENT(OUT)       :: ierr
  END SUBROUTINE dquad

  FUNCTION dwf(dx, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: dx
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION dwf
END INTERFACE

REAL (dp)  :: dxm(mc*ncapm+mp), dwm(mc*ncapm+mp)
INTEGER    :: i, ierr, im1tn, incap, k, mtncap

! The arrays  dxp,dyp  are assumed to have dimension  mp  if mp > 0,

IF(idelta <= 0) idelta = 1
IF(n < 1) THEN
  ierrd = -1
  RETURN
END IF

incap = 1
kount = -1
ierrd = 0
dbeta(1:n) = 0._dp
ncap = (2*n-1) / idelta
20 dbe(1:n) = dbeta(1:n)
kount = kount + 1
IF(kount > 1) incap = 2**(kount/5)*n
ncap = ncap + incap
IF(ncap > ncapm) THEN
  ierrd = ncapm
  RETURN
END IF
mtncap = mc*ncap
DO  i=1,mc
  im1tn = (i-1)*ncap
  IF(iq == 1) THEN
    CALL dquad(ncap, dx, dw, i, ierr)
    IF(ierr /= 0) THEN
      ierrd = i
      RETURN
    END IF
  ELSE
    CALL dqgp(ncap, dx, dw, i, mc, finld, finrd, dendl, dendr, dwf)
  END IF

  DO  k=1,ncap
    dxm(im1tn+k) = dx(k)
    dwm(im1tn+k) = dw(k)
  END DO
END DO
IF(mp /= 0) THEN
  DO  k=1,mp
    dxm(mtncap+k) = dxp(k)
    dwm(mtncap+k) = dyp(k)
  END DO
END IF
IF(irout == 1) THEN
  CALL dsti(n, mtncap+mp, dxm, dwm, dalpha, dbeta, ied)
ELSE
  CALL dlancz(n, mtncap+mp, dxm, dwm, dalpha, dbeta, ied)
END IF
DO  k=1,n
  IF(ABS(dbeta(k)-dbe(k)) > deps*ABS(dbeta(k))) GO TO 20
END DO

RETURN
END SUBROUTINE dmcdis



SUBROUTINE dqgp(n, dx, dw, i, mcd, finld, finrd, dendl, dendr, dwf)

! N.B. Arguments DXFER, DWFER & IERR have been removed.
 
! N.B. The name dwf of the user's function has been added as an extra argument.

! This is a general-purpose discretization routine that can be used
! as an alternative to the routine  quad  in the multiple-component
! discretization procedure  mcdis.  It takes no account of the special
! nature of the weight function involved and hence may result in slow
! convergence of the discretization procedure.  This routine, therefore,
! should be used only as a last resort, when no better, more natural
! discretization can be found.

! It is assumed that there are  mc >= 1  disjoint component intervals.
! The discretization is effected by the Fejer quadrature rule,
! suitably transformed to the respective interval.  An interval that
! extends to minus infinity has to be indexed by  1; one that extends
! to plus infinity has to be indexed by  mc.

! The output variable  ierr  is given the value  0.  Additional input
! parameters and working space used by this routine are as follows:

!          mc      - the number of component intervals; type integer
!          finl    - a logical variable to be set .true. if the extreme
!                    left interval is finite, and .false. otherwise
!          finr    - a logical variable to be set .true. if the extreme
!                    right interval is finite, and .false. otherwise
!          endl    - an array of dimension  mc  containing the left
!                    endpoints of the component intervals; if the
!                    first of these extends to minus infinity,  endl(1)
!                    can be set to an arbitrary value
!          endr    - an array of dimension  mc  containing the right
!                    endpoints of the component intervals; if the
!                    last of these extends to plus infinity,  endr(mc)
!                    can be set to an arbitrary value
!          xfer,wfer-working arrays holding the Fejer nodes and
!                    weights, respectively, for the interval [-1,1].

! The user has to supply the routine

!              REAL (dp) function dwf(dx,i),

! which evaluates the weight function in REAL (dp) at the
! point  dx  on the i-th component interval.

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: dx(:)
REAL (dp), INTENT(OUT)  :: dw(:)
INTEGER, INTENT(IN)     :: i
INTEGER, INTENT(IN)     :: mcd
LOGICAL, INTENT(IN)     :: finld
LOGICAL, INTENT(IN)     :: finrd
REAL (dp), INTENT(IN)   :: dendl(:)
REAL (dp), INTENT(IN)   :: dendr(:)

INTERFACE
  FUNCTION dwf(dx, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: dx
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION dwf
END INTERFACE

REAL (dp) :: dphi, dphi1, dxfer(n), dwfer(n)
INTEGER   :: k

IF(i == 1) CALL dfejer(n, dxfer, dwfer)
IF(i > 1 .AND. i < mcd) GO TO 60
IF(mcd == 1) THEN
  IF(finld .AND. finrd) GO TO 60
  IF(finld) GO TO 20
  IF(finrd) GO TO 40
  DO  k=1,n
    CALL dsymtr(dxfer(k), dphi, dphi1)
    dx(k) = dphi
    dw(k) = dwfer(k)*dwf(dphi, i)*dphi1
  END DO
  RETURN
ELSE
  IF((i == 1 .AND. finld) .OR. (i == mcd .AND. finrd)) GO TO 60
  IF(i == 1) GO TO 40
END IF

20 DO  k=1,n
  CALL dtr(dxfer(k), dphi, dphi1)
  dx(k) = dendl(mcd) + dphi
  dw(k) = dwfer(k)*dwf(dx(k), mcd)*dphi1
END DO
RETURN

40 DO  k=1,n
  CALL dtr(-dxfer(k), dphi, dphi1)
  dx(k) = dendr(1) - dphi
  dw(k) = dwfer(k)*dwf(dx(k), 1)*dphi1
END DO
RETURN

60 DO  k=1,n
  dx(k) = .5_dp*((dendr(i)-dendl(i))*dxfer(k) + dendr(i) + dendl(i))
  dw(k) = .5_dp*(dendr(i)-dendl(i))*dwfer(k)*dwf(dx(k),i)
END DO

RETURN
END SUBROUTINE dqgp



SUBROUTINE dsymtr(dt, dphi, dphi1)

! This implements a particular transformation  x=phi(t)  mapping
! the t-interval [-1,1] to the x-interval [-oo,oo].

!        input:   t
!        output:  phi=phi(t)
!                 phi1=derivative of phi(t)

REAL (dp), INTENT(IN)   :: dt
REAL (dp), INTENT(OUT)  :: dphi
REAL (dp), INTENT(OUT)  :: dphi1

REAL (dp) :: dt2

dt2 = dt*dt
dphi = dt/(1.-dt2)
dphi1 = (dt2 + 1._dp)/(dt2 - 1._dp)**2

RETURN
END SUBROUTINE dsymtr



SUBROUTINE dtr(dt, dphi, dphi1)

! This implements a particular transformation  x=phi(t)  mapping
! the t-interval [-1,1] to the x-interval [0,oo].

!         input:   t
!         output:  phi=phi(t)
!                  phi1=derivative of phi(t)

REAL (dp), INTENT(IN)   :: dt
REAL (dp), INTENT(OUT)  :: dphi
REAL (dp), INTENT(OUT)  :: dphi1

dphi = (1._dp + dt)/(1._dp - dt)
dphi1 = 2._dp/(dt - 1._dp)**2

RETURN
END SUBROUTINE dtr



SUBROUTINE dfejer(n, dx, dw)

! This routine generates the n-point Fejer quadrature rule.

!         input:   n   - the number of quadrature nodes
!         output:  x,w - arrays of dimension  n  holding the quadrature
!                        nodes and weights, respectively; the nodes
!                        are ordered increasingly

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: dx(:)
REAL (dp), INTENT(OUT)  :: dw(:)

REAL (dp) :: dpi, dn, dc1, dc0, dt, dsum, dc2
INTEGER   :: k, m, nh, np1h

dpi = 4._dp*ATAN(1._dp)
nh = n/2
np1h = (n+1)/2
dn = n
DO  k=1,nh
  dx(n+1-k) = COS(.5_dp*(2*k-1)*dpi/dn)
  dx(k) = -dx(n+1-k)
END DO
IF(2*nh /= n) dx(np1h) = 0._dp
DO  k=1,np1h
  dc1 = 1._dp
  dc0 = 2._dp*dx(k)*dx(k) - 1._dp
  dt = 2._dp*dc0
  dsum = dc0/3._dp
  DO  m=2,nh
    dc2 = dc1
    dc1 = dc0
    dc0 = dt*dc1 - dc2
    dsum = dsum + dc0/DBLE(4*m*m-1)
  END DO
  dw(k) = 2._dp*(1._dp-2._dp*dsum)/dn
  dw(n+1-k) = dw(k)
END DO

RETURN
END SUBROUTINE dfejer



SUBROUTINE dmcheb(n, ncapm, mcd, mp, dxp, dyp, dquad, deps, iq, idelta,  &
                  finld, finrd, dendl, dendr, dxfer, dwfer, da, db, dnu, &
                  dalpha, dbeta, ncap, kount, ierrd, dx, dw, ds, dwf)
 
! N.B. The name dwf of the user's function has been added as an extra argument.
! N.B. Arguments DBE, DXM, DWM, DS0, DS1 & DS2 have been removed.

! This is a multiple-component discretized modified Chebyshev
! algorithm, basically a modified Chebyshev algorithm in which the
! modified moments are discretized in the same manner as the inner
! product in the discretization procedure  mcdis.  The input and
! output parameters are as in  mcdis.  In addition, the arrays  a,b
! must be filled with the recursion coefficients  a(k-1),b(k-1),
! k=1,2,...,2*n-1, defining the modified moments.  The arrays
! be,x,w,xm,wm,s,s0,s1,s2  are used for working space.  The routine
! calls upon the subroutine  cheb.  The routine exits immediately with
! ierr=-1  if  n  is not in range.

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: ncapm
INTEGER, INTENT(IN)        :: mcd
INTEGER, INTENT(IN)        :: mp
REAL (dp), INTENT(IN)      :: dxp(:)
REAL (dp), INTENT(IN)      :: dyp(:)
REAL (dp), INTENT(IN)      :: deps
INTEGER, INTENT(IN OUT)    :: iq
INTEGER, INTENT(OUT)       :: idelta
LOGICAL, INTENT(IN OUT)    :: finld
LOGICAL, INTENT(IN OUT)    :: finrd
REAL (dp), INTENT(IN OUT)  :: dendl(:)
REAL (dp), INTENT(IN OUT)  :: dendr(:)
REAL (dp), INTENT(IN OUT)  :: dxfer(:)
REAL (dp), INTENT(IN OUT)  :: dwfer(:)
REAL (dp), INTENT(IN)      :: da(:)
REAL (dp), INTENT(IN)      :: db(:)
REAL (dp), INTENT(OUT)     :: dnu(:)
REAL (dp), INTENT(IN OUT)  :: dalpha(:)
REAL (dp), INTENT(OUT)     :: dbeta(:)
INTEGER, INTENT(OUT)       :: ncap
INTEGER, INTENT(OUT)       :: kount
INTEGER, INTENT(OUT)       :: ierrd
REAL (dp), INTENT(OUT)     :: dx(:)
REAL (dp), INTENT(OUT)     :: dw(:)
REAL (dp), INTENT(IN OUT)  :: ds(:)

INTERFACE
  SUBROUTINE dquad(n, dx, dw, i, ierr)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN OUT)  :: dx(:), dw(:)
    INTEGER, INTENT(IN)        :: i
    INTEGER, INTENT(OUT)       :: ierr
  END SUBROUTINE dquad

  FUNCTION dwf(dx, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: dx
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION dwf
END INTERFACE

REAL (dp) :: dbe(n), dxm(mcd*ncapm+mp), dwm(mcd*ncapm+mp)
REAL (dp) :: dsum, dp1, dpp, dpm1
INTEGER   :: i, ierr, incap, im1tn, k, km1, l, mtncap, mtnpmp, nd

! The arrays  dxp,dyp  are assumed to have dimension  mp  if mp > 0,
! the arrays  da,db  dimension 2*n-1, the arrays  dnu,ds0,ds1,ds2
! dimension  2*n, and the arrays  dxm,dwm  dimension  mc*ncapm+mp.

nd = 2*n
IF(idelta <= 0) idelta = 1
IF(n < 1) THEN
  ierrd = -1
  RETURN
END IF

incap = 1
kount = -1
ierrd = 0
dbeta(1:n) = 0._dp
ncap = (nd-1)/idelta
20 dbe(1:n) = dbeta(1:n)
kount = kount + 1
IF(kount > 1) incap = 2**(kount/5)*n
ncap = ncap + incap
IF(ncap > ncapm) THEN
  ierrd = ncapm
  RETURN
END IF

mtncap = mcd*ncap
DO  i=1,mcd
  im1tn = (i-1)*ncap
  IF(iq == 1) THEN
    CALL dquad(ncap, dx, dw, i, ierr)
    IF(ierr /= 0) THEN
      ierrd = i
      RETURN
    END IF
  ELSE
    CALL dqgp(ncap, dx, dw, i, mcd, finld, finrd, dendl, dendr, dwf)
  END IF

  DO  k=1,ncap
    dxm(im1tn+k) = dx(k)
    dwm(im1tn+k) = dw(k)
  END DO
END DO
IF(mp /= 0) THEN
  DO  k=1,mp
    dxm(mtncap+k) = dxp(k)
    dwm(mtncap+k) = dyp(k)
  END DO
END IF
mtnpmp = mtncap + mp
DO  k=1,nd
  km1 = k-1
  dsum = 0._dp
  DO  i=1,mtnpmp
    dp1 = 0._dp
    dpp = 1._dp
    IF(k > 1) THEN
      DO  l=1,km1
        dpm1 = dp1
        dp1 = dpp
        dpp = (dxm(i)-da(l))*dp1 - db(l)*dpm1
      END DO
    END IF
    dsum = dsum + dwm(i)*dpp
  END DO
  dnu(k) = dsum
END DO

CALL dcheb(n, da, db, dnu, dalpha, dbeta, ds, ierr)
DO  k=1,n
  IF(ABS(dbeta(k)-dbe(k)) > deps*ABS(dbeta(k))) GO TO 20
END DO

RETURN
END SUBROUTINE dmcheb



SUBROUTINE dchri(n, iopt, da, db, dx, dy, dhr, dhi, dalpha, dbeta, ierr)
 
! This subroutine implements the Christoffel or generalized Christoffel
! theorem. In all cases except  iopt=7, it uses nonlinear recurrence
! algorithms described in W. Gautschi,``An algorithmic implementation
! of the generalized Christoffel theorem'', Numerical Integration
! (G. Haemmerlin, ed.), Birkhaeuser, Basel, 1982, pp. 89-106. The case
! iopt=7  incorporates a QR step with shift  x  in the manner of
! J. Kautsky and G.H. Golub, ``On the calculation of Jacobi matrices'',
! Linear Algebra Appl. 52/53, 1983, 439-455, using the algorithm of
! Eq. (67.11) on p. 567 in J.H. Wilkinson,``The Algebraic Eigenvalue
! Problem'', Clarendon Press, Oxford, 1965. Given the recursion
! coefficients  a(k),b(k), k=0,1,...,n, for the (monic) orthogonal
! polynomials with respect to some measure  dlambda(t), it generates the
! recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the measure

!              (t-x)dlambda(t)               if  iopt=1
!              [(t-x)**2+y**2]dlambda(t)     if  iopt=2
!              (t**2+y**2)dlambda(t) with    if  iopt=3
!                dlambda(t) and supp(dlambda)
!                symmetric  with respect to
!                the origin
!              dlambda(t)/(t-x)              if  iopt=4
!              dlambda(t)/[(t-x)**2+y**2]    if  iopt=5
!              dlambda(t)/(t**2+y**2) with   if  iopt=6
!                dlambda(t) and supp(dlambda)
!                symmetric with respect to
!                the origin
!              [(t-x)**2]dlambda(t)          if  iopt=7


!      Input:   n  - - - the number of recurrence coefficients
!                        desired; type integer
!               iopt - - an integer selecting the desired weight distribution
!               a,b  - - arrays of dimension  n+1  containing the
!                        recursion coefficients a(k-1),b(k-1),k=1,2,
!                        ...,n+1, of the polynomials orthogonal with
!                        respect to the given measure  dlambda(t)
!               x,y  - - real parameters defining the linear and
!                        quadratic factors, or divisors, of  dlambda(t)
!               hr,hi  - the real and imaginary part, respectively, of
!                        the integral of dlambda(t)/(z-t), where z=x+iy;
!                        the parameter  hr  is used only if  iopt=4 or
!                        5, the parameter  hi  only if  iopt=5 or 6

!      Output:  alpha,beta - - arrays of dimension  n  containing the
!                         desired recursion coefficients  alpha(k-1),
!                         beta(k-1), k=1,2,...,n

! It is assumed that  n  is larger than or equal to 2. Otherwise, the
! routine exits immediately with the error flag  ierr  set equal to 1.
! If  iopt  is not between 1 and 7, the routine exits with  ierr=2.

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: iopt
REAL (dp), INTENT(IN)   :: da(:)
REAL (dp), INTENT(IN)   :: db(:)
REAL (dp), INTENT(IN)   :: dx
REAL (dp), INTENT(IN)   :: dy
REAL (dp), INTENT(IN)   :: dhr
REAL (dp), INTENT(IN)   :: dhi
REAL (dp), INTENT(OUT)  :: dalpha(:)
REAL (dp), INTENT(OUT)  :: dbeta(:)
INTEGER, INTENT(OUT)    :: ierr

REAL (dp) :: deps, de, dq, ds, dt, deio, dd, der, dei, deroo, deioo, dso,  &
             dero, deoo, deo, du, dc, dc0, dgam, dcm1, dp2
INTEGER   :: k, nm1

! The arrays  da, db  are assumed to have dimension  n+1.

deps = 5._dp*EPSILON(0.0_dp)
ierr = 0
IF(n < 2) THEN
  ierr = 1
  RETURN
END IF

IF(iopt == 1) THEN
  de = 0._dp
  DO  k=1,n
    dq = da(k) - de - dx
    dbeta(k) = dq*de
    de = db(k+1)/dq
    dalpha(k) = dx + dq + de
  END DO
  dbeta(1) = db(1)*(da(1) - dx)
ELSE IF(iopt == 2) THEN
  ds = dx - da(1)
  dt = dy
  deio = 0._dp
  DO  k=1,n
    dd = ds*ds + dt*dt
    der = -db(k+1)*ds/dd
    dei = db(k+1)*dt/dd
    ds = dx + der - da(k+1)
    dt = dy + dei
    dalpha(k) = dx + dt*der/dei - ds*dei/dt
    dbeta(k) = dt*deio*(1._dp + (der/dei)**2)
    deio = dei
  END DO
  dbeta(1) = db(1)*(db(2) + (da(1)-dx)**2 + dy*dy)
ELSE IF(iopt == 3) THEN
  dt = dy
  deio = 0._dp
  DO  k=1,n
    dei = db(k+1)/dt
    dt = dy + dei
    dalpha(k) = 0._dp
    dbeta(k) = dt*deio
    deio = dei
  END DO
  dbeta(1) = db(1)*(db(2) + dy*dy)
ELSE IF(iopt == 4) THEN
  dalpha(1) = dx - db(1)/dhr
  dbeta(1) = -dhr
  dq = -db(1)/dhr
  DO  k=2,n
    de = da(k-1) - dx - dq
    dbeta(k) = dq*de
    dq = db(k)/de
    dalpha(k) = dq + de + dx
  END DO
ELSE IF(iopt == 5) THEN
  nm1 = n-1
  dd = dhr*dhr + dhi*dhi
  deroo = da(1) - dx + db(1)*dhr/dd
  deioo = -db(1)*dhi/dd - dy
  dalpha(1) = dx + dhr*dy/dhi
  dbeta(1) = -dhi/dy
  dalpha(2) = dx - db(1)*dhi*deroo/(dd*deioo) + dhr*deioo/dhi
  dbeta(2) = dy*deioo*(1._dp + (dhr/dhi)**2)
  IF(n == 2) RETURN
  dso = db(2)/(deroo**2 + deioo**2)
  dero = da(2) - dx - dso*deroo
  deio = dso*deioo - dy
  dalpha(3) = dx + deroo*deio/deioo + dso*deioo*dero/deio
  dbeta(3) = -db(1)*dhi*deio*(1._dp + (deroo/deioo)**2)/dd
  IF(n == 3) RETURN
  DO  k=3,nm1
    ds = db(k)/(dero**2 + deio**2)
    der = da(k) - dx - ds*dero
    dei = ds*deio - dy
    dalpha(k+1) = dx + dero*dei/deio + ds*deio*der/dei
    dbeta(k+1) = dso*deioo*dei*(1._dp + (dero/deio)**2)
    deroo = dero
    deioo = deio
    dero = der
    deio = dei
    dso = ds
  END DO
ELSE IF(iopt == 6) THEN
  nm1 = n-1
  deoo = -db(1)/dhi-dy
  deo = db(2)/deoo-dy
  dalpha(1) = 0._dp
  dbeta(1) = -dhi/dy
  dalpha(2) = 0._dp
  dbeta(2) = dy*deoo
  IF(n == 2) RETURN
  dalpha(3) = 0._dp
  dbeta(3) = -db(1)*deo/dhi
  IF(n == 3) RETURN
  DO  k=3,nm1
    de = db(k)/deo-dy
    dbeta(k+1) = db(k-1)*de/deoo
    dalpha(k+1) = 0._dp
    deoo = deo
    deo = de
  END DO
ELSE IF(iopt == 7) THEN
  du = 0._dp
  dc = 1._dp
  dc0 = 0._dp
  DO  k=1,n
    dgam = da(k)-dx-du
    dcm1 = dc0
    dc0 = dc
    IF(ABS(dc0) > deps) THEN
      dp2 = (dgam**2)/dc0
    ELSE
      dp2 = dcm1*db(k)
    END IF
    IF(k > 1) dbeta(k) = ds*(dp2+db(k+1))
    ds = db(k+1)/(dp2+db(k+1))
    dc = dp2/(dp2+db(k+1))
    du = ds*(dgam+da(k+1)-dx)
    dalpha(k) = dgam+du+dx
  END DO
  dbeta(1) = db(1)*(db(2) + (dx-da(1))**2)
ELSE
  ierr=2
END IF

RETURN
END SUBROUTINE dchri



SUBROUTINE dknum(n, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, nu, ierr)

! N.B. Arguments DROLDR & DROLDI have been removed.
 
! This routine generates

!   rho(k)(z)=integral pi(k)(t)dlambda(t)/(z-t), k=0,1,2,...,n,

! where  pi(k)(t)  is the (monic) k-th degree orthogonal polynomial
! with respect to the measure  dlambda(t), and the integral is extended
! over the support of  dlambda. It is assumed that  z  is a complex
! number outside the smallest interval containing the support of
! dlambda. The quantities  rho(k)(z)  are computed as the first  n+1
! members of the minimal solution of the basic three-term recurrence
! relation

!      y(k+1)(z)=(z-a(k))y(k)(z)-b(k)y(k-1)(z), k=0,1,2,...,

! satisfied by the orthogonal polynomials  pi(k)(z).

!   Input:  n  - -  the largest integer  k  for which  rho(k)  is
!                   desired
!           nu0  -  an estimate of the starting backward recurrence
!                   index; if no better estimate is known, set
!                   nu0 = 3*n/2; for Jacobi, Laguerre and Hermite
!                   weight functions, estimates of  nu0  are generated
!                   respectively by the routines  nu0jac,nu0lag  and
!                   nu0her
!           numax - an integer larger than  n  cutting off backward
!                   recursion in case of nonconvergence; if  nu0
!                   exceeds  numax, then the routine aborts with the
!                   error flag  ierr  set equal to  nu0
!           z - - - the variable in  rho(k)(z); type complex
!           eps - - the relative accuracy to which the  rho(k)  are
!                   desired
!           a,b - - arrays of dimension  numax  to be supplied with the
!                   recurrence coefficients  a(k-1), b(k-1), k=1,2,...,
!                   numax.

!   Output: rho - - an array of dimension  n+1  containing the results
!                   rho(k)=rho(k-1)(z), k=1,2,...,n+1; type complex
!           nu  - - the starting backward recurrence index that yields
!                   convergence
!           ierr  - an error flag equal to zero on normal return, equal
!                   to  nu0  if  nu0 > numax, and equal to  numax in
!                   case of nonconvergence.

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN OUT)  :: nu0
INTEGER, INTENT(IN)      :: numax
REAL (dp), INTENT(IN)    :: dx
REAL (dp), INTENT(IN)    :: dy
REAL (dp), INTENT(IN)    :: deps
REAL (dp), INTENT(IN)    :: da(:)
REAL (dp), INTENT(IN)    :: db(:)
REAL (dp), INTENT(OUT)   :: drhor(:)
REAL (dp), INTENT(OUT)   :: drhoi(:)
INTEGER, INTENT(OUT)     :: nu
INTEGER, INTENT(OUT)     :: ierr

REAL (dp) :: drr, dri, dden, dt, droldr(n+1), droldi(n+1)
INTEGER   :: j, j1, k, np1

! The arrays  drhor,drhoi,droldr,droldi  are assumed to have dimension  n+1.

ierr = 0
np1 = n + 1
IF(nu0 > numax) THEN
  ierr = nu0
  RETURN
END IF

IF(nu0 < np1) nu0 = np1
nu = nu0 - 5
DO  k=1,np1
  drhor(k) = 0._dp
  drhoi(k) = 0._dp
END DO
20 nu = nu + 5
IF(nu > numax) THEN
  ierr = numax
  GO TO 60
END IF
DO  k=1,np1
  droldr(k) = drhor(k)
  droldi(k) = drhoi(k)
END DO
drr = 0._dp
dri = 0._dp
DO  j=1,nu
  j1=nu-j+1
  dden = (dx-da(j1)-drr)**2 + (dy-dri)**2
  drr = db(j1)*(dx-da(j1)-drr)/dden
  dri = -db(j1)*(dy-dri)/dden
  IF(j1 <= np1) THEN
    drhor(j1) = drr
    drhoi(j1) = dri
  END IF
END DO
DO  k=1,np1
  IF((drhor(k)-droldr(k))**2+(drhoi(k)-droldi(k))**2 >  &
      (deps**2)*(drhor(k)**2+drhoi(k)**2)) GO TO 20
END DO

60 IF(n == 0) RETURN
DO  k=2,np1
  dt = drhor(k)*drhor(k-1)-drhoi(k)*drhoi(k-1)
  drhoi(k) = drhor(k)*drhoi(k-1)+drhoi(k)*drhor(k-1)
  drhor(k) = dt
END DO

RETURN
END SUBROUTINE dknum



SUBROUTINE dkern(n,nu0,numax,dx,dy,deps,da,db,dkerr,dkeri,  &
                 nu,ierr,droldr,droldi)
 
! This routine generates the kernels in the Gauss quadrature remainder
! term, namely

!           K(k)(z)=rho(k)(z)/pi(k)(z), k=0,1,2,...,n,

! where  rho(k)  are the output quantities of the routine  knum, and
! pi(k)  the (monic) orthogonal polynomials.  The results are returned
! in the array  ker  as ker(k)=K(k-1)(z), k=1,2,...,n+1.  All the other
! input and output parameters have the same meaning as in the routine knum.

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN OUT)    :: nu0
INTEGER, INTENT(IN OUT)    :: numax
REAL (dp), INTENT(IN OUT)  :: dx
REAL (dp), INTENT(IN)      :: dy
REAL (dp), INTENT(IN OUT)  :: deps
REAL (dp), INTENT(IN)      :: da(:)
REAL (dp), INTENT(IN)      :: db(:)
REAL (dp), INTENT(OUT)     :: dkerr(:)
REAL (dp), INTENT(OUT)     :: dkeri(:)
INTEGER, INTENT(IN OUT)    :: nu
INTEGER, INTENT(IN OUT)    :: ierr
REAL (dp), INTENT(IN OUT)  :: droldr(:)
REAL (dp), INTENT(IN OUT)  :: droldi(:)

REAL (dp) :: dp0r, dp0i, dpr, dpi, dpm1r, dpm1i, dden, dt
INTEGER   :: k

! The arrays  dkerr,dkeri,droldr,droldi  are assumed to have dimension  n+1.

CALL dknum(n, nu0, numax, dx, dy, deps, da, db, dkerr, dkeri, nu, ierr)
IF(ierr /= 0) RETURN
dp0r=0._dp
dp0i=0._dp
dpr=1._dp
dpi=0._dp
DO  k=1,n
  dpm1r=dp0r
  dpm1i=dp0i
  dp0r=dpr
  dp0i=dpi
  dpr=(dx-da(k))*dp0r-dy*dp0i-db(k)*dpm1r
  dpi=(dx-da(k))*dp0i+dy*dp0r-db(k)*dpm1i
  dden=dpr**2+dpi**2
  dt=(dkerr(k+1)*dpr+dkeri(k+1)*dpi)/dden
  dkeri(k+1)=(dkeri(k+1)*dpr-dkerr(k+1)*dpi)/dden
  dkerr(k+1)=dt
END DO
RETURN
END SUBROUTINE dkern



SUBROUTINE dgchri(n, iopt, nu0, numax, deps, da, db, dx, dy, dalpha, dbeta,  &
                  nu, ierr, ierrc, dnu, drhor, drhoi, droldr, droldi)

! N.B. Arguments DS, DS0, DS1 & DS2 have been removed.
 
! This routine implements the generalized Christoffel theorem, using
! the method of modified moments (cf. Section 4 of W. Gautschi,
! ``Minimal solutions of three-term recurrence relations and orthogonal
! polynomials'', Math. Comp. 36, 1981, 547-554). Given the recursion
! coefficients  a(k), b(k), k=0,1,...n, for the (monic) orthogonal
! polynomials with respect to some measure  dlambda(t), it generates
! the recursion coefficients  alpha(k), beta(k), k=0,1,2,...,n-1 for
! the measure

!         dlambda(t)/(t-x)        if iopt=1
!         dlambda(t)/{(t-x)**2+y**2} if iopt=2

!   Input:  n   - - the number of recurrence coefficients desired;
!                   type integer
!           iopt  - an integer selecting the desired weight distribution
!           nu0   - an integer estimating the starting backward
!                   recurrence index; in the absence of any better
!                   choice, take  nu0 = 3*n
!           numax - an integer controlling termination of backward
!                   recursion in case of nonconvergence; a conservative
!                   choice is  numax = 500
!           eps - - a relative error tolerance; type real
!           a,b - - arrays of dimension numax to be supplied with the
!                   recursion coefficients a(k)=alpha(k-1),b(k)=beta(k),
!                   k=1,2,...,numax, for the measure  dlambda
!           x,y - - real parameters defining the linear and quadratic
!                   divisors of  dlambda

!   Output: alpha,beta - arrays of dimension  n  containing the desired
!                   recursion coefficients  alpha(k-1), beta(k-1), k=1,
!                   2,...,n
!           nu  - - the backward recurrence index yielding convergence;
!                   in case of nonconvergence,  nu  will have the value numax
!           ierr  - an error flag, where
!                   ierr=0     on normal return
!                   ierr=1     if  iopt  is neither 1 nor 2
!                   ierr=nu0   if  nu0 > numax
!                   ierr=numax if the backward recurrence algorithm does
!                              not converge
!                   ierr=-1    if  n  is not in range
!           ierrc - an error flag inherited from the routine  cheb

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: iopt
INTEGER, INTENT(IN OUT)    :: nu0
INTEGER, INTENT(IN)        :: numax
REAL (dp), INTENT(IN)      :: deps
REAL (dp), INTENT(IN)      :: da(:)
REAL (dp), INTENT(IN)      :: db(:)
REAL (dp), INTENT(IN OUT)  :: dx
REAL (dp), INTENT(IN OUT)  :: dy
REAL (dp), INTENT(OUT)     :: dalpha(:)
REAL (dp), INTENT(OUT)     :: dbeta(:)
INTEGER, INTENT(IN OUT)    :: nu
INTEGER, INTENT(OUT)       :: ierr
INTEGER, INTENT(IN OUT)    :: ierrc
REAL (dp), INTENT(OUT)     :: dnu(:)
REAL (dp), INTENT(OUT)     :: drhor(:)
REAL (dp), INTENT(OUT)     :: drhoi(:)
REAL (dp), INTENT(IN OUT)  :: droldr(:)
REAL (dp), INTENT(IN OUT)  :: droldi(:)

INTEGER    :: nd, ndm1
REAL (dp)  :: ds(n)

! The arrays  dnu,drhor,drhoi,droldr,droldi,ds0,ds1,ds2  are assumed
! to have dimension  2*n.

IF(n < 1) THEN
  ierr=-1
  RETURN
END IF
ierr = 0
nd = 2*n
ndm1 = nd - 1
IF(iopt == 1) THEN
  CALL dknum(ndm1, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, nu, ierr)
  dnu(1:nd) = -drhor(1:nd)
  CALL dcheb(n, da, db, dnu, dalpha, dbeta, ds, ierrc)
ELSE IF(iopt == 2) THEN
  dy=ABS(dy)
  CALL dknum(ndm1, nu0, numax, dx, dy, deps, da, db, drhor, drhoi, nu, ierr)
  dnu(1:nd) = -drhoi(1:nd)/dy
  CALL dcheb(n, da, db, dnu, dalpha, dbeta, ds, ierrc)
ELSE
  ierr=1
END IF

RETURN
END SUBROUTINE dgchri



SUBROUTINE dgauss(n, dalpha, dbeta, deps, dzero, dweigh, ierr)

! N.B. Argument DE has been removed.
 
! Given  n  and a measure  dlambda, this routine generates the n-point
! Gaussian quadrature formula

!     integral over supp(dlambda) of f(x)dlambda(x)

!        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).

! The nodes are returned as  zero(k)=x(k) and the weights as
! weight(k)=w(k), k=1,2,...,n.  The user has to supply the recursion
! coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure
! dlambda. The routine computes the nodes as eigenvalues, and the
! weights in term of the first component of the respective normalized
! eigenvectors of the n-th order Jacobi matrix associated with  dlambda.
! It uses a translation and adaptation of the algol procedure  imtql2,
! Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified
! by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for
! Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack
! routine  imtql2.

!        Input:  n - - the number of points in the Gaussian quadrature
!                      formula; type integer
!                alpha,beta - - arrays of dimension  n  to be filled with
!                      the values of  alpha(k-1), beta(k-1), k=1,2, ...,n
!                eps - the relative accuracy desired in the nodes and weights

!        Output: zero- array of dimension  n  containing the Gaussian nodes
!                      (in increasing order)  zero(k)=x(k), k=1,2, ...,n
!                weight - array of dimension  n  containing the
!                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
!                ierr- an error flag equal to  0  on normal return,
!                      equal to  i  if the QR algorithm does not
!                      converge within 30 iterations on evaluating the
!                      i-th eigenvalue, equal to  -1  if  n  is not in range,
!                      and equal to  -2  if one of the beta's is negative.

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: dalpha(:)
REAL (dp), INTENT(IN)   :: dbeta(:)
REAL (dp), INTENT(IN)   :: deps
REAL (dp), INTENT(OUT)  :: dzero(:)
REAL (dp), INTENT(OUT)  :: dweigh(:)
INTEGER, INTENT(OUT)    :: ierr

REAL (dp) :: dpp, dg, dr, ds, dc, df, db, de(n)
INTEGER   :: i, ii, j, k, l, m, mml

IF(n < 1) THEN
  ierr = -1
  RETURN
END IF
ierr = 0
dzero(1) = dalpha(1)
IF(dbeta(1) < 0._dp) THEN
  ierr = -2
  RETURN
END IF

dweigh(1) = dbeta(1)
IF (n == 1) RETURN
dweigh(1) = 1._dp
de(n) = 0._dp
DO  k=2,n
  dzero(k) = dalpha(k)
  IF(dbeta(k) < 0._dp) THEN
    ierr = -2
    RETURN
  END IF
  de(k-1) = SQRT(dbeta(k))
  dweigh(k) = 0._dp
END DO
DO  l=1,n
  j = 0
  105 DO  m=l,n
    IF(m == n) EXIT
    IF(ABS(de(m)) <= deps*(ABS(dzero(m)) + ABS(dzero(m+1)))) EXIT
  END DO
  dpp = dzero(l)
  IF(m == l) CYCLE
  IF(j == 30) GO TO 400
  j = j + 1
  dg = (dzero(l+1) - dpp) / (2._dp*de(l))
  dr = SQRT(dg*dg + 1._dp)
  dg = dzero(m) - dpp + de(l)/(dg + SIGN(dr,dg))
  ds = 1._dp
  dc = 1._dp
  dpp = 0._dp
  mml = m - l
  DO  ii=1,mml
    i = m-ii
    df = ds*de(i)
    db = dc*de(i)
    IF(ABS(df) < ABS(dg)) GO TO 150
    dc = dg/df
    dr = SQRT(dc*dc+1._dp)
    de(i+1) = df*dr
    ds = 1._dp/dr
    dc = dc*ds
    GO TO 160

    150 ds = df/dg
    dr = SQRT(ds*ds+1._dp)
    de(i+1) = dg*dr
    dc = 1._dp/dr
    ds = ds*dc

    160 dg = dzero(i+1)-dpp
    dr = (dzero(i)-dg)*ds+2._dp*dc*db
    dpp = ds*dr
    dzero(i+1) = dg+dpp
    dg = dc*dr-db
    df = dweigh(i+1)
    dweigh(i+1) = ds*dweigh(i)+dc*df
    dweigh(i) = dc*dweigh(i)-ds*df
  END DO
  dzero(l) = dzero(l)-dpp
  de(l) = dg
  de(m) = 0._dp
  GO TO 105
END DO
DO  ii=2,n
  i = ii-1
  k = i
  dpp = dzero(i)
  DO  j=ii,n
    IF(dzero(j) >= dpp) CYCLE
    k = j
    dpp = dzero(j)
  END DO
  IF(k == i) CYCLE
  dzero(k) = dzero(i)
  dzero(i) = dpp
  dpp = dweigh(i)
  dweigh(i) = dweigh(k)
  dweigh(k) = dpp
END DO
DO  k=1,n
  dweigh(k) = dbeta(1)*dweigh(k)*dweigh(k)
END DO
RETURN

400 ierr = l
RETURN
END SUBROUTINE dgauss



SUBROUTINE dradau(n, dalpha, dbeta, dend, dzero, dweigh, ierr)

! N.B. Arguments DE, DA & DB have been removed.
 
! Given  n  and a measure  dlambda, this routine generates the
! (n+1)-point Gauss-Radau quadrature formula

!   integral over supp(dlambda) of f(t)dlambda(t)

!     = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).

! The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
! =w(k), k=0,1,2,...,n. The user has to supply the recursion
! coefficients  alpha(k), beta(k), k=0,1,2,...,n, for the measure
! dlambda. The nodes and weights are computed as eigenvalues and
! in terms of the first component of the respective normalized
! eigenvectors of a slightly modified Jacobi matrix of order  n+1.
! To do this, the routine calls upon the subroutine  gauss.

!    Input:  n - -  the number of interior points in the Gauss-Radau
!                   formula; type integer
!            alpha,beta - arrays of dimension  n+1  to be supplied with
!                   the recursion coefficients  alpha(k-1), beta(k-1),
!                   k=1,2,...,n+1; the coefficient  alpha(n+1)  is not
!                   used by the routine
!            end -  the prescribed endpoint  x(0)  of the Gauss-Radau
!                   formula; type real

!    Output: zero - array of dimension  n+1  containing the nodes (in
!                   increasing order)  zero(k)=x(k), k=0,1,2,...,n
!            weight-array of dimension  n+1  containing the weights
!                   weight(k)=w(k), k=0,1,2,...,n
!            ierr - an error flag inherited from the routine  gauss

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: dalpha(:)
REAL (dp), INTENT(IN)      :: dbeta(:)
REAL (dp), INTENT(IN)      :: dend
REAL (dp), INTENT(IN OUT)  :: dzero(:)
REAL (dp), INTENT(IN OUT)  :: dweigh(:)
INTEGER, INTENT(IN OUT)    :: ierr

REAL (dp) :: depsma, dp0, dp1, dpm1, da(n+1), db(n+1)
INTEGER   :: k, np1

! The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
! dimension  n+1.

depsma = EPSILON(0.0_dp)

! depsma is the machine REAL (dp).

np1=n+1
DO  k=1,np1
  da(k)=dalpha(k)
  db(k)=dbeta(k)
END DO
dp0=0._dp
dp1=1._dp
DO  k=1,n
  dpm1=dp0
  dp0=dp1
  dp1=(dend-da(k))*dp0-db(k)*dpm1
END DO
da(np1) = dend - db(np1)*dp0/dp1
CALL dgauss(np1, da, db, depsma, dzero, dweigh, ierr)

RETURN
END SUBROUTINE dradau



SUBROUTINE dlob(n, dalpha, dbeta, dleft, dright, dzero, dweigh, ierr)

! N.B. Arguments DE, DA & DB have been removed.
 
! Given  n  and a measure  dlambda, this routine generates the
! (n+2)-point Gauss-Lobatto quadrature formula

!   integral over supp(dlambda) of f(x)dlambda(x)

!      = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k))

!              + w(n+1)f(x(n+1)) + R(n;f).

! The nodes are returned as  zero(k)=x(k), the weights as  weight(k)
! =w(k), k=0,1,...,n,n+1. The user has to supply the recursion
! coefficients  alpha(k), beta(k), k=0,1,...,n,n+1, for the measure
! dlambda. The nodes and weights are computed in terms of the
! eigenvalues and first component of the normalized eigenvectors of
! a slightly modified Jacobi matrix of order  n+2. The routine calls
! upon the subroutine  gauss  and the function subroutine  r1mach.

!   Input:  n - -  the number of interior points in the Gauss-Lobatto
!                  formula; type integer
!           alpha,beta - arrays of dimension  n+2  to be supplied with
!                  the recursion coefficients  alpha(k-1), beta(k-1),
!                  k=1,2,...,n+2, of the underlying measure; the
!                  routine does not use  alpha(n+2), beta(n+2)
!           aleft,right - the prescribed left and right endpoints
!                  x(0)  and  x(n+1)  of the Gauss-Lobatto formula

!   Output: zero - an array of dimension  n+2  containing the nodes (in
!                  increasing order)  zero(k)=x(k), k=0,1,...,n,n+1
!           weight-an array of dimension  n+2  containing the weights
!                  weight(k)=w(k), k=0,1,...,n,n+1
!           ierr - an error flag inherited from the routine  gauss

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: dalpha(:)
REAL (dp), INTENT(IN)      :: dbeta(:)
REAL (dp), INTENT(IN OUT)  :: dleft
REAL (dp), INTENT(IN)      :: dright
REAL (dp), INTENT(IN OUT)  :: dzero(:)
REAL (dp), INTENT(IN OUT)  :: dweigh(:)
INTEGER, INTENT(IN OUT)    :: ierr

REAL (dp) :: da(n+2), db(n+2)
REAL (dp) :: depsma, dp0l, dp0r, dp1l, dp1r, dpm1l, dpm1r, ddet
INTEGER   :: k, np1, np2

! The arrays  dalpha,dbeta,dzero,dweigh,de,da,db  are assumed to have
! dimension  n+2.

depsma = EPSILON(0.0_dp)

! depsma is the machine REAL (dp).

np1 = n+1
np2 = n+2
DO  k=1,np2
  da(k) = dalpha(k)
  db(k) = dbeta(k)
END DO
dp0l = 0._dp
dp0r = 0._dp
dp1l = 1._dp
dp1r = 1._dp
DO  k=1,np1
  dpm1l = dp0l
  dp0l = dp1l
  dpm1r = dp0r
  dp0r = dp1r
  dp1l = (dleft-da(k))*dp0l-db(k)*dpm1l
  dp1r = (dright-da(k))*dp0r-db(k)*dpm1r
END DO
ddet = dp1l*dp0r-dp1r*dp0l
da(np2) = (dleft*dp1l*dp0r-dright*dp1r*dp0l)/ddet
db(np2) = (dright-dleft)*dp1l*dp1r/ddet
CALL dgauss(np2,da,db,depsma,dzero,dweigh,ierr)

RETURN
END SUBROUTINE dlob



FUNCTION nu0jac(n, z, eps) RESULT(fn_val)
 
! This is an auxiliary function routine providing a starting backward
! recurrence index for the Jacobi measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.

INTEGER, INTENT(IN)  :: n
COMPLEX, INTENT(IN)  :: z
REAL, INTENT(IN)     :: eps
REAL                 :: fn_val

REAL  :: angle, pi, r, x, x2, y, y2

pi = 4.*ATAN(1.)
x = REAL(z)
y = ABS(AIMAG(z))
IF(x < 1.) THEN
  IF(x < -1.) angle = .5*(2.*pi + ATAN(y/(x-1.)) + ATAN(y/(x+1.)))
  IF(x == -1.) angle = .5*(1.5*pi - ATAN(.5*y))
  IF(x > -1.) angle = .5*(pi+ATAN(y/(x-1.)) + ATAN(y/(x+1.)))
ELSE
  IF(x == 1.) angle = .5*(.5*pi+ATAN(.5*y))
  IF(x > 1.) angle = .5*(ATAN(y/(x-1.)) + ATAN(y/(x+1.)))
END IF
x2 = x*x
y2 = y*y
r = ((x2-y2-1.)**2 + 4.*x2*y2)**.25
r = SQRT((x + r*COS(angle))**2 + (y + r*SIN(angle))**2)
fn_val = REAL(n+1) + .5*LOG(1./eps)/LOG(r)

RETURN
END FUNCTION nu0jac



FUNCTION nu0lag(n, z, al, eps) RESULT(fn_val)
 
! This is an auxiliary function routine providing a starting backward
! recurrence index for the Laguerre measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.

INTEGER, INTENT(IN)    :: n
COMPLEX, INTENT(IN)    :: z
REAL (dp), INTENT(IN)  :: al
REAL (dp), INTENT(IN)  :: eps
REAL (dp)              :: fn_val

REAL (dp)  :: phi, pi, x, y

pi = 4.*ATAN(1.)
x = REAL(z)
y = AIMAG(z)
phi = .5*pi
IF(y < 0.) phi = 1.5*pi
IF(x == 0.) GO TO 10
phi = ATAN(y/x)
IF(y > 0. .AND. x > 0.) GO TO 10
phi = phi+pi
IF(x < 0.) GO TO 10
phi = phi + pi
10 fn_val = (SQRT(REAL(n+1) + .5*(al+1.)) + LOG(1./eps)/(4.*(x*x+  &
            y*y)**.25*COS(.5*(phi-pi))))**2 - .5*(al+1.)
RETURN
END FUNCTION nu0lag



FUNCTION nu0her(n, z, eps) RESULT(fn_val)
 
! This is an auxiliary function routine providing a starting backward
! recurrence index for the Hermite measure that can be used in place
! of  nu0  in the routines  knum  and  dknum.

INTEGER, INTENT(IN)    :: n
COMPLEX, INTENT(IN)    :: z
REAL (dp), INTENT(IN)  :: eps
REAL (dp)              :: fn_val

fn_val = 2.*(SQRT(.5*REAL(n+1)) + .25*LOG(1./eps)/ ABS(AIMAG(z)))**2
RETURN
END FUNCTION nu0her

END MODULE orthpol
