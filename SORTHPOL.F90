MODULE s_orthpol
! Some of the single precision routines from orthpol

IMPLICIT NONE


CONTAINS


SUBROUTINE mcdis(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,irout,finl,finr,  &
                 endl,endr,xfer,wfer,alpha,beta,ncap,kount,ierr,ie,  &
                 be,x,w,xm,wm,p0,p1,p2)

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-03  Time: 16:33:39

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
!                     ncapm=500  will usually be satisfactory; type integer
!            mc  - -  the number of disjoint intervals in the
!                     continuous part of the spectrum; type integer
!            mp  - -  the number of points in the discrete part of
!                     the spectrum; type integer.  If there is no
!                     point spectrum, set  mp=0.
!            xp  - -  an array of dimension  mp  containing the
!                     abscissas of the point spectrum
!            yp  - -  an array of dimension  mp  containing the jumps
!                     of the point spectrum
!            quad  -  a subroutine determining the discretization of
!                     the inner product on each component interval,
!                     or a dummy routine if  iq  is not equal to  1 (see below)
!            eps  - - the desired relative accuracy of the nonzero
!                     recursion coefficients; type real
!            iq   - - an integer selecting a user-supplied quadrature routine
!                     quad  if  iq=1  or the ORTHPOL routine qgp  otherwise
!            idelta - a nonzero integer, typically  1  or  2, inducing fast
!                     convergence in the case of special quadrature routines
!            irout  - an integer selecting the routine for generating the
!                     recursion coefficients from the discrete inner product.
!                     Specifically,  irout=1  selects the routine  sti,
!                     whereas any other value selects the routine  lancz

! The logical variables  finl,finr, the arrays  endl,endr  of
! dimension  mc, and the arrays  xfer,wfer  of dimension  ncapm  are
! input variables to the subroutine  qgp  and are used (and hence need
! to be properly dimensioned) only if  iq  is not equal to  1.

!    Output:  alpha,beta - arrays of dimension n, holding as k-th element
!                     alpha(k-1), beta(k-1), k=1,2,...,n, respectively
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
! for working space. The routine calls upon the subroutine  sti  or
! lancz, depending on the choice of  irout.

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: ncapm
INTEGER, INTENT(IN)      :: mc
INTEGER, INTENT(IN)      :: mp
REAL, INTENT(IN)         :: xp(:)
REAL, INTENT(IN)         :: yp(:)
REAL, INTENT(IN)         :: eps
INTEGER, INTENT(IN)      :: iq
INTEGER, INTENT(IN OUT)  :: idelta
INTEGER, INTENT(IN)      :: irout
LOGICAL, INTENT(IN)      :: finl
LOGICAL, INTENT(IN)      :: finr
REAL, INTENT(IN)         :: endl(:)
REAL, INTENT(IN)         :: endr(:)
REAL, INTENT(IN OUT)     :: xfer(:)
REAL, INTENT(IN OUT)     :: wfer(:)
REAL, INTENT(IN OUT)     :: alpha(:)
REAL, INTENT(OUT)        :: beta(:)
INTEGER, INTENT(OUT)     :: ncap
INTEGER, INTENT(OUT)     :: kount
INTEGER, INTENT(OUT)     :: ierr
INTEGER, INTENT(IN OUT)  :: ie
REAL, INTENT(OUT)        :: be(:)
REAL, INTENT(IN OUT)     :: x(:)
REAL, INTENT(IN OUT)     :: w(:)
REAL, INTENT(OUT)        :: xm(:)
REAL, INTENT(OUT)        :: wm(:)
REAL, INTENT(IN OUT)     :: p0(:)
REAL, INTENT(IN OUT)     :: p1(:)
REAL, INTENT(IN OUT)     :: p2(:)

! The arrays  xp,yp  are assumed to have dimension  mp  if mp > 0, and
! the arrays  xm,wm,p0,p1,p2  dimension  mc*ncapm+mp.

INTERFACE
  SUBROUTINE quad(n,x,w,i,ierr)
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: n
    REAL, INTENT(IN OUT)  :: x(:)
    REAL, INTENT(IN OUT)  :: w(:)
    INTEGER, INTENT(IN)   :: i
    INTEGER, INTENT(OUT)  :: ierr
  END SUBROUTINE quad
END INTERFACE

INTEGER  :: i, incap, im1tn, k, mtncap

IF (idelta <= 0) idelta = 1
IF (n < 1) THEN
  ierr = -1
  RETURN
END IF

! Initialization

incap = 1
kount = -1
ierr = 0
beta(1:n) = 0.0
ncap = (2*n-1) / idelta
20 be(1:n) = beta(1:n)
kount = kount + 1
IF (kount > 1) incap = 2 ** (kount/5) * n
ncap = ncap + incap
IF (ncap > ncapm) THEN
  ierr = ncapm
  RETURN
END IF

! Discretization of the inner product

mtncap = mc * ncap
DO  i = 1, mc
  im1tn = (i-1) * ncap
  IF (iq == 1) THEN
    CALL quad(ncap,x,w,i,ierr)
  ELSE
    CALL qgp(ncap,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)
  END IF
  IF (ierr /= 0) THEN
    ierr = i
    RETURN
  END IF
  DO  k = 1, ncap
    xm(im1tn+k) = x(k)
    wm(im1tn+k) = w(k)
  END DO
END DO
IF (mp /= 0) THEN
  DO  k = 1, mp
    xm(mtncap+k) = xp(k)
    wm(mtncap+k) = yp(k)
  END DO
END IF

! Computation of the desired recursion coefficients

IF (irout == 1) THEN
  CALL sti(n,mtncap+mp,xm,wm,alpha,beta,ie,p0,p1,p2)
ELSE
  CALL lancz(n,mtncap+mp,xm,wm,alpha,beta,ie,p0,p1)
END IF

! In the following statement, the absolute value of the beta's is
! used to guard against failure in cases where the routine is applied
! to variable-sign weight functions and hence the positivity of the
! beta's is not guaranteed.

DO  k = 1, n
  IF (ABS(beta(k)-be(k)) > eps*ABS(beta(k))) GO TO 20
END DO

RETURN
END SUBROUTINE mcdis



SUBROUTINE gauss(n,alpha,beta,eps,zero,weight,ierr,e)

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-03  Time: 22:44:07

! Given  n  and a measure  dlambda, this routine generates the n-point
! Gaussian quadrature formula

!     integral over supp(dlambda) of f(x)dlambda(x)

!        = sum from k=1 to k=n of w(k)f(x(k)) + R(n;f).

! The nodes are returned as  zero(k)=x(k) and the weights as
! weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion
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
!                alpha,beta - - arrays of dimension  n  to be filled
!                      with the values of  alpha(k-1), beta(k-1), k=1,2,
!                      ...,n
!                eps - the relative accuracy desired in the nodes
!                      and weights

!        Output: zero- array of dimension  n  containing the Gaussian
!                      nodes (in increasing order)  zero(k)=x(k), k=1,2,
!                      ...,n
!                weight - array of dimension  n  containing the
!                      Gaussian weights  weight(k)=w(k), k=1,2,...,n
!                ierr- an error flag equal to  0  on normal return,
!                      equal to  i  if the QR algorithm does not
!                      converge within 30 iterations on evaluating the
!                      i-th eigenvalue, equal to  -1  if  n  is not in
!                      range, and equal to  -2  if one of the beta's is
!                      negative.

! The array  e  is needed for working space.

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: alpha(n)
REAL, INTENT(IN)      :: beta(n)
REAL, INTENT(IN)      :: eps
REAL, INTENT(OUT)     :: zero(n)
REAL, INTENT(OUT)     :: weight(n)
INTEGER, INTENT(OUT)  :: ierr
REAL, INTENT(OUT)     :: e(n)

INTEGER  :: i, ii, j, k, l, m, mml
REAL     :: b, c, f, g, p, r, s

IF (n < 1) THEN
  ierr = -1
  RETURN
END IF
ierr = 0
zero(1) = alpha(1)
IF (beta(1) < 0.) THEN
  ierr = -2
  RETURN
END IF
weight(1) = beta(1)
IF (n == 1) RETURN
weight(1) = 1.
e(n) = 0.
DO  k = 2, n
  zero(k) = alpha(k)
  IF (beta(k) < 0.) THEN
    ierr = -2
    RETURN
  END IF
  e(k-1) = SQRT(beta(k))
  weight(k) = 0.
END DO
DO  l = 1, n
  j = 0

! Look for a small subdiagonal element.

  20 DO  m = l, n
    IF (m == n) EXIT
    IF (ABS(e(m)) <= eps*(ABS(zero(m))+ABS(zero(m+1)))) EXIT
  END DO
  p = zero(l)
  IF (m /= l) THEN
    IF (j == 30) GO TO 100
    j = j + 1

! Form shift.

    g = (zero(l+1)-p) / (2.*e(l))
    r = SQRT(g*g+1.)
    g = zero(m) - p + e(l) / (g+SIGN(r,g))
    s = 1.
    c = 1.
    p = 0.
    mml = m - l

! For i=m-1 step -1 until l do ...

    DO  ii = 1, mml
      i = m - ii
      f = s * e(i)
      b = c * e(i)
      IF (ABS(f) >= ABS(g)) THEN
        c = g / f
        r = SQRT(c*c+1.)
        e(i+1) = f * r
        s = 1. / r
        c = c * s
      ELSE
        s = f / g
        r = SQRT(s*s+1.)
        e(i+1) = g * r
        c = 1. / r
        s = s * c
      END IF
      g = zero(i+1) - p
      r = (zero(i)-g) * s + 2. * c * b
      p = s * r
      zero(i+1) = g + p
      g = c * r - b

! Form first component of vector.

      f = weight(i+1)
      weight(i+1) = s * weight(i) + c * f
      weight(i) = c * weight(i) - s * f
    END DO
    zero(l) = zero(l) - p
    e(l) = g
    e(m) = 0.
    GO TO 20
  END IF
END DO

! Order eigenvalues and eigenvectors.

DO  ii = 2, n
  i = ii - 1
  k = i
  p = zero(i)
  DO  j = ii, n
    IF (zero(j) < p) THEN
      k = j
      p = zero(j)
    END IF
  END DO
  IF (k /= i) THEN
    zero(k) = zero(i)
    zero(i) = p
    p = weight(i)
    weight(i) = weight(k)
    weight(k) = p
  END IF
END DO
DO  k = 1, n
  weight(k) = beta(1) * weight(k) * weight(k)
END DO
RETURN

! Set error - no convergence to an eigenvalue after 30 iterations.

100 ierr = l
RETURN
END SUBROUTINE gauss



SUBROUTINE sti(n,ncap,x,w,alpha,beta,ierr,p0,p1,p2)

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-03  Time: 22:44:12

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

! The routine uses the function subroutine  r1mach.

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: ncap
REAL, INTENT(IN)      :: x(ncap)
REAL, INTENT(IN)      :: w(ncap)
REAL, INTENT(OUT)     :: alpha(n)
REAL, INTENT(OUT)     :: beta(n)
INTEGER, INTENT(OUT)  :: ierr
REAL, INTENT(OUT)     :: p0(ncap)
REAL, INTENT(OUT)     :: p1(ncap)
REAL, INTENT(OUT)     :: p2(ncap)

INTEGER  :: k, m, nm1
REAL     :: hhuge, sum0, sum1, sum2, t, ttiny

ttiny = 10. * TINY(0.0)
hhuge = .1 * HUGE(0.0)
ierr = 0
IF (n <= 0 .OR. n > ncap) THEN
  ierr = 1
  RETURN
END IF
nm1 = n - 1

! Compute the first alpha- and beta-coefficient.

sum0 = 0.
sum1 = 0.
DO  m = 1, ncap
  sum0 = sum0 + w(m)
  sum1 = sum1 + w(m) * x(m)
END DO
alpha(1) = sum1 / sum0
beta(1) = sum0
IF (n == 1) RETURN

! Compute the remaining alpha- and beta-coefficients.

DO  m = 1, ncap
  p1(m) = 0.
  p2(m) = 1.
END DO
DO  k = 1, nm1
  sum1 = 0.
  sum2 = 0.
  DO  m = 1, ncap

! The following statement is designed to avoid an overflow condition
! in the computation of  p2(m)  when the weights  w(m)  go to zero
! faster (and underflow) than the  p2(m)  grow.

    IF (w(m) /= 0.) THEN
      p0(m) = p1(m)
      p1(m) = p2(m)
      p2(m) = (x(m)-alpha(k)) * p1(m) - beta(k) * p0(m)

! Check for impending overflow.

      IF (ABS(p2(m)) > hhuge .OR. ABS(sum2) > hhuge) THEN
        ierr = k
        RETURN
      END IF
      t = w(m) * p2(m) * p2(m)
      sum1 = sum1 + t
      sum2 = sum2 + t * x(m)
    END IF
  END DO

! Check for impending underflow.

  IF (ABS(sum1) < ttiny) THEN
    ierr = -k
    RETURN
  END IF
  alpha(k+1) = sum2 / sum1
  beta(k+1) = sum1 / sum0
  sum0 = sum1
END DO
RETURN
END SUBROUTINE sti



SUBROUTINE lancz(n,ncap,x,w,alpha,beta,ierr,p0,p1)

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-03  Time: 22:44:20

! This routine carries out the same task as the routine  sti, but
! uses the more stable Lanczos method. The meaning of the input
! and output parameters is the same as in the routine  sti. (This
! routine is adapted from the routine RKPW in W.B. Gragg and
! W.J. Harrod,``The numerically stable reconstruction of Jacobi
! matrices from spectral data'', Numer. Math. 44, 1984, 317-335.)


INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: ncap
REAL, INTENT(IN)      :: x(ncap)
REAL, INTENT(IN)      :: w(ncap)
REAL, INTENT(OUT)     :: alpha(n)
REAL, INTENT(OUT)     :: beta(n)
INTEGER, INTENT(OUT)  :: ierr
REAL, INTENT(OUT)     :: p0(ncap)
REAL, INTENT(OUT)     :: p1(ncap)

INTEGER  :: i, k
REAL     :: gam, pi, rho, sig, t, tk, tmp, tsig, xlam

IF (n <= 0 .OR. n > ncap) THEN
  ierr = 1
  RETURN
ELSE
  ierr = 0
END IF
DO  i = 1, ncap
  p0(i) = x(i)
  p1(i) = 0.
END DO
p1(1) = w(1)
DO  i = 1, ncap - 1
  pi = w(i+1)
  gam = 1.
  sig = 0.
  t = 0.
  xlam = x(i+1)
  DO  k = 1, i + 1
    rho = p1(k) + pi
    tmp = gam * rho
    tsig = sig
    IF (rho <= 0.) THEN
      gam = 1.
      sig = 0.
    ELSE
      gam = p1(k) / rho
      sig = pi / rho
    END IF
    tk = sig * (p0(k)-xlam) - gam * t
    p0(k) = p0(k) - (tk-t)
    t = tk
    IF (sig <= 0.) THEN
      pi = tsig * p1(k)
    ELSE
      pi = (t**2) / sig
    END IF
    tsig = sig
    p1(k) = tmp
  END DO
END DO

DO  k = 1, n
  alpha(k) = p0(k)
  beta(k) = p1(k)
END DO

RETURN
END SUBROUTINE lancz



SUBROUTINE qgp(n,x,w,i,ierr,mc,finl,finr,endl,endr,xfer,wfer)

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-05  Time: 11:17:33

! This is a general-purpose discretization routine that can be used
! as an alternative to the routine  quad  in the multiple-component
! discretization procedure  mcdis. It takes no account of the special
! nature of the weight function involved and hence may result in slow
! convergence of the discretization procedure. This routine, therefore,
! should be used only as a last resort, when no better, more natural
! discretization can be found.

! It is assumed that there are  mc.ge.1  disjoint component intervals.
! The discretization is effected by the Fejer quadrature rule,
! suitably transformed to the respective interval. An interval that
! extends to minus infinity has to be indexed by  1; one that extends
! to plus infinity has to be indexed by  mc.

! The output variable  ierr  is given the value  0. Additional input
! parameters and working space used by this routine are as follows:

!          mc      - the number of component intervals; type integer
!          finl    - a logical variable to be set .true. if the
!                    extreme left interval is finite, and .false.
!                    otherwise
!          finr    - a logical variable to be set .true. if the
!                    extreme right interval is finite, and .false.
!                    otherwise
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

!                     function wf(x,i),

! which evaluates the weight function at the point  x  on the i-th
! component interval. The routine also uses the subroutines  fejer,
! symtr  and  tr, which are appended.

INTEGER, INTENT(IN)   :: n
REAL, INTENT(OUT)     :: x(n)
REAL, INTENT(OUT)     :: w(n)
INTEGER, INTENT(IN)   :: i
INTEGER, INTENT(OUT)  :: ierr
INTEGER, INTENT(IN)   :: mc
LOGICAL, INTENT(IN)   :: finl
LOGICAL, INTENT(IN)   :: finr
REAL, INTENT(IN)      :: endl(mc)
REAL, INTENT(IN)      :: endr(mc)
REAL, INTENT(OUT)     :: xfer(:)
REAL, INTENT(OUT)     :: wfer(:)

! The arrays  xfer,wfer  are dimensioned in the routine  mcdis.

INTERFACE
  FUNCTION wf(x, i) RESULT(fn_val)
    IMPLICIT NONE
    REAL, INTENT(IN)     :: x
    INTEGER, INTENT(IN)  :: i
    REAL                 :: fn_val
  END FUNCTION wf
END INTERFACE

INTEGER  :: k
REAL     :: phi, phi1

ierr = 0
IF (i == 1) CALL fejer(n,xfer,wfer)
IF (i <= 1 .OR. i >= mc) THEN
  IF (mc == 1) THEN
    IF (finl .AND. finr) GO TO 60
    IF (finl) GO TO 20
    IF (finr) GO TO 40
    DO  k = 1, n
      CALL symtr(xfer(k),phi,phi1)
      x(k) = phi
      w(k) = wfer(k) * wf(phi,i) * phi1
    END DO
    RETURN
  ELSE
    IF ((i == 1 .AND. finl) .OR. (i == mc .AND. finr)) GO TO 60
    IF (i == 1) GO TO 40
  END IF
  20 DO  k = 1, n
    CALL tr(xfer(k),phi,phi1)
    x(k) = endl(mc) + phi
    w(k) = wfer(k) * wf(x(k),mc) * phi1
  END DO
  RETURN
  40 DO  k = 1, n
    CALL tr(-xfer(k),phi,phi1)
    x(k) = endr(1) - phi
    w(k) = wfer(k) * wf(x(k),1) * phi1
  END DO
  RETURN
END IF

60 DO  k = 1, n
  x(k) = .5 * ((endr(i)-endl(i))*xfer(k) + endr(i) + endl(i))
  w(k) = .5 * (endr(i)-endl(i)) * wfer(k) * wf(x(k),i)
END DO
RETURN
END SUBROUTINE qgp



SUBROUTINE symtr(t,phi,phi1)

! This implements a particular transformation  x=phi(t)  mapping
! the t-interval [-1,1] to the x-interval [-oo,oo].

!        input:   t
!        output:  phi=phi(t)
!                 phi1=derivative of phi(t)

REAL, INTENT(IN)   :: t
REAL, INTENT(OUT)  :: phi, phi1

REAL  :: t2

t2 = t * t
phi = t / (1.-t2)
phi1 = (t2+1.) / (t2-1.) ** 2
RETURN
END SUBROUTINE symtr



SUBROUTINE tr(t,phi,phi1)

! This implements a particular transformation  x=phi(t)  mapping
! the t-interval [-1,1] to the x-interval [0,oo].

!         input:   t
!         output:  phi=phi(t)
!                  phi1=derivative of phi(t)

REAL, INTENT(IN)   :: t
REAL, INTENT(OUT)  :: phi, phi1

phi = (1.+t) / (1.-t)
phi1 = 2. / (t-1.) ** 2
RETURN
END SUBROUTINE tr



SUBROUTINE fejer(n,x,w)

! This routine generates the n-point Fejer quadrature rule.

!   input:   n   - the number of quadrature nodes
!   output:  x,w - arrays of dimension  n  holding the quadrature nodes and
!                  weights, respectively; the nodes are ordered increasingly

INTEGER, INTENT(IN)  :: n
REAL, INTENT(OUT)    :: x(n)
REAL, INTENT(OUT)    :: w(n)

INTEGER  :: k, m, nh, np1h
REAL     :: c0, c1, c2, fn, pi, sum, t

pi = 4. * ATAN(1.)
nh = n / 2
np1h = (n+1) / 2
fn = n
DO  k = 1, nh
  x(n+1-k) = COS(.5*REAL(2*k-1)*pi/fn)
  x(k) = -x(n+1-k)
END DO
IF (2*nh /= n) x(np1h) = 0.
DO  k = 1, np1h
  c1 = 1.
  c0 = 2. * x(k) * x(k) - 1.
  t = 2. * c0
  sum = c0 / 3.
  DO  m = 2, nh
    c2 = c1
    c1 = c0
    c0 = t * c1 - c2
    sum = sum + c0 / REAL(4*m*m-1)
  END DO
  w(k) = 2. * (1.-2.*sum) / fn
  w(n+1-k) = w(k)
END DO
RETURN
END SUBROUTINE fejer



SUBROUTINE recur(n,ipoly,al,be,a,b,ierr)

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-01  Time: 23:13:41

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
!                     7=generalized Laguerre polynomial with parameter  al
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

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: ipoly
REAL, INTENT(IN)      :: al
REAL, INTENT(IN)      :: be
REAL, INTENT(OUT)     :: a(:)
REAL, INTENT(OUT)     :: b(:)
INTEGER, INTENT(OUT)  :: ierr

INTEGER  :: k
REAL     :: almach, alpbe, al2, be2, fkm1, t

IF(n < 1) THEN
  ierr=3
  RETURN
END IF
almach=LOG(HUGE(0.0))
ierr=0
a(1:n)=0.
IF(ipoly == 1) THEN
  b(1)=2.
  IF (n == 1) RETURN
  DO  k=2,n
    fkm1=k-1
    b(k)=1./(4.-1./(fkm1*fkm1))
  END DO
  RETURN
ELSE IF (ipoly == 2) THEN
  a(1)=.5
  b(1)=1.
  IF(n == 1) RETURN
  DO  k=2,n
    a(k)=.5
    fkm1=REAL(k-1)
    b(k)=.25/(4.-1./(fkm1*fkm1))
  END DO
  RETURN
ELSE IF(ipoly == 3) THEN
  b(1)=4.*ATAN(1.)
  IF(n == 1) RETURN
  b(2)=.5
  IF(n == 2) RETURN
  DO  k=3,n
    b(k)=.25
  END DO
  RETURN
ELSE IF(ipoly == 4) THEN
  b(1)=2.*ATAN(1.)
  IF(n == 1) RETURN
  DO  k=2,n
    b(k)=.25
  END DO
  RETURN
ELSE IF(ipoly == 5) THEN
  b(1)=4.*ATAN(1.)
  a(1)=.5
  IF(n == 1) RETURN
  DO  k=2,n
    b(k)=.25
  END DO
  RETURN
ELSE IF(ipoly == 6) THEN
  IF(al <= -1. .OR. be <= -1.) THEN
    ierr=1
    RETURN
  ELSE
    alpbe=al+be
    a(1)=(be-al)/(alpbe+2.)
    t=(alpbe+1.)*LOG(2.) + alga(al+1.) + alga(be+1.) - alga(alpbe+2.)
    IF(t > almach) THEN
      ierr=2
      b(1)=HUGE(0.0)
    ELSE
      b(1)=EXP(t)
    END IF
    IF(n == 1) RETURN
    al2=al*al
    be2=be*be
    a(2)=(be2-al2)/((alpbe+2.)*(alpbe+4.))
    b(2)=4.*(al+1.)*(be+1.)/((alpbe+3.)*(alpbe+2.)**2)
    IF(n == 2) RETURN
    DO  k=3,n
      fkm1=REAL(k-1)
      a(k)=.25*(be2-al2)/(fkm1*fkm1*(1.+.5*alpbe/fkm1)*  &
          (1.+.5*(alpbe+2.)/fkm1))
      b(k)=.25*(1.+al/fkm1)*(1.+be/fkm1)*(1.+alpbe/fkm1)/  &
          ((1.+.5*(alpbe+1.)/fkm1)*(1.+.5*(alpbe-1.)/fkm1)*(1.+.5*alpbe/fkm1)**2)
    END DO
    RETURN
  END IF
ELSE IF(ipoly == 7) THEN
  IF(al <= -1.) THEN
    ierr=1
    RETURN
  ELSE
    a(1)=al+1.
    b(1)=gamma(al+1.)
    IF(ierr == 2) b(1)=HUGE(0.0)
    IF(n == 1) RETURN
    DO  k=2,n
      fkm1=REAL(k-1)
      a(k)=2.*fkm1+al+1.
      b(k)=fkm1*(fkm1+al)
    END DO
    RETURN
  END IF
ELSE IF(ipoly == 8) THEN
  b(1)=SQRT(4.*ATAN(1.))
  IF(n == 1) RETURN
  DO  k=2,n
    b(k)=.5*REAL(k-1)
  END DO
  RETURN
ELSE
  ierr=4
END IF

RETURN
END SUBROUTINE recur



FUNCTION alga(x) RESULT(fn_val)

! This is an auxiliary function subroutine (not optimized in any
! sense) evaluating the logarithm of the gamma function for positive
! arguments  x. It is called by the subroutine  gamma. The integer  m0
! in the first executable statement is the smallest integer  m  such
! that  1*3*5* ... *(2*m+1)/(2**m)  is greater than or equal to the
! largest machine-representable number. The routine is based on a
! rational approximation valid on [.5,1.5] due to W.J. Cody and
! K.E. Hillstrom; see Math. Comp. 21, 1967, 198-203, in particular the
! case  n=7  in Table II. For the computation of  m0  it calls upon the
! function subroutines  t  and  r1mach. The former, appended below,
! evaluates the inverse function  t = t(y)  of  y = t ln t.

REAL, INTENT(IN)  :: x
REAL              :: fn_val

REAL, PARAMETER  :: cnum(8) = (/   &
       4.120843185847770, 85.68982062831317, 243.175243524421,  &
    -261.7218583856145, -922.2613728801522, -517.6383498023218,  &
     -77.41064071332953,  -2.208843997216182 /),  &
    cden(8) = (/1., 45.64677187585908, 377.8372484823942, 951.323597679706,  &
                846.0755362020782, 262.3083470269460, 24.43519662506312,   &
                0.4097792921092615 /)

INTEGER  :: k, m, m0, mm1
REAL     :: p, sden, snum, xe, xi

! The constants in the statement below are  exp(1.)  and  .5*LOG(8.).

m0=2.71828*t((LOG(HUGE(0.0))-1.03972)/2.71828)
xi=x
IF(x-xi > .5) xi=xi + 1.
m=xi - 1

! Computation of log gamma on the standard interval (1/2, 3/2]

xe=x - m
snum=cnum(1)
sden=cden(1)
DO  k=2,8
  snum=xe*snum + cnum(k)
  sden=xe*sden + cden(k)
END DO
fn_val=(xe-1.)*snum/sden

! Computation of log gamma on (0,1/2]

IF(m == -1) THEN
  fn_val=fn_val - LOG(x)
  RETURN
ELSE IF(m == 0) THEN
  RETURN
ELSE

! Computation of log gamma on (3/2,5/2]

  p=xe
  IF(m == 1) THEN
    fn_val=fn_val + LOG(p)
    RETURN
  ELSE

! Computation of log gamma for arguments larger than 5/2

    mm1=m-1

! The else-clause in the next statement is designed to avoid possible
! overflow in the computation of  p  in the if-clause, at the expense
! of computing many logarithms.

    IF(m < m0) THEN
      DO  k=1,mm1
        p=(xe + k)*p
      END DO
      fn_val=fn_val + LOG(p)
      RETURN
    ELSE
      fn_val=fn_val + LOG(xe)
      DO  k=1,mm1
        fn_val=fn_val + LOG(xe + k)
      END DO
      RETURN
    END IF
  END IF
END IF

RETURN
END FUNCTION alga



FUNCTION gamma(x) RESULT(fn_val)

! This evaluates the gamma function for real positive  x, using the
! function subroutines  alga  and  r1mach. In case of overflow, the
! routine returns the largest machine-representable number and the
! error flag  ierr=2.

REAL, INTENT(IN)      :: x
REAL                  :: fn_val

REAL     :: almach, t
INTEGER  :: ierr

almach=LOG(HUGE(0.0))
ierr=0
t=alga(x)
IF(t >= almach) THEN
  ierr=2
  fn_val=HUGE(0.0)
ELSE
  fn_val=EXP(t)
END IF

IF (ierr /= 0) WRITE(*, *) ' ** Overflow in function GAMMA **'

RETURN
END FUNCTION gamma



FUNCTION t(y) RESULT(fn_val)

! This evaluates the inverse function  t = t(y)  of y = t ln t  for
! nonnegative  y  to an accuracy of about one percent. For the
! approximation used, see pp. 51-52 in W. Gautschi,``Computational
! aspects of three-term recurrence relations'', SIAM Rev. 9, 1967, 24-82.

REAL, INTENT(IN)  :: y
REAL              :: fn_val

REAL  :: p, z

IF(y <= 10.) THEN
  p=.000057941*y-.00176148
  p=y*p+.0208645
  p=y*p-.129013
  p=y*p+.85777
  fn_val=y*p+1.0125
ELSE
  z=LOG(y)-.775
  p=(.775-LOG(z))/(1.+z)
  p=1./(1.+p)
  fn_val=y*p/z
END IF
RETURN
END FUNCTION t

END MODULE s_orthpol
