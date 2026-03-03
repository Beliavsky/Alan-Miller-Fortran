PROGRAM test2
 
! This test generates the first n recursion coefficients for the
! orthogonal polynomials relative to the weight function

!        (x**sigma)*ln(1/x)  on (0,1],  sigma = -.5, 0, .5,

! where n=100 when using modified (Legendre) moments, and n=12 when
! using ordinary moments. It prints the double-precision values of the
! coefficients as well as the relative errors of the respective single-
! precision values and the maximum relative errors.

! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50
! Latest revision - 7 June 2001

USE orthpol
IMPLICIT NONE

REAL       :: a(199), b(199), fnu(200), alpha(100), beta(100), s(100),  &
              s0(200), s1(200), s2(200), sigma
REAL (dp)  :: dsigma, da(199), db(199), dnu(200), dalpha(100), dbeta(100), &
              ds(100), eamax, ebmax, erra, errb
LOGICAL    :: modmom, intexp
INTEGER    :: iderr, ierr, is, k, kamax, kbmax, km1, n, ndm1

modmom  =.true.

! Generate the recursion coefficients for the polynomials defining the
! modified resp. ordinary moments.

10 IF(modmom) THEN
  n = 100
  ndm1 = 2*n-1
  CALL recur(ndm1, 2, 0., 0., a, b, ierr)
  CALL drecur(ndm1, 2, 0.d0, 0.d0, da, db, iderr)
ELSE
  n = 12
  ndm1 = 2*n-1
  DO  k=1,ndm1
    a(k) = 0.
    b(k) = 0.
    da(k) = 0.d0
    db(k) = 0.d0
  END DO
END IF
DO  is=1,3
  dsigma = -.5D0 + .5D0*DBLE(is-1)
  sigma = dsigma
  IF(is == 2) THEN
    intexp = .true.
  ELSE
    intexp = .false.
  END IF
  
! Compute the modified resp. ordinary moments using Eqs. (3.12) and
! (3.11) of the companion paper. On machines with limited exponent
! range, some of the high-order modified moments may underflow, without
! this having any deteriorating effect on the accuracy.
  
  CALL fmm(n, modmom, intexp, sigma, fnu)
  CALL dmm(n, modmom, intexp, dsigma, dnu)
  
! Compute the desired recursion coefficients by means of the modified
! Chebyshev algorithm; for the latter, see, e.g., Section 2.4 of
! W. Gautschi, ``On generating orthogonal polynomials'', SIAM J. Sci.
! Statist. Comput. 3, 1982, 289-317.
  
  CALL cheb(n, a, b, fnu, alpha, beta, s, ierr, s0, s1, s2)
  
! On machines with limited single-precision exponent range, the routine
! cheb  may generate an underflow exception, which however is harmless
! and can be ignored.
  
  CALL dcheb(n, da, db, dnu, dalpha, dbeta, ds, iderr)
  WRITE(*,1) ierr, iderr
  1 FORMAT(/t7, ' ierr in cheb = ', i4, '  iderr in dcheb = ', i4/)
  
! Compute and print the relative errors and their maxima.
  
  eamax = 0.
  ebmax = 0.
  WRITE(*,2) sigma
  2 FORMAT(/t33, 'sigma = ', f5.1)
  WRITE(*,3)
  3 FORMAT(/t4, 'k', t17, 'dalpha(k)', t35, 'dbeta(k)'/)
  DO  k=1,n
    km1 = k-1
    IF(iderr == 0 .OR. km1 < ABS(iderr)) THEN
      WRITE(*,4) km1, dalpha(k), dbeta(k)
      4 FORMAT(' ', i3, 2g24.16)
      IF(ierr == 0 .OR. km1 < ABS(ierr)) THEN
        erra = ABS(alpha(k)-dalpha(k))/dalpha(k)
        errb = ABS(beta(k)-dbeta(k))/dbeta(k)
        WRITE(*,5) erra, errb
        5 FORMAT(t5, e12.4, t39, e12.4)
        IF(erra > eamax) THEN
          eamax = erra
          kamax = km1
        END IF
        IF(errb > ebmax) THEN
          ebmax = errb
          kbmax = km1
        END IF
      END IF
    END IF
  END DO
  WRITE(*,6) eamax, kamax, ebmax, kbmax
  6 FORMAT(/t7, 'eamax =', e11.4, ' at', i3, '         ebmax =', e11.4,  &
           ' at', i3//)
END DO

IF (modmom) THEN
  modmom = .FALSE.
  WRITE(*, '(//)')
  WRITE(*, *) '   output of test2 with modmom=.false.:'
  WRITE(*, *)
  GO TO 10
END IF

STOP


CONTAINS


SUBROUTINE fmm(n, modmom, intexp, sigma, fnu)

! This generates the first  2*n  modified moments (if modmom=.true.)
! relative to shifted monic Legendre polynomials, using Eq. (3.12) of
! the companion paper, and the first  2*n  ordinary moments (if modmom
! =.false.) by Eq. (3.11), of the weight function

!          (x**sigma)*ln(1/x)  on (0,1],   sigma > -1,

! for sigma an integer (if intexp=.true.) or a real number (if intexp
! =.false.).  In either case, the input variable  sigma  is of type real.

INTEGER, INTENT(IN)      :: n
LOGICAL, INTENT(IN)      :: modmom
LOGICAL, INTENT(IN OUT)  :: intexp
REAL, INTENT(IN)         :: sigma
REAL, INTENT(OUT)        :: fnu(:)

REAL     :: sigp1, c, fk, p, s, fi, q, fiq, fkm1
INTEGER  :: i, iq, isigma, isigp1, isigp2, isigp3, k, kmax, km1, nd

! The array  fnu  is assumed to have dimension  2*n.

nd = 2*n
sigp1 = sigma + 1.0
IF(modmom) THEN
  isigma = dsigma
  isigp1 = isigma + 1
  isigp2 = isigma + 2
  isigp3 = isigma + 3
  IF(intexp .AND. isigp1 < nd) THEN
    kmax = isigp1
  ELSE
    kmax = nd
  END IF
  c = 1.d0
  DO  k=1,kmax
    km1 = k-1
    fk = k
    p = 1.0
    s = 1.0/sigp1
    IF(kmax > 1) THEN
      DO  i=1,km1
        fi = i
        p = (sigp1 - fi)*p/(sigp1 + fi)
        s = s + 1.0/(sigp1 + fi) - 1.0/(sigp1 - fi)
      END DO
    END IF
    fnu(k) = c*s*p/sigp1
    c = fk*c/(4.0*fk - 2.0)
  END DO
  IF(.NOT.intexp .OR. isigp1 >= nd) RETURN
  q = -0.5
  IF(isigma > 0) THEN
    DO  iq=1,isigma
      fiq = iq
      q = fiq*fiq*q/((2.0*fiq + 1.0)*(2.0*fiq + 2.0))
    END DO
  END IF
  fnu(isigp2) = c*q
  IF(isigp2 == nd) RETURN
  DO  k=isigp3,nd
    km1 = k-1
    fkm1 = km1
    fnu(k) = -fkm1*(fkm1-sigp1)*fnu(km1)/((4.0*fkm1 - 2.0)* (fkm1 + sigp1))
  END DO
  RETURN
ELSE
  DO  k=1,nd
    fkm1 = k-1
    fnu(k) = (1.0/(sigp1 + fkm1))**2
  END DO
END IF

RETURN
END SUBROUTINE fmm



SUBROUTINE dmm(n, modmom, intexp, dsigma, dnu)

! This is a double-precision version of the routine  fmm.

INTEGER, INTENT(IN)      :: n
LOGICAL, INTENT(IN OUT)  :: modmom
LOGICAL, INTENT(IN OUT)  :: intexp
REAL (dp), INTENT(IN)    :: dsigma
REAL (dp), INTENT(OUT)   :: dnu(:)

REAL (dp) :: dsigp1, dc, dk, dpp, ds, di, dq, diq, dkm1
INTEGER   :: i, iq, isigma, isigp1, isigp2, isigp3, k, kmax, km1, nd

! The array  dnu  is assumed to have dimension  2*n.

nd = 2*n
dsigp1 = dsigma + 1.d0
IF(modmom) THEN
  isigma = dsigma
  isigp1 = isigma + 1
  isigp2 = isigma + 2
  isigp3 = isigma + 3
  IF(intexp .AND. isigp1 < nd) THEN
    kmax = isigp1
  ELSE
    kmax = nd
  END IF
  dc = 1.d0
  DO  k=1,kmax
    km1 = k-1
    dk = k
    dpp = 1.d0
    ds = 1.d0/dsigp1
    IF(kmax > 1) THEN
      DO  i=1,km1
        di = i
        dpp = (dsigp1-di)*dpp/(dsigp1+di)
        ds = ds + 1.d0/(dsigp1+di) - 1.d0/(dsigp1-di)
      END DO
    END IF
    dnu(k) = dc*ds*dpp/dsigp1
    dc = dk*dc/(4.d0*dk-2.d0)
  END DO
  IF(.NOT.intexp .OR. isigp1 >= nd) RETURN
  dq = -.5D0
  IF(isigma > 0) THEN
    DO  iq=1,isigma
      diq = iq
      dq = diq*diq*dq/((2.d0*diq+1.d0)*(2.d0*diq+2.d0))
    END DO
  END IF
  dnu(isigp2) = dc*dq
  IF(isigp2 == nd) RETURN
  DO  k=isigp3,nd
    km1 = k-1
    dkm1 = km1
    dnu(k) = -dkm1*(dkm1-dsigp1)*dnu(km1)/((4.d0*dkm1-2.d0)*(dkm1+dsigp1))
  END DO
  RETURN
ELSE
  DO  k=1,nd
    dkm1 = k-1
    dnu(k) = (1.d0/(dsigp1+dkm1))**2
  END DO
END IF

RETURN
END SUBROUTINE dmm



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



SUBROUTINE cheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)

! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-01  Time: 23:13:48

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
!               fnu-- array of dimension  2*n  to be filled with the
!                     values of the modified moments  fnu(k-1), k=1,2,
!                     ...,2*n
!     Output:   alpha,beta-- arrays containing, respectively, the
!                     recursion coefficients  alpha(k-1),beta(k-1),
!                     k=1,2,...,n, where  beta(0)  is the total mass.
!               s - - array containing the normalization factors
!                     s(k)=integral [pi(k)(x)]**2 dlambda(x), k=0,1,
!                     2,...,n-1.
!               ierr- an error flag, equal to  0  on normal return,
!                     equal to  1  if  abs(fnu(0))  is less than the
!                     machine zero, equal to  2  if  n  is out of range,
!                     equal to  -k  if  s(k), k=0,1,2,...,n-1, is about
!                     to underflow, and equal to  +k  if it is about to
!                     overflow.

! The arrays  s0,s1,s2  are needed for working space.

! On machines with limited exponent range, the occurrence of underflow
! [overflow] in the computation of the  alpha's  and  beta's  can often
! be avoided by multiplying all modified moments by a sufficiently large
! [small] scaling factor and dividing the new  beta(0)  by the same
! scaling factor.

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN)      :: a(:)
REAL, INTENT(IN OUT)  :: b(:)
REAL, INTENT(IN)      :: fnu(:)
REAL, INTENT(OUT)     :: alpha(:)
REAL, INTENT(OUT)     :: beta(:)
REAL, INTENT(OUT)     :: s(:)
INTEGER, INTENT(OUT)  :: ierr
REAL, INTENT(OUT)     :: s0(:)
REAL, INTENT(OUT)     :: s1(:)
REAL, INTENT(OUT)     :: s2(:)

! The arrays  a,b  are assumed to have dimension  2*n-1, the arrays
! fnu,s0,s1,s2  dimension  2*n.

INTEGER  :: k, l, lk, nd
REAL     :: hhuge, ttiny

nd=2*n
ttiny=10.*TINY(0.0)
hhuge=.1*HUGE(0.0)
ierr=0
IF(ABS(fnu(1)) < ttiny) THEN
  ierr=1
  RETURN
END IF
IF(n < 1) THEN
  ierr=2
  RETURN
END IF

! Initialization

alpha(1)=a(1) + fnu(2)/fnu(1)
beta(1)=fnu(1)
IF(n == 1) RETURN
s(1)=fnu(1)
DO  l=1,nd
  s0(l)=0.
  s1(l)=fnu(l)
END DO

! Continuation

DO  k=2,n
  lk=nd-k+1
  DO  l=k,lk

! The quantities  s2(l)  for l > k are auxiliary quantities which may
! be zero or may become so small as to underflow, without however
! causing any harm.

    s2(l)=s1(l+1) - (alpha(k-1) - a(l))*s1(l) - beta(k-1)*s0(l) + b(l)*s1(l-1)
    IF(l == k) s(k)=s2(k)

! Check impending underflow or overflow

    IF(ABS(s(k)) < ttiny) THEN
      ierr=-(k-1)
      RETURN
    ELSE IF(ABS(s(k)) > hhuge) THEN
      ierr=k-1
      RETURN
    END IF
  END DO

! Compute the alpha- and beta-coefficient

  alpha(k)=a(k) + (s2(k+1)/s2(k)) - (s1(k)/s1(k-1))
  beta(k)=s2(k)/s1(k-1)
  DO  l=k,lk
    s0(l)=s1(l)
    s1(l)=s2(l)
  END DO
END DO
RETURN
END SUBROUTINE cheb

END PROGRAM test2

