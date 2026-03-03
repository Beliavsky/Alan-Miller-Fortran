PROGRAM test1
USE orthpol, ONLY: dp, drecur, dcheb
IMPLICIT NONE

! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50
! Latest revision - 7 June 2001

REAL       :: fnu(160),f(80),f0(80),rr(80),a(159),b(159),alpha(80),  &
              beta(80),s(80),s0(160),s1(160),s2(160)
REAL (dp)  :: deps,dom2,dnu(160),d(80),d0(80),drr(80),da(159),db(159),  &
              dalpha(80),dbeta(80),ds(80)
REAL       :: eps, errb, om2
INTEGER    :: iderr, ierr, iom, k, km1, n, ndm1
LOGICAL    :: modmom
REAL (dp)  :: doom2(7) = (/ .1D0,.3D0,.5D0,.7D0,.9D0,.99D0,.999D0 /)

! This test generates the first n beta-coefficients in the recurrence
! relation for the orthogonal polynomials relative to the weight function:

!         ((1-om2*x**2)*(1-x**2))**(-1/2)  on (-1,1)

! for om2=.1(.2).9,.99,.999, both in single and double precision,
! using modified moments if  modmom=.true.  and ordinary moments
! otherwise. In the former case, n=80, in the latter, n=20. Printed
! are the double-precision values of the coefficients along with the
! relative errors of the single-precision values.

WRITE(*,1)
1 FORMAT(/)
WRITE(*, '(t6, a)') 'output of test1 with modmom=.true.:'
WRITE(*,1)
modmom=.TRUE.
eps=EPSILON(0.0)
deps=EPSILON(0.0_dp)

10 IF(modmom) THEN
  n=80
ELSE
  n=20
END IF
ndm1=2*n-1
DO  iom=1,7
  dom2=doom2(iom)
  om2=dom2
  
! Compute the modified resp. ordinary moments using Eqs. (3.7) and (3.9)
! of the companion paper. On machines with limited exponent range, some
! of the high-order modified moments may underflow, without this having
! any deteriorating effect on the accuracy.
  
  CALL fmm(n,eps,modmom,om2,fnu,ierr,f,f0,rr)
  CALL dmm(n,deps,modmom,dom2,dnu,iderr,d,d0,drr)
  IF(ierr /= 0 .OR. iderr /= 0) THEN
    WRITE(*,2) ierr,iderr,om2
    2 FORMAT(/t6, 'ierr in fmm = ', i1, '  iderr in dmm = ', i1,  &
             '  for om2 = ', f8.4/)
    CYCLE
  END IF
  
! Generate the recursion coefficients for the polynomials defining the
! modified resp. ordinary moments.
  
  IF(modmom) THEN
    CALL recur(ndm1,3,0.,0.,a,b,ierr)
    CALL drecur(ndm1,3,0.d0,0.d0,da,db,iderr)
  ELSE
    DO  k=1,ndm1
      a(k)=0.
      b(k)=0.
      da(k)=0.d0
      db(k)=0.d0
    END DO
  END IF
  
! Compute the desired recursion coefficients by means of the modified
! Chebyshev algorithm; for the latter, see, e.g., Section 2.4 of
! W. Gautschi, ``On generating orthogonal polynomials'', SIAM J. Sci.
! Statist. Comput. 3, 1982, 289-317.
  
  CALL cheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
  
! On machines with limited single-precision exponent range, the routine
! cheb  may generate an underflow exception, which however is harmless
! and can be ignored.
  
  CALL dcheb(n,da,db,dnu,dalpha,dbeta,ds,iderr)

  WRITE(*,3) ierr,iderr
  3 FORMAT(/t6, 'ierr in cheb = ', i3, '  iderr in dcheb = ', i3/)
  WRITE(*,4)
  4 FORMAT(/t6, 'k', t21, 'dbeta(k)'/)
  DO  k=1,n
    km1=k-1
    IF(iderr == 0 .OR. km1 < ABS(iderr)) THEN
      IF(ierr == 0 .OR. km1 < ABS(ierr)) THEN
        errb=ABS(beta(k)-dbeta(k))/dbeta(k)
        IF(k == 1) THEN
          WRITE(*,5) km1,dbeta(k),errb,om2
          5 FORMAT(' ', i5, g24.16, e12.4, '   om2 =', f6.3)
        ELSE
          WRITE(*,6) km1,dbeta(k),errb
          6 FORMAT(' ', i5, g24.16, e12.4)
        END IF
      ELSE
        WRITE(*,7) km1,dbeta(k)
        7 FORMAT(' ', i5, g24.16)
      END IF
    END IF
  END DO
  WRITE(*,1)
END DO

IF (modmom) THEN
  modmom = .FALSE.
  WRITE (*, 1)
  WRITE(*, '(t6, a)') 'output of test1 with modmom=.false.:'
  WRITE (*, 1)
  GO TO 10
END IF

STOP

CONTAINS


SUBROUTINE fmm(n,eps,modmom,om2,fnu,ierr,f,f0,rr)

! This routine generates the modified (Chebyshev) resp. ordinary
! moments of the weight function

!          ((1-om2*x**2)*(1-x**2))**(-1/2)  on (-1,1)

! using Eqs. (3.7) resp. (3.9) of the companion paper.

INTEGER, INTENT(IN)      :: n
REAL, INTENT(IN)         :: eps
LOGICAL, INTENT(IN OUT)  :: modmom
REAL, INTENT(IN)         :: om2
REAL, INTENT(OUT)        :: fnu(:)
INTEGER, INTENT(OUT)     :: ierr
REAL, INTENT(OUT)        :: f(:)
REAL, INTENT(OUT)        :: f0(:)
REAL, INTENT(OUT)        :: rr(:)

! The array  fnu  is assumed to have dimension  2*n.

INTEGER  :: i, k, k1, k1m1, n1, nd, ndm1, nu
REAL     :: c, c0, c1, fn1, pi, q, q1, r, s, sum

ierr=0
nd=2*n
ndm1=nd-1
pi=4.*ATAN(1.)

! Compute the Fourier coefficients of ((1-om2*sin(theta)**2))**(-1/2)
! as minimal solution of a three-term recurrence relation as described
! on pp.310-311 of W. Gautschi,``On generating orthogonal polynomials'',
! SIAM J. Sci. Statist. Comput. 3, 1982, 289-317.

q=om2/(2. - om2 + 2.*SQRT(1.-om2))
q1=(1. + q*q)/q
f(1:n)=0.
nu=nd

20 nu=nu+10
f0(1:n)=f(1:n)
IF(nu > 500) THEN
  ierr=1
  RETURN
END IF
r=0.
s=0.
DO  k=1,nu
  n1=nu-k+1
  fn1=n1
  r=-(fn1-.5)/(fn1*q1 + (fn1+.5)*r)
  s=r*(2.+s)
  IF(n1 <= n) rr(n1)=r
END DO
c0=1./(1.+s)
f(1)=rr(1)*c0
IF(n > 1) THEN
  DO  k=2,n
    f(k)=rr(k)*f(k-1)
  END DO
END IF
DO  k=1,n
  IF(ABS(f(k)-f0(k)) > eps*ABS(f(k))) GO TO 20
END DO

! Compute the desired modified resp. ordinary moments in term of
! the above Fourier coefficients.

fnu(1)=pi*c0
IF(n == 1) RETURN
fnu(2)=0.
IF(n == 2) RETURN
IF(modmom) THEN
  c=2.*pi
  DO  k=3,ndm1,2
    k1=(k-1)/2
    c=-.25*c
    fnu(k)=c*f(k1)
    fnu(k+1)=0.
  END DO
ELSE
  c=.5*pi
  fnu(3)=c*(c0-f(1))
  fnu(4)=0.
  c=-c
  DO  k=5,ndm1,2
    k1=(k-1)/2
    k1m1=k1-1
    c=-.25*c
    c1=1.
    sum=f(k1)
    DO  i=1,k1m1
      c1=-c1*REAL(2*k1-i+1)/REAL(i)
      sum=sum + c1*f(k1-i)
    END DO
    c1=-c1*REAL(k1+1)/REAL(2*k1)
    sum=sum + c1*c0
    fnu(k)=c*sum
    fnu(k+1)=0.
  END DO
END IF

RETURN
END SUBROUTINE fmm



SUBROUTINE dmm(n,deps,modmom,dom2,dnu,ierrd,d,d0,drr)

! This is a double-precision version of the routine  fmm.

INTEGER, INTENT(IN)      :: n
REAL (dp), INTENT(IN)    :: deps
LOGICAL, INTENT(IN)      :: modmom
REAL (dp), INTENT(IN)    :: dom2
REAL (dp), INTENT(OUT)   :: dnu(:)
INTEGER, INTENT(OUT)     :: ierrd
REAL (dp), INTENT(OUT)   :: d(:)
REAL (dp), INTENT(OUT)   :: d0(:)
REAL (dp), INTENT(OUT)   :: drr(:)

! The array  dnu  is assumed to have dimension  2*n.

REAL (dp) :: dpi,dq, dq1,dr,ds,dn1,dc0,dc,dc1,dsum
INTEGER   :: i, k, k1, k1m1, n1, nd, ndm1, nud

ierrd=0
nd=2*n
ndm1=nd-1
dpi=4.d0*ATAN(1.d0)
dq=dom2/(2.d0 - dom2 + 2.d0*SQRT(1.d0-dom2))
dq1=(1.d0 + dq*dq)/dq
d(1:n)=0.d0
nud=nd

20 nud=nud+10
d0(1:n)=d(1:n)
IF(nud > 1000) THEN
  ierrd=1
  RETURN
END IF
dr=0.d0
ds=0.d0
DO  k=1,nud
  n1=nud-k+1
  dn1=DBLE(n1)
  dr=-(dn1-.5D0)/(dn1*dq1 + (dn1+.5D0)*dr)
  ds=dr*(2.d0 + ds)
  IF(n1 <= n) drr(n1)=dr
END DO
dc0=1.d0/(1.d0+ds)
d(1)=drr(1)*dc0
IF(n > 1) THEN
  DO  k=2,n
    d(k)=drr(k)*d(k-1)
  END DO
END IF
DO  k=1,n
  IF(ABS(d(k)-d0(k)) > deps*ABS(d(k))) GO TO 20
END DO

dnu(1)=dpi*dc0
IF(n == 1) RETURN
dnu(2)=0.d0
IF(n == 2) RETURN
IF(modmom) THEN
  dc=2.d0*dpi
  DO  k=3,ndm1,2
    k1=(k-1)/2
    dc=-.25D0*dc
    dnu(k)=dc*d(k1)
    dnu(k+1)=0.d0
  END DO
ELSE
  dc=.5D0*dpi
  dnu(3)=dc*(dc0-d(1))
  dnu(4)=0.d0
  dc=-dc
  DO  k=5,ndm1,2
    k1=(k-1)/2
    k1m1=k1-1
    dc=-.25D0*dc
    dc1=1.d0
    dsum=d(k1)
    DO  i=1,k1m1
      dc1=-dc1*DBLE(2*k1-i+1)/DBLE(i)
      dsum=dsum+dc1*d(k1-i)
    END DO
    dc1=-dc1*DBLE(k1+1)/DBLE(2*k1)
    dsum=dsum+dc1*dc0
    dnu(k)=dc*dsum
    dnu(k+1)=0.d0
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

END PROGRAM test1

