PROGRAM test10
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50

USE orthpol
USE s_orthpol
IMPLICIT NONE

REAL      :: a(31),b(31),alpha(31),beta(31),z(11),w(11),e(11),  &
             a1(31),b1(31),betap(20,12),erram(12),errbm(12),erra,errb,  &
             epsma
REAL (dp) :: depsma,da(31),db(31),dalpha(31),dbeta(31),  &
             dz(11),dw(11),da1(31),db1(31)
INTEGER   :: ierr, im, ip, ip4, k, km1, m, n, npm

! epsma and depsma are the machine single and REAL (dp).

epsma=EPSILON(0.0)
depsma=EPSILON(0.0_dp)

! This test applies the routines  indp  and  dindp  to generate the
! first 20 recursion coefficients of the induced Legendre polynomials
! pind(k,m)(.), m=0,1,2,...,11, that is, of the polynomials orthogonal
! relative to the weight function

!                 [p(m)(x)]**2    on [-1,1],

! where  p(m)(.)  is the (monic) Legendre polynomial of degree m.
! (When m=0, then  pind(k,0)(.)=p(k)(.).)  The routine also prints the
! absolute and relative errors, respectively, of the alpha- and beta-
! coefficients.

n=20
DO  im=1,12
  m=im-1
  npm=n+m
  
! Generate the Legendre recurrence coefficients required in the
! routines  indp  and dindp.
  
  CALL recur(npm,1,0.,0.,a,b,ierr)
  CALL drecur(npm,1,0.d0,0.d0,da,db,ierr)
  
! Compute the desired recursion coefficients.
  
  CALL indp(n,m,a,b,epsma,alpha,beta,ierr,z,w,e,a1,b1)
  CALL dindp(n,m,da,db,depsma,dalpha,dbeta,ierr,dz,dw,da1,db1)
  
! Compute and print the respective errors.
  
  erram(im)=0.
  errbm(im)=0.
  DO  k=1,n
    erra=ABS(alpha(k)-dalpha(k))
    errb=ABS((beta(k)-dbeta(k))/dbeta(k))
    IF(erra > erram(im)) erram(im)=erra
    IF(errb > errbm(im)) errbm(im)=errb
    betap(k,im)=dbeta(k)
  END DO
END DO
DO  ip=1,3
  ip4=1+4*(ip-1)
  WRITE(*,1) ip4-1, ip4, ip4+1, ip4+2
  1 FORMAT('     k  m=', i1, '  beta(k)  m=', i1, '  beta(k)  m=', i2,  &
           ' beta(k)  m=', i2, ' beta(k)'/)
  DO  k=1,n
    km1=k-1
    WRITE(*,2) km1, betap(k,ip4), betap(k,ip4+1), betap(k,ip4+2), betap(k,ip4+3)
    2 FORMAT(' ', i5, 4F14.10)
  END DO
  WRITE(*,3) erram(ip4),erram(ip4+1),erram(ip4+2),erram(ip4+3)
  3 FORMAT(/'    erra', e12.4, 3E14.4)
  WRITE(*,4) errbm(ip4), errbm(ip4+1), errbm(ip4+2), errbm(ip4+3)
  4 FORMAT('    errb', e12.4, 3E14.4//)
END DO
STOP


CONTAINS


SUBROUTINE indp(n,m,a,b,eps,alpha,beta,ierr,z,w,e,a1,b1)

! If  p(m)(.)  denotes the (monic) orthogonal polynomial of degree  m
! relative to the weight function  w(x), then the corresponding m-th
! induced orthogonal polynomials  pind(k,m)(.), k=0,1,2,..., are those
! orthogonal with respect to the weight function

!                   (p(m)(x)**2)*w(x).

! (For background on induced orthogonal polynomials, including an
! algorithm for generating their recursion coefficients, see W. Gautschi
! and S. Li,``A set of orthogonal polynomials induced by a given
! orthogonal polynomial'', Aequationes Math., to appear.) This routine
! obtains the first n recurrence coefficients of the m-th induced
! orthogonal polynomials by an m-fold application of the routine  chri
! with  iopt=7, the shifts taken being, in succession, the zeros of
! p(m)(.).

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: m
REAL, INTENT(IN)      :: a(:)
REAL, INTENT(IN)      :: b(:)
REAL, INTENT(IN OUT)  :: eps
REAL, INTENT(OUT)     :: alpha(:)
REAL, INTENT(OUT)     :: beta(:)
INTEGER, INTENT(OUT)  :: ierr
REAL, INTENT(OUT)     :: z(m)
REAL, INTENT(IN OUT)  :: w(m)
REAL, INTENT(IN OUT)  :: e(m)
REAL, INTENT(OUT)     :: a1(:)
REAL, INTENT(OUT)     :: b1(:)

! The arrays  a,b,alpha,beta,a1,b1  are assumed to have dimension  n+m.

INTEGER  :: ierrc, imu, k, mi, npm
REAL     :: x

npm=n+m
DO  k=1,npm
  alpha(k)=a(k)
  beta(k)=b(k)
END DO
IF(m == 0) RETURN
CALL gauss(m,a,b,eps,z,w,ierr,e)
DO  imu=1,m
  mi=npm-imu
  DO  k=1,mi+1
    a1(k)=alpha(k)
    b1(k)=beta(k)
  END DO
  x=z(imu)
  CALL chri(mi,7,a1,b1,x,0.,0.,0.,alpha,beta,ierrc)
END DO
RETURN
END SUBROUTINE indp



SUBROUTINE dindp(n,m,da,db,deps,dalpha,dbeta,ierr,dz,dw,da1,db1)

! N.B. Argument DE has been removed.

! This is a double-precision version of the routine  indp.

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN)      :: da(:)
REAL (dp), INTENT(IN)      :: db(:)
REAL (dp), INTENT(IN OUT)  :: deps
REAL (dp), INTENT(OUT)     :: dalpha(:)
REAL (dp), INTENT(OUT)     :: dbeta(:)
INTEGER, INTENT(OUT)       :: ierr
REAL (dp), INTENT(OUT)     :: dz(:)
REAL (dp), INTENT(IN OUT)  :: dw(:)
REAL (dp), INTENT(OUT)     :: da1(:)
REAL (dp), INTENT(OUT)     :: db1(:)

REAL (dp)  :: dx
INTEGER    :: ierrc, imu, k, mi, npm

! The arrays  da,db,dalpha,dbeta,da1,db1  are assumed to have dimension  n+m.

npm=n+m
DO  k=1,npm
  dalpha(k)=da(k)
  dbeta(k)=db(k)
END DO
IF(m == 0) RETURN

CALL dgauss(m,da,db,deps,dz,dw,ierr)
DO  imu=1,m
  mi=npm-imu
  DO  k=1,mi+1
    da1(k)=dalpha(k)
    db1(k)=dbeta(k)
  END DO
  dx=dz(imu)
  CALL dchri(mi,7,da1,db1,dx,0.d0,0.d0,0.d0,dalpha,dbeta,ierrc)
END DO
RETURN
END SUBROUTINE dindp



SUBROUTINE chri(n,iopt,a,b,x,y,hr,hi,alpha,beta,ierr)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-03  Time: 20:15:05

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
! polynomials with respect to some measure  dlambda(t), it generates
! the recursion coefficients  alpha(k),beta(k), k=0,1,...,n-1, for the
! measure

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
!               iopt - - an integer selecting the desired weight
!                        distribution
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

! The routine uses the function subroutine  r1mach  to evaluate the
! constant  eps, which is used only if  iopt=7.

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: iopt
REAL, INTENT(IN)      :: a(:)
REAL, INTENT(IN)      :: b(:)
REAL, INTENT(IN)      :: x
REAL, INTENT(IN)      :: y
REAL, INTENT(IN)      :: hr
REAL, INTENT(IN)      :: hi
REAL, INTENT(OUT)     :: alpha(n)
REAL, INTENT(OUT)     :: beta(n)
INTEGER, INTENT(OUT)  :: ierr

! The arrays  a,b  are assumed to have dimension  n+1.

INTEGER  :: k, nm1
REAL     :: c, c0, cm1, d, e, ei, eio, eioo, eo, eoo, eps, er, ero, eroo,  &
            gamma, p2, q, s, so, t, u

eps = 5. * EPSILON(0.0)

! The quantity  eps  is a constant slightly larger than the machine
! precision.

ierr = 0
IF (n < 2) THEN
  ierr = 1
  RETURN
END IF

! What follows implements Eq. (3.7) of W. Gautschi, op. cit.

IF (iopt == 1) THEN
  e = 0.
  DO  k = 1, n
    q = a(k) - e - x
    beta(k) = q * e
    e = b(k+1) / q
    alpha(k) = x + q + e
  END DO
  
! Set the first beta-coefficient as discussed in Section 5.1 of the
! companion paper.
  
  beta(1) = b(1) * (a(1)-x)
  RETURN
  
! What follows implements Eq. (4.7) of W. Gautschi, op. cit.
  
ELSE IF (iopt == 2) THEN
  s = x - a(1)
  t = y
  eio = 0.
  DO  k = 1, n
    d = s * s + t * t
    er = -b(k+1) * s / d
    ei = b(k+1) * t / d
    s = x + er - a(k+1)
    t = y + ei
    alpha(k) = x + t * er / ei - s * ei / t
    beta(k) = t * eio * (1.+(er/ei)**2)
    eio = ei
  END DO
  
! Set the first beta-coefficient.
  
  beta(1) = b(1) * (b(2)+(a(1)-x)**2+y*y)
  RETURN
  
! What follows implements Eq. (4.8) of W. Gautschi, op. cit.
  
ELSE IF (iopt == 3) THEN
  t = y
  eio = 0.
  DO  k = 1, n
    ei = b(k+1) / t
    t = y + ei
    alpha(k) = 0.
    beta(k) = t * eio
    eio = ei
  END DO
  
! Set the first beta-coefficient.
  
  beta(1) = b(1) * (b(2)+y*y)
  RETURN
  
! What follows implements Eqs. (5.1),(5.2) of W. Gautschi, op. cit.
  
ELSE IF (iopt == 4) THEN
  alpha(1) = x - b(1) / hr
  beta(1) = -hr
  q = -b(1) / hr
  DO  k = 2, n
    e = a(k-1) - x - q
    beta(k) = q * e
    q = b(k) / e
    alpha(k) = q + e + x
  END DO
  RETURN
  
! What follows implements Eq. (5.8) of W. Gautschi, op. cit.
  
ELSE IF (iopt == 5) THEN
  nm1 = n - 1
  d = hr * hr + hi * hi
  eroo = a(1) - x + b(1) * hr / d
  eioo = -b(1) * hi / d - y
  alpha(1) = x + hr * y / hi
  beta(1) = -hi / y
  alpha(2) = x - b(1) * hi * eroo / (d*eioo) + hr * eioo / hi
  beta(2) = y * eioo * (1.+(hr/hi)**2)
  IF (n == 2) RETURN
  so = b(2) / (eroo**2+eioo**2)
  ero = a(2) - x - so * eroo
  eio = so * eioo - y
  alpha(3) = x + eroo * eio / eioo + so * eioo * ero / eio
  beta(3) = -b(1) * hi * eio * (1.+(eroo/eioo)**2) / d
  IF (n == 3) RETURN
  DO  k = 3, nm1
    s = b(k) / (ero**2+eio**2)
    er = a(k) - x - s * ero
    ei = s * eio - y
    alpha(k+1) = x + ero * ei / eio + s * eio * er / ei
    beta(k+1) = so * eioo * ei * (1.+(ero/eio)**2)
    eroo = ero
    eioo = eio
    ero = er
    eio = ei
    so = s
  END DO
  RETURN
  
! What follows implements Eq. (5.9) of W. Gautschi, op. cit.
  
ELSE IF (iopt == 6) THEN
  nm1 = n - 1
  eoo = -b(1) / hi - y
  eo = b(2) / eoo - y
  alpha(1) = 0.
  beta(1) = -hi / y
  alpha(2) = 0.
  beta(2) = y * eoo
  IF (n == 2) RETURN
  alpha(3) = 0.
  beta(3) = -b(1) * eo / hi
  IF (n == 3) RETURN
  DO  k = 3, nm1
    e = b(k) / eo - y
    beta(k+1) = b(k-1) * e / eoo
    alpha(k+1) = 0.
    eoo = eo
    eo = e
  END DO
  RETURN
  
! What follows implements a QR step with shift  x.
  
ELSE IF (iopt == 7) THEN
  u = 0.
  c = 1.
  c0 = 0.
  DO  k = 1, n
    gamma = a(k) - x - u
    cm1 = c0
    c0 = c
    IF (ABS(c0) > eps) THEN
      p2 = (gamma**2) / c0
    ELSE
      p2 = cm1 * b(k)
    END IF
    IF (k > 1) beta(k) = s * (p2+b(k+1))
    s = b(k+1) / (p2+b(k+1))
    c = p2 / (p2+b(k+1))
    u = s * (gamma+a(k+1)-x)
    alpha(k) = gamma + u + x
  END DO
  beta(1) = b(1) * (b(2)+(x-a(1))**2)
ELSE
  ierr = 2
END IF

RETURN
END SUBROUTINE chri

END PROGRAM test10
