PROGRAM test11
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:51

USE orthpol
USE s_orthpol
IMPLICIT NONE

COMPLEX   :: rho(80), rold(80), z, e
REAL      :: a(500), b(500), alpha(40), beta(40), alphc(40), betc(40), &
             fnu(80), s(40), s0(80), s1(80), s2(80), alphr(40), betr(40), &
             alphcr(40), betcr(40), acmv, acrmv, agmv, al, amv, bcmv,  &
             bcrmv, be, bgmv, bmv, dacrmv, damv, dbcrmv, dbmv, edacr,  &
             edacrm, edbcr, edbcrm, eps, epsd, epsma, eracr, eracrm, erbcr, &
             erbcrm, erra, errac, errb, errbc, erracm, errag, errbg,  &
             erragm, erram, errbcm, errbgm, errbm, errda, errdam, errdb,  &
             errdbm, fndiv, fndm1, hi, hr, pi, r, theta, x, y
REAL (dp) :: depsma, deps, dal, dbe, dhi, da(800), db(800), dx, dy,  &
             dalpha(40), dbeta(40), dnu(80), drhor(80), drhoi(80),  &
             droldr(80), droldi(80), dhr, dalphc(40), dbetc(40), dalphr(40), &
             dbetr(40), dalcr(40), dbecr(40)
INTEGER   :: ial, ibe, ibemax, ierr, ierrc, ierrcd, ierrd, ierrg, ierrgd, &
             ipoly, ir, ith, ix, k, n, nd, ndiv, ndivm1, ndm1, nm1, &
             nu, nud, nu0, nu0d, nu0dv, nu0v, numax, numaxd
REAL  :: xx(5) = (/ 1.001, 1.01, 1.04, 1.07, 1.1 /)
REAL  :: rr(5) = (/ 1.05, 1.1625, 1.275, 1.3875, 1.5 /)

! This test is to illustrate the dissimilar performance of the routines
! chri  and  gchri  in the case of division of the Jacobi weight
! function  w(t;alj,bej)  with parameters  alj,bej  by either a linear
! divisor  t-x  or a quadratic divisor  (t-x)**2 + y**2 . In either
! case, the parameters selected are  alj=-.8(.4).8, bej=alj(.4).8.  In
! the former case, x = -1.001, -1.01, -1.04, -1.07 and -1.1, whereas in
! the latter case,  x and  y  are chosen to lie, regularly spaced, on
! the upper half of an ellipse with foci at +1 and -1 and sum of the
! semiaxes equal to  rho = 1.05, 1.1625, 1.275, 1.3875 and 1.5. The
! routines are run in both single and REAL (dp) with n=40, the
! results of the latter being used to calculate, and print, the maximum
! absolute and relative error of the single-precision alpha- and beta-
! coefficients, respectively. Also printed are the starting recurrence
! indexes required in the backward recurrence schemes of  gchri,dgchri
! to achieve single- resp. double-precision accuracy. This information
! is contained in the first line of each 3-line block of the output,
! where in the case of quadratic divisors only average values (averaged
! over the upper half of the respective ellipse) are shown. The second
! and third line of each 3-line block display the maximum
! ``reconstruction error'', that is, the maximum errors in the alpha's
! and beta's if the coefficients produced by  gchri,chri  and  dgchri,
! dchri  are fed back to the routines  chri  and  dchri  with  iopt=1
! to recover the original recursion coefficients in single and double
! precision.

WRITE(*,1)
1 FORMAT(/)
epsma=EPSILON(0.0)
depsma=EPSILON(0.0_dp)

! epsma and depsma are the machine single and REAL (dp).

n=40
nm1=n-1
nd=2*n
ndm1=nd-1
numax=500
numaxd=800
eps=10.*epsma
deps=100.d0*depsma
epsd=deps
ipoly=6
DO  ial=1,5
  al=-.8 + .4*(ial-1)
  dal=al
  ibemax=6-ial
  DO  ibe=1,ibemax
    be=al + .4*(ibe-1)
    dbe=be
    WRITE(*,2) al,be
    2 FORMAT(///' al = ', f6.2, '  be = ', f6.2//)
    hi=0.
    dhi=0.d0
    
! Generate the Jacobi recurrence coefficients to be used in the
! backward recurrence algorithm of the routines  gchri  and  dgchri.
    
    CALL recur(numax,ipoly,al,be,a,b,ierr)
    CALL drecur(numaxd,ipoly,dal,dbe,da,db,ierrd)
    WRITE(*,3)
    3 FORMAT(t31, 'gchri', t56, 'chri')
    WRITE(*,4)
    4 FORMAT(t6, 'x     nu0  nud0    erra        errb          erra',  &
             '       errb'/)
    DO  ix=1,5
      x=-xx(ix)
      dx=x
      y=0.
      dy=0.d0
      z=CMPLX(x,y)
      
! Compute the starting index for backward recurrence.
      
      nu0=nu0jac(ndm1,z,eps)
      nu0d=nu0jac(ndm1,z,epsd)
      
! Generate the recurrence coefficients for the Jacobi weight function
! divided by a linear divisor, using the routines  gchri,dgchri.
      
      CALL gchri(n,1,nu0,numax,eps,a,b,x,y,alpha,beta,nu,ierrg,  &
                 ierrc,fnu,rho,rold,s,s0,s1,s2)
      
! On machines with limited single-precision exponent range, the routine
! cheb  used in  gchri  may have generated an underflow exception,
! which however is harmless and can be ignored.
      
      CALL dgchri(n, 1, nu0d, numaxd, deps, da, db, dx, dy, dalpha, dbeta,  &
                  nud, ierrgd, ierrcd, dnu, drhor, drhoi, droldr, droldi)
      IF(ierrg /= 0 .OR. ierrc /= 0 .OR. ierrgd /= 0 .OR.ierrcd  /= 0) THEN
        WRITE(*,5) ierrg,ierrgd,al,be,x
        5 FORMAT(/' ierrg in gchri = ',i4,' ierrg in dgchri = ',  &
            i4,' for al = ',f6.2,' be = ',f6.2,' x = ',f7.4)
        WRITE(*,6) ierrc,ierrcd,al,be,x
        6 FORMAT(' ierrc in gchri = ',i4,' ierrc in dgchri = ',  &
                 i4,'  or al = ',f6.2,' be = ',f6.2,' x = ',f7.4/)
        CYCLE
      END IF
      
! Generate the recurrence coefficients for the Jacobi weight function
! divided by a linear divisor, using the routines  chri,dchri.
      
      hr=rho(1)
      dhr=drhor(1)
      CALL chri(n,4,a,b,x,y,hr,hi,alphc,betc,ierr)
      CALL dchri(n,4,da,db,dx,dy,dhr,dhi,dalphc,dbetc,ierr)
      
! Do the reconstruction.
      
      CALL chri(nm1,1,alpha,beta,x,y,0.,0.,alphr,betr,ierr)
      CALL dchri(nm1,1,dalpha,dbeta,dx,dy,0.d0,0.d0,dalphr,dbetr,ierr)
      CALL chri(nm1,1,alphc,betc,x,y,0.,0.,alphcr,betcr,ierr)
      CALL dchri(nm1,1,dalphc,dbetc,dx,dy,0.d0,0.d0,dalcr,dbecr,ierr)
      
! Compute and print the maximum errors.
      
      erragm=0.
      errbgm=0.
      erracm=0.
      errbcm=0.
      erram=0.
      errbm=0.
      errdam=0.
      errdbm=0.
      eracrm=0.
      erbcrm=0.
      edacrm=0.
      edbcrm=0.
      DO  k=1,n
        errag=ABS(alpha(k)-dalpha(k))
        errbg=ABS((beta(k)-dbeta(k))/dbeta(k))
        errac=ABS(alphc(k)-dalphc(k))
        errbc=ABS((betc(k)-dbetc(k))/dbetc(k))
        IF(k < n) THEN
          erra=ABS(alphr(k)-a(k))
          errb=ABS((betr(k)-b(k))/b(k))
          errda=ABS(dalphr(k)-da(k))
          errdb=ABS((dbetr(k)-db(k))/db(k))
          eracr=ABS(alphcr(k)-a(k))
          erbcr=ABS((betcr(k)-b(k))/b(k))
          edacr=ABS(dalcr(k)-da(k))
          edbcr=ABS((dbecr(k)-db(k))/db(k))
          IF(erra > erram) erram=erra
          IF(errb > errbm) errbm=errb
          IF(errda > errdam) errdam=errda
          IF(errdb > errdbm) errdbm=errdb
          IF(eracr > eracrm) eracrm=eracr
          IF(erbcr > erbcrm) erbcrm=erbcr
          IF(edacr > edacrm) edacrm=edacr
          IF(edbcr > edbcrm) edbcrm=edbcr
        END IF
        IF(errag > erragm) erragm=errag
        IF(errbg > errbgm) errbgm=errbg
        IF(errac > erracm) erracm=errac
        IF(errbc > errbcm) errbcm=errbc
      END DO
      WRITE(*,7) x,nu0,nu0d,erragm,errbgm,erracm,errbcm
      7 FORMAT(/' ',f7.4,2I6,2E12.4, '  ', 2E12.4)
      IF(ix == 1) THEN
        WRITE(*,8) erram,errbm,errdam,errdbm
        8 FORMAT(t12,'reconstr.',2E12.4, '  ', 2E12.4)
        WRITE(*,9) eracrm,erbcrm,edacrm,edbcrm
        9 FORMAT(t13,'errors', '  ', 2E12.4, '  ', 2E12.4)
      ELSE
        WRITE(*,11) erram,errbm,errdam,errdbm
        11 FORMAT(t21,2E12.4, '  ', 2E12.4)
        WRITE(*,11) eracrm,erbcrm,edacrm,edbcrm
      END IF
    END DO
    WRITE(*,1)
    ndiv=20
    ndivm1=ndiv-1
    fndiv=ndiv
    fndm1=ndivm1
    pi=4.*ATAN(1.)
    WRITE(*,3)
    WRITE(*,12)
    12 FORMAT(t5, 'rho    nu0  nud0    erra        errb          erra',  &
              '       errb'/)
    DO  ir=1,5
      r=rr(ir)
      agmv=0.
      bgmv=0.
      acmv=0.
      bcmv=0.
      amv=0.
      bmv=0.
      damv=0.
      dbmv=0.
      acrmv=0.
      bcrmv=0.
      dacrmv=0.
      dbcrmv=0.
      nu0v=0
      nu0dv=0
      DO  ith=1,ndivm1
        
! Generate the points on the ellipse.
        
        theta=pi*REAL(ith)/fndiv
        e=CMPLX(COS(theta),SIN(theta))
        z=.5*(r*e+1./(r*e))
        x=REAL(z)
        y=AIMAG(z)
        dx=x
        dy=y
        
! Compute the starting index for backward recurrence.
        
        nu0=nu0jac(ndm1,z,eps)
        nu0d=nu0jac(ndm1,z,epsd)
        
! Generate the recurrence coefficients for the Jacobi weight function
! divided by a quadratic divisor, using the routines  gchri,dgchri.
        
        CALL gchri(n,2,nu0,numax,eps,a,b,x,y,alpha,beta,nu,ierrg,  &
                   ierrc,fnu,rho,rold,s,s0,s1,s2)
        CALL dgchri(n,2,nu0d,numaxd,deps,da,db,dx,dy,dalpha,dbeta,  &
                    nud,ierrgd,ierrcd,dnu,drhor,drhoi,droldr,droldi)
        IF(ierrg /= 0 .OR.ierrc /= 0 .OR. ierrgd /= 0 .OR. ierrcd  /= 0) THEN
          WRITE(*,5) ierrg,ierrgd,al,be,x
          WRITE(*,6) ierrc,ierrcd,al,be,x
          CYCLE
        END IF
        nu0v=nu0v+nu0
        nu0dv=nu0dv+nu0d
        
! Generate the recurrence coefficients for the Jacobi weight function
! divided by a quadratic divisor, using the routines  chri,dchri.
        
        hr=REAL(rho(1))
        hi=AIMAG(rho(1))
        dhr=drhor(1)
        dhi=drhoi(1)
        CALL chri(n,5,a,b,x,y,hr,hi,alphc,betc,ierr)
        CALL dchri(n,5,da,db,dx,dy,dhr,dhi,dalphc,dbetc,ierr)
        
! Do the reconstruction.
        
        CALL chri(nm1,2,alpha,beta,x,y,0.,0.,alphr,betr,ierr)
        CALL dchri(nm1,2,dalpha,dbeta,dx,dy,0.d0,0.d0,dalphr, dbetr,ierr)
        CALL chri(nm1,2,alphc,betc,x,y,0.,0.,alphcr,betcr,ierr)
        CALL dchri(nm1,2,dalphc,dbetc,dx,dy,0.d0,0.d0,dalcr, dbecr,ierr)
        
! Compute and print the maximum average errors.
        
        erragm=0.
        errbgm=0.
        erracm=0.
        errbcm=0.
        erram=0.
        errbm=0.
        errdam=0.
        errdbm=0.
        eracrm=0.
        erbcrm=0.
        edacrm=0.
        edbcrm=0.
        DO  k=1,n
          errag=ABS(alpha(k)-dalpha(k))
          errbg=ABS((beta(k)-dbeta(k))/dbeta(k))
          errac=ABS(alphc(k)-dalphc(k))
          errbc=ABS((betc(k)-dbetc(k))/dbetc(k))
          IF(k < n) THEN
            erra=ABS(alphr(k)-a(k))
            errb=ABS((betr(k)-b(k))/b(k))
            errda=ABS(dalphr(k)-da(k))
            errdb=ABS((dbetr(k)-db(k))/db(k))
            eracr=ABS(alphcr(k)-a(k))
            erbcr=ABS((betcr(k)-b(k))/b(k))
            edacr=ABS(dalcr(k)-da(k))
            edbcr=ABS((dbecr(k)-db(k))/db(k))
            IF(erra > erram) erram=erra
            IF(errb > errbm) errbm=errb
            IF(errda > errdam) errdam=errda
            IF(errdb > errdbm) errdbm=errdb
            IF(eracr > eracrm) eracrm=eracr
            IF(erbcr > erbcrm) erbcrm=erbcr
            IF(edacr > edacrm) edacrm=edacr
            IF(edbcr > edbcrm) edbcrm=edbcr
          END IF
          IF(errag > erragm) erragm=errag
          IF(errbg > errbgm) errbgm=errbg
          IF(errac > erracm) erracm=errac
          IF(errbc > errbcm) errbcm=errbc
        END DO
        agmv=agmv+erragm
        bgmv=bgmv+errbgm
        acmv=acmv+erracm
        bcmv=bcmv+errbcm
        amv=amv+erram
        bmv=bmv+errbm
        damv=damv+errdam
        dbmv=dbmv+errdbm
        acrmv=acrmv+eracrm
        bcrmv=bcrmv+erbcrm
        dacrmv=dacrmv+edacrm
        dbcrmv=dbcrmv+edbcrm
      END DO
      nu0=nu0v/fndm1
      nu0d=nu0dv/fndm1
      erragm=agmv/fndm1
      errbgm=bgmv/fndm1
      erracm=acmv/fndm1
      errbcm=bcmv/fndm1
      erram=amv/fndm1
      errbm=bmv/fndm1
      errdam=damv/fndm1
      errdbm=dbmv/fndm1
      eracrm=acrmv/fndm1
      erbcrm=bcrmv/fndm1
      edacrm=dacrmv/fndm1
      edbcrm=dbcrmv/fndm1
      IF(ir == 1) THEN
        WRITE(*,7) r,nu0,nu0d,erragm,errbgm,erracm,errbcm
        WRITE(*,8) erram,errbm,errdam,errdbm
        WRITE(*,9) eracrm,erbcrm,edacrm,edbcrm
      ELSE
        WRITE(*,7) r,nu0,nu0d,erragm,errbgm,erracm,errbcm
        WRITE(*,11) erram,errbm,errdam,errdbm
        WRITE(*,11) eracrm,erbcrm,edacrm,edbcrm
      END IF
    END DO
  END DO
END DO
STOP


CONTAINS


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
REAL, INTENT(OUT)     :: alpha(:)
REAL, INTENT(OUT)     :: beta(:)
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



SUBROUTINE gchri(n,iopt,nu0,numax,eps,a,b,x,y,alpha,beta,nu,ierr,  &
                 ierrc,fnu,rho,rold,s,s0,s1,s2)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-07  Time: 08:53:13

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
!                   in case of nonconvergence,  nu  will have the value
!                   numax
!           ierr  - an error flag, where
!                   ierr=0     on normal return
!                   ierr=1     if  iopt  is neither 1 nor 2
!                   ierr=nu0   if  nu0 > numax
!                   ierr=numax if the backward recurrence algorithm does
!                              not converge
!                   ierr=-1    if  n  is not in range
!           ierrc - an error flag inherited from the routine  cheb

! The arrays  fnu,s,s0,s1,s2  are working space. The routine calls
! upon the routines  knum  and  cheb.

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN)      :: iopt
INTEGER, INTENT(IN OUT)  :: nu0
INTEGER, INTENT(IN)      :: numax
REAL, INTENT(IN OUT)     :: eps
REAL, INTENT(IN OUT)     :: a(numax)
REAL, INTENT(IN OUT)     :: b(numax)
REAL, INTENT(IN OUT)     :: x
REAL, INTENT(OUT)        :: y
REAL, INTENT(IN OUT)     :: alpha(n)
REAL, INTENT(IN OUT)     :: beta(n)
INTEGER, INTENT(IN OUT)  :: nu
INTEGER, INTENT(OUT)     :: ierr
INTEGER, INTENT(IN OUT)  :: ierrc
REAL, INTENT(OUT)        :: fnu(:)
COMPLEX, INTENT(IN OUT)  :: rho(:)
COMPLEX, INTENT(IN OUT)  :: rold(:)
REAL, INTENT(IN OUT)     :: s(n)
REAL, INTENT(IN OUT)     :: s0(:)
REAL, INTENT(IN OUT)     :: s1(:)
REAL, INTENT(IN OUT)     :: s2(:)

COMPLEX :: z

! The arrays  fnu,rho,rold,s0,s1,s2  are assumed to have dimension  2*n.

IF (n < 1) THEN
  ierr = -1
  RETURN
END IF
ierr = 0
nd = 2 * n
ndm1 = nd - 1

! Linear divisor

IF (iopt == 1) THEN
  
! Generate the modified moments of  dlambda.
  
  z = CMPLX(x,0.)
  CALL knum(ndm1,nu0,numax,z,eps,a,b,rho,nu,ierr)
  DO  k = 1, nd
    fnu(k) = -REAL(rho(k))
  END DO
  
! Compute the desired recursion coefficients by means of the modified
! Chebyshev algorithm.
  
  CALL cheb(n,a,b,fnu,alpha,beta,s,ierrc,s0,s1,s2)
  RETURN
  
! Quadratic divisor
  
ELSE IF (iopt == 2) THEN
  
! Generate the modified moments of  dlambda.
  
  y = ABS(y)
  z = CMPLX(x,y)
  CALL knum(ndm1,nu0,numax,z,eps,a,b,rho,nu,ierr)
  DO  k = 1, nd
    fnu(k) = -AIMAG(rho(k)) / y
  END DO
  
! Compute the desired recursion coefficients by means of the modified
! Chebyshev algorithm.
  
  CALL cheb(n,a,b,fnu,alpha,beta,s,ierrc,s0,s1,s2)
ELSE
  ierr = 1
END IF

RETURN
END SUBROUTINE gchri



SUBROUTINE cheb(n,a,b,fnu,alpha,beta,s,ierr,s0,s1,s2)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-07  Time: 09:43:27

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
REAL, INTENT(OUT)     :: alpha(n)
REAL, INTENT(OUT)     :: beta(n)
REAL, INTENT(OUT)     :: s(n)
INTEGER, INTENT(OUT)  :: ierr
REAL, INTENT(OUT)     :: s0(:)
REAL, INTENT(OUT)     :: s1(:)
REAL, INTENT(OUT)     :: s2(:)

! The arrays  a,b  are assumed to have dimension  2*n-1, the arrays
! fnu,s0,s1,s2  dimension  2*n.

INTEGER  :: k, l, lk, nd
REAL     :: hhuge, ttiny

nd = 2 * n
ttiny = 10. * TINY(0.0)
hhuge = .1 * HUGE(0.0)
ierr = 0
IF (ABS(fnu(1)) < ttiny) THEN
  ierr = 1
  RETURN
END IF
IF (n < 1) THEN
  ierr = 2
  RETURN
END IF

! Initialization

alpha(1) = a(1) + fnu(2) / fnu(1)
beta(1) = fnu(1)
IF (n == 1) RETURN
s(1) = fnu(1)
DO  l = 1, nd
  s0(l) = 0.
  s1(l) = fnu(l)
END DO

! Continuation

DO  k = 2, n
  lk = nd - k + 1
  DO  l = k, lk
    
! The quantities  s2(l)  for l > k are auxiliary quantities which may
! be zero or may become so small as to underflow, without however
! causing any harm.
    
    s2(l) = s1(l+1) - (alpha(k-1)-a(l)) * s1(l) - beta(k-1) *  &
            s0(l) + b(l) * s1(l-1)
    IF (l == k) s(k) = s2(k)
    
! Check impending underflow or overflow
    
    IF (ABS(s(k)) < ttiny) THEN
      ierr = -(k-1)
      RETURN
    ELSE IF (ABS(s(k)) > hhuge) THEN
      ierr = k - 1
      RETURN
    END IF
  END DO
  
! Compute the alpha- and beta-coefficient
  
  alpha(k) = a(k) + (s2(k+1)/s2(k)) - (s1(k)/s1(k-1))
  beta(k) = s2(k) / s1(k-1)
  DO  l = k, lk
    s0(l) = s1(l)
    s1(l) = s2(l)
  END DO
END DO
RETURN
END SUBROUTINE cheb



SUBROUTINE knum(n,nu0,numax,z,eps,a,b,rho,nu,ierr)

! N.B. Argument ROLD has been removed.
 
! Code converted using TO_F90 by Alan Miller
! Date: 2001-06-07  Time: 09:43:34

! This routine generates

!   rho(k)(z) = integral pi(k)(t)dlambda(t)/(z-t), k=0,1,2,...,n,

! where  pi(k)(t)  is the (monic) k-th degree orthogonal polynomial
! with respect to the measure  dlambda(t), and the integral is extended
! over the support of  dlambda. It is assumed that  z  is a complex
! number outside the smallest interval containing the support of dlambda.
! The quantities  rho(k)(z)  are computed as the first  n+1  members of the
! minimal solution of the basic three-term recurrence relation

!      y(k+1)(z)=(z-a(k))y(k)(z)-b(k)y(k-1)(z), k=0,1,2,...,

! satisfied by the orthogonal polynomials  pi(k)(z).

!   Input:  n  - -  the largest integer  k  for which  rho(k)  is desired
!           nu0  -  an estimate of the starting backward recurrence
!                   index; if no better estimate is known, set
!                   nu0 = 3*n/2; for Jacobi, Laguerre and Hermite
!                   weight functions, estimates of  nu0  are generated
!                   respectively by the routines  nu0jac,nu0lag  and nu0her
!           numax - an integer larger than  n  cutting off backward
!                   recursion in case of nonconvergence; if  nu0
!                   exceeds  numax, then the routine aborts with the
!                   error flag  ierr  set equal to  nu0
!           z - - - the variable in  rho(k)(z); type complex
!           eps - - the relative accuracy to which the  rho(k)  are desired
!           a,b - - arrays of dimension  numax  to be supplied with the
!                   recurrence coefficients  a(k-1), b(k-1), k=1,2,..., numax.

!   Output: rho - - an array of dimension  n+1  containing the results
!                   rho(k)=rho(k-1)(z), k=1,2,...,n+1; type complex
!           nu  - - the starting backward recurrence index that yields
!                   convergence
!           ierr  - an error flag equal to zero on normal return, equal to  nu0
!                   if  nu0 > numax, and equal to  numax in case of
!                   nonconvergence.

! The complex array  rold  of dimension  n+1  is used for working space.

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN OUT)  :: nu0
INTEGER, INTENT(IN)      :: numax
COMPLEX, INTENT(IN OUT)  :: z
REAL, INTENT(IN)         :: eps
REAL, INTENT(IN OUT)     :: a(numax)
REAL, INTENT(IN OUT)     :: b(numax)
COMPLEX, INTENT(OUT)     :: rho(:)
INTEGER, INTENT(OUT)     :: nu
INTEGER, INTENT(OUT)     :: ierr

COMPLEX  :: r, rold(n+1)

! The arrays  rho,rold  are assumed to have dimension  n+1.

INTEGER  :: j, j1, k, np1

ierr = 0
np1 = n + 1
IF (nu0 > numax) THEN
  ierr = nu0
  RETURN
END IF
IF (nu0 < np1) nu0 = np1
nu = nu0 - 5
rho(1:np1) = (0.,0.)

20 nu = nu + 5
IF (nu > numax) THEN
  ierr = numax
  GO TO 60
END IF
rold(1:np1) = rho(1:np1)
r = (0.,0.)
DO  j = 1, nu
  j1 = nu - j + 1
  r = CMPLX(b(j1),0.) / (z-CMPLX(a(j1),0.)-r)
  IF (j1 <= np1) rho(j1) = r
END DO
DO  k = 1, np1
  IF (ABS(rho(k)-rold(k)) > eps*ABS(rho(k))) GO TO 20
END DO

60 IF (n == 0) RETURN
DO  k = 2, np1
  rho(k) = rho(k) * rho(k-1)
END DO

RETURN
END SUBROUTINE knum

END PROGRAM test11
