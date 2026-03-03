MODULE common_test6
USE orthpol, ONLY: dp
IMPLICIT NONE

! COMMON/s/ a, b, e, epsma

REAL (dp), SAVE :: da(100), db(100), de(100), depsma
REAL, SAVE      :: a(100), b(100), e(100), epsma

END MODULE common_test6



PROGRAM test6
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50

! This test generates in single and REAL (dp) the first 40
! recursion coefficients of the orthogonal polynomials belonging to
! the logistic density function

!              exp(-x)/((1 + exp(-x))**2)  on (-oo, oo).

! It prints the double-precision beta-coefficients (the alpha's being
! all zero) along with the absolute and relative errors of the alpha-
! coefficients resp. beta-coefficients.

USE orthpol
USE s_orthpol
USE common_test6
IMPLICIT NONE

! EXTERNAL qlag,dqlag

INTERFACE
  SUBROUTINE qlag(n,x,w,i,ierr)
    USE common_test6
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: n
    REAL, INTENT(IN OUT)  :: x(:)
    REAL, INTENT(IN OUT)  :: w(:)
    INTEGER, INTENT(IN)   :: i
    INTEGER, INTENT(OUT)  :: ierr
  END SUBROUTINE qlag

  SUBROUTINE dqlag(n,dx,dw,i,ierr)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN OUT)  :: dx(:)
    REAL (dp), INTENT(IN OUT)  :: dw(:)
    INTEGER, INTENT(IN)        :: i
    INTEGER, INTENT(OUT)       :: ierr
  END SUBROUTINE dqlag

  FUNCTION wf(x, i) RESULT(fn_val)
    IMPLICIT NONE
    REAL, INTENT(IN)     :: x
    INTEGER, INTENT(IN)  :: i
    REAL                 :: fn_val
  END FUNCTION wf

  FUNCTION dwf(dx, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: dx
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION dwf
END INTERFACE

REAL      :: xp(1) = 0.0, yp(1) = 0.0, endl(1) = 0.0, endr(1) = 0.0, &
             xfer(1), wfer(1), alpha(40), beta(40), be(40), x(100), w(100), &
             xm(200), wm(200), p0(200), p1(200), p2(200), eps, erra, errb
REAL (dp) :: dxp(1) = 0.0_dp, dyp(1) = 0.0_dp, dendl(1) = 0.0_dp,  &
             dendr(1) = 0.0_dp, dalpha(40), dbeta(40), dbe(40), dx(400), &
             dw(400), deps
LOGICAL   :: finl = .TRUE., finr = .TRUE., finld = .TRUE., finrd = .TRUE.
INTEGER   :: idelta, ie, ied, ierr, ierrd, iq, irout, k, km1, kount, kountd, &
             mc, mp, n, ncap, ncapd, ncapm, ncpmd

! COMMON/s/a,b,e,epsma
! COMMON/d/da,db,de,depsma

WRITE(*,1)
1 FORMAT(/)

! epsma and depsma are the machine single and REAL (dp).

iq=1
idelta=1
irout=1
epsma=EPSILON(0.0)
depsma=EPSILON(0.0_dp)
n=40
mc=2
mp=0
ncapm=100
ncpmd=400
eps=5000.*epsma
deps=1000.*depsma

! Compute the desired coefficients.  On machines with limited exponent range,
! some of the weights in the Gauss-Laguerre quadrature rule may underflow.

CALL mcdis(n,ncapm,mc,mp,xp,yp,qlag,eps,iq,idelta,irout,finl,finr,endl,endr, &
           xfer,wfer,alpha,beta,ncap,kount,ierr,ie,be,x,w,xm,wm,p0,p1,p2)
CALL dmcdis(n,ncpmd,mc,mp,dxp,dyp,dqlag,deps,iq,idelta,irout,finld,finrd,  &
            dendl,dendr,dalpha,dbeta,ncapd,kountd, ierrd,ied,dbe,dx,dw,dwf)
WRITE(*,2) ncap, kount, ierr, ie
2 FORMAT(/' ncap = ', i3, ' kount = ', i2, ' ierr = ', i3, ' ie = ', i5)
WRITE(*,3) ncapd, kountd, ierrd, ied
3 FORMAT(' ncapd= ', i3, ' kountd= ', i2, ' ierrd= ', i3, ' ied= ', i5/)

! Print the results.

WRITE(*,4)
4 FORMAT(/t6, 'k        dbeta(k)             erra        errb'/)
DO  k=1,n
  km1=k-1
  erra=ABS(alpha(k))
  errb=ABS((beta(k)-dbeta(k))/dbeta(k))
  IF(ied == 0 .OR. km1 < ied) THEN
    IF(ie == 0 .OR. km1 < ie) THEN
      WRITE(*,5) km1,dbeta(k),erra,errb
      5 FORMAT(' ', i5, g24.16, 2E12.4)
    ELSE
      WRITE(*,6) km1,dbeta(k)
      6 FORMAT(' ', i5, g24.16)
    END IF
  END IF
END DO

STOP

END PROGRAM test6




SUBROUTINE dqlag(n,dx,dw,i,ierr)
USE common_test6
USE orthpol
IMPLICIT NONE

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: dx(:)
REAL (dp), INTENT(IN OUT)  :: dw(:)
INTEGER, INTENT(IN)        :: i
INTEGER, INTENT(OUT)       :: ierr

! REAL (dp) :: da(300),db(300),de(300),depsma
! COMMON/d/da,db,de,depsma

INTEGER  :: k

CALL drecur(n,7,0.d0,0.d0,da,db,ierr)
CALL dgauss(n,da,db,depsma,dx,dw,ierr)
DO  k=1,n
  dw(k)=dw(k)/((1.d0 + EXP(-dx(k)))**2)
  IF(i == 1) dx(k)=-dx(k)
END DO
RETURN
END SUBROUTINE dqlag



FUNCTION wf(x,i) RESULT(fn_val)

! This is a dummy function. It is never called, since the routine
! qgp  in  mcdis  which requires  wf  is not activated in this test.

IMPLICIT NONE
REAL, INTENT(IN)     :: x
INTEGER, INTENT(IN)  :: i
REAL                 :: fn_val

fn_val=0.
RETURN
END FUNCTION wf



FUNCTION dwf(dx,i) RESULT(fn_val)

! This is a dummy function. It is never called, since the routine
! dqgp  in  dmcdis  which requires  dwf  is not activated in this test.

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: dx
INTEGER, INTENT(IN)    :: i
REAL (dp)              :: fn_val

fn_val=0.0_dp
RETURN
END FUNCTION dwf



SUBROUTINE qlag(n,x,w,i,ierr)
USE common_test6
USE s_orthpol
IMPLICIT NONE

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: x(:)
REAL, INTENT(IN OUT)  :: w(:)
INTEGER, INTENT(IN)   :: i
INTEGER, INTENT(OUT)  :: ierr

! DIMENSION  a(100),b(100),e(100)
! COMMON/s/a,b,e,epsma

INTEGER  :: k

CALL recur(n,7,0.,0.,a,b,ierr)
CALL gauss(n,a,b,epsma,x,w,ierr,e)
DO  k=1,n
  w(k)=w(k)/((1.0 + EXP(-x(k)))**2)
  IF(i == 1) x(k)=-x(k)
END DO
RETURN
END SUBROUTINE qlag

