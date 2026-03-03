PROGRAM test7
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50

USE orthpol
USE s_orthpol
IMPLICIT NONE

! EXTERNAL quad,dquad
REAL      :: endl(4),endr(4),xp(1)=0.0,yp(1)=0.0,xfer(100),wfer(100),alpha(40)=0.0, &
             beta(40)=0.0,be(40),x(100),w(100),xm(400),wm(400),p0(400),  &
             p1(400),p2(400),eps,epsma,erra, errb,fi
REAL (dp) :: depsma, dendl(4), dendr(4), dxp(1)=0.0_dp, dyp(1)=0.0_dp, &
             dalpha(40)=0.0, dbeta(40)=0.0,dbe(40), dx(250), dw(250), di, deps
LOGICAL   :: finl,finr,finld,finrd
INTEGER   :: i, idelta, ie, ied, ierr, ierrd, iq, irout, k, km1, kount=0,  &
             kountd=0, mc, mcd, mp, n, ncap=0, ncapd=0, ncapm, ncpmd

INTERFACE
  SUBROUTINE quad(n,x,w,i,ierr)
    USE common_test6
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: n
    REAL, INTENT(IN OUT)  :: x(:)
    REAL, INTENT(IN OUT)  :: w(:)
    INTEGER, INTENT(IN)   :: i
    INTEGER, INTENT(OUT)  :: ierr
  END SUBROUTINE quad

  SUBROUTINE dquad(n, dx, dw, i, ierr)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)    :: n
    REAL (dp), INTENT(IN)  :: dx(:)
    REAL (dp), INTENT(IN)  :: dw(:)
    INTEGER, INTENT(IN)    :: i
    INTEGER, INTENT(IN)    :: ierr
  END SUBROUTINE dquad

  FUNCTION dwf(dx, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: dx
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION dwf
END INTERFACE

! This test generates in single and REAL (dp) the first 40
! recurrence coefficients of the orthogonal polynomials for the half-
! range Hermite weight function

!                   exp(-x**2)  on  (0,oo).

! Printed are the double-precision values of the alpha- and beta-
! coefficients along with the respective relative errors of the single-
! precision values.

finl=.true.
finr=.false.
finld=.true.
finrd=.false.
epsma=EPSILON(0.0)
depsma=EPSILON(0.0_dp)

! epsma and depsma are the machine single and REAL (dp).

iq=2
idelta=1
irout=1
n=40
mc=4
mcd=4
mp=0
ncapm=100
ncpmd=250

! Set up the partition for the discretization of the inner product.

DO  i=1,4
  fi=i
  di=fi
  endl(i)=3.*(fi-1.)
  endr(i)=3.*fi
  dendl(i)=3.d0*(di-1.d0)
  dendr(i)=3.d0*di
END DO
eps=50.*epsma
deps=1000.*depsma

! Compute the desired recursion coefficients by the multiple-component
! discretization procedure. On the third and fourth subinterval of the
! partition, the quadrature weights produced by  qgp  resp.  dqgp  may
! underflow, on the fourth subinterval even on machines with large
! exponent range.

CALL mcdis(n,ncapm,mc,mp,xp,yp,quad,eps,iq,idelta,irout,finl,finr,endl,endr, &
           xfer,wfer,alpha,beta,ncap,kount,ierr,ie,be,x,w,xm,wm,p0,p1,p2)
CALL dmcdis(n,ncpmd,mcd,mp,dxp,dyp,dquad,deps,iq,idelta,irout,finld,finrd,  &
            dendl,dendr,dalpha,dbeta,ncapd,kountd, ierrd,ied,dbe,dx,dw,dwf)
WRITE(*,1) ncap,kount,ierr,ie
1 FORMAT(/' ncap = ',i3,' kount = ',i2,' ierr = ',i3,' ie = ',i3)
WRITE(*,2) ncapd,kountd,ierrd,ied
2 FORMAT(' ncapd =',i3,' kountd =',i2,' ierrd =',i3,' ied =',i3/)

! Print the results.

WRITE(*,3)
3 FORMAT(/t6, 'k         dalpha(k)                dbeta(k)')
WRITE(*,4)
4 FORMAT(t12, 'erra               errb')
DO  k=1,n
  km1=k-1
  erra=ABS((alpha(k) - dalpha(k))/dalpha(k))
  errb=ABS((beta(k) - dbeta(k))/dbeta(k))
  WRITE(*,5) km1, dalpha(k), dbeta(k)
  5 FORMAT(' ', i5, 2g24.16)
  WRITE(*,6) erra,errb
  6 FORMAT(t7, e12.4, t40, e12.4)
END DO

STOP

END PROGRAM test7



SUBROUTINE quad(n, x, w, i, ierr)
IMPLICIT NONE

INTEGER, INTENT(IN)   :: n
REAL, INTENT(IN OUT)  :: x(:)
REAL, INTENT(IN OUT)  :: w(:)
INTEGER, INTENT(IN)   :: i
INTEGER, INTENT(OUT)  :: ierr

WRITE(*, *) ' User has selected the wrong SP-quadrature routine.'
ierr = 1
STOP
END SUBROUTINE quad



SUBROUTINE dquad(n, dx, dw, i, ierr)
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: dx(:)
REAL (dp), INTENT(IN)  :: dw(:)
INTEGER, INTENT(IN)    :: i
INTEGER, INTENT(IN)    :: ierr

WRITE(*, *) ' User has selected the wrong DP-quadrature routine.'
STOP
END SUBROUTINE dquad



FUNCTION wf(x,i) RESULT(fn_val)
IMPLICIT NONE
REAL, INTENT(IN)     :: x
INTEGER, INTENT(IN)  :: i
REAL                 :: fn_val

fn_val=EXP(-x*x)
RETURN
END FUNCTION wf



FUNCTION dwf(dx, i) RESULT(fn_val)
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
REAL (dp), INTENT(IN)  :: dx
INTEGER, INTENT(IN)    :: i
REAL (dp)              :: fn_val

fn_val = EXP(-dx*dx)

RETURN
END FUNCTION dwf
