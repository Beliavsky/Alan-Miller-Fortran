MODULE common_test5
USE orthpol, ONLY: dp
IMPLICIT NONE

! COMMON/s/ a, b, e, epsma

REAL (dp), SAVE :: a(41), b(41), e(41), epsma

END MODULE common_test5



PROGRAM test5
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50

! This test generates the first 40 recursion coefficients for
! polynomials orthogonal with respect to the Jacobi weight function
! with parameters  alj = -.8(.2)1., bej = -.8(.2)1.  and an added mass
! point of strength  y = .5, 1, 2, 4, 8  at the left end point.  It
! also computes the maximum relative errors (absolute errors for alpha-
! coefficients near zero) of the computed coefficients by comparing
! them against the exact coefficients known analytically.

USE orthpol
USE common_test5
IMPLICIT NONE

! EXTERNAL qjac
INTERFACE
  SUBROUTINE qjac(n, x, w, i, ierr)
    IMPLICIT NONE
    INTEGER, PARAMETER      :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: x(:)
    REAL (dp), INTENT(OUT)  :: w(:)
    INTEGER, INTENT(IN)     :: i
    INTEGER, INTENT(OUT)    :: ierr
  END SUBROUTINE qjac

  FUNCTION dwf(dx, i) RESULT(fn_val)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    REAL (dp), INTENT(IN)  :: dx
    INTEGER, INTENT(IN)    :: i
    REAL (dp)              :: fn_val
  END FUNCTION dwf
END INTERFACE

REAL (dp)  :: xp(1), yp(1), endl(1), endr(1), alpha(40), beta(40), be(40), &
              x(41), w(41), eps, erra, errb, erram, alpham, errbm, betam
REAL (dp)  :: dy, dalj, dbej, dalpbe, daex(40), dbex(40), dnum, dd, dkm1, dden
LOGICAL    :: finl, finr
INTEGER    :: ia, ib, idelta, iderr, ie, ierr, iq, irout, iy, k, kam, kbm, &
              km1, kount, mc, mp, n, ncap, ncapm

! COMMON/s/a,b,e,epsma

WRITE(*,*)

! epsma is the machine double precision

epsma = EPSILON(0.0_dp)
iq = 1
irout = 1
idelta = 2
n = 40
mc = 1
mp = 1
xp(1) = -1.
ncapm = 41
eps = 5000.*epsma
dy = .25_dp

! finl, finr, endl & endr are not needed when iq = 1
finl = .TRUE.
finr = .TRUE.
endl(1) = 0.0_dp
endr = 0.0_dp

DO  iy=1,5
  dy = 2.*dy
  yp(1) = dy
  WRITE(*,2) dy
  2 FORMAT(/' y = ', f6.2/)
  WRITE(*,3)
  3 FORMAT('  alj  bej     erra       errb      alpha',  &
           '    beta    ka  kb ierr ie it'/)
  DO  ia=1,10
    dalj = -1.0_dp + .2_dp*ia
    DO  ib=1,10
      dbej = -1.0_dp + .2_dp*ib
      dalpbe = dalj + dbej
      
! Generate the Jacobi recurrence coefficients.
      
!      CALL recur(ncapm,6,alj,bej,a,b,ierr)
      CALL drecur(ncapm, 6, dalj, dbej, a, b, iderr)
      
! Compute the desired recursion coefficients.
      
      CALL dmcdis(n, ncapm, mc, mp, xp, yp, qjac, eps, iq, idelta, irout,  &
                  finl, finr, endl, endr, alpha, beta, ncap,  &
                  kount, ierr, ie, be, x, w, dwf)
      
! Compute the exact coefficients by Eqs. (4.19)-(4.21) of the companion
! paper along with the relative errors (absolute errors for alpha-
! coefficients close to zero).
      
      daex(1) = (a(1) - dy)/(1.d0 + dy)
      dbex(1) = 1.d0 + dy
      erra = ABS(alpha(1) - daex(1))
      IF(ABS(daex(1)) > eps) erra = erra/ABS(daex(1))
      errb = ABS((beta(1) - dbex(1)) / dbex(1))
      erram = erra
      alpham = alpha(1)
      errbm = errb
      betam = beta(1)
      kam = 0
      kbm = 0
      dnum = 1.d0 + dy
      dd = 1.d0
      DO  k=2,n
        km1 = k-1
        dkm1 = DBLE(km1)
        dden = dnum
        IF(k > 2) dd = (dbej + dkm1)*(dalpbe + dkm1)*dd/((dalj + dkm1-1.d0)  &
                       *(dkm1 - 1.d0))
        dnum = (1.d0 + (dbej + dkm1 + 1.d0)*(dalpbe + dkm1 + 1.d0)*dy*dd/  &
               (dkm1*(dalj + dkm1)))/(1.d0 + dy*dd)
        daex(k) = a(k) + 2.d0*dkm1*(dkm1 + dalj)*(dnum - 1.d0)/  &
                  ((dalpbe + 2.d0*dkm1)*(dalpbe + 2.d0*dkm1 + 1.d0))  &
                  + 2.d0*(dbej + dkm1 + 1.d0)*(dalpbe + dkm1 + 1.d0)*((1.d0/dnum)  &
                  - 1.d0)/((dalpbe + 2.d0*dkm1 + 1.d0)*(dalpbe + 2.d0*dkm1 + 2.d0))
        dbex(k) = dnum*b(k)/dden
        erra = ABS(alpha(k) - daex(k))
        IF(ABS(daex(k)) > eps) erra = erra/ABS(daex(k))
        errb = ABS((beta(k) - dbex(k)) / dbex(k))
        IF(erra > erram) THEN
          erram = erra
          alpham = alpha(k)
          kam = km1
        END IF
        IF(errb > errbm) THEN
          errbm = errb
          betam = beta(k)
          kbm = km1
        END IF
      END DO
      
! Print the results.
      
      WRITE(*,4) dalj, dbej, erram, errbm, alpham, betam, kam, kbm, ierr,  &
                 ie, kount
      4 FORMAT(' ', 2F5.2, 2E11.4, 2F9.6, 4I4, i2)
    END DO
    WRITE(*,*)
  END DO
END DO

STOP
END PROGRAM test5



SUBROUTINE qjac(n, x, w, i, ierr)
USE orthpol
USE common_test5
IMPLICIT NONE

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: x(:)
REAL (dp), INTENT(OUT)  :: w(:)
INTEGER, INTENT(IN)     :: i
INTEGER, INTENT(OUT)    :: ierr

! DIMENSION  a(41),b(41),e(41)
! COMMON/s/a,b,e,epsma

CALL dgauss(n, a, b, epsma, x, w, ierr)
w(1:n) = w(1:n) / b(1)

RETURN
END SUBROUTINE qjac



FUNCTION dwf(x, i) RESULT(fn_val)
USE orthpol, ONLY: dp
IMPLICIT NONE

REAL (dp), INTENT(IN)  :: x
INTEGER, INTENT(IN)    :: i
REAL (dp)              :: fn_val

! This is a dummy function. It is never called, since the routine
! qgp  in  mcdis  which requires  wf  is not activated in this test.

fn_val = 0.0_dp

RETURN
END FUNCTION dwf
