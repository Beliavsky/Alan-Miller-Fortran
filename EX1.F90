PROGRAM ex1
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-12  Time: 20:30:16

!       ROSENBROCK'S BANANA

USE Equation_Solver
IMPLICIT NONE

INTEGER, PARAMETER  :: mgll = 10, n = 2, nunit = 10
REAL (dp)  :: alpha, confac, delfac, delta, epsmch, etafac, fcnnew, fdtolj, &
              ftol, jac(n,n), lam0, mstpf, nsttol, omega, ratiof, sigma,   &
              stptol
INTEGER    :: acptcr, itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns,  &
              maxqns, minqns, narmij, niejev, njacch, njetot, output,   &
              qnupdm, stopcr, supprs, trmcod, trupdm
REAL (dp)  :: a(n,n), boundl(n), boundu(n), delf(n), fsave(n),  &
              ftrack(0:mgll-1), fvec(n), fvecc(n), h(n,n), hhpi(n),  &
              plee(n,n), rdiag(n), s(n), sbar(n), scalef(n), scalex(n),  &
              sn(n), ssdhat(n), strack(0:mgll-1), vhat(n), wv1(n), wv2(n),  &
              wv3(n), wv4(n), xc(n), xplus(n), xsave(n)
LOGICAL    :: absnew, cauchy, deuflh, geoms, linesr, newton, overch
CHARACTER (LEN=6) :: help

! EXTERNAL fcn, jacob
! COMMON /nnes_1/ matsup
! COMMON /nnes_2/ wrnsup
! COMMON /nnes_3/ bypass
! COMMON /nnes_4/ nfetot

INTERFACE
  SUBROUTINE fcn (overfl, n, fvec, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: fvec(n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE fcn

  SUBROUTINE jacob (overfl, n, jac, xc)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    LOGICAL, INTENT(OUT)    :: overfl
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: jac(n,n)
    REAL (dp), INTENT(IN)   :: xc(n)
  END SUBROUTINE jacob
END INTERFACE

OPEN (UNIT=nunit, FILE='EX1.OUT', STATUS='UNKNOWN')
xc(1) = -1.2D0
xc(2) = 1.0D0
CALL setup(absnew, cauchy, deuflh, geoms, linesr, newton, overch, acptcr,  &
           itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns, maxqns,   &
           minqns, n, narmij, niejev, njacch, output, qnupdm, stopcr, supprs, &
           trupdm, alpha, confac, delta, delfac, epsmch, etafac, fdtolj,    &
           ftol, lam0, mstpf, nsttol, omega, ratiof, sigma, stptol, boundl, &
           boundu, scalef, scalex, help)
CALL nnes(absnew, cauchy, deuflh, geoms, linesr, newton, overch, acptcr,  &
          itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns, maxqns, mgll,  &
          minqns, n, narmij, niejev, njacch, njetot, nunit, output, qnupdm,  &
          stopcr, supprs, trmcod, trupdm, alpha, confac, delta, delfac,   &
          epsmch, etafac, fcnnew, fdtolj, ftol, lam0, mstpf, nsttol, omega,  &
          ratiof, sigma, stptol, a, boundl, boundu, delf, fsave, ftrack,   &
          fvec, fvecc, h, hhpi, jac, plee, rdiag, s, sbar, scalef, scalex,  &
          sn, ssdhat, strack, vhat, wv1, wv2, wv3, wv4, xc, xplus, xsave,   &
          help, fcn, jacob)
STOP
END PROGRAM ex1



SUBROUTINE fcn(overfl, n, fvec, xc)

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: fvec(n)
REAL (dp), INTENT(IN)   :: xc(n)

! COMMON /nnes_4/ nfetot

overfl = .false.
! nfetot = nfetot + 1

!       Rosenbrock Banana Function

!       f1 = 10(x2 - x1**2)
!       f2 = 1 - x1

!       start:  (-1.2,  1)
!               (6.39, -0.221)

!       sol'n:  (1,1)

fvec(1) = 10.0D0 * (xc(2) - xc(1)**2)
fvec(2) = 1.0D0 - xc(1)
RETURN
END SUBROUTINE fcn



SUBROUTINE jacob(overfl, n, jac, xc)

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: jac(n,n)
REAL (dp), INTENT(IN)   :: xc(n)

overfl = .false.
jac = 0.0_dp
RETURN
END SUBROUTINE jacob
