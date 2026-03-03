PROGRAM ex2
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-12  Time: 20:30:20

USE Equation_Solver
IMPLICIT NONE

INTEGER, PARAMETER  :: mgll = 10, n = 6, nunit = 10
REAL (dp)  :: jac(n,n), lam0, mstpf, nsttol
INTEGER    :: acptcr, itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns,  &
              maxqns, minqns, narmij, niejev, njacch, njetot, output,   &
              qnupdm, stopcr, supprs, trmcod, trupdm
REAL (dp)  :: a(n,n), alpha, boundl(n), boundu(n), confac, delf(n), delfac, &
              delta, epsmch, etafac, fcnnew, fdtolj, fsave(n), ftol,  &
              ftrack(0:mgll-1), fvec(n), fvecc(n), h(n,n), hhpi(n),  &
              omega, plee(n,n), ratiof, rdiag(n), s(n), sbar(n), scalef(n), &
              scalex(n), sigma, sn(n), ssdhat(n), stptol, strack(0:mgll-1), &
              vhat(n), xc(n), xplus(n), xsave(n)
LOGICAL    :: absnew, cauchy, deuflh, geoms, linesr, newton, overch
CHARACTER (LEN=6)  :: help

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

OPEN (UNIT=nunit, FILE='EX2.OUT', STATUS='UNKNOWN')
xc(1:n) = 1.0D0
CALL setup(absnew, cauchy, deuflh, geoms, linesr, newton, overch, acptcr,  &
           itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns, maxqns,   &
           minqns, n, narmij, niejev, njacch, output, qnupdm, stopcr,  &
           supprs, trupdm, alpha, confac, delta, delfac, epsmch, etafac,  &
           fdtolj, ftol, lam0, mstpf, nsttol, omega, ratiof, sigma, stptol, &
           boundl, boundu, scalef, scalex, help)
output = 5
boundl(3) = 0.0D0
boundl(4) = 0.0D0
boundl(6) = 0.0D0

!       STOPS PREMATURELY UNLESS THIS ARE SET VERY STIFF

nsttol = 1.0D-16
stptol = 1.0D-16

!       THIS PARAMETER AFFECTS ILL-CONDITIONED JACOBIAN MECHANISM

!       THIS WILL WORK, =.FALSE. WILL FAIL

bypass = .true.
CALL nnes(absnew, cauchy, deuflh, geoms, linesr, newton, overch, acptcr,  &
          itsclf, itsclx, jactyp, jupdm, maxexp, maxit, maxns, maxqns, mgll, &
          minqns, n, narmij, niejev, njacch, njetot, nunit, output, qnupdm,  &
          stopcr, supprs, trmcod, trupdm, alpha, confac, delta, delfac,   &
          epsmch, etafac, fcnnew, fdtolj, ftol, lam0, mstpf, nsttol, omega,  &
          ratiof, sigma, stptol, a, boundl, boundu, delf, fsave, ftrack,   &
          fvec, fvecc, h, hhpi, jac, plee, rdiag, s, sbar, scalef, scalex,  &
          sn, ssdhat, strack, vhat, xc, xplus, xsave, help, fcn, jacob)
STOP
END PROGRAM ex2



SUBROUTINE fcn(overfl, n, fvec, xc)

USE Equation_Solver, ONLY: nfetot
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

LOGICAL, INTENT(OUT)    :: overfl
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: fvec(n)
REAL (dp), INTENT(IN)   :: xc(n)

! COMMON /nnes_4/ nfetot

overfl = .false.
nfetot = nfetot + 1

!    Hiebert's 2nd Chemical Engineering Problem

!    source: Hiebert; Sandia Technical Report #SAND80-0181
!            Sandia National Laboratories, Albuquerque, NM (1980)

!    f1 = X1 + X2 + X4 - .001
!    f2 = X5 + X6 -55
!    f3 = X1 + X2 + X3 + 2X5 + X6 - 110.001
!    f4 = X1 - 0.1X2
!    f5 = X1 - 10000 X3 X4
!    f6 = X5 - 5.5e15 X3 X6

!    start:  (0, 0, 0, 0, 0, 0)
!            (1, 1, 1, 1, 1, 1)
!            (1e-4, 1e-3, 0, 1e-4, 55, 1e-4)
!            (10, 10, 10, 10, 10, 10, 10)

!    sol'n   (8.264e-5, 8.264e-4, 9.091e-5, 9.091e-5, 55, 1.1e-10)

fvec(1) = xc(1) + xc(2) + xc(4) - 0.001D0
fvec(2) = xc(5) + xc(6) - 55.0D0
fvec(3) = xc(1) + xc(2) + xc(3) + 2.0D0 * xc(5) + xc(6) - 110.001D0
fvec(4) = xc(1) - 0.1D0 * xc(2)
fvec(5) = xc(1) - 1.0D04 * xc(3) * xc(4)
fvec(6) = xc(5) - 5.5D15 * xc(3) * xc(6)
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
