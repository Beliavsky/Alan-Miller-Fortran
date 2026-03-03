PROGRAM test3
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-02-04  Time: 20:39:50
! Latest revision - 7 June 2001

! This test applies both the Stieltjes procedure (cf. Section 2.1 of
! W. Gautschi, ``On generating orthogonal polynomials'', SIAM J. Sci.
! Statist. Comput. 3, 1982, 289-317) and the Lanczos algorithm (cf.
! W.B. Gragg and W.J. Harrod, ``The numerically stable reconstruction of
! Jacobi matrices from spectral data'', Numer. Math. 44, 1984, 317-335)
! to generate the  N  recursion coefficients, N = 40, 80, 160, 320,
! for the (monic) orthogonal polynomials relative to the discrete inner
! product supported on  N  equally spaced points on [-1,1] (including
! the end points) and having equal weights 2/N.  The routine prints the
! absolute errors in the alpha-coefficients and the relative errors in
! the beta-coefficients, these being computed using known formulae for
! the coefficients.  The maxima of these errors are also printed.

USE orthpol
IMPLICIT NONE

REAL (dp) :: x(320), w(320), bexact(320), als(320), bes(320), all(320), &
             bel(320), erral, erralm, erras, errbl, errblm, errbs, fkm1, &
             fncap, fncm1
INTEGER   :: ierrl, ierrs, in, incap, inm1, k, ncap

ncap = 20
DO  incap=1,4
  ncap = 2*ncap
  fncap = ncap
  fncm1 = ncap-1
  
! Generate the abscissae and weights of the discrete inner product and
! the exact beta-coefficients.
  
  x(1) = -1.
  w(1) = 2./fncap
  bexact(1) = 2.
  DO  k=2,ncap
    fkm1 = k-1
    x(k) = -1. + 2.*fkm1/fncm1
    w(k) = w(1)
    bexact(k) = (1. + 1./fncm1)**2*(1. - (fkm1/fncap)**2)/ (4. - (1./fkm1)**2)
  END DO
  
! Compute the desired coefficients, first by the Stieltjes procedure,
! and then by the Lanczos algorithm.  Indicate via the error flag  ierrs
! whether a critical underflow condition has arisen in Stieltjes's
! procedure.  (There may, in addition, occur harmless underflow, which
! the routine  sti  does not test for.)
  
  CALL dsti(ncap, ncap, x, w, als, bes, ierrs)
  CALL dlancz(ncap, ncap, x, w, all, bel, ierrl)
  WRITE(*,1) ierrs, ierrl
  1 FORMAT(/'     ierr in sti = ', i4, '           ierr in lancz = ', i3/)
  
! Compute and print the absolute errors of the alpha-coefficients and
! the relative errors of the beta-coefficients as well as the maximum
! respective errors.
  
  erralm = 0.
  errblm = 0.
  WRITE(*,2)
  2 FORMAT(t6, 'k    erra        errb          erra        errb'/)
  DO  in=1,ncap
    inm1 = in-1
    erras = ABS(als(in))
    errbs = ABS((bes(in) - bexact(in))/bexact(in))
    erral = ABS(all(in))
    errbl = ABS((bel(in) - bexact(in))/bexact(in))
    IF(erral > erralm) erralm = erral
    IF(errbl > errblm) errblm = errbl
    IF(ierrs == 0 .OR. inm1 < ABS(ierrs)) THEN
      IF(in == 1) THEN
        WRITE(*,3) inm1, erras, errbs, erral, errbl, ncap
        3 FORMAT(' ', i5, 2E12.4, '  ', 2E12.4, '   N =', i4)
      ELSE
        WRITE(*,4) inm1, erras, errbs, erral, errbl
        4 FORMAT(' ', i5, 2E12.4, '  ', 2E12.4)
      END IF
    ELSE
      IF(in == 1) THEN
        WRITE(*,5) inm1, erral, errbl, ncap
        5 FORMAT(' ', i5, t34, 2E12.4, '   N =', i4)
      ELSE
        WRITE(*,6) inm1, erral, errbl
        6 FORMAT(' ', i5, t34, 2E12.4)
      END IF
    END IF
  END DO
  WRITE(*,7) erralm, errblm
  7 FORMAT(/t33, 2E12.4//)
END DO
STOP
END PROGRAM test3
