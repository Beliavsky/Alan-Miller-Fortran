MODULE tron

! ************************************************************************* !
!                                                                           !
!            COPYRIGHT NOTIFICATION                                         !
!                                                                           !
! This program discloses material protectable under copyright laws of       !
! the United States. Permission to copy and modify this software and its    !
! documentation for internal research use is hereby granted, provided       !
! that this notice is retained thereon and on all copies or modifications.  !
! The University of Chicago makes no representations as to the suitability  !
! and operability of this software for any purpose.                         !
! It is provided "as is" without express or implied warranty.               !
!                                                                           !
! Use of this software for commercial purposes is expressly prohibited      !
! without contacting                                                        !
!                                                                           !
!    Jorge J. More'                                                         !
!    Mathematics and Computer Science Division                              !
!    Argonne National Laboratory                                            !
!    9700 S. Cass Ave.                                                      !
!    Argonne, Illinois 60439-4844                                           !
!    e-mail: more@mcs.anl.gov                                               !
!                                                                           !
! Argonne National Laboratory with facilities in the states of              !
! Illinois and Idaho, is owned by The United States Government, and         !
! operated by the University of Chicago under provision of a contract       !
! with the Department of Energy.                                            !
!                                                                           !
! ************************************************************************* !
!                                                                           !
!            ADDITIONAL INFORMATION                                         !
!                                                                           !
! Chih-Jen Lin and Jorge J. More',                                          !
! Newton's method for large bound-constrained optimization problems,        !
! Argonne National Laboratory,                                              !
! Mathematics and Computer Science Division,                                !
! Preprint ANL/MCS-P724-0898,                                               !
! August 1998 (Revised March 1999).                                         !
!                                                                           !
! http://www.mcs.anl.gov/~more/papers/tron.ps.gz                            !
!                                                                           !
! The Fortran 77 code can be downloaded from:                               !
!                                                                           !
! http://www-unix.mcs.anl.gov/~more/tron                                    !
!                                                                           !
! ************************************************************************* !
!                                                                           !
! Last modification (Fortran 77 version): April 29, 1999                    !


IMPLICIT NONE
INTEGER, PARAMETER                     :: dp = SELECTED_REAL_KIND(14, 60)
REAL (dp), PARAMETER, PRIVATE          :: zero = 0.0_dp, one = 1.0_dp
REAL (dp), ALLOCATABLE, SAVE, PRIVATE  :: xc(:), s(:), dsave(:), wa(:)
INTEGER, ALLOCATABLE, SAVE, PRIVATE    :: indfree(:), isave(:)

CONTAINS


SUBROUTINE dtron(n, x, xl, xu, f, g, a, adiag, acol_ptr, arow_ind,  &
                 frtol, fatol, fmin, cgtol, itermax, delta, task,  &
                 b, bdiag, bcol_ptr, brow_ind,   &
                 l, ldiag, lcol_ptr, lrow_ind, iterscg)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32
! Latest revision - 5 July 1999

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: xl(:)
REAL (dp), INTENT(IN)      :: xu(:)
REAL (dp), INTENT(IN OUT)  :: f
REAL (dp), INTENT(IN)      :: g(:)
REAL (dp), INTENT(IN)      :: a(:)
REAL (dp), INTENT(IN)      :: adiag(:)
INTEGER, INTENT(IN)        :: acol_ptr(:)    ! acol_ptr(n+1)
INTEGER, INTENT(IN)        :: arow_ind(:)
REAL (dp), INTENT(IN)      :: frtol
REAL (dp), INTENT(IN)      :: fatol
REAL (dp), INTENT(IN)      :: fmin
REAL (dp), INTENT(IN)      :: cgtol
INTEGER, INTENT(IN)        :: itermax
REAL (dp), INTENT(IN OUT)  :: delta
CHARACTER (LEN=*), INTENT(IN OUT) :: task
REAL (dp), INTENT(OUT)     :: b(:)
REAL (dp), INTENT(OUT)     :: bdiag(:)
INTEGER, INTENT(OUT)       :: bcol_ptr(:)    ! bcol_ptr(n+1)
INTEGER, INTENT(OUT)       :: brow_ind(:)
REAL (dp), INTENT(OUT)     :: l(:)
REAL (dp), INTENT(OUT)     :: ldiag(:)
INTEGER, INTENT(OUT)       :: lcol_ptr(:)    ! lcol_ptr(n+1)
INTEGER, INTENT(OUT)       :: lrow_ind(:)
INTEGER, INTENT(OUT)       :: iterscg

!  *********

!  Subroutine dtron

!  This subroutine implements a trust region Newton method for the
!  solution of large bound-constrained optimization problems

!        min { f(x) : xl <= x <= xu }

!  where the Hessian matrix is sparse. The user must evaluate the
!  function, gradient, and the Hessian matrix.

!  This subroutine uses reverse communication.
!  The user must choose an initial approximation x to the minimizer,
!  and make an initial call with task set to 'START'.
!  On exit task indicates the required action.

!  A typical invocation has the following outline:

!  Compute a starting vector x.
!  Compute the sparsity pattern of the Hessian matrix and
!  store in compressed column storage in (acol_ptr,arow_ind).

!  task = 'START'
!  do while (search)

!     if (task .eq. 'F' .or. task .eq. 'START') then
!        Evaluate the function at x and store in f.
!     end if
!     if (task .eq. 'GH' .or. task .eq. 'START') then
!        Evaluate the gradient at x and store in g.
!        Evaluate the Hessian at x and store in compressed
!        column storage in (a, adiag, acol_ptr, arow_ind)
!     end if

!     call dtron(n, x, xl, xu, f, g, a, adiag, acol_ptr, arow_ind,
!                frtol, fatol, fmin, cgtol, itermax, delta, task,
!                b, bdiag, bcol_ptr, brow_ind,
!                l, ldiag, lcol_ptr, lrow_ind, iterscg)

!     if (task(1:4) .eq. 'CONV') search = .false.

!   end do

!  NOTE: The user must not alter work arrays between calls.

!  The subroutine statement is

!    subroutine dtron(n, x, xl, xu, f, g, a, adiag, acol_ptr, arow_ind,
!                     frtol, fatol, fmin, cgtol, itermax, delta, task,
!                     b, bdiag, bcol_ptr, brow_ind,
!                     l, ldiag, lcol_ptr, lrow_ind, isave)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is the final minimizer.

!    xl is a REAL (dp) array of dimension n.
!      On entry xl is the vector of lower bounds.
!      On exit xl is unchanged.

!    xu is a REAL (dp) array of dimension n.
!      On entry xu is the vector of upper bounds.
!      On exit xu is unchanged.

!    f is a REAL (dp) variable.
!      On entry f must contain the function at x.
!      On exit f is unchanged.   This is NOT true.

!    g is a REAL (dp) array of dimension n.
!      On entry g must contain the gradient at x.
!      On exit g is unchanged.

!    a is a REAL (dp) array of dimension nnz.
!      On entry a must contain the strict lower triangular part
!         of A in compressed column storage.
!      On exit a is unchanged.

!    adiag is a REAL (dp) array of dimension n.
!      On entry adiag must contain the diagonal elements of A.
!      On exit adiag is unchanged.

!    acol_ptr is an integer array of dimension n + 1.
!      On entry acol_ptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         acol_ptr(j), ... , acol_ptr(j+1) - 1.
!      On exit acol_ptr is unchanged.

!    arow_ind is an integer array of dimension nnz.
!      On entry arow_ind must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit arow_ind is unchanged.

!    frtol is a REAL (dp) variable.
!      On entry frtol specifies the relative error desired in the function.
!         Convergence occurs if the estimate of the relative error between f(x)
!         and f(xsol), where xsol is a local minimizer, is less than frtol.
!      On exit frtol is unchanged.

!    fatol is a REAL (dp) variable.
!      On entry fatol specifies the absolute error desired in the function.
!         Convergence occurs if the estimate of the absolute error between f(x)
!         and f(xsol), where xsol is a local minimizer, is less than fatol.
!      On exit fatol is unchanged.

!    fmin is a REAL (dp) variable.
!      On entry fmin specifies a lower bound for the function.
!         The subroutine exits with a warning if f < fmin.
!      On exit fmin is unchanged.

!    cgtol is a REAL (dp) variable.
!      On entry cgtol specifies the convergence criteria for
!         the conjugate gradient method.
!      On exit cgtol is unchanged.

!    itermax is an integer variable.
!      On entry itermax specifies the limit on the number of
!         conjugate gradient iterations.
!      On exit itermax is unchanged.

!    delta is a REAL (dp) variable.
!      On entry delta is the trust region bound.
!      On exit delta is unchanged.  This is NOT true.

!    task is a character variable of length at least 60.
!      On initial entry task must be set to 'START'.
!      On exit task indicates the required action:

!         If task(1:1) = 'F' then evaluate the function at x.

!         If task(1:2) = 'GH' then evaluate the gradient and the
!         Hessian matrix at x.

!         If task(1:4) = 'CONV' then the search is successful.

!         If task(1:4) = 'WARN' then the subroutine is not able
!         to satisfy the convergence conditions. The exit value
!         of x contains the best approximation found.

!    bdiag is a REAL (dp) array of dimension n.
!      On entry bdiag need not be specified.
!      On exit bdiag contains the diagonal elements of B.

!    bcol_ptr is an integer array of dimension n + 1.
!      On entry bcol_ptr need not be specified
!      On exit bcol_ptr contains pointers to the columns of B.
!         The nonzeros in column j of B are in the
!         bcol_ptr(j), ... , bcol_ptr(j+1) - 1 positions of b.

!    brow_ind is an integer array of dimension nnz.
!      On entry brow_ind need not be specified.
!      On exit brow_ind contains row indices for the strict lower
!         triangular part of B in compressed column storage.

!    l is a REAL (dp) array of dimension nnz + n*p.
!      On entry l need not be specified.
!      On exit l contains the strict lower triangular part
!         of L in compressed column storage.

!    ldiag is a REAL (dp) array of dimension n.
!      On entry ldiag need not be specified.
!      On exit ldiag contains the diagonal elements of L.

!    lcol_ptr is an integer array of dimension n + 1.
!      On entry lcol_ptr need not be specified.
!      On exit lcol_ptr contains pointers to the columns of L.
!         The nonzeros in column j of L are in the
!         lcol_ptr(j), ... , lcol_ptr(j+1) - 1 positions of l.

!    lrow_ind is an integer array of dimension nnz + n*p.
!      On entry lrow_ind need not be specified.
!      On exit lrow_ind contains row indices for the strict lower
!         triangular part of L in compressed column storage.

!  Subprograms called

!    MINPACK-2  ......  dcauchy, dspcg, dssyax

!    Level 1 BLAS  ...  dcopy

!  MINPACK-2 Project. May 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp

!     Parameters for updating the iterates.

REAL (dp), PARAMETER :: eta0=1D-4, eta1=0.25_dp, eta2=0.75_dp

!     Parameters for updating the trust region size delta.

REAL (dp), PARAMETER :: sigma1=0.25_dp, sigma2=0.5_dp, sigma3=4.0_dp

LOGICAL   :: search
INTEGER   :: info, iter, iters
REAL (dp) :: alphac, fc, prered, actred, snorm, gs, alpha
CHARACTER (LEN=60) :: work

! REAL (dp) :: ddot, dnrm2
! EXTERNAL dcauchy, dspcg, dssyax
! EXTERNAL dcopy, ddot, dnrm2

!     Initialization section.

IF (task(1:5) == 'START') THEN
  
!        Initialize local variables.
  
  iter = 1
  iterscg = 0
  alphac = one
  work = 'COMPUTE'
  IF (ALLOCATED(xc)) DEALLOCATE( xc, s, indfree, dsave, wa, isave )
  ALLOCATE( xc(n), s(n), indfree(n), dsave(n), wa(n), isave(n) )
  
ELSE
  
!        Restore local variables.
  
  IF (isave(1) == 1) THEN
    work = 'COMPUTE'
  ELSE IF (isave(1) == 2) THEN
    work = 'EVALUATE'
  END IF
  iter = isave(2)
  iterscg = isave(3)
  fc = dsave(1)
  alphac = dsave(2)
END IF

!     Search for a lower function value.

search = .true.
DO WHILE (search)
  
!        Compute a step and evaluate the function at the trial point.
  
  IF (work == 'COMPUTE') THEN
    
!           Save the best function value and the best x.
    
    fc = f
    xc(1:n) = x(1:n)
    
!           Compute the Cauchy step and store in s.
    
    CALL dcauchy(n, x, xl, xu, a, adiag, acol_ptr, arow_ind, g, delta,  &
                 alphac, s)
    
!           Compute the projected Newton step.
    
    CALL dspcg(n, x, xl, xu, a, adiag, acol_ptr, arow_ind, g, delta,  &
               cgtol, s, 5, itermax, iters, info, b, bdiag, bcol_ptr, brow_ind,  &
               l, ldiag, lcol_ptr, lrow_ind)
    
    iterscg = iterscg + iters
    task = 'F'
  END IF
  
!        Evaluate the step and determine if the step is successful.
  
  IF (work == 'EVALUATE') THEN
    
!           Compute the predicted reduction.
    
    CALL dssyax(n, a, adiag, acol_ptr, arow_ind, s, wa)
    prered = DOT_PRODUCT( s(1:n), p5*wa(1:n) - g(1:n) )
    
!           Compute the actual reduction.
    
    actred =  fc - f
    
!           On the first iteration, adjust the initial step bound.
    
    snorm = dnrm2(n, s, 1)
    IF (iter == 1)  delta = MIN(delta, snorm)
    
!           Compute prediction alpha*snorm of the step.
    
    gs = DOT_PRODUCT( g(1:n), s(1:n) )
    IF (f-fc-gs <= zero) THEN
      alpha = sigma3
    ELSE
      alpha = MAX(sigma1, -p5*(gs/(f-fc-gs)))
    END IF
    
!           Update the trust region bound according to the ratio
!           of actual to predicted reduction.
    
    IF (actred < eta0*prered) THEN
      
!              Reduce delta.  Step is not successful.
      
      delta = MIN(MAX(alpha, sigma1)*snorm, sigma2*delta)
      
    ELSE IF (actred < eta1*prered) THEN
      
!              Reduce delta.  Step is not sufficiently successful.
      
      delta = MAX(sigma1*delta, MIN(alpha*snorm, sigma2*delta))
      
    ELSE IF (actred < eta2*prered) THEN
      
!              The ratio of actual to predicted reduction is in
!              the interval (eta1, eta2).  We are allowed to either
!              increase or decrease delta.
      
      delta = MAX(sigma1*delta, MIN(alpha*snorm, sigma3*delta))
      
    ELSE
      
!              The ratio of actual to predicted reduction exceeds eta2.
!              Do not decrease delta.
      
      delta = MAX(delta, MIN(alpha*snorm, sigma3*delta))
      
    END IF
    
!           Update the iterate.
    
    IF (actred > eta0*prered) THEN
      
!              Successful iterate.
      
      task = 'GH'
      iter = iter + 1
      
    ELSE
      
!              Unsuccessful iterate.
      
      task = 'NFGH'
      x(1:n) = xc(1:n)
      f = fc
      
    END IF
    
!           Test for convergence.
    
    IF (f < fmin) task = 'WARNING: F < FMIN'
    IF (ABS(actred) <= fatol .AND. prered <= fatol) task =  &
        'CONVERGENCE: FATOL TEST SATISFIED'
    IF (ABS(actred) <= frtol*ABS(f) .AND. prered <= frtol*ABS(f)) task =  &
        'CONVERGENCE: FRTOL TEST SATISFIED'
    
  END IF
  
!        Determine what needs to be done on the next call.
  
  IF (task == 'F') THEN
    work = 'EVALUATE'
  ELSE
    work = 'COMPUTE'
  END IF
  
!        Continue the search if the step is not successful.
  
  IF (task /= 'NFGH') search = .false.
  
END DO

!     Save local variables.

IF (work == 'COMPUTE') THEN
  isave(1) = 1
ELSE IF (work == 'EVALUATE') THEN
  isave(1) = 2
END IF

isave(2) = iter
isave(3) = iterscg
dsave(1) = fc
dsave(2) = alphac

RETURN
END SUBROUTINE dtron



SUBROUTINE dcauchy(n, x, xl, xu, a, diag, col_ptr, row_ind, g, delta,  &
                   alpha, s)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:33

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN)      :: xl(:)
REAL (dp), INTENT(IN)      :: xu(:)
REAL (dp), INTENT(IN)      :: a(:)
REAL (dp), INTENT(IN)      :: diag(:)
INTEGER, INTENT(IN)        :: col_ptr(:)    ! col_ptr(n+1)
INTEGER, INTENT(IN)        :: row_ind(:)
REAL (dp), INTENT(IN)      :: g(:)
REAL (dp), INTENT(IN)      :: delta
REAL (dp), INTENT(IN OUT)  :: alpha
REAL (dp), INTENT(OUT)     :: s(:)


!  **********

!  Subroutine dcauchy

!  This subroutine computes a Cauchy step that satisfies a trust
!  region constraint and a sufficient decrease condition.

!  The Cauchy step is computed for the quadratic

!        q(s) = 0.5*s'*A*s + g'*s,

!  where A is a symmetric matrix in compressed row storage, and
!  g is a vector. Given a parameter alpha, the Cauchy step is

!        s[alpha] = P[x - alpha*g] - x,

!  with P the projection onto the n-dimensional interval [xl,xu].
!  The Cauchy step satisfies the trust region constraint and the
!  sufficient decrease condition

!        || s || <= delta,      q(s) <= mu_0*(g'*s),

!  where mu_0 is a constant in (0,1).

!  The subroutine statement is

!    subroutine dcauchy(n, x, xl, xu, a, diag, col_ptr, row_ind, g, delta,
!                       alpha, s)

!  where


!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is unchanged.

!    xl is a REAL (dp) array of dimension n.
!      On entry xl is the vector of lower bounds.
!      On exit xl is unchanged.

!    xu is a REAL (dp) array of dimension n.
!      On entry xu is the vector of upper bounds.
!      On exit xu is unchanged.

!    a is a REAL (dp) array of dimension nnz.
!      On entry a must contain the strict lower triangular part
!         of A in compressed column storage.
!      On exit a is unchanged.

!    diag is a REAL (dp) array of dimension n.
!      On entry diag must contain the diagonal elements of A.
!      On exit diag is unchanged.

!    col_ptr is an integer array of dimension n + 1.
!      On entry col_ptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         col_ptr(j), ... , col_ptr(j+1) - 1.
!      On exit col_ptr is unchanged.

!    row_ind is an integer array of dimension nnz.
!      On entry row_ind must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit row_ind is unchanged.

!    g is a REAL (dp) array of dimension n.
!      On entry g specifies the gradient g.
!      On exit g is unchanged.

!    delta is a REAL (dp) variable.
!      On entry delta is the trust region size.
!      On exit delta is unchanged.

!    alpha is a REAL (dp) variable.
!      On entry alpha is the current estimate of the step.
!      On exit alpha defines the Cauchy step s[alpha].

!    s is a REAL (dp) array of dimension n.
!      On entry s need not be specified.
!      On exit s is the Cauchy step s[alpha].

!  Subprograms called

!    MINPACK-2  ......  dbreakpt, dgpstep, dssyax

!    Level 1 BLAS  ...  ddot, dnrm2

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp

!     Constant that defines sufficient decrease.

REAL (dp), PARAMETER :: mu0=0.01_dp

!     Interpolation and extrapolation factors.

REAL (dp), PARAMETER :: interpf=0.1_dp, extrapf=10.0_dp

LOGICAL   :: search, interp
INTEGER   :: nbrpt, nsteps
REAL (dp) :: alphas, brptmax, brptmin, gts, q
REAL (dp) :: wa(n)

! REAL (dp) :: dnrm2, ddot
! EXTERNAL dnrm2, ddot

!     Find the minimal and maximal break-point on x - alpha*g.

wa(1:n) = g(1:n)
wa(1:n) = - wa(1:n)
CALL dbreakpt(n, x, xl, xu, g, nbrpt, brptmin, brptmax)

!     Evaluate the initial alpha and decide if the algorithm
!     must interpolate or extrapolate.

CALL dgpstep(n, x, xl, xu, -alpha, g, s)
IF (dnrm2(n, s, 1) > delta) THEN
  interp = .true.
ELSE
  CALL dssyax(n, a, diag, col_ptr, row_ind, s, wa)
  gts = DOT_PRODUCT( g(1:n), s(1:n) )
  q = p5*DOT_PRODUCT( s(1:n), wa(1:n) ) + gts
  interp = (q >= mu0*gts)
END IF

nsteps = 1

!     Either interpolate or extrapolate to find a successful step.

IF (interp) THEN
  
!        Reduce alpha until a successful step is found.
  
  search = .true.
  DO WHILE (search)
    
!           This is a crude interpolation procedure that
!           will be replaced in future versions of the code.
    
    nsteps = nsteps + 1
    alpha = interpf*alpha
    CALL dgpstep(n, x, xl, xu, -alpha, g, s)
    IF (dnrm2(n, s, 1) <= delta) THEN
      CALL dssyax(n, a, diag, col_ptr, row_ind, s, wa)
      gts = DOT_PRODUCT( g(1:n), s(1:n) )
      q = p5*DOT_PRODUCT( s(1:n), wa(1:n) ) + gts
      search = (q > mu0*gts)
    END IF
  END DO
  
ELSE
  
!        Increase alpha until a successful step is found.
  
  search = .true.
  DO WHILE (search .AND. alpha <= brptmax)
    
!           This is a crude extrapolation procedure that
!           will be replaced in future versions of the code.
    
    nsteps = nsteps + 1
    alphas = alpha
    alpha = extrapf*alpha
    CALL dgpstep(n, x, xl, xu, -alpha, g, s)
    IF (dnrm2(n, s, 1) <= delta) THEN
      CALL dssyax(n, a, diag, col_ptr, row_ind, s, wa)
      gts = DOT_PRODUCT( g(1:n), s(1:n) )
      q = p5*DOT_PRODUCT( s(1:n), wa(1:n) ) + gts
      search = (q < mu0*gts)
    ELSE
      search = .false.
    END IF
  END DO
  
!        Recover the last successful step.
  
  alpha = alphas
  CALL dgpstep(n, x, xl, xu, -alpha, g, s)
END IF

RETURN

END SUBROUTINE dcauchy



SUBROUTINE dspcg(n, x, xl, xu, a, adiag, acol_ptr, arow_ind, g, delta,  &
                 rtol, s, nv, itermax, iters, info,  &
                 b, bdiag, bcol_ptr, brow_ind,  &
                 l, ldiag, lcol_ptr, lrow_ind)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: xl(:)
REAL (dp), INTENT(IN)      :: xu(:)
REAL (dp), INTENT(IN)      :: a(:)
REAL (dp), INTENT(IN)      :: adiag(:)
INTEGER, INTENT(IN)        :: acol_ptr(:)    ! acol_ptr(n+1)
INTEGER, INTENT(IN)        :: arow_ind(:)
REAL (dp), INTENT(IN)      :: g(:)
REAL (dp), INTENT(IN)      :: delta
REAL (dp), INTENT(IN)      :: rtol
REAL (dp), INTENT(IN OUT)  :: s(:)
INTEGER, INTENT(IN)        :: nv
INTEGER, INTENT(IN)        :: itermax
INTEGER, INTENT(OUT)       :: iters
INTEGER, INTENT(OUT)       :: info
REAL (dp), INTENT(OUT)     :: b(:)
REAL (dp), INTENT(OUT)     :: bdiag(:)
INTEGER, INTENT(OUT)       :: bcol_ptr(:)    ! bcol_ptr(n+1)
INTEGER, INTENT(OUT)       :: brow_ind(:)
REAL (dp), INTENT(OUT)     :: l(:)
REAL (dp), INTENT(OUT)     :: ldiag(:)
INTEGER, INTENT(OUT)       :: lcol_ptr(:)    ! lcol_ptr(n+1)
INTEGER, INTENT(OUT)       :: lrow_ind(:)


!  *********

!  Subroutine dspcg

!  This subroutine generates a sequence of approximate minimizers
!  for the subproblem

!        min { q(x) : xl <= x <= xu }.

!  The quadratic is defined by

!        q(x[0]+s) = 0.5*s'*A*s + g'*s,

!  where x[0] is a base point provided by the user, A is a symmetric
!  matrix in compressed column storage, and g is a vector.

!  At each stage we have an approximate minimizer x[k], and generate a direction
!  p[k] by using a preconditioned conjugate gradient method on the subproblem

!        min { q(x[k]+p) : || L'*p || <= delta, s(fixed) = 0 },

!  where fixed is the set of variables fixed at x[k], delta is the trust region
!  bound, and L is an incomplete Cholesky factorization of the submatrix

!        B = A(free:free),

!  where free is the set of free variables at x[k]. Given p[k],
!  the next minimizer x[k+1] is generated by a projected search.

!  The starting point for this subroutine is x[1] = x[0] + s, where
!  x[0] is a base point and s is the Cauchy step.

!  The subroutine converges when the step s satisfies

!        || (g + A*s)[free] || <= rtol*|| g[free] ||

!  In this case the final x is an approximate minimizer in the
!  face defined by the free variables.

!  The subroutine terminates when the trust region bound does
!  not allow further progress, that is, || L'*p[k] || = delta.
!  In this case the final x satisfies q(x) < q(x[k]).

!  The subroutine statement is

!    subroutine dspcg(n, x, xl, xu, a, adiag, acol_ptr, arow_ind, g, delta,
!                     rtol, s, nv, itermax, iters, info,
!                     b, bdiag, bcol_ptr, brow_ind,
!                     l, ldiag, lcol_ptr, lrow_ind)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is the final minimizer.

!    xl is a REAL (dp) array of dimension n.
!      On entry xl is the vector of lower bounds.
!      On exit xl is unchanged.

!    xu is a REAL (dp) array of dimension n.
!      On entry xu is the vector of upper bounds.
!      On exit xu is unchanged.

!    a is a REAL (dp) array of dimension nnz.
!      On entry a must contain the strict lower triangular part
!         of A in compressed column storage.
!      On exit a is unchanged.

!    adiag is a REAL (dp) array of dimension n.
!      On entry adiag must contain the diagonal elements of A.
!      On exit adiag is unchanged.

!    acol_ptr is an integer array of dimension n + 1.
!      On entry acol_ptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         acol_ptr(j), ... , acol_ptr(j+1) - 1.
!      On exit acol_ptr is unchanged.

!    arow_ind is an integer array of dimension nnz.
!      On entry arow_ind must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit arow_ind is unchanged.

!    g is a REAL (dp) array of dimension n.
!      On entry g must contain the vector g.
!      On exit g is unchanged.

!    delta is a REAL (dp) variable.
!      On entry delta is the trust region size.
!      On exit delta is unchanged.

!    rtol is a REAL (dp) variable.
!      On entry rtol specifies the accuracy of the final minimizer.
!      On exit rtol is unchanged.

!    s is a REAL (dp) array of dimension n.
!      On entry s is the Cauchy step.
!      On exit s contain the final step.

!    nv is an integer variable.
!      On entry nv specifies the amount of memory available for the
!         incomplete Cholesky factorization.
!      On exit nv is unchanged.

!    itermax is an integer variable.
!      On entry itermax specifies the limit on the number of
!         conjugate gradient iterations.
!      On exit itermax is unchanged.

!    iters is an integer variable.
!      On entry iters need not be specified.
!      On exit iters is set to the number of conjugate
!         gradient iterations.

!    info is an integer variable.
!      On entry info need not be specified.
!      On exit info is set as follows:

!          info = 1  Convergence. The final step s satisfies
!                    || (g + A*s)[free] || <= rtol*|| g[free] ||,
!                    and the final x is an approximate minimizer
!                    in the face defined by the free variables.

!          info = 2  Termination. The trust region bound does
!                    not allow further progress.

!          info = 3  Failure to converge within itermax iterations.

!    b is a REAL (dp) array of dimension nnz.
!      On entry b need not be specified.
!      On exit b contains the strict lower triangular part
!         of B in compressed column storage.

!    bdiag is a REAL (dp) array of dimension n.
!      On entry bdiag need not be specified.
!      On exit bdiag contains the diagonal elements of B.

!    bcol_ptr is an integer array of dimension n + 1.
!      On entry bcol_ptr need not be specified
!      On exit bcol_ptr contains pointers to the columns of B.
!         The nonzeros in column j of B are in the
!         bcol_ptr(j), ... , bcol_ptr(j+1) - 1 positions of b.

!    brow_ind is an integer array of dimension nnz.
!      On entry brow_ind need not be specified.
!      On exit brow_ind contains row indices for the strict lower
!         triangular part of B in compressed column storage.

!    l is a REAL (dp) array of dimension nnz + n*nv.
!      On entry l need not be specified.
!      On exit l contains the strict lower triangular part
!         of L in compressed column storage.

!    ldiag is a REAL (dp) array of dimension n.
!      On entry ldiag need not be specified.
!      On exit ldiag contains the diagonal elements of L.

!    lcol_ptr is an integer array of dimension n + 1.
!      On entry lcol_ptr need not be specified.
!      On exit lcol_ptr contains pointers to the columns of L.
!         The nonzeros in column j of L are in the
!         lcol_ptr(j), ... , lcol_ptr(j+1) - 1 positions of l.

!    lrow_ind is an integer array of dimension nnz + n*nv.
!      On entry lrow_ind need not be specified.
!      On exit lrow_ind contains row indices for the strict lower
!         triangular part of L in compressed column storage.

!  Subprograms called

!    MINPACK-2  ......  dmid, dprsrch, dtrpcg,
!                       dssyax, dstrsol

!    Level 1 BLAS  ...  daxpy, dnrm2

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER   :: indfree(n), iwa(3*n)
REAL (dp) :: gfree(n), w(n), wa(3*n)
INTEGER   :: infotr, ip, j, jfree, itertr, nfaces, nfree, nnz
REAL (dp) :: alpha, gfnorm, gfnormf
REAL (dp) :: tol, stol

! REAL (dp) :: dnrm2
! EXTERNAL dmid, dprsrch, dtrpcg, dssyax, dstrsol
! EXTERNAL daxpy, dnrm2

!     Compute A*(x[1] - x[0]) and store in w.

CALL dssyax(n, a, adiag, acol_ptr, arow_ind, s, w)

!     Compute the Cauchy point.

x(1:n) = x(1:n) + s(1:n)
CALL dmid(n, x, xl, xu)

!     Start the main iteration loop.
!     There are at most n iterations because at each iteration
!     at least one variable becomes active.

iters = 0
DO nfaces = 1, n
  
!        Determine the free variables at the current minimizer.
!        The indices of the free variables are stored in the first
!        n free positions of the array indfree.
!        The array iwa is used to detect free variables by setting
!        iwa(i) = 1 if the ith variable is free, otherwise iwa(i) = 0.
  
  nfree  = 0
  DO j = 1, n
    IF (xl(j) < x(j) .AND. x(j) < xu(j)) THEN
      nfree = nfree + 1
      indfree(nfree) = j
      iwa(j) = nfree
    ELSE
      iwa(j) = 0
    END IF
  END DO
  
!        Exit if there are no free constraints.
  
  IF (nfree == 0) THEN
    info = 1
    RETURN
  END IF
  
!        Obtain the submatrix of A for the free variables.
!        Recall that iwa allows the detection of free variables.
  
  bcol_ptr(1) = 1
  nnz = 0
  DO j = 1, nfree
    jfree = indfree(j)
    bdiag(j) = adiag(jfree)
    DO ip = acol_ptr(jfree), acol_ptr(jfree+1)-1
      IF (iwa(arow_ind(ip)) > 0) THEN
        nnz = nnz + 1
        brow_ind(nnz) = iwa(arow_ind(ip))
        b(nnz) = a(ip)
      END IF
    END DO
    bcol_ptr(j+1) = nnz + 1
  END DO
  
!        Compute the incomplete Cholesky factorization.
!        The variable nnz is used to specify an upper bound on the
!        length of b, so we need to make sure that nnz >= 1.
  
  alpha = zero
  nnz = MAX(nnz,1)
  CALL dicfs(nfree, nnz, b, bdiag, bcol_ptr, brow_ind,  &
             l, ldiag, lcol_ptr, lrow_ind, nv, alpha)
  
!        Compute the gradient grad q(x[k]) = g + A*(x[k] - x[0]),
!        of q at x[k] for the free variables.
!        Recall that w contains  A*(x[k] - x[0]).
!        Compute the norm of the reduced gradient Z'*g.
  
  DO j = 1, nfree
    gfree(j) = w(indfree(j)) + g(indfree(j))
    wa(j) = g(indfree(j))
  END DO
  gfnorm = dnrm2(nfree, wa, 1)
  
!        Solve the trust region subproblem in the free variables
!        to generate a direction p[k]. Store p[k] in the array w.
  
  tol = rtol*gfnorm
  stol = zero
  
  CALL dtrpcg(nfree, b, bdiag, bcol_ptr, brow_ind, gfree, delta,  &
              l, ldiag, lcol_ptr, lrow_ind, tol, stol, itermax, w, itertr,  &
              infotr)
  
  iters = iters + itertr
  CALL dstrsol(nfree, l, ldiag, lcol_ptr, lrow_ind, w, 'T')
  
!        Use a projected search to obtain the next iterate.
!        The projected search algorithm stores s[k] in w.
  
  DO j = 1, nfree
    wa(j) = x(indfree(j))
    wa(n+j) = xl(indfree(j))
    wa(2*n+j) = xu(indfree(j))
  END DO
  
  CALL dprsrch(nfree, wa, wa(n+1:), wa(2*n+1:),  &
               b, bdiag, bcol_ptr, brow_ind, gfree, w)
  
!        Update the minimizer and the step.
!        Note that s now contains x[k+1] - x[0].
  
  DO j = 1, nfree
    x(indfree(j)) = wa(j)
    s(indfree(j)) = s(indfree(j)) + w(j)
  END DO
  
!        Compute A*(x[k+1] - x[0]) and store in w.
  
  CALL dssyax(n, a, adiag, acol_ptr, arow_ind, s, w)
  
!        Compute the gradient grad q(x[k+1]) = g + A*(x[k+1] - x[0])
!        of q at x[k+1] for the free variables.
  
  DO j = 1, nfree
    gfree(j) = w(indfree(j)) + g(indfree(j))
  END DO
  gfnormf = dnrm2(nfree, gfree, 1)
  
!        Convergence and termination test.
!        We terminate if the preconditioned conjugate gradient method
!        encounters a direction of negative curvature, or
!        if the step is at the trust region bound.
  
  IF (gfnormf <= rtol*gfnorm) THEN
    info = 1
    RETURN
  ELSE IF (infotr == 3 .OR. infotr == 4) THEN
    info = 2
    RETURN
  ELSE IF (iters > itermax) THEN
    info = 3
    RETURN
  END IF
  
END DO

RETURN

END SUBROUTINE dspcg



SUBROUTINE dbreakpt(n, x, xl, xu, w, nbrpt, brptmin, brptmax)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(IN)   :: xl(:)
REAL (dp), INTENT(IN)   :: xu(:)
REAL (dp), INTENT(IN)   :: w(:)
INTEGER, INTENT(OUT)    :: nbrpt
REAL (dp), INTENT(OUT)  :: brptmin
REAL (dp), INTENT(OUT)  :: brptmax


!  **********

!  Subroutine dbreakpt

!  This subroutine computes the number of break-points, and
!  the minimal and maximal break-points of the projection of
!  x + alpha*w on the n-dimensional interval [xl,xu].

!  The subroutine statement is

!    subroutine dbreakpt(n, x, xl, xu, w, nbrpt, brptmin, brptmax)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is unchanged.

!    xl is a REAL (dp) array of dimension n.
!      On entry xl is the vector of lower bounds.
!      On exit xl is unchanged.

!    xu is a REAL (dp) array of dimension n.
!      On entry xu is the vector of upper bounds.
!      On exit xu is unchanged.

!    w is a REAL (dp) array of dimension n.
!      On entry w specifies the vector w.
!      On exit w is unchanged.

!    nbrpt is an integer variable.
!      On entry nbrpt need not be specified.
!      On exit nbrpt is the number of break points.

!    brptmin is a REAL (dp) variable
!      On entry brptmin need not be specified.
!      On exit brptmin is minimal break-point.

!    brptmax is a REAL (dp) variable
!      On entry brptmax need not be specified.
!      On exit brptmax is maximal break-point.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********


INTEGER   :: i
REAL (dp) :: brpt

nbrpt = 0
DO i = 1, n
  IF (x(i) < xu(i) .AND. w(i) > zero) THEN
    nbrpt = nbrpt + 1
    brpt =  (xu(i) - x(i))/w(i)
    IF (nbrpt == 1) THEN
      brptmin = brpt
      brptmax = brpt
    ELSE
      brptmin = MIN(brpt,brptmin)
      brptmax = MAX(brpt,brptmax)
    END IF
  ELSE IF (x(i) > xl(i) .AND. w(i) < zero) THEN
    nbrpt = nbrpt + 1
    brpt = (xl(i) - x(i))/w(i)
    IF (nbrpt == 1) THEN
      brptmin = brpt
      brptmax = brpt
    ELSE
      brptmin = MIN(brpt,brptmin)
      brptmax = MAX(brpt,brptmax)
    END IF
  END IF
END DO

!     Handle the exceptional case.

IF (nbrpt == 0) THEN
  brptmin = zero
  brptmax = zero
END IF

RETURN

END SUBROUTINE dbreakpt



SUBROUTINE dgpstep(n, x, xl, xu, alpha, w, s)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(IN)   :: xl(:)
REAL (dp), INTENT(IN)   :: xu(:)
REAL (dp), INTENT(IN)   :: alpha
REAL (dp), INTENT(IN)   :: w(:)
REAL (dp), INTENT(OUT)  :: s(:)

!  **********

!  Subroutine dgpstep

!  This subroutine computes the gradient projection step

!        s = P[x + alpha*w] - x,

!  where P is the projection on the n-dimensional interval [xl,xu].

!  The subroutine statement is

!    subroutine  dgpstep(n,x,xl,xu,alpha,w,s)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is unchanged.

!    xl is a REAL (dp) array of dimension n.
!      On entry xl is the vector of lower bounds.
!      On exit xl is unchanged.

!    xu is a REAL (dp) array of dimension n.
!      On entry xu is the vector of upper bounds.
!      On exit xu is unchanged.

!    alpha is a REAL (dp) variable.
!      On entry alpha specifies the scalar alpha.
!      On exit alpha is unchanged.

!    w is a REAL (dp) array of dimension n.
!      On entry w specifies the vector w.
!      On exit w is unchanged.

!    s is a REAL (dp) array of dimension n.
!      On entry s need not be specified.
!      On exit s contains the gradient projection step.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER   :: i

!     This computation of the gradient projection step avoids
!     rounding errors for the components that are feasible.

DO i = 1, n
  IF (x(i) + alpha*w(i) < xl(i)) THEN
    s(i) = xl(i) - x(i)
  ELSE IF (x(i) + alpha*w(i) > xu(i)) THEN
    s(i) = xu(i) - x(i)
  ELSE
    s(i) = alpha*w(i)
  END IF
END DO

RETURN

END SUBROUTINE dgpstep



SUBROUTINE dmid(n, x, xl, xu)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: xl(:)
REAL (dp), INTENT(IN)      :: xu(:)

!  **********

!  Subroutine dmid

!  This subroutine computes the projection of x
!  on the n-dimensional interval [xl,xu].

!  The subroutine statement is

!    subroutine dmid(n, x, xl, xu)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is the projection of x on [xl,xu].

!    xl is a REAL (dp) array of dimension n.
!      On entry xl is the vector of lower bounds.
!      On exit xl is unchanged.

!    xu is a REAL (dp) array of dimension n.
!      On entry xu is the vector of upper bounds.
!      On exit xu is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER   :: i

DO i = 1, n
  x(i) = MAX(xl(i), MIN(x(i), xu(i)))
END DO

RETURN
END SUBROUTINE dmid



SUBROUTINE dprsrch(n, x, xl, xu, a, diag, col_ptr, row_ind, g, w)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: xl(:)
REAL (dp), INTENT(IN)      :: xu(:)
REAL (dp), INTENT(IN)      :: a(:)
REAL (dp), INTENT(IN)      :: diag(:)
INTEGER, INTENT(IN)        :: col_ptr(:)    ! col_ptr(n+1)
INTEGER, INTENT(IN)        :: row_ind(:)
REAL (dp), INTENT(IN)      :: g(:)
REAL (dp), INTENT(IN OUT)  :: w(:)

!  **********

!  Subroutine dprsrch

!  This subroutine uses a projected search to compute a step
!  that satisfies a sufficient decrease condition for the quadratic

!        q(s) = 0.5*s'*A*s + g'*s,

!  where A is a symmetric matrix in compressed column storage,
!  and g is a vector. Given the parameter alpha, the step is

!        s[alpha] = P[x + alpha*w] - x,

!  where w is the search direction and P the projection onto the
!  n-dimensional interval [xl,xu]. The final step s = s[alpha]
!  satisfies the sufficient decrease condition

!        q(s) <= mu_0*(g'*s),

!  where mu_0 is a constant in (0,1).

!  The search direction w must be a descent direction for the
!  quadratic q at x such that the quadratic is decreasing
!  in the ray  x + alpha*w for 0 <= alpha <= 1.

!  The subroutine statement is

!    subroutine dprsrch(n, x, xl, xu, a, diag, col_ptr, row_ind, g, w)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is set to the final point P[x + alpha*w].

!    xl is a REAL (dp) array of dimension n.
!      On entry xl is the vector of lower bounds.
!      On exit xl is unchanged.

!    xu is a REAL (dp) array of dimension n.
!      On entry xu is the vector of upper bounds.
!      On exit xu is unchanged.

!    a is a REAL (dp) array of dimension nnz.
!      On entry a must contain the strict lower triangular part
!         of A in compressed column storage.
!      On exit a is unchanged.

!    diag is a REAL (dp) array of dimension n.
!      On entry diag must contain the diagonal elements of A.
!      On exit diag is unchanged.

!    col_ptr is an integer array of dimension n + 1.
!      On entry col_ptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         col_ptr(j), ... , col_ptr(j+1) - 1.
!      On exit col_ptr is unchanged.

!    row_ind is an integer array of dimension nnz.
!      On entry row_ind must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit row_ind is unchanged.

!    g is a REAL (dp) array of dimension n.
!      On entry g specifies the vector g.
!      On exit g is unchanged.

!    w is a double prevision array of dimension n.
!      On entry w specifies the search direction.
!      On exit w is the step s[alpha].

!  Subprograms called

!    MINPACK-2  ......  dbreakpt, dgpstep, dmid, dssyax

!    Level 1 BLAS  ...  daxpy, dcopy, ddot

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp

!     Constant that defines sufficient decrease.

REAL (dp), PARAMETER :: mu0=0.01_dp

!     Interpolation factor.

REAL (dp), PARAMETER :: interpf=0.5_dp

REAL (dp) :: wa1(n), wa2(n)
LOGICAL   :: search
INTEGER   :: nbrpt, nsteps
REAL (dp) :: alpha, brptmin, brptmax, gts, q

! REAL (dp) :: ddot
! EXTERNAL daxpy, dcopy, ddot
! EXTERNAL dbreakpt, dgpstep, dmid, dssyax

!     Set the initial alpha = 1 because the quadratic function is
!     decreasing in the ray x + alpha*w for 0 <= alpha <= 1.

alpha = one
nsteps = 0

!     Find the smallest break-point on the ray x + alpha*w.

CALL dbreakpt(n, x, xl, xu, w, nbrpt, brptmin, brptmax)

!     Reduce alpha until the sufficient decrease condition is
!     satisfied or x + alpha*w is feasible.

search = .true.
DO WHILE (search .AND. alpha > brptmin)
  
!        Calculate P[x + alpha*w] - x and check the sufficient
!        decrease condition.
  
  nsteps = nsteps + 1
  CALL dgpstep(n, x, xl, xu, alpha, w, wa1)
  CALL dssyax(n, a, diag, col_ptr, row_ind, wa1, wa2)
  gts = DOT_PRODUCT( g(1:n), wa1(1:n) )
  q = p5*DOT_PRODUCT( wa1(1:n), wa2(1:n) ) + gts
  IF (q <= mu0*gts) THEN
    search = .false.
  ELSE
    
!           This is a crude interpolation procedure that
!           will be replaced in future versions of the code.
    
    alpha = interpf*alpha
    
  END IF
END DO

!     Force at least one more constraint to be added to the active
!     set if alpha < brptmin and the full step is not successful.
!     There is sufficient decrease because the quadratic function
!     is decreasing in the ray x + alpha*w for 0 <= alpha <= 1.

IF (alpha < one .AND. alpha < brptmin) alpha = brptmin

!     Compute the final iterate and step.

CALL dgpstep(n, x, xl, xu, alpha, w, wa1)
x(1:n) = x(1:n) + alpha * w(1:n)
CALL dmid(n, x, xl, xu)
w(1:n) = wa1(1:n)

RETURN

END SUBROUTINE dprsrch



SUBROUTINE dtrpcg(n, a, adiag, acol_ptr, arow_ind, g, delta,  &
                  l, ldiag, lcol_ptr, lrow_ind,   &
                  tol, stol, itermax, w, iters, info)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:33

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: a(:)
REAL (dp), INTENT(IN)   :: adiag(:)
INTEGER, INTENT(IN)     :: acol_ptr(:)    ! acol_ptr(n+1)
INTEGER, INTENT(IN)     :: arow_ind(:)
REAL (dp), INTENT(IN)   :: g(:)
REAL (dp), INTENT(IN)   :: delta
REAL (dp), INTENT(IN)   :: l(:)
REAL (dp), INTENT(IN)   :: ldiag(:)
INTEGER, INTENT(IN)     :: lcol_ptr(:)    ! lcol_ptr(n+1)
INTEGER, INTENT(IN)     :: lrow_ind(:)
REAL (dp), INTENT(IN)   :: tol
REAL (dp), INTENT(IN)   :: stol
INTEGER, INTENT(IN)     :: itermax
REAL (dp), INTENT(OUT)  :: w(:)
INTEGER, INTENT(OUT)    :: iters
INTEGER, INTENT(OUT)    :: info

!  *********

!  Subroutine dtrpcg

!  Given a sparse symmetric matrix A in compressed column storage,
!  this subroutine uses a preconditioned conjugate gradient method
!  to find an approximate minimizer of the trust region subproblem

!        min { q(s) : || L'*s || <= delta }.

!  where q is the quadratic

!        q(s) = 0.5*s'*A*s + g'*s,

!  A is a symmetric matrix in compressed column storage, L is a lower
!  triangular matrix in compressed column storage, and g is a vector.

!  This subroutine generates the conjugate gradient iterates for
!  the equivalent problem

!        min { Q(w) : || w || <= delta }.

!  where Q is the quadratic defined by

!        Q(w) = q(s),      w = L'*s.

!  Termination occurs if the conjugate gradient iterates leave
!  the trust region, a negative curvature direction is generated,
!  or one of the following two convergence tests is satisfied.

!  Convergence in the original variables:

!        || grad q(s) || <= tol

!  Convergence in the scaled variables:

!        || grad Q(w) || <= stol

!  Note that if w = L'*s, then L*grad Q(w) = grad q(s).

!  The subroutine statement is

!    subroutine dtrcg(n, a, adiag, acol_ptr, arow_ind, g, delta,
!                    l, ldiag, lcol_ptr, lrow_ind,
!                    tol, stol, itermax, w, iters, info)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    a is a REAL (dp) array of dimension nnz.
!      On entry a must contain the strict lower triangular part
!         of A in compressed column storage.
!      On exit a is unchanged.

!    adiag is a REAL (dp) array of dimension n.
!      On entry adiag must contain the diagonal elements of A.
!      On exit adiag is unchanged.

!    acol_ptr is an integer array of dimension n + 1.
!      On entry acol_ptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         acol_ptr(j), ... , acol_ptr(j+1) - 1.
!      On exit acol_ptr is unchanged.

!    arow_ind is an integer array of dimension nnz.
!      On entry arow_ind must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit arow_ind is unchanged.

!    g is a REAL (dp) array of dimension n.
!      On entry g must contain the vector g.
!      On exit g is unchanged.

!    delta is a REAL (dp) variable.
!      On entry delta is the trust region size.
!      On exit delta is unchanged.

! The descriptions below for l, ldiag, lcol_ptr & lrow_ind appear to be wrong.
! All of these quantities must be INPUT; they are not changed.

!    l is a REAL (dp) array of dimension nnz+n*p.
!      On entry l need not be specified.
!      On exit l contains the strict lower triangular part
!         of L in compressed column storage.

!    ldiag is a REAL (dp) array of dimension n.
!      On entry ldiag need not be specified.
!      On exit ldiag contains the diagonal elements of L.

!    lcol_ptr is an integer array of dimension n + 1.
!      On entry lcol_ptr need not be specified.
!      On exit lcol_ptr contains pointers to the columns of L.
!         The nonzeros in column j of L are in the
!         lcol_ptr(j), ... , lcol_ptr(j+1) - 1 positions of l.

!    lrow_ind is an integer array of dimension nnz+n*p.
!      On entry lrow_ind need not be specified.
!      On exit lrow_ind contains row indices for the strict lower
!         triangular part of L in compressed column storage.

!    tol is a REAL (dp) variable.
!      On entry tol specifies the convergence test
!         in the un-scaled variables.
!      On exit tol is unchanged

!    stol is a REAL (dp) variable.
!      On entry stol specifies the convergence test
!         in the scaled variables.
!      On exit stol is unchanged

!    itermax is an integer variable.
!      On entry itermax specifies the limit on the number of
!         conjugate gradient iterations.
!      On exit itermax is unchanged.

!    w is a REAL (dp) array of dimension n.
!      On entry w need not be specified.
!      On exit w contains the final conjugate gradient iterate.

!    iters is an integer variable.
!      On entry iters need not be specified.
!      On exit iters is set to the number of conjugate
!         gradient iterations.

!    info is an integer variable.
!      On entry info need not be specified.
!      On exit info is set as follows:

!          info = 1  Convergence in the original variables.
!                    || grad q(s) || <= tol

!          info = 2  Convergence in the scaled variables.
!                    || grad Q(w) || <= stol

!          info = 3  Negative curvature direction generated.
!                    In this case || w || = delta and a direction
!                    of negative curvature w can be recovered by
!                    solving L'*w = p.

!          info = 4  Conjugate gradient iterates exit the
!                    trust region. In this case || w || = delta.

!          info = 5  Failure to converge within itermax iterations.

!  Subprograms called

!    MINPACK-2  ......  dtrqsol, dstrsol, dssyax

!    Level 1 BLAS  ...  daxpy, dcopy, ddot

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

REAL (dp) :: p(n), q(n), r(n), t(n), z(n)
REAL (dp) :: alpha, beta, ptq, rho, rtr, sigma
REAL (dp) :: rnorm, rnorm0, tnorm

! REAL (dp) :: ddot
! EXTERNAL dtrqsol, dstrsol, dssyax
! EXTERNAL daxpy, dcopy, ddot

!     Initialize the iterate w and the residual r.

w(1:n) = zero

!     Initialize the residual t of grad q to -g.
!     Initialize the residual r of grad Q by solving L*r = -g.
!     Note that t = L*r.

t(1:n) = g(1:n)
t(1:n) = - t(1:n)
r(1:n) = t(1:n)
CALL dstrsol(n, l, ldiag, lcol_ptr, lrow_ind, r, 'N')

!     Initialize the direction p.

p(1:n) = r(1:n)

!     Initialize rho and the norms of r and t.

rho = SUM( r(1:n)**2 )
rnorm0 = SQRT(rho)

!     Exit if g = 0.

IF (rnorm0 == zero) THEN
  info = 1
  RETURN
END IF

iters = 0
DO iters = 1, itermax
  
!        Compute z by solving L'*z = p.
  
  z(1:n) = p(1:n)
  CALL dstrsol(n, l, ldiag, lcol_ptr, lrow_ind, z, 'T')
  
!        Compute q by solving L*q = A*z and save L*q for
!        use in updating the residual t.
  
  CALL dssyax(n, a, adiag, acol_ptr, arow_ind, z, q)
  z(1:n) = q(1:n)
  CALL dstrsol(n, l, ldiag, lcol_ptr, lrow_ind, q, 'N')
  
!        Compute alpha and determine sigma such that the trust region
!        constraint || w + sigma*p || = delta is satisfied.
  
  ptq = DOT_PRODUCT( p(1:n), q(1:n) )
  IF (ptq > zero) THEN
    alpha = rho/ptq
  ELSE
    alpha = zero
  END IF
  CALL dtrqsol(n, w, p, delta, sigma)
  
!        Exit if there is negative curvature or if the
!        iterates exit the trust region.
  
  IF (ptq <= zero .OR. alpha >= sigma) THEN
    w(1:n) = w(1:n) + sigma * p(1:n)
    IF (ptq <= zero) THEN
      info = 3
    ELSE
      info = 4
    END IF
    
    RETURN
    
  END IF
  
!        Update w and the residuals r and t.
!        Note that t = L*r.
  
  w(1:n) = w(1:n) + alpha * p(1:n)
  r(1:n) = r(1:n) - alpha * q(1:n)
  t(1:n) = t(1:n) - alpha * z(1:n)
  
!        Exit if the residual convergence test is satisfied.
  
  rtr = SUM( r(1:n)**2 )
  rnorm = SQRT(rtr)
  tnorm = dnrm2(n, t, 1)
  
  IF (tnorm <= tol) THEN
    info = 1
    RETURN
  END IF
  
  IF (rnorm <= stol) THEN
    info = 2
    RETURN
  END IF
  
!        Compute p = r + beta*p and update rho.
  
  beta = rtr/rho
  p(1:n) = beta * p(1:n)
  p(1:n) = p(1:n) + r(1:n)
  rho = rtr
  
END DO

info = 5

RETURN

END SUBROUTINE dtrpcg



SUBROUTINE dtrqsol(n, x, p, delta, sigma)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(IN)   :: p(:)
REAL (dp), INTENT(IN)   :: delta
REAL (dp), INTENT(OUT)  :: sigma

!  **********

!  Subroutine dtrqsol

!  This subroutine computes the largest (non-negative) solution
!  of the quadratic trust region equation

!        ||x + sigma*p|| = delta.

!  The code is only guaranteed to produce a non-negative solution
!  if ||x|| <= delta, and p != 0. If the trust region equation has
!  no solution, sigma = 0.

!  The subroutine statement ix

!    dtrqsol(n, x, p, delta, sigma)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x must contain the vector x.
!      On exit x is unchanged.

!    p is a REAL (dp) array of dimension n.
!      On entry p must contain the vector p.
!      On exit p is unchanged.

!    delta is a REAL (dp) variable.
!      On entry delta specifies the scalar delta.
!      On exit delta is unchanged.

!    sigma is a REAL (dp) variable.
!      On entry sigma need not be specified.
!      On exit sigma contains the non-negative solution.

!  Subprograms called

!    Level 1 BLAS  ...  ddot

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********


REAL (dp) :: dsq, ptp, ptx, rad, xtx

! REAL (dp) :: ddot

ptx = DOT_PRODUCT( p(1:n), x(1:n) )
ptp = SUM( p(1:n)**2 )
xtx = SUM( x(1:n)**2 )
dsq = delta**2

!     Guard against abnormal cases.

rad = ptx**2 + ptp*(dsq - xtx)
rad = SQRT(MAX(rad, zero))

IF (ptx > zero) THEN
  sigma = (dsq - xtx)/(ptx + rad)
ELSE IF (rad > zero) THEN
  sigma = (rad - ptx)/ptp
ELSE
  sigma = zero
END IF

RETURN

END SUBROUTINE dtrqsol



SUBROUTINE dssyax(n, a, adiag, jptr, indr, x, y)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 14:03:11

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: a(:)
REAL (dp), INTENT(IN)   :: adiag(:)
INTEGER, INTENT(IN)     :: jptr(:)    ! jptr(n+1)
INTEGER, INTENT(IN)     :: indr(:)
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(OUT)  :: y(:)

!  **********

!  Subroutine dssyax

!  This subroutine computes the matrix-vector product y = A*x,
!  where A is a symmetric matrix with the strict lower triangular
!  part in compressed column storage.

!  The subroutine statement is

!    subroutine dssyax(n, a, adiag, jptr, indr, x, y)

!  where

!    n is an integer variable.
!      On entry n is the order of A.
!      On exit n is unchanged.

!    a is a REAL (dp) array of dimension *.
!      On entry a must contain the strict lower triangular part
!         of A in compressed column storage.
!      On exit a is unchanged.

!    adiag is a REAL (dp) array of dimension n.
!      On entry adiag must contain the diagonal elements of A.
!      On exit adiag is unchanged.

!    jptr is an integer array of dimension n + 1.
!      On entry jptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         jptr(j), ... , jptr(j+1) - 1.
!      On exit jptr is unchanged.

!    indr is an integer array of dimension *.
!      On entry indr must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit indr is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x must contain the vector x.
!      On exit x is unchanged.

!    y is a REAL (dp) array of dimension n.
!      On entry y need not be specified.
!      On exit y contains the product A*x.

!  MINPACK-2 Project. May 1998.
!  Argonne National Laboratory.

!  **********


INTEGER   :: i, j, k
REAL (dp) :: sum

y(1:n) = adiag(1:n)*x(1:n)

DO j = 1, n
  sum = zero
  DO i = jptr(j), jptr(j+1) - 1
    k = indr(i)
    sum = sum + a(i)*x(k)
    y(k) = y(k) + a(i)*x(j)
  END DO
  y(j) = y(j) + sum
END DO

RETURN

END SUBROUTINE dssyax



FUNCTION dnrm2 ( n, x, incx) RESULT(fn_val)

!  Euclidean norm of the n-vector stored in x() with storage increment incx .
!  if n <= 0 return with result = 0.
!  if n >= 1 then incx must be >= 1

!  c.l.lawson, 1978 jan 08
!  modified to correct failure to update ix, 1/25/92.
!  modified 3/93 to return if incx <= 0.
!  This version by Alan.Miller @ vic.cmis.csiro.au
!  Latest revision - 22 January 1999

!  four phase method using two built-in constants that are
!  hopefully applicable to all machines.
!      cutlo = maximum of  SQRT(u/eps)  over all known machines.
!      cuthi = minimum of  SQRT(v)      over all known machines.
!  where
!      eps = smallest no. such that eps + 1. > 1.
!      u   = smallest positive no.   (underflow limit)
!      v   = largest  no.            (overflow  limit)

!  brief outline of algorithm..

!  phase 1    scans zero components.
!  move to phase 2 when a component is nonzero and <= cutlo
!  move to phase 3 when a component is > cutlo
!  move to phase 4 when a component is >= cuthi/m
!  where m = n for x() real and m = 2*n for complex.

INTEGER, INTENT(IN)   :: n, incx
REAL (dp), INTENT(IN) :: x(:)
REAL (dp)             :: fn_val

! Local variables
INTEGER     :: i, ix, j, next
REAL (dp)   :: cuthi, cutlo, hitest, sum, xmax

IF(n <= 0 .OR. incx <= 0) THEN
  fn_val = zero
  RETURN
END IF

! Set machine-dependent constants

cutlo = SQRT( TINY(one) / EPSILON(one) )
cuthi = SQRT( HUGE(one) )

next = 1
sum = zero
i = 1
ix = 1
!                                                 begin main loop
20 SELECT CASE (next)
  CASE (1)
     IF( ABS(x(i)) > cutlo) GO TO 85
     next = 2
     xmax = zero
     GO TO 20

  CASE (2)
!                   phase 1.  sum is zero

     IF( x(i) == zero) GO TO 200
     IF( ABS(x(i)) > cutlo) GO TO 85

!                                prepare for phase 2.   x(i) is very small.
     next = 3
     GO TO 105

  CASE (3)
!                   phase 2.  sum is small.
!                             scale to avoid destructive underflow.

     IF( ABS(x(i)) > cutlo ) THEN
!                  prepare for phase 3.

       sum = (sum * xmax) * xmax
       GO TO 85
     END IF

  CASE (4)
     GO TO 110
END SELECT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                     common code for phases 2 and 4.
!                     in phase 4 sum is large.  scale to avoid overflow.

110 IF( ABS(x(i)) <= xmax ) GO TO 115
sum = one + sum * (xmax / x(i))**2
xmax = ABS(x(i))
GO TO 200

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!                   phase 3.  sum is mid-range.  no scaling.

!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)

85 hitest = cuthi / REAL( n, dp )

DO j = ix, n
  IF(ABS(x(i)) >= hitest) GO TO 100
  sum = sum + x(i)**2
  i = i + incx
END DO
fn_val = SQRT( sum )
RETURN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!                                prepare for phase 4.
!                                ABS(x(i)) is very large
100 ix = j
next = 4
sum = (sum / x(i)) / x(i)
!                                Set xmax; large if next = 4, small if next = 3
105 xmax = ABS(x(i))

115 sum = sum + (x(i)/xmax)**2

200 ix = ix + 1
i = i + incx
IF( ix <= n ) GO TO 20

!              end of main loop.

!              compute square root and adjust for scaling.

fn_val = xmax * SQRT(sum)

RETURN
END FUNCTION dnrm2



SUBROUTINE dicfs(n, nnz, a, adiag, acol_ptr, arow_ind,  &
                 l, ldiag, lcol_ptr, lrow_ind,  &
                 p, alpha)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 14:02:36

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nnz
REAL (dp), INTENT(IN)      :: a(:)
REAL (dp), INTENT(IN)      :: adiag(:)
INTEGER, INTENT(IN)        :: acol_ptr(:)   ! acol_ptr(n+1)
INTEGER, INTENT(IN)        :: arow_ind(:)
REAL (dp), INTENT(OUT)     :: l(:)          ! l(nnz+n*p)
REAL (dp), INTENT(OUT)     :: ldiag(:)
INTEGER, INTENT(OUT)       :: lcol_ptr(:)   ! lcol_ptr(n+1)
INTEGER, INTENT(OUT)       :: lrow_ind(:)   ! lrow_ind(nnz+n*p)
INTEGER, INTENT(IN)        :: p
REAL (dp), INTENT(IN OUT)  :: alpha

!  *********

!  Subroutine dicfs

!  Given a symmetric matrix A in compressed column storage, this
!  subroutine computes an incomplete Cholesky factor of A + alpha*D,
!  where alpha is a shift and D is the diagonal matrix with entries
!  set to the l2 norms of the columns of A.

!  The subroutine statement is

!    subroutine dicfs(n, nnz, a, adiag, acol_ptr, arow_ind,
!                     l, ldiag, lcol_ptr, lrow_ind,
!                     p, alpha)

!  where

!    n is an integer variable.
!      On entry n is the order of A.
!      On exit n is unchanged.

!    nnz is an integer variable.
!      On entry nnz is the number of nonzeros in the strict lower
!         triangular part of A.
!      On exit nnz is unchanged.

!    a is a REAL (dp) array of dimension nnz.
!      On entry a must contain the strict lower triangular part
!         of A in compressed column storage.
!      On exit a is unchanged.

!    adiag is a REAL (dp) array of dimension n.
!      On entry adiag must contain the diagonal elements of A.
!      On exit adiag is unchanged.

!    acol_ptr is an integer array of dimension n + 1.
!      On entry acol_ptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         acol_ptr(j), ... , acol_ptr(j+1) - 1.
!      On exit acol_ptr is unchanged.

!    arow_ind is an integer array of dimension nnz.
!      On entry arow_ind must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit arow_ind is unchanged.

!    l is a REAL (dp) array of dimension nnz+n*p.
!      On entry l need not be specified.
!      On exit l contains the strict lower triangular part
!         of L in compressed column storage.

!    ldiag is a REAL (dp) array of dimension n.
!      On entry ldiag need not be specified.
!      On exit ldiag contains the diagonal elements of L.

!    lcol_ptr is an integer array of dimension n + 1.
!      On entry lcol_ptr need not be specified.
!      On exit lcol_ptr contains pointers to the columns of L.
!         The nonzeros in column j of L are in the
!         lcol_ptr(j), ... , lcol_ptr(j+1) - 1 positions of l.

!    lrow_ind is an integer array of dimension nnz+n*p.
!      On entry lrow_ind need not be specified.
!      On exit lrow_ind contains row indices for the strict lower
!         triangular part of L in compressed column storage.

!    p is an integer variable.
!      On entry p specifes the amount of memory available for the
!         incomplete Cholesky factorization.
!      On exit p is unchanged.

!    alpha is a REAL (dp) variable.
!      On entry alpha is the initial guess of the shift.
!      On exit alpha is final shift

!  Subprograms called

!    MINPACK-2  ......  dicf

!  MINPACK-2 Project. October 1998.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER, PARAMETER :: nbmax=3

REAL (dp), PARAMETER :: alpham=1.0D-3, two=2.0_dp
REAL (dp), PARAMETER :: nbfactor=512

REAL (dp) :: wa1(n), wa2(n)
INTEGER   :: i, info, j, k, nb
REAL (dp) :: alphas

! EXTERNAL dicf

!     Compute the l2 norms of the columns of A.

wa1(1:n) = adiag(1:n)**2
DO j = 1, n
  DO i = acol_ptr(j), acol_ptr(j+1)-1
    k = arow_ind(i)
    wa1(j) = wa1(j) + a(i)**2
    wa1(k) = wa1(k) + a(i)**2
  END DO
END DO
wa1(1:n) = SQRT(wa1(1:n))

!     Compute the scaling matrix D.

DO i = 1, n
  IF (wa1(i) > zero) THEN
    wa2(i) = one/SQRT(wa1(i))
  ELSE
    wa2(i) = one
  END IF
END DO

!     Determine a lower bound for the step.

IF (alpha <= zero) THEN
  alphas = alpham
ELSE
  alphas = alpha
END IF

!     Compute the initial shift.

alpha = zero
DO i = 1, n
  IF (adiag(i) == zero) THEN
    alpha = alphas
  ELSE
    alpha = MAX(alpha, -adiag(i)*(wa2(i)**2))
  END IF
END DO
IF (alpha > zero) alpha = MAX(alpha,alphas)

!     Search for an acceptable shift. During the search we decrease
!     the lower bound alphas until we determine a lower bound that
!     is not acceptable. We then increase the shift.
!     The lower bound is decreased by nbfactor at most nbmax times.

nb = 1
DO
  
!        Copy the sparsity structure of A into L.
  
  lcol_ptr(1:n+1) = acol_ptr(1:n+1)
  lrow_ind(1:nnz) = arow_ind(1:nnz)
  
!        Scale A and store in the lower triangular matrix L.
  
  ldiag(1:n) = adiag(1:n)*(wa2(1:n)**2) + alpha
  DO j = 1, n
    DO i = acol_ptr(j), acol_ptr(j+1)-1
      l(i) = a(i)*wa2(j)*wa2(arow_ind(i))
    END DO
  END DO
  
!        Attempt an incomplete factorization.
  
  CALL dicf(n, nnz, l, ldiag, lcol_ptr, lrow_ind, p, info)
  
!        If the factorization exists, then test for termination.
!        Otherwise increment the shift.
  
  IF (info >= 0) THEN
    
!           If the shift is at the lower bound, reduce the shift.
!           Otherwise undo the scaling of L and exit.
    
    IF (alpha == alphas .AND. nb < nbmax) THEN
      alphas = alphas/nbfactor
      alpha = alphas
      nb = nb + 1
    ELSE
      ldiag(1:n) = ldiag(1:n)/wa2(1:n)
      DO j = 1, lcol_ptr(n+1)-1
        l(j) = l(j)/wa2(lrow_ind(j))
      END DO
      RETURN
    END IF
  ELSE
    alpha = MAX(two*alpha, alphas)
  END IF
END DO

RETURN

END SUBROUTINE dicfs



SUBROUTINE dicf(n, nnz, a, diag, col_ptr, row_ind, p, info)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 14:02:25

INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: nnz
REAL (dp), INTENT(IN OUT)  :: a(:)
REAL (dp), INTENT(IN OUT)  :: diag(:)
INTEGER, INTENT(IN OUT)    :: col_ptr(:)    ! col_ptr(n+1)
INTEGER, INTENT(IN OUT)    :: row_ind(:)
INTEGER, INTENT(IN)        :: p
INTEGER, INTENT(OUT)       :: info

!  *********

!  Subroutine dicf

!  Given a sparse symmetric matrix A in compressed row storage,
!  this subroutine computes an incomplete Cholesky factorization.

!  Implementation of dicf is based on the Jones-Plassmann code.
!  Arrays indf and list define the data structure.
!  At the beginning of the computation of the j-th column,

!    For k < j, indf(k) is the index of a for the first
!    nonzero l(i,k) in the k-th column with i >= j.

!    For k < j, list(i) is a pointer to a linked list of column
!    indices k with i = row_ind(indf(k)).

!  For the computation of the j-th column, the array indr records
!  the row indices. Hence, if nlj is the number of nonzeros in the
!  j-th column, then indr(1),...,indr(nlj) are the row indices.
!  Also, for i > j, indf(i) marks the row indices in the j-th
!  column so that indf(i) = 1 if l(i,j) is not zero.

!  The subroutine statement is

!    subroutine dicf(n, nnz, a, diag, col_ptr, row_ind, p, info)

!  where

!    n is an integer variable.
!      On entry n is the order of A.
!      On exit n is unchanged.

!    nnz is an integer variable.
!      On entry nnz is the number of nonzeros in the strict lower
!         triangular part of A.
!      On exit nnz is unchanged.

!    a is a REAL (dp) array of dimension nnz+n*p.
!      On entry the first nnz entries of a must contain the strict
!         lower triangular part of A in compressed column storage.
!      On exit a contains the strict lower triangular part
!         of L in compressed column storage.

!    diag is a REAL (dp) array of dimension n.
!      On entry diag must contain the diagonal elements of A.
!      On exit diag contains the diagonal elements of L.

!    col_ptr is an integer array of dimension n + 1.
!      On entry col_ptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         col_ptr(j), ... , col_ptr(j+1) - 1.
!      On exit col_ptr contains pointers to the columns of L.
!         The nonzeros in column j of L are in the
!         col_ptr(j), ... , col_ptr(j+1) - 1 positions of l.

!    row_ind is an integer array of dimension nnz+n*p.
!      On entry row_ind must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit row_ind contains row indices for the strict lower
!         triangular part of L in compressed column storage.

!    p is an integer variable.
!      On entry p specifes the amount of memory available for the
!         incomplete Cholesky factorization.
!      On exit p is unchanged.

!    info is an integer variable.
!      On entry info need not be specified.
!      On exit info = 0 if the factorization succeeds, and
!         info < 0 if the -info pivot is not positive.

!  Subprograms called

!    MINPACK-2  ......  dsel2, ihsort, insort

!    Level 1 BLAS  ...  daxpy, dcopy, ddot, dnrm2

!  MINPACK-2 Project. May 1998.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER, PARAMETER :: insortf=20

INTEGER   :: i, ip, j, k, kth, l, nlj, newk, np, mlj, npj
INTEGER   :: isj, iej, isk, iek, newisj, newiej
REAL (dp) :: lval
INTEGER   :: indr(n), indf(n), list(n)
REAL (dp) :: w(n)

! EXTERNAL dsel2, ihsort, insort

info = 0
DO j = 1, n
  indf(j) = 0
  list(j) = 0
END DO

!     Make room for L by moving A to the last n*p positions in a.

np = n*p
col_ptr(1:n+1) = col_ptr(1:n+1) + np
npj = np + nnz
DO j = nnz, 1, -1
  row_ind(npj) = row_ind(j)
  a(npj) = a(j)
  npj = npj - 1
END DO

!     Compute the incomplete Cholesky factorization.

isj = col_ptr(1)
col_ptr(1) = 1
DO j = 1, n
  
!        Load column j into the array w. The first and last elements
!        of the j-th column of A are a(isj) and a(iej).
  
  nlj = 0
  iej = col_ptr(j+1) - 1
  DO ip = isj, iej
    i = row_ind(ip)
    w(i) = a(ip)
    nlj = nlj + 1
    indr(nlj) = i
    indf(i) = 1
  END DO
  
!        Exit if the current pivot is not positive.
  
  IF (diag(j) <= zero) THEN
    info = -j
    RETURN
  END IF
  diag(j) = SQRT(diag(j))
  
!        Update column j using the previous columns.
  
  k = list(j)
  DO WHILE (k /= 0)
    isk = indf(k)
    iek = col_ptr(k+1) - 1
    
!           Set lval to l(j,k).
    
    lval = a(isk)
    
!           Update indf and list.
    
    newk = list(k)
    isk = isk + 1
    IF (isk < iek) THEN
      indf(k) = isk
      list(k) = list(row_ind(isk))
      list(row_ind(isk)) = k
    END IF
    k = newk
    
!           Compute the update a(i,i) <- a(i,j) - l(i,k)*l(j,k).
!           In this loop we pick up l(i,k) for k < j and i > j.
    
    DO ip = isk, iek
      i = row_ind(ip)
      IF (indf(i) /= 0) THEN
        w(i) = w(i) - lval*a(ip)
      ELSE
        indf(i) = 1
        nlj = nlj + 1
        indr(nlj) = i
        w(i) = - lval*a(ip)
      END IF
    END DO
  END DO
  
!        Compute the j-th column of L.
  
  DO k = 1, nlj
    w(indr(k)) = w(indr(k))/diag(j)
  END DO
  
!        Set mlj to the number of nonzeros to be retained.
  
  mlj = MIN(iej-isj+1+p, nlj)
  kth = nlj - mlj + 1
  
  IF (nlj >= 1) THEN
    
!           Determine the kth smallest elements in the current
!           column, and hence, the largest mlj elements.
    
    CALL dsel2(nlj, w, indr, kth)
    
!           Sort the row indices of the selected elements.  Insertion
!           sort is used for small arrays, and heap sort for larger
!           arrays. The sorting of the row indices is required so that
!           we can retrieve l(i,k) with i > k from indf(k).
    
    IF (mlj <= insortf) THEN
      CALL insort(mlj, indr(kth:))
    ELSE
      CALL ihsort(mlj, indr(kth:))
    END IF
  END IF
  
!        Store the largest elements in L. The first and last elements
!        of the j-th column of L are a(newisj) and a(newiej).
  
  newisj = col_ptr(j)
  newiej = newisj + mlj -1
  DO k = newisj, newiej
    a(k) = w(indr(k-newisj+kth))
    row_ind(k) = indr(k-newisj+kth)
  END DO
  
!        Update the diagonal elements.
  
  DO k = kth, nlj
    l = indr(k)
    diag(l) = diag(l) - w(l)**2
  END DO
  
!        Update indf and list for the j-th column.
  
  IF (newisj < newiej) THEN
    indf(j) = newisj
    list(j) = list(row_ind(newisj))
    list(row_ind(newisj)) = j
  END IF
  
!        Clear out elements j+1,...,n of the array indf.
  
  DO k = 1, nlj
    indf(indr(k)) = 0
  END DO
  
!        Update isj and col_ptr.
  
  isj = col_ptr(j+1)
  col_ptr(j+1) = newiej + 1
  
END DO

RETURN
END SUBROUTINE dicf



SUBROUTINE dstrsol(n, l, ldiag, jptr, indr, r, task)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 14:03:19

INTEGER, INTENT(IN)            :: n
REAL (dp), INTENT(IN)          :: l(:)
REAL (dp), INTENT(IN)          :: ldiag(:)
INTEGER, INTENT(IN)            :: jptr(:)    ! jptr(n+1)
INTEGER, INTENT(IN)            :: indr(:)
REAL (dp), INTENT(IN OUT)      :: r(:)
CHARACTER (LEN=*), INTENT(IN)  :: task

!  **********

!  Subroutine dstrsol

!  This subroutine solves the triangular systems L*x = r or L'*x = r.

!  The subroutine statement is

!    subroutine dstrsol(n, l, ldiag, jptr, indr, r, task)

!  where

!    n is an integer variable.
!      On entry n is the order of L.
!      On exit n is unchanged.

!    l is a REAL (dp) array of dimension *.
!      On entry l must contain the nonzeros in the strict lower
!         triangular part of L in compressed column storage.
!      On exit l is unchanged.

!    ldiag is a REAL (dp) array of dimension n.
!      On entry ldiag must contain the diagonal elements of L.
!      On exit ldiag is unchanged.

!    jptr is an integer array of dimension n + 1.
!      On entry jptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         jptr(j), ... , jptr(j+1) - 1.
!      On exit jptr is unchanged.

!    indr is an integer array of dimension *.
!      On entry indr must contain row indices for the strict
!         lower triangular part of L in compressed column storage.
!      On exit indr is unchanged.

!    r is a REAL (dp) array of dimension n.
!      On entry r must contain the vector r.
!      On exit r contains the solution vector x.

!    task is a character variable of length 60.
!      On entry
!         task(1:1) = 'N' if we need to solve L*x = r
!         task(1:1) = 'T' if we need to solve L'*x = r
!      On exit task is unchanged.

!  MINPACK-2 Project. May 1998.
!  Argonne National Laboratory.

!  **********

INTEGER   :: i, j, k
REAL (dp) :: temp

!     Solve L*x =r and store the result in r.

IF (task(1:1) == 'N') THEN
  
  DO j = 1, n
    temp = r(j)/ldiag(j)
    DO k = jptr(j), jptr(j+1) - 1
      i = indr(k)
      r(i) = r(i) - l(k)*temp
    END DO
    r(j) = temp
  END DO
  
  RETURN
  
END IF

!     Solve L'*x =r and store the result in r.

IF (task(1:1) == 'T') THEN
  
  r(n) = r(n)/ldiag(n)
  DO j = n - 1, 1, -1
    temp = zero
    DO k = jptr(j), jptr(j+1) - 1
      temp = temp + l(k)*r(indr(k))
    END DO
    r(j) = (r(j) - temp)/ldiag(j)
  END DO
  
  RETURN
  
END IF

RETURN
END SUBROUTINE dstrsol



SUBROUTINE dsel2(n, x, keys, k)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 14:02:58

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
INTEGER, INTENT(IN OUT)    :: keys(:)
INTEGER, INTENT(IN)        :: k

!  **********

!  Subroutine dsel2

!  Given an array x of length n, this subroutine permutes
!  the elements of the array keys so that

!    abs(x(keys(i))) <= abs(x(keys(k))),  1 <= i <= k,
!    abs(x(keys(k))) <= abs(x(keys(i))),  k <= i <= n.

!  In other words, the smallest k elements of x in absolute value are
!  x(keys(i)), i = 1,...,k, and x(keys(k)) is the kth smallest element.

!  The subroutine statement is

!    subroutine dsel2(n,x,keys,k)

!  where

!    n is an integer variable.
!      On entry n is the number of keys.
!      On exit n is unchanged.

!    x is a REAL (dp) array of length n.
!      On entry x is the array to be sorted.
!      On exit x is unchanged.

!    keys is an integer array of length n.
!      On entry keys is the array of indices for x.
!      On exit keys is permuted so that the smallest k elements
!         of x in absolute value are x(keys(i)), i = 1,...,k, and
!         x(keys(k)) is the kth smallest element.

!    k is an integer.
!      On entry k specifes the kth largest element.
!      On exit k is unchanged.

!  MINPACK-2 Project. March 1998.
!  Argonne National Laboratory.
!  William D. Kastak, Chih-Jen Lin, and Jorge J. More'.

!  **********

INTEGER   :: i, l, lc, lp, m, p, p1, p2, p3, u
INTEGER   :: swap
REAL (dp) :: abskeys

IF (n <= 1 .OR. k <= 0 .OR. k > n) RETURN

u = n
l = 1
lc = n
lp = 2*n

!     Start of iteration loop.

DO WHILE (l < u)
  
!        Choose the partition as the median of the elements in
!        positions l+s*(u-l) for s = 0, 0.25, 0.5, 0.75, 1.
!        Move the partition element into position l.
  
  p1 = (u+3*l)/4
  p2 = (u+l)/2
  p3 = (3*u+l)/4
  
!        Order the elements in positions l and p1.
  
  IF (ABS(x(keys(l))) > ABS(x(keys(p1)))) THEN
    swap = keys(l)
    keys(l) = keys(p1)
    keys(p1) = swap
  END IF
  
!        Order the elements in positions p2 and p3.
  
  IF (ABS(x(keys(p2))) > ABS(x(keys(p3)))) THEN
    swap = keys(p2)
    keys(p2) = keys(p3)
    keys(p3) = swap
  END IF
  
!        Swap the larger of the elements in positions p1
!        and p3, with the element in position u, and reorder
!        the first two pairs of elements as necessary.
  
  IF (ABS(x(keys(p3))) > ABS(x(keys(p1)))) THEN
    swap = keys(p3)
    keys(p3) = keys(u)
    keys(u) = swap
    IF (ABS(x(keys(p2))) > ABS(x(keys(p3)))) THEN
      swap = keys(p2)
      keys(p2) = keys(p3)
      keys(p3) = swap
    END IF
  ELSE
    swap = keys(p1)
    keys(p1) = keys(u)
    keys(u) = swap
    IF (ABS(x(keys(l))) > ABS(x(keys(p1)))) THEN
      swap = keys(l)
      keys(l) = keys(p1)
      keys(p1) = swap
    END IF
  END IF
  
!        If we define a(i) = abs(x(keys(i)) for i = 1,..., n, we have
!        permuted keys so that
  
!          a(l) <= a(p1), a(p2) <= a(p3), max(a(p1),a(p3)) <= a(u).
  
!        Find the third largest element of the four remaining
!        elements (the median), and place in position l.
  
  IF (ABS(x(keys(p1))) > ABS(x(keys(p3)))) THEN
    IF (ABS(x(keys(l))) <= ABS(x(keys(p3)))) THEN
      swap = keys(l)
      keys(l) = keys(p3)
      keys(p3) = swap
    END IF
  ELSE
    IF (ABS(x(keys(p2))) <= ABS(x(keys(p1)))) THEN
      swap = keys(l)
      keys(l) = keys(p1)
      keys(p1) = swap
    ELSE
      swap = keys(l)
      keys(l) = keys(p2)
      keys(p2) = swap
    END IF
  END IF
  
!        Partition the array about the element in position l.
  
  m = l
  abskeys = ABS(x(keys(l)))
  DO i = l+1, u
    IF (ABS(x(keys(i))) < abskeys) THEN
      m = m + 1
      swap = keys(m)
      keys(m) = keys(i)
      keys(i) = swap
    END IF
  END DO
  
!        Move the partition element into position m.
  
  swap = keys(l)
  keys(l) = keys(m)
  keys(m) = swap
  
!        Adjust the values of l and u.
  
  IF (k >= m) l = m + 1
  IF (k <= m) u = m - 1
  
!        Check for multiple medians if the length of the subarray
!        has not decreased by 1/3 after two consecutive iterations.
  
  IF (3*(u-l) > 2*lp .AND. k > m) THEN
    
!           Partition the remaining elements into those elements
!           equal to x(m), and those greater than x(m).  Adjust
!           the values of l and u.
    
    p = m
    DO i = m+1, u
      IF (ABS(x(keys(i))) == ABS(x(keys(m)))) THEN
        p = p + 1
        swap = keys(p)
        keys(p) = keys(i)
        keys(i) = swap
      END IF
    END DO
    l = p + 1
    IF (k <= p) u = p - 1
  END IF
  
!        Update the length indicators for the subarray.
  
  lp = lc
  lc = u-l
  
END DO

RETURN

END SUBROUTINE dsel2



SUBROUTINE ihsort(n, keys)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 14:03:58

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN OUT)  :: keys(:)

!  **********

!  Subroutine ihsort

!  Given an integer array keys of length n, this subroutine uses
!  a heap sort to sort the keys in increasing order.

!  This subroutine is a minor modification of code written by
!  Mark Jones and Paul Plassmann.

!  The subroutine statement is

!    subroutine ihsort(n, keys)

!  where

!    n is an integer variable.
!      On entry n is the number of keys.
!      On exit n is unchanged.

!    keys is an integer array of length n.
!      On entry keys is the array to be sorted.
!      On exit keys is permuted to increasing order.

!  MINPACK-2 Project. March 1998.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER :: k, m, lheap, rheap, mid
INTEGER :: x

IF (n <= 1) RETURN

!     Build the heap.

mid = n/2
DO k = mid, 1, -1
  x = keys(k)
  lheap = k
  rheap = n
  m = lheap*2
  DO WHILE (m <= rheap)
    IF (m < rheap) THEN
      IF (keys(m) < keys(m+1)) m = m + 1
    END IF
    IF (x >= keys(m)) THEN
      m = rheap + 1
    ELSE
      keys(lheap) = keys(m)
      lheap = m
      m = 2*lheap
    END IF
  END DO
  keys(lheap) = x
END DO

!     Sort the heap.

DO k = n, 2, -1
  x = keys(k)
  keys(k) = keys(1)
  lheap = 1
  rheap = k-1
  m = 2
  DO WHILE (m <= rheap)
    IF (m < rheap) THEN
      IF (keys(m) < keys(m+1)) m = m+1
    END IF
    IF (x >= keys(m)) THEN
      m = rheap + 1
    ELSE
      keys(lheap) = keys(m)
      lheap = m
      m = 2*lheap
    END IF
  END DO
  keys(lheap) = x
END DO

RETURN

END SUBROUTINE ihsort



SUBROUTINE insort(n, keys)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 14:04:03

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN OUT)  :: keys(:)

!  **********

!  Subroutine insort

!  Given an integer array keys of length n, this subroutine uses
!  an insertion sort to sort the keys in increasing order.

!  The subroutine statement is

!    subroutine insort(n, keys)

!  where

!    n is an integer variable.
!      On entry n is the number of keys.
!      On exit n is unchanged.

!    keys is an integer array of length n.
!      On entry keys is the array to be sorted.
!      On exit keys is permuted to increasing order.

!  MINPACK-2 Project. March 1998.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER :: i, j, ind

DO j = 2, n
  ind = keys(j)
  i = j - 1
  DO WHILE (i > 0 .AND. keys(i) > ind)
    keys(i+1) = keys(i)
    i = i - 1
  END DO
  keys(i+1) = ind
END DO

RETURN
END SUBROUTINE insort

END MODULE tron
