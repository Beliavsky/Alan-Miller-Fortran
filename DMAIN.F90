PROGRAM dmain
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:45:25
! Latest version - 5 July 1999

USE tron
USE coloring
IMPLICIT NONE

INTEGER, PARAMETER  :: nmax=90000, nnzmax=10*nmax
INTEGER, PARAMETER  :: nread=1, nwrite=2
!     **********

!     Driver for bound-constrained problems.

!     Subprograms called

!       USER ........... dminfg, dminhs, dminsp, dminxb

!       MINPACK-2 ...... dtron, dgpnrm2, dsphesd, dmid, dtimer

!       Level 1 BLAS ... dnrm2

!     MINPACK-2 Project. March 1999.
!     Argonne National Laboratory.
!     Chih-Jen Lin and Jorge J. More'.

!     **********
REAL (dp), PARAMETER :: zero=0.0_dp, one=1.0_dp

CHARACTER (LEN=60) :: task

LOGICAL   :: search
INTEGER   :: arow_ind(nnzmax), acol_ind(nnzmax), brow_ind(nnzmax)
INTEGER   :: arow_ptr(nmax+1), acol_ptr(nmax+1), bcol_ptr(nmax+1),  &
             lrow_ind(nnzmax), lcol_ptr(nmax+1)
INTEGER   :: n, itermax
REAL (dp) :: f, delta
REAL (dp) :: x(nmax), xl(nmax), xu(nmax), g(nmax), hs(nmax)
REAL (dp) :: adiag(nmax), bdiag(nmax), ldiag(nmax)
REAL (dp) :: a(nnzmax), b(nnzmax), l(nnzmax)

!     Tolerances.

REAL (dp) :: cgtol, frtol, fatol, fmin

!     Summary information.

INTEGER   :: iterscg, nbind, nfev, nfree, ngev, nhev, nhsev

!     Evaluation of the Hessian matrix.

INTEGER   :: maxgrp, numgrp
INTEGER   :: listp(nmax), ngrp(nmax)
REAL (dp) :: y(nmax), eta(nmax)

!     Test problems.

CHARACTER (LEN=6) :: prob
CHARACTER (LEN=2) :: ch
INTEGER           :: i, info, j, maxfev, nnz, nx, ny
REAL (dp)         :: gnorm, gnorm0, par, gtol

!     Timing.

REAL :: ttimes, ttimef, times, timef
REAL :: fgtime, htime, ttime

OPEN (nread, FILE='tron.dat', STATUS='OLD')
OPEN (nwrite, FILE='tron.inf')

DO
  
  READ (nread,*) prob, n, nx, ny, par
  WRITE (*,*)    prob, n, nx, ny, par
  
  IF (prob(1:4) == 'STOP') THEN
    CLOSE(nread)
    STOP
  END IF
  
!        Generate the initial point and project into [xl,xu].
  
  ch = 'XS'
  CALL dminfg(n, nx, ny, x, f, g, ch, prob, par)
  CALL dminxb(n, nx, ny, xl, xu, prob)
  CALL dmid(n, x, xl, xu)
  
!        Initialize variables.
  
  nfev = 0
  ngev = 0
  nhev = 0
  nhsev = 0
  fgtime = zero
  htime = zero
  
!        Set parameters.
  
  itermax = n
  maxfev = 1000
  fatol = zero
  frtol = 1.d-12
  fmin = -1.0D+32
  cgtol = 0.1_dp
  gtol = 1.0D-5
  
  CALL dtimer(ttimes)
  
!        Calculate the sparsity pattern.
  
  CALL dminsp(n, nx, ny, nnz, arow_ind, acol_ind, prob)
  CALL dsetsp(n, nnz, arow_ind, acol_ind, acol_ptr, arow_ptr, 1,  &
              info, listp, ngrp, maxgrp)
  IF (info <= 0) THEN
    WRITE (nwrite,*) 'ERROR: INFO IN SUBROUTINE SETSP IS ', info
    STOP
  END IF
  
!        Start the iteration.
  
  task = 'START'
  search = .TRUE.
  DO WHILE (search)
    
!           Function evaluation.
    
    IF (task == 'F' .OR. task == 'START') THEN
      CALL dtimer(times)
      ch = 'F'
      CALL dminfg(n, nx, ny, x, f, g, ch, prob, par)
      nfev = nfev + 1
      CALL dtimer(timef)
      fgtime = fgtime + (timef - times)
    END IF
    
!           Evaluate the gradient and the Hessian matrix.
    
    IF (task == 'GH' .OR. task == 'START') THEN
      CALL dtimer(times)
      ch = 'G'
      CALL dminfg(n, nx, ny, x, f, g, ch, prob, par)
      ngev = ngev + 1
      CALL dtimer(timef)
      fgtime = fgtime + (timef - times)
      
!              Evaluate the Hessian matrix.
      
      CALL dtimer(times)
      DO i = 1, n
        hs(i) = zero
        eta(i) = one
      END DO
      DO numgrp = 1, maxgrp
        DO j = 1, n
          IF (ngrp(j) == numgrp) hs(j) = one
        END DO
        CALL dminhs(n, nx, ny, x, hs, y, prob, par)
        CALL dsphesd(n, arow_ind, acol_ind, arow_ptr, acol_ptr,  &
                     listp, ngrp, maxgrp, numgrp, eta, y, a, adiag)
        DO j = 1, n
          IF (ngrp(j) == numgrp) hs(j) = zero
        END DO
      END DO
      nhev = nhev + 1
      nhsev = nhsev + maxgrp
      CALL dtimer(timef)
      htime = htime + (timef - times)
      
    END IF
    
!           Initialize the trust region bound.
    
    IF (task == 'START') THEN
      gnorm0 = dnrm2(n, g, 1)
      delta = dnrm2(n, g, 1)
    END IF
    
!           Stopping criteria.
    
    IF (task == 'GH' .OR. task == 'START') THEN
      gnorm = dgpnrm2(n, x, xl, xu, g)
      IF (gnorm <= gtol*gnorm0) THEN
        search = .false.
        task = 'CONVERGENCE: GTOL TEST SATISFIED'
      END IF
    END IF
    
    IF (nfev > maxfev) THEN
      search = .false.
      task = 'ERROR: NFEV > MAXFEV'
    END IF
    
!           Call the optimizer.
    
    IF (search) THEN
      CALL dtron(n, x, xl, xu, f, g, a, adiag, acol_ptr, arow_ind,  &
                 frtol, fatol, fmin, cgtol, itermax, delta, task,   &
                 b, bdiag, bcol_ptr, brow_ind,  &
                 l, ldiag, lcol_ptr, lrow_ind, iterscg)
    END IF
    
!           Exit search if the algorithm has converged.
    
    IF (task(1:4) == 'CONV') search = .false.
    
  END DO
  
  CALL dtimer(ttimef)
  
  WRITE (*, *) task
  
!        Summary information.
  
  nfree = 0
  nbind = 0
  DO i = 1, n
    IF (xl(i) < x(i) .AND. x(i) < xu(i)) THEN
      nfree = nfree + 1
    ELSE IF ((x(i) == xl(i) .AND. g(i) >= zero) .OR.  &
             (x(i) == xu(i) .AND. g(i) <= zero) .OR. (xl(i) == xu(i))) THEN
      nbind = nbind + 1
    END IF
  END DO
  
  ch = 'FG'
  CALL dminfg(n, nx, ny, x, f, g, ch, prob, par)
  gnorm = dgpnrm2(n, x, xl, xu, g)
  WRITE (nwrite,1000) prob, n, maxgrp, DBLE(nnz)/n, nfree, n-nfree, nbind,  &
                      nfev, ngev, nhev, nhsev, iterscg, f, gnorm
  
!        Timing information.
  
  ttime = ttimef - ttimes
  fgtime = 100. * fgtime/ttime
  htime  = 100. * htime/ttime
  
  WRITE (nwrite,2000) ttime, fgtime, htime, task
  
END DO
STOP

1000 FORMAT (' Problem ',  a6,                                    //  &
             ' Number of variables                         ', i12/   &
             ' Number of coloring groups                   ', i12/   &
             ' Average number of nonzeros in the strictly  '    /   &
             ' lower triangular part of the Hessian matrix ', f12.2/   &
             ' Number of free variables                    ', i12/   &
             ' Number of active variables                  ', i12/   &
             ' Number of binding variables                 ', i12/   &
             ' Number of function evaluations              ', i12/   &
             ' Number of gradient evaluations              ', i12/   &
             ' Number of Hessian evaluations               ', i12/   &
             ' Number of Hessian-vector evaluations        ', i12/   &
             ' Number of conjugate gradient iterations     ', i12 //   &
             ' Function value at final iterate          '   , g15.8/   &
             ' Projected gradient at final iterate      '   , g15.3 /)

2000 FORMAT (' Total execution time                        ', f12.2/   &
             ' Percentage in function evaluations          ', f12.1/   &
             ' Percentage in Hessian evaluations           ', f12.1 //   &
             ' Exit message     '                           , a60 /)

CONTAINS


SUBROUTINE dtimer(time)

REAL, INTENT(OUT) :: time

CALL CPU_TIME(time)
RETURN

END SUBROUTINE dtimer



SUBROUTINE dminfg(n, nx, ny, x, f, g, task, prob, par)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)                :: n
INTEGER, INTENT(IN)                :: nx
INTEGER, INTENT(IN)                :: ny
REAL (dp), INTENT(IN OUT)          :: x(:)
REAL (dp), INTENT(OUT)             :: f
REAL (dp), INTENT(OUT)             :: g(:)
CHARACTER (LEN=*), INTENT(IN OUT)  :: task
CHARACTER (LEN=*), INTENT(IN OUT)  :: prob
REAL (dp), INTENT(IN)              :: par

!  *********

!  Subroutine dminfg

!  This subroutine computes the function and gradient for
!  the minimization problem from the MINPACK-2 test problem
!  collection specified by the character variable prob.

!  The subroutine statement is

!    subroutine dminfg(n, nx, ny, x, f, g, task, prob, par)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first coordinate
!         direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction. If the problem is formulated in
!         one spatial dimension, ny = 1.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is unchanged.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    g is a REAL (dp) array of dimension n.
!      On entry g need not be specified.
!      On exit g contains the gradient evaluated at x if task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at x.
!          'G'     Evaluate the gradient vector at x.
!          'FG'    Evaluate the function and the gradient at x.
!          'XS'    Set x to the standard starting point xs.

!      On exit task may be changed.
!         task can be changed in routine dljcfg.

!    prob is a character*6 variable.
!      On entry prob specifies the problem.
!      On exit prob is set to 'ERROR' if prob is not an
!         acceptable problem name.  Otherwise prob is unchanged.

!    par is a REAL (dp) variable.
!      On entry par specifies a probem-dependent parameter.
!      On exit par is unchanged.

!  Subprograms called

!    MINPACK-2 ... deptfg, dgl1fg, dgl2fg, dmsafg, dmsabc,
!                  dljcfg, dodcfg, dpjbfg, dsscfg

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Brett M. Averick and Jorge J. More'.

!  **********

REAL (dp), PARAMETER   :: ten=10.0_dp
REAL (dp), ALLOCATABLE :: bottom(:), top(:), left(:), right(:)

! EXTERNAL deptfg, dgl1fg, dgl2fg, dmsafg, dmsabc, dljcfg, dodcfg, dpjbfg,  &
!          dsscfg

!     Select a problem.

SELECT CASE ( prob(1:4) )
  CASE ( 'DEPT' )
    CALL deptfg(nx, ny, x, f, g, task, par)
  CASE( 'DPJB' )
    CALL dpjbfg(nx, ny, x, f, g, task, par, ten)
  CASE( 'DMSA' )
    ALLOCATE( bottom(nx+2), top(nx+2), left(ny+2), right(ny+2) )
    CALL dmsabc(nx, ny, bottom, top, left, right)
    CALL dmsafg(nx, ny, x, f, g, task, bottom, top, left, right)
    DEALLOCATE( bottom, top, left, right )
  CASE( 'DODC' )
    CALL dodcfg(nx, ny, x, f, g, task, par)
  CASE( 'DSSC' )
    CALL dsscfg(nx, ny, x, f, g, task, par)
  CASE( 'DGL1' )
    CALL dgl1fg(n, x, f, g, task, par)
  CASE( 'DGL2' )
    CALL dgl2fg(nx, ny, x, f, g, task, INT(par))
  CASE( 'DLJ2' )
    CALL dljcfg(n, x, f, g, task, 2, n/2)
  CASE( 'DLJ3' )
    CALL dljcfg(n, x, f, g, task, 3, n/3)
  CASE DEFAULT
    prob = 'ERROR'
END SELECT

RETURN
END SUBROUTINE dminfg



SUBROUTINE dminhs(n, nx, ny, x, s, hs, prob, par)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)                :: n
INTEGER, INTENT(IN)                :: nx
INTEGER, INTENT(IN)                :: ny
REAL (dp), INTENT(IN OUT)          :: x(:)
REAL (dp), INTENT(IN OUT)          :: s(:)
REAL (dp), INTENT(OUT)             :: hs(:)
CHARACTER (LEN=*), INTENT(IN OUT)  :: prob
REAL (dp), INTENT(IN OUT)          :: par

!  **********

!  Subroutine dminhs

!  This subroutine computes the Hessian-vector product for
!  the minimization problem from the MINPACK-2 test problem
!  collection specified by the character variable prob.

!  The subroutine statement is

!    subroutine dminhs(n, nx, ny, x, s, hs, prob, par)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction. If the problem is formulated in
!         one spatial dimension, ny = 1.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x may be changed.
!         Routine dgl2hs can change x & s.

!    s is a REAL (dp) array of dimension n.
!      On entry s specifies a vector s.
!      On exit s may be changed.

!    hs is a REAL (dp) array of dimension n.
!      On entry hs need not be specified.
!      On exit hs contains the product H*s where H is the Hessian matrix at x.

!    prob is a character*6 variable.
!      On entry prob specifies the problem.
!      On exit prob is set to 'ERROR' if prob is not an acceptable problem name.
!         Otherwise prob is unchanged.

!    par is a REAL (dp) variable.
!      On entry par specifies a problem-dependent parameter.
!      On exit par is unchanged.

!  Subprograms called

!    MINPACK-2 ... depths, dgl1hs, dgl2hs, dmsahs, dmsabc,
!                  dodchs, dpjbhs, dsschs

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Brett M. Averick and Jorge J. More'.

!  **********

REAL (dp), PARAMETER   :: b=10.0_dp
REAL (dp), ALLOCATABLE :: bottom(:), top(:), left(:), right(:)

! EXTERNAL depths, dgl1hs, dgl2hs, dmsahs, dmsabc, dodchs, dpjbhs, dsschs

!     Select a problem.

SELECT CASE ( prob(1:4) )
  CASE ( 'DEPT' )
    CALL depths(nx, ny, s, hs)
  CASE( 'DGL1' )
    CALL dgl1hs(n, x, s, hs, par)
  CASE( 'DGL2' )
    CALL dgl2hs(nx, ny, x, s, hs, INT(par))
  CASE( 'DMSA' )
    ALLOCATE( bottom(nx+2), top(nx+2), left(ny+2), right(ny+2) )
    CALL dmsabc(nx, ny, bottom, top, left, right)
    CALL dmsahs(nx, ny, x, s, hs, bottom, top, left, right)
    DEALLOCATE( bottom, top, left, right )
  CASE( 'DODC' )
    CALL dodchs(nx, ny, x, s, hs, par)
  CASE( 'DPJB' )
    CALL dpjbhs(nx, ny, s, hs, par, b)
  CASE( 'DSSC' )
    CALL dsschs(nx, ny, x, s, hs, par)
  CASE DEFAULT
    prob = 'ERROR'
END SELECT

RETURN
END SUBROUTINE dminhs



SUBROUTINE dminsp(n, nx, ny, nnz, indrow, indcol, prob)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

INTEGER, INTENT(IN)                :: n
INTEGER, INTENT(IN)                :: nx
INTEGER, INTENT(IN)                :: ny
INTEGER, INTENT(OUT)               :: nnz
INTEGER, INTENT(OUT)               :: indrow(:)
INTEGER, INTENT(OUT)               :: indcol(:)
CHARACTER (LEN=*), INTENT(IN OUT)  :: prob

!  *********

!  Subroutine dminsp

!  This subroutine computes the sparsity pattern for
!  the minimization problem from the MINPACK-2 test problem
!  collection specified by the character variable prob.

!  The subroutine statement is

!    dminsp(n, nx, ny, nnz, indrow, indcol, prob)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction. If the problem is formulated in
!         one spatial dimension, ny = 1.
!      On exit ny is unchanged.

!    nnz is an integer variable.
!      On entry nnz need not be specified.
!      On exit nnz is set to the number of nonzeros in the
!         lower triangle of the Hessian matrix.

!    indrow is an integer array of dimension at least nnz.
!      On entry indrow need not be specified.
!      On exit indrow contains the row indices of the nonzeros
!        in the lower triangle of the Hessian matrix.

!    indcol is an integer array of dimension at least nnz.
!      On entry indcol need not be specified.
!      On exit indcol contains the column indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!    prob is a character*6 variable.
!      On entry prob specifies the problem.
!      On exit prob is set to 'ERROR' if prob is not an
!         acceptable problem name. Otherwise prob is unchanged.

!  Subprograms called

!    MINPACK-2 ... deptsp, dgl1sp, dgl2sp, dmsasp, dmsabc,
!                  dodcsp, dpjbsp, dsscsp

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Brett M. Averick and Jorge J. More'.

!  **********

! EXTERNAL deptsp, dgl1sp, dgl2sp, dmsasp, dodcsp, dpjbsp, dsscsp

!     Select a problem.

SELECT CASE ( prob(1:4) )
  CASE ( 'DEPT' )
    CALL deptsp(nx, ny, nnz, indrow, indcol)
  CASE( 'DGL1' )
    CALL dgl1sp(n, nnz, indrow, indcol)
  CASE( 'DGL2' )
    CALL dgl2sp(nx, ny, nnz, indrow, indcol)
  CASE( 'DMSA' )
    CALL dmsasp(nx, ny, nnz, indrow, indcol)
  CASE( 'DODC' )
    CALL dodcsp(nx, ny, nnz, indrow, indcol)
  CASE( 'DPJB' )
    CALL dpjbsp(nx, ny, nnz, indrow, indcol)
  CASE( 'DSSC' )
    CALL dsscsp(nx, ny, nnz, indrow, indcol)
  CASE DEFAULT
    prob = 'ERROR'
END SELECT

RETURN
END SUBROUTINE dminsp



SUBROUTINE dminxb(n, nx, ny, xl, xu, prob)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:11

INTEGER, INTENT(IN)            :: n
INTEGER, INTENT(IN)            :: nx
INTEGER, INTENT(IN)            :: ny
REAL (dp), INTENT(OUT)         :: xl(:)
REAL (dp), INTENT(OUT)         :: xu(:)
CHARACTER (LEN=*), INTENT(IN)  :: prob

!  *********

!  Subroutine dminxb

!  This subroutine computes the lower and upper bounds for
!  the minimization problem from the MINPACK-2 test problem
!  collection specified by the character variable prob.

!  The subroutine statement is

!    dminxb(n, nx, ny, xl, xu, prob)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction. If the problem is formulated in
!         one spatial dimension, ny = 1.
!      On exit ny is unchanged.

!    xl is a REAL (dp) array of dimension n.
!      On entry xl need not be specified.
!      On exit xl is the vector of lower bounds.

!    xu is a REAL (dp) array of dimension n.
!      On entry xu need not be specified.
!      On exit xu is the vector of upper bounds.

!    prob is a character*6 variable.
!      On entry prob specifies the problem.
!      On exit prob is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Brett M. Averick and Jorge J. More'.

!  **********

REAL (dp), PARAMETER :: p0001=1.0D-4, p001=1.0D-3, p01=1.0D-2, p1=1.0D-1

REAL (dp), PARAMETER :: oned2=1.0D2, xbmax=1.0D20

INTEGER   :: i, j, k
REAL (dp) :: hx, hy
REAL (dp) :: tmp

!     Select a problem.

IF (prob(1:4) == 'DEPT') THEN
  hx = one/DBLE(nx+1)
  hy = one/DBLE(ny+1)
  DO j = 1, ny
    tmp = DBLE(MIN(j,ny-j+1))*hy
    DO i = 1, nx
      k = nx*(j-1) + i
      xu(k) = MIN(DBLE(MIN(i,nx-i+1))*hx,tmp)
      xl(k) = -xu(k)
    END DO
  END DO
ELSE IF (prob(1:5) == 'DMSA1') THEN
  DO i = 1, n
    xl(i) = -4*p1
    xu(i) =  4*p1
  END DO
ELSE IF (prob(1:5) == 'DMSA2') THEN
  DO i = 1, n
    xl(i) = -2*p1
    xu(i) =  2*p1
  END DO
ELSE IF (prob(1:5) == 'DMSA3') THEN
  DO i = 1, n
    xl(i) = -p1
    xu(i) =  p1
  END DO
ELSE IF (prob(1:4) == 'DPJB') THEN
  DO i = 1, n
    xl(i) = zero
    xu(i) = oned2
  END DO
ELSE IF (prob(1:5) == 'DSSC1') THEN
  DO i = 1, n
    xl(i) = p1
    xu(i) = one
  END DO
ELSE IF (prob(1:5) == 'DSSC2') THEN
  DO i = 1, n
    xl(i) = p01
    xu(i) = one
  END DO
ELSE IF (prob(1:5) == 'DSSC3') THEN
  DO i = 1, n
    xl(i) = p001
    xu(i) = one
  END DO
ELSE IF (prob(1:5) == 'DSSC4') THEN
  DO i = 1, n
    xl(i) = p0001
    xu(i) = one
  END DO
ELSE
  DO i = 1, n
    xl(i) = -xbmax
    xu(i) =  xbmax
  END DO
END IF

RETURN
END SUBROUTINE dminxb



FUNCTION dgpnrm2(n, x, xl, xu, g) RESULT(norm)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(:)
REAL (dp), INTENT(IN)  :: xl(:)
REAL (dp), INTENT(IN)  :: xu(:)
REAL (dp), INTENT(IN)  :: g(:)
REAL (dp)              :: norm

!  **********

!  Function dgpnrm2

!  This function computes the Euclidean norm of the projected gradient at x.

!  The function statement is

!    function dgpnrm2(n, x, xl, xu, g)

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

!    g is a REAL (dp) array of dimension n.
!      On entry g specifies the gradient g.
!      On exit g is unchanged.

!  MINPACK-2 Project. May 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER :: i

norm = zero
DO i = 1, n
  IF (xl(i) /= xu(i)) THEN
    IF (x(i) == xl(i)) THEN
      norm = norm + MIN(g(i), zero)**2
    ELSE IF (x(i) == xu(i)) THEN
      norm = norm + MAX(g(i), zero)**2
    ELSE
      norm = norm + g(i)**2
    END IF
  END IF
END DO
norm = SQRT(norm)

RETURN

END FUNCTION dgpnrm2



SUBROUTINE dsphesd(n, row_ind, col_ind, row_ptr, col_ptr, listp, ngrp,  &
                   maxgrp, numgrp, eta, fhesd, fhes, diag)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(IN OUT)  :: row_ind(:)
INTEGER, INTENT(IN OUT)  :: col_ind(:)
INTEGER, INTENT(IN)      :: row_ptr(:)    ! row_ptr(n+1)
INTEGER, INTENT(IN)      :: col_ptr(:)    ! col_ptr(n+1)
INTEGER, INTENT(OUT)     :: listp(:)
INTEGER, INTENT(IN)      :: ngrp(:)
INTEGER, INTENT(IN)      :: maxgrp
INTEGER, INTENT(IN)      :: numgrp
REAL (dp), INTENT(IN)    :: eta(:)
REAL (dp), INTENT(IN)    :: fhesd(:)
REAL (dp), INTENT(OUT)   :: fhes(:)
REAL (dp), INTENT(OUT)   :: diag(:)

!  **********

!  Subroutine dsphesd

!  This subroutine computes an approximation to the (symmetric) Hessian matrix
!  of a function by a substitution method.
!  The lower triangular part of the approximation is stored in compressed
!  column storage.

!  This subroutine requires a symmetric permutation of the Hessian matrix and
!  a partition of the columns of the Hessian matrix consistent with the
!  determination of the Hessian matrix by a lower triangular substitution
!  method.   This information can be provided by subroutine dssm.

!  The symmetric permutation of the Hessian matrix is defined by the array
!  listp.  This array is only used internally.

!  The partition of the Hessian matrix is defined by the array ngrp by setting
!  ngrp(j) to the group number of column j.
!  The user must provide an approximation to the columns of the Hessian matrix
!  in each group by specifying a difference parameter vector eta and an
!  approximation to A*d where A is the Hessian matrix and the vector d is
!  defined by the following section of code.

!        do j = 1, n
!           d(j) = 0.0
!           if (ngrp(j) .eq. numgrp) d(j) = eta(j)
!        end do

!  In the above code numgrp is a group number and eta(j) is the
!  difference parameter used to approximate column j of the
!  Hessian matrix. Suitable values for eta(j) must be provided.

!  As mentioned above, an approximation to A*d must be provided.
!  For example, if grad f(x) is the gradient of the function at x, then

!        grad f(x+d) - grad f(x)

!  corresponds to the forward difference approximation.

!  The lower triangular substitution method requires that the approximations
!  to A*d for all the groups be stored in special locations of the array fhes.
!  This is done by calls with numgrp = 1, 2, ... ,maxgrp. On the call with
!  numgrp = maxgrp, the array fhes is overwritten with the approximation to the
!  lower triangular part of the Hessian matrix.

!  The subroutine statement is

!    subroutine dsphesd(n, row_ind, col_ind, row_ptr, col_ptr, listp, ngrp,
!                       maxgrp, numgrp, eta, fhesd, fhes, diag)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    row_ind is an integer array of dimension nnz.
!      On entry row_ind must contain row indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit row_ind is unchanged.  NOT TRUE!

!    col_ind is an integer array of dimension nnz.
!      On entry col_ind must contain column indices for the strict
!         lower triangular part of A in compressed column storage.
!      On exit col_ind is unchanged.  NOT TRUE!

!    row_ptr is an integer array of dimension n + 1.
!      On entry row_ptr must contain pointers to the rows of A.
!         The nonzeros in row j of A must be in positions
!         row_ptr(j), ... , row_ptr(j+1) - 1 of col_ind.
!      On exit row_ptr is unchanged.

!    col_ptr is an integer array of dimension n + 1.
!      On entry col_ptr must contain pointers to the columns of A.
!         The nonzeros in column j of A must be in positions
!         col_ptr(j), ... , col_ptr(j+1) - 1 of row_ind.
!      On exit col_ptr is unchanged.

!    listp is an integer array of length n.
!      On input listp need not be specified.
!      On output listp specifies a permutation of the matrix.
!         Element (i,j) of the matrix is the (listp(i), listp(j))
!         element of the permuted matrix.

!    ngrp is an integer array of length n.
!      On entry ngrp need not be specified.
!      On exit ngrp specifies the partition of the columns of A.
!         Column j belongs to group ngrp(j).

!    maxgrp is an integer variable.
!      On entry maxgrp need not be specified.  WRONG.
!      On exit maxgrp specifies the number of groups in
!         the partition of the columns of A.

!    numgrp is an integer variable.
!      On input numgrp must be set to a group number.
!      On output numgrp is unchanged.

!    eta is a REAL (dp) variable.
!      On input eta is the difference parameter vector.
!      On output eta is unchanged.

!    fhesd is a REAL (dp) array of length n.
!      On input fhesd contains an approximation to A*d, where A
!        is the Hessian matrix and d is the difference vector for group numgrp.
!      On output fhesd is unchanged.

!    fhes is a REAL (dp) array of length nnz.
!      On input fhes need not be specified.
!      On output fhes is overwritten. When numgrp = maxgrp the
!         array fhes contains an approximation to the Hessian matrix
!         in compressed column storage. The elements in column j of
!         the strict lower triangular part of the Hessian matrix are

!            fhes(k), k = col_ptr(j),...,col_ptr(j+1)-1,

!         and the row indices for these elements are

!            row_ind(k), k = col_ptr(j),...,col_ptr(j+1)-1.

!    diag is a REAL (dp) array of length n.
!      On input diag need not be specified.
!      On output diag is overwritten.  When numgrp = maxgrp the array diag
!         contains the diagonal elements of an approximation to the Hessian
!         matrix.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER   :: i, ip, irow, j, jp, k, l, numg, numl
INTEGER   :: iwk(n)
REAL (dp) :: sum

!     Store the i-th element of gradient difference fhesd corresponding to
!     group numgrp if there is a position (i,j) such that ngrp(j) = numgrp and
!     (i,j) is mapped onto the lower triangular part of the permuted matrix.

DO j = 1, n
  IF (ngrp(j) == numgrp) THEN
    diag(j) = fhesd(j)/eta(j)
    numl = listp(j)
    DO ip = row_ptr(j), row_ptr(j+1)-1
      i = col_ind(ip)
      IF (listp(i) > numl) THEN
        DO jp = col_ptr(i), col_ptr(i+1)-1
          IF (row_ind(jp) == j) THEN
            fhes(jp) = fhesd(i)
            EXIT
          END IF
        END DO
      END IF
    END DO
    DO jp = col_ptr(j), col_ptr(j+1)-1
      i = row_ind(jp)
      IF (listp(i) >= numl) fhes(jp) = fhesd(i)
    END DO
  END IF
END DO

!     Exit if this is not the last group.

IF (numgrp < maxgrp) RETURN

!     Mark all column indices j such that (i,j) is mapped onto
!     the lower triangular part of the permuted matrix.

DO i = 1, n
  numl = listp(i)
  DO ip = row_ptr(i), row_ptr(i+1)-1
    j = col_ind(ip)
    IF (numl >= listp(j)) col_ind(ip) = -col_ind(ip)
  END DO
  DO jp = col_ptr(i), col_ptr(i+1)-1
    j = row_ind(jp)
    IF (numl > listp(j)) row_ind(jp) = -row_ind(jp)
  END DO
END DO

!     Invert the array listp.

DO j = 1, n
  iwk(listp(j)) = j
END DO
listp(1:n) = iwk(1:n)

!     Determine the lower triangular part of the original matrix.

DO irow = n, 1, -1
  i = listp(irow)
  
!        Find the positions of the elements in the i-th row of the lower
!        triangular part of the original matrix that have already been
!        determined.
  
  DO ip = row_ptr(i), row_ptr(i+1)-1
    j = col_ind(ip)
    IF (j > 0) THEN
      DO jp = col_ptr(j), col_ptr(j+1)-1
        IF (row_ind(jp) == i) THEN
          iwk(j) = jp
          EXIT
        END IF
      END DO
    END IF
  END DO
  
!        Determine the elements in the i-th row of the lower
!        triangular part of the original matrix which get mapped
!        onto the lower triangular part of the permuted matrix.
  
  DO k = row_ptr(i), row_ptr(i+1)-1
    j = -col_ind(k)
    IF (j > 0) THEN
      col_ind(k) = j
      
!              Determine the (i,j) element.
      
      numg = ngrp(j)
      sum = zero
      DO ip = row_ptr(i), row_ptr(i+1)-1
        l = ABS(col_ind(ip))
        IF (ngrp(l) == numg .AND. l /= j) sum = sum + fhes(iwk(l))*eta(l)
      END DO
      DO jp = col_ptr(i), col_ptr(i+1)-1
        l = ABS(row_ind(jp))
        IF (ngrp(l) == numg .AND. l /= j) sum = sum + fhes(jp)*eta(l)
      END DO
      
!              Store the (i,j) element.
      
      DO jp = col_ptr(j), col_ptr(j+1)-1
        IF (row_ind(jp) == i) THEN
          fhes(jp) = (fhes(jp) - sum)/eta(j)
          EXIT
        END IF
      END DO
    END IF
  END DO
  
  
!        Determine the elements in the i-th row of the strict upper
!        triangular part of the original matrix which get mapped
!        onto the lower triangular part of the permuted matrix.
  
  DO k = col_ptr(i), col_ptr(i+1)-1
    j = -row_ind(k)
    IF (j > 0) THEN
      row_ind(k) = j
      
!              Determine the (i,j) element.
      
      numg = ngrp(j)
      sum = zero
      DO ip = row_ptr(i), row_ptr(i+1)-1
        l = ABS(col_ind(ip))
        IF (ngrp(l) == numg) sum = sum + fhes(iwk(l))*eta(l)
      END DO
      DO jp = col_ptr(i), col_ptr(i+1)-1
        l = ABS(row_ind(jp))
        IF (ngrp(l) == numg .AND. l /= j) sum = sum + fhes(jp)*eta(l)
      END DO
      
!              Store the (i,j) element.
      
      fhes(k) = (fhes(k) - sum)/eta(j)
    END IF
  END DO
END DO

!     Re-invert the array listp.

DO j = 1, n
  iwk(listp(j)) = j
END DO
listp(1:n) = iwk(1:n)

RETURN
END SUBROUTINE dsphesd



SUBROUTINE dsetsp(n, nnz, row_ind, col_ind, col_ptr, row_ptr,  &
                  method, info, listp, ngrp, maxgrp)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 11:18:32

INTEGER, INTENT(IN)      :: n
INTEGER, INTENT(OUT)     :: nnz
INTEGER, INTENT(OUT)     :: row_ind(:)
INTEGER, INTENT(IN OUT)  :: col_ind(:)
INTEGER, INTENT(IN OUT)  :: col_ptr(:)    ! col_ptr(n+1)
INTEGER, INTENT(IN OUT)  :: row_ptr(:)    ! row_ptr(n+1)
INTEGER, INTENT(IN)      :: method
INTEGER, INTENT(OUT)     :: info
INTEGER, INTENT(OUT)     :: listp(:)
INTEGER, INTENT(OUT)     :: ngrp(:)
INTEGER, INTENT(OUT)     :: maxgrp

!  **********

!  Subroutine dsetsp

!  Given the non-zero elements of a symmetric matrix A in coordinate
!  format, this subroutine computes the coloring information for
!  determining A from a lower triangular substitution method.

!  The sparsity pattern of the matrix A is specified by the
!  arrays row_ind and col_ind.  On input the indices for the
!  non-zero elements in the lower triangular part of A are

!        (row_ind(k),col_ind(k)), k = 1,2,...,nnz.

!  The (row_ind(k),col_ind(k)) pairs may be specified in any order.
!  Duplicate input pairs are permitted, but they are eliminated.
!  Diagonal elements must be part of the sparsity pattern.
!  Any pair (row_ind(k),col_ind(k)), where row_ind(k) is less than
!  col_ind(k), is replaced by the pair (col_ind(k),row_ind(k)).

!  The input coordinate format is changed to a storage format.
!  On output the strict lower triangular part of A is stored in both
!  compressed column storage and compressed row storage.

!  The information required to determine the matrix from a lower
!  triangular substitution method is obtained from subroutine dssm.

!  The subroutine statement is

!    subroutine dsetsp(n, nnz, row_ind, col_ind, col_ptr, row_ptr,
!                      method, info, listp, ngrp, maxgrp, iwa)

!  where

!    n is an integer variable.
!      On entry n is the order of the matrix.
!      On exit n is unchanged.

!    nnz is an integer variable.
!      On entry nnz is the number of non-zeros entries in the
!        coordinate format.
!      On exit nnz is the number of non-zeroes in the strict
!        lower triangular part of the matrix A.

!    row_ind is an integer array of length nnz.
!      On entry row_ind must contain the row indices of the non-zero
!         elements of A in coordinate format.
!      On exit row_ind contains row indices for the strict
!         lower triangular part of A in compressed column storage.

!    col_ind is an integer array of length nnz.
!      On entry col_ind must contain the column indices of the
!         non-zero elements of A in coordinate format.
!      On exit col_ind contains column indices for the strict
!         lower triangular part of A in compressed row storage.

!    row_ptr is an integer array of length n + 1.
!      On entry row_ptr need not be specified.
!      On exit row_ptr must contain pointers to the rows of A.
!         The non-zeros in row j of A must be in positions
!         row_ptr(j), ... , row_ptr(j+1) - 1 of col_ind.

!    col_ptr is an integer array of length n + 1.
!      On entry col_ptr need not be specified.
!      On exit col_ptr must contain pointers to the columns of A.
!         The non-zeros in column j of A must be in positions
!         col_ptr(j), ... , col_ptr(j+1) - 1 of row_ind.

!    method is an integer variable.
!      On input with method = 1, the direct method is used to
!        determine the partition and symmetric permutation.
!        Otherwise, the indirect method is used.

!    info is an integer variable.
!      On input info need not be specified.
!      On output info is set as follows.

!         info = 1  Normal termination.

!         info = 0  Input n or nnz are not positive.

!         info < 0  There is an error in the sparsity pattern.
!                   If k = -info then row_ind(k) or col_ind(k) is
!                   not an integer between 1 and n, or the k-th
!                   diagonal element is not in the sparsity pattern.

!    listp is an integer array of length n.
!      On input listp need not be specified.
!      On output listp specifies a permutation of the matrix.
!         Element (i,j) of the matrix is the (listp(i),listp(j))
!         element of the permuted matrix.

!    ngrp is an integer array of length n.
!      On entry ngrp need not be specified.
!      On exit ngrp specifies the partition of the columns
!         of A. Column j belongs to group ngrp(j).

!    maxgrp is an integer variable.
!      On entry maxgrp need not be specified.
!      On exit maxgrp specifies the number of groups in
!         the partition of the columns of A.

!  Subprograms called

!    MINPACK-2  ......  dssm

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Chih-Jen Lin and Jorge J. More'.

!  **********

INTEGER :: i, j, mingrp, ndiag

!     Subroutine dssm first checks the sparsity structure.
!     If there are no errors, the input format is changed to
!     compressed column storage. Finally, the information required
!     to determine the matrix from matrix-vector products is obtained.

CALL dssm(n, nnz, row_ind, col_ind, method, listp,  &
          ngrp, maxgrp, mingrp, info, row_ptr, col_ptr)

!     Exit if there are errors on input.

IF (info <= 0) RETURN

!     Change the sparsity structure to exclude the diagonal entries.

ndiag = 0
DO j = 1, n
  DO i = col_ptr(j), col_ptr(j+1) - 1
    IF (row_ind(i) == j) THEN
      ndiag = ndiag + 1
    ELSE
      row_ind(i-ndiag) =  row_ind(i)
    END IF
  END DO
  col_ptr(j) = col_ptr(j) - (j - 1)
END DO
col_ptr(n+1) = col_ptr(n+1) - n
nnz = nnz - n

RETURN
END SUBROUTINE dsetsp



SUBROUTINE deptfg(nx, ny, x, f, fgrad, task, c)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:07

INTEGER, INTENT(IN)            :: nx
INTEGER, INTENT(IN)            :: ny
REAL (dp), INTENT(IN OUT)      :: x(:)    ! x(nx*ny)
REAL (dp), INTENT(OUT)         :: f
REAL (dp), INTENT(OUT)         :: fgrad(:)    ! fgrad(nx*ny)
CHARACTER (LEN=*), INTENT(IN)  :: task
REAL (dp), INTENT(IN)          :: c

!  **********

!  Subroutine deptfg

!  This subroutine computes the function and gradient of the
!  elastic-plastic torsion problem.

!  The subroutine statement is

!    subroutine deptfg(nx, ny, x, f, fgrad, task, c)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first coordinate
!         direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'.  Otherwise
!         x is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    fgrad is a REAL (dp) array of dimension nx*ny.
!      On entry fgrad need not be specified.
!      On exit fgrad contains the gradient evaluated at x if
!         task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at x.
!          'G'     Evaluate the gradient at x.
!          'FG'    Evaluate the function and the gradient at x.
!          'XS'    Set x to the standard starting point xs.
!          'XL'    Set x to the lower bound xl.
!          'XU'    Set x to the upper bound xu.

!      On exit task is unchanged.

!    c is a REAL (dp) variable.
!      On entry c is the angle of twist per unit length.
!      On exit c is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick and Jorge J. More'.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp, three=3.0_dp

LOGICAL   :: feval, geval
INTEGER   :: i, j, k
REAL (dp) :: area, cdiv3, dvdx, dvdy, flin, fquad, hx, hy,  &
             temp, temp1, v, vb, vl, vr, vt

hx = one/DBLE(nx+1)
hy = one/DBLE(ny+1)
area = p5*hx*hy
cdiv3 = c/three

!     Compute a lower bound for x if task = 'XL' or
!     an upper bound if task = 'XU'.

IF (task(1:2) == 'XL' .OR. task(1:2) == 'XU') THEN
  IF (task(1:2) == 'XL') temp1 = -one
  IF (task(1:2) == 'XU') temp1 = one
  DO  j = 1, ny
    temp = DBLE(MIN(j,ny-j+1))*hy
    DO  i = 1, nx
      k = nx*(j-1) + i
      x(k) = SIGN(MIN(DBLE(MIN(i,nx-i+1))*hx,temp),temp1)
    END DO
  END DO
  
  RETURN
  
END IF

!     Compute the standard starting point if task = 'XS'.

IF (task(1:2) == 'XS') THEN
  DO  j = 1, ny
    temp = DBLE(MIN(j,ny-j+1))*hy
    DO  i = 1, nx
      k = nx*(j-1) + i
      x(k) = MIN(DBLE(MIN(i,nx-i+1))*hx,temp)
    END DO
  END DO
  
  RETURN
  
END IF

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  feval = .true.
ELSE
  feval = .false.
END IF
IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') THEN
  geval = .true.
ELSE
  geval = .false.
END IF

!     Evaluate the function if task = 'F', the gradient if task = 'G',
!     or both if task = 'FG'.

IF (feval) THEN
  fquad = zero
  flin = zero
END IF
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = zero
  END DO
END IF

!     Computation of the function and the gradient over the lower
!     triangular elements.

DO  j = 0, ny
  DO  i = 0, nx
    k = nx*(j-1) + i
    v = zero
    vr = zero
    vt = zero
    IF (i >= 1 .AND. j >= 1) v = x(k)
    IF (i < nx .AND. j > 0) vr = x(k+1)
    IF (i > 0 .AND. j < ny) vt = x(k+nx)
    dvdx = (vr-v)/hx
    dvdy = (vt-v)/hy
    IF (feval) THEN
      fquad = fquad + dvdx**2 + dvdy**2
      flin = flin - cdiv3*(v+vr+vt)
    END IF
    IF (geval) THEN
      IF (i /= 0 .AND. j /= 0)  &
          fgrad(k) = fgrad(k) - dvdx/hx - dvdy/hy - cdiv3
      IF (i /= nx .AND. j /= 0) fgrad(k+1) = fgrad(k+1) + dvdx/hx - cdiv3
      IF (i /= 0 .AND. j /= ny) fgrad(k+nx) = fgrad(k+nx) + dvdy/hy - cdiv3
    END IF
  END DO
END DO

!     Computation of the function and the gradient over the upper
!     triangular elements.

DO  j = 1, ny + 1
  DO  i = 1, nx + 1
    k = nx*(j-1) + i
    vb = zero
    vl = zero
    v = zero
    IF (i <= nx .AND. j > 1) vb = x(k-nx)
    IF (i > 1 .AND. j <= ny) vl = x(k-1)
    IF (i <= nx .AND. j <= ny) v = x(k)
    dvdx = (v-vl)/hx
    dvdy = (v-vb)/hy
    IF (feval) THEN
      fquad = fquad + dvdx**2 + dvdy**2
      flin = flin - cdiv3*(vb+vl+v)
    END IF
    IF (geval) THEN
      IF (i /= nx+1 .AND. j /= 1) fgrad(k-nx) = fgrad(k-nx) - dvdy/hy - cdiv3
      IF (i /= 1 .AND. j /= ny+1) fgrad(k-1) = fgrad(k-1) - dvdx/hx - cdiv3
      IF (i /= nx+1 .AND. j /= ny+1)  &
          fgrad(k) = fgrad(k) + dvdx/hx + dvdy/hy - cdiv3
    END IF
  END DO
END DO

!     Scale the result.

IF (feval) f = area*(p5*fquad+flin)
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = area*fgrad(k)
  END DO
END IF

RETURN
END SUBROUTINE deptfg



SUBROUTINE dpjbfg(nx, ny, x, f, fgrad, task, ecc, b)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:08

INTEGER, INTENT(IN)            :: nx
INTEGER, INTENT(IN)            :: ny
REAL (dp), INTENT(OUT)         :: x(:)    ! x(nx*ny)
REAL (dp), INTENT(OUT)         :: f
REAL (dp), INTENT(OUT)         :: fgrad(:)    ! fgrad(nx*ny)
CHARACTER (LEN=*), INTENT(IN)  :: task
REAL (dp), INTENT(IN)          :: ecc
REAL (dp), INTENT(IN)          :: b

!  **********

!  Subroutine dpjbfg

!  This subroutine computes the function and gradient of the
!  pressure distribution in a journal bearing problem.

!  The subroutine statement is

!    subroutine dpjbfg(nx, ny, x, f, fgrad, task, ecc, b)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
!         x is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    fgrad is a REAL (dp) array of dimension nx*ny.
!      On entry fgrad need not be specified.
!      On exit fgrad contains the gradient evaluated at x if task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at x.
!          'G'     Evaluate the gradient at x.
!          'FG'    Evaluate the function and the gradient at x.
!          'XS'    Set x to the standard starting point xs.
!          'XL'    Set x to the lower bound xl.

!      On exit task is unchanged.

!    ecc is a REAL (dp) variable
!      On entry ecc is the eccentricity in (0,1).
!      On exit ecc is unchanged

!    b is a REAL (dp) variable
!      On entry b defines the domain as D = (0,2*pi) X (0,2*b).
!      On exit b is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick and Jorge J. More'.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp
REAL (dp), PARAMETER :: two=2.0_dp
REAL (dp), PARAMETER :: four=4.0_dp
REAL (dp), PARAMETER :: six=6.0_dp

LOGICAL   :: feval, geval
INTEGER   :: i, j, k
REAL (dp) :: dvdx, dvdy, ehxhy, flin, fquad, hx, hxhy, hy, pi,  &
             temp, trule, v, vb, vl, vr, vt, xi

! REAL (dp) :: p

! p(xi) = (1 + ecc*COS(xi))**3

!     Initialization.

pi = four*ATAN(one)
hx = two*pi/DBLE(nx+1)
hy = two*b/DBLE(ny+1)
hxhy = hx*hy
ehxhy = ecc*hxhy

!     Compute the lower bound xl for x if task = 'XL'.

IF (task(1:2) == 'XL') THEN
  DO  k = 1, nx*ny
    x(k) = zero
  END DO
  
  RETURN
  
END IF

!     Compute the standard starting point if task = 'XS'.

IF (task(1:2) == 'XS') THEN
  DO  i = 1, nx
    temp = MAX(SIN(DBLE(i)*hx), zero)
    DO  j = 1, ny
      k = nx*(j-1) + i
      x(k) = temp
    END DO
  END DO
  
  RETURN
  
END IF

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  feval = .true.
ELSE
  feval = .false.
END IF
IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') THEN
  geval = .true.
ELSE
  geval = .false.
END IF

!     Compute the function if task = 'F', the gradient if task = 'G', or
!     both if task = 'FG'.

IF (feval) THEN
  fquad = zero
  flin = zero
END IF
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = zero
  END DO
END IF

!     Computation of the quadratic part of the function and corresponding
!     components of the gradient over the lower triangular elements.

DO  i = 0, nx
  xi = DBLE(i)*hx
  trule = hxhy*(p(xi, ecc) + p(xi+hx, ecc) + p(xi, ecc))/six
  DO  j = 0, ny
    k = nx*(j-1) + i
    v = zero
    vr = zero
    vt = zero
    IF (i /= 0 .AND. j /= 0) v = x(k)
    IF (i /= nx .AND. j /= 0) vr = x(k+1)
    IF (i /= 0 .AND. j /= ny) vt = x(k+nx)
    dvdx = (vr-v)/hx
    dvdy = (vt-v)/hy
    IF (feval) fquad = fquad + trule*(dvdx**2+dvdy**2)
    IF (geval) THEN
      IF (i /= 0 .AND. j /= 0) fgrad(k) = fgrad(k) - trule*(dvdx/hx+dvdy/hy)
      IF (i /= nx .AND. j /= 0) fgrad(k+1) = fgrad(k+1) + trule*dvdx/hx
      IF (i /= 0 .AND. j /= ny) fgrad(k+nx) = fgrad(k+nx) + trule*dvdy/hy
    END IF
  END DO
END DO

!     Computation of the quadratic part of the function and corresponding
!     components of the gradient over the upper triangular elements.

DO  i = 1, nx + 1
  xi = DBLE(i)*hx
  trule = hxhy*(p(xi, ecc) + p(xi-hx, ecc) + p(xi, ecc))/six
  DO  j = 1, ny + 1
    k = nx*(j-1) + i
    vb = zero
    vl = zero
    v = zero
    IF (i /= nx+1 .AND. j /= 1) vb = x(k-nx)
    IF (i /= 1 .AND. j /= ny+1) vl = x(k-1)
    IF (i /= nx+1 .AND. j /= ny+1) v = x(k)
    dvdx = (v-vl)/hx
    dvdy = (v-vb)/hy
    IF (feval) fquad = fquad + trule*(dvdx**2+dvdy**2)
    IF (geval) THEN
      IF (i <= nx .AND. j > 1) fgrad(k-nx) = fgrad(k-nx) - trule*dvdy/hy
      IF (i > 1 .AND. j <= ny) fgrad(k-1) = fgrad(k-1) - trule*dvdx/hx
      IF (i <= nx .AND. j <= ny) fgrad(k) = fgrad(k) + trule*(dvdx/hx+dvdy/hy)
    END IF
  END DO
END DO

!     Computation of the linear part of the function and
!     corresponding components of the gradient.

DO  i = 1, nx
  temp = SIN(DBLE(i)*hx)
  IF (feval) THEN
    DO  j = 1, ny
      k = nx*(j-1) + i
      flin = flin + temp*x(k)
    END DO
  END IF
  IF (geval) THEN
    DO  j = 1, ny
      k = nx*(j-1) + i
      fgrad(k) = fgrad(k) - ehxhy*temp
    END DO
  END IF
END DO

!     Finish off the function.

IF (feval) f = p5*fquad - ehxhy*flin

RETURN
END SUBROUTINE dpjbfg




FUNCTION p(xi, ecc) RESULT(fn_val)
! This is the statement function which was part of routine DPJBFG.

REAL (dp), INTENT(IN) :: xi, ecc
REAL (dp)             :: fn_val

fn_val = (1 + ecc*COS(xi))**3
RETURN
END FUNCTION p




SUBROUTINE dmsabc(nx, ny, bottom, top, left, right)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

INTEGER, INTENT(IN)     :: nx
INTEGER, INTENT(IN)     :: ny
REAL (dp), INTENT(OUT)  :: bottom(:)    ! bottom(nx+2)
REAL (dp), INTENT(OUT)  :: top(:)       ! top(nx+2)
REAL (dp), INTENT(OUT)  :: left(:)      ! left(ny+2)
REAL (dp), INTENT(OUT)  :: right(:)     ! right(ny+2)

!  **********

!  Subroutine dmsabc

!  This subroutine computes Enneper's boundary conditions for the minimal
!  surface area problem on the unit square centered at the origin.

!  The subroutine statement is

!    subroutine dmsabc(nx, ny, bottom, top, left, right)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first coordinate
!         direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    bottom is a REAL (dp) array of dimension nx + 2.
!      On entry bottom need not be specified.
!      On exit bottom contains boundary values for the bottom
!         boundary of the domain.

!    top is a REAL (dp) array of dimension nx + 2.
!      On entry top need not be specified.
!      On exit top contains boundary values for the top boundary of the domain.
!    left is a REAL (dp) array of dimension ny + 2.
!      On entry left need not be specified.
!      On exit left contains boundary values for the left boundary
!         of the domain.

!    right is a REAL (dp) array of dimension ny + 2.
!      On entry right need not be specified.
!      On exit right contains boundary values for the right boundary
!         of the domain.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: two=2.0_dp, three=3.0_dp, tol=1.0D-10, b=-.50_dp, &
                        t=.50_dp, l=-.50_dp, r=.50_dp
INTEGER, PARAMETER   :: maxit=5

INTEGER   :: i, j, k, limit
REAL (dp) :: det, fnorm, hx, hy, xt, yt
REAL (dp) :: nf(2), njac(2,2), u(2)

!     Compute Enneper's boundary conditions: bottom, top, left, then
!     right.  Enneper's boundary values are obtained by defining
!     bv(x,y) = u**2 - v**2 where u and v are the unique solutions of
!     x = u + u*(v**2) - (u**3)/3, y = -v - (u**2)*v + (v**3)/3.

hx = (r-l)/DBLE(nx+1)
hy = (t-b)/DBLE(ny+1)

DO  j = 1, 4
  IF (j == 1) THEN
    yt = b
    xt = l
    limit = nx + 2
  ELSE IF (j == 2) THEN
    yt = t
    xt = l
    limit = nx + 2
  ELSE IF (j == 3) THEN
    yt = b
    xt = l
    limit = ny + 2
  ELSE IF (j == 4) THEN
    yt = b
    xt = r
    limit = ny + 2
  END IF
  
!        Use Newton's method to solve xt = u + u*(v**2) - (u**3)/3,
!        yt = -v - (u**2)*v + (v**3)/3.
  
  DO  i = 1, limit
    u(1) = xt
    u(2) = -yt
    DO  k = 1, maxit
      nf(1) = u(1) + u(1)*u(2)**2 - u(1)**3/three - xt
      nf(2) = -u(2) - u(1)**2*u(2) + u(2)**3/three - yt
      fnorm = SQRT(nf(1)*nf(1)+nf(2)*nf(2))
      IF (fnorm <= tol) EXIT
      njac(1,1) = one + u(2)**2 - u(1)**2
      njac(1,2) = two*u(1)*u(2)
      njac(2,1) = -two*u(1)*u(2)
      njac(2,2) = -one - u(1)**2 + u(2)**2
      det = njac(1,1)*njac(2,2) - njac(1,2)*njac(2,1)
      u(1) = u(1) - (njac(2,2)*nf(1)-njac(1,2)*nf(2))/det
      u(2) = u(2) - (njac(1,1)*nf(2)-njac(2,1)*nf(1))/det
    END DO
    
    IF (j == 1) THEN
      bottom(i) = u(1)*u(1) - u(2)*u(2)
      xt = xt + hx
    ELSE IF (j == 2) THEN
      top(i) = u(1)*u(1) - u(2)*u(2)
      xt = xt + hx
    ELSE IF (j == 3) THEN
      left(i) = u(1)*u(1) - u(2)*u(2)
      yt = yt + hy
    ELSE IF (j == 4) THEN
      right(i) = u(1)*u(1) - u(2)*u(2)
      yt = yt + hy
    END IF
  END DO
END DO

RETURN
END SUBROUTINE dmsabc



SUBROUTINE dmsafg(nx, ny, x, f, fgrad, task, bottom, top, left, right)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:08

INTEGER, INTENT(IN)            :: nx
INTEGER, INTENT(IN)            :: ny
REAL (dp), INTENT(OUT)         :: x(:)        ! x(nx*ny)
REAL (dp), INTENT(OUT)         :: f
REAL (dp), INTENT(OUT)         :: fgrad(:)    ! fgrad(nx*ny)
CHARACTER (LEN=*), INTENT(IN)  :: task
REAL (dp), INTENT(IN)          :: bottom(:)   ! bottom(nx+2)
REAL (dp), INTENT(IN)          :: top(:)      ! top(nx+2)
REAL (dp), INTENT(IN)          :: left(:)     ! left(ny+2)
REAL (dp), INTENT(IN)          :: right(:)    ! right(ny+2)

!  **********

!  Subroutine dmsafg

!  This subroutine computes the function and gradient of the
!  minimal surface area problem.

!  The subroutine statement is

!    subroutine dmsafg(nx, ny, x, f, fgrad, task, bottom, top, left, right)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
!         x is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    fgrad is a REAL (dp) array of dimension nx*ny.
!      On entry fgrad need not be specified.
!      On exit fgrad contains the gradient evaluated at x if
!         task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at x.
!          'G'     Evaluate the gradient at x.
!          'FG'    Evaluate the function and the gradient at x.
!          'XS'    Set x to the standard starting point xs.

!      On exit task is unchanged.

!    bottom is a REAL (dp) array of dimension nx + 2.
!      On entry bottom must contain boundary data beginning
!         with the lower left corner of the domain.
!      On exit bottom is unchanged.

!    top is a REAL (dp) array of dimension nx + 2.
!      On entry top must contain boundary data beginning with
!         the upper left corner of the domain.
!      On exit top is unchanged.

!    left is a REAL (dp) array of dimension ny + 2.
!      On entry left must contain boundary data beginning with
!         the lower left corner of the domain.
!      On exit left is unchanged.

!    right is a REAL (dp) array of dimension ny + 2.
!      On entry right must contain boundary data beginning with
!         the lower right corner of the domain.
!      On exit right is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp, two=2.0_dp

LOGICAL   :: feval, geval
INTEGER   :: i, j, k
REAL (dp) :: alphaj, area, betai, dvdx, dvdy, fl, fu, hx, hy,  &
             v, vb, vl, vr, vt, xline, yline

!     Initialize.

hx = one/DBLE(nx+1)
hy = one/DBLE(ny+1)
area = p5*hx*hy

!     Compute the standard starting point if task = 'XS'.

IF (task(1:2) == 'XS') THEN
  DO  j = 1, ny
    alphaj = DBLE(j)*hy
    DO  i = 1, nx
      k = nx*(j-1) + i
      betai = DBLE(i)*hx
      yline = alphaj*top(i+1) + (one-alphaj)*bottom(i+1)
      xline = betai*right(j+1) + (one-betai)*left(j+1)
      x(k) = (yline+xline)/two
    END DO
  END DO
  
  RETURN
  
END IF

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  feval = .true.
ELSE
  feval = .false.
END IF
IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') THEN
  geval = .true.
ELSE
  geval = .false.
END IF

!     Evaluate the function if task = 'F', the gradient if task = 'G',
!     or both if task = 'FG'.

IF (feval) f = zero
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = zero
  END DO
END IF

!     Computation of the function and gradient over the lower
!     triangular elements.

DO  j = 0, ny
  DO  i = 0, nx
    k = nx*(j-1) + i
    IF (i >= 1 .AND. j >= 1) THEN
      v = x(k)
    ELSE
      IF (j == 0) v = bottom(i+1)
      IF (i == 0) v = left(j+1)
    END IF
    IF (i < nx .AND. j > 0) THEN
      vr = x(k+1)
    ELSE
      IF (i == nx) vr = right(j+1)
      IF (j == 0) vr = bottom(i+2)
    END IF
    IF (i > 0 .AND. j < ny) THEN
      vt = x(k+nx)
    ELSE
      IF (i == 0) vt = left(j+2)
      IF (j == ny) vt = top(i+1)
    END IF
    dvdx = (vr-v)/hx
    dvdy = (vt-v)/hy
    fl = SQRT(one+dvdx**2+dvdy**2)
    IF (feval) f = f + fl
    IF (geval) THEN
      IF (i >= 1 .AND. j >= 1) fgrad(k) = fgrad(k) - (dvdx/hx+dvdy/hy)/fl
      IF (i < nx .AND. j > 0) fgrad(k+1) = fgrad(k+1) + (dvdx/hx)/fl
      IF (i > 0 .AND. j < ny) fgrad(k+nx) = fgrad(k+nx) + (dvdy/hy)/fl
    END IF
  END DO
END DO

!     Computation of the function and the gradient over the upper
!     triangular elements.

DO  j = 1, ny + 1
  DO  i = 1, nx + 1
    k = nx*(j-1) + i
    IF (i <= nx .AND. j > 1) THEN
      vb = x(k-nx)
    ELSE
      IF (j == 1) vb = bottom(i+1)
      IF (i == nx+1) vb = right(j)
    END IF
    IF (i > 1 .AND. j <= ny) THEN
      vl = x(k-1)
    ELSE
      IF (j == ny+1) vl = top(i)
      IF (i == 1) vl = left(j+1)
    END IF
    IF (i <= nx .AND. j <= ny) THEN
      v = x(k)
    ELSE
      IF (i == nx+1) v = right(j+1)
      IF (j == ny+1) v = top(i+1)
    END IF
    dvdx = (v-vl)/hx
    dvdy = (v-vb)/hy
    fu = SQRT(one+dvdx**2+dvdy**2)
    IF (feval) f = f + fu
    IF (geval) THEN
      IF (i <= nx .AND. j > 1) fgrad(k-nx) = fgrad(k-nx) - (dvdy/hy)/fu
      IF (i > 1 .AND. j <= ny) fgrad(k-1) = fgrad(k-1) - (dvdx/hx)/fu
      IF (i <= nx .AND. j <= ny) fgrad(k) = fgrad(k) + (dvdx/hx+dvdy/hy)/fu
    END IF
  END DO
END DO

!     Scale the function and the gradient.

IF (feval) f = area*f
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = area*fgrad(k)
  END DO
END IF

RETURN
END SUBROUTINE dmsafg



SUBROUTINE dodcfg(nx, ny, x, f, fgrad, task, lambda)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:08

INTEGER, INTENT(IN)            :: nx
INTEGER, INTENT(IN)            :: ny
REAL (dp), INTENT(OUT)         :: x(:)    ! x(nx*ny)
REAL (dp), INTENT(OUT)         :: f
REAL (dp), INTENT(OUT)         :: fgrad(:)    ! fgrad(nx*ny)
CHARACTER (LEN=*), INTENT(IN)  :: task
REAL (dp), INTENT(IN)          :: lambda

!  **********

!  Subroutine dodcfg

!  This subroutine computes the function and gradient of the
!  optimal design with composite materials problem.

!  The subroutine statement is

!    subroutine dodcfg(nx, ny, x, f, fgrad, task, lambda)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
!         x is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    fgrad is a REAL (dp) array of dimension nx*ny.
!      On entry fgrad need not be specified.
!      On exit fgrad contains the gradient evaluated at x if
!         task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at x.
!          'G'     Evaluate the gradient at x.
!          'FG'    Evaluate the function and the gradient at x.
!          'XS'    Set x to the standard starting point xs.

!      On exit task is unchanged.

!    lambda is a REAL (dp) variable.
!      On entry lambda is the Lagrange multiplier.
!      On exit lambda is unchanged.

!  Subprograms called

!    MINPACK-supplied   ...   dodcps

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp, two=2.0_dp, mu1=one, mu2=two

LOGICAL   :: feval, geval
INTEGER   :: i, j, k
REAL (dp) :: area, dpsi, dpsip, dvdx, dvdy, gradv, hx, hxhy,  &
             hy, temp, t1, t2, v, vb, vl, vr, vt

! EXTERNAL dodcps

!     Initialization.

hx = one/DBLE(nx+1)
hy = one/DBLE(ny+1)
hxhy = hx*hy
area = p5*hxhy

!     Compute the break points.

t1 = SQRT(two*lambda*mu1/mu2)
t2 = SQRT(two*lambda*mu2/mu1)

!     Compute the standard starting point if task = 'XS'.

IF (task(1:2) == 'XS') THEN
  DO  j = 1, ny
    temp = DBLE(MIN(j,ny-j+1))*hy
    DO  i = 1, nx
      k = nx*(j-1) + i
      x(k) = -(MIN(DBLE(MIN(i,nx-i+1))*hx,temp))**2
    END DO
  END DO
  
  RETURN
  
END IF

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  feval = .true.
ELSE
  feval = .false.
END IF
IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') THEN
  geval = .true.
ELSE
  geval = .false.
END IF

!     Evaluate the function if task = 'F', the gradient if task = 'G',
!     or both if task = 'FG'.

IF (feval) f = zero
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = zero
  END DO
END IF

!     Computation of the function and the gradient over the lower
!     triangular elements.

DO  j = 0, ny
  DO  i = 0, nx
    k = nx*(j-1) + i
    v = zero
    vr = zero
    vt = zero
    IF (j >= 1 .AND. i >= 1) v = x(k)
    IF (i < nx .AND. j > 0) vr = x(k+1)
    IF (i > 0 .AND. j < ny) vt = x(k+nx)
    dvdx = (vr-v)/hx
    dvdy = (vt-v)/hy
    gradv = dvdx**2 + dvdy**2
    IF (feval) THEN
      CALL dodcps(gradv, mu1, mu2, t1, t2, dpsi, 0, lambda)
      f = f + dpsi
    END IF
    IF (geval) THEN
      CALL dodcps(gradv, mu1, mu2, t1, t2, dpsip, 1, lambda)
      IF (i >= 1 .AND. j >= 1)  &
          fgrad(k) = fgrad(k) - two*(dvdx/hx+dvdy/hy)*dpsip
      IF (i < nx .AND. j > 0) fgrad(k+1) = fgrad(k+1) + two*(dvdx/hx)*dpsip
      IF (i > 0 .AND. j < ny) fgrad(k+nx) = fgrad(k+nx) + two*(dvdy/hy)*dpsip
    END IF
  END DO
END DO

!     Computation of the function and the gradient over the upper
!     triangular elements.

DO  j = 1, ny + 1
  DO  i = 1, nx + 1
    k = nx*(j-1) + i
    vb = zero
    vl = zero
    v = zero
    IF (i <= nx .AND. j > 1) vb = x(k-nx)
    IF (i > 1 .AND. j <= ny) vl = x(k-1)
    IF (i <= nx .AND. j <= ny) v = x(k)
    dvdx = (v-vl)/hx
    dvdy = (v-vb)/hy
    gradv = dvdx**2 + dvdy**2
    IF (feval) THEN
      CALL dodcps(gradv,mu1,mu2,t1,t2,dpsi,0,lambda)
      f = f + dpsi
    END IF
    IF (geval) THEN
      CALL dodcps(gradv,mu1,mu2,t1,t2,dpsip,1,lambda)
      IF (i <= nx .AND. j > 1)  &
          fgrad(k-nx) = fgrad(k-nx) - two*(dvdy/hy)*dpsip
      IF (i > 1 .AND. j <= ny) fgrad(k-1) = fgrad(k-1) - two*(dvdx/hx)*dpsip
      IF (i <= nx .AND. j <= ny)  &
          fgrad(k) = fgrad(k) + two*(dvdx/hx+dvdy/hy)*dpsip
    END IF
  END DO
END DO

!     Scale the function.

IF (feval) f = area*f

!     Integrate v over the domain.

IF (feval) THEN
  temp = zero
  DO  k = 1, nx*ny
    temp = temp + x(k)
  END DO
  f = f + hxhy*temp
END IF
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = area*fgrad(k) + hxhy
  END DO
END IF

RETURN
END SUBROUTINE dodcfg



SUBROUTINE dsscfg(nx, ny, x, f, fgrad, task, lambda)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)            :: nx
INTEGER, INTENT(IN)            :: ny
REAL (dp), INTENT(OUT)         :: x(:)    ! x(nx*ny)
REAL (dp), INTENT(OUT)         :: f
REAL (dp), INTENT(OUT)         :: fgrad(:)    ! fgrad(nx*ny)
CHARACTER (LEN=*), INTENT(IN)  :: task
REAL (dp), INTENT(IN)          :: lambda

!  **********

!  Subroutine dsscfg

!  This subroutine computes the function and gradient of the
!  steady state combustion problem.

!  The subroutine statement is

!    subroutine dsscfg(nx, ny, x, f, fgrad, task, lambda)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'.  Otherwise
!         x is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    fgrad is a REAL (dp) array of dimension nx*ny.
!      On entry fgrad need not be specified.
!      On exit fgrad contains the gradient evaluated at x if
!         task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!         'F'      Evaluate the function at x.
!         'G'      Evaluate the gradient at x.
!         'FG'     Evaluate the function and the gradient at x.
!         'XS'     Set x to the standard starting point xs.
!         'XL'     Set x to the lower bound xl.
!         'XU'     Set x to the upper bound xu.

!      On exit task is unchanged.

!    lambda is a REAL (dp) variable.
!      On entry lambda is a nonnegative Frank-Kamenetski parameter.
!      On exit lambda is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp, three=3.0_dp

LOGICAL   :: feval, geval
INTEGER   :: i, j, k
REAL (dp) :: area, dvdx, dvdy, expv, expvb, expvl, expvr,  &
             expvt, fexp, fquad, hx, hy, v, vb, vl, vr, vt, temp, temp1

!     Initialization.

hx = one/DBLE(nx+1)
hy = one/DBLE(ny+1)
area = p5*hx*hy

!     Compute the standard starting point if task = 'XS'.

IF (task(1:2) == 'XS') THEN
  
  temp1 = lambda/(lambda+one)
  DO  j = 1, ny
    temp = DBLE(MIN(j,ny-j+1))*hy
    DO  i = 1, nx
      k = nx*(j-1) + i
      x(k) = temp1*SQRT(MIN(DBLE(MIN(i,nx-i+1))*hx,temp))
    END DO
  END DO
  
  RETURN
  
END IF

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  feval = .true.
ELSE
  feval = .false.
END IF
IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') THEN
  geval = .true.
ELSE
  geval = .false.
END IF

!     Compute the function if task = 'F', the gradient if task = 'G', or
!     both if task = 'FG'.

IF (feval) THEN
  fquad = zero
  fexp = zero
END IF
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = zero
  END DO
END IF

!     Computation of the function and the gradient over the lower
!     triangular elements.  The trapezoidal rule is used to estimate
!     the integral of the exponential term.

DO  j = 0, ny
  DO  i = 0, nx
    k = nx*(j-1) + i
    v = zero
    vr = zero
    vt = zero
    IF (i /= 0 .AND. j /= 0) v = x(k)
    IF (i /= nx .AND. j /= 0) vr = x(k+1)
    IF (i /= 0 .AND. j /= ny) vt = x(k+nx)
    dvdx = (vr-v)/hx
    dvdy = (vt-v)/hy
    expv = EXP(v)
    expvr = EXP(vr)
    expvt = EXP(vt)
    IF (feval) THEN
      fquad = fquad + dvdx**2 + dvdy**2
      fexp = fexp - lambda*(expv+expvr+expvt)/three
    END IF
    IF (geval) THEN
      IF (i /= 0 .AND. j /= 0) fgrad(k) = fgrad(k) -  &
          dvdx/hx - dvdy/hy - lambda*expv/three
      IF (i /= nx .AND. j /= 0) fgrad(k+1) = fgrad(k+1) +  &
          dvdx/hx - lambda*expvr/three
      IF (i /= 0 .AND. j /= ny) fgrad(k+nx) = fgrad(k+nx) +  &
          dvdy/hy - lambda*expvt/three
    END IF
  END DO
END DO

!     Computation of the function and the gradient over the upper
!     triangular elements.  The trapezoidal rule is used to estimate
!     the integral of the exponential term.

DO  j = 1, ny + 1
  DO  i = 1, nx + 1
    k = nx*(j-1) + i
    vb = zero
    vl = zero
    v = zero
    IF (i /= nx+1 .AND. j /= 1) vb = x(k-nx)
    IF (i /= 1 .AND. j /= ny+1) vl = x(k-1)
    IF (i /= nx+1 .AND. j /= ny+1) v = x(k)
    dvdx = (v-vl)/hx
    dvdy = (v-vb)/hy
    expvb = EXP(vb)
    expvl = EXP(vl)
    expv = EXP(v)
    IF (feval) THEN
      fquad = fquad + dvdx**2 + dvdy**2
      fexp = fexp - lambda*(expvb+expvl+expv)/three
    END IF
    IF (geval) THEN
      IF (i /= nx+1 .AND. j /= 1) fgrad(k-nx) = fgrad(k-nx) - dvdy/hy -  &
          lambda*expvb/three
      IF (i /= 1 .AND. j /= ny+1) fgrad(k-1) = fgrad(k-1) -  &
          dvdx/hx - lambda*expvl/three
      IF (i /= nx+1 .AND. j /= ny+1) fgrad(k) = fgrad(k) +  &
          dvdx/hx + dvdy/hy - lambda*expv/three
    END IF
  END DO
END DO

!     Scale the result.

IF (feval) f = area*(p5*fquad+fexp)
IF (geval) THEN
  DO  k = 1, nx*ny
    fgrad(k) = area*fgrad(k)
  END DO
END IF

RETURN
END SUBROUTINE dsscfg



SUBROUTINE dgl1fg(n, x, f, fgrad, task, t)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:08

INTEGER, INTENT(IN)            :: n
REAL (dp), INTENT(IN OUT)      :: x(:)
REAL (dp), INTENT(OUT)         :: f
REAL (dp), INTENT(OUT)         :: fgrad(:)
CHARACTER (LEN=*), INTENT(IN)  :: task
REAL (dp), INTENT(IN)          :: t

!  **********

!  Subroutine dgl1fg

!  This subroutine computes the function and gradient of the
!  Ginzburg-Landau (1-dimensional) problem.

!  The subroutine statement is

!    subroutine dgl1fg(n, x, f, fgrad, task, t)

!  where

!    n is an integer variable.
!      On entry n is the number of subintervals in the domain.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
!         x is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F'
!         or 'FG'.

!    fgrad is a REAL (dp) array of dimension n.
!      On entry fgrad need not be specified.
!      On exit fgrad contains the gradient evaluated at x if
!         task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at x.
!          'G'     Evaluate the gradient at x.
!          'FG'    Evaluate the function and the gradient at x.
!          'XS'    Set x to the standard starting point xs.

!      On exit task is unchanged.

!    t is a REAL (dp) variable.
!      On entry t is a temperature in (3.73,7.32).
!      On exit t is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick, Cheryl Hile, and Jorge J. More'.

!  **********

REAL (dp), PARAMETER :: two=2.0_dp, three=3.0_dp, four=4.0_dp, ten=10.0_dp, &
                        sxteen=16.0_dp

INTEGER   :: i, n1, n2
REAL (dp) :: alphan, alphas, betan, betas, c, dn, ds, ec, em,  &
             f1, f2, f3, fac, gamma, h1, h2, hbar, hcn, hcs, penn, pens, pi, &
             tcn, tcs, temp

!     Initialization.

!     Set electron mass (grams), speed of light (cm/sec), and
!     electronic charge (esu).

em = 9.11D-28
c = 2.99D+10
ec = 4.80D-10

!     Set length of a half-layer of lead and tin (10**3-angstroms),
!     d = ds + dn.

ds = 1.0_dp
dn = 2.2_dp

!     Set critical temperature for lead and tin (Kelvin).

tcs = 7.32_dp
tcn = 3.73_dp

!     Set critical magnetic field for lead and tin at zero
!     temperature (gauss).

hcs = 803.0_dp
hcn = 309.0_dp

!     Set penetration depth for lead and tin at zero temperature
!     (cm).

pens = 3.7D-6
penn = 3.4D-6

!     Compute pi.

pi = four*ATAN(one)

!     Set initial values for temperature dependent constants alphas,
!     alphan (ergs), and betas, betan (ergs-cm**3).

alphas = -two*((ec/c)**2/em)*(hcs**2)*(pens**2)
alphan = -two*((ec/c)**2/em)*(hcn**2)*(penn**2)
betas = sxteen*pi*(((ec/c)**2/em)**2)*(hcs**2)*(pens**4)
betan = sxteen*pi*(((ec/c)**2/em)**2)*(hcn**2)*(penn**4)

alphas = alphas*((one-(t/tcs)**2)/(one+(t/tcs)**2))
alphan = alphan*((one-(t/tcn)**2)/(one+(t/tcn)**2))
betas = betas/((one+(t/tcs)**2)**2)
betan = betan/((one+(t/tcn)**2)**2)

!     Set Planck's constant (erg-sec).

hbar = 1.05459D-27

!     Set temperature dependent constant gamma (erg-cm**2).

gamma = hbar**2/(four*em)

!     Scale temperature dependent constants to the same units.
!     This makes the order parameter dimensionless.

fac = 1.0D6
alphas = alphas*(fac**3)
alphan = alphan*(fac**3)
betas = betas*(fac**6)
betan = betan*(fac**6)
gamma = gamma*(fac**5)

!     Compute the number of subintervals in (-d,-ds), in (-ds,ds),
!     and in (ds,d).

n1 = n/4
n2 = n - 2*n1
h1 = dn/DBLE(n1)
h2 = (two*ds)/DBLE(n2)

!     Compute the standard starting point if task = 'XS'.

IF (task(1:2) == 'XS') THEN
  temp = SQRT((betas+betan)/(two*(ABS(alphas)+ABS(alphan))))
  x(1:n) = temp
  RETURN
END IF

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  f1 = zero
  f2 = zero
  f3 = zero
END IF
IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') fgrad(1:n) = zero

!     Evaluate the function over the intervals (-d, -ds), (-ds, ds),
!     and (ds, d) if task = 'F' or task = 'FG'.

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  DO  i = 1, n1
    f1 = f1 + (alphan/three)*(x(i+1)**2+x(i+1)*x(i)+x(i)**2) +  &
        (betan/ten)*(x(i+1)**4+x(i+1)**3*x(i)+  &
        x(i+1)**2*x(i)**2+x(i+1)*x(i)**3+x(i)**4) + gamma*((x(i+1)-x(i))/h1)**2
  END DO
  DO  i = n1 + 1, n1 + n2
    f2 = f2 + (alphas/three)*(x(i+1)**2+x(i+1)*x(i)+x(i)**2) +  &
        (betas/ten)*(x(i+1)**4+x(i+1)**3*x(i)+  &
        x(i+1)**2*x(i)**2+x(i+1)*x(i)**3+x(i)**4) + gamma*((x(i+1)-x(i))/h2)**2
  END DO
  DO  i = n1 + n2 + 1, n - 1
    f3 = f3 + (alphan/three)*(x(i+1)**2+x(i+1)*x(i)+x(i)**2) +  &
        (betan/ten)*(x(i+1)**4+x(i+1)**3*x(i)+  &
        x(i+1)**2*x(i)**2+x(i+1)*x(i)**3+x(i)**4) + gamma*((x(i+1)-x(i))/h1)**2
  END DO
  
!        Special case for the right subinterval where x(n+1) = x(1).
  
  f3 = f3 + (alphan/three)*(x(1)**2+x(1)*x(n)+x(n)**2) +  &
      (betan/ten)*(x(1)**4+x(1)**3*x(n)+x(1)**2*x(n)**2+  &
      x(1)*x(n)**3+x(n)**4) + gamma*((x(1)-x(n))/h1)**2
  f = h1*f1 + h2*f2 + h1*f3
END IF

!     Evaluate the gradient over the intervals (-d, -ds), (-ds, ds),
!     and (ds, d).

IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') THEN
  DO  i = 1, n1
    fgrad(i) = fgrad(i) + h1*((alphan/three)*(x(i+1)+two*x(i))+  &
        (betan/ten)*(x(i+1)**3+two*(x(i+1)**2)*x(i)+  &
        three*x(i+1)*(x(i)**2)+four*x(i)**3)- gamma*(two/h1)*((x(i+1)-x(i))/h1))
    fgrad(i+1) = fgrad(i+1) + h1* ((alphan/three)*(two*x(i+1)+x(i))+  &
        (betan/ten)*(four*x(i+1)**3+three*(x(i+  &
        1)**2)*x(i)+two*x(i+1)*(x(i)**2)+x(i)**3)+  &
        gamma*(two/h1)*((x(i+1)-x(i))/h1))
  END DO
  DO  i = n1 + 1, n1 + n2
    fgrad(i) = fgrad(i) + h2*((alphas/three)*(x(i+1)+two*x(i))+  &
        (betas/ten)*(x(i+1)**3+two*(x(i+1)**2)*x(i)+  &
        three*x(i+1)*(x(i)**2)+four*x(i)**3)- gamma*(two/h2)*((x(i+1)-x(i))/h2))
    fgrad(i+1) = fgrad(i+1) + h2* ((alphas/three)*(two*x(i+1)+x(i))+  &
        (betas/ten)*(four*x(i+1)**3+three*(x(i+  &
        1)**2)*x(i)+two*x(i+1)*(x(i)**2)+x(i)**3)+  &
        gamma*(two/h2)*((x(i+1)-x(i))/h2))
  END DO
  DO  i = n1 + n2 + 1, n - 1
    fgrad(i) = fgrad(i) + h1*((alphan/three)*(x(i+1)+two*x(i))+  &
        (betan/ten)*(x(i+1)**3+two*(x(i+1)**2)*x(i)+  &
        three*x(i+1)*(x(i)**2)+four*x(i)**3)- gamma*(two/h1)*((x(i+1)-x(i))/h1))
    fgrad(i+1) = fgrad(i+1) + h1* ((alphan/three)*(two*x(i+1)+x(i))+  &
        (betan/ten)*(four*x(i+1)**3+three*(x(i+  &
        1)**2)*x(i)+two*x(i+1)*(x(i)**2)+x(i)**3)+  &
        gamma*(two/h1)*((x(i+1)-x(i))/h1))
  END DO
  
!        Special case for the right subinterval where x(n+1) = x(1).
  
  fgrad(n) = fgrad(n) + h1*((alphan/three)*(x(1)+two*x(n))+  &
      (betan/ten)*(x(1)**3+two*(x(1)**2)*x(n)+  &
      three*x(1)*(x(n)**2)+four*x(n)**3)- gamma*(two/h1)*((x(1)-x(n))/h1))
  fgrad(1) = fgrad(1) + h1*((alphan/three)*(two*x(1)+x(n))+  &
      (betan/ten)*(four*x(1)**3+three*(x(1)**2)*x(n)+  &
      two*x(1)*(x(n)**2)+x(n)**3)+ gamma*(two/h1)*((x(1)-x(n))/h1))
END IF

RETURN
END SUBROUTINE dgl1fg



SUBROUTINE dgl2fg(nx, ny, x, f, fgrad, task, vornum)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:08

INTEGER, INTENT(IN)            :: nx
INTEGER, INTENT(IN)            :: ny
REAL (dp), INTENT(IN OUT)      :: x(:)    ! x(4*nx*ny)
REAL (dp), INTENT(OUT)         :: f
REAL (dp), INTENT(OUT)         :: fgrad(:)    ! fgrad(4*nx*ny)
CHARACTER (LEN=*), INTENT(IN)  :: task
INTEGER, INTENT(IN)            :: vornum

!  **********

!  Subroutine dgl2fg

!  This subroutine computes the function and gradient of the
!  Ginzburg-Landau (2-dimensional) superconductivity problem.

!  The subroutine statement is

!    subroutine dgl2fg(nx, ny, x, f, fgrad, task, w, vornum)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension 4*nx*ny.
!      On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
!         x is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    fgrad is a REAL (dp) array of dimension 4*nx*ny.
!      On entry fgrad need not be specified.
!      On exit fgrad contains the gradient evaluated at x if
!         task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at x.
!          'G'     Evaluate the gradient at x.
!          'FG'    Evaluate the function and the gradient at x.
!          'XS'    Set x to the standard starting point xs.

!      On exit task is unchanged.

!    vornum is an integer variable.
!      On entry vornum specifies the number of vortices.
!      On exit vornum is unchanged.

!  Subprograms called

!    MINPACK-supplied   ...   dgl2fc

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick, Paul L. Plassmann, and Stephen J. Wright.

!  **********

REAL (dp) :: wx(nx+1,ny+1), wy(nx+1,ny+1), vpotx(nx+1,ny+1),  &
             vpoty(nx+1,ny+1), gradx(nx,ny), grady(nx,ny), gradax(nx,ny),  &
             graday(nx,ny)
INTEGER   :: j, k1, k2, k3, k4

! Copy elements of x to arrays wx, wy, vpotx & vpoty

k1 = 0
k2 = nx*ny
k3 = k2 + nx*ny
k4 = k3 + nx*ny
DO j = 1, ny
  wx(1:nx,j) = x(k1+1:k1+nx)
  wy(1:nx,j) = x(k2+1:k2+nx)
  vpotx(1:nx,j) = x(k3+1:k3+nx)
  vpoty(1:nx,j) = x(k4+1:k4+nx)
  k1 = k1 + nx
  k2 = k2 + nx
  k3 = k3 + nx
  k4 = k4 + nx
END DO

wx(:,ny+1) = zero
wx(nx+1,:) = zero
wy(:,ny+1) = zero
wy(nx+1,:) = zero
vpotx(:,ny+1) = zero
vpotx(nx+1,:) = zero
vpoty(:,ny+1) = zero
vpoty(nx+1,:) = zero

CALL dgl2fc(nx, ny, wx, wy, vpotx, vpoty, f, gradx, grady, gradax,  &
                  graday, task, vornum)

! Copy elements back into x & fgrad.

k1 = 0
k2 = nx*ny
k3 = k2 + nx*ny
k4 = k3 + nx*ny
DO j = 1, ny
  x(k1+1:k1+nx) = wx(1:nx,j)
  x(k2+1:k2+nx) = wy(1:nx,j)
  x(k3+1:k3+nx) = vpotx(1:nx,j)
  x(k4+1:k4+nx) = vpoty(1:nx,j)
  fgrad(k1+1:k1+nx) = gradx(1:nx,j)
  fgrad(k2+1:k2+nx) = grady(1:nx,j)
  fgrad(k3+1:k3+nx) = gradax(1:nx,j)
  fgrad(k4+1:k4+nx) = graday(1:nx,j)
  k1 = k1 + nx
  k2 = k2 + nx
  k3 = k3 + nx
  k4 = k4 + nx
END DO

RETURN
END SUBROUTINE dgl2fg



SUBROUTINE dgl2fc(nx, ny, x, y, vpotx, vpoty, f, gradx, grady, gradax,  &
                  graday, task, vornum)

INTEGER, INTENT(IN)            :: nx
INTEGER, INTENT(IN)            :: ny
REAL (dp), INTENT(IN OUT)      :: x(:,:)        ! x(nx+1,ny+1)
REAL (dp), INTENT(IN OUT)      :: y(:,:)        ! y(nx+1,ny+1)
REAL (dp), INTENT(OUT)         :: vpotx(:,:)    ! vpotx(nx+1,ny+1)
REAL (dp), INTENT(OUT)         :: vpoty(:,:)    ! vpoty(nx+1,ny+1)
REAL (dp), INTENT(OUT)         :: f
REAL (dp), INTENT(OUT)         :: gradx(:,:)    ! gradx(nx,ny)
REAL (dp), INTENT(OUT)         :: grady(:,:)    ! grady(nx,ny)
REAL (dp), INTENT(OUT)         :: gradax(:,:)   ! gradax(nx,ny)
REAL (dp), INTENT(OUT)         :: graday(:,:)   ! graday(nx,ny)
CHARACTER (LEN=*), INTENT(IN)  :: task
INTEGER, INTENT(IN)            :: vornum

!  **********

!  Subroutine dgl2fc

!  This subroutine computes the function and gradient of the
!  Ginzburg-Landau (2-dimensional) superconductivity problem.

!  The subroutine statement is

!    subroutine dgl2fc(nx, ny, x, y, vpotx, vpoty, f,
! +                    gradx, grady, gradax, graday, task, vornum)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the real part of the order parameter
!         if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
!         x is set according to task.

!    y is a REAL (dp) array of dimension nx*ny.
!      On entry y specifies the imaginary part of the order parameter
!         if task = 'F', 'G', or 'FG'.
!         Otherwise y need not be specified.
!      On exit y is unchanged if task = 'F', 'G', or 'FG'. Otherwise
!         y is set according to task.

!    vpotx is a REAL (dp) array of dimension nx*ny.
!      On entry vpotx specifies the x component of the vector
!         potential if task = 'F', 'G', or 'FG'.
!         Otherwise vpotx need not be specified.
!      On exit vpotx is unchanged if task = 'F', 'G', or 'FG'.
!         Otherwise vpotx is set according to task.

!    vpoty is a REAL (dp) array of dimension nx*ny.
!      On entry vpoty specifies the y component of the vector
!         potential if task = 'F', 'G', or 'FG'.
!         Otherwise vpoty need not be specified.
!      On exit vpoty is unchanged if task = 'F', 'G', or 'FG'.
!         Otherwise vpoty is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    gradx is a REAL (dp) array of dimension nx*ny.
!      On entry gradx need not be specified.
!      On exit gradx contains the gradient with respect to x
!         of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.

!    grady is a REAL (dp) array of dimension nx*ny.
!      On entry grady need not be specified.
!      On exit grady contains the gradient with respect to y
!         of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
!         if task = 'G' or 'FG'.

!    gradax is a REAL (dp) array of dimension nx*ny.
!      On entry gradax need not be specified.
!      On exit gradax contains the gradient with respect to vpotx
!         of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
!         if task = 'G' or 'FG'.

!    graday is a REAL (dp) array of dimension nx*ny.
!      On entry graday need not be specified.
!      On exit graday contains the gradient with respect to vpoty
!         of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
!         if task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at (x,y,vpotx,vpoty).
!          'G'     Evaluate the gradient at (x,y,vpotx,vpoty).
!          'FG'    Evaluate the function and the gradient at
!                      (x,y,vpotx,vpoty).
!          'XS'    Set (x,y,vpotx,vpoty) to the standard starting
!                      point xs.

!      On exit task is unchanged.

!    vornum is an integer variable.
!      On entry vornum is the number of vortices.
!      On exit vornum is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick, Paul L. Plassmann, and Stephen J. Wright.

!  **********

REAL (dp), PARAMETER :: two=2.0_dp, three=3.0_dp, four=4.0_dp, five=5.0_dp

INTEGER   :: i, j
REAL (dp) :: arg, bave, cfac, delsq, fcond, ffield, fkin,  &
             fkinx1, fkinx2, fkiny1, fkiny2, hx, hy, pi, sfac,  &
             sqn, sqrtv, tkappa, x1, x2, xpt, xy, ypt

!     Initialize.

tkappa = five
hx = SQRT(vornum/two)*three/DBLE(nx)
hy = SQRT(vornum/two)*three*SQRT(three)/DBLE(ny)
sqn = DBLE(nx*ny)
pi = four*ATAN(one)
bave = two*pi*vornum*tkappa/(sqn*hx*hy)
sqrtv = SQRT(DBLE(vornum))*pi

IF (task(1:2) == 'XS') THEN
  
!        Initial Order Parameter.
  
  DO  j = 1, ny + 1
    ypt = (DBLE(j)-one)*hy
    DO  i = 1, nx + 1
      xpt = (DBLE(i)-one)*hx
      x(i,j) = one - (SIN(sqrtv*xpt/(two*three))*  &
          SIN(sqrtv*ypt/(two*SQRT(three)*three)))**2
      y(i,j) = zero
    END DO
  END DO
  
!        Initial Vector Potential.
  
  DO  j = 1, ny + 1
    DO  i = 1, nx + 1
      xpt = (DBLE(i)-one)*hx
      vpotx(i,j) = zero
      vpoty(i,j) = bave*xpt/tkappa
    END DO
  END DO
  
  RETURN
  
END IF

!     Enforce vortex constraint and boundary conditions.

!     Right face for order parameter and vector potential.

DO  j = 1, ny + 1
  arg = two*pi*vornum*(DBLE(j)-one)/DBLE(ny)
  x(nx+1,j) = x(1,j)*COS(arg) - y(1,j)*SIN(arg)
  y(nx+1,j) = x(1,j)*SIN(arg) + y(1,j)*COS(arg)
  vpotx(nx+1,j) = vpotx(1,j)
  vpoty(nx+1,j) = vpoty(1,j) + two*pi*vornum/(DBLE(ny)*hy)
END DO

!     Top face for order parameter and vector potential.

DO  i = 1, nx + 1
  x(i,ny+1) = x(i,1)
  y(i,ny+1) = y(i,1)
  vpotx(i,ny+1) = vpotx(i,1)
  vpoty(i,ny+1) = vpoty(i,1)
END DO

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  
!        Compute the Condensation Energy Density
  
  fcond = zero
  DO  i = 1, nx
    DO  j = 1, ny
      delsq = x(i,j)**2 + y(i,j)**2
      fcond = fcond - delsq + (delsq**2)/two
    END DO
  END DO
  fcond = fcond/sqn
  
!        Compute the Kinetic Energy Density.
  
  fkin = zero
  DO  i = 1, nx
    DO  j = 1, ny
      x1 = x(i+1,j) - x(i,j)*COS(hx*vpotx(i,j)) + y(i,j)*SIN(hx*vpotx(i,j))
      x2 = y(i+1,j) - y(i,j)*COS(hx*vpotx(i,j)) - x(i,j)*SIN(hx*vpotx(i,j))
      fkin = fkin + (x1**2+x2**2)/(hx**2)
      x1 = x(i,j+1) - x(i,j)*COS(hy*vpoty(i,j)) + y(i,j)*SIN(hy*vpoty(i,j))
      x2 = y(i,j+1) - y(i,j)*COS(hy*vpoty(i,j)) - x(i,j)*SIN(hy*vpoty(i,j))
      fkin = fkin + (x1**2+x2**2)/(hy**2)
    END DO
  END DO
  fkin = fkin/sqn
  
!        Compute the Magnetic Field Energy Density.
  
  ffield = zero
  DO  i = 1, nx
    DO  j = 1, ny
      xy = (vpoty(i+1,j)-vpoty(i,j))/hx - (vpotx(i,j+1)-vpotx(i,j))/hy
      ffield = ffield + xy**2
    END DO
  END DO
  ffield = ffield*(tkappa**2)/sqn
  f = fcond + fkin + ffield
END IF

IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') THEN
  
  DO  j = 1, ny
    DO  i = 1, nx
      gradx(i,j) = x(i,j)*(-one+x(i,j)**2+y(i,j)**2)
      gradx(i,j) = gradx(i,j)*two/sqn
      grady(i,j) = y(i,j)*(-one+x(i,j)**2+y(i,j)**2)
      grady(i,j) = grady(i,j)*two/sqn
      gradax(i,j) = zero
      graday(i,j) = zero
    END DO
  END DO
  
!        Kinetic Energy Part, Interior Points
  
  DO  i = 2, nx
    DO  j = 2, ny
      fkinx1 = (two/(hx*hx*sqn))*(x(i+1,j)- x(i,j)*COS(hx*vpotx(i,j))+  &
          y(i,j)*SIN(hx*vpotx(i,j)))
      fkinx2 = (two/(hx*hx*sqn))*(y(i+1,j)- y(i,j)*COS(hx*vpotx(i,j))-  &
          x(i,j)*SIN(hx*vpotx(i,j)))
      fkiny1 = (two/(hy*hy*sqn))*(x(i,j+1)- x(i,j)*COS(hy*vpoty(i,j))+  &
          y(i,j)*SIN(hy*vpoty(i,j)))
      fkiny2 = (two/(hy*hy*sqn))*(y(i,j+1)- y(i,j)*COS(hy*vpoty(i,j))-  &
          x(i,j)*SIN(hy*vpoty(i,j)))
      ffield = (vpotx(i,j)-vpotx(i,j+1))/hy + (vpoty(i+1,j)-vpoty(i,j))/hx
      ffield = (two*(tkappa**2)/sqn)*ffield
      gradx(i,j) = gradx(i,j) - COS(hx*vpotx(i,j))*fkinx1 -  &
          SIN(hx*vpotx(i,j))*fkinx2 - COS(hy*vpoty(i,j))*fkiny1 -  &
          SIN(hy*vpoty(i,j))*fkiny2
      grady(i,j) = grady(i,j) + SIN(hx*vpotx(i,j))*fkinx1 -  &
          COS(hx*vpotx(i,j))*fkinx2 + SIN(hy*vpoty(i,j))*fkiny1 -  &
          COS(hy*vpoty(i,j))*fkiny2
      gradax(i,j) = gradax(i,j) + ffield/hy +  &
          fkinx1*(hx*x(i,j)*SIN(hx*vpotx(i,j))+ hx*y(i,j)*COS(hx*vpotx(i,j))) +  &
          fkinx2*(hx*y(i,j)*SIN(hx*vpotx(i,j))- hx*x(i,j)*COS(hx*vpotx(i,j)))
      graday(i,j) = graday(i,j) - ffield/hx +  &
          fkiny1*(hy*x(i,j)*SIN(hy*vpoty(i,j))+ hy*y(i,j)*COS(hy*vpoty(i,j))) +  &
          fkiny2*(hy*y(i,j)*SIN(hy*vpoty(i,j))- hy*x(i,j)*COS(hy*vpoty(i,j)))
      fkinx1 = (two/(hx*hx*sqn))*(x(i,j)- x(i-1,j)*COS(hx*vpotx(i-1,j))+  &
          y(i-1,j)*SIN(hx*vpotx(i-1,j)))
      fkinx2 = (two/(hx*hx*sqn))*(y(i,j)- y(i-1,j)*COS(hx*vpotx(i-1,j))-  &
          x(i-1,j)*SIN(hx*vpotx(i-1,j)))
      fkiny1 = (two/(hy*hy*sqn))*(x(i,j)- x(i,j-1)*COS(hy*vpoty(i,j-1))+  &
          y(i,j-1)*SIN(hy*vpoty(i,j-1)))
      fkiny2 = (two/(hy*hy*sqn))*(y(i,j)- y(i,j-1)*COS(hy*vpoty(i,j-1))-  &
          x(i,j-1)*SIN(hy*vpoty(i,j-1)))
      gradx(i,j) = gradx(i,j) + fkinx1 + fkiny1
      grady(i,j) = grady(i,j) + fkinx2 + fkiny2
      ffield = (vpotx(i,j-1)-vpotx(i,j))/hy + (vpoty(i+1,j-1)-vpoty(i,j-1))/hx
      ffield = (two*(tkappa**2)/sqn)*ffield
      gradax(i,j) = gradax(i,j) - ffield/hy
      ffield = (vpotx(i-1,j)-vpotx(i-1,j+1))/hy + (vpoty(i,j)-vpoty(i-1,j))/hx
      ffield = (two*(tkappa**2)/sqn)*ffield
      graday(i,j) = graday(i,j) + ffield/hx
    END DO
  END DO
  
!        Kinetic Energy Part, Boundary Points.
  
!        Bottom J = 1
  
  DO  i = 2, nx
    fkinx1 = (two/(hx*hx*sqn))*(x(i+1,1)- x(i,1)*COS(hx*vpotx(i,1))+  &
        y(i,1)*SIN(hx*vpotx(i,1)))
    fkinx2 = (two/(hx*hx*sqn))*(y(i+1,1)- y(i,1)*COS(hx*vpotx(i,1))-  &
        x(i,1)*SIN(hx*vpotx(i,1)))
    fkiny1 = (two/(hy*hy*sqn))*(x(i,2)- x(i,1)*COS(hy*vpoty(i,1))+  &
        y(i,1)*SIN(hy*vpoty(i,1)))
    fkiny2 = (two/(hy*hy*sqn))*(y(i,2)- y(i,1)*COS(hy*vpoty(i,1))-  &
        x(i,1)*SIN(hy*vpoty(i,1)))
    ffield = (vpotx(i,1)-vpotx(i,2))/hy + (vpoty(i+1,1)-vpoty(i,1))/hx
    ffield = (two*(tkappa**2)/sqn)*ffield
    gradx(i,1) = gradx(i,1) - COS(hx*vpotx(i,1))*fkinx1 -  &
        SIN(hx*vpotx(i,1))*fkinx2 - COS(hy*vpoty(i,1))*fkiny1 -  &
        SIN(hy*vpoty(i,1))*fkiny2
    grady(i,1) = grady(i,1) + SIN(hx*vpotx(i,1))*fkinx1 -  &
        COS(hx*vpotx(i,1))*fkinx2 + SIN(hy*vpoty(i,1))*fkiny1 -  &
        COS(hy*vpoty(i,1))*fkiny2
    gradax(i,1) = gradax(i,1) + ffield/hy +  &
        fkinx1*(hx*x(i,1)*SIN(hx*vpotx(i,1))+ hx*y(i,1)*COS(hx*vpotx(i,1))) +  &
        fkinx2*(hx*y(i,1)*SIN(hx*vpotx(i,1))- hx*x(i,1)*COS(hx*vpotx(i,1)))
    graday(i,1) = graday(i,1) - ffield/hx +  &
        fkiny1*(hy*x(i,1)*SIN(hy*vpoty(i,1))+ hy*y(i,1)*COS(hy*vpoty(i,1))) +  &
        fkiny2*(hy*y(i,1)*SIN(hy*vpoty(i,1))- hy*x(i,1)*COS(hy*vpoty(i,1)))
    fkinx1 = (two/(hx*hx*sqn))*(x(i,1)- x(i-1,1)*COS(hx*vpotx(i-1,1))+  &
        y(i-1,1)*SIN(hx*vpotx(i-1,1)))
    fkinx2 = (two/(hx*hx*sqn))*(y(i,1)- y(i-1,1)*COS(hx*vpotx(i-1,1))-  &
        x(i-1,1)*SIN(hx*vpotx(i-1,1)))
    fkiny1 = (two/(hy*hy*sqn))*(x(i,ny+1)- x(i,ny)*COS(hy*vpoty(i,ny))+  &
        y(i,ny)*SIN(hy*vpoty(i,ny)))
    fkiny2 = (two/(hy*hy*sqn))*(y(i,ny+1)- y(i,ny)*COS(hy*vpoty(i,ny))-  &
        x(i,ny)*SIN(hy*vpoty(i,ny)))
    gradx(i,1) = gradx(i,1) + fkinx1 + fkiny1
    grady(i,1) = grady(i,1) + fkinx2 + fkiny2
    ffield = (vpotx(i,ny)-vpotx(i,ny+1))/hy + (vpoty(i+1,ny)-vpoty(i,ny))/hx
    ffield = (two*(tkappa**2)/sqn)*ffield
    gradax(i,1) = gradax(i,1) - ffield/hy
    ffield = (vpotx(i-1,1)-vpotx(i-1,2))/hy + (vpoty(i,1)-vpoty(i-1,1))/hx
    ffield = (two*(tkappa**2)/sqn)*ffield
    graday(i,1) = graday(i,1) + ffield/hx
  END DO
  
!        Left I = 1.
  
  DO  j = 2, ny
    fkinx1 = (two/(hx*hx*sqn))*(x(2,j)- x(1,j)*COS(hx*vpotx(1,j))+  &
        y(1,j)*SIN(hx*vpotx(1,j)))
    fkinx2 = (two/(hx*hx*sqn))*(y(2,j)- y(1,j)*COS(hx*vpotx(1,j))-  &
        x(1,j)*SIN(hx*vpotx(1,j)))
    fkiny1 = (two/(hy*hy*sqn))*(x(1,j+1)- x(1,j)*COS(hy*vpoty(1,j))+  &
        y(1,j)*SIN(hy*vpoty(1,j)))
    fkiny2 = (two/(hy*hy*sqn))*(y(1,j+1)- y(1,j)*COS(hy*vpoty(1,j))-  &
        x(1,j)*SIN(hy*vpoty(1,j)))
    ffield = (vpotx(1,j)-vpotx(1,j+1))/hy + (vpoty(2,j)-vpoty(1,j))/hx
    ffield = (two*(tkappa**2)/sqn)*ffield
    gradx(1,j) = gradx(1,j) - COS(hx*vpotx(1,j))*fkinx1 -  &
        SIN(hx*vpotx(1,j))*fkinx2 - COS(hy*vpoty(1,j))*fkiny1 -  &
        SIN(hy*vpoty(1,j))*fkiny2
    grady(1,j) = grady(1,j) + SIN(hx*vpotx(1,j))*fkinx1 -  &
        COS(hx*vpotx(1,j))*fkinx2 + SIN(hy*vpoty(1,j))*fkiny1 -  &
        COS(hy*vpoty(1,j))*fkiny2
    gradax(1,j) = gradax(1,j) + ffield/hy +  &
        fkinx1*(hx*x(1,j)*SIN(hx*vpotx(1,j))+ hx*y(1,j)*COS(hx*vpotx(1,j))) +  &
        fkinx2*(hx*y(1,j)*SIN(hx*vpotx(1,j))- hx*x(1,j)*COS(hx*vpotx(1,j)))
    graday(1,j) = graday(1,j) - ffield/hx +  &
        fkiny1*(hy*x(1,j)*SIN(hy*vpoty(1,j))+ hy*y(1,j)*COS(hy*vpoty(1,j))) +  &
        fkiny2*(hy*y(1,j)*SIN(hy*vpoty(1,j))- hy*x(1,j)*COS(hy*vpoty(1,j)))
    fkinx1 = (two/(hx*hx*sqn))*(x(nx+1,j)- x(nx,j)*COS(hx*vpotx(nx,j))+  &
        y(nx,j)*SIN(hx*vpotx(nx,j)))
    fkinx2 = (two/(hx*hx*sqn))*(y(nx+1,j)- y(nx,j)*COS(hx*vpotx(nx,j))-  &
        x(nx,j)*SIN(hx*vpotx(nx,j)))
    fkiny1 = (two/(hy*hy*sqn))*(x(1,j)- x(1,j-1)*COS(hy*vpoty(1,j-1))+  &
        y(1,j-1)*SIN(hy*vpoty(1,j-1)))
    fkiny2 = (two/(hy*hy*sqn))*(y(1,j)- y(1,j-1)*COS(hy*vpoty(1,j-1))-  &
        x(1,j-1)*SIN(hy*vpoty(1,j-1)))
    sfac = SIN(two*pi*vornum*(j-one)/DBLE(ny))
    cfac = COS(two*pi*vornum*(j-one)/DBLE(ny))
    gradx(1,j) = gradx(1,j) + cfac*fkinx1 + sfac*fkinx2 + fkiny1
    grady(1,j) = grady(1,j) - sfac*fkinx1 + cfac*fkinx2 + fkiny2
    ffield = (vpotx(1,j-1)-vpotx(1,j))/hy + (vpoty(2,j-1)-vpoty(1,j-1))/hx
    ffield = (two*(tkappa**2)/sqn)*ffield
    gradax(1,j) = gradax(1,j) - ffield/hy
    ffield = (vpotx(nx,j)-vpotx(nx,j+1))/hy + (vpoty(nx+1,j)-vpoty(nx,j))/hx
    ffield = (two*(tkappa**2)/sqn)*ffield
    graday(1,j) = graday(1,j) + ffield/hx
  END DO
  
!        Kinetic Energy Part, at origin (only needed in zero field).
  
  fkinx1 = (two/(hx*hx*sqn))*(x(2,1)-x(1,1)*COS(hx*vpotx(1,1))+  &
      y(1,1)*SIN(hx*vpotx(1,1)))
  fkinx2 = (two/(hx*hx*sqn))*(y(2,1)-y(1,1)*COS(hx*vpotx(1,1))-  &
      x(1,1)*SIN(hx*vpotx(1,1)))
  fkiny1 = (two/(hy*hy*sqn))*(x(1,2)-x(1,1)*COS(hy*vpoty(1,1))+  &
      y(1,1)*SIN(hy*vpoty(1,1)))
  fkiny2 = (two/(hy*hy*sqn))*(y(1,2)-y(1,1)*COS(hy*vpoty(1,1))-  &
      x(1,1)*SIN(hy*vpoty(1,1)))
  ffield = (vpotx(1,1)-vpotx(1,2))/hy + (vpoty(2,1)-vpoty(1,1))/hx
  ffield = (two*(tkappa**2)/sqn)*ffield
  gradx(1,1) = gradx(1,1) - COS(hx*vpotx(1,1))*fkinx1 -  &
      SIN(hx*vpotx(1,1))*fkinx2 - COS(hy*vpoty(1,1))*fkiny1 -  &
      SIN(hy*vpoty(1,1))*fkiny2
  grady(1,1) = grady(1,1) + SIN(hx*vpotx(1,1))*fkinx1 -  &
      COS(hx*vpotx(1,1))*fkinx2 + SIN(hy*vpoty(1,1))*fkiny1 -  &
      COS(hy*vpoty(1,1))*fkiny2
  gradax(1,1) = gradax(1,1) + ffield/hy +  &
      fkinx1*(hx*x(1,1)*SIN(hx*vpotx(1,1))+ hx*y(1,1)*COS(hx*vpotx(1,1))) +  &
      fkinx2*(hx*y(1,1)*SIN(hx*vpotx(1,1))- hx*x(1,1)*COS(hx*vpotx(1,1)))
  graday(1,1) = graday(1,1) - ffield/hx +  &
      fkiny1*(hy*x(1,1)*SIN(hy*vpoty(1,1))+ hy*y(1,1)*COS(hy*vpoty(1,1))) +  &
      fkiny2*(hy*y(1,1)*SIN(hy*vpoty(1,1))- hy*x(1,1)*COS(hy*vpoty(1,1)))
  fkinx1 = (two/(hx*hx*sqn))*(x(nx+1,1)- x(nx,1)*COS(hx*vpotx(nx,1))+  &
      y(nx,1)*SIN(hx*vpotx(nx,1)))
  fkinx2 = (two/(hx*hx*sqn))*(y(nx+1,1)- y(nx,1)*COS(hx*vpotx(nx,1))-  &
      x(nx,1)*SIN(hx*vpotx(nx,1)))
  fkiny1 = (two/(hy*hy*sqn))*(x(1,ny+1)- x(1,ny)*COS(hy*vpoty(1,ny))+  &
      y(1,ny)*SIN(hy*vpoty(1,ny)))
  fkiny2 = (two/(hy*hy*sqn))*(y(1,ny+1)- y(1,ny)*COS(hy*vpoty(1,ny))-  &
      x(1,ny)*SIN(hy*vpoty(1,ny)))
  gradx(1,1) = gradx(1,1) + fkinx1 + fkiny1
  grady(1,1) = grady(1,1) + fkinx2 + fkiny2
  ffield = (vpotx(1,ny)-vpotx(1,ny+1))/hy + (vpoty(2,ny)-vpoty(1,ny))/hx
  ffield = (two*(tkappa**2)/sqn)*ffield
  gradax(1,1) = gradax(1,1) - ffield/hy
  ffield = (vpotx(nx,1)-vpotx(nx,2))/hy + (vpoty(nx+1,1)-vpoty(nx,1))/hx
  ffield = (two*(tkappa**2)/sqn)*ffield
  graday(1,1) = graday(1,1) + ffield/hx
END IF

RETURN
END SUBROUTINE dgl2fc



SUBROUTINE dljcfg(n, x, f, fgrad, task, ndim, natoms)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:08

INTEGER, INTENT(IN)                :: n
REAL (dp), INTENT(OUT)             :: x(:)
REAL (dp), INTENT(OUT)             :: f
REAL (dp), INTENT(OUT)             :: fgrad(:)
CHARACTER (LEN=*), INTENT(IN OUT)  :: task
INTEGER, INTENT(IN)                :: ndim
INTEGER, INTENT(IN)                :: natoms

!  **********

!  Subroutine dljcfg

!  This subroutine computes the function and gradient of the
!  Leonard-Jones clusters (molecular conformation) problem.

!  The subroutine statement is

!    subroutine dljcfg(n, x, f, fgrad, task, natoms, ndim)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!         For the 2-dimensional problem n = 2*natoms.
!         For the 3-dimensional problem n = 3*natoms.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
!         Otherwise x need not be specified.
!      On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
!         x is set according to task.

!    f is a REAL (dp) variable.
!      On entry f need not be specified.
!      On exit f is set to the function evaluated at x if task = 'F' or 'FG'.

!    fgrad is a REAL (dp) array of dimension n.
!      On entry fgrad need not be specified.
!      On exit fgrad contains the gradient evaluated at x if
!         task = 'G' or 'FG'.

!    task is a character*60 variable.
!      On entry task specifies the action of the subroutine:

!         task               action
!         ----               ------
!          'F'     Evaluate the function at x.
!          'G'     Evaluate the gradient at x.
!          'FG'    Evaluate the function and the gradient at x.
!          'XS'    Set x to the standard starting point xs.

!      On exit task is unchanged.  WRONG!

!    natoms is an integer variable.
!      On entry natoms is the number of atoms in the cluster.
!      On exit natoms is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick, R. S. Maier and G. L. Xue

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp, two=2.0_dp, three=3.0_dp, six=6.0_dp

INTEGER   :: ctr, i, icrtn, il, ileft, isqrtn, j, jl, k
REAL (dp) :: rij, temp, xx, yy, zz

!     Check input for errors.

IF (ndim == 2 .AND. n /= 2*natoms) task = 'ERROR: N MUST .EQ. 2*NATOMS'
IF (ndim == 3 .AND. n /= 3*natoms) task = 'ERROR: N MUST .EQ. 3*NATOMS'
IF (ndim /= 2 .AND. ndim /= 3) task = 'ERROR: NDIM MUST .EQ. 2 OR 3'
IF (task(1:5) == 'ERROR') RETURN

!     Compute the standard starting point if task = 'XS'

IF (task(1:2) == 'XS') THEN
  IF (ndim == 2) THEN
    isqrtn = INT(SQRT(DBLE(natoms)))
    ileft = natoms - isqrtn*isqrtn
    xx = zero
    yy = zero
    DO  j = 1, isqrtn + ileft/isqrtn + 1
      DO  i = 1, MIN(isqrtn,natoms-(j-1)*isqrtn)
        ctr = (j-1)*isqrtn + i
        x(2*ctr-1) = xx
        x(2*ctr) = yy
        xx = xx + one
      END DO
      yy = yy + one
      xx = zero
    END DO
  ELSE IF (ndim == 3) THEN
    icrtn = INT((natoms+p5)**(one/three))
    ileft = natoms - icrtn*icrtn*icrtn
    xx = zero
    yy = zero
    zz = zero
    DO  k = 1, icrtn + ileft/(icrtn*icrtn) + 1
      jl = MIN(icrtn,(natoms-(k-1)*icrtn*icrtn)/icrtn+1)
      DO  j = 1, jl
        il = MIN(icrtn,natoms-(k-1)*icrtn*icrtn-(j-1)*icrtn)
        DO  i = 1, il
          ctr = (k-1)*icrtn*icrtn + (j-1)*icrtn + i
          x(3*ctr-2) = xx
          x(3*ctr-1) = yy
          x(3*ctr) = zz
          xx = xx + one
        END DO
        yy = yy + one
        xx = zero
      END DO
      yy = zero
      zz = zz + one
    END DO
  END IF
  
  RETURN
  
END IF

!     Evaluate the function if task = 'F', the gradient if task = 'G',
!     or both if task = 'FG'.

IF (task(1:1) == 'F' .OR. task(1:2) == 'FG') THEN
  f = zero
  IF (ndim == 2) THEN
    DO  j = 2, natoms
      DO  i = 1, j - 1
        xx = x(2*j-1) - x(2*i-1)
        yy = x(2*j) - x(2*i)
        rij = xx**2 + yy**2
        temp = one/rij/rij/rij
        f = f + temp*(temp-two)
      END DO
    END DO
  ELSE IF (ndim == 3) THEN
    DO  j = 2, natoms
      DO  i = 1, j - 1
        xx = x(3*j-2) - x(3*i-2)
        yy = x(3*j-1) - x(3*i-1)
        zz = x(3*j) - x(3*i)
        rij = xx*xx + yy*yy + zz*zz
        temp = one/rij/rij/rij
        f = f + temp*(temp-two)
      END DO
    END DO
  END IF
END IF

!     Compute the gradient

IF (task(1:1) == 'G' .OR. task(1:2) == 'FG') THEN
  DO  i = 1, n
    fgrad(i) = zero
  END DO
  IF (ndim == 2) THEN
    DO  j = 2, natoms
      DO  i = 1, j - 1
        xx = x(2*j-1) - x(2*i-1)
        yy = x(2*j) - x(2*i)
        rij = xx**2 + yy**2
        temp = one/rij/rij/rij
        temp = six*temp*(two*temp-two)/rij
        fgrad(2*j-1) = fgrad(2*j-1) - xx*temp
        fgrad(2*j) = fgrad(2*j) - yy*temp
        fgrad(2*i-1) = fgrad(2*i-1) + xx*temp
        fgrad(2*i) = fgrad(2*i) + yy*temp
      END DO
    END DO
  ELSE IF (ndim == 3) THEN
    DO  j = 2, natoms
      DO  i = 1, j - 1
        xx = x(3*j-2) - x(3*i-2)
        yy = x(3*j-1) - x(3*i-1)
        zz = x(3*j) - x(3*i)
        rij = xx*xx + yy*yy + zz*zz
        temp = one/rij/rij/rij
        temp = six*temp*(two*temp-two)/rij
        fgrad(3*j-2) = fgrad(3*j-2) - xx*temp
        fgrad(3*j-1) = fgrad(3*j-1) - yy*temp
        fgrad(3*j) = fgrad(3*j) - zz*temp
        fgrad(3*i-2) = fgrad(3*i-2) + xx*temp
        fgrad(3*i-1) = fgrad(3*i-1) + yy*temp
        fgrad(3*i) = fgrad(3*i) + zz*temp
      END DO
    END DO
  END IF
END IF

RETURN
END SUBROUTINE dljcfg



SUBROUTINE depths(nx, ny, s, y)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)     :: nx
INTEGER, INTENT(IN)     :: ny
REAL (dp), INTENT(IN)   :: s(:)    ! s(nx*ny)
REAL (dp), INTENT(OUT)  :: y(:)    ! y(nx*ny)

!  **********

!  Subroutine depths

!  This subroutine computes the product H*s = y, where H is the
!  Hessian matrix for the elastic plastic torsion problem.

!  The subroutine statement is

!    subroutine depths(nx, ny, s, y)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    s is a REAL (dp) array of dimension nx*ny.
!      On entry s contains the vector s.
!      On exit s is unchanged.

!    y is a REAL (dp) array of dimension nx*ny.
!      On entry out need not be specified.
!      On exit y contains H*s.


!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp

INTEGER   :: i, j, k, km1, kmnx, kp1, kpnx
REAL (dp) :: area, hx, hxhx, hy, hyhy, v, vb, vl, vr, vt

hx = one/DBLE(nx+1)
hy = one/DBLE(ny+1)
hxhx = one/(hx*hx)
hyhy = one/(hy*hy)
area = p5*hx*hy

y(1:nx*ny) = zero

!     Computation of H*s over the lower triangular elements.

DO  j = 0, ny
  k = nx*(j-1)
  DO  i = 0, nx
    v = zero
    vr = zero
    vt = zero
    IF (i /= 0 .AND. j /= 0) v = s(k)
    IF (i /= nx .AND. j /= 0) THEN
      kp1 = k + 1
      vr = s(kp1)
      y(kp1) = y(kp1) + hxhx*(vr-v)
    END IF
    IF (i /= 0 .AND. j /= ny) THEN
      kpnx = k + nx
      vt = s(kpnx)
      y(kpnx) = y(kpnx) + hyhy*(vt-v)
    END IF
    IF (i /= 0 .AND. j /= 0) y(k) = y(k) + hxhx*(v-vr) + hyhy*(v-vt)
    k = k + 1
  END DO
END DO

!     Computation of H*s over the upper triangular elements.

DO  j = 1, ny + 1
  k = nx*(j-1)
  DO  i = 1, nx + 1
    k = k + 1
    v = zero
    vl = zero
    vb = zero
    IF (i /= nx+1 .AND. j /= ny+1) v = s(k)
    IF (i /= nx+1 .AND. j /= 1) THEN
      kmnx = k - nx
      vb = s(kmnx)
      y(kmnx) = y(kmnx) + hyhy*(vb-v)
    END IF
    IF (i /= 1 .AND. j /= ny+1) THEN
      km1 = k - 1
      vl = s(km1)
      y(km1) = y(km1) + hxhx*(vl-v)
    END IF
    IF (i /= nx+1 .AND. j /= ny+1) y(k) = y(k) + hxhx*(v-vl) + hyhy*(v-vb)
  END DO
END DO

!     Scale the result.

DO  k = 1, nx*ny
  y(k) = area*y(k)
END DO

RETURN
END SUBROUTINE depths



SUBROUTINE dgl1hs(n, x, s, y, t)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(IN)   :: s(:)
REAL (dp), INTENT(OUT)  :: y(:)
REAL (dp), INTENT(IN)   :: t

!  **********

!  Subroutine dglhs

!  This subroutine computes the product f''(x)*s = y, where
!  f''(x) is the Hessian matrix for the Ginzburg-Landau
!  (1-dimensional) problem.

!  The subroutine statement is

!    subroutine dgl1hs(n, x, s, y, t)

!  where

!    n is an integer variable.
!      On entry n is the number of variables.
!      On exit n is unchanged.

!    x is a REAL (dp) array of dimension n.
!      On entry x specifies the vector x.
!      On exit x is unchanged.

!    y is a REAL (dp) array of dimension n.
!      On entry out need not be specified.
!      On exit y contains f''(x)*s.

!    s is a REAL (dp) array of dimension n.
!      On entry s contains the vector s.
!      On exit s is unchanged.

!    t is a REAL (dp) variable.
!      On entry t is a temperature in (3.73, 7.32).
!      On exit t is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick and Jorge J. More'.

!  **********

REAL (dp), PARAMETER :: two=2.0_dp, three=3.0_dp, four=4.0_dp, six=6.0_dp, &
                        ten=10.0_dp, twelve=12.0_dp, sxteen=16.0_dp

INTEGER   :: i, n1, n2
REAL (dp) :: alphan, alphas, betan, betas, c, dn, ds, ec, em,  &
             fac, gamma, h1, h2, hbar, hcn, hcs, penn, pens, pi, tcn, tcs

!     Initialization.

!     Set electron mass (grams), speed of light (cm/sec), and
!     electronic charge (esu).

em = 9.11D-28
c = 2.99D+10
ec = 4.80D-10

!     Set length of a half-layer of lead and tin (10**3-angstroms),
!     d = ds + dn.

ds = 1.0_dp
dn = 2.2_dp

!     Set critical temperature for lead and tin (Kelvin).

tcs = 7.32_dp
tcn = 3.73_dp

!     Set critical magnetic field for lead and tin at zero
!     temperature (gauss).

hcs = 803.0_dp
hcn = 309.0_dp

!     Set penetration depth for lead and tin at zero temperature
!     (cm).

pens = 3.7D-6
penn = 3.4D-6

!     Compute pi.

pi = four*ATAN(one)

!     Set initial values for temperature dependent constants alphas,
!     alphan (ergs), and betas, betan (ergs-cm**3).

alphas = -two*((ec/c)**2/em)*(hcs**2)*(pens**2)
alphan = -two*((ec/c)**2/em)*(hcn**2)*(penn**2)
betas = sxteen*pi*(((ec/c)**2/em)**2)*(hcs**2)*(pens**4)
betan = sxteen*pi*(((ec/c)**2/em)**2)*(hcn**2)*(penn**4)

alphas = alphas*((one-(t/tcs)**2)/(one+(t/tcs)**2))
alphan = alphan*((one-(t/tcn)**2)/(one+(t/tcn)**2))
betas = betas/((one+(t/tcs)**2)**2)
betan = betan/((one+(t/tcn)**2)**2)

!     Set Planck's constant (erg-sec).

hbar = 1.05459D-27

!     Set temperature dependent constant gamma (erg-cm**2).

gamma = hbar**2/(four*em)

!     Scale temperature dependent constants to the same units.
!     This makes the order parameter dimensionless.

fac = 1.0D6
alphas = alphas*(fac**3)
alphan = alphan*(fac**3)
betas = betas*(fac**6)
betan = betan*(fac**6)
gamma = gamma*(fac**5)

!     Compute the number of subintervals in (-d,-ds), in (-ds,ds),
!     and in (ds,d).

n1 = n/4
n2 = n - 2*n1
h1 = dn/DBLE(n1)
h2 = (two*ds)/DBLE(n2)

y(1:n) = zero

!     Evaluate f''(x)*s over the intervals (-d, -ds), (-ds, ds),
!     and (ds, d).

DO  i = 1, n1
  y(i) = y(i) + h1*(two*alphan/three+  &
      (betan/ten)*(two*x(i+1)**2+six*x(i+1)*x(i)+  &
      twelve*x(i)**2)+two*gamma/h1/h1)*s(i) + h1*(alphan/three+(betan/ten)*  &
      (three*x(i+1)**2+four*x(i+1)*x(i)+three*x(i)**2)- two*gamma/h1/h1)*s(i+1)
  y(i+1) = y(i+1) + h1*(alphan/three+  &
      (betan/ten)*(three*x(i+1)**2+four*x(i+1)*x(i)+  &
      three*x(i)**2)-two*gamma/h1/h1)*s(i) + h1*(two*alphan/three+(betan/ten)*  &
      (twelve*x(i+1)**2+six*x(i+1)*x(i)+two*x(i)**2)+ two*gamma/h1/h1)*s(i+1)
END DO
DO  i = n1 + 1, n1 + n2
  y(i) = y(i) + h2*(two*alphas/three+  &
      (betas/ten)*(two*x(i+1)**2+six*x(i+1)*x(i)+  &
      twelve*x(i)**2)+two*gamma/h2/h2)*s(i) + h2*(alphas/three+(betas/ten)*  &
      (three*x(i+1)**2+four*x(i+1)*x(i)+three*x(i)**2)- two*gamma/h2/h2)*s(i+1)
  y(i+1) = y(i+1) + h2*(alphas/three+  &
      (betas/ten)*(three*x(i+1)**2+four*x(i+1)*x(i)+  &
      three*x(i)**2)-two*gamma/h2/h2)*s(i) + h2*(two*alphas/three+(betas/ten)*  &
      (twelve*x(i+1)**2+six*x(i+1)*x(i)+two*x(i)**2)+ two*gamma/h2/h2)*s(i+1)
END DO
DO  i = n1 + n2 + 1, n - 1
  y(i) = y(i) + h1*(two*alphan/three+  &
      (betan/ten)*(two*x(i+1)**2+six*x(i+1)*x(i)+  &
      twelve*x(i)**2)+two*gamma/h1/h1)*s(i) + h1*(alphan/three+(betan/ten)*  &
      (three*x(i+1)**2+four*x(i+1)*x(i)+three*x(i)**2)- two*gamma/h1/h1)*s(i+1)
  y(i+1) = y(i+1) + h1*(alphan/three+  &
      (betan/ten)*(three*x(i+1)**2+four*x(i+1)*x(i)+  &
      three*x(i)**2)-two*gamma/h1/h1)*s(i) + h1*(two*alphan/three+(betan/ten)*  &
      (twelve*x(i+1)**2+six*x(i+1)*x(i)+two*x(i)**2)+ two*gamma/h1/h1)*s(i+1)
END DO

!     Special case for the right subinterval where x(n+1) = x(1).

y(1) = y(1) + h1*(alphan/three+(betan/ten)*  &
    (three*x(1)**2+four*x(1)*x(n)+three*x(n)**2)-  &
    two*gamma/h1/h1)*s(n) + h1*(two*alphan/three+  &
    (betan/ten)*(twelve*x(1)**2+six*x(1)*x(n)+two*x(n)**2)+ two*gamma/h1/h1)*s(1)
y(n) = y(n) + h1*(two*alphan/three+  &
    (betan/ten)*(two*x(1)**2+six*x(1)*x(n)+twelve*x(n)**2)+  &
    two*gamma/h1/h1)*s(n) + h1*(alphan/three+  &
    (betan/ten)*(three*x(1)**2+four*x(1)*x(n)+three*x(n)**2)-  &
    two*gamma/h1/h1)*s(1)

RETURN
END SUBROUTINE dgl1hs



SUBROUTINE dgl2hs(nx, ny, x, s, y, vornum)

! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)        :: nx
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN OUT)  :: x(:)    ! x(4*nx*ny)
REAL (dp), INTENT(IN OUT)  :: s(:)    ! s(4*nx*ny)
REAL (dp), INTENT(IN OUT)  :: y(:)    ! y(4*nx*ny)
INTEGER, INTENT(IN)        :: vornum

!  **********

!  Subroutine dgl2hs

!  This subroutine computes the product f''(x)*s = y, where f''(x) is the
!  Hessian matrix for the Ginzburg-Landau (2-dimensional) problem evaluted at x.

!  The subroutine statement is

!    subroutine dgl2hs(nx, ny, x, s, y, task, wa1, wa2, vornum)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension 4*nx*ny.
!      On entry x specifies the vector x.
!      On exit x is unchanged.   WRONG!

!    s is a REAL (dp) array of dimension 4*nx*ny.
!      On entry s contains the vector s.
!      On exit s is unchanged.   WRONG!

!    y is a REAL (dp) array of dimension 4*nx*ny.
!      On entry y need not be specified.
!      On exit y contains f''(x)*s.

!    vornum is an integer variable.
!      On entry vornum specifies the number of vortices.
!      On exit vornum is unchanged.

!  Subprograms called

!    MINPACK-2  ...  dgl2co

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Brett M. Averick

!  **********

REAL (dp) :: xx(nx+1,ny+1), sx(1,nx+1,ny+1), yw(nx+1,ny+1), sy(1,nx+1,ny+1), &
             vpotx(nx+1,ny+1), svpotx(1,nx+1,ny+1), vpoty(nx+1,ny+1),  &
             svpoty(1,nx+1,ny+1), yx(1,nx,ny), yy(1,nx,ny), yvpotx(1,nx,ny), &
             yvpoty(1,nx,ny)
INTEGER   :: j, k1, k2, k3, k4

! Initialize arrays

k1 = 0
k2 = nx*ny
k3 = k2 + nx*ny
k4 = k3 + nx*ny
DO j = 1, ny
  xx(1:nx,j)   = x(k1+1:k1+nx)
  sx(1,1:nx,j) = s(k1+1:k1+nx)
  yw(1:nx,j)   = x(k2+1:k2+nx)
  sy(1,1:nx,j) = s(k2+1:k2+nx)
  vpotx(1:nx,j)    = x(k3+1:k3+nx)
  svpotx(1,1:nx,j) = s(k3+1:k3+nx)
  vpoty(1:nx,j)    = x(k4+1:k4+nx)
  svpoty(1,1:nx,j) = s(k4+1:k4+nx)
  yx(1,1:nx,j) = y(k1+1:k1+nx)
  yy(1,1:nx,j) = y(k2+1:k2+nx)
  yvpotx(1,1:nx,j) = y(k3+1:k3+nx)
  yvpoty(1,1:nx,j) = y(k4+1:k4+nx)
  k1 = k1 + nx
  k2 = k2 + nx
  k3 = k3 + nx
  k4 = k4 + nx
END DO

xx(nx+1,:) = zero
xx(:,ny+1) = zero
sx(1,nx+1,:) = zero
sx(1,:,ny+1) = zero
yw(nx+1,:) = zero
yw(:,ny+1) = zero
sy(1,nx+1,:) = zero
sy(1,:,ny+1) = zero
vpotx(nx+1,:) = zero
vpotx(:,ny+1) = zero
svpotx(1,nx+1,:) = zero
svpotx(1,:,ny+1) = zero
vpoty(nx+1,:) = zero
vpoty(:,ny+1) = zero
svpoty(1,nx+1,:) = zero
svpoty(1,:,ny+1) = zero

CALL dgl2co(1, nx, ny, xx, sx, 1, yw, sy, 1, vpotx, svpotx,  &
            1, vpoty, svpoty, 1, yx, 1, yy, 1,  &
            yvpotx, 1, yvpoty, 1, vornum)

! Copy results back into 1-D arrays

k1 = 0
k2 = nx*ny
k3 = k2 + nx*ny
k4 = k3 + nx*ny
DO j = 1, ny
  x(k1+1:k1+nx) = xx(1:nx,j)
  s(k1+1:k1+nx) = sx(1,1:nx,j)
  x(k2+1:k2+nx) = yw(1:nx,j)
  s(k2+1:k2+nx) = sy(1,1:nx,j)
  x(k3+1:k3+nx) = vpotx(1:nx,j)
  s(k3+1:k3+nx) = svpotx(1,1:nx,j)
  x(k4+1:k4+nx) = vpoty(1:nx,j)
  s(k4+1:k4+nx) = svpoty(1,1:nx,j)
  y(k1+1:k1+nx) = yx(1,1:nx,j)
  y(k2+1:k2+nx) = yy(1,1:nx,j)
  y(k3+1:k3+nx) = yvpotx(1,1:nx,j)
  y(k4+1:k4+nx) = yvpoty(1,1:nx,j)
  k1 = k1 + nx
  k2 = k2 + nx
  k3 = k3 + nx
  k4 = k4 + nx
END DO

RETURN
END SUBROUTINE dgl2hs



SUBROUTINE dmsahs(nx, ny, x, s, y, bottom, top, left, right)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)     :: nx
INTEGER, INTENT(IN)     :: ny
REAL (dp), INTENT(IN)   :: x(:)    ! x(nx*ny)
REAL (dp), INTENT(IN)   :: s(:)    ! s(nx*ny)
REAL (dp), INTENT(OUT)  :: y(:)    ! y(nx*ny)
REAL (dp), INTENT(IN)   :: bottom(:)
REAL (dp), INTENT(IN)   :: top(:)
REAL (dp), INTENT(IN)   :: left(:)    ! left(ny+2)
REAL (dp), INTENT(IN)   :: right(:)   ! right(ny+2)

!  **********

!  Subroutine dmsahs

!  This subroutine computes the product f''(x)*s = y, where
!  f''(x) is the Hessian matrix for the minimal surface area
!  problem exaluted at x.

!  The subroutine statement is

!    subroutine dmsahs(nx, ny, x, s, y, bottom, top, left, right)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the vector x.
!      On exit x is unchanged.

!    s is a REAL (dp) array of dimension nx*ny.
!      On entry s contains the vector s.
!      On exit s is unchanged.

!    y is a REAL (dp) array of dimension nx*ny.
!      On entry out need not be specified.
!      On exit y contains f''(x)*s.

!    bottom is a REAL (dp) array of dimension nx + 2.
!      On entry bottom must contain appropriate boundary data.
!      On exit bottom is unchanged.

!    top is a REAL (dp) work array of dimension nx + 2.
!      On entry top must contain appropriate boundary data.
!      On exit top is unchanged.

!    left is a REAL (dp) work array of dimension ny + 2.
!      On entry left must contain appropriate boundary data.
!      On exit left is unchanged.

!    right is a REAL (dp) work array of dimension ny + 2.
!      On entry right must contain appropriate boundary data.
!      On exit right is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp

INTEGER   :: i, j, k
REAL (dp) :: area, dvdx, dvdxhx, dvdy, dvdyhy, dzdx, dzdxhx,  &
             dzdy, dzdyhy, fl, fl3, fu, fu3, hx, hy, tl, tu,  &
             v, vb, vl, vr, vt, z, zb, zl, zr, zt

hx = one/DBLE(nx+1)
hy = one/DBLE(ny+1)
area = p5*hx*hy

DO  k = 1, nx*ny
  y(k) = zero
END DO

!     Computation of f''(x)*s over the lower triangular elements.

DO  j = 0, ny
  DO  i = 0, nx
    k = nx*(j-1) + i
    IF (i /= 0 .AND. j /= 0) THEN
      v = x(k)
      z = s(k)
    ELSE
      IF (j == 0) v = bottom(i+1)
      IF (i == 0) v = left(j+1)
      z = zero
    END IF
    IF (i /= nx .AND. j /= 0) THEN
      vr = x(k+1)
      zr = s(k+1)
    ELSE
      IF (i == nx) vr = right(j+1)
      IF (j == 0) vr = bottom(i+2)
      zr = zero
    END IF
    IF (i /= 0 .AND. j /= ny) THEN
      vt = x(k+nx)
      zt = s(k+nx)
    ELSE
      IF (i == 0) vt = left(j+2)
      IF (j == ny) vt = top(i+1)
      zt = zero
    END IF
    dvdx = (vr-v)/hx
    dvdy = (vt-v)/hy
    dzdx = (zr-z)/hx
    dzdy = (zt-z)/hy
    dvdxhx = dvdx/hx
    dvdyhy = dvdy/hy
    dzdxhx = dzdx/hx
    dzdyhy = dzdy/hy
    tl = one + dvdx**2 + dvdy**2
    fl = SQRT(tl)
    fl3 = fl*tl
    IF (i /= 0 .AND. j /= 0) y(k) = y(k) +  &
        (dvdx*dzdx+dvdy*dzdy)*(dvdxhx+dvdyhy)/fl3 - (dzdxhx+dzdyhy)/fl
    IF (i /= nx .AND. j /= 0) y(k+1) = y(k+1) -  &
        (dvdx*dzdx+dvdy*dzdy)*dvdxhx/fl3 + dzdxhx/fl
    IF (i /= 0 .AND. j /= ny) y(k+nx) = y(k+nx) -  &
        (dvdx*dzdx+dvdy*dzdy)*dvdyhy/fl3 + dzdyhy/fl
  END DO
END DO

!     Computation of f''(x)*s over the upper triangular elements.

DO  j = 1, ny + 1
  DO  i = 1, nx + 1
    k = nx*(j-1) + i
    IF (i /= nx+1 .AND. j /= 1) THEN
      vb = x(k-nx)
      zb = s(k-nx)
    ELSE
      IF (j == 1) vb = bottom(i+1)
      IF (i == nx+1) vb = right(j)
      zb = zero
    END IF
    IF (i /= 1 .AND. j /= ny+1) THEN
      vl = x(k-1)
      zl = s(k-1)
    ELSE
      IF (j == ny+1) vl = top(i)
      IF (i == 1) vl = left(j+1)
      zl = zero
    END IF
    IF (i /= nx+1 .AND. j /= ny+1) THEN
      v = x(k)
      z = s(k)
    ELSE
      IF (i == nx+1) v = right(j+1)
      IF (j == ny+1) v = top(i+1)
      z = zero
    END IF
    dvdx = (v-vl)/hx
    dvdy = (v-vb)/hy
    dzdx = (z-zl)/hx
    dzdy = (z-zb)/hy
    dvdxhx = dvdx/hx
    dvdyhy = dvdy/hy
    dzdxhx = dzdx/hx
    dzdyhy = dzdy/hy
    tu = one + dvdx**2 + dvdy**2
    fu = SQRT(tu)
    fu3 = fu*tu
    IF (i /= nx+1 .AND. j /= ny+1) y(k) = y(k) -  &
        (dvdx*dzdx+dvdy*dzdy)*(dvdxhx+dvdyhy)/fu3 + (dzdxhx+dzdyhy)/fu
    IF (i /= 1 .AND. j /= ny+1) y(k-1) = y(k-1) +  &
        (dvdx*dzdx+dvdy*dzdy)*dvdxhx/fu3 - dzdxhx/fu
    IF (i /= nx+1 .AND. j /= 1) y(k-nx) = y(k-nx) +  &
        (dvdx*dzdx+dvdy*dzdy)*dvdyhy/fu3 - dzdyhy/fu
  END DO
END DO

!     Scale the result by the area.

DO  k = 1, nx*ny
  y(k) = area*y(k)
END DO

RETURN
END SUBROUTINE dmsahs



SUBROUTINE dodchs(nx, ny, x, s, y, lambda)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)     :: nx
INTEGER, INTENT(IN)     :: ny
REAL (dp), INTENT(IN)   :: x(:)    ! x(nx*ny)
REAL (dp), INTENT(IN)   :: s(:)    ! s(nx*ny)
REAL (dp), INTENT(OUT)  :: y(:)    ! y(nx*ny)
REAL (dp), INTENT(IN)   :: lambda

!  **********

!  Subroutine dodchs

!  This subroutine computes the productc f''(x)*s = y, where
!  f''(x) is the Hessian matrix for the optimal design with
!  composites problem evaluated at x.

!  The subroutine statement is

!    subroutine dodchs(nx, ny, x, s, y, lambda)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the vector x.
!      On exit x is unchanged.

!    s is a REAL (dp) array of dimension nx*ny.
!      On entry s contains the vector s.
!      On exit s is unchanged.

!    y is a REAL (dp) array of dimension nx*ny.
!      On entry out need not be specified.
!      On exit y contains f''(x)*s.

!    lambda is a REAL (dp) variable.
!      On entry lambda is the Lagrange multiplier.
!      On exit lambda is unchanged.

!  Subprograms called

!    MINPACK-supplied   ...   dodcps

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp, two=2.0_dp, four=4.0_dp

REAL (dp), PARAMETER :: mu1=1.0_dp, mu2=2.0_dp

INTEGER   :: i, j, k
REAL (dp) :: area, dpsip, dpsipp, dvdx, dvdxhx, dvdy, dvdyhy,  &
             dzdx, dzdxhx, dzdy, dzdyhy, gradv, hx, hy, t1,  &
             t2, v, vb, vl, vr, vt, z, zb, zl, zr, zt

! EXTERNAL dodcps

hx = one/DBLE(nx+1)
hy = one/DBLE(ny+1)
area = p5*hx*hy
t1 = SQRT(two*lambda*mu1/mu2)
t2 = SQRT(two*lambda*mu2/mu1)

DO  k = 1, nx*ny
  y(k) = zero
END DO

!     Computation of f''(x)*s over the lower triangular elements.

DO  j = 0, ny
  DO  i = 0, nx
    k = nx*(j-1) + i
    v = zero
    vr = zero
    vt = zero
    z = zero
    zr = zero
    zt = zero
    IF (i /= 0 .AND. j /= 0) THEN
      v = x(k)
      z = s(k)
    END IF
    IF (i /= nx .AND. j /= 0) THEN
      vr = x(k+1)
      zr = s(k+1)
    END IF
    IF (i /= 0 .AND. j /= ny) THEN
      vt = x(k+nx)
      zt = s(k+nx)
    END IF
    dvdx = (vr-v)/hx
    dvdy = (vt-v)/hy
    dzdx = (zr-z)/hx
    dzdy = (zt-z)/hy
    dvdxhx = dvdx/hx
    dvdyhy = dvdy/hy
    dzdxhx = dzdx/hx
    dzdyhy = dzdy/hy
    gradv = dvdx**2 + dvdy**2
    CALL dodcps(gradv, mu1, mu2, t1, t2, dpsip, 1, lambda)
    CALL dodcps(gradv, mu1, mu2, t1, t2, dpsipp, 2, lambda)
    IF (i /= nx .AND. j /= 0) y(k+1) = y(k+1) +  &
        four*dpsipp*dvdxhx*(dvdx*dzdx+dvdy*dzdy) + two*dpsip*dzdxhx
    IF (i /= 0 .AND. j /= ny) y(k+nx) = y(k+nx) +  &
        four*dpsipp*dvdyhy*(dvdx*dzdx+dvdy*dzdy) + two*dpsip*dzdyhy
    IF (i /= 0 .AND. j /= 0) y(k) = y(k) -  &
        four*dpsipp*(dvdxhx+dvdyhy)*(dvdx*dzdx+dvdy*dzdy) -  &
        two*dpsip*(dzdxhx+dzdyhy)
  END DO
END DO

!     Computation of f''(x)*s over the upper triangular elements.

DO  j = 1, ny + 1
  DO  i = 1, nx + 1
    k = nx*(j-1) + i
    vb = zero
    vl = zero
    v = zero
    zb = zero
    zl = zero
    z = zero
    IF (i /= nx+1 .AND. j /= 1) THEN
      vb = x(k-nx)
      zb = s(k-nx)
    END IF
    IF (i /= 1 .AND. j /= ny+1) THEN
      vl = x(k-1)
      zl = s(k-1)
    END IF
    IF (i /= nx+1 .AND. j /= ny+1) THEN
      v = x(k)
      z = s(k)
    END IF
    dvdx = (v-vl)/hx
    dvdy = (v-vb)/hy
    dzdx = (z-zl)/hx
    dzdy = (z-zb)/hy
    dvdxhx = dvdx/hx
    dvdyhy = dvdy/hy
    dzdxhx = dzdx/hx
    dzdyhy = dzdy/hy
    gradv = dvdx**2 + dvdy**2
    CALL dodcps(gradv, mu1, mu2, t1, t2, dpsip, 1, lambda)
    CALL dodcps(gradv, mu1, mu2, t1, t2, dpsipp, 2, lambda)
    IF (i /= nx+1 .AND. j /= 1) y(k-nx) = y(k-nx) -  &
        four*dpsipp*dvdyhy*(dvdx*dzdx+dvdy*dzdy) - two*dpsip*dzdyhy
    IF (i /= 1 .AND. j /= ny+1) y(k-1) = y(k-1) -  &
        four*dpsipp*dvdxhx*(dvdx*dzdx+dvdy*dzdy) - two*dpsip*dzdxhx
    IF (i /= nx+1 .AND. j /= ny+1) y(k) = y(k) +  &
        four*dpsipp*(dvdxhx+dvdyhy)*(dvdx*dzdx+dvdy*dzdy) +  &
        two*dpsip*(dzdxhx+dzdyhy)
  END DO
END DO

!     Scale the result by the area.

DO  k = 1, nx*ny
  y(k) = area*y(k)
END DO

RETURN
END SUBROUTINE dodchs



SUBROUTINE dpjbhs(nx, ny, s, y, ecc, b)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)     :: nx
INTEGER, INTENT(IN)     :: ny
REAL (dp), INTENT(IN)   :: s(:)    ! s(nx*ny)
REAL (dp), INTENT(OUT)  :: y(:)    ! y(nx*ny)
REAL (dp), INTENT(IN)   :: ecc
REAL (dp), INTENT(IN)   :: b

!  **********

!  Subroutine dpjbhs

!  This subroutine computes the product H*s = y, where H is the
!  Hessian matrix for the pressure distribution in a journal
!  bearing problem.

!  The subroutine statement is

!    subroutine dpjbhs(nx, ny, s, y, ecc, b)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    s is a REAL (dp) array of dimension nx*ny.
!      On entry s contains the vector s.
!      On exit s is unchanged.

!    y is a REAL (dp) array of dimension nx*ny.
!      On entry out need not be specified.
!      On exit y contains H*s.

!    ecc is a REAL (dp) variable
!      On entry ecc is the eccentricity in (0,1).
!      On exit ecc is unchanged

!    b is a REAL (dp) variable
!      On entry b defines the domain as D = (0,2*pi) X (0,2*b).
!      On exit b is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: two=2.0_dp, four=4.0_dp, six=6.0_dp

INTEGER   :: i, j, k
REAL (dp) :: hx, hxhx, hxhy, hy, hyhy, pi, trule, v, vb, vl, vr, vt, xi

! p(xi) = (one + ecc*COS(xi))**3

pi = four*ATAN(one)
hx = two*pi/DBLE(nx+1)
hy = two*b/DBLE(ny+1)
hxhy = hx*hy
hxhx = one/(hx*hx)
hyhy = one/(hy*hy)

DO  k = 1, nx*ny
  y(k) = zero
END DO

!     Computation of H*s over the lower triangular elements.

DO  i = 0, nx
  xi = DBLE(i)*hx
  trule = hxhy*(p(xi, ecc) + p(xi+hx, ecc) + p(xi, ecc))/six
  DO  j = 0, ny
    k = nx*(j-1) + i
    v = zero
    vr = zero
    vt = zero
    IF (i /= 0 .AND. j /= 0) v = s(k)
    IF (i /= nx .AND. j /= 0) THEN
      vr = s(k+1)
      y(k+1) = y(k+1) + trule*hxhx*(vr-v)
    END IF
    IF (i /= 0 .AND. j /= ny) THEN
      vt = s(k+nx)
      y(k+nx) = y(k+nx) + trule*hyhy*(vt-v)
    END IF
    IF (i /= 0 .AND. j /= 0) y(k) = y(k) + trule*(hxhx*(v-vr)+hyhy*(v-vt))
  END DO
END DO

!     Computation of H*s over the upper triangular elements.

DO  i = 1, nx + 1
  xi = DBLE(i)*hx
  trule = hxhy*(p(xi, ecc) + p(xi-hx, ecc) + p(xi, ecc))/six
  DO  j = 1, ny + 1
    k = nx*(j-1) + i
    vb = zero
    vl = zero
    v = zero
    IF (i /= nx+1 .AND. j /= ny+1) v = s(k)
    IF (i /= nx+1 .AND. j /= 1) THEN
      vb = s(k-nx)
      y(k-nx) = y(k-nx) + trule*hyhy*(vb-v)
    END IF
    IF (i /= 1 .AND. j /= ny+1) THEN
      vl = s(k-1)
      y(k-1) = y(k-1) + trule*hxhx*(vl-v)
    END IF
    IF (i /= nx+1 .AND. j /= ny+1)  &
        y(k) = y(k) + trule*(hyhy*(v-vb) + hxhx*(v-vl))
  END DO
END DO

RETURN
END SUBROUTINE dpjbhs



SUBROUTINE dsschs(nx, ny, x, s, y, lambda)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)     :: nx
INTEGER, INTENT(IN)     :: ny
REAL (dp), INTENT(IN)   :: x(:)    ! x(nx*ny)
REAL (dp), INTENT(IN)   :: s(:)    ! s(nx*ny)
REAL (dp), INTENT(OUT)  :: y(:)    ! y(nx*ny)
REAL (dp), INTENT(IN)   :: lambda

!  **********

!  Subroutine dsschs

!  This subroutine computes the product f''(x)*s = y, where
!  f''(x) is the Hessian matrix of the steady state combustion
!  problem exaluated at x.

!  The subroutine statement is

!    subroutine dsschs(nx, ny, x, s, y, lambda)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    x is a REAL (dp) array of dimension nx*ny.
!      On entry x specifies the vector x.
!      On exit x is unchanged.

!    s is a REAL (dp) array of dimension nx*ny.
!      On entry s contains the vector s.
!      On exit s is unchanged.

!    y is a REAL (dp) array of dimension nx*ny.
!      On entry out need not be specified.
!      On exit y contains f''(x)*s.

!    lambda is a REAL (dp) variable.
!      On entry lambda is the multiplier from the Bratu problem.
!      On exit lambda is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: p5=0.5_dp, three=3.0_dp

INTEGER   :: i, j, k
REAL (dp) :: area, hx, hxhx, hy, hyhy, v, vb, vl, vr, vt

hx = one/DBLE(nx+1)
hy = one/DBLE(ny+1)
hxhx = one/(hx*hx)
hyhy = one/(hy*hy)
area = p5*hx*hy

DO  k = 1, nx*ny
  y(k) = zero
END DO

!     Computation of f''(x)*s over the lower triangular elements.

DO  j = 0, ny
  DO  i = 0, nx
    k = nx*(j-1) + i
    v = zero
    vr = zero
    vt = zero
    IF (i /= 0 .AND. j /= 0) v = s(k)
    IF (i /= nx .AND. j /= 0) THEN
      vr = s(k+1)
      y(k+1) = y(k+1) + hxhx*(vr-v) - lambda*EXP(x(k+1))*vr/three
    END IF
    IF (i /= 0 .AND. j /= ny) THEN
      vt = s(k+nx)
      y(k+nx) = y(k+nx) + hyhy*(vt-v) - lambda*EXP(x(k+nx))*vt/three
    END IF
    IF (i /= 0 .AND. j /= 0) y(k) = y(k) + hxhx*(v-vr) +  &
                                    hyhy*(v-vt) - lambda*EXP(x(k))*v/three
  END DO
END DO

!     Computation of f''(x)*s over the upper triangular elements.

DO  j = 1, ny + 1
  DO  i = 1, nx + 1
    k = nx*(j-1) + i
    vb = zero
    vl = zero
    v = zero
    IF (i /= nx+1 .AND. j /= ny+1) v = s(k)
    IF (i /= nx+1 .AND. j /= 1) THEN
      vb = s(k-nx)
      y(k-nx) = y(k-nx) + hyhy*(vb-v) - lambda*EXP(x(k-nx))*vb/three
    END IF
    IF (i /= 1 .AND. j /= ny+1) THEN
      vl = s(k-1)
      y(k-1) = y(k-1) + hxhx*(vl-v) - lambda*EXP(x(k-1))*vl/three
    END IF
    IF (i /= nx+1 .AND. j /= ny+1) y(k) = y(k) + hxhx*(v-vl) + hyhy*(v-vb) - &
                                          lambda*EXP(x(k))*v/three
  END DO
END DO

DO  k = 1, nx*ny
  y(k) = area*y(k)
END DO

RETURN
END SUBROUTINE dsschs



SUBROUTINE deptsp(nx, ny, nnz, indrow, indcol)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:09

INTEGER, INTENT(IN)   :: nx
INTEGER, INTENT(IN)   :: ny
INTEGER, INTENT(OUT)  :: nnz
INTEGER, INTENT(OUT)  :: indrow(:)
INTEGER, INTENT(OUT)  :: indcol(:)

!  **********

!  Subroutine deptsp

!  This subroutine defines the sparsity structure of the Hessian
!  matrix for the elastic plastic torsion problem.

!  The subroutine statement is

!    subroutine deptsp(nx, ny, nnz, indrow, indcol)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    nnz is an integer variable.
!      On entry nnz need not be specified.
!      On exit nnz is set to the number of nonzeros in the
!         lower triangle of the Hessian matrix.

!    indrow is an integer array of dimension at least nnz.
!      On entry indrow need not be specified.
!      On exit indrow contains the row indices of the nonzeros
!        in the lower triangle of the Hessian matrix.

!    indcol is an integer array of dimension at least nnz.
!      On entry indcol need not be specified.
!      On exit indcol contains the column indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

INTEGER :: i, j

!     Compute the sparsity structure.

nnz = 0
DO  j = 1, ny
  DO  i = 1, nx
    nnz = nnz + 1
    indrow(nnz) = (j-1)*nx + i
    indcol(nnz) = (j-1)*nx + i
    IF (i /= nx) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + 1
      indcol(nnz) = (j-1)*nx + i
    END IF
    IF (j /= ny) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + nx
      indcol(nnz) = (j-1)*nx + i
      IF (i /= 1) THEN
        nnz = nnz + 1
        indrow(nnz) = (j-1)*nx + i + nx - 1
        indcol(nnz) = (j-1)*nx + i
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE deptsp



SUBROUTINE dgl1sp(n, nnz, indrow, indcol)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(OUT)  :: nnz
INTEGER, INTENT(OUT)  :: indrow(:)
INTEGER, INTENT(OUT)  :: indcol(:)

!  **********

!  Subroutine dgl1sp

!  This subroutine defines the sparsity structure of the Hessian
!  matrix of the Ginzburg-Landau (1-dimensional) problem.

!  The subroutine statement is

!    subroutine dgl1sp(n,nnz,indrow,indcol)

!  where

!    n is an integer variable.
!      On entry n is the number of grid points.
!      On exit n is unchanged.

!    nnz is an integer variable.
!      On entry nnz need not be specified
!      On exit nnz is set to the number of nonzero index pairs
!         in the sparsity structure. Redundancy is permitted.

!    indrow is an integer array of dimension at least nnz.
!      On entry indrow need not be specified.
!      On exit indrow contains the row indices of the nonzeros
!         in the sparsity structure of the Hessian matrix.

!    indcol is an integer array of dimension at least nnz.
!      On entry indcol need not be specified.
!      On exit indcol contains the column indices of the nonzeros
!         in the sparsity structure of the Hessian matrix.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

INTEGER :: j

!     Compute the sparsity structure.

nnz = 0
DO  j = 1, n
  nnz = nnz + 1
  indrow(nnz) = j
  indcol(nnz) = j
END DO
DO  j = 1, n - 1
  nnz = nnz + 1
  indrow(nnz) = j + 1
  indcol(nnz) = j
END DO
nnz = nnz + 1
indrow(nnz) = n
indcol(nnz) = 1

RETURN
END SUBROUTINE dgl1sp



SUBROUTINE dgl2sp(nx, ny, nnz, indrow, indcol)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

INTEGER, INTENT(IN)   :: nx
INTEGER, INTENT(IN)   :: ny
INTEGER, INTENT(OUT)  :: nnz
INTEGER, INTENT(OUT)  :: indrow(:)
INTEGER, INTENT(OUT)  :: indcol(:)

!  **********

!  Subroutine dgl2sp

!  This subroutine defines the sparsity structure of the Hessian
!  matrix for the Ginzburg-Landau (2-dimensional) problem.

!  The subroutine statement is

!    subroutine dgl2sp(nx, ny, nnz, indrow, indcol)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    nnz is an integer variable.
!      On entry nnz need not be specified.
!      On exit nnz is set to the number of nonzeros in the
!         lower triangle of the Hessian matrix.

!    indrow is an integer array of dimension at least nnz.
!      On entry indrow need not be specified.
!      On exit indrow contains the row indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!    indcol is an integer array of dimension at least nnz.
!      On entry indcol need not be specified.
!      On exit indcol contains the column indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Richard Carter, Paul Plassmann, and Steve Wright.

!  **********

INTEGER :: dof, i, ix, iy, j, k, n

!     Compute the diagonal

n = 4*nx*ny
DO  i = 1, n
  indrow(i) = i
  indcol(i) = i
END DO

!     Order each of the four nodal degrees of freedom sequentially for
!     each vertex (ix,iy), ix = 1, ..., nx , iy = 1, ..., ny.

nnz = n
j = 0
DO  iy = 1, ny
  DO  ix = 1, nx
    DO  dof = 1, 4
      j = j + 1
      DO  k = dof + 1, 4
        nnz = nnz + 1
        indrow(nnz) = 4*(nx*(iy-1)+ix-1) + k
        indcol(nnz) = j
      END DO
      
!              East vertex.
      
      IF (ix /= nx) THEN
        k = 4*(nx*(iy-1)+ix) + 1
        IF (dof <= 3) THEN
          nnz = nnz + 2
          indrow(nnz-1) = k
          indcol(nnz-1) = j
          indrow(nnz) = k + 1
          indcol(nnz) = j
          IF (dof == 3) THEN
            nnz = nnz + 1
            indrow(nnz) = k + 3
            indcol(nnz) = j
          END IF
        ELSE
          nnz = nnz + 1
          indrow(nnz) = k + 3
          indcol(nnz) = j
        END IF
      END IF
      
!              West column of grid.
      
      IF (ix == 1) THEN
        k = 4*(nx*(iy-1)+nx-1) + 1
        IF (dof <= 2) THEN
          nnz = nnz + 3
          indrow(nnz-2) = k
          indcol(nnz-2) = j
          indrow(nnz-1) = k + 1
          indcol(nnz-1) = j
          indrow(nnz) = k + 2
          indcol(nnz) = j
        ELSE IF (dof == 4) THEN
          nnz = nnz + 2
          indrow(nnz-1) = k + 2
          indcol(nnz-1) = j
          indrow(nnz) = k + 3
          indcol(nnz) = j
        END IF
      END IF
      
!              North-west vertex (if not on western column).
      
      IF ((iy /= ny) .AND. (ix /= 1)) THEN
        k = 4*(nx*iy+ix-2) + 1
        IF (dof == 4) THEN
          nnz = nnz + 1
          indrow(nnz) = k + 2
          indcol(nnz) = j
        END IF
      END IF
      
!              North vertex.
      
      IF (iy /= ny) THEN
        k = 4*(nx*iy+ix-1) + 1
        IF (dof <= 2) THEN
          nnz = nnz + 2
          indrow(nnz-1) = k
          indcol(nnz-1) = j
          indrow(nnz) = k + 1
          indcol(nnz) = j
        ELSE IF (dof == 3) THEN
          nnz = nnz + 1
          indrow(nnz) = k + 2
          indcol(nnz) = j
        ELSE
          nnz = nnz + 3
          indrow(nnz-2) = k
          indcol(nnz-2) = j
          indrow(nnz-1) = k + 1
          indcol(nnz-1) = j
          indrow(nnz) = k + 2
          indcol(nnz) = j
        END IF
      END IF
      
!              North-west vertex (if on western column).
      
      IF ((iy /= ny) .AND. (ix == 1)) THEN
        k = 4*(nx*iy+nx-1) + 1
        IF (dof == 4) THEN
          nnz = nnz + 1
          indrow(nnz) = k + 2
          indcol(nnz) = j
        END IF
      END IF
      
!              South-east vertex (if on south-east corner).
      
      IF ((iy == 1) .AND. (ix == nx)) THEN
        k = 4*nx*(ny-1) + 1
        IF (dof == 3) THEN
          nnz = nnz + 1
          indrow(nnz) = k + 3
          indcol(nnz) = j
        END IF
      END IF
      
!              South row of grid.
      
      IF (iy == 1) THEN
        k = 4*(nx*(ny-1)+ix-1) + 1
        IF (dof <= 2) THEN
          nnz = nnz + 3
          indrow(nnz-2) = k
          indcol(nnz-2) = j
          indrow(nnz-1) = k + 1
          indcol(nnz-1) = j
          indrow(nnz) = k + 3
          indcol(nnz) = j
        ELSE IF (dof == 3) THEN
          nnz = nnz + 2
          indrow(nnz-1) = k + 2
          indcol(nnz-1) = j
          indrow(nnz) = k + 3
          indcol(nnz) = j
        END IF
      END IF
      
!              South-east vertex (if on south border, but not on
!              south-east corner).
      
      IF ((iy == 1) .AND. (ix /= nx)) THEN
        k = 4*(nx*(ny-1)+ix) + 1
        IF (dof == 3) THEN
          nnz = nnz + 1
          indrow(nnz) = k + 3
          indcol(nnz) = j
        END IF
      END IF
      
    END DO
  END DO
END DO

!     Reorder by degree of freedom.

DO  k = 1, nnz
  j = MOD(indrow(k),4) - 1
  IF (j == -1) j = 3
  indrow(k) = j*nx*ny + (indrow(k)-1)/4 + 1
  j = MOD(indcol(k),4) - 1
  IF (j == -1) j = 3
  indcol(k) = j*nx*ny + (indcol(k)-1)/4 + 1
END DO

!     Make sure all elements are in the lower half.

DO  k = 1, nnz
  IF (indcol(k) > indrow(k)) THEN
    i = indrow(k)
    indrow(k) = indcol(k)
    indcol(k) = i
  END IF
END DO

RETURN
END SUBROUTINE dgl2sp



SUBROUTINE dmsasp(nx, ny, nnz, indrow, indcol)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

INTEGER, INTENT(IN)   :: nx
INTEGER, INTENT(IN)   :: ny
INTEGER, INTENT(OUT)  :: nnz
INTEGER, INTENT(OUT)  :: indrow(:)
INTEGER, INTENT(OUT)  :: indcol(:)

!  **********

!  Subroutine dmsasp

!  This subroutine defines the sparsity structure of the Hessian
!  matrix for the minimal surface area problem.

!  The subroutine statement is

!    subroutine dmsasp(nx, ny, nnz, indrow, indcol)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    nnz is an integer variable.
!      On entry nnz need not be specified.
!      On exit nnz is set to the number of nonzeros in the
!         lower triangle of the Hessian matrix.

!    indrow is an integer array of dimension at least nnz.
!      On entry indrow need not be specified.
!      On exit indrow contains the row indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!    indcol is an integer array of dimension at least nnz.
!      On entry indcol need not be specified.
!      On exit indcol contains the column indices of the nonzeros
!        in the lower triangle of the Hessian matrix.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

INTEGER :: i, j

!     Compute the sparsity structure.

nnz = 0
DO  j = 1, ny
  DO  i = 1, nx
    nnz = nnz + 1
    indrow(nnz) = (j-1)*nx + i
    indcol(nnz) = (j-1)*nx + i
    IF (i /= nx) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + 1
      indcol(nnz) = (j-1)*nx + i
    END IF
    IF (j /= ny) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + nx
      indcol(nnz) = (j-1)*nx + i
      IF (i /= 1) THEN
        nnz = nnz + 1
        indrow(nnz) = (j-1)*nx + i + nx - 1
        indcol(nnz) = (j-1)*nx + i
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE dmsasp



SUBROUTINE dodcsp(nx, ny, nnz, indrow, indcol)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

INTEGER, INTENT(IN)   :: nx
INTEGER, INTENT(IN)   :: ny
INTEGER, INTENT(OUT)  :: nnz
INTEGER, INTENT(OUT)  :: indrow(:)
INTEGER, INTENT(OUT)  :: indcol(:)

!  **********

!  Subroutine dodcsp

!  This subroutine defines the sparsity structure of the Hessian
!  matrix for the optimal design with composites problem.

!  The subroutine statement is

!    subroutine dodcsp(nx, ny, nnz, indrow, indcol)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    nnz is an integer variable.
!      On entry nnz need not be specified.
!      On exit nnz is set to the number of nonzeros in the
!         lower triangle of the Hessian matrix.

!    indrow is an integer array of dimension at least nnz.
!      On entry indrow need not be specified.
!      On exit indrow contains the row indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!    indcol is an integer array of dimension at least nnz.
!      On entry indcol need not be specified.
!      On exit indcol contains the column indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

INTEGER :: i, j

!     Compute the sparsity structure.

nnz = 0
DO  j = 1, ny
  DO  i = 1, nx
    nnz = nnz + 1
    indrow(nnz) = (j-1)*nx + i
    indcol(nnz) = (j-1)*nx + i
    IF (i /= nx) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + 1
      indcol(nnz) = (j-1)*nx + i
    END IF
    IF (j /= ny) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + nx
      indcol(nnz) = (j-1)*nx + i
      IF (i /= 1) THEN
        nnz = nnz + 1
        indrow(nnz) = (j-1)*nx + i + nx - 1
        indcol(nnz) = (j-1)*nx + i
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE dodcsp



SUBROUTINE dpjbsp(nx, ny, nnz, indrow, indcol)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

INTEGER, INTENT(IN)   :: nx
INTEGER, INTENT(IN)   :: ny
INTEGER, INTENT(OUT)  :: nnz
INTEGER, INTENT(OUT)  :: indrow(:)
INTEGER, INTENT(OUT)  :: indcol(:)

!  **********

!  Subroutine dpjbsp

!  This subroutine defines the sparsity structure of the Hessian
!  matrix for the pressure distribution in a journal bearing problem.

!  The subroutine statement is

!    subroutine dpjbsp(nx, ny, nnz, indrow, indcol)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    nnz is an integer variable.
!      On entry nnz need not be specified.
!      On exit nnz is set to the number of nonzero in the
!         lower triangle of the Hessian matrix.

!    indrow is an integer array of dimension at least nnz.
!      On entry indrow need not be specified.
!      On exit indrow contains the row indices of the nonzeros
!         in the lower triangle of the Hessian matrix

!    indcol is an integer array of dimension at least nnz.
!      On entry indcol need not be specified.
!      On exit indcol contains the column indices of the nonzeros
!        in the lower triangle of the Hessian matrix

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

INTEGER :: i, j

!     Compute the sparsity structure.

nnz = 0
DO  j = 1, ny
  DO  i = 1, nx
    nnz = nnz + 1
    indrow(nnz) = (j-1)*nx + i
    indcol(nnz) = (j-1)*nx + i
    IF (i /= nx) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + 1
      indcol(nnz) = (j-1)*nx + i
    END IF
    IF (j /= ny) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + nx
      indcol(nnz) = (j-1)*nx + i
      IF (i /= 1) THEN
        nnz = nnz + 1
        indrow(nnz) = (j-1)*nx + i + nx - 1
        indcol(nnz) = (j-1)*nx + i
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE dpjbsp



SUBROUTINE dsscsp(nx, ny, nnz, indrow, indcol)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

INTEGER, INTENT(IN)   :: nx
INTEGER, INTENT(IN)   :: ny
INTEGER, INTENT(OUT)  :: nnz
INTEGER, INTENT(OUT)  :: indrow(:)
INTEGER, INTENT(OUT)  :: indcol(:)

!  **********

!  Subroutine dsscsp

!  This subroutine defines the sparsity structure of the Hessian
!  matrix for the steady state combustion problem.

!  The subroutine statement is

!    subroutine dsscsp(nx, ny, nnz, indrow, indcol)

!  where

!    nx is an integer variable.
!      On entry nx is the number of grid points in the first
!         coordinate direction.
!      On exit nx is unchanged.

!    ny is an integer variable.
!      On entry ny is the number of grid points in the second
!         coordinate direction.
!      On exit ny is unchanged.

!    nnz is an integer variable.
!      On entry nnz need not be specified.
!      On exit nnz is set to the number of nonzeros in the
!         lower triangle of the Hessian matrix.

!    indrow is an integer array of dimension at least nnz.
!      On entry indrow need not be specified.
!      On exit indrow contains the row indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!    indcol is an integer array of dimension at least nnz.
!      On entry indcol need not be specified.
!      On exit indcol contains the column indices of the nonzeros
!         in the lower triangle of the Hessian matrix.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

INTEGER :: i, j

!     Compute the sparsity structure.

nnz = 0
DO  j = 1, ny
  DO  i = 1, nx
    nnz = nnz + 1
    indrow(nnz) = (j-1)*nx + i
    indcol(nnz) = (j-1)*nx + i
    IF (i /= nx) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + 1
      indcol(nnz) = (j-1)*nx + i
    END IF
    IF (j /= ny) THEN
      nnz = nnz + 1
      indrow(nnz) = (j-1)*nx + i + nx
      indcol(nnz) = (j-1)*nx + i
      IF (i /= 1) THEN
        nnz = nnz + 1
        indrow(nnz) = (j-1)*nx + i + nx - 1
        indcol(nnz) = (j-1)*nx + i
      END IF
    END IF
  END DO
END DO

RETURN
END SUBROUTINE dsscsp



SUBROUTINE dodcps(t, mu1, mu2, t1, t2, result, option, lambda)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:10

REAL (dp), INTENT(IN)   :: t
REAL (dp), INTENT(IN)   :: mu1
REAL (dp), INTENT(IN)   :: mu2
REAL (dp), INTENT(IN)   :: t1
REAL (dp), INTENT(IN)   :: t2
REAL (dp), INTENT(OUT)  :: result
INTEGER, INTENT(IN)     :: option
REAL (dp), INTENT(IN)   :: lambda

!  **********

!  This subroutine computes the function psi(t) and the scaled
!  functions psi'(t)/t and psi''(t)/t for the optimal design
!  with composite materials problem.

!  The subroutine statement is

!    subroutine dodcps(t, mu1, mu2, t1, t2, result, option, lambda)

!  where

!    t is a REAL (dp) variable.
!      On entry t is the variable t
!      On exit t is unchanged

!    mu1 is a REAL (dp) variable.
!      On entry mu1 is the reciprocal shear modulus of material 1.
!      On exit mu1 is unchanged.

!    mu2 is a REAL (dp) variable.
!      On entry mu2 is the reciprocal shear modulus of material 2.
!      On exit mu2 is unchanged.

!    t1 is a REAL (dp) variable.
!      On entry t1 is the first breakpoint.
!      On exit t1 is unchanged.

!    t2 is a REAL (dp) variable.
!      On entry t2 is the second breakpoint.
!      On exit t2 is unchanged.

!    result is a REAL (dp) variable.
!      On entry result need not be specified.
!      On exit result is set according to task.

!    option is an integer variable.
!      On entry option specifies the action of the subroutine:

!         if option = 0 then evaluate the function psi(t).
!         if option = 1 then evaluate the scaled function psi'(t)/t.
!         if option = 2 then evaluate the scaled function psi''(t)/t.

!     On option task is unchanged.

!    lambda is a REAL (dp) variable
!      On entry lambda is the Lagrange multiplier.
!      On exit lambda is unchanged.

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory and University of Minnesota.
!  Brett M. Averick.

!  **********

REAL (dp), PARAMETER :: p25=0.25_dp, p5=0.5_dp

REAL (dp) :: sqrtt

sqrtt = SQRT(t)

IF (option == 0) THEN
  IF (sqrtt <= t1) THEN
    result = p5*mu2*t
  ELSE IF (sqrtt > t1 .AND. sqrtt < t2) THEN
    result = mu2*t1*sqrtt - lambda*mu1
  ELSE IF (sqrtt >= t2) THEN
    result = p5*mu1*t + lambda*(mu2-mu1)
  END IF
ELSE IF (option == 1) THEN
  IF (sqrtt <= t1) THEN
    result = p5*mu2
  ELSE IF (sqrtt > t1 .AND. sqrtt < t2) THEN
    result = p5*mu2*t1/sqrtt
  ELSE IF (sqrtt >= t2) THEN
    result = p5*mu1
  END IF
ELSE IF (option == 2) THEN
  IF (sqrtt <= t1) THEN
    result = zero
  ELSE IF (sqrtt > t1 .AND. sqrtt < t2) THEN
    result = -p25*mu2*t1/(sqrtt*t)
  ELSE IF (sqrtt >= t2) THEN
    result = zero
  END IF
END IF

RETURN
END SUBROUTINE dodcps



SUBROUTINE dgl2co(gp, nx, ny, x, sx, ldsx, y, sy, ldsy, vpotx, svpotx,  &
                  ldsvpx, vpoty, svpoty, ldsvpy, yx, ldyx, yy, ldyy,  &
                  yvpotx, ldyvpx, yvpoty, ldyvpy, vornum)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-06-29  Time: 15:49:11

INTEGER, INTENT(IN)        :: gp
INTEGER, INTENT(IN)        :: nx
INTEGER, INTENT(IN)        :: ny
REAL (dp), INTENT(IN OUT)  :: x(:,:)         ! x(nx+1,ny+1)
REAL (dp), INTENT(IN OUT)  :: sx(:,:,:)      ! sx(ldsx,nx+1,ny+1)
INTEGER, INTENT(IN)        :: ldsx
REAL (dp), INTENT(IN OUT)  :: y(:,:)         ! y(nx+1,ny+1)
REAL (dp), INTENT(IN OUT)  :: sy(:,:,:)      ! sy(ldsy,nx+1,ny+1)
INTEGER, INTENT(IN)        :: ldsy
REAL (dp), INTENT(IN OUT)  :: vpotx(:,:)     ! vpotx(nx+1,ny+1)
REAL (dp), INTENT(IN OUT)  :: svpotx(:,:,:)  ! svpotx(ldsvpx,nx+1,ny+1)
INTEGER, INTENT(IN)        :: ldsvpx
REAL (dp), INTENT(IN OUT)  :: vpoty(:,:)     ! vpoty(nx+1,ny+1)
REAL (dp), INTENT(IN OUT)  :: svpoty(:,:,:)  ! svpoty(ldsvpy,nx+1,ny+1)
INTEGER, INTENT(IN)        :: ldsvpy
REAL (dp), INTENT(OUT)     :: yx(:,:,:)      ! yx(ldyx,nx,ny)
INTEGER, INTENT(IN)        :: ldyx
REAL (dp), INTENT(OUT)     :: yy(:,:,:)      ! yy(ldyy,nx,ny)
INTEGER, INTENT(IN)        :: ldyy
REAL (dp), INTENT(OUT)     :: yvpotx(:,:,:)  ! yvpotx(ldyvpx,nx,ny)
INTEGER, INTENT(IN)        :: ldyvpx
REAL (dp), INTENT(OUT)     :: yvpoty(:,:,:)  ! yvpoty(ldyvpy,nx,ny)
INTEGER, INTENT(IN)        :: ldyvpy
INTEGER, INTENT(IN)        :: vornum

!  **********

!  Subroutine dgl2co

!  This subroutine computes the product f''(x)*s = y, where f''(x) is the
!  Hessian matrix for the Ginzburg-Landau (2-dimensional) problem evaluted at x.

!  This subroutine was obtained by running dgl2fg.f through
!  ADIFOR, then through nag_polish, then edited for readibility.

!  The subroutine statement is

!    subroutine dgl2co(gp, nx, ny, x, sx, ldsx, y, sy, ldsy, vpotx, svpotx,
!                      ldsvpx, vpoty, svpoty, ldsvpy, yx, ldyx, yy, ldyy,
!                      yvpotx, ldyvpx, yvpoty, ldyvpy, vornum)

!  MINPACK-2 Project. March 1999.
!  Argonne National Laboratory.
!  Brett M. Averick.

!  **********

INTEGER, PARAMETER   :: gpmax=20
REAL (dp), PARAMETER :: two=2.0_dp, three=3.0_dp, four=4.0_dp, five=5.0_dp

INTEGER   :: i, j
REAL (dp) :: arg, cfac, d1, d1bar, ffbar, d2, d2bar, d3, d3bar, d4, d4bar, &
             d5, d5bar, d6, d7, d7bar, d8, d8bar, d9, d9bar, d10, d10bar,  &
             d11, d12, d13bar, d14, d14bar, d15, d17, d17bar, d18, d18bar, &
             d19, d19bar, d20, d20bar, d21, d23, d23bar, d24, d25bar, d26, &
             d27, d29, ffield, fkinx1, fkinx2, fkiny1, fkiny2, hx, hy, pi, &
             sfac, sqn, tkappa
REAL (dp) :: gffld(gpmax), gfknx1(gpmax), gfknx2(gpmax),  &
             gfkny1(gpmax), gfkny2(gpmax)

!     Initialize.

tkappa = five
hx = SQRT(vornum/two)*three/DBLE(nx)
hy = SQRT(vornum/two)*three*SQRT(three)/DBLE(ny)
sqn = DBLE(nx*ny)
pi = four*ATAN(one)

!     Enforce vortex constraint and boundary conditions.

!     Right face for order parameter and vector potential.

DO  j = 1, ny + 1
  arg = two*pi*vornum*(DBLE(j)-one)/DBLE(ny)
  d3bar = -SIN(arg)
  d1bar = COS(arg)
  sx(1:gp,nx+1,j) = d1bar*sx(1:gp,1,j) + d3bar*sy(1:gp,1,j)
  x(nx+1,j) = x(1,j)*COS(arg) - y(1,j)*SIN(arg)
  d3bar = COS(arg)
  d1bar = SIN(arg)
  sy(1:gp,nx+1,j) = d1bar*sx(1:gp,1,j) + d3bar*sy(1:gp,1,j)
  y(nx+1,j) = x(1,j)*SIN(arg) + y(1,j)*COS(arg)
  svpotx(1:gp,nx+1,j) = svpotx(1:gp,1,j)
  vpotx(nx+1,j) = vpotx(1,j)
  svpoty(1:gp,nx+1,j) = svpoty(1:gp,1,j)
  vpoty(nx+1,j) = vpoty(1,j) + two*pi*vornum/(DBLE(ny)*hy)
END DO

!     Top face for order parameter and vector potential.

DO  i = 1, nx + 1
  sx(1:gp,i,ny+1) = sx(1:gp,i,1)
  x(i,ny+1) = x(i,1)
  sy(1:gp,i,ny+1) = sy(1:gp,i,1)
  y(i,ny+1) = y(i,1)
  svpotx(1:gp,i,ny+1) = svpotx(1:gp,i,1)
  vpotx(i,ny+1) = vpotx(i,1)
  svpoty(1:gp,i,ny+1) = svpoty(1:gp,i,1)
  vpoty(i,ny+1) = vpoty(i,1)
END DO

DO  j = 1, ny
  DO  i = 1, nx
    d1 = x(i,j)
    d2 = x(i,j)
    d3 = d2**2
    d5 = y(i,j)
    d6 = d5**2
    d7 = -one + d3 + d6
    IF (d5 /= zero) THEN
      d5bar = 2*(d1*(d6/d5))
    ELSE
      d5bar = zero
    END IF
    IF (d2 /= zero) THEN
      d2bar = 2*(d1*(d3/d2))
    ELSE
      d2bar = zero
    END IF
    yx(1:gp,i,j) = d7*sx(1:gp,i,j) + d2bar*sx(1:gp,i,j) + d5bar*sy(1:gp,i,j)
    d1bar = (one/sqn)*two
    yx(1:gp,i,j) = d1bar*yx(1:gp,i,j)
    d1 = y(i,j)
    d2 = x(i,j)
    d3 = d2**2
    d5 = y(i,j)
    d6 = d5**2
    d7 = -one + d3 + d6
    IF (d5 /= zero) THEN
      d5bar = 2*(d1*(d6/d5))
    ELSE
      d5bar = zero
    END IF
    IF (d2 /= zero) THEN
      d2bar = 2*(d1*(d3/d2))
    ELSE
      d2bar = zero
    END IF
    yy(1:gp,i,j) = d7*sy(1:gp,i,j) + d2bar*sx(1:gp,i,j) + d5bar*sy(1:gp,i,j)
    d1bar = (one/sqn)*two
    yy(1:gp,i,j) = d1bar*yy(1:gp,i,j)
    yvpotx(1:gp,i,j) = zero
    yvpoty(1:gp,i,j) = zero
  END DO
END DO

!     Kinetic energy part, interior points

DO  i = 2, nx
  DO  j = 2, ny
    d2 = x(i,j)
    d4 = hx*vpotx(i,j)
    d5 = COS(d4)
    d8 = y(i,j)
    d10 = hx*vpotx(i,j)
    d11 = SIN(d10)
    d13bar = (two/(hx*hx*sqn))
    d8bar = d13bar*d11
    d9bar = COS(d10)*(d13bar*d8)*hx
    d2bar = -d13bar*d5
    d3bar = (-SIN(d4)*(-d13bar*d2))*hx
    gfknx1(1:gp) = d13bar*sx(1:gp,i+1,j) + d2bar*sx(1:gp,i,j) +  &
                   d3bar*svpotx(1:gp,i,j) + d8bar*sy(1:gp,i,j) + &
                   d9bar*svpotx(1:gp,i,j)
    fkinx1 = (two/(hx*hx*sqn))*(x(i+1,j)-d2*d5+d8*d11)
    d2 = y(i,j)
    d4 = hx*vpotx(i,j)
    d5 = COS(d4)
    d8 = x(i,j)
    d10 = hx*vpotx(i,j)
    d11 = SIN(d10)
    d13bar = (two/(hx*hx*sqn))
    d8bar = -d13bar*d11
    d9bar = COS(d10)*(-d13bar*d8)*hx
    d2bar = -d13bar*d5
    d3bar = (-SIN(d4)*(-d13bar*d2))*hx
    gfknx2(1:gp) = d13bar*sy(1:gp,i+1,j) + d2bar*sy(1:gp,i,j) +  &
                   d3bar*svpotx(1:gp,i,j) + d8bar*sx(1:gp,i,j) + &
                   d9bar*svpotx(1:gp,i,j)
    fkinx2 = (two/(hx*hx*sqn))*(y(i+1,j)-d2*d5-d8*d11)
    d2 = x(i,j)
    d4 = hy*vpoty(i,j)
    d5 = COS(d4)
    d8 = y(i,j)
    d10 = hy*vpoty(i,j)
    d11 = SIN(d10)
    d13bar = (two/(hy*hy*sqn))
    d8bar = d13bar*d11
    d9bar = COS(d10)*(d13bar*d8)*hy
    d2bar = -d13bar*d5
    d3bar = (-SIN(d4)*(-d13bar*d2))*hy
    gfkny1(1:gp) = d13bar*sx(1:gp,i,j+1) + d2bar*sx(1:gp,i,j) +  &
                   d3bar*svpoty(1:gp,i,j) + d8bar*sy(1:gp,i,j) + &
                   d9bar*svpoty(1:gp,i,j)
    fkiny1 = (two/(hy*hy*sqn))*(x(i,j+1)-d2*d5+d8*d11)
    d2 = y(i,j)
    d4 = hy*vpoty(i,j)
    d5 = COS(d4)
    d8 = x(i,j)
    d10 = hy*vpoty(i,j)
    d11 = SIN(d10)
    d13bar = (two/(hy*hy*sqn))
    d8bar = -d13bar*d11
    d9bar = COS(d10)*(-d13bar*d8)*hy
    d2bar = -d13bar*d5
    d3bar = (-SIN(d4)*(-d13bar*d2))*hy
    gfkny2(1:gp) = d13bar*sy(1:gp,i,j+1) + d2bar*sy(1:gp,i,j) +  &
                   d3bar*svpoty(1:gp,i,j) + d8bar*sx(1:gp,i,j) + &
                   d9bar*svpoty(1:gp,i,j)
    fkiny2 = (two/(hy*hy*sqn))*(y(i,j+1)-d2*d5-d8*d11)
    d7bar = (one/hx)
    d3bar = (one/hy)
    gffld(1:gp) = d3bar*svpotx(1:gp,i,j) + (-d3bar*svpotx(1:gp,i,j+1)) +  &
                  d7bar*svpoty(1:gp,i+1,j) + (-d7bar*svpoty(1:gp,i,j))
    ffield = (vpotx(i,j)-vpotx(i,j+1))/hy + (vpoty(i+1,j)-vpoty(i,j))/hx
    ffbar = (two*(tkappa**2)/sqn)
    gffld(1:gp) = ffbar*gffld(1:gp)
    ffield = (two*(tkappa**2)/sqn)*ffield
    d3 = hx*vpotx(i,j)
    d5 = -COS(d3)
    d9 = hx*vpotx(i,j)
    d11 = -SIN(d9)
    d15 = hy*vpoty(i,j)
    d17 = -COS(d15)
    d21 = hy*vpoty(i,j)
    d23 = -SIN(d21)
    d20bar = COS(d21)*(-fkiny2)*hy
    d14bar = (-SIN(d15)*(-fkiny1))*hy
    d8bar = COS(d9)*(-fkinx2)*hx
    d2bar = (-SIN(d3)*(-fkinx1))*hx
    yx(1:gp,i,j) = d5*gfknx1(1:gp) + d11*gfknx2(1:gp) +  &
                   d17*gfkny1(1:gp) + d23*gfkny2(1:gp) + yx(1:gp,i,j) + &
                   d2bar*svpotx(1:gp,i,j) + d8bar*svpotx(1:gp,i,j) + &
                   d14bar*svpoty(1:gp,i,j) + d20bar*svpoty(1:gp,i,j)
    d3 = hx*vpotx(i,j)
    d4 = SIN(d3)
    d8 = hx*vpotx(i,j)
    d10 = -COS(d8)
    d14 = hy*vpoty(i,j)
    d15 = SIN(d14)
    d19 = hy*vpoty(i,j)
    d21 = -COS(d19)
    d18bar = (-SIN(d19)*(-fkiny2))*hy
    d13bar = COS(d14)*fkiny1*hy
    d7bar = (-SIN(d8)*(-fkinx2))*hx
    d2bar = COS(d3)*fkinx1*hx
    yy(1:gp,i,j) = d4*gfknx1(1:gp) + d10*gfknx2(1:gp) +  &
          d15*gfkny1(1:gp) + d21*gfkny2(1:gp) + yy(1:gp,i,j) + &
          d2bar*svpotx(1:gp,i,j) + d7bar*svpotx(1:gp,i,j) + &
          d13bar*svpoty(1:gp,i,j) + d18bar*svpoty(1:gp,i,j)
    d3 = hx*x(i,j)
    d5 = hx*vpotx(i,j)
    d6 = SIN(d5)
    d9 = hx*y(i,j)
    d11 = hx*vpotx(i,j)
    d12 = COS(d11)
    d14 = d3*d6 + d9*d12
    d18 = hx*y(i,j)
    d20 = hx*vpotx(i,j)
    d21 = SIN(d20)
    d24 = hx*x(i,j)
    d26 = hx*vpotx(i,j)
    d27 = COS(d26)
    d29 = d18*d21 - d24*d27
    ffbar = (one/hy)
    d25bar = (-SIN(d26)*(-fkinx2*d24))*hx
    d23bar = -fkinx2*d27*hx
    d19bar = COS(d20)*(fkinx2*d18)*hx
    d17bar = fkinx2*d21*hx
    d10bar = (-SIN(d11)*(fkinx1*d9))*hx
    d8bar = fkinx1*d12*hx
    d4bar = COS(d5)*(fkinx1*d3)*hx
    d2bar = fkinx1*d6*hx
    yvpotx(1:gp,i,j) = d14*gfknx1(1:gp) + d29*gfknx2(1:gp) +  &
          ffbar*gffld(1:gp) + yvpotx(1:gp,i,j) + d2bar*sx(1:gp,i,j) +  &
          d4bar*svpotx(1:gp,i,j) + d8bar*sy(1:gp,i,j) +  &
          d10bar*svpotx(1:gp,i,j) + d17bar*sy(1:gp,i,j) +  &
          d19bar*svpotx(1:gp,i,j) + d23bar*sx(1:gp,i,j) +  &
          d25bar*svpotx(1:gp,i,j)
    d3 = hy*x(i,j)
    d5 = hy*vpoty(i,j)
    d6 = SIN(d5)
    d9 = hy*y(i,j)
    d11 = hy*vpoty(i,j)
    d12 = COS(d11)
    d14 = d3*d6 + d9*d12
    d18 = hy*y(i,j)
    d20 = hy*vpoty(i,j)
    d21 = SIN(d20)
    d24 = hy*x(i,j)
    d26 = hy*vpoty(i,j)
    d27 = COS(d26)
    d29 = d18*d21 - d24*d27
    ffbar = -(one/hx)
    d25bar = (-SIN(d26)*(-fkiny2*d24))*hy
    d23bar = -fkiny2*d27*hy
    d19bar = COS(d20)*(fkiny2*d18)*hy
    d17bar = fkiny2*d21*hy
    d10bar = (-SIN(d11)*(fkiny1*d9))*hy
    d8bar = fkiny1*d12*hy
    d4bar = COS(d5)*(fkiny1*d3)*hy
    d2bar = fkiny1*d6*hy
    yvpoty(1:gp,i,j) = d14*gfkny1(1:gp) + d29*gfkny2(1:gp) +  &
          ffbar*gffld(1:gp) + yvpoty(1:gp,i,j) + d2bar*sx(1:gp,i,j) +  &
          d4bar*svpoty(1:gp,i,j) + d8bar*sy(1:gp,i,j) +  &
          d10bar*svpoty(1:gp,i,j) + d17bar*sy(1:gp,i,j) +  &
          d19bar*svpoty(1:gp,i,j) + d23bar*sx(1:gp,i,j) +  &
          d25bar*svpoty(1:gp,i,j)
    d2 = x(i-1,j)
    d4 = hx*vpotx(i-1,j)
    d5 = COS(d4)
    d8 = y(i-1,j)
    d10 = hx*vpotx(i-1,j)
    d11 = SIN(d10)
    d13bar = (two/(hx*hx*sqn))
    d8bar = d13bar*d11
    d9bar = COS(d10)*(d13bar*d8)*hx
    d2bar = -d13bar*d5
    d3bar = (-SIN(d4)*(-d13bar*d2))*hx
    gfknx1(1:gp) = d13bar*sx(1:gp,i,j) + d2bar*sx(1:gp,i-1,j) +  &
                   d3bar*svpotx(1:gp,i-1,j) + d8bar*sy(1:gp,i-1,j) + &
                   d9bar*svpotx(1:gp,i-1,j)
    fkinx1 = (two/(hx*hx*sqn))*(x(i,j)-d2*d5+d8*d11)
    d2 = y(i-1,j)
    d4 = hx*vpotx(i-1,j)
    d5 = COS(d4)
    d8 = x(i-1,j)
    d10 = hx*vpotx(i-1,j)
    d11 = SIN(d10)
    d13bar = (two/(hx*hx*sqn))
    d8bar = -d13bar*d11
    d9bar = COS(d10)*(-d13bar*d8)*hx
    d2bar = -d13bar*d5
    d3bar = (-SIN(d4)*(-d13bar*d2))*hx
    gfknx2(1:gp) = d13bar*sy(1:gp,i,j) + d2bar*sy(1:gp,i-1,j) +  &
                   d3bar*svpotx(1:gp,i-1,j) + d8bar*sx(1:gp,i-1,j) + &
                   d9bar*svpotx(1:gp,i-1,j)
    fkinx2 = (two/(hx*hx*sqn))*(y(i,j)-d2*d5-d8*d11)
    d2 = x(i,j-1)
    d4 = hy*vpoty(i,j-1)
    d5 = COS(d4)
    d8 = y(i,j-1)
    d10 = hy*vpoty(i,j-1)
    d11 = SIN(d10)
    d13bar = (two/(hy*hy*sqn))
    d8bar = d13bar*d11
    d9bar = COS(d10)*(d13bar*d8)*hy
    d2bar = -d13bar*d5
    d3bar = (-SIN(d4)*(-d13bar*d2))*hy
    gfkny1(1:gp) = d13bar*sx(1:gp,i,j) + d2bar*sx(1:gp,i,j-1) +  &
                   d3bar*svpoty(1:gp,i,j-1) + d8bar*sy(1:gp,i,j-1) + &
                   d9bar*svpoty(1:gp,i,j-1)
    fkiny1 = (two/(hy*hy*sqn))*(x(i,j)-d2*d5+d8*d11)
    d2 = y(i,j-1)
    d4 = hy*vpoty(i,j-1)
    d5 = COS(d4)
    d8 = x(i,j-1)
    d10 = hy*vpoty(i,j-1)
    d11 = SIN(d10)
    d13bar = (two/(hy*hy*sqn))
    d8bar = -d13bar*d11
    d9bar = COS(d10)*(-d13bar*d8)*hy
    d2bar = -d13bar*d5
    d3bar = (-SIN(d4)*(-d13bar*d2))*hy
    gfkny2(1:gp) = d13bar*sy(1:gp,i,j) + d2bar*sy(1:gp,i,j-1) +  &
                   d3bar*svpoty(1:gp,i,j-1) + d8bar*sx(1:gp,i,j-1) + &
                   d9bar*svpoty(1:gp,i,j-1)
    fkiny2 = (two/(hy*hy*sqn))*(y(i,j)-d2*d5-d8*d11)
    yx(1:gp,i,j) = gfknx1(1:gp) + gfkny1(1:gp) + yx(1:gp,i,j)
    yy(1:gp,i,j) = gfknx2(1:gp) + gfkny2(1:gp) + yy(1:gp,i,j)
    d7bar = (one/hx)
    d3bar = (one/hy)
    gffld(1:gp) = d3bar*svpotx(1:gp,i,j-1) + (-d3bar*svpotx(1:gp,i,j)) +  &
                  d7bar*svpoty(1:gp,i+1,j-1) + (-d7bar*svpoty(1:gp,i,j-1))
    ffield = (vpotx(i,j-1)-vpotx(i,j))/hy + (vpoty(i+1,j-1)-vpoty(i,j-1))/hx
    ffbar = (two*(tkappa**2)/sqn)
    gffld(1:gp) = ffbar*gffld(1:gp)
    ffield = (two*(tkappa**2)/sqn)*ffield
    ffbar = -(one/hy)
    yvpotx(1:gp,i,j) = ffbar*gffld(1:gp) + yvpotx(1:gp,i,j)
    d7bar = (one/hx)
    d3bar = (one/hy)
    gffld(1:gp) = d3bar*svpotx(1:gp,i-1,j) + (-d3bar*svpotx(1:gp,i-1,j+1)) +  &
                  d7bar*svpoty(1:gp,i,j) + (-d7bar*svpoty(1:gp,i-1,j))
    ffield = (vpotx(i-1,j)-vpotx(i-1,j+1))/hy + (vpoty(i,j)-vpoty(i-1,j))/hx
    ffbar = (two*(tkappa**2)/sqn)
    gffld(1:gp) = ffbar*gffld(1:gp)
    ffield = (two*(tkappa**2)/sqn)*ffield
    ffbar = (one/hx)
    yvpoty(1:gp,i,j) = ffbar*gffld(1:gp) + yvpoty(1:gp,i,j)
  END DO
END DO

!     Kinetic energy part, boundary points.

!     Bottom j = 1

DO  i = 2, nx
  d2 = x(i,1)
  d4 = hx*vpotx(i,1)
  d5 = COS(d4)
  d8 = y(i,1)
  d10 = hx*vpotx(i,1)
  d11 = SIN(d10)
  d13bar = (two/(hx*hx*sqn))
  d8bar = d13bar*d11
  d9bar = COS(d10)*(d13bar*d8)*hx
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hx
  gfknx1(1:gp) = d13bar*sx(1:gp,i+1,1) + d2bar*sx(1:gp,i,1) +  &
                 d3bar*svpotx(1:gp,i,1) + d8bar*sy(1:gp,i,1) + &
                 d9bar*svpotx(1:gp,i,1)
  fkinx1 = (two/(hx*hx*sqn))*(x(i+1,1)-d2*d5+d8*d11)
  d2 = y(i,1)
  d4 = hx*vpotx(i,1)
  d5 = COS(d4)
  d8 = x(i,1)
  d10 = hx*vpotx(i,1)
  d11 = SIN(d10)
  d13bar = (two/(hx*hx*sqn))
  d8bar = -d13bar*d11
  d9bar = COS(d10)*(-d13bar*d8)*hx
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hx
  gfknx2(1:gp) = d13bar*sy(1:gp,i+1,1) + d2bar*sy(1:gp,i,1) +  &
                 d3bar*svpotx(1:gp,i,1) + d8bar*sx(1:gp,i,1) + &
                 d9bar*svpotx(1:gp,i,1)
  fkinx2 = (two/(hx*hx*sqn))*(y(i+1,1)-d2*d5-d8*d11)
  d2 = x(i,1)
  d4 = hy*vpoty(i,1)
  d5 = COS(d4)
  d8 = y(i,1)
  d10 = hy*vpoty(i,1)
  d11 = SIN(d10)
  d13bar = (two/(hy*hy*sqn))
  d8bar = d13bar*d11
  d9bar = COS(d10)*(d13bar*d8)*hy
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hy
  gfkny1(1:gp) = d13bar*sx(1:gp,i,2) + d2bar*sx(1:gp,i,1) +  &
                 d3bar*svpoty(1:gp,i,1) + d8bar*sy(1:gp,i,1) + &
                 d9bar*svpoty(1:gp,i,1)
  fkiny1 = (two/(hy*hy*sqn))*(x(i,2)-d2*d5+d8*d11)
  d2 = y(i,1)
  d4 = hy*vpoty(i,1)
  d5 = COS(d4)
  d8 = x(i,1)
  d10 = hy*vpoty(i,1)
  d11 = SIN(d10)
  d13bar = (two/(hy*hy*sqn))
  d8bar = -d13bar*d11
  d9bar = COS(d10)*(-d13bar*d8)*hy
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hy
  gfkny2(1:gp) = d13bar*sy(1:gp,i,2) + d2bar*sy(1:gp,i,1) +  &
                 d3bar*svpoty(1:gp,i,1) + d8bar*sx(1:gp,i,1) + &
                 d9bar*svpoty(1:gp,i,1)
  fkiny2 = (two/(hy*hy*sqn))*(y(i,2)-d2*d5-d8*d11)
  d7bar = (one/hx)
  d3bar = (one/hy)
  gffld(1:gp) = d3bar*svpotx(1:gp,i,1) + (-d3bar*svpotx(1:gp,i,2)) +  &
                d7bar*svpoty(1:gp,i+1,1) + (-d7bar*svpoty(1:gp,i,1))
  ffield = (vpotx(i,1)-vpotx(i,2))/hy + (vpoty(i+1,1)-vpoty(i,1))/hx
  ffbar = (two*(tkappa**2)/sqn)
  gffld(1:gp) = ffbar*gffld(1:gp)
  ffield = (two*(tkappa**2)/sqn)*ffield
  d3 = hx*vpotx(i,1)
  d5 = -COS(d3)
  d9 = hx*vpotx(i,1)
  d11 = -SIN(d9)
  d15 = hy*vpoty(i,1)
  d17 = -COS(d15)
  d21 = hy*vpoty(i,1)
  d23 = -SIN(d21)
  d20bar = COS(d21)*(-fkiny2)*hy
  d14bar = (-SIN(d15)*(-fkiny1))*hy
  d8bar = COS(d9)*(-fkinx2)*hx
  d2bar = (-SIN(d3)*(-fkinx1))*hx
  yx(1:gp,i,1) = d5*gfknx1(1:gp) + d11*gfknx2(1:gp) +  &
                 d17*gfkny1(1:gp) + d23*gfkny2(1:gp) + yx(1:gp,i,1) +  &
                 d2bar*svpotx(1:gp,i,1) + d8bar*svpotx(1:gp,i,1) +  &
                 d14bar*svpoty(1:gp,i,1) + d20bar*svpoty(1:gp,i,1)
  d3 = hx*vpotx(i,1)
  d4 = SIN(d3)
  d8 = hx*vpotx(i,1)
  d10 = -COS(d8)
  d14 = hy*vpoty(i,1)
  d15 = SIN(d14)
  d19 = hy*vpoty(i,1)
  d21 = -COS(d19)
  d18bar = (-SIN(d19)*(-fkiny2))*hy
  d13bar = COS(d14)*fkiny1*hy
  d7bar = (-SIN(d8)*(-fkinx2))*hx
  d2bar = COS(d3)*fkinx1*hx
  yy(1:gp,i,1) = d4*gfknx1(1:gp) + d10*gfknx2(1:gp) +  &
                 d15*gfkny1(1:gp) + d21*gfkny2(1:gp) + yy(1:gp,i,1) +  &
                 d2bar*svpotx(1:gp,i,1) + d7bar*svpotx(1:gp,i,1) +  &
                 d13bar*svpoty(1:gp,i,1) + d18bar*svpoty(1:gp,i,1)
  d3 = hx*x(i,1)
  d5 = hx*vpotx(i,1)
  d6 = SIN(d5)
  d9 = hx*y(i,1)
  d11 = hx*vpotx(i,1)
  d12 = COS(d11)
  d14 = d3*d6 + d9*d12
  d18 = hx*y(i,1)
  d20 = hx*vpotx(i,1)
  d21 = SIN(d20)
  d24 = hx*x(i,1)
  d26 = hx*vpotx(i,1)
  d27 = COS(d26)
  d29 = d18*d21 - d24*d27
  ffbar = (one/hy)
  d25bar = (-SIN(d26)*(-fkinx2*d24))*hx
  d23bar = -fkinx2*d27*hx
  d19bar = COS(d20)*(fkinx2*d18)*hx
  d17bar = fkinx2*d21*hx
  d10bar = (-SIN(d11)*(fkinx1*d9))*hx
  d8bar = fkinx1*d12*hx
  d4bar = COS(d5)*(fkinx1*d3)*hx
  d2bar = fkinx1*d6*hx
  yvpotx(1:gp,i,1) = d14*gfknx1(1:gp) + d29*gfknx2(1:gp) +  &
                     ffbar*gffld(1:gp) + yvpotx(1:gp,i,1) +  &
                     d2bar*sx(1:gp,i,1) + d4bar*svpotx(1:gp,i,1) +  &
                     d8bar*sy(1:gp,i,1) + d10bar*svpotx(1:gp,i,1) +  &
                     d17bar*sy(1:gp,i,1) + d19bar*svpotx(1:gp,i,1) +  &
                     d23bar*sx(1:gp,i,1) + d25bar*svpotx(1:gp,i,1)
  d3 = hy*x(i,1)
  d5 = hy*vpoty(i,1)
  d6 = SIN(d5)
  d9 = hy*y(i,1)
  d11 = hy*vpoty(i,1)
  d12 = COS(d11)
  d14 = d3*d6 + d9*d12
  d18 = hy*y(i,1)
  d20 = hy*vpoty(i,1)
  d21 = SIN(d20)
  d24 = hy*x(i,1)
  d26 = hy*vpoty(i,1)
  d27 = COS(d26)
  d29 = d18*d21 - d24*d27
  ffbar = -(one/hx)
  d25bar = (-SIN(d26)*(-fkiny2*d24))*hy
  d23bar = -fkiny2*d27*hy
  d19bar = COS(d20)*(fkiny2*d18)*hy
  d17bar = fkiny2*d21*hy
  d10bar = (-SIN(d11)*(fkiny1*d9))*hy
  d8bar = fkiny1*d12*hy
  d4bar = COS(d5)*(fkiny1*d3)*hy
  d2bar = fkiny1*d6*hy
  yvpoty(1:gp,i,1) = d14*gfkny1(1:gp) + d29*gfkny2(1:gp) +  &
                     ffbar*gffld(1:gp) + yvpoty(1:gp,i,1) +  &
                     d2bar*sx(1:gp,i,1) + d4bar*svpoty(1:gp,i,1) +  &
                     d8bar*sy(1:gp,i,1) + d10bar*svpoty(1:gp,i,1) +  &
                     d17bar*sy(1:gp,i,1) + d19bar*svpoty(1:gp,i,1) +  &
                     d23bar*sx(1:gp,i,1) + d25bar*svpoty(1:gp,i,1)
  d2 = x(i-1,1)
  d4 = hx*vpotx(i-1,1)
  d5 = COS(d4)
  d8 = y(i-1,1)
  d10 = hx*vpotx(i-1,1)
  d11 = SIN(d10)
  d13bar = (two/(hx*hx*sqn))
  d8bar = d13bar*d11
  d9bar = COS(d10)*(d13bar*d8)*hx
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hx
  gfknx1(1:gp) = d13bar*sx(1:gp,i,1) + d2bar*sx(1:gp,i-1,1) +  &
                 d3bar*svpotx(1:gp,i-1,1) + d8bar*sy(1:gp,i-1,1) + &
                 d9bar*svpotx(1:gp,i-1,1)
  fkinx1 = (two/(hx*hx*sqn))*(x(i,1)-d2*d5+d8*d11)
  d2 = y(i-1,1)
  d4 = hx*vpotx(i-1,1)
  d5 = COS(d4)
  d8 = x(i-1,1)
  d10 = hx*vpotx(i-1,1)
  d11 = SIN(d10)
  d13bar = (two/(hx*hx*sqn))
  d8bar = -d13bar*d11
  d9bar = COS(d10)*(-d13bar*d8)*hx
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hx
  gfknx2(1:gp) = d13bar*sy(1:gp,i,1) + d2bar*sy(1:gp,i-1,1) +  &
                 d3bar*svpotx(1:gp,i-1,1) + d8bar*sx(1:gp,i-1,1) + &
                 d9bar*svpotx(1:gp,i-1,1)
  fkinx2 = (two/(hx*hx*sqn))*(y(i,1)-d2*d5-d8*d11)
  d2 = x(i,ny)
  d4 = hy*vpoty(i,ny)
  d5 = COS(d4)
  d8 = y(i,ny)
  d10 = hy*vpoty(i,ny)
  d11 = SIN(d10)
  d13bar = (two/(hy*hy*sqn))
  d8bar = d13bar*d11
  d9bar = COS(d10)*(d13bar*d8)*hy
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hy
  gfkny1(1:gp) = d13bar*sx(1:gp,i,ny+1) + d2bar*sx(1:gp,i,ny) +  &
                 d3bar*svpoty(1:gp,i,ny) + d8bar*sy(1:gp,i,ny) + &
                 d9bar*svpoty(1:gp,i,ny)
  fkiny1 = (two/(hy*hy*sqn))*(x(i,ny+1)-d2*d5+d8*d11)
  d2 = y(i,ny)
  d4 = hy*vpoty(i,ny)
  d5 = COS(d4)
  d8 = x(i,ny)
  d10 = hy*vpoty(i,ny)
  d11 = SIN(d10)
  d13bar = (two/(hy*hy*sqn))
  d8bar = -d13bar*d11
  d9bar = COS(d10)*(-d13bar*d8)*hy
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hy
  gfkny2(1:gp) = d13bar*sy(1:gp,i,ny+1) + d2bar*sy(1:gp,i,ny) +  &
                 d3bar*svpoty(1:gp,i,ny) + d8bar*sx(1:gp,i,ny) + &
                 d9bar*svpoty(1:gp,i,ny)
  fkiny2 = (two/(hy*hy*sqn))*(y(i,ny+1)-d2*d5-d8*d11)
  yx(1:gp,i,1) = gfknx1(1:gp) + gfkny1(1:gp) + yx(1:gp,i,1)
  yy(1:gp,i,1) = gfknx2(1:gp) + gfkny2(1:gp) + yy(1:gp,i,1)
  d7bar = (one/hx)
  d3bar = (one/hy)
  gffld(1:gp) = d3bar*svpotx(1:gp,i,ny) + (-d3bar*svpotx(1:gp,i,ny+1)) +  &
                d7bar*svpoty(1:gp,i+1,ny) + (-d7bar*svpoty(1:gp,i,ny))
  ffield = (vpotx(i,ny)-vpotx(i,ny+1))/hy + (vpoty(i+1,ny)-vpoty(i,ny))/hx
  ffbar = (two*(tkappa**2)/sqn)
  gffld(1:gp) = ffbar*gffld(1:gp)
  ffield = (two*(tkappa**2)/sqn)*ffield
  ffbar = -(one/hy)
  yvpotx(1:gp,i,1) = ffbar*gffld(1:gp) + yvpotx(1:gp,i,1)
  d7bar = (one/hx)
  d3bar = (one/hy)
  gffld(1:gp) = d3bar*svpotx(1:gp,i-1,1) + (-d3bar*svpotx(1:gp,i-1,2)) +  &
                d7bar*svpoty(1:gp,i,1) + (-d7bar*svpoty(1:gp,i-1,1))
  ffield = (vpotx(i-1,1)-vpotx(i-1,2))/hy + (vpoty(i,1)-vpoty(i-1,1))/hx
  ffbar = (two*(tkappa**2)/sqn)
  gffld(1:gp) = ffbar*gffld(1:gp)
  ffield = (two*(tkappa**2)/sqn)*ffield
  ffbar = (one/hx)
  yvpoty(1:gp,i,1) = ffbar*gffld(1:gp) + yvpoty(1:gp,i,1)
END DO

!     Left i = 1.

DO  j = 2, ny
  d2 = x(1,j)
  d4 = hx*vpotx(1,j)
  d5 = COS(d4)
  d8 = y(1,j)
  d10 = hx*vpotx(1,j)
  d11 = SIN(d10)
  d13bar = (two/(hx*hx*sqn))
  d8bar = d13bar*d11
  d9bar = COS(d10)*(d13bar*d8)*hx
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hx
  gfknx1(1:gp) = d13bar*sx(1:gp,2,j) + d2bar*sx(1:gp,1,j) +  &
                 d3bar*svpotx(1:gp,1,j) + d8bar*sy(1:gp,1,j) + &
                 d9bar*svpotx(1:gp,1,j)
  fkinx1 = (two/(hx*hx*sqn))*(x(2,j)-d2*d5+d8*d11)
  d2 = y(1,j)
  d4 = hx*vpotx(1,j)
  d5 = COS(d4)
  d8 = x(1,j)
  d10 = hx*vpotx(1,j)
  d11 = SIN(d10)
  d13bar = (two/(hx*hx*sqn))
  d8bar = -d13bar*d11
  d9bar = COS(d10)*(-d13bar*d8)*hx
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hx
  gfknx2(1:gp) = d13bar*sy(1:gp,2,j) + d2bar*sy(1:gp,1,j) +  &
                 d3bar*svpotx(1:gp,1,j) + d8bar*sx(1:gp,1,j) + &
                 d9bar*svpotx(1:gp,1,j)
  fkinx2 = (two/(hx*hx*sqn))*(y(2,j)-d2*d5-d8*d11)
  d2 = x(1,j)
  d4 = hy*vpoty(1,j)
  d5 = COS(d4)
  d8 = y(1,j)
  d10 = hy*vpoty(1,j)
  d11 = SIN(d10)
  d13bar = (two/(hy*hy*sqn))
  d8bar = d13bar*d11
  d9bar = COS(d10)*(d13bar*d8)*hy
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hy
  gfkny1(1:gp) = d13bar*sx(1:gp,1,j+1) + d2bar*sx(1:gp,1,j) +  &
                 d3bar*svpoty(1:gp,1,j) + d8bar*sy(1:gp,1,j) + &
                 d9bar*svpoty(1:gp,1,j)
  fkiny1 = (two/(hy*hy*sqn))*(x(1,j+1)-d2*d5+d8*d11)
  d2 = y(1,j)
  d4 = hy*vpoty(1,j)
  d5 = COS(d4)
  d8 = x(1,j)
  d10 = hy*vpoty(1,j)
  d11 = SIN(d10)
  d13bar = (two/(hy*hy*sqn))
  d8bar = -d13bar*d11
  d9bar = COS(d10)*(-d13bar*d8)*hy
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hy
  gfkny2(1:gp) = d13bar*sy(1:gp,1,j+1) + d2bar*sy(1:gp,1,j) +  &
                 d3bar*svpoty(1:gp,1,j) + d8bar*sx(1:gp,1,j) + &
                 d9bar*svpoty(1:gp,1,j)
  fkiny2 = (two/(hy*hy*sqn))*(y(1,j+1)-d2*d5-d8*d11)
  d7bar = (one/hx)
  d3bar = (one/hy)
  gffld(1:gp) = d3bar*svpotx(1:gp,1,j) + (-d3bar*svpotx(1:gp,1,j+1)) +  &
                d7bar*svpoty(1:gp,2,j) + (-d7bar*svpoty(1:gp,1,j))
  ffield = (vpotx(1,j)-vpotx(1,j+1))/hy + (vpoty(2,j)-vpoty(1,j))/hx
  ffbar = (two*(tkappa**2)/sqn)
  gffld(1:gp) = ffbar*gffld(1:gp)
  ffield = (two*(tkappa**2)/sqn)*ffield
  d3 = hx*vpotx(1,j)
  d5 = -COS(d3)
  d9 = hx*vpotx(1,j)
  d11 = -SIN(d9)
  d15 = hy*vpoty(1,j)
  d17 = -COS(d15)
  d21 = hy*vpoty(1,j)
  d23 = -SIN(d21)
  d20bar = COS(d21)*(-fkiny2)*hy
  d14bar = (-SIN(d15)*(-fkiny1))*hy
  d8bar = COS(d9)*(-fkinx2)*hx
  d2bar = (-SIN(d3)*(-fkinx1))*hx
  yx(1:gp,1,j) = d5*gfknx1(1:gp) + d11*gfknx2(1:gp) +  &
                 d17*gfkny1(1:gp) + d23*gfkny2(1:gp) + yx(1:gp,1,j) +  &
                 d2bar*svpotx(1:gp,1,j) + d8bar*svpotx(1:gp,1,j) +  &
                 d14bar*svpoty(1:gp,1,j) + d20bar*svpoty(1:gp,1,j)
  d3 = hx*vpotx(1,j)
  d4 = SIN(d3)
  d8 = hx*vpotx(1,j)
  d10 = -COS(d8)
  d14 = hy*vpoty(1,j)
  d15 = SIN(d14)
  d19 = hy*vpoty(1,j)
  d21 = -COS(d19)
  d18bar = (-SIN(d19)*(-fkiny2))*hy
  d13bar = COS(d14)*fkiny1*hy
  d7bar = (-SIN(d8)*(-fkinx2))*hx
  d2bar = COS(d3)*fkinx1*hx
  yy(1:gp,1,j) = d4*gfknx1(1:gp) + d10*gfknx2(1:gp) +  &
                 d15*gfkny1(1:gp) + d21*gfkny2(1:gp) + yy(1:gp,1,j) +  &
                 d2bar*svpotx(1:gp,1,j) + d7bar*svpotx(1:gp,1,j) +  &
                 d13bar*svpoty(1:gp,1,j) + d18bar*svpoty(1:gp,1,j)
  d3 = hx*x(1,j)
  d5 = hx*vpotx(1,j)
  d6 = SIN(d5)
  d9 = hx*y(1,j)
  d11 = hx*vpotx(1,j)
  d12 = COS(d11)
  d14 = d3*d6 + d9*d12
  d18 = hx*y(1,j)
  d20 = hx*vpotx(1,j)
  d21 = SIN(d20)
  d24 = hx*x(1,j)
  d26 = hx*vpotx(1,j)
  d27 = COS(d26)
  d29 = d18*d21 - d24*d27
  ffbar = (one/hy)
  d25bar = (-SIN(d26)*(-fkinx2*d24))*hx
  d23bar = -fkinx2*d27*hx
  d19bar = COS(d20)*(fkinx2*d18)*hx
  d17bar = fkinx2*d21*hx
  d10bar = (-SIN(d11)*(fkinx1*d9))*hx
  d8bar = fkinx1*d12*hx
  d4bar = COS(d5)*(fkinx1*d3)*hx
  d2bar = fkinx1*d6*hx
  yvpotx(1:gp,1,j) = d14*gfknx1(1:gp) + d29*gfknx2(1:gp) +  &
                     ffbar*gffld(1:gp) + yvpotx(1:gp,1,j) +  &
                     d2bar*sx(1:gp,1,j) + d4bar*svpotx(1:gp,1,j) +  &
                     d8bar*sy(1:gp,1,j) + d10bar*svpotx(1:gp,1,j) +  &
                     d17bar*sy(1:gp,1,j) + d19bar*svpotx(1:gp,1,j) +  &
                     d23bar*sx(1:gp,1,j) + d25bar*svpotx(1:gp,1,j)
  d3 = hy*x(1,j)
  d5 = hy*vpoty(1,j)
  d6 = SIN(d5)
  d9 = hy*y(1,j)
  d11 = hy*vpoty(1,j)
  d12 = COS(d11)
  d14 = d3*d6 + d9*d12
  d18 = hy*y(1,j)
  d20 = hy*vpoty(1,j)
  d21 = SIN(d20)
  d24 = hy*x(1,j)
  d26 = hy*vpoty(1,j)
  d27 = COS(d26)
  d29 = d18*d21 - d24*d27
  ffbar = -(one/hx)
  d25bar = (-SIN(d26)*(-fkiny2*d24))*hy
  d23bar = -fkiny2*d27*hy
  d19bar = COS(d20)*(fkiny2*d18)*hy
  d17bar = fkiny2*d21*hy
  d10bar = (-SIN(d11)*(fkiny1*d9))*hy
  d8bar = fkiny1*d12*hy
  d4bar = COS(d5)*(fkiny1*d3)*hy
  d2bar = fkiny1*d6*hy
  yvpoty(1:gp,1,j) = d14*gfkny1(1:gp) + d29*gfkny2(1:gp) +  &
                     ffbar*gffld(1:gp) + yvpoty(1:gp,1,j) +  &
                     d2bar*sx(1:gp,1,j) + d4bar*svpoty(1:gp,1,j) +  &
                     d8bar*sy(1:gp,1,j) + d10bar*svpoty(1:gp,1,j) +  &
                     d17bar*sy(1:gp,1,j) + d19bar*svpoty(1:gp,1,j) +  &
                     d23bar*sx(1:gp,1,j) + d25bar*svpoty(1:gp,1,j)
  d2 = x(nx,j)
  d4 = hx*vpotx(nx,j)
  d5 = COS(d4)
  d8 = y(nx,j)
  d10 = hx*vpotx(nx,j)
  d11 = SIN(d10)
  d13bar = (two/(hx*hx*sqn))
  d8bar = d13bar*d11
  d9bar = COS(d10)*(d13bar*d8)*hx
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hx
  gfknx1(1:gp) = d13bar*sx(1:gp,nx+1,j) + d2bar*sx(1:gp,nx,j) +  &
                 d3bar*svpotx(1:gp,nx,j) + d8bar*sy(1:gp,nx,j) + &
                 d9bar*svpotx(1:gp,nx,j)
  fkinx1 = (two/(hx*hx*sqn))*(x(nx+1,j)-d2*d5+d8*d11)
  d2 = y(nx,j)
  d4 = hx*vpotx(nx,j)
  d5 = COS(d4)
  d8 = x(nx,j)
  d10 = hx*vpotx(nx,j)
  d11 = SIN(d10)
  d13bar = (two/(hx*hx*sqn))
  d8bar = -d13bar*d11
  d9bar = COS(d10)*(-d13bar*d8)*hx
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hx
  gfknx2(1:gp) = d13bar*sy(1:gp,nx+1,j) + d2bar*sy(1:gp,nx,j) +  &
                 d3bar*svpotx(1:gp,nx,j) + d8bar*sx(1:gp,nx,j) + &
                 d9bar*svpotx(1:gp,nx,j)
  fkinx2 = (two/(hx*hx*sqn))*(y(nx+1,j)-d2*d5-d8*d11)
  d2 = x(1,j-1)
  d4 = hy*vpoty(1,j-1)
  d5 = COS(d4)
  d8 = y(1,j-1)
  d10 = hy*vpoty(1,j-1)
  d11 = SIN(d10)
  d13bar = (two/(hy*hy*sqn))
  d8bar = d13bar*d11
  d9bar = COS(d10)*(d13bar*d8)*hy
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hy
  gfkny1(1:gp) = d13bar*sx(1:gp,1,j) + d2bar*sx(1:gp,1,j-1) +  &
                 d3bar*svpoty(1:gp,1,j-1) + d8bar*sy(1:gp,1,j-1) + &
                 d9bar*svpoty(1:gp,1,j-1)
  fkiny1 = (two/(hy*hy*sqn))*(x(1,j)-d2*d5+d8*d11)
  d2 = y(1,j-1)
  d4 = hy*vpoty(1,j-1)
  d5 = COS(d4)
  d8 = x(1,j-1)
  d10 = hy*vpoty(1,j-1)
  d11 = SIN(d10)
  d13bar = (two/(hy*hy*sqn))
  d8bar = -d13bar*d11
  d9bar = COS(d10)*(-d13bar*d8)*hy
  d2bar = -d13bar*d5
  d3bar = (-SIN(d4)*(-d13bar*d2))*hy
  gfkny2(1:gp) = d13bar*sy(1:gp,1,j) + d2bar*sy(1:gp,1,j-1) +  &
                 d3bar*svpoty(1:gp,1,j-1) + d8bar*sx(1:gp,1,j-1) + &
                 d9bar*svpoty(1:gp,1,j-1)
  fkiny2 = (two/(hy*hy*sqn))*(y(1,j)-d2*d5-d8*d11)
  sfac = SIN(two*pi*vornum*(j-one)/DBLE(ny))
  cfac = COS(two*pi*vornum*(j-one)/DBLE(ny))
  yx(1:gp,1,j) = cfac*gfknx1(1:gp) + sfac*gfknx2(1:gp) + gfkny1(1:gp) + yx(1:gp,1,j)
  yy(1:gp,1,j) = -sfac*gfknx1(1:gp) + cfac*gfknx2(1:gp) + gfkny2(1:gp) + yy(1:gp,1,j)
  d7bar = (one/hx)
  d3bar = (one/hy)
  gffld(1:gp) = d3bar*svpotx(1:gp,1,j-1) + (-d3bar*svpotx(1:gp,1,j)) +  &
                d7bar*svpoty(1:gp,2,j-1) + (-d7bar*svpoty(1:gp,1,j-1))
  ffield = (vpotx(1,j-1)-vpotx(1,j))/hy + (vpoty(2,j-1)-vpoty(1,j-1))/hx
  ffbar = (two*(tkappa**2)/sqn)
  gffld(1:gp) = ffbar*gffld(1:gp)
  ffield = (two*(tkappa**2)/sqn)*ffield
  ffbar = -(one/hy)
  yvpotx(1:gp,1,j) = ffbar*gffld(1:gp) + yvpotx(1:gp,1,j)
  d7bar = (one/hx)
  d3bar = (one/hy)
  gffld(1:gp) = d3bar*svpotx(1:gp,nx,j) + (-d3bar*svpotx(1:gp,nx,j+1)) +  &
                d7bar*svpoty(1:gp,nx+1,j) + (-d7bar*svpoty(1:gp,nx,j))
  ffield = (vpotx(nx,j)-vpotx(nx,j+1))/hy + (vpoty(nx+1,j)-vpoty(nx,j))/hx
  ffbar = (two*(tkappa**2)/sqn)
  gffld(1:gp) = ffbar*gffld(1:gp)
  ffield = (two*(tkappa**2)/sqn)*ffield
  ffbar = (one/hx)
  yvpoty(1:gp,1,j) = ffbar*gffld(1:gp) + yvpoty(1:gp,1,j)
END DO

!     Kinetic energy part, at origin (only needed in zero field).

d2 = x(1,1)
d4 = hx*vpotx(1,1)
d5 = COS(d4)
d8 = y(1,1)
d10 = hx*vpotx(1,1)
d11 = SIN(d10)
d13bar = (two/(hx*hx*sqn))
d8bar = d13bar*d11
d9bar = COS(d10)*(d13bar*d8)*hx
d2bar = -d13bar*d5
d3bar = (-SIN(d4)*(-d13bar*d2))*hx
gfknx1(1:gp) = d13bar*sx(1:gp,2,1) + d2bar*sx(1:gp,1,1) +  &
               d3bar*svpotx(1:gp,1,1) + d8bar*sy(1:gp,1,1) + &
               d9bar*svpotx(1:gp,1,1)
fkinx1 = (two/(hx*hx*sqn))*(x(2,1)-d2*d5+d8*d11)
d2 = y(1,1)
d4 = hx*vpotx(1,1)
d5 = COS(d4)
d8 = x(1,1)
d10 = hx*vpotx(1,1)
d11 = SIN(d10)
d13bar = (two/(hx*hx*sqn))
d8bar = -d13bar*d11
d9bar = COS(d10)*(-d13bar*d8)*hx
d2bar = -d13bar*d5
d3bar = (-SIN(d4)*(-d13bar*d2))*hx
gfknx2(1:gp) = d13bar*sy(1:gp,2,1) + d2bar*sy(1:gp,1,1) +  &
               d3bar*svpotx(1:gp,1,1) + d8bar*sx(1:gp,1,1) + &
               d9bar*svpotx(1:gp,1,1)
fkinx2 = (two/(hx*hx*sqn))*(y(2,1)-d2*d5-d8*d11)
d2 = x(1,1)
d4 = hy*vpoty(1,1)
d5 = COS(d4)
d8 = y(1,1)
d10 = hy*vpoty(1,1)
d11 = SIN(d10)
d13bar = (two/(hy*hy*sqn))
d8bar = d13bar*d11
d9bar = COS(d10)*(d13bar*d8)*hy
d2bar = -d13bar*d5
d3bar = (-SIN(d4)*(-d13bar*d2))*hy
gfkny1(1:gp) = d13bar*sx(1:gp,1,2) + d2bar*sx(1:gp,1,1) +  &
               d3bar*svpoty(1:gp,1,1) + d8bar*sy(1:gp,1,1) + d9bar*svpoty(1:gp,1,1)
fkiny1 = (two/(hy*hy*sqn))*(x(1,2)-d2*d5+d8*d11)
d2 = y(1,1)
d4 = hy*vpoty(1,1)
d5 = COS(d4)
d8 = x(1,1)
d10 = hy*vpoty(1,1)
d11 = SIN(d10)
d13bar = (two/(hy*hy*sqn))
d8bar = -d13bar*d11
d9bar = COS(d10)*(-d13bar*d8)*hy
d2bar = -d13bar*d5
d3bar = (-SIN(d4)*(-d13bar*d2))*hy
gfkny2(1:gp) = d13bar*sy(1:gp,1,2) + d2bar*sy(1:gp,1,1) +  &
               d3bar*svpoty(1:gp,1,1) + d8bar*sx(1:gp,1,1) + &
               d9bar*svpoty(1:gp,1,1)
fkiny2 = (two/(hy*hy*sqn))*(y(1,2)-d2*d5-d8*d11)
d7bar = (one/hx)
d3bar = (one/hy)
gffld(1:gp) = d3bar*svpotx(1:gp,1,1) + (-d3bar*svpotx(1:gp,1,2)) +  &
              d7bar*svpoty(1:gp,2,1) + (-d7bar*svpoty(1:gp,1,1))
ffield = (vpotx(1,1)-vpotx(1,2))/hy + (vpoty(2,1)-vpoty(1,1))/hx
ffbar = (two*(tkappa**2)/sqn)
gffld(1:gp) = ffbar*gffld(1:gp)
ffield = (two*(tkappa**2)/sqn)*ffield
d3 = hx*vpotx(1,1)
d5 = -COS(d3)
d9 = hx*vpotx(1,1)
d11 = -SIN(d9)
d15 = hy*vpoty(1,1)
d17 = -COS(d15)
d21 = hy*vpoty(1,1)
d23 = -SIN(d21)
d20bar = COS(d21)*(-fkiny2)*hy
d14bar = (-SIN(d15)*(-fkiny1))*hy
d8bar = COS(d9)*(-fkinx2)*hx
d2bar = (-SIN(d3)*(-fkinx1))*hx
yx(1:gp,1,1) = d5*gfknx1(1:gp) + d11*gfknx2(1:gp) + d17*gfkny1(1:gp) +  &
               d23*gfkny2(1:gp) + yx(1:gp,1,1) +  &
               d2bar*svpotx(1:gp,1,1) + d8bar*svpotx(1:gp,1,1) +  &
               d14bar*svpoty(1:gp,1,1) + d20bar*svpoty(1:gp,1,1)
d3 = hx*vpotx(1,1)
d4 = SIN(d3)
d8 = hx*vpotx(1,1)
d10 = -COS(d8)
d14 = hy*vpoty(1,1)
d15 = SIN(d14)
d19 = hy*vpoty(1,1)
d21 = -COS(d19)
d18bar = (-SIN(d19)*(-fkiny2))*hy
d13bar = COS(d14)*fkiny1*hy
d7bar = (-SIN(d8)*(-fkinx2))*hx
d2bar = COS(d3)*fkinx1*hx
yy(1:gp,1,1) = d4*gfknx1(1:gp) + d10*gfknx2(1:gp) + d15*gfkny1(1:gp) +  &
               d21*gfkny2(1:gp) + yy(1:gp,1,1) +  &
               d2bar*svpotx(1:gp,1,1) + d7bar*svpotx(1:gp,1,1) +  &
               d13bar*svpoty(1:gp,1,1) + d18bar*svpoty(1:gp,1,1)
d3 = hx*x(1,1)
d5 = hx*vpotx(1,1)
d6 = SIN(d5)
d9 = hx*y(1,1)
d11 = hx*vpotx(1,1)
d12 = COS(d11)
d14 = d3*d6 + d9*d12
d18 = hx*y(1,1)
d20 = hx*vpotx(1,1)
d21 = SIN(d20)
d24 = hx*x(1,1)
d26 = hx*vpotx(1,1)
d27 = COS(d26)
d29 = d18*d21 - d24*d27
ffbar = (one/hy)
d25bar = (-SIN(d26)*(-fkinx2*d24))*hx
d23bar = -fkinx2*d27*hx
d19bar = COS(d20)*(fkinx2*d18)*hx
d17bar = fkinx2*d21*hx
d10bar = (-SIN(d11)*(fkinx1*d9))*hx
d8bar = fkinx1*d12*hx
d4bar = COS(d5)*(fkinx1*d3)*hx
d2bar = fkinx1*d6*hx
yvpotx(1:gp,1,1) = d14*gfknx1(1:gp) + d29*gfknx2(1:gp) +  &
                   ffbar*gffld(1:gp) + yvpotx(1:gp,1,1) +  &
                   d2bar*sx(1:gp,1,1) + d4bar*svpotx(1:gp,1,1) +  &
                   d8bar*sy(1:gp,1,1) + d10bar*svpotx(1:gp,1,1) +  &
                   d17bar*sy(1:gp,1,1) + d19bar*svpotx(1:gp,1,1) +  &
                   d23bar*sx(1:gp,1,1) + d25bar*svpotx(1:gp,1,1)
d3 = hy*x(1,1)
d5 = hy*vpoty(1,1)
d6 = SIN(d5)
d9 = hy*y(1,1)
d11 = hy*vpoty(1,1)
d12 = COS(d11)
d14 = d3*d6 + d9*d12
d18 = hy*y(1,1)
d20 = hy*vpoty(1,1)
d21 = SIN(d20)
d24 = hy*x(1,1)
d26 = hy*vpoty(1,1)
d27 = COS(d26)
d29 = d18*d21 - d24*d27
ffbar = -(one/hx)
d25bar = (-SIN(d26)*(-fkiny2*d24))*hy
d23bar = -fkiny2*d27*hy
d19bar = COS(d20)*(fkiny2*d18)*hy
d17bar = fkiny2*d21*hy
d10bar = (-SIN(d11)*(fkiny1*d9))*hy
d8bar = fkiny1*d12*hy
d4bar = COS(d5)*(fkiny1*d3)*hy
d2bar = fkiny1*d6*hy
yvpoty(1:gp,1,1) = d14*gfkny1(1:gp) + d29*gfkny2(1:gp) +  &
                   ffbar*gffld(1:gp) + yvpoty(1:gp,1,1) +  &
                   d2bar*sx(1:gp,1,1) + d4bar*svpoty(1:gp,1,1) +  &
                   d8bar*sy(1:gp,1,1) + d10bar*svpoty(1:gp,1,1) +  &
                   d17bar*sy(1:gp,1,1) + d19bar*svpoty(1:gp,1,1) +  &
                   d23bar*sx(1:gp,1,1) + d25bar*svpoty(1:gp,1,1)
d2 = x(nx,1)
d4 = hx*vpotx(nx,1)
d5 = COS(d4)
d8 = y(nx,1)
d10 = hx*vpotx(nx,1)
d11 = SIN(d10)
d13bar = (two/(hx*hx*sqn))
d8bar = d13bar*d11
d9bar = COS(d10)*(d13bar*d8)*hx
d2bar = -d13bar*d5
d3bar = (-SIN(d4)*(-d13bar*d2))*hx
gfknx1(1:gp) = d13bar*sx(1:gp,nx+1,1) + d2bar*sx(1:gp,nx,1) +  &
               d3bar*svpotx(1:gp,nx,1) + d8bar*sy(1:gp,nx,1) + &
               d9bar*svpotx(1:gp,nx,1)
fkinx1 = (two/(hx*hx*sqn))*(x(nx+1,1)-d2*d5+d8*d11)
d2 = y(nx,1)
d4 = hx*vpotx(nx,1)
d5 = COS(d4)
d8 = x(nx,1)
d10 = hx*vpotx(nx,1)
d11 = SIN(d10)
d13bar = (two/(hx*hx*sqn))
d8bar = -d13bar*d11
d9bar = COS(d10)*(-d13bar*d8)*hx
d2bar = -d13bar*d5
d3bar = (-SIN(d4)*(-d13bar*d2))*hx
gfknx2(1:gp) = d13bar*sy(1:gp,nx+1,1) + d2bar*sy(1:gp,nx,1) +  &
               d3bar*svpotx(1:gp,nx,1) + d8bar*sx(1:gp,nx,1) + &
               d9bar*svpotx(1:gp,nx,1)
fkinx2 = (two/(hx*hx*sqn))*(y(nx+1,1)-d2*d5-d8*d11)
d2 = x(1,ny)
d4 = hy*vpoty(1,ny)
d5 = COS(d4)
d8 = y(1,ny)
d10 = hy*vpoty(1,ny)
d11 = SIN(d10)
d13bar = (two/(hy*hy*sqn))
d8bar = d13bar*d11
d9bar = COS(d10)*(d13bar*d8)*hy
d2bar = -d13bar*d5
d3bar = (-SIN(d4)*(-d13bar*d2))*hy
gfkny1(1:gp) = d13bar*sx(1:gp,1,ny+1) + d2bar*sx(1:gp,1,ny) +  &
               d3bar*svpoty(1:gp,1,ny) + d8bar*sy(1:gp,1,ny) + &
               d9bar*svpoty(1:gp,1,ny)
fkiny1 = (two/(hy*hy*sqn))*(x(1,ny+1)-d2*d5+d8*d11)
d2 = y(1,ny)
d4 = hy*vpoty(1,ny)
d5 = COS(d4)
d8 = x(1,ny)
d10 = hy*vpoty(1,ny)
d11 = SIN(d10)
d13bar = (two/(hy*hy*sqn))
d8bar = -d13bar*d11
d9bar = COS(d10)*(-d13bar*d8)*hy
d2bar = -d13bar*d5
d3bar = (-SIN(d4)*(-d13bar*d2))*hy
gfkny2(1:gp) = d13bar*sy(1:gp,1,ny+1) + d2bar*sy(1:gp,1,ny) +  &
               d3bar*svpoty(1:gp,1,ny) + d8bar*sx(1:gp,1,ny) + &
               d9bar*svpoty(1:gp,1,ny)
fkiny2 = (two/(hy*hy*sqn))*(y(1,ny+1)-d2*d5-d8*d11)
yx(1:gp,1,1) = gfknx1(1:gp) + gfkny1(1:gp) + yx(1:gp,1,1)
yy(1:gp,1,1) = gfknx2(1:gp) + gfkny2(1:gp) + yy(1:gp,1,1)
d7bar = (one/hx)
d3bar = (one/hy)
gffld(1:gp) = d3bar*svpotx(1:gp,1,ny) + (-d3bar*svpotx(1:gp,1,ny+1)) +  &
              d7bar*svpoty(1:gp,2,ny) + (-d7bar*svpoty(1:gp,1,ny))
ffield = (vpotx(1,ny)-vpotx(1,ny+1))/hy + (vpoty(2,ny)-vpoty(1,ny))/hx
ffbar = (two*(tkappa**2)/sqn)
gffld(1:gp) = ffbar*gffld(1:gp)
ffield = (two*(tkappa**2)/sqn)*ffield
ffbar = -(one/hy)
yvpotx(1:gp,1,1) = ffbar*gffld(1:gp) + yvpotx(1:gp,1,1)
d7bar = (one/hx)
d3bar = (one/hy)
gffld(1:gp) = d3bar*svpotx(1:gp,nx,1) + (-d3bar*svpotx(1:gp,nx,2)) +  &
              d7bar*svpoty(1:gp,nx+1,1) + (-d7bar*svpoty(1:gp,nx,1))
ffield = (vpotx(nx,1)-vpotx(nx,2))/hy + (vpoty(nx+1,1)-vpoty(nx,1))/hx
ffbar = (two*(tkappa**2)/sqn)
gffld(1:gp) = ffbar*gffld(1:gp)
ffield = (two*(tkappa**2)/sqn)*ffield
ffbar = (one/hx)
yvpoty(1:gp,1,1) = ffbar*gffld(1:gp) + yvpoty(1:gp,1,1)

RETURN
END SUBROUTINE dgl2co

END PROGRAM dmain
