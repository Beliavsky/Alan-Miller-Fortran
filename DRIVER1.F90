MODULE Common_Bounds
IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

!     PARAMETERS
INTEGER, PARAMETER  :: nmax = 100

!     COMMON ARRAYS
REAL (dp), SAVE  :: l(nmax), u(nmax)

END MODULE Common_Bounds



PROGRAM spgma1
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-04  Time: 17:43:49

!     SPG Test driver. This master file uses SPG to solve
!     an easy (very easy) bound constrained problem.

!     This version 02 FEB 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.
!     Final revision 03 JUL 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.

USE Spectral_Projected_Grad
USE Common_Bounds
IMPLICIT NONE

!     LOCAL SCALARS
REAL (dp)  :: pginfn, pgtwon, eps, eps2, f
INTEGER    :: fcnt, flag, gcnt, i, iter, m, maxfc, maxit, n
LOGICAL    :: output

!     LOCAL ARRAYS
REAL (dp)  :: x(nmax)

!     COMMON BLOCKS
! COMMON /bounds/ l, u

!     DEFINE DIMENSION OF THE PROBLEM

n = 10

!     DEFINE BOUNDS

DO i = 1, n
  l(i) = -100.0D0
  u(i) = 50.0D0
END DO

!     DEFINE INITIAL POINT

x(1:n) = 60.0_dp

!     SET UP OPTIMIZER PARAMETERS

output = .true.
maxit = 1000
maxfc = 2000
eps = 0.0D0
eps2 = 1.0D-6
m = 10

!     CALL THE OPTIMIZER

CALL spg(n, x, m, eps, eps2, maxit, maxfc, output, f, pginfn, pgtwon, iter,  &
         fcnt, gcnt, flag)

!     WRITE STATISTICS

WRITE (*,FMT=5000) f, pginfn, SQRT(pgtwon), flag
WRITE (*,FMT=5100) iter, fcnt, gcnt

STOP

5000 FORMAT (/' F = ', g17.10 / ' PGINFNORM = ', g16.10 /  &
             ' PGTWONORM =  ', g16.10 / ' FLAG = ', i1)
5100 FORMAT (/' ITER = ', i10 / ' FCNT = ', i10 / ' GCNT = ', i10)

END PROGRAM spgma1



SUBROUTINE evalf(n, x, f, inform)

!  This subroutine computes the objective function.

!  On Entry:

!  n     integer,
!        size of the problem,

!  x     REAL (dp) x(n),
!        point at which the function will be evaluated.

!  On Return

!  f     REAL (dp),
!        function value at x,

!  inform integer,
!        termination parameter:
!        0 = the function was successfully evaluated,
!        1 = some error occurs in the function evaluation.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(n)
REAL (dp), INTENT(OUT)  :: f
INTEGER, INTENT(OUT)    :: inform

inform = 0

f = SUM( x(1:n)**2 )

RETURN
END SUBROUTINE evalf



SUBROUTINE evalg(n, x, g, inform)

!  This subroutine computes the gradient of the objective function.

!  On Entry:

!  n     integer,
!        size of the problem,

!  x     REAL (dp) x(n),
!        point at which the gradient will be evaluated.

!  On Return

!  g     REAL (dp) g(n),
!        gradient vector at x,

!  inform integer,
!        termination parameter:
!        0 = the gradient was successfully evaluated,
!        1 = some error occurs in the gradient evaluation.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(n)
REAL (dp), INTENT(OUT)  :: g(n)
INTEGER, INTENT(OUT)    :: inform

inform = 0

g(1:n) = 2.0_dp * x(1:n)

RETURN
END SUBROUTINE evalg



SUBROUTINE proj(n, x, inform)

!  This subroutine computes the projection of an arbitrary point onto the
!  feasible set.  Since the feasible set can be described in many ways, its
!  information is in a common block.

!  On Entry:

!  n     integer,
!        size of the problem,

!  x     REAL (dp) x(n),
!        point that will be projected.

!  On Return

!  x     REAL (dp) x(n),
!        projected point,

!  inform integer,
!        termination parameter:
!        0 = the projection was successfully done,
!        1 = some error occurs in the projection.

!     COMMON ARRAYS
! REAL (dp)  :: l(nmax), u(nmax)
!     COMMON BLOCKS
! COMMON /bounds/ l, u

USE Common_Bounds
IMPLICIT NONE

!     PARAMETERS

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
INTEGER, INTENT(OUT)       :: inform

!     LOCAL SCALARS
INTEGER  :: i

inform = 0

DO i = 1, n
  x(i) = MAX(l(i), MIN(x(i),u(i)))
END DO

RETURN
END SUBROUTINE proj
