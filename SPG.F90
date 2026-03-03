MODULE Spectral_Projected_Grad

!   ALGORITHM 813, COLLECTED ALGORITHMS FROM ACM.
!   THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!   VOL. 27,NO. 3, September, 2001, P.  340--349.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

! Interfaces to user-supplied functions

INTERFACE
  SUBROUTINE evalf(n, x, f, inform)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: f
    INTEGER, INTENT(OUT)    :: inform
  END SUBROUTINE evalf

  SUBROUTINE evalg(n, x, g, inform)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: x(:)
    REAL (dp), INTENT(OUT)  :: g(:)
    INTEGER, INTENT(OUT)    :: inform
  END SUBROUTINE evalg

  SUBROUTINE proj(n, x, inform)
    IMPLICIT NONE
    INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN OUT)  :: x(:)
    INTEGER, INTENT(OUT)       :: inform
  END SUBROUTINE proj
END INTERFACE

PRIVATE
PUBLIC  :: spg

CONTAINS


SUBROUTINE spg(n, x, m, eps, eps2, maxit, maxfc, output, f, pginfn, pgtwon,  &
               iter, fcnt, gcnt, flag)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2002-05-04  Time: 17:43:43

!  Subroutine SPG implements the Spectral Projected Gradient Method (Version 2:
!  "continuous projected gradient direction") to find the local minimizers of
!  a given function with convex constraints, described in

!  E. G. Birgin, J. M. Martinez, and M. Raydan, "Nonmonotone spectral
!  projected gradient methods on convex sets", SIAM Journal on Optimization
!  10, pp. 1196-1211, 2000.

!  and

!  E. G. Birgin, J. M. Martinez, and M. Raydan, "Algorithm 813: SPG: software
!  for convex-constrained optimization", ACM Transactions on Mathematical
!  Software, Vol.27(3), pp. 340-9, 2001.

!  The user must supply the external subroutines evalf, evalg and proj to
!  evaluate the objective function and its gradient and to project an
!  arbitrary point onto the feasible region.

!  This version 17 JAN 2000 by E.G.Birgin, J.M.Martinez and M.Raydan.
!  Reformatted 03 OCT 2000 by Tim Hopkins.
!  Final revision 03 JUL 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.

!  On Entry:

!  n     integer,
!        size of the problem,

!  x     REAL (dp) x(n),
!        initial guess,

!  m     integer,
!        number of previous function values to be considered
!        in the nonmonotone line search,

!  eps   REAL (dp),
!        stopping criterion: ||projected grad||_inf < eps,

!  eps2  REAL (dp),
!        stopping criterion: ||projected grad||_2 < eps2,

!  maxit integer,
!        maximum number of iterations,

!  maxfc integer,
!        maximum number of function evaluations,

!  output logical,
!        true: print some information at each iteration,
!        false: no print.

!  On Return:

!  x     REAL (dp) x(n),
!        approximation to the local minimizer,

!  f     REAL (dp),
!        function value at the approximation to the local minimizer,

!  pginfn REAL (dp),
!        ||projected grad||_inf at the final iteration,

!  pgtwon REAL (dp),
!        ||projected grad||_2^2 at the final iteration,

!  iter  integer,
!        number of iterations,

!  fcnt  integer,
!        number of function evaluations,

!  gcnt  integer,
!        number of gradient evaluations,

!  flag  integer,
!        termination parameter:
!        0= convergence with projected gradient infinite-norm,
!        1= convergence with projected gradient 2-norm,
!        2= too many iterations,
!        3= too many function evaluations,
!        4= error in proj subroutine,
!        5= error in evalf subroutine,
!        6= error in evalg subroutine.

!  PARAMETERS

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
INTEGER, INTENT(IN)        :: m
REAL (dp), INTENT(IN)      :: eps
REAL (dp), INTENT(IN)      :: eps2
INTEGER, INTENT(IN)        :: maxit
INTEGER, INTENT(IN)        :: maxfc
LOGICAL, INTENT(IN)        :: output
REAL (dp), INTENT(OUT)     :: f
REAL (dp), INTENT(OUT)     :: pginfn
REAL (dp), INTENT(OUT)     :: pgtwon
INTEGER, INTENT(OUT)       :: iter
INTEGER, INTENT(OUT)       :: fcnt
INTEGER, INTENT(OUT)       :: gcnt
INTEGER, INTENT(OUT)       :: flag

REAL (dp), PARAMETER  :: lmin = 1.0D-30, lmax = 1.0D+30

!     LOCAL SCALARS
REAL (dp)  :: fbest, fnew, gtd, lambda, sts, sty
INTEGER    :: i, lsflag, inform

!     LOCAL
REAL (dp)  :: pg(n), g(n), gnew(n), s(n), y(n),  &
              d(n), xbest(n), xnew(n), lastfv(0:m-1)

!     EXTERNAL SUBROUTINES
! EXTERNAL ls, evalf, evalg, proj

!     INITIALIZATION

iter = 0
fcnt = 0
gcnt = 0

lastfv(0:m-1) = -1.0D+99

!     PROJECT INITIAL GUESS

CALL proj(n, x, inform)

IF (inform /= 0) THEN
! ERROR IN PROJ SUBROUTINE, STOP
  flag = 4
  GO TO 20
END IF

!     INITIALIZE BEST SOLUTION

xbest(1:n) = x(1:n)

!     EVALUATE FUNCTION AND GRADIENT

CALL evalf(n, x, f, inform)
fcnt = fcnt + 1

IF (inform /= 0) THEN
! ERROR IN EVALF SUBROUTINE, STOP
  flag = 5
  GO TO 20
END IF

CALL evalg(n, x, g, inform)
gcnt = gcnt + 1

IF (inform /= 0) THEN
! ERROR IN EVALG SUBROUTINE, STOP
  flag = 6
  GO TO 20
END IF

!     STORE FUNCTION VALUE FOR THE NONMONOTONE LINE SEARCH

lastfv(0) = f

!     INITIALIZE BEST FUNCTION VALUE

fbest = f

!     COMPUTE CONTINUOUS PROJECTED GRADIENT (AND ITS NORMS)

pg(1:n) = x(1:n) - g(1:n)

CALL proj(n, pg, inform)

IF (inform /= 0) THEN
! ERROR IN PROJ SUBROUTINE, STOP
  flag = 4
  GO TO 20
END IF

pgtwon = 0.0D0
pginfn = 0.0D0
DO i = 1, n
  pg(i) = pg(i) - x(i)
  pgtwon = pgtwon + pg(i) ** 2
  pginfn = MAX(pginfn, ABS(pg(i)))
END DO

!     PRINT ITERATION INFORMATION

IF (output) THEN
  WRITE (*,FMT=5000) iter, f, pginfn
END IF

!     DEFINE INITIAL SPECTRAL STEPLENGTH

IF (pginfn /= 0.0D0) THEN
  lambda = MIN(lmax, MAX(lmin, 1.0D0/pginfn))
END IF

!     MAIN LOOP

!     TEST STOPPING CRITERIA

10 IF (pginfn <= eps) THEN
! GRADIENT INFINITE-NORM STOPPING CRITERION SATISFIED, STOP
  flag = 0
  GO TO 20
END IF

IF (pgtwon <= eps2**2) THEN
! GRADIENT 2-NORM STOPPING CRITERION SATISFIED, STOP
  flag = 1
  GO TO 20
END IF

IF (iter > maxit) THEN
! MAXIMUM NUMBER OF ITERATIONS EXCEEDED, STOP
  flag = 2
  GO TO 20
END IF

IF (fcnt > maxfc) THEN
! MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED, STOP
  flag = 3
  GO TO 20
END IF

!     DO AN ITERATION

iter = iter + 1

!     COMPUTE THE SPECTRAL PROJECTED GRADIENT DIRECTION
!     AND <G,D>

d(1:n) = x(1:n) - lambda * g(1:n)

CALL proj(n, d, inform)

IF (inform /= 0) THEN
! ERROR IN PROJ SUBROUTINE, STOP
  flag = 4
  GO TO 20
END IF

gtd = 0.0D0
DO i = 1, n
  d(i) = d(i) - x(i)
  gtd = gtd + g(i) * d(i)
END DO

!     NONMONOTONE LINE SEARCH

CALL ls(n, x, f, d, gtd, m, lastfv, maxfc, fcnt, fnew, xnew, lsflag)

IF (lsflag == 3) THEN
! THE NUMBER OF FUNCTION EVALUATIONS WAS EXCEEDED
! INSIDE  THE LINE SEARCH, STOP
  flag = 3
  GO TO 20
END IF

!     SET NEW FUNCTION VALUE AND SAVE IT FOR THE NONMONOTONE
!     LINE SEARCH

f = fnew
lastfv(MOD(iter,m)) = f

!     COMPARE THE NEW FUNCTION VALUE AGAINST THE BEST FUNCTION VALUE AND,
!     IF SMALLER, UPDATE THE BEST FUNCTION VALUE AND THE CORRESPONDING
!     BEST POINT

IF (f < fbest) THEN
  fbest = f
  xbest(1:n) = xnew(1:n)
END IF

!     EVALUATE THE GRADIENT AT THE NEW ITERATE

CALL evalg(n, xnew, gnew, inform)
gcnt = gcnt + 1

IF (inform /= 0) THEN
! ERROR IN EVALG SUBROUTINE, STOP
  flag = 6
  GO TO 20
END IF

!     COMPUTE S = XNEW - X, Y = GNEW - G, <S,S>, <S,Y>,
!     THE CONTINUOUS PROJECTED GRADIENT AND ITS NORMS

sts = 0.0D0
sty = 0.0D0
DO i = 1, n
  s(i) = xnew(i) - x(i)
  y(i) = gnew(i) - g(i)
  sts = sts + s(i) * s(i)
  sty = sty + s(i) * y(i)
  x(i) = xnew(i)
  g(i) = gnew(i)
  pg(i) = x(i) - g(i)
END DO

CALL proj(n, pg, inform)

IF (inform /= 0) THEN
! ERROR IN PROJ SUBROUTINE, STOP
  flag = 4
  GO TO 20
END IF

pgtwon = 0.0D0
pginfn = 0.0D0
DO i = 1, n
  pg(i) = pg(i) - x(i)
  pgtwon = pgtwon + pg(i) ** 2
  pginfn = MAX(pginfn, ABS(pg(i)))
END DO

!     PRINT ITERATION INFORMATION

IF (output) THEN
  WRITE (*,FMT=5000) iter, f, pginfn
END IF

!     COMPUTE SPECTRAL STEPLENGTH

IF (sty <= 0.0D0) THEN
  lambda = lmax
  
ELSE
  lambda = MIN(lmax, MAX(lmin,sts/sty))
END IF

!     FINISH THE ITERATION

GO TO 10

!     STOP

!     SET X AND F WITH THE BEST SOLUTION AND ITS FUNCTION VALUE

20 f = fbest
x(1:n) = xbest(1:n)
RETURN

5000 FORMAT (' ITER= ', i10,' F= ', g17.10, '  PGINFNORM= ', g17.10)
END SUBROUTINE spg



SUBROUTINE ls(n, x, f, d, gtd, m, lastfv, maxfc, fcnt, fnew, xnew, flag)

!  Subroutine LS implements a nonmonotone line search with
!  safeguarded quadratic interpolation.

!  This version 17 JAN 2000 by E.G.Birgin, J.M.Martinez and M.Raydan.
!  Reformatted 03 OCT 2000 by Tim Hopkins.
!  Final revision 03 JUL 2001 by E.G.Birgin, J.M.Martinez and M.Raydan.

!  On Entry:

!  n     integer,
!        size of the problem,

!  x     REAL (dp) x(n),
!        initial guess,

!  f     REAL (dp),
!        function value at the actual point,

!  d     REAL (dp) d(n),
!        search direction,

!  gtd   REAL (dp),
!        internal product <g,d>, where g is the gradient at x,

!  m     integer,
!        number of previous function values to be considered
!        in the nonmonotone line search,

!  lastfv REAL (dp) lastfv(m),
!        last m function values,

!  maxfc integer,
!        maximum number of function evaluations,

!  fcnt  integer,
!        actual number of function evaluations.

!  On Return:

!  fcnt  integer,
!        actual number of function evaluations,

!  fnew  REAL (dp),
!        function value at the new point,

!  xnew  REAL (dp) xnew(n),
!        new point,

!  flag  integer,
!        0= convergence with nonmonotone Armijo-like criterion,
!        3= too many function evaluations,
!        5= error in evalf subroutine.

!     PARAMETERS

INTEGER, INTENT(IN)      :: n
REAL (dp), INTENT(IN)    :: x(n)
REAL (dp), INTENT(IN)    :: f
REAL (dp), INTENT(IN)    :: d(n)
REAL (dp), INTENT(IN)    :: gtd
INTEGER, INTENT(IN)      :: m
REAL (dp), INTENT(IN)    :: lastfv(0:m-1)
INTEGER, INTENT(IN)      :: maxfc
INTEGER, INTENT(IN OUT)  :: fcnt
REAL (dp), INTENT(OUT)   :: fnew
REAL (dp), INTENT(OUT)   :: xnew(n)
INTEGER, INTENT(OUT)     :: flag

REAL (dp), PARAMETER  :: gamma = 1.0D-04

!     LOCAL SCALARS
REAL (dp)  :: alpha, atemp, fmax
INTEGER    :: i, inform

!     EXTERNAL SUBROUTINES
! EXTERNAL evalf

!     INITIALIZATION

!     COMPUTE THE MAXIMUM FUNCTIONAL VALUE OF THE LAST M ITERATIONS

fmax = lastfv(0)
DO i = 1, m - 1
  fmax = MAX(fmax, lastfv(i))
END DO

!     COMPUTE FIRST TRIAL

alpha = 1.0D0

xnew(1:n) = x(1:n) + d(1:n)

CALL evalf(n, xnew, fnew, inform)
fcnt = fcnt + 1

IF (inform /= 0) THEN
! ERROR IN EVALF SUBROUTINE, STOP
  flag = 5
  GO TO 20
END IF

!     MAIN LOOP

!     TEST STOPPING CRITERIA

10 IF (fnew <= fmax + gamma*alpha*gtd) THEN
! NONMONOTONE ARMIJO-LIKE STOPPING CRITERION SATISFIED, STOP
  flag = 0
  GO TO 20
END IF

IF (fcnt >= maxfc) THEN
! MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED, STOP
  flag = 3
  GO TO 20
END IF

!     DO AN ITERATION

!     SAFEGUARDED QUADRATIC INTERPOLATION

IF (alpha <= 0.1D0) THEN
  alpha = alpha / 2.0D0
  
ELSE
  atemp = (-gtd*alpha**2) / (2.0D0*(fnew - f - alpha*gtd))
  IF (atemp < 0.1D0 .OR. atemp > 0.9D0*alpha) THEN
    atemp = alpha / 2.0D0
  END IF
  alpha = atemp
END IF

!     COMPUTE TRIAL POINT

xnew(1:n) = x(1:n) + alpha * d(1:n)

!     EVALUATE FUNCTION

CALL evalf(n, xnew, fnew, inform)
fcnt = fcnt + 1

IF (inform /= 0) THEN
! ERROR IN EVALF SUBROUTINE, STOP
  flag = 5
  GO TO 20
END IF

!     ITERATE

GO TO 10

!     STOP

20 RETURN

END SUBROUTINE ls

END MODULE Spectral_Projected_Grad
