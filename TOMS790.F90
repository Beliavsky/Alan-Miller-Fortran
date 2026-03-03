MODULE Cubic_Shepard

!      ALGORITHM 790, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 25,NO. 1,     March, 1999, P.   70-- 73.

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

CONTAINS


SUBROUTINE cshep2 (n, x, y, f, nc, nw, nr,  lcell, lnext, xmin,  &
                   ymin, dx, dy, rmax, rw, a, ier)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-08-10  Time: 20:55:57

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN)      :: y(:)
REAL (dp), INTENT(IN)      :: f(:)
INTEGER, INTENT(IN)        :: nc
INTEGER, INTENT(IN)        :: nw
INTEGER, INTENT(IN)        :: nr
INTEGER, INTENT(IN OUT)    :: lcell(:,:)     ! lcell(nr,nr)
INTEGER, INTENT(OUT)       :: lnext(:)
REAL (dp), INTENT(OUT)     :: xmin
REAL (dp), INTENT(OUT)     :: ymin
REAL (dp), INTENT(OUT)     :: dx
REAL (dp), INTENT(OUT)     :: dy
REAL (dp), INTENT(OUT)     :: rmax
REAL (dp), INTENT(OUT)     :: rw(:)
REAL (dp), INTENT(IN OUT)  :: a(:,:)     ! a(9,n)
INTEGER, INTENT(OUT)       :: ier

!***********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/13/97

!   This subroutine computes a set of parameters defining a C2 (twice
! continuously differentiable) bivariate function C(X,Y) which interpolates
! data values F at a set of N arbitrarily distributed points (X,Y) in the
! plane (nodes).
! The interpolant C may be evaluated at an arbitrary point by function CS2VAL,
! and its first partial derivatives are computed by Subroutine CS2GRD.

!   The interpolation scheme is a modified Cubic Shepard method:

! C = [W(1)*C(1)+W(2)*C(2)+..+W(N)*C(N)]/[W(1)+W(2)+..+W(N)]

! for bivariate functions W(k) and C(k).  The nodal functions are given by

!  C(k)(x,y) = A(1,k)*(x-X(k))**3 +
!              A(2,k)*(x-X(k))**2*(y-Y(k)) +
!              A(3,k)*(x-X(k))*(y-Y(k))**2 +
!              A(4,k)*(y-Y(k))**3 + A(5,k)*(x-X(k))**2 +
!              A(6,k)*(x-X(k))*(y-Y(k)) + A(7,k)*(y-Y(k))**2
!              + A(8,k)*(x-X(k)) + A(9,k)*(y-Y(k)) + F(k) .

! Thus, C(k) is a cubic function which interpolates the data value at node k.
! Its coefficients A(,k) are obtained by a weighted least squares fit to the
! closest NC data points with weights similar to W(k).
! Note that the radius of influence for the least squares fit is fixed for
! each k, but varies with k.

! The weights are taken to be

!   W(k)(x,y) = ( (R(k)-D(k))+ / R(k)*D(k) )**3 ,

! where (R(k)-D(k))+ = 0 if R(k) < D(k), and D(k)(x,y) is the Euclidean
! distance between (x,y) and (X(k),Y(k)).  The radius of influence R(k) varies
! with k and is chosen so that NW nodes are within the radius.
! Note that W(k) is not defined at node (X(k),Y(k)), but C(x,y) has limit F(k)
! as (x,y) approaches (X(k),Y(k)).

! On input:

!       N = Number of nodes and data values.  N >= 10.

!       X,Y = Arrays of length N containing the Cartesian coordinates of
!             the nodes.

!       F = Array of length N containing the data values in one-to-one
!           correspondence with the nodes.

!       NC = Number of data points to be used in the least squares fit for
!            coefficients defining the nodal functions C(k).  Values found to
!            be optimal for test data sets ranged from 11 to 25.
!            A recommended value for general data sets is NC = 17.
!            For nodes lying on (or close to) a rectangular grid, the
!            recommended value is NC = 11.  In any case, NC must be in the
!            range 9 to Min(40,N-1).

!       NW = Number of nodes within (and defining) the radii of influence R(k)
!            which enter into the weights W(k).
!            For N sufficiently large, a recommended value is NW = 30.
!            In general, NW should be about 1.5*NC.  1 <= NW <= Min(40,N-1).

!       NR = Number of rows and columns in the cell grid defined in Subroutine
!            STORE2.  A rectangle containing the nodes is partitioned into
!            cells in order to increase search efficiency.  NR = Sqrt(N/3)
!            is recommended.  NR >= 1.

! The above parameters are not altered by this routine.

!       LCELL = Array of length >= NR**2.

!       LNEXT = Array of length >= N.

!       RW = Array of length >= N.

!       A = Array of length >= 9N.

! On output:

!       LCELL = NR by NR array of nodal indexes associated with cells.
!               Refer to Subroutine STORE2.

!       LNEXT = Array of length N containing next-node indexes.
!               Refer to Subroutine STORE2.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell dimensions.
!                         Refer to Subroutine STORE2.

!       RMAX = Largest element in RW -- maximum radius R(k).

!       RW = Array containing the the radii R(k) which enter into the weights
!            W(k).

!       A = 9 by N array containing the coefficients for cubic nodal
!           function C(k) in column k.

!   Note that the output parameters described above are not defined unless
!   IER = 0.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NC, NW, or NR is outside its valid range.
!             IER = 2 if duplicate nodes were encountered.
!             IER = 3 if all nodes are collinear.

! Modules required by CSHEP2:  GETNP2, GIVENS, ROTATE, SETUP2, STORE2

! Intrinsic functions called by CSHEP2:  ABS, DBLE, MAX, MIN, SQRT

!***********************************************************


INTEGER, PARAMETER :: lmx=40
INTEGER   :: i, ierr, ip1, irm1, irow, j, jp1, k, lmax,  &
             lnp, neq, nn, nnc, nnr, nnw, np, npts(lmx), ncwmax
REAL (dp) :: b(10,10), c, ddx, ddy, dmin, fk, rc, rs, rsmx, rsold, &
             rws, s, sf, sfc, sfs, stf, sum, t, xk, xmn, yk, ymn

REAL (dp), PARAMETER :: dtol = 0.01_dp, rtol = 1.0e-05_dp

! Local parameters:

! B =          Transpose of the augmented regression matrix
! C =          First component of the plane rotation used to zero the lower
!                triangle of B**T -- computed by Subroutine GIVENS
! DDX,DDY =    Local variables for DX and DY
! DMIN =       Minimum of the magnitudes of the diagonal elements of the
!                regression matrix after zeros are introduced below the
!                diagonal
! DTOL =       Tolerance for detecting an ill-conditioned system.
!                The system is accepted when DMIN*RC >= DTOL.
! FK =         Data value at mode K -- F(K)
! I =          Index for A, B, and NPTS
! IERR =       Error flag for the call to Subroutine STORE2
! IP1 =        I+1
! IRM1 =       IROW-1
! IROW =       Row index for B
! J =          Index for A and B
! JP1 =        J+1
! K =          Nodal function index and column index for A
! LMAX =       Maximum number of NPTS elements
! LMX =        Maximum value of LMAX
! LNP =        Current length of NPTS
! NEQ =        Number of equations in the least squares fit
! NN,NNC,NNR = Local copies of N, NC, and NR
! NNW =        Local copy of NW
! NP =         NPTS element
! NPTS =       Array containing the indexes of a sequence of nodes to be used
!                in the least squares fit or to compute RW.
!                The nodes are ordered by distance from K, and the last
!                element (usually indexed by LNP) is used only to determine RC,
!                or RW(K) if NW > NC.
! NCWMAX =     Max(NC,NW)
! RC =         Radius of influence which enters into the weights for C(K)
!                (see Subroutine SETUP2)
! RS =         Squared distance between K and NPTS(LNP) -- used to compute
!                RC and RW(K)
! RSMX =       Maximum squared RW element encountered
! RSOLD =      Squared distance between K and NPTS(LNP-1) -- used to compute a
!                relative change in RS between succeeding NPTS elements
! RTOL =       Tolerance for detecting a sufficiently large relative change
!                in RS.  If the change is not greater than RTOL, the nodes are
!                treated as being the same distance from K
! RWS =        Current squared value of RW(K)
! S =          Second component of the plane rotation determined by subroutine
!                GIVENS
! SF =         Scale factor for the linear terms (columns 8 and 9) in the least
!                squares fit -- inverse of the root-mean-square distance
!                between K and the nodes (other than K) in the least squares fit
! SFS =        Scale factor for the quadratic terms (columns
!                5, 6, and 7) in the least squares fit -- SF*SF
! SFC =        Scale factor for the cubic terms (first 4 columns) in the least
!                squares fit -- SF**3
! STF =        Marquardt stabilization factor used to damp out the first 4
!                solution components (third partials of the cubic) when the
!                system is ill-conditioned.  As STF increases, the fitting
!                function approaches a quadratic polynomial.
! SUM =        Sum of squared Euclidean distances between node K and the nodes
!                used in the least squares fit (unless additional nodes are
!                added for stability)
! T =          Temporary variable for accumulating a scalar product in the
!                back solve
! XK,YK =      Coordinates of node K -- X(K), Y(K)
! XMN,YMN =    Local variables for XMIN and YMIN

nn = n
nnc = nc
nnw = nw
nnr = nr
ncwmax = MAX(nnc,nnw)
lmax = MIN(lmx,nn-1)
IF (nnc < 9 .OR. nnw < 1 .OR. ncwmax > lmax .OR. nnr < 1) GO TO 21

! Create the cell data structure, and initialize RSMX.

CALL store2 (nn, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, ierr)
IF (ierr /= 0) GO TO 23
rsmx = 0.

! Outer loop on node K:

DO  k = 1,nn
  xk = x(k)
  yk = y(k)
  fk = f(k)
  
! Mark node K to exclude it from the search for nearest neighbors.
  
  lnext(k) = -lnext(k)
  
! Initialize for loop on NPTS.
  
  rs = 0.
  sum = 0.
  rws = 0.
  rc = 0.
  lnp = 0
  
! Compute NPTS, LNP, RWS, NEQ, RC, and SFS.
  
  1   sum = sum + rs
  IF (lnp == lmax) GO TO 2
  lnp = lnp + 1
  rsold = rs
  CALL getnp2 (xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, np, rs)
  IF (rs == 0.) GO TO 22
  npts(lnp) = np
  IF ( (rs-rsold)/rs < rtol ) GO TO 1
  IF (rws == 0. .AND. lnp > nnw) rws = rs
  IF (rc == 0. .AND. lnp > nnc) THEN
    
!   RC = 0 (not yet computed) and LNP > NC.  RC = Sqrt(RS) is sufficiently
!     large to (strictly) include NC nodes.
!     The least squares fit will include NEQ = LNP - 1 equations for
!     9 <= NC <= NEQ < LMAX <= N-1.
    
    neq = lnp - 1
    rc = SQRT(rs)
    sfs = DBLE(neq)/sum
  END IF
  
!   Bottom of loop -- test for termination.
  
  IF (lnp > ncwmax) GO TO 3
  GO TO 1
  
! All LMAX nodes are included in NPTS.  RWS and/or RC**2 is (arbitrarily)
!   taken to be 10 percent larger than the distance RS to the last node
!   included.
  
  2   IF (rws == 0.) rws = 1.1*rs
  IF (rc == 0.) THEN
    neq = lmax
    rc = SQRT(1.1*rs)
    sfs = DBLE(neq)/sum
  END IF
  
! Store RW(K), update RSMX if necessary, and compute SF and SFC.
  
  3   rw(k) = SQRT(rws)
  IF (rws > rsmx) rsmx = rws
  sf = SQRT(sfs)
  sfc = sf*sfs
  
! A Q-R decomposition is used to solve the least squares system.
!   The transpose of the augmented regression matrix is stored in B with
!   columns (rows of B) defined as follows:  1-4 are the cubic terms,
!   5-7 are the quadratic terms, 8 and 9 are the linear terms,
!   and the last column is the right hand side.
  
! Set up the equations and zero out the lower triangle with Givens rotations.
  
  i = 0
  4 i = i + 1
  np = npts(i)
  irow = MIN(i,10)
  CALL setup2 (xk, yk, fk, x(np), y(np), f(np), sf, sfs,  sfc, rc, b(1,irow))
  IF (i == 1) GO TO 4
  irm1 = irow-1
  DO  j = 1,irm1
    jp1 = j + 1
    CALL givens (b(j,j), b(j,irow), c, s)
    CALL rotate (10-j, c, s, b(jp1:,j), b(jp1:, irow))
  END DO
  IF (i < neq) GO TO 4
  
! Test the system for ill-conditioning.
  
  dmin = MIN( ABS(b(1,1)), ABS(b(2,2)), ABS(b(3,3)), ABS(b(4,4)), ABS(b(5,5)), &
              ABS(b(6,6)), ABS(b(7,7)), ABS(b(8,8)), ABS(b(9,9)) )
  IF (dmin*rc >= dtol) GO TO 11
  IF (neq == lmax) GO TO 7
  
! Increase RC and add another equation to the system to improve the
!   conditioning.  The number of NPTS elements is also increased if necessary.
  
  6   rsold = rs
  neq = neq + 1
  IF (neq == lmax) THEN
    rc = SQRT(1.1*rs)
    GO TO 4
  END IF
  IF (neq < lnp) THEN
    
!   NEQ < LNP.
    
    np = npts(neq+1)
    rs = (x(np)-xk)**2 + (y(np)-yk)**2
    IF ( (rs-rsold)/rs < rtol ) GO TO 6
    rc = SQRT(rs)
    GO TO 4
  END IF
  
!   NEQ = LNP.  Add an element to NPTS.
  
  lnp = lnp + 1
  CALL getnp2 (xk, yk, x, y, nnr, lcell, lnext, xmn, ymn, ddx, ddy, np, rs)
  IF (np == 0) GO TO 22
  npts(lnp) = np
  IF ( (rs-rsold)/rs < rtol ) GO TO 6
  rc = SQRT(rs)
  GO TO 4
  
! Stabilize the system by damping third partials -- add
!   multiples of the first four unit vectors to the first
!   four equations.
  
  7 stf = 1.0/rc
  DO  i = 1,4
    b(i,10) = stf
    ip1 = i + 1
    DO  j = ip1,10
      b(j,10) = 0.
    END DO
    DO  j = i,9
      jp1 = j + 1
      CALL givens (b(j,j), b(j,10), c, s)
      CALL rotate (10-j, c, s, b(jp1:,j), b(jp1:,10))
    END DO
  END DO
  
! Test the damped system for ill-conditioning.
  
  dmin = MIN( ABS(b(5,5)), ABS(b(6,6)), ABS(b(7,7)), ABS(b(8,8)), ABS(b(9,9)) )
  IF (dmin*rc < dtol) GO TO 23
  
! Solve the 9 by 9 triangular system for the coefficients.
  
  11 DO  i = 9,1,-1
    t = 0.
    IF (i /= 9) THEN
      ip1 = i + 1
      DO  j = ip1,9
        t = t + b(j,i)*a(j,k)
      END DO
    END IF
    a(i,k) = (b(10,i)-t)/b(i,i)
  END DO
  
! Scale the coefficients to adjust for the column scaling.
  
  a(1:4,k) = a(1:4,k)*sfc
  a(5,k) = a(5,k)*sfs
  a(6,k) = a(6,k)*sfs
  a(7,k) = a(7,k)*sfs
  a(8,k) = a(8,k)*sf
  a(9,k) = a(9,k)*sf
  
! Unmark K and the elements of NPTS.
  
  lnext(k) = -lnext(k)
  DO  i = 1,lnp
    np = npts(i)
    lnext(np) = -lnext(np)
  END DO
END DO

! No errors encountered.

xmin = xmn
ymin = ymn
dx = ddx
dy = ddy
rmax = SQRT(rsmx)
ier = 0
RETURN

! N, NC, NW, or NR is outside its valid range.

21 ier = 1
RETURN

! Duplicate nodes were encountered by GETNP2.

22 ier = 2
RETURN

! No unique solution due to collinear nodes.

23 xmin = xmn
ymin = ymn
dx = ddx
dy = ddy
ier = 3
RETURN
END SUBROUTINE cshep2



FUNCTION cs2val (px, py, n, x, y, f, nr, lcell, lnext, xmin, ymin, dx, dy, &
                 rmax, rw, a) RESULT(fn_val)

REAL (dp), INTENT(IN)      :: px
REAL (dp), INTENT(IN)      :: py
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN)      :: y(:)
REAL (dp), INTENT(IN)      :: f(:)
INTEGER, INTENT(IN)        :: nr
INTEGER, INTENT(IN)        :: lcell(:,:)    ! lcell(nr,nr)
INTEGER, INTENT(IN)        :: lnext(:)
REAL (dp), INTENT(IN)      :: xmin
REAL (dp), INTENT(IN)      :: ymin
REAL (dp), INTENT(IN)      :: dx
REAL (dp), INTENT(IN)      :: dy
REAL (dp), INTENT(IN)      :: rmax
REAL (dp), INTENT(IN)      :: rw(:)
REAL (dp), INTENT(IN)      :: a(:,:)     ! a(9,n)
REAL (dp)                  :: fn_val

!***********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   This function returns the value C(PX,PY), where C is the weighted sum of
! cubic nodal functions defined in Subroutine CSHEP2.  CS2GRD may be called to
! compute a gradient of C along with the value, and/or to test for errors.
! CS2HES may be called to compute a value, first partial derivatives, and
! second partial derivatives at a point.

! On input:

!       PX,PY = Cartesian coordinates of the point P at
!               which C is to be evaluated.

!       N = Number of nodes and data values defining C.
!           N >= 10.

!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.

!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR >= 1.

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.

!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell dimensions.
!                         DX and DY must be positive.
!                         Refer to Subroutine STORE2.

!       RMAX = Largest element in RW -- maximum radius R(k).

!       RW = Array containing the the radii R(k) which enter
!            into the weights W(k) defining C.

!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.

!   Input parameters are not altered by this function.  The parameters other
! than PX and PY should be input unaltered from their values on output from
! CSHEP2.  This function should not be called if a nonzero error flag was
! returned by CSHEP2.

! On output:

!       CS2VAL = Function value C(PX,PY) unless N, NR, DX, DY, or RMAX is
!                invalid, in which case no value is returned.

! Modules required by CS2VAL:  NONE

! Intrinsic functions called by CS2VAL:  INT, SQRT

!***********************************************************

INTEGER   :: i, imax, imin, j, jmax, jmin, k, kp
REAL (dp) :: d, delx, dely, r, sw, swc, w, xp, yp

! Local parameters:

! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells intersected by a disk of
!               radius RMAX centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells intersected by a disk
!               of radius RMAX centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! XP,YP =     Local copies of PX and PY -- coordinates of P

xp = px
yp = py
IF (n < 10 .OR. nr < 1 .OR. dx <= 0. .OR. dy <= 0. .OR. rmax < 0.) RETURN

! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining the range of the
!   search for nodes whose radii include P.  The cells which must be searched
!   are those intersected by (or contained in) a circle of radius RMAX
!   centered at P.

imin = INT((xp-xmin-rmax)/dx) + 1
imax = INT((xp-xmin+rmax)/dx) + 1
IF (imin < 1) imin = 1
IF (imax > nr) imax = nr
jmin = INT((yp-ymin-rmax)/dy) + 1
jmax = INT((yp-ymin+rmax)/dy) + 1
IF (jmin < 1) jmin = 1
IF (jmax > nr) jmax = nr

! The following is a test for no cells within the circle of radius RMAX.

IF (imin > imax .OR. jmin > jmax) GO TO 6

! Accumulate weight values in SW and weighted nodal function
!   values in SWC.  The weights are W(K) = ((R-D)+/(R*D))**3
!   for R = RW(K) and D = distance between P and node K.

sw = 0.
swc = 0.

! Outer loop on cells (I,J).

DO  j = jmin,jmax
  DO  i = imin,imax
    k = lcell(i,j)
    IF (k == 0) CYCLE
    
! Inner loop on nodes K.
    
    1 delx = xp - x(k)
    dely = yp - y(k)
    d = SQRT(delx*delx + dely*dely)
    r = rw(k)
    IF (d >= r) GO TO 2
    IF (d == 0.) GO TO 5
    w = (1.0/d - 1.0/r)**3
    sw = sw + w
    swc = swc + w*( ( (a(1,k)*delx + a(2,k)*dely + a(5,k))*delx +  &
                    (a(3,k)*dely + a(6,k))*dely + a(8,k) )*delx +  &
                    ( (a(4,k)*dely + a(7,k))*dely + a(9,k) )*dely + f(k) )
    
! Bottom of loop on nodes in cell (I,J).
    
    2 kp = k
    k = lnext(kp)
    IF (k /= kp) GO TO 1
  END DO
END DO

! SW = 0 iff P is not within the radius R(K) for any node K.

IF (sw == 0.) GO TO 6
fn_val = swc/sw
RETURN

! (PX,PY) = (X(K),Y(K)).

5 fn_val = f(k)
RETURN

! All weights are 0 at P.

6 fn_val = 0.
RETURN
END FUNCTION cs2val



SUBROUTINE cs2grd (px, py, n, x, y, f, nr, lcell, lnext, xmin,  &
                   ymin, dx, dy, rmax, rw, a, c, cx, cy, ier)

REAL (dp), INTENT(IN)       :: px
REAL (dp), INTENT(IN)       :: py
INTEGER, INTENT(IN)         :: n
REAL (dp), INTENT(IN)       :: x(:)
REAL (dp), INTENT(IN)       :: y(:)
REAL (dp), INTENT(IN)       :: f(:)
INTEGER, INTENT(IN)         :: nr
INTEGER, INTENT(IN)         :: lcell(:,:)    ! lcell(nr,nr)
INTEGER, INTENT(IN)         :: lnext(:)
REAL (dp), INTENT(IN)       :: xmin
REAL (dp), INTENT(IN)       :: ymin
REAL (dp), INTENT(IN)       :: dx
REAL (dp), INTENT(IN OUT)   :: dy
REAL (dp), INTENT(IN OUT)   :: rmax
REAL (dp), INTENT(IN)       :: rw(:)
REAL (dp), INTENT(IN)       :: a(:,:)       ! a(9,n)
REAL (dp), INTENT(OUT)      :: c
REAL (dp), INTENT(OUT)      :: cx
REAL (dp), INTENT(OUT)      :: cy
INTEGER, INTENT(OUT)        :: ier

!***********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   This subroutine computes the value and gradient at P = (PX,PY) of the
! interpolatory function C defined in Subroutine CSHEP2.
! C is a weighted sum of cubic nodal functions.

! On input:

!       PX,PY = Cartesian coordinates of the point P at which C and
!               its partial derivatives are to be evaluated.

!       N = Number of nodes and data values defining C.
!           N >= 10.

!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.

!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR >= 1.

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.

!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell dimensions.
!                         DX and DY must be positive.
!                         Refer to Subroutine STORE2.

!       RMAX = Largest element in RW -- maximum radius R(k).

!       RW = Array of length N containing the the radii R(k)
!            which enter into the weights W(k) defining C.

!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.

!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input unaltered from their
! values on output from CSHEP2.  This subroutine should not be called if a
! nonzero error flag was returned by CSHEP2.

! On output:

!       C = Value of C at (PX,PY) unless IER .EQ. 1, in
!           which case no values are returned.

!       CX,CY = First partial derivatives of C at (PX,PY) unless IER .EQ. 1.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but
!                     (PX,PY) is not within the radius R(k)
!                     for any node k (and thus C=CX=CY=0).

! Modules required by CS2GRD:  None

! Intrinsic functions called by CS2GRD:  INT, SQRT

!***********************************************************

INTEGER   :: i, imax, imin, j, jmax, jmin, k, kp
REAL (dp) :: ck, ckx, cky, d, delx, dely, r, sw,  &
             swc, swcx, swcy, sws, swx, swy, t, w, wx, wy, xp, yp

! Local parameters:

! CK =        Value of cubic nodal function C(K) at P
! CKX,CKY =   Partial derivatives of C(K) with respect to X and Y, respectively
! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells intersected by a disk
!               of radius RMAX centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells intersected by a disk
!               of radius RMAX centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! SWCX,SWCY = Partial derivatives of SWC with respect to X and Y, respectively
! SWS =       SW**2
! SWX,SWY =   Partial derivatives of SW with respect to X and Y, respectively
! T =         Temporary variable
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! WX,WY =     Partial derivatives of W with respect to X and Y, respectively
! XP,YP =     Local copies of PX and PY -- coordinates of P

xp = px
yp = py
IF (n < 10 .OR. nr < 1 .OR. dx <= 0. .OR. dy <= 0. .OR. rmax < 0.) GO TO 6

! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining the range of the
!   search for nodes whose radii include P.
!   The cells which must be searched are those intersected by (or contained in)
!   a circle of radius RMAX centered at P.

imin = INT((xp-xmin-rmax)/dx) + 1
imax = INT((xp-xmin+rmax)/dx) + 1
IF (imin < 1) imin = 1
IF (imax > nr) imax = nr
jmin = INT((yp-ymin-rmax)/dy) + 1
jmax = INT((yp-ymin+rmax)/dy) + 1
IF (jmin < 1) jmin = 1
IF (jmax > nr) jmax = nr

! The following is a test for no cells within the circle of radius RMAX.

IF (imin > imax .OR. jmin > jmax) GO TO 7

! C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is from K = 1 to N,
!     C(K) is the cubic nodal function value,
!     and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and distance D(K).  Thus

!        CX = (SWCX*SW - SWC*SWX)/SW**2  and
!        CY = (SWCY*SW - SWC*SWY)/SW**2

!   where SWCX and SWX are partial derivatives with respect to X of SWC and SW,
!   respectively.  SWCY and SWY are defined similarly.

sw = 0.
swx = 0.
swy = 0.
swc = 0.
swcx = 0.
swcy = 0.

! Outer loop on cells (I,J).

DO  j = jmin,jmax
  DO  i = imin,imax
    k = lcell(i,j)
    IF (k == 0) CYCLE
    
! Inner loop on nodes K.
    
    1 delx = xp - x(k)
    dely = yp - y(k)
    d = SQRT(delx*delx + dely*dely)
    r = rw(k)
    IF (d >= r) GO TO 2
    IF (d == 0.) GO TO 5
    t = (1.0/d - 1.0/r)
    w = t**3
    t = -3.0*t*t/(d**3)
    wx = delx*t
    wy = dely*t
    t = a(2,k)*delx + a(3,k)*dely + a(6,k)
    cky = ( 3.0*a(4,k)*dely + a(3,k)*delx + 2.0*a(7,k) )*dely + t*delx + a(9,k)
    t = t*dely + a(8,k)
    ckx = ( 3.0*a(1,k)*delx + a(2,k)*dely + 2.0*a(5,k) )*delx + t
    ck = ( (a(1,k)*delx + a(5,k))*delx + t )*delx +  &
         ( (a(4,k)*dely + a(7,k))*dely + a(9,k) )*dely + f(k)
    sw = sw + w
    swx = swx + wx
    swy = swy + wy
    swc = swc + w*ck
    swcx = swcx + wx*ck + w*ckx
    swcy = swcy + wy*ck + w*cky
    
! Bottom of loop on nodes in cell (I,J).
    
    2 kp = k
    k = lnext(kp)
    IF (k /= kp) GO TO 1
  END DO
END DO

! SW = 0 iff P is not within the radius R(K) for any node K.

IF (sw == 0.) GO TO 7
c = swc/sw
sws = sw*sw
cx = (swcx*sw - swc*swx)/sws
cy = (swcy*sw - swc*swy)/sws
ier = 0
RETURN

! (PX,PY) = (X(K),Y(K)).

5 c = f(k)
cx = a(8,k)
cy = a(9,k)
ier = 0
RETURN

! Invalid input parameter.

6 ier = 1
RETURN

! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D >= RW(K) for all K.

7 c = 0.
cx = 0.
cy = 0.
ier = 2
RETURN
END SUBROUTINE cs2grd



SUBROUTINE cs2hes (px, py, n, x, y, f, nr, lcell, lnext, xmin,  &
                   ymin, dx, dy, rmax, rw, a, c, cx, cy, cxx, cxy, cyy, ier)

REAL (dp), INTENT(IN)      :: px
REAL (dp), INTENT(IN)      :: py
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: x(:)
REAL (dp), INTENT(IN)      :: y(:)
REAL (dp), INTENT(IN)      :: f(:)
INTEGER, INTENT(IN)        :: nr
INTEGER, INTENT(IN)        :: lcell(:,:)    ! lcell(nr,nr)
INTEGER, INTENT(IN)        :: lnext(:)
REAL (dp), INTENT(IN)      :: xmin
REAL (dp), INTENT(IN)      :: ymin
REAL (dp), INTENT(IN)      :: dx
REAL (dp), INTENT(IN OUT)  :: dy
REAL (dp), INTENT(IN OUT)  :: rmax
REAL (dp), INTENT(IN)      :: rw(:)
REAL (dp), INTENT(IN)      :: a(:,:)      ! a(9,n)
REAL (dp), INTENT(OUT)     :: c
REAL (dp), INTENT(OUT)     :: cx
REAL (dp), INTENT(OUT)     :: cy
REAL (dp), INTENT(OUT)     :: cxx
REAL (dp), INTENT(OUT)     :: cxy
REAL (dp), INTENT(OUT)     :: cyy
INTEGER, INTENT(OUT)       :: ier

!***********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   This subroutine computes the value, gradient, and Hessian at P = (PX,PY)
! of the interpolatory function C defined in Subroutine CSHEP2.
! C is a weighted sum of cubic nodal functions.

! On input:

!       PX,PY = Cartesian coordinates of the point P at which C and its
!               partial derivatives are to be evaluated.

!       N = Number of nodes and data values defining C.
!           N >= 10.

!       X,Y,F = Arrays of length N containing the nodes and
!               data values interpolated by C.

!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR >= 1.

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.

!       LNEXT = Array of length N containing next-node
!               indexes.  Refer to Subroutine STORE2.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell dimensions.
!                         DX and DY must be positive.
!                         Refer to Subroutine STORE2.

!       RMAX = Largest element in RW -- maximum radius R(k).

!       RW = Array of length N containing the the radii R(k)
!            which enter into the weights W(k) defining C.

!       A = 9 by N array containing the coefficients for
!           cubic nodal function C(k) in column k.

!   Input parameters are not altered by this subroutine.
! The parameters other than PX and PY should be input unaltered from their
! values on output from CSHEP2.  This subroutine should not be called if a
! nonzero error flag was returned by CSHEP2.

! On output:

!       C = Value of C at (PX,PY) unless IER .EQ. 1, in
!           which case no values are returned.

!       CX,CY = First partial derivatives of C at (PX,PY)
!               unless IER .EQ. 1.

!       CXX,CXY,CYY = Second partial derivatives of C at
!                     (PX,PY) unless IER .EQ. 1.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N, NR, DX, DY or RMAX is invalid.
!             IER = 2 if no errors were encountered but (PX,PY) is not within
!                     the radius R(k) for any node k (and thus C = 0).

! Modules required by CS2HES:  None

! Intrinsic functions called by CS2HES:  INT, SQRT

!***********************************************************

INTEGER   :: i, imax, imin, j, jmax, jmin, k, kp
REAL (dp) :: ck, ckx, ckxx, ckxy, cky, ckyy, d, delx, dely, dxsq, dysq, r, &
             sw, swc, swcx, swcxx, swcxy, swcy, swcyy, sws, swx, swxx, swxy,  &
             swy, swyy, t1, t2, t3, t4, w, wx, wxx, wxy, wy, wyy, xp, yp

! Local parameters:

! CK =        Value of cubic nodal function C(K) at P
! CKX,CKY =   Partial derivatives of C(K) with respect to X and Y, respectively
! CKXX,CKXY,CKYY = Second partial derivatives of CK
! D =         Distance between P and node K
! DELX =      XP - X(K)
! DELY =      YP - Y(K)
! DXSQ,DYSQ = DELX**2, DELY**2
! I =         Cell row index in the range IMIN to IMAX
! IMIN,IMAX = Range of cell row indexes of the cells intersected by a disk
!               of radius RMAX centered at P
! J =         Cell column index in the range JMIN to JMAX
! JMIN,JMAX = Range of cell column indexes of the cells intersected by a disk
!               of radius RMAX centered at P
! K =         Index of a node in cell (I,J)
! KP =        Previous value of K in the sequence of nodes in cell (I,J)
! R =         Radius of influence for node K
! SW =        Sum of weights W(K)
! SWC =       Sum of weighted nodal function values at P
! SWCX,SWCY = Partial derivatives of SWC with respect to X and Y, respectively
! SWCXX,SWCXY,SWCYY = Second partial derivatives of SWC
! SWS =       SW**2
! SWX,SWY =   Partial derivatives of SW with respect to X and Y, respectively
! SWXX,SWXY,SWYY = Second partial derivatives of SW
! T1,T2,T3,T4 = Temporary variables
! W =         Weight W(K) value at P:  ((R-D)+/(R*D))**3,
!               where (R-D)+ = 0 if R < D
! WX,WY =     Partial derivatives of W with respect to X and Y, respectively
! WXX,WXY,WYY = Second partial derivatives of W
! XP,YP =     Local copies of PX and PY -- coordinates of P

xp = px
yp = py
IF (n < 10 .OR. nr < 1 .OR. dx <= 0. .OR. dy <= 0. .OR. rmax < 0.) GO TO 6

! Set IMIN, IMAX, JMIN, and JMAX to cell indexes defining the range of the
!   search for nodes whose radii include P.
!   The cells which must be searched are those intersected by (or contained in)
!   a circle of radius RMAX centered at P.

imin = INT((xp-xmin-rmax)/dx) + 1
imax = INT((xp-xmin+rmax)/dx) + 1
IF (imin < 1) imin = 1
IF (imax > nr) imax = nr
jmin = INT((yp-ymin-rmax)/dy) + 1
jmax = INT((yp-ymin+rmax)/dy) + 1
IF (jmin < 1) jmin = 1
IF (jmax > nr) jmax = nr

! The following is a test for no cells within the circle of radius RMAX.

IF (imin > imax .OR. jmin > jmax) GO TO 7

! C = SWC/SW = Sum(W(K)*C(K))/Sum(W(K)), where the sum is from K = 1 to N,
!     C(K) is the cubic nodal function value,
!     and W(K) = ((R-D)+/(R*D))**3 for radius R(K) and distance D(K).  Thus

!        CX = (SWCX*SW - SWC*SWX)/SW**2  and
!        CY = (SWCY*SW - SWC*SWY)/SW**2

!   where SWCX and SWX are partial derivatives with respect to x of SWC and SW,
!   respectively.  SWCY and SWY are defined similarly.  The second partials are

!        CXX = ( SW*(SWCXX -    2*SWX*CX) - SWC*SWXX )/SW**2
!        CXY = ( SW*(SWCXY-SWX*CY-SWY*CX) - SWC*SWXY )/SW**2
!        CYY = ( SW*(SWCYY -    2*SWY*CY) - SWC*SWYY )/SW**2

!   where SWCXX and SWXX are second partials with respect to x, SWCXY and SWXY
!   are mixed partials, and SWCYY and SWYY are second partials with respect
!   to y.

sw = 0.
swx = 0.
swy = 0.
swxx = 0.
swxy = 0.
swyy = 0.
swc = 0.
swcx = 0.
swcy = 0.
swcxx = 0.
swcxy = 0.
swcyy = 0.

! Outer loop on cells (I,J).

DO  j = jmin,jmax
  DO  i = imin,imax
    k = lcell(i,j)
    IF (k == 0) CYCLE
    
! Inner loop on nodes K.
    
    1 delx = xp - x(k)
    dely = yp - y(k)
    dxsq = delx*delx
    dysq = dely*dely
    d = SQRT(dxsq + dysq)
    r = rw(k)
    IF (d >= r) GO TO 2
    IF (d == 0.) GO TO 5
    t1 = (1.0/d - 1.0/r)
    w = t1**3
    t2 = -3.0*t1*t1/(d**3)
    wx = delx*t2
    wy = dely*t2
    t1 = 3.0*t1*(2.0+3.0*d*t1)/(d**6)
    wxx = t1*dxsq + t2
    wxy = t1*delx*dely
    wyy = t1*dysq + t2
    t1 = a(1,k)*delx + a(2,k)*dely + a(5,k)
    t2 = t1 + t1 + a(1,k)*delx
    t3 = a(4,k)*dely + a(3,k)*delx + a(7,k)
    t4 = t3 + t3 + a(4,k)*dely
    ck = (t1*delx + a(6,k)*dely + a(8,k))*delx +  &
         (t3*dely + a(9,k))*dely + f(k)
    ckx = t2*delx + (a(3,k)*dely+a(6,k))*dely + a(8,k)
    cky = t4*dely + (a(2,k)*delx+a(6,k))*delx + a(9,k)
    ckxx = t2 + 3.0*a(1,k)*delx
    ckxy = 2.0*(a(2,k)*delx + a(3,k)*dely) + a(6,k)
    ckyy = t4 + 3.0*a(4,k)*dely
    sw = sw + w
    swx = swx + wx
    swy = swy + wy
    swxx = swxx + wxx
    swxy = swxy + wxy
    swyy = swyy + wyy
    swc = swc + w*ck
    swcx = swcx + wx*ck + w*ckx
    swcy = swcy + wy*ck + w*cky
    swcxx = swcxx + w*ckxx + 2.0*wx*ckx + ck*wxx
    swcxy = swcxy + w*ckxy + wx*cky + wy*ckx + ck*wxy
    swcyy = swcyy + w*ckyy + 2.0*wy*cky + ck*wyy
    
! Bottom of loop on nodes in cell (I,J).
    
    2 kp = k
    k = lnext(kp)
    IF (k /= kp) GO TO 1
  END DO
END DO

! SW = 0 iff P is not within the radius R(K) for any node K.

IF (sw == 0.) GO TO 7
c = swc/sw
sws = sw*sw
cx = (swcx*sw - swc*swx)/sws
cy = (swcy*sw - swc*swy)/sws
cxx = (sw*(swcxx-2.0*swx*cx) - swc*swxx)/sws
cxy = (sw*(swcxy-swy*cx-swx*cy) - swc*swxy)/sws
cyy = (sw*(swcyy-2.0*swy*cy) - swc*swyy)/sws
ier = 0
RETURN

! (PX,PY) = (X(K),Y(K)).

5 c = f(k)
cx = a(8,k)
cy = a(9,k)
cxx = 2.0*a(5,k)
cxy = a(6,k)
cyy = 2.0*a(7,k)
ier = 0
RETURN

! Invalid input parameter.

6 ier = 1
RETURN

! No cells contain a point within RMAX of P, or
!   SW = 0 and thus D >= RW(K) for all K.

7 c = 0.
cx = 0.
cy = 0.
cxx = 0.
cxy = 0.
cyy = 0.
ier = 2
RETURN
END SUBROUTINE cs2hes



SUBROUTINE getnp2 (px, py, x, y, nr, lcell, lnext, xmin, ymin, dx, dy, np, dsq)

REAL (dp), INTENT(IN)    :: px
REAL (dp), INTENT(IN)    :: py
REAL (dp), INTENT(IN)    :: x(:)
REAL (dp), INTENT(IN)    :: y(:)
INTEGER, INTENT(IN)      :: nr
INTEGER, INTENT(IN)      :: lcell(:,:)    ! lcell(nr,nr)
INTEGER, INTENT(IN OUT)  :: lnext(:)
REAL (dp), INTENT(IN)    :: xmin
REAL (dp), INTENT(IN)    :: ymin
REAL (dp), INTENT(IN)    :: dx
REAL (dp), INTENT(IN)    :: dy
INTEGER, INTENT(OUT)     :: np
REAL (dp), INTENT(OUT)   :: dsq

!***********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   Given a set of N nodes and the data structure defined in
! Subroutine STORE2, this subroutine uses the cell method to
! find the closest unmarked node NP to a specified point P.
! NP is then marked by setting LNEXT(NP) to -LNEXT(NP).  (A
! node is marked if and only if the corresponding LNEXT ele-
! ment is negative.  The absolute values of LNEXT elements,
! however, must be preserved.)  Thus, the closest M nodes to
! P may be determined by a sequence of M calls to this rou-
! tine.  Note that if the nearest neighbor to node K is to
! be determined (PX = X(K) and PY = Y(K)), then K should be
! marked before the call to this routine.

!   The search is begun in the cell containing (or closest
! to) P and proceeds outward in rectangular layers until all
! cells which contain points within distance R of P have
! been searched, where R is the distance from P to the first
! unmarked node encountered (infinite if no unmarked nodes are present).

!   This code is essentially unaltered from the subroutine
! of the same name in QSHEP2D.

! On input:

!       PX,PY = Cartesian coordinates of the point P whose
!               nearest unmarked neighbor is to be found.

!       X,Y = Arrays of length N, for N >= 2, containing
!             the Cartesian coordinates of the nodes.

!       NR = Number of rows and columns in the cell grid.
!            Refer to Subroutine STORE2.  NR >= 1.

!       LCELL = NR by NR array of nodal indexes associated
!               with cells.  Refer to Subroutine STORE2.

!       LNEXT = Array of length N containing next-node
!               indexes (or their negatives).  Refer to
!               Subroutine STORE2.

!       XMIN,YMIN,DX,DY = Minimum nodal coordinates and cell
!                         dimensions.  DX and DY must be
!                         positive.  Refer to Subroutine
!                         STORE2.

!   Input parameters other than LNEXT are not altered by
! this routine.  With the exception of (PX,PY) and the signs
! of LNEXT elements, these parameters should be unaltered
! from their values on output from Subroutine STORE2.

! On output:

!       NP = Index (for X and Y) of the nearest unmarked
!            node to P, or 0 if all nodes are marked or NR
!            < 1 or DX <= 0 or DY <= 0.  LNEXT(NP)
!            < 0 IF NP .NE. 0.

!       DSQ = Squared Euclidean distance between P and node
!             NP, or 0 if NP = 0.

! Modules required by GETNP2:  None

! Intrinsic functions called by GETNP2:  ABS, INT, SQRT

!***********************************************************

INTEGER   :: i, i0, i1, i2, imax, imin, j, j0, j1, j2, jmax, jmin, l, lmin, ln
LOGICAL   :: first
REAL (dp) :: delx, dely, r, rsmin, rsq, xp, yp

! Local parameters:

! DELX,DELY =   PX-XMIN, PY-YMIN
! FIRST =       Logical variable with value TRUE iff the
!                 first unmarked node has yet to be encountered
! I,J =         Cell indexes in the range [I1,I2] X [J1,J2]
! I0,J0 =       Indexes of the cell containing or closest to P
! I1,I2,J1,J2 = Range of cell indexes defining the layer whose intersection
!                 with the range [IMIN,IMAX] X [JMIN,JMAX] is currently
!                 being searched
! IMIN,IMAX =   Cell row indexes defining the range of the search
! JMIN,JMAX =   Cell column indexes defining the range of the search
! L,LN =        Indexes of nodes in cell (I,J)
! LMIN =        Current candidate for NP
! R =           Distance from P to node LMIN
! RSMIN =       Squared distance from P to node LMIN
! RSQ =         Squared distance from P to node L
! XP,YP =       Local copy of PX,PY -- coordinates of P

xp = px
yp = py

! Test for invalid input parameters.

IF (nr < 1 .OR. dx <= 0. .OR. dy <= 0.) GO TO 9

! Initialize parameters.

first = .true.
imin = 1
imax = nr
jmin = 1
jmax = nr
delx = xp - xmin
dely = yp - ymin
i0 = INT(delx/dx) + 1
IF (i0 < 1) i0 = 1
IF (i0 > nr) i0 = nr
j0 = INT(dely/dy) + 1
IF (j0 < 1) j0 = 1
IF (j0 > nr) j0 = nr
i1 = i0
i2 = i0
j1 = j0
j2 = j0

! Outer loop on layers, inner loop on layer cells, excluding
!   those outside the range [IMIN,IMAX] X [JMIN,JMAX].

1 loop6:  DO  j = j1,j2
  IF (j > jmax) GO TO 7
  IF (j < jmin) CYCLE loop6
  DO  i = i1,i2
    IF (i > imax) CYCLE loop6
    IF (i < imin) CYCLE
    IF (j /= j1 .AND. j /= j2 .AND. i /= i1 .AND.  i /= i2) CYCLE
    
! Search cell (I,J) for unmarked nodes L.
    
    l = lcell(i,j)
    IF (l == 0) CYCLE
    
!   Loop on nodes in cell (I,J).
    
    2 ln = lnext(l)
    IF (ln < 0) GO TO 4
    
!   Node L is not marked.
    
    rsq = (x(l)-xp)**2 + (y(l)-yp)**2
    IF (.NOT. first) GO TO 3
    
!   Node L is the first unmarked neighbor of P encountered.
!     Initialize LMIN to the current candidate for NP, and
!     RSMIN to the squared distance from P to LMIN.  IMIN,
!     IMAX, JMIN, and JMAX are updated to define the smal-
!     lest rectangle containing a circle of radius R =
!     Sqrt(RSMIN) centered at P, and contained in [1,NR] X
!     [1,NR] (except that, if P is outside the rectangle
!     defined by the nodes, it is possible that IMIN > NR,
!     IMAX < 1, JMIN > NR, or JMAX < 1).  FIRST is reset to FALSE.
    
    lmin = l
    rsmin = rsq
    r = SQRT(rsmin)
    imin = INT((delx-r)/dx) + 1
    IF (imin < 1) imin = 1
    imax = INT((delx+r)/dx) + 1
    IF (imax > nr) imax = nr
    jmin = INT((dely-r)/dy) + 1
    IF (jmin < 1) jmin = 1
    jmax = INT((dely+r)/dy) + 1
    IF (jmax > nr) jmax = nr
    first = .false.
    GO TO 4
    
!   Test for node L closer than LMIN to P.
    
    3 IF (rsq >= rsmin) GO TO 4
    
!   Update LMIN and RSMIN.
    
    lmin = l
    rsmin = rsq
    
!   Test for termination of loop on nodes in cell (I,J).
    
    4 IF (ABS(ln) == l) CYCLE
    l = ABS(ln)
    GO TO 2
  END DO
END DO loop6

! Test for termination of loop on cell layers.

7 IF (i1 <= imin .AND. i2 >= imax .AND. j1 <= jmin .AND. j2 >= jmax) GO TO 8
i1 = i1 - 1
i2 = i2 + 1
j1 = j1 - 1
j2 = j2 + 1
GO TO 1

! Unless no unmarked nodes were encountered, LMIN is the closest unmarked node
! to P.

8 IF (first) GO TO 9
np = lmin
dsq = rsmin
lnext(lmin) = -lnext(lmin)
RETURN

! Error:  NR, DX, or DY is invalid or all nodes are marked.

9 np = 0
dsq = 0.
RETURN
END SUBROUTINE getnp2



SUBROUTINE givens ( a,b, c,s)

REAL (dp), INTENT(IN OUT)  :: a
REAL (dp), INTENT(IN OUT)  :: b
REAL (dp), INTENT(OUT)     :: c
REAL (dp), INTENT(OUT)     :: s

!***********************************************************

!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!   This subroutine constructs the Givens plane rotation,

!           ( C  S)
!       G = (     ) , where C*C + S*S = 1,
!           (-S  C)

! which zeros the second component of the vector (A,B)**T (transposed).
! Subroutine ROTATE may be called to apply the transformation to a
! 2 by N matrix.

!   This routine is identical to subroutine SROTG from the
! LINPACK BLAS (Basic Linear Algebra Subroutines).

! On input:

!       A,B = Components of the vector defining the rotation.
!             These are overwritten by values R and Z
!             (described below) which define C and S.

! On output:

!       A = Signed Euclidean norm R of the input vector:
!           R = +/-SQRT(A*A + B*B)

!       B = Value Z such that:
!             C = SQRT(1-Z*Z) and S=Z if ABS(Z) <= 1, and
!             C = 1/Z and S = SQRT(1-C*C) if ABS(Z) > 1.

!       C = +/-(A/R) or 1 if R = 0.

!       S = +/-(B/R) or 0 if R = 0.

! Modules required by GIVENS:  None

! Intrinsic functions called by GIVENS:  ABS, SQRT

!***********************************************************

REAL (dp) :: aa, bb, r, u, v

! Local parameters:

! AA,BB = Local copies of A and B
! R =     C*A + S*B = +/-SQRT(A*A+B*B)
! U,V =   Variables used to scale A and B for computing R

aa = a
bb = b
IF (ABS(aa) <= ABS(bb)) GO TO 1

! ABS(A) > ABS(B).

u = aa + aa
v = bb/u
r = SQRT(.25 + v*v) * u
c = aa/r
s = v * (c + c)

! Note that R has the sign of A, C > 0, and S has SIGN(A)*SIGN(B).

b = s
a = r
RETURN

! ABS(A) <= ABS(B).

1 IF (bb == 0.) GO TO 2
u = bb + bb
v = aa/u

! Store R in A.

a = SQRT(.25 + v*v) * u
s = bb/a
c = v * (s + s)

! Note that R has the sign of B, S > 0, and C has SIGN(A)*SIGN(B).

b = 1.
IF (c /= 0.) b = 1./c
RETURN

! A = B = 0.

2 c = 1.
s = 0.
RETURN
END SUBROUTINE givens



SUBROUTINE rotate (n,c,s, x,y )

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: c
REAL (dp), INTENT(IN)      :: s
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN OUT)  :: y(:)

!***********************************************************

!                                               From SRFPACK
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   09/01/88

!                                                ( C  S)
!   This subroutine applies the Givens rotation  (     )  to
!                                                (-S  C)
!                    (X(1) ... X(N))
! the 2 by N matrix  (             ) .
!                    (Y(1) ... Y(N))

!   This routine is identical to subroutine SROT from the
! LINPACK BLAS (Basic Linear Algebra Subroutines).

! On input:

!       N = Number of columns to be rotated.

!       C,S = Elements of the Givens rotation.  Refer to
!             subroutine GIVENS.

! The above parameters are not altered by this routine.

!       X,Y = Arrays of length >= N containing the components
!             of the vectors to be rotated.

! On output:

!       X,Y = Arrays containing the rotated vectors (not altered if N < 1).

! Modules required by ROTATE:  None

!***********************************************************

INTEGER   :: i
REAL (dp) :: xi, yi

DO  i = 1,n
  xi = x(i)
  yi = y(i)
  x(i) = c*xi + s*yi
  y(i) = -s*xi + c*yi
END DO
RETURN
END SUBROUTINE rotate



SUBROUTINE setup2 (xk, yk, zk, xi, yi, zi, s1, s2, s3, r, row)

REAL (dp), INTENT(IN)   :: xk
REAL (dp), INTENT(IN)   :: yk
REAL (dp), INTENT(IN)   :: zk
REAL (dp), INTENT(IN)   :: xi
REAL (dp), INTENT(IN)   :: yi
REAL (dp), INTENT(IN)   :: zi
REAL (dp), INTENT(IN)   :: s1
REAL (dp), INTENT(IN)   :: s2
REAL (dp), INTENT(IN)   :: s3
REAL (dp), INTENT(IN)   :: r
REAL (dp), INTENT(OUT)  :: row(10)

!***********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   02/03/97

!   This subroutine sets up the I-th row of an augmented regression matrix
! for a weighted least squares fit of a cubic function f(x,y) to a set of data
! values z, where f(XK,YK) = ZK.  The first four columns (cubic terms) are
! scaled by S3, the next three columns (quadratic terms) are scaled by S2,
! and the eighth and ninth columns (linear terms) are scaled by S1.

! On input:

!       XK,YK = Coordinates of node K.

!       ZK = Data value at node K to be interpolated by f.

!       XI,YI,ZI = Coordinates and data value at node I.

!       S1,S2,S3 = Scale factors.

!       R = Radius of influence about node K defining the
!           weight.

! The above parameters are not altered by this routine.

!       ROW = Array of length 10.

! On output:

!       ROW = Array containing a row of the augmented regression matrix.

! Modules required by SETUP2:  None

! Intrinsic function called by SETUP2:  SQRT

!***********************************************************

REAL (dp) :: d, dx, dxsq, dy, dysq, w, w1, w2, w3

! Local parameters:

! D =    Distance between nodes K and I
! DX =   XI - XK
! DXSQ = DX*DX
! DY =   YI - YK
! DYSQ = DY*DY
! I =    DO-loop index
! W =    Weight associated with the row:  (R-D)/(R*D)
!          (0 if D = 0 or D > R)
! W1 =   S1*W
! W2 =   S2*W
! W3 =   W3*W

dx = xi - xk
dy = yi - yk
dxsq = dx*dx
dysq = dy*dy
d = SQRT(dxsq + dysq)
IF (d <= 0. .OR. d >= r) GO TO 1
w = (r-d)/r/d
w1 = s1*w
w2 = s2*w
w3 = s3*w
row(1) = dxsq*dx*w3
row(2) = dxsq*dy*w3
row(3) = dx*dysq*w3
row(4) = dysq*dy*w3
row(5) = dxsq*w2
row(6) = dx*dy*w2
row(7) = dysq*w2
row(8) = dx*w1
row(9) = dy*w1
row(10) = (zi - zk)*w
RETURN

! Nodes K and I coincide or node I is outside of the radius of influence.
!   Set ROW to the zero vector.

1 row(1:10) = 0.

RETURN
END SUBROUTINE setup2



SUBROUTINE store2 (n, x, y, nr, lcell, lnext, xmin, ymin, dx, dy, ier)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: x(:)
REAL (dp), INTENT(IN)   :: y(:)
INTEGER, INTENT(IN)     :: nr
INTEGER, INTENT(OUT)    :: lcell(:,:)     ! lcell(nr,nr)
INTEGER, INTENT(OUT)    :: lnext(:)
REAL (dp), INTENT(OUT)  :: xmin
REAL (dp), INTENT(OUT)  :: ymin
REAL (dp), INTENT(OUT)  :: dx
REAL (dp), INTENT(OUT)  :: dy
INTEGER, INTENT(OUT)    :: ier

!***********************************************************

!                                               From CSHEP2D
!                                            Robert J. Renka
!                                  Dept. of Computer Science
!                                       Univ. of North Texas
!                                           renka@cs.unt.edu
!                                                   03/28/97

!   Given a set of N arbitrarily distributed nodes in the
! plane, this subroutine creates a data structure for a
! cell-based method of solving closest-point problems.  The
! smallest rectangle containing the nodes is partitioned
! into an NR by NR uniform grid of cells, and nodes are as-
! sociated with cells.  In particular, the data structure
! stores the indexes of the nodes contained in each cell.
! For a uniform random distribution of nodes, the nearest
! node to an arbitrary point can be determined in constant
! expected time.

!   This code is essentially unaltered from the subroutine
! of the same name in QSHEP2D.

! On input:

!       N = Number of nodes.  N >= 2.

!       X,Y = Arrays of length N containing the Cartesian
!             coordinates of the nodes.

!       NR = Number of rows and columns in the grid.  The
!            cell density (average number of nodes per cell)
!            is D = N/(NR**2).  A recommended value, based
!            on empirical evidence, is D = 3 -- NR =
!            Sqrt(N/3).  NR >= 1.

! The above parameters are not altered by this routine.

!       LCELL = Array of length >= NR**2.

!       LNEXT = Array of length >= N.

! On output:

!       LCELL = NR by NR cell array such that LCELL(I,J)
!               contains the index (for X and Y) of the
!               first node (node with smallest index) in
!               cell (I,J), or LCELL(I,J) = 0 if no nodes
!               are contained in the cell.  The upper right
!               corner of cell (I,J) has coordinates (XMIN+
!               I*DX,YMIN+J*DY).  LCELL is not defined if IER .NE. 0.

!       LNEXT = Array of next-node indexes such that
!               LNEXT(K) contains the index of the next node
!               in the cell which contains node K, or
!               LNEXT(K) = K if K is the last node in the
!               cell for K = 1,...,N.  (The nodes contained
!               in a cell are ordered by their indexes.)
!               If, for example, cell (I,J) contains nodes
!               2, 3, and 5 (and no others), then LCELL(I,J)
!               = 2, LNEXT(2) = 3, LNEXT(3) = 5, and
!               LNEXT(5) = 5.  LNEXT is not defined if IER .NE. 0.

!       XMIN,YMIN = Cartesian coordinates of the lower left
!                   corner of the rectangle defined by the
!                   nodes (smallest nodal coordinates) un-
!                   less IER = 1.  The upper right corner is
!                   (XMAX,YMAX) for XMAX = XMIN + NR*DX and
!                   YMAX = YMIN + NR*DY.

!       DX,DY = Dimensions of the cells unless IER = 1.  DX
!               = (XMAX-XMIN)/NR and DY = (YMAX-YMIN)/NR,
!               where XMIN, XMAX, YMIN, and YMAX are the
!               extrema of X and Y.

!       IER = Error indicator:
!             IER = 0 if no errors were encountered.
!             IER = 1 if N < 2 or NR < 1.
!             IER = 2 if DX = 0 or DY = 0.

! Modules required by STORE2:  None

! Intrinsic functions called by STORE2:  DBLE, INT

!***********************************************************

INTEGER   :: i, j, k, l, nn, nnr
REAL (dp) :: delx, dely, xmn, xmx, ymn, ymx

! Local parameters:

! DELX,DELY = Components of the cell dimensions -- local copies of DX,DY
! I,J =       Cell indexes
! K =         Nodal index
! L =         Index of a node in cell (I,J)
! NN =        Local copy of N
! NNR =       Local copy of NR
! XMN,XMX =   Range of nodal X coordinates
! YMN,YMX =   Range of nodal Y coordinates

nn = n
nnr = nr
IF (nn < 2 .OR. nnr < 1) GO TO 5

! Compute the dimensions of the rectangle containing the nodes.

xmn = x(1)
xmx = xmn
ymn = y(1)
ymx = ymn
DO  k = 2,nn
  IF (x(k) < xmn) xmn = x(k)
  IF (x(k) > xmx) xmx = x(k)
  IF (y(k) < ymn) ymn = y(k)
  IF (y(k) > ymx) ymx = y(k)
END DO
xmin = xmn
ymin = ymn

! Compute cell dimensions and test for zero area.

delx = (xmx-xmn)/DBLE(nnr)
dely = (ymx-ymn)/DBLE(nnr)
dx = delx
dy = dely
IF (delx == 0. .OR. dely == 0.) GO TO 6

! Initialize LCELL.

lcell(1:nnr,1:nnr) = 0

! Loop on nodes, storing indexes in LCELL and LNEXT.

DO  k = nn,1,-1
  i = INT((x(k)-xmn)/delx) + 1
  IF (i > nnr) i = nnr
  j = INT((y(k)-ymn)/dely) + 1
  IF (j > nnr) j = nnr
  l = lcell(i,j)
  lnext(k) = l
  IF (l == 0) lnext(k) = k
  lcell(i,j) = k
END DO

! No errors encountered.

ier = 0
RETURN

! Invalid input parameter.

5 ier = 1
RETURN

! DX = 0 or DY = 0.

6 ier = 2
RETURN
END SUBROUTINE store2

END MODULE Cubic_Shepard
