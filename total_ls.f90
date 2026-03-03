MODULE Total_LS

!     Total Least Squares

!     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!     REVISIONS: 1988, February 15.

! From the vanhuffel directory at netlib
! Code converted using TO_F90 by Alan Miller
! Date: 2001-04-27  Time: 17:51:50


IMPLICIT NONE
INTEGER, PARAMETER, PRIVATE  :: dp = SELECTED_REAL_KIND(12, 60)

CONTAINS


SUBROUTINE ptls(c,ldc,m,n,l,rank,theta,x,ldx,q,inul,tol1,tol2,ierr,iwarn)

! N.B. Arguments WRK, IWRK & LWRK have been removed.

!   Subroutine PTLS solves a set of linear equations by a Total Least
!   Squares Approximation, based on the Partial SVD.

!   The TLS problem assumes an overdetermined set of linear equations
!   AX = B, where both the data matrix A as well as the observation
!   matrix B are inaccurate. If the perturbations D on the data [A;B]
!   have zero mean and if their covariance matrix E(D'D) with E the expected
!   value operator, equals the identity matrix up to an unknown
!   scaling factor (e.g. when all the errors are independent and equally
!   sized), then the routine computes a strongly consistent estimate
!   of the true solution of the corresponding unperturbed set.
!   The routine also solves determined and underdetermined sets of
!   equations by computing the minimum norm solution.
!   It is assumed that all preprocessing measures (scaling, coordinate
!   transformations, whitening, ... ) of the data have been performed
!   in advance.

!   C - DOUBLE PRECISION array of DIMENSION (LDC,N+L).
!       The leading M x (N + L) part of this array contains A and B
!       as follows:
!       The first N columns contain the data matrix A, and the last
!       L columns the observation matrix B (right-hand sides).
!       NOTE that this array is overwritten.
!
!   LDC - INTEGER.
!       The declared leading dimension of the array C.
!       LDC >= max(M,(N+L)).
!
!   M - INTEGER.
!       The number of rows in the data matrix A and the observation
!       matrix B.
!
!   N - INTEGER.
!       The number of columns in the data matrix A.
!
!   L - INTEGER.
!       The number of columns in the observation matrix B.
!
!   RANK - INTEGER.
!       The rank of the TLS approximation [A+DA;B+DB].
!       On entry there are two possibilities:
!       i)  RANK <  0: the rank is to be computed by the routine.
!       ii) RANK >= 0: specifies the given rank.
!           RANK <= min(M,N).
!           NOTE that RANK may be overwritten, if C has multiple
!                singular values or if the upper triangular matrix
!                F is singular (see 6 METHOD DESCRIPTION).
!
!   THETA - DOUBLE PRECISION.
!       On entry, there are two possibilities (depending on RANK):
!       i)  RANK <  0: the rank of the TLS approximation [A+DA;B+DB]
!           is computed from THETA and equals min(M,N+L) - d, where
!           d equals the number of singular values of [A;B] <= THETA.
!           THETA >= 0.
!       ii) RANK >= 0: THETA is an initial estimate for computing a
!           lower bound of the RANK largest singular values of [A;B].
!           If not available, assign a negative value (< 0) to THETA.
!           NOTE that THETA is overwritten.
!
!   ARGUMENTS OUT
!
!   C - DOUBLE PRECISION array of DIMENSION (LDC,N+L).
!       The first N + L components of the columns of C whose
!       index i corresponds with an INUL(i) = .TRUE., are the
!       N + L - RANK base vectors of the right singular subspace
!       corresponding to the singular values of C which are <= THETA.
!       The TLS solution is computed from these vectors.
!
!   RANK - INTEGER.
!       RANK specifies the rank of [A+DA;B+DB].
!       If not specified by the user, then RANK is computed by the
!       routine.
!       If specified by the user, then the specified RANK is changed
!       by the routine if the RANK-th and the (RANK+1)-th singular
!       value of [A;B] coincide up to the precision given by TOL1,
!       or if the upper triangular matrix F (see 6 METHOD DESCRIP-
!       TION) is singular.
!
!   THETA - DOUBLE PRECISION.
!       If RANK >= 0, then THETA specifies the computed bound such
!       that precisely RANK singular values of [A;B] are greater
!       than THETA + TOL1.
!
!   X - DOUBLE PRECISION array of DIMENSION (LDX,L).
!       The leading N x L part of this array contains the solutions
!       to the TLS problems specified by A and B.
!
!   LDX - INTEGER.
!       The declared leading dimension of the array X .
!       LDX >= N.
!
!   Q - DOUBLE PRECISION array of DIMENSION (min(M,N+L)+min(M+1,N+L))
!       Returns the partially diagonalized bidiagonal computed from
!       C, at the moment that the desired singular subspace has been
!       found.
!       The first p = min(M,N+L) entries of Q contain the diagonal
!       elements q(1),...,q(p) and the entries Q(p+2),...,Q(p+s)
!       (with s = min(M+1,N+L)) contain the superdiagonal elements
!       e(2),...,e(s). Q(p+1) = 0.
!
!   INUL - LOGICAL array of DIMENSION (N+L).
!       The indices of the elements of INUL with value .TRUE. indicate
!       the columns in C containing the base vectors of the right singular
!       subspace of C from which the TLS solution has been computed.

!   WORK SPACE
!
!   WRK - DOUBLE PRECISION array of DIMENSION (t).
!       The dimension t = M + N + L + (N+L)*(N+L-1)/2 if M >= N + L;
!                         M + N + L + M*(N+L-(M+1)/2) if M <  N + L.
!
!   IWRK - INTEGER array of DIMENSION (N+L).
!
!   LWRK - LOGICAL array of DIMENSION (N+L).

!   TOLERANCES
!
!   TOL1 - DOUBLE PRECISION.
!       This parameter defines the multiplicity of singular values by
!       considering all singular values within an interval of length
!       TOL1 as coinciding. TOL1 is used in checking how many singu-
!       lar values are <= THETA. Also in computing an appropriate
!       upper bound THETA by a bisection method, TOL1 is used as stop
!       criterion defining the minimal subinterval length.
!       According to the numerical properties of the SVD, TOL1 must
!       be >= !!C!! x EPS where !!C!! denotes the L2-norm and EPS
!       is the machine precision.
!
!   TOL2 - DOUBLE PRECISION.
!       Working precision for the computations in PTLS.  This parameter
!       specifies that elements of matrices used in the computation,
!       which are <= TOL2 in absolute value, are considered to be zero.

!   ERROR INDICATORS
!
!   IERR - INTEGER.
!       On return, IERR contains 0 unless the subroutine has failed.
!
!   IWARN - INTEGER.
!       On return, IWARN contains 0 unless RANK has been lowered by the
!       routine.

!   ERROR INDICATORS and WARNINGS:
!
!   Errors detected by the routine:
!
!       IERR = 0: successful completion.
!              1: number M of rows of array C = [A;B] smaller than 1.
!              2: number N of columns of matrix A smaller than 1.
!              3: number L of columns of matrix B smaller than 1.
!              4: leading dimension LDC of array C is smaller than
!                 max(M,N+L).
!              5: leading dimension LDX of array X is smaller than N.
!              6: rank of the TLS approximation [A+DA;B+DB] is larger
!                 than min(M,N).
!              7: the parameters RANK and THETA are both negative (< 0).
!              8: tolerance TOL1 is negative.
!              9: tolerance TOL2 is negative.
!             10: maximum number of QR/QL iteration steps (50) exceeded.

!   Warnings given by the routine:
!
!       IWARN = 0: no warnings.
!               1: the rank of matrix C has been lowered because a
!                  singular value of multiplicity > 1 has been found.
!               2: the rank of matrix C has been lowered because a
!                  singular upper triangular matrix F has been found.

REAL (dp), INTENT(OUT)     :: c(:,:)   ! c(ldc,*)
INTEGER, INTENT(IN)        :: ldc
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: l
INTEGER, INTENT(IN OUT)    :: rank
REAL (dp), INTENT(IN OUT)  :: theta
REAL (dp), INTENT(OUT)     :: x(:,:)   ! x(ldx,*)
INTEGER, INTENT(IN)        :: ldx
REAL (dp), INTENT(IN OUT)  :: q(:)
LOGICAL, INTENT(OUT)       :: inul(n+l)
REAL (dp), INTENT(IN)      :: tol1
REAL (dp), INTENT(IN)      :: tol2
INTEGER, INTENT(OUT)       :: ierr
INTEGER, INTENT(OUT)       :: iwarn

!     .. External Subroutines/Functions ..
! REAL (dp) :: ddot
! EXTERNAL daxpy, ddot, dcopy, dqrdc, bidiag, cancel, housh, init, qrql
!     .. Intrinsic Functions ..
! INTRINSIC ABS, MAX, MIN
!     .. Local Scalars ..
INTEGER    :: i, j, p, nl, pp1, j1, mc, nj, k, ii, i1, n1, mnl1
REAL (dp)  :: temp, hh, inprod
LOGICAL    :: zero
!     .. Local Arrays ..
INTEGER    :: dumar1(1)
REAL (dp)  :: dumar2(1,1)
!     .. Workspaces ..
REAL (dp), ALLOCATABLE  :: wrk(:)
INTEGER  :: iwrk(n+l)
LOGICAL  :: lwrk(n+l)
!     .. Executable Statements ..

ierr = 0
iwarn = 0
IF (m < 1) ierr = 1
IF (n < 1) ierr = 2
IF (l < 1) ierr = 3
IF (ldc < MAX(m,n+l)) ierr = 4
IF (ldx < n) ierr = 5
IF (rank > MIN(m,n)) ierr = 6
IF (rank < 0 .AND. theta < 0.0D0) ierr = 7
IF (tol1 < 0.0D0) ierr = 8
IF (tol2 < 0.0D0) ierr = 9
IF (ierr /= 0) RETURN

!     Initializations.

nl = n + l
IF (m >= nl) THEN
  ALLOCATE( wrk(m + nl + nl*(nl-1)/2) )
ELSE
  ALLOCATE( wrk(m + nl + m*(nl-(m+1)/2)) )
END IF
mnl1 = m + nl + 1
p = MIN(m,nl)
DO  i = 1, p
  inul(i) = .false.
  lwrk(i) = .false.
END DO

pp1 = p + 1
DO  i = pp1, nl
  inul(i) = .true.
  lwrk(i) = .false.
END DO

x(1:n,1:l) = 0.0D0

!     Step 1: Bidiagonalization phase
!             -----------------------
!     1.a): If M >= 5*(N+L)/3, transform C into upper triangular form R
!           by Householder transformations.

mc = m
IF (3*m >= 5*nl) THEN
  CALL dqrdc(c,m,nl,q,dumar1,0)
  DO  j = 1, nl
    j1 = j + 1
    c(j1:nl,j) = 0.0D0
  END DO
  mc = nl
END IF

!     1.b): Transform C (or R) into bidiagonal form Q using Householder
!           transformations.

CALL bidiag(c,mc,nl,q,q(p+1:))

!     Store the Householder transformations performed onto the rows of C
!     in the last storage locations of the work array WRK.

mc = MIN(nl-2,m)
IF (mc > 0) THEN
  k = mnl1
  DO  ii = 1, mc
    j = mc - ii + 1
    nj = nl - j
    wrk(k:k+nj-1) = c(j,j+1:j+nj)
    k = k + nj
  END DO
END IF

!     1.c): Initialize the right singular base matrix V with the identity
!           matrix (V overwrites C).

CALL init(c,nl,nl)

!     1.d): If M < N+L, bring the bidiagonal Q to M by M by cancelling
!           its last superdiagonal element using Givens rotations.

pp1 = p
IF (m < nl) THEN
  pp1 = p + 1
  CALL cancel(dumar2,c,q,q(pp1:),p,pp1,pp1,pp1,tol2,.false.,.true.)
END IF

!     REPEAT

!        Compute the Householder matrix Q and matrices F and Y such that
!        F is nonsingular.

!        Step 2: Partial diagonalization phase.
!                -----------------------------
!        Diagonalize the bidiagonal Q partially until convergence to
!        the desired right singular subspace.

80 CALL qrql(dumar2,c,p,pp1,rank,theta,q,q(p+1:),inul,tol1,tol2,  &
             .false.,.true.,ierr,iwarn)

IF (ierr /= 0) RETURN

!        Step 3: Back transformation phase.
!                -------------------------
!        Apply the Householder transformations (stored in WRK) performed
!        onto the rows of C during the bidiagonalization phase, to the
!        selected base vectors in the right singular base matrix of C.

DO  i = 1, nl
  IF (inul(i) .AND. (.NOT.lwrk(i))) THEN
    k = mnl1
    DO  ii = 1, mc
      j = mc - ii + 1
      nj = nl - j
      j1 = j + 1
      IF (ABS(wrk(k)) > tol2) THEN
        temp = -DOT_PRODUCT( wrk(k:k+nj-1), c(j1:j1+nj-1,i) ) / wrk(k)
        c(j1:j1+nj-1,i) = c(j1:j1+nj-1,i) + temp * wrk(k:k+nj-1)
        k = k + nj
      END IF
    END DO
    lwrk(i) = .true.
  END IF
END DO
IF (rank <= 0) RETURN

!        Step 4: Compute matrices F and Y using Householder transf. Q.
!                ------------------------
k = 0
DO  i = 1, nl
  IF (inul(i)) THEN
    k = k + 1
    iwrk(k) = i
  END IF
END DO

IF (k < l) THEN
  
!           Rank TLS approximation is larger than min(M,N).
  
  ierr = 3
  RETURN
END IF

n1 = n + 1
zero = .false.
i = nl

!        WHILE ((K > 1) .AND. (I > N) .AND. (.NOT.ZERO)) DO
120 IF (k > 1 .AND. i > n .AND. (.NOT.zero)) THEN
  DO  j = 1, k
    wrk(j) = c(i,iwrk(j))
  END DO
  
!           Compute Householder transformation.
  
  CALL housh(wrk,k,k,tol2,zero,temp)
  IF (.NOT.zero) THEN
    
!             Apply Householder transformation onto the selected base vectors.
    
    DO  i1 = 1, i
      inprod = 0.0_dp
      DO  j = 1, k
        inprod = inprod + wrk(j) * c(i1,iwrk(j))
      END DO
      hh = inprod * temp
      DO  j = 1, k
        j1 = iwrk(j)
        c(i1,j1) = c(i1,j1) - wrk(j) * hh
      END DO
    END DO
    
    k = k - 1
  END IF
  i = i - 1
  GO TO 120
END IF
!        END WHILE 10

IF (.NOT.zero) k = n1 - rank

!        If F singular, lower the rank of the TLS approximation .

IF (ABS(c(n1,iwrk(k))) <= tol2) THEN
  rank = rank - 1
  iwarn = 2
  theta = -1.0D0
  GO TO 80
END IF
!     UNTIL ((.NOT.ZERO) .AND. (F nonsingular))

!     Step 5: Compute TLS solution
!             --------------------
!     Solve X F = -Y  by forward elimination  (F is upper triangular).

nj = iwrk(k)
x(1:n,1) = x(1:n,1) - c(1:n,nj) / c(n1,nj)

DO  j = 2, l
  j1 = j - 1
  nj = iwrk(k+j1)
  temp = c(n+j,nj)
  DO  i = 1, n
    x(i,j) = -(c(i,nj) + DOT_PRODUCT( c(n1:n1+j1-1,nj), x(i,1:j1) )) / temp
  END DO
END DO
RETURN
! *** Last line of PTLS ***********************************************
END SUBROUTINE ptls



SUBROUTINE qlstep(u,v,q,e,m,n,i,k,shift,wantu,wantv)

! N.B. Arguments LDU & LDV have been removed.
 
!  PURPOSE:

!  The subroutine QLSTEP performs one QL iteration step onto the
!  unreduced subbidiagonal Jk:

!           !Q(i) E(i+1)  0  ...    0  !
!           ! 0   Q(i+1) E(i+2)     .  !
!      Jk = ! .                     .  !
!           ! .                        !
!           ! .                    E(k)!
!           ! 0   ...              Q(k)!

!  with k <= p and i >= 1, p = min(M,N), of the bidiagonal J:

!           !Q(1) E(2)  0    ...   0  !
!           ! 0   Q(2) E(3)        .  !
!       J = ! .                    .  !
!           ! .                   E(p)!
!           ! 0   ...             Q(p)!

!  Hereby, Jk is transformed to  S'Jk T with S and T products of
!  Givens rotations. These Givens rotations S (resp.,T) will be post-
!  multiplied into U (resp.,V), if WANTU (resp.,WANTV) = .TRUE.

!  ARGUMENT LIST:

!  U - REAL (dp) array of DIMENSION (LDU,min(M,N)).
!      On entry, U may contain the M by p (p=min(M,N)) left transformation
!      matrix.
!      On return, if WANTU = .TRUE., the Givens rotations S on the
!      left have been postmultiplied into U.
!      NOTE: U is not referenced if WANTU = .FALSE.
!  LDU - INTEGER.
!      LDU is the leading dimension of the array U (LDU >= M).
!  V - REAL (dp) array of DIMENSION (LDV,min(M,N)).
!      On entry, V may contain the N by p (p=min(M,N)) right transformation
!      matrix.
!      On return, if WANTV = .TRUE., the Givens rotations T on the
!      right have been postmultiplied into V.
!      NOTE: V is not referenced if WANTV is .false.
!  LDV - INTEGER.
!      LDV is the leading dimension of the array V (LDV >= N).
!  Q - REAL (dp) array of DIMENSION (min(M,N)).
!      On entry, Q contains the diagonal entries of the bidiagonal J.
!      On return, Q contains the diagonal entries of the transformed
!      matrix S' J T.
!  E - REAL (dp) array of DIMENSION (min(M,N)).
!      On entry, E contains the superdiagonal entries of J.
!      On return, E contains the superdiagonal entries of the transformed
!      matrix S' J T. E(k+1) = 0 if k < min(M,N).
!  M - INTEGER.
!      M is the number of rows of the matrix U.
!  N - INTEGER.
!      N is the number of rows of the matrix V.
!  I - INTEGER.
!      I is the index of the first diagonal entry of the considered
!      unreduced subbidiagonal Jk of J.
!  K - INTEGER.
!      K is the index of the last diagonal entry of the considered
!      unreduced subbidiagonal Jk of J.
!  SHIFT - REAL (dp).
!      Value of the shift used in the QL iteration step.
!  WANTU - LOGICAL.
!      WANTU = .TRUE. if the Givens rotations S must be postmultiplied
!      on the left into U, else .FALSE.
!  WANTV - LOGICAL.
!      WANTV = .TRUE. if the Givens rotations T must be postmultiplied
!      on the left into V, else .FALSE.

!  EXTERNAL SUBROUTINES AND FUNCTIONS:

!  DROT from BLAS.

!  METHOD DESCRIPTION:

!  QL iterations diagonalize the bidiagonal by zeroing the super-
!  diagonal elements of Jk from top to bottom.
!  The routine QLSTEP overwrites Jk with the bidiagonal matrix
!  S' Jk T where S and T are Givens rotations.
!  T is essentially the orthogonal matrix that would be obtained by
!  applying one implicit symmetric shift QL step onto the matrix
!  Jk'Jk. This step factors the matrix (Jk'Jk - shift*I) into a
!  product of an orthogonal matrix T and a lower triangular matrix.
!  See [1,Sec.8.2-8.3] and [2] for more details.

!  REFERENCES:
!  [1] G.H. Golub and C.F. Van Loan, Matrix Computations. The Johns
!      Hopkins University Press, Baltimore,Maryland (1983).
!  [2] H. Bowdler, R.S. Martin and J.H. Wilkinson, The QR and QL algorithms
!      for symmetric matrices. Numer. Math., 11 (1968), 293-306.

!  CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!  REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN OUT)  :: u(:,:)   ! u(ldu,*)
REAL (dp), INTENT(IN OUT)  :: v(:,:)   ! v(ldv,*)
REAL (dp), INTENT(IN OUT)  :: q(:)
REAL (dp), INTENT(IN OUT)  :: e(:)
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: i
INTEGER, INTENT(IN)        :: k
REAL (dp), INTENT(IN)      :: shift
LOGICAL, INTENT(IN)        :: wantu
LOGICAL, INTENT(IN)        :: wantv

!     .. External Subroutines/Functions ..
! EXTERNAL drot
!     .. Intrinsic Functions ..
! INTRINSIC MIN, SQRT
!     .. Local Scalars ..
INTEGER    :: ik, j, ll, jj
REAL (dp)  :: f, g, h, c, s, x, y, z
!     .. Executable Statements ..

x = q(k)
g = shift
f = (x-g) * (x+g) / x
c = 1.0D0
s = 1.0D0
ik = k - i
DO  jj = 1, ik
  j = k - jj
  ll = j + 1
  g = e(ll)
  y = q(j)
  h = s * g
  g = c * g
  z = SQRT(f**2 + h**2)
  IF (jj > 1) e(ll+1) = z
  c = f / z
  s = h / z
  f = x * c + g * s
  g = -x * s + g * c
  h = y * s
  y = y * c
  IF (wantu) CALL drot(m,u(1:,ll),1,u(1:,j),1,c,s)
  z = SQRT(f**2 + h**2)
  q(ll) = z
  c = f / z
  s = h / z
  f = c * g + s * y
  x = -s * g + c * y
  IF (wantv) CALL drot(n,v(1:,ll),1,v(1:,j),1,c,s)
END DO
e(i+1) = f
q(i) = x
IF (k < MIN(m,n)) e(k+1) = 0.0D0
RETURN
! Last line of QLSTEP *************************************************
END SUBROUTINE qlstep



SUBROUTINE qrql(u,v,m,n,rank,theta,q,e,inul,tol1,tol2,wantu,wantv,ierr,iwarn)

! N.B. Arguments LDU & LDV have been removed.
 
!     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!     REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN OUT)  :: u(:,:)   ! u(ldu,*)
REAL (dp), INTENT(IN OUT)  :: v(:,:)   ! v(ldv,*)
INTEGER, INTENT(IN OUT)    :: m
INTEGER, INTENT(IN OUT)    :: n
INTEGER, INTENT(IN OUT)    :: rank
REAL (dp), INTENT(IN OUT)  :: theta
REAL (dp), INTENT(IN OUT)  :: q(:)
REAL (dp), INTENT(IN OUT)  :: e(:)
LOGICAL, INTENT(IN OUT)    :: inul(:)
REAL (dp), INTENT(IN)      :: tol1
REAL (dp), INTENT(IN)      :: tol2
LOGICAL, INTENT(IN)        :: wantu
LOGICAL, INTENT(IN)        :: wantv
INTEGER, INTENT(OUT)       :: ierr
INTEGER, INTENT(IN OUT)    :: iwarn

!     .. External Subroutines/Functions ..
! INTEGER :: nsingv
! EXTERNAL cancel, estim, qrstep, qlstep, nsingv
!     .. Intrinsic Functions ..
! INTRINSIC ABS, MIN
!     .. Local Scalars ..
INTEGER    :: i1, k, j, i, numeig, iter, maxit, p, r
REAL (dp)  :: x, shift
LOGICAL    :: noc12
!     .. Executable Statements ..

ierr = 0
p = MIN(m,n)

!     Estimate THETA (if not fixed by the user), and set R.

IF (rank >= 0) THEN
  j = p - rank
  CALL estim(q,e,p,j,theta,tol1,tol2,iwarn)
  IF (j <= 0) RETURN
  r = p - j
ELSE
  r = 0
END IF

maxit = 50
rank = p
DO  i = 1, p
  IF (inul(i)) rank = rank - 1
END DO
e(1) = 0.0D0

!     From now K is the smallest known index such that the subbidia-
!     gonals with indices > K belong to C1 or C2.
!     RANK = P - SUM(dimensions of known elements of C2).

k = p
!     WHILE (C3 NOT EMPTY) DO
20 IF (rank > r .AND. k > 0) THEN
!        WHILE (K > 0 .AND INUL(K)) DO
  
!        Search for the rightmost index of a subbidiagonal of C1 or C3.
  
  30 IF (k > 0) THEN
    IF (inul(k)) THEN
      k = k - 1
      GO TO 30
    END IF
  END IF
!        END WHILE 3
  
  IF (k == 0) RETURN
  
  iter = 0
  noc12 = .true.
!        WHILE ((ITER < MAXIT) .AND. (No element of C1 or C2 found)) DO
  40 IF (iter < maxit .AND. noc12) THEN
    
!           Search for negligible Q(I) or E(I).
    
    i = k
    x = ABS(q(i))
    shift = x
!           WHILE (ABS(Q(I)) > TOL2 .AND. ABS(E(I)) > TOL2) DO
    50 IF (x > tol2 .AND. ABS(e(i)) > tol2) THEN
      i = i - 1
      x = ABS(q(i))
      IF (x < shift) shift = x
      GO TO 50
    END IF
!           END WHILE 5
    
!           Classify the subbidiagonal found.
    
    j = k - i + 1
    IF (x <= tol2 .OR. k == i) THEN
      noc12 = .false.
    ELSE
      numeig = nsingv(q(i:),e(i:),j,theta,tol1,tol2)
      IF (numeig >= j .OR. numeig <= 0) noc12 = .false.
    END IF
    IF (noc12) THEN
      IF (shift >= theta) shift = 0.0D0
      IF (ABS(q(k)) <= ABS(q(i))) THEN
        CALL qrstep(u,v,q,e,m,n,i,k,shift,wantu,wantv)
      ELSE
        CALL qlstep(u,v,q,e,m,n,i,k,shift,wantu,wantv)
      END IF
      iter = iter + 1
    END IF
    GO TO 40
  END IF
!        END WHILE 4
  
  IF (iter == maxit) THEN
    ierr = 10
    RETURN
  END IF
  
  IF (x <= tol2) THEN
    
!           Split at negligible diagonal element abs(Q(I)) <= TOL2.
    
    CALL cancel(u,v,q,e,m,n,i,k,tol2,wantu,wantv)
    inul(i) = .true.
    rank = rank - 1
  ELSE
    
!           A negligible superdiagonal element abs(E(I)) <= TOL2 has
!           been found, the corresponding subbidiagonal belongs to
!           C1 or C2. Treat this subbidiagonal.
    
    IF (j >= 2) THEN
      IF (numeig == j) THEN
        DO  i1 = i, k
          inul(i1) = .true.
        END DO
        rank = rank - j
        k = k - j
      ELSE
        k = i - 1
      END IF
    ELSE
      IF (x <= theta + tol1) THEN
        inul(i) = .true.
        rank = rank - 1
      END IF
      k = k - 1
    END IF
  END IF
  GO TO 20
END IF
!     END WHILE 2
RETURN
! *** Last line of QRQL ***********************************************
END SUBROUTINE qrql



SUBROUTINE qrstep(u,v,q,e,m,n,i,k,shift,wantu,wantv)

! N.B. Arguments LDU & LDV have been removed.
 
!  PURPOSE:

!  The subroutine QRSTEP performs one QR iteration step onto the
!  unreduced subbidiagonal Jk:

!           !Q(i) E(i+1)  0  ...    0  !
!           ! 0   Q(i+1) E(i+2)     .  !
!      Jk = ! .                     .  !
!           ! .                        !
!           ! .                    E(k)!
!           ! 0   ...              Q(k)!

!  with k <= p and i >= 1, p = min(M,N), of the bidiagonal J:

!           !Q(1) E(2)  0    ...   0  !
!           ! 0   Q(2) E(3)        .  !
!       J = ! .                    .  !
!           ! .                   E(p)!
!           ! 0   ...             Q(p)!

!  Hereby, Jk is transformed to  S'Jk T with S and T products of
!  Givens rotations. These Givens rotations S (resp.,T) will be post-
!  multiplied into U (resp.,V), if WANTU (resp.,WANTV) = .TRUE.

!  ARGUMENT LIST:

!  U - REAL (dp) array of DIMENSION (LDU,min(M,N)).
!      On entry, U may contain the M by p (p=min(M,N)) left transformation
!      matrix.
!      On return, if WANTU = .TRUE., the Givens rotations S on the
!      left have been postmultiplied into U.
!      NOTE: U is not referenced if WANTU = .FALSE.
!  LDU - INTEGER.
!      LDU is the leading dimension of the array U (LDU >= M).
!  V - REAL (dp) array of DIMENSION (LDV,min(M,N)).
!      On entry, V may contain the N by p (p=min(M,N)) right transformation
!      matrix.
!      On return, if WANTV = .TRUE., the Givens rotations T on the
!      right have been postmultiplied into V.
!      NOTE: V is not referenced if WANTV is .false.
!  LDV - INTEGER.
!      LDV is the leading dimension of the array V (LDV >= N).
!  Q - REAL (dp) array of DIMENSION (min(M,N)).
!      On entry, Q contains the diagonal entries of the bidiagonal J.
!      On return, Q contains the diagonal entries of the transformed
!      matrix S' J T.
!  E - REAL (dp) array of DIMENSION (min(M,N)).
!      On entry, E contains the superdiagonal entries of J.
!      On return, E contains the superdiagonal entries of the transformed
!      matrix S' J T. E(i) = 0.
!  M - INTEGER.
!      M is the number of rows of the matrix U.
!  N - INTEGER.
!      N is the number of rows of the matrix V.
!  I - INTEGER.
!      I is the index of the first diagonal entry of the considered
!      unreduced subbidiagonal Jk of J.
!  K - INTEGER.
!      K is the index of the last diagonal entry of the considered
!      unreduced subbidiagonal Jk of J.
!  SHIFT - REAL (dp).
!      Value of the shift used in the QR iteration step.
!  WANTU - LOGICAL.
!      WANTU = .TRUE. if the Givens rotations S must be postmultiplied
!      on the left into U, else .FALSE.
!  WANTV - LOGICAL.
!      WANTV = .TRUE. if the Givens rotations T must be postmultiplied
!      on the left into V, else .FALSE.

!  EXTERNAL SUBROUTINES AND FUNCTIONS:

!  DROT from BLAS.

!  METHOD DESCRIPTION:

!  QR iterations diagonalize the bidiagonal by zeroing the superdiagonal
!  elements of Jk from bottom to top.
!  The routine QRSTEP overwrites Jk with the bidiagonal matrix
!  S' Jk T where S and T are Givens rotations.
!  T is essentially the orthogonal matrix that would be obtained by
!  applying one implicit symmetric shift QR step onto the matrix
!  Jk'Jk. This step factors the matrix (Jk'Jk - shift*I) into a
!  product of an orthogonal matrix T and an upper triangular matrix.
!  See [1,Sec.8.2-8.3] for more details.

!  REFERENCES:
!  [1] G.H. Golub and C.F. Van Loan, Matrix Computations. The Johns
!      Hopkins University Press, Baltimore,Maryland (1983).

!  CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!  REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN OUT)  :: u(:,:)   ! u(ldu,*)
REAL (dp), INTENT(IN OUT)  :: v(:,:)   ! v(ldv,*)
REAL (dp), INTENT(IN OUT)  :: q(:)
REAL (dp), INTENT(IN OUT)  :: e(:)
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: i
INTEGER, INTENT(IN)        :: k
REAL (dp), INTENT(IN)      :: shift
LOGICAL, INTENT(IN)        :: wantu
LOGICAL, INTENT(IN)        :: wantv

!     .. External Subroutines/Functions ..
! EXTERNAL drot
!     .. Intrinsic Functions ..
! INTRINSIC SQRT
!     .. Local Scalars ..
INTEGER    :: i1, j, ll
REAL (dp)  :: f, g, h, c, s, x, y, z
!     .. Executable Statements ..

x = q(i)
g = shift
f = (x-g) * (x+g) / x
c = 1.0D0
s = 1.0D0
i1 = i + 1
DO  j = i1, k
  ll = j - 1
  g = e(j)
  y = q(j)
  h = s * g
  g = c * g
  z = SQRT(f**2 + h**2)
  e(ll) = z
  c = f / z
  s = h / z
  f = x * c + g * s
  g = -x * s + g * c
  h = y * s
  y = y * c
  IF (wantv) CALL drot(n,v(1:,ll),1,v(1:,j),1,c,s)
  z = SQRT(f**2 + h**2)
  q(ll) = z
  c = f / z
  s = h / z
  f = c * g + s * y
  x = -s * g + c * y
  IF (wantu) CALL drot(m,u(1:,ll),1,u(1:,j),1,c,s)
END DO
e(k) = f
q(k) = x
e(i) = 0.0D0
RETURN
! *** Last line of QRSTEP *********************************************
END SUBROUTINE qrstep



SUBROUTINE bidiag(x,n,p,q,e)

! N.B. Arguments LDX & WRK have been removed.
 
!     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!     REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN OUT)  :: x(:,:)   ! x(ldx,*)
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: p
REAL (dp), INTENT(OUT)     :: q(:)
REAL (dp), INTENT(OUT)     :: e(:)

!     .. External Subroutines/Functions ..
! REAL (dp) :: dnrm2, ddot
! EXTERNAL dscal, dnrm2, ddot, daxpy
!     .. Intrinsic Functions ..
! INTRINSIC MAX, MIN, SIGN
!     .. Local Scalars ..
INTEGER    :: nct, nrt, lu, l, lp1, j, i, npp, pp1, m
REAL (dp)  :: t, wrk(n+p)
!     .. Executable Statements ..

!     Reduce X to bidiagonal form, storing the diagonal elements in Q
!     and the superdiagonal elements in E.

npp = n + p
pp1 = p + 1
nct = MIN(n-1,p)
nrt = MAX(0,MIN(p-2,n))
lu = MAX(nct,nrt)
DO  l = 1, lu
  lp1 = l + 1
  IF (l <= nct) THEN
    
!           Compute the transformation for the L-th column and place the
!           L-th diagonal in Q(L).
    
    q(l) = dnrm2(n-l+1,x(l:,l),1)
    IF (q(l) /= 0.0D0) THEN
      IF (x(l,l) /= 0.0D0) q(l) = SIGN(q(l), x(l,l))
      x(l:n,l) = x(l:n,l) / q(l)
      x(l,l) = 1.0D0 + x(l,l)
    END IF
    q(l) = -q(l)
  END IF
  
  DO  j = lp1, p
    IF (l <= nct) THEN
      IF (q(l) /= 0.0D0) THEN
        
!              Apply the transformation.
        
        t = -DOT_PRODUCT( x(l:n,l), x(l:n,j) ) / x(l,l)
        x(l:n,j) = x(l:n,j) + t * x(l:n,l)
      END IF
    END IF
    
!           Place the L-th row of X into WRK for the subsequent
!           calculation of the row transformation.
    
    wrk(j) = x(l,j)
  END DO
  IF (l <= nrt) THEN
    
!           Compute the L-th row transformation and place the L-th
!           superdiagonal in E(L).
    
    wrk(l) = dnrm2(p-l, wrk(lp1:), 1)
    IF (wrk(l) /= 0.0D0) THEN
      IF (wrk(lp1) /= 0.0D0) wrk(l) = SIGN(wrk(l),wrk(lp1))
      wrk(lp1:lp1+p-l-1) = wrk(lp1:lp1+p-l-1) / wrk(l)
      wrk(lp1) = 1.0D0 + wrk(lp1)
    END IF
    wrk(l) = -wrk(l)
    e(lp1) = wrk(l)
    IF (lp1 <= n .AND. wrk(l) /= 0.0D0) THEN
      
!              Apply the transformation.
      
      DO  i = pp1, npp
        wrk(i) = 0.0D0
      END DO
      DO  j = lp1, p
        wrk(pp1:pp1+n-l-1) = wrk(pp1:pp1+n-l-1) + wrk(j) * x(lp1:lp1+n-l-1,j)
      END DO
      DO  j = lp1, p
        x(lp1:lp1+n-l-1,j) = x(lp1:lp1+n-l-1,j) -  &
                             (wrk(j)/wrk(lp1)) * wrk(pp1:pp1+n-l-1)
      END DO
    END IF
    x(l,lp1:p) = wrk(lp1:p)
  END IF
END DO

!     Set up the final bidiagonal matrix elements, if necessary.

e(1) = 0.0D0
m = MIN(p-1,n)
IF (nct < p) q(n) = x(n,n)
IF (nrt < m) e(p) = x(p-1,p)
RETURN
! *** Last line of BIDIAG *********************************************
END SUBROUTINE bidiag



SUBROUTINE cancel(u,v,q,e,m,n,i,k,tol,wantu,wantv)

! N.B. Arguments LDU & LDV have been removed.
 
!  PURPOSE:

!  Either, subroutine CANCEL separates a zero singular value of a
!  subbidiagonal matrix of order k, k <= p, of the bidiagonal

!            !Q(1) E(2)  0    ...   0  !
!            ! 0   Q(2) E(3)        .  !
!        J = ! .                    .  !
!            ! .                   E(p)!
!            ! 0   ...             Q(p)!

!  with p = min(M,N), by annihilating one or two superdiagonal
!  elements E(i) and/or E(i+1).
!  Or, CANCEL annihilates the element E(M+1) of the bidiagonal matrix

!            !Q(1) E(2)  0    ...   0     0   !
!            ! 0   Q(2) E(3)        .     .   !
!        J = ! .                    .     .   !
!            ! .                   E(M)   .   !
!            ! 0   ...             Q(M) E(M+1)!

!  ARGUMENT LIST:

!  U - REAL (dp) array of DIMENSION (LDU,p).
!      On entry, U contains the M by p (p = min(M,N)) left transformation
!      matrix.
!      On return, the Givens rotations S on the left, annihilating
!      E(i+1), have been postmultiplied into U.
!      NOTE: U is not referenced if WANTU = .FALSE. .
!  LDU - INTEGER.
!      LDU is the leading dimension of the array U (LDU >= M).
!  V - REAL (dp) array of DIMENSION (LDV,s).
!      On entry, V contains the N by s (s = min(M+1,N)) right transformation
!      matrix.
!      On return, the Givens rotations T on the right, annihilating
!      E(i), have been postmultiplied into V.
!      NOTE: V is not referenced if WANTV = .FALSE. .
!  LDV - INTEGER.
!      LDV is the leading dimension of the array V (LDV >= N).
!  Q - REAL (dp) array of DIMENSION (p).
!      On entry, Q contains the diagonal entries of the bidiagonal J.
!      p = min(M,N).
!      On return, Q contains the diagonal elements of the transformed
!      bidiagonal S' J T.
!  E - REAL (dp) array of DIMENSION (s).
!      On entry, E(i), i=2,...,s, contain the superdiagonal entries
!      of the bidiagonal J.  s = min(M+1,N), E(1) = 0.0D0.
!      On return, E contains the superdiagonal elements of the transformed
!      bidiagonal S' J T.
!  M - INTEGER.
!      M is the number of rows of the matrix U.
!  N - INTEGER.
!      N is the number of rows of the matrix V.
!  I - INTEGER.
!      Either, I is the index of the negligible diagonal entry Q(I)
!      of the bidiagonal J, i,e. abs(Q(I)) <= TOL, I <= p.
!      Or, I = M + 1 if E(M+1) is to be annihilated.
!  K - INTEGER.
!      Either, K is the index of the last diagonal entry of the considered
!      subbidiagonal of J, i.e. abs(E(K+1)) <= TOL, K <= p.
!      Or, K = M + 1 if E(M+1) is to be annihilated.
!  TOL - REAL (dp).
!      Specifies that matrix elements Q(i), which are <= TOL in
!      absolute value, are considered to be zero.
!  WANTU - LOGICAL.
!      Logical indicating the need for postmultiplying the Givens
!      rotations S on the left into U.
!  WANTV - LOGICAL.
!      Logical indicating the need for postmultiplying the Givens
!      rotations T on the right into V.

!  EXTERNAL SUBROUTINES and FUNCTIONS:

!  DROT from BLAS.

!  METHOD DESCRIPTION:

!  Let the considered subbidiagonal be

!            !Q(1) E(2)  0                    ...   0  !
!            ! 0   Q(2) E(3)                  ...      !
!            ! .                              ...      !
!            !             Q(i-1) E(i)              .  !
!       Jk = !                    Q(i) E(i+1)       .  !
!            !                         Q(i+1) .        !
!            ! .                              ..       !
!            ! .                                   E(k)!
!            ! 0    ...                       ...  Q(k)!

!  A zero singular value of Jk manifests itself by a zero diagonal
!  entry Q(i) or in practice, a negligible value of Q(i).
!  We call Q(i) negligible if abs(Q(i)) <= TOL.
!  When such a negligible diagonal element Q(i) in Jk is present,
!  the subbidiagonal Jk is splitted by the routine CANCEL into 2 or
!  3 unreduced subbidiagonals by annihilating E(i+1) (if i<k) using
!  Givens rotations S on the left and by annihilating E(i) (if i>1)
!  using Givens rotations T on the right until Jk is reduced to the form :

!            !Q(1) E(2)  0                ...   0  !
!            ! 0         .                ...      !
!            ! .                          ...      !
!            !         Q(i-1) 0                 .  !
!  S' Jk T = !                0   0             .  !
!            !                   Q(i+1)   .        !
!            ! .                          ..       !
!            ! .                               E(k)!
!            ! 0    ...                   ...  Q(k)!

!  For more details, see [1, pp.11.12-11.14].
!  The case of the annihilation of E(M+1) can be treated by the same
!  process. This may be seen by augmenting the matrix J with an extra
!  row of zeros, i.e. by introducing Q(M+1) = 0.

!  REFERENCES:

!  [1] J.J. Dongarra, J.R. Bunch, C.B. Moler and G.W. Stewart,
!      LINPACK User's Guide. SIAM, Philadelphia (1979).

!  CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!  REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN OUT)  :: u(:,:)   ! u(ldu,*)
REAL (dp), INTENT(IN OUT)  :: v(:,:)   ! v(ldv,*)
REAL (dp), INTENT(IN OUT)  :: q(:)
REAL (dp), INTENT(IN OUT)  :: e(:)
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: i
INTEGER, INTENT(IN)        :: k
REAL (dp), INTENT(IN)      :: tol
LOGICAL, INTENT(IN)        :: wantu
LOGICAL, INTENT(IN)        :: wantv

!     .. External Subroutines/Functions ..
! EXTERNAL drot
!     .. Intrinsic Functions ..
! INTRINSIC ABS, SQRT
!     .. Local Scalars ..
INTEGER    :: i1, l, l1
REAL (dp)  :: c, s, f, g, h
!     .. Executable Statements ..

IF (i <= m) q(i) = 0.0D0

!     Annihilate E(I+1) (if I < K).

IF (i < k) THEN
  c = 0.0D0
  s = 1.0D0
  i1 = i + 1
  DO  l = i1, k
    g = e(l)
    f = s * g
    e(l) = c * g
    IF (ABS(f) <= tol) GO TO 20
    g = q(l)
    h = SQRT(f**2 + g**2)
    q(l) = h
    c = g / h
    s = -f / h
    IF (wantu) CALL drot(m,u(1:,i),1,u(1:,l),1,c,s)
  END DO
END IF

!     Annihilate E(I) (if I > 1).

20 IF (i > 1) THEN
  i1 = i - 1
  f = e(i)
  e(i) = 0.0D0
  DO  l1 = 1, i1
    IF (ABS(f) <= tol) RETURN
    l = i - l1
    g = q(l)
    IF (ABS(g) <= tol) THEN
      g = 0.0D0
      h = ABS(f)
    ELSE
      h = SQRT(f**2 + g**2)
    END IF
    q(l) = h
    c = g / h
    s = -f / h
    g = e(l)
    f = s * g
    e(l) = c * g
    IF (wantv) CALL drot(n,v(1:,i),1,v(1:,l),1,c,s)
  END DO
  e(1) = 0.0D0
END IF
RETURN
! *** Last line of CANCEL *********************************************
END SUBROUTINE cancel



FUNCTION damin(x,nx,incx) RESULT(fn_val)
 
!  PURPOSE:

!  The function DAMIN computes the absolute minimal value of NX
!  elements in an array.
!  The function returns the value zero if NX < 1.

!  ARGUMENT LIST:

!  X - REAL (dp) array of DIMENSION (NX x INCX).
!      X is the one-dimensional array of which the absolute minimal
!      value of the elements is to be computed.
!  NX - INTEGER.
!      NX is the number of elements in X to be examined.
!  INCX - INTEGER.
!      INCX is the increment to be taken in the array X, defining
!      the distance between two consecutive elements.
!      INCX = 1, if all elements are contiguous in memory.

!  CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!  REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN)  :: x(:)   ! x(nx*incx)
INTEGER, INTENT(IN)    :: nx
INTEGER, INTENT(IN)    :: incx
REAL (dp)              :: fn_val

!     .. Intrinsic Functions ..
! INTRINSIC ABS
!     .. Local Scalars ..
INTEGER    :: i, ix, nincx
REAL (dp)  :: dx
!     .. Executable Statements ..

fn_val = 0.0D0
IF (nx < 1) RETURN

fn_val = ABS(x(1))
IF (nx == 1) RETURN

ix = 1 + incx
nincx = nx * incx
DO  i = ix, nincx, incx
  dx = ABS(x(i))
  IF (dx < fn_val) fn_val = dx
END DO
RETURN
! *** Last line of DAMIN **********************************************
END FUNCTION damin



SUBROUTINE estim(q,e,n,l,theta,tol1,tol2,iwarn)
 
!     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!     REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN)      :: q(:)
REAL (dp), INTENT(IN)      :: e(:)
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN OUT)    :: l
REAL (dp), INTENT(IN OUT)  :: theta
REAL (dp), INTENT(IN)      :: tol1
REAL (dp), INTENT(IN)      :: tol2
INTEGER, INTENT(OUT)       :: iwarn

!     .. External Subroutines/Functions ..
! INTEGER :: nsingv
! REAL (dp) :: damin
! EXTERNAL damin, nsingv
!     .. Intrinsic Functions ..
! INTRINSIC ABS, MAX
!     .. Local Scalars ..
INTEGER    :: num, i, numz
REAL (dp)  :: th, h1, h2, y, z, sumz
!     .. Executable Statements ..

iwarn = 0
IF (l < 0 .OR. l > n) RETURN

!     Step 1: Initialization of THETA.
!             -----------------------
IF (l == 0) theta = 0.0D0
IF (theta < 0.0D0) THEN
  IF (l == 1) THEN
    
!           An upper bound which is close if S(N-1) >> S(N):
    
    theta = damin(q,n,1)
  ELSE
    
!           An experimentally established estimate which is good if
!           S(N-L) >> S(N-L+1):
    
    theta = ABS(q(n-l+1))
  END IF
END IF

!     Step 2: Check quality initial estimate THETA.
!             ------------------------------------
num = nsingv(q,e,n,theta,tol1,tol2)
IF (num == l) RETURN

!     Step 3: Initialization starting values for bisection method.
!             ---------------------------------------------------
!     Let S(i), i=1,...,N, be the singular values of J in decreasing order.
!     Then, the computed Y and Z will be such that
!        (number of S(i) <= Y) < L < (number of S(i) <= Z).

IF (num < l) THEN
  y = theta
  th = ABS(q(1))
  z = th
  numz = n
  DO  i = 2, n
    h1 = ABS(e(i))
    h2 = ABS(q(i))
    sumz = MAX(th+h1, h2+h1)
    IF (sumz > z) z = sumz
    th = h2
  END DO
ELSE
  z = theta
  y = 0.0D0
  numz = num
END IF

!     Step 4: Bisection method for finding the upper bound on the L
!             smallest singular values of the bidiagonal.
!             ------------------------------------------
!     A sequence of subintervals [Y,Z] is produced such that
!            (number of S(i) <= Y) < L < (number of S(i) <= Z).
!     NUM : number of S(i) <= TH,
!     NUMZ: number of S(i) <= Z.

!     WHILE ((NUM .NE. L) .AND. (Z-Y) > TOL1) DO
20 IF (num /= l .AND. (z-y) > tol1) THEN
  th = (y+z) / 2.0D0
  num = nsingv(q,e,n,th,tol1,tol2)
  IF (num < l) THEN
    y = th
  ELSE
    z = th
    numz = num
  END IF
  GO TO 20
END IF
!     END WHILE 2

!     If (Z - Y) <= TOL1, then at least two singular values of J lie in
!     the interval [Y,Z] within a distance < TOL1 from each other.
!     S(N-L) ans S(N-L+1) are then assumed to coincide.
!     L is increased, and a warning is given.

IF (num /= l) THEN
  theta = z
  l = numz
  iwarn = 1
ELSE
  theta = th
END IF
RETURN
! *** Last line of ESTIM **********************************************
END SUBROUTINE estim



SUBROUTINE housh(dummy,k,j,tol,zero,s)
 
!  PURPOSE:

!  The subroutine HOUSH computes a Householder transformation H,
!  H = I - S x UU', that 'mirrors' a vector DUMMY(1,...,K) to the
!  J-th unit vector.

!  ARGUMENT LIST :

!  DUMMY - REAL (dp) array of DIMENSION (K).
!      A row or column vector of a matrix that has to be mirrored
!      to the corresponding unit vector EJ = (0,...,1,0,...,0).
!      On return, DUMMY contains the U-vector of the transformation
!      matrix H = I - S x UU'.
!  K - INTEGER.
!      The dimension of DUMMY.
!  J - INTEGER.
!      The transformation preserves the J-th element of DUMMY to
!      become zero.  All the other elements are transformed to zero.
!  TOL - REAL (dp).
!      If on entry, norm(DUMMY) < TOL, ZERO is put equal to .TRUE.
!  ZERO - LOGICAL.
!      See the description of TOL.
!  S - REAL (dp).
!      On return, S contains the scalar S of the transformation matrix H.

!  REFERENCES:

!  [1] A. Emami-Naeine and P. Van Dooren, Computation of Zeros of
!      Linear Multivariable Systems.
!      Automatica, 18,No.4 (1982), 415-430.

!  CONTRIBUTOR: P. Van Dooren (Philips Res. Laboratory, Brussels).

!  REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN OUT)  :: dummy(:)
INTEGER, INTENT(IN)        :: k
INTEGER, INTENT(IN)        :: j
REAL (dp), INTENT(IN)      :: tol
LOGICAL, INTENT(OUT)       :: zero
REAL (dp), INTENT(OUT)     :: s

!     .. Intrinsic Functions ..
! INTRINSIC SQRT
!     .. Local Scalars ..
REAL (dp)  :: alfa, dum1
!     .. Executable Statements ..

zero = .true.
s = SUM( dummy(1:k)**2 )
alfa = SQRT(s)
IF (alfa <= tol) RETURN
zero = .false.
dum1 = dummy(j)
IF (dum1 > 0.0D0) alfa = -alfa
dummy(j) = dum1 - alfa
s = 1.0D0 / (s - alfa*dum1)
RETURN
! *** Last line of HOUSH **********************************************
END SUBROUTINE housh



SUBROUTINE init(x,m,n)

! N.B. Argument LDX has been removed.
 
!  PURPOSE:

!  The subroutine INIT initializes an M by N matrix X with the M by N
!  identity matrix, characterized by unit diagonal entries and zero
!  off-diagonal elements.

!  ARGUMENT LIST:

!  X - REAL (dp) array of DIMENSION (LDX,N)
!      On return, X contains the M by N identity matrix.
!  LDX - INTEGER
!      LDX is the leading dimension of the array X (LDX >= M).
!  M - INTEGER
!      M is the number of rows of the matrix X.
!  N - INTEGER
!      N is the number of columns of the matrix X.

!  CONTRIBUTOR: S. Van Huffel, (ESAT Laboratory, KU Leuven).

!  REVISIONS: 1988, February 15.

REAL (dp), INTENT(OUT)  :: x(:,:)
INTEGER, INTENT(IN)     :: m
INTEGER, INTENT(IN)     :: n

!     .. Local Scalars ..
INTEGER  :: j
!     .. Executable Statements ..

x(1:m,1:n) = 0.0D0
DO j = 1, n
  x(j,j) = 1.0D0
END DO
RETURN
! *** Last line of INIT ***********************************************
END SUBROUTINE init



FUNCTION nsingv(q,e,k,theta,tol1,tol2) RESULT(fn_val)
 
!     CONTRIBUTOR: S. Van Huffel (ESAT Laboratory, KU Leuven).

!     REVISIONS: 1988, February 15.

REAL (dp), INTENT(IN)  :: q(:)
REAL (dp), INTENT(IN)  :: e(:)
INTEGER, INTENT(IN)    :: k
REAL (dp), INTENT(IN)  :: theta
REAL (dp), INTENT(IN)  :: tol1
REAL (dp), INTENT(IN)  :: tol2
INTEGER                :: fn_val

!     .. Intrinsic Functions ..
! INTRINSIC ABS
!     .. Local Scalars ..
INTEGER    :: j, numeig
REAL (dp)  :: r, t
!     .. Executable Statements ..

IF (theta < 0.0D0) THEN
  fn_val = 0
  RETURN
END IF
t = -theta - tol1
numeig = k
IF (ABS(q(1)) <= tol2) THEN
  r = t
ELSE
  r = t - q(1) * (q(1)/t)
  IF (r > 0.0D0) numeig = numeig - 1
END IF

DO  j = 2, k
  IF (ABS(e(j)) <= tol2) THEN
    r = t
  ELSE
    r = t - e(j) * (e(j)/r)
    IF (r > 0.0D0) numeig = numeig - 1
  END IF
  IF (ABS(q(j)) <= tol2) THEN
    r = t
  ELSE
    r = t - q(j) * (q(j)/r)
    IF (r > 0.0D0) numeig = numeig - 1
  END IF
END DO

fn_val = numeig
RETURN
! *** Last line of NSINGV *********************************************
END FUNCTION nsingv



SUBROUTINE dqrdc(x,n,p,qraux,jpvt,job)

! Arguments LDX & WORK have been removed.
 
REAL (dp), INTENT(IN OUT)  :: x(:,:)   ! x(ldx,1)
INTEGER, INTENT(IN)        :: n
INTEGER, INTENT(IN)        :: p
REAL (dp), INTENT(OUT)     :: qraux(:)
INTEGER, INTENT(IN OUT)    :: jpvt(:)
INTEGER, INTENT(IN)        :: job

!  dqrdc uses householder transformations to compute the qr factorization
!  of an n by p matrix x.  Column pivoting based on the 2-norms of the
!  reduced columns may be performed at the users option.

!  on entry

!     x       REAL (dp)(ldx,p), where ldx >= n.
!             x contains the matrix whose decomposition is to be computed.

!     ldx     integer.
!             ldx is the leading dimension of the array x.

!     n       integer.
!             n is the number of rows of the matrix x.

!     p       integer.
!             p is the number of columns of the matrix x.

!     jpvt    integer(p).
!             jpvt contains integers that control the selection of the
!             pivot columns.  The k-th column x(k) of x is placed in one
!             of three classes according to the value of jpvt(k).

!                if jpvt(k) > 0, then x(k) is an initial column.

!                if jpvt(k) = 0, then x(k) is a free column.

!                if jpvt(k) < 0, then x(k) is a final column.

!             Before the decomposition is computed, initial columns are moved
!             to the beginning of the array x and final columns to the end.
!             Both initial and final columns are frozen in place during the
!             computation and only free columns are moved.  At the k-th stage
!             of the reduction, if x(k) is occupied by a free column it is
!             interchanged with the free column of largest reduced norm.
!             jpvt is not referenced if job = 0.

!     work    REAL (dp)(p).
!             work is a work array.  work is not referenced if job = 0.

!     job     integer.
!             job is an integer that initiates column pivoting.
!             if job = 0, no pivoting is done.
!             if job /= 0, pivoting is done.

!  on return

!     x       x contains in its upper triangle the upper triangular matrix R
!             of the QR factorization.   Below its diagonal x contains
!             information from which the orthogonal part of the decomposition
!             can be recovered.  Note that if pivoting has been requested,
!             the decomposition is not that of the original matrix X but that
!             of X with its columns permuted as described by jpvt.

!     qraux   REAL (dp)(p).
!             qraux contains further information required to recover
!             the orthogonal part of the decomposition.

!     jpvt    jpvt(k) contains the index of the column of the original matrix
!             that has been interchanged into the k-th column,
!             if pivoting was requested.

!  linpack. this version dated 08/14/78 .
!  g.w. stewart, university of maryland, argonne national lab.

!  dqrdc uses the following functions and subprograms.

!  blas daxpy,ddot,dscal,dswap,dnrm2
!  fortran abs,max,min,sqrt

!     internal variables

INTEGER    :: j, jj, jp, l, lp1, lup, maxj, pl, pu
REAL (dp)  :: maxnrm, tt, nrmxl, t, work(p)
LOGICAL    :: negj, swapj

pl = 1
pu = 0
IF (job == 0) GO TO 60

!        Pivoting has been requested.  Rearrange the columns
!        according to jpvt.

DO  j = 1, p
  swapj = jpvt(j) > 0
  negj = jpvt(j) < 0
  jpvt(j) = j
  IF (negj) jpvt(j) = -j
  IF (swapj) THEN
    IF (j /= pl) CALL dswap(n,x(1:,pl),1,x(1:,j),1)
    jpvt(j) = jpvt(pl)
    jpvt(pl) = j
    pl = pl + 1
  END IF
END DO
pu = p
DO  jj = 1, p
  j = p - jj + 1
  IF (jpvt(j) < 0) THEN
    jpvt(j) = -jpvt(j)
    IF (j /= pu) THEN
      CALL dswap(n,x(1:,pu),1,x(1:,j),1)
      jp = jpvt(pu)
      jpvt(pu) = jpvt(j)
      jpvt(j) = jp
    END IF
    pu = pu - 1
  END IF
END DO

!     compute the norms of the free columns.

60 IF (pu < pl) GO TO 80
DO  j = pl, pu
  qraux(j) = dnrm2(n,x(1:,j),1)
  work(j) = qraux(j)
END DO

!     perform the householder reduction of x.

80 lup = MIN(n,p)
DO  l = 1, lup
  IF (l < pl .OR. l >= pu) GO TO 120
  
!           locate the column of largest norm and bring it
!           into the pivot position.
  
  maxnrm = 0.0D0
  maxj = l
  DO  j = l, pu
    IF (qraux(j) > maxnrm) THEN
      maxnrm = qraux(j)
      maxj = j
    END IF
  END DO
  IF (maxj /= l) THEN
    CALL dswap(n,x(1:,l),1,x(1:,maxj),1)
    qraux(maxj) = qraux(l)
    work(maxj) = work(l)
    jp = jpvt(maxj)
    jpvt(maxj) = jpvt(l)
    jpvt(l) = jp
  END IF

  120 qraux(l) = 0.0D0
  IF (l == n) CYCLE
  
!           compute the householder transformation for column l.
  
  nrmxl = dnrm2(n-l+1,x(l:,l),1)
  IF (nrmxl == 0.0D0) CYCLE
  IF (x(l,l) /= 0.0D0) nrmxl = SIGN(nrmxl,x(l,l))
  x(l:n,l) = x(l:n,l) / nrmxl
  x(l,l) = 1.0D0 + x(l,l)
  
!              apply the transformation to the remaining columns,
!              updating the norms.
  
  lp1 = l + 1
  DO  j = lp1, p
    t = -DOT_PRODUCT( x(l:n,l), x(l:n,j) ) / x(l,l)
    x(l:n,j) = x(l:n,j) + t * x(l:n,l)
    IF (j < pl .OR. j > pu) CYCLE
    IF (qraux(j) == 0.0D0) CYCLE
    tt = 1.0D0 - (ABS(x(l,j))/qraux(j))**2
    tt = MAX(tt,0.0D0)
    t = tt
    tt = 1.0D0 + 0.05D0*tt*(qraux(j)/work(j))**2
    IF (tt /= 1.0D0) THEN
      qraux(j) = qraux(j)*SQRT(t)
    ELSE
      qraux(j) = dnrm2(n-l, x(l+1:,j), 1)
      work(j) = qraux(j)
    END IF
  END DO
  
!              save the transformation.
  
  qraux(l) = x(l,l)
  x(l,l) = -nrmxl
END DO

RETURN
END SUBROUTINE dqrdc



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
INTEGER               :: i, ix, j, next
REAL (dp)             :: cuthi, cutlo, hitest, sum, xmax
REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

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
IF( i <= n ) GO TO 20

!              end of main loop.

!              compute square root and adjust for scaling.

fn_val = xmax * SQRT(sum)

RETURN
END FUNCTION dnrm2



SUBROUTINE drot (n, x, incx, y, incy, c, s)

!     applies a plane rotation.
!     jack dongarra, linpack, 3/11/78.

INTEGER, INTENT(IN)       :: n, incx, incy
REAL (dp), INTENT(IN OUT) :: x(:), y(:)
REAL (dp), INTENT(IN)     :: c, s

! Local variables
REAL (dp) :: dtemp
INTEGER   :: i, ix, iy

IF(n <= 0) RETURN
IF(incx == 1 .AND. incy == 1) GO TO 20

!       code for unequal increments or equal increments not equal
!         to 1

ix = 1
iy = 1
IF(incx < 0) ix = (-n+1)*incx + 1
IF(incy < 0) iy = (-n+1)*incy + 1
DO i = 1,n
  dtemp = c*x(ix) + s*y(iy)
  y(iy) = c*y(iy) - s*x(ix)
  x(ix) = dtemp
  ix = ix + incx
  iy = iy + incy
END DO
RETURN

!       code for both increments equal to 1

20 DO i = 1, n
  dtemp = c*x(i) + s*y(i)
  y(i) = c*y(i) - s*x(i)
  x(i) = dtemp
END DO
RETURN
END SUBROUTINE drot



SUBROUTINE dswap (n, x, incx, y, incy)

!     interchanges two vectors.

INTEGER, INTENT(IN)       :: n, incx, incy
REAL (dp), INTENT(IN OUT) :: x(:), y(:)

! Local variables
REAL (dp) :: temp(n)

IF(n <= 0) RETURN
IF(incx == 1 .AND. incy == 1) THEN
  temp = x(:n)
  x(:n) = y(:n)
  y(:n) = temp
  RETURN
END IF

temp = x(:n*incx:incx)
x(:n*incx:incx) = y(:n*incy:incy)
y(:n*incy:incy) = temp

RETURN
END SUBROUTINE dswap

END MODULE Total_LS
