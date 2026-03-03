MODULE ridge
! Regularization (ridge regression), including SVD calculation,
! using the LSQ module.   LSQ must have already been used to have formed
! a QR-factorization.   The ridge regression / regularization options are
! only usable if there are more cases than predictors, and a non-zero
! estimate of residual variance is available.

! The methods used here for regularization all use the SVD.   This is not
! the most efficient way of performing regularization, but it does preserve
! the QR-factorization for subsequent use.

! The alternative way of performing ridge regression / regularization is by
! adding extra artificial lines of data of the kind:
! Variable: X1  X2  ...  Xi  ...  Xk  Y
! Value:     0   0        1        0  0
! with weight = lambda, for each i in turn.   This has the effect of adding a
! diagonal matrix to the X'X-matrix with elements lambda on the diagonal.
! This can be done using module LSQ (or almost any other least-squares code).
! Using LSQ, it will increase the number of observations (nobs) for each row
! added.  If nobs is needed, it should be reduced accordingly.

! Author: Alan Miller
! This code is ELF90 compatible
! Latest revision - 1 August 1997

USE lsq
IMPLICIT NONE

REAL (lsq_kind), ALLOCATABLE, SAVE :: a(:,:), v(:,:), sigma(:), b(:,:)
REAL (lsq_kind), PARAMETER         :: zero = 0.0_lsq_kind, one = 1.0_lsq_kind
LOGICAL, SAVE                      :: svd_initialized = .FALSE.


CONTAINS


SUBROUTINE initialize_svd(nreq)
! Initialization prior to the calculation of the SVD for the first nreq rows
! and columns of the QR factorization.

INTEGER, INTENT(IN) :: nreq

IF (ALLOCATED( a )) DEALLOCATE( a )
IF (ALLOCATED( v )) DEALLOCATE( v )
IF (ALLOCATED( sigma )) DEALLOCATE( sigma )
IF (ALLOCATED( b )) DEALLOCATE( b )
ALLOCATE( a(nreq,nreq), v(nreq,nreq), sigma(nreq), b(nreq,1) )

IF (nreq > ncol) THEN
  WRITE(*, *) 'Error in call to initialize_svd'
  WRITE(*, '(a, i10, a, i6)') ' nreq = ', nreq, ' > ncol = ', ncol
  RETURN
END IF

CALL qr_to_2d(nreq)
svd_initialized = .TRUE.

RETURN
END SUBROUTINE initialize_svd



SUBROUTINE qr_to_2d(nreq)
! The elements of R are scaled by the row multipliers and stored in a 2D
! array (a).   The right-hand side projections are scaled and copied into b.

INTEGER, INTENT(IN) :: nreq

! Local variables
INTEGER         :: pos, inc, row, col
REAL (lsq_kind) :: scale

a = zero
pos = 1
inc = ncol - nreq
DO row = 1, nreq
  scale = SQRT( d(row) )
  b(row,1) = rhs(row) * scale
  a(row,row) = scale
  DO col = row+1, nreq
    a(row,col) = r(pos) * scale
    pos = pos + 1
  END DO
  pos = pos + inc
END DO

RETURN
END SUBROUTINE qr_to_2d



SUBROUTINE grsvd(nu, nv, nb, m, n, w, matu, u, matv, v, b, irhs, ierr)

!  This subroutine is a translation of the Algol procedure svd,
!  Num. Math. 14, 403-420 (1970) by Golub and Reinsch.
!  Handbook for automatic computation, vol ii-linear algebra, 134-151 (1971).
!
!  This subroutine determines the singular value decomposition
!       t
!  a=usv  of a real m by n rectangular matrix.  Householder
!  bidiagonalization and a variant of the qr algorithm are used.
!  grsvd assumes that a copy of the matrix a is in the array u. it
!  also assumes m >= n.  if m < n, then compute the singular
!                          t       t    t            t
!  value decomposition of a .  If a =uwv , then a=vwu .
!
!  grsvd can also be used to compute the minimal length least squares
!  solution to the overdetermined linear system a*x=b.
!
!  on input-
!
!     nu must be set to the row dimension of two-dimensional array u
!       as declared in the calling program dimension statement.
!       note that nu must be at least as large as m,
!
!     nv must be set to the row dimension of the two-dimensional
!       array parameter v as declared in the calling program
!       dimension statement.  nv must be at least as large as n,
!
!     nb must be set to the row dimension of two-dimensional array b
!       as declared in the calling program dimension statement.
!       Note that nb must be at least as large as m,
!
!     m is the number of rows of a (and u),
!
!     n is the number of columns of a (and u) and the order of v,
!
!     u contains the rectangular input matrix (a) to be decomposed,
!
!     b contains the irhs right-hand-sides of the overdetermined
!       linear system a*x=b.  if irhs > 0,  then on output,
!       the first n components of these irhs columns of b
!                     t
!       will contain u b.  Thus, to compute the minimal length least
!                                             +
!       squares solution, one must compute v*w  times the columns of
!                 +                        +
!       b, where w  is a diagonal matrix, w (i)=0 if w(i) is
!       negligible, otherwise is 1/w(i).  If irhs=0, b may coincide
!       with a or u and will not be referenced,
!
!     irhs is the number of right-hand-sides of the overdetermined
!       system a*x=b.  irhs should be set to zero if only the singular
!       value decomposition of a is desired,
!
!     matu should be set to .true. if the u matrix in the
!       decomposition is desired, and to .false. otherwise,
!
!     matv should be set to .true. if the v matrix in the
!       decomposition is desired, and to .false. otherwise.
!
!  on output-
!
!     w contains the n (non-negative) singular values of a (the
!       diagonal elements of s).  they are unordered.  if an
!       error exit is made, the singular values should be correct
!       for indices ierr+1,ierr+2,...,n,
!
!     u contains the matrix u (orthogonal column vectors) of the
!       decomposition if matu has been set to .true.  otherwise
!       u is used as a temporary array.
!       if an error exit is made, the columns of u corresponding
!       to indices of correct singular values should be correct,
!
!     v contains the matrix v (orthogonal) of the decomposition if
!       matv has been set to .true.  otherwise v is not referenced.
!       if an error exit is made, the columns of v corresponding to
!       indices of correct singular values should be correct,
!
!     ierr is set to
!       zero       for normal return,
!       k          if the k-th singular value has not been
!                  determined after 30 iterations,
!       -1         if irhs < 0 ,
!       -2         if m < n ,
!       -3         if nu < m .or. nb < m,
!       -4         if nv < n .
!
!     The original subroutine was checked by the pfort verifier
!     (ryder, b.g. 'the pfort verifier', software - practice and
!     experience, vol.4, 359-377, 1974) for adherence to a large,
!     carefully defined, portable subset of american national standar
!     fortran called pfort.
!
!     original version of this code is subroutine svd in release 2 of eispack.
!
!     modified by Tony F. Chan,
!                 Comp. Sci. Dept, Yale Univ.,
!                 Box 2158, Yale station,
!                 CT 06520, USA
!     last modified (fortran 77 version): january, 1982.
!
!     N.B. The last argument (rv1) has been dropped from the arguments;
!          the array was used as work space.
!
!  ------------------------------------------------------------------
!

IMPLICIT NONE

INTEGER, INTENT(IN)             :: nu, nv, nb, m, n, irhs
INTEGER, INTENT(OUT)            :: ierr
REAL (lsq_kind), INTENT(IN OUT) :: u(:,:), b(:,:)
REAL (lsq_kind), INTENT(OUT)    :: w(:), v(:,:)
LOGICAL, INTENT(IN)             :: matu, matv

!     Local variables
INTEGER         :: i, j, k, l, id, ii, i1, kk, k1, ll, l1, mn, its, krhs
REAL (lsq_kind) :: c, f, g, h, rv1(n), s, t, x, y, z, eps, scale

ierr = 0
IF (irhs < 0) THEN
  ierr = -1
  RETURN
END IF
IF (m < n) THEN
  ierr = -2
  RETURN
END IF
IF (nu < m .OR. nb < m) THEN
  ierr = -3
  RETURN
END IF
IF (nv < n) THEN
  ierr = -4
  RETURN
END IF

!     ********** Householder reduction to bidiagonal form **********
g = zero
scale = zero
x = zero

DO i=1,n
  l = i + 1
  rv1(i) = scale*g
  g = zero
  s = zero
  scale = zero

!     Compute left transformations that zero the subdiagonal elements
!     of the i-th column.

  DO k=i,m
    scale = scale + ABS(u(k,i))
  END DO

  IF (scale == zero) GO TO 160

  DO k=i,m
    u(k,i) = u(k,i)/scale
    s = s + u(k,i)**2
  END DO

  f = u(i,i)
  g = -SIGN(SQRT(s),f)
  h = f*g - s
  u(i,i) = f - g
  IF (i == n) GO TO 100

!     Apply left transformations to remaining columns of a.

  DO j=l,n
    s = DOT_PRODUCT( u(i:m,i), u(i:m,j) )
    f = s/h
    u(i:m,j) = u(i:m,j) + f*u(i:m,i)
  END DO
!
!      Apply left transformations to the columns of b if irhs > 0
!
  100 DO j=1,irhs
    s = DOT_PRODUCT( u(i:m,i), b(i:m,j) )
    f = s/h
    b(i:m,j) = b(i:m,j) + f*u(i:m,i)
  END DO

!     Compute right transformations.

  u(i:m,i) = scale*u(i:m,i)

  160 w(i) = scale*g
  g = zero
  s = zero
  scale = zero
  IF (i > m .OR. i == n) GO TO 250

  DO k=l,n
    scale = scale + ABS(u(i,k))
  END DO

  IF (scale == zero) GO TO 250

  DO k=l,n
    u(i,k) = u(i,k)/scale
    s = s + u(i,k)**2
  END DO

  f = u(i,l)
  g = -SIGN(SQRT(s),f)
  h = f*g - s
  u(i,l) = f - g

  rv1(l:n) = u(i,l:n)/h

  IF (i == m) GO TO 230

  DO j=l,m
    s = DOT_PRODUCT( u(j,l:n), u(i,l:n) )
    u(j,l:n) = u(j,l:n) + s*rv1(l:n)
  END DO

  230 u(i,l:n) = scale*u(i,l:n)

  250 x = MAX(x, ABS(w(i))+ABS(rv1(i)))
END DO
!     ********** accumulation of right-hand transformations **********
IF (.NOT.matv) GO TO 350
!     ********** for i=n step -1 until 1 do -- **********
DO ii=1,n
  i = n + 1 - ii
  IF (i == n) GO TO 330
  IF (g == zero) GO TO 310

!     ********** double division avoids possible underflow **********
  v(l:n,i) = (u(i,l:n) / u(i,l)) / g

  DO j=l,n
    s = DOT_PRODUCT( u(i,l:n), v(l:n,j) )
    v(l:n,j) = v(l:n,j) + s*v(l:n,i)
  END DO

  310 v(i,l:n) = zero
  v(l:n,i) = zero

  330 v(i,i) = one
  g = rv1(i)
  l = i
END DO
!     ********** accumulation of left-hand transformations **********
350 IF (.NOT.matu) GO TO 470
!     **********for i=min(m,n) step -1 until 1 do -- **********
mn = n
IF (m < n) mn = m

DO ii=1,mn
  i = mn + 1 - ii
  l = i + 1
  g = w(i)
  IF (i == n) GO TO 370

  u(i,l:n) = zero

  370 IF (g == zero) GO TO 430
  IF (i == mn) GO TO 410

  DO j=l,n
    s = DOT_PRODUCT( u(l:n,i), u(l:n,j) )
!     ********** double division avoids possible underflow **********
    f = (s/u(i,i)) / g
    u(i:m,j) = u(i:m,j) + f*u(i:m,i)
  END DO

  410 u(i:m,i) = u(i:m,i)/g

  GO TO 450

  430 u(i:m,i) = zero

  450 u(i,i) = u(i,i) + one
END DO
!     ********** diagonalization of the bidiagonal form **********
470 eps = EPSILON(one) * x
!     ********** for k=n step -1 until 1 do -- **********
DO kk=1,n
  k1 = n - kk
  k = k1 + 1
  its = 0
!     ********** test for splitting.
!                for l=k step -1 until 1 do -- **********
  480 DO ll=1,k
    l1 = k - ll
    l = l1 + 1
    IF (ABS(rv1(l)) <= eps) GO TO 550
!     ********** rv1(1) is always zero, so there is no exit
!                through the bottom of the loop **********
    IF (ABS(w(l1)) <= eps) GO TO 500
  END DO
!     ********** cancellation of rv1(l) if l greater than 1 **********
  500 c = zero
  s = one

  DO i=l,k
    f = s*rv1(i)
    rv1(i) = c*rv1(i)
    IF (ABS(f) <= eps) GO TO 550
    g = w(i)
    h = SQRT(f*f + g*g)
    w(i) = h
    c = g/h
    s = -f/h

!     apply left transformations to b if irhs > 0

    DO j=1,irhs
      y = b(l1,j)
      z = b(i,j)
      b(l1,j) = y*c + z*s
      b(i,j) = -y*s + z*c
    END DO

    IF (.NOT.matu) CYCLE

    DO j=1,m
      y = u(j,l1)
      z = u(j,i)
      u(j,l1) = y*c + z*s
      u(j,i) = -y*s + z*c
    END DO

  END DO
!     ********** test for convergence **********
  550 z = w(k)
  IF (l == k) GO TO 630
!     ********** shift from bottom 2 by 2 minor **********
  IF (its == 30) GO TO 660
  its = its + 1
  x = w(l)
  y = w(k1)
  g = rv1(k1)
  h = rv1(k)
  f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2.0*h*y)
  g = SQRT(f*f + one)
  f = ((x-z)*(x+z) + h*(y/(f+SIGN(g,f))-h))/x
!     ********** next qr transformation **********
  c = one
  s = one

  DO i1=l,k1
    i = i1 + 1
    g = rv1(i)
    y = w(i)
    h = s*g
    g = c*g
    z = SQRT(f*f+h*h)
    rv1(i1) = z
    c = f/z
    s = h/z
    f = x*c + g*s
    g = -x*s + g*c
    h = y*s
    y = y*c
    IF (.NOT.matv) GO TO 570

    DO j=1,n
      x = v(j,i1)
      z = v(j,i)
      v(j,i1) = x*c + z*s
      v(j,i) = -x*s + z*c
    END DO

    570 z = SQRT(f*f+h*h)
    w(i1) = z
!     ********** rotation can be arbitrary if z is zero **********
    IF (z == zero) GO TO 580
    c = f/z
    s = h/z
    580 f = c*g + s*y
    x = -s*g + c*y

!     apply left transformations to b if irhs > 0

    DO j=1, irhs
      y = b(i1,j)
      z = b(i,j)
      b(i1,j) = y*c + z*s
      b(i,j) = -y*s + z*c
    END DO

    IF (.NOT.matu) CYCLE

    DO j=1,m
      y = u(j,i1)
      z = u(j,i)
      u(j,i1) = y*c + z*s
      u(j,i) = -y*s + z*c
    END DO

  END DO

  rv1(l) = zero
  rv1(k) = f
  w(k) = x
  GO TO 480
!     ********** convergence **********
  630 IF (z >= zero) CYCLE
!     ********** w(k) is made non-negative **********
  w(k) = -z
  IF (.NOT.matv) CYCLE

  v(1:n,k) = -v(1:n,k)

END DO

! Re-order the singular values.
! Code from TOMS algorithm 581 by Tony Chan

!     Sort singular values and exchange columns of u and v accordingly.
!     Selection sort minimizes swapping of u and v.

DO i=1, n-1
!...    find index of maximum singular value
  id = i
  DO j=i+1, n
    IF (w(j) > w(id)) id = j
  END DO
  IF (id == i) CYCLE
!...    swap singular values and vectors
  t = w(i)
  w(i) = w(id)
  w(id) = t
  IF (matv) CALL dswap(n, v(1:,i), 1, v(1:,id), 1)
  IF (matu) CALL dswap(m, u(1:,i), 1, u(1:,id), 1)
  IF (irhs < 1) CYCLE
  DO krhs=1, irhs
    t = b(i,krhs)
    b(i,krhs) = b(id,krhs)
    b(id,krhs) = t
  END DO
END DO

RETURN
!     ********** set error -- no convergence to a
!                singular value after 30 iterations **********
660 ierr = k
RETURN
!     ********** last card of grsvd **********

END SUBROUTINE grsvd


SUBROUTINE dswap(n, dx, incx, dy, incy)

!     interchanges two vectors.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.

INTEGER, PARAMETER        :: dp = SELECTED_REAL_KIND(14, 100)
REAL (dp), INTENT(IN OUT) :: dx(:), dy(:)
INTEGER, INTENT(IN)       :: incx, incy, n

!     Local variables
INTEGER   :: i, ix, iy, m, mp1
REAL (dp) :: dtemp

IF (n <= 0) RETURN
IF (incx == 1 .AND. incy == 1) GO TO 20

!       code for unequal increments or equal increments not equal to 1

ix = 1
iy = 1
IF (incx < 0) ix = (-n+1)*incx + 1
IF (incy < 0) iy = (-n+1)*incy + 1
DO i=1,n
  dtemp = dx(ix)
  dx(ix) = dy(iy)
  dy(iy) = dtemp
  ix = ix + incx
  iy = iy + incy
END DO
RETURN

!       code for both increments equal to 1
!
!
!       clean-up loop

20 m = MOD(n,3)
IF (m == 0) GO TO 40
DO i = 1, m
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
END DO
IF (n < 3) RETURN
40 mp1 = m + 1
DO i = mp1,n,3
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
  dtemp = dx(i+1)
  dx(i+1) = dy(i+1)
  dy(i+1) = dtemp
  dtemp = dx(i+2)
  dx(i+2) = dy(i+2)
  dy(i+2) = dtemp
END DO

RETURN
END SUBROUTINE dswap



SUBROUTINE regularize(nreq, fit_const, option, beta, ier, lambda, col_scale, &
                      bias, covar)
! Ridge regression / regularization.

! Arguments:
! nreq      = number of variables to be used (include the constant if one is
!             being fitted).   Must be <= the total number of variables used
!             in forming the QR-factorization (ncol).   If nreq < ncol, the
!             first nreq variables, in their current order, are used.
! fit_const = .TRUE. if a constant is being fitted, in which case no shrinkage
!             is applied to the first orthogonal direction (the centroid).
!           = .FALSE. otherwise
! option    = the kind of regularization to be performed
!           = 0 if the user is supplying an array, lambda(), of elements to
!               augment the diagonal.   In this case, both of the optional
!               arguments lambda and col_scale must be supplied.
!           = 1 for the original Hoerl & Kennard method
!           = 2 for Hemmerle's iterated method
!           = 3 for Bayesian estimate
!           = 4 for a method due to Browne & Rock and also Hemmerle & Brantle
! beta      = the array of ridge regression coefficients
! ier       = error indicator
!           = 0 if successful
!           = 1 if nreq > ncol
!           = 2 if option < 0 or option > 4
!           = 3 if the dimension of beta in the calling program < nreq
!           = 4 if there are no degrees of freedom for estimating the residual
!               variance
!           = 5 if option = 0 and either lambda or col_scale is absent.
! lambda (optional)
!           = array of elements to be added to the diagonal of the X'X matrix.
! col_scale (optional)
!           = logical variable which indicates whether the elements of lambda
!             are to be scaled.   If col_scale = .TRUE., lambda(i) is multiplied
!             by the sum of squares of values of the i-th variable.   If
!             fit_const = .TRUE. as well, the mean is first subtracted from
!             each element in column i.   This is equivalent to the common
!             practice of applying ridge regression to correlation matrices,
!             that is to variables which have been standardized.
!             There is no scaling if col_scale = .FALSE.
! bias (optional)
!           = estimated vector of biases in beta introduced by regularization
! covar (optional)
!           = estimated covariance matrix for beta

! References (the numbers correspond to the option numbers):
! 1. Hoerl, A.E. & Kennard, R.W. (1970) `Ridge regression: biased estimation
!           for non-orthogonal problems', Technometrics, v.12, 55-67.
! 2. Hemmerle, W.J. (1975) `An explicit solution for generalized ridge
!           regression, Technometrics, v.17, 309-314.
! 3. Troskie, C.G. & Chalton, D.O. (1996) `A Bayesian estimate for the constants
!           in ridge regression', S.African Statist. J., v.30, 119-137.
! 4a. Browne, M.W. & Rock, D.A. (1978) `The choice of additive constants in
!           ridge regression', S.African Statist. J., v.12, 57-74.
! 4b. Hemmerle, W.J. & Brantle, T.F. (1978) 'Explicit and constrained
!           generalized ridge estimation', Technometrics, v.20, 109-120.

INTEGER, INTENT(IN)                    :: nreq, option
LOGICAL, INTENT(IN)                    :: fit_const
REAL (lsq_kind), INTENT(IN), OPTIONAL  :: lambda(:)
REAL (lsq_kind), INTENT(OUT)           :: beta(:)
INTEGER, INTENT(OUT)                   :: ier
REAL (lsq_kind), INTENT(OUT), OPTIONAL :: bias(:), covar(:,:)
LOGICAL, INTENT(IN), OPTIONAL          :: col_scale

! Local variables
INTEGER         :: irhs = 1, i1, i
REAL (lsq_kind) :: stdev, std_proj, two = 2.0_lsq_kind, half = 0.5_lsq_kind, &
                   quart = 0.25_lsq_kind
REAL (lsq_kind), ALLOCATABLE :: temp(:), shrinkage(:), col_norm(:)

IF (nreq > ncol) THEN
  ier = 1
ELSE IF (option < 0 .OR. option > 4) THEN
  ier = 2
ELSE IF (SIZE(beta) < nreq) THEN
  ier = 3
ELSE IF (nobs <= nreq) THEN
  ier = 4
ELSE IF (option == 0 .AND.   &
        (.NOT. PRESENT(lambda) .OR. .NOT. PRESENT(col_scale))) THEN
  ier = 5
ELSE 
  ier = 0
END IF
IF (ier > 0) RETURN

CALL initialize_svd(nreq)
ALLOCATE( shrinkage(nreq) )

IF (fit_const) THEN
  i1 = 2
  shrinkage(1) = one
ELSE
  i1 = 1
END IF

! If option = 0 & col_scale = .TRUE., scale the columns of the QR factorization.
! The sums of squares of column elements will be stored in array col_norm.

IF (option == 0 .AND. col_scale) THEN
  ALLOCATE( col_norm(nreq) )
  DO i = i1, nreq
    col_norm(i) = SQRT( SUM( a(i1:i,i)**2 ) )
    a(i1:i,i) = a(i1:i,i) / col_norm(i)
  END DO
END IF

CALL grsvd(nreq, nreq, nreq, nreq, nreq, sigma, .FALSE., a, .TRUE., v, b,  &
           irhs, ier)

CALL ss()
stdev = SQRT( rss(nreq) / (nobs - nreq) )

! Shrink the SVD projections in array b.
SELECT CASE (option)
  CASE (0)
    shrinkage(i1:nreq) = sigma(i1:nreq)**2 /  &
                         (sigma(i1:nreq)**2 + lambda(i1:nreq))
    b(i1:nreq,1) = b(i1:nreq,1) * shrinkage(i1:nreq)
  CASE (1)
    DO i = i1, nreq
      std_proj = b(i,1) / stdev
      shrinkage(i) = std_proj**2 / (one + std_proj**2)
      b(i, 1) = b(i,1) * shrinkage(i)
    END DO

  CASE (2)
    DO i = i1, nreq
      std_proj = b(i,1) / stdev
      IF (ABS(std_proj) < two) THEN
        shrinkage(i) = zero
        b(i,1) = zero
      ELSE
        shrinkage(i) = half + SQRT( quart - one / std_proj**2 )
        b(i,1) = b(i,1) * shrinkage(i)
      END IF
    END DO

  CASE (3)
    DO i = i1, nreq
      std_proj = b(i,1) / stdev
      shrinkage(i) = (one + std_proj**2) / (two + std_proj**2)
      b(i, 1) = b(i,1) * shrinkage(i)
    END DO

  CASE (4)
    DO i = i1, nreq
      std_proj = b(i,1) / stdev
      IF (ABS(std_proj) < one) THEN
        shrinkage(i) = zero
        b(i,1) = zero
      ELSE
        shrinkage(i) = one - one / std_proj**2
        b(i,1) = b(i,1) * shrinkage(i)
      END IF
    END DO
END SELECT

! Calculate regression coefficients from the shrunken projections
ALLOCATE( temp(nreq) )
temp = b(:,1) / sigma
beta(1:nreq) = MATMUL(v, temp)
IF (ALLOCATED(col_norm)) beta(i1:nreq) = beta(i1:nreq) / col_norm(i1:nreq)

! Estimate bias, if argument is present
IF (PRESENT(bias)) THEN
  temp = temp * (one - shrinkage)
  bias = MATMUL(v, temp)
  IF (ALLOCATED(col_norm)) bias(i1:nreq) = bias(i1:nreq) / col_norm(i1:nreq)
END IF

! Estimate covariance, if argument is present
IF (PRESENT(covar)) THEN
  covar = zero
  DO i = 1, nreq
    covar(i,i) = stdev * shrinkage(i) / sigma(i)
  END DO
  covar = MATMUL(v, covar)
  covar = MATMUL(covar, TRANSPOSE(covar))
  IF (ALLOCATED(col_norm)) THEN
    DO i = i1, nreq
      covar(i1:nreq, i) = covar(i1:nreq, i) * col_norm(i)
    END DO
    DO i = i1, nreq
      covar(i, i1:nreq) = covar(i, i1:nreq) * col_norm(i)
    END DO
  END IF
END IF

DEALLOCATE( temp, shrinkage )
IF (ALLOCATED(col_norm)) DEALLOCATE( col_norm )

RETURN
END SUBROUTINE regularize


END MODULE ridge
