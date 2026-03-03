SUBROUTINE syminv(a, n, nn, c, nullty, ifault)

!    Algorithm AS7, Applied Statistics, vol.17, 1968, p.198.

! N.B. Argument W has been removed.

!    Forms in c( ) as lower triangle, a generalised inverse
!    of the positive semi-definite symmetric matrix a( )
!    order n, stored as lower triangle.

!    arguments:-
!    a()     = input, the symmetric matrix to be inverted, stored in
!              lower triangular form
!    n       = input, order of the matrix
!    nn      = input, the size of the a and c arrays >=  n*(n+1)/2
!    c()     = output, the inverse of a (a generalized inverse if c is
!              singular), also stored in lower triangular.
!              c and a may occupy the same locations.
!    nullty  = output, the rank deficiency of a.
!    ifault  = output, error indicator
!                    = 1 if n < 1
!                    = 2 if a is not +ve semi-definite
!                    = 3 if nn < n*(n+1)/2
!                    = 0 otherwise

!***************************************************************************

IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
REAL (dp), INTENT(IN)  :: a(:)
INTEGER, INTENT(IN)    :: n, nn
REAL (dp), INTENT(OUT) :: c(:)
INTEGER, INTENT(OUT)   :: nullty, ifault

INTERFACE
  SUBROUTINE chol (a, n, nn, u, nullty, ifault)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
    REAL (dp), INTENT(IN)  :: a(:)
    INTEGER, INTENT(IN)    :: n, nn
    REAL (dp), INTENT(OUT) :: u(:)
    INTEGER, INTENT(OUT)   :: nullty, ifault
  END SUBROUTINE chol
END INTERFACE

! Local variables
REAL (dp)              :: w(n), x
REAL (dp), PARAMETER   :: zero = 0.0_dp, one = 1.0_dp
INTEGER                :: i, icol, irow, j, jcol, k, l, mdiag, ndiag

!       Cholesky factorization of a, result in c

CALL chol(a, n, nn, c, nullty, ifault)
IF (ifault /= 0) RETURN

!       Invert c & form the product (cinv)'*cinv, where cinv is the inverse
!       of c, row by row starting with the last row.
!       irow = the row number, ndiag = location of last element in the row.

ndiag = nn
DO irow = n, 1, -1
  l = ndiag
  IF (c(ndiag) == zero) THEN
    DO j = irow, n
      c(l) = zero
      l = l + j
    END DO
    GO TO 80
  END IF

  DO i = irow, n
    w(i) = c(l)
    l = l + i
  END DO

  icol = n
  jcol = nn
  mdiag = nn
  DO
    l = jcol
    x = zero
    IF(icol == irow) x = one / w(irow)
    k = n
    DO
      IF(k == irow) EXIT
      x = x - w(k)*c(l)
      k = k - 1
      l = l - 1
      IF(l > mdiag) l = l - k + 1
    END DO

    c(l) = x/w(irow)
    IF(icol == irow) EXIT

    mdiag = mdiag - icol
    icol = icol - 1
    jcol = jcol - 1
  END DO

  80 ndiag = ndiag - irow
END DO

RETURN
END SUBROUTINE syminv




SUBROUTINE subinv(a, nm, b, n, c, nullty, ifault, det)

!       Remark asr 44.1   Appl. Statist. (1982) Vol.31, No.3

!       A revised and enhanced version of
!       algorithm as7, applied statistics, vol.17, 1968.

! N.B. Arguments W & NDIM have been removed, and DET is OPTIONAL.

!       Forms in c() as lower triangle, a generalised inverse
!       of a sub-matrix of the positive semi-definite symmetric
!       matrix a() of order n, stored as lower triangle.
!       c() may co-incide with a().  nullty is returned as the nullity
!       of a().  ifault is returned as 1 if n.lt.1, 2 if the
!       submatrix is not positive semi-definite, 3 if the elements
!       of b are inadmissible, otherwise zero.

!       arguments:-
!       a()     = input, the symmetric matrix to be inverted, stored in
!                 lower triangular form
!       nm      = input, the order of a
!       n       = input, order of the sub-matrix to be inverted
!       b       = input, an array containing the numbers of rows and
!                 columns of a() to be included.
!       c()     = output, the inverse of a (a generalized inverse if c is
!                 singular), also stored in lower triangular.
!                 c and a may occupy the same locations.
!       nullty  = output, the rank deficiency of a.
!       ifault  = output, error indicator
!                       = 1 if n < 1
!                       = 2 if the submatrix of a is not +ve semi-definite
!                       = 3 if elements of b are inadmissible
!                       = 0 otherwise

!     Auxiliary routine required: SUBCHL from ASR44 (see AS6)

!***************************************************************************

IMPLICIT NONE
INTEGER, PARAMETER               :: dp = SELECTED_REAL_KIND(15, 60)
REAL (dp), INTENT(IN)            :: a(:)
INTEGER, INTENT(IN)              :: nm, b(:), n
REAL (dp), INTENT(OUT)           :: c(:)
REAL (dp), INTENT(OUT), OPTIONAL :: det
INTEGER, INTENT(OUT)             :: nullty, ifault

INTERFACE
  SUBROUTINE subchl (a, b, n, u, nullty, ifault, det)
    IMPLICIT NONE
    INTEGER, PARAMETER               :: dp = SELECTED_REAL_KIND(15, 60)
    REAL (dp), INTENT(IN)            :: a(:)
    REAL (dp), INTENT(OUT)           :: u(:)
    REAL (dp), INTENT(OUT), OPTIONAL :: det
    INTEGER, INTENT(IN)              :: b(:), n
    INTEGER, INTENT(OUT)             :: nullty, ifault
  END SUBROUTINE subchl
END INTERFACE

! Local variables
REAL (dp)              :: w(n), x
REAL (dp), PARAMETER   :: zero = 0.0_dp, one = 1.0_dp
INTEGER                :: i, icol, irow, j, jcol, k, l, mdiag, ndiag, nn, nrow

ifault = 3
IF (n > nm) RETURN
IF (b(1) < 1 .OR. b(1) > nm - n + 1) RETURN
IF (n == 1) GO TO 19
DO i = 2, n
  IF (b(i) <= b(i-1) .OR. b(i) > nm - n + i) RETURN
END DO
19 nrow = n
ifault = 1
IF(nrow <= 0) GO TO 100
ifault = 0

!       Cholesky factorization of a, result in c

IF (PRESENT(det)) THEN
  CALL subchl(a, b, nrow, c, nullty, ifault, det)
ELSE
  CALL subchl(a, b, nrow, c, nullty, ifault )
END IF
IF(ifault /= 0) GO TO 100

!       invert c & form the product (cinv)'*cinv, where cinv is the inverse
!       of c, row by row starting with the last row.
!       irow = the row number, ndiag = location of last element in the row.

nn = nrow*(nrow + 1)/2
irow = nrow
ndiag = nn
DO irow = nrow, 1, -1
  IF(c(ndiag) == zero) GO TO 60
  l = ndiag
  DO i = irow, nrow
    w(i) = c(l)
    l = l + i
  END DO
  icol = nrow
  jcol = nn
  mdiag = nn

  DO
    l = jcol
    x = zero
    IF(icol == irow) x = one / w(irow)
    k = nrow
    DO
      IF(k == irow) EXIT
      x = x - w(k)*c(l)
      k = k - 1
      l = l - 1
      IF(l > mdiag) l = l - k + 1
    END DO

    c(l) = x / w(irow)
    IF(icol == irow) GO TO 80
    mdiag = mdiag - icol
    icol = icol - 1
    jcol = jcol - 1
  END DO

  60 l = ndiag
  DO j = irow, nrow
    c(l) = zero
    l = l + j
  END DO

  80 ndiag = ndiag - irow
END DO

100 RETURN
END SUBROUTINE subinv



PROGRAM test_as6_and_as7
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
REAL (dp), ALLOCATABLE :: a(:), c(:), aa(:), u(:)
INTEGER, ALLOCATABLE   :: b(:)
INTEGER                :: n, nn, nullty, ifault, status, row, col, pos1, pos2, &
                          i, k, end_row, destn
REAL (dp)              :: total, det, maxerr

INTERFACE
  SUBROUTINE syminv(a, n, nn, c, nullty, ifault)
    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
    REAL (dp), INTENT(IN)  :: a(:)
    INTEGER, INTENT(IN)    :: n, nn
    REAL (dp), INTENT(OUT) :: c(:)
    INTEGER, INTENT(OUT)   :: nullty, ifault
  END SUBROUTINE syminv

  SUBROUTINE subinv(a, nm, b, n, c, nullty, ifault, det)
    IMPLICIT NONE
    INTEGER, PARAMETER               :: dp = SELECTED_REAL_KIND(15, 60)
    REAL (dp), INTENT(IN)            :: a(:)
    INTEGER, INTENT(IN)              :: nm, b(:), n
    REAL (dp), INTENT(OUT)           :: c(:)
    REAL (dp), INTENT(OUT), OPTIONAL :: det
    INTEGER, INTENT(OUT)             :: nullty, ifault
  END SUBROUTINE subinv
END INTERFACE

WRITE(*, *) 'Enter number of rows for matrix: '
READ(*, *) n
nn = n*(n+1)/2
ALLOCATE( a(nn), c(nn), aa(nn), u(nn), b(n), STAT=status )
IF (status /= 0) THEN
  WRITE(*, *) 'Unable to allocate space for arrays'
  STOP
END IF

! Generate random lower-triangular matrix.
CALL RANDOM_NUMBER( c )

! aa = c.c'
! This code could be used if aa & c occupied the same locations.

det = 1.0_dp
end_row = nn
destn = nn
DO row = n, 1, -1
  pos1 = end_row
  DO col = row, 1, -1
    pos2 = destn
    IF (col == row) det = det * c(destn)
    total = c(pos1)*c(pos2)
    DO k = col-1, 1, -1
      pos1 = pos1 - 1
      pos2 = pos2 - 1
      total = total + c(pos1)*c(pos2)
    END DO
    aa(destn) = total
    destn = destn - 1
    pos1 = pos1 - 1
  END DO
  end_row = destn
END DO

WRITE(*, '(a, g13.5)') ' Determinant of original Cholesky factor = ', det
WRITE(*, '(a, g13.5)') ' Expected determinant of inv(A''A)        = ',  &
                       1.0_dp / det**2

a = aa

! Form the Cholesky factorization, then invert & calculate inv(A).
CALL syminv(a, n, nn, u, nullty, ifault)
IF (nullty /= 0) WRITE(*, *) 'NULLTY = ', nullty
IF (ifault /= 0) WRITE(*, *) 'IFAULT = ', ifault

! Reverse the process
DO i = 1, n
  b(i) = i
END DO
! Apart from signs, u should be identical to the original c matrix.
CALL subinv(u, n, b, n, a, nullty, ifault, det)
IF (nullty /= 0) WRITE(*, *) 'NULLTY = ', nullty
IF (ifault /= 0) WRITE(*, *) 'IFAULT = ', ifault
WRITE(*, '(a, g13.5)') ' Determinant now = ', det

! Calculate maximum absolute difference between elements of A & AA.

maxerr = 0.0_dp
DO i = 1, nn
  maxerr = MAX( maxerr, ABS( a(i) - aa(i) ) )
END DO
WRITE(*, '(a, g12.4)') ' Max. abs. error = ', maxerr

DEALLOCATE( a, c, aa, u, b )
STOP
END PROGRAM test_as6_and_as7
