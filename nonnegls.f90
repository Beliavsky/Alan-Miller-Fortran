SUBROUTINE nonnegls(nin, b)
! Non-negative least squares fitting after INCLUD has been used to form
! an orthogonal reduction of the data.

! Author: Alan.Miller @ vic.cmis.csiro.au
! Latest version - 19 October 1998

USE lsq
IMPLICIT NONE

INTEGER, INTENT(OUT)         :: nin
REAL (lsq_kind), INTENT(OUT) :: b(:)

! Local variables

REAL (lsq_kind)            :: alpha, sxy(ncol), sxx(ncol), temp, wmax, z(ncol)
INTEGER                    :: i, ier, t
REAL (lsq_kind), PARAMETER :: zero = 0.0_lsq_kind

! The step numbers below are those for Lawson & Hanson's NNLS algorithm
! from page 161 of their book:
! Lawson, C.L. & Hanson, R.J. `Solving least squares problems', Prentice-Hall,
! 1974, republished in SIAM Classics in Applied Mathematics 1995.

! Step 1

nin = 0
b(:ncol) = zero

! Step 2

2 CALL sums_of_sq_and_prod(nin+1, ncol, sxy, sxx)

! Step 4
! We find the max of (sxy)^2 / sxx conditional upon sxy > 0.
! This is better than L & H's (sxy - b'.sxx) which is dependent upon the size
! of the X-variables and so tends to pick X-variables which have large values
! ahead of X-variables with small values.

t = 0
wmax = zero
DO i = nin+1, ncol
  IF (sxy(i) > zero) THEN
    temp = sxy(i)**2 / sxx(i)
    IF (temp > wmax) THEN
      wmax = temp
      t = i
    END IF
  END IF
END DO

! Step 3

IF (t == 0) RETURN

! Step 5, add the variable from position t to the model.

nin = nin + 1
IF (t > nin) CALL vmove(t, nin, ier)

! Step 6, calculate the regression coefficients.

6 CALL regcf(z, nin, ier)

! Step 7, if all coefficients > 0, set b = z, then go back to step 2, unless
! all the variables have been included.

IF (ALL(z(:nin) > zero)) THEN
  b(:nin) = z(:nin)
  IF (nin == ncol) RETURN
  GO TO 2
END IF

! Steps 8-11, one or more negative regression coefficients.
! Steps 8-9, calculate shrinkage factor alpha.

t = 0
alpha = 1.0
nin = nin - 1
DO i = 1, nin
  IF (z(i) <= zero) THEN
    temp = b(i) / (b(i) - z(i))
    IF (temp < alpha) THEN
      alpha = temp
      t = i
    END IF
  END IF
END DO

! Step 10, shrink regression coefficients

IF (t > 0) THEN
  b(1:nin) = b(1:nin) + alpha * (z(1:nin) - b(1:nin))
ELSE
  GO TO 2                    ! This should never happen
END IF

! Step 11, remove variables with coefficients <= 0.

CALL vmove(t, nin, ier)
nin = nin - 1
GO TO 6


CONTAINS


SUBROUTINE sums_of_sq_and_prod(first, last, sxy, sxx)
! Calculate the sums of squares, and cross-products, of those parts of the
! X and Y variables which are orthogonal to the X-variables already in the
! model (if any) in positions 1 .. first-1.  Ignore any variables after last.

INTEGER, INTENT(IN)          :: first, last
REAL (lsq_kind), INTENT(OUT) :: sxy(:), sxx(:)

!     Local variables

INTEGER         :: inc, pos, row, col
REAL (lsq_kind) :: zero = 0.0, diag, dy

!     Accumulate sums of squares & products from row FIRST.
!     The orthogonal reduction gave:
!     X = Q.sqrt(D).R   where Q is orthonormal and D is diagonal and such
!     that the diagonal elements of R are 1's.
!     The sums of squares and products are given by  R'DR
!     where only rows & columns FIRST .. LAST of R & D are used.

sxx(first:last) = zero
sxy(first:last) = zero
inc = ncol - last
pos = row_ptr(first)
DO row = first, last
  diag = d(row)
  dy = diag * rhs(row)
  sxx(row) = sxx(row) + diag
  sxy(row) = sxy(row) + dy
  DO col = row+1, last
    sxx(col) = sxx(col) + diag * r(pos)**2
    sxy(col) = sxy(col) + dy * r(pos)
    pos = pos + 1
  END DO
  pos = pos + inc
END DO

RETURN
END SUBROUTINE sums_of_sq_and_prod

END SUBROUTINE nonnegls
