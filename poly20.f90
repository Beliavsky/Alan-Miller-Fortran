PROGRAM poly20
! Generate the coefficients of the polynomial:
!    (x-1)(x-2)(x-3) ... (x-20)
! and then solve it.

USE poly_zeroes
IMPLICIT NONE

INTEGER           :: i, j, n, nitmax, iter
REAL (KIND=dp)    :: p(0:20), new_root, eps, big, small, radius(20), pr
COMPLEX (KIND=dp) :: root(20), poly(0:20)
LOGICAL           :: err(21)

INTERFACE
  SUBROUTINE sort(n, x, y, e)
    USE poly_zeroes
    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: n
    COMPLEX (KIND=dp), INTENT(IN OUT) :: x(:)
    REAL (KIND=dp), INTENT(IN OUT)    :: y(:)
    LOGICAL, INTENT(IN OUT)           :: e(:)
  END SUBROUTINE sort
END INTERFACE

eps    = EPSILON(1.0D0)
small  = TINY(1.0D0)
big    = HUGE(1.0D0)
nitmax = 30

p(0) = -1.0_dp
p(1) =  1.0_dp
n = 20
DO i = 2, n
  new_root = DBLE(i)
  p(i) = p(i-1)
  DO j = i-1, 1, -1
    p(j) = p(j-1) - new_root * p(j)
  END DO
  p(0) = - new_root * p(0)
END DO

poly = CMPLX(p, KIND=dp)

CALL polzeros (n, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL sort(n, root, radius, err)
OPEN(UNIT=2, FILE='poly20.out', STATUS='UNKNOWN', ACTION='WRITE')
DO i = 1, n
  pr = ABS(root(i))
  IF(pr /= 0.) pr=radius(i)/pr
  IF(radius(i) <= 0) pr = -1
  WRITE(*,300) i, err(i), root(i), radius(i), pr
  WRITE(2,300) i, err(i), root(i), radius(i), pr
END DO
WRITE(*,*)' ITERATIONS =', iter
WRITE(2,*)' ITERATIONS =', iter

STOP
300 FORMAT(' ',i4,' ERR=',l1,' X=',e21.15,',',e21.15,' R=',e7.1,' RR=',e7.1)

END PROGRAM poly20



!************************************************************************
!                             SUBROUTINE SORT                           *
!************************************************************************
!   SORT  the vector X, according to nonincreasing real parts,          *
!   the same permutation is applied to vectors Y and E.                 *
!************************************************************************
SUBROUTINE sort(n, x, y, e)
USE poly_zeroes
IMPLICIT NONE
INTEGER, INTENT(IN)               :: n
COMPLEX (KIND=dp), INTENT(IN OUT) :: x(:)
REAL (KIND=dp), INTENT(IN OUT)    :: y(:)
LOGICAL, INTENT(IN OUT)           :: e(:)

! Local variables
INTEGER           :: k, i, imax
COMPLEX (KIND=dp) :: temp
REAL (KIND=dp)    :: yt, amax, rxi
LOGICAL           :: et

DO k = 1, n-1
  amax = REAL(x(k))
  imax = k
  DO i = k+1,n
    rxi = REAL(x(i))
    IF (amax < rxi) THEN
      amax = rxi
      imax = i
    END IF
  END DO
  temp = x(k)
  x(k) = x(imax)
  x(imax) = temp
  yt = y(k)
  et = e(k)
  y(k) = y(imax)
  y(imax) = yt
  e(k) = e(imax)
  e(imax) = et
END DO
RETURN
END SUBROUTINE sort

