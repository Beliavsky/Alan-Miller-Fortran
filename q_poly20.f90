PROGRAM poly20
! Generate the coefficients of the polynomial:
!    (x-1)(x-2)(x-3) ... (x-20)
! and then solve it.

USE quadruple_precision
USE quad_prec_complex
USE poly_zeroes
IMPLICIT NONE

INTEGER            :: i, j, n, nitmax, iter, ier
REAL (dp)          :: eps, big, small
TYPE (quad)        :: p(0:20), new_root, radius(20), pr
TYPE (qc)          :: root(20), poly(0:20)
LOGICAL            :: err(21)
CHARACTER (LEN=36) :: string

eps    = 1.d-29
small  = TINY(1.0_dp) * 1.e+16_dp
big    = HUGE(1.0_dp)
nitmax = 50

p(0) = quad(-1.0_dp, 0.0_dp)
p(1) = quad( 1.0_dp, 0.0_dp)
n = 20
DO i = 2, n
  new_root = quad( DBLE(i), 0.0_dp )
  p(i) = p(i-1)
  DO j = i-1, 1, -1
    p(j) = p(j-1) - new_root * p(j)
  END DO
  p(0) = - new_root * p(0)
END DO

OPEN(UNIT=2, FILE='poly20.out', STATUS='UNKNOWN', ACTION='WRITE')
WRITE(2, '(a)') 'Coefficients of input polynomial'
DO i = 0, n
  CALL quad_string(p(i), string, ier)
  IF (ier == 0) WRITE(2, '(i4, "  ", a36)') i, string
END DO
WRITE(2, *)

DO i = 0, n
  poly(i) = CMPLX(p(i), quad(0.0_dp, 0.0_dp) )
END DO

CALL polzeros (n, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL sort(n, root, radius, err)
WRITE(2, '(a)') 'Solution:'
WRITE(*, '(t16, a, t40, a)') ' Real part', 'Imaginary part'
WRITE(2, '(t16, a, t40, a)') ' Real part', 'Imaginary part'
DO i = 1, n
  pr = ABS(root(i))
  IF(pr%hi /= 0.0_dp) pr = radius(i) / pr
  IF(radius(i)%hi <= 0.0_dp) pr = quad( -1.0_dp, 0.0_dp )
  WRITE(*,300) i, err(i), root(i)%qr%hi, root(i)%qi%hi, radius(i)%hi, pr%hi
  WRITE(2,300) i, err(i), root(i)%qr%hi, root(i)%qi%hi, radius(i)%hi, pr%hi
END DO
WRITE(*,*)' ITERATIONS =', iter
WRITE(2,*)' ITERATIONS =', iter

STOP
300 FORMAT(' ', i3, ' ERR=', l1,' X=', e21.15, ', ', e21.15, '  R=', e7.1,  &
           ' RR=', e7.1)


CONTAINS

!************************************************************************
!                             SUBROUTINE SORT                           *
!************************************************************************
!   SORT  the vector X, according to nonincreasing real parts,          *
!   the same permutation is applied to vectors Y and E.                 *
!************************************************************************
SUBROUTINE sort(n, x, y, e)

INTEGER, INTENT(IN)         :: n
TYPE (qc), INTENT(IN OUT)   :: x(:)
TYPE (quad), INTENT(IN OUT) :: y(:)
LOGICAL, INTENT(IN OUT)     :: e(:)

! Local variables
INTEGER     :: k, i, imax
TYPE (qc)   :: temp
TYPE (quad) :: yt, amax, rxi
LOGICAL     :: et

DO k = 1, n-1
  amax = x(k)%qr
  imax = k
  DO i = k+1,n
    rxi = x(i)%qr
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

END PROGRAM poly20

